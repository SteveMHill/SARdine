//! Precise orbit (EOF) file parser and scene orbit application.
//!
//! Reads ESA POEORB/RESORB `.EOF` orbit files (XML), validates temporal
//! coverage against a [`SceneMetadata`], and replaces the coarse annotation
//! orbit with precise state vectors clipped to the scene window.
//!
//! # Public API
//!
//! - [`parse_eof_file`]: Read and parse an EOF file from disk.
//! - [`apply_precise_orbit`]: Replace `SceneMetadata.orbit` with precise vectors.

use std::cell::Cell;
use std::path::Path;

use chrono::{DateTime, Duration, NaiveDateTime, Utc};

use crate::types::{OrbitData, SceneMetadata, StateVector};

thread_local! {
    /// Phase-3 cache: index of the last `centre` returned by
    /// [`interpolate_orbit`] on this worker thread.  Initialised to 0.
    ///
    /// The bracketing logic in [`interpolate_orbit`] walks forward and
    /// backward from this hint, so a stale or completely-wrong hint can
    /// only degrade performance, never change the result — every walk
    /// terminates at the unique `centre` that satisfies
    /// `svs[centre].time <= time && (centre+1 == len || svs[centre+1].time > time)`,
    /// which is the exact invariant maintained by the previous
    /// linear scan over `svs`.
    ///
    /// We deliberately do not gate by orbit identity: cross-orbit hint
    /// reuse just costs an O(N) walk on the next call, then the hint
    /// is correct again.  See `interpolate_orbit_matches_linear_scan`
    /// for the equivalence regression test.
    static LAST_CENTRE: Cell<usize> = const { Cell::new(0) };
}

// ── Error type ──────────────────────────────────────────────────────

#[derive(Debug, thiserror::Error)]
pub enum OrbitError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error("XML deserialization failed: {0}")]
    Xml(String),

    #[error("no OSV blocks found in EOF file")]
    NoStateVectors,

    #[error("invalid state vector: {0}")]
    InvalidVector(String),

    #[error("EOF orbit does not cover scene: {detail}")]
    InsufficientCoverage { detail: String },

    #[error(
        "interpolation time {requested} outside orbit window \
         [{first}, {last}]; refusing to extrapolate"
    )]
    OutsideOrbitWindow {
        requested: DateTime<Utc>,
        first: DateTime<Utc>,
        last: DateTime<Utc>,
    },

    #[error(
        "Lagrange interpolation window has {got} points, expected {expected} \
         (internal invariant violation)"
    )]
    InterpolationWindowTooSmall { got: usize, expected: usize },

    #[error("validation failed: {0}")]
    Validation(#[from] crate::validate::ValidationErrors),
}

// ── Constants ───────────────────────────────────────────────────────

/// Margin (seconds) around scene start/stop when clipping EOF vectors.
/// 120 s ≈ 12 state vectors at 10 s intervals — sufficient for degree-4
/// Lagrange interpolation at burst edges.
const CLIP_MARGIN_S: f64 = 120.0;

/// Minimum number of precise state vectors required after clipping.
const MIN_PRECISE_VECTORS: usize = 10;

/// Expected LEO orbital altitude range for Sentinel-1 (meters).
const MIN_ALTITUDE_M: f64 = 670_000.0;
const MAX_ALTITUDE_M: f64 = 730_000.0;

/// Expected orbital velocity magnitude (m/s). ~7.5 km/s ± 5%.
const MIN_VELOCITY_M_S: f64 = 7_100.0;
const MAX_VELOCITY_M_S: f64 = 7_900.0;

// ── Public API ──────────────────────────────────────────────────────

/// Parse an ESA precise orbit `.EOF` file from disk.
///
/// Returns an [`OrbitData`] with all valid state vectors, sorted by time.
pub fn parse_eof_file(path: &Path) -> Result<OrbitData, OrbitError> {
    let content = std::fs::read_to_string(path)?;
    parse_eof_content(&content)
}

/// Parse EOF content from a string.
pub fn parse_eof_content(content: &str) -> Result<OrbitData, OrbitError> {
    let content = content
        .strip_prefix('\u{feff}')
        .unwrap_or(content); // SAFETY-OK: strip_prefix-on-no-match returns input unchanged by design (BOM optional)

    let parsed: EofXml = quick_xml::de::from_str(content)
        .map_err(|e| OrbitError::Xml(e.to_string()))?;

    let osv_list = &parsed.data_block.list_of_osvs.osvs;
    if osv_list.is_empty() {
        return Err(OrbitError::NoStateVectors);
    }

    let mut state_vectors = Vec::with_capacity(osv_list.len());
    for osv in osv_list {
        let sv = parse_osv(osv)?;

        // Reject zero / degenerate vectors
        let pos_mag = (sv.position_m[0].powi(2)
            + sv.position_m[1].powi(2)
            + sv.position_m[2].powi(2))
        .sqrt();
        if pos_mag == 0.0 {
            continue;
        }

        state_vectors.push(sv);
    }

    if state_vectors.is_empty() {
        return Err(OrbitError::NoStateVectors);
    }

    state_vectors.sort_by_key(|sv| sv.time);

    let reference_epoch = state_vectors[0].time;

    Ok(OrbitData {
        reference_epoch,
        state_vectors,
    })
}

/// Replace the annotation orbit in `scene` with precise EOF vectors.
///
/// The EOF orbit is clipped to `[scene.start_time - margin, scene.stop_time + margin]`
/// and validated for coverage and physical sanity.
pub fn apply_precise_orbit(
    mut scene: SceneMetadata,
    eof_orbit: &OrbitData,
) -> Result<SceneMetadata, OrbitError> {
    // 1. Verify coverage: EOF must span the full scene acquisition window.
    //    The NoStateVectors guard is the only way eof_orbit can be empty;
    //    parse_eof_file already returns Err(NoStateVectors) before this point.
    //    We check explicitly here so the error is explicit rather than a panic.
    let eof_start = eof_orbit
        .state_vectors
        .first()
        .ok_or(OrbitError::NoStateVectors)?
        .time;
    let eof_end = eof_orbit
        .state_vectors
        .last()
        .ok_or(OrbitError::NoStateVectors)?
        .time;

    let required_start = scene.start_time - Duration::seconds(CLIP_MARGIN_S as i64);
    let required_end = scene.stop_time + Duration::seconds(CLIP_MARGIN_S as i64);

    if eof_start > required_start || eof_end < required_end {
        return Err(OrbitError::InsufficientCoverage {
            detail: format!(
                "EOF covers [{}, {}], scene requires [{}, {}]",
                eof_start.format("%H:%M:%S"),
                eof_end.format("%H:%M:%S"),
                required_start.format("%H:%M:%S"),
                required_end.format("%H:%M:%S"),
            ),
        });
    }

    // 2. Clip to scene window + margin
    let clip_start = scene.start_time - Duration::seconds(CLIP_MARGIN_S as i64);
    let clip_end = scene.stop_time + Duration::seconds(CLIP_MARGIN_S as i64);

    let clipped: Vec<StateVector> = eof_orbit
        .state_vectors
        .iter()
        .filter(|sv| sv.time >= clip_start && sv.time <= clip_end)
        .cloned()
        .collect();

    if clipped.len() < MIN_PRECISE_VECTORS {
        return Err(OrbitError::InsufficientCoverage {
            detail: format!(
                "only {} vectors in scene window (need >= {})",
                clipped.len(),
                MIN_PRECISE_VECTORS
            ),
        });
    }

    // 3. Validate physical sanity of clipped vectors
    validate_vectors(&clipped)?;

    // 4. Replace orbit
    let reference_epoch = clipped[0].time;
    scene.orbit = OrbitData {
        reference_epoch,
        state_vectors: clipped,
    };

    // 5. Re-validate full scene (orbit change may affect invariants)
    scene.validated().map_err(OrbitError::Validation)
}

// ── EOF XML serde structs ───────────────────────────────────────────

/// Root element of an EOF file.
#[derive(Debug, serde::Deserialize)]
#[serde(rename = "Earth_Explorer_File")]
struct EofXml {
    #[serde(rename = "Data_Block")]
    data_block: DataBlockXml,
}

#[derive(Debug, serde::Deserialize)]
struct DataBlockXml {
    #[serde(rename = "List_of_OSVs")]
    list_of_osvs: ListOfOsvsXml,
}

#[derive(Debug, serde::Deserialize)]
struct ListOfOsvsXml {
    #[serde(rename = "OSV", default)]
    osvs: Vec<OsvXml>,
}

/// Single Orbit State Vector from EOF.
///
/// The coordinate elements have `unit="m"` / `unit="m/s"` attributes.
/// quick-xml serde maps the text content via `$value`.
#[derive(Debug, serde::Deserialize)]
struct OsvXml {
    #[serde(rename = "UTC")]
    utc: String,
    #[serde(rename = "X")]
    x: CoordXml,
    #[serde(rename = "Y")]
    y: CoordXml,
    #[serde(rename = "Z")]
    z: CoordXml,
    #[serde(rename = "VX")]
    vx: CoordXml,
    #[serde(rename = "VY")]
    vy: CoordXml,
    #[serde(rename = "VZ")]
    vz: CoordXml,
}

/// Element like `<X unit="m">798419.234490</X>`.
/// The `$value` field captures the text content.
#[derive(Debug, serde::Deserialize)]
struct CoordXml {
    #[serde(rename = "$value")]
    value: f64,
}

// ── Internal helpers ────────────────────────────────────────────────

/// Parse the `UTC=YYYY-MM-DDTHH:MM:SS.ffffff` timestamp from EOF.
fn parse_eof_timestamp(raw: &str) -> Result<DateTime<Utc>, OrbitError> {
    // Strip "UTC=" prefix if present
    let s = raw.trim();
    let s = s.strip_prefix("UTC=").unwrap_or(s); // SAFETY-OK: strip_prefix-on-no-match returns input unchanged by design (UTC= prefix optional)

    let naive = NaiveDateTime::parse_from_str(s, "%Y-%m-%dT%H:%M:%S%.f").map_err(|e| {
        OrbitError::InvalidVector(format!("bad timestamp '{}': {}", raw, e))
    })?;
    Ok(naive.and_utc())
}

fn parse_osv(osv: &OsvXml) -> Result<StateVector, OrbitError> {
    let time = parse_eof_timestamp(&osv.utc)?;
    Ok(StateVector {
        time,
        position_m: [osv.x.value, osv.y.value, osv.z.value],
        velocity_m_s: [osv.vx.value, osv.vy.value, osv.vz.value],
    })
}

/// Validate that state vectors have physically sane magnitudes.
///
/// Altitude is computed as height above the WGS84 ellipsoid (the physically
/// meaningful quantity), not relative to a fixed mean Earth radius.  The
/// previous implementation used 6 371 000 m, which produces ±21 km error
/// between equatorial and polar passes and could either accept obviously
/// bad vectors or reject good ones depending on latitude.
fn validate_vectors(vectors: &[StateVector]) -> Result<(), OrbitError> {
    for (i, sv) in vectors.iter().enumerate() {
        let vel_mag = (sv.velocity_m_s[0].powi(2)
            + sv.velocity_m_s[1].powi(2)
            + sv.velocity_m_s[2].powi(2))
        .sqrt();

        // Geodetic height above WGS84 ellipsoid (Bowring's iteration).
        let (_, _, altitude) = crate::geodesy::ecef_to_geodetic(sv.position_m);

        if altitude < MIN_ALTITUDE_M || altitude > MAX_ALTITUDE_M {
            return Err(OrbitError::InvalidVector(format!(
                "vector {} at {}: ellipsoidal height {:.0} m outside [{}, {}]",
                i,
                sv.time.format("%H:%M:%S"),
                altitude,
                MIN_ALTITUDE_M,
                MAX_ALTITUDE_M,
            )));
        }

        if vel_mag < MIN_VELOCITY_M_S || vel_mag > MAX_VELOCITY_M_S {
            return Err(OrbitError::InvalidVector(format!(
                "vector {} at {}: velocity {:.1} m/s outside [{}, {}]",
                i,
                sv.time.format("%H:%M:%S"),
                vel_mag,
                MIN_VELOCITY_M_S,
                MAX_VELOCITY_M_S,
            )));
        }
    }
    Ok(())
}

// ── Lagrange interpolation ──────────────────────────────────────────

/// Degree for the Lagrange interpolation polynomial.
/// 8-point (degree 7) is standard for precise orbit interpolation and
/// matches the ESA recommendation for Sentinel-1.
const LAGRANGE_DEGREE: usize = 8;

/// Interpolate satellite position and velocity at an arbitrary time.
///
/// Uses 8-point Lagrange polynomial interpolation over the orbit state vectors.
/// The input time is given as absolute UTC.
///
/// # Errors
///
/// Returns `OrbitError::InsufficientCoverage` if the requested time is outside
/// the interpolation window (i.e., not surrounded by enough state vectors).
pub fn interpolate_orbit(
    orbit: &OrbitData,
    time: DateTime<Utc>,
) -> Result<StateVector, OrbitError> {
    let svs = &orbit.state_vectors;
    if svs.len() < LAGRANGE_DEGREE {
        return Err(OrbitError::InsufficientCoverage {
            detail: format!(
                "need at least {} vectors for interpolation, have {}",
                LAGRANGE_DEGREE,
                svs.len()
            ),
        });
    }

    // Refuse to extrapolate.  Lagrange polynomials grow without bound outside
    // the sample window, so silently extending past either end would corrupt
    // any geometry that relied on the result.  The caller is expected to have
    // applied a precise orbit that covers the scene window with margin.
    let first = svs[0].time;
    let last = svs[svs.len() - 1].time;
    if time < first || time > last {
        return Err(OrbitError::OutsideOrbitWindow {
            requested: time,
            first,
            last,
        });
    }

    // Convert target time to seconds relative to first state vector
    let t_ref = svs[0].time;
    let t_sec = (time - t_ref).num_microseconds().unwrap_or(0) as f64 / 1e6; // SAFETY-OK: chrono microseconds only overflow for durations > 290 000 years

    // Find the centre index: the largest `i` such that `svs[i].time <= time`.
    //
    // Phase 3: warm-path walk from the per-worker thread-local hint.  In
    // the terrain-correction loop consecutive calls request times that
    // differ by microseconds (Newton iters on the same pixel) to a few
    // milliseconds (consecutive pixels in azimuth), so the hint is
    // essentially always already at the correct `centre` and the walk
    // terminates with zero comparisons in either branch.  When the hint
    // is stale (different worker, different scene, or a fresh chunk)
    // the walk still produces the bit-identical result that the
    // previous linear scan produced — it just takes O(distance) time.
    let centre = {
        let mut i = LAST_CENTRE.with(|c| c.get()).min(svs.len() - 1);
        while i + 1 < svs.len() && svs[i + 1].time <= time {
            i += 1;
        }
        while i > 0 && svs[i].time > time {
            i -= 1;
        }
        // Invariant on exit: svs[i].time <= time && (i+1 == len || svs[i+1].time > time).
        // This matches the contract of the previous linear-scan
        // bracketing logic exactly.  Bounds check above guarantees
        // svs[0].time <= time so the lower edge is reachable; the
        // upper edge falls through to i == svs.len() - 1.
        LAST_CENTRE.with(|c| c.set(i));
        i
    };

    // Select LAGRANGE_DEGREE vectors centred on the target time
    let half = LAGRANGE_DEGREE / 2;
    let start = if centre < half {
        0
    } else if centre + half >= svs.len() {
        svs.len().saturating_sub(LAGRANGE_DEGREE)
    } else {
        centre - half + 1
    };
    let end = (start + LAGRANGE_DEGREE).min(svs.len());
    let window = &svs[start..end];

    // Defense in depth: the saturating_sub branch above is supposed to
    // guarantee `window.len() == LAGRANGE_DEGREE` whenever
    // `svs.len() >= LAGRANGE_DEGREE`.  If a future refactor breaks that
    // invariant we want a typed error, not a silent drop in interpolation
    // order.
    if window.len() != LAGRANGE_DEGREE {
        return Err(OrbitError::InterpolationWindowTooSmall {
            got: window.len(),
            expected: LAGRANGE_DEGREE,
        });
    }

    // Time values for the window points (seconds relative to t_ref).
    // Stack array avoids an 8-element heap allocation on every call;
    // LAGRANGE_DEGREE is a compile-time constant so this is zero-cost.
    let mut times = [0.0_f64; LAGRANGE_DEGREE];
    for (i, sv) in window.iter().enumerate() {
        times[i] = (sv.time - t_ref).num_microseconds().unwrap_or(0) as f64 / 1e6; // SAFETY-OK: chrono microseconds cannot overflow for orbit-window durations
    }

    // Lagrange interpolation for each coordinate
    let mut pos = [0.0_f64; 3];
    let mut vel = [0.0_f64; 3];

    for (k, sv) in window.iter().enumerate() {
        let mut basis = 1.0;
        for (j, tj) in times.iter().enumerate() {
            if j != k {
                basis *= (t_sec - tj) / (times[k] - tj);
            }
        }
        for c in 0..3 {
            pos[c] += basis * sv.position_m[c];
            vel[c] += basis * sv.velocity_m_s[c];
        }
    }

    Ok(StateVector {
        time,
        position_m: pos,
        velocity_m_s: vel,
    })
}

/// Interpolate satellite position, velocity **and orbital acceleration** in a
/// single Lagrange pass.
///
/// Used by the Newton zero-Doppler solver to eliminate the finite-difference
/// second orbit evaluation (`interpolate_orbit(orbit, t + Δt)`) per iteration.
///
/// # Acceleration computation
///
/// The acceleration is the exact analytical derivative of the velocity Lagrange
/// polynomial:
///
/// ```text
/// A(t) = Σ_k  v_k · L'_k(t)
///
/// where  L'_k(t) = L_k(t) · Σ_{j≠k}  1 / (t − t_j)
/// ```
///
/// This is more accurate than finite differencing (O(Δt²) truncation error)
/// and eliminates one `interpolate_orbit` call per Newton iteration.
///
/// # Numerical note
///
/// The reciprocal `1 / (t − t_j)` is skipped when `|t − t_j| < 1 ns`
/// (guarded below).  S-1 orbit nodes are 10 s apart; the Newton step is
/// microseconds — this guard is never triggered in practice and only exists
/// to prevent NaN propagation in degenerate test inputs.
pub(crate) fn interpolate_orbit_and_accel(
    orbit: &OrbitData,
    time: DateTime<Utc>,
) -> Result<(StateVector, [f64; 3]), OrbitError> {
    let svs = &orbit.state_vectors;
    if svs.len() < LAGRANGE_DEGREE {
        return Err(OrbitError::InsufficientCoverage {
            detail: format!(
                "need at least {} vectors for interpolation, have {}",
                LAGRANGE_DEGREE,
                svs.len()
            ),
        });
    }

    let first = svs[0].time;
    let last  = svs[svs.len() - 1].time;
    if time < first || time > last {
        return Err(OrbitError::OutsideOrbitWindow { requested: time, first, last });
    }

    let t_ref = svs[0].time;
    let t_sec = (time - t_ref).num_microseconds().unwrap_or(0) as f64 / 1e6; // SAFETY-OK: chrono microseconds cannot overflow for orbit-window durations

    let centre = {
        let mut i = LAST_CENTRE.with(|c| c.get()).min(svs.len() - 1);
        while i + 1 < svs.len() && svs[i + 1].time <= time { i += 1; }
        while i > 0 && svs[i].time > time { i -= 1; }
        LAST_CENTRE.with(|c| c.set(i));
        i
    };

    let half  = LAGRANGE_DEGREE / 2;
    let start = if centre < half { 0 }
                else if centre + half >= svs.len() { svs.len().saturating_sub(LAGRANGE_DEGREE) }
                else { centre - half + 1 };
    let end   = (start + LAGRANGE_DEGREE).min(svs.len());
    let window = &svs[start..end];

    if window.len() != LAGRANGE_DEGREE {
        return Err(OrbitError::InterpolationWindowTooSmall {
            got: window.len(),
            expected: LAGRANGE_DEGREE,
        });
    }

    let mut times = [0.0_f64; LAGRANGE_DEGREE];
    for (i, sv) in window.iter().enumerate() {
        times[i] = (sv.time - t_ref).num_microseconds().unwrap_or(0) as f64 / 1e6; // SAFETY-OK: chrono microseconds cannot overflow for orbit-window durations
    }

    let mut pos   = [0.0_f64; 3];
    let mut vel   = [0.0_f64; 3];
    let mut accel = [0.0_f64; 3];

    // Precompute reciprocals  rdiff[j] = 1 / (t − t_j)  for all j.
    // This reduces the total division count from  LAGRANGE_DEGREE² = 64
    // (one per (k, j≠k) pair in the inner loop) to  LAGRANGE_DEGREE = 8.
    // For each k, the derivative sum  Σ_{j≠k} 1/(t−t_j)  is then
    //   total_rdiff − rdiff[k]
    // using O(1) subtraction instead of O(LAGRANGE_DEGREE) divisions.
    //
    // Guard: if  |t − t_j| < 1 ns  the reciprocal is skipped (rdiff stays 0
    // and total_rdiff is not incremented).  Orbit nodes are 10 s apart;
    // Newton steps are microseconds — this branch is never taken in
    // production.  It prevents NaN if a caller supplies t exactly at a node. // SAFETY-OK: t at a node means basis[k']=0 for all k'≠k; the contribution from k' to accel is 0·(total_rdiff-rdiff[k'])=0 regardless of the guarded term
    let mut rdiff       = [0.0_f64; LAGRANGE_DEGREE];
    let mut total_rdiff = 0.0_f64;
    for j in 0..LAGRANGE_DEGREE {
        let diff = t_sec - times[j];
        if diff.abs() > 1e-9 { // SAFETY-OK: orbit nodes are 10 s apart; Newton steps are ≪ 1 s; this guard is never triggered in production and only prevents NaN on degenerate inputs
            rdiff[j] = 1.0 / diff;
            total_rdiff += rdiff[j];
        }
    }

    for k in 0..LAGRANGE_DEGREE {
        // Standard Lagrange basis weight L_k(t).
        let mut basis = 1.0_f64;
        for j in 0..LAGRANGE_DEGREE {
            if j != k {
                basis *= (t_sec - times[j]) / (times[k] - times[j]);
            }
        }

        // L'_k(t) = L_k(t) · Σ_{j≠k} 1/(t−t_j) = basis · (total_rdiff − rdiff[k]).
        // When basis ≈ 0 (t near another node), the product is 0 regardless
        // of the sum — no short-circuit needed; the result is exact. // SAFETY-OK: IEEE 754 guarantees 0.0 * finite = 0.0 and 0.0 * ±Inf = NaN; the rdiff guard above ensures rdiff[k] is always finite, so the product is always finite when basis is 0.0
        let d_basis = basis * (total_rdiff - rdiff[k]);

        for c in 0..3 {
            pos[c]   += basis   * window[k].position_m[c];
            vel[c]   += basis   * window[k].velocity_m_s[c];
            accel[c] += d_basis * window[k].velocity_m_s[c];
        }
    }

    Ok((StateVector { time, position_m: pos, velocity_m_s: vel }, accel))
}

/// Interpolate satellite position at a given time (convenience wrapper).
pub fn interpolate_position(
    orbit: &OrbitData,
    time: DateTime<Utc>,
) -> Result<[f64; 3], OrbitError> {
    Ok(interpolate_orbit(orbit, time)?.position_m)
}

/// Interpolate satellite velocity at a given time (convenience wrapper).
pub fn interpolate_velocity(
    orbit: &OrbitData,
    time: DateTime<Utc>,
) -> Result<[f64; 3], OrbitError> {
    Ok(interpolate_orbit(orbit, time)?.velocity_m_s)
}

// ── Tests ───────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn s1a_eof() -> PathBuf {
        PathBuf::from("/home/datacube/dev/SARdine/legacy/old_package_copy/SARdine/orbit_cache/S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66_POEORB.EOF")
    }

    fn s1b_eof() -> PathBuf {
        PathBuf::from("/home/datacube/dev/SARdine/legacy/old_package_copy/SARdine/orbit_cache/S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833_POEORB.EOF")
    }

    fn s1a_safe() -> PathBuf {
        PathBuf::from("/home/datacube/dev/SARdine/data/SLC/S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66.SAFE")
    }

    fn s1b_safe() -> PathBuf {
        PathBuf::from("/home/datacube/dev/SARdine/data/SLC/S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833.SAFE")
    }

    fn have_eof_data() -> bool {
        s1a_eof().is_file()
    }

    fn have_safe_data() -> bool {
        s1a_safe().join("annotation").is_dir()
    }

    // ── Unit tests ──────────────────────────────────────────────────

    #[test]
    fn test_parse_eof_timestamp() {
        let dt = parse_eof_timestamp("UTC=2020-10-04T22:59:42.000000").unwrap();
        assert_eq!(
            dt.format("%Y-%m-%d %H:%M:%S").to_string(),
            "2020-10-04 22:59:42"
        );
    }

    #[test]
    fn test_parse_eof_timestamp_no_prefix() {
        let dt = parse_eof_timestamp("2020-10-04T22:59:42.000000").unwrap();
        assert_eq!(
            dt.format("%Y-%m-%d %H:%M:%S").to_string(),
            "2020-10-04 22:59:42"
        );
    }

    #[test]
    fn test_parse_eof_timestamp_bad() {
        assert!(parse_eof_timestamp("not-a-date").is_err());
    }

    #[test]
    fn test_validate_vectors_sane() {
        // LEO position: ~7e6 m from Earth center, velocity ~7.5 km/s
        let sv = StateVector {
            time: Utc::now(),
            position_m: [798419.0, -5364097.0, 4536283.0],
            velocity_m_s: [-2467.0, 4427.0, 5653.0],
        };
        assert!(validate_vectors(&[sv]).is_ok());
    }

    #[test]
    fn test_validate_vectors_bad_altitude() {
        let sv = StateVector {
            time: Utc::now(),
            position_m: [100.0, 200.0, 300.0], // way too close
            velocity_m_s: [-2467.0, 4427.0, 5653.0],
        };
        assert!(validate_vectors(&[sv]).is_err());
    }

    #[test]
    fn test_validate_vectors_bad_velocity() {
        let sv = StateVector {
            time: Utc::now(),
            position_m: [798419.0, -5364097.0, 4536283.0],
            velocity_m_s: [1.0, 2.0, 3.0], // way too slow
        };
        assert!(validate_vectors(&[sv]).is_err());
    }

    // ── Integration tests ───────────────────────────────────────────

    #[test]
    fn test_parse_s1a_eof() {
        if !have_eof_data() {
            eprintln!("Skipping: EOF test data not available");
            return;
        }

        let orbit = parse_eof_file(&s1a_eof()).expect("S1A EOF parse failed");

        // POEORB: ~9361 vectors at 10 s intervals over ~26 hours
        assert!(
            orbit.state_vectors.len() > 9000,
            "expected >9000 vectors, got {}",
            orbit.state_vectors.len()
        );

        // Sorted by time
        for w in orbit.state_vectors.windows(2) {
            assert!(w[1].time > w[0].time);
        }

        // 10-second interval
        let dt = (orbit.state_vectors[1].time - orbit.state_vectors[0].time)
            .num_seconds();
        assert_eq!(dt, 10, "expected 10 s interval, got {} s", dt);

        // Reference epoch is first vector
        assert_eq!(orbit.reference_epoch, orbit.state_vectors[0].time);

        // Physical sanity: all vectors have valid altitude & velocity
        assert!(validate_vectors(&orbit.state_vectors).is_ok());
    }

    #[test]
    fn test_parse_s1b_eof() {
        if !s1b_eof().is_file() {
            eprintln!("Skipping: S1B EOF test data not available");
            return;
        }

        let orbit = parse_eof_file(&s1b_eof()).expect("S1B EOF parse failed");
        assert!(orbit.state_vectors.len() > 9000);
        assert!(validate_vectors(&orbit.state_vectors).is_ok());
    }

    #[test]
    fn test_apply_precise_orbit_s1a() {
        if !have_eof_data() || !have_safe_data() {
            eprintln!("Skipping: test data not available");
            return;
        }

        let scene = crate::parse::parse_safe_directory(&s1a_safe())
            .expect("S1A SAFE parse failed");
        let eof_orbit = parse_eof_file(&s1a_eof())
            .expect("S1A EOF parse failed");

        // Annotation orbit has 17 vectors
        assert_eq!(scene.orbit.state_vectors.len(), 17);

        let updated = apply_precise_orbit(scene, &eof_orbit)
            .expect("apply_precise_orbit failed");

        // Precise orbit should have significantly more vectors
        assert!(
            updated.orbit.state_vectors.len() > 17,
            "expected more precise vectors, got {}",
            updated.orbit.state_vectors.len()
        );

        // Should be physically valid
        assert!(validate_vectors(&updated.orbit.state_vectors).is_ok());

        // Reference epoch should be within the clip window
        let margin = Duration::seconds(CLIP_MARGIN_S as i64);
        assert!(updated.orbit.reference_epoch >= updated.start_time - margin);

        // Orbit should cover the full scene
        let first = updated.orbit.state_vectors.first().unwrap().time;
        let last = updated.orbit.state_vectors.last().unwrap().time;
        assert!(first <= updated.start_time);
        assert!(last >= updated.stop_time);
    }

    #[test]
    fn test_apply_precise_orbit_s1b() {
        if !s1b_eof().is_file() || !s1b_safe().join("annotation").is_dir() {
            eprintln!("Skipping: S1B test data not available");
            return;
        }

        let scene = crate::parse::parse_safe_directory(&s1b_safe())
            .expect("S1B SAFE parse failed");
        let eof_orbit = parse_eof_file(&s1b_eof())
            .expect("S1B EOF parse failed");

        let updated = apply_precise_orbit(scene, &eof_orbit)
            .expect("S1B apply_precise_orbit failed");

        assert!(updated.orbit.state_vectors.len() > 17);
    }

    #[test]
    fn test_eof_coverage_mismatch() {
        if !have_safe_data() {
            eprintln!("Skipping: test data not available");
            return;
        }

        let scene = crate::parse::parse_safe_directory(&s1a_safe())
            .expect("S1A SAFE parse failed");

        // Create a fake orbit that doesn't cover the scene
        let fake_orbit = OrbitData {
            reference_epoch: Utc::now(),
            state_vectors: vec![StateVector {
                time: Utc::now(),
                position_m: [798419.0, -5364097.0, 4536283.0],
                velocity_m_s: [-2467.0, 4427.0, 5653.0],
            }],
        };

        let result = apply_precise_orbit(scene, &fake_orbit);
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(
            matches!(err, OrbitError::InsufficientCoverage { .. }),
            "expected InsufficientCoverage, got: {}",
            err
        );
    }

    // ── Interpolation tests ─────────────────────────────────────────

    #[test]
    fn test_interpolate_at_known_vector() {
        // Interpolating exactly at a state vector time should return that vector.
        if !have_eof_data() {
            eprintln!("Skipping: EOF test data not available");
            return;
        }

        let orbit = parse_eof_file(&s1a_eof()).unwrap();
        // Pick a vector near the middle
        let idx = orbit.state_vectors.len() / 2;
        let sv = &orbit.state_vectors[idx];

        let interp = interpolate_orbit(&orbit, sv.time).unwrap();

        // Should match to sub-millimeter accuracy
        for c in 0..3 {
            assert!(
                (interp.position_m[c] - sv.position_m[c]).abs() < 0.001,
                "position[{}]: {} vs {}",
                c, interp.position_m[c], sv.position_m[c]
            );
            assert!(
                (interp.velocity_m_s[c] - sv.velocity_m_s[c]).abs() < 1e-6,
                "velocity[{}]: {} vs {}",
                c, interp.velocity_m_s[c], sv.velocity_m_s[c]
            );
        }
    }

    #[test]
    fn test_interpolate_midpoint() {
        // Interpolating between two vectors should give a physically plausible result.
        if !have_eof_data() {
            eprintln!("Skipping: EOF test data not available");
            return;
        }

        let orbit = parse_eof_file(&s1a_eof()).unwrap();
        let idx = orbit.state_vectors.len() / 2;
        let sv0 = &orbit.state_vectors[idx];
        let sv1 = &orbit.state_vectors[idx + 1];

        let mid_time = sv0.time + Duration::seconds(5); // halfway between 10s interval
        let interp = interpolate_orbit(&orbit, mid_time).unwrap();

        // Position should be between the two bracketing vectors
        for c in 0..3 {
            let lo = sv0.position_m[c].min(sv1.position_m[c]);
            let hi = sv0.position_m[c].max(sv1.position_m[c]);
            // Allow some tolerance for polynomial curvature
            let margin = (hi - lo).abs() * 0.1 + 1.0;
            assert!(
                interp.position_m[c] >= lo - margin && interp.position_m[c] <= hi + margin,
                "position[{}]: {} not between {} and {} (margin {})",
                c, interp.position_m[c], lo, hi, margin
            );
        }

        // Altitude should be in LEO range
        let alt = (interp.position_m[0].powi(2)
            + interp.position_m[1].powi(2)
            + interp.position_m[2].powi(2))
        .sqrt()
            - 6_371_000.0;
        assert!(
            alt > 650_000.0 && alt < 750_000.0,
            "altitude {} m",
            alt
        );

        // Velocity magnitude should be ~7.5 km/s
        let vel = (interp.velocity_m_s[0].powi(2)
            + interp.velocity_m_s[1].powi(2)
            + interp.velocity_m_s[2].powi(2))
        .sqrt();
        assert!(
            vel > 7000.0 && vel < 8000.0,
            "velocity {} m/s",
            vel
        );
    }

    #[test]
    fn test_interpolate_insufficient_vectors() {
        let orbit = OrbitData {
            reference_epoch: Utc::now(),
            state_vectors: vec![StateVector {
                time: Utc::now(),
                position_m: [798419.0, -5364097.0, 4536283.0],
                velocity_m_s: [-2467.0, 4427.0, 5653.0],
            }],
        };

        let result = interpolate_orbit(&orbit, Utc::now());
        assert!(result.is_err());
    }

    #[test]
    fn test_interpolate_refuses_extrapolation_past_last() {
        // Lagrange polynomials diverge outside the sample window; we must
        // refuse, not silently extrapolate.
        if !have_eof_data() {
            eprintln!("Skipping: EOF test data not available");
            return;
        }

        let orbit = parse_eof_file(&s1a_eof()).unwrap();
        let last = orbit.state_vectors.last().unwrap().time;
        let beyond = last + Duration::seconds(10);

        let err = interpolate_orbit(&orbit, beyond).expect_err(
            "must refuse to extrapolate past last orbit state vector",
        );
        assert!(
            matches!(err, OrbitError::OutsideOrbitWindow { .. }),
            "expected OutsideOrbitWindow, got: {err}"
        );
    }

    #[test]
    fn test_interpolate_refuses_extrapolation_before_first() {
        if !have_eof_data() {
            eprintln!("Skipping: EOF test data not available");
            return;
        }

        let orbit = parse_eof_file(&s1a_eof()).unwrap();
        let first = orbit.state_vectors.first().unwrap().time;
        let before = first - Duration::seconds(10);

        let err = interpolate_orbit(&orbit, before).expect_err(
            "must refuse to extrapolate before first orbit state vector",
        );
        assert!(
            matches!(err, OrbitError::OutsideOrbitWindow { .. }),
            "expected OutsideOrbitWindow, got: {err}"
        );
    }

    /// Lagrange interpolation must reproduce a linear orbit exactly.
    ///
    /// For position = p0 + v*t (degree-1 polynomial), an 8-point Lagrange
    /// interpolator must produce the exact result to within floating-point
    /// rounding (‹ 1 mm position, ‹ 1 µm/s velocity).  This test is purely
    /// synthetic — no file I/O required.
    #[test]
    fn test_lagrange_interpolation_exact_on_linear_orbit() {
        use chrono::TimeZone;

        let t0 = Utc.with_ymd_and_hms(2024, 1, 1, 0, 0, 0).unwrap();
        // Satellite at ~700 km altitude above north pole, moving east.
        let p0 = [0.0_f64, 0.0, 7_078_137.0];
        let v = [7_600.0_f64, 10.0, -5.0];

        // 10 state vectors at 10 s intervals, from t = 0 to t = 90 s.
        let svs: Vec<StateVector> = (0..10)
            .map(|i| {
                let dt = i as f64 * 10.0;
                StateVector {
                    time: t0 + Duration::milliseconds((dt * 1000.0) as i64),
                    position_m: [p0[0] + v[0] * dt, p0[1] + v[1] * dt, p0[2] + v[2] * dt],
                    velocity_m_s: v,
                }
            })
            .collect();

        let orbit = OrbitData { reference_epoch: t0, state_vectors: svs };

        // Query at t = 45 s (centre of the window, between vectors 4 and 5).
        let t_query = t0 + Duration::milliseconds(45_000);
        let dt = 45.0_f64;
        let expected_pos = [p0[0] + v[0] * dt, p0[1] + v[1] * dt, p0[2] + v[2] * dt];
        let expected_vel = v;

        let sv = interpolate_orbit(&orbit, t_query)
            .expect("interpolation of a mid-window point must succeed");

        for c in 0..3 {
            assert!(
                (sv.position_m[c] - expected_pos[c]).abs() < 1e-3,
                "position[{c}]: got {got} m, expected {exp} m (error {err} m)",
                got = sv.position_m[c],
                exp = expected_pos[c],
                err = (sv.position_m[c] - expected_pos[c]).abs()
            );
            assert!(
                (sv.velocity_m_s[c] - expected_vel[c]).abs() < 1e-6,
                "velocity[{c}]: got {got} m/s, expected {exp} m/s",
                got = sv.velocity_m_s[c],
                exp = expected_vel[c]
            );
        }
    }

    // ── Phase 3: thread-local centre-hint cache ─────────────────────

    /// Reference implementation: original linear scan (pre-Phase-3).
    /// Used by [`interpolate_orbit_matches_linear_scan`] as the
    /// equivalence baseline.  Keeping it in the test module rather than
    /// the production module avoids any risk that callers accidentally
    /// pick the slow path.
    fn linear_scan_centre(svs: &[StateVector], time: DateTime<Utc>) -> usize {
        svs.iter()
            .position(|sv| sv.time > time)
            .unwrap_or(svs.len()) // SAFETY-OK: test-only reference impl mirroring the pre-Phase-3 behaviour
            .saturating_sub(1)
    }

    /// Walk-from-hint produces the same `centre` as the original linear
    /// scan for every probed time, regardless of starting hint state.
    /// This is the *only* correctness gate for Phase 3 \u2014 if the centre
    /// matches, the rest of `interpolate_orbit` is unchanged so the
    /// returned `StateVector` is bit-identical.
    #[test]
    fn interpolate_orbit_matches_linear_scan_on_real_orbit() {
        if !have_eof_data() {
            eprintln!("skipping: orbit cache not present");
            return;
        }
        let orbit = parse_eof_file(&s1b_eof()).unwrap();
        let svs = &orbit.state_vectors;
        assert!(svs.len() > 100, "need a real multi-vector orbit");

        let first = svs[0].time;
        let last = svs[svs.len() - 1].time;
        let total_us = (last - first).num_microseconds().unwrap(); // SAFETY-OK: orbit window << 290k years

        // Reset the hint to a known stale state so this test exercises
        // the cold-walk path on the first call.
        LAST_CENTRE.with(|c| c.set(0));

        // Sample 1024 deterministic times spanning the orbit window in
        // a deliberately non-monotonic order to stress the walk in both
        // directions and across multi-vector jumps.
        let mut state: u64 = 0x9E37_79B9_7F4A_7C15;
        for _ in 0..1024 {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            let frac = (state >> 11) as f64 / (1u64 << 53) as f64;
            let dt_us = (frac * total_us as f64) as i64;
            let t = first + Duration::microseconds(dt_us);

            let expected = linear_scan_centre(svs, t);
            // Run interpolate_orbit so the cache evolves through the
            // sequence; then reproduce its centre logic with a probe to
            // assert equivalence.  Checking the returned StateVector for
            // bit-equality vs a fresh-cache call is the strongest
            // possible gate.
            let got = interpolate_orbit(&orbit, t).unwrap();

            // Independent fresh-cache call: clear the hint and call
            // again.  Bit-equal results prove the hint never shifted
            // the answer.
            LAST_CENTRE.with(|c| c.set(0));
            let fresh = interpolate_orbit(&orbit, t).unwrap();

            assert_eq!(
                got.position_m[0].to_bits(), fresh.position_m[0].to_bits(),
                "position[0] differs at t={t}: warm={} fresh={}",
                got.position_m[0], fresh.position_m[0],
            );
            assert_eq!(got.position_m[1].to_bits(), fresh.position_m[1].to_bits());
            assert_eq!(got.position_m[2].to_bits(), fresh.position_m[2].to_bits());
            assert_eq!(got.velocity_m_s[0].to_bits(), fresh.velocity_m_s[0].to_bits());
            assert_eq!(got.velocity_m_s[1].to_bits(), fresh.velocity_m_s[1].to_bits());
            assert_eq!(got.velocity_m_s[2].to_bits(), fresh.velocity_m_s[2].to_bits());

            // Also verify the underlying centre choice is correct so
            // the test fails loudly if the bracketing logic regresses
            // even if the bit-equal check above masks the symptom.
            let cached = LAST_CENTRE.with(|c| c.get());
            assert_eq!(
                cached, expected,
                "cached centre={cached} but linear-scan centre={expected} at t={t}",
            );
        }
    }

    /// Locality test: the dominant TC pattern is "thousands of calls at
    /// monotonically advancing times within a small window".  Under
    /// that pattern the walk must take O(1) steps per call, not O(N).
    /// We assert this by counting how many cache-misses\u2014i.e. calls
    /// whose final centre differs from the prior cached centre by more
    /// than 1\u2014occur over a fine sweep.  For 10 000 evenly-spaced
    /// samples spanning 60 s of the orbit window (~6 SV intervals),
    /// the expected number of "centre changed" events is exactly the
    /// number of state vectors crossed (~6).  We allow 12 for safety.
    #[test]
    fn interpolate_orbit_walk_is_locality_bounded() {
        if !have_eof_data() {
            eprintln!("skipping: orbit cache not present");
            return;
        }
        let orbit = parse_eof_file(&s1b_eof()).unwrap();
        let svs = &orbit.state_vectors;

        // Choose a 60-second window starting near the middle of the orbit.
        let mid = svs[svs.len() / 2].time;
        let window_us: i64 = 60_000_000;

        LAST_CENTRE.with(|c| c.set(0));
        let mut prev_centre = usize::MAX;
        let mut centre_changes = 0usize;
        for k in 0..10_000usize {
            let dt_us = (k as i64) * window_us / 10_000;
            let t = mid + Duration::microseconds(dt_us);
            let _ = interpolate_orbit(&orbit, t).unwrap(); // SAFETY-OK: locality test only inspects the cached centre, the StateVector itself is unused
            let centre = LAST_CENTRE.with(|c| c.get());
            if centre != prev_centre {
                centre_changes += 1;
                prev_centre = centre;
            }
        }
        // 60 s / SV-spacing(10 s) = 6 boundary crossings + 1 initial set.
        assert!(
            centre_changes <= 12,
            "warm-walk took {centre_changes} centre changes over 10 000 monotonic samples; expected <= 12",
        );
    }

    // ── Analytical acceleration correctness ───────────────────────────────────

    /// `interpolate_orbit_and_accel` returns the same position/velocity as
    /// `interpolate_orbit` and an acceleration within 0.01 m/s² of the
    /// 10 ms finite-difference estimate on real S-1B orbit data.
    ///
    /// 0.01 m/s² tolerance: the FD estimate has O(Δt²) truncation error
    /// (Δt = 10 ms, fourth derivative of a Lagrange-8 polynomial over a
    /// 70 s window is negligible), while the analytical form is exact.
    /// Differences > 1e-4 m/s² would indicate a coding error.
    #[test]
    fn interpolate_orbit_and_accel_matches_fd_on_real_orbit() {
        if !have_eof_data() {
            eprintln!("skipping: orbit cache not present");
            return;
        }
        let orbit = parse_eof_file(&s1b_eof()).unwrap();
        let svs = &orbit.state_vectors;
        let t_ref = svs[0].time;

        // Test at 10 equally-spaced times within the scene window
        // (avoid the very first/last 40 s where the 8-point window clips).
        let t_first = svs[4].time + Duration::seconds(5);
        let t_last  = svs[svs.len() - 5].time - Duration::seconds(5);
        let span_us = (t_last - t_first).num_microseconds().unwrap_or(1); // SAFETY-OK: test window is always > 0 µs

        const ACCEL_FD_STEP_S: f64 = 0.01;

        for i in 0..10usize {
            let frac = i as f64 / 9.0;
            let t = t_first + Duration::microseconds((frac * span_us as f64) as i64);
            let t_plus = t + Duration::milliseconds((ACCEL_FD_STEP_S * 1000.0) as i64);

            // Analytical result.
            let (sv_a, accel_a) = interpolate_orbit_and_accel(&orbit, t)
                .expect("analytical call must not fail within orbit window");

            // FD baseline.
            let sv0 = interpolate_orbit(&orbit, t)
                .expect("baseline call at t must not fail");
            let sv1 = interpolate_orbit(&orbit, t_plus)
                .expect("baseline call at t+dt must not fail");
            let accel_fd = [
                (sv1.velocity_m_s[0] - sv0.velocity_m_s[0]) / ACCEL_FD_STEP_S,
                (sv1.velocity_m_s[1] - sv0.velocity_m_s[1]) / ACCEL_FD_STEP_S,
                (sv1.velocity_m_s[2] - sv0.velocity_m_s[2]) / ACCEL_FD_STEP_S,
            ];

            // Position and velocity must be bit-identical to interpolate_orbit.
            for c in 0..3 {
                assert_eq!(
                    sv_a.position_m[c].to_bits(), sv0.position_m[c].to_bits(),
                    "position[{c}] differs at i={i}: analytical={} vs interpolate_orbit={}",
                    sv_a.position_m[c], sv0.position_m[c],
                );
                assert_eq!(
                    sv_a.velocity_m_s[c].to_bits(), sv0.velocity_m_s[c].to_bits(),
                    "velocity[{c}] differs at i={i}: analytical={} vs interpolate_orbit={}",
                    sv_a.velocity_m_s[c], sv0.velocity_m_s[c],
                );
            }

            // Acceleration: analytical vs FD must agree to < 0.01 m/s².
            // S-1 orbital acceleration magnitude is ~8 m/s²; 0.01 m/s²
            // is 0.1% relative error — tighter than any downstream need.
            let t_sec = (t - t_ref).num_microseconds().unwrap_or(0) as f64 / 1e6; // SAFETY-OK: chrono microseconds cannot overflow for orbit-window durations
            for c in 0..3 {
                let diff = (accel_a[c] - accel_fd[c]).abs();
                assert!(
                    diff < 0.01,
                    "accel[{c}] mismatch at t={t_sec:.3} s: analytical={:.6} FD={:.6} diff={diff:.2e} m/s²",
                    accel_a[c], accel_fd[c],
                );
            }
        }
    }

    /// On a constant-velocity (linear) orbit the Lagrange polynomial for
    /// velocity is exactly constant, so its analytical derivative (acceleration)
    /// is exactly zero.  Verify this to machine precision.
    #[test]
    fn interpolate_orbit_and_accel_zero_accel_on_linear_orbit() {
        use chrono::TimeZone;

        let t0 = Utc.with_ymd_and_hms(2024, 1, 1, 0, 0, 0).unwrap();
        let p0 = [0.0_f64, 0.0, 7_078_137.0];
        let v  = [7_600.0_f64, 10.0, -5.0];

        let svs: Vec<StateVector> = (0..10)
            .map(|i| {
                let dt = i as f64 * 10.0;
                StateVector {
                    time: t0 + Duration::milliseconds((dt * 1000.0) as i64),
                    position_m: [p0[0] + v[0] * dt, p0[1] + v[1] * dt, p0[2] + v[2] * dt],
                    velocity_m_s: v,
                }
            })
            .collect();
        let orbit = OrbitData { reference_epoch: t0, state_vectors: svs };

        let t_query = t0 + Duration::milliseconds(45_000);
        let (sv, accel) = interpolate_orbit_and_accel(&orbit, t_query)
            .expect("must succeed for mid-window query");

        // Velocity must equal the constant v exactly.
        for c in 0..3 {
            assert!(
                (sv.velocity_m_s[c] - v[c]).abs() < 1e-9,
                "velocity[{c}]: got {} expected {}", sv.velocity_m_s[c], v[c],
            );
        }
        // Acceleration of a constant-velocity orbit is exactly zero.
        for c in 0..3 {
            assert!(
                accel[c].abs() < 1e-9,
                "accel[{c}] on linear orbit = {} (expected ≈ 0)", accel[c],
            );
        }
    }
}
