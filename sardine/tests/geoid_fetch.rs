/// Integration test: verify that EGM96 geoid data loads correctly.
///
/// This test is gated on:
/// - `--features geoid-fetch` (requires the `ureq` HTTP dependency)
/// - the cached `us_nga_egm96_15.tif` being present in `~/.sardine/geoid/`
///   OR network access to download it
///
/// Run with:
///   cargo test --features geoid-fetch --test geoid_fetch
#[cfg(feature = "geoid-fetch")]
mod tests {
    use sardine::geoid_fetch::fetch_egm96;
    use sardine::geoid::GeoidModel;

    /// Verify the full fetch-and-convert path produces a usable geoid model.
    ///
    /// Checks:
    /// 1. `fetch_egm96()` succeeds and returns an `Egm96Grid`.
    /// 2. The `egm96_2p5deg.bin` cache file is written to disk.
    /// 3. The geoid model returns physically plausible undulation values at two
    ///    reference points (Central Europe and the equatorial Pacific).
    /// 4. A second call to `fetch_egm96()` hits the fast path (bin already exists).
    #[test]
    fn egm96_round_trip() {
        let grid = fetch_egm96()
            .expect("fetch_egm96 should succeed (cached tif present or network available)");

        // The .bin file must now exist.
        let home = std::env::var("HOME").expect("HOME must be set");
        let bin_path = std::path::PathBuf::from(&home)
            .join(".sardine")
            .join("geoid")
            .join("egm96_2p5deg.bin");
        assert!(
            bin_path.exists(),
            "egm96_2p5deg.bin was not written to {}",
            bin_path.display()
        );

        // The bin file must be the expected size: 73 rows × 144 cols × 4 bytes.
        let meta = std::fs::metadata(&bin_path).expect("bin file metadata");
        assert_eq!(
            meta.len(),
            73 * 144 * 4,
            "egm96_2p5deg.bin has unexpected size"
        );

        // Wrap in GeoidModel and spot-check undulation values.
        let model = GeoidModel::Egm96(grid);

        // Central Europe (Frankfurt ~48°N, ~8°E): EGM96 undulation ≈ +47 m
        let u_europe = model.undulation_m(48.0, 8.0);
        assert!(
            u_europe > 40.0 && u_europe < 55.0,
            "Europe undulation out of plausible range: {u_europe:.1} m"
        );

        // Equatorial Pacific (0°N, 180°E): EGM96 undulation varies; roughly −30 to +30 m
        let u_pacific = model.undulation_m(0.0, 180.0);
        assert!(
            u_pacific > -50.0 && u_pacific < 50.0,
            "Pacific undulation out of plausible range: {u_pacific:.1} m"
        );

        // Fast path: second call should succeed immediately (bin already written).
        let grid2 = fetch_egm96().expect("second fetch_egm96 should hit fast path");
        let model2 = GeoidModel::Egm96(grid2);
        let u2 = model2.undulation_m(48.0, 8.0);
        assert!(
            (u_europe - u2).abs() < 0.5,
            "fast-path result differs from first call: {u_europe:.2} vs {u2:.2}"
        );
    }
}
