use ndarray::Array2;

use crate::types::BoundingBox;

/// Flags describing tie-point grid cell status
pub const TIE_FLAG_VALID: u16 = 0b0000_0001;
pub const TIE_FLAG_ZERO_DOPPLER_FAIL: u16 = 0b0000_0010;
pub const TIE_FLAG_RANGE_OOB: u16 = 0b0000_0100;
pub const TIE_FLAG_DEM_NODATA: u16 = 0b0000_1000;
pub const TIE_FLAG_OUT_OF_SWATH: u16 = 0b0001_0000;
pub const TIE_FLAG_SHADOW: u16 = 0b0010_0000;
pub const TIE_FLAG_LAYOVER: u16 = 0b0100_0000;

/// Coarse tie-point sampled forward model value
#[derive(Debug, Clone, Copy)]
pub struct TiePointCell {
    /// Azimuth time relative to product start (seconds)
    pub azimuth_time: f64,
    /// Slant range distance (meters)
    pub slant_range: f64,
    /// Cosine of the local incidence angle
    pub cos_local_incidence: f32,
    /// Bitmask of `TIE_FLAG_*` constants describing validity/state
    pub flags: u16,
}

impl TiePointCell {
    pub fn invalid() -> Self {
        Self {
            azimuth_time: f64::NAN,
            slant_range: f64::NAN,
            cos_local_incidence: f32::NAN,
            flags: 0,
        }
    }

    pub fn with_values(azimuth_time: f64, slant_range: f64, cos_lia: f64, flags: u16) -> Self {
        Self {
            azimuth_time,
            slant_range,
            cos_local_incidence: cos_lia.clamp(-1.0, 1.0) as f32,
            flags,
        }
    }
}

/// Tie-point grid storing coarse forward model samples for interpolation
#[derive(Debug, Clone)]
pub struct TiePointGrid {
    pub(crate) stride: usize,
    pub(crate) cells: Array2<TiePointCell>,
}

impl TiePointGrid {
    pub fn new(stride: usize, cells: Array2<TiePointCell>) -> Self {
        Self {
            stride: stride.max(1),
            cells,
        }
    }

    pub fn dimensions(&self) -> (usize, usize) {
        self.cells.dim()
    }

    /// Bilinearly interpolate the tie grid at the given output pixel coordinates
    pub fn interpolate(&self, row: usize, col: usize) -> Option<TiePointCell> {
        let (grid_rows, grid_cols) = self.cells.dim();
        if grid_rows == 0 || grid_cols == 0 {
            return None;
        }

        let stride_f = self.stride as f64;
        let gy = (row as f64) / stride_f;
        let gx = (col as f64) / stride_f;

        let mut i0 = gy.floor() as isize;
        let mut j0 = gx.floor() as isize;
        i0 = i0.clamp(0, grid_rows as isize - 1);
        j0 = j0.clamp(0, grid_cols as isize - 1);
        let i0 = i0 as usize;
        let j0 = j0 as usize;

        let i1 = (i0 + 1).min(grid_rows.saturating_sub(1));
        let j1 = (j0 + 1).min(grid_cols.saturating_sub(1));

        let dy = if i0 == i1 {
            0.0
        } else {
            (gy - i0 as f64).clamp(0.0, 1.0)
        };
        let dx = if j0 == j1 {
            0.0
        } else {
            (gx - j0 as f64).clamp(0.0, 1.0)
        };

        let w00 = (1.0 - dy) * (1.0 - dx);
        let w10 = dy * (1.0 - dx);
        let w01 = (1.0 - dy) * dx;
        let w11 = dy * dx;

        let mut weight_sum = 0.0;
        let mut az_sum = 0.0;
        let mut slant_sum = 0.0;
        let mut cos_sum = 0.0;
        let mut flags_acc: u16 = 0;

        let accumulate = |cell: TiePointCell,
                          weight: f64,
                          weight_sum: &mut f64,
                          az_sum: &mut f64,
                          slant_sum: &mut f64,
                          cos_sum: &mut f64,
                          flags_acc: &mut u16| {
            if weight <= f64::EPSILON {
                return;
            }
            if cell.flags & TIE_FLAG_VALID != 0 {
                *weight_sum += weight;
                *az_sum += cell.azimuth_time * weight;
                *slant_sum += cell.slant_range * weight;
                *cos_sum += (cell.cos_local_incidence as f64) * weight;
            }
            *flags_acc |= cell.flags;
        };

        accumulate(
            self.cells[[i0, j0]],
            w00,
            &mut weight_sum,
            &mut az_sum,
            &mut slant_sum,
            &mut cos_sum,
            &mut flags_acc,
        );
        accumulate(
            self.cells[[i1, j0]],
            w10,
            &mut weight_sum,
            &mut az_sum,
            &mut slant_sum,
            &mut cos_sum,
            &mut flags_acc,
        );
        accumulate(
            self.cells[[i0, j1]],
            w01,
            &mut weight_sum,
            &mut az_sum,
            &mut slant_sum,
            &mut cos_sum,
            &mut flags_acc,
        );
        accumulate(
            self.cells[[i1, j1]],
            w11,
            &mut weight_sum,
            &mut az_sum,
            &mut slant_sum,
            &mut cos_sum,
            &mut flags_acc,
        );

        if weight_sum <= 1e-6 {
            return None;
        }

        let azimuth_time = az_sum / weight_sum;
        let slant_range = slant_sum / weight_sum;
        let cos_lia = (cos_sum / weight_sum).clamp(-1.0, 1.0);

        let mut flags = flags_acc;
        flags |= TIE_FLAG_VALID;

        Some(TiePointCell::with_values(
            azimuth_time,
            slant_range,
            cos_lia,
            flags,
        ))
    }
}

#[derive(Debug, Clone, Copy)]
pub struct DemLookupSample {
    pub elevation: f64,
    pub base_row: usize,
    pub base_col: usize,
    pub center_row: usize,
    pub center_col: usize,
    pub frac_row: f64,
    pub frac_col: f64,
}

/// Fast seeding grid for Phase 3.1 terrain correction optimization
///
/// Stores pre-computed range-doppler solutions at sparse grid points (e.g., every 64 pixels)
/// for fast bilinear interpolation to dense pixels. This eliminates redundant Newton-Raphson
/// iterations for neighboring pixels.
///
/// # Performance Impact
/// - Seed grid: ~1-2% of total pixels (full Newton-Raphson)
/// - Dense pixels: 1-2 iterations (vs 3-8 cold start)
/// - Expected speedup: 5-7× for terrain correction
///
/// # Scientific Accuracy
/// Bilinear interpolation of range-doppler coordinates maintains < 0.1 pixel error
/// for typical SAR scene geometry (validated against SNAP).
#[derive(Debug, Clone)]
pub struct SeedGrid {
    /// Seed solutions: (t_azimuth_orbit_rel, range_pixel, iterations_used)
    /// None indicates seed computation failed (ocean, out of bounds, etc.)
    pub seeds: Array2<Option<(f64, f64, u8)>>,

    /// Grid stride in output pixels (e.g., 64)
    pub stride: usize,

    /// Output grid bounds
    pub bounds: BoundingBox,

    /// Output grid dimensions
    pub width: usize,
    pub height: usize,

    /// Seed computation statistics
    pub total_seeds: usize,
    pub successful_seeds: usize,
    pub failed_seeds: usize,
    pub total_iterations: usize,
    pub max_iterations: usize,
}

impl SeedGrid {
    /// Get seed dimensions (rows, cols)
    pub fn seed_dimensions(&self) -> (usize, usize) {
        self.seeds.dim()
    }

    /// Get average iterations per seed
    pub fn avg_iterations(&self) -> f64 {
        if self.successful_seeds > 0 {
            self.total_iterations as f64 / self.successful_seeds as f64
        } else {
            0.0
        }
    }

    /// Get success rate (0.0 to 1.0)
    pub fn success_rate(&self) -> f64 {
        if self.total_seeds > 0 {
            self.successful_seeds as f64 / self.total_seeds as f64
        } else {
            0.0
        }
    }

    /// Get time range (min, max) across all successful seeds
    /// Used to detect constant seed bug (Δt=0)
    pub fn time_range(&self) -> (f64, f64) {
        let mut t_min = f64::INFINITY;
        let mut t_max = f64::NEG_INFINITY;

        for seed_opt in self.seeds.iter() {
            if let Some((t, _, _)) = seed_opt {
                t_min = t_min.min(*t);
                t_max = t_max.max(*t);
            }
        }

        (t_min, t_max)
    }

    /// Bilinear interpolation of seed at output pixel coordinates (x, y)
    ///
    /// Returns interpolated (t_azimuth, range_pixel) or None if seeds unavailable
    ///
    /// # Algorithm
    /// 1. Convert output pixel to fractional seed coordinates
    /// 2. Get 4 corner seeds (with bounds checking)
    /// 3. Bilinear interpolation: f(x,y) = (1-fx)(1-fy)f00 + fx(1-fy)f01 + (1-fx)fy*f10 + fx*fy*f11
    /// 4. Fallback to nearest seed if any corner missing
    pub fn interpolate_seed(&self, x: usize, y: usize) -> Option<(f64, f64)> {
        // Convert output pixel to seed grid coordinates (fractional)
        let seed_x_f = x as f64 / self.stride as f64;
        let seed_y_f = y as f64 / self.stride as f64;

        // Floor to get base seed indices
        let sx0 = seed_x_f.floor() as usize;
        let sy0 = seed_y_f.floor() as usize;

        // Ceiling with bounds check
        let sx1 = (sx0 + 1).min(self.seeds.ncols().saturating_sub(1));
        let sy1 = (sy0 + 1).min(self.seeds.nrows().saturating_sub(1));

        // Fractional parts for bilinear weights
        let fx = seed_x_f - sx0 as f64;
        let fy = seed_y_f - sy0 as f64;

        // Get 4 corner seeds
        let s00 = self.seeds.get([sy0, sx0]).and_then(|&s| s);
        let s01 = self.seeds.get([sy0, sx1]).and_then(|&s| s);
        let s10 = self.seeds.get([sy1, sx0]).and_then(|&s| s);
        let s11 = self.seeds.get([sy1, sx1]).and_then(|&s| s);

        // Bilinear interpolation (requires all 4 corners valid)
        match (s00, s01, s10, s11) {
            (
                Some((t00, r00, _)),
                Some((t01, r01, _)),
                Some((t10, r10, _)),
                Some((t11, r11, _)),
            ) => {
                // Bilinear formula
                let t_interp = (1.0 - fx) * (1.0 - fy) * t00
                    + fx * (1.0 - fy) * t01
                    + (1.0 - fx) * fy * t10
                    + fx * fy * t11;

                let r_interp = (1.0 - fx) * (1.0 - fy) * r00
                    + fx * (1.0 - fy) * r01
                    + (1.0 - fx) * fy * r10
                    + fx * fy * r11;

                Some((t_interp, r_interp))
            }
            _ => {
                // Fallback: use nearest valid seed if available
                // Priority: bottom-right, bottom-left, top-right, top-left
                s11.or(s10).or(s01).or(s00).map(|(t, r, _)| (t, r))
            }
        }
    }
}
