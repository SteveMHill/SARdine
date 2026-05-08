use crate::types::{SarError, SarResult};

/// Speckle filtering parameters
#[derive(Debug, Clone)]
pub struct SpeckleFilterParams {
    pub window_size: usize,
    pub num_looks: f32,
    pub edge_threshold: f32,
    pub damping_factor: f32,
    pub cv_threshold: f32,
    pub tile_size: usize,
}

impl Default for SpeckleFilterParams {
    fn default() -> Self {
        Self {
            window_size: 7,
            num_looks: 1.0,
            edge_threshold: 0.5,
            damping_factor: 1.0,
            cv_threshold: 0.5,
            tile_size: 0,
        }
    }
}

/// Tile size configuration for cache optimization
#[derive(Debug, Clone, Copy)]
pub enum TileSize {
    Small = 64,
    Medium = 128,
    Large = 256,
    Auto = 0,
}

impl SpeckleFilterParams {
    pub fn validate_window(&self) -> SarResult<()> {
        if self.window_size == 0 || self.window_size % 2 == 0 {
            return Err(SarError::Processing(
                "Window size must be odd and >= 1".to_string(),
            ));
        }
        Ok(())
    }
}
