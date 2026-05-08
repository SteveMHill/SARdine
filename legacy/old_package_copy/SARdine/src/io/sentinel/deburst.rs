#![allow(dead_code, unused_variables)]
use super::*;
use crate::core::deburst::DeburstConfig;

impl SlcReader {
    /// Shared deburst preparation: reads SLC, extracts bursts, and builds a TopSarDeburstProcessor
    fn prepare_deburst(
        &mut self,
        pol: Polarization,
    ) -> SarResult<(SarImage, crate::core::deburst::TopSarDeburstProcessor)> {
        use crate::core::deburst::{DeburstProcessor, TopSarDeburstProcessor};

        log::info!("Starting deburst preparation for polarization {:?}", pol);

        let slc_data = self.read_slc_data(pol)?;
        log::debug!("Read SLC data with dimensions: {:?}", slc_data.dim());

        let annotation_files = self.find_all_annotation_files()?;
        let files = annotation_files.get(&pol).ok_or_else(|| {
            SarError::Processing(format!("No annotation files found for {:?}", pol))
        })?;
        let first = files
            .first()
            .ok_or_else(|| SarError::Processing(format!("Empty annotation list for {:?}", pol)))?;

        let annotation_content = self.read_file_as_string(first)?;

        let subswath_id = first
            .to_lowercase()
            .split('-')
            .find(|s| s.starts_with("iw"))
            .map(|s| s.to_uppercase());

        // Use DeburstProcessor's static method for burst extraction
        let burst_info = DeburstProcessor::extract_burst_info_from_annotation_with_subswath(
            &annotation_content,
            slc_data.dim().0,
            slc_data.dim().1,
            None,
        )?;

        log::info!("Burst info extracted: {} bursts", burst_info.len(),);

        let cache_dir = std::env::var("SARDINE_ORBIT_CACHE")
            .map(std::path::PathBuf::from)
            .map_err(|_| {
                SarError::Processing(
                    "SARDINE_ORBIT_CACHE not set. Required for orbit-based velocity in deburst."
                        .to_string(),
                )
            })?;
        let orbit_data = self.rget_orbit_data(Some(cache_dir.as_path()))?;
        if orbit_data.state_vectors.is_empty() {
            return Err(SarError::Processing(
                "No orbit state vectors available for velocity calculation".to_string(),
            ));
        }
        let v = &orbit_data.state_vectors[0].velocity;
        let satellite_velocity = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        // Tightened bounds: 7200-7800 m/s (typical S1 is 7500-7700 m/s)
        if satellite_velocity < 7200.0 || satellite_velocity > 7800.0 {
            return Err(SarError::Processing(format!(
                "Invalid satellite velocity: {:.1} m/s (expected 7200-7800 m/s)",
                satellite_velocity
            )));
        }

        let config = DeburstConfig::default();
        let processor = TopSarDeburstProcessor::new(burst_info, config, satellite_velocity);
        Ok((slc_data, processor))
    }

    /// Deburst SLC data for a specific polarization
    pub fn deburst_slc(&mut self, pol: Polarization) -> SarResult<SarImage> {
        use crate::core::deburst::{DeburstConfig, TopSarDeburstProcessor};

        log::info!("Starting deburst for polarization {:?}", pol);

        let slc_data = self.read_slc_data(pol)?;
        log::debug!("Read SLC data with dimensions: {:?}", slc_data.dim());

        let annotation_files = self.find_all_annotation_files()?;
        let files = annotation_files.get(&pol).ok_or_else(|| {
            SarError::Processing(format!("No annotation files found for {:?}", pol))
        })?;
        let first = files
            .first()
            .ok_or_else(|| SarError::Processing(format!("Empty annotation list for {:?}", pol)))?;

        let annotation_content = self.read_file_as_string(first)?;

        // Use DeburstProcessor's static method for burst extraction
        use crate::core::deburst::DeburstProcessor;
        let burst_info = DeburstProcessor::extract_burst_info_from_annotation_with_subswath(
            &annotation_content,
            slc_data.dim().0,
            slc_data.dim().1,
            None,
        )?;

        log::info!("Burst info extracted: {} bursts", burst_info.len());

        let cache_dir = std::env::var("SARDINE_ORBIT_CACHE")
            .map(std::path::PathBuf::from)
            .map_err(|_| {
                SarError::Processing(
                    "SARDINE_ORBIT_CACHE not set. Required for orbit-based velocity in deburst."
                        .to_string(),
                )
            })?;
        let orbit_data = self.rget_orbit_data(Some(cache_dir.as_path()))?;
        if orbit_data.state_vectors.is_empty() {
            return Err(SarError::Processing(
                "No orbit state vectors available for velocity calculation".to_string(),
            ));
        }
        let v = &orbit_data.state_vectors[0].velocity;
        let satellite_velocity = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();

        // Build optimized config
        let mut config = DeburstConfig::default();
        config.blend_overlap = true;
        // FIX: Increased blend_lines from 10 to 50 for smoother burst boundaries
        // This eliminates visible stripes between bursts (typical overlaps are 50-200 lines)
        config.blend_lines = 50;
        config.apply_deramp = true;
        config.preserve_phase = true;
        // Range-dependent deramp is always enabled (required for TOPS IW)
        config.use_annotation_timing = true;
        config.enable_bilinear_interp = false;
        config.enable_hit_count_mask = true;
        config.power_preservation_check = true;
        config.enable_row_equalization = true;
        config.fill_small_gaps = true;
        config.enable_loop_unrolling = true;
        config.enable_simd = true;
        config.enable_tiled_processing = false;
        config.tile_lines = 2048;
        config.tile_samples = 4096;

        let processor = TopSarDeburstProcessor::new(burst_info, config, satellite_velocity);
        let deburst_data = processor.deburst_topsar(&slc_data)?;

        log::info!(
            "Deburst completed. Output dimensions: {:?}",
            deburst_data.dim()
        );
        Ok(deburst_data)
    }

    /// Deburst SLC data with custom configuration for TOPSAR processing
    pub fn deburst_slc_with_config(
        &mut self,
        pol: Polarization,
        config: DeburstConfig,
    ) -> SarResult<SarImage> {
        use crate::core::deburst::{DeburstProcessor, TopSarDeburstProcessor};

        let slc_data = self.read_slc_data(pol)?;
        let annotation_files = self.find_all_annotation_files()?;
        let files = annotation_files.get(&pol).ok_or_else(|| {
            SarError::Processing(format!("No annotation files found for {:?}", pol))
        })?;
        let first = files
            .first()
            .ok_or_else(|| SarError::Processing(format!("Empty annotation list for {:?}", pol)))?;

        let annotation_content = self.read_file_as_string(first)?;

        let burst_info = DeburstProcessor::extract_burst_info_from_annotation_with_subswath(
            &annotation_content,
            slc_data.dim().0,
            slc_data.dim().1,
            None,
        )?;

        let cache_dir = std::env::var("SARDINE_ORBIT_CACHE")
            .map(std::path::PathBuf::from)
            .map_err(|_| {
                SarError::Processing(
                    "SARDINE_ORBIT_CACHE not set. Required for orbit-based velocity in deburst."
                        .to_string(),
                )
            })?;
        let orbit_data = self.rget_orbit_data(Some(cache_dir.as_path()))?;
        let v = &orbit_data.state_vectors[0].velocity;
        let satellite_velocity = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();

        let processor = TopSarDeburstProcessor::new(burst_info, config, satellite_velocity);
        let deburst_data = processor.deburst_topsar(&slc_data)?;

        log::info!(
            "TOPSAR deburst completed. Output dimensions: {:?}",
            deburst_data.dim()
        );
        Ok(deburst_data)
    }
}
