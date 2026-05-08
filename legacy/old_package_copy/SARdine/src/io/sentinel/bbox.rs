use super::*;

impl SlcReader {
    /// Get bounding box for specific subswaths only (e.g., ["IW1"] or ["IW1", "IW2"])
    /// This returns a bbox that covers ONLY the requested subswaths, not the full product.
    pub fn get_subswath_bounding_box(
        &mut self,
        subswaths: &[String],
        pol: Polarization,
    ) -> SarResult<crate::types::BoundingBox> {
        log::info!(
            "🌐 Extracting bounding box for specific subswaths: {:?}",
            subswaths
        );

        let all_annotation_files = self.find_all_annotation_files()?;
        let pol_files = all_annotation_files.get(&pol).ok_or_else(|| {
            SarError::Processing(format!(
                "No annotation files found for polarization {}",
                pol
            ))
        })?;

        let mut merged_min_lat = f64::INFINITY;
        let mut merged_max_lat = f64::NEG_INFINITY;
        let mut merged_min_lon = f64::INFINITY;
        let mut merged_max_lon = f64::NEG_INFINITY;
        let mut processed_count = 0;

        let subswaths_lower: Vec<String> = subswaths.iter().map(|s| s.to_lowercase()).collect();

        for annotation_file in pol_files {
            let file_lower = annotation_file.to_lowercase();
            let is_requested = subswaths_lower.iter().any(|sw| file_lower.contains(sw));

            if !is_requested {
                log::debug!(
                    "Skipping annotation file {} (not in requested subswaths)",
                    annotation_file
                );
                continue;
            }

            match self.read_annotation_file_raw(annotation_file) {
                Ok(sar_metadata) => {
                    let bbox = &sar_metadata.bounding_box;

                    if bbox.min_lat < bbox.max_lat
                        && bbox.min_lon < bbox.max_lon
                        && bbox.min_lat != 0.0
                        && bbox.max_lat != 0.0
                    {
                        merged_min_lat = merged_min_lat.min(bbox.min_lat);
                        merged_max_lat = merged_max_lat.max(bbox.max_lat);
                        merged_min_lon = merged_min_lon.min(bbox.min_lon);
                        merged_max_lon = merged_max_lon.max(bbox.max_lon);
                        processed_count += 1;

                        log::info!(
                            "✅ Subswath {} bbox: [{:.6}, {:.6}, {:.6}, {:.6}]",
                            annotation_file,
                            bbox.min_lon,
                            bbox.min_lat,
                            bbox.max_lon,
                            bbox.max_lat
                        );
                    }
                }
                Err(e) => {
                    log::warn!(
                        "⚠️  Failed to read annotation file {}: {}",
                        annotation_file,
                        e
                    );
                }
            }
        }

        if processed_count == 0 {
            return Err(SarError::Processing(format!(
                "No valid bounding boxes found for subswaths {:?}",
                subswaths
            )));
        }

        let merged_bbox = crate::types::BoundingBox {
            min_lat: merged_min_lat,
            max_lat: merged_max_lat,
            min_lon: merged_min_lon,
            max_lon: merged_max_lon,
        };

        log::info!(
            "🎯 SUBSWATH-SPECIFIC BOUNDING BOX for {:?}: [{:.6}, {:.6}, {:.6}, {:.6}]",
            subswaths,
            merged_bbox.min_lon,
            merged_bbox.min_lat,
            merged_bbox.max_lon,
            merged_bbox.max_lat
        );
        log::info!(
            "   📏 Coverage: {:.3}deg x {:.3}deg ({:.1}km x {:.1}km)",
            merged_bbox.max_lon - merged_bbox.min_lon,
            merged_bbox.max_lat - merged_bbox.min_lat,
            (merged_bbox.max_lon - merged_bbox.min_lon) * 111.0,
            (merged_bbox.max_lat - merged_bbox.min_lat) * 111.0
        );

        Ok(merged_bbox)
    }

    /// Extract merged bounding box from ALL annotation files for complete scene coverage
    pub fn extract_merged_bounding_box_all_subswaths(
        &mut self,
        pol: Polarization,
    ) -> SarResult<crate::types::BoundingBox> {
        log::debug!(
            "🌐 Starting merged bounding box extraction for polarization {}",
            pol
        );
        log::info!("🌐 Extracting merged bounding box from ALL subswaths for complete coverage");

        log::debug!("📁 Getting all annotation files...");
        let all_annotation_files = match self.find_all_annotation_files() {
            Ok(files) => files,
            Err(e) => {
                log::error!("Failed to find annotation files: {}", e);
                return Err(e);
            }
        };

        log::debug!(
            "📁 Found annotation files for {} polarizations",
            all_annotation_files.len()
        );

        let pol_files = all_annotation_files.get(&pol).ok_or_else(|| {
            log::error!(
                "No annotation files found for polarization {}",
                pol
            );
            SarError::Processing(format!(
                "No annotation files found for polarization {}",
                pol
            ))
        })?;

        log::debug!(
            "📁 Found {} annotation files for {}",
            pol_files.len(),
            pol
        );
        log::info!(
            "📁 Found {} annotation files for {}: merging bounding boxes",
            pol_files.len(),
            pol
        );

        let mut merged_min_lat = f64::INFINITY;
        let mut merged_max_lat = f64::NEG_INFINITY;
        let mut merged_min_lon = f64::INFINITY;
        let mut merged_max_lon = f64::NEG_INFINITY;

        let mut processed_count = 0;

        for annotation_file in pol_files {
            log::debug!("📄 Processing annotation file: {}", annotation_file);

            match self.read_annotation_file_raw(annotation_file) {
                Ok(sar_metadata) => {
                    let bbox = &sar_metadata.bounding_box;

                    if bbox.min_lat < bbox.max_lat
                        && bbox.min_lon < bbox.max_lon
                        && bbox.min_lat != 0.0
                        && bbox.max_lat != 0.0
                        && bbox.min_lon != 0.0
                        && bbox.max_lon != 0.0
                    {
                        merged_min_lat = merged_min_lat.min(bbox.min_lat);
                        merged_max_lat = merged_max_lat.max(bbox.max_lat);
                        merged_min_lon = merged_min_lon.min(bbox.min_lon);
                        merged_max_lon = merged_max_lon.max(bbox.max_lon);

                        processed_count += 1;

                        log::info!(
                            "✅ Subswath {}: [{:.6}, {:.6}, {:.6}, {:.6}]",
                            processed_count,
                            bbox.min_lon,
                            bbox.min_lat,
                            bbox.max_lon,
                            bbox.max_lat
                        );
                    } else {
                        log::warn!("⚠️  Invalid bounding box in {}, skipping", annotation_file);
                    }
                }
                Err(e) => {
                    log::warn!(
                        "⚠️  Failed to read annotation file {}: {}",
                        annotation_file,
                        e
                    );
                }
            }
        }

        if processed_count == 0 {
            return Err(SarError::Processing(
                "No valid bounding boxes found in any annotation files".to_string(),
            ));
        }

        let merged_bbox = crate::types::BoundingBox {
            min_lat: merged_min_lat,
            max_lat: merged_max_lat,
            min_lon: merged_min_lon,
            max_lon: merged_max_lon,
        };

        log::info!("🎯 MERGED BOUNDING BOX from {} subswaths:", processed_count);
        log::info!(
            "   📍 [{:.6}, {:.6}, {:.6}, {:.6}]",
            merged_bbox.min_lon,
            merged_bbox.min_lat,
            merged_bbox.max_lon,
            merged_bbox.max_lat
        );
        log::info!(
            "   📏 Coverage: {:.3}deg x {:.3}deg ({:.1}km x {:.1}km)",
            merged_bbox.max_lon - merged_bbox.min_lon,
            merged_bbox.max_lat - merged_bbox.min_lat,
            (merged_bbox.max_lon - merged_bbox.min_lon) * 111.0,
            (merged_bbox.max_lat - merged_bbox.min_lat) * 111.0
        );

        Ok(merged_bbox)
    }
}
