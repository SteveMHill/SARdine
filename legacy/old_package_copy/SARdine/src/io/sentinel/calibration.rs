use super::*;

impl SlcReader {
    /// Find ALL calibration files for each polarization in the ZIP
    pub fn find_all_calibration_files(&mut self) -> SarResult<HashMap<Polarization, Vec<String>>> {
        log::debug!("Finding all calibration files (multi-subswath support)");

        let files = self.list_files()?;
        let mut calibration_files: HashMap<Polarization, Vec<String>> = HashMap::new();

        for file in files {
            // Only include actual calibration files, not noise files
            if Self::is_calibration_xml(&file) {
                if let Some(pol) = Self::extract_polarization(&file) {
                    calibration_files
                        .entry(pol)
                        .or_insert_with(Vec::new)
                        .push(file);
                }
            }
        }

        log::info!(
            "Found calibration files: {:?}",
            calibration_files
                .iter()
                .map(|(pol, files)| (pol, files.len()))
                .collect::<Vec<_>>()
        );
        Ok(calibration_files)
    }

    /// Find calibration files for each polarization in the ZIP (legacy - single file per pol)
    pub fn find_calibration_files(&mut self) -> SarResult<HashMap<Polarization, String>> {
        log::debug!("Finding calibration files (legacy single-file mode)");

        let files = self.list_files()?;
        let mut calibration_files = HashMap::new();

        for file in files {
            if Self::is_calibration_xml(&file) {
                if let Some(pol) = Self::extract_polarization(&file) {
                    if !calibration_files.contains_key(&pol) {
                        calibration_files.insert(pol, file);
                    }
                }
            }
        }

        log::info!("Found {} calibration files", calibration_files.len());
        Ok(calibration_files)
    }

    /// Read calibration data for a specific polarization (handles multiple subswaths)
    pub fn read_calibration_data(
        &mut self,
        pol: Polarization,
    ) -> SarResult<crate::core::calibration::CalibrationCoefficients> {
        log::debug!("Reading calibration data for polarization {:?}", pol);

        let all_calibration_files = self.find_all_calibration_files()?;
        let file_paths = all_calibration_files.get(&pol).ok_or_else(|| {
            SarError::Processing(format!(
                "No calibration files found for polarization {:?}",
                pol
            ))
        })?;

        log::info!(
            "Found {} calibration files for {:?}: {:?}",
            file_paths.len(),
            pol,
            file_paths
        );

        let mut merged_calibration = None;

        for file_path in file_paths {
            let xml_content = self.read_file_as_string(file_path)?;
            let calibration_data =
                crate::core::calibration::parse_calibration_from_xml(&xml_content)?;

            log::info!(
                "Parsed {} calibration vectors from {} (swath: {})",
                calibration_data.vectors.len(),
                file_path,
                calibration_data.swath
            );

            match merged_calibration {
                None => {
                    let mut base_calibration = calibration_data;
                    base_calibration.swath = "IW".to_string();
                    merged_calibration = Some(base_calibration);
                }
                Some(ref mut merged) => {
                    merged.vectors.extend(calibration_data.vectors);
                    log::info!(
                        "Merged calibration vectors - total now: {}",
                        merged.vectors.len()
                    );
                }
            }
        }

        let final_calibration = merged_calibration.ok_or_else(|| {
            SarError::Processing(format!(
                "No valid calibration data found for polarization {:?}",
                pol
            ))
        })?;

        log::info!(
            "Successfully merged calibration data for {:?} - total {} vectors from {} files",
            pol,
            final_calibration.vectors.len(),
            file_paths.len()
        );
        Ok(final_calibration)
    }

    /// Read calibration data with caching
    pub fn read_calibration_data_cached(
        &mut self,
        pol: Polarization,
    ) -> SarResult<CalibrationCoefficients> {
        let cache_key = format!("{}", pol);

        if let Some(cached) = self.calibration_cache.get(&cache_key) {
            log::debug!("Using cached calibration data for {}", pol);
            return Ok(cached.clone());
        }

        log::debug!("Reading calibration data for {} (not cached)", pol);
        let cal_data = self.read_calibration_data(pol)?;
        self.calibration_cache.insert(cache_key, cal_data.clone());

        Ok(cal_data)
    }

    /// Read calibration data for all available polarizations
    pub fn read_all_calibration_data(
        &mut self,
    ) -> SarResult<HashMap<Polarization, crate::core::calibration::CalibrationCoefficients>> {
        let mut all_cal_data = HashMap::new();

        let available_pols = self
            .find_all_annotation_files()?
            .keys()
            .cloned()
            .collect::<Vec<_>>();
        for pol in available_pols {
            match self.read_calibration_data(pol) {
                Ok(cal_data) => {
                    all_cal_data.insert(pol, cal_data);
                }
                Err(e) => {
                    log::warn!("Failed to read calibration data for {:?}: {}", pol, e);
                }
            }
        }

        Ok(all_cal_data)
    }
}
