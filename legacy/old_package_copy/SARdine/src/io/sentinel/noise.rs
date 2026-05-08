use super::*;

impl SlcReader {
    /// Find noise files for all available polarizations
    pub fn find_noise_files(&mut self) -> SarResult<HashMap<Polarization, String>> {
        log::debug!("Finding noise files");

        let files = self.list_files()?;
        let mut noise_files = HashMap::new();

        for file in files {
            if file.contains("annotation/calibration/")
                && file.ends_with(".xml")
                && file.contains("noise-")
            {
                if let Some(pol) = Self::extract_polarization(&file) {
                    if !noise_files.contains_key(&pol) {
                        noise_files.insert(pol, file);
                    }
                }
            }
        }

        log::info!("Found {} noise files", noise_files.len());
        Ok(noise_files)
    }

    /// Read thermal noise data for a specific polarization
    pub fn read_noise_data(
        &mut self,
        pol: Polarization,
    ) -> SarResult<crate::core::calibration::NoiseCoefficients> {
        log::debug!("Reading noise data for polarization {:?}", pol);

        let noise_files = self.find_noise_files()?;
        let file_path = noise_files.get(&pol).ok_or_else(|| {
            SarError::Processing(format!("No noise file found for polarization {:?}", pol))
        })?;

        let xml_content = self.read_file_as_string(file_path)?;
        let noise_data = crate::core::calibration::parse_noise_from_xml(&xml_content)?;

        log::info!("Successfully parsed noise data for polarization {:?}", pol);
        Ok(noise_data)
    }
}
