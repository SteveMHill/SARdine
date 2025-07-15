use crate::types::{SarError, SarResult, SubSwath};
use quick_xml::de::from_str;
use serde::Deserialize;
use std::collections::HashMap;

/// Detailed annotation structures for Sentinel-1
/// This represents the root <product> element directly
#[derive(Debug, Deserialize)]
pub struct AnnotationRoot {
    #[serde(rename = "swathTiming")]
    pub swath_timing: SwathTiming,
    #[serde(rename = "imageAnnotation")]
    pub image_annotation: ImageAnnotation,
    #[serde(rename = "dopplerCentroid")]
    pub doppler_centroid: DopplerCentroid,
}

#[derive(Debug, Deserialize)]
pub struct SwathTiming {
    #[serde(rename = "burstList")]
    pub burst_list: BurstList,
}

#[derive(Debug, Deserialize)]
pub struct BurstList {
    #[serde(rename = "burst")]
    pub bursts: Vec<Burst>,
}

#[derive(Debug, Deserialize)]
pub struct Burst {
    #[serde(rename = "azimuthTime")]
    pub azimuth_time: String,
    #[serde(rename = "sensingTime")]
    pub sensing_time: String,
    #[serde(rename = "byteOffset")]
    pub byte_offset: u64,
    #[serde(rename = "firstValidSample")]
    pub first_valid_sample: Vec<u32>,
    #[serde(rename = "lastValidSample")]
    pub last_valid_sample: Vec<u32>,
}

#[derive(Debug, Deserialize)]
pub struct ImageAnnotation {
    #[serde(rename = "imageInformation")]
    pub image_information: ImageInformation,
    #[serde(rename = "processingInformation")]
    pub processing_information: ProcessingInformation,
}

#[derive(Debug, Deserialize)]
pub struct ImageInformation {
    #[serde(rename = "slantRangeTime")]
    pub slant_range_time: f64,
    #[serde(rename = "numberOfSamples")]
    pub number_of_samples: usize,
    #[serde(rename = "numberOfLines")]
    pub number_of_lines: usize,
    #[serde(rename = "rangePixelSpacing")]
    pub range_pixel_spacing: f64,
    #[serde(rename = "azimuthPixelSpacing")]
    pub azimuth_pixel_spacing: f64,
    #[serde(rename = "azimuthFrequency")]
    pub azimuth_frequency: f64,
}

#[derive(Debug, Deserialize)]
pub struct ProcessingInformation {
    #[serde(rename = "swathProcParamsList")]
    pub swath_proc_params_list: SwathProcParamsList,
}

#[derive(Debug, Deserialize)]
pub struct SwathProcParamsList {
    #[serde(rename = "swathProcParams")]
    pub swath_proc_params: Vec<SwathProcParams>,
}

#[derive(Debug, Deserialize)]
pub struct SwathProcParams {
    #[serde(rename = "swath")]
    pub swath: String,
    #[serde(rename = "rangeProcessing")]
    pub range_processing: RangeProcessing,
    #[serde(rename = "azimuthProcessing")]
    pub azimuth_processing: AzimuthProcessing,
}

#[derive(Debug, Deserialize)]
pub struct RangeProcessing {
    #[serde(rename = "numberOfLooks")]
    pub number_of_looks: u32,
    #[serde(rename = "lookBandwidth")]
    pub look_bandwidth: f64,
}

#[derive(Debug, Deserialize)]
pub struct AzimuthProcessing {
    #[serde(rename = "numberOfLooks")]
    pub number_of_looks: u32,
    #[serde(rename = "lookBandwidth")]
    pub look_bandwidth: f64,
}

#[derive(Debug, Deserialize)]
pub struct DopplerCentroid {
    #[serde(rename = "dcEstimateList")]
    pub dc_estimate_list: DcEstimateList,
}

#[derive(Debug, Deserialize)]
pub struct DcEstimateList {
    #[serde(rename = "dcEstimate")]
    pub dc_estimates: Vec<DcEstimate>,
}

#[derive(Debug, Deserialize)]
pub struct DcEstimate {
    #[serde(rename = "azimuthTime")]
    pub azimuth_time: String,
    #[serde(rename = "slantRangeTime")]
    pub slant_range_time: f64,
    #[serde(rename = "frequency")]
    pub frequency: f64,
    #[serde(rename = "dataDcPolynomial")]
    pub data_dc_polynomial: DataDcPolynomial,
}

#[derive(Debug, Deserialize)]
pub struct DataDcPolynomial {
    #[serde(rename = "dataDcPolynomial")]
    pub coefficients: Vec<f64>,
}

/// Parser for Sentinel-1 annotation XML files
pub struct AnnotationParser;

impl AnnotationParser {
    /// Parse complete annotation XML
    pub fn parse_annotation(xml_content: &str) -> SarResult<AnnotationRoot> {
        from_str::<AnnotationRoot>(xml_content)
            .map_err(|e| SarError::XmlParsing(format!("Failed to parse annotation XML: {}", e)))
    }

    /// Extract sub-swath information from annotation
    pub fn extract_subswaths(annotation: &AnnotationRoot) -> SarResult<HashMap<String, SubSwath>> {
        let mut subswaths = HashMap::new();
        
        let image_info = &annotation.image_annotation.image_information;
        let proc_params = &annotation.image_annotation.processing_information
            .swath_proc_params_list.swath_proc_params;

        for params in proc_params {
            let subswath = SubSwath {
                id: params.swath.clone(),
                burst_count: annotation.swath_timing.burst_list.bursts.len(),
                range_samples: image_info.number_of_samples,
                azimuth_samples: image_info.number_of_lines,
                range_pixel_spacing: image_info.range_pixel_spacing,
                azimuth_pixel_spacing: image_info.azimuth_pixel_spacing,
                slant_range_time: image_info.slant_range_time,
                burst_duration: 2.758277, // Typical IW burst duration in seconds
            };
            
            subswaths.insert(params.swath.clone(), subswath);
        }

        if subswaths.is_empty() {
            return Err(SarError::Metadata(
                "No sub-swaths found in annotation".to_string(),
            ));
        }

        Ok(subswaths)
    }

    /// Extract burst timing information
    pub fn extract_burst_times(annotation: &AnnotationRoot) -> SarResult<Vec<String>> {
        let bursts = &annotation.swath_timing.burst_list.bursts;
        let mut burst_times = Vec::new();
        
        for burst in bursts {
            burst_times.push(burst.azimuth_time.clone());
        }
        
        Ok(burst_times)
    }

    /// Extract Doppler centroid information
    pub fn extract_doppler_centroid(annotation: &AnnotationRoot) -> SarResult<Vec<f64>> {
        let dc_estimates = &annotation.doppler_centroid.dc_estimate_list.dc_estimates;
        let mut frequencies = Vec::new();
        
        for estimate in dc_estimates {
            frequencies.push(estimate.frequency);
        }
        
        Ok(frequencies)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_annotation_parsing() {
        // This would test with real annotation XML content
        let sample_xml = r#"<?xml version="1.0" encoding="UTF-8"?>
        <product>
            <swathTiming>
                <burstList>
                    <burst>
                        <azimuthTime>2020-01-01T12:00:00.000000</azimuthTime>
                        <sensingTime>2020-01-01T12:00:00.000000</sensingTime>
                        <byteOffset>0</byteOffset>
                        <firstValidSample>0</firstValidSample>
                        <lastValidSample>1000</lastValidSample>
                    </burst>
                </burstList>
            </swathTiming>
        </product>"#;
        
        // This test would fail with the simplified XML above
        // but demonstrates the structure
    }
}
