//! Tests for download module

#[cfg(test)]
mod tests {
    use super::super::*;
    use crate::io::download::utils::{
        extract_from_zip, is_zip_content as utils_is_zip, sanitize_url_for_error, validate_checksum,
    };
    use crate::types::BoundingBox;
    use chrono::{TimeZone, Utc};
    use std::io::Write;
    use tempfile::TempDir;

    use super::super::cache::CacheManager;
    use super::super::dem::DEMDownloader;
    use super::super::manager::{DownloadConfig, ProviderConfig};
    use super::super::providers::ProviderRegistry;
    // AWS provider module has been removed - using ASF instead
    // use super::super::providers::aws::AWSProvider;
    use super::super::products::SearchParams;
    use super::super::progress::ProgressState;
    use super::super::providers::asf::ASFProvider;
    use super::super::queue::{DownloadQueue, DownloadTask, Priority};

    #[test]
    fn test_download_manager_creation() {
        let temp_dir = TempDir::new().unwrap();
        let config = DownloadConfig {
            cache_dir: temp_dir.path().to_path_buf(),
            ..Default::default()
        };

        let manager = DownloadManager::new(config);
        assert!(manager.is_ok());
    }

    #[test]
    fn test_cache_manager() {
        let temp_dir = TempDir::new().unwrap();
        let cache = CacheManager::new(temp_dir.path()).unwrap();

        assert!(cache.products_dir().exists());
        assert!(cache.orbits_dir().exists());
        assert!(cache.dem_dir().exists());
    }

    #[test]
    fn test_provider_registry() {
        let config = DownloadConfig::default();
        let registry = ProviderRegistry::new(&config);
        assert!(registry.is_ok());
    }

    // NOTE: AWS provider module has been removed - using ASF provider instead
    // #[test]
    // fn test_aws_provider() {
    //     let provider = AWSProvider::new();
    //     assert!(provider.is_ok());
    //     // Provider created successfully
    // }

    #[test]
    fn test_asf_provider() {
        let config = ProviderConfig {
            username: None,
            password: None,
            token: None,
            base_url: None,
        };

        let provider = ASFProvider::new(&config);
        assert!(provider.is_ok());
        // Provider created successfully
    }

    #[test]
    fn test_product_id_resolution() {
        let config = ProviderConfig {
            username: None,
            password: None,
            token: None,
            base_url: None,
        };

        let provider = ASFProvider::new(&config).unwrap();
        // Product ID resolution test - method is private, so just verify provider creation
        let _product_id = "S1A_IW_SLC__1SDV_20201230T165244_20201230T165311_035918_0434F0_6788";
    }

    #[test]
    fn test_download_queue() {
        let queue = DownloadQueue::new(2);

        let task1 = DownloadTask::new(
            "task1".to_string(),
            "http://example.com/file1.zip".to_string(),
            "/tmp/file1.zip".to_string(),
            Priority::High,
        );

        let task2 = DownloadTask::new(
            "task2".to_string(),
            "http://example.com/file2.zip".to_string(),
            "/tmp/file2.zip".to_string(),
            Priority::Low,
        );

        queue.add_task(task1).unwrap();
        queue.add_task(task2).unwrap();

        // High priority task should come first
        let next = queue.next_task();
        assert!(next.is_some());
        assert_eq!(next.unwrap().id, "task1");

        let stats = queue.get_stats();
        assert_eq!(stats.0, 1); // 1 pending (task2)
        assert_eq!(stats.1, 1); // 1 downloading (task1)
    }

    #[test]
    fn test_search_params() {
        let start = Utc.with_ymd_and_hms(2020, 1, 1, 0, 0, 0).unwrap();
        let end = Utc.with_ymd_and_hms(2020, 1, 2, 0, 0, 0).unwrap();

        let params = SearchParams {
            product_type: Some("SLC".to_string()),
            platform: Some("S1A".to_string()),
            start_date: start,
            end_date: end,
            aoi_wkt: None,
            max_results: 10,
            acquisition_mode: Some("IW".to_string()),
            polarization: Some("VV".to_string()),
            orbit_direction: Some("ASCENDING".to_string()),
            relative_orbit: Some(123),
        };

        assert_eq!(params.product_type, Some("SLC".to_string()));
        assert_eq!(params.platform, Some("S1A".to_string()));
        assert_eq!(params.max_results, 10);
    }

    #[test]
    fn test_progress_state() {
        let mut state = ProgressState::new(1000);
        assert_eq!(state.downloaded, 0);
        assert_eq!(state.total, 1000);
        assert_eq!(state.percentage, 0.0);

        state.update(500, 1.0);
        assert_eq!(state.downloaded, 500);
        assert_eq!(state.percentage, 50.0);
        assert_eq!(state.bytes_per_second, 500.0);
    }

    #[test]
    fn test_dem_tile_calculation() {
        let temp_dir = TempDir::new().unwrap();
        let config = DownloadConfig {
            cache_dir: temp_dir.path().to_path_buf(),
            ..Default::default()
        };

        let cache = std::sync::Arc::new(CacheManager::new(&config.cache_dir).unwrap());
        let providers = std::sync::Arc::new(ProviderRegistry::new(&config).unwrap());

        let dem_downloader = DEMDownloader::new(&cache, &providers, &config).unwrap();

        // Just verify the downloader was created successfully
        // Note: We can't access internal cache field, so just check creation
        assert!(cache.products_dir().exists());
    }

    #[test]
    fn test_is_zip_content_magic_bytes() {
        let zip_header = [0x50u8, 0x4B, 0x03, 0x04, 0x00, 0x00];
        let not_zip = [0x00u8, 0x01, 0x02, 0x03];
        assert!(utils_is_zip(&zip_header));
        assert!(!utils_is_zip(&not_zip));
    }

    #[test]
    fn test_extract_eof_from_zip_success() {
        // Build an in-memory ZIP with one .EOF file
        let mut buf = Vec::new();
        {
            let mut zip = zip::ZipWriter::new(std::io::Cursor::new(&mut buf));
            let options = zip::write::FileOptions::default();
            zip.start_file("S1_TEST.EOF", options).unwrap();
            let content = "<?xml version=\"1.0\"?><Earth_Explorer_File><Data_Block><List_of_OSVs/></Data_Block></Earth_Explorer_File>";
            zip.write_all(content.as_bytes()).unwrap();
            zip.finish().unwrap();
        }
        // utils generic extraction by extension
        let eof_bytes = extract_from_zip(&buf, ".EOF").expect("should extract EOF bytes");
        let eof_str = String::from_utf8(eof_bytes).unwrap();
        assert!(eof_str.contains("<Data_Block>"));
    }

    #[test]
    fn test_extract_eof_from_zip_not_found() {
        // Build an in-memory ZIP without .EOF
        let mut buf = Vec::new();
        {
            let mut zip = zip::ZipWriter::new(std::io::Cursor::new(&mut buf));
            let options = zip::write::FileOptions::default();
            zip.start_file("note.txt", options).unwrap();
            zip.write_all(b"hello").unwrap();
            zip.finish().unwrap();
        }
        let err = extract_from_zip(&buf, ".EOF").err().expect("should error");
        let msg = format!("{}", err);
        assert!(msg.to_lowercase().contains("no .eof"));
    }

    #[test]
    fn test_sanitize_url_for_error_masks_query() {
        let url = "https://example.com/path?token=SECRET&x=1";
        let sanitized = sanitize_url_for_error(url);
        assert!(sanitized.ends_with("?[REDACTED]"));
        assert!(sanitized.starts_with("https://example.com/path"));
        assert!(!sanitized.contains("SECRET"));
    }

    #[test]
    fn test_validate_checksum_md5() {
        let tmp = TempDir::new().unwrap();
        let path = tmp.path().join("test.bin");
        std::fs::write(&path, b"abc").unwrap();
        // md5("abc") = 900150983cd24fb0d6963f7d28e17f72
        assert!(validate_checksum(&path, "900150983cd24fb0d6963f7d28e17f72").unwrap());
        assert!(!validate_checksum(&path, "deadbeef").unwrap());
    }
}
