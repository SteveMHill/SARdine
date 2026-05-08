#![allow(dead_code, unused_variables)]
//! Download provider system - ASF DAAC provider for Sentinel-1 products

pub mod asf;

use super::manager::{DownloadConfig, ProviderConfig};
use super::products::{ProductMetadata, SearchParams};
use super::progress::ProgressCallback;
use crate::types::{SarError, SarResult};
use std::path::{Path, PathBuf};

/// Provider trait for different data sources
pub trait DownloadProvider: Send + Sync {
    /// Get provider name
    fn name(&self) -> &str;

    /// Check if provider supports product search
    fn supports_search(&self) -> bool;

    /// Check if provider requires authentication
    fn requires_auth(&self) -> bool;

    /// Check if provider is available (e.g., credentials valid)
    fn is_available(&self) -> bool;

    /// Resolve product ID to download URL
    /// Returns None if provider doesn't support direct ID resolution
    fn resolve_product_url(&self, product_id: &str) -> SarResult<Option<String>> {
        // Default implementation: no direct resolution
        Ok(None)
    }

    /// Search for products using catalog API
    fn search_products(&self, params: &SearchParams) -> SarResult<Vec<ProductMetadata>> {
        if !self.supports_search() {
            return Err(SarError::Processing(format!(
                "Provider {} does not support product search",
                self.name()
            )));
        }
        Err(SarError::Processing(format!(
            "Product search not yet implemented for provider {}",
            self.name()
        )))
    }

    /// Download a product from a URL with authentication
    fn download_product(
        &self,
        url: &str,
        output_path: &Path,
        progress: Option<&ProgressCallback>,
    ) -> SarResult<PathBuf> {
        // Default implementation uses generic download
        use crate::io::download::utils::download_from_url;
        download_from_url(url, Some(output_path), 3600)?;
        Ok(output_path.to_path_buf())
    }
}

/// Provider registry - manages available providers (ASF only)
pub struct ProviderRegistry {
    providers: Vec<Box<dyn DownloadProvider>>,
    enabled: Vec<String>,
}

impl ProviderRegistry {
    pub fn new(config: &DownloadConfig) -> SarResult<Self> {
        let mut providers: Vec<Box<dyn DownloadProvider>> = Vec::new();

        // ASF is the only supported provider
        let enabled = vec!["asf".to_string()];

        // Initialize ASF provider
        if let Some(provider_config) = config.provider_configs.get("asf") {
            providers.push(Box::new(asf::ASFProvider::new(provider_config)?));
        } else {
            // Try to get from environment
            let username = std::env::var("SARDINE_ASF_USERNAME").ok();
            let password = std::env::var("SARDINE_ASF_PASSWORD").ok();
            let token = std::env::var("SARDINE_ASF_TOKEN")
                .ok()
                .or_else(|| std::env::var("SARDINE_EDL_TOKEN").ok());
            let provider_config = ProviderConfig {
                username,
                password,
                token,
                base_url: None,
            };
            providers.push(Box::new(asf::ASFProvider::new(&provider_config)?));
        }

        Ok(Self { providers, enabled })
    }

    pub fn get_provider(&self, name: &str) -> Option<&dyn DownloadProvider> {
        self.providers
            .iter()
            .find(|p| p.name() == name)
            .map(|p| p.as_ref())
    }

    pub fn get_providers(&self) -> &[Box<dyn DownloadProvider>] {
        &self.providers
    }

    pub fn enabled(&self) -> &[String] {
        &self.enabled
    }
}
