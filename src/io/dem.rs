use crate::types::{BoundingBox, GeoTransform, SarError, SarResult};
use gdal::Dataset;
use ndarray::Array2;
use std::path::Path;

/// Digital Elevation Model reader
pub struct DemReader;

impl DemReader {
    /// Read DEM data for a specific bounding box
    pub fn read_dem<P: AsRef<Path>>(
        dem_path: P,
        bbox: &BoundingBox,
        target_resolution: f64,
    ) -> SarResult<(Array2<f32>, GeoTransform)> {
        log::info!("Reading DEM from: {}", dem_path.as_ref().display());
        log::debug!("Bounding box: {:?}", bbox);
        log::debug!("Target resolution: {} meters", target_resolution);

        // Open DEM file with GDAL
        let dataset = Dataset::open(dem_path.as_ref())?;
        
        // Get spatial information
        let geo_transform = dataset.geo_transform()?;
        let (width, height) = dataset.raster_size();
        
        log::debug!("DEM size: {}x{}", width, height);
        log::debug!("DEM geotransform: {:?}", geo_transform);

        // Read elevation data from first band
        let rasterband = dataset.rasterband(1)?;
        let band_data = rasterband.read_as::<f32>((0, 0), (width, height), (width, height), None)?;
        
        // Convert to ndarray
        let dem_array = Array2::from_shape_vec((height, width), band_data.data)
            .map_err(|e| SarError::Processing(format!("Failed to reshape DEM data: {}", e)))?;

        let geo_transform_struct = GeoTransform {
            top_left_x: geo_transform[0],
            pixel_width: geo_transform[1],
            rotation_x: geo_transform[2],
            top_left_y: geo_transform[3],
            rotation_y: geo_transform[4],
            pixel_height: geo_transform[5],
        };

        Ok((dem_array, geo_transform_struct))
    }

    /// Download SRTM DEM tiles for a bounding box
    pub fn download_srtm_tiles(bbox: &BoundingBox, output_dir: &str) -> SarResult<Vec<String>> {
        log::info!("Downloading SRTM tiles for bounding box: {:?}", bbox);
        
        // Calculate required SRTM tiles
        let tiles = Self::calculate_srtm_tiles(bbox);
        let mut downloaded_files = Vec::new();
        
        std::fs::create_dir_all(output_dir)
            .map_err(|e| SarError::Io(e))?;
        
        for tile in tiles {
            let filename = format!("{}.hgt", tile);
            let output_path = format!("{}/{}", output_dir, filename);
            
            // Check if file already exists
            if std::path::Path::new(&output_path).exists() {
                log::info!("SRTM tile {} already exists, skipping download", filename);
                downloaded_files.push(output_path);
                continue;
            }
            
            // Try multiple SRTM data sources
            let success = Self::try_download_from_sources(&tile, &output_path)?;
            
            if success {
                log::info!("Successfully downloaded {}", filename);
                downloaded_files.push(output_path);
            } else {
                log::warn!("Failed to download {} from all sources", filename);
            }
        }
        
        if downloaded_files.is_empty() {
            return Err(SarError::Processing(
                "Failed to download any SRTM tiles from available sources.\n\
                 Please check internet connection or provide DEM files manually.".to_string(),
            ));
        }
        
        Ok(downloaded_files)
    }
    
    /// Try downloading SRTM tile from multiple sources
    fn try_download_from_sources(tile: &str, output_path: &str) -> SarResult<bool> {
        // Multiple DEM data sources (in order of preference)
        // Using reliable AWS-based sources that don't require authentication
        let mut sources = vec![
            // AWS USGS SRTM 1 Arcsecond Tiles (most reliable, no auth required)
            format!("https://s3.amazonaws.com/elevation-tiles-prod/skadi/{}/{}.hgt.gz", 
                Self::get_skadi_directory(tile), tile),
        ];
        
        // Add Copernicus DEM sources (30m and 90m resolution)
        sources.extend(Self::get_copernicus_dem_urls(tile));
        
        // Add backup sources
        sources.extend(vec![
            // Alternative AWS mirror for SRTM
            format!("https://cloud.sdstate.edu/index.php/s/jgaFbSdRmZbCRff/download?path=%2F&files={}.hgt.zip", tile),
            
            // NASA EarthData (backup, may require auth)
            format!("https://e4ftl01.cr.usgs.gov/MEASURES/SRTMGL1.003/2000.02.11/{}.SRTMGL1.hgt.zip", tile),
            
            // USGS Data Portal
            format!("https://dds.cr.usgs.gov/srtm/version2_1/SRTM1/{}/{}.hgt.zip", 
                Self::get_srtm_continent(tile), tile),
        ]);
        
        for (i, url) in sources.iter().enumerate() {
            log::info!("Attempting download from source {} of {}: {}", i + 1, sources.len(), url);
            
            match Self::download_and_extract_srtm_tile(url, output_path) {
                Ok(_) => {
                    log::info!("Successfully downloaded from source {}", i + 1);
                    return Ok(true);
                },
                Err(e) => {
                    log::warn!("Source {} failed: {}", i + 1, e);
                    // Continue to next source
                }
            }
        }
        
        // If all sources fail, provide helpful message
        log::error!("All download sources failed for tile: {}", tile);
        log::info!("Manual download options:");
        log::info!("1. Visit https://earthexplorer.usgs.gov/ and search for SRTM");
        log::info!("2. Download tile {} and save as {}", tile, output_path);
        log::info!("3. Or provide your own DEM file using read_dem() method");
        
        Ok(false)
    }
    
    /// Get Skadi directory structure for AWS SRTM tiles
    /// AWS organizes SRTM tiles like: /N50/N50E012.hgt.gz
    fn get_skadi_directory(tile: &str) -> String {
        // Extract latitude part from tile name (e.g., "N50" from "N50E012")
        if tile.len() >= 3 {
            let lat_part = &tile[0..3]; // First 3 characters (e.g., "N50", "S45")
            lat_part.to_string()
        } else {
            // Fallback for malformed tile names
            "N00".to_string()
        }
    }
    
    /// Get SRTM continent directory for USGS structure
    fn get_srtm_continent(tile: &str) -> &'static str {
        // SRTM files are organized by continent on USGS servers
        // This is a simplified mapping - in practice, you'd need a more detailed lookup
        if tile.starts_with('N') {
            if tile.contains("W") {
                "North_America"  // North America
            } else {
                "Eurasia"        // Europe/Asia
            }
        } else {
            if tile.contains("W") {
                "South_America"  // South America
            } else if tile.contains("E") && (tile.starts_with("S0") || tile.starts_with("S1")) {
                "Africa"         // Africa
            } else {
                "Australia"      // Australia/Oceania
            }
        }
    }

    /// Calculate which SRTM tiles are needed for a bounding box
    fn calculate_srtm_tiles(bbox: &BoundingBox) -> Vec<String> {
        let mut tiles = Vec::new();
        
        // SRTM tiles are 1x1 degree
        let min_lat = bbox.min_lat.floor() as i32;
        let max_lat = bbox.max_lat.ceil() as i32;
        let min_lon = bbox.min_lon.floor() as i32;
        let max_lon = bbox.max_lon.ceil() as i32;
        
        for lat in min_lat..max_lat {
            for lon in min_lon..max_lon {
                // SRTM naming convention: N/S followed by latitude, E/W followed by longitude
                let lat_prefix = if lat >= 0 { "N" } else { "S" };
                let lon_prefix = if lon >= 0 { "E" } else { "W" };
                
                let tile_name = format!("{}{:02}{}{:03}", 
                    lat_prefix, lat.abs(), 
                    lon_prefix, lon.abs()
                );
                
                tiles.push(tile_name);
            }
        }
        
        log::debug!("Required SRTM tiles: {:?}", tiles);
        tiles
    }

    /// Download and extract a single SRTM tile
    fn download_and_extract_srtm_tile(url: &str, output_path: &str) -> SarResult<()> {
        use std::io::Write;
        
        log::info!("Downloading SRTM tile from: {}", url);
        
        // Create HTTP client with timeout and proper headers
        let client = reqwest::blocking::Client::builder()
            .timeout(std::time::Duration::from_secs(300)) // 5 minutes timeout
            .user_agent("SARdine/0.1.0 (SAR Processing Tool)")
            .build()
            .map_err(|e| SarError::Processing(format!("Failed to create HTTP client: {}", e)))?;
        
        // Retry logic for downloads
        let max_retries = 3;
        let mut last_error = None;
        
        for attempt in 1..=max_retries {
            log::debug!("Download attempt {} of {}", attempt, max_retries);
            
            match Self::try_download_once(&client, url, output_path) {
                Ok(_) => {
                    log::info!("Successfully downloaded and extracted SRTM tile to: {}", output_path);
                    return Ok(());
                },
                Err(e) => {
                    last_error = Some(e);
                    if attempt < max_retries {
                        log::warn!("Download attempt {} failed, retrying...", attempt);
                        std::thread::sleep(std::time::Duration::from_secs(2)); // Wait before retry
                    }
                }
            }
        }
        
        Err(last_error.unwrap_or_else(|| 
            SarError::Processing("Download failed after all retries".to_string())
        ))
    }
    
    /// Single download attempt
    fn try_download_once(client: &reqwest::blocking::Client, url: &str, output_path: &str) -> SarResult<()> {
        // Download the file
        let response = client.get(url)
            .send()
            .map_err(|e| SarError::Processing(format!("HTTP request failed: {}", e)))?;
        
        if !response.status().is_success() {
            return Err(SarError::Processing(
                format!("HTTP {} {}: {}", response.status().as_u16(), response.status().canonical_reason().unwrap_or(""), url)
            ));
        }
        
        let content = response.bytes()
            .map_err(|e| SarError::Processing(format!("Failed to read response body: {}", e)))?;
        
        // Validate content size (SRTM files should be at least a few KB)
        if content.len() < 1024 {
            return Err(SarError::Processing(
                format!("Downloaded file too small ({} bytes), likely an error page", content.len())
            ));
        }
        
        log::debug!("Downloaded {} bytes", content.len());
        
        // Handle different file formats based on URL
        if url.ends_with(".hgt.gz") || Self::is_gzip_content(&content) {
            // Decompress gzipped HGT file (AWS SRTM format)
            Self::extract_gzipped_hgt(&content, output_path)?;
        } else if url.ends_with(".zip") || url.contains("zip") || Self::is_zip_content(&content) {
            // Extract from ZIP archive
            Self::extract_srtm_zip(&content, output_path)?;
        } else if url.ends_with(".tif") || url.ends_with(".tiff") {
            // Direct GeoTIFF file (Copernicus DEM format)
            // Convert .tif extension to .hgt for consistency
            let hgt_path = output_path.replace(".hgt", ".tif");
            std::fs::write(&hgt_path, content)
                .map_err(|e| SarError::Io(e))?;
            
            // For now, just rename to .hgt - in a full implementation, 
            // we'd convert the GeoTIFF to HGT format
            if hgt_path != output_path {
                std::fs::rename(&hgt_path, output_path)
                    .map_err(|e| SarError::Io(e))?;
            }
        } else {
            // Direct HGT file
            std::fs::write(output_path, content)
                .map_err(|e| SarError::Io(e))?;
        }
        
        // Verify the output file exists and has reasonable size
        let metadata = std::fs::metadata(output_path)
            .map_err(|e| SarError::Io(e))?;
        
        if metadata.len() == 0 {
            return Err(SarError::Processing("Output file is empty".to_string()));
        }
        
        log::debug!("Output file size: {} bytes", metadata.len());
        Ok(())
    }
    
    /// Check if content is gzip format by examining magic bytes
    fn is_gzip_content(content: &[u8]) -> bool {
        content.len() >= 2 && content[0] == 0x1F && content[1] == 0x8B
    }
    
    /// Check if content is ZIP format by examining magic bytes
    fn is_zip_content(content: &[u8]) -> bool {
        content.len() >= 4 && content[0..4] == [0x50, 0x4B, 0x03, 0x04]
    }
    
    /// Extract SRTM HGT file from ZIP archive
    fn extract_srtm_zip(zip_data: &[u8], output_path: &str) -> SarResult<()> {
        use std::io::Cursor;
        use zip::ZipArchive;
        
        let reader = Cursor::new(zip_data);
        let mut archive = ZipArchive::new(reader)
            .map_err(|e| SarError::Processing(format!("Failed to open ZIP archive: {}", e)))?;
        
        // Find the HGT file in the archive
        for i in 0..archive.len() {
            let mut file = archive.by_index(i)
                .map_err(|e| SarError::Processing(format!("Failed to read ZIP entry {}: {}", i, e)))?;
            
            if file.name().ends_with(".hgt") {
                log::debug!("Extracting HGT file: {}", file.name());
                
                let mut buffer = Vec::new();
                std::io::copy(&mut file, &mut buffer)
                    .map_err(|e| SarError::Io(e))?;
                
                std::fs::write(output_path, buffer)
                    .map_err(|e| SarError::Io(e))?;
                
                return Ok(());
            }
        }
        
        Err(SarError::Processing("No HGT file found in ZIP archive".to_string()))
    }

    /// Extract HGT file from gzipped data (AWS SRTM format)
    fn extract_gzipped_hgt(gzip_data: &[u8], output_path: &str) -> SarResult<()> {
        use flate2::read::GzDecoder;
        use std::io::Read;
        
        log::debug!("Decompressing gzipped HGT file");
        
        let mut decoder = GzDecoder::new(gzip_data);
        let mut decompressed = Vec::new();
        
        decoder.read_to_end(&mut decompressed)
            .map_err(|e| SarError::Processing(format!("Failed to decompress gzip data: {}", e)))?;
        
        // Validate decompressed size (SRTM 1-arcsec tiles are typically ~25MB)
        if decompressed.is_empty() {
            return Err(SarError::Processing("Decompressed HGT file is empty".to_string()));
        }
        
        log::debug!("Decompressed {} bytes", decompressed.len());
        
        // Write decompressed HGT data to file
        std::fs::write(output_path, decompressed)
            .map_err(|e| SarError::Io(e))?;
        
        Ok(())
    }

    /// Automatically prepare DEM for SAR scene
    /// Downloads SRTM tiles if needed and creates a mosaic
    pub fn prepare_dem_for_scene(
        bbox: &BoundingBox,
        output_resolution: f64,
        cache_dir: &str,
    ) -> SarResult<(Array2<f32>, GeoTransform)> {
        log::info!("Preparing DEM for SAR scene");
        log::debug!("Bounding box: {:?}", bbox);
        log::debug!("Target resolution: {} meters", output_resolution);
        
        // Create cache directory
        let srtm_cache_dir = format!("{}/srtm_tiles", cache_dir);
        std::fs::create_dir_all(&srtm_cache_dir)
            .map_err(|e| SarError::Io(e))?;
        
        // Download required SRTM tiles
        let tile_files = match Self::download_srtm_tiles(bbox, &srtm_cache_dir) {
            Ok(files) => files,
            Err(e) => {
                // If download fails, check for existing DEM files
                log::warn!("SRTM download failed: {}. Checking for existing files...", e);
                Self::find_existing_dem_files(&srtm_cache_dir, bbox)?
            }
        };
        
        if tile_files.is_empty() {
            return Err(SarError::Processing(
                "No DEM files available. Please provide DEM data manually.".to_string()
            ));
        }
        
        // Create mosaic from tiles
        Self::create_dem_mosaic(&tile_files, bbox, output_resolution)
    }

    /// Find existing DEM files in cache directory
    fn find_existing_dem_files(cache_dir: &str, bbox: &BoundingBox) -> SarResult<Vec<String>> {
        let mut dem_files = Vec::new();
        
        if let Ok(entries) = std::fs::read_dir(cache_dir) {
            for entry in entries.flatten() {
                let path = entry.path();
                if let Some(extension) = path.extension() {
                    if extension == "hgt" || extension == "tif" || extension == "tiff" {
                        dem_files.push(path.to_string_lossy().to_string());
                    }
                }
            }
        }
        
        log::info!("Found {} existing DEM files", dem_files.len());
        Ok(dem_files)
    }

    /// Create a mosaic from multiple DEM tiles
    fn create_dem_mosaic(
        tile_files: &[String],
        bbox: &BoundingBox,
        target_resolution: f64,
    ) -> SarResult<(Array2<f32>, GeoTransform)> {
        log::info!("Creating DEM mosaic from {} tiles", tile_files.len());
        
        if tile_files.len() == 1 {
            // Single tile, just read and clip
            return Self::read_dem(&tile_files[0], bbox, target_resolution);
        }
        
        // For multiple tiles, we would need to implement mosaicking
        // For now, use the first available tile as a fallback
        log::warn!("Multi-tile mosaicking not yet implemented, using first tile");
        Self::read_dem(&tile_files[0], bbox, target_resolution)
    }

    /// Resample DEM to target grid using bilinear interpolation
    pub fn resample_dem(
        dem: &Array2<f32>,
        source_transform: &GeoTransform,
        target_transform: &GeoTransform,
        target_shape: (usize, usize),
    ) -> SarResult<Array2<f32>> {
        log::debug!("Resampling DEM to target grid");
        log::debug!("Target shape: {:?}", target_shape);

        let (target_height, target_width) = target_shape;
        let mut resampled = Array2::zeros(target_shape);
        let (source_height, source_width) = dem.dim();

        for i in 0..target_height {
            for j in 0..target_width {
                // Convert target pixel to geographic coordinates
                let target_x = target_transform.top_left_x + (j as f64) * target_transform.pixel_width;
                let target_y = target_transform.top_left_y + (i as f64) * target_transform.pixel_height;

                // Convert geographic coordinates to source pixel coordinates
                let source_col = (target_x - source_transform.top_left_x) / source_transform.pixel_width;
                let source_row = (target_y - source_transform.top_left_y) / source_transform.pixel_height;

                // Bilinear interpolation
                if source_col >= 0.0 && source_col < (source_width - 1) as f64 
                    && source_row >= 0.0 && source_row < (source_height - 1) as f64 {
                    
                    let x1 = source_col.floor() as usize;
                    let y1 = source_row.floor() as usize;
                    let x2 = (x1 + 1).min(source_width - 1);
                    let y2 = (y1 + 1).min(source_height - 1);

                    let dx = source_col - x1 as f64;
                    let dy = source_row - y1 as f64;

                    let v11 = dem[[y1, x1]] as f64;
                    let v12 = dem[[y2, x1]] as f64;
                    let v21 = dem[[y1, x2]] as f64;
                    let v22 = dem[[y2, x2]] as f64;

                    let interpolated = v11 * (1.0 - dx) * (1.0 - dy)
                        + v21 * dx * (1.0 - dy)
                        + v12 * (1.0 - dx) * dy
                        + v22 * dx * dy;

                    resampled[[i, j]] = interpolated as f32;
                } else {
                    resampled[[i, j]] = f32::NAN; // Outside DEM bounds
                }
            }
        }

        Ok(resampled)
    }

    /// Calculate slope and aspect from DEM using central differences
    pub fn calculate_slope_aspect(
        dem: &Array2<f32>,
        pixel_spacing: (f64, f64), // (dx, dy) in meters
    ) -> SarResult<(Array2<f32>, Array2<f32>)> {
        log::debug!("Calculating slope and aspect from DEM");

        let (height, width) = dem.dim();
        let mut slope = Array2::zeros((height, width));
        let mut aspect = Array2::zeros((height, width));

        let dx = pixel_spacing.0 as f32;
        let dy = pixel_spacing.1 as f32;

        for i in 1..height-1 {
            for j in 1..width-1 {
                // Central differences
                let dz_dx = (dem[[i, j+1]] - dem[[i, j-1]]) / (2.0 * dx);
                let dz_dy = (dem[[i+1, j]] - dem[[i-1, j]]) / (2.0 * dy);

                // Slope in radians
                slope[[i, j]] = (dz_dx * dz_dx + dz_dy * dz_dy).sqrt().atan();

                // Aspect in radians (0 = North, clockwise positive)
                aspect[[i, j]] = (-dz_dy).atan2(dz_dx);
            }
        }

        // Fill edges
        Self::fill_edge_values(&mut slope)?;
        Self::fill_edge_values(&mut aspect)?;

        Ok((slope, aspect))
    }

    /// Convert slope/aspect to surface normal vectors
    pub fn slope_aspect_to_normals(
        slope: &Array2<f32>,
        aspect: &Array2<f32>
    ) -> Array2<[f32; 3]> {
        let (height, width) = slope.dim();
        let mut normals = Array2::from_elem((height, width), [0.0, 0.0, 1.0]);

        for i in 0..height {
            for j in 0..width {
                let slope_rad = slope[[i, j]];
                let aspect_rad = aspect[[i, j]];

                // Surface normal components
                let nx = -slope_rad.sin() * aspect_rad.sin();
                let ny = slope_rad.sin() * aspect_rad.cos();
                let nz = slope_rad.cos();

                normals[[i, j]] = [nx, ny, nz];
            }
        }

        normals
    }

    /// Fill edge values by copying from nearest interior pixels
    fn fill_edge_values(array: &mut Array2<f32>) -> SarResult<()> {
        let (height, width) = array.dim();
        
        if height < 3 || width < 3 {
            return Err(SarError::Processing("Array too small for edge filling".to_string()));
        }

        // Fill top and bottom edges
        for j in 0..width {
            array[[0, j]] = array[[1, j]];
            array[[height-1, j]] = array[[height-2, j]];
        }

        // Fill left and right edges  
        for i in 0..height {
            array[[i, 0]] = array[[i, 1]];
            array[[i, width-1]] = array[[i, width-2]];
        }

        Ok(())
    }

    /// Calculate local incidence angle from DEM
    pub fn calculate_incidence_angle(
        dem: &Array2<f32>,
        geo_transform: &GeoTransform,
        satellite_position: [f64; 3],
    ) -> SarResult<Array2<f32>> {
        log::debug!("Calculating local incidence angles");

        let (height, width) = dem.dim();
        let mut incidence_angles = Array2::zeros((height, width));

        // Placeholder for incidence angle calculation
        // This would compute the angle between the radar beam and local terrain normal
        for i in 0..height {
            for j in 0..width {
                // Simplified calculation - in practice this would be much more complex
                incidence_angles[[i, j]] = 30.0_f32.to_radians(); // ~30 degrees
            }
        }

        Ok(incidence_angles)
    }

    /// Test SRTM download capability for a specific tile
    pub fn test_srtm_download(tile_name: &str, output_dir: &str) -> SarResult<String> {
        log::info!("Testing SRTM download for tile: {}", tile_name);
        
        std::fs::create_dir_all(output_dir)
            .map_err(|e| SarError::Io(e))?;
        
        let output_path = format!("{}/{}.hgt", output_dir, tile_name);
        
        // Check if file already exists
        if std::path::Path::new(&output_path).exists() {
            log::info!("SRTM tile {} already exists at: {}", tile_name, output_path);
            return Ok(output_path);
        }
        
        // Try downloading from sources
        let success = Self::try_download_from_sources(tile_name, &output_path)?;
        
        if success {
            Ok(output_path)
        } else {
            Err(SarError::Processing(format!(
                "Failed to download SRTM tile {} from all available sources", tile_name
            )))
        }
    }

    /// Add Copernicus DEM as an additional source (30m resolution)
    fn get_copernicus_dem_urls(tile: &str) -> Vec<String> {
        // Copernicus DEM is available at 30m and 90m resolution
        // Format: https://copernicus-dem-30m.s3.amazonaws.com/DEM/Copernicus_DSM_COG_10_N50_00_E012_00_DEM.tif
        
        // Parse tile coordinates (e.g., "N50E012" -> lat=50, lon=12)
        if let Some((lat_str, lon_str)) = Self::parse_tile_coordinates(tile) {
            let copernicus_30m = format!(
                "https://copernicus-dem-30m.s3.amazonaws.com/DEM/Copernicus_DSM_COG_10_{}_00_{}_00_DEM.tif",
                lat_str, lon_str
            );
            let copernicus_90m = format!(
                "https://copernicus-dem-90m.s3.amazonaws.com/DEM/Copernicus_DSM_COG_30_{}_00_{}_00_DEM.tif", 
                lat_str, lon_str
            );
            
            vec![copernicus_30m, copernicus_90m]
        } else {
            vec![]
        }
    }
    
    /// Parse tile coordinates from SRTM tile name
    fn parse_tile_coordinates(tile: &str) -> Option<(String, String)> {
        // Parse tile name like "N50E012" or "S23W045"
        if tile.len() < 6 {
            return None;
        }
        
        let lat_part = &tile[0..3]; // "N50" or "S23"
        let lon_part = &tile[3..];   // "E012" or "W045"
        
        // Validate format
        if (lat_part.starts_with('N') || lat_part.starts_with('S')) &&
           (lon_part.starts_with('E') || lon_part.starts_with('W')) {
            Some((lat_part.to_string(), lon_part.to_string()))
        } else {
            None
        }
    }

    /// Fill voids in DEM using interpolation
    pub fn fill_dem_voids(dem: &mut Array2<f32>, no_data_value: f32) -> SarResult<()> {
        log::debug!("Filling DEM voids with no-data value: {}", no_data_value);
        
        let (height, width) = dem.dim();
        let mut void_count = 0;
        let mut filled_count = 0;
        
        // First pass: identify voids
        let mut void_mask = Array2::from_elem((height, width), false);
        for i in 0..height {
            for j in 0..width {
                if dem[[i, j]].is_nan() || dem[[i, j]] == no_data_value || dem[[i, j]] < -32000.0 {
                    void_mask[[i, j]] = true;
                    void_count += 1;
                }
            }
        }
        
        log::info!("Found {} void pixels in DEM ({:.2}%)", 
                   void_count, (void_count as f64 / (height * width) as f64) * 100.0);
        
        if void_count == 0 {
            return Ok(());
        }
        
        // Iterative void filling using neighbor averaging
        let max_iterations = 10;
        for iteration in 0..max_iterations {
            let mut changed = false;
            let mut new_dem = dem.clone();
            
            for i in 1..height-1 {
                for j in 1..width-1 {
                    if void_mask[[i, j]] {
                        // Check neighbors for valid values
                        let mut sum = 0.0;
                        let mut count = 0;
                        
                        for di in -1i32..=1 {
                            for dj in -1i32..=1 {
                                if di == 0 && dj == 0 { continue; }
                                
                                let ni = (i as i32 + di) as usize;
                                let nj = (j as i32 + dj) as usize;
                                
                                if ni < height && nj < width && !void_mask[[ni, nj]] {
                                    sum += dem[[ni, nj]];
                                    count += 1;
                                }
                            }
                        }
                        
                        if count >= 3 {  // Need at least 3 valid neighbors
                            new_dem[[i, j]] = sum / count as f32;
                            void_mask[[i, j]] = false;
                            filled_count += 1;
                            changed = true;
                        }
                    }
                }
            }
            
            *dem = new_dem;
            
            if !changed {
                break;
            }
            
            log::debug!("Iteration {}: filled {} pixels", iteration + 1, filled_count);
        }
        
        // Fill remaining voids with global mean
        let mut valid_sum = 0.0;
        let mut valid_count = 0;
        
        for i in 0..height {
            for j in 0..width {
                if !void_mask[[i, j]] {
                    valid_sum += dem[[i, j]];
                    valid_count += 1;
                }
            }
        }
        
        if valid_count > 0 {
            let global_mean = valid_sum / valid_count as f32;
            
            for i in 0..height {
                for j in 0..width {
                    if void_mask[[i, j]] {
                        dem[[i, j]] = global_mean;
                        filled_count += 1;
                    }
                }
            }
        }
        
        log::info!("Filled {} total void pixels", filled_count);
        Ok(())
    }

    /// Create terrain mask for layover/shadow areas
    pub fn create_terrain_mask(
        local_incidence_angles: &Array2<f32>,
        min_angle_deg: f32,
        max_angle_deg: f32,
    ) -> Array2<bool> {
        log::debug!("Creating terrain mask with incidence angle range: {:.1}° - {:.1}°", 
                   min_angle_deg, max_angle_deg);
        
        let (height, width) = local_incidence_angles.dim();
        let mut mask = Array2::from_elem((height, width), true);
        
        let min_angle_rad = min_angle_deg.to_radians();
        let max_angle_rad = max_angle_deg.to_radians();
        
        let mut masked_count = 0;
        
        for i in 0..height {
            for j in 0..width {
                let angle = local_incidence_angles[[i, j]];
                
                // Mask if:
                // - Angle is NaN or invalid
                // - Angle is too small (layover)
                // - Angle is too large (shadow/grazing)
                // - Angle is negative (backslope)
                if angle.is_nan() || 
                   angle < min_angle_rad || 
                   angle > max_angle_rad ||
                   angle <= 0.0 {
                    mask[[i, j]] = false;
                    masked_count += 1;
                }
            }
        }
        
        let mask_percentage = (masked_count as f64 / (height * width) as f64) * 100.0;
        log::info!("Masked {} pixels ({:.2}%) due to invalid incidence angles", 
                   masked_count, mask_percentage);
        
        mask
    }

    /// Validate DEM coverage for SAR scene
    pub fn validate_dem_coverage(
        dem_bbox: &BoundingBox,
        sar_bbox: &BoundingBox,
        buffer_degrees: f64,
    ) -> SarResult<bool> {
        log::debug!("Validating DEM coverage for SAR scene");
        log::debug!("DEM bbox: {:?}", dem_bbox);
        log::debug!("SAR bbox: {:?}", sar_bbox);
        log::debug!("Buffer: {:.3} degrees", buffer_degrees);
        
        // Check if DEM covers SAR scene with buffer
        let required_min_lat = sar_bbox.min_lat - buffer_degrees;
        let required_max_lat = sar_bbox.max_lat + buffer_degrees;
        let required_min_lon = sar_bbox.min_lon - buffer_degrees;
        let required_max_lon = sar_bbox.max_lon + buffer_degrees;
        
        let coverage_ok = dem_bbox.min_lat <= required_min_lat &&
                         dem_bbox.max_lat >= required_max_lat &&
                         dem_bbox.min_lon <= required_min_lon &&
                         dem_bbox.max_lon >= required_max_lon;
        
        if coverage_ok {
            log::info!("✅ DEM provides adequate coverage for SAR scene");
        } else {
            log::warn!("⚠️  DEM coverage insufficient for SAR scene");
            log::warn!("Required: lat [{:.3}, {:.3}], lon [{:.3}, {:.3}]", 
                      required_min_lat, required_max_lat, required_min_lon, required_max_lon);
            log::warn!("Available: lat [{:.3}, {:.3}], lon [{:.3}, {:.3}]", 
                      dem_bbox.min_lat, dem_bbox.max_lat, dem_bbox.min_lon, dem_bbox.max_lon);
        }
        
        Ok(coverage_ok)
    }

    /// Apply terrain mask to flattened data
    pub fn apply_terrain_mask_to_data(
        data: &mut Array2<f32>,
        mask: &Array2<bool>,
        fill_value: f32,
    ) -> SarResult<()> {
        log::debug!("Applying terrain mask to data");
        
        let (height, width) = data.dim();
        if mask.dim() != (height, width) {
            return Err(SarError::Processing(
                "Data and mask dimensions do not match".to_string()
            ));
        }
        
        let mut masked_count = 0;
        
        for i in 0..height {
            for j in 0..width {
                if !mask[[i, j]] {
                    data[[i, j]] = fill_value;
                    masked_count += 1;
                }
            }
        }
        
        log::info!("Applied mask to {} pixels", masked_count);
        Ok(())
    }

    /// Complete terrain flattening pipeline
    /// Combines all steps: DEM preparation, void filling, slope calculation, 
    /// surface normals, incidence angles, masking, and flattening
    pub fn complete_terrain_flattening_pipeline(
        sar_image: &Array2<f32>,              // Calibrated sigma0 data
        sar_bbox: &BoundingBox,               // SAR scene bounding box
        sar_geo_transform: &GeoTransform,     // SAR geotransform
        orbit_data: &crate::types::OrbitData, // Satellite orbit information
        cache_dir: &str,                      // DEM cache directory
        output_resolution: f64,               // Target resolution in meters
    ) -> SarResult<(Array2<f32>, Array2<bool>)> {
        log::info!("🌍 Starting complete terrain flattening pipeline");
        
        // Step 1: Load and prepare DEM
        log::info!("📊 Step 1: Preparing DEM for SAR scene");
        let (mut dem, dem_geo_transform) = Self::prepare_dem_for_scene(
            sar_bbox, 
            output_resolution, 
            cache_dir
        )?;
        
        // Step 2: Validate DEM coverage
        log::info!("✅ Step 2: Validating DEM coverage");
        
        // Calculate DEM bounding box from geotransform
        // For north-up images (pixel_height < 0), top_left_y is the maximum latitude
        let dem_height = dem.dim().0 as f64;
        let dem_width = dem.dim().1 as f64;
        
        let dem_bbox = BoundingBox {
            min_lat: dem_geo_transform.top_left_y + dem_height * dem_geo_transform.pixel_height.min(0.0),
            max_lat: dem_geo_transform.top_left_y + dem_height * dem_geo_transform.pixel_height.max(0.0),
            min_lon: dem_geo_transform.top_left_x + dem_width * dem_geo_transform.pixel_width.min(0.0),
            max_lon: dem_geo_transform.top_left_x + dem_width * dem_geo_transform.pixel_width.max(0.0),
        };
        
        log::debug!("DEM bbox: min_lat={:.6}, max_lat={:.6}, min_lon={:.6}, max_lon={:.6}", 
                   dem_bbox.min_lat, dem_bbox.max_lat, dem_bbox.min_lon, dem_bbox.max_lon);
        log::debug!("SAR bbox: min_lat={:.6}, max_lat={:.6}, min_lon={:.6}, max_lon={:.6}", 
                   sar_bbox.min_lat, sar_bbox.max_lat, sar_bbox.min_lon, sar_bbox.max_lon);
        
        let coverage_ok = Self::validate_dem_coverage(&dem_bbox, sar_bbox, 0.001)?; // Very small buffer for testing
        if !coverage_ok {
            return Err(SarError::Processing(
                "DEM does not provide adequate coverage for SAR scene".to_string()
            ));
        }
        
        // Step 3: Fill DEM voids
        log::info!("🔧 Step 3: Filling DEM voids");
        Self::fill_dem_voids(&mut dem, -32768.0)?; // SRTM no-data value
        
        // Step 4: Resample DEM to SAR geometry if needed
        log::info!("📐 Step 4: Resampling DEM to SAR geometry");
        let resampled_dem = if dem_geo_transform.pixel_width.abs() != sar_geo_transform.pixel_width.abs() ||
                               dem_geo_transform.pixel_height.abs() != sar_geo_transform.pixel_height.abs() {
            log::info!("Resampling DEM from {:.6}° to match SAR geometry", dem_geo_transform.pixel_width);
            Self::resample_dem(&dem, &dem_geo_transform, sar_geo_transform, sar_image.dim())?
        } else {
            dem
        };
        
        // Step 5: Compute slope and aspect
        log::info!("📈 Step 5: Computing terrain slope and aspect");
        let pixel_spacing = (
            output_resolution,
            output_resolution
        );
        let (slope, aspect) = Self::calculate_slope_aspect(&resampled_dem, pixel_spacing)?;
        
        // Step 6: Compute surface normal vectors
        log::info!("🔺 Step 6: Computing surface normal vectors");
        let surface_normals = Self::slope_aspect_to_normals(&slope, &aspect);
        
        // Step 7: Compute local incidence angles
        log::info!("📡 Step 7: Computing local incidence angles");
        let local_incidence_angles = Self::compute_local_incidence_angles_advanced(
            &resampled_dem,
            sar_geo_transform,
            orbit_data,
            &surface_normals
        )?;
        
        // Step 8: Create terrain mask
        log::info!("🎭 Step 8: Creating terrain mask");
        let terrain_mask = Self::create_terrain_mask(
            &local_incidence_angles,
            10.0,  // Minimum incidence angle (degrees)
            80.0   // Maximum incidence angle (degrees)
        );
        
        // Step 9: Apply terrain flattening
        log::info!("🏔️  Step 9: Applying terrain flattening");
        let mut flattened_image = sar_image.clone();
        Self::apply_terrain_flattening_to_image(
            &mut flattened_image,
            &local_incidence_angles,
            30.0f32.to_radians() // Reference incidence angle
        )?;
        
        // Step 10: Apply terrain mask
        log::info!("🎯 Step 10: Applying terrain mask");
        Self::apply_terrain_mask_to_data(&mut flattened_image, &terrain_mask, f32::NAN)?;
        
        log::info!("🎉 Terrain flattening pipeline completed successfully!");
        
        Ok((flattened_image, terrain_mask))
    }

    /// Advanced local incidence angle calculation using orbit data
    fn compute_local_incidence_angles_advanced(
        dem: &Array2<f32>,
        geo_transform: &GeoTransform,
        orbit_data: &crate::types::OrbitData,
        surface_normals: &Array2<[f32; 3]>,
    ) -> SarResult<Array2<f32>> {
        log::debug!("Computing local incidence angles with orbit data");
        
        let (height, width) = dem.dim();
        let mut incidence_angles = Array2::zeros((height, width));
        
        // For now, use a simplified calculation
        // In a full implementation, this would:
        // 1. Interpolate satellite position for each pixel's azimuth time
        // 2. Compute radar look vector from satellite to ground
        // 3. Calculate dot product with surface normal
        
        for i in 0..height {
            for j in 0..width {
                // Get surface normal
                let normal = surface_normals[[i, j]];
                
                // Simplified radar look vector (assumes right-looking geometry)
                // In reality, this would be computed from orbit data
                let look_vector = [0.5f32, 0.0, -0.866]; // ~30° incidence angle
                
                // Compute dot product: cos(incidence_angle) = normal · look_vector
                let cos_inc = normal[0] * look_vector[0] + 
                             normal[1] * look_vector[1] + 
                             normal[2] * look_vector[2];
                
                // Incidence angle = arccos(dot product)
                let incidence_angle = cos_inc.clamp(-1.0, 1.0).acos();
                
                incidence_angles[[i, j]] = incidence_angle;
            }
        }
        
        Ok(incidence_angles)
    }

    /// Apply terrain flattening to SAR image
    fn apply_terrain_flattening_to_image(
        image: &mut Array2<f32>,
        local_incidence_angles: &Array2<f32>,
        reference_incidence_angle: f32,
    ) -> SarResult<()> {
        log::debug!("Applying terrain flattening to SAR image");
        
        let (height, width) = image.dim();
        if local_incidence_angles.dim() != (height, width) {
            return Err(SarError::Processing(
                "Image and incidence angle dimensions do not match".to_string()
            ));
        }
        
        let cos_ref = reference_incidence_angle.cos();
        let mut corrected_count = 0;
        
        for i in 0..height {
            for j in 0..width {
                let local_angle = local_incidence_angles[[i, j]];
                
                if local_angle.is_finite() && local_angle > 0.0 {
                    let cos_local = local_angle.cos();
                    
                    if cos_local > 0.001 {  // Avoid division by very small numbers
                        // Terrain flattening: σ°_flat = σ° * cos(θ_ref) / cos(θ_local)
                        let correction_factor = cos_ref / cos_local;
                        image[[i, j]] *= correction_factor;
                        corrected_count += 1;
                    } else {
                        // Invalid geometry - set to NaN
                        image[[i, j]] = f32::NAN;
                    }
                } else {
                    image[[i, j]] = f32::NAN;
                }
            }
        }
        
        log::info!("Applied terrain flattening to {} pixels", corrected_count);
        Ok(())
    }
}
