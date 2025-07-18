# Terrain Correction (Geocoding) Implementation

## Overview

Terrain correction, also known as geocoding or orthorectification, transforms SAR data from slant range geometry to map-projected coordinates. This process corrects for topographic distortions and provides geographically accurate SAR products.

## 🎯 What is Terrain Correction?

SAR sensors acquire data in **slant range geometry**, where distances are measured along the radar beam path. This geometry contains distortions due to:

- **Foreshortening**: Steep slopes appear compressed
- **Layover**: Mountaintops appear in front of their bases
- **Shadow**: Areas hidden behind terrain features
- **Range migration**: Curved Earth and orbital effects

Terrain correction transforms this data into **map coordinates** (UTM, WGS84, etc.) for:
- Geographic accuracy
- GIS integration
- Multi-temporal analysis
- Scientific measurements

## 🔄 Implementation Approach

### Range-Doppler Terrain Correction

SARdine implements Range-Doppler terrain correction, the same method used by ESA SNAP:

1. **Range Equation**: `R = c·t/2`
   - Relates slant range to travel time
   - Accounts for platform motion

2. **Doppler Equation**: `f_d = (2·v·cos(θ))/λ`
   - Relates azimuth position to Doppler frequency
   - Determines along-track location

3. **DEM Intersection**: 
   - Solves for ground point elevation
   - Corrects topographic distortions

## 🏗️ Core Components

### TerrainCorrector

Main processor for terrain correction operations:

```rust
pub struct TerrainCorrector {
    dem: Array2<f32>,           // Digital elevation model
    dem_transform: GeoTransform, // DEM spatial reference
    dem_nodata: f32,            // No-data fill value
    output_crs: u32,            // Target coordinate system (EPSG)
    output_spacing: f64,        // Output pixel spacing (meters)
}
```

### RangeDopplerParams

Radar system parameters for geometric calculations:

```rust
pub struct RangeDopplerParams {
    range_pixel_spacing: f64,    // Range resolution (meters)
    azimuth_pixel_spacing: f64,  // Azimuth resolution (meters)
    slant_range_time: f64,       // Time to first pixel (seconds)
    prf: f64,                    // Pulse repetition frequency (Hz)
    wavelength: f64,             // Radar wavelength (meters)
    speed_of_light: f64,         // Light speed constant
}
```

## 📐 Processing Algorithm

### 1. Backward Geocoding

For each output map pixel:

```
map_coordinates → geographic_coordinates → elevation → sar_pixel
```

### 2. Coordinate Transformations

#### Geographic to ECEF
```rust
fn latlon_to_ecef(lat: f64, lon: f64, elevation: f64) -> [f64; 3] {
    let a = 6_378_137.0; // WGS84 semi-major axis
    let e2 = 0.00669437999014; // First eccentricity squared
    
    let lat_rad = lat.to_radians();
    let lon_rad = lon.to_radians();
    let n = a / (1.0 - e2 * lat_rad.sin().powi(2)).sqrt();
    
    [
        (n + elevation) * lat_rad.cos() * lon_rad.cos(),
        (n + elevation) * lat_rad.cos() * lon_rad.sin(),
        (n * (1.0 - e2) + elevation) * lat_rad.sin(),
    ]
}
```

#### SAR Geometry Calculation
```rust
fn latlon_to_sar_pixel(
    lat: f64, lon: f64, elevation: f64,
    orbit_data: &OrbitData,
    params: &RangeDopplerParams,
) -> Option<(usize, usize)> {
    // Convert to ECEF coordinates
    let ground_point = latlon_to_ecef(lat, lon, elevation);
    
    // Find satellite closest approach time
    let azimuth_time = find_azimuth_time(&ground_point, orbit_data);
    
    // Interpolate satellite state
    let sat_state = interpolate_satellite_state(orbit_data, azimuth_time);
    
    // Calculate slant range
    let range_vector = [
        ground_point[0] - sat_state.position[0],
        ground_point[1] - sat_state.position[1],
        ground_point[2] - sat_state.position[2],
    ];
    let slant_range = range_vector.magnitude();
    
    // Convert to pixel coordinates
    let range_pixel = slant_range_to_pixel(slant_range, params);
    let azimuth_pixel = azimuth_time_to_pixel(azimuth_time, params);
    
    Some((range_pixel, azimuth_pixel))
}
```

### 3. DEM Operations

#### Elevation Interpolation
```rust
fn get_elevation_at_latlon(lat: f64, lon: f64) -> Option<f64> {
    // Convert to DEM pixel coordinates
    let col = (lon - dem_transform.top_left_x) / dem_transform.pixel_width;
    let row = (lat - dem_transform.top_left_y) / dem_transform.pixel_height;
    
    // Bilinear interpolation
    bilinear_interpolate(&dem, col, row)
}
```

### 4. Output Grid Creation

```rust
fn create_output_grid(bounds: &BoundingBox) -> (usize, usize, GeoTransform) {
    // Calculate pixel size in target CRS
    let pixel_size_lat = output_spacing / meters_per_degree_lat;
    let pixel_size_lon = output_spacing / meters_per_degree_lon;
    
    // Grid dimensions
    let width = ((bounds.max_lon - bounds.min_lon) / pixel_size_lon).ceil();
    let height = ((bounds.max_lat - bounds.min_lat) / pixel_size_lat).ceil();
    
    // Geotransform
    let transform = GeoTransform {
        top_left_x: bounds.min_lon,
        pixel_width: pixel_size_lon,
        top_left_y: bounds.max_lat,
        pixel_height: -pixel_size_lat, // North-up
        // ... other fields
    };
    
    (width, height, transform)
}
```

## 🐍 Python API

### Complete Pipeline
```python
import sardine

# Perform terrain correction
sardine.terrain_correction(
    sar_image=[[...]], # 2D list of SAR values
    dem_path="/path/to/dem.tif",
    orbit_data=[(time, position, velocity), ...],
    sar_bbox=(min_lon, min_lat, max_lon, max_lat),
    output_path="/path/to/output.tif",
    output_crs=32610,  # UTM Zone 10N
    output_spacing=10.0 # 10m pixels
)
```

### Coordinate Utilities
```python
# Convert coordinates
ecef = sardine.latlon_to_ecef(lat=37.7749, lon=-122.4194, elevation=50.0)
print(f"ECEF: {ecef}")  # [X, Y, Z] in meters
```

## 🖥️ CLI Interface

### Basic Geocoding
```bash
# Perform terrain correction
sardine geocode sar_image.npy dem.tif orbit.EOF output.tif \
  --bbox "37.0,38.0,-122.5,-121.5" \
  --output-crs 4326 \
  --output-spacing 10.0
```

### DEM Testing
```bash
# Test DEM loading
sardine test-dem dem.tif --output-crs 4326 --output-spacing 10.0
```

## 📊 Supported Formats

### Input Data
- **SAR Images**: NumPy arrays (.npy), GeoTIFF
- **DEMs**: GeoTIFF (SRTM, Copernicus, custom)
- **Orbits**: Sentinel-1 precise orbit files (.EOF)

### Output Formats
- **GeoTIFF**: Georeferenced raster with CRS information
- **Compression**: LZW, DEFLATE, or uncompressed
- **Data Types**: Float32 for backscatter values

## 🎯 Coordinate Reference Systems

### Supported CRS
- **EPSG:4326**: WGS84 Geographic (lat/lon)
- **UTM Zones**: Local projected coordinates
- **Custom**: Any EPSG-defined CRS

### CRS Selection Guidelines
- **Global analysis**: EPSG:4326 (WGS84)
- **Local studies**: UTM zones for area
- **Polar regions**: Polar stereographic projections

## ⚡ Performance Optimizations

### Backward Geocoding
- Processes output pixels → SAR pixels
- Reduces geometric distortions
- Better sampling control

### Bilinear Interpolation
- Smooth value transitions
- Preserves radiometric accuracy
- Efficient implementation

### Orbit Interpolation
- Linear interpolation between state vectors
- Cached satellite positions
- Optimized distance calculations

## 🧪 Testing and Validation

### Synthetic Data Tests
```python
# Create test data
from examples.complete_terrain_correction_workflow import *

# Run workflow
test_dir = perform_terrain_correction_workflow()
print(f"Test data in: {test_dir}")
```

### Real Data Validation
1. Use known ground control points
2. Compare with SNAP geocoding results
3. Validate geometric accuracy
4. Check radiometric preservation

## 🚀 Example Workflow

### Complete Processing Chain
```python
import sardine

# Load Sentinel-1 product
reader = sardine.SlcReader("S1A_IW_SLC__*.zip")

# Get product information
info = reader.get_product_info()
bbox = info.bounding_box

# Download DEM and orbit
dem_files = sardine.download_srtm_tiles(bbox, "./dem_cache")
orbit_data = sardine.load_orbit_file("precise_orbit.EOF")

# Process SAR data
calibrated = reader.calibrate_and_multilook("VV", "sigma0", 4, 1)

# Perform terrain correction
geocoded = sardine.terrain_correction(
    calibrated, dem_files[0], orbit_data, bbox,
    "geocoded_sigma0.tif", 32610, 10.0
)
```

## 🔧 Configuration Options

### Quality vs Speed Trade-offs
- **High Quality**: Small output pixels, bilinear interpolation
- **Fast Processing**: Larger pixels, nearest neighbor
- **Memory Efficient**: Process in chunks

### Output Customization
- **Pixel spacing**: 1m-1000m depending on application
- **Coordinate system**: Match local standards
- **Compression**: Balance file size vs access speed

## 📚 References

1. **ESA Sentinel-1 Toolbox**: SNAP Range-Doppler Terrain Correction
2. **Zebker, H.A. & Goldstein, R.M.**: Topographic mapping from SAR observations
3. **Bamler, R. & Hartl, P.**: Synthetic aperture radar interferometry
4. **Small, D.**: Flattening gamma SAR imagery

## ✅ Summary

SARdine's terrain correction implementation provides:

- **Accurate geocoding** using Range-Doppler equations
- **Automatic DEM integration** with SRTM/Copernicus data
- **Flexible output projections** for any EPSG coordinate system
- **High-performance processing** with Rust backend
- **Complete Python/CLI APIs** for integration

The implementation follows established SAR processing principles while providing modern performance and usability improvements over traditional tools like SNAP.

---

*For implementation details, see the source code in `src/core/terrain_correction.rs`*
