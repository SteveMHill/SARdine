use sardine::io::orbit::OrbitReader;
use tempfile::NamedTempFile;
use std::io::Write;

#[test]
fn test_xml_orbit_format_parsing() {
    // Initialize logging
    env_logger::init();
    
    println!("=== Testing XML Orbit Format Parsing ===");
    
    // Create a sample orbit file in XML format
    let sample_xml = r#"<?xml version="1.0" encoding="UTF-8"?>
<Earth_Explorer_File>
  <Earth_Fixed_Coordinate_System>ITRF_2008</Earth_Fixed_Coordinate_System>
  <Data_Block type="xml">
    <List_of_OSVs count="3">
      <OSV>
        <UTC>UTC=2020-01-03T17:00:00.000000</UTC>
        <TAI>TAI=2020-01-03T17:00:37.000000</TAI>
        <X unit="m">X=1234567.890</X>
        <Y unit="m">Y=2345678.901</Y>
        <Z unit="m">Z=3456789.012</Z>
        <VX unit="m/s">VX=1234.5678</VX>
        <VY unit="m/s">VY=2345.6789</VY>
        <VZ unit="m/s">VZ=3456.7890</VZ>
      </OSV>
      <OSV>
        <UTC>UTC=2020-01-03T17:01:00.000000</UTC>
        <TAI>TAI=2020-01-03T17:01:37.000000</TAI>
        <X unit="m">X=1234567.890</X>
        <Y unit="m">Y=2345678.901</Y>
        <Z unit="m">Z=3456789.012</Z>
        <VX unit="m/s">VX=1234.5678</VX>
        <VY unit="m/s">VY=2345.6789</VY>
        <VZ unit="m/s">VZ=3456.7890</VZ>
      </OSV>
      <OSV>
        <UTC>UTC=2020-01-03T17:02:00.000000</UTC>
        <TAI>TAI=2020-01-03T17:02:37.000000</TAI>
        <X unit="m">X=1234567.890</X>
        <Y unit="m">Y=2345678.901</Y>
        <Z unit="m">Z=3456789.012</Z>
        <VX unit="m/s">VX=1234.5678</VX>
        <VY unit="m/s">VY=2345.6789</VY>
        <VZ unit="m/s">VZ=3456.7890</VZ>
      </OSV>
    </List_of_OSVs>
  </Data_Block>
</Earth_Explorer_File>"#;
    
    // Write to temporary file
    let mut temp_file = NamedTempFile::new().expect("Failed to create temp file");
    temp_file.write_all(sample_xml.as_bytes()).expect("Failed to write sample data");
    
    // Test parsing
    match OrbitReader::read_orbit_file(temp_file.path()) {
        Ok(orbit_data) => {
            println!("✅ Successfully parsed XML orbit file!");
            println!("  - State vectors: {}", orbit_data.state_vectors.len());
            
            assert_eq!(orbit_data.state_vectors.len(), 3, "Should have parsed 3 state vectors");
            
            let first_sv = &orbit_data.state_vectors[0];
            println!("  - First state vector time: {}", first_sv.time);
            println!("  - Position: [{:.3}, {:.3}, {:.3}] m", 
                    first_sv.position[0], first_sv.position[1], first_sv.position[2]);
            println!("  - Velocity: [{:.6}, {:.6}, {:.6}] m/s", 
                    first_sv.velocity[0], first_sv.velocity[1], first_sv.velocity[2]);
            
            // Verify the values
            assert!((first_sv.position[0] - 1234567.890).abs() < 1e-3);
            assert!((first_sv.position[1] - 2345678.901).abs() < 1e-3);
            assert!((first_sv.position[2] - 3456789.012).abs() < 1e-3);
            assert!((first_sv.velocity[0] - 1234.5678).abs() < 1e-6);
            assert!((first_sv.velocity[1] - 2345.6789).abs() < 1e-6);
            assert!((first_sv.velocity[2] - 3456.7890).abs() < 1e-6);
            
            println!("✅ All values parsed correctly!");
        },
        Err(e) => {
            panic!("❌ Failed to parse XML orbit file: {}", e);
        }
    }
}
