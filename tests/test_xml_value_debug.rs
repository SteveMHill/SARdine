use sardine::io::orbit::OrbitReader;

#[test]
fn test_xml_value_extraction() {
    // Test the XML value extraction directly
    println!("=== Testing XML Value Extraction ===");
    
    let test_lines = vec![
        r#"<UTC>UTC=2020-01-03T17:00:00.000000</UTC>"#,
        r#"<X unit="m">X=1234567.890</X>"#,
        r#"<VX unit="m/s">VX=1234.5678</VX>"#,
    ];
    
    for line in test_lines {
        println!("Testing line: {}", line);
        
        // Try to extract UTC value
        if line.contains("<UTC>") {
            if let Some(value) = extract_xml_value_test(line, "UTC") {
                println!("  UTC value: '{}'", value);
            }
        }
        
        // Try to extract X value
        if line.contains("<X ") {
            if let Some(value) = extract_xml_value_test(line, "X") {
                println!("  X value: '{}'", value);
            }
        }
        
        // Try to extract VX value
        if line.contains("<VX ") {
            if let Some(value) = extract_xml_value_test(line, "VX") {
                println!("  VX value: '{}'", value);
            }
        }
    }
}

/// Copy of the extract_xml_value function for testing
fn extract_xml_value_test(line: &str, tag: &str) -> Option<String> {
    // Look for opening tag (with or without attributes)
    let start_pattern = format!("<{}", tag);
    let end_tag = format!("</{}>", tag);
    
    println!("    Looking for start_pattern: '{}' and end_tag: '{}'", start_pattern, end_tag);
    
    if let Some(start) = line.find(&start_pattern) {
        println!("    Found start pattern at position {}", start);
        // Find the end of the opening tag
        if let Some(tag_end) = line[start..].find('>') {
            let content_start = start + tag_end + 1;
            println!("    Found opening tag end at relative position {}, content starts at {}", tag_end, content_start);
            if let Some(end) = line[content_start..].find(&end_tag) {
                let result = line[content_start..content_start + end].to_string();
                println!("    Found end tag at relative position {}, extracted: '{}'", end, result);
                return Some(result);
            } else {
                println!("    End tag not found in remaining string: '{}'", &line[content_start..]);
            }
        } else {
            println!("    Opening tag end '>' not found");
        }
    } else {
        println!("    Start pattern not found in line");
    }
    None
}
