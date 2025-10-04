/// Type-safe coordinate frames to prevent index confusion and violations
/// 
/// This module implements the "Coordinate frames" refactor from the practical refactors:
/// Newtypes/structs for BurstIdx, SubswathIdx, StitchedIdx; explicit converters.

use std::fmt;
use crate::types::SarResult;

// ============================================================================
// COORDINATE FRAME INDICES - Type-safe indexing
// ============================================================================

/// Type-safe burst-relative index
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct BurstIdx(pub usize);

/// Type-safe subswath-relative index
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SubswathIdx(pub usize);

/// Type-safe stitched/merged image index
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct StitchedIdx(pub usize);

/// Type-safe line index in any coordinate frame
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct LineIdx<Frame>(pub usize, std::marker::PhantomData<Frame>);

/// Type-safe pixel/sample index in any coordinate frame
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct PixelIdx<Frame>(pub usize, std::marker::PhantomData<Frame>);

// Coordinate frame markers
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct BurstFrame;
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct SubswathFrame;
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct StitchedFrame;

// Type aliases for common coordinate pairs
pub type BurstLineIdx = LineIdx<BurstFrame>;
pub type BurstPixelIdx = PixelIdx<BurstFrame>;
pub type SubswathLineIdx = LineIdx<SubswathFrame>;
pub type SubswathPixelIdx = PixelIdx<SubswathFrame>;
pub type StitchedLineIdx = LineIdx<StitchedFrame>;
pub type StitchedPixelIdx = PixelIdx<StitchedFrame>;

// ============================================================================
// COORDINATE CONVERSIONS - Explicit transformations between frames
// ============================================================================

/// Burst to subswath coordinate converter
#[derive(Debug, Clone)]
pub struct BurstToSubswathConverter {
    /// Offset of this burst within the subswath
    pub burst_line_offset: usize,
    pub burst_pixel_offset: usize,
    /// Size of this burst
    pub burst_lines: usize,
    pub burst_pixels: usize,
    /// Subswath bounds for validation
    pub subswath_lines: usize,
    pub subswath_pixels: usize,
}

/// Subswath to stitched coordinate converter
#[derive(Debug, Clone)]
pub struct SubswathToStitchedConverter {
    /// Global offset of this subswath in the stitched image
    pub subswath_line_offset: usize,
    pub subswath_pixel_offset: usize,
    /// Size of this subswath
    pub subswath_lines: usize,
    pub subswath_pixels: usize,
    /// Stitched image bounds for validation
    pub stitched_lines: usize,
    pub stitched_pixels: usize,
}

/// Complete burst to stitched coordinate converter (composition)
#[derive(Debug, Clone)]
pub struct BurstToStitchedConverter {
    pub burst_to_subswath: BurstToSubswathConverter,
    pub subswath_to_stitched: SubswathToStitchedConverter,
}

// ============================================================================
// INDEX IMPLEMENTATIONS
// ============================================================================

impl BurstIdx {
    pub fn new(value: usize) -> Self {
        Self(value)
    }
    
    pub fn value(&self) -> usize {
        self.0
    }
}

impl SubswathIdx {
    pub fn new(value: usize) -> Self {
        Self(value)
    }
    
    pub fn value(&self) -> usize {
        self.0
    }
}

impl StitchedIdx {
    pub fn new(value: usize) -> Self {
        Self(value)
    }
    
    pub fn value(&self) -> usize {
        self.0
    }
}

impl<Frame> LineIdx<Frame> {
    pub fn new(value: usize) -> Self {
        Self(value, std::marker::PhantomData)
    }
    
    pub fn value(&self) -> usize {
        self.0
    }
    
    pub fn checked_add(&self, offset: usize) -> Option<Self> {
        self.0.checked_add(offset).map(|v| Self::new(v))
    }
    
    pub fn checked_sub(&self, offset: usize) -> Option<Self> {
        self.0.checked_sub(offset).map(|v| Self::new(v))
    }
}

impl<Frame> PixelIdx<Frame> {
    pub fn new(value: usize) -> Self {
        Self(value, std::marker::PhantomData)
    }
    
    pub fn value(&self) -> usize {
        self.0
    }
    
    pub fn checked_add(&self, offset: usize) -> Option<Self> {
        self.0.checked_add(offset).map(|v| Self::new(v))
    }
    
    pub fn checked_sub(&self, offset: usize) -> Option<Self> {
        self.0.checked_sub(offset).map(|v| Self::new(v))
    }
}

// ============================================================================
// COORDINATE POSITION STRUCTURES
// ============================================================================

/// Type-safe coordinate position in a specific frame
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CoordinatePosition<Frame> {
    pub line: LineIdx<Frame>,
    pub pixel: PixelIdx<Frame>,
}

pub type BurstPosition = CoordinatePosition<BurstFrame>;
pub type SubswathPosition = CoordinatePosition<SubswathFrame>;
pub type StitchedPosition = CoordinatePosition<StitchedFrame>;

impl<Frame> CoordinatePosition<Frame> {
    pub fn new(line: usize, pixel: usize) -> Self {
        Self {
            line: LineIdx::new(line),
            pixel: PixelIdx::new(pixel),
        }
    }
    
    pub fn from_indices(line: LineIdx<Frame>, pixel: PixelIdx<Frame>) -> Self {
        Self { line, pixel }
    }
}

// ============================================================================
// CONVERTER IMPLEMENTATIONS
// ============================================================================

impl BurstToSubswathConverter {
    pub fn new(
        burst_line_offset: usize,
        burst_pixel_offset: usize,
        burst_lines: usize,
        burst_pixels: usize,
        subswath_lines: usize,
        subswath_pixels: usize,
    ) -> Self {
        Self {
            burst_line_offset,
            burst_pixel_offset,
            burst_lines,
            burst_pixels,
            subswath_lines,
            subswath_pixels,
        }
    }
    
    /// Convert burst coordinates to subswath coordinates
    pub fn convert_position(&self, burst_pos: BurstPosition) -> SarResult<SubswathPosition> {
        // Validate burst coordinates are within bounds
        if burst_pos.line.value() >= self.burst_lines {
            return Err(crate::types::SarError::InvalidParameter(format!(
                "Burst line index {} exceeds burst bounds {}", 
                burst_pos.line.value(), self.burst_lines
            )));
        }
        
        if burst_pos.pixel.value() >= self.burst_pixels {
            return Err(crate::types::SarError::InvalidParameter(format!(
                "Burst pixel index {} exceeds burst bounds {}", 
                burst_pos.pixel.value(), self.burst_pixels
            )));
        }
        
        // Convert to subswath coordinates
        let subswath_line = burst_pos.line.value() + self.burst_line_offset;
        let subswath_pixel = burst_pos.pixel.value() + self.burst_pixel_offset;
        
        // Validate converted coordinates are within subswath bounds
        if subswath_line >= self.subswath_lines {
            return Err(crate::types::SarError::InvalidParameter(format!(
                "Converted subswath line {} exceeds subswath bounds {}", 
                subswath_line, self.subswath_lines
            )));
        }
        
        if subswath_pixel >= self.subswath_pixels {
            return Err(crate::types::SarError::InvalidParameter(format!(
                "Converted subswath pixel {} exceeds subswath bounds {}", 
                subswath_pixel, self.subswath_pixels
            )));
        }
        
        Ok(SubswathPosition::new(subswath_line, subswath_pixel))
    }
    
    /// Convert individual line index
    pub fn convert_line(&self, burst_line: BurstLineIdx) -> SarResult<SubswathLineIdx> {
        if burst_line.value() >= self.burst_lines {
            return Err(crate::types::SarError::InvalidParameter(format!(
                "Burst line index {} exceeds burst bounds {}", 
                burst_line.value(), self.burst_lines
            )));
        }
        
        let subswath_line = burst_line.value() + self.burst_line_offset;
        if subswath_line >= self.subswath_lines {
            return Err(crate::types::SarError::InvalidParameter(format!(
                "Converted subswath line {} exceeds subswath bounds {}", 
                subswath_line, self.subswath_lines
            )));
        }
        
        Ok(SubswathLineIdx::new(subswath_line))
    }
    
    /// Convert individual pixel index
    pub fn convert_pixel(&self, burst_pixel: BurstPixelIdx) -> SarResult<SubswathPixelIdx> {
        if burst_pixel.value() >= self.burst_pixels {
            return Err(crate::types::SarError::InvalidParameter(format!(
                "Burst pixel index {} exceeds burst bounds {}", 
                burst_pixel.value(), self.burst_pixels
            )));
        }
        
        let subswath_pixel = burst_pixel.value() + self.burst_pixel_offset;
        if subswath_pixel >= self.subswath_pixels {
            return Err(crate::types::SarError::InvalidParameter(format!(
                "Converted subswath pixel {} exceeds subswath bounds {}", 
                subswath_pixel, self.subswath_pixels
            )));
        }
        
        Ok(SubswathPixelIdx::new(subswath_pixel))
    }
}

impl SubswathToStitchedConverter {
    pub fn new(
        subswath_line_offset: usize,
        subswath_pixel_offset: usize,
        subswath_lines: usize,
        subswath_pixels: usize,
        stitched_lines: usize,
        stitched_pixels: usize,
    ) -> Self {
        Self {
            subswath_line_offset,
            subswath_pixel_offset,
            subswath_lines,
            subswath_pixels,
            stitched_lines,
            stitched_pixels,
        }
    }
    
    /// Convert subswath coordinates to stitched coordinates
    pub fn convert_position(&self, subswath_pos: SubswathPosition) -> SarResult<StitchedPosition> {
        // Validate subswath coordinates are within bounds
        if subswath_pos.line.value() >= self.subswath_lines {
            return Err(crate::types::SarError::InvalidParameter(format!(
                "Subswath line index {} exceeds subswath bounds {}", 
                subswath_pos.line.value(), self.subswath_lines
            )));
        }
        
        if subswath_pos.pixel.value() >= self.subswath_pixels {
            return Err(crate::types::SarError::InvalidParameter(format!(
                "Subswath pixel index {} exceeds subswath bounds {}", 
                subswath_pos.pixel.value(), self.subswath_pixels
            )));
        }
        
        // Convert to stitched coordinates
        let stitched_line = subswath_pos.line.value() + self.subswath_line_offset;
        let stitched_pixel = subswath_pos.pixel.value() + self.subswath_pixel_offset;
        
        // Validate converted coordinates are within stitched bounds
        if stitched_line >= self.stitched_lines {
            return Err(crate::types::SarError::InvalidParameter(format!(
                "Converted stitched line {} exceeds stitched bounds {}", 
                stitched_line, self.stitched_lines
            )));
        }
        
        if stitched_pixel >= self.stitched_pixels {
            return Err(crate::types::SarError::InvalidParameter(format!(
                "Converted stitched pixel {} exceeds stitched bounds {}", 
                stitched_pixel, self.stitched_pixels
            )));
        }
        
        Ok(StitchedPosition::new(stitched_line, stitched_pixel))
    }
    
    /// Convert individual line index
    pub fn convert_line(&self, subswath_line: SubswathLineIdx) -> SarResult<StitchedLineIdx> {
        if subswath_line.value() >= self.subswath_lines {
            return Err(crate::types::SarError::InvalidParameter(format!(
                "Subswath line index {} exceeds subswath bounds {}", 
                subswath_line.value(), self.subswath_lines
            )));
        }
        
        let stitched_line = subswath_line.value() + self.subswath_line_offset;
        if stitched_line >= self.stitched_lines {
            return Err(crate::types::SarError::InvalidParameter(format!(
                "Converted stitched line {} exceeds stitched bounds {}", 
                stitched_line, self.stitched_lines
            )));
        }
        
        Ok(StitchedLineIdx::new(stitched_line))
    }
    
    /// Convert individual pixel index
    pub fn convert_pixel(&self, subswath_pixel: SubswathPixelIdx) -> SarResult<StitchedPixelIdx> {
        if subswath_pixel.value() >= self.subswath_pixels {
            return Err(crate::types::SarError::InvalidParameter(format!(
                "Subswath pixel index {} exceeds subswath bounds {}", 
                subswath_pixel.value(), self.subswath_pixels
            )));
        }
        
        let stitched_pixel = subswath_pixel.value() + self.subswath_pixel_offset;
        if stitched_pixel >= self.stitched_pixels {
            return Err(crate::types::SarError::InvalidParameter(format!(
                "Converted stitched pixel {} exceeds stitched bounds {}", 
                stitched_pixel, self.stitched_pixels
            )));
        }
        
        Ok(StitchedPixelIdx::new(stitched_pixel))
    }
}

impl BurstToStitchedConverter {
    pub fn new(
        burst_to_subswath: BurstToSubswathConverter,
        subswath_to_stitched: SubswathToStitchedConverter,
    ) -> Self {
        Self {
            burst_to_subswath,
            subswath_to_stitched,
        }
    }
    
    /// Convert burst coordinates directly to stitched coordinates
    pub fn convert_position(&self, burst_pos: BurstPosition) -> SarResult<StitchedPosition> {
        let subswath_pos = self.burst_to_subswath.convert_position(burst_pos)?;
        self.subswath_to_stitched.convert_position(subswath_pos)
    }
    
    /// Convert individual line index
    pub fn convert_line(&self, burst_line: BurstLineIdx) -> SarResult<StitchedLineIdx> {
        let subswath_line = self.burst_to_subswath.convert_line(burst_line)?;
        self.subswath_to_stitched.convert_line(subswath_line)
    }
    
    /// Convert individual pixel index
    pub fn convert_pixel(&self, burst_pixel: BurstPixelIdx) -> SarResult<StitchedPixelIdx> {
        let subswath_pixel = self.burst_to_subswath.convert_pixel(burst_pixel)?;
        self.subswath_to_stitched.convert_pixel(subswath_pixel)
    }
}

// ============================================================================
// DISPLAY IMPLEMENTATIONS
// ============================================================================

impl fmt::Display for BurstIdx {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "B{}", self.0)
    }
}

impl fmt::Display for SubswathIdx {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "IW{}", self.0 + 1) // IW subswaths are 1-indexed
    }
}

impl fmt::Display for StitchedIdx {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "S{}", self.0)
    }
}

impl<Frame> fmt::Display for LineIdx<Frame> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "L{}", self.0)
    }
}

impl<Frame> fmt::Display for PixelIdx<Frame> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "P{}", self.0)
    }
}

impl<Frame> fmt::Display for CoordinatePosition<Frame> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "({}, {})", self.line, self.pixel)
    }
}

// ============================================================================
// UTILITY FUNCTIONS - Safe coordinate operations
// ============================================================================

/// Coordinate frame utilities
pub mod coord_utils {
    use super::*;
    use crate::types::SarError;
    
    /// Create burst to subswath converter from dimensions
    pub fn create_burst_to_subswath_converter(
        burst_lines: u32,
        burst_samples: u32,
        subswath_lines: u32,
        subswath_samples: u32,
    ) -> SarResult<BurstToSubswathConverter> {
        // Validate burst is within subswath bounds
        if burst_lines > subswath_lines {
            return Err(SarError::InvalidParameter(format!(
                "Burst extends beyond subswath azimuth bounds: {} > {}",
                burst_lines, subswath_lines
            )));
        }
        
        if burst_samples > subswath_samples {
            return Err(SarError::InvalidParameter(format!(
                "Burst extends beyond subswath range bounds: {} > {}",
                burst_samples, subswath_samples
            )));
        }
        
        Ok(BurstToSubswathConverter {
            burst_line_offset: 0,
            burst_pixel_offset: 0,
            burst_lines: burst_lines as usize,
            burst_pixels: burst_samples as usize,
            subswath_lines: subswath_lines as usize,
            subswath_pixels: subswath_samples as usize,
        })
    }
    
    /// Create subswath to stitched converter from dimensions
    pub fn create_subswath_to_stitched_converter(
        subswath_lines: usize,
        subswath_pixels: usize,
        stitched_lines: usize,
        stitched_pixels: usize,
        line_offset: usize,
        pixel_offset: usize,
    ) -> SarResult<SubswathToStitchedConverter> {
        // Validate subswath fits within stitched bounds
        let subswath_end_line = line_offset + subswath_lines;
        let subswath_end_pixel = pixel_offset + subswath_pixels;
        
        if subswath_end_line > stitched_lines {
            return Err(SarError::InvalidParameter(format!(
                "Subswath extends beyond stitched azimuth bounds: {} > {}",
                subswath_end_line, stitched_lines
            )));
        }
        
        if subswath_end_pixel > stitched_pixels {
            return Err(SarError::InvalidParameter(format!(
                "Subswath extends beyond stitched range bounds: {} > {}",
                subswath_end_pixel, stitched_pixels
            )));
        }
        
        Ok(SubswathToStitchedConverter::new(
            line_offset,
            pixel_offset,
            subswath_lines,
            subswath_pixels,
            stitched_lines,
            stitched_pixels,
        ))
    }
    
    /// Check if position is within bounds for a given frame
    pub fn is_position_valid<Frame>(
        pos: &CoordinatePosition<Frame>,
        max_lines: usize,
        max_pixels: usize,
    ) -> bool {
        pos.line.value() < max_lines && pos.pixel.value() < max_pixels
    }
    
    /// Safe array access using type-safe coordinates
    pub fn safe_array_access<'a, T, Frame>(
        array: &'a ndarray::Array2<T>,
        pos: &CoordinatePosition<Frame>,
    ) -> Option<&'a T> {
        let (rows, cols) = array.dim();
        if is_position_valid(pos, rows, cols) {
            Some(&array[[pos.line.value(), pos.pixel.value()]])
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_coordinate_index_creation() {
        let burst_idx = BurstIdx::new(5);
        assert_eq!(burst_idx.value(), 5);
        assert_eq!(format!("{}", burst_idx), "B5");
        
        let line_idx = BurstLineIdx::new(100);
        assert_eq!(line_idx.value(), 100);
        assert_eq!(format!("{}", line_idx), "L100");
        
        let pos = BurstPosition::new(100, 200);
        assert_eq!(pos.line.value(), 100);
        assert_eq!(pos.pixel.value(), 200);
        assert_eq!(format!("{}", pos), "(L100, P200)");
    }
    
    #[test]
    fn test_burst_to_subswath_conversion() {
        let converter = BurstToSubswathConverter::new(
            50,    // burst line offset
            100,   // burst pixel offset
            1000,  // burst lines
            2000,  // burst pixels
            5000,  // subswath lines
            8000,  // subswath pixels
        );
        
        let burst_pos = BurstPosition::new(500, 1000);
        let subswath_pos = converter.convert_position(burst_pos).unwrap();
        
        assert_eq!(subswath_pos.line.value(), 550); // 500 + 50
        assert_eq!(subswath_pos.pixel.value(), 1100); // 1000 + 100
    }
    
    #[test]
    fn test_subswath_to_stitched_conversion() {
        let converter = SubswathToStitchedConverter::new(
            1000,  // subswath line offset
            2000,  // subswath pixel offset
            3000,  // subswath lines
            4000,  // subswath pixels
            10000, // stitched lines
            15000, // stitched pixels
        );
        
        let subswath_pos = SubswathPosition::new(1500, 2500);
        let stitched_pos = converter.convert_position(subswath_pos).unwrap();
        
        assert_eq!(stitched_pos.line.value(), 2500); // 1500 + 1000
        assert_eq!(stitched_pos.pixel.value(), 4500); // 2500 + 2000
    }
    
    #[test]
    fn test_bounds_checking() {
        let converter = BurstToSubswathConverter::new(
            0, 0, 100, 200, 1000, 2000
        );
        
        // Valid conversion
        let valid_pos = BurstPosition::new(50, 100);
        assert!(converter.convert_position(valid_pos).is_ok());
        
        // Invalid burst coordinates
        let invalid_pos = BurstPosition::new(150, 100); // Line 150 >= 100 burst lines
        assert!(converter.convert_position(invalid_pos).is_err());
        
        let invalid_pos = BurstPosition::new(50, 250); // Pixel 250 >= 200 burst pixels
        assert!(converter.convert_position(invalid_pos).is_err());
    }
    
    #[test]
    fn test_type_safety() {
        // These should compile and work correctly
        let burst_line = BurstLineIdx::new(100);
        let subswath_line = SubswathLineIdx::new(200);
        let stitched_line = StitchedLineIdx::new(300);
        
        // Type safety prevents mixing incompatible coordinate frames
        // These would NOT compile:
        // let bad = burst_line + subswath_line;    // Cannot mix coordinate frames
        // let bad = burst_pos == subswath_pos;     // Cannot compare different frames
        
        // But operations within the same frame work fine
        let burst_pos1 = BurstPosition::new(100, 200);
        let burst_pos2 = BurstPosition::new(100, 200);
        assert_eq!(burst_pos1, burst_pos2);
        
        // And explicit conversions work
        assert_eq!(burst_line.value(), 100);
        assert_eq!(subswath_line.value(), 200);
        assert_eq!(stitched_line.value(), 300);
    }
    
    #[test]
    fn test_complete_conversion_chain() {
        // Set up complete conversion chain: burst -> subswath -> stitched
        let burst_to_subswath = BurstToSubswathConverter::new(
            10, 20,   // burst offsets
            100, 200, // burst size
            1000, 2000 // subswath size
        );
        
        let subswath_to_stitched = SubswathToStitchedConverter::new(
            100, 300,  // subswath offsets
            1000, 2000, // subswath size
            5000, 10000 // stitched size
        );
        
        let complete_converter = BurstToStitchedConverter::new(
            burst_to_subswath,
            subswath_to_stitched,
        );
        
        // Test complete conversion
        let burst_pos = BurstPosition::new(50, 100);
        let stitched_pos = complete_converter.convert_position(burst_pos).unwrap();
        
        // burst(50,100) -> subswath(60,120) -> stitched(160,420)
        assert_eq!(stitched_pos.line.value(), 160);  // 50 + 10 + 100
        assert_eq!(stitched_pos.pixel.value(), 420); // 100 + 20 + 300
    }
}