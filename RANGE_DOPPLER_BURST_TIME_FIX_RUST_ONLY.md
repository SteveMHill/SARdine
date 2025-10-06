# Range-Doppler Time Base Fix - REVISED IMPLEMENTATION PLAN

## Root Cause (Confirmed)

**The burst timing information already exists in Rust but isn't being used by the Range-Doppler solver.**

### What We Have:
✅ Burst `azimuth_time` extracted in `annotation.rs` (line 1145)
✅ Burst-aware time calculation in `topsar_merge.rs` (`AzimuthTimingModel`)
✅ All annotation data parsed in Rust, no Python parsing needed

### What's Missing:
❌ Burst timing not passed to `RangeDopplerParams`
❌ Newton-Raphson solver uses `product_start_rel_s` (manifest time) instead of `first_burst_time`
❌ Result: 100% Range-Doppler failures because solver searches wrong time window

## Simplified Fix (Rust-Only, No Python Changes)

### Step 1: Add First Burst Time to RangeDopplerParams

**File**: `SARdine/src/core/terrain_correction.rs`

Add to `RangeDopplerParams` structure (around line 1107):

```rust
pub struct RangeDopplerParams {
    // ... existing fields ...
    
    /// CRITICAL FIX: First burst azimuth time (seconds relative to orbit_ref_epoch)
    /// This is the ACTUAL data start time for TOPS/IW data
    /// DO NOT use product_start_rel_s which is the manifest time (wrong!)
    pub first_burst_time_rel_orbit: Option<f64>,
    
    /// Optional: All burst times for per-line burst-aware geocoding
    pub burst_times: Option<Vec<f64>>, // Each burst start time, orbit-relative
}
```

### Step 2: Extract First Burst Time When Creating RangeDopplerParams

**File**: `SARdine/src/lib.rs` (around line 3220)

Modify where `RangeDopplerParams` is created:

```rust
// CRITICAL FIX: Extract first burst azimuth time from annotation
let first_burst_time_rel_orbit = if let Some(annotation_data) = real_metadata.get("annotation") {
    // Parse burst list from annotation if available
    extract_first_burst_time(annotation_data, orbit_ref_epoch_utc)
} else {
    // Fallback: use product_start_rel_s (less accurate)
    log::warn!("⚠️  No burst timing available, using product_start_rel_s as fallback");
    None
};

let rd_params = crate::core::terrain_correction::RangeDopplerParams {
    // ... existing fields ...
    first_burst_time_rel_orbit,
    burst_times: None, // Optional: extract all burst times later
};
```

### Step 3: Modify Newton-Raphson Initial Guess

**File**: `SARdine/src/core/terrain_correction.rs` (around line 2480)

Change initial guess calculation:

```rust
// CRITICAL FIX: Use first burst time instead of product start time
let initial_time_base = params.first_burst_time_rel_orbit
    .unwrap_or(params.product_start_rel_s);
    
// Initial guess: mid-burst time
let initial_guess = initial_time_base + (params.product_duration / 2.0);

log::info!(
    "🎯 BURST-AWARE TIME BASE:\n\
     ├─ product_start_rel_s = {:.6}s (manifest time)\n\
     ├─ first_burst_time_rel_orbit = {:.6}s (actual data start)\n\
     ├─ Δt (burst vs manifest) = {:.6}s\n\
     └─ Using burst time for Range-Doppler solver",
    params.product_start_rel_s,
    initial_time_base,
    initial_time_base - params.product_start_rel_s
);
```

### Step 4: Update Azimuth Time Calculation for Per-Pixel

**File**: `SARdine/src/core/terrain_correction.rs` (around line 2315)

Change `azimuth_time_from_start` calculation:

```rust
// CRITICAL FIX: Compute azimuth time relative to FIRST BURST, not product manifest
let time_base = params.first_burst_time_rel_orbit.unwrap_or(params.product_start_rel_s);
let azimuth_time_from_start = azimuth_time_rel_orbit - time_base;

// Update validation to use burst-based expected time
let expected_max_time = params.product_duration; // Time from first burst to end
```

### Step 5: Helper Function to Extract First Burst Time

**File**: `SARdine/src/lib.rs` (new function)

```rust
/// Extract first burst azimuth time from annotation data
/// Returns time in seconds relative to orbit_ref_epoch
fn extract_first_burst_time(
    annotation: &crate::io::annotation::AnnotationRoot,
    orbit_ref_epoch_utc: f64
) -> Option<f64> {
    let swath_timing = annotation.swath_timing.as_ref()?;
    let burst_list = swath_timing.burst_list.as_ref()?;
    let bursts = burst_list.bursts.as_ref()?;
    
    if bursts.is_empty() {
        return None;
    }
    
    // Get first burst azimuth time
    let first_burst = &bursts[0];
    if let Some(azimuth_time_str) = &first_burst.azimuth_time {
        // Parse datetime string to UTC seconds
        use chrono::{DateTime, Utc};
        if let Ok(dt) = DateTime::parse_from_rfc3339(azimuth_time_str) {
            let burst_time_utc = dt.timestamp() as f64 
                + (dt.timestamp_subsec_nanos() as f64) * 1e-9;
            let burst_time_rel = burst_time_utc - orbit_ref_epoch_utc;
            
            log::info!(
                "✅ Extracted first burst time:\n\
                 ├─ Burst azimuth time (UTC): {:.6}s\n\
                 ├─ Orbit ref epoch (UTC): {:.6}s\n\
                 └─ Burst time (orbit-relative): {:.6}s",
                burst_time_utc,
                orbit_ref_epoch_utc,
                burst_time_rel
            );
            
            return Some(burst_time_rel);
        }
    }
    
    None
}
```

## Expected Results After Fix

### Immediate Changes:
- Newton-Raphson solver will search around correct time window
- Doppler residual will bracket zero (sign change detected)
- First pixels will geocode successfully
- Range-Doppler failure rate drops from 100% to <10%

### Diagnostic Output:
```
🎯 BURST-AWARE TIME BASE:
├─ product_start_rel_s = 67983.46s (manifest time)
├─ first_burst_time_rel_orbit = 67985.12s (actual data start)  
├─ Δt (burst vs manifest) = 1.66s
└─ Using burst time for Range-Doppler solver

🔍 DOPPLER RESIDUAL BRACKETING SCAN:
  k=-10, t=67984.5s, f_d=-125.3 Hz
  k=  0, t=67985.5s, f_d=-15.7 Hz
  k= 10, t=67986.5s, f_d=95.8 Hz
  ✅ Sign change detected between t=67985.45s and t=67985.50s
  ✅ Doppler residual brackets zero - Newton-Raphson should converge
```

## Why No Python Changes Needed

1. **Annotation parsing is in Rust** (`annotation.rs`)
2. **Burst data structures exist** (`struct Burst` with `azimuth_time`)
3. **Just need to connect the data** to `RangeDopplerParams`
4. **Python only passes through metadata** that Rust already extracted

## Implementation Order

1. ✅ Add `first_burst_time_rel_orbit` field to `RangeDopplerParams`
2. ✅ Create `extract_first_burst_time()` helper function
3. ✅ Populate field when creating `RangeDopplerParams` in `lib.rs`
4. ✅ Update Newton-Raphson initial guess to use burst time
5. ✅ Update per-pixel azimuth time calculation
6. ✅ Test and verify Doppler residual brackets zero

## Testing

```bash
# Build with new changes
cd SARdine && cargo build --release

# Test - should see immediate improvement
python3 -m sardine.cli backscatter \
  data/S1A_IW_SLC__1SDV_20201230T165244_20201230T165311_035918_0434F0_6788.SAFE \
  pipeline_output \
  --resolution 20 \
  --polarization VV

# Look for:
# - "✅ Extracted first burst time"
# - "🎯 BURST-AWARE TIME BASE"  
# - "✅ Sign change detected"
# - Range-Doppler failure rate < 10% (was 100%)
```

## Key Insight

**The fix is entirely in Rust plumbing** - connecting existing burst timing data to where the Range-Doppler solver needs it. No XML parsing in Python needed because Rust already does it!
