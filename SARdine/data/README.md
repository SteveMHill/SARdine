# Test Data

This directory is for Sentinel-1 SLC test data files.

## Getting Test Data

To test SARdine, you'll need Sentinel-1 SLC data. You can download test data from:

1. **Copernicus Open Access Hub**: https://scihub.copernicus.eu/
2. **Alaska Satellite Facility (ASF)**: https://search.asf.alaska.edu/
3. **ESA Science Hub**: https://scihub.copernicus.eu/dhus/

## Recommended Test Files

For initial testing, we recommend downloading a small Sentinel-1 IW SLC scene:
- Product type: SLC (Single Look Complex)
- Mode: IW (Interferometric Wide swath)
- Size: Usually 4-8 GB per scene

## File Organization

Place downloaded `.zip` files directly in this directory:
```
data/
├── README.md
└── S1A_IW_SLC_*.zip   # Your Sentinel-1 scenes
```

## Note

Large test data files (`.zip`, `.SAFE/`) are automatically ignored by git to keep the repository lightweight.
