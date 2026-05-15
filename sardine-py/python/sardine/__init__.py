"""SARdine — Sentinel-1 SAR backscatter processing.

Public Python API (1:1 with the ``sardine`` CLI):

    sardine.process(safe, dem, output, geoid, **kwargs)
    sardine.grd(safe, output, **kwargs)

Optional, feature-gated helpers (only available when the wheel was built
with the corresponding Cargo feature — inspect with ``sardine.features()``):

    sardine.fetch_orbit(safe, cache_dir)              # 'orbit-fetch'
    sardine.download_slc(product_id, output_dir, …)   # 'slc-fetch'
    sardine.fetch_geoid()                             # 'geoid-fetch'

All Rust calls release the GIL and raise ``RuntimeError`` (with the full
chained Rust diagnostic) on failure.  See each function's docstring for
the complete argument list and defaults — these match the corresponding
CLI flags exactly.

Example
-------
>>> import sardine
>>> sardine.features()
{'geoid_fetch': True, 'orbit_fetch': True, 'slc_fetch': False}
>>> sardine.process(
...     safe="/path/to/S1A_IW_SLC__1SDV_….SAFE",
...     dem="/path/to/srtm1/",
...     output="/tmp/out.tif",
...     geoid="auto",        # requires geoid-fetch feature
...     crs="auto",          # pick UTM zone from scene centre
...     pixel_spacing_m=10.0,
...     write_lia=True,
... )
"""

from ._sardine import (
    __version__,
    download_slc,
    features,
    fetch_geoid,
    fetch_orbit,
    grd,
    insar,
    polsar,
    process,
)

__all__ = [
    "__version__",
    "download_slc",
    "features",
    "fetch_geoid",
    "fetch_orbit",
    "grd",
    "insar",
    "polsar",
    "process",
]
