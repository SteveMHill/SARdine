"""Smoke tests for the `sardine` Python extension.

These are intentionally cheap: they only verify that the module loads,
the public API is the one we expect, and that obviously bad inputs
surface as `RuntimeError` (rather than panicking the interpreter, hanging,
or silently succeeding).

They do **not** exercise the full pipeline — that requires real SAFE +
DEM + orbit data and is covered by the Rust integration suite.

Run with::

    cd new_pipeline/sardine-py
    maturin develop --release             # or --features fetch-all
    pytest -q
"""

from __future__ import annotations

import importlib

import pytest

sardine = importlib.import_module("sardine")


# --- module surface ---------------------------------------------------------


def test_version_is_nonempty_string() -> None:
    assert isinstance(sardine.__version__, str)
    assert sardine.__version__  # non-empty


@pytest.mark.parametrize(
    "name",
    ["process", "grd", "fetch_orbit", "download_slc", "fetch_geoid", "features"],
)
def test_public_symbol_is_callable(name: str) -> None:
    assert callable(getattr(sardine, name)), f"sardine.{name} should be callable"


def test_features_returns_dict_with_expected_keys() -> None:
    feats = sardine.features()
    assert isinstance(feats, dict)
    for key in ("orbit_fetch", "slc_fetch", "geoid_fetch"):
        assert key in feats, f"features() missing {key!r}"
        assert isinstance(feats[key], bool)


# --- error paths ------------------------------------------------------------


def test_process_with_bogus_paths_raises_runtime_error(tmp_path) -> None:
    """A non-existent SAFE must surface as a RuntimeError, not a panic."""
    out = tmp_path / "out.tif"
    with pytest.raises(RuntimeError):
        sardine.process(
            str(tmp_path / "no_such.SAFE"),
            str(tmp_path / "no_such_dem"),
            str(out),
            "zero",
        )
    assert not out.exists(), "process must not write output on failure"


def test_grd_with_bogus_safe_raises_runtime_error(tmp_path) -> None:
    out = tmp_path / "grd.tif"
    with pytest.raises(RuntimeError):
        sardine.grd(str(tmp_path / "no_such.SAFE"), str(out))
    assert not out.exists()


def test_fetch_orbit_with_bogus_safe_raises_runtime_error(tmp_path) -> None:
    """If `orbit-fetch` is compiled in, parsing the bogus SAFE fails;
    if not, the missing-feature error is raised. Either way: RuntimeError."""
    with pytest.raises(RuntimeError):
        sardine.fetch_orbit(str(tmp_path / "no_such.SAFE"), str(tmp_path / "cache"))


# --- new kwargs accepted / rejected -----------------------------------------


def test_process_accepts_valid_mode_kwarg(tmp_path) -> None:
    """mode="nrb" must be accepted (pipeline fails later on missing files, not at parse)."""
    with pytest.raises(RuntimeError):
        sardine.process(
            str(tmp_path / "no_such.SAFE"),
            str(tmp_path / "no_dem"),
            str(tmp_path / "out.tif"),
            "zero",
            mode="nrb",
        )


def test_process_rejects_invalid_mode(tmp_path) -> None:
    """An unrecognised mode string must raise RuntimeError, not panic."""
    with pytest.raises(RuntimeError, match="(?i)mode|invalid|unknown"):
        sardine.process(
            str(tmp_path / "no_such.SAFE"),
            str(tmp_path / "no_dem"),
            str(tmp_path / "out.tif"),
            "zero",
            mode="beta0",
        )


def test_process_accepts_valid_speckle_order(tmp_path) -> None:
    with pytest.raises(RuntimeError):
        sardine.process(
            str(tmp_path / "no_such.SAFE"),
            str(tmp_path / "no_dem"),
            str(tmp_path / "out.tif"),
            "zero",
            speckle_order="before",
        )


def test_process_rejects_invalid_speckle_order(tmp_path) -> None:
    with pytest.raises(RuntimeError):
        sardine.process(
            str(tmp_path / "no_such.SAFE"),
            str(tmp_path / "no_dem"),
            str(tmp_path / "out.tif"),
            "zero",
            speckle_order="during",
        )


def test_process_accepts_valid_iw(tmp_path) -> None:
    with pytest.raises(RuntimeError):
        sardine.process(
            str(tmp_path / "no_such.SAFE"),
            str(tmp_path / "no_dem"),
            str(tmp_path / "out.tif"),
            "zero",
            iw="IW1",
        )


def test_process_rejects_invalid_iw(tmp_path) -> None:
    with pytest.raises(RuntimeError):
        sardine.process(
            str(tmp_path / "no_such.SAFE"),
            str(tmp_path / "no_dem"),
            str(tmp_path / "out.tif"),
            "zero",
            iw="IW9",
        )


def test_grd_accepts_valid_iw(tmp_path) -> None:
    with pytest.raises(RuntimeError):
        sardine.grd(str(tmp_path / "no_such.SAFE"), str(tmp_path / "grd.tif"), iw="IW2")


def test_grd_rejects_invalid_iw(tmp_path) -> None:
    with pytest.raises(RuntimeError):
        sardine.grd(str(tmp_path / "no_such.SAFE"), str(tmp_path / "grd.tif"), iw="IW9")


def test_download_slc_with_bogus_token_raises_runtime_error(tmp_path) -> None:
    """If `slc-fetch` is compiled in, the bogus token fails the auth path;
    if not, the missing-feature error is raised. Either way: RuntimeError."""
    with pytest.raises(RuntimeError):
        sardine.download_slc(
            "S1A_IW_SLC__1SDV_INVALID_PRODUCT_ID",
            str(tmp_path / "dl"),
            token="not-a-real-token",
        )
