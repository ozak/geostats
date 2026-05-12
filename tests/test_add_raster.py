"""
Tests for add_raster() — registering user-provided rasters as first-class measures.
"""
import os
import numpy as np
import pytest
import rasterio
from rasterio.transform import from_bounds
import geopandas as gpd
from shapely.geometry import box

from geostats.main import (
    add_raster,
    pathmeasures, wgs84_measures, cea_measures, main_measures,
    namemeasures, _user_single_files, _user_registered,
)

GEOSTATS_DATA = os.path.join(os.path.expanduser('~'), 'geostats-data')
ELEVATION_DATA = os.path.join(GEOSTATS_DATA, 'GLOBE', 'tifs')


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_tif(path, epsg, fill=5.0, bounds=(-20.0, -20.0, 20.0, 20.0), size=(40, 40)):
    """Write a constant-value single-band GeoTIFF with an EPSG CRS."""
    west, south, east, north = bounds
    transform = from_bounds(west, south, east, north, size[1], size[0])
    data = np.full(size, fill, dtype=np.float32)
    with rasterio.open(
        path, 'w', driver='GTiff',
        height=size[0], width=size[1],
        count=1, dtype='float32',
        crs=f'EPSG:{epsg}',
        transform=transform,
    ) as dst:
        dst.write(data, 1)


def _write_tif_crs(path, crs_str, fill=5.0, bounds=(-1e6, -1e6, 1e6, 1e6), size=(40, 40)):
    """Write a constant-value single-band GeoTIFF with an arbitrary CRS string."""
    west, south, east, north = bounds
    transform = from_bounds(west, south, east, north, size[1], size[0])
    data = np.full(size, fill, dtype=np.float32)
    with rasterio.open(
        path, 'w', driver='GTiff',
        height=size[0], width=size[1],
        count=1, dtype='float32',
        crs=crs_str,
        transform=transform,
    ) as dst:
        dst.write(data, 1)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def wgs84_tif(tmp_path):
    p = tmp_path / 'ndvi.tif'
    _write_tif(str(p), epsg=4326, fill=5.0)
    return str(p)


@pytest.fixture
def cea_tif(tmp_path):
    """A raster in Lambert CEA (ESRI:54034)."""
    p = tmp_path / 'terrain_cea.tif'
    _write_tif_crs(str(p), crs_str='ESRI:54034', fill=3.0)
    return str(p)


@pytest.fixture
def non_wgs84_tif(tmp_path):
    """A raster in Web Mercator (EPSG:3857) — unsupported CRS."""
    p = tmp_path / 'ruggedness.tif'
    _write_tif(str(p), epsg=3857, fill=3.0, bounds=(-2e6, -2e6, 2e6, 2e6))
    return str(p)


@pytest.fixture
def wgs84_tif_dir(tmp_path):
    """Directory with two WGS84 rasters, fill values 10 and 20."""
    d = tmp_path / 'precip'
    d.mkdir()
    _write_tif(str(d / 'precip_2000.tif'), epsg=4326, fill=10.0)
    _write_tif(str(d / 'precip_2001.tif'), epsg=4326, fill=20.0)
    return str(d)


@pytest.fixture
def mixed_unsupported_dir(tmp_path):
    """Directory with one WGS84 and one Web Mercator (unsupported) raster."""
    d = tmp_path / 'mixed_bad'
    d.mkdir()
    _write_tif(str(d / 'file_a.tif'), epsg=4326, fill=5.0)
    _write_tif(str(d / 'file_b.tif'), epsg=3857, fill=5.0,
               bounds=(-2e6, -2e6, 2e6, 2e6))
    return str(d)


@pytest.fixture
def mixed_wgs84_cea_dir(tmp_path):
    """Directory with one WGS84 and one CEA raster — valid CRSes but mixed."""
    d = tmp_path / 'mixed_proj'
    d.mkdir()
    _write_tif(str(d / 'file_wgs.tif'), epsg=4326, fill=5.0)
    _write_tif_crs(str(d / 'file_cea.tif'), crs_str='ESRI:54034', fill=3.0)
    return str(d)


@pytest.fixture
def simple_gdf():
    """One polygon centred on the equator, well inside the test raster extent."""
    return gpd.GeoDataFrame(
        {'id': [1]},
        geometry=[box(-10.0, -10.0, 10.0, 10.0)],
        crs='EPSG:4326',
    )


@pytest.fixture(autouse=True)
def _cleanup(request):
    """
    Remove any measures registered during a test so tests don't bleed into
    each other (add_raster mutates module-level dicts).
    """
    names_before = set(pathmeasures.keys())
    yield
    for name in list(pathmeasures.keys()):
        if name not in names_before:
            pathmeasures.pop(name, None)
            namemeasures.pop(name, None)
            _user_single_files.pop(name, None)
            _user_registered.pop(name, None)
            for lst in (wgs84_measures, cea_measures, main_measures):
                if name in lst:
                    lst.remove(name)


# ---------------------------------------------------------------------------
# Registration tests — basic path/list wiring
# ---------------------------------------------------------------------------

def test_single_file_wgs84_goes_to_wgs84_measures(wgs84_tif):
    add_raster(wgs84_tif, name='ndvi_test')
    assert 'ndvi_test' in wgs84_measures
    assert 'ndvi_test' not in cea_measures
    assert 'ndvi_test' in main_measures


def test_single_file_path_stored_as_parent_dir(wgs84_tif):
    add_raster(wgs84_tif, name='ndvi_test')
    expected_dir = os.path.dirname(os.path.abspath(wgs84_tif)) + os.sep
    assert pathmeasures['ndvi_test'] == expected_dir


def test_single_file_basename_stored_in_user_single_files(wgs84_tif):
    add_raster(wgs84_tif, name='ndvi_test')
    assert _user_single_files['ndvi_test'] == 'ndvi.tif'


def test_single_file_namemeasures_sentinel_is_none(wgs84_tif):
    add_raster(wgs84_tif, name='ndvi_test')
    assert namemeasures['ndvi_test'] is None


def test_directory_goes_to_wgs84_measures(wgs84_tif_dir):
    add_raster(wgs84_tif_dir, name='Precip_test')
    assert 'Precip_test' in wgs84_measures


def test_directory_namemeasures_is_minus4(wgs84_tif_dir):
    add_raster(wgs84_tif_dir, name='Precip_test')
    assert namemeasures['Precip_test'] == -4


def test_directory_not_in_user_single_files(wgs84_tif_dir):
    add_raster(wgs84_tif_dir, name='Precip_test')
    assert 'Precip_test' not in _user_single_files


def test_empty_directory_raises(tmp_path):
    with pytest.raises(ValueError, match='No .tif files'):
        add_raster(str(tmp_path), name='bad')


# ---------------------------------------------------------------------------
# CRS validation tests
# ---------------------------------------------------------------------------

def test_non_wgs84_non_cea_single_file_raises(non_wgs84_tif):
    """Auto-detection must raise for CRS that is neither WGS84 nor CEA."""
    with pytest.raises(ValueError, match='neither WGS84'):
        add_raster(non_wgs84_tif, name='rug_test')


def test_directory_with_unsupported_file_raises(mixed_unsupported_dir):
    """Directory containing a non-WGS84/CEA file must raise on auto-detect."""
    with pytest.raises(ValueError, match='neither WGS84'):
        add_raster(mixed_unsupported_dir, name='mixed_test')


def test_directory_mixed_wgs84_and_cea_raises(mixed_wgs84_cea_dir):
    """Directory mixing WGS84 and CEA rasters must raise on auto-detect."""
    with pytest.raises(ValueError, match='mix of WGS84 and CEA'):
        add_raster(mixed_wgs84_cea_dir, name='mixed_test')


def test_explicit_wgs84_on_wgs84_file(wgs84_tif):
    """Explicit crs='wgs84' on a WGS84 file is valid."""
    add_raster(wgs84_tif, name='ndvi_test', crs='wgs84')
    assert 'ndvi_test' in wgs84_measures
    assert 'ndvi_test' not in cea_measures


def test_explicit_cea_on_cea_file(cea_tif):
    """Explicit crs='cea' on a CEA file is valid."""
    add_raster(cea_tif, name='terrain_test', crs='cea')
    assert 'terrain_test' in cea_measures
    assert 'terrain_test' not in wgs84_measures


def test_explicit_wgs84_on_non_wgs84_file_raises(non_wgs84_tif):
    """Explicit crs='wgs84' on a fully-unsupported CRS file must raise.

    EPSG:3857 is not WGS84 or CEA, so _classify_raster_crs raises
    'neither WGS84' before we even reach the mismatch check.
    """
    with pytest.raises(ValueError, match='neither WGS84'):
        add_raster(non_wgs84_tif, name='rug_test', crs='wgs84')


def test_explicit_cea_on_wgs84_file_raises(wgs84_tif):
    """Explicit crs='cea' on a WGS84 file must raise."""
    with pytest.raises(ValueError, match='CRS mismatch'):
        add_raster(wgs84_tif, name='ndvi_test', crs='cea')


def test_explicit_wgs84_on_mixed_dir_raises(mixed_wgs84_cea_dir):
    """Explicit crs='wgs84' on a dir containing a CEA file must raise."""
    with pytest.raises(ValueError, match='CRS mismatch'):
        add_raster(mixed_wgs84_cea_dir, name='mixed_test', crs='wgs84')


# ---------------------------------------------------------------------------
# Computation tests
# ---------------------------------------------------------------------------

def test_single_file_produces_named_columns(wgs84_tif, simple_gdf):
    from geostats.main import geostats as GeoStats
    add_raster(wgs84_tif, name='ndvi_test')
    G = GeoStats(simple_gdf, measures=['ndvi_test'],
                 stats=['mean', 'count'], add_stats={}, adds=True)
    G.geostats()
    assert 'ndvi_testmean' in G.df.columns
    assert 'ndvi_testcount' in G.df.columns


def test_single_file_mean_matches_fill_value(wgs84_tif, simple_gdf):
    from geostats.main import geostats as GeoStats
    add_raster(wgs84_tif, name='ndvi_test')
    G = GeoStats(simple_gdf, measures=['ndvi_test'],
                 stats=['mean'], add_stats={}, adds=True)
    G.geostats()
    assert G.df['ndvi_testmean'].iloc[0] == pytest.approx(5.0, abs=0.01)


def test_directory_produces_per_file_columns(wgs84_tif_dir, simple_gdf):
    from geostats.main import geostats as GeoStats
    add_raster(wgs84_tif_dir, name='Precip_test')
    G = GeoStats(simple_gdf, measures=['Precip_test'],
                 stats=['mean'], add_stats={}, adds=True)
    G.geostats()
    assert 'precip_2000mean' in G.df.columns
    assert 'precip_2001mean' in G.df.columns


def test_directory_mean_matches_fill_values(wgs84_tif_dir, simple_gdf):
    from geostats.main import geostats as GeoStats
    add_raster(wgs84_tif_dir, name='Precip_test')
    G = GeoStats(simple_gdf, measures=['Precip_test'],
                 stats=['mean'], add_stats={}, adds=True)
    G.geostats()
    assert G.df['precip_2000mean'].iloc[0] == pytest.approx(10.0, abs=0.01)
    assert G.df['precip_2001mean'].iloc[0] == pytest.approx(20.0, abs=0.01)


def test_two_custom_measures_coexist(wgs84_tif, tmp_path, simple_gdf):
    """Two independently registered custom rasters must both appear in output."""
    from geostats.main import geostats as GeoStats
    second = tmp_path / 'albedo.tif'
    _write_tif(str(second), epsg=4326, fill=9.0)

    add_raster(wgs84_tif, name='ndvi_test')
    add_raster(str(second), name='albedo_test')

    G = GeoStats(simple_gdf, measures=['ndvi_test', 'albedo_test'],
                 stats=['mean'], add_stats={}, adds=True)
    G.geostats()
    assert 'ndvi_testmean' in G.df.columns
    assert 'albedo_testmean' in G.df.columns
    assert G.df['ndvi_testmean'].iloc[0] == pytest.approx(5.0, abs=0.01)
    assert G.df['albedo_testmean'].iloc[0] == pytest.approx(9.0, abs=0.01)


@pytest.mark.skipif(
    not os.path.isdir(ELEVATION_DATA),
    reason='~/geostats-data/GLOBE/tifs/ not present — skipping custom+builtin coexistence test',
)
def test_custom_and_builtin_measures_coexist(wgs84_tif, simple_gdf):
    """A custom raster alongside a built-in measure must both produce columns."""
    from geostats.main import geostats as GeoStats
    add_raster(wgs84_tif, name='ndvi_test')
    G = GeoStats(simple_gdf, measures=['ndvi_test', 'Elevation'],
                 stats=['mean'], add_stats={}, adds=True)
    G.geostats()
    assert 'ndvi_testmean' in G.df.columns
    assert 'globecylmean' in G.df.columns


# ---------------------------------------------------------------------------
# Robustness tests
# ---------------------------------------------------------------------------

def test_builtin_name_collision_raises(wgs84_tif):
    """Registering a name that matches a built-in measure must raise ValueError."""
    with pytest.raises(ValueError, match="built-in measure"):
        add_raster(wgs84_tif, name='Suitability')


def test_invalid_crs_string_raises(wgs84_tif):
    """Passing an unrecognised crs value must raise ValueError."""
    with pytest.raises(ValueError, match="crs must be"):
        add_raster(wgs84_tif, name='ndvi_test', crs='mercator')


def test_double_registration_same_params_is_idempotent(wgs84_tif):
    """Registering the same name + path + crs twice should be a no-op (no error)."""
    add_raster(wgs84_tif, name='ndvi_test')
    add_raster(wgs84_tif, name='ndvi_test')  # identical — must not raise
    assert wgs84_measures.count('ndvi_test') == 1
    assert main_measures.count('ndvi_test') == 1


def test_double_registration_different_path_raises(wgs84_tif, tmp_path):
    """Re-registering the same name with a different path must raise ValueError."""
    other = tmp_path / 'other.tif'
    _write_tif(str(other), epsg=4326, fill=1.0)
    add_raster(wgs84_tif, name='ndvi_test')
    with pytest.raises(ValueError, match="already registered"):
        add_raster(str(other), name='ndvi_test')


def test_double_registration_different_crs_raises(cea_tif, tmp_path):
    """Re-registering the same name with a different crs must raise ValueError."""
    wgs = tmp_path / 'wgs.tif'
    _write_tif(str(wgs), epsg=4326, fill=1.0)
    add_raster(str(wgs), name='my_test', crs='wgs84')
    with pytest.raises(ValueError, match="already registered"):
        add_raster(cea_tif, name='my_test', crs='cea')


def test_no_duplicate_in_measure_lists(wgs84_tif):
    """Multiple add_raster calls with same args must not duplicate list entries."""
    add_raster(wgs84_tif, name='ndvi_test')
    add_raster(wgs84_tif, name='ndvi_test')
    assert wgs84_measures.count('ndvi_test') == 1
    assert main_measures.count('ndvi_test') == 1
