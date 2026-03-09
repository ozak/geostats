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
    namemeasures, _user_single_files,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_tif(path, epsg, fill=5.0, bounds=(-20.0, -20.0, 20.0, 20.0), size=(40, 40)):
    """Write a constant-value single-band GeoTIFF."""
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


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def wgs84_tif(tmp_path):
    p = tmp_path / 'ndvi.tif'
    _write_tif(str(p), epsg=4326, fill=5.0)
    return str(p)


@pytest.fixture
def non_wgs84_tif(tmp_path):
    """A raster in Web Mercator (EPSG:3857) — will be classified as 'cea'."""
    p = tmp_path / 'ruggedness.tif'
    _write_tif(str(p), epsg=3857, fill=3.0,
               bounds=(-2e6, -2e6, 2e6, 2e6))
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
            for lst in (wgs84_measures, cea_measures, main_measures):
                if name in lst:
                    lst.remove(name)


# ---------------------------------------------------------------------------
# Registration tests
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


def test_non_wgs84_single_file_goes_to_cea_measures(non_wgs84_tif):
    add_raster(non_wgs84_tif, name='rug_test')
    assert 'rug_test' in cea_measures
    assert 'rug_test' not in wgs84_measures


def test_directory_goes_to_wgs84_measures(wgs84_tif_dir):
    add_raster(wgs84_tif_dir, name='Precip_test')
    assert 'Precip_test' in wgs84_measures


def test_directory_namemeasures_is_minus4(wgs84_tif_dir):
    add_raster(wgs84_tif_dir, name='Precip_test')
    assert namemeasures['Precip_test'] == -4


def test_directory_not_in_user_single_files(wgs84_tif_dir):
    add_raster(wgs84_tif_dir, name='Precip_test')
    assert 'Precip_test' not in _user_single_files


def test_explicit_crs_overrides_autodetect(wgs84_tif):
    """Passing crs='cea' should override the WGS84 autodetect."""
    add_raster(wgs84_tif, name='ndvi_cea', crs='cea')
    assert 'ndvi_cea' in cea_measures
    assert 'ndvi_cea' not in wgs84_measures


def test_empty_directory_raises(tmp_path):
    with pytest.raises(ValueError, match='No .tif files'):
        add_raster(str(tmp_path), name='bad')


# ---------------------------------------------------------------------------
# Computation tests
# ---------------------------------------------------------------------------

def test_single_file_produces_named_columns(wgs84_tif, simple_gdf):
    from geostats.main import geostats as GeoStats
    add_raster(wgs84_tif, name='ndvi_test')
    G = GeoStats(simple_gdf, measures=['ndvi_test'],
                 stats=['mean', 'count'], add_stats={}, adds=True)
    G.geostats()
    # rasterstats returns stat keys without separator: myvar + 'mean' -> 'ndvi_testmean'
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


def test_custom_and_builtin_measures_coexist(wgs84_tif, simple_gdf):
    """A custom raster in the measures list shouldn't break built-in measure lookup."""
    from geostats.main import geostats as GeoStats
    add_raster(wgs84_tif, name='ndvi_test')
    G = GeoStats(simple_gdf, measures=['ndvi_test'],
                 stats=['mean'], add_stats={}, adds=True)
    G.geostats()
    assert 'ndvi_testmean' in G.df.columns
