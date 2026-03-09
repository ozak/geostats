"""
Regression tests using real geographic data (Natural Earth countries).

 - test_builtin_measure_*  : verify existing built-in pipeline is intact.
   Skipped automatically when ~/geostats-data/ is not present.

 - test_add_raster_*       : verify add_raster() works with real country
   geometries (Natural Earth) and a synthetic raster over the same extent.
"""
import os
import io
import zipfile
import numpy as np
import pytest
import requests
import rasterio
from rasterio.transform import from_bounds
import geopandas as gpd

from geostats.main import (
    add_raster,
    pathmeasures, wgs84_measures, cea_measures, main_measures,
    namemeasures, _user_single_files,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

GEOSTATS_DATA = os.path.join(os.path.expanduser('~'), 'geostats-data')

NATURALEARTH_URL = (
    'https://naciscdn.org/naturalearth/110m/cultural/'
    'ne_110m_admin_0_countries.zip'
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope='session')
def ne_countries(tmp_path_factory):
    """Download Natural Earth 110m countries once per session."""
    dest = tmp_path_factory.mktemp('ne_data')
    try:
        r = requests.get(NATURALEARTH_URL, timeout=30)
        r.raise_for_status()
    except Exception as exc:
        pytest.skip(f'Could not download Natural Earth data: {exc}')
    with zipfile.ZipFile(io.BytesIO(r.content)) as zf:
        zf.extractall(dest)
    shp = next(dest.glob('*.shp'), None)
    if shp is None:
        pytest.skip('Natural Earth shapefile not found after extraction')
    return gpd.read_file(str(shp))


@pytest.fixture(scope='session')
def global_wgs84_tif(tmp_path_factory):
    """A global WGS84 raster (fill=7.0) at 1-degree resolution."""
    p = tmp_path_factory.mktemp('raster') / 'custom_global.tif'
    transform = from_bounds(-180, -90, 180, 90, 360, 180)
    data = np.full((180, 360), 7.0, dtype=np.float32)
    with rasterio.open(
        str(p), 'w', driver='GTiff',
        height=180, width=360,
        count=1, dtype='float32',
        crs='EPSG:4326',
        transform=transform,
    ) as dst:
        dst.write(data, 1)
    return str(p)


@pytest.fixture(autouse=True)
def _cleanup():
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
# Regression tests — built-in measures
# ---------------------------------------------------------------------------

@pytest.mark.skipif(
    not os.path.isdir(GEOSTATS_DATA),
    reason='~/geostats-data/ not present — skipping built-in measure regression test',
)
def test_builtin_elevation_produces_columns(ne_countries):
    """Elevation (CEA raster) should produce stat columns for all countries."""
    from geostats.main import geostats as GeoStats
    # Reset index: main.py uses range(len(df)) with .loc, so index must be 0-based
    subset = ne_countries[ne_countries['CONTINENT'] == 'South America'].copy()
    G = GeoStats(subset, measures=['Elevation'],
                 stats=['mean', 'min', 'max'], add_stats={}, adds=True)
    G.geostats()
    for col in ['globecylmean', 'globecylmin', 'globecylmax']:
        assert col in G.df.columns, f"Missing column: {col}"


@pytest.mark.skipif(
    not os.path.isdir(GEOSTATS_DATA),
    reason='~/geostats-data/ not present — skipping built-in measure regression test',
)
def test_builtin_elevation_values_are_plausible(ne_countries):
    """Country mean elevations should be non-negative finite numbers."""
    from geostats.main import geostats as GeoStats
    subset = ne_countries[ne_countries['CONTINENT'] == 'South America'].copy()
    G = GeoStats(subset, measures=['Elevation'],
                 stats=['mean'], add_stats={}, adds=True)
    G.geostats()
    col = 'globecylmean'
    vals = G.df[col].dropna()
    assert len(vals) > 0
    assert (vals >= 0).all(), 'Negative elevation mean — unexpected'
    assert np.isfinite(vals).all()


# ---------------------------------------------------------------------------
# add_raster tests with real country geometries
# ---------------------------------------------------------------------------

def test_add_raster_real_countries_columns(ne_countries, global_wgs84_tif):
    """add_raster() should produce one column per stat for every country row."""
    from geostats.main import geostats as GeoStats
    subset = ne_countries[ne_countries['CONTINENT'] == 'Africa'].copy()
    add_raster(global_wgs84_tif, name='custom')
    G = GeoStats(subset, measures=['custom'],
                 stats=['mean', 'count'], add_stats={}, adds=True)
    G.geostats()
    assert 'custommean' in G.df.columns
    assert 'customcount' in G.df.columns
    assert len(G.df) == len(subset)


def test_add_raster_real_countries_mean_value(ne_countries, global_wgs84_tif):
    """Every country should get mean == fill value (7.0) from the global raster."""
    from geostats.main import geostats as GeoStats
    subset = ne_countries[ne_countries['CONTINENT'] == 'Africa'].copy()
    add_raster(global_wgs84_tif, name='custom')
    G = GeoStats(subset, measures=['custom'],
                 stats=['mean'], add_stats={}, adds=True)
    G.geostats()
    vals = G.df['custommean'].dropna()
    assert len(vals) > 0
    assert np.allclose(vals.values, 7.0, atol=0.01)


def test_add_raster_real_countries_count_positive(ne_countries, global_wgs84_tif):
    """Every country polygon should intersect at least one raster cell."""
    from geostats.main import geostats as GeoStats
    subset = ne_countries[ne_countries['CONTINENT'] == 'Africa'].copy()
    add_raster(global_wgs84_tif, name='custom')
    G = GeoStats(subset, measures=['custom'],
                 stats=['count'], add_stats={}, adds=True)
    G.geostats()
    assert (G.df['customcount'].fillna(0) > 0).all()
