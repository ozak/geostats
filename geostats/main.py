#!/usr/bin/env python
# coding: utf-8
'''
======================================================
Author:  Ömer Özak, 2014 (ozak at smu.edu)
Website: http://omerozak.com
GitHub:  https://github.com/ozak/
======================================================
Program to compute all geographical statistics for regions given a shapefile
'''
from __future__ import division, print_function
import sys, os
os.environ['USE_PYGEOS'] = '0'
import pandas as pd
from rasterstats import zonal_stats
import numpy as np
from rtree import index
from shapely.ops import unary_union, linemerge
import isounidecode
import geopandas as gp
from geopandas.tools import sjoin, overlay
import pyproj
from pyproj import CRS
import re
from multiprocessing import Pool
from multiprocessing import set_start_method
from rasterstats import zonal_stats
from affine import Affine
import warnings

# Multiprocessing options
try:
    set_start_method('fork')
except:
    print('context has already been set')
n_procs = os.cpu_count()

# Fix paths for data sources
# We include the WGS84 (EPSG:4326) and Lanbert Cylindrical Equal Area (ESRI:54034) projections of the data.
# In the future we may expand to have on-the-fly reprojection (part of to-do list)
# Paths
path = os.getenv('HOME') + '/geostats-data/'
pathmeasures = {'Suitability' : path + '/Ramankutty/tifs/',
                'Suitability2' : path + '/Ramankutty/tifs/',
                'Lights' : path + '/Lights/tifs/',
                'Lights2' : path + '/Lights/tifs/',
                'Elevation' : path + '/GLOBE/tifs/',
                'Elevation2' : path + '/GLOBE/tifs/',
                'Elevation3' : path + '/ETOPO/tifs/',
                'Elevation4' : path + '/ETOPO/tifs/',                
                'RIX' : path + '/RIX/tifs/',
                'RIX2' : path + '/RIX/tifs/',
                'CRU' : path + '/CRU/tifs/',
                'CRU2' : path + '/CRU/tifs/',
                'CropsYield' : path + '/CropsYield/tifs/',
                'CSICrops' : path + '/CSICrops/tifs/',
                'CSI' : path + '/CSI/tifs/',
                'CSI2' : path + '/CSI/tifs/',
                'CSICycle' : path + '/CSICycle/tifs/',
                'CSICycle2' : path + '/CSICycle/tifs/',
                'CSICycleExtra' : path + '/CSICycle/tifs/',
                'CSICycleExtra2' : path + '/CSICycle/tifs/',
                'CSIPlow' : path + '/CSIPlow/tifs/',
                'CSIPlow2' : path + '/CSIPlow/tifs/',
                'Malaria' : path + '/Malaria/tifs/',
                'Malaria2' : path + '/Malaria/tifs/',
                'HMI' : path + '/HMI/tifs/',
                'PopDensAfrica' : path + '/UNEP/tifs/',
                'Popdens' : path + '/GPW/v3/popdens/',
                'Population' : path + '/GPW/v3/population/',
                'Population-GPWv4' : path + '/GPW/v4/population/',
                'PopDensity-GPWv4' : path + '/GPW/v4/popdens/',
                'HYDE' : path + '/HYDE/tifs/',
                'Tsetse' : path + 'Tsetse/tifs/',
                'Ecodiversity' : path + '/Ecodiversity/',
                'EcodiversityWWF' : path + '/Ecodiversity/WWF/',
                'EcodiversityLGM' : path + '/EcodiversityLGM/',
                'Sea100' : path + '/Waters-Coasts/',
                'Coast' : path + '/Waters-Coasts/',
                'Inwater' : path + '/Waters-Coasts/',
                'PerInwater' : path + '/Waters-Coasts/',
                'FluctInwater' : path + '/Waters-Coasts/',
                'Landscan' : path + '/Landscan/tifs/',
                'HLD' : path + '/HLD/tifs/',
                }

# Paths to measures available for computations (if adding new measures one needs to change all places where the lists/dicts are created)
# Grouop measures by type, projection, etc.

# Define and sort main measures
main_measures = [
    'Suitability', 'Suitability2', 'Lights', 'Lights2', 'Elevation','Elevation3', 'RIX', 'Elevation2', 'Elevation4', 'RIX2', 
    'CRU', 'CRU2', 'CSI', 'CSI2', 'CSICycle2', 'CSICycleExtra2', 'CSIPlow2', 
    'CSICrops', 'CropsYield', 'CSICycle', 'CSICycleExtra', 'CSIPlow', 'Malaria', 'Malaria2', 'HMI', 
    'PopDensAfrica', 'Popdens', 'Population', 'Population-GPWv4', 'PopDensity-GPWv4', 'HYDE', 'Tsetse', 'Ecodiversity', 'EcodiversityLGM', 
    'Sea100', 'Coast', 'Inwater', 'PerInwater', 'FluctInwater', 'Landscan', 'HLD'
]
main_measures.sort()

# Define and sort WGS84 measures
wgs84_measures = [
    'Suitability2', 'CRU2', 'CSI2', 'CSICycle2', 'CSICycleExtra2', 'CSIPlow2', 'Lights2', 'Elevation2', 'Elevation4', 'RIX2', 
    'Malaria2', 'PopDensAfrica', 'Popdens', 'Population', 'HYDE', 
    'Population-GPWv4', 'PopDensity-GPWv4', 'Landscan', 'CSICrops', 'CropsYield', 'HLD'
]
wgs84_measures.sort()

# Define and sort shapefile measures
shp_measures = [
    'Ecodiversity', 'EcodiversityLGM', 'Sea100', 'Coast', 'Inwater', 'PerInwater', 'FluctInwater'
]
shp_measures.sort()

# Define and sort ecodiversity measures
ecodiv_measures = ['Ecodiversity', 'EcodiversityLGM']
ecodiv_measures.sort()

cea_measures = list(set(main_measures).difference(set(wgs84_measures)).difference(set(shp_measures)))
cea_measures.sort()

noecodiv_measures = list(set(cea_measures).difference(set(ecodiv_measures)))
noecodiv_measures.sort()

# Identify how many characters need to be adjusted for correct name in each source
namemeasures = {'Suitability' : -4,
                'Suitability2' : -4,
                'Lights' : 7,
                'Lights2' : 7,
                'Elevation' : -4,
                'Elevation2' : -4,
                'Elevation3' : -4,
                'Elevation4' : -4,
                'RIX' : -4,
                'RIX2' : -4,
                'CRU' : -4,
                'CRU2' : -4,
                'CropsYield' : -4,
                'CSICrops' : -4,
                'CSI' : -4,
                'CSI2' : -4,
                'CSICycle' : -4,
                'CSICycle2' : -4,
                'CSICycleExtra' : -4,
                'CSICycleExtra2' : -4,
                'CSIPlow' : -4,
                'CSIPlow2' : -4,
                'Malaria' : -4,
                'Malaria2' : -4,
                'HMI' : -4,
                'PopDensAfrica' : -4,
                'Popdens' : -4,
                'Population' : -4,
                'Population-GPWv4':-4, 
                'PopDensity-GPWv4':-4,
                'HYDE' : -4,
                'Tsetse' : -7,
                'Ecodiversity' : 0,
                'EcodiversityWWF' : 0,
                'EcodiversityLGM' : 0,
                'Sea100' : 0,
                'Coast' : 0,
                'Inwater' : 0,
                'PerInwater' : 0,
                'FluctInwater' : 0,
                'Landscan' : -4,
                'HLD' : 0,
                }

# Functions to perform various operations
# Create vector of polygons and iso-codes
def GISData(myshp, adds=False):
    '''
    This function imports the shape file and converts everyhting to Shapely, and WKT and saves it in a
    pandas dataframe
    Additionally, it computes areas and boundary lengths if adds=True
    Usage:
    GISData(shp=file,[adds=adds])
    '''
    # Shape file
    if isinstance(myshp, str):
        data = gp.GeoDataFrame.from_file(myshp)
    elif isinstance(myshp, gp.GeoDataFrame):
        data = myshp.copy()
    return data

def gini(self):
    '''
    gini(geo)

    Return computed Gini coefficient.
    '''
    if self.count()>1:
        xsort = sorted(self.data[self.mask==False].flatten()) # increasing order
        y = np.cumsum(xsort)
        B = sum(y) / (y[-1] * len(xsort))
        return 1 +  1./len(xsort) - 2*B
    else:
        return 1

# The DataFrames to construct and keep all the additional data
mystats=['min', 'max', 'median', 'mean', 'majority', 'sum','std','count']
addstats={'gini':gini}

# Convert project to equal cylindrical area
p=pyproj.Proj(proj='cea',ellps='WGS84')

#wgs84 = {u'datum': u'WGS84', u'no_defs': True, u'proj': u'longlat'}
#cea = {u'datum': u'WGS84', u'lat_ts': 0, u'lon_0': 0, u'no_defs': True, u'proj': u'cea',u'units': u'm',  u'x_0': 0, u'y_0': 0, 'over':True}
wgs84 = CRS("EPSG:4326")
cea = CRS("ESRI:54034")

# Multiprocessing functions
# Sea
def seaintersects(row):
    if row['sea']!=[]:
        area=row['geometry'].intersection(row['sea']).area/1e6
    else:
        area=0
    return area

def compute_seas(sea100cyl):
    sea100cyl['bbox'] = sea100cyl.geometry.apply(lambda x: x.bounds)
    idxsea = index.Index()
    [idxsea.insert(i, sea100cyl.bbox[i]) for i in sea100cyl.index]
    return idxsea

def process_geometry_seas_df(dfin):
    dfin['bbox'] = dfin.geometry.apply(lambda x: x.bounds)
    dfin['seas'] = dfin.bbox.apply(lambda x: list(idxsea.intersection(x)))
    dfin['sea'] = dfin.seas.apply(lambda x: unary_union(sea100cyl.geometry.iloc[x]))
    dfin['sea100']=dfin.apply(lambda row: seaintersects(row),axis=1)
    dfin['sea100pct']=dfin['sea100']/dfin['area']
    return dfin[['sea100', 'sea100pct']]

def process_geometry_seas(row2):
    row = row2.copy()
    row['bbox'] = row.geometry.bounds
    row['seas'] = list(idxsea.intersection(row.bbox))
    row['sea'] = unary_union(sea100cyl.geometry.iloc[row.seas])
    if row['sea']!=[]:
        row['sea100'] = row['geometry'].intersection(row['sea']).area / 1e6
    else:
        row['sea100'] = 0
    row['sea100pct'] = row['sea100'] / row['area']
    return row[['sea100', 'sea100pct']]

# Coast
# Define a function to calculate coastline length
def coastintersects(row):
    if row['coast'] != []:
        area = row['geometry'].intersection(row['coast']).length / 1e3
    else:
        area = 0
    return area

def process_geometry_coasts_df(dfin):
    dfin['bbox'] = dfin.geometry.apply(lambda x: x.bounds)
    dfin['coasts'] = dfin.bbox.apply(lambda x: list(idxcoast.intersection(x)))
    dfin['coast'] = dfin.coasts.apply(lambda x: unary_union(coastcyl.iloc[x].geometry.values))
    # Apply the coastintersects function to calculate coastline length for each row
    dfin['coastlen'] = dfin.apply(lambda row: coastintersects(row), axis=1)
    return dfin[['coastlen']]


# Ecodiversity
# Parallelization
def ecodiv_parallel(args):
    row, idxeco = args
    if row['ecos'] != []:
        tarea = row['geometry'].area
        polys = ecological.iloc[row.ecos].dissolve(by=['OBJECTID', 'WWF_MHTNAM'])
        share = polys.apply(lambda x: row.geometry.buffer(0).intersection(x.geometry.buffer(0)).area, axis=1) / tarea
        sharesq = share**2
        sharesq2 = share * (1 - 2 * share)**2
        sharesq = sharesq.sum()
        sharesq2 = sharesq2.sum()
        ecodiversity = 1 - sharesq
        ecopolarization = 1 - sharesq2
    else:
        ecodiversity = 0
        ecopolarization = 0
    return [ecodiversity, ecopolarization]

def ecodivwwf_parallel(args):
    row, idxeco = args
    if row['ecos'] != []:
        tarea = row['geometry'].area
        polys = ecological.iloc[row.ecos].dissolve(by=['ECO_NAME', 'ECO_NUM', 'ECO_ID', 'ECO_SYM', 'eco_code'])
        share = polys.apply(lambda x: row.geometry.buffer(0).intersection(x.geometry.buffer(0)).area, axis=1) / tarea
        sharesq = share**2
        sharesq2 = share * (1 - 2 * share)**2
        sharesq = sharesq.sum()
        sharesq2 = sharesq2.sum()
        ecodiversity = 1 - sharesq
        ecopolarization = 1 - sharesq2
    else:
        ecodiversity = 0
        ecopolarization = 0
    return [ecodiversity, ecopolarization]

def ecodivLGM_parallel(args):
    row, idxeco = args
    if row['ecos'] != []:
        tarea = row['geometry'].area
        polys = ecologicalLGM.iloc[row.ecos]
        share = polys.apply(lambda x: row.geometry.buffer(0).intersection(x.geometry.buffer(0)).area, axis=1) / tarea
        sharesq = share**2
        sharesq2 = share * (1 - 2 * share)**2
        sharesq = sharesq.sum()
        sharesq2 = sharesq2.sum()
        ecodiversity = 1 - sharesq
        ecopolarization = 1 - sharesq2
    else:
        ecodiversity = 0
        ecopolarization = 0
    return [ecodiversity, ecopolarization]


# Geostats rasters
def zs_parallel2(x):
    [n, shape, rpathi, stats, copy_properties, add_stats, all_touched] = x
    return zonal_stats(shape, rpathi, stats=stats, copy_properties=copy_properties,
                    add_stats=add_stats, all_touched=all_touched)[0]

def geostats_MP(shapefile, measures = ['All'], stats=mystats, copy_properties=True,
             add_stats=addstats, all_touched=True, adds=True, correct_cea=True, correct_wgs=True):
    '''Compute Zonalstats for measures using MP'''
    df = GISData(myshp = shapefile, adds=adds)
    if df.crs!=wgs84:
        dfnocyl = df.copy().to_crs(wgs84)
    else:
        dfnocyl = df.copy()
    if df.crs!=cea:
        df = df.to_crs(cea)
    if correct_cea:
        df.loc[df.is_valid==False, 'geometry'] = df.loc[df.is_valid==False, 'geometry'].apply(lambda x: x.buffer(0))
    if correct_wgs:
        dfnocyl.loc[dfnocyl.is_valid==False, 'geometry'] = dfnocyl.loc[dfnocyl.is_valid==False, 'geometry'].apply(lambda x: x.buffer(0))
    if adds:
        if df.geom_type.values[0].find('Polygon')!=-1:
            df['area']=df.area/1e6
            df['boundary']=df.boundary.length/1e3
        if df.geom_type.values[0].find('Point')!=-1:
            df['boundary']=df.boundary.length/1e3
    # Lat lon
    df['X']=df.centroid.apply(lambda c: c.coords.xy[0][0])
    df['Y']=df.centroid.apply(lambda c: c.coords.xy[1][0])
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning)
        df['lon']=dfnocyl.centroid.apply(lambda c: c.coords.xy[0][0])
        df['lat']=dfnocyl.centroid.apply(lambda c: c.coords.xy[1][0])

    if np.sum([m.upper()=='ALL' for m in measures])>0:
        measures = list(main_measures)
    else:
        measures = list(measures)
    pass

    for measure in measures:
        print('Computing Stats for:',measure)
        rpath = pathmeasures[measure]
        if measure in cea_measures:
            dfin = df
            dfin.crs = df.crs
            mytiffiles=[mytiffile for mytiffile in os.listdir(rpath) if mytiffile.endswith('cyl.tif')]
            if measure == 'HMI':
                mytiffiles=[mytiffile for mytiffile in os.listdir(rpath) if mytiffile.endswith('.tif')]
            elif measure == 'CSICycle':
                mytiffiles=[mytiffile for mytiffile in os.listdir(rpath) if mytiffile.endswith('cyl.tif') and (mytiffile.find('1500')!=-1 or mytiffile.find('dif')!=-1)]
            elif measure == 'CSICycleExtra':
                mytiffiles=[mytiffile for mytiffile in os.listdir(rpath) if mytiffile.endswith('cyl.tif') and (mytiffile.find('1500')==-1 and mytiffile.find('dif')==-1)]
        elif measure in wgs84_measures:
            dfin = dfnocyl
            dfin.crs = dfnocyl.crs
            mytiffiles=[mytiffile for mytiffile in os.listdir(rpath) if mytiffile.endswith('.tif') and mytiffile.find('cyl.tif')==-1]
            if measure=='HLD':
                mytiffiles=[mytiffile for mytiffile in mytiffiles if mytiffile.find('Harmonized')!=-1]
            elif measure == 'CSICycle2':
                mytiffiles=[mytiffile for mytiffile in os.listdir(rpath) if mytiffile.endswith('cyl.tif')==False and (mytiffile.find('1500')!=-1 or mytiffile.find('dif')!=-1)]
            elif measure == 'CSICycleExtra2':
                mytiffiles=[mytiffile for mytiffile in os.listdir(rpath) if mytiffile.endswith('cyl.tif')==False and (mytiffile.find('1500')==-1 and mytiffile.find('dif')==-1)]
        else:
            dfin = df
            dfin.crs = df.crs
        if measure not in shp_measures:
            mytiffiles.sort()
            for i in mytiffiles:
                myvar=i[:namemeasures[measure]]
                if measure=='CSICycleWater' or measure=='CSIWater':
                    myvar = 'w' +  myvar
                elif measure=='CSICrops':
                    myvar = 'CSI' + myvar
                elif measure=='Lights':
                    myvar = myvar + 'cyl'
                elif measure=='HLD':
                    myvar = 'HLD' + re.findall(r'\d+', i)[0]
                # Multiprocessing of stats
                print(myvar)
                n = range(0, len(dfin))
                mypool = [[n, dfin.loc[n, ['geometry']], rpath + i, stats, copy_properties, add_stats, all_touched] for n in range(0, len(dfin))]
                pool = Pool(n_procs)
                results = pool.map(zs_parallel2, mypool)
                pool.close()
                df2 = pd.DataFrame(data=results)
                mypstats = [myvar + j for j in df2.columns.values]
                df2.columns = mypstats
                df2.fillna(0,inplace=True)
                df = pd.merge(df, df2, left_index=True, right_index=True, how='outer')
    return df

class geostats(object):
    """
    Implements the geostats Class, which will be used to construct the computations, etc. and hold everything in a unique object.

    Usage:
        Create geostats object:
            A = geostats(shapefile, measures = ['All'], stats=mystats, copy_properties=True,
                         add_stats=addstats, all_touched=True, adds=True, correct_cea=True, correct_wgs=True)

    where
            shapefile: string with the location of the shapefile to use or a geo-pandas data frame
            measures: List of measures to compute among those available, default ['All']
                      Possible measures :
            stats: List of statistics to be computed, determined by raster_stats module, see help(raster_stats)
                   default ['min', 'max', 'median', 'mean', 'majority', 'sum','std','count']
            copy_properties: Boolean, include or not all properties from the dataframe/shapefile
                             default = True
            add_stats: Dictionary of additional statistics to be computed, see help(raster_stats)
            all_touched: Boolean, default True, whether to include all cells touched by polygon in computations
            adds: Boolean, default True, compute area, boundary, centroid
    """
    # Initialize Object
    def __init__(self, shapefile, measures = ['All'], stats=mystats, copy_properties=True,
                 add_stats=addstats, all_touched=True, adds=True, correct_cea=True, correct_wgs=True, **kwargs):
        '''
        Initialize
        '''
        super(geostats, self).__init__()
        self.shapefile = shapefile
        self.df = GISData(myshp = shapefile, adds=adds)
        if self.df.crs!=wgs84:
            self.dfnocyl = self.df.copy().to_crs(wgs84)
        else:
            self.dfnocyl = self.df.copy()
        if self.df.crs!=cea:
            self.df = self.df.to_crs(cea)
        if correct_cea:
            self.df.loc[self.df.is_valid==False, 'geometry'] = self.df.loc[self.df.is_valid==False, 'geometry'].apply(lambda x: x.buffer(0))
        if correct_wgs:
            self.dfnocyl.loc[self.dfnocyl.is_valid==False, 'geometry'] = self.dfnocyl.loc[self.dfnocyl.is_valid==False, 'geometry'].apply(lambda x: x.buffer(0))
        if adds:
            if self.df.geom_type.values[0].find('Polygon')!=-1:
                self.df['area']=self.df.area/1e6
                self.df['boundary']=self.df.boundary.length/1e3
            if self.df.geom_type.values[0].find('Point')!=-1:
                self.df['boundary']=self.df.boundary.length/1e3
        # Lat lon
        self.df['X']=self.df.centroid.apply(lambda c: c.coords.xy[0][0])
        self.df['Y']=self.df.centroid.apply(lambda c: c.coords.xy[1][0])
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning)
            self.df['lon']=self.dfnocyl.centroid.apply(lambda c: c.coords.xy[0][0])
            self.df['lat']=self.dfnocyl.centroid.apply(lambda c: c.coords.xy[1][0])
        self.stats = stats
        self.copy_properties = copy_properties
        self.add_stats = add_stats
        self.all_touched = all_touched
        if np.sum([m.upper()=='ALL' for m in measures])>0:
            self.measures = list(main_measures)
        else:
            self.measures = list(measures)
        pass

    def geostats(self, **kwargs):
        '''
        Compute the zonal stats
        '''
        for measure in self.measures:
            print('Computing Stats for:',measure)
            rpath = pathmeasures[measure]
            if measure in cea_measures:
                dfin = self.df
                dfin.crs = self.df.crs
                mytiffiles=[mytiffile for mytiffile in os.listdir(rpath) if mytiffile.endswith('cyl.tif')]
                if measure == 'HMI':
                    mytiffiles=[mytiffile for mytiffile in os.listdir(rpath) if mytiffile.endswith('.tif')]
                elif measure == 'CSI':
                    mytiffiles=[mytiffile for mytiffile in os.listdir(rpath) if mytiffile.endswith('cyl.tif') and (mytiffile.find('1500')!=-1 or mytiffile.find('dif')!=-1)]
                elif measure == 'CSICycle':
                    mytiffiles=[mytiffile for mytiffile in os.listdir(rpath) if mytiffile.endswith('cyl.tif') and (mytiffile.find('1500')!=-1 or mytiffile.find('dif')!=-1)]
                elif measure == 'CSICycleExtra':
                    mytiffiles=[mytiffile for mytiffile in os.listdir(rpath) if mytiffile.endswith('cyl.tif') and (mytiffile.find('1500')==-1 and mytiffile.find('dif')==-1)]
            elif measure in wgs84_measures:
                dfin = self.dfnocyl
                dfin.crs = self.dfnocyl.crs
                mytiffiles=[mytiffile for mytiffile in os.listdir(rpath) if mytiffile.endswith('.tif') and mytiffile.find('cyl.tif')==-1]
                if measure=='HLD':
                    mytiffiles=[mytiffile for mytiffile in mytiffiles if mytiffile.find('Harmonized')!=-1]
                elif measure == 'CSI2':
                    mytiffiles=[mytiffile for mytiffile in os.listdir(rpath) if mytiffile.endswith('cyl.tif')==False and (mytiffile.find('1500')!=-1 or mytiffile.find('dif')!=-1)]
                elif measure == 'CSICycle2':
                    mytiffiles=[mytiffile for mytiffile in os.listdir(rpath) if mytiffile.endswith('cyl.tif')==False and (mytiffile.find('1500')!=-1 or mytiffile.find('dif')!=-1)]
                elif measure == 'CSICycleExtra2':
                    mytiffiles=[mytiffile for mytiffile in os.listdir(rpath) if mytiffile.endswith('cyl.tif')==False and (mytiffile.find('1500')==-1 and mytiffile.find('dif')==-1)]
            else:
                dfin = self.df
                dfin.crs = self.df.crs
            if measure not in shp_measures:
                mytiffiles.sort()
                for i in mytiffiles:
                    myvar=i[:namemeasures[measure]]
                    if measure=='CSICycleWater' or measure=='CSIWater':
                        myvar = 'w' +  myvar
                    elif measure=='CSICrops':
                        myvar = 'CSI' + myvar
                    elif measure=='Lights':
                        myvar = myvar + 'cyl'
                    elif measure=='GEcon':
                        myvar = 'GE' + myvar
                    elif measure=='Population-GPWv4':
                        myvar = myvar.replace('population_count_', 'popc').replace('adjusted_to_2015_unwpp_country_totals_', 'adj').replace('rev11_', '').replace('_30_sec', '').replace('_', '')
                    elif measure=='PopDensity-GPWv4':
                        myvar = myvar.replace('population_density_', 'popd').replace('adjusted_to_2015_unwpp_country_totals_', 'adj').replace('rev11_', '').replace('_30_sec', '').replace('_', '')
                    elif measure=='Landscan':
                        myvar = myvar.replace('landscan', 'ls').replace('-global-', '')
                    elif measure=='HLD':
                        myvar = 'HLD' + re.findall(r'\d+', i)[0]
                    # Multiprocessing of stats
                    print(myvar)
                    n = range(0, len(dfin))
                    mypool = [[n, dfin.loc[n, ['geometry']], rpath + i, self.stats, self.copy_properties, self.add_stats, self.all_touched] for n in range(0, len(dfin))]
                    pool = Pool(n_procs)
                    results = pool.map(zs_parallel2, mypool)
                    pool.close()
                    df2 = pd.DataFrame(data=results)
                    mypstats = [myvar + j for j in df2.columns.values]
                    df2.columns = mypstats
                    df2.fillna(0,inplace=True)
                    self.df = pd.merge(self.df, df2, left_index=True, right_index=True, how='outer')
            else:
                if measure in ['Sea100', 'Inwater', 'PerInwater', 'FluctInwater']:
                    if measure=='Sea100':
                        seafile = 'coastlcyl_buffer100.shp'
                        seavar = 'sea100'
                    elif measure=='Inwater':
                        seafile = 'inwateracyl100.shp'
                        seavar = 'inwater100'
                    elif measure=='PerInwater':
                        seafile = 'per_inwateracyl100.shp'
                        seavar = 'inwater100per'
                    elif measure=='FluctInwater':
                        seafile = 'fluct_inwateracyl100.shp'
                        seavar = 'inwater100fluct'

                    # Load GIS data from shapefile into the sea100cyl DataFrame
                    global sea100cyl
                    sea100cyl = GISData(myshp=pathmeasures[measure] + seafile)

                    # Calculate bounding boxes for the sea100cyl geometry and add them as a column
                    sea100cyl['bbox'] = sea100cyl.geometry.apply(lambda x: x.bounds)

                    # Copy geometry and area columns from self.df to a new DataFrame called dfin
                    dfin = self.df[['geometry', 'area']].copy()

                    # Set the CRS of dfin to match self.df
                    dfin.crs = self.df.crs

                    # Calculate bounding boxes for dfin's geometry and add them as a column
                    #dfin['bbox'] = dfin.geometry.apply(lambda x: x.bounds)

                    # Create a spatial index for sea100cyl
                    global idxsea
                    idxsea = index.Index()
                    [idxsea.insert(sea100cyl.index[i], sea100cyl.bbox[i]) for i in range(sea100cyl.shape[0])]

                    # Determine which seas each dfin bounding box intersects with
                    #dfin['seas'] = dfin.bbox.apply(lambda x: list(idxsea.intersection(x)))

                    # Combine the geometries of intersecting seas using unary union
                    #dfin['sea'] = dfin.seas.apply(lambda x: unary_union(sea100cyl.geometry.iloc[x]))

                    with Pool(processes=n_procs) as pool:
                        # Process the geometry using parallel processing
                        dfin_chunks = np.array_split(dfin, n_procs)
                        processed_dfin = pool.map(process_geometry_seas_df, [chunk for chunk in dfin_chunks])

                    # Combine the results if needed
                    processed_dfin = pd.concat(processed_dfin, ignore_index=True)
                    processed_dfin.columns = [seavar, seavar + 'pct']

                    ## Apply the processed results to dfin
                    #processed_dfin = pd.DataFrame(processed_rows, columns=['sea100', 'sea100pct'], index=dfin.index)
                    dfin = pd.concat([dfin, processed_dfin], axis=1)

                    # Merge the sea100 and sea100pct columns from dfin into self.df
                    self.df = pd.merge(self.df, dfin[[seavar, seavar + 'pct']], left_index=True, right_index=True, how='outer')
                # Coast
                elif measure == 'Coast':
                    # Load GIS data from shapefile into the coastcyl DataFrame
                    global coastcyl
                    coastcyl = GISData(myshp=pathmeasures[measure] + 'coastcyl.shp')

                    # Calculate bounding boxes for the coastcyl geometry and add them as a column
                    coastcyl['bbox'] = coastcyl.geometry.apply(lambda x: x.bounds)

                    # Copy geometry and area columns from self.df to a new DataFrame called dfin
                    dfin = self.df[['geometry', 'area']].copy()

                    # Set the CRS of dfin to match self.df
                    dfin.crs = self.df.crs

                    # Calculate bounding boxes for dfin's geometry and add them as a column
                    #dfin['bbox'] = dfin.geometry.apply(lambda x: x.bounds)

                    # Create a spatial index for coastcyl
                    global idxcoast
                    idxcoast = index.Index()
                    [idxcoast.insert(coastcyl.index[i], coastcyl.bbox[i]) for i in range(coastcyl.shape[0])]

                    # Determine which coasts each dfin bounding box intersects with
                    #dfin['coasts'] = dfin.bbox.apply(lambda x: list(idxcoast.intersection(x)))

                    # Combine the geometries of intersecting coasts using unary union
                    #dfin['coast'] = dfin.coasts.apply(lambda x: unary_union(coastcyl.iloc[x].geometry.values))

                    with Pool(processes=n_procs) as pool:
                        # Process the geometry using parallel processing
                        dfin_chunks = np.array_split(dfin, n_procs)
                        processed_dfin = pool.map(process_geometry_coasts_df, [chunk for chunk in dfin_chunks])

                    # Combine the results if needed
                    processed_dfin = pd.concat(processed_dfin, ignore_index=True)
                    processed_dfin.columns = ['coastlen']

                    ## Apply the processed results to dfin
                    #processed_dfin = pd.DataFrame(processed_rows, columns=['sea100', 'sea100pct'], index=dfin.index)
                    dfin = pd.concat([dfin, processed_dfin], axis=1)

                    # Merge the coastlen column from dfin into self.df
                    self.df = pd.merge(self.df, dfin[['coastlen']], left_index=True, right_index=True, how='outer')

                elif measure == 'Ecodiversity':
                    # Load ecological GIS data from shapefile into the ecological DataFrame
                    global ecological
                    ecological = GISData(myshp=pathmeasures[measure] + 'Biomes_Worldcyl.shp')

                    # Convert Multipolygon into polygons for speeding up computations
                    ecological = ecological.explode(index_parts=True).reset_index()

                    # Calculate bounding boxes for ecological geometry and add them as a column
                    ecological['bbox'] = ecological.geometry.apply(lambda x: x.bounds)

                    # Create spatial index for ecological geometries
                    idxeco = index.Index()
                    [idxeco.insert(ecological.index[i], ecological.bbox[i]) for i in range(ecological.shape[0])]

                    # Copy geometry from self.df to dfin
                    dfin = self.df[['geometry']].copy()

                    # Set the CRS of dfin to match self.df
                    dfin.crs = self.df.crs

                    # Calculate bounding boxes for dfin's geometry and add them as a column
                    dfin['bbox'] = dfin.geometry.apply(lambda x: x.bounds)

                    # Determine which ecological geometries each dfin bounding box intersects with
                    dfin['ecos'] = dfin.bbox.apply(lambda x: list(idxeco.intersection(x)))

                    with Pool(processes=n_procs) as pool:
                        processed_rows = pool.map(ecodiv_parallel, [(row, idxeco) for idx, row in dfin.iterrows()])

                    # Apply the processed results to dfin
                    processed_dfin = pd.DataFrame(processed_rows, columns=['ecodiversity', 'ecopolarization'], index=dfin.index)
                    dfin = pd.concat([dfin, processed_dfin], axis=1)

                    # Merge ecodiversity and ecopolarization columns from dfin into self.df
                    self.df = pd.merge(self.df, dfin[['ecodiversity', 'ecopolarization']], left_index=True, right_index=True, how='outer')
                elif measure == 'EcodiversityLGM':
                    # Load ecological GIS data from shapefile into the ecological DataFrame
                    global ecologicalLGM
                    ecologicalLGM = GISData(myshp=pathmeasures[measure] + 'world_cutcyl.shp')

                    # Calculate bounding boxes for ecological geometry and add them as a column
                    ecologicalLGM['bbox'] = ecologicalLGM.geometry.apply(lambda x: x.bounds)

                    # Create spatial index for ecological geometries
                    idxeco = index.Index()
                    [idxeco.insert(ecologicalLGM.index[i], ecologicalLGM.bbox[i]) for i in range(ecologicalLGM.shape[0])]

                    # Copy geometry from self.df to dfin
                    dfin = self.df[['geometry']].copy()

                    # Set the CRS of dfin to match self.df
                    dfin.crs = self.df.crs

                    # Calculate bounding boxes for dfin's geometry and add them as a column
                    dfin['bbox'] = dfin.geometry.apply(lambda x: x.bounds)

                    # Determine which ecological geometries each dfin bounding box intersects with
                    dfin['ecos'] = dfin.bbox.apply(lambda x: list(idxeco.intersection(x)))

                    with Pool(processes=n_procs) as pool:
                        processed_rows = pool.map(ecodivLGM_parallel, [(row, idxeco) for idx, row in dfin.iterrows()])

                    # Apply the processed results to dfin
                    processed_dfin = pd.DataFrame(processed_rows, columns=['ecodiversityLGM', 'ecopolarizationLGM'], index=dfin.index)
                    dfin = pd.concat([dfin, processed_dfin], axis=1)

                    self.df = pd.merge(self.df, dfin[['ecodiversityLGM','ecopolarizationLGM']], left_index=True, right_index=True, how='outer')
        pass
