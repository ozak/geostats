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
#from osgeo import gdal, gdalnumeric, ogr, osr
#from gdalconst import *
#from PIL import Image, ImageDraw
#import tarfile
#import gzip
#import shutil, glob
import pandas as pd
#from pyGDsandbox.dataIO import df2dbf, dbf2df
from rasterstats import zonal_stats
#import pysal as ps
#import shapely
#from shapely.wkt import loads, dumps
#from pysal.contrib import shapely_ext
import numpy as np
from rtree import index
from shapely.ops import unary_union, linemerge
import isounidecode
import geopandas as gp
from geopandas.tools import sjoin, overlay
import pyproj
from pyproj import CRS
#import hmi
#import georasters as gr
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

# Determine drive where other data is located
try:
    os.listdir("/Volumes/Mirror RAID")
    drv="/Volumes/Backup/MacPro/Mirror RAID/"
except:
    drv="/Users/bizcocho/Desktop/"

# Fix paths of Ethnic, country and coastline Shapefiles
# Although we include the WGS84 versions, we will only use the projected ones, in order to have everything
if sys.platform=='darwin':
    path = os.getenv('HOME') + '/LatexMeGD/CulturalDistance/'
    pathcyl=path + '/data/GIS/cyl/'
    pathcntry=path + 'data/GIS/GMI/'
    pathcntrycyl=path + 'data/GIS/GMI/cyl/'
    pathout=path + 'data/GIS/data/'
    pathgis=path + 'data/GIS/'
    pathmeasures = {'Suitability' : pathgis + 'Ramankutty/',
                    'Lights' : '/Volumes/My Book/GIS/Lights/tifs/',
                    'Lights2' : '/Volumes/My Book/GIS/Lights/tifs/',
                    'Elevation' : drv + '/Geographical_Index/GIS/GLOBE/',
                    'Elevation2' : drv + '/Geographical_Index/GIS/GLOBE/',
                    'RIX' : drv + '/Geographical_Index/GIS/RIX/',
                    'RIX2' : drv + '/Geographical_Index/GIS/RIX/',
                    'CRU' : drv + '/Geographical_Index/DataOthers/CRU/tif/',
                    'CropsYield' : drv + '/Geographical_Index/DataOthers/FAO/GAEZ/Agroclimatic/yield/rain_fed/',
                    'CSICrops' : drv + '/Geographical_Index/DataOthers/FAO/GAEZ/Agroclimatic/energyyield/rain_fed/',
                    'CSI' : drv + '/Geographical_Index/DataOthers/FAO/GAEZ/Agroclimatic/bestcropcal/CSI/',
                    'CSIWater' : drv + '/Geographical_Index/DataOthers/FAO/GAEZ/Agroclimatic/bestcropcal/CSI_water/',
                    'CSICycle' : drv + '/Geographical_Index/DataOthers/FAO/GAEZ/Agroclimatic/bestcropcal/latest/',
                    'CSICycleExtra' : drv + '/Geographical_Index/DataOthers/FAO/GAEZ/Agroclimatic/bestcropcal/latest/',
                    'CSICycleWater' : drv + '/Geographical_Index/DataOthers/FAO/GAEZ/Agroclimatic/bestcropcal/updated_water/',
                    'CSIPlow' : drv + '/Geographical_Index/DataOthers/FAO/GAEZ/Agroclimatic/bestcropcal/plow/',
                    'Malaria' : pathgis + 'Malaria/',
                    'HMI' : drv + '/Geographical_Index/GIS/HMI_Diff/',
                    'PopDensAfrica' : '/Volumes/Macintosh HD 2/GIS/UNEP/',
                    'Popdens' : '/Volumes/My Book/GIS/GPW/popdens/',
                    'Population' : '/Volumes/My Book/GIS/GPW/population/',
                    'HYDE' : '/Volumes/Backup/MacPro/Macintosh HD 2/GIS/HYDE/hyde31_final/tif/',
                    'Tsetse' : pathgis + 'Tsetse/',
                    'Ecodiversity' : pathgis + 'Ecological_Zones/',
                    'EcodiversityLGM' : pathgis + 'PaleoVegetation/',
                    'Sea100' : pathcntrycyl,
                    'Coast' : pathcntrycyl,
                    'Inwater' : pathcntrycyl,
                    'PerInwater' : pathcntrycyl,
                    'FluctInwater' : pathcntrycyl,
                    'Landscan' : '/Volumes/Backup/MacPro/Macintosh HD 2/GIS/Landscan/',
                    'GEcon' : '/Volumes/Backup/MacPro/Macintosh HD 2/GIS/GEconRaster/',
                    'HLD' : '/Volumes/My Book/GIS/Harmonized Light Data/light/',
                    }
elif sys.platform.startswith('linux'):
    path = os.getenv('HOME') + '/LatexMe/CulturalDistance/'
    pathcyl=path + '/data/GIS/cyl/'
    pathcntry=path + 'data/GIS/GMI/'
    pathcntrycyl=path + 'data/GIS/GMI/cyl/'
    pathout=path + 'data/GIS/data/'
    pathgis=path + 'data/GIS/'
    drv = os.getenv('HOME') + '/DataOthers/'
    pathmeasures = {'Suitability' : pathgis + 'Ramankutty/',
                    'Lights' : drv + '/Lights/tifs/',
                    'Lights2' : drv + '/Lights/tifs/',
                    'Elevation' : drv + '/GLOBE/',
                    'Elevation2' : drv + '/GLOBE/',
                    'RIX' : drv + '/RIX/',
                    'RIX2' : drv + '/RIX/',
                    'CRU' : drv + '/CRU/tif/',
                    'CropsYield' : drv + '/FAO/GAEZ/Agroclimatic/yield/rain_fed/',
                    'CSICrops' : drv + '/FAO/GAEZ/Agroclimatic/energyyield/rain_fed/',
                    'CSI' : drv + '/FAO/GAEZ/Agroclimatic/bestcropcal/CSI/',
                    'CSIWater' : drv + '/FAO/GAEZ/Agroclimatic/bestcropcal/CSI_water/',
                    'CSICycle' : drv + '/FAO/GAEZ/Agroclimatic/bestcropcal/latest/',
                    'CSICycleExtra' : drv + '/FAO/GAEZ/Agroclimatic/bestcropcal/latest/',
                    'CSICycleWater' : drv + '/FAO/GAEZ/Agroclimatic/bestcropcal/updated_water/',
                    'CSIPlow' : drv + '/FAO/GAEZ/Agroclimatic/bestcropcal/plow/',
                    'Malaria' : pathgis + 'Malaria/',
                    'HMI' : drv + '/HMI_Diff/',
                    'PopDensAfrica' : drv + '/UNEP/',
                    'Popdens' : drv + '/GPW/popdens/',
                    'Population' : drv + '/GPW/population/',
                    'HYDE' : drv + '/HYDE/hyde31_final/tif/',
                    'Tsetse' : pathgis + 'Tsetse/',
                    'Ecodiversity' : pathgis + 'Ecological_Zones/',
                    'EcodiversityLGM' : pathgis + 'PaleoVegetation/',
                    'Sea100' : pathcntrycyl,
                    'Coast' : pathcntrycyl,
                    'Inwater' : pathcntrycyl,
                    'PerInwater' : pathcntrycyl,
                    'FluctInwater' : pathcntrycyl,
                    'Landscan' : drv +  '/Landscan/',
                    'GEcon' : os.getenv('HOME')+'/LatexMe/CountryStability/data/GEconRaster/',
                    'HLD' : drv + '/Harmonized Light Data/light/',
}

# Paths to measures available for computations (if adding new measures one needs to change all places where the lists/dicts are created)
mainmeasures = ['Suitability','Lights','Lights2','Elevation','RIX','Elevation2','RIX2','CRU','CSI', 'CSICrops', 'CropsYield',
                'CSICycle', 'CSICycleExtra','CSIPlow','Malaria','HMI','PopDensAfrica','Popdens','Population',
                'HYDE','Tsetse','Ecodiversity','EcodiversityLGM','Sea100','Coast',
                'Inwater','PerInwater','FluctInwater','Landscan', 'GEcon', 'HLD']
mainmeasures.sort()

mainmeasuresw = ['Suitability','Lights','Lights2','Elevation','RIX','Elevation2','RIX2','CRU','CSI', 'CSICrops', 'CropsYield', 'CSIWater',
                'CSICycle', 'CSICycleExtra', 'CSICycleWater','CSIPlow','Malaria','HMI','PopDensAfrica','Popdens','Population',
                'HYDE','Tsetse','Ecodiversity','EcodiversityLGM','Sea100','Coast',
                'Inwater','PerInwater','FluctInwater','Landscan', 'GEcon', 'HLD']
mainmeasuresw.sort()

wgs84measures = ['Lights2','Elevation2','RIX2','PopDensAfrica','Popdens','Population','HYDE','Landscan','CSICrops','CropsYield', 'GEcon', 'HLD']
wgs84measures.sort()

shpmeasures = ['Ecodiversity','EcodiversityLGM','Sea100','Coast','Inwater','PerInwater','FluctInwater',]
shpmeasures.sort()

ecodivmeasures = ['Ecodiversity','EcodiversityLGM']
ecodivmeasures.sort()

ceameasures = list(set(mainmeasuresw).difference(set(wgs84measures)).difference(set(shpmeasures)))
ceameasures.sort()

noecodivmeasures = list(set(ceameasures).difference(set(ecodivmeasures)))
noecodivmeasures.sort()

namemeasures = {'Suitability' : -7,
                'Lights' : 7,
                'Lights2' : 7,
                'Elevation' : -4,
                'Elevation2' : -4,
                'RIX' : -4,
                'RIX2' : -4,
                'CRU' : -7,
                'CropsYield' : -4,
                'CSICrops' : -4,
                'CSI' : -7,
                'CSIWater' : -7,
                'CSICycle' : -7,
                'CSICycleExtra' : -7,
                'CSICycleWater' : -7,
                'CSIPlow' : -7,
                'Malaria' : -7,
                'HMI' : -4,
                'PopDensAfrica' : -4,
                'Popdens' : -4,
                'Population' : -4,
                'HYDE' : -4,
                'Tsetse' : -7,
                'Ecodiversity' : 0,
                'EcodiversityLGM' : 0,
                'Sea100' : 0,
                'Coast' : 0,
                'Inwater' : 0,
                'PerInwater' : 0,
                'FluctInwater' : 0,
                'Landscan' : -4,
                'GEcon' : -4,
                'HLD' : 0,
                }

# Functions to perform various operations
# Create vector of polygons and iso-codes
def GISData(myshp=pathcntrycyl + 'DCW_countriescyl.shp', adds=False):
    '''
    This function imports the shape file and converts everyhting to Shapely, and WKT and saves it in a
    pandas dataframe
    Additionally, it computes areas and boundary lengths if adds=True
    Usage:
    GISData(shp=file,[adds=adds])
    '''
    # Shape file
    if isinstance(myshp,str):
        data = gp.GeoDataFrame.from_file(myshp)
    elif isinstance(myshp,gp.GeoDataFrame):
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
        measures = list(mainmeasures)
    else:
        measures = list(measures)
    pass

    for measure in measures:
        print('Computing Stats for:',measure)
        rpath = pathmeasures[measure]
        if measure in ceameasures:
            dfin = df
            dfin.crs = df.crs
            mytiffiles=[mytiffile for mytiffile in os.listdir(rpath) if mytiffile.endswith('cyl.tif')]
            if measure == 'HMI':
                mytiffiles=[mytiffile for mytiffile in os.listdir(rpath) if mytiffile.endswith('.tif')]
            elif measure == 'CSICycle':
                mytiffiles=[mytiffile for mytiffile in os.listdir(rpath) if mytiffile.endswith('cyl.tif') and (mytiffile.find('1500')!=-1 or mytiffile.find('dif')!=-1)]
            elif measure == 'CSICycleExtra':
                mytiffiles=[mytiffile for mytiffile in os.listdir(rpath) if mytiffile.endswith('cyl.tif') and (mytiffile.find('1500')==-1 and mytiffile.find('dif')==-1)]
        elif measure in wgs84measures:
            dfin = dfnocyl
            dfin.crs = dfnocyl.crs
            mytiffiles=[mytiffile for mytiffile in os.listdir(rpath) if mytiffile.endswith('.tif') and mytiffile.find('cyl.tif')==-1]
            if measure=='HLD':
                mytiffiles=[mytiffile for mytiffile in mytiffiles if mytiffile.find('Harmonized')!=-1]
        else:
            dfin = df
            dfin.crs = df.crs
        if measure not in shpmeasures:
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
        #self.df.geometry = self.df.geometry.apply(lambda x: x.buffer(0))
        #self.dfnocyl.geometry = self.dfnocyl.geometry.apply(lambda x: x.buffer(0))
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
            self.measures = list(mainmeasures)
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
            if measure in ceameasures:
                dfin = self.df
                dfin.crs = self.df.crs
                mytiffiles=[mytiffile for mytiffile in os.listdir(rpath) if mytiffile.endswith('cyl.tif')]
                if measure == 'HMI':
                    mytiffiles=[mytiffile for mytiffile in os.listdir(rpath) if mytiffile.endswith('.tif')]
                elif measure == 'CSICycle':
                    mytiffiles=[mytiffile for mytiffile in os.listdir(rpath) if mytiffile.endswith('cyl.tif') and (mytiffile.find('1500')!=-1 or mytiffile.find('dif')!=-1)]
                elif measure == 'CSICycleExtra':
                    mytiffiles=[mytiffile for mytiffile in os.listdir(rpath) if mytiffile.endswith('cyl.tif') and (mytiffile.find('1500')==-1 and mytiffile.find('dif')==-1)]
            elif measure in wgs84measures:
                dfin = self.dfnocyl
                dfin.crs = self.dfnocyl.crs
                mytiffiles=[mytiffile for mytiffile in os.listdir(rpath) if mytiffile.endswith('.tif') and mytiffile.find('cyl.tif')==-1]
                if measure=='HLD':
                    mytiffiles=[mytiffile for mytiffile in mytiffiles if mytiffile.find('Harmonized')!=-1]
            else:
                dfin = self.df
                dfin.crs = self.df.crs
            if measure not in shpmeasures:
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
