# -*- coding: utf-8 -*-
"""

This is a Python code that defines a function named get_netCDF_vals. The function takes in a
shapefile of a fire perimeter with associated information, path to McElhinny NetCDF files,
FWI code, resolution, and a factor multiplier as input parameters. It returns the maximum
value in the fire, which is either the closest point to the convex hull of the fire or a sum
of the points inside the fire. The function creates a regular grid of points inside the
bounding box of the input shapefile and extracts weather data from NetCDF files for each point
inside the fire. It then calculates the mean, median, and maximum values of the weather data
for the points inside the fire. The function handles cases where the shapefile has multiple
polygons, not continuous, that make up the shape.


"""

#Here we are importing the packages we need. 
import geopandas as gpd
import pandas as pd 
from geopandas.tools import sjoin
from shapely.geometry import LineString
from shapely.geometry import Point
from shapely.geometry import Polygon
from shapely.geometry import shape
from descartes import PolygonPatch
import time
import math
import scipy.stats as stats
import numpy as np
import os, sys
from pyproj import CRS, Transformer
import fiona

import statsmodels.api as sm
import statsmodels.formula.api as smf

import matplotlib.pyplot as plt
import matplotlib as mpl
from math import floor

from shapely.ops import unary_union

import warnings
warnings.filterwarnings('ignore')

from osgeo import ogr, gdal,osr

from datetime import datetime
from datetime import timedelta
from pytz import timezone
import xarray as xr

def get_netCDF_vals(fire_shapefile, path, code_type, res,factor,print_map=False,first_21_days=False):
    '''This is a function to get the value inside the fire.
    We will use to calculate the mean, median, max value for a fire.
    
    Parameters
    ----------
        fire_shapefile : GeoDataFrame
            fire perimeter + associated information
        path : string
            path to McElhinny NetCDF files on drive
        code_type : string 
            FWI code
        res : float
            specified resolution for regular grid of points inside fire
        factor : int
            multiplier for number of points considered
        print_map : bool
            if True, print a map of the points inside the fire and their values
            by default False
        first_21_days : bool
            if True, take average / mean / median of first 21 days of the fire
            if False, just take the metrics on the report date
            by default False
            
    Returns
    ----------
        float
            - mean, median, maximum, p90_mean, p90_med, p90_max  
    '''


    bounds = gpd.GeoDataFrame(geometry=[fire_shapefile['geometry']]).bounds 
    repdat = pd.to_datetime(fire_shapefile['REP_DATE'])
    
    if first_21_days: 
        edate = pd.to_datetime(fire_shapefile['REP_DATE']) + pd.DateOffset(days=21)
        all_dat = pd.date_range(repdat,edate-timedelta(days=1),freq='d')
    else: 
        all_dat = [repdat]
    xmax = np.nanmax(bounds['maxx'])
    xmin = np.nanmin(bounds['minx'])
    ymax = np.nanmax(bounds['maxy'])
    ymin = np.nanmin(bounds['miny'])

    # Calculate the number of rows cols to fill the bounding box at that resolution
    num_col = int((xmax - xmin) / res)
    num_row = int((ymax - ymin) / res)

    # Add the bounding box coords to the dataset so we can extrapolate the interpolation to cover whole area
    yProj_extent = [bounds['maxy'], bounds['miny']]
    xProj_extent = [bounds['maxx'], bounds['minx']]

    # Get the value for lat lon in each cell we just made
    Yi = np.linspace(np.min(yProj_extent), np.max(yProj_extent), num_row*factor)
    Xi = np.linspace(np.min(xProj_extent), np.max(xProj_extent), num_col*factor)

    # This line creates a meshgrid of x and y coordinates.
    Xi, Yi = np.meshgrid(Xi, Yi)
    # These lines flatten the meshgrid and convert it to a list.
    concat = np.array((Xi.flatten(), Yi.flatten())).T
    send_to_list = concat.tolist()

    meshPoints = [Point(item) for item in send_to_list]
    gdf = gpd.GeoDataFrame(geometry=meshPoints)

    DF = fire_shapefile
    try: #If there is a single polygon in the shapefile
        DF = unary_union(Polygon(DF['geometry'])) #Multipolygon --> Polygon
        poly_define = gpd.GeoDataFrame(geometry=[DF])
        # Get points falling in fire 
        within_fire = gdf[gdf.geometry.within(poly_define['geometry'][0])]
    # Catch the case where multiple polygons not continuous make up the shp
    except (NotImplementedError,AttributeError,TypeError) as e:
        DF = [unary_union(Polygon(geom)) for geom in list(DF['geometry'])]
        poly_define = gpd.GeoDataFrame(geometry=DF)
        # Left spatial join 
        within_fire = sjoin(gdf, poly_define, how='left',op='within')
        # Drop points that are not in the fire 
        within_fire = within_fire[~np.isnan(within_fire['index_right'])]


    inside_fire = []
    lon = []
    lat = [] 
    listP = within_fire


    
    for idx,p in listP.iterrows():
        
        mx,my=np.array(p['geometry'].coords.xy[0])[0], np.array(p['geometry'].coords.xy[1])[0]

        weather_val = extract(path,code_type,my,mx,all_dat)
        inside_fire.append(weather_val)
        lon.append(mx)
        lat.append(my)
        
    
    within_fire['max_val'] = inside_fire
    within_fire['lon'] = lon
    within_fire['lat'] = lat



    if len(within_fire) > 0:
        mean = np.nanmean([x[0] for x in inside_fire])  # get the mean val inside the fire for median in 21 days (or rep date)
        median = np.nanmedian([x[0] for x in inside_fire])
        maximum = np.nanmax([x[0] for x in inside_fire])

        p90_mean = np.nanmean([x[1] for x in inside_fire])  # get the mean val inside the fire for median in 21 days (or rep date)
        p90_med = np.nanmedian([x[1] for x in inside_fire])
        p90_max = np.nanmax([x[1] for x in inside_fire])


    else: 
        mean = np.nan
        median = np.nan
        maximum = np.nan
        p90_mean = np.nan
        p90_med = np.nan
        p90_max = np.nan

    if print_map: 

        fig, ax = plt.subplots(figsize=(15, 15))
        fire = gpd.GeoDataFrame(geometry=[fire_shapefile['geometry']])
        vals = [x[0] for x in inside_fire]
        sc = plt.scatter(within_fire['lon'],within_fire['lat'],
                             c=vals,cmap='Spectral_r',edgecolors='None',
                             vmin=min(vals),vmax=max(vals),s=15)

        
        fire.plot(ax=ax, facecolor="none", edgecolor='k',linewidth=1, zorder=14, alpha=1)
        plt.show()


    return

if __name__ == "__main__":

    dirname = ''

    var_name = 'FFMC' # Change variable name of interest here
    var_name_extended = 'fine_fuel_moisture_code'

    path_ffmc = 'ffmc/' # Path to stored files

    years = list(range(1985,2015+1))

    for y in years: 

        year = y 

        shp = gpd.read_file('doriana_fires_4326.shp')
        shp_type = shp[shp['SIZE_HA'] >= 10]
        shp_type = shp_type[shp_type['YEAR'] == year]

        ffmc_list = []

        for index, fire in shp_type.iterrows():
            year = int(fire['YEAR'])
            ave,med,max_val,ave90,med90,max_val90 = get_netCDF_vals(fire,path_ffmc+'/'+var_name_extended+'_'+str(year)+'.nc',var_name,0.001,1) #0.001

            ffmc_list.append([ave,med,max_val,ave90,med90,max_val90])

        shp_type[var_name+'_mean'] = [x[0] for x in ffmc_list]
        shp_type[var_name+'_median'] = [x[1] for x in ffmc_list]
        shp_type[var_name+'_max'] = [x[2] for x in ffmc_list]

        shp_type[var_name+'_p90mean'] = [x[3] for x in ffmc_list]
        shp_type[var_name+'_p90median'] = [x[4] for x in ffmc_list]
        shp_type[var_name+'_p90max'] = [x[5] for x in ffmc_list]

        shp_type.to_csv(var_name+'_in_cnfdb_fires_'+str(year)+'.txt',sep=',')
