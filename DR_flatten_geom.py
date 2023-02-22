#coding: utf-8

import geopandas as gpd
import pandas as pd 
from geopandas.tools import sjoin
from shapely.geometry import LineString
from shapely.geometry import Point, Polygon
from shapely.geometry import shape
from shapely.geometry.multipolygon import MultiPolygon
from descartes import PolygonPatch
import time
import math
import scipy.stats as stats
import numpy as np
import os, sys
from pyproj import CRS, Transformer
import fiona
import statsmodels.api as sm

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib as mpl
from osgeo import ogr, gdal,osr
from math import floor

def convert_3D_2D(geometry):
    '''
    Takes a GeoSeries of 3D Multi/Polygons (has_z) 
    # and returns a list of 2D Multi/Polygons
    '''
    new_geo = []
    for p in geometry:
        if p.has_z:
            if p.geom_type == 'Polygon':
                lines = [xy[:2] for xy in list(p.exterior.coords)]
                new_p = Polygon(lines)
                new_geo.append(new_p)
            elif p.geom_type == 'MultiPolygon':
                new_multi_p = []
                for ap in p:
                    lines = [xy[:2] for xy in list(ap.exterior.coords)]
                    new_p = Polygon(lines)
                    new_multi_p.append(new_p)
                new_geo.append(MultiPolygon(new_multi_p))
    return new_geo

if __name__ == "__main__":

    # Initiate a plot. 

    fig, ax = plt.subplots(figsize=(15, 15))

    # Read in shapefile with Z geometry values. 
    na_map = gpd.read_file('fire_severity_selected_v003.shp')
    
    crs = {'init': 'esri:102001'}
    
    na_map.plot(ax=ax, facecolor="none", edgecolor='k',linewidth=1, zorder=14, alpha=1)
    
    # Makes sure map looks correct. 
    plt.show()

    # Get rid of the Z values. 
    ng = convert_3D_2D(na_map['geometry'])

    # Make a new geometry column for the GeoDataFrame. 
    na_map['geometry'] = ng

    # Export to new ESRI shapefile. 
    na_map.to_file('DR_fires_flatten.shp', driver = 'ESRI Shapefile')
