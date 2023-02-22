# -*- coding: utf-8 -*-
"""

Code to extract either (a) pre-NBR, (b) post-NBR, or (c) calculate RBR for
values within certain specified fire perimeters in BC for DR's thesis project.

"""


#Here we are importing the packages we need.
#The above code imports necessary libraries such as
#Geopandas, Shapely, Descartes, PyProj, Fiona, etc.
#These libraries are used for geographic data manipulation,
#statistical analysis, plotting, etc.

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

def get_val_in_fire_canlabs(fire_shapefile, shapefile,data_surface,\
                            transform,size,srcds,year_fire,fire_years_raster,\
                            data_pre,data_post,factor,threshold=0.8,show=False,return_var = 'RBR'):
    '''This is a function to get the value inside the fire.
    We will use to calculate the mean, median, max, & p90 value for a fire.
    
    Parameters
    ----------
        fire_shapefile : string
            path to the fire shapefile 
        shapefile : string
            path to the study area shapefile
        data_surface : ndarray
            an array of values in the study area (dNBR, Canlabs) 
        transform : list 
            list describing GeoTransform of raster 
        size : list 
            pixel dimensions
        srcds : GDAL object
            read in raster
        year_fire : int
            year of fire from CNFDB
        fire_years_raster : ndarray 
            fire years raster from canlabs in ndarray format
        data_pre : ndarray
            an array of values in the study area (pre-NBR, Canlabs)
        data_post : ndarray
            an array of values in the study area (post-NBR, Canlabs)
        factor : int
            factor by which to multiply number of points in each pixel 
            higher = more accurate 
        threshold : float
            how much cloud mask to tolerate      
            by default it is 0.8
        show : bool
            whether to show a map of the generated points inside the fire perimeter and their values
            by default it is False
        return_var : str
            which value to return for fire perimeter: pre-NBR, post-NBR, RBR
            by default it is RBR, others are 'pre' or 'post'
            
    Returns
    ----------
        float
            average value inside fire
            median value inside fire
            maximum value inside fire
            90th percentile value inside fire
            cause of no data inside fire 
    '''

    #Flip the arrays upside down so that they can match the orientation of the shapefile


    fire_years_raster = np.flipud(fire_years_raster)
    data_surface = np.flipud(data_surface)
    data_pre = np.flipud(data_pre)

    #Get the bounding box of the shapefile
    bounds = gpd.GeoDataFrame(geometry=[fire_shapefile['geometry']]).bounds 
    xmax = np.nanmax(bounds['maxx'])
    xmin = np.nanmin(bounds['minx'])
    ymax = np.nanmax(bounds['maxy'])
    ymin = np.nanmin(bounds['miny'])

    #Calculate the origin, max, width, and height for the area to be filled with points

    xOrigin = transform[0]
    yOrigin = transform[3]
    xMax = xOrigin + transform[1] * size[0]
    yMin = yOrigin + transform[5] * size[1]
    pixelWidth = transform[1]
    pixelHeight = -transform[5]

    # Calculate the number of rows cols to fill the bounding box at that resolution
    num_col = int((xmax - xmin) / pixelHeight)+1
    num_row = int((ymax - ymin) / pixelWidth)+1

    # Add the bounding box coords to the dataset so we can extrapolate the interpolation to cover whole area
    yProj_extent = [bounds['maxy'], bounds['miny']]
    xProj_extent = [bounds['maxx'], bounds['minx']]

    # Get the value for lat lon in each cell we just made
    # Factor refers to how many points per pixel 
    Yi = np.linspace(np.min(yProj_extent), np.max(yProj_extent), num_row*factor)
    Xi = np.linspace(np.min(xProj_extent), np.max(xProj_extent), num_col*factor)

    Xi, Yi = np.meshgrid(Xi, Yi)
    # Because we are not using the lookup file, send in X,Y order
    # Concatenate the generated points and send them to a list
    concat = np.array((Xi.flatten(), Yi.flatten())).T
    send_to_list = concat.tolist()

    meshPoints = [Point(item) for item in send_to_list]

    gdf = gpd.GeoDataFrame(geometry=meshPoints)

    DF = fire_shapefile

    # Separate the shapefile polygons and check if they are a single polygon or multiple polygons
    try: #If there is a single polygon in the shapefile

        DF = unary_union(DF['geometry'])
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


    #Get the data values at the points that are within the fire polygon


    inside_fire = []
    years_inside_fire = [] 
    inside_fire_pre = []
    inside_fire_post = [] 
    lon = []
    lat = [] 
    listP = within_fire

    for idx,p in listP.iterrows():
        mx,my=np.array(p['geometry'].coords.xy[0])[0], np.array(p['geometry'].coords.xy[1])[0]

        col = int(np.floor((mx - xOrigin)) / pixelHeight)
        row = int(np.floor((my - yOrigin)) / pixelWidth)
        
        sev = data_surface[row-1][col]
            
##        if sev != -32768.0: #No data value
        inside_fire.append(data_surface[row-1][col])
        inside_fire_pre.append(data_pre[row-1][col])
        inside_fire_post.append(data_post[row-1][col])
        years_inside_fire.append(int(fire_years_raster[row-1][col]))
            
##        else:
##            inside_fire.append(np.nan)
##            inside_fire_pre.append(np.nan)
##            years_inside_fire.append(np.nan)
        
        lon.append(mx)
        lat.append(my)

    # Populate the within_fire dataframe with values that have been extracted above
        
    
    within_fire['dNBR'] = inside_fire
    within_fire['pre'] = inside_fire_pre
    within_fire['post'] = inside_fire_post

    #Make sure the no data values from Canlabs are reassigned to nan values
    
    within_fire.loc[within_fire['pre'] == -32768.0,'pre'] = np.nan
    within_fire.loc[within_fire['dNBR'] == -32768.0,'dNBR'] = np.nan
    within_fire.loc[within_fire['post'] == -32768.0,'post'] = np.nan

    # Add lat lon year information 
    within_fire['lon'] = lon
    within_fire['lat'] = lat
    within_fire['years'] = years_inside_fire

    #Calculate RBR
    within_fire['process_step'] = (within_fire['pre']/1000)+1.001
    within_fire['RBR'] = within_fire['dNBR']/within_fire['process_step']


    if show: # If specified, then show a map of the generated points inside the fire perimeter

        if return_var == 'RBR':
            ras = data_surface
        elif return_var == 'pre':
            ras = data_pre
        else:
            ras = data_post 

        
        fig, ax = plt.subplots(figsize=(15, 15))
        gk = plt.imshow(data_surface,extent=(
                    xOrigin,xMax,yOrigin,
                    yMin),vmin=0,vmax=2015,cmap='Spectral')
            
        sc = plt.scatter(mx,my,
                             c=ras[row-1][col],cmap='Spectral',edgecolors='k',
                             vmin=0,vmax=2015,s=38)
        gpd.GeoDataFrame(geometry=
                             [fire_shapefile['geometry']]).plot(ax=ax,
                                                                facecolor='None',
                                                                edgecolor='k')
        c = fig.colorbar(gk)
        c.set_label(return_var, rotation=270)
        ax.set_title(fire_shapefile['FIRE_ID'])
        ax.set_ylim([my-50000, my+50000])
        ax.set_yticks([])
        ax.set_xlim([mx-50000,mx+50000])
        ax.set_xticks([])
        plt.show()

    # Proceed if the fire year from the CNFDB fire matches the information from Canlabs

    if year_fire in list(within_fire['years']):

        # Get rid of all pixels where fire year does not match
        # This is really important in the case where the area got burned twice 

        within_fire = within_fire[within_fire['years'] == year_fire]

        # Proceed if pixels are left within the fire 

        if len(within_fire) > 0:

            # Check the cloud mask 
            numnan = len(within_fire[~np.isnan(within_fire[return_var])]) / len(within_fire[return_var])


            # Ensure that enough data is present (clouds do not cover significant amount of fire) 

            if numnan >= threshold:
                mean = np.nanmean(list(within_fire[return_var]))  # get the mean val inside the fire
                median = np.nanmedian(list(within_fire[return_var]))
                maximum = np.nanmax(list(within_fire[return_var]))

                no_nan = within_fire[within_fire[return_var].notnull()]
                p90 = np.percentile(list(no_nan[return_var]), 90)
                nd_c = 'data_present'

            # If not enough, assigned variables to be returned to nan and save information about why 
            else:
                mean = np.nan
                median = np.nan
                maximum = np.nan
                p90 = np.nan
                nd_c = 'threshold_not_met'

        # This is the case where there are no Canlabs pixels within the fire... such as if clouds cover all of it. 
                
        else:
            mean = np.nan
            median = np.nan
            maximum = np.nan
            p90 = np.nan
            nd_c = 'no_pixel_in_fire'

    # The case where the fire year does not match or there is no year data available. 

    else: 

        mean = np.nan
        median = np.nan
        maximum = np.nan      
        p90 = np.nan
        nd_c = 'no_year_data'

    # Return the information from the Canlabs data inside the fire perimeter 

    return mean,median,maximum,p90,nd_c


if __name__ == "__main__":

    # The following code is reading in various geospatial datasets and performing data manipulation on them.

    # The first two lines set up some variables, dirname and rv.
    # The next several lines use the geopandas library to read in a shapefile of Canada in the Lambert Conformal
    # Conic projection (bc_canada_lambert.shp), and another shapefile of fire severity data (fire_severity_selected_v003.shp).
    # The fire_severity_selected_v003.shp file is filtered to only include records where the fire size is greater than
    # or equal to 10 hectares, the source agency is 'BC', the year is between 1985 and 2015, and the cause is 'L' (lightning).
    # The filtered data is then projected to the same Lambert Conformal Conic projection as the bc_canada_lambert.shp file.
    # The remaining lines read in various GeoTIFF rasters (clipped_year.tif, clipped_post_NBR.tif, clipped_pre_NBR.tif, clipped_dNBR.tif)
    # using the gdal library, extract data from the raster bands as arrays, and store them in variables (ydata, postdata, predata, ddata).

    
    dirname = ''

    rv = 'RBR'

    sa = gpd.read_file('bc_canada_lambert.shp')
    shp = gpd.read_file('fire_severity_selected_v003.shp')
    shp = shp[shp['SIZE_HA'] >= 10]
    shp_prov = shp[shp['SRC_AGENCY'].isin(['BC'])]
    shp_year = shp_prov[shp_prov['YEAR'] >= 1985]
    shp_year = shp_year[shp_year['YEAR'] <= 2015]
    shp_type = shp_year[shp_year['CAUSE'] == 'L'].to_crs('PROJCRS["NAD_1983_Canada_Lambert",BASEGEOGCRS["NAD83",DATUM["North American Datum 1983",ELLIPSOID["GRS 1980",6378137,298.257222101004,LENGTHUNIT["metre",1]]],PRIMEM["Greenwich",0,ANGLEUNIT["degree",0.0174532925199433]],ID["EPSG",4269]],CONVERSION["Lambert Conic Conformal (2SP)",METHOD["Lambert Conic Conformal (2SP)",ID["EPSG",9802]],PARAMETER["Latitude of false origin",0,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8821]],PARAMETER["Longitude of false origin",-95,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8822]],PARAMETER["Latitude of 1st standard parallel",49,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8823]],PARAMETER["Latitude of 2nd standard parallel",77,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8824]],PARAMETER["Easting at false origin",0,LENGTHUNIT["metre",1],ID["EPSG",8826]],PARAMETER["Northing at false origin",0,LENGTHUNIT["metre",1],ID["EPSG",8827]]],CS[Cartesian,2],AXIS["easting",east,ORDER[1],LENGTHUNIT["metre",1,ID["EPSG",9001]]],AXIS["northing",north,ORDER[2],LENGTHUNIT["metre",1,ID["EPSG",9001]]]]') #esri:102001


    fire_sev_list = []
    y_ds1 = gdal.Open('clipped_year.tif')
    yb2=y_ds1.GetRasterBand(1)
    ycols = y_ds1.RasterXSize
    yrows = y_ds1.RasterYSize
    ydata = yb2.ReadAsArray(0, 0, ycols, yrows)

    src_ds1 = gdal.Open('clipped_post_NBR.tif')
    rb2=src_ds1.GetRasterBand(1)
    transform=src_ds1.GetGeoTransform()
    cols = src_ds1.RasterXSize
    rows = src_ds1.RasterYSize
    postdata = rb2.ReadAsArray(0, 0, cols, rows)

    src_ds3 = gdal.Open('clipped_pre_NBR.tif')
    rb3=src_ds3.GetRasterBand(1)
    cols = src_ds1.RasterXSize #same
    rows = src_ds1.RasterYSize
    predata = rb3.ReadAsArray(0, 0, cols, rows)

    src_ds3 = gdal.Open('clipped_dNBR.tif')
    rb3=src_ds3.GetRasterBand(1)
    cols = src_ds1.RasterXSize #same
    rows = src_ds1.RasterYSize
    ddata = rb3.ReadAsArray(0, 0, cols, rows)

    # This following code performs the following tasks:

    # It loops over each row of the pandas DataFrame shp_type using the .iterrows() method and extracts the
    # FIRE_ID and YEAR values for each row.
    # For each row, it calls the get_val_in_fire_canlabs() function with a series of arguments including fire,
    # shp, ddata, transform, cols, rows, src_ds1, year, ydata, predata, postdata, 1, threshold=0.8, and return_var='RBR'.
    # This function likely calculates the fire severity metrics for the fire associated with the current row of shp_type.

    # The results of the get_val_in_fire_canlabs() function are appended to the fire_sev_list list.

    # The results in fire_sev_list are then added as new columns to the shp_type DataFrame, with column names based on the
    # value of rv. Specifically, the mean, median, maximum, and 90th percentile values are added as columns, along with a
    # new column indicating the cause of the fire (rv+'_nd_cause').

    # The resulting shp_type DataFrame is saved as a CSV file named 'fire_'+rv+'_in_fires.csv'.

    for index, fire in shp_type.iterrows():
        print(fire['FIRE_ID'])
        year = int(fire['YEAR'])

        ave,med,max_val,p90,nd = get_val_in_fire_canlabs(fire, shp,ddata,transform,\
                                                  (cols,rows,),src_ds1,year,\
                                                  ydata,predata,postdata,1,threshold=0.8,return_var='RBR') #0.8
        fire_sev_list.append([ave,med,max_val,p90,nd])

    shp_type['fire_'+rv+'_mean'] = [x[0] for x in fire_sev_list]
    shp_type['fire_'+rv+'_median'] = [x[1] for x in fire_sev_list]
    shp_type['fire_'+rv+'_max'] = [x[2] for x in fire_sev_list]
    shp_type['fire_'+rv+'_p90'] = [x[3] for x in fire_sev_list]
    shp_type[rv+'_nd_cause'] = [x[4] for x in fire_sev_list]

    shp_type.to_csv('fire_'+rv+'_in_fires.csv',sep=',')

    # This follwing code performs several operations on the GeoDataFrame object named shp_type:

    # It selects rows where the 'fire_' + rv + '_median' column does not contain a NaN value and assigns the result to shp_type.
    # It extracts a subset of columns from shp_type and saves the result to a CSV file named 'fire_' + rv + '_in_cnfdb_fires_nonnan.csv'.
    # It creates a new GeoDataFrame with the same rows and columns as shp_type and assigns it to shp_type.
    # It sets the coordinate reference system (CRS) of the new GeoDataFrame to 'PROJCRS["NAD_1983_Canada_Lambert", ... ]' and assigns the result to shp_type.
    # It saves the new GeoDataFrame to a shapefile named 'fire_' + rv + '_selected'.

    shp_type = shp_type[~np.isnan(shp_type['fire_'+rv+'_median'])]
    shp_type[['SRC_AGENCY','FIRE_ID','FIRENAME','YEAR','MONTH','DAY','REP_DATE',\
              'DATE_TYPE','OUT_DATE','DECADE','SIZE_HA','CALC_HA','CAUSE','ACQ_DATE',\
              'fire_'+rv+'_mean','fire_'+rv+'_max','fire_'+rv+'_median','fire_'+rv+'_p90']].to_csv('fire_'+rv+'_in_cnfdb_fires_nonnan.csv',sep=',')
    shp_type= gpd.GeoDataFrame(shp_type,crs='PROJCRS["NAD_1983_Canada_Lambert",BASEGEOGCRS["NAD83",DATUM["North American Datum 1983",ELLIPSOID\
["GRS 1980",6378137,298.257222101004,LENGTHUNIT["metre",1]]],PRIMEM["Greenwich",0,ANGLEUNIT["degree",0.0174532925199433]]\
,ID["EPSG",4269]],CONVERSION["Lambert Conic Conformal (2SP)",METHOD["Lambert Conic Conformal (2SP)",ID["EPSG",9802]],\
PARAMETER["Latitude of false origin",0,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8821]],PARAMETER["Longitude of false origin",\
-95,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8822]],PARAMETER["Latitude of 1st standard parallel",49,ANGLEUNIT\
["degree",0.0174532925199433],ID["EPSG",8823]],PARAMETER["Latitude of 2nd standard parallel",77,ANGLEUNIT\
["degree",0.0174532925199433],ID["EPSG",8824]],PARAMETER["Easting at false origin",0,LENGTHUNIT["metre",1],ID["EPSG",8826]],\
PARAMETER["Northing at false origin",0,LENGTHUNIT["metre",1],ID["EPSG",8827]]],CS[Cartesian,2],AXIS["easting",east,ORDER[1],\
LENGTHUNIT["metre",1,ID["EPSG",9001]]],AXIS["northing",north,ORDER[2],LENGTHUNIT["metre",1,ID["EPSG",9001]]]]',geometry=shp_type['geometry'])
    shp_type.to_file(driver = 'ESRI Shapefile',filename='fire_'+rv+'_selected')
