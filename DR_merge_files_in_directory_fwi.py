# -*- coding: utf-8 -*-
"""

Code to merge all FWI files for individual years in specified directory into one csv file. 
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

if __name__ == "__main__":

    dirname = ''

    # FWI value of interest

    var_name = 'DMC'

    tracker = []

    # Get list of years data was originally generated for 

    years = list(range(1985,2015+1))

    # For each file in specified directory (syear) 

    for f in os.listdir('syear/'):

        # Make sure file is a text file, as there are zip files in the directory as well
        # Also make sure that the FWI metric name is in the file name 

        if var_name in f and f.endswith('.txt'):
            
            # Read the text file using Pandas 

            df = pd.read_csv('syear/'+f)

            # Append the dataframe to the 'tracker' list 

            tracker.append(df[['FIRE_ID',var_name+'_mean',var_name+'_median',\
                              var_name+'_max',var_name+'_p90mean']])

    # Merge the dataframes in the 'tracker' list and reset the index

    merger = pd.concat(tracker).reset_index()

    # Rename the file according to specifications of project
    merger[var_name+'_p90'] = merger[var_name+'_p90mean']

    # Export the merged file to a comma-delimited csv 
    merger[['FIRE_ID',\
            var_name+'_mean',var_name+'_median',var_name+'_max',var_name+'_p90']].to_csv(var_name+'_in_cnfdb_fires_merged_6col.csv',sep=',')
