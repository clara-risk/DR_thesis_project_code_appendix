# Master's Thesis Project Data Preparation Appendix - Doriana Romualdi 

The following repository contains the Python code that was used to prepare the data for a Master's thesis project as well as the data in final form.  

# 1: DR_canlabs_extract_preNBR_postNBR_RBR.py 

This script was used to extract the preNBR, postNBR, and calculate RBR from the preNBR and dNBR values from Canlabs, found here: https://ftp.maps.canada.ca/pub/nrcan_rncan/Forest-fires_Incendie-de-foret/CanLaBS-Burned_Severity-Severite_des_feux/ 

The final product of this script is a csv file containing the maximum, median, mean, and 90th percentile values inside each selected fire perimeter from the Canadian National Fire Database (CNFDB) in British Columbia. 

# 2: DR_extract_netCDF_fwi_data.py 

This script was used to extract the FWI netCDF data originating from this publication: https://essd.copernicus.org/articles/12/1823/2020/essd-12-1823-2020-assets.html and retrieved from records stored on Zenodo, for example: https://zenodo.org/record/3540924#.Y_gu3HbMJPY. 

The final product of this script is a csv file containing the maximum, median, mean, and 90th percentile values for the FWI values inside each selected fire perimeter from the Canadian National Fire Database (CNFDB) in British Columbia on the fire report date. The script contains an option to calculate these metrics for the first 21 days of the fire if necessary as well. 

The repository also contains basic Python scripts used to perform intermediate data processing tasks: 

# 1: DR_flatten_geom.py 

This script was used to convert the CNFDB shapefile geometry data into 2D data (only containing latitude and longitude) from 3D data (containing latitude, longitude, and elevation data, i.e. "Polygon Z" or "MultiPolygon Z" data). This was done to standardize the data and ensure that it could be used to effectively extract raster data pixels within a fire perimeter. 

# 2: DR_merge_files_in_directory_fwi.py 

This script was used to merge all text files for a specified FWI variable (DC, DMC, etc.) containing data for different years (1985-2015) into a single dataframe (exported as a csv file). 
