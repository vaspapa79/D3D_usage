# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 10:41:35 2023

@author: vaspapa
"""

import xarray as xr
import os

# Directory containing the original NetCDF files
input_directory = r'C:\Python_scripts\WRF\original'

# Directory to save the modified NetCDF files
output_directory = r'C:\Python_scripts\WRF\modified'

# List of files in the input directory
files = [f for f in os.listdir(input_directory) if f.endswith('.nc')]

# List of variables you want to keep
variables_to_keep = ['LandUse','T_2m', 'p_sfc', 'u_10m_gr','wd_10m','ws_10m', 'precip_c', 'v_10m_gr', 'rh_2m', 'SW_d']

# Loop through each file
for file in files:
    # Create the full path to the input file
    input_file_path = os.path.join(input_directory, file)

    # Open the NetCDF file
    ds = xr.open_dataset(input_file_path)

    # Select only the desired variables
    ds_subset = ds[variables_to_keep]

    # Create a new file name for the subset in the output directory
    output_file_path = os.path.join(output_directory, 'sub_' + file)

    # Save the subset to a new NetCDF file
    ds_subset.to_netcdf(output_file_path)

    # Close the current dataset
    ds.close()

from pyproj import CRS, Transformer

# Directory containing the original NetCDF files
original_directory = r'C:\Python_scripts\WRF\modified'

# Directory to save the Cartesian NetCDF files
cartesian_directory = r'C:\Python_scripts\WRF\cartesian'

# List of files in the original directory
files = [f for f in os.listdir(original_directory) if f.endswith('.nc')]

# Define the projection for conversion (WGS 84)
original_crs = CRS.from_string('+proj=latlong +datum=WGS84')

# Define the projection for the Cartesian coordinates (UTM Zone 34N)
cartesian_crs = CRS.from_epsg(32634)

# Create a transformer for the coordinate conversion
transformer = Transformer.from_crs(original_crs, cartesian_crs, always_xy=True)

# Loop through each file
for file in files:
    # Create the full path to the original file
    original_file_path = os.path.join(original_directory, file)

    # Open the original NetCDF file
    ds = xr.open_dataset(original_file_path)

    # Convert 'lat' and 'lon' to Cartesian coordinates
    lon, lat = xr.broadcast(ds['lon'], ds['lat'])  # Ensure lon and lat have the same shape
    x, y = transformer.transform(lon.values.flatten(), lat.values.flatten())

    # Update 'lat' and 'lon' variables with Cartesian coordinates
    ds['lat'].values = y.reshape(ds['lat'].shape)
    ds['lon'].values = x.reshape(ds['lon'].shape)

    # Create a new file name for the Cartesian dataset
    cartesian_file_path = os.path.join(cartesian_directory, 'cs_' + file)

    # Save the Cartesian dataset to a new NetCDF file in the Cartesian directory
    ds.to_netcdf(cartesian_file_path)

    # Close the current dataset
    ds.close()

# Directory containing the Cartesian NetCDF files
cartesian_directory = r'C:\Python_scripts\WRF\cartesian'

# List of files in the Cartesian directory
files = [f for f in os.listdir(cartesian_directory) if f.endswith('.nc')]

# Loop through each file
for file in files:
    # Create the full path to the Cartesian file
    cartesian_file_path = os.path.join(cartesian_directory, file)

    # Open the Cartesian NetCDF file
    ds = xr.open_dataset(cartesian_file_path)

    # Print information about the file
    print(f"File: {file}")
    print("Variables:")
    print(ds.data_vars)
    print("Coordinates:")
    print(ds.coords)
    print("Attributes:")
    print(ds.attrs)
    print("\n")

    # Plot the 'lat' and 'lon' variables if you want to visualize them
    #ds['T_2m'].plot()

    # Close the current dataset
    ds.close()
