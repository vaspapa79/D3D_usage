
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 09:17:33 2023

@author: vaspapa
"""

import netCDF4 as nc
import numpy as np
import os
import csv
from datetime import datetime

def haversine(lat1, lon1, lat2, lon2):
    """
    Calculate the Haversine distance between two points.
    """
    R = 6371  # Radius of the Earth in kilometers
    dlat = np.radians(lat2 - lat1)
    dlon = np.radians(lon2 - lon1)
    a = np.sin(dlat/2) * np.sin(dlat/2) + np.cos(np.radians(lat1)) * np.cos(np.radians(lat2)) * np.sin(dlon/2) * np.sin(dlon/2)
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))
    distance = R * c
    return distance

def circular_mean(angles, weights=None):
    """
    Calculate the circular mean of a set of angles.

    Parameters:
    - angles: numpy array of angles in degrees.
    - weights: optional weights for each angle.

    Returns:
    - Circular mean in degrees.
    """
    if weights is None:
        weights = np.ones_like(angles)

    x = np.sum(weights * np.cos(np.radians(angles)))
    y = np.sum(weights * np.sin(np.radians(angles)))

    mean_angle = np.degrees(np.arctan2(y, x))
    mean_angle = (mean_angle + 360) % 360  # Ensure the result is in [0, 360) range

    return mean_angle

def mean_t_2m_in_cycle(dataset, center_lat, center_lon, radius_km):
    """
    Calculate the mean value of 'T_2m' variable inside a specific cycle around a given center point.
    """
    lat = dataset.variables['lat'][:]
    lon = dataset.variables['lon'][:]
    t_2m = dataset.variables['T_2m'][:]

    # Find indices of the grid points within the specified radius around the center point
    indices_within_radius = np.where(haversine(center_lat, center_lon, lat, lon) <= radius_km)

    # Extract 'T_2m' values within the specified cycle
    t_2m_in_cycle = t_2m[:, indices_within_radius[0], indices_within_radius[1]]

    # Calculate the mean value
    mean_t_2m = np.mean(t_2m_in_cycle)

    return mean_t_2m

#def mean_cldfra_in_cycle(dataset, center_lat, center_lon, radius_km):
#    """
#    Calculate the mean value of 'T_2m' variable inside a specific cycle around a given center point.
#    """
#    lat = dataset.variables['lat'][:]
#    lon = dataset.variables['lon'][:]
#    cldfra = dataset.variables['CLDFRA'][:]
#
#    # Find indices of the grid points within the specified radius around the center point
#    indices_within_radius = np.where(haversine(center_lat, center_lon, lat, lon) <= radius_km)
#
#    # Extract 'T_2m' values within the specified cycle
#    cldfra_in_cycle = cldfra[:, indices_within_radius[0], indices_within_radius[1]]
#
#    # Calculate the mean value
#    mean_cldfra = np.mean(cldfra_in_cycle)
#
#    return mean_cldfra)
#

def mean_sw_d_in_cycle(dataset, center_lat, center_lon, radius_km):
    """
    Calculate the mean value of 'T_2m' variable inside a specific cycle around a given center point.
    """
    lat = dataset.variables['lat'][:]
    lon = dataset.variables['lon'][:]
    sw_d = dataset.variables['SW_d'][:]

    # Find indices of the grid points within the specified radius around the center point
    indices_within_radius = np.where(haversine(center_lat, center_lon, lat, lon) <= radius_km)

    # Extract 'T_2m' values within the specified cycle
    sw_d_in_cycle = sw_d[:, indices_within_radius[0], indices_within_radius[1]]

    # Calculate the mean value
    mean_sw_d = np.mean(sw_d_in_cycle)

    return mean_sw_d

def mean_rh_in_cycle(dataset, center_lat, center_lon, radius_km):
    """
    Calculate the mean value of 'RH' variable inside a specific cycle around a given center point.
    """
    lat = dataset.variables['lat'][:]
    lon = dataset.variables['lon'][:]
    rh = dataset.variables['rh_2m'][:]

    # Find indices of the grid points within the specified radius around the center point
    indices_within_radius = np.where(haversine(center_lat, center_lon, lat, lon) <= radius_km)

    # Extract 'RH' values within the specified cycle
    rh_in_cycle = rh[:, indices_within_radius[0], indices_within_radius[1]]

    # Calculate the mean value
    mean_rh = np.mean(rh_in_cycle)

    return mean_rh

def mean_wd_in_cycle(dataset, center_lat, center_lon, radius_km):
    """
    Calculate the mean value of 'WD' variable inside a specific cycle around a given center point.
    """
    lat = dataset.variables['lat'][:]
    lon = dataset.variables['lon'][:]
    wd = dataset.variables['wd_10m'][:]

    # Find indices of the grid points within the specified radius around the center point
    indices_within_radius = np.where(haversine(center_lat, center_lon, lat, lon) <= radius_km)

    # Extract 'WD' values within the specified cycle
    wd_in_cycle = wd[:, indices_within_radius[0], indices_within_radius[1]]

    # Calculate the mean value using circular_mean function
    mean_wd = circular_mean(wd_in_cycle)

    return mean_wd

def mean_ws_in_cycle(dataset, center_lat, center_lon, radius_km):
    """
    Calculate the mean value of 'WS' variable inside a specific cycle around a given center point.
    """
    lat = dataset.variables['lat'][:]
    lon = dataset.variables['lon'][:]
    ws = dataset.variables['ws_10m'][:]

    # Find indices of the grid points within the specified radius around the center point
    indices_within_radius = np.where(haversine(center_lat, center_lon, lat, lon) <= radius_km)

    # Extract 'WS' values within the specified cycle
    ws_in_cycle = ws[:, indices_within_radius[0], indices_within_radius[1]]

    # Calculate the mean value
    mean_ws = np.mean(ws_in_cycle)

    return mean_ws

def values_in_cycle(dataset, variable_name, center_lat, center_lon, radius_km):
    """
    Extract values of a variable within a specific cycle around a given center point.
    """
    lat = dataset.variables['lat'][:]
    lon = dataset.variables['lon'][:]
    variable = dataset.variables[variable_name][:]

    # Find indices of the grid points within the specified radius around the center point
    indices_within_radius = np.where(haversine(center_lat, center_lon, lat, lon) <= radius_km)

    # Extract variable values within the specified cycle
    values_in_cycle = variable[:, indices_within_radius[0], indices_within_radius[1]]

    return values_in_cycle
# Directory containing the NetCDF files
directory = 'C:\\Python_scripts\\WRF\\original\\'

# List all files in the directory
files = [f for f in os.listdir(directory) if f.startswith('wrfout_d04_') and f.endswith('_sn.nc')]

# Specify the center point and radius of the cycle
center_lat = 40.183720 # kesaria: 40.183720, velventos: 40.255576
center_lon = 21.859620 # kesaria: 21.859620, velventos: 22.056100
radius_km = 0.2 # grid resolution: 200m  

# Loop through each file
for file_name in files:
    file_path = os.path.join(directory, file_name)

    # Load the NetCDF file
    dataset = nc.Dataset(file_path)

    # Extract values within the specified cycle
    t_2m_values_in_cycle = values_in_cycle(dataset, 'T_2m', center_lat, center_lon, radius_km)
    rh_values_in_cycle = values_in_cycle(dataset, 'rh_2m', center_lat, center_lon, radius_km)
    wd_values_in_cycle = values_in_cycle(dataset, 'wd_10m', center_lat, center_lon, radius_km)
    ws_values_in_cycle = values_in_cycle(dataset, 'ws_10m', center_lat, center_lon, radius_km)
    sw_d_values_in_cycle = values_in_cycle(dataset, 'SW_d', center_lat, center_lon, radius_km)
#    cldfra_values_in_cycle = values_in_cycle(dataset, 'CLDFRA', center_lat, center_lon, radius_km)
    
    # Calculate the mean values in the specified cycle
    mean_t_2m = mean_t_2m_in_cycle(dataset, center_lat, center_lon, radius_km)
    mean_rh = mean_rh_in_cycle(dataset, center_lat, center_lon, radius_km)
    mean_wd_10m = mean_wd_in_cycle(dataset, center_lat, center_lon, radius_km)
    mean_ws_10m = mean_ws_in_cycle(dataset, center_lat, center_lon, radius_km)
    mean_sw_d = mean_sw_d_in_cycle(dataset, center_lat, center_lon, radius_km)
#    mean_cldfra = mean_cldfra_in_cycle(dataset, center_lat, center_lon, radius_km)
    
    print(f"\nFile: {file_name}")
    print("Mean T_2m in the specified cycle:", mean_t_2m)
    print("Mean RH_2m in the specified cycle:", mean_rh)
    print("Mean wd_10m in the specified cycle:", mean_wd_10m)
    print("Mean ws_10m in the specified cycle:", mean_ws_10m)
    print("Mean SW_d in the specified cycle:", mean_sw_d)
#    print("Mean CLDFRA in the specified cycle:", mean_cldfra)
  
    # Print or analyze the values
    print(f"\nFile: {file_name}")
    
    print("T_2m values in the specified cycle:")
    print(t_2m_values_in_cycle)

    print("\nRH_2m values in the specified cycle:")
    print(rh_values_in_cycle)

    print("\nWind Direction values in the specified cycle:")
    print(wd_values_in_cycle)

    print("\nWind Speed values in the specified cycle:")
    print(ws_values_in_cycle)
   
    print("\nSW_d values in the specified cycle:")
    print(sw_d_values_in_cycle)

#    print("\nCLDFRA values in the specified cycle:")
#    print(cldfra_values_in_cycle)
    
# Close the NetCDF dataset
    dataset.close()


