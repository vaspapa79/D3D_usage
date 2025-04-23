# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 13:11:32 2023

@author: vaspapa
"""

import xarray as xr
import numpy as np
import os

# Function to convert NetCDF files to Delft3D format
def convert_to_delft3d(nc_files, output_file):
    # Read the first NetCDF file to get dimensions and time information
    ds_first = xr.open_dataset(nc_files[0])
    sw_d_first = ds_first['SW_d'].values
    time, south_north, west_east = sw_d_first.shape

    # Open output file for writing
    with open(output_file, 'w') as f:
        # Write adjusted header
        f.write("FileVersion     = 1.03\n")
        f.write("filetype        = meteo_on_equidistant_grid\n")
        f.write("NODATA_value    = -999.000\n")
        f.write(f"n_cols          = {west_east}\n")
        f.write(f"n_rows          = {south_north}\n")
        f.write("grid_unit       = m\n")  # Adjusted grid unit to meters
        f.write("x_llcenter      = 560791.9\n")  # Adjusted x_llcenter value
        f.write("y_llcenter      = 4432953.0\n")  # Adjusted y_llcenter value
        f.write("dx              = 206.5\n")
        f.write("dy              = 206.5\n")
        f.write("n_quantity      = 1\n")
        f.write("quantity1       = sw_radiation_flux\n")
        f.write("unit1           = W/m2\n")

        # Initialize the counter to start at the desired reference time (example 2009-11-01 10:00:00 +00:00)
        reference_time = np.datetime64("2018-07-00T00:00:00Z")
        counter_hours = 0.0

        # Iterate through each NetCDF file
        for nc_file in nc_files:
            # Read NetCDF file
            ds = xr.open_dataset(nc_file)
            sw_d = ds['SW_d'].values

            # Get the time dimension
            times = ds.time.values

            # Iterate over each time step and write data
            for i, time_val in enumerate(times):
                # Calculate the time difference from the reference time in hours
                time_difference_hours = (time_val - reference_time) / np.timedelta64(1, 'h')

                # Format time as desired
                formatted_time = f"TIME = {time_difference_hours:.1f} hours since {str(reference_time)[:-6].replace('T', ' ')}:00:00 +00:00\n"
                f.write(formatted_time)

                # Update the counter by 1 hour for the next time step
                counter_hours += 1.0

                # Write data with three decimal places for each time step
                formatted_values = ['%07.2f' % value if value >= 0 else '%.2f' % value for value in sw_d[i].flatten()]

                # Write values, starting a new line every 200 values
                for j, val in enumerate(formatted_values, start=1):
                    f.write(val)
                    if j % 200 == 0 or j == len(formatted_values):
                        f.write('\n')
                    else:
                        f.write('\t')

            print(f"Data from {nc_file} written to {output_file}")

    print(f"Conversion completed. Delft3D file saved as: {output_file}")

# Directory containing NetCDF files
nc_dir = r'C:\Python_scripts\WRF\cartesian'

# List of NetCDF files
nc_files = [os.path.join(nc_dir, file) for file in os.listdir(nc_dir) if file.endswith('.nc')]

# Output Delft3D file with a unique name
output_file = os.path.join(nc_dir, 'sw_radiation_flux.ams')

# Convert multiple NetCDF files to a single Delft3D file
convert_to_delft3d(nc_files, output_file)


# Function to convert NetCDF files to Delft3D format
def convert_to_delft3d(nc_files, output_file):
    # Read the first NetCDF file to get dimensions and time information
    ds_first = xr.open_dataset(nc_files[0])
    y_wind_first = ds_first['v_10m_gr'].values
    time, south_north, west_east = y_wind_first.shape

    # Open output file for writing
    with open(output_file, 'w') as f:
        # Write adjusted header
        f.write("FileVersion     = 1.03\n")
        f.write("filetype        = meteo_on_equidistant_grid\n")
        f.write("NODATA_value    = -999.000\n")
        f.write(f"n_cols          = {west_east}\n")
        f.write(f"n_rows          = {south_north}\n")
        f.write("grid_unit       = m\n")  # Adjusted grid unit to meters
        f.write("x_llcenter      = 560791.9\n")  # Adjusted x_llcenter value
        f.write("y_llcenter      = 4432953.0\n")  # Adjusted y_llcenter value
        f.write("dx              = 206.5\n")
        f.write("dy              = 206.5\n")
        f.write("n_quantity      = 1\n")
        f.write("quantity1       = y_wind\n")
        f.write("unit1           = m s-1\n")

        # Initialize the counter to start at the desired reference time (2009-11-01 10:00:00 +00:00)
        reference_time = np.datetime64("2018-07-00T00:00:00Z")
        counter_hours = 0.0

        # Iterate through each NetCDF file
        for nc_file in nc_files:
            # Read NetCDF file
            ds = xr.open_dataset(nc_file)
            y_wind = ds['v_10m_gr'].values

            # Get the time dimension
            times = ds.time.values

            # Iterate over each time step and write data
            for i, time_val in enumerate(times):
                # Calculate the time difference from the reference time in hours
                time_difference_hours = (time_val - reference_time) / np.timedelta64(1, 'h')

                # Format time as desired
                formatted_time = f"TIME = {time_difference_hours:.1f} hours since {str(reference_time)[:-6].replace('T', ' ')}:00:00 +00:00\n"
                f.write(formatted_time)

                # Update the counter by 1 hour for the next time step
                counter_hours += 1.0

                # Write data with three decimal places for each time step
                formatted_values = ['%.3f' % value if value >= 0 else '%.2f' % value for value in y_wind[i].flatten()]

                # Write values, starting a new line every 200 values
                for j, val in enumerate(formatted_values, start=1):
                    f.write(val)
                    if j % 200 == 0 or j == len(formatted_values):
                        f.write('\n')
                    else:
                        f.write('\t')

            print(f"Data from {nc_file} written to {output_file}")

    print(f"Conversion completed. Delft3D file saved as: {output_file}")

# Directory containing NetCDF files
nc_dir = r'C:\Python_scripts\WRF\cartesian'

# List of NetCDF files
nc_files = [os.path.join(nc_dir, file) for file in os.listdir(nc_dir) if file.endswith('.nc')]

# Output Delft3D file with a unique name
output_file = os.path.join(nc_dir, 'y_wind.amv')

# Convert multiple NetCDF files to a single Delft3D file
convert_to_delft3d(nc_files, output_file)

# Function to convert NetCDF files to Delft3D format
def convert_to_delft3d(nc_files, output_file):
    # Read the first NetCDF file to get dimensions and time information
    ds_first = xr.open_dataset(nc_files[0])
    x_wind_first = ds_first['u_10m_gr'].values
    time, south_north, west_east = x_wind_first.shape

    # Open output file for writing
    with open(output_file, 'w') as f:
        # Write adjusted header
        f.write("FileVersion     = 1.03\n")
        f.write("filetype        = meteo_on_equidistant_grid\n")
        f.write("NODATA_value    = -999.000\n")
        f.write(f"n_cols          = {west_east}\n")
        f.write(f"n_rows          = {south_north}\n")
        f.write("grid_unit       = m\n")  # Adjusted grid unit to meters
        f.write("x_llcenter      = 560791.9\n")  # Adjusted x_llcenter value
        f.write("y_llcenter      = 4432953.0\n")  # Adjusted y_llcenter value
        f.write("dx              = 206.5\n")
        f.write("dy              = 206.5\n")
        f.write("n_quantity      = 1\n")
        f.write("quantity1       = x_wind\n")
        f.write("unit1           = m s-1\n")

        # Initialize the counter to start at the desired reference time (2009-11-01 10:00:00 +00:00)
        reference_time = np.datetime64("2018-07-00T00:00:00Z")
        counter_hours = 0.0

        # Iterate through each NetCDF file
        for nc_file in nc_files:
            # Read NetCDF file
            ds = xr.open_dataset(nc_file)
            x_wind = ds['u_10m_gr'].values

            # Get the time dimension
            times = ds.time.values

            # Iterate over each time step and write data
            for i, time_val in enumerate(times):
                # Calculate the time difference from the reference time in hours
                time_difference_hours = (time_val - reference_time) / np.timedelta64(1, 'h')

                # Format time as desired
                formatted_time = f"TIME = {time_difference_hours:.1f} hours since {str(reference_time)[:-6].replace('T', ' ')}:00:00 +00:00\n"
                f.write(formatted_time)

                # Update the counter by 1 hour for the next time step
                counter_hours += 1.0

                # Write data with three decimal places for each time step
                formatted_values = ['%.3f' % value if value >= 0 else '%.2f' % value for value in x_wind[i].flatten()]

                # Write values, starting a new line every 200 values
                for j, val in enumerate(formatted_values, start=1):
                    f.write(val)
                    if j % 200 == 0 or j == len(formatted_values):
                        f.write('\n')
                    else:
                        f.write('\t')

            print(f"Data from {nc_file} written to {output_file}")

    print(f"Conversion completed. Delft3D file saved as: {output_file}")

# Directory containing NetCDF files
nc_dir = r'C:\Python_scripts\WRF\cartesian'

# List of NetCDF files
nc_files = [os.path.join(nc_dir, file) for file in os.listdir(nc_dir) if file.endswith('.nc')]

# Output Delft3D file with a unique name
output_file = os.path.join(nc_dir, 'x_wind.amu')

# Convert multiple NetCDF files to a single Delft3D file
convert_to_delft3d(nc_files, output_file)

# Function to convert NetCDF files to Delft3D format
def convert_to_delft3d(nc_files, output_file):
    # Read the first NetCDF file to get dimensions and time information
    ds_first = xr.open_dataset(nc_files[0])
    rh_2m_first = ds_first['rh_2m'].values
    time, south_north, west_east = rh_2m_first.shape

    # Open output file for writing
    with open(output_file, 'w') as f:
        # Write adjusted header
        f.write("FileVersion     = 1.03\n")
        f.write("filetype        = meteo_on_equidistant_grid\n")
        f.write("NODATA_value    = -999.000\n")
        f.write(f"n_cols          = {west_east}\n")
        f.write(f"n_rows          = {south_north}\n")
        f.write("grid_unit       = m\n")  # Adjusted grid unit to meters
        f.write("x_llcenter      = 560791.9\n")  # Adjusted x_llcenter value
        f.write("y_llcenter      = 4432953.0\n")  # Adjusted y_llcenter value
        f.write("dx              = 206.5\n")
        f.write("dy              = 206.5\n")
        f.write("n_quantity      = 1\n")
        f.write("quantity1       = relative_humidity\n")
        f.write("unit1           = %\n")

        # Initialize the counter to start at the desired reference time (2009-11-01 10:00:00 +00:00)
        reference_time = np.datetime64("2018-07-00T00:00:00Z")
        counter_hours = 0.0

        # Iterate through each NetCDF file
        for nc_file in nc_files:
            # Read NetCDF file
            ds = xr.open_dataset(nc_file)
            rh_2m = ds['rh_2m'].values

            # Get the time dimension
            times = ds.time.values

            # Iterate over each time step and write data
            for i, time_val in enumerate(times):
                # Calculate the time difference from the reference time in hours
                time_difference_hours = (time_val - reference_time) / np.timedelta64(1, 'h')

                # Format time as desired
                formatted_time = f"TIME = {time_difference_hours:.1f} hours since {str(reference_time)[:-6].replace('T', ' ')}:00:00 +00:00\n"
                f.write(formatted_time)

                # Update the counter by 1 hour for the next time step
                counter_hours += 1.0

                # Write data with three decimal places for each time step
                formatted_values = ['%.3f' % value if value >= 0 else '%.2f' % value for value in rh_2m[i].flatten()]

                # Write values, starting a new line every 200 values
                for j, val in enumerate(formatted_values, start=1):
                    f.write(val)
                    if j % 200 == 0 or j == len(formatted_values):
                        f.write('\n')
                    else:
                        f.write('\t')

            print(f"Data from {nc_file} written to {output_file}")

    print(f"Conversion completed. Delft3D file saved as: {output_file}")

# Directory containing NetCDF files
nc_dir = r'C:\Python_scripts\WRF\cartesian'

# List of NetCDF files
nc_files = [os.path.join(nc_dir, file) for file in os.listdir(nc_dir) if file.endswith('.nc')]

# Output Delft3D file with a unique name
output_file = os.path.join(nc_dir, 'rel_humidity.amr')

# Convert multiple NetCDF files to a single Delft3D file
convert_to_delft3d(nc_files, output_file)

# Function to convert NetCDF files to Delft3D format
def convert_to_delft3d(nc_files, output_file):
    # Read the first NetCDF file to get dimensions and time information
    ds_first = xr.open_dataset(nc_files[0])
    precip_c_first = ds_first['precip_c'].values
    time, south_north, west_east = precip_c_first.shape

    # Open output file for writing
    with open(output_file, 'w') as f:
        # Write adjusted header
        f.write("FileVersion     = 1.03\n")
        f.write("filetype        = meteo_on_equidistant_grid\n")
        f.write("NODATA_value    = -999.000\n")
        f.write(f"n_cols          = {west_east}\n")
        f.write(f"n_rows          = {south_north}\n")
        f.write("grid_unit       = m\n")  # Adjusted grid unit to meters
        f.write("x_llcenter      = 560791.9\n")  # Adjusted x_llcenter value
        f.write("y_llcenter      = 4432953.0\n")  # Adjusted y_llcenter value
        f.write("dx              = 206.5\n")
        f.write("dy              = 206.5\n")
        f.write("n_quantity      = 1\n")
        f.write("quantity1       = precipitation\n")
        f.write("unit1           = mm/h\n")

        # Initialize the counter to start at the desired reference time (2009-11-01 10:00:00 +00:00)
        reference_time = np.datetime64("2018-07-00T00:00:00Z")
        counter_hours = 0.0

        # Iterate through each NetCDF file
        for nc_file in nc_files:
            # Read NetCDF file
            ds = xr.open_dataset(nc_file)
            precip_c = ds['precip_c'].values

            # Get the time dimension
            times = ds.time.values

            # Iterate over each time step and write data
            for i, time_val in enumerate(times):
                # Calculate the time difference from the reference time in hours
                time_difference_hours = (time_val - reference_time) / np.timedelta64(1, 'h')

                # Format time as desired
                formatted_time = f"TIME = {time_difference_hours:.1f} hours since {str(reference_time)[:-6].replace('T', ' ')}:00:00 +00:00\n"
                f.write(formatted_time)

                # Update the counter by 1 hour for the next time step
                counter_hours += 1.0

                # Write data with three decimal places for each time step
                formatted_values = ['%.3f' % value if value >= 0 else '%.2f' % value for value in precip_c[i].flatten()]

                # Write values, starting a new line every 200 values
                for j, val in enumerate(formatted_values, start=1):
                    f.write(val)
                    if j % 200 == 0 or j == len(formatted_values):
                        f.write('\n')
                    else:
                        f.write('\t')

            print(f"Data from {nc_file} written to {output_file}")

    print(f"Conversion completed. Delft3D file saved as: {output_file}")

# Directory containing NetCDF files
nc_dir = r'C:\Python_scripts\WRF\cartesian'

# List of NetCDF files
nc_files = [os.path.join(nc_dir, file) for file in os.listdir(nc_dir) if file.endswith('.nc')]

# Output Delft3D file with a unique name
output_file = os.path.join(nc_dir, 'precipitation.ampr')

# Convert multiple NetCDF files to a single Delft3D file
convert_to_delft3d(nc_files, output_file)

# Function to convert NetCDF files to Delft3D format
def convert_to_delft3d(nc_files, output_file):
    # Read the first NetCDF file to get dimensions and time information
    ds_first = xr.open_dataset(nc_files[0])
    t_2m_first = ds_first['T_2m'].values
    time, south_north, west_east = t_2m_first.shape

    # Open output file for writing
    with open(output_file, 'w') as f:
        # Write adjusted header
        f.write("FileVersion     = 1.03\n")
        f.write("filetype        = meteo_on_equidistant_grid\n")
        f.write("NODATA_value    = -999.000\n")
        f.write(f"n_cols          = {west_east}\n")
        f.write(f"n_rows          = {south_north}\n")
        f.write("grid_unit       = m\n")  # Adjusted grid unit to meters
        f.write("x_llcenter      = 560791.9\n")  # Adjusted x_llcenter value
        f.write("y_llcenter      = 4432953.0\n")  # Adjusted y_llcenter value
        f.write("dx              = 206.5\n")
        f.write("dy              = 206.5\n")
        f.write("n_quantity      = 1\n")
        f.write("quantity1       = air_temperature\n")
        f.write("unit1           = Celsius\n")
        
        # Initialize the counter to start at the desired reference time (2009-11-01 10:00:00 +00:00)
        reference_time = np.datetime64("2018-07-02T00:00:00Z")
        counter_hours = 0.0

        # Iterate through each NetCDF file
        for nc_file in nc_files:
            # Read NetCDF file
            ds = xr.open_dataset(nc_file)
            t_2m = ds['T_2m'].values

            # Get the time dimension
            times = ds.time.values

            # Iterate over each time step and write data
            for i, time_val in enumerate(times):
                # Calculate the time difference from the reference time in hours
                time_difference_hours = (time_val - reference_time) / np.timedelta64(1, 'h')

                # Format time as desired
                formatted_time = f"TIME = {time_difference_hours:.1f} hours since {str(reference_time)[:-6].replace('T', ' ')}:00:00 +00:00\n"
                f.write(formatted_time)

                # Update the counter by 1 hour for the next time step
                counter_hours += 1.0

                # Write data with three decimal places for each time step
                formatted_values = ['%.3f' % value if value >= 0 else '%.2f' % value for value in t_2m[i].flatten()]

                # Write values, starting a new line every 200 values
                for j, val in enumerate(formatted_values, start=1):
                    f.write(val)
                    if j % 200 == 0 or j == len(formatted_values):
                        f.write('\n')
                    else:
                        f.write('\t')

            print(f"Data from {nc_file} written to {output_file}")

    print(f"Conversion completed. Delft3D file saved as: {output_file}")

# Directory containing NetCDF files
nc_dir = r'C:\Python_scripts\WRF\cartesian'

# List of NetCDF files
nc_files = [os.path.join(nc_dir, file) for file in os.listdir(nc_dir) if file.endswith('.nc')]

# Output Delft3D file with a unique name
output_file = os.path.join(nc_dir, 'air_temperature.amt')

# Convert multiple NetCDF files to a single Delft3D file
convert_to_delft3d(nc_files, output_file)

# Function to convert NetCDF files to Delft3D format
def convert_to_delft3d(nc_files, output_file):
    # Read the first NetCDF file to get dimensions and time information
    ds_first = xr.open_dataset(nc_files[0])
    p_sfc_first = ds_first['p_sfc'].values
    time, south_north, west_east = p_sfc_first.shape

    # Open output file for writing
    with open(output_file, 'w') as f:
        # Write adjusted header
        f.write("FileVersion     = 1.03\n")
        f.write("filetype        = meteo_on_equidistant_grid\n")
        f.write("NODATA_value    = -999.000\n")
        f.write(f"n_cols          = {west_east}\n")
        f.write(f"n_rows          = {south_north}\n")
        f.write("grid_unit       = m\n")  # Adjusted grid unit to meters
        f.write("x_llcenter      = 560791.9\n")  # Adjusted x_llcenter value
        f.write("y_llcenter      = 4432953.0\n")  # Adjusted y_llcenter value
        f.write("dx              = 206.5\n")
        f.write("dy              = 206.5\n")
        f.write("n_quantity      = 1\n")
        f.write("quantity1       = air_pressure\n")
        f.write("unit1           = mbar\n")
        
        # Initialize the counter to start at the desired reference time (2009-11-01 10:00:00 +00:00)
        reference_time = np.datetime64("2018-07-00T00:00:00Z")
        counter_hours = 0.0

        # Iterate through each NetCDF file
        for nc_file in nc_files:
            # Read NetCDF file
            ds = xr.open_dataset(nc_file)
            p_sfc = ds['p_sfc'].values

            # Get the time dimension
            times = ds.time.values

            # Iterate over each time step and write data
            for i, time_val in enumerate(times):
                # Calculate the time difference from the reference time in hours
                time_difference_hours = (time_val - reference_time) / np.timedelta64(1, 'h')

                # Format time as desired
                formatted_time = f"TIME = {time_difference_hours:.1f} hours since {str(reference_time)[:-6].replace('T', ' ')}:00:00 +00:00\n"
                f.write(formatted_time)

                # Update the counter by 1 hour for the next time step
                counter_hours += 1.0

                # Write data with three decimal places for each time step
                formatted_values = ['%.3f' % value if value >= 0 else '%.2f' % value for value in p_sfc[i].flatten()]

                # Write values, starting a new line every 200 values
                for j, val in enumerate(formatted_values, start=1):
                    f.write(val)
                    if j % 200 == 0 or j == len(formatted_values):
                        f.write('\n')
                    else:
                        f.write('\t')

            print(f"Data from {nc_file} written to {output_file}")

    print(f"Conversion completed. Delft3D file saved as: {output_file}")

# Directory containing NetCDF files
nc_dir = r'C:\Python_scripts\WRF\cartesian'

# List of NetCDF files
nc_files = [os.path.join(nc_dir, file) for file in os.listdir(nc_dir) if file.endswith('.nc')]

# Output Delft3D file with a unique name
output_file = os.path.join(nc_dir, 'air_pressure.amp')

# Convert multiple NetCDF files to a single Delft3D file
convert_to_delft3d(nc_files, output_file)