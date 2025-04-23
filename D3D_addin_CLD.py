# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 13:11:32 2023

@author: vaspapa
"""

import os
import numpy as np
import matplotlib.pyplot as plt

# Specify the directory containing your NPZ files
npz_directory = 'C:/Python_scripts/WRF/original/'

# Specify the specific NPZ file to plot
specific_npz_filename = 'CLDFRA_2July2018_1400.npz'
specific_npz_file_path = os.path.join(npz_directory, specific_npz_filename)

# Check if the file exists
if os.path.exists(specific_npz_file_path):
    # Load data from the specified NPZ file
    data_npz = np.load(specific_npz_file_path)

    # Print the keys in the NPZ file
    print(f"Keys in {specific_npz_filename}: {list(data_npz.keys())}")

    # Determine the correct key to use
    key_to_plot = 'CLDFRA' if 'CLDFRA' in data_npz else 'arr_0'

    # Plot the data using the correct key
    plt.imshow(data_npz[key_to_plot], cmap='viridis', origin='lower')
    plt.colorbar(label=key_to_plot)
    plt.title(f'Plot of {key_to_plot} from {specific_npz_filename}')
    plt.xlabel('X-axis label')  # Add X-axis label if applicable
    plt.ylabel('Y-axis label')  # Add Y-axis label if applicable
    plt.show()
else:
    print(f"File {specific_npz_filename} not found in the specified directory.")

import os
import numpy as np
import xarray as xr

def convert_to_delft3d(npz_files, output_file):
    # Read the first NPZ file to get dimensions and time information
    data_npz_first = np.load(npz_files[0])
    cloudiness_first = data_npz_first['CLDFRA']
    
    # Get dimensions
    if len(cloudiness_first.shape) == 2:
        time, south_north, west_east = 1, *cloudiness_first.shape
    else:
        time, south_north, west_east = cloudiness_first.shape
    
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
        f.write("quantity1       = cloudiness\n")
        f.write("unit1           = %\n")

        # Initialize the counter to start at the desired reference time (2009-11-01 10:00:00 +00:00)
        reference_time = np.datetime64("2018-07-00T00:00:00Z")
        counter_hours = 0.0

        # Iterate through each NPZ file
        for npz_file in npz_files:
            # Read NPZ file
            data_npz = np.load(npz_file)
            cloudiness = data_npz['CLDFRA']

            # Get the time dimension (assuming a single time step per file)
            times = np.array([reference_time + np.timedelta64(int(counter_hours), 'h')])

            # Iterate over each time step and write data
            for i, time_val in enumerate(times):
                # Calculate the time difference from the reference time in hours
                time_difference_hours = (time_val - reference_time) / np.timedelta64(1, 'h')

                # Format time as desired
                formatted_time = f"TIME = {time_difference_hours:.1f} hours since {str(reference_time)[:-6].replace('T', ' ')}:00:00 +00:00\n"
                f.write(formatted_time)

                # Update the counter by 1 hour for the next time step
                counter_hours += 1.0

                # Multiply values by 100%
                cloudiness_percent = cloudiness.flatten() * 100.0

                # Write data with three decimal places for each time step
                formatted_values = ['%.3f' % val for val in cloudiness_percent]

                # Write values, starting a new line every 200 values
                for j, val in enumerate(formatted_values, start=1):
                    f.write(val)
                    if j % 200 == 0 or j == len(formatted_values):
                        f.write('\n')
                    else:
                        f.write('\t')

            print(f"Data from {npz_file} written to {output_file}")

    print(f"Conversion completed. Delft3D file saved as: {output_file}")

# Directory containing NPZ files
npz_dir = 'C:/Python_scripts/WRF/original/'

# List of NPZ files (assumes that the NPZ files have been renamed to end with '_renamed.npz')
npz_files = [os.path.join(npz_dir, file) for file in os.listdir(npz_dir) if file.endswith('_renamed.npz')]

# Output Delft3D file with a unique name
output_file = os.path.join(npz_dir, 'cloudiness.amc')

# Convert multiple NPZ files to a single Delft3D file
convert_to_delft3d(npz_files, output_file)
