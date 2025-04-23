# WRF to Delft3D Conversion Tools

This repository includes a collection of Python scripts to:
- Convert cloudiness, wind, radiation, and other WRF model outputs into Delft3D format
- Perform spatial averaging and point extraction
- Transform NetCDF coordinates from geographic to UTM

## Scripts
- `D3D_addin_CLD.py`: Converts cloudiness `.npz` files into Delft3D format
- `D3D_addin_parameters.py`: Converts multiple variables (e.g., radiation, humidity) to `.amc`, `.amt`, etc.
- `wrf_degrees2cartesian.py`: Converts WRF `.nc` files from lat/lon to Cartesian coordinates
- `extract_values_wrf.py`: Extracts and averages weather parameters from `.nc` files within defined areas

## Data Availability
The data used for these scripts is **not included** in the repository due to size and licensing.  
**Please contact `your_email@example.com` to request access.**

## Requirements
These scripts require:
- `xarray`
- `numpy`
- `matplotlib`
- `netCDF4`
- `pyproj`

## License
GNU General Public License v2.0
