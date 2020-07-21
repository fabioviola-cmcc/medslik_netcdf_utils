# NetCDF Utils

This is a simple set of utilities written and maintained by Fabio Viola to deal with NetCDF and Medslik files.

## Coordinate conversion

This script converts latitude and longitude from degrees (e.g. 40.452213, 17.191939) to degrees, minutes (e.g. 40 deg 27.13 min, 17 deg 11.52 min).

```bash
$ python coords_conv.py 40.452213 17.191939
```

## Crop NetCDF files

To crop a NetCDF file we need the input and output filenames, the name for lat and lon dimensions and variables, the bounding box and the variables to crop. Example:

```bash
python crop_nc.py --inputFile=20190805_U.nc --outputFile=cropped.nc --latVar=nav_lat --lonVar=nav_lon --cropVars4="vozocrtx" --bbox="35,10,45,15" --latDim=y --lonDim=x
```