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

## Visualise winds

To visualise winds, the user must specify the input file, lat and lon variables, u10 and v10 variables, and the desired timestep. Optionally, he can require to save data to a file specifying its filename (and also to only save to file with `--onlySave`). `scale` and `arrowsize` are optional and allow users to customise the output.

```bash
python view_nc_winds.py --inputFile=samples/winds_20190805.nc --latVariable=lat --lonVariable=lon --plotVariables=U10M,V10M --timestep=1 --scale=3 --arrowSize=300 --outputFile=output.png
```