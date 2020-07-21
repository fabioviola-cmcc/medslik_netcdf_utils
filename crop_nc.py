#!/usr/bin/python3

# requirements
import sys
import getopt
import logging
import numpy as np
from netCDF4 import Dataset


# help
def printHelp():
    print("Mandatory parameters:")
    print(" --inputFile=  to specify the input file")
    print(" --outputFile=  to specify the output file")
    print(" --latDim=  to specify the dimension for the latitude")
    print(" --lonDim=  to specify the dimension fro the longitude")
    print(" --latVar=  to specify the variable holding the latitude")
    print(" --lonVar=  to specify the variable holding the longitude")
    print(" --cropVars3=  to specify a comma-separated list of the 3-dimensional variables to crop")
    print(" --cropVars4=  to specify a comma-separated list of the 4-dimensional variables to crop")
    print(" --bbox= to specify lat1, lon1, lat2, lon2 of the bounding box")

    
# main
if __name__ == "__main__":

    # get a logger
    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger('NetCDF_cropper')
    logger.setLevel(logging.DEBUG)

    # init vars
    cropVars3 = cropVars4 = []
    
    # read input and output files
    try:
        options, rem = getopt.getopt(sys.argv[1:], 'h', ['inputFile=','outputFile=', 'latVar=', 'lonVar=', 'cropVars4=', 'cropVars3=', 'bbox=', 'latDim=', 'lonDim=', 'help'])
        for opt, arg in options:
            if opt == '--inputFile':
                inputFile = arg
            elif opt == '--outputFile':
                outputFile = arg
            elif opt == '--latDim':
                latDim = arg
            elif opt == '--lonDim':
                lonDim = arg
            elif opt == '--latVar':
                latVar = arg
            elif opt == '--lonVar':
                lonVar = arg
            elif opt == '--cropVars4':
                cropVars4 = arg.split(",")
            elif opt == '--cropVars3':
                cropVars3 = arg.split(",")
            elif opt == '--bbox':
                lat1,lon1,lat2,lon2 = arg.split(",")
            elif opt in ('-h', '--help'):
                printHelp(logger)
                sys.exit(0)
    
    except getopt.GetoptError:
        logger.error("Wrong arguments!")
        printHelp(logger)
        sys.exit(1)
        
    # open the input and output files
    logger.debug("Opening input and output files")
    in_ds = Dataset(inputFile, "r")
    out_ds = Dataset(outputFile, "w")

    # read the input file and prepare cropped area
    latitudes = in_ds.variables[latVar][:]    
    ywindow = np.logical_and(latitudes>=float(lat1), latitudes<=float(lat2))
    longitudes = in_ds.variables[lonVar][:]
    xwindow = np.logical_and(longitudes>=float(lon1), longitudes<=float(lon2))
    
    # create dimensions
    for dd in in_ds.dimensions:
        logger.debug("Creating dimension %s" % dd)

        # check if lat or lon
        if dd == latDim:
            out_ds.createDimension(dd, latitudes[ywindow].size)
        elif dd == lonDim:
            out_ds.createDimension(dd, longitudes[xwindow].size)
        else:
            out_ds.createDimension(dd, None)
    
    # create variables
    for vv in in_ds.variables:
        
        # ...if not belonging to these groups
        if (not(vv in cropVars4) and not(vv in cropVars3) and not(vv in [latVar, lonVar])):

            # debug print
            logger.debug("Creating variable %s" % vv)            
            
            # create variable
            out_ds.createVariable(vv, in_ds.variables[vv].datatype, in_ds.variables[vv].dimensions)

            # copy values for all the variables except lat, lon and cropVars3 and 4
            out_ds.variables[vv][:] = in_ds.variables[vv][:].data

    # create variables for lat and lon
    logger.debug("Creating variable %s" % latVar)            
    lv1 = out_ds.createVariable(latVar, in_ds.variables[latVar].datatype, (latDim))    
    lv1[:] = latitudes[ywindow].data
    logger.debug("Creating variable %s" % lonVar)            
    lv2 = out_ds.createVariable(lonVar, in_ds.variables[lonVar].datatype, (lonDim))    
    lv2[:] = longitudes[xwindow].data

    # copy values for 4-dim variables
    for vv in cropVars4:        
        logger.debug("Creating variable %s" % vv)        
        vvv = out_ds.createVariable(vv, in_ds.variables[vv].datatype, in_ds.variables[vv].dimensions)        
        vvv[:] = in_ds.variables[vv][:,:,ywindow,xwindow]
        
    # crop 3-dimensional variables
    for vv in cropVars3:        
        logger.debug("Creating variable %s" % vv)        
        vvv = out_ds.createVariable(vv, in_ds.variables[vv].datatype, in_ds.variables[vv].dimensions)        
        vvv[:] = in_ds.variables[vv][:,ywindow,xwindow]

    # write output file
    logger.debug("Closing input and output files")
    in_ds.close()
    out_ds.close()
    
