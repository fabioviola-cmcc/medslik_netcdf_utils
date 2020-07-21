#!/usr/bin/python

# requirements
import sys
import getopt
import logging
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
import warnings


# help
def printHelp():
    pass


# getRangeIndexes
def getRangeIndexes(arr, var_min, var_max):
    return np.where((arr >= var_min) & (arr <= var_max))[0]



# main
if __name__ == "__main__":

    # disable warnings
    warnings.filterwarnings("ignore")

    #############################################################
    #
    # Configure logger
    #
    #############################################################

    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger('wind_viewer')
    logger.setLevel(logging.DEBUG)
    logging.getLogger("matplotlib").setLevel(logging.CRITICAL)


    #############################################################
    #
    # read input params
    #
    #############################################################

    logger.debug("Reading input parameters")

    # initialize variables
    scale = None
    outputFile = None
    onlySave = False
    arrowSize=100

    # read params
    try:
        options, rem = getopt.getopt(sys.argv[1:], 'i:hf:p:t:', ['inputFile=', 'help', 'latVariable=', 'lonVariable=', 'plotVariables=', 'timestep=', 'scale=', 'arrowSize=', 'onlySave', 'outputFile='])
    
        for opt, arg in options:
            if opt in ('-i', '--inputFile'):
                inputFile = arg
            elif opt in ('--latVariable'):                
                latVar = arg
            elif opt in ('--lonVariable'):                
                lonVar = arg
            elif opt in ('-p', '--plotVariables'):
                plotVar1, plotVar2 = arg.split(",")
            elif opt in ('-t', '--timestep'):
                ts = int(arg)
            elif opt == "--scale":
                scale = int(arg)
            elif opt in "--outputFile":
                outputFile = arg
            elif opt in "--arrowSize":
                arrowSize = int(arg)
            elif opt in "--onlySave":
                onlySave = True
            elif opt in ('-h', '--help'):
                showHelp(logger)
                sys.exit(0)
    
    except getopt.GetoptError:
        showHelp(logger)
        sys.exit(1)

        
    #############################################################
    #
    # read input file and extract data
    #
    #############################################################

    # debug print
    logger.debug("Reading input file and extracting data")
    
    # open input file
    ds = Dataset(inputFile)

    # select the timestep
    u = ds.variables[plotVar1][ts,:,:]
    v = ds.variables[plotVar2][ts,:,:]

    # read u and v
    umax = ds.variables[plotVar1][ts,:,:]
    vmax = ds.variables[plotVar2][ts,:,:]

    # read lat and lon
    lons = ds.variables[lonVar][:]
    lats = ds.variables[latVar][:]
    lat_indexes = getRangeIndexes(lats, lats.min(), lats.max())
    lon_indexes = getRangeIndexes(lons, lons.min(), lons.max())    
    lats_sel = lats[lat_indexes]
    lons_sel = lons[lon_indexes]
    xx, yy = np.meshgrid(lons_sel, lats_sel)

    #############################################################
    #
    # plotting data
    #
    #############################################################
    
    # plot the map
    m = Basemap(llcrnrlat=lats.min(), urcrnrlat=lats.max(),
                llcrnrlon=lons.min(), urcrnrlon=lons.max(),                    
                resolution='l')

    # add coastlines, states, and country boundaries
    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()
    
    # add color
    m.fillcontinents(color='coral',lake_color='aqua')

    # set the title
    plt.title(inputFile)
    
    # color the sea
    lon, lat = np.meshgrid(lons, lats)
    xi, yi = m(lon, lat)
    cs = m.pcolor(xi, yi, np.squeeze(umax))
    
    # draw arrows
    X = lons[::scale]
    Y = lats[::scale]
    UU = umax[::scale,::scale]
    VV = vmax[::scale,::scale]
    m.quiver(X, Y, UU, VV, scale=arrowSize)

    # save or show    
    if outputFile:
        logger.debug("Saving output to: %s" % outputFile)
        plt.savefig(outputFile)
    if not(onlySave):
        plt.show()
