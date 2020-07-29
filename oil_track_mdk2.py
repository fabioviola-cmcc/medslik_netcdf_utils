# requirements
import os
import pdb
import sys
import json
import numpy
import imageio
import netCDF4
import warnings
import traceback
import numpy as np
from datetime import  *
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

# muting warnings
warnings.filterwarnings("ignore")


"""
Application:
    This script has been developed to plot oil concentrations found at the 
    sea surface based on MEDSLIK II outputs. The outputs are png figures.
"""

def oil_track(sNcFile, ds, time_index):

    # load ncfile
    ncfile = netCDF4.Dataset(sNcFile,'r')

    # load start position of the spill
    y0 = ncfile.variables['non_evaporative_volume'].initial_position_y
    x0 = ncfile.variables['non_evaporative_volume'].initial_position_x
    print ("Spill location = " + str(x0) + "W ::::: " + str(y0) + "N ")

    # variable extraction
    lats = ncfile.variables['latitude'][time_index,:]
    lons = ncfile.variables['longitude'][time_index,:]
        
    # generate output grid
    grid_min_longitude = numpy.min(lons)-ds
    grid_min_latitude = numpy.min(lats)-ds
    grid_max_longitude = numpy.max(lons)+ds
    grid_max_latitude = numpy.max(lats)+ds

    try:
        x_points = numpy.arange(grid_min_longitude+ds/2,grid_max_longitude,ds)
    except:
        pdb.set_trace()
    y_points = numpy.arange(grid_min_latitude+ds/2,grid_max_latitude,ds)
    box = numpy.zeros((len(y_points),len(x_points)))
    area = (ds*110)**2

    # conversion factor - barrel to tonnes
    oil_density = ncfile.variables['non_evaporative_volume'].oil_density
    #parcel_volume = ncfile.variables['non_evaporative_volume'].volume_of_parcel
    rbm3=0.158987
    barrel2tonnes=1/(rbm3*(oil_density/1000))

    # extract variables of interest
    particle_status = ncfile.variables['particle_status'][time_index,:]
    evaporative_volume = ncfile.variables['evaporative_volume'][time_index,:]
    non_evaporative_volume = ncfile.variables['non_evaporative_volume'][time_index,:]

    # Particle status guide
#        is=0 parcel not released
#        is=1 in the spreading surface slick
#        is=2 on surface but not spreading
#        is=3 dispersed into water column
#        is=-nsg beached on shore segment number nsg

    iNoise=numpy.logical_or(particle_status <= 0, particle_status > 2).nonzero()[0]
    lats = numpy.delete(lats, (iNoise), axis=0)
    lons = numpy.delete(lons, (iNoise), axis=0)
    evaporative_volume = numpy.delete(evaporative_volume, (iNoise), axis=0)
    non_evaporative_volume = numpy.delete(non_evaporative_volume, (iNoise), axis=0)


    xp=numpy.round((lons-numpy.min(x_points))/ds)
    yp=numpy.round((lats-numpy.min(y_points))/ds)
    total_volume=(evaporative_volume+non_evaporative_volume)/barrel2tonnes

    for aa in range(0,len(xp)):
        box[yp[aa].astype(int),xp[aa].astype(int)] = box[yp[aa].astype(int),xp[aa].astype(int)] + total_volume[aa]/area

    return x0, y0,x_points, y_points, box


################################################################################
# USER INPUTS
################################################################################



# read the configuration file
#
# we need to set:
# - sNcFile: the file containing MEDSLIK II netcdf outputs
# - time_line: time steps of interest (hours by default -- Python counting starts from 0).
#   It may be a single number e.g. [146] or a list of numbers e.g. np.arange(0,15)
# - grid_reolution: used to estimate concentrations (in degrees)
# - output folder: where .png files will be placed

config = None
sNcFile= None
grid_resolution = 0.15/110
time_line = None
output_folder = None

# parse configuration file, if any
if len(sys.argv) > 1:

        configFile = sys.argv[1]
        print("Parsing configuration file %s" % configFile)
        
        with open(configFile) as f:
                try:
                        config = json.load(f)
                        sNcFile = config["oil_track"]["spill_prop_file"]
                        output_folder = config["oil_track"]["output_folder"]
                        st = config["oil_track"]["time_line"]["start"]
                        en = config["oil_track"]["time_line"]["end"]
                        inc = config["oil_track"]["time_line"]["increment"]
                        time_line = numpy.arange(st, en, inc)
                except KeyError:
                        print("Missing entry in configuration file! Abort!")
                        sys.exit()
        print("Done!")

else:
    print("You need to provide the config file!")
    sys.exit(0)
        

# create output folder, if needed
if not os.path.exists(output_folder):
    os.mkdir(output_folder)


################################################################################
# USER INPUTS - OVER!
################################################################################
# From here onwards, the script should do everything pretty much automatic
# bugs/errors are expected and in case you unfortunate enough to find out one,
# feel free to send us comments/corrections.

# plotting loop
cc= 0.

# initialize animation
images = []

# loop over the timesteps
for ii in time_line:

    print("[oil_track_mdk2] -- Timestep: %s" % ii)
    
    # extract values
    print("[oil_track_mdk2] -- Extracting values")
    x0, y0, x_points, y_points, box = oil_track(sNcFile, grid_resolution, ii)

    if cc == 0:

        dt = numpy.max(time_line)
        max_ds = (dt*.8)/110.
        grid_min_longitude = x0 - max_ds
        grid_min_latitude = y0 - max_ds
        grid_max_longitude = x0 + max_ds
        grid_max_latitude = y0 + max_ds

        # initiate coastline
        print("[oil_track_mdk2] -- Initiating Basemap")
        m = Basemap(llcrnrlon=grid_min_longitude-.1,llcrnrlat=grid_min_latitude-.1,\
            urcrnrlon=grid_max_longitude+.101,urcrnrlat=grid_max_latitude+.101,\
            rsphere=(6378137.00,6356752.3142),\
            resolution='f',projection='merc',\
            lat_0=(grid_max_latitude + grid_min_latitude)/2,\
			lon_0=(grid_max_longitude + grid_min_longitude)/2)
        cc = cc + 1        

    # Plot trajectory - set map
    plt.figure(figsize=(7, 7)) # This increases resolution
    
    # Plot trajectory - use only non zero areas
    oiled_grid_points=numpy.argwhere(box!=0)
    X=x_points[oiled_grid_points[:,1]]
    Y=y_points[oiled_grid_points[:,0]]
    x,y=m(X,Y)
    x0,y0=m(x0,y0)
    box_plot=box[oiled_grid_points[:,0],oiled_grid_points[:,1]]

    # Plot coastline
    m.drawcoastlines(linewidth=0.05)
    m.fillcontinents(alpha=0.1)
    
    m.drawmeridians(numpy.arange(grid_min_longitude-.1,grid_max_longitude+.101, \
    	((grid_max_longitude+.101)-(grid_min_longitude-.1))/3),\
    	labels=[0,0,0,1],color='black',linewidth=0.03) # draw parallels
    m.drawparallels(numpy.arange(grid_min_latitude-.1,grid_max_latitude+.101, \
    	((grid_max_latitude+.101)-(grid_min_latitude-.1))/3), \
    	labels=[1,0,0,0],color='black',linewidth=0.03) # draw meridians
        
    # Plot trajectory
    m.scatter(x, y, s=[1.], c=numpy.log(box_plot), vmin=numpy.log(0.05),vmax=numpy.log(1.0), edgecolor='')
    m.plot(x0,y0,'k+',markersize=5)
    
    # Colorbar setup
    ticks = numpy.log([0.05, 0.1, 0.5, 1.0])
    cbar = plt.colorbar(ticks=ticks, format='$%.2f$', orientation='horizontal')
    cbar.ax.set_xticklabels(['0.01','0.05','0.1','0.5','1.0'])    
    cbar.set_label('tons/km2')
    plt.title('Surface oil concentrations for '+ '%03d' % (ii+1) + 'h')
    plt.savefig(output_folder + '/surface_oil_' + '%03d' % (ii+1) + 'h.png',dpi=600,bbox_inches='tight')
    # plt.show()
    plt.close('all')
    images.append(imageio.imread(output_folder + '/surface_oil_' + '%03d' % (ii+1) + 'h.png'))

# generate output gif
print("Generating output gif in %s/anim.gif" % output_folder)
imageio.mimsave(output_folder + '/anim.gif', images, fps=2)
