########################################################################################################################
# This script is to be used with realtime sticknet obs
# The purpose of this script is to create the station plots and to fill out the html file
# that drives the observation table on the website (hence "grid").
#
# This script should be run every 1 MINUTE during a project
#
# Initially developed by A. Hill, Overhauled by J. McDonald 2021
########################################################################################################################
import matplotlib
matplotlib.use('Agg')

import cartopy.crs as ccrs 
import cartopy.feature as cfeature
from metpy.plots import USCOUNTIES
from metpy.plots import  StationPlot

import datetime as dt
import matplotlib.pyplot as plt


from netCDF4 import num2date
import numpy as np
import glob
import SNmods as snmods

# IMPORTANT
from info import probe_locs, savedir

import pandas as pd
from functions import calc_dewpoint,calc_thetae,calc_thetav,C_to_F,calc_mslp,convert_wind, parse_currtime
from functions_plotting import scale_bar, plot_logo
from collections import OrderedDict

# # read in cartopy information    
crs = ccrs.PlateCarree()
# Get data to plot state and province boundaries
states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lakes',
        scale='10m',
        facecolor='none')


# create the station plot and write to html

endtime = dt.datetime(2017,4,30,19,9) #dt.datetime.utcnow()
starttime = endtime 

# find all active StickNets
probes = list(probe_locs.keys())


# open (create) html file for writing
# html_filepath = '/Users/severe/Research/VORTEXSE/data.csv'
# html_file = open(html_filepath,'w')

# write to html AND gret the data we want to plot
# probe ID, date, pressure, temperature, wind speed, wind direction, relative humidity,
data_df = snmods.get_sticknet_data(starttime,endtime, dataset='latest',probes=probes)
         #,html=html_file)
    
if not data_df.empty: # only works if you have data

    lats =       data_df['Lat'].values
    lons =       data_df['Lon'].values
    elevations = data_df['Elevation'].values
    IDs =        [ids[1:] for ids in data_df.index.values]

    dew_plot = calc_dewpoint(data_df['T'].values,data_df['RH'].values)
    mslp = calc_mslp(data_df['T'].values, data_df['P'].values, elevations)
    u,v = convert_wind(data_df['WS'].values,data_df['WD'].values)

    obtime = data_df['date'][0]

    ###### STATION PLOT

    ### find extent of plot
    # exact center of plot
    clat = np.amin(lats)+ abs(np.amax(lats) - np.amin(lats))/2
    clon = np.amin(lons) + abs(np.amin(lons) - np.amax(lons))/2

    # NOTE: change these hardcoded values if you want to change relative domain size
    # larger (smaller) numbers = larger (smaller) domain
    dlat = 0.75 * abs(np.amax(lats) - np.amin(lats))
    dlon = 0.75 * abs(np.amin(lons) - np.amax(lons))
    
    if dlon < .1:
        dlon = .1
    if dlat < .1:
        dlat = .1

    # find corners using the center and the buffers
    north_lat, south_lat = clat+dlat, clat-dlat
    west_lon, east_lon = clon-dlon, clon+dlon


    ### initialize figure
    fig = plt.figure(figsize = [10,10])
    ax = fig.add_subplot(1,1,1, projection=crs)
    ax.set_extent([west_lon, east_lon, north_lat,south_lat], crs )
    ax.add_feature(states_provinces, edgecolor='k', alpha=0.25, linewidth=1)
    ax.add_feature(USCOUNTIES.with_scale('500k'), alpha=0.4, linewidth=0.2)

    # mark locations of SN
    ax.plot(lons,lats,marker='s',color='0.4',markersize=5, linewidth=0)

    # use metpy to plot T, Td, MSLP (coded), and the 4-letter identifiers
    stationplot = StationPlot(ax, lons, lats,clip_on=True, transform=crs, fontsize=10)
    stationplot.plot_parameter((-1.5,1), C_to_F(data_df['T'].values), color='#b30000', formatter='0.1f')
    stationplot.plot_parameter((-1.5,-1), C_to_F(dew_plot), color='darkgreen', formatter='0.1f')
    stationplot.plot_parameter((1.5,1), mslp, formatter=lambda v: format(10 * v, '.0f')[-3:])
    stationplot.plot_text((1.5, -.9), IDs, fontsize=9, weight='bold')

    # Add wind barbs manually do to centering issue with stationplot
    ax.barbs(lons, lats, u, v, length=7.5,sizes={'emptybarb':.18}, lw=0.8)


    ### title
    ax.set_title('Observations at {}'.format(obtime.strftime('%D %H:%M UTC')),
                 fontsize=18, y=1.01, weight='bold', color='0.3')

    ### Plot the TTU logo, have it update location based on shape of base map
    plot_logo(fig, ax)


    ### Scale Bar
    center = np.average([ax.get_position().x0, ax.get_position().x1])
    scale_len = np.ceil((dlon*10)/3.5)*10
    if scale_len < 10: scale_len = 10
    scale_bar(fig, ax, length=int(scale_len), location=(center, ax.get_position().y0),fontsize=10)

    # save figure 
    nametime = parse_currtime()[1] # note that this may be different from obs time on plot
                                   # THATS OKAY... it will make data drops more apparent!
    plt.savefig('{0}station_{1}.png'.format(savedir, nametime),dpi=300,bbox_inches = 'tight')
    plt.close()
    
else:
    print('No data at {}'.format(analysis_time.strftime('%D %H:%M UTC')))
