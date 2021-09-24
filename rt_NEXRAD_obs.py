########################################################################################################################
# This script is to be used with realtime sticknet obs
# The purpose of this script is to create reflectivity and velocity sticknet station plots, 
# and reflectivity and thetaV and thetaE objective analysis plots
#
# This script should be run every 10 MINUTES during a project
#
# Initially developed by A. Hill, Overhauled by J. McDonald 2021
########################################################################################################################

import matplotlib
matplotlib.use('Agg')

import cartopy.crs as ccrs 
import cartopy.feature as cfeature
import cartopy
from metpy.plots import USCOUNTIES
from metpy.plots import  StationPlot
from metpy.interpolate import interpolate_to_grid

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.axes as maxes
import matplotlib.patheffects as path_effects
import cmocean


import pandas as pd
from netCDF4 import num2date
import numpy as np
import pyart;
import glob
import pytz
import datetime as dt
from collections import OrderedDict

# IMPORTANT
from info import TZ, savedir, probe_locs

import SNmods as snmods
from functions import calc_dewpoint,calc_thetae,calc_thetav,C_to_F,calc_mslp,convert_wind,parse_currtime
from functions_plotting import plot_logo, scale_bar, drop_nans_var, homeyer_rainbow_ramp
from functions_radar import get_radar_from_aws,_nearestDate, get_radar_ID

# get refl colormap
homeyer_rainbow_ramp()

# read in cartopy information    
crs = ccrs.PlateCarree()
# Get data to plot state and province boundaries
states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lakes',
        scale='10m',
        facecolor='none')


def add_radar(ax, radar, field):
    ''' 
    returns information for title
    '''
    if field == 'reflectivity':
        vmin = 0
        vmax = 70
        sweep = 0
        cbar_label = 'Reflectivity (dBZ)'
        cmap = 'pyart_HomeyerRainbow_ramp'
        
    elif field == 'velocity':
        vmin  = -45
        vmax  = 45 
        sweep = 1
        cbar_label = 'Radial Velocity (m/s)'
        cmap = cmocean.cm.balance
        
        
    display = pyart.graph.RadarMapDisplay(radar)
    lat_0 = display.loc[0]
    lon_0 = display.loc[1]

    ### Grab indeces and times from radar object to make nice title
    index_at_start = radar.sweep_start_ray_index['data'][0]
    time_at_start_of_radar = num2date(radar.time['data'][index_at_start], 
                                      radar.time['units'])
    index_at_start = radar.sweep_start_ray_index['data'][0]
    ts = num2date(radar.time['data'][index_at_start], radar.time['units'])

    rdate = dt.datetime(ts.year, ts.month, ts.day, ts.hour, ts.minute, ts.second)
    timezone = pytz.timezone(TZ) #if this errors, replace TZ with "US/Central" or "US/Eastern", etc
    local_time = timezone.fromutc(rdate)
    fancy_date_string = local_time.strftime('%A, %B %d at %I:%M %p %Z')
    fancy_date_string_utc = time_at_start_of_radar.strftime('%Y-%m-%d %H:%M UTC')
      
    display.plot_ppi_map(field, sweep, cmap=cmap, ax=ax, vmin=vmin, vmax=vmax, 
                         colorbar_flag=False, title_flag=False,
                        lat_lines=[0], lon_lines=[0])


    ### gives stylistic control over colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.15, axes_class=maxes.Axes)
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    cbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), pad=0.02, cax=cax, orientation='vertical')
    cbar.set_label(label=cbar_label, size=16, weight='bold', labelpad=12)
    cbar.ax.tick_params(labelsize=15)

    ### Mark the radar
    ax.scatter(lon_0, lat_0, marker='o',s=30, c='r', ec='k', label=radar.metadata['instrument_name'])
    
    return fancy_date_string, fancy_date_string_utc, radar.metadata['instrument_name']


def radar_and_sticknet(lats, lons, data_df, mslp, dew_plot, IDs, u, v, obtime, radar, product='reflectivity'):

    ### find extent of plot
    # exact center of plot
    clat = np.amin(lats)+ abs(np.amax(lats) - np.amin(lats))/2
    clon = np.amin(lons) + abs(np.amin(lons) - np.amax(lons))/2

    # NOTE: change these hardcoded values if you want to change relative domain size
    # larger (smaller) numbers = larger (smaller) domain
    dlat = 0.78 * abs(np.amax(lats) - np.amin(lats))
    dlon = 0.78 * abs(np.amin(lons) - np.amax(lons))
    
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
    ax.add_feature(USCOUNTIES.with_scale('20m'), alpha=0.4, linewidth=0.2)

    #### RADAR
    rtime, rtime_utc, radarname = add_radar(ax, radar, product)

    ##### STICKNET
    ax.plot(lons,lats,marker='s',color='k',markersize=5, linewidth=0)
    
    pe= [path_effects.Stroke(linewidth=3, foreground='white'), path_effects.Normal()]

    # use metpy to plot T, Td, MSLP (coded), and the 4-letter identifiers
    stationplot = StationPlot(ax, lons, lats,clip_on=True, transform=crs, fontsize=9, weight='bold')
    
    stationplot.plot_parameter((-1.5,1), data_df['T'].values, color='#b30000', 
                               formatter='0.1f', path_effects=pe)
    stationplot.plot_parameter((-1.5,-1), dew_plot, color='darkgreen', 
                               formatter='0.1f', path_effects=pe)
    stationplot.plot_parameter((1.5,1), mslp, formatter=lambda v: format(10 * v, '.0f')[-3:],
                               path_effects=pe)
    stationplot.plot_text((1.5, -.9), IDs, fontsize=8, weight='bold', path_effects=pe)

    # Add wind barbs, dropping bad vals, cuz WS can == -999
    ax.barbs(lons, lats, u, v, length=7.5,sizes={'emptybarb':.18}, lw=1.2,zorder=20)

    ### title
    ax.set_title('{0} {1}\n{2} ({3})\nStickNet Observations at {4}'.format(radarname,
                                product.capitalize(), rtime, rtime_utc, obtime.strftime('%D %H:%M UTC')),
                 fontsize=16, y=1.01)

    ### Plot the TTU logo, have it update location based on shape of base map
    plt.draw()
    plot_logo(fig, ax,alpha=0.7)


    ### Scale Bar
    center = np.average([ax.get_position().x0, ax.get_position().x1])
    scale_len = np.ceil((dlon*10)/3.5)*10
    if scale_len < 10: scale_len = 10
    scale_bar(fig, ax, length=int(scale_len), location=(center, ax.get_position().y0),fontsize=10)
    
    ### legend
    ax.legend(loc=1, fontsize=11,handletextpad=0.01,borderpad=0.3)
    
    # save figure 
    nametime = parse_currtime()[1] # note that this may be different from last time on plot
                                 # THATS OKAY... it will make data drops more apparent!
    plt.savefig('{0}{1}_station_{2}.png'.format(savedir,product,nametime),dpi=300,bbox_inches = 'tight')
    plt.close()

    
def radar_and_oban(lats, lons, variable, name, savename, obtime, radar, line_color='C0',dK=1, product='reflectivity'):

    ### find extent of plot
    # exact center of plot
    clat = np.amin(lats)+ abs(np.amax(lats) - np.amin(lats))/2
    clon = np.amin(lons) + abs(np.amin(lons) - np.amax(lons))/2

    # NOTE: change these hardcoded values if you want to change relative domain size
    # larger (smaller) numbers = larger (smaller) domain
    dlat = 0.78 * abs(np.amax(lats) - np.amin(lats))
    dlon = 0.78 * abs(np.amin(lons) - np.amax(lons))
    
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
    ax.add_feature(USCOUNTIES.with_scale('20m'), alpha=0.4, linewidth=0.2)

    #### RADAR
    rtime, rtime_utc, radarname = add_radar(ax, radar, product)
    
    ### OBAN
    ### interpolate data. Hres is 0.05 degrees 
    varn, latsn, lonsn, FLAG = drop_nans_var(variable, lats, lons)
    gx,gy,img1 = interpolate_to_grid(np.asarray(lonsn),np.asarray(latsn),np.asarray(varn),\
                         interp_type='natural_neighbor',hres=.05)
    levels = np.arange(int(np.floor(np.amin(variable))),int(np.ceil(np.amax(variable)))+1,dK)
    ax.contour(gx, gy, img1, colors='w', levels=levels, linewidths=3)
    c=ax.contour(gx, gy, img1, colors=line_color, levels=levels, linewidths=1.5)
    
    # labels
    cl = ax.clabel(c, fontsize=10,inline=1, inline_spacing=2,fmt='%i', 
              rightside_up=True, use_clabeltext=True)
    for t in cl:
        t.set_weight('bold')
        t.set_path_effects([path_effects.Stroke(linewidth=3, foreground='white'), path_effects.Normal()])
        
        
        
    # add box around the grid for reference
    x0, x1 = np.min(gx), np.max(gx)
    y0, y1 = np.min(gy), np.max(gy)
    
    ax.plot([x0,x0,x1,x1,x0],[y0,y1, y1, y0,y0], color='k', alpha=0.2, linewidth=3, linestyle='--')

    ### title
    ax.set_title('{0} {1}\n{2} ({3})\nStickNet {4} Analysis at {5}'.format(radarname,
                                product.capitalize(), rtime, rtime_utc, name, obtime.strftime('%D %H:%M UTC')),
                 fontsize=16, y=1.01)

    ### Plot the TTU logo, have it update location based on shape of base map
    plt.draw()
    plot_logo(fig, ax, alpha=0.7)


    ### Scale Bar
    center = np.average([ax.get_position().x0, ax.get_position().x1])
    scale_len = np.ceil((dlon*10)/3.5)*10
    if scale_len < 10: scale_len = 10
    scale_bar(fig, ax, length=int(scale_len), location=(center, ax.get_position().y0),fontsize=10)
    
    ### legend
    ax.legend(loc=1, fontsize=11,handletextpad=0.01,borderpad=0.3)
    
    ## warn if variable data is missing!
    if FLAG != 0:
        plt.draw()
        ax.text(center,0.98,  f'{FLAG} StickNets missing, plot may be inaccurate', 
                fontsize=10,color='r',ha='center', va='center', transform=ax.transAxes)

    # save figure 
    nametime = parse_currtime()[1] # note that this may be different from last time on plot
                                 # THATS OKAY... it will make data drops more apparent!
    plt.savefig('{0}{1}_{2}_{3}.png'.format(savedir,product,savename,nametime),dpi=300,bbox_inches = 'tight')
    plt.close()


################################    
# pull in data

analysis_time = dt.datetime(2017,4,30,19,9) #dt.datetime.utcnow()
probes = list(probe_locs.keys())

# Grab the data from files in filedir we want to plot
data_df = snmods.get_sticknet_data(analysis_time,analysis_time,dataset='latest',probes=probes)

if not data_df.empty: # only works if you have data

    # pull sticknet locations and time
    lats =       data_df['Lat'].values
    lons =       data_df['Lon'].values
    IDs =        [ids[1:] for ids in data_df.index.values]
    obtime = data_df['date'][0]

    # calculate variables
    dewp_c = calc_dewpoint(data_df['T'].values,data_df['RH'].values)
    thetae = calc_thetae(data_df['T'].values,dewp_c, data_df['P'].values)
    thetav = calc_thetav(data_df['T'].values,dewp_c, data_df['P'].values)
    mslp = calc_mslp(data_df['T'].values, data_df['P'].values, data_df['Elevation'].values)
    u,v = convert_wind(data_df['WS'].values,data_df['WD'].values)

    # find appropriate radar, and pull in data
    clat = np.amin(lats)+ abs(np.amax(lats) - np.amin(lats))/2
    clon = np.amin(lons) + abs(np.amin(lons) - np.amax(lons))/2
    station, city = get_radar_ID(clat, clon)
    print('Radar: {}, from {}'.format(station, city))

    radar_namelist, radar_list = get_radar_from_aws(station, analysis_time, analysis_time)
    radar = radar_list[0]


    # radar products and station plots - this will be of Lubbock if no data is available
    radar_and_sticknet(lats, lons, data_df, mslp, dewp_c, IDs, u, v, obtime, radar, product='reflectivity')
    radar_and_sticknet(lats, lons, data_df, mslp, dewp_c, IDs, u, v, obtime, radar, product='velocity')

    # reflectivity and oban plots
    if len(data_df) >= 4: # at least 4 pts needed for interpolation
        radar_and_oban(lats, lons, thetav, r'${\theta_v}$','TV', obtime, radar, line_color='#002775', dK=1) 
        radar_and_oban(lats, lons, thetae,r'${\theta_e}$', 'TE', obtime, radar, line_color='#014522', dK=1)

else:
    print('No data at {}'.format(analysis_time.strftime('%D %H:%M UTC')))
    