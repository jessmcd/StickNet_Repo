#####################################################
# Functions used in scripts for the VORTEX SE project
#####################################################

import numpy as np
import datetime as dt
import matplotlib as mpl
import matplotlib.pyplot as plt
#import probe_info
import pandas as pd
import subprocess as sp
import os

def get_winddir_string(ws):
    if ws >= 0.0 and ws < 22.5:
        return 'N'
    elif ws >= 22.5 and ws < 67.5:
        return 'NE'
    elif ws >= 67.5 and ws < 112.5:
        return 'E'
    elif ws >= 112.5 and ws < 157.5:
        return 'SE'
    elif ws >= 157.5 and ws < 202.5:
        return 'S'
    elif ws >= 202.5 and ws < 247.5:
        return 'SW'
    elif ws >= 247.5 and ws < 292.5:
        return 'W'
    elif ws >= 292.5 and ws < 337.5:
        return 'NW'
    else:
        return 'N'


def calc_thetae(T,Td,P):
    """ Calculate equivalent potential temperature from Bolton (1980)
    Inputs: Temperature (T, celcius), Dewpoint (Td, celcius) and Station Pressure (P, hPa)
    Output: Equivalent Potential Temperature (theta_e, Kelvin) """
    e = 6.11*(10**((7.5*Td)/(237.3+Td)))                              # vapor pressure, uses degrees C
    w = 0.622 * e/(P-e)                                               # mixing ratio, uses hPa for pressure variables
    T_K,Td_K = T + 273.15,Td + 273.15                                 # convert T and Td to Kelvin
    Tl = 1.0/(1.0/(Td_K-56.0) + np.log(T_K/Td_K)/800.0) + 56.0        # approximated temperature at LCL (Kelvin)
    theta_l = T_K * ((1000.0/(P-e))**0.2854) * (T_K/Tl)**(0.28*w)     # dry potential temperature at LCL (Kelvin)
    theta_e = theta_l * np.exp(((3036.0/Tl)-1.78)*w*(1.0+0.448*w))    # equivalent potential temp (Kelvin)
    return theta_e

def calc_thetav(T,Td,P):
    """ Calculate theta v from Bolton (1980)
    theta_v = theta (1+0.61w)
    Inputs: Temperature (T, celcius), Dewpoint (Td, celcius), Station Pressure (P, hPa)"""
    e = 6.11*(10**((7.5*Td)/(237.3+Td)))
    w = 0.622 * e/(P-e)
    kappa = 2/7.
    theta = (T+273.15)*((1000/P)**kappa)
    theta_v = theta*(1+0.61*w)
    return theta_v

def convert_wind(ws,dir):
    """ convert wind speed to u and v components (in knots) for plotting wind barbs """
    new_dir = 270-dir
    u = (ws*1.94384)*np.cos(new_dir * np.pi/180)
    v = (ws*1.94384)*np.sin(new_dir * np.pi/180)
    return u,v

# calc station pressure from MSLP, T, H
def calc_station_pressure(p_slp,t,h):
    """ p_slp in mb, t in K, h in meters """
    return p_slp * np.exp(-h/(t+29.263))

def alt_to_mb(mm):
    return mm*33.8637526

def calc_dewpoint(T,RH):
    RH = np.ma.masked_values(RH,-999.9)
    num = np.log(RH/100) + (17.625*T)/(243.04+T)
    denom = 17.625 - np.log(RH/100) - (17.625*T)/(243.04+T)
    return 243.04*num/denom

def calc_windchill(T,V):
    """ Temperature in fahrenheit and wind speed V in miles per hour"""
    return 35.74 + 0.6215*T - 32.75*(V**0.16) + 0.4275*T*(V**0.16)

def calc_heatindex(T,RH):
    """ Temperature in fahrenheit and relative humidity in percent """
    line1 = -42.379 + (2.04901523*T) + (10.14333127*RH)
    line2 = (0.22475541 * T * RH) + (6.83783 * 10**-3 * T**2)
    line3 = (5.481717*10**-1 * RH**2) - (1.22874 * 10**-3 * T**2 * RH)
    line4 = (8.5282*10**-4 * T * RH**2) - (1.99*10**-6 * T**2 * RH**2)
    return line1-line2-line3+line4

def calc_mslp(T,P,h):
    return P*(1-(0.0065*h)/(T+0.0065*h+273.15))**(-5.257)

def C_to_F(temp):
    return np.round(temp*1.8 + 32,decimals=1)

def parse_currtime():
    """ parse the current datetime, in datetime format, to a specific string format for titles """
    currtime = dt.datetime.utcnow()
    #year = currtime.strftime("%Y")
    #month = currtime.strftime("%m")
    #day = currtime.strftime("%d")
    #hour = currtime.strftime("%H")
    #min = currtime.strftime("%M")
    #return "{0}-{1}-{2} {3}:{4} UTC".format(year,month,day,hour,min),currtime.strftime("%Y%m%d%H%M")
    return currtime.strftime("%Y-%m-%d %H:%M UTC"),currtime.strftime("%Y%m%d_%H%M")


# Make meteogram plot
def plot_meteogram(met, probe_id, savedir):
    ''' 
    Plots the meteogram of a single sticknet. Accepts any length of time,
    but is designed to look best for 24 hour periods. 
    -----------
    Inputs:
    met       - a pandas dictionary containing Temp ['T'], RH ['RH'], Pressure ['P'],
                windspeed ['WS'], 3 second wind gust ['WSMAX'], and wind dir ['WD'] 
                for a single StickNet.
    elevation - the elevation in meters of the Sticknet. Used to calculated 
                MSLP, ThetaE, and ThetaV
    probe_id  - 4 letter/number station identifier.
    -----------
    outputs: N/A. Saves meteogram to hardcoded location with hardcoded date in name. 
             Edit this for post-processing.  

    Upgraded J. McDonald 2021. 
    '''

    # pull and format data
    elevation = met.attrs['elevation']
    tempf     = C_to_F(met['T'].values)
    dewp      = calc_dewpoint(met['T'].values,met['RH'].values)
    dewpf     = C_to_F(dewp)
    thetae    = calc_thetae(met['T'].values, dewp, met['P'].values)
    thetav    = calc_thetav(met['T'].values, dewp, met['P'].values)
    mslp      = calc_mslp(met['T'].values,met['P'].values, elevation)
    dates     = pd.to_datetime(met.index)
    ws        = met['WS']*1.94384     # convert from m/s to kt
    ws3sec    = met['WSMAX']*1.94384  # convert from m/s to kt
    wd        = met['WD']
    RH        = met['RH']


    # set up fontsizes 
    plt.rc('xtick',labelsize=11)
    plt.rc('ytick',labelsize=13)
    
    # set up figure
    fig, axes = plt.subplots(5, figsize = (12,18))
    top_anchor = 1.15 # changes title/legend relation to their plot
    lfs = 13 # legend fontsize
    ylfs = 14 # ylable fontsize
    ec = '0.5' # color of legend outline. default is 0.8


    # plot wind direction and speed

    axes[0].set_ylim(0,max(ws3sec)+2)
    axes[0].plot(dates,ws,color='#57858c',alpha = 0.5, linewidth=1)
    ln1 = axes[0].fill_between(dates,ws,plt.ylim()[0],color='#57858c',alpha = 0.5,label='Wind Speed')
    ln2 = axes[0].scatter(dates, ws3sec, marker='s', color = '#655978', 
                          s=1.2, alpha=0.5, label='3-sec Gust', zorder=0)
    axes[0].set_ylabel('Wind Speed\n(knots)', multialignment='center',fontsize = ylfs)

    ax0 = axes[0].twinx()
    ln3 = ax0.plot(dates,wd,'.k', alpha = 0.8, markersize=3, label='Wind Direction')
    ax0.set_ylim(0,360)
    ax0.set_ylabel('Wind Direction\n(degrees)', multialignment='center',fontsize=ylfs)
    ax0.set_yticks(np.arange(45,405,90))
    ax0.set_yticklabels(['NE','SE','SW','NW'])

    lns = [ln1,ln2,ln3[0]] 
    labs = [l.get_label() for l in lns]
    ax0.legend(lns, labs, loc='upper left', markerscale=5,
               bbox_to_anchor=(0.005,top_anchor),ncol=3,prop={'size':lfs},framealpha=1,edgecolor=ec)


    # plot temp and dewpoint and RH

    ymax,ymin = (max(tempf)+3,min(dewpf)-3)
    axes[1].set_ylim(ymin,ymax)
    axes[1].set_ylabel('Temperature\n(F)', multialignment='center',fontsize = ylfs)
    ln1=axes[1].fill_between(dates,tempf,dewpf,color = '#a83939',alpha = 0.6, label = 'Temperature') 
    ln2=axes[1].fill_between(dates,dewpf,plt.ylim()[0],color = '#207849', alpha = 0.7,label = 'Dewpoint')

    ax1 = axes[1].twinx()
    ln3 = ax1.plot(dates,RH,color='#2b0000', linewidth=3,alpha=0.5, label='Relative Humidity')
    ax1.set_ylim(0,103)
    ax1.set_ylabel('Relative Humidity\n(%)',multialignment='center',fontsize = ylfs)
    lns = [ln1,ln2, ln3[0]] 
    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs, loc='upper left', 
               bbox_to_anchor=(0.005,top_anchor),ncol=3,prop={'size':lfs},framealpha=1,edgecolor=ec)


    # MSLP

    axes[2].plot(dates,mslp,color = '#E37609',linewidth=3,alpha = 0.8,label='Mean Sea Level Pressure')
    ymax,ymin = np.round(np.max(mslp))+1,np.round(np.min(mslp))-1
    axes[2].set_ylim(ymin,ymax)
    axes[2].set_ylabel('Mean Sea Level\nPressure (hPa)', multialignment='center',fontsize = ylfs)
    axes[2].legend(loc='upper left', bbox_to_anchor=(0.005,top_anchor),ncol=1,prop={'size':lfs},framealpha=1,
                  edgecolor=ec)


    # plot theta_v and theta_e

    ymax,ymin = (max(thetav)+2,min(thetav)-2)
    axes[3].set_ylim(ymin,ymax)
    ln1=axes[3].plot(dates,thetav, color = '#084b69',alpha = 0.7, linewidth=4,
                     label = r"${\theta_v}$")
    axes[3].set_ylabel(r'${\theta_v}$ (K)', multialignment='center',fontsize = ylfs)

    ax3 = axes[3].twinx()
    ymax,ymin = (max(thetae)+2,min(thetae)-2)
    ax3.set_ylim(ymin,ymax)
    ln2=ax3.plot(dates,thetae, color = '#508a50',alpha = 0.6,linewidth=3,
                 label = r"${\theta_e}$")
    ax3.set_ylabel(r'${\theta_e}$ (K)', multialignment='center',fontsize = ylfs)
    lns= [ln1[0],ln2[0]]
    labs=[l.get_label() for l in lns]
    ax3.legend(lns,labs,loc='upper left', bbox_to_anchor=(0.005,top_anchor),
               ncol=3,prop={'size':lfs},framealpha=1,edgecolor=ec)

    # theta_v gradient

    gradient = thetav[1:] - thetav[0:-1]

    axes[4].plot(dates[1:], gradient, color='#204e63',alpha=0.9,label=r"1-min ${\Delta\theta_v}$")
    axes[4].set_ylim(gradient.min()-1, gradient.max()+1)
    axes[4].set_ylabel(r'${\Delta\theta_v}$ (K/min)', multialignment='center',fontsize = ylfs)
    axes[4].legend(loc='upper left', bbox_to_anchor=(0.005,top_anchor),ncol=3,prop={'size':lfs},framealpha=1,
                  edgecolor=ec)
    axes[4].text(0.02, 0.1, 'Largest Decrease: {} K/min at {}'.format(np.round(min(gradient),2),
                dates[np.argmin(gradient)+1].strftime('%d/%H:%M UTC')), fontsize=14, color ='#204e63',
                 transform=axes[4].transAxes)

    #28264d the pretty purple color
    
    # number of hours rounded up
    length = np.ceil((dates[-1]-dates[0]).seconds/(60*60) + (dates[-1]-dates[0]).days*24)

    for ax in axes.flatten():
        
        # assign proper date formatting based on number of hours plotted
        if length > 8:
            ax.xaxis.set_major_formatter(mpl.dates.DateFormatter('%d/%H UTC'))
        else:
            ax.xaxis.set_major_formatter(mpl.dates.DateFormatter('%d/%H%M UTC'))
            if length==5: # things look bad at this length for some reason
                ax.tick_params(labelsize=10)
                
            
        # set grids and things
        ax.grid(color='k', linestyle=':', linewidth=0.3, alpha=0.7)
        ax.set_xlim(dates.min(),dates.max()) # removes white space buffer

    # title              
    plt.suptitle('{}-Hour Meteogram for Station {}'.format(int(length), probe_id), fontsize = 18, 
                  x=0.135, y = .94, weight='bold', ha='left')
    axes[0].set_title('Latest Observation: {}'.format(dates[-1].strftime('%D %H:%M UTC')),
                      ha='left',x=0.015, y=1.21, fontsize=16)
    fig.subplots_adjust(hspace=0.35)


    # Plot the TTU logo
    im = plt.imread("TTU_Logo.tif")
    newax = fig.add_axes([0.04,0.87,0.07,0.07],anchor='NW',zorder=10)
    newax.imshow(im, alpha=0.4)
    newax.axis('off')

    # save figure 
    obtime = parse_currtime()[1] # note that this may be different from last time on plot
                                 # THATS OKAY... it will make data drops more apparent!
    plt.savefig('{0}{1}_meteogram_{2}.png'.format(savedir,probe_id,obtime),dpi=300,bbox_inches = 'tight')
    plt.close()
    
    
    