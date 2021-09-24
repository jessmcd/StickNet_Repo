# created by J. McDonald, 2019-2021


import numpy as np
import cartopy.crs as ccrs 
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import matplotlib.patches as patches
import matplotlib
import cmocean
import pyart;

def rgb(hexcol):
    return matplotlib.colors.to_rgba(hexcol)

def make_cmap(colors, n_bin=50):
    cols = []
    for c in colors:
        cols.append(rgb(c))
    return LinearSegmentedColormap.from_list('newcmap', colors, N=n_bin)

# truncate an existing colormap
def min_col(cmap, minval=0.0, maxval=1.0, n=15):
    new_cmap = LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

def homeyer_rainbow_ramp():
    color_list = [pyart.graph.cm_colorblind.HomeyerRainbow(i) for i in np.arange(0, 1, .01)]
    ucmap = LinearSegmentedColormap.from_list('c', color_list, N=100)
    color_array = ucmap(np.arange(0, 1, .01))
    # change alpha values
    color_array[:12,-1] = [0,0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9, 1]
    # create a colormap object
    map_object = LinearSegmentedColormap.from_list(name='pyart_HomeyerRainbow_ramp',colors=color_array)
    # register this new colormap with matplotlib
    plt.register_cmap(cmap=map_object)

def drop_nans_wind(u,v,lats,lons):
    ''' drops nan values, for use in interpolation'''
    idx = np.where(np.isnan(u)==False)[0]
    if len(idx) < len(u):
        return u[idx],v[idx],lats[idx],lons[idx]
    else:
        return u,v,lats,lons
    
def drop_nans_var(var,lats,lons):
    ''' drops nan values, for use in interpolation
        also returns the N of missing values'''
    idx = np.where(np.isnan(var)==False)[0]
    if len(idx) < len(var):
        n_sn_missing = len(var)-len(idx)
        return var[idx],lats[idx],lons[idx], abs(n_sn_missing)
    else:
        return var,lats,lons, 0

def no_data_check(data_df):
    ''' if there is no data, return empty data centered on Lubbock'''
    
    if data_df.empty: # no sticknets have been deployed! Center on Lubbock
    
        fake_data = {'Lat':[33.5779], 'Lon':[-101.8552], 'T':[np.nan], 
                     'RH':[np.nan], 'P':[np.nan], 'WS':[np.nan],
                  'WD':[np.nan], 'WSMAX':[np.nan], 'BATT':[np.nan], 
                     'date':[dt.datetime.utcnow()], 'Elevation':[np.nan]}

        data_df = pd.DataFrame(fake_data, index=['0NAN'])
    return data_df


def plot_logo(fig, ax, height=0.09, alpha=0.4, bumpx=False, bumpy=False):
    '''
    plots the ttu logo in the upper left corner.
    fig is the plt.figure object, and ax is the axis object.
    Height specifies the height AND width of the logo, 
    alpha is the transparency.
    
    If for some reason, the logo is plotting above the plot, set bumpy = True
    If the logo is plotting too far to the left, set bumpx = True
    
    **** IMPORTANT: you need to have plt.draw() ahead of using this function
    
    if its not plotting in the right place, set bump=True
    '''
    fix=0
    if bumpx:
        fix=(height/2)+0.005
        
    fixy=0
    if bumpy:
        fixy = (height)+0.01
       
    
    im = plt.imread("TTU_Logo.tif")
    # get corner
    p1 = ax.get_position()
    # put logo in corner
    ax2 = fig.add_axes([p1.x0+fix+0.01, p1.y1-height-0.01-fixy,height,height],anchor='NW',zorder=20) # only the last two nums matter
    ax2.imshow(im, alpha=alpha)
    ax2.axis('off')

    
    
    
def scale_bar(fig, ax, length=None, location=(0.5, 0.05),N=4,
              height=1, fontsize=10, text_pad=1, edge_width=0.5,
              edge_color='k', text_color='k', color1='w', color2='0.7', 
              unit_label='km', top=True, ticks=True, tick_len=.8):
    """
    Adds a scale bar to cartopy plots. Only use even numbers in km for length!
    ***you MUST use plt.draw() before calling this function, or placement will be wrong***
   
    INPUTS:
        fig         - the matplotlib figure, for sizing the scale bar appropriately
        ax          - the axis object to draw the scale bar on.
        length      - the length of the scale bar in km. Only use even numbers!
        location    - (x,y) the center of the scalebar in *figure coords*
                      (i.e., (0.5, 0.5) is the middle of the plot). Will be different
                      from ax coords if ax is projected
        N           - The number of alternating bars. Only works for even numbers!
        height      - the height of the scale bar (default=1)
        fontsize    - the size of the text (default=10)
        text_pad    - controls the spacing between the text and the scale bar (default=1)
        edge_width  - the width of the scale bar edge
        edge_color  - the color of the edge (default='k')
        text_color  - the color of the text (default='k')
        color1      - the first alternating color of the scale (default='w')
        color2      - the second alternating color of the scale (default='0.7')
        unit_label  - adds unit to the scale. Set to '' if you don't want it (default='km')
        top         - if True, numbers are located above the scale, otherwise they are below (default=True)
        ticks       - if True, adds ticks in line with the numbers (default=True)
        tick_len    - the length of the ticks, valid only if ticks=True (default=.8)

    Dependencies: requires import matplotlib.patches as patches and cartopy.crs as ccrs 
    
    Edited by J. McDonald 2021. Based on GitHub - https://github.com/SciTools/cartopy/issues/490
    """

    #Get the limits of the axis in lat long
    llx0, llx1, lly0, lly1 = ax.get_extent(ccrs.PlateCarree())
    #Make Transverse Mercator proj (tmc) centered on the middle of the map,
    sbllx = (llx1 + llx0) / 2
    sblly = (lly1 - lly0) / 2  
    tmc = ccrs.TransverseMercator(sbllx, sblly)
    
    #Get the extent of the plotted area in coordinates in meters
    x0, x1, y0, y1 = ax.get_extent(tmc)
 
    # if no length added, make scale ~1/5th of the plot width
    if length == None: 
        length= (int((x1 - x0)*.2/1000)- (int((x1 - x0)*.2/1000)%10))*1000 # meters
    else:
        length = length*1000 # convert to meters

    # get the width in ax-relative coords, then convert to fig-relative coords
    ax_width = length/(x1 - x0)
    fig_width = ax_width*(ax.get_position().x1-ax.get_position().x0)

    # add new axes of the correct width (the same as the length desired)
    ax2 = fig.add_axes([location[0]-fig_width/2, location[1] ,fig_width,.25], anchor='C', zorder=0)
    ax2.set_xlim(0,1)
    ax2.axis('off')

    # horizontal center of scale within ax2, keeps numbers from touching bottom when
    # using ax.get_position().y0 for the y location
    y = .09
    fc = [color1, color2]

    # add in the alternating colors, (x,y), width, height
    for i in np.arange(0, N, 1):
        r = patches.Rectangle((i/N  , y-.035*height/2), 1/N, .035*height, lw=edge_width,ec=edge_color, 
                              fc=fc[i%2], transform=ax2.transAxes, clip_on=False)
        ax2.add_patch(r)


    buffer = .035*height / 2

    # if ticks, let tick height > scale width
    tick_h = 1
    if not ticks: tick_h = 0; tick_len=0

    # put numbers on top of scale
    # the units have a "tick_len" adjustment, so that changing the text_pad for longer
    # tick_lens won't result in a unit spacing that's stupidly large
    if top: 

        ticks_sign = 1
        ax2.text(0.5, y-buffer-text_pad/100-.01+tick_len/100,'km', transform=ax2.transAxes,
                            ha ='center', va='top', color=text_color, fontsize=fontsize)

        for i in np.arange(0, N-1, 1):
            ax2.text(i/2, y+buffer+text_pad/100 , str(int((i/2000)*length)), transform=ax2.transAxes,
                            ha ='center', va='bottom', color=text_color, fontsize=fontsize)
            # ticks
            ax2.plot([i/2, i/2],[y, y+(.035*height*tick_len)*tick_h], lw=edge_width, 
                     c=edge_color,transform=ax2.transAxes, clip_on=False)

    # put numbers on bottom of scale    
    if not top: 
        ticks_sign = -1
        # km
        ax2.text(0.5, y+buffer+text_pad/100-tick_len/100,'km', transform=ax2.transAxes,
                            ha ='center', va='bottom', color=text_color, fontsize=fontsize)

        for i in np.arange(0, N-1, 1):
            ax2.text(i/2, y-buffer-text_pad/100-0.01 , str(int((i/2000)*length)), transform=ax2.transAxes,
                            ha ='center', va='top', color=text_color, fontsize=fontsize)
             # ticks
            ax2.plot([i/2, i/2],[y, y-(.035*height*tick_len)*tick_h], lw=edge_width, 
                     c=edge_color,transform=ax2.transAxes, clip_on=False)
            