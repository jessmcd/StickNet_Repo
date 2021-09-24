#####################################################
# Functions used in scripts for the VORTEX SE project
#####################################################

### modified by Jessica M. McDonald 2017/2018

import numpy as np
import datetime as dt
from boto.s3.connection import S3Connection
import tempfile
import pyart
import pandas as pd

# find nearest radar so you don't have to hardcode it into the scripts
def get_radar_ID(lat, lon):
    
    ''' 
    Lat and lon are the location of the center of your domain (average sticknet loc?). 
    Returns the 4-letter identifier (e.g., 'KHTX' for hytop) needed to get radar data using AWS, 
    and the name of the city that the radar is located in for a sanity check.
    
    Created 2021 by J. McDonald
    '''
    
    RadarLocs = pd.read_csv('NEXRAD_Radar_Locations.csv')

    dx = (lon-RadarLocs['LON'])*40000*np.cos((lat+RadarLocs['LAT'])*np.pi/360)/360
    dy = (lat-RadarLocs['LAT'])*40000/360
    ds = np.sqrt(dx**2+dy**2) # km
    
    return  RadarLocs.iloc[np.argmin(ds)]['ID'], RadarLocs.iloc[np.argmin(ds)]['CITY']

# Helper function for the radar search
def _nearestDate(dates, pivot):
    return min(dates, key=lambda x: abs(x - pivot))

# Fetch radar from AWS using boto S3connection
def get_radar_from_aws(site, datetime_t, datetime_te):
    """
    Get the closest volume of NEXRAD data to a particular datetime.
    Parameters
    ----------
    site : string
        four letter radar designation
    datetime_t : datetime
        desired start date time
    datetime_te : datetime
        desired end date time
    Returns
    -------
    radar : Py-ART Radar Object
        Radar closest to the queried datetime
    """

    # First create the query string for the bucket knowing
    # how NOAA and AWS store the data
    my_pref = datetime_t.strftime('%Y/%m/%d/') + site

    # Connect to the bucket
    conn = S3Connection(anon = True)
    bucket = conn.get_bucket('noaa-nexrad-level2')

    # Get a list of files
    bucket_list = list(bucket.list(prefix = my_pref))

    # Create a list of keys and datetimes to allow easy searching
    keys = []
    datetimes = []

    # Populate the list
    for i in range(len(bucket_list)):
        this_str = str(bucket_list[i].key)
        # Ensure that correct datetime data is grabbed from 
        # the two possible filename formats
        if 'gz' in this_str:
            endme = this_str[-22:-4] # Grabs datetime data from list of radar files
            #endme = this_str[20:38]# if getting mismatched format error, try this instead
            fmt = '%Y%m%d_%H%M%S_V0'
            dts = dt.datetime.strptime(endme, fmt)
            datetimes.append(dts)
            keys.append(bucket_list[i])

        if this_str[-3::] == 'V06':
            endme = this_str[-19::] # Grabs datetime data from list of radar files
            fmt = '%Y%m%d_%H%M%S_V06'
            dts = dt.datetime.strptime(endme, fmt)
            datetimes.append(dts)
            keys.append(bucket_list[i])

    # Find the closest available radar to the start and end datetimes
    closest_datetime_b = _nearestDate(datetimes, datetime_t)
    closest_datetime_e = _nearestDate(datetimes, datetime_te)
    
    # Locate the indices of start/end times within datetimes list
    index_b = datetimes.index(closest_datetime_b)
    index_e = datetimes.index(closest_datetime_e)
    print(index_b, index_e)

    # Grab all radar files within the start and end times
    radar_namelist = keys[index_b:index_e+1]
    radar_list=[]
    for i in range(np.shape(radar_namelist)[0]):
        # Creating a temp file provides temporary storage area for radar data
        localfile = tempfile.NamedTemporaryFile() 
        # Put radar filename within the namelist
        radar_namelist[i].get_contents_to_filename(localfile.name)
        # Append Pyart radar data to radar list
        radar_list.append(pyart.io.read(localfile.name))
        
    return radar_namelist,radar_list

### 8/15/2018 - Added by Jessica M. McDonald
### Read in ARMOR data 

import glob 
import re

# sort files so that they're in order... because they aren't for some reason

def sorted_nicely(list):
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(list, key=alphanum_key)

def get_ARMOR_data(filepath, starttime, endtime):
    
    ''' give filepath to ARMOR data, starttime and endtime should be datetime objects.
    Used to Return Pyart.Radar object, but ran into memory issue. now only returns the actual
    file paths. Use pyart.io.uf.read_uf(file) to open in your own script, one at a time'''
    
    bucket_list= glob.glob(filepath+'*.uf')
    bucket_list = sorted_nicely(bucket_list)

    datetimes = []
    for i in range(len(bucket_list)):
        date_str = bucket_list[i][-21:-7]
        fmt = '%Y%m%d%H%M%S'
        dts = dt.datetime.strptime(date_str, fmt)
        datetimes.append(dts)
        
     #Find the closest available radar to the start and end datetimes
    closest_datetime_b = _nearestDate(datetimes, starttime)
    closest_datetime_e = _nearestDate(datetimes, endtime)

    # Locate the indices of start/end times within datetimes list
    index_b = datetimes.index(closest_datetime_b)
    index_e = datetimes.index(closest_datetime_e)

    print(index_b, index_e)
    
    radar_namelist = bucket_list[index_b:index_e+1]
    #radar_list=[]
    #for i in range(np.shape(radar_namelist)[0]):
        # Append Pyart radar data to radar list
        #radar_list.append(pyart.io.uf.read_uf('{0}'.format(radar_namelist[0])))
        
    return radar_namelist
