########################################################################################################################
# This script is to be used with realtime sticknet obs
# The purpose of this script is to create 24 hour meteograms 
#
# This script should be run every 5 MINUTES during a project
#
# Initially developed by A. Hill, Overhauled by J. McDonald 2021
########################################################################################################################

import matplotlib
matplotlib.use('Agg')

import SNmods as snmods
from info import probe_locs
import datetime as dt


endtime = dt.datetime(2017,4,30,19,9) #dt.datetime.utcnow()
starttime = endtime - dt.timedelta(hours=24)

probes = list(probe_locs.keys())[0:3]

# get available data and plot!
data_df = snmods.get_sticknet_data(starttime, endtime,probes=probes,
                                   dataset='subset', 
                                   plotmeteograms=True) # plot meteograms!