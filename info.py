######################
# This script will contain all hardcoded variables 
#
# Fill this out as is needed. Contains:
#
# CURRENT lat, lon and elevation of probes in field, "probe_locs"
#
# TIME ZONE that StickNets are in, "TZ"
#
# FILE PATH to realtime files, "filedir"
#
# FILE PATH to directory to save plots to, "savedir"
#
#####################



# thoughts on handling these during the field...
# enter probe IDs that you know you're going to deploy that day, enter 0s for everything
# this will allow the information to start populating in the obs table on the website, but plots won't work.
# OR, if you know where they're going to be put (like VSE), enter in estimated lat/lons, and update when you can

              # ID.    LAT.     LON.      Elev.
probe_locs = {'0102A':[34.29590,-87.58710,282],
              '0103A':[34.85510,-86.00170,513],
              '0104A':[35.33910,-87.03240,239],
              '0105A':[34.62240,-86.08000,182],
              '0106A':[34.21490,-87.16190,245],
              '0107A':[34.19439,-86.79825,231],
              '0108A':[35.32360,-86.63470,230],
              '0109A':[35.30450,-87.51850,315],
              '0110A':[34.16320,-86.33410,280],
              '0111A':[34.55040,-86.55820,174],
              '0112A':[34.72540,-87.46250,174],
              '0213A':[34.90160,-86.53860,239],
              '0214A':[34.93090,-86.97640,256],
              '0215A':[35.03800,-87.47190,237],
              '0216A':[35.37190,-86.09980,305],
              '0217A':[34.61760,-87.10140,192]}

TZ = "US/Central" # "US/Mountain", "US/Eastern, "US/Central"

filedir = '/Users/jessmcd/Documents/GitHub_repos/StickNetRepo/RT_tests_data/'

savedir = '/Users/jessmcd/Documents/GitHub_repos/StickNetRepo/RT_tests/'