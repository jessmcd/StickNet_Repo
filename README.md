# "sticknet" : Repository containing code to parse and analyze StickNet data


### Important information:
  all scripts rely on the info.py script. This script contains hardcoded variables, such as the StickNet locations, 
  the time zone, the directory that holds the realtime data, and the directory that holds all the created figures. 


### Dependencies:

  - boto3
  - pandas
  - datetime
  - py-art
  - metpy
  - scipy 
  - cartopy
  - netCDF4
  - cmocean
  - StickNet object class
  
  Everything else is basic python (numpy, matplotlib, math, etc. ) 
  


### Other Info:
  All the .ipynb scripts are for testing purposes. The .py scripts are considered the final production code. Do not do any major     changes to the .py scripts without testing in the .ipynb scripts first!
