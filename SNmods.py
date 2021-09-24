import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import datetime as dt
import glob
#import SNobject as sn
#import probe_info
from info import probe_locs, savedir, filedir
from functions import convert_wind,parse_currtime,calc_dewpoint,plot_meteogram,get_winddir_string,C_to_F,calc_mslp

# Hard Coded Stesonet locations (year can change)
#probe_locs = probe_info.probe_locs_2017

# rows for StickNet data units for Vars are C, %, hPa, m/s, m/s, and degrees
col_names = ['probe', 'Lat', 'Lon', 'YYMMDD', 'HHMM','T', 'RH', 'P', 'WS','WSMAX', 'WD', 'BATT']

class SNFile(object):
    def __init__(self, filename):
        """ initialize a SN object based on a filename given """
        self.year = int(filename[-17:-13])
        self.month = int(filename[-13:-11])
        self.day = int(filename[-11:-9])
        self.hour = int(filename[-8:-6])
        self.minute = int(filename[-6:-4])
        self.datetime = dt.datetime(self.year,self.month,self.day,self.hour,self.minute)
        self.filename = filename
        self.probe = filename[-23:-18]

    def read_realtime(self):
        """ read the one line of the RT sticknet text file and split apart.
       """
        f = open(self.filename)
        line = f.readline().split(',')
        return line

#     def read_daq(self,identifier='stesonet'):
#         """ read the DAQ box output text file and parse into columns. Identify with 
#         call input whether the filetype is from StesoNet probe or finescale, as this determines 
#         number of columns to read """
#         if identifier == 'stesonet':
#             sn = pd.read_csv(self.filename, names=['time','T','RH','P','windsp','winddir','batt'], 
#                 dtype={'time':np.str},header=1,parse_dates=[0],date_parser=self.parse,error_bad_lines=False)
#         elif identifier == 'finescale':
#         # for rapid probes, that don't have the extra data point at the end of the line (i.e. battery voltage)
#             sn = pd.read_csv(self.filename, names=['time','T','RH','P','windsp','winddir'],
#                          dtype={'time':np.str},header=1,parse_dates=[0],date_parser=self.parse,error_bad_lines=False)
         
#         sn['ID'] = pd.Series(self.probe,index=sn.index)
#         sn.index = sn['time']
#         return sn

#     def read_daq_qced(self,identifier):
#         """ read the QC'd DAQ box output text file and parse into columns. Identify with 
#         call input whether the filetype is from StesoNet probe or finescale, as this determines 
#         number of columns to read """        
#         if identifier == 'stesonet':
#             sn = pd.read_csv(self.filename, names=['time','T','RH','P','windsp','winddir','batt','tflag','wflag'], 
#                 dtype={'time':np.str},header=1,parse_dates=[0],date_parser=self.parse,error_bad_lines=False)
#         elif identifier == 'finescale':
#         # for rapid probes, that don't have the extra data point at the end of the line (i.e. battery voltage)
#             sn = pd.read_csv(self.filename, names=['time','T','RH','P','windsp','winddir','tflag','wflag'],
#                          dtype={'time':np.str},header=1,parse_dates=[0],date_parser=self.parse,error_bad_lines=False)
         
#         sn['ID'] = pd.Series(self.probe,index=sn.index)
#         sn.index = sn['time']
#         return sn

def get_files(filedir, probe_id, starttime, endtime):
    ''' this seems like a janky way to do it, but is actually 3x faster than
        making a loop of files.'''

    day_str = np.array([(starttime+dt.timedelta(days=i)).strftime("%Y%m%d") 
                     for i in range((endtime-starttime).days+1)])


    # find all SN files with days between starttime and endtime
    files = [glob.glob(filedir+probe_id+'_'+string+'_*.txt') for string in day_str]
    files = [f for subf in files for f in subf] # flatten list in case of multiple days
    
    if files:
         # sort files by date, then find nearest indices for all the dates, and loop over that
        sorted_files = sorted(files,key=lambda f: dt.datetime.strptime(f[-17:-4],'%Y%m%d_%H%M'))
        fdates = [dt.datetime.strptime(f[-17:-4],'%Y%m%d_%H%M') for f in sorted_files] 
        _, idx1 = min((abs(val-starttime), idx) for (idx, val) in enumerate(fdates))
        _, idx2 = min((abs(val-endtime),   idx) for (idx, val) in enumerate(fdates))
        files= sorted_files[idx1:idx2+1]
        
    return files


def format_ID(n):
    if n > 12:
        probe_id = "02{0}A".format("%02d"%n)
    else:
        probe_id = "01{0}A".format("%02d"%n)
        
    return probe_id
    
def write_to_html(html,filedir, probe_id, endtime, d):
    ''' all inputs are created in other function
    Writes to html file for website Observation Table.
    Writes out City name, 4-letter city identifer, probe id, current T, Td,
    RH, WS, WD (string format), wind gust, 24 hr max T, 24 hour min T, 24 hour Max wind gust,
    index of max wind gust (idk why this is needed), batter, MSLP.
    Units for temps are F, wind is kt, pressure is hPa, and battery is Volts. '''
    
    # 24 hour data 
    s = endtime - dt.timedelta(hours=24)
    files = get_files(filedir, probe_id, s, endtime)

    data, dates = [],[] # to put data in for faster looping
    for f in files:
        sn_data = SNFile(f)
        data.append(sn_data.read_realtime())
        dates.append(sn_data.datetime)

    met24 = pd.DataFrame(data, index=dates, columns=col_names)
    met24 = met24[s:endtime] # just in case 
    met24.drop(columns=['probe'], inplace=True)
    met24 = met24.astype(float)


    Tmax = np.round(C_to_F(np.amax(met24['T'].values)),1)     # F
    Tmin = np.round(C_to_F(np.amin(met24['T'].values)),1)      # F 
    WSgustmax = np.round(np.amax(met24['WSMAX'].values)* 1.94384, 1)# kts
    wsindx = np.argmax(met24['WSMAX'].values) # not sure why I need this...

    # current vals
    T = C_to_F(d['T']).values[0]
    Td = C_to_F(calc_dewpoint(d['T'].values, d['RH'].values))[0]
    RH = d['RH'].values[0]
    WS = np.round(d['WS']* 1.94384, 1).values[0]
    WSgust = np.round(d['WSMAX']* 1.94384, 1).values[0]
    WD = get_winddir_string(d['WD'].values)
    batt = d['BATT'].values[0]
    MSLP = np.round(calc_mslp(d['T'].values,d['P'].values,d['Elevation'].values),1)[0]
    city = probe_locs[probe_id][4]
    id_name = probe_locs[probe_id][3]
    time = sn_data.datetime

    text1 = f'{city},{id_name},{probe_id},{time},{T},{Td},{RH},{WS},{WD},{WSgust:4.1f},'
    text2 = f'{Tmax},{Tmin},{WSgustmax:4.1f},{wsindx},{batt},{MSLP}'
    text = text1+text2

    print(text)
    #html.write(text)

    
def get_sticknet_data(starttime, endtime, probes=[], dataset="subset",
                      plotmeteograms=False, returndata=False, html=None,
                      filedir=filedir, savedir=savedir):
    
    '''
    This is a somewhat complicated function that performs multiple tasks depending on the inputs.
    It is designed to be helpful for the website, but should also be used for post processing and analysis. 
    
    INPUTS:
   
    starttime - datetime object of start time of data you want. Only important when dataset="subset"
    endtime   - datetime object of end time of data you want. 
    probes    - int number associated with the SNs you want (1 = "0101A", 24 = "0224A") OR
                the full id string. Must be in list format
                
    dataset   - can be "subset" or "latest"
    
        if "subset":  Get a range of data. If starttime = endtime, only one minute of data is returned.
                      Used for plotting 1 or more meteograms, or analysis on a single SN
    
          plotmeteograms - default False. Set to True if you want meteograms to be plotted in the data 
                           range you selected. Set returndata to False, or you'll waste some computing power.
          returndata     - default False. set to True if you want the pandas dataframe of the data 
                           returned. If multiple probes are chosen in the "probes" var, only the 
                           last is returned. Latitude, longitude, and elevation are set as attributes
                           for the dataframe.(i.e, df.attrs['latitude'] will give you the latitude)
                           
        if "latest":   Get the a single line of data closest to "endtime". If starttime != endtime, 
                       functionality is the same, but it will be slower.
                       Used for making station plots, oban plots, and for creating files useful for 
                       the website
        
          html           - default None. Set html = an opened text file, and it will write out the
                           data needed for the Obs Table on the website to that file.
                           
    filepaths: These should be updated in info.py, otherwise set them to what you want (useful for analysis)
    filedir   - location of the data in real time form (each file only contains 1 min of data for 1 SN)
    savedir   - where to save meteograms. 
                           
    OUTPUTS:
        Returns a Pandas dataframe that will be full ("latest" OR ("subset"&returndata)) 
        or empty ("subset").
        
        All returned dataframes will include the attribute "units", which outputs a reminder
        of the units of each variable.
        
        "latest" returns a dataframe with probe_id as index, with dates as a column.
        "subset" & returndata=True returns a dataframe with dates as the index. Probe ID is 
                 contained in the dataframe attrs
    
    '''

    # create empty dataframe
    met = pd.DataFrame(columns=col_names)
 
    # loop over all probes
    for i in probes:
        
        # check for type of probe input
        if type(i) == str:
            probe_id = i
        else:
            probe_id = format_ID(i)

        files = get_files(filedir, probe_id, starttime, endtime)

        if files: # not empty 
            
            print(probe_id)

            ### get a subset of data
            if dataset == 'subset':

                    data, dates = [],[] # to put data in for faster looping

                    for f in files: # append data to lists, LOTS faster than appending with pandas
                        sn_data = SNFile(f)
                        data.append(sn_data.read_realtime())
                        dates.append(sn_data.datetime) 

                    met = pd.DataFrame(data, index=dates, columns=col_names)

                    # in case the "nearest time" is really not that near, use Pandas to make sure
                    # times are within what you requested. Clean up dataframe
                    met.drop(columns=['Lat', 'Lon', 'YYMMDD', 'HHMM','probe'], inplace = True)
                    met.sort_index(ascending=True, inplace=True)
                    met = met[starttime:endtime]
                    met = met.astype(float)

                    # set attributes                    
                    met.attrs = {'probe':probe_id, 
                                 'latitude':probe_locs[probe_id][0],
                                 'longitude':probe_locs[probe_id][1], 
                                 'elevation':int(probe_locs[probe_id][2]),
                                 'units':'T(C), RH(%), P(hPa), WS(m/s), WSMAX(m/s), WD(deg), elevation(m)'}

                    if plotmeteograms:
                        plot_meteogram(met, probe_id, savedir) 

                    if returndata:
                        if i == probes[-1]: # only do for the last item, or it will kill the loop
                            return met.replace(-999.9,np.nan)

                    met = pd.DataFrame(columns=col_names) # start over for next SN

            ### get the time closest to end time
            ### will find the most recent file available IF endtime is utc.now()
            
            if dataset == 'latest': # get file closest to endtime

                    if len(files) > 1:
                        files = [files[-1]] # need the extra brackets so the next line works

                    sn_data = SNFile(files[0])
                    d = sn_data.read_realtime()
                    
                    # have date be a column, and probe_id be an index
                    d.insert(0, sn_data.datetime)
                    d = pd.DataFrame([d], index=[sn_data.probe], columns=['date']+col_names).replace(-999.9,np.nan)
                    d.drop(columns=['probe'], inplace=True)
                    d[col_names[1:]] = d[col_names[1:]].apply(pd.to_numeric) 
                    # add more info
                    d['Elevation'] = int(probe_locs[probe_id][2])
                    d['Lat']       = probe_locs[probe_id][0] # overwrite with potentially/ 
                    d['Lon']       = probe_locs[probe_id][1] # more accurate location from info


                    # data is less than 5 minutes old
                    if abs(sn_data.datetime - endtime) < dt.timedelta(minutes=5):
                        met = met.append(d)  # only have to append a few times, so it's okay to use pandas
                        # Write out info for Obs Table on website
                        if html:
                            write_to_html(html,filedir, probe_id, endtime, d) 

                    else: # no data
                        continue


    # return dataframe. Will be empty if subset % returndata=False
    if (dataset == 'latest'):
        if met.empty == False:
            met.drop(columns=['YYMMDD', 'HHMM', 'probe'], inplace=True)
            met.attrs = {'units': 'T(C), RH(%), P(hPa), WS(m/s), WSMAX(m/s), WD(deg), elevation(m)'}
            # change to floats
            to_nums = col_names[5:]+['Lat','Lon','Elevation']
            met[to_nums] = met[to_nums].apply(pd.to_numeric) 
            met.sort_index(ascending=True, inplace=True)
            
        return met.replace(-999.9,np.nan)