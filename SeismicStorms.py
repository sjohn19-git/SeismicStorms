#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 20:37:44 2024
@author: sebin john
email:sjohn19@alaska.edu
"""
#-------------------------------------------------------------
"""
This code generates a map which displays significant 
wave height data and microseismic data from different 
station specified by the ./station_list.txt file.
Significant wave height is obtained from NOAA wave watch
3 model.
Seismic data is obtained from the AEC database
This code always generates a map 24 hours behind the current time
This is because the availability of wave model is constrained
in the NOAA server
output location of the map is:./hourly maps/{datetime}.png
This code will generate a wave directory, PSDs directory
RESPs directory to store the data and metadata.
It reads files from the ./grids and ./cpts for plotting. 
If these folders are not available code will fail
Alaska_network_station_location.csv is another file this
code needs to plot the map.
This code employs median filtering for data processing
so medfilt.py should be located in the same directory of the 
code
"""

from obspy import UTCDateTime
import os
from os.path import join as jn
import pickle
import importlib.util
# loading updated version of wf2obspy
module_path = '/usr/local/aec/'+os.environ['ANT_VER']+'/data/python/wf2obspy.py'
spec = importlib.util.spec_from_file_location('wf2obspy', module_path)
wf2obspy = importlib.util.module_from_spec(spec)
spec.loader.exec_module(wf2obspy)
import pandas as pd
from datetime import datetime,timedelta
import scipy.signal as sig
from obspy.signal.invsim import evalresp
from obspy.clients.iris import Client
import numpy as np
import sqlite3
import sys
# Get the directory of the script
if getattr(sys, 'frozen', False):  # if the script is frozen (e.g., by PyInstaller)
    script_dir = os.path.dirname(sys.executable)
elif '__file__' in globals():  # if __file__ is defined (i.e., not in interactive mode)
    script_dir = os.path.dirname(os.path.abspath(__file__))
else:  # if __file__ is not defined (e.g., in interactive mode or during import)
    script_dir = os.getcwd()
root_path=script_dir
os.chdir(root_path)
from medfilt import medfilt #import a custom code to median filter
import wget
import glob
import pytz
import pygrib
import xarray as xr
import pygmt
from scipy.interpolate import griddata


def get_folder_size(folder_path):
    '''This function returns the size of a folder
        
        :param folder_path:(str) Path of the folder
    '''
    
    #function to calculate the size of a folder
    total_size = 0

    for dirpath, dirnames, filenames in os.walk(folder_path):
        for filename in filenames:
            file_path = os.path.join(dirpath, filename)
            total_size += os.path.getsize(file_path)

    # Convert bytes to megabytes
    total_size_mb = total_size / (1024 * 1024)
    return total_size_mb

root_path=os.getcwd()
class SeismicStorms: 
    '''This class deals with downloading and 
    processing of seismic data.'''
    
    def __init__(self,st):
        self.starttime=UTCDateTime(st)
        self.endtime=UTCDateTime(st)+3600
        self.compute_psds()
        
    def wel(self,stream,nfft,windlap):
        ''' This function returns the Power spectra of 
        seismic data along with period bins.
        For more info refer to scipy.signal.welch method
        This function reads response files from ./RESPs
        directory. If not available attempts download
        from IRIS
        
        :param stream: obspy stream 
        :param nfft: (int) no of datapoints per segment
        :param windlap: (int) overlap window size
        '''
        
        #this dunction calculates power spectra using welch algorithm
        try:    
            tr=stream[0].copy()
        except:
            power=np.ones((1,8192))*np.nan
            peri=np.nan
            return power,peri
        if not os.path.exists(jn(root_path,"RESPS")):
            os.makedirs(jn(root_path,"./RESPS"))
        client=Client()
        if not os.path.exists(jn(root_path,"RESPS/"+tr.stats.station+"RESP")): 
            #check for response files in local directory
            # if not exists get response from iris
            resp=client.resp("AK", tr.stats.station,'*','BHZ',self.starttime,self.endtime,filename=jn(root_path,"RESPS/"+tr.stats.station+"RESP"))
        '''Following line gets the instrument response from a SEED RESP-file.
        For more info refer to obspy.signal.invsim.evalresp'''
        resp1 = evalresp(t_samp = tr.stats.delta, nfft=nfft, 
                  filename=jn(root_path,"RESPS/"+tr.stats.station+"RESP"), date = tr.stats.starttime,
                  station=tr.stats.station, channel=tr.stats.channel,
                  locid=tr.stats.location, network=tr.stats.network,
                  units="ACC")
        data=tr.data
        power=[]
        f,Pxx=sig.welch(data,fs=50,nperseg=nfft,noverlap=nfft*windlap)
        dbPxx=10.*np.log10(Pxx[1:]/(np.abs(resp1[1:])**2))
        power.append(dbPxx)
        power=np.array(power)
        return power,1/f[1:]
    
    def compute_psds(self):
        '''
        Returns the mean PSD in 5-10s band for all stations
        Downloads seismic waveform data from the AEC database 
        for stations specified in ./station_list.txt file 
        After downloading data, this function passes data
        to wel() to calculate psds
        '''
        
        stations=[]
        with open("./station_list.txt", 'r') as f:  # loading stations to plot temperature map
            content = f.read()
            sta = content.replace('\n', '').split(',')
            stations.extend(sta)
        self.stations=stations
        meansm=[]
        for sta in stations:
            try:# handling day boundary exception in wf2obspy
                if self.endtime.hour==0:
                    self.endtime=self.endtime-1
                    stream = wf2obspy.get_waveforms(["AK","AV"],sta, "*", "BHZ", self.starttime, self.endtime)
                    stream[0].data=np.ma.concatenate([stream[0].data,np.ma.masked_array(np.zeros(50))])
                else:
                    stream = wf2obspy.get_waveforms(["AK","AV"],sta, "*", "BHZ", self.starttime, self.endtime)
                nfft=2**14
                windlap=0.5
                psdWel,periWel=self.wel(stream,nfft,windlap)
                msm=np.mean(psdWel[0][31:65]) #extratcting 5-10s band energy from power spectra
                meansm.append(msm)
            except Exception as e:
                print(f"Error processing station {sta}: {e}")
                msm=np.nan
                meansm.append(msm)
        self.meansm=meansm
        return meansm 
    # Function to create the database and table if they don't exist
    
    def create_database(self): 
        '''Creates a ./PSDs/seismic_database.db file for saving mean seismic power 
           between 5-10 s.
           Checks if database file exists creates one if not exists.
           This database only has one table containing station names,
           time and mean seismic power (5-10s)
           '''
           
        if not os.path.exists(jn(root_path,"PSDs/")):
            os.makedirs(jn(root_path,"PSDs/"))
        connection = sqlite3.connect(jn(root_path,"PSDs/seismic_database.db"))
        cursor = connection.cursor()
    
        # Create the PSD table if it doesn't exist
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS PSD (
                station TEXT,
                time TEXT,  -- Store UTCDateTime as string
                power REAL,
                PRIMARY KEY (station, time)
            )
        ''')
    
        connection.commit()
        connection.close()
    
    def insert_psd_data(self,station, time, power):
        ''' This function ensures that only 7 data points
        are stored in the database and also inserts the data
        in the right column
        This function modifies ./PSDs/seismic_database.db
        
        :param stations: (str) station name
        :param time: (UTCDateTime) time of data
        :param power: (float) mean seismic power in 5-10s
        '''
        
        connection = sqlite3.connect(jn(root_path,"PSDs/seismic_database.db"))
        cursor = connection.cursor()
    
        # Convert UTCDateTime to string format
        time_str = str(time)
    
        # Check if a record with the same station and time already exists
        cursor.execute('SELECT * FROM PSD WHERE station=? AND time=?', (station, time_str))
        existing_record = cursor.fetchone()
    
        if existing_record:
            # Update the existing record
            cursor.execute('UPDATE PSD SET power=? WHERE station=? AND time=?', (power, station, time_str))
        else:
            # Insert the new PSD data
            cursor.execute('INSERT INTO PSD (station, time, power) VALUES (?, ?, ?)', (station, time_str, power))
            cursor.execute('SELECT * FROM PSD WHERE station=?', (station,))
            records = cursor.fetchall()
            if len(records) >= 8:
                # Delete the oldest entry for the station
                oldest_time = min(records, key=lambda x: datetime.strptime(x[1], '%Y-%m-%dT%H:%M:%S.%fZ'))[1]
                cursor.execute('DELETE FROM PSD WHERE station=? AND time=?', (station, oldest_time))
        connection.commit()
        connection.close()
        
    def add_value(self,timecheck):
        ''' This function does mostly time-keeping
        Every time this function runs the last processed time
        is stored in a file called previous_runtime.pkl
        This function loads previous_runtime and creates it if
        it doesn't exist.
        This function also checks if the current processing time is 
        within two hours of the previous run.
        If not, reinitiates the .PSDs/seismic_database.db array
        If yes, updates the .PSDs/seismic_database.db
        
        :param timecheck: 1/0 (flag to turn off time check
                               helps in debugging)
        '''
        
        try:#load info about previous run time
            with (open(jn(root_path,"PSDs/previous_runtime.pkl"), "rb")) as openfile:
                pr_time=pickle.load(openfile)
            with open(jn(root_path,"PSDs/previous_runtime.pkl"), "wb") as fp:   #Pickling
                pickle.dump(self.starttime, fp)#write current run time
        except:# if previous run time info not available default it to a value
            with open(jn(root_path,"PSDs/previous_runtime.pkl"), "wb") as fp:   #Pickling
                pickle.dump(UTCDateTime(1960,1,1), fp)
            pr_time=UTCDateTime(1960,1,1)
        print("previous runtime is",pr_time)
        if (self.starttime-pr_time>3600*2 and timecheck==1): #checking if current run time is with in 2 hour of previous run time
        #if run time is not continues reinitiate the database
            if os.path.exists(jn(root_path,"PSDs/seismic_database.db")):
                os.remove(jn(root_path,"PSDs/seismic_database.db"))
            print("time not continuous reinitiating")
            self.create_database()
            for jj in range(len(self.stations)):
                t=self.starttime+(self.endtime-self.starttime)/2
                power=self.meansm[jj]
                sta=self.stations[jj]
                self.insert_psd_data(sta, t, power)
        else:
            if os.path.exists(jn(root_path,"PSDs/seismic_database.db")):
                print("connecting to database")
                for jj in range(len(self.stations)):
                    t=self.starttime+(self.endtime-self.starttime)/2
                    power=self.meansm[jj]
                    sta=self.stations[jj]
                    self.insert_psd_data(sta, t, power)
            else: #if database not found
                print("database not found reinitiating")
                self.create_database()
                for jj in range(len(self.stations)):
                    t=self.starttime+(self.endtime-self.starttime)/2
                    power=self.meansm[jj]
                    sta=self.stations[jj]
                    self.insert_psd_data(sta, t, power)
    
    def filter_seis(self,windw): # this function filters seismic data
        '''
        Returns filtered array.
        This function filters the seismic data in 
        database using median filtering in the 
        specified windw size
        uses .PSDs/seismic_database.db
        
        :param windw: (int) window size'''
        
        connection = sqlite3.connect(jn(root_path, "PSDs/seismic_database.db"))
        cursor = connection.cursor()
    
        # Select all data from the PSD table
        cursor.execute('SELECT * FROM PSD')
        records = cursor.fetchall()
        stations,times,powers=[],[],[]
        for record in records:
            stations.append(record[0])
            times.append(record[1])
            powers.append(record[2])
        data = {'station': stations, 'time': times, 'power': powers}
        df = pd.DataFrame(data)
        df['time'] = pd.to_datetime(df['time'])
        df_pivoted = df.pivot(index='station', columns='time', values='power')              
        ms_arr=np.array(df_pivoted[:])
        ms_arr_fl=np.empty((ms_arr.shape[0]))
        if windw>ms_arr.shape[1]: #generally windw is 7 hr but if the database is newly reinitated 
        #then we dont have enough data to filter in 7 hr window
        #following lines adjusts windw accordingly
            if ms_arr.shape[1]% 2 != 0:
                windw=ms_arr.shape[1]#length of database is 7 so this ensures the filtering window is 7 hr
            else:
                windw=ms_arr.shape[1]-1 #satisfying conditon that window length must be odd
        for i in range(ms_arr.shape[0]):
            if sum(~np.isnan(ms_arr[i,:]))>=4:# conditon atleast 4 datapoints should exist for result to consider valid
                ms_arr_fl[i]=medfilt(ms_arr[i,:],windw)[round(windw/2)]
            else:
                ms_arr_fl[i]=np.nan 
        self.ms_arr_fl=ms_arr_fl


class Wave: # this class downloads appropriate wave height data
    '''This class handles all processing associated with waves'''
    
    def __init__(self,st):
        self.starttime=UTCDateTime(st)
        self.endtime=UTCDateTime(st)+3600
    
    def str_maker(self):
        '''Returns a list of parameters that defines the 
        significant wave height file to be 
        downloaded in the format [year, month, day, hour, forecast duration]
        This function identifies according to time
        of the day which significant wave height forecast should be downloaded
        '''
        
        t=self.starttime
        h=t.hour
        yr=str((t).year)
        mnt=str((t).month).zfill(2)
        day=str((t).day).zfill(2)
        # The following section handles which wave height data to be downloaded
        #Forecasts are available for wave height data
        # we use a combination of forecast and nowcast to acquire the most relevant wave height data
        if (h<=1):
            yr=str((t-24*3600).year)
            mnt=str((t-24*3600).month).zfill(2)
            day=str((t-24*3600).day).zfill(2)
            hr="18"
            fr="06"
        elif (h>1 and h<=3):
            yr=str((t-24*3600).year)
            mnt=str((t-24*3600).month).zfill(2)
            day=str((t-24*3600).day).zfill(2)
            hr="18"
            fr="09"
        elif (h==4):
            hr="00"
            fr="03"
        elif (h>4 and h<=7):
            hr="00"
            fr="06"
        elif (h>7 and h<=9):
            hr="00"
            fr="09"
        elif (h==10):
            hr="06"
            fr="03"
        elif (h>10 and h<=13):
            hr="06"
            fr="06"
        elif (h>13 and h<=15):
            hr="06"
            fr="09"
        elif (h==16):
            hr="12"
            fr="03"
        elif (h>16 and h<=19):
            hr="12"
            fr="06"
        elif (h>19 and h<=21):
            hr="12"
            fr="09"
        elif (h==22):
            hr="18"
            fr="03"
        else:
            hr="18"
            fr="06"
        ls=[yr,mnt,day,hr,fr]
        wave_t=UTCDateTime(int(ls[0]),int(ls[1]),int(ls[2]),int(ls[3]))
        self.wv_name=str(wave_t)[:13]+"_fr_"+fr+".grib2"
        return ls
    
    def download(self):
        '''This function downloads waveheight data
        using the prameters defined by str_maker()
        This function saves the 
        significant wave height data to the ./wave
        directory and make sure the total size of this
        directory stays below 1 GB'''
        
        ls=self.str_maker()
        yr,mnt,day,hr,fr=ls
        ##creating url
        str1="https://nomads.ncep.noaa.gov/cgi-bin/filter_gefs_wave_0p25.pl?file=gefs.wave.t"
        str2="{hr}z.c00.global.0p25.f0{fr}."
        str3="grib2&lev_surface=on&var_HTSGW=on&subregion=&leftlon=150&rightlon=260&toplat=78&bottomlat=40&dir=%2Fgefs."
        str4="{yr}{mnt}{day}%2F{hr}%2Fwave%2Fgridded"
        url=str1+str2.format(hr=hr,fr=fr)+str3+str4.format(yr=yr,mnt=mnt,day=day,hr=hr)
        print("downloading",self.wv_name)
        if not os.path.isdir(jn(root_path,"wave/")):
            os.makedirs(jn(root_path,"wave/"))
        os.chdir(jn(root_path,"wave/"))
        wget.download(url) #downloading
        if get_folder_size(jn(root_path,"wave"))>1024: #ensuring the size of wave height directory is below 1 gb
            open_g=[]
            for root,dire,file in os.walk(jn(root_path,"wave/")):
                open_g.append(file)
            
            def custom_sort(file_name):
                parts = file_name.split('_')
                date_part = parts[0]
                fr_part = parts[2].split('.')[0]
                # Convert date and time to datetime object
                dt = datetime.strptime(f"{date_part}", "%Y-%m-%dT%H")
                # Add fr_part to hours and handle day increment
                dt += timedelta(hours=int(fr_part))
                return dt
            
            sorted_list = sorted(open_g[0], key=custom_sort)
            os.remove((os.path.join(root,sorted_list[0])))
        for file in glob.glob(jn(root_path,"wave/gefs.*.grib2")):
            os.rename(file,jn(root_path,"wave/",self.wv_name))
        os.chdir(root_path)
        return str(url)
        

def plot_map(title,stations,ms_arr_fl,wv_name):
    '''This function plots the final map and save
    it to the folder ./hourly_maps
    Function make use of ./alaska_network_station_location.csv
    to find station location
    Uses files in ./cpts and ./grids folder
    Uses of wave height file from ./wave directory
    
    :param title:(datetime) title of the plot 
    :param stations:(list/array) stations to be plotted
    :param ms_arr_fl:(int) filteres secondary microseisms (5-10s)
    :param wv_name:(str) name of significant wave height data
    to be plotted'''
    
    #function to plot pygmt map
    root=jn(root_path,"wave/")
    grbs=pygrib.open(os.path.join(root,wv_name))# loading wave data
    grb=grbs[1]
    data_wave = grb.values
    latg, long = grb.latlons()
    x_s=np.arange(150,260,0.15)
    y_s=np.arange(30,78,0.15) 
    grid_x, grid_y = np.meshgrid(x_s,y_s)
    lon=long.flatten()
    lat=latg.flatten()
    for i in range(len(lon.flatten())):
        if lon[i]<0:
            lon[i]=lon[i]+360
    points=np.array((lon,lat)).T
    values=data_wave.filled(fill_value=np.nan).flatten()
    #interpolating wave height data so that the plot looks smooth
    grid_z0 = griddata(points, values, (grid_x, grid_y), method='linear')
    
    grid = xr.DataArray(
        grid_z0, dims=["lat", "lon"], coords={"lat": y_s, "lon": x_s}
        )
    stat=stations
    stat_loc=pd.read_csv("./Alaska_network_station_location.csv")
    lon_psd=np.array([])
    lat_psd=np.array([])
    for ele in stat:
        lon_psd=np.append(lon_psd,(stat_loc[stat_loc["Station Code"]==ele]["Longitude"].iloc[0]))
        lat_psd=np.append(lat_psd,(stat_loc[stat_loc["Station Code"]==ele]["Latitude"].iloc[0]))
    
    lon_psd[np.where(lon_psd[lon_psd<0])]=lon_psd[np.where(lon_psd[lon_psd<0])]+360
    ms_arr_fl[np.isnan(ms_arr_fl)]=-180
    #grid1=pygmt.datasets.load_earth_relief(resolution='02m', region=[150,260, 30, 78])
    fig=pygmt.Figure()
    reg="175/40/257/72r"
    proj="L-159/35/33/45/10c"
    with pygmt.config(MAP_FRAME_TYPE="plain",MAP_FRAME_PEN="1p,black"):
        fig.basemap(region=reg, projection=proj,frame="lrbt")
    fig.grdimage(grid="./grids/dem.nc",region=reg,projection=proj, cmap="./cpts/topo.cpt",nan_transparent=True)
    fig.grdimage(grid=grid,region=reg,projection=proj, cmap="cpts/expll1.cpt",nan_transparent=True)
    fig.coast(region=reg, projection=proj,shorelines="0.05p,black",borders=["1/0.1p,black","2/0.05p,gray30"],area_thresh='600')
    with pygmt.config(MAP_FRAME_TYPE="plain",MAP_FRAME_PEN="0.5p,black",FONT="10p,black"):
        fig.colorbar(projection=proj, cmap="cpts/expll1.cpt",frame=["x+locean\twave\theight", "y+lm"],position="n0.55/-0.05+w4c/0.3c+h")
        fig.colorbar(cmap="cpts/psd.cpt",frame=["x+lseismic\tpower", "y+ldb"],projection=proj,position="n0.02/-0.05+w4c/0.3c+h")  
    fig.plot(x=lon_psd,y=lat_psd,style="c0.2c",fill=ms_arr_fl,cmap="cpts/psd.cpt",pen="black", projection=proj)
    fig.text(text=str(title)[0:16],x=-160,y=74,projection=proj)
    if not os.path.exists(jn(root_path,"hourly_maps")):
            os.makedirs(jn(root_path,"hourly_maps"))
    fig.savefig(jn(root_path,"hourly_maps/")+str(title)[0:16]+".jpg")
    if get_folder_size(jn(root_path,"hourly_maps"))>500: #ensuring the size of the hourly_maps folder stays below 500 mb
        open_g=[]
        for root,dire,file in os.walk(jn(root_path,"hourly_maps")):
            for name in file:
                open_g.append(name)
        open_g.sort()
        os.remove((os.path.join(root,open_g[0])))
    return fig 



print("using default st and et")
timec=datetime.utcnow().replace(microsecond=0, second=0, minute=0)
time=UTCDateTime(timec)-(24*3600) #setting start time as 24 hr behind current time
print("processing ",time)
AKDT= pytz.timezone('America/Anchorage')
titime=datetime.now(AKDT).replace(microsecond=0, second=0, minute=0)
title=titime-timedelta(hours=23,minutes=30)
ms=SeismicStorms(time)
ms.add_value(1)# 1 here is a flag to turn on time check this should be on normally but turning this off can help in trouble shooting
ms.filter_seis(7)
wv=Wave(time)
ls=wv.str_maker()
url=wv.download()
plot_map(title,ms.stations,ms.ms_arr_fl,wv.wv_name)
os.chdir(root_path)
# with open("log.txt", "a") as file:
#     current_time=datetime.now(AKDT).strftime("%Y-%m-%d %H:%M:%S")
#     file.write(current_time + "\n")

    



     
   
       
# connection = sqlite3.connect(jn(root_path,"PSDs/seismic_database.db"))
# cursor = connection.cursor()       
# cursor.execute('SELECT * FROM PSD')


# # Fetch all rows
# rows = cursor.fetchall()       
# print("station | time                 | power")
# print("-" * 40)

# # Print each row
# for row in rows:
#     # Assuming the time is stored as a floating-point number in the database
#     print(f"{row[0]:<8} | {row[1]:19} | {str(row[2])[:7]}")

# # Close the database connection
# connection.close()


