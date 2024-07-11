# Data Visuals Repository
This repository is where all code created for the [Data Visuals Project](https://app.asana.com/0/1205732312106349/1205732312106361) will be stored.  

*README updated :: 20 October, 2023 by Gabe Paris*

### Author
Sebin John

## Command Line Usage

python ./SeismicStorms.py 

## Description

SeismicStorms.py generates a map that displays significant wave height data and microseismic data from different 
stations specified by the station list file. Significant wave height is obtained from the NOAA wave watch 3 models.
Seismic data is obtained from the AEC database. SeismicStorms.py always generates a map 24 hours behind the current time.
This is because the availability of wave models is constrained in the NOAA data server.

output location of the map is ./hourly maps/{datetime}.png
SeismicStorms.py will generate a ./wave/, ./PSDs/, ./RESPs/ to store the data and metadata.

Wave directory: Stores wave height data. The size of the wave directory is limited to 1 GB.

PSDs directory: stores local database of PSDs (microseismic) from seismic stations. This SQL database inside this directory named
"seismic_database.db" only has one table. 

RESPS directory: stores response files for seismic stations. If not available checks IRIS to download.

## Required Files

For this code to run there should be these files in the same directory as the code.

station list file - contains comma-delimited station names to be plotted

grids directory - contains topography grid and shading grid for the map

cpts directory - Specifies color map to be used. (This is hardcoded)

Alaska_network_station_location.csv - Locations of the stations in the AK network

medfilt.py - median filtering Python module

## Debugging

./error_log.txt file stores the standard output and errors for debugging

## Dependencies

SeismicStorms.py depends on the following:  

1. Conda python  
2. Matplotlib  
3. Cartopy
4. tqdm
5. numpy
6. obspy
7. pickle
8. pandas
9. scipy
10. sqlite3
11. wget
12. pygmt
13. xarray
14. pygrib
15. pytz
16. globe

## Caveats

The user doesn't have control over the time plotted, which is a serious limitation. Attempts to address this can be made with feature requests.

This code will only run in systems connected to the AEC database using Antelope

## Example

python ./SeismicStorms.py 

```processing  2024-04-23T19:00:00.000000Z
Error processing station S12K: No response data available (HTTPError: HTTP Error 404: )
Error processing station KINC: No response data available (HTTPError: HTTP Error 404: )
Error processing station LSSA: No response data available (HTTPError: HTTP Error 404: )
previous runtime is 2024-04-23T19:00:00.000000Z
connecting to database
downloading 2024-04-23T12_fr_06.grib2
100% [..............................................................................] 40756 / 40756
```

* IF no response is available please remove such station from the station list or add the response file manually to the RESP directory
