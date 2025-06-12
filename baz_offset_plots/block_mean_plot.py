import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
import obspy
from obspy.geodetics import gps2dist_azimuth
from obspy import read, Stream, UTCDateTime,read_events
from obspy.core.event import Origin
from obspy.clients.fdsn import Client
import os
import glob as glob

import cartopy.crs as ccrs
import sys
from netCDF4 import Dataset
import circle as cir_robin
from importlib import reload
reload(cir_robin)
import requests
import glob as glob
import imageio.v2 as imageio
from IPython import get_ipython
import cartopy.crs as ccrs


######

def get_aperture(stlat,stlong):
    stlat= [i for i in stlat]
    stlong= [i for i in stlong]
    max=[]
    for i in range(len(stlat)):
        lat,long=stlat[i],stlong[i]
        dist=[]
        for j in range(len(stlat)):

            d,v,b=gps2dist_azimuth(lat, long, stlat[j], stlong[j])
            dist.append(d)

        max.append(np.max(dist))
    return np.max(max),max

# Function to extract the numeric part from the filenames
def extract_number(filename):
    match = re.search(r'(\d+)', filename)
    return int(match.group(0)) if match else 0
##
# Function to create GIF from images
def create_gif(image_files, output_gif, duration=2):
    # Sort the files based on the numeric part of their names
    image_files.sort(key=extract_number)

    images = []
    for filename in image_files:
        images.append(imageio.imread(filename))
    imageio.mimsave(output_gif, images, duration=duration)

#
def read_station_coordinates(file_name):
    coordinates = []
    with open(file_name, 'r') as file:
        for line in file:
            parts = line.split(':')
            if len(parts) > 2:
                lat, lon = map(float, parts[2].strip().split()[:2])
                coordinates.append((lat, lon))
    return coordinates
##
# not using this fucntion anymnore. made a list of lat longs
def get_all_AK_st():
    # file_path = 'sac_files/230702_102743_PA/'
    st_all=Stream()

    for filename in glob.glob("sac_files/230702_102743_PA/*.sac"):
        st=read(filename,format='sac')
        st_all.extend(st)

    print(len(st_all))
    latlist=[]
    lonlist=[]

    latlist = [tr.stats.sac['stla'] for tr in st_all]
    lonlist = [tr.stats.sac['stlo'] for tr in st_all]
    return latlist,lonlist

###
def plot_baz_circle(m,file,marker):
    for line in open(file,'r'):
        line=line.split()
        if float(line[2])-1.25*float(line[3]) < float(line[1]) < float(line[2])+1.25*float(line[3]):
            if float(line[1]) in [-1,-0.5,0,1]:
                m.scatter(float(line[5]),float(line[4]), latlon=True,s=165, marker=marker, facecolor='none', edgecolor='black',alpha=.75,linewidth=.75)
            if float(line[1]) > 1.1:
                m.scatter(float(line[5]),float(line[4]), latlon=True,s=abs(float(line[1]))*165, marker=marker, facecolor='darkred', edgecolor='white',alpha=.125)
            if float(line[1]) < -1.1:
                m.scatter(float(line[5]),float(line[4]), latlon=True,s=abs(float(line[1]))*165, marker=marker, facecolor='darkblue', edgecolor='white',alpha=.125)
        else:
            print('Grid:{} max not reprensentative of vespa'.format(float(line[0])))

    ##

# coordinates = read_station_coordinates(file_name)
# if not coordinates:
#     # return None
#     print('')
# lats, lons = zip(*coordinates)
####
# latlist_all,lonlist_all=get_all_AK_st()
##
get_ipython().magic('reset -sf')

latlist_all=[]
lonlist_all=[]
for line in open('AK_4arrays.txt','r'):
    line=line.split()
    latlist_all.append(float(line[0]))
    lonlist_all.append(float(line[1]))

plt.ion()
fig = plt.figure(figsize=(8, 8))
m = Basemap(projection='lcc', resolution='i',width=1.85E6, height=1.85E6,lat_0=63, lon_0=-152)
# m.etopo(scale=1, alpha=0.9,cmap='Grays')
# m.shadedrelief(scale=5, alpha=0.75)#, cmap='Grays')
m.drawmapboundary(fill_color='aliceblue')
m.drawcoastlines(color='dimgray',linewidth=.8)
m.drawcountries(color='gray',linewidth=.7)

m.fillcontinents(color="white", lake_color='aliceblue')

m.drawstates(color='gray')
m.drawmeridians(range(190, 238, 3), color='k', linewidth=.25, dashes=[4, 4], labels=[0, 0, 0, 1])
m.drawparallels(range(53, 74, 1), color='k', linewidth=0.25, dashes=[4, 4], labels=[1, 0, 0, 0])
m.scatter(-146.886597,64.7714, latlon=True,s=55, marker='v',facecolor='none', edgecolor='xkcd:muted blue',alpha=.1,zorder=5)


m.scatter(lonlist_all, latlist_all, latlon=True,s=35, marker='^', facecolor='white', edgecolor='xkcd:muted blue',alpha=.85)

#SA_eq_2021102629_maxVals_stmax.txt
# plot_baz_circle(m,'SA_eq_202111281052_maxVals_stmax.txt','p')
#
# folder_pattern_pa = "vals_zmax/*.txt"
# folder_pattern_sa = "vals_zmax/SA/*.txt"
# matching_files_pa = glob.glob(folder_pattern_pa)
# matching_files_sa = glob.glob(folder_pattern_sa)
#
# for file in matching_files_sa:
#     plot_baz_circle(m,file,'o')

for line in open('block_mean_vals/block_mean_mean_low_slow_PA.txt','r'):
    line=line.split()
    if float(line[2]) in [-1,-0.5,0,1]:
        # print(' ')
        m.scatter(float(line[0]),float(line[1]), latlon=True,s=abs(float(line[2]))*165, marker='o', facecolor='none', edgecolor='black',alpha=.75,linewidth=.75)
    elif float(line[2]) > 1.1:
        m.scatter(float(line[0]),float(line[1]), latlon=True,s=abs(float(line[2]))*165, marker='o', facecolor='darkred', edgecolor='white',alpha=.25)
        m.scatter(float(line[0]),float(line[1]), latlon=True,s=abs(float(line[3]))*165, marker='|', facecolor='darkred',alpha=.75)

    elif float(line[2]) < 1.1:
        m.scatter(float(line[0]),float(line[1]), latlon=True,s=abs(float(line[2]))*165, marker='o', facecolor='darkblue', edgecolor='white',alpha=.25)
        m.scatter(float(line[0]),float(line[1]), latlon=True,s=abs(float(line[3]))*165, marker='|', facecolor='darkblue',alpha=.75)


# plot_baz_circle(m,'eq_202112291825_maxVals_1.75.txt','d')#'*',
# plot_baz_circle(m,'eq_2023721027_maxVals_stmax.txt','>')

x,y=m(-141,56)
x2,y2 = m(-141,58.5)
circle1 = plt.Circle((x, y), y2-y, color='beige',fill=True)
# circle = Circle(xy=m(-146.886597,64.7714),radius=5, fill=False)
plt.gca().add_patch(circle1)

m.drawgreatcircle(-141, 56, 174.90, -21.13, del_s=100.0) # 2021 10 02 06:29
m.drawgreatcircle(-141, 56, 127.58, -7.55,  del_s=100.0) # 2021 12 29 18:25
m.drawgreatcircle(-141, 56, 146.50, -6.29, del_s=100.0)
m.drawgreatcircle(-141, 56, 130.00, -7.06, del_s=100.0)
m.drawgreatcircle(-141, 56, -174.95, -17.88, del_s=100.0)
m.drawgreatcircle(-141, 56, 110.70, -5.60, del_s=100.0)
m.drawgreatcircle(-141, 56, 123.48, -6.69, del_s=100.0)

# ### SA
# m.drawgreatcircle(-141, 56, -76.81, -4.45, del_s=100.0) # 2021 11 28 10:52
# m.drawgreatcircle(-141, 56, -66.65, -23.50, del_s=100.0) # 2022 05 10 23:06
# m.drawgreatcircle(-141, 56, -70.20, -14.89, del_s=100.0) # 2022 05 26 12:02
# m.drawgreatcircle(-141, 56, -63.10, -26.75, del_s=100.0) # 2023 01 20 22:09

plt.savefig('bl_mean_mean_low_slow_PA.png',dpi=200,bbox_inches='tight', pad_inches=0.1)

# plt.close()
# for line in open('eq_202111281052_maxVals_stmax.txt','r'):
#     line=line.split()
#     if float(line[1]) == 0:
#         m.scatter(float(line[5]),float(line[4]), latlon=True,s=165, marker='d', facecolor='none', edgecolor='black',alpha=.75)
#     if float(line[1]) > 0:
#         m.scatter(float(line[5]),float(line[4]), latlon=True,s=abs(float(line[1]))*165, marker='d', facecolor='darkred', edgecolor='white',alpha=.15)
#     if float(line[1]) < 0:
#         m.scatter(float(line[5]),float(line[4]), latlon=True,s=abs(float(line[1]))*165, marker='d', facecolor='darkblue', edgecolor='white',alpha=.15)

#SA event
# -4.45 -76.81 126.00 -- 20217 11 28 10:52
# -23.50 -66.65 220.00 2022 05 10 23:06
# -14.89 -70.20 252.00 2022 05 26 12:02
# -26.75 -63.10 596.70 2023 01 20 22:09
#PA
### -21.13 174.90 527.00 -- 20217 10 02 06:29
## -7.55 127.58 164.30 -- 20217 12 29 18.25
## -6.29 146.50 116.00 -- 2022 09 10 23:47
## -7.06 130.00 105.20 -- 2023 01 09 17:47
## -17.88 -174.95 229.00 -- 2023 07 02 10-27
# -5.60 110.70 538.70 2020 07 06 22:54
# -6.69 123.48 627.80 2020 08 21 04:09

sys.exit()

#
