import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
import obspy
from obspy.geodetics import gps2dist_azimuth
from obspy import read, Stream, UTCDateTime,read_events
from obspy.core.event import Origin
from obspy.clients.fdsn import Client
import os
import cartopy.crs as ccrs
import sys
from netCDF4 import Dataset
import circle as cir_robin
from importlib import reload
reload(cir_robin)
import requests
import glob as glob
import imageio.v2 as imageio
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


# sys.exit()

def plot_stations(file_name):
    ##
    coordinates = read_station_coordinates(file_name)
    if not coordinates:
        # return None
        print('')
    lats, lons = zip(*coordinates)
    ####
    # latlist_all,lonlist_all=get_all_AK_st()
    ##
    latlist_all=[]
    lonlist_all=[]
    for line in open('AK_4arrays.txt','r'):
        line=line.split()
        latlist_all.append(float(line[0]))
        lonlist_all.append(float(line[1]))

    fig = plt.figure(figsize=(8, 8))
    m = Basemap(projection='lcc', resolution='i',width=1.85E6, height=1.85E6,lat_0=63, lon_0=-152)
    # m.etopo(scale=1, alpha=0.9,cmap='Grays')
    # m.shadedrelief(scale=5, alpha=0.75)#, cmap='Grays')
    m.drawcoastlines(color='dimgray',linewidth=.8)
    m.drawcountries(color='gray',linewidth=.7)
    m.fillcontinents(color="seashell")#, lake_color='#DDEEFF')
    m.drawstates(color='gray')
    m.scatter(lonlist_all, latlist_all, latlon=True,s=55, marker='^', facecolor='white', edgecolor='xkcd:muted blue',alpha=.85)
    ######
    x, y = m(lons, lats)

    apt,a=get_aperture(lats,lons)
    m.scatter(x, y, marker='^', c='xkcd:light burgundy', s=65,alpha=.75)


    m.scatter(-146.886597,64.7714, latlon=True,s=55, marker='v',facecolor='none', edgecolor='xkcd:muted blue',alpha=.1,zorder=5)
    m.drawmeridians(range(190, 238, 3), color='k', linewidth=.25, dashes=[4, 4], labels=[0, 0, 0, 1])
    m.drawparallels(range(53, 74, 1), color='k', linewidth=0.25, dashes=[4, 4], labels=[1, 0, 0, 0])
    num_stations = len(coordinates)
    plt.text(0.6, 0.2, f'Number of Stations: {num_stations}',transform=plt.gca().transAxes,fontsize=11, c='xkcd:light burgundy', verticalalignment='top')
    plt.text(0.6, 0.15, f'Aperture: {apt/1000:.2f} km', transform=plt.gca().transAxes, fontsize=11, c='xkcd:muted blue', verticalalignment='top')

    output_file = file_name.replace('.txt', '.png')
    output_file=output_file.replace('st_files','st_imgs')
    plt.savefig(output_file,dpi=300,bbox_inches='tight', pad_inches=0.1)
    plt.close(fig)
    return output_file
#
def calculate_aperture_all_files(parent_folder, pattern_stfiles):
    # Initialize a dictionary to store the extracted numbers by folder
    # folder_numbers = {}
    apt_all=[]
    # Loop through each folder in the parent folder
    for root, dirs, files in os.walk(parent_folder):
        # Find files matching the pattern
        # numbers = []
        for file in glob.glob(os.path.join(root, pattern_stfiles)):
            file_path = os.path.join(parent_folder, file)
            coordinates = read_station_coordinates(file_path)
            if not coordinates:
                # return None
                print('')
            lats, lons = zip(*coordinates)
            apt,_=get_aperture(lats,lons)
            apt_all.append(apt)

    return np.array(apt_all)
###
# image_files = []
# folder_path='/Users/keyser/Research/AK_all_stations/sac_files/221109_093842_PA_inc2_r2.5/st_files'
# for file_name in os.listdir(folder_path):
#     if file_name.endswith('.txt'):  # Adjust the extension based on your files
#         file_path = os.path.join(folder_path, file_name)
#         output_file = plot_stations(file_path)
#         # sys.exit()
#         if output_file:
#             image_files.append(output_file)

###
pattern_stfiles = 'STA_DISTANCE_LOC_gridnumber*.txt'
parent_folder='/Users/keyser/Research/AK_all_stations/sac_files'
apt_all=calculate_aperture_all_files(parent_folder, pattern_stfiles)
apt=apt_all/1000
ave_aperture=np.mean(apt)

lower, upper = np.percentile(np.array(apt), [5, 95])
print(f"90% of the values lie between {lower:.2f} and {upper:.2f}")

sys.exit()

apt_all=[]
for root, dirs, files in os.walk(parent_folder):
    print(root)
    for file_name in os.listdir(folder_path):

        if file_name.endswith('.txt'):  # Adjust the extension based on your files
            file_path = os.path.join(folder_path, file_name)
            coordinates = read_station_coordinates(file_path)
            if not coordinates:
                # return None
                print('')
            lats, lons = zip(*coordinates)
            apt,_=get_aperture(lats,lons)
            apt_all.append(apt)


# Create GIF

###### do following three steps on terminal after making st_imgs folder in eq directory
# mv -f STA_D* /Users/keyser/Research/AK_all_stations/sac_files/230615_180628_PA/st_imgs/
# /Users/keyser/Research/AK_all_stations/sac_files/230615_180628_PA/st_imgs/
# magick -delay 80 -quiet -quality 50 -loop 0 *.png out.gif

create_gif(image_files, 'stations_animation_4.gif')

sys.exit()
