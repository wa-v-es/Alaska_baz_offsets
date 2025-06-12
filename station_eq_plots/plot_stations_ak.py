import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
import obspy
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

# boundaries = requests.get("https://raw.githubusercontent.com/fraxen/tectonicplates/master/GeoJSON/PB2002_boundaries.json").json()
def add_bathymetry(m):
    etopo_data = Dataset('etopo1.nc')
    lons = etopo_data.variables['X'][:]
    lats = etopo_data.variables['Y'][:]
    topoin = etopo_data.variables['Z'][:]

    topoin, lons = shiftgrid(180., topoin, lons, start=False)
    lons, lats = np.meshgrid(lons, lats)
    x, y = m(lons, lats)

    ocean_topo = maskoceans(lons, lats, topoin, inlands=False)
    m.contourf(x, y, ocean_topo, cmap='gray', alpha=0.5)
####
station='ILAR'
sta_lat= 64.768097
sta_long=-146.783203

station='AL_NW'
sta_lat= 70
sta_long=-157

station='AL_SE'
sta_lat_= 59.5
sta_long_=-136.5

client = Client("IRIS")

###########
# catfile = 'events_ak_2024_6.7.xml'
#
# catfile_NW = 'events_ak_NW_2024_6.7.xml'
# catfile_SE = 'events_ak_SE_2024_6.7.xml'
### get eq data
#earthquakes
starttime= UTCDateTime('2020-11-01T00:00:01')
endtime= UTCDateTime('2024-05-14T00:00:01')


inventory = client.get_stations(network="AK,CN,AV,AT",minlatitude=55,maxlongitude=-137,
minlongitude=-168,starttime=starttime,endtime=endtime)#,level='response')
# fig=inventory.plot(show=False,color='indianred')
# catalog.plot(fig=fig)
# fig = plt.figure(figsize=(8, 8))
network_data = {'AK': [], 'CN': [], 'AV': [], 'AT': []}

# Extract and group by network
for network in inventory:
    for station in network:
        lon = station.longitude
        lat = station.latitude
        name = station.code
        net = network.code
        if net in network_data:
            network_data[net].append((lon, lat, name))

# Write each network's data to a file
for net, records in network_data.items():
    with open(f"{net}_stations.txt", "w") as f:
        for lon, lat, name in records:
            f.write(f"{lon:.5f} {lat:.5f} {name}\n")

# inventory.plot('local',show=False,resolution='h',water_fill_color='white',color_per_network={"AK": 'darkseagreen', "CN": "xkcd:heather","AV": "xkcd:muted blue","AT": "xkcd:light burgundy"},label=False,marker='^',size=22, fontsize=45,edgecolors='b',linewidths=3)

# plt.show()
sys.exit()

# plt.savefig('AK_all.pdf',bbox_inches='tight', dpi=500,pad_inches=0.1)
########
# method #2
## ran the next bit from ipython interface

# data = os.path.join('sac_files/230702_102743_PA/', '')
file_path = 'sac_files/230702_102743_PA/'
st_all=Stream()

for filename in glob.glob("sac_files/230702_102743_PA/*.sac"):
    st=read(filename,format='sac')
    st_all.extend(st)

print(len(st_all))
latlist=[]
lonlist=[]

latlist = [tr.stats.sac['stla'] for tr in st_all]
lonlist = [tr.stats.sac['stlo'] for tr in st_all]
plt.ion()

fig = plt.figure(figsize=(8, 8))
m = Basemap(projection='lcc', resolution='i',width=1.85E6, height=1.85E6,lat_0=63, lon_0=-152)
# m.etopo(scale=1, alpha=0.9,cmap='Grays')
# m.shadedrelief(scale=5, alpha=0.75)#, cmap='Grays')
m.drawcoastlines(color='dimgray',linewidth=.8)
m.drawcountries(color='gray',linewidth=.7)
m.fillcontinents(color="seashell")#, lake_color='#DDEEFF')
m.drawstates(color='gray')
m.scatter(lonlist, latlist, latlon=True,s=55, marker='^', facecolor='white', edgecolor='xkcd:muted blue',alpha=.85)
m.scatter(-146.886597,64.7714, latlon=True,s=55, marker='v',facecolor='none', edgecolor='xkcd:muted blue',alpha=.1,zorder=5)
m.drawmeridians(range(190, 238, 3), color='k', linewidth=.25, dashes=[4, 4], labels=[0, 0, 0, 1])
m.drawparallels(range(53, 74, 1), color='k', linewidth=0.25, dashes=[4, 4], labels=[1, 0, 0, 0])
# plt.savefig('AK_230702_102743_PA_tr.png',dpi=500,bbox_inches='tight', pad_inches=0.1)
# fig.show()
sys.exit()
