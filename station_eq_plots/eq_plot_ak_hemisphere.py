import matplotlib.pyplot as plt
import numpy as np
import obspy
from obspy import read, Stream, UTCDateTime,read_events
from obspy.core.event import Origin
from obspy.clients.fdsn import Client
import os
import cartopy.crs as ccrs
import sys
import cartopy.feature as cfeature
import circle as cir_robin
from importlib import reload
reload(cir_robin)
import requests
import cartopy.io.img_tiles as cimgt
# from cartopy.io.img_tiles import StamenTerrain


# tiler = cimgt.Stamen('toner-lite')  # Greyscale tile set
# tiler = StamenTerrain()
# mercator = tiler.crs

boundaries = requests.get("https://raw.githubusercontent.com/fraxen/tectonicplates/master/GeoJSON/PB2002_boundaries.json").json()

####
station='ILAR'
sta_lat= 64.768097
sta_long=-146.783203

station='AL'
sta_lat= 65
sta_long=-148

client = Client("IRIS")

###########
catfile = '../xmls/events_6.7_2024mw_AK_allmag.xml'
# catfile = 'events_ZS_2019_6.7.xml'
### get eq data
#earthquakes
starttime= UTCDateTime('2020-11-01T00:00:01')
endtime= UTCDateTime('2024-05-01T00:00:01')

# starttime= UTCDateTime('2017-10-01T00:00:01') # for ZS
# endtime= UTCDateTime('2019-12-31T00:00:01') # for ZS
catalog=read_events(catfile)

lats, longs = [], []
mags = []
azi=[]
#
for event in catalog:
    lats.append(event.origins[0].latitude)
    longs.append(event.origins[0].longitude)
    mags.append(event.magnitudes[0].mag)
#
fig, ax=plt.subplots(figsize=(4,3))
plt.axis('off')

# ax = plt.axes(projection=ccrs.Mollweide(central_longitude=sta_long))
ax = plt.axes(projection=ccrs.LambertConformal(central_longitude=-155, central_latitude=0,cutoff=-35))
# ax = plt.axes(projection=ccrs.AzimuthalEquidistant(central_longitude=sta_long,central_latitude=sta_lat))
ax.set_global()
# ax.stock_img()
ax.set_facecolor('none')
ax.add_feature(cfeature.OCEAN.with_scale('110m'), facecolor='darkgrey', zorder=0)
ax.add_feature(cfeature.LAND.with_scale('110m'), facecolor='lavenderblush', edgecolor='black', linewidth=0.2, zorder=1)
# ax.add_image(tiler, 4)
# ax.coastlines(color='black', linewidth=.15,)
ax.plot(sta_long, sta_lat, color='indigo', marker='^', markersize=5, transform=ccrs.Geodetic())
# ax.plot(sta_long_, sta_lat_, color='indigo', marker='^', markersize=1, transform=ccrs.Geodetic())

min_marker_size = .5
for i in range(len(lats)): #plot eqs
    # x,y = eq_map(lon, lat)
    msize = mags[i] * min_marker_size
    # marker_string = get_marker_color(mag)
    ax.plot(longs[i], lats[i],color='darkred',marker='o',markersize=msize,alpha=.8,transform=ccrs.Geodetic())

# X,Y=cir_robin.equi(sta_long, sta_lat, 7215)
# X1,Y1=cir_robin.equi(sta_long, sta_lat, 11100)
#
# plt.plot(X,Y,transform=ccrs.Geodetic(),lw=.7,alpha=.6,linestyle='--',c='maroon')
# plt.plot(X1,Y1,transform=ccrs.Geodetic(),lw=.7,alpha=.6,linestyle='--',c='maroon')

#plot gcp
for event in catalog:
    ori= event.preferred_origin()
    plt.plot([sta_long, ori.longitude],[sta_lat, ori.latitude],  transform=ccrs.Geodetic(),color='black',lw=.25,linestyle='dotted',alpha=.85)

# Plot boundaries.
for f in boundaries["features"]:
    c = np.array(f["geometry"]["coordinates"])
    lng, lat = c[:, 0], c[:, 1]
    x, y = lng, lat
    mask = np.bitwise_or(np.abs(x) > 1e15, np.abs(y) > 1e15)
    x = np.ma.array(x)
    y = np.ma.array(y)
    x.mask = mask
    y.mask = mask
    plt.plot(x, y, color="Navy", lw=.35,transform=ccrs.Geodetic())

# ax.text(sta_long+5, sta_lat-3, 'AK_NW',fontsize=8,fontfamily='serif', color='indigo',transform=ccrs.Geodetic())
# ax.text(sta_long_+5, sta_lat_-3, 'AK_SE',fontsize=8,fontfamily='serif', color='indigo',transform=ccrs.Geodetic())

# ax.text(-145, -8, '65°',fontsize=10,fontfamily='serif', color='maroon',transform=ccrs.Geodetic())
# ax.text(-145, -42, '100°',fontsize=10,fontfamily='serif', color='maroon',transform=ccrs.Geodetic())
# ax.axis('off')
# ax.set_axis_off()
# ax.set_frame_on(False)

for pos in ['right', 'top', 'bottom', 'left']:
    plt.gca().spines[pos].set_visible(False)

# plt.show()


plt.savefig('eq_AK_all_conf.png',dpi=300,bbox_inches='tight', pad_inches=0.05,transparent=True)
