###
import numpy as np
# from numpy import load
import obspy
import miller_alaskamoho_srl2018 as alaskamoho
import cartopy
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.feature as cfeature
import xarray as xr
from cmcrameri import cm
import matplotlib.colors as colors

import os.path as path
import stripy as stripy
import sys
from matplotlib.colors import ListedColormap

#######
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
ds = xr.open_dataset('sediment/sedmap.grd')
stations = '/Users/keyser/Research/TOMOGRAD-main/STA_DISTANCE_LOC_gridnumber97.txt'
coordinates = read_station_coordinates(stations)
if not coordinates:
    # return None
    print('')
lats_st, lons_st = zip(*coordinates)

# xmin, xmax = ds['x_range'].values
# ymin, ymax = ds['y_range'].values
# dx, dy = ds['spacing'].values
#
# lons = np.arange(xmin, xmax + dx, dx)
# lats = np.arange(ymin, ymax + dy, dy)
#
# z_flat = ds['z'].values
#
# nx = len(lons)
# ny = len(lats)
#
# Z = z_flat.reshape(ny, nx)
#
#
# # Create meshgrid
# lon2d, lat2d = np.meshgrid(lons, lats)
lon, lat, Z = np.loadtxt("sediment/sedmap.gmt", unpack=True)

ulon = np.unique(lon)
ulat = np.unique(lat)

Zgrid = Z.reshape(len(ulat), len(ulon))

lon2d, lat2d = np.meshgrid(ulon, ulat)

# Create new colormap
# cmapA = ListedColormap(colA)
# cmapA = cmap
proj = ccrs.Stereographic(central_longitude=-154, central_latitude=90, true_scale_latitude=60)
# plt.clf()
plt.ion()
fig = plt.figure(figsize=(15, 8), facecolor=None)
ax1 = plt.subplot(1, 1, 1, projection=proj)
plt.rcParams.update({'font.size': 14})


# ax1.gridlines(draw_labels=False, linewidth=0.5, color='gray', alpha=0.5)

ax1.set_extent([-165,-138,55,68.5], crs=ccrs.PlateCarree())

# grat = cartopy.feature.NaturalEarthFeature(category="physical", scale="10m", name="graticules_5")
# ax1.add_feature(grat, linewidth=0.5,linestyle="--",edgecolor="#000000",facecolor="None", zorder=2)

ax1.coastlines(resolution="10m",color="#111111", linewidth=0.5, zorder=99)
ax1.add_feature(cfeature.BORDERS.with_scale('10m'), linestyle=':')
ax1.add_feature(cfeature.OCEAN.with_scale('10m'))#,alpha=0.2,facecolor='xkcd:dusty blue')

gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=.8, color='gray', alpha=0.5, linestyle='--',rotate_labels=False,
        x_inline=False, y_inline=False)

gl.xlocator = mticker.FixedLocator([-160,-150,-140])
gl.ylocator = mticker.FixedLocator([55,60,65,70])

# gl.xlines = True
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.top_labels = False
gl.left_labels = False
gl.xlabel_style = {'size': 15}
gl.ylabel_style = {'size': 15}

cmap_c=cm.batlowW_r
img = ax1.pcolormesh(
    lon2d, lat2d, Zgrid,
    cmap=cmap_c,
    vmax=.7,
    vmin=0,
    transform=ccrs.PlateCarree(),
    shading='auto'
)#norm=colors.Normalize(vmin=-5, vmax=5),
# sys.exit()
ax1.plot(lons_st, lats_st, marker='^',markersize=8, linestyle='None', markerfacecolor='none', markeredgecolor='navy', transform=ccrs.PlateCarree())
# ax1.plot(np.mean(lons_st), np.mean(lats_st), marker='^',markersize=13, linestyle='None', markerfacecolor='royalblue', markeredgecolor='navy', transform=ccrs.PlateCarree())

cbar = plt.colorbar(img, orientation='horizontal',location='top',ax=ax1,extend='max', shrink=0.3, pad=0.01)
cbar.set_label('Sediment thickness (km)',fontsize=15)
# fig.savefig('sedi_crust1.jpg', dpi=400,bbox_inches='tight', pad_inches=0.1)
plt.show()
##
