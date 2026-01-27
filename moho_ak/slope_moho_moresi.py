
import obspy
import miller_alaskamoho_srl2018 as alaskamoho
import os.path as path
import stripy as stripy
import numpy as np

import cartopy
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from shapely.geometry import Point
import shapely.geometry as sgeom
from shapely.geometry import Polygon
from cartopy.feature import ShapelyFeature
from cartopy.io.shapereader import Reader
from cartopy.io import shapereader
# import shapely.geometry as sgeom
from matplotlib.path import Path
from matplotlib.patches import PathPatch
import glob as glob
from cmcrameri import cm
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.pyplot as plt
import matplotlib.path as mpath
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import sys
import pygmt
import xarray as xr
import matplotlib.colors as colors

####
#
def plot_baz_circle(ax1,file,marker):
    for line in open(file,'r'):
        line=line.split()
        # if float(line[2])-1.5*float(line[3]) < float(line[1]) < float(line[2])+1.5*float(line[3]) and float(line[3]) < 3.5:
        if float(line[3]) < 3:

            if -.6 < float(line[1]) < .6:
                ax1.scatter(float(line[5]),float(line[4]),s=abs(float(line[1]))*165, marker='o', facecolor='none', edgecolor='black',alpha=.25,linewidth=.75,zorder=10,transform=ccrs.PlateCarree())

            if float(line[1]) > .6:
                ax1.scatter(float(line[5]),float(line[4]),s=abs(float(line[1]))*165, marker='o', facecolor='rebeccapurple', alpha=.125,zorder=10,transform=ccrs.PlateCarree())
                ax1.scatter(float(line[5]),float(line[4]),s=abs(float(line[1]))*165, marker='o', facecolor='None', edgecolor='white',alpha=.95,linewidth=.35,zorder=10,transform=ccrs.PlateCarree())

            if float(line[1]) < -.6:
                ax1.scatter(float(line[5]),float(line[4]),s=abs(float(line[1]))*165, marker='o', facecolor='seagreen', alpha=.125,zorder=20,transform=ccrs.PlateCarree())
                ax1.scatter(float(line[5]),float(line[4]),s=abs(float(line[1]))*165, marker='o', facecolor='None', edgecolor='white',alpha=.95,linewidth=.35,zorder=10,transform=ccrs.PlateCarree())

        else:
            print('Grid:{} max not reprensentative of vespa'.format(float(line[0])))

#
def plot_baz_circle_mean(ax1,file,marker):
    for line in open(file,'r'):
        line=line.split()
        # if float(line[2])-1.5*float(line[3]) < float(line[1]) < float(line[2])+1.5*float(line[3]) and float(line[3]) < 3.5:
        if float(line[3]) < 3:

            if -.6 < float(line[2]) < .6:
                ax1.scatter(float(line[5]),float(line[4]),s=abs(float(line[2]))*265, marker='o', facecolor='none', edgecolor='black',alpha=.75,linewidth=.75,zorder=10,transform=ccrs.PlateCarree())

            if float(line[2]) > .6:
                ax1.scatter(float(line[5]),float(line[4]),s=abs(float(line[2]))*165, marker='o', facecolor='rebeccapurple', alpha=.125,zorder=10,transform=ccrs.PlateCarree())
                ax1.scatter(float(line[5]),float(line[4]),s=abs(float(line[2]))*165, marker='o', facecolor='None', edgecolor='white',alpha=.95,linewidth=.25,zorder=10,transform=ccrs.PlateCarree())

            if float(line[2]) < -.6:
                ax1.scatter(float(line[5]),float(line[4]),s=abs(float(line[2]))*165, marker='o', facecolor='seagreen', alpha=.125,zorder=20,transform=ccrs.PlateCarree())
                ax1.scatter(float(line[5]),float(line[4]),s=abs(float(line[2]))*165, marker='o', facecolor='None', edgecolor='white',alpha=.95,linewidth=.25,zorder=10,transform=ccrs.PlateCarree())

        else:
            print('Grid:{} max not reprensentative of vespa'.format(float(line[0])))

####
def add_topo(ax,exx=1.5):

    # Load ETOPO1 topography data
    topogrd_path = '/Users/keyser/Documents/mirage/etopo/ETOPO1_Ice_g_gmt5.grd'
    ds = xr.open_dataset(topogrd_path)

    # Extract variables (coordinate names might be 'x', 'y', or 'lon', 'lat')
    elevation = ds['z']
    lon = ds['x']  # or ds['lon'] depending on file
    lat = ds['y']  # or ds['lat']

    # Subset the data for Alaska region (make sure lat slicing is in descending order!)
    lon_min, lon_max = -165, -135
    lat_min, lat_max = 55, 70
    subset = elevation.sel(x=slice(lon_min, lon_max), y=slice(lat_max, lat_min))  # y descending

    # Optional exaggeration
    exaggeration_factor = exx  # Set your value, e.g., 10
    topo_exaggerated = subset * exaggeration_factor

    # Meshgrid for pcolormesh (must match data dims: [y,x])
    lon2d, lat2d = np.meshgrid(subset['x'], subset['y'])
    cmap = plt.get_cmap('Greys')
    # norm = colors.Normalize(vmin=float(topo_exaggerated.min()), vmax=float(topo_exaggerated.max()))
    img = ax.pcolormesh(lon2d, lat2d, topo_exaggerated, cmap=cmap, transform=ccrs.PlateCarree())#norm=norm,


msmoho_opt  = alaskamoho.MohoModel_opt
#MohoModel_opt is the ensemble of models that lie within 1.25% of the model that best predicts the test data
# and 1.25% of the model that best predicts all the data
model=msmoho_opt
plot_range=[25, 45]
show_bg_image=False
show_coastline=True
show_bg_image=True
raw_data_points=None
cmap=None

goodgrid  = model.gridF
grid_data = model.surface
grad_data = model.surface_slope
quality   = model.quality

for_TA=False
if cmap is None:
    cmap = plt.cm.Greys

# Transparent colours
from matplotlib.colors import ListedColormap

colA = cmap(np.arange(cmap.N))
colA[:,-1] = np.linspace(0.0, 1, cmap.N) #replaces the last column (alpha channel) in colA with values from 0.0 to 1, creating a gradient in transparency across the color map.

# Create new colormap
cmapA = ListedColormap(colA)

# try:
#     from osgeo import gdal
#     globalsrelief       = gdal.Open("ShadedRelief/GRAY_HR_SR_OB.tif")
#     globalsrelief_img   = globalsrelief.ReadAsArray()/255.0  # .transpose(1,2,0)
#     globalsrelief_img_q = globalsrelief_img[0:globalsrelief_img.shape[0]//4, 0:globalsrelief_img.shape[1]//4]
#
# except ImportError:
#     show_bg_image = False
proj = ccrs.Stereographic(central_longitude=-154, central_latitude=90, true_scale_latitude=60)
if for_TA:
    proj = ccrs.Stereographic(central_longitude=-90, central_latitude=90, true_scale_latitude=37)

plt.rcParams.update({'font.size': 19})

plt.ion()
fig = plt.figure(figsize=(15, 9), facecolor=None)
ax1 = plt.subplot(1, 1, 1, projection=proj)


# ax1.gridlines(draw_labels=False, linewidth=0.5, color='gray', alpha=0.5)

# ax1  = plt.subplot(111, projection=ccrs.LambertConformal(central_longitude=-154, central_latitude=50))#,
                                                   # standard_parallels=(55,65) ))
# ax1.spines['geo'].set_visible(False)
ax1.set_extent([-165,-138,55,70.5], crs=ccrs.PlateCarree())
if for_TA:
    ax1.set_extent([-95,-85,30,47.5], crs=ccrs.PlateCarree())

# if(show_bg_image):
#     ax1.imshow(globalsrelief_img_q**0.5, origin='upper', transform=ccrs.PlateCarree(),
#        extent=[-180,-90,45,90], zorder=0, alpha=1.0, cmap=plt.cm.Greys_r, vmin=0.0, vmax=1.0)

# add_topo(ax1,1.5)

ax1.coastlines(resolution="10m",color="#111111", linewidth=0.55)

ax1.add_feature(cfeature.OCEAN.with_scale('10m'),alpha=0.2,facecolor='xkcd:dusty blue')
ax1.add_feature(cfeature.BORDERS.with_scale('10m'), linestyle=':')
# ax1.add_feature(cfeature.LAKES.with_scale('50m'), alpha=0.5,facecolor='xkcd:teal blue')
# ax1.add_feature(cfeature.RIVERS.with_scale('10m'))

gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=.8, color='gray', alpha=0.5, linestyle='--',rotate_labels=False,
        x_inline=False, y_inline=False)
gl.xlocator = mticker.FixedLocator([ -160, -150, -140])
gl.ylocator = mticker.FixedLocator([55,60,65,70])
if for_TA:
    gl.xlocator = mticker.FixedLocator([-95,-90,-85,-80])
    gl.ylocator = mticker.FixedLocator([30,35,40,45])
# gl.xlines = True
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 19}
gl.ylabel_style = {'size': 19}

# plt.show()
# sys.exit()

lons = np.degrees(goodgrid.lons)%360.0
lats = np.degrees(goodgrid.lats)

gdata2 = grad_data.copy()
gdata2[ quality == 0] = -1000000.0


# grad_grd=pygmt.load_dataarray('grad_.25_posit_mask.grd')
# grad_grd.plot(ax=ax1,cmap=plt.cm.Greys,add_colorbar=True,vmin=0,mouseover=True,transform=ccrs.PlateCarree()) #plt.cm.Greys

# cnt0=ax1.tricontourf(lons, lats, goodgrid.simplices, gdata2,
#                cmap=cmapA,
#                levels=np.linspace(0.0,0.1,11),
#                extend="max",
#                transform=ccrs.PlateCarree(),
#                      zorder=1)

centroids_q, err = goodgrid.interpolate(goodgrid.lons[goodgrid.simplices].mean(axis=1),
                                     goodgrid.lats[goodgrid.simplices].mean(axis=1), quality)


mask = centroids_q == 0.0

# cnt=ax1.tricontour(lons, lats, goodgrid.simplices, grid_data,
#                alpha=1.0, colors="#000000", mask=mask,
#                levels=[25, 29, 33, 37, 41, 45],
#                linewidths=[0.5,0.5, 0.75, 0.75, 1.0, 1.0],
#                linestyles=[":","-"],
#                transform=ccrs.PlateCarree(),
#                    zorder=5)

# plt.colorbar(ax=ax1, mappable=cnt0, shrink=0.4,pad=.012, extend='max', drawedges=False)

## Locations of raw data collection

if raw_data_points is not None:

    m = ax1.scatter(raw_data_points['lon'], raw_data_points['lat'],  color="Black",
                   edgecolor="Black", linewidth=0.5,
                   marker="+", s=25.0, transform=ccrs.Geodetic(), zorder=25)


# fig.savefig('grad_.25_posit_mask_60.png', dpi=400,bbox_inches='tight', pad_inches=0.1)

# plt.show()
# sys.exit()
#plot individual POINTS
# 230304_064122_PA_inc2_r2.5
folder_pattern_pa = "/Users/keyser/Research/AK_all_stations/sac_files/max_vals_coherence/*maxVals_low_slow.txt"

if for_TA:
    folder_pattern_pa = "/Users/keyser/Research/TA_arrays/sac_files/TA_max_vals_coherence/*maxVals_low_slow.txt"

# folder_pattern_pa = "/Users/keyser/Research/AK_all_stations/sac_files_with_P/max_vals_coherence/*PA_maxVals_low_slow.txt"
# folder_pattern_sa = "/Users/keyser/Research/AK_all_stations/sac_files/max_vals_coherence/*SA_maxVals_low_slow.txt"

matching_files_pa = glob.glob(folder_pattern_pa)
# matching_files_sa = glob.glob(folder_pattern_sa)
#plots individual max baz offsets
#
# for file in matching_files_pa:
#     plot_baz_circle(ax1,file,'o')

# plots individual mean baz offsets
for file in matching_files_pa:
    plot_baz_circle_mean(ax1,file,'o')

# plt.show()
mean_max='mean'
low_high='low'
filename='MohoSlope_indi_{}_{}_slow_PA.png'.format(mean_max,low_high)
##
legend_elements = [Line2D([0], [0], color='white',marker='o', lw=4, label='1$^\circ$ offset',
                    markerfacecolor='rebeccapurple',markersize=20,alpha=.75),
                   Line2D([0], [1], marker='o', color='white', label='-2$^\circ$ offset',
                          markerfacecolor='seagreen', markersize=30,alpha=.75),]
                #   Patch(facecolor='black', edgecolor='k',alpha=.5,
                #         label='Moho grad.')]
# fig.subplots_adjust(top=0.8)
#
if for_TA:
    loca='lower right'
else:
    loca='upper right'
ax1.legend(handles=legend_elements,loc=loca, labelspacing=1.5,handletextpad=1.5,
    borderaxespad=.75,bbox_transform=fig.transFigure,fontsize=17 )
plt.show()
fig.savefig('without_P_offset_mean_PA_new.png', dpi=400,bbox_inches='tight', pad_inches=0.1)
#without_P_offset_mean_PA_new.png
sys.exit()

### add baz offset block mean
mean_max='mean'
low_high='low'
filename='MohoSlope_block_{}_{}_slow_PA_new.png'.format(mean_max,low_high)
block_file='/Users/keyser/Research/AK_all_stations/baz_offset_plots/block_mean_vals/block_{}_{}_slow_PA.txt'.format(mean_max,low_high)
for line in open(block_file,'r'):
    line=line.split()
    if -.1 < float(line[2]) < .1:
        # print(' ')
        ax1.scatter(float(line[0]),float(line[1]),s=165, marker='o', facecolor='none', edgecolor='black',alpha=.75,linewidth=.75,zorder=10,transform=ccrs.PlateCarree())
    elif float(line[2]) > .1:
        ax1.scatter(float(line[0]),float(line[1]),s=abs(float(line[2]))*165, marker='o', facecolor='rebeccapurple', alpha=.25,zorder=10,transform=ccrs.PlateCarree())
        ax1.scatter(float(line[0]),float(line[1]),s=abs(float(line[2]))*165, marker='o', facecolor='None', edgecolor='white',alpha=.95,linewidth=1.25,zorder=10,transform=ccrs.PlateCarree())

        ax1.scatter(float(line[0]),float(line[1]),s=abs(float(line[3]))*165, marker='|', facecolor='rebeccapurple',alpha=.75,zorder=10,linewidth=1.75,transform=ccrs.PlateCarree())

    elif float(line[2]) < -.1:
        ax1.scatter(float(line[0]),float(line[1]),s=abs(float(line[2]))*165, marker='o', facecolor='seagreen', alpha=.25,zorder=10,transform=ccrs.PlateCarree())
        ax1.scatter(float(line[0]),float(line[1]),s=abs(float(line[2]))*165, marker='o', facecolor='None', edgecolor='white',alpha=.95,linewidth=1.25,zorder=10,transform=ccrs.PlateCarree())

        ax1.scatter(float(line[0]),float(line[1]),s=abs(float(line[3]))*165, marker='|', facecolor='seagreen',alpha=.75,zorder=10,linewidth=1.75,transform=ccrs.PlateCarree())


plt.show()
# fig.savefig(filename, dpi=400,bbox_inches='tight', pad_inches=0.1)
#######
