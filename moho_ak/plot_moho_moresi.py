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

import os.path as path
import stripy as stripy
import sys
#######

moho_point='Models/AlaskaMoho.npz'
AlaskaMoho_point=np.load(moho_point,allow_pickle = True)
AlaskaMoho_point['alaska_moho']

moho_finegrid='Models/AlaskaMohoErrs-AlaskaMohoFineGrid.npz'
# AlaskaMoho_finegrid=np.load(moho_finegrid,allow_pickle = True)
# AlaskaMoho_finegrid['alaska_moho']

# print(AlaskaMoho_finegrid.files)
mohoraw = alaskamoho.MohoErr
msmoho_opt  = alaskamoho.MohoModel_opt
msmoho_min  = alaskamoho.MohoModel_min
msmoho_minj = alaskamoho.MohoModel_minj

# sys.exit()
####
filename="AlaskaMohoOpt1pct_pt.png"

model=msmoho_opt
print(msmoho_opt.description)
print('----------------------')


plot_range=[25, 45]
show_coastline=True
show_bg_image=True
raw_data_points=alaskamoho.MohoErr
cmap=None
####
goodgrid = model.gridF
quality = model.quality
grid_data = model.surface

try:
    import gdal
    globalsrelief       = gdal.Open("ShadedRelief/GRAY_HR_SR_OB.tif")
    globalsrelief_img   = globalsrelief.ReadAsArray()/255.0  # .transpose(1,2,0)
    globalsrelief_img_q = globalsrelief_img[0:globalsrelief_img.shape[0]//4, 0:globalsrelief_img.shape[1]//4]

except ImportError:
    show_bg_image = False

if cmap == None:
    cmap = plt.cm.RdYlBu

# Transparent colours
from matplotlib.colors import ListedColormap
###
colA = cmap(np.arange(cmap.N))
colA[:,-1] = 0.25 + 0.5 * np.linspace(-1.0, 1.0, cmap.N)**2.0
#adjusts the opacity based on the quadratic curve, setting values between 0.25 and 0.75.
##

# Create new colormap
cmapA = ListedColormap(colA)
# cmapA = cmap
proj = ccrs.Stereographic(central_longitude=-154, central_latitude=90, true_scale_latitude=60)
# plt.clf()
plt.ion()
fig = plt.figure(figsize=(15, 8), facecolor=None)
ax1 = plt.subplot(1, 1, 1, projection=proj)


# ax1.gridlines(draw_labels=False, linewidth=0.5, color='gray', alpha=0.5)

ax1.set_extent([-165,-138,55,70.5], crs=ccrs.PlateCarree())

# grat = cartopy.feature.NaturalEarthFeature(category="physical", scale="10m", name="graticules_5")
# ax1.add_feature(grat, linewidth=0.5,linestyle="--",edgecolor="#000000",facecolor="None", zorder=2)

ax1.coastlines(resolution="10m",color="#111111", linewidth=0.5, zorder=99)
ax1.add_feature(cfeature.BORDERS.with_scale('10m'), linestyle=':')
ax1.add_feature(cfeature.OCEAN.with_scale('10m'),alpha=0.2,facecolor='xkcd:dusty blue')

gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=.8, color='gray', alpha=0.5, linestyle='--',rotate_labels=False,
        x_inline=False, y_inline=False)

gl.xlocator = mticker.FixedLocator([-160,-150,-140])
gl.ylocator = mticker.FixedLocator([55,60,65,70])

# gl.xlines = True
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 14}
gl.ylabel_style = {'size': 14}

# sys.exit()

lons = np.degrees(goodgrid.lons)%360.0 #
lats = np.degrees(goodgrid.lats)

gdata2 = grid_data.copy()
gdata2[ quality == 0.0 ] = -10000000

cnt0=ax1.tricontourf(lons, lats, goodgrid.simplices, gdata2,
               cmap=cmapA,levels=np.linspace(plot_range[0], plot_range[1], 11),
               extend="max",transform=ccrs.PlateCarree(),zorder=10)

gdata2 = grid_data.copy()
gdata2[ quality < 0.05 ] = -10000000

cnt0=ax1.tricontourf(lons, lats, goodgrid.simplices, gdata2,cmap=cmapA,
               levels=np.linspace(plot_range[0], plot_range[1], 11),extend="max",
               transform=ccrs.PlateCarree(),
                     # alpha=0.5,
                     zorder=11)
# raw_data_points= None
if raw_data_points is not None:

    m = ax1.scatter(raw_data_points['lon'], raw_data_points['lat'],  color="Black",
                   edgecolor="Black", linewidth=0.5,
                   marker="+", s=25.0, transform=ccrs.Geodetic(),alpha=.75, zorder=25)
###


cbar=plt.colorbar(ax=ax1, mappable=cnt0, shrink=0.5, extend='max', drawedges=False, pad=0.02 )
cbar.set_label("Moho depth (km)",fontsize=13)
# fig.savefig(filename, dpi=600,bbox_inches='tight', pad_inches=0.1)
plt.show()
##
