###
import numpy as np
# from numpy import load
import obspy
import miller_alaskamoho_srl2018 as alaskamoho
import cartopy
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

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

sys.exit()
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

fig = plt.figure(figsize=(15, 8), facecolor=None)
ax1  = plt.subplot(111, projection=ccrs.AlbersEqualArea(central_longitude=-154, central_latitude=50,
                                                   standard_parallels=(55,65) ))
ax1.set_extent([-170,-130,50,70])

grat = cartopy.feature.NaturalEarthFeature(category="physical", scale="50m", name="graticules_5")
ax1.add_feature(grat, linewidth=0.5,linestyle="--",edgecolor="#000000",facecolor="None", zorder=2)

if(show_bg_image):
    ax1.imshow(globalsrelief_img_q**0.5, origin='upper', transform=ccrs.PlateCarree(),
               extent=[-180,-90,45,90], zorder=0, alpha=1.0, cmap=plt.cm.Greys_r, vmin=0.0, vmax=1.0)

ax1.coastlines(resolution="50m",color="#111111", linewidth=0.5, zorder=99)


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

if raw_data_points is not None:

    m = ax1.scatter(raw_data_points['lon'], raw_data_points['lat'],  color="Black",
                   edgecolor="Black", linewidth=0.5,
                   marker="+", s=25.0, transform=ccrs.Geodetic(), zorder=25)


plt.colorbar(ax=ax1, mappable=cnt0, shrink=0.5, extend='max', drawedges=False )

# fig.savefig(filename, dpi=600,bbox_inches='tight', pad_inches=0.1)
plt.show()
##
