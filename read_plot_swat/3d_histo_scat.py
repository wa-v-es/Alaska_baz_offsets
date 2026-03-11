# using swat output, finds the raw count of number of scatterers in each cell in area of interest
# not equal area.
import sys
import numpy as np
import pandas as pd
from obspy.taup import TauPyModel
from obspy.geodetics.base import gps2dist_azimuth, kilometers2degrees
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

csv_path = "/Users/keyser/Research/sct_wat/scattererwhereartthou/examples/swat_230402_180411_all_grids.csv"
df = pd.read_csv(csv_path)
evt_lat = float(df["evtlat"].iloc[0]);  evt_lon = float(df["evtlon"].iloc[0]);  evt_z = float(df["evtdepth"].iloc[0])
sta_lat = float(df["stalat"].iloc[0]);  sta_lon = float(df["stalon"].iloc[0]);  sta_z = 0.0

model = TauPyModel(model="iasp91")
lat = df["scatlat"].to_numpy()
lon = df["scatlon"].to_numpy()
dep = df["scatdepth"].to_numpy()

lon360 = np.mod(lon, 360.0)


###create bins
LAT_MIN, LAT_MAX = -30.0, 50.0
LON_MIN_360, LON_MAX_360 = 90.0, 300.0
DLAT = 0.2
DLON = 0.2

Z_MIN, Z_MAX = 50.0, 2850.0
DZ = 100.0

##
m = ((lat >= LAT_MIN) & (lat <= LAT_MAX) &
    (lon360 >= LON_MIN_360) & (lon360 <= LON_MAX_360) &
    (dep >= Z_MIN) & (dep <= Z_MAX))

lat_f = lat[m]
lon_f = lon360[m]
dep_f = dep[m]

print(f"Kept {lat_f.size:,} / {lat.size:,} scatterers in window")

# -------------------------
# BIN EDGES
# -------------------------
lat_edges = np.arange(LAT_MIN, LAT_MAX + DLAT, DLAT)
lon_edges = np.arange(LON_MIN_360, LON_MAX_360 + DLON, DLON)
# depth bins: 50–2850 every 100 km
dep_edges = np.arange(Z_MIN, Z_MAX + DZ + 1e-9, DZ)

samples = np.column_stack([lat_f, lon_f, dep_f])
counts, edges = np.histogramdd(samples, bins=[lat_edges, lon_edges, dep_edges])

counts = counts.astype(np.uint32)

print("counts shape (nLat, nLon, nDep):", counts.shape)
print("total counted:", int(counts.sum()))
###########
#### CHOOSE DEPTH HERE..
z_target = 550.0
for z_target in np.arange(50,1350,100):
# for z_target in [250]:

    k = np.searchsorted(dep_edges, z_target, side="right") - 1
    if k < 0 or k >= counts.shape[2]:
        raise ValueError("Requested depth is outside dep_edges range.")
    z_lo, z_hi = dep_edges[k], dep_edges[k+1]
    print(f"Plotting bin k={k} covering depth [{z_lo:.0f}, {z_hi:.0f}) km")

    ###
    counts_2d = counts[:, :, k]
    lon_edges_plot = lon_edges - 180.0
    # plot on map
    proj = ccrs.Robinson(central_longitude=180)
    pc_pacific = ccrs.PlateCarree(central_longitude=180)
    plt.ion()
    fig = plt.figure(figsize=(8,8))
    ax = plt.axes(projection=proj)

    ax.add_feature(cfeature.LAND, facecolor="0.92", zorder=0)
    ax.add_feature(cfeature.OCEAN, facecolor="1.0", zorder=0)
    ax.set_extent((-90, 120, -20, 75), crs=pc_pacific)
    ax.scatter(evt_lon,evt_lat,s=30,c='maroon',marker='*',transform=ccrs.PlateCarree(),zorder=20)
    ax.scatter(sta_lon,sta_lat,s=30,c='black',marker='v',transform=ccrs.PlateCarree(),zorder=20)

    # ax.set_extent((LON_MIN_360, LON_MAX_360, LAT_MIN, LAT_MAX), crs=pc_pacific)

    #  plot counts as a gridded field
    # pcolormesh expects edges; counts_2d shape is (len(lat_edges)-1, len(lon_edges)-1)
    counts_plot = np.ma.masked_less(counts_2d, 1)
    vmax = int(counts_2d.max()) if counts_2d.size else 1
    m = ax.pcolormesh(lon_edges_plot, lat_edges, counts_2d,transform=pc_pacific,cmap="GnBu",shading="auto",norm=LogNorm(vmin=1, vmax=vmax))
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5, zorder=10,color='white')

    cb = plt.colorbar(m, ax=ax, orientation="horizontal", pad=0.04, fraction=0.06)
    cb.set_label("Raw scatterer count per 0.2°×0.2°×100 km bin")

    ax.set_title(f"Scatterer density at ~{z_target:.0f} km")
    plt.savefig("figs_samples_depth/230402_180411_{}km.png".format(z_target), dpi=300, bbox_inches='tight', pad_inches=0.1)
    plt.close()

plt.show()
