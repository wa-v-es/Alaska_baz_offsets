# plots output on swat (cvs file) in surface map view.
import numpy as np
import matplotlib.pyplot as plt

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.colors import ListedColormap, BoundaryNorm
import pandas as pd
import sys
from obspy.taup import TauPyModel
from obspy.geodetics.base import gps2dist_azimuth, kilometers2degrees


EARTH_R_KM = 6371.0

def wrap180(a):
    return (a + 180.0) % 360.0 - 180.0

def delta_deg(lat1, lon1, lat2, lon2):
    dist_m, _, _ = gps2dist_azimuth(lat1, lon1, lat2, lon2)
    return kilometers2degrees(dist_m / 1000.0)

def get_rp_for_leg(model, phase, src_depth_km, delta_deg_val,rcv_depth_km):
    rps = model.get_ray_paths(
        source_depth_in_km=float(src_depth_km),
        distance_in_degree=float(delta_deg_val),
        phase_list=[phase],receiver_depth_in_km=rcv_depth_km)
    if not rps:
        return None
    return rps[0]

plt.ion()
csv_path = "/Users/keyser/Research/sct_wat/scattererwhereartthou/examples/swat_230402_180411.csv"
df = pd.read_csv(csv_path)
# single_phase=["P","Ped"]
bounce=["pP","PP"]

# df = df.sample(n=5000, random_state=0)
df["n_bounces"] = df["evt_scat_phase"].isin(bounce).astype(int) + df["sta_scat_phase"].isin(bounce).astype(int)

#####
evt_lat = float(df["evtlat"].iloc[0]);  evt_lon = float(df["evtlon"].iloc[0]);  evt_z = float(df["evtdepth"].iloc[0])
sta_lat = float(df["stalat"].iloc[0]);  sta_lon = float(df["stalon"].iloc[0]);  sta_z = 0.0

model = TauPyModel(model="iasp91")

delta_deg_val= delta_deg(evt_lat, evt_lon, sta_lat, sta_lon)

rp = get_rp_for_leg(model, "P", evt_z, delta_deg_val,0)
rp_PP=get_rp_for_leg(model, "PP", evt_z, delta_deg_val,0)

ray_p_P=np.round(rp.ray_param_sec_degree,3)
GCP_time=np.round(rp.time,3)

lat = df["scatlat"].to_numpy()
lon = df["scatlon"].to_numpy()
depth = df["scatdepth"].to_numpy()

dbaz = wrap180(df["scatbaz"].to_numpy() - df["baz_GCP"].to_numpy())
pval = df["sta_scat_p"].to_numpy()
del_time = df["scat_time"].to_numpy() - float(GCP_time)
N_bounce = df["n_bounces"].to_numpy().astype(int)

#########
proj = ccrs.Robinson(central_longitude=180)
pc = ccrs.PlateCarree()
PACIFIC_EXTENT = (120, 280, -50, 70)

fig, axes = plt.subplots(2, 3, figsize=(20, 13),subplot_kw=dict(projection=proj),constrained_layout=True)
axes = axes.ravel()
def setup_ax(ax, title):
    ax.set_title(title, fontsize=11)
    ax.add_feature(cfeature.LAND, facecolor="honeydew",alpha=.5, zorder=0)
    ax.add_feature(cfeature.OCEAN, facecolor="slategray",alpha=.3, zorder=0)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.1, zorder=1)
    # ax.set_global()
    ax.set_extent(PACIFIC_EXTENT, crs=pc)

#  scatter kwargs
scatter_kw = dict(s=0.4,transform=pc,rasterized=True,alpha=.5)#linewidths=0.1,edgecolors="white",

# 1) Time delay
ax = axes[0]
setup_ax(ax, "Δt (s)")
m = ax.scatter(lon, lat, c=del_time, cmap="viridis", **scatter_kw)#vmin=dt_vmin, vmax=dt_vmax,
plt.colorbar(m, ax=ax, orientation="horizontal", pad=0.03, fraction=0.05)

# 2) Δbaz
ax = axes[1]
setup_ax(ax, "Δbaz (°)")
m = ax.scatter(lon, lat, c=dbaz, cmap="RdBu_r", **scatter_kw) #vmin=baz_vmin, vmax=baz_vmax
plt.colorbar(m, ax=ax, orientation="horizontal", pad=0.03, fraction=0.05)

# 3) Slowness / ray parameter
ax = axes[2]
setup_ax(ax, "Slow at station (s/°)")
m = ax.scatter(lon, lat, c=pval, cmap="cividis", **scatter_kw) #vmin=p_vmin, vmax=p_vmax,
plt.colorbar(m, ax=ax, orientation="horizontal", pad=0.03, fraction=0.05)

# 4) Depth
ax = axes[3]
setup_ax(ax, "Scatterer depth (km)")
m = ax.scatter(lon, lat, c=depth, cmap="plasma", **scatter_kw)#vmin=z_vmin, vmax=z_vmax,
plt.colorbar(m, ax=ax, orientation="horizontal", pad=0.03, fraction=0.05)

# 5) Bounce points (discrete 0/1/2)
ax = axes[4]
setup_ax(ax, "# surface bounces")

bounce_cmap = ListedColormap(["#440154", "#21918c", "#fde725"])  # 0,1,2
bounce_norm = BoundaryNorm(boundaries=[-0.5, 0.5, 1.5, 2.5], ncolors=3)

m = ax.scatter(lon, lat, c=N_bounce, cmap=bounce_cmap, norm=bounce_norm, **scatter_kw)
cb = plt.colorbar(m, ax=ax, orientation="horizontal", pad=0.03, fraction=0.05, ticks=[0, 1, 2])
cb.set_ticklabels(["0", "1", "2"])

plt.show()
plt.savefig("230402_180411_grid_65_-157.png", dpi=700, bbox_inches='tight', pad_inches=0.1)
