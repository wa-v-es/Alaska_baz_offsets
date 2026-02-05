####
import numpy as np
import math
import plotly.graph_objects as go

from obspy.taup import TauPyModel
from obspy.geodetics.base import gps2dist_azimuth, kilometers2degrees

EARTH_R_KM = 6371.0

######
src_lat, src_lon, src_depth_km = 10.0, 20.0, 300.0
rcv_lat, rcv_lon = 35.0, 70.0

# -----------------------
# fxs
# -----------------------
def latlon_to_unit(lat_deg, lon_deg):
    lat = math.radians(lat_deg)
    lon = math.radians(lon_deg)
    return np.array([math.cos(lat)*math.cos(lon),
                     math.cos(lat)*math.sin(lon),
                     math.sin(lat)], float)

def rotate_about_axis(v, axis, ang_rad):
    axis = axis / np.linalg.norm(axis)
    c = math.cos(ang_rad)
    s = math.sin(ang_rad)
    return v*c + np.cross(axis, v)*s + axis*np.dot(axis, v)*(1-c)

def ecef_from_latlon_depth(lat_deg, lon_deg, depth_km):
    r = EARTH_R_KM - depth_km
    lat = math.radians(lat_deg)
    lon = math.radians(lon_deg)
    x = r * math.cos(lat) * math.cos(lon)
    y = r * math.cos(lat) * math.sin(lon)
    z = r * math.sin(lat)
    return x, y, z

# -----------------------
# 1) Read scatterers
# -----------------------
fname = "tube_scatter_candidates.txt"
# file columns: scat_lat scat_lon scat_depth_km d1_deg dt_s dP_sdeg dbaz_deg
data = np.loadtxt(fname, skiprows=4)

scat_lat = data[:, 0]
scat_lon = data[:, 1]
scat_dep = data[:, 2]
scat_dts = data[:, 4]
scat_slow = data[:, 5]
scat_baz = data[:, 6]


# Convert to XYZ (km)
sx, sy, sz = [], [], []
for la, lo, de in zip(scat_lat, scat_lon, scat_dep):
    x, y, z = ecef_from_latlon_depth(la, lo, de)
    sx.append(x); sy.append(y); sz.append(z)
sx = np.array(sx); sy = np.array(sy); sz = np.array(sz)

# -----------------------
# 2) Compute direct P ray path (TauP) and convert to XYZ
# -----------------------
model = TauPyModel(model="iasp91")

dist_m, az_sr, baz_rs = gps2dist_azimuth(src_lat, src_lon, rcv_lat, rcv_lon)
delta_sr_deg = kilometers2degrees(dist_m / 1000.0)

rps = model.get_ray_paths(
    source_depth_in_km=src_depth_km,
    distance_in_degree=delta_sr_deg,
    phase_list=["P"])
if not rps:
    raise RuntimeError("No direct P ray path returned by TauP for this geometry.")
rp = rps[0]

# TauP raypath distances are in radians -> degrees
path_dist_deg = np.degrees(np.asarray(rp.path["dist"], float))
path_depth_km = np.asarray(rp.path["depth"], float)

# Build great-circle rotation axis (spherical)
s_hat = latlon_to_unit(src_lat, src_lon)
r_hat = latlon_to_unit(rcv_lat, rcv_lon)
axis = np.cross(s_hat, r_hat)
axis_norm = np.linalg.norm(axis)
if axis_norm < 1e-12:
    raise RuntimeError("Source/receiver geometry degenerate (axis ill-defined).")
axis = axis / axis_norm

# Convert ray path (dist, depth) to XYZ by rotating along great-circle
rx, ry, rz = [], [], []
for ddeg, zkm in zip(path_dist_deg, path_depth_km):
    u_surf = rotate_about_axis(s_hat, axis, math.radians(float(ddeg)))
    r = EARTH_R_KM - float(zkm)
    p = r * u_surf
    rx.append(p[0]); ry.append(p[1]); rz.append(p[2])
rx = np.array(rx); ry = np.array(ry); rz = np.array(rz)

# Source and receiver points
src_x, src_y, src_z = ecef_from_latlon_depth(src_lat, src_lon, src_depth_km)
rcv_x, rcv_y, rcv_z = ecef_from_latlon_depth(rcv_lat, rcv_lon, 0.0)

# -----------------------
# 3) Plotly Scatter3d
# -----------------------
fig = go.Figure()

# Scatterers
fig.add_trace(go.Scatter3d(
    x=sx, y=sy, z=sz,
    mode="markers",
    marker=dict(
        size=2.5,
        color=scat_dts,           # color by depth
        colorscale="Viridis",
        colorbar=dict(title="Time (s)",x=-0.12,
        xanchor="left",
        len=0.75),
        opacity=0.7,),
    name="wrt time",
    text=[f"lat={la:.3f}, lon={lo:.3f}, z={de:.1f} km" for la, lo, de in zip(scat_lat, scat_lon, scat_dep)],
    hoverinfo="text"))

fig.add_trace(go.Scatter3d(
    x=sx, y=sy, z=sz,
    mode="markers",
    marker=dict(
        size=2.5,
        color=scat_baz,           # color by depth
        colorscale="PRGn",
        colorbar=dict(title="Baz (deg)",x=-0.12,
        xanchor="left",
        len=0.75),
        opacity=0.7,),
    name="wrt baz",
    text=[f"lat={la:.3f}, lon={lo:.3f}, z={de:.1f} km" for la, lo, de in zip(scat_lat, scat_lon, scat_dep)],
    hoverinfo="text"))

fig.add_trace(go.Scatter3d(
    x=sx, y=sy, z=sz,
    mode="markers",
    marker=dict(
        size=2.5,
        color=scat_slow,           # color by depth
        colorscale="cividis",
        colorbar=dict(title="Slow (s/deg)",x=-0.12,
        xanchor="left",
        len=0.75),
        opacity=0.7,),
    name="wrt slow",
    text=[f"lat={la:.3f}, lon={lo:.3f}, z={de:.1f} km" for la, lo, de in zip(scat_lat, scat_lon, scat_dep)],
    hoverinfo="text"))

# Ray path
fig.add_trace(go.Scatter3d(
    x=rx, y=ry, z=rz,
    mode="lines",
    line=dict(width=5),
    name="Direct P ray path"))

# Source / Receiver markers
fig.add_trace(go.Scatter3d(
    x=[src_x], y=[src_y], z=[src_z],
    mode="markers+text",
    marker=dict(size=6),
    text=["SRC"],
    textposition="top center",
    name="Source"))

fig.add_trace(go.Scatter3d(
    x=[rcv_x], y=[rcv_y], z=[rcv_z],
    mode="markers+text",
    marker=dict(size=6),
    text=["RCV"],
    textposition="top center",
    name="Receiver"))

fig.update_layout(
    title=f"Tube scatterers (N={len(scat_dep)}) + direct P ray (iasp91) | Δ≈{delta_sr_deg:.2f}°",
    scene=dict(
        xaxis_title="X (km)",
        yaxis_title="Y (km)",
        zaxis_title="Z (km)",
        aspectmode="data",
    ),legend=dict(itemsizing="constant"))

fig.show()
