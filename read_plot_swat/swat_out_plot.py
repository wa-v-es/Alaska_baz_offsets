# reads output in csv from swat and makes 3d plots using plotly
# Feb 25 2026
import numpy as np
import math
import plotly.graph_objects as go
import sys
from obspy.taup import TauPyModel
from obspy.geodetics.base import gps2dist_azimuth, kilometers2degrees
import pandas as pd
import sys

EARTH_R_KM = 6371.0

def wrap180(a):
    return (a + 180.0) % 360.0 - 180.0

def ecef_from_latlon_depth(lat_deg, lon_deg, depth_km):
    r = EARTH_R_KM - depth_km
    lat = np.deg2rad(lat_deg)
    lon = np.deg2rad(lon_deg)
    x = r * np.cos(lat) * np.cos(lon)
    y = r * np.cos(lat) * np.sin(lon)
    z = r * np.sin(lat)
    return x, y, z
def latlon_to_unit(lat_deg, lon_deg):
    lat = math.radians(lat_deg); lon = math.radians(lon_deg)
    return np.array([math.cos(lat)*math.cos(lon),
                     math.cos(lat)*math.sin(lon),
                     math.sin(lat)], float)

def rotate_about_axis(v, axis, ang_rad):
    axis = axis / np.linalg.norm(axis)
    c = math.cos(ang_rad); s = math.sin(ang_rad)
    return v*c + np.cross(axis, v)*s + axis*np.dot(axis, v)*(1-c)
#
def get_xyz_gcp_ray(rp,evt_lat, evt_lon,sta_lat, sta_lon):
    path_dist_deg = np.degrees(np.asarray(rp.path["dist"], float))
    path_depth_km = np.asarray(rp.path["depth"], float)

    # Great-circle rotation axis (spherical)
    s_hat = latlon_to_unit(evt_lat, evt_lon)
    r_hat = latlon_to_unit(sta_lat, sta_lon)
    axis = np.cross(s_hat, r_hat)
    axis = axis / np.linalg.norm(axis)

    rx, ry, rz = [], [], []
    for ddeg, zkm in zip(path_dist_deg, path_depth_km):
        u_surf = rotate_about_axis(s_hat, axis, math.radians(float(ddeg)))
        r = EARTH_R_KM - float(zkm)
        p = r * u_surf
        rx.append(p[0]); ry.append(p[1]); rz.append(p[2])

    return rx,ry,rz

csv_path = "/Users/keyser/Research/sct_wat/scattererwhereartthou/examples/swat_try.csv"
df = pd.read_csv(csv_path)

###
evt_lat = float(df["evtlat"].iloc[0]);  evt_lon = float(df["evtlon"].iloc[0]);  evt_z = float(df["evtdepth"].iloc[0])
sta_lat = float(df["stalat"].iloc[0]);  sta_lon = float(df["stalon"].iloc[0]);  sta_z = 0.0

evt_x, evt_y, evt_z3 = ecef_from_latlon_depth(evt_lat, evt_lon, evt_z)
sta_x, sta_y, sta_z3 = ecef_from_latlon_depth(sta_lat, sta_lon, sta_z)

model = TauPyModel(model="iasp91")

dist_m, az_sr, baz_rs = gps2dist_azimuth(evt_lat, evt_lon, sta_lat, sta_lon)
delta_sr_deg = kilometers2degrees(dist_m / 1000.0)

rps = model.get_ray_paths(source_depth_in_km=evt_z, distance_in_degree=delta_sr_deg, phase_list=["P","PP"])
rp = rps[0]
rp_PP=rps[1]
ray_p_P=np.round(rp.ray_param_sec_degree,3)
GCP_time=np.round(rp.time,3)
# sys.exit()
# TauP raypath distances are in radians -> degrees
rx,ry,rz=get_xyz_gcp_ray(rp,evt_lat, evt_lon,sta_lat, sta_lon)
rx_PP,ry_PP,rz_PP=get_xyz_gcp_ray(rp_PP,evt_lat, evt_lon,sta_lat, sta_lon)

# Compute XYZ for scatterers
x, y, z = ecef_from_latlon_depth(df["scatlat"].to_numpy(),
                                 df["scatlon"].to_numpy(),
                                 df["scatdepth"].to_numpy())

# ["scatlat", "scatlon", "scatdepth",
# "scatdistdeg", "scatbaz", "sta_scat_p", "scat_time",
# 'sta_scat_phase','evt_scat_phase'
# "evtlat", "evtlon", "evtdepth",
# "stalat", "stalon",'baz_GCP']

# color by: del baz, absolute P, and del time..

dbaz = wrap180(df["scatbaz"].to_numpy() - df["baz_GCP"].to_numpy())
pval = df["sta_scat_p"].to_numpy()
del_time=df["scat_time"].to_numpy() - GCP_time

# Hover text
hover = (
    "scatlat=" + df["scatlat"].round(3).astype(str) +
    ", scatlon=" + df["scatlon"].round(3).astype(str) +
    ", z=" + df["scatdepth"].round(1).astype(str) + " km" +
    "<br>dbaz=" + np.round(dbaz, 2).astype(str) + " deg" +
    "<br>p=" + df["sta_scat_p"].round(3).astype(str) +
    "<br>dt? scat_time=" + df["scat_time"].round(2).astype(str) +
    "<br>phases: " + df["evt_scat_phase"].astype(str) + " / " + df["sta_scat_phase"].astype(str))

fig = go.Figure()
fig.add_trace(go.Scatter3d(
    x=rx, y=ry, z=rz,
    mode="lines",
    line=dict(width=4,color='cyan'),
    name="Direct P"
))
fig.add_trace(go.Scatter3d(
    x=rx_PP, y=ry_PP, z=rz_PP,
    mode="lines",
    line=dict(width=4,color='magenta'),
    name="PP ray"
))
fig.add_trace(go.Scatter3d(
    x=[evt_x], y=[evt_y], z=[evt_z3],
    mode="markers+text",
    marker=dict(size=8, symbol="cross"),
    text=["EVT"],
    textposition="top center",
    name="Event"
))

fig.add_trace(go.Scatter3d(
    x=[sta_x], y=[sta_y], z=[sta_z3],
    mode="markers+text",
    marker=dict(size=9, symbol="square"),
    text=["STA"],
    textposition="top center",
    name="Station"
))
# color by dbaz
fig.add_trace(go.Scatter3d(
    x=x, y=y, z=z,
    mode="markers",
    marker=dict(
        size=2.5,
        color=dbaz,
        colorscale="RdBu",
        colorbar=dict(title="Δ baz (°)", x=-0.12, xanchor="left", len=0.75),
        opacity=0.7,
    ),
    name="Color: dbaz",
    text=hover,
    hoverinfo="text",
    visible=True
))

# color by sta_scat_p
fig.add_trace(go.Scatter3d(
    x=x, y=y, z=z,
    mode="markers",
    marker=dict(
        size=2.5,
        color=pval,
        colorscale="Viridis",
        colorbar=dict(title="Slow at station (s/°)", x=-0.12, xanchor="left", len=0.75),
        opacity=0.7,
    ),
    name="Color: sta_scat_p",
    text=hover,
    hoverinfo="text",
    visible=True
))
# color by time
fig.add_trace(go.Scatter3d(
    x=x, y=y, z=z,
    mode="markers",
    marker=dict(
        size=2.5,
        color=pval,
        colorscale="plasma",
        colorbar=dict(title="Δ time (s)", x=-0.12, xanchor="left", len=0.75),
        opacity=0.7,
    ),
    name="Color: dt",
    text=hover,
    hoverinfo="text",
    visible=True
))

fig.update_layout(
    title=f"Scatterers between sP & PP. P slow={ray_p_P:.2f}",
    # scene controls 3D axes, aspect, camera, etc.
    scene=dict(xaxis_title="X (km)", yaxis_title="Y (km)", zaxis_title="Z (km)",
        aspectmode="data"),
    margin=dict(l=120, r=180, t=60, b=40),
    legend=dict(x=1.02, y=1.0, xanchor="left", yanchor="top"),
    updatemenus=[dict(
        type="dropdown",
        x=1.02, y=0.8,
        buttons=[
            dict(label="Color: scatbaz - baz_GCP",
                 method="update",
                 args=[{"visible": [True,True, True, True, True, False, False]},
                       {"title": "Scatterers in 3D (color = scatbaz - baz_GCP)"}]),
            dict(label="Color: sta_scat_p",
                 method="update",
                 args=[{"visible": [True,True, True, True, False, True, False]},
                       {"title": "Scatterers in 3D (color = sta_scat_p)"}]),
            dict(label="Color: dt",
                 method="update",
                 args=[{"visible": [True,True, True, True, False, False, True]},
                       {"title": "Scatterers in 3D (color = dt)"}]),
        ],
    )]
)

fig.show()
