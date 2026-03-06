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

def delta_deg(lat1, lon1, lat2, lon2):
    dist_m, _, _ = gps2dist_azimuth(lat1, lon1, lat2, lon2)
    return kilometers2degrees(dist_m / 1000.0)

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

def get_rp_for_leg(model, phase, src_depth_km, delta_deg_val,rcv_depth_km):
    rps = model.get_ray_paths(
        source_depth_in_km=float(src_depth_km),
        distance_in_degree=float(delta_deg_val),
        phase_list=[phase],receiver_depth_in_km=rcv_depth_km)
    if not rps:
        return None
    return rps[0]

def raypaths_for_row(model, row):

    evt_lat, evt_lon, evt_z = float(row.evtlat), float(row.evtlon), float(row.evtdepth)
    sta_lat, sta_lon = float(row.stalat), float(row.stalon)
    scat_lat, scat_lon, scat_z = float(row.scatlat), float(row.scatlon), float(row.scatdepth)

    phase1 = str(row.evt_scat_phase).strip()  # evt to scat
    phase2 = str(row.sta_scat_phase).strip()  # station to scat

    # PHASE_MAP = {"Ped": "P"}  # changes all ped to P
    # phase1 = PHASE_MAP.get(phase1, phase1)
    # phase2 = PHASE_MAP.get(phase2, phase2)

    # distances
    d_evt_scat = delta_deg(evt_lat, evt_lon, scat_lat, scat_lon)
    d_scat_sta = delta_deg(scat_lat, scat_lon, sta_lat, sta_lon)

    # TauP rays
    rp1 = get_rp_for_leg(model, phase1, evt_z, d_evt_scat,scat_z)   # event depth
    rp2 = get_rp_for_leg(model, phase2, 0 , d_scat_sta,scat_z)  # scatter depth treated as receiver depth

    if rp1 is None:
        print("didn't work for leg 1",row)
        return None, None

    if rp2 is None:
        print("didn't work for leg 2",row)
        return None, None

    # Convert to XYZ
    # leg1 geometry uses the evt to scat gcp
    rx1, ry1, rz1 = get_xyz_gcp_ray(rp1, evt_lat, evt_lon, scat_lat, scat_lon)

    # leg2 geometry uses the station to scat gcp
    rx2, ry2, rz2 = get_xyz_gcp_ray(rp2, sta_lat, sta_lon, scat_lat, scat_lon)

    return (rx1, ry1, rz1), (rx2, ry2, rz2)

csv_path = "/Users/keyser/Research/sct_wat/scattererwhereartthou/examples/swat_230402_180411.csv"
df = pd.read_csv(csv_path)
single_phase=["P","Ped"]
bounce=["pP","PP"]

# df = df.sample(n=500, random_state=0)
df["n_bounces"] = df["evt_scat_phase"].isin(bounce).astype(int) + df["sta_scat_phase"].isin(bounce).astype(int)

# df= df[df["n_bounces"] == 2]
# sys.exit()
###
evt_lat = float(df["evtlat"].iloc[0]);  evt_lon = float(df["evtlon"].iloc[0]);  evt_z = float(df["evtdepth"].iloc[0])
sta_lat = float(df["stalat"].iloc[0]);  sta_lon = float(df["stalon"].iloc[0]);  sta_z = 0.0

evt_x, evt_y, evt_z3 = ecef_from_latlon_depth(evt_lat, evt_lon, evt_z)
sta_x, sta_y, sta_z3 = ecef_from_latlon_depth(sta_lat, sta_lon, sta_z)

model = TauPyModel(model="iasp91")

delta_deg_val= delta_deg(evt_lat, evt_lon, sta_lat, sta_lon)

# rps = model.get_ray_paths(source_depth_in_km=evt_z, distance_in_degree=delta_sr_deg, phase_list=["P","PP"])
rp = get_rp_for_leg(model, "P", evt_z, delta_deg_val,0)
rp_PP=get_rp_for_leg(model, "PP", evt_z, delta_deg_val,0)

ray_p_P=np.round(rp.ray_param_sec_degree,3)
GCP_time=np.round(rp.time,3)
# sys.exit()
# TauP raypath distances are in radians -> degrees
rx_P,ry_P,rz_P=get_xyz_gcp_ray(rp,evt_lat, evt_lon,sta_lat, sta_lon)
rx_PP,ry_PP,rz_PP=get_xyz_gcp_ray(rp_PP,evt_lat, evt_lon,sta_lat, sta_lon)

# Compute XYZ for scatterers
x, y, z = ecef_from_latlon_depth(df["scatlat"].to_numpy(),
                                 df["scatlon"].to_numpy(),
                                 df["scatdepth"].to_numpy())

# ["scatlat", "scatlon", "scatdepth", "scatdistdeg", "scatbaz", "sta_scat_p", "scat_time",
# 'sta_scat_phase','evt_scat_phase', "evtlat", "evtlon", "evtdepth", "stalat", "stalon",'baz_GCP']

# color by: del baz, absolute P, del time, # bounces ..

raylegs = []  # list of (leg1_xyz, leg2_xyz)

for row in df.itertuples(index=False):
    leg1, leg2 = raypaths_for_row(model, row)
    if leg1 is None:
        # break
        continue
    raylegs.append((leg1, leg2))

print("Computed raypaths for rows:", len(raylegs))

xs1, ys1, zs1 = [], [], []
xs2, ys2, zs2 = [], [], []

for (leg1, leg2) in raylegs:
    rx, ry, rz = leg1
    xs1.extend(rx); ys1.extend(ry); zs1.extend(rz)
    xs1.append(None); ys1.append(None); zs1.append(None)

    rx, ry, rz = leg2
    xs2.extend(rx); ys2.extend(ry); zs2.extend(rz)
    xs2.append(None); ys2.append(None); zs2.append(None)


dbaz = wrap180(df["scatbaz"].to_numpy() - df["baz_GCP"].to_numpy())
pval = df["sta_scat_p"].to_numpy()
del_time=df["scat_time"].to_numpy() - GCP_time
N_bounce=df["n_bounces"].to_numpy()

# Hover text
hover = (
    "scatlat=" + df["scatlat"].round(3).astype(str) +
    ", scatlon=" + df["scatlon"].round(3).astype(str) +
    ", z=" + df["scatdepth"].round(1).astype(str) + " km" +
    "<br>dbaz=" + np.round(dbaz, 2).astype(str) + " deg" +
    "<br>p=" + df["sta_scat_p"].round(3).astype(str) +
    "<br>dt? scat_time=" + df["scat_time"].round(2).astype(str) +
    "<br>phases: " + df["evt_scat_phase"].astype(str) + " / " + df["sta_scat_phase"].astype(str))

# sys.exit()
fig = go.Figure()
fig.add_trace(go.Scatter3d(
    x=rx_P, y=ry_P, z=rz_P,
    mode="lines",
    line=dict(width=4,color='black'),
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
    name="Δ baz",
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
    name="Slow",
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
        color=del_time,
        colorscale="plasma",
        colorbar=dict(title="Δ time (s)", x=-0.12, xanchor="left", len=0.75),
        opacity=0.7,
    ),
    name="Δt",
    text=hover,
    hoverinfo="text",
    visible=True
))

fig.add_trace(go.Scatter3d(
    x=x, y=y, z=z,
    mode="markers",
    marker=dict(
        size=2.5,
        color=N_bounce,
        cmin=0,
        cmax=2,
        colorscale=[
        [0.00, "yellow"], [0.333333, "yellow"],
        [0.333334, "darkorange"], [0.666666, "darkorange"],
        [0.666667, "navy"], [1.00, "navy"],   ],
        colorbar=dict(title="# surface bounces", x=-0.12, xanchor="left", len=0.5,tickmode="array",
        tickvals=[0, 1, 2],
        ticktext=["0", "1", "2"],),
        opacity=0.5,
    ),
    name="#",
    text=hover,
    hoverinfo="text",
    visible=True
))

# for (leg1, leg2) in raylegs:
#     rx, ry, rz = leg1
#     fig.add_trace(go.Scatter3d(x=rx, y=ry, z=rz, mode="lines",line=dict(width=2.5,color='cadetblue'),name="rays",opacity=0.65, showlegend=False,visible=True))
#     rx2, ry2, rz2 = leg2
#     fig.add_trace(go.Scatter3d(x=rx2, y=ry2, z=rz2, mode="lines",line=dict(width=2.5,color='indianred'),name="rays",opacity=0.65, showlegend=False,visible=True))


fig.add_trace(go.Scatter3d(
    x=xs1, y=ys1, z=zs1,
    mode="lines",
    line=dict(width=1.5, color="cadetblue"),
    name="ray:evt2scat",
    opacity=0.35,
    showlegend=False,
    visible=False
))
fig.add_trace(go.Scatter3d(
    x=xs2, y=ys2, z=zs2,
    mode="lines",
    line=dict(width=1.5, color="indianred"),
    name="ray:sta2scat",
    opacity=0.35,
    showlegend=False,
    visible=False
))

fig.update_layout(
    #title=f"Scatterers between sP & PP. P slow={ray_p_P:.2f}",
    # scene controls 3D axes, aspect, camera, etc.
    scene=dict(xaxis_title="X (km)", yaxis_title="Y (km)", zaxis_title="Z (km)",
        aspectmode="data"),
    margin=dict(l=120, r=180, t=60, b=40),
    legend=dict(x=1.02, y=1.0, xanchor="left", yanchor="top"),
    updatemenus=[dict(
        type="dropdown",
        x=1.02, y=0.8,
        buttons=[
            dict(label="Δ baz",
                 method="update",
                 args=[{"visible": [True,True, True, True, True, False, False,False, False,False]},
                       {"title": f"Scatterers between sP & PP. P slow={ray_p_P:.2f}"}]),
            dict(label="Slow",
                 method="update",
                 args=[{"visible": [True,True, True, True, False, True, False,False, False,False]},
                       {"title": f"Scatterers between sP & PP. P slow={ray_p_P:.2f}"}]),
            dict(label="Δt",
                 method="update",
                 args=[{"visible": [True,True, True, True, False, False, True,False, False,False]},
                       {"title": f"Scatterers between sP & PP. P slow={ray_p_P:.2f}"}]),
            dict(label="#",
                 method="update",
                 args=[{"visible": [True,True, True, True, False, False, False,True,False, False]},
                       {"title": f"Scatterers between sP & PP. P slow={ray_p_P:.2f}"}]),
            dict(label="rays",
                 method="update",
                 args=[{"visible": [True,True, True, True, False, False, False,True,True,True]},
                       {"title": f"Scatterers between sP & PP. P slow={ray_p_P:.2f}"}]),
            ],
    )]
)


fig.update_layout(
    title=f"Scatterers between sP & PP. P slow={ray_p_P:.2f}")
fig.show()
