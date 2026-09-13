#!/usr/bin/env python
# reads scatterers picked by hand (json file) in vespagrams and finds scatterer locations using SWAT.
taup_path="~/Research/sct_wat/TauP/build/install/TauP/bin/taup"

import csv
import taup
from scattererwhereartthou import SWAT, mapplot, sliceplot
import sys,re,os
import glob as glob
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
import pandas as pd
sys.path.append("../")
from swat_out_plot import read_swat_plotly
import json

##
def create_panda(swatList):
    rows = []
    for SctDist in swatList:
        for sc in SctDist.scatterers:
            rows.append({
                "scatlat": sc.scat.lat,
                "scatlon": sc.scat.lon,
                "scatdepth": sc.scat.depth,
                "scatdistdeg": sc.scat.distdeg,
                "scatbaz": sc.scat_baz,
                "sta_scat_p": sc.sta_scat_rayparam,
                "scat_time": sc.scat.time + sc.evt_scat.time,
                "sta_scat_phase": sc.sta_scat_phase,
                "evt_scat_phase": sc.evt_scat.phase,
                "evtlat": ans.evtlat,
                "evtlon": ans.evtlon,
                "evtdepth": ans.evtdepth,
                "stalat": ans.stalat,
                "stalon": ans.stalon,
                "baz_GCP": ans.esbaz,
                "del_baz":sc.scat_baz-ans.esbaz
            })

    return pd.DataFrame(rows)

def plot_3d_locations(points,figname=None):
    # Plot sctrs in 3D as latitude, longitude and depth.
    points = np.asarray(points)
    lat = points[:, 0]
    lon = points[:, 1]
    depth = points[:, 2]
    ### lat lon to xyx..
    lat0 = np.mean(lat)
    lon0 = np.mean(lon)

    R = 6371.0  # Earth radius in km

    x = np.radians(lon - lon0) * R * np.cos(np.radians(lat0))
    y = np.radians(lat - lat0) * R
    z = depth

    xyz = np.column_stack((x, y, z))
    hull = ConvexHull(xyz)

    print(f"Convex hull volume: {hull.volume:.2f} km³ / {hull.volume/(111.32**3):.2f} degree³")

    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(111, projection="3d")

    ###
    for simplex in hull.simplices:
        simplex = np.append(simplex, simplex[0])

        ax.plot(x[simplex],y[simplex],z[simplex],"k-",linewidth=0.8,alpha=0.25)
    sc = ax.scatter(x,y,z,c=depth,cmap="viridis",s=50,edgecolor="k")

    ax.set_xlabel("Distance E-W (km)")
    ax.set_ylabel("Distance N-S (km)")
    ax.set_zlabel("Depth (km)")

    ax.invert_zaxis()
    ax.view_init(elev=-20, azim=-35,roll=10)
    cbar = fig.colorbar(sc, ax=ax, pad=0.1,shrink=0.5,fraction=.15)
    cbar.set_label("Depth (km)")
    ax.set_title(f"Volume: {hull.volume:.2f} km³ / {hull.volume/(111.32**3):.2f} deg³. # sct:{len(points)}")
    plt.tight_layout()
    if figname:
        plt.savefig(figname,dpi=300,bbox_inches='tight', pad_inches=0.1)
    plt.show()

    return hull
##

sct_json="/Users/keyser/Research/AK_all_stations/sac_files_.1slow/220914_110406_PA_inc2_r2.5/py_picks/grid_num_109_2022914114_PICKS.json"

with open(sct_json, "r") as file:
    scatterers = json.load(file)

# sys.exit()

with taup.TauPServer(taup_path=taup_path) as taupserver:

    model="iasp91"
    phase="P"   # reference phase
    max_dist_step=2.0 # max separation between path scatterers in degrees, default is 2 deg
    min_dist_step=0.05
    # ith=2
    evt=(scatterers['SRC_LAT'] ,scatterers['SRC_LON'])
    eventdepth=(scatterers['SRC_DEP'])
    sta=(scatterers['REC_LAT'] ,scatterers['REC_LON'])
    for i, scat in enumerate(scatterers['sct']):

        slow_sct=scat['SCAT_slow_max']
        time_sct=round(scat['SCAT_time_max'],2)
        bazoffset=scat['SCAT_baz_max']
        sc_time_delta=round(max(scat['SCAT_sl_time_5_delta'],scat['SCAT_bz_time_5_delta'],3),2)
        sc_slow_delta=round(max(scat['SCAT_slow_5_delta'], .1),2)
        sc_baz_delta= max(scat['SCAT_baz_5_delta'], 1)
        print(f"Delta time/slow/baz used: {sc_time_delta}sec, {sc_slow_delta}sec/deg, {sc_baz_delta}deg")
        ###
        bazdelta=sc_baz_delta/2

        sta_scat_revphase="P,Ped,PP,PS" ###
        # evt_scat_phase="p,s,P,S,Ped,Sed,pP,sP,pS,sS,PP,SS,SP,PS"

        sta_scat_revphase='P,Ped,PP'
        evt_scat_phase='p,P,Ped'

        # for d in distazResult.distances:
        #     km = f"Km: {d.km}" if d.km is not None else ""
        #     print(f"{d.disttype.type} from {sta} to {evt}: Dist: {d.deg} Az: {d.az} Baz: {d.baz} {km}")

        swatList = []
        swat = SWAT(taupserver, eventdepth, model=model,
            sta_scat_revphase=sta_scat_revphase,
            evt_scat_phase=evt_scat_phase)
        swat.event(*evt)
        swat.station(*sta)
        swat.max_dist_step = max_dist_step
        swat.min_dist_step = min_dist_step

        baz_GCP=swat.es_baz
        # ans = swat.find_via_path(5.25, 949.65, bazoffset=7.5, bazdelta=.2)

        slow_list=[slow_sct-sc_slow_delta/2,slow_sct,slow_sct+sc_slow_delta/2]
        time_list=[time_sct-sc_time_delta/2,time_sct,time_sct+sc_time_delta/2]
        print(f"slow: {slow_list}, traveltimes: {time_list}, bazOff:{bazoffset}, bazdelta:{bazdelta}")
        # for i,sl in enumerate(slow_list):
        ans = swat.find_via_path(slow_list, time_list, bazoffset=bazoffset, bazdelta=bazdelta)
        print(f"Length of sct: {len(ans.scatterers)}")#", for sl:{sl}, time:{time}")
        swatList.append(ans)

        # print(f"bazoff:{swatList[0].bazoffset}, bazdel:{swatList[0].bazdelta}, esbaz:{swatList[0].esbaz}")
        len_all=0
        print(f"\n ....Output.... \n")
        sct_loc=[]
        for SctDist in swatList:
            len_all+=len(SctDist.scatterers)
            for sct in SctDist.scatterers:
                # print(f"")
                # print(f"slow:{sct.sta_scat_rayparam}, total_time:{sct.scat.time+sct.evt_scat.time:.2f}, baz: {sct.scat_baz-baz_GCP:.2f}")
                # print(f"Phase: {sct.evt_scat.phase} & {sct.sta_scat_phase}. Lat, Long, depth:{sct.scat.lat:.4f}, {sct.scat.lon:.4f}, {sct.scat.depth:.4f}")
                # print("--------------------------------------------------------------------------------")
                sct_loc.append((sct.scat.lat,sct.scat.lon,sct.scat.depth))
            #
        print(f"Length of all sct: {len_all}")
        df= create_panda(swatList)
        # read_swat_plotly(taupserver,csv_path=None,data_swat=df,plotrays=True)

        hull_convex=plot_3d_locations(sct_loc,'220914_109_{}_P_min_.05.png'.format(i))
        # break

# print(f"Volume of potential sct: {hull_convex.volume/(111.32**3):.2f} degree³")
