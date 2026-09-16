#!/usr/bin/env python
# reads scatterers picked by hand (json file) in vespagrams and finds scatterer locations using SWAT.

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
from scipy.spatial import Delaunay

##
def create_panda(swatList,ans):
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

def get_hull_volume(points):
    # Plot sctrs in 3D as x y z.
    points = np.asarray(points)
    lat = points[:, 0]
    lon = points[:, 1]
    depth = points[:, 2]
    # converting to earth centered cartesian..
    R = 6371.0  # km
    r = R - depth

    lat_rad = np.radians(lat)
    lon_rad = np.radians(lon)

    X = r * np.cos(lat_rad) * np.cos(lon_rad)
    Y = r * np.cos(lat_rad) * np.sin(lon_rad)
    Z = r * np.sin(lat_rad)

    xyz = np.column_stack((X, Y, Z))
    hull = ConvexHull(xyz)
    print(f"Convex hull volume: {hull.volume:.2f} km³ / {hull.volume/(111.32**3):.2f} degree³")

    return xyz,hull

def plot_hull_3d(xyz,hull,figname=None):
    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(111, projection="3d")

    for simplex in hull.simplices:
        simplex = np.append(simplex, simplex[0])

        ax.plot(xyz[simplex, 0],xyz[simplex, 1],xyz[simplex, 2],"k-",linewidth=0.8,alpha=0.25)
    sc = ax.scatter(xyz[:, 0], xyz[:, 1], xyz[:, 2],c=depth,cmap="viridis",s=50,edgecolor="k")

    ax.set_xlabel("X (km)")
    ax.set_ylabel("Y (km)")
    ax.set_zlabel("Z (km)")

    ax.invert_zaxis()
    ax.view_init(elev=-20, azim=-35,roll=10)
    cbar = fig.colorbar(sc, ax=ax, pad=0.1,shrink=0.5,fraction=.15)
    cbar.set_label("Depth (km)")
    ax.set_title(f"Volume: {hull.volume:.2f} km³ / {hull.volume/(111.32**3):.2f} deg³. # sct:{len(points)}")
    plt.tight_layout()
    if figname:
        plt.savefig(figname,dpi=300,bbox_inches='tight', pad_inches=0.1)
    plt.show()

def swat_sct_volume(taupserver,scatterers,i,scat,figname=None):
    """
    read one sct at a time. using swat, finds potential sctrs.
    finds delta time/baz/slow from picked vals.
    for the scatteres, finds a convex volume.
    the volume caluclate is done in cartesian coordinates.
    around lat0 and lon0 (mean lat long of all scat locations.)
    """
    model="iasp91"
    phase="P"   # reference phase
    max_dist_step=2.0 # max separation between path scatterers in degrees, default is 2 deg
    min_dist_step=0.05
    evt=(scatterers['SRC_LAT'] ,scatterers['SRC_LON'])
    eventdepth=(scatterers['SRC_DEP'])
    sta=(scatterers['REC_LAT'] ,scatterers['REC_LON'])
    ######
    slow_sct=scat['SCAT_slow_max']
    time_sct=round(scat['SCAT_time_max'],2)
    bazoffset=scat['SCAT_baz_max']
    sc_time_delta=round(max(scat['SCAT_sl_time_5_delta'],scat['SCAT_bz_time_5_delta'],3),2)
    sc_slow_delta=round(max(scat['SCAT_slow_5_delta'], .1),2)
    sc_baz_delta= max(scat['SCAT_baz_5_delta'], 1)
    #
    print("Scatterer props..")
    print(f"Delta time/slow/baz used: {sc_time_delta}sec, {sc_slow_delta}sec/deg, {sc_baz_delta}deg")
    ###
    bazdelta=sc_baz_delta/2

    sta_scat_revphase="P,Ped,PP,PS" ###
    # evt_scat_phase="p,s,P,S,Ped,Sed,pP,sP,pS,sS,PP,SS,SP,PS"

    sta_scat_revphase='P,Ped,PP'
    evt_scat_phase='p,P,Ped'

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
    sct_loc=[]
    for SctDist in swatList:
        len_all+=len(SctDist.scatterers)
        for sct in SctDist.scatterers:
            # print(f"")
            # print(f"slow:{sct.sta_scat_rayparam}, total_time:{sct.scat.time+sct.evt_scat.time:.2f}, baz: {sct.scat_baz-baz_GCP:.2f}")
            # print(f"Phase: {sct.evt_scat.phase} & {sct.sta_scat_phase}. Lat, Long, depth:{sct.scat.lat:.4f}, {sct.scat.lon:.4f}, {sct.scat.depth:.4f}")
            sct_loc.append((sct.scat.lat,sct.scat.lon,sct.scat.depth))
        #
    return swatList,sct_loc

def fibonacci_sphere(number_points):
    #https://stackoverflow.com/questions/9600801/evenly-distributing-n-points-on-a-sphere
    phi = np.pi * (np.sqrt(5.) - 1.)

    i = np.arange(number_points)
    y = 1 - 2 * i / (number_points - 1)
    radius = np.sqrt(1 - y**2)
    theta = phi * i

    x = np.cos(theta) * radius
    z = np.sin(theta) * radius

    return np.column_stack((x, y, z))

def find_weights_Scatterer(hull_convex,fib_grid):
    # weights, inside = hull_to_bin_weights(hull_convex,lat0,lon0)

    weights = np.zeros(len(fib_grid))

    inside = np.all(fib_grid @ hull_convex.equations[:, :-1].T+ hull_convex.equations[:, -1] <= 1e-8,axis=1)

    n_inside=inside.sum()
    print(f"Number of cells inside hull: {n_inside}")

    if n_inside>0:
        weights[inside] += 1/n_inside
    return weights
    # counts, edges = np.histogramdd(samples,bins=[lat_edges, lon_edges, dep_edges])

#
def create_Fib_grid(delta_deg=1,depth_delta=100):
    """
    for a delta_deg2 area, creates a fibonacci_sphere for each depth (depth_delta).
    returns fib_grid.
    the number of points at each depth slice change such that the area is conserved.
    """
    R = 6371.0
    #area per point at surface
    area_point = (np.radians(delta_deg) * R)**2

    radii = np.arange(2900, 6370, depth_delta)
    n_points = np.round(4 * np.pi * radii**2 / area_point).astype(int)

    fib_grid = []
    for r, n in zip(radii, n_points):
        fib = fibonacci_sphere(n)
        fib_grid.append(fib * r)

    fib_grid = np.vstack(fib_grid)

    return fib_grid

def loadScatterers(sct_json):
    with open(sct_json, "r") as file:
        scatterers = json.load(file)
        return scatterers

def justOne():
    taup_path="~/Research/sct_wat/TauP/build/install/TauP/bin/taup"
    sct_json="/Users/keyser/Research/AK_all_stations/sac_files_.1slow/220914_110406_PA_inc2_r2.5/py_picks/grid_num_109_2022914114_PICKS.json"
    ### bin edges..

    scatterers = loadScatterers(sct_json)
    i=2
    scat = scatterers['sct'][i]
    with taup.TauPServer(taup_path=taup_path) as taupserver:
        figname='220914_109_{}_P_min_.05.png'.format(i)
        hull_convex,swatList=swat_sct_volume(taupserver,scatterers,i, scat)
        fib_grid=create_Fib_grid(delta_deg=1,depth_delta=100)
        # weights_all = np.zeros(len(fib_grid))
        weights=find_weights_Scatterer(hull_convex,fib_grid)

    return weights

# def main():
taup_path="~/Research/sct_wat/TauP/build/install/TauP/bin/taup"
sct_json="/Users/keyser/Research/AK_all_stations/sac_files_.1slow/220914_110406_PA_inc2_r2.5/py_picks/grid_num_109_2022914114_PICKS.json"
### bin edges..

scatterers = loadScatterers(sct_json)
with taup.TauPServer(taup_path=taup_path) as taupserver:
    for i, scat in enumerate(scatterers['sct']):
        figname='220914_109_{}_P_min_.05.png'.format(i)
        swatList,sct_loc=swat_sct_volume(taupserver,scatterers,i, scat)
        # df= create_panda(swatList,ans)
        # read_swat_plotly(taupserver,csv_path=None,data_swat=df,plotrays=True)
        xyz,hull_convex=get_hull_volume(sct_loc)
        # plot_hull_3d(xyz,hull,figname=None)
        fib_grid=create_Fib_grid(delta_deg=.5,depth_delta=100)
        weights=find_weights_Scatterer(hull_convex,fib_grid)

        print(f'Scatterer#{i+1} done')
        print("------------------\n")
        # break

# hull_convex,lat0,lon0,swatList=swat_sct_volume(taupserver,scatterers,i, scat)
# if __name__ == '__main__':
#     main()
