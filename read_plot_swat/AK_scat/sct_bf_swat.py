#!/usr/bin/env python
# reads scatterers picked by hand in vespagrams and finds scatterer locations
taup_path="~/Research/sct_wat/TauP/build/install/TauP/bin/taup"

import csv
import taup
from scattererwhereartthou import SWAT, mapplot, sliceplot
import sys,re,os
import glob as glob
import numpy as np
import matplotlib.pyplot as plt
##
def plot_3d_locations(points):
    # Plot sctrs in 3D as latitude, longitude and depth.
    points = np.asarray(points)
    lat = points[:, 0]
    lon = points[:, 1]
    depth = points[:, 2]

    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(111, projection="3d")

    sc = ax.scatter(lon,lat,depth,c=depth,cmap="viridis",s=50,edgecolor="k")

    ax.set_xlabel("Longitude (°)")
    ax.set_ylabel("Latitude (°)")
    ax.set_zlabel("Depth (km)")

    ax.invert_zaxis()

    cbar = fig.colorbar(sc, ax=ax, pad=0.1,shrink=0.5,fraction=.15)
    cbar.set_label("Depth (km)")
    plt.tight_layout()
    plt.show()
##
file="/Users/keyser/Research/AK_all_stations/sac_files_.1slow/220914_110406_PA_inc2_r2.5/py_picks/grid_num_109_2022914114_AK_PICKS_amp_f_3.dat"
#C1-'SRC_LAT' C2-'SRC_LON' C3-'SRC_DEP' C4-'REC_LAT' C5-'REC_LON' C6-'DIST' C7-'BAZ' C8-'SCAT_TIME' C9-'SCAT_SLOW' C10-'SCAT_BAZ' C11-'ABS_BAZ' C12-'SNR_BEAM'
with open(file, 'r') as f:
    l = [line.split() for line in f]

evt=(float(l[0][0]),float(l[0][1]))
eventdepth=(float(l[0][2]))
sta=(float(l[0][3]),float(l[0][4]))

time,slow,baz=[],[],[]
for sct in l:
    time.append(float(sct[7]))
    slow.append(float(sct[8]))
    baz.append(float(sct[9]))

print(f'Backazi from text: {float(l[0][6])}')
# sys.exit()

model="iasp91"
phase="P"   # reference phase
max_dist_step=2.0 # max separation between path scatterers in degrees, default is 2 deg
bazoffset=baz[0]
bazdelta=1
sta_scat_revphase="P,Ped,PP,PS" ###
# evt_scat_phase="p,s,P,S,Ped,Sed,pP,sP,pS,sS,PP,SS,SP,PS"

sta_scat_revphase='P,Ped,PP'
evt_scat_phase='p,P,Ped'


with taup.TauPServer(taup_path=taup_path) as taupserver:

    # for d in distazResult.distances:
    #     km = f"Km: {d.km}" if d.km is not None else ""
    #     print(f"{d.disttype.type} from {sta} to {evt}: Dist: {d.deg} Az: {d.az} Baz: {d.baz} {km}")

    swatList = []
    swat = SWAT(taupserver, eventdepth, model=model,
        sta_scat_revphase=sta_scat_revphase,
        evt_scat_phase=evt_scat_phase)
    swat.event(*evt)
    swat.station(*sta)
    swat.dist_step = max_dist_step
    baz_GCP=swat.es_baz

    slow_list=[slow[0]-.25,slow[0],slow[0]+.25]
    time_list=[time[0]-2,time[0],time[0]+2]
    print(f"slow: {slow_list}, traveltimes: {time[0]}, bazOff:{bazoffset}")
    for i,sl in enumerate(slow_list):
        ans = swat.find_via_path(sl, time_list, bazoffset=bazoffset, bazdelta=bazdelta)
        print(f"Length of sct: {len(ans.scatterers)}, for sl:{sl}, time:{time[0]}")
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
        # print(f"Phase: {sct.evt_scat.phase} & {sct.sta_scat_phase}. Lat, Long, depth:{sct.scat.lat:.2f}, {sct.scat.lon:.2f}, {sct.scat.depth:.2f}")
        # print("--------------------------------------------------------------------------------")
        sct_loc.append((sct.scat.lat,sct.scat.lon,sct.scat.depth))
    #
print(f"Length of all sct: {len_all}")

plot_3d_locations(sct_loc)
# print("NEED TO MAKE A FUNCTION TO GET LAT LON DEPTH AND A FUNCTION TO PLOT IT!!!")
# for sct in swatList[0].scatterers:
#     bazdiff=sct.scat_baz-baz_GCP
#     print(f'baz sct {sct.scat_baz} - ori {baz_GCP}: {bazdiff}')
# ###
