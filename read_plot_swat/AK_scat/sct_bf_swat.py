#!/usr/bin/env python
# reads scatterers picked by hand in vespagrams and finds scatterer locations
taup_path="~/Research/sct_wat/TauP/build/install/TauP/bin/taup"

import csv
import taup
from scattererwhereartthou import SWAT, mapplot, sliceplot
import sys,re,os
import glob as glob
import numpy as np
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

    params = taup.DistazQuery()
    params.geodist(["spherical"])
    # params.geodist(["spherical", "geocentric", "geodetic"])
    params.event(*evt)
    params.station(*sta)
    distazResult = params.calc(taupserver)
    baz_GCP=distazResult.distances[0].baz

    for d in distazResult.distances:
        km = f"Km: {d.km}" if d.km is not None else ""
        print(f"{d.disttype.type} from {sta} to {evt}: Dist: {d.deg} Az: {d.az} Baz: {d.baz} {km}")
    swatList = []
    swat = SWAT(taupserver, eventdepth, model=model,
        sta_scat_revphase=sta_scat_revphase,
        evt_scat_phase=evt_scat_phase)
    swat.event(*evt)
    swat.station(*sta)
    swat.dist_step = max_dist_step

    # for a in timeResult.arrivals:
    #     print(f"Arrival: {a}")
    #     # traveltimes = [a.time+delay for delay in delaytimes] # used when using delay..
    # traveltimes = time # used for absolute
    print(f"slow: {slow[0]} traveltimes: {time[0]}")
    # for i,sl in enumerate(slow):
    ans = swat.find_via_path(slow[0], time[0], bazoffset=bazoffset, bazdelta=bazdelta)
    # print(f"Length of sct: {len(ans.scatterers)}")

    swatList.append(ans)
###
print(f"Length of sct: {len(swatList[0].scatterers)}")
print(f"bazoff:{swatList[0].bazoffset}, bazdel:{swatList[0].bazdelta}, esbaz:{swatList[0].esbaz}")
for sct in swatList[0].scatterers:
    bazdiff=sct.scat_baz-baz_GCP
    print(f'baz sct {sct.scat_baz} - ori {baz_GCP}: {bazdiff}')
###
