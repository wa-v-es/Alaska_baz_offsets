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

# sys.exit()

model="iasp91"
phase="P"   # reference phase
max_dist_step=1.0 # max separation between path scatterers in degrees, default is 2 deg
bazoffset=0
bazdelta=1
sta_scat_revphase="P,Ped,PP,PS" ###
evt_scat_phase="p,s,P,S,Ped,Sed,pP,sP,pS,sS,PP,SS,SP,PS"

sta_scat_revphase='P,Ped,PP'
evt_scat_phase='p,P,Ped'


with taup.TauPServer(taup_path=taup_path) as taupserver:

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
    traveltimes = time # used for absolute
    print(f"slow: {slow} traveltimes: {traveltimes}")
    ans = swat.find_via_path(slow[0], time[0], bazoffset=baz[0], bazdelta=bazdelta)
    swatList.append(ans)
###
