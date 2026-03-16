#compares amplitude of phases wrt direct P"
import csv
import taup
from scattererwhereartthou import SWAT, mapplot, sliceplot
import sys,re,os
import glob as glob
import numpy as np
# from obspy.taup import TauPyModel
import taup
import requests
from collections import defaultdict
import math
import matplotlib.pyplot as plt
###
INPUT_FILE = "amps_230402_180411.txt"   # change to your filename

def parse_file(path):
    """
    Parse file, returning a dict:
      grouped[(col0, col1)][phase] = max_abs_amplitude_for_that_phase
    """
    grouped = defaultdict(lambda: dict())

    with open(path, 'r') as f:
        for raw in f:
            line = raw.strip()
            parts = line.split()

            phase = parts[2]
            amp_str = parts[-2]
            amp = float(amp_str)


            # grouping key for a "record" (change if you want a different grouping)
            record_key = (parts[0], parts[1])

            # keep the largest absolute amplitude for each phase within the same record
            prev = grouped[record_key].get(phase)
            if prev is None or abs(amp) > abs(prev):
                grouped[record_key][phase] = amp

    return grouped

taup_path="~/Research/sct_wat/TauP-3.2.0-SNAPSHOT6/bin/taup"

sta=(64.67, -155.88)
evt=(-4.33,143.16 )   # eq lat, lon -4.33 143.16 70.00
eventdepth=70  # eq depth
model="iasp91"

phases=['PP','sP','pP','PPP','SP']
phases=['PP']

with taup.TauPServer( taup_path=taup_path) as taupserver:
    params = taup.TimeQuery()
    params.model(model)
    params.event(*evt)
    params.sourcedepth(eventdepth)
    params.station(*sta)
    params.phase('P')
    timeResult = params.calc(taupserver)
##
   82.02    70.0   P        731.39     5.216     22.54    15.79      0.0    82.02   = P       4.2e-05  0.0e+00
   82.02    70.0   pP       749.39     5.253    157.29    15.90      0.0    82.02   = pP     -3.5e-05  0.0e+00
   82.02    70.0   sP       757.33     5.245    167.68    15.88      0.0    82.02   = sP      5.3e-05  0.0e+00
   82.02    70.0   PP       919.67     8.214     37.14    25.37      0.0    82.02   = PP     -2.1e-05  0.0e+00
   82.02    70.0   PPP     1030.12     8.975     41.28    27.92      0.0    82.02   = PPP     8.9e-06  0.0e+00
   82.02    70.0   PPP     1038.35     9.851     46.39    30.92      0.0    82.02   = PPP     5.2e-06  0.0e+00
   82.02    70.0   PPP     1039.38     9.671     45.30    30.29      0.0    82.02   = PPP     3.4e-06  0.0e+00
   82.02    70.0   PPP     1129.59    13.545     84.59    44.95      0.0    82.02   = PPP     8.4e-07  0.0e+00
   82.02    70.0   PPP     1131.99    13.407     80.20    44.37      0.0    82.02   = PPP     9.0e-07  0.0e+00
   82.02    70.0   SP      1386.64    12.309     30.04    39.94      0.0    82.02   = SP      4.0e-05  0.0e+00
   82.02    70.0   SP      1394.96    13.545     33.43    44.95      0.0    82.02   = SP      1.4e-05  0.0e+00
   82.02    70.0   SP      1395.81    13.402     33.03    44.35      0.0    82.02   = SP      1.7e-05  0.0e+00
   82.02    70.0   PSP     1421.77    12.540     67.17    40.85      0.0    82.02   = PSP     1.6e-07  0.0e+00
   82.02    70.0   PSP     1422.33    12.950     72.14    42.49      0.0    82.02   = PSP     2.8e-05  0.0e+00
   82.02    70.0   PPP     2366.78     4.596     19.74    13.87      0.0   277.98   = PPP     2.7e-06  0.0e+00
    # for phase in phases:
    #     params.phase(phase)
    #     timeResult = params.calc(taupserver)
