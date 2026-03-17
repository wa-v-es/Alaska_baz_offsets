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
import sys
###

def parse_file(path):
    """
    Parse file, returning a dict:
      grouped[(col0, col1)][phase] = max_abs_amplitude_for_that_phase
    """
    # grouped = defaultdict(lambda: dict())
    grouped = defaultdict()

    with open(path, 'r') as f:
        lines = f.readlines()
        for raw in lines[3:]:
            line = raw.strip()
            parts = line.split()

            phase = parts[2]
            amp_str = parts[-2]
            amp = float(amp_str)

            # grouping key
            # record_key = (parts[0], parts[1])

            # keep the largest absolute amplitude for each phase within the same record
            prev = grouped.get(phase)
            if prev is None or abs(amp) > abs(prev):
                grouped[phase] = amp

    return grouped

amp_file = "amps_230402_180411.txt"
grouped = parse_file(amp_file)
phase_ratios = defaultdict(list)
amp_ref_abs = abs(grouped['P'])
for phase, amp in grouped.items():

    ratio = np.round(abs(amp) / amp_ref_abs,3)
    phase_ratios[phase].append(ratio)

phases = dict(
    sorted(phase_ratios.items(), key=lambda x: x[1][0], reverse=True))

print(phase_ratios)

plt.figure(figsize=(10, 6))

for i, p in enumerate(phases):
    y = phase_ratios[p]
    x = [i + 1] * len(y)  # align points vertically per phase

    color = 'rosybrown' if 'S' in p else 'olivedrab'

    plt.scatter(x, y, marker='D', alpha=0.8,s=55, color=color)

plt.xticks(range(1, len(phases) + 1), phases)
# plt.yscale('log')
plt.ylabel("Amplitude P/ phase")
plt.xlabel("phase")
plt.title("Phase amplitude ratios")

plt.grid(which='both', linestyle='--', linewidth=0.5, alpha=0.6)
plt.tight_layout()
plt.show()

sys.exit()
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
    # for phase in phases:
    #     params.phase(phase)
    #     timeResult = params.calc(taupserver)
