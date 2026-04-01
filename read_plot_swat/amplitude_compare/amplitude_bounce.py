#compares amplitude of phases wrt direct P"
import csv
import taup
from scattererwhereartthou import SWAT, mapplot, sliceplot
import sys,re,os
import glob as glob
import numpy as np
# from obspy.taup import TauPyModel
# import taup
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

def get_dict_amps(TimeResult):
    """
    readfs taup_py output and gets pahse names and amplitudes.
    """
    grouped = defaultdict()
    for a in TimeResult.arrivals:
        phase=a.phase
        amp=float(a.amp.factorpsv)
        grouped[phase] = amp

    return grouped


# amp_file = "amps_230402_180411.txt"
# taup_path="~/Research/sct_wat/TauP-3.2.0-SNAPSHOT6/bin/taup"
taup_path="~/Research/sct_wat/TauP/build/install/TauP/bin/taup"

sta=(64.67, -155.88)
evt=(-4.33,143.16 )   # eq lat, lon -4.33 143.16 70.00
eventdepth=70  # eq depth
model="ak135fcont"

phases=['S','P','PP','SS','SP','PS','sP','pP','sS','pS','PPP']
# phases=['PP']
SDR=[45, 90, 0]
plt.rcParams.update({'font.size': 15})

plt.figure(figsize=(12, 6))

with taup.TauPServer(taup_path=taup_path) as taupserver:
    for s in [0,45,90]:
        for d in [0,45,90]:
            for r in [0,45,90]:

                params = taup.TimeQuery()
                params.model(model)
                params.event(*evt)
                params.sourcedepth(eventdepth)
                params.station(*sta)
                params.phase(phases)
                params.amp(True)
                params.mw(7)
                params.strikediprake([s,d,r])
                # params.az(20)
                TimeResult = params.calc(taupserver)

                # print()
                grouped = get_dict_amps(TimeResult)

                # sys.exit()
                phase_ratios = defaultdict(list)
                amp_ref_abs = abs(grouped['P'])

                for phase, amp in grouped.items():

                    ratio = np.round(abs(amp) / amp_ref_abs,3)
                    phase_ratios[phase].append(ratio)

                # phases = dict(
                #     sorted(phase_ratios.items(), key=lambda x: x[1][0], reverse=True))

                # print(phase_ratios)
                for i, p in enumerate(phases):
                    y = phase_ratios[p]
                    x = [i + 1] * len(y)  # align points vertically per phase

                    color = 'darkseagreen' if ('S' in p or 's' in p) else 'palevioletred'

                    plt.scatter(x, y, marker='D', alpha=0.7,s=55, color=color)

plt.xticks(range(1, len(phases) + 1), phases)
# plt.yscale('log')
plt.ylabel("Amplitude P/ phase")
plt.xlabel("phase")
plt.title(f"Phase amp ratios for SDR in 0, 45, 90")

plt.grid(which='both', linestyle='--', linewidth=0.5, alpha=0.6)
plt.tight_layout()
# plt.show()
plt.savefig("230402_180411_amps_0_45_90.png", dpi=500, bbox_inches='tight', pad_inches=0.1)
