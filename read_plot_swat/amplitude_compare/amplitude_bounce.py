#compares amplitude of phases wrt direct P"
import csv
import taup
from scattererwhereartthou import SWAT, mapplot, sliceplot
import sys,re,os
import glob as glob
import numpy as np
# from obspy.taup import TauPyModel
import requests
from collections import defaultdict
import math
import matplotlib.pyplot as plt
import sys
import pandas as pd
import plotly.express as px
import sys
###

def parse_file(path):
    """
    wrote when I copied the result ffrom taup_web :]
    """
    grouped = defaultdict()
    with open(path, 'r') as f:
        lines = f.readlines()
        for raw in lines[3:]:
            line = raw.strip()
            parts = line.split()
            phase = parts[2]
            amp_str = parts[-2]
            amp = float(amp_str)
            # keep the largest absolute amplitude for each phase within the same record
            prev = grouped.get(phase)
            if prev is None or abs(amp) > abs(prev):
                grouped[phase] = amp
    return grouped
###
def get_dict_amps(TimeResult):
    """
    readfs taup_py output and gets phase names and amplitudes.
    """
    grouped = defaultdict()
    for a in TimeResult.arrivals:
        phase=a.phase
        amp=float(a.amp.factorpsv)
        prev = grouped.get(phase)
        if prev is None or abs(amp) > abs(prev):
            if prev != None: #debug step
                pass
                # print(f"prev:{prev}; phase:{phase}; abs amp{abs(amp)}")
            grouped[phase] = amp

    return grouped

def get_phase_ratio(taupserver,model,evt,eventdepth,sta,phases,s,d,r,):
    params = taup.TimeQuery()
    params.model(model)
    params.event(*evt)
    params.sourcedepth(eventdepth)
    params.station(*sta)
    params.phase('P')
    params.amp(True)
    params.mw(7)
    TimeResult = params.calc(taupserver)
    # get amplitude of P for explosion
    grouped = get_dict_amps(TimeResult)
    amp_ref_exp_abs = abs(grouped['P'])

    # now get amplitudfe of all phases for sdr vals
    params.phase(phases)
    params.strikediprake([s,d,r])
    # params.az(20)
    TimeResult = params.calc(taupserver)
    grouped = get_dict_amps(TimeResult)
    amp_ref_abs = abs(grouped['P'])

    phase_ratios = defaultdict(list)
    phase_ratios_exp = defaultdict(list)

    for phase, amp in grouped.items():

        ratio = np.round(abs(amp) / amp_ref_abs,3)
        ratio_exp = np.round(abs(amp) / amp_ref_exp_abs,3)

        phase_ratios[phase].append(ratio)
        phase_ratios_exp[phase].append(ratio_exp)


    return grouped,phase_ratios,phase_ratios_exp,amp_ref_exp_abs

####
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
#rake −180° and 180°
plt.rcParams.update({'font.size': 15})

plt.figure(figsize=(12, 6))
all_phaseR = defaultdict(list)

##
with taup.TauPServer(taup_path=taup_path) as taupserver:

    for s in np.arange(0, 181, 45):
        for d in np.arange(0, 91, 45):
            for r in np.arange(-180, 181, 45):
                amps,phase_ratios,phase_ratios_exp,amp_P_exp_abs= get_phase_ratio(taupserver,model,evt,eventdepth,sta,phases,s,d,r)
                #### change here for Explosion P vs normal P.
                for ph, vals in phase_ratios_exp.items():
                    all_phaseR[ph].extend(vals)
                for i, p in enumerate(phases):
                    y = phase_ratios_exp[p]
                    x = [i + 1] * len(y)
                    color = 'darkseagreen' if ('S' in p or 's' in p) else 'palevioletred'
                    plt.scatter(x, y, marker='X', alpha=0.2,s=55, color=color)

# keep uniques vals in
for ph in all_phaseR:
    all_phaseR[ph] = list(set(all_phaseR[ph]))

min_phaseR = {ph: (float(np.min(v)) if v else None) for ph, v in all_phaseR.items()}
max_phaseR = {ph: (float(np.max(v)) if v else None) for ph, v in all_phaseR.items()}
mean_phaseR = {ph: (float(np.round(np.mean(v),3)) if v else None) for ph, v in all_phaseR.items()}
median_phaseR = {ph: (float(np.round(np.median(v),3)) if v else None) for ph, v in all_phaseR.items()}

print('mean Phase Ratios:\n', mean_phaseR)
print('median Phase Ratios:\n', median_phaseR)


plt.xticks(range(1, len(phases) + 1), phases)
# plt.yscale('log')
plt.ylabel("Amplitude P/ phase")
plt.xlabel("phase")
# plt.title(f"Phase amp P exp")
plt.title(f"Phase amp ratios wrt P exp; Δ20")

plt.grid(which='both', linestyle='--', linewidth=0.5, alpha=0.6)
plt.tight_layout()
# plt.savefig("230402_180411_amps_expl.png", dpi=500, bbox_inches='tight', pad_inches=0.1)
plt.savefig("230402_180411_amps_wrt_Pexp_20.png", dpi=500, bbox_inches='tight', pad_inches=0.1)
# plt.show()
plt.close()
sys.exit()

#########
## bit for getting actaul phases from a swat output
csv_path='../230402_180411_S_10baz.csv'
# csv_path='../230402_180411_reso_dup.csv'


df = pd.read_csv(csv_path)
bounce=["pP","PP",'sP','SP','SS','sS']
s_phases=['sP','SP','SS','sS','s','S','Sed']
print('Len of read csv:',len(df),'\n')

# df = df.sample(n=500, random_state=0)
df["n_bounces"] = df["evt_scat_phase"].isin(bounce).astype(int) + df["sta_scat_phase"].isin(bounce).astype(int)
df["n_Sphases"] = df["evt_scat_phase"].isin(s_phases).astype(int) + df["sta_scat_phase"].isin(s_phases).astype(int)
##
# dups_all = df[df.duplicated(keep=False)]

df = df.drop_duplicates().reset_index(drop=True)
print('Len of unique scats:',len(df),'\n')
##
#### count pairs
pair_counts = (df[["evt_scat_phase", "sta_scat_phase"]]
      .astype("string")
      .apply(lambda c: c.str.strip())
      .groupby(["evt_scat_phase", "sta_scat_phase"])
      .size()
      .sort_values(ascending=False))

# print(pair_counts)
pairs_df = pair_counts.reset_index(name="count")
print(pairs_df)

# matrix of all combinations
mat = pairs_df.pivot(index="evt_scat_phase", columns="sta_scat_phase", values="count").fillna(0)
row_order = mat.sum(axis=1).sort_values(ascending=False).index
col_order = mat.sum(axis=0).sort_values(ascending=False).index
mat2 = mat.loc[row_order, col_order]
mat_int = mat2.astype(int)

### for p, Ped, s, sed
# map aliases for these
phase_alias = {"p": "P", "Ped": "P",
            "s": "S", "Sed": "S"}

def alias(ph: str) -> str:
    return phase_alias.get(ph, ph)

# vector for row/col phase amps
row_w = np.array([float(median_phaseR.get(alias(ph), 1.0)) for ph in mat2.index], dtype=float)      # (nrow,)
col_w = np.array([float(median_phaseR.get(alias(ph), 1.0)) for ph in mat2.columns], dtype=float)   # (ncol,)

mat_color = row_w[:, None] * col_w[None, :]

color_by_amp=True

# sys.exit()
if color_by_amp:
    fig = px.imshow(
        mat_color, #mat2 #np.log10(mat + 1)
        x=mat2.columns,
        y=mat2.index,
        labels=dict(x="sta_scat_phase", y="evt_scat_phase", color="Median Amp*"),
        aspect="auto",
        color_continuous_scale="PuBu",)#"PuBu"
else:
    fig = px.imshow(
        mat2, #mat2 #np.log10(mat + 1)
        x=mat2.columns,
        y=mat2.index,
        labels=dict(x="sta_scat_phase", y="evt_scat_phase", color="#"),
        aspect="auto",
        color_continuous_scale="GnBu",)
fig.update_traces(text=mat_int.values, texttemplate="%{text:d}")

fig.update_layout(
    xaxis_title_font=dict(size=18),
    yaxis_title_font=dict(size=18),
    xaxis_tickfont=dict(size=14),
    yaxis_tickfont=dict(size=14),)
fig.update_layout(width=400, height=800)
fig.update_layout(
    title="Median Amp*amp wrt P expl",
    xaxis_side="bottom",)
# fig.write_image("freq_phases", format="png")
fig.show()
