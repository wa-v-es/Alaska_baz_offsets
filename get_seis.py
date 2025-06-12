import matplotlib.pyplot as plt
import obspy
from obspy import read, Stream, UTCDateTime,read_events
from obspy.core.event import Origin, Catalog
from obspy.taup import TauPyModel
from obspy.taup.taup_geo import calc_dist,calc_dist_azi
from obspy.geodetics import gps2dist_azimuth
from obspy.clients.fdsn import Client
import numpy as np
from obspy.clients.fdsn import Client as FDSN_Client
from obspy.clients.earthworm import Client as EW_Client
from obspy.clients.fdsn.header import FDSNNoDataException
import datetime
import os
import sys
from rf import read_rf, rfstats,RFStream
from rf import iter_event_data, IterMultipleComponents
from tqdm import tqdm
from obspy.io.sac.util import obspy_to_sac_header
import warnings


###
# get taup time using
#bin/taup find --amp -h 137 --deg 89 --max 2 --time 780 1050 --exclude 20 210 35 | awk 'NR > 3 {if ($3 ~ /[pP]$/) print $0}'
st=read('220914_110406_AK.SCM*')
st.filter('bandpass',freqmin=0.05,freqmax=.5)
end=st[0].stats.endtime
startime=end-450
endtime=end-50
# sys.exit()

sttt=st.trim(startime,endtime)
# sttt.plot()
end_v=end-200
P=759.21
PP=970.52
pP=794.13
sP=808.82
S660P=809.69
pPP=1001.38
S660PP=1027
fig = plt.figure(figsize=(14,2))
plt.tick_params(axis='both', left=False, top=False, right=False, bottom=False, labelleft=False, labeltop=False, labelright=False, labelbottom=False)
sttt.plot(color='cadetblue',bgcolor='ivory', fig=fig, orientation="horizontal")
plt.title('EQ: 2022:09:14:11:04; New Hebrides event')
plt.axvline(end_v,.1,.9,lw=1.5,ls='--',c='navy')
plt.axvline(end_v-PP+P,.1,.9,ls='--',lw=1.5,c='navy')
plt.axvline(end_v-PP+pP,.1,.9,ls='--',lw=1.25,c='indianred')
plt.axvline(end_v-PP+sP,.1,.9,ls='--',lw=1.25,c='indianred')
plt.axvline(end_v-PP+pPP,.1,.9,ls='--',lw=1.25,c='indianred')

# plt.axvline(end_v-PP+S660P,.1,.9,ls='-.',lw=1.5,c='seagreen')
# plt.axvline(end_v-PP+S660PP,.1,.9,ls='-.',lw=1.5,c='seagreen')



plt.show()
# fig_name=os.path.join(folder_path,'eq_{}_{}.png'.format(time_string,i))
plt.savefig('220914_scm_marked.jpg', dpi=300, bbox_inches='tight', pad_inches=0.1)
plt.close()
