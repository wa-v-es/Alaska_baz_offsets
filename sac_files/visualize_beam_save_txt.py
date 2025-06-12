import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from mpl_point_clicker import clicker
from mpl_interactions import zoom_factory, panhandler
import pygmt
from matplotlib import ticker
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from cmcrameri import cm
import glob as glob
import re
from obspy.core import UTCDateTime as utc
from IPython import get_ipython
from obspy.taup.taup_geo import calc_dist,calc_dist_azi
from obspy.taup import TauPyModel
import subprocess

####
def extract_gridnumber(filename):
    match = re.search(r'gridnum(\d+)_', filename)
    if match:
        return int(match.group(1))
    return None
####
def extract_datapackfile(grid_number,folder_datapack):
    for file_name in os.listdir(folder_datapack):
        if file_name.endswith('.txt'):
            grid_num=extract_gridnumber(file_name)
            if grid_num==grid_number:
                print(file_name)
                return file_name
###
def extract_region_from_grid(grid_number,grid_folder,type):
    file_pattern = os.path.join(grid_folder, 'xf_{}_grid_{}*.grd'.format(type,grid_number))

    for filename in glob.glob(file_pattern):
        # print(filename)
        grd_file = pygmt.load_dataarray(filename)
        grd_info=pygmt.grdinfo(filename,per_column=True)
        grd_info=grd_info.split()
        # print(grd_info)
        region_type=[float(grd_info[0])+10, float(grd_info[1])-10, float(grd_info[2]), float(grd_info[3]),float(grd_info[4]),float(grd_info[5])]
        return(grd_file,region_type)
######
def calc_tt(eq_lat,eq_long,st_lat,st_long,eq_depth):
    model = TauPyModel(model="ak135")
    dist=calc_dist(eq_lat,eq_long,st_lat,st_long,6400,0)
    try:
        arr_PP = model.get_travel_times(source_depth_in_km=eq_depth,distance_in_degree=dist,phase_list=["PP"])[0]
        arr_pP=model.get_travel_times(source_depth_in_km=eq_depth,distance_in_degree=dist,phase_list=["pP"])[0]
        arr_sP=model.get_travel_times(source_depth_in_km=eq_depth,distance_in_degree=dist,phase_list=["sP"])[0]
    # arr_P = model.get_travel_times(source_depth_in_km=ori.depth/1000,distance_in_degree=dist,phase_list=["P"])[0]
    except:
        print('some phase didnt arrive')
    return(arr_PP,arr_pP,arr_sP)

def extract_grid_list(grid_folder):
    file_pattern = os.path.join(grid_folder, 'xf_slow_grid_*.grd')
    files = glob.glob(file_pattern)

    # List to store the extracted grid numbers
    grid_numbers = []

    # Regular expression to extract numbers after 'grid_'
    pattern = re.compile(r'grid_(\d+)_')

    for filename in files:
        match = pattern.search(filename)
        if match:
            grid_number = int(match.group(1))  # Convert to integer if needed
            grid_numbers.append(grid_number)
    grid_numbers.sort()
    return(grid_numbers)

################# CHANGE HERE FOR DIFF GRID #'s and folders
# %reset -f

get_ipython().magic('reset -sf')
#103-107, 125, 127, 150,168 #210121_122304_PA

# 6, 76,77, 100-07, 109, 125,127-8, 150 # 210818_101005_PA
# 76, 77, 100-107, 147-150, 168-69# 230615_180628_PA
grid_number=57
plot_amp_factor=3.5
print('plot_amp_factor=',plot_amp_factor)
print('grid_number=',grid_number)

main_folder='/Users/keyser/Research/AK_all_stations/sac_files/211002_062917_PA_inc2_r2.5/'
folder_datapack=main_folder+'data_pack/'
grid_folder=main_folder+'grid_folder'
pick_folder=main_folder+'py_picks/'

grid_list=extract_grid_list(grid_folder)
# sys.exit()
#################
beam_deets=folder_datapack+extract_datapackfile(grid_number,folder_datapack)

#beam_deets '/Users/keyser/Research/sub_array_alaska/sac_files/230702_102743_PA/Datapack_20230702_1027_.05_.5Hz_60samps_Zcomp_WRHbase_gridnum2_num8_PP_Y_N_0.0_Y_-1.txt'
patterns = {
    "Origin": re.compile(r"Origin: (\d+) (\d+) (\d+) (\d+):(\d+)"),
    "ArrCen": re.compile(r"ArrCen la/lo/elv: (\d+\.\d+) (-?\d+\.\d+) (\d+)"),
    "ArrBaseStn": re.compile(r"ArrBaseStn: (\w+), grid la/lp (\d+), (-?\d+)"),
    "Event": re.compile(r"Event la/lo/dp: (-?\d+\.\d+) (-?\d+\.\d+) (\d+\.\d+)"),
    "Dist": re.compile(r"Dist: (\d+\.\d+)"),
    "Baz": re.compile(r"Baz \(Arr-Evt\): (\d+\.\d+)"),
    "Frequencies": re.compile(r"Frequencies: (\.\d+) - (\.\d+) Hz"),
    "TrcesSNR": re.compile(r"TrcesSNR mn,SD,min,max: (\d+\.\d+) (\d+\.\d+) (\d+\.\d+) (\d+\.\d+)"),
    "PredPP": re.compile(r"Pred PP \(prem\) time/U: (\d+\.\d+) (\d+\.\d+)")
}

# Initialize a dictionary to store the extracted values
deets = {
    "Origin": [],
    "ArrCen": [],
    "ArrBaseStn": [],
    "Event": [],
    "Dist": [],
    "Baz": [],
    "Frequencies": [],
    "TrcesSNR": [],
    "PredPP": []
}

# Read the file and match lines with the defined patterns
with open(beam_deets, 'r') as file:
    for line in file:
        for key, pattern in patterns.items():
            match = pattern.search(line)
            if match:
                # Handle ArrBaseStn separately to exclude non-numeric values
                if key == "ArrBaseStn":
                    # Convert only numeric values, excluding the first group (station name)
                    deets[key].extend(match.groups()[1:])
                else:
                    deets[key].extend(map(float, match.groups()))

# Convert numeric strings to floats for ArrBaseStn
deets["ArrBaseStn"] = [float(x) for x in deets["ArrBaseStn"]]

# Print the extracted values
for key, values in deets.items():
    print(f"{key}: {values}")
#############

#------------------------
# gmt blockmean $beam_slow_xyz -I$interp_slow $range_slow | gmt surface -G$beam_slow_grid -I$interp_slow $range_slow
# gmt grdimage -X$xoffset0 -Y$yoffset0 $beam_slow_grid $range_slow $frame -C$temps"TEMP_BEAM.cpt_$$" -B$time_border/2:"Slowness (s/deg)":Wsen:."Beam - $phase": -K -P > $outps
# gmt makecpt -C/Users/keyser/Documents/ScientificColourMaps8/lipari/lipari.cpt --COLOR_FOREGROUND=white -T0/$grid_max_cpt/$xf_grid_inc -I -Z > $temps"TEMP_XF.cpt_$$"
# gmt makecpt -Cpolar -T-$beam_grid_max/$beam_grid_max/$beam_grid_inc -Z -I > $temps"TEMP_BEAM.cpt_$$"


interp_slow=[0.1,0.05]
# interp_baz="0.1/0.5"

slow_grd,region_slow=extract_region_from_grid(grid_number,grid_folder,'slow')
print(region_slow)
#####
baz_grd,region_baz=extract_region_from_grid(grid_number,grid_folder,'baz')
print(region_baz)
####
for line in open(beam_deets,'r'):
    line=line.split()
    if line[6]=='Pred':
        PP_t_s=[float(line[10]),float(line[11])]

##### xarray Plotting
# sys.exit()
try:
    arr_PP,arr_pP,arr_sP=calc_tt(deets['Event'][0],deets['Event'][1],deets['ArrCen'][0],deets['ArrCen'][1],deets['Event'][2])
except:
    print('one of phases didnt arrive')
plt.ion()

fig = plt.figure(figsize=(16, 6))


ax1 = fig.add_axes([0.05, 0.05, 0.38, 0.85]) #[left, bottom, width, height]
ax2 = fig.add_axes([0.53, 0.05, 0.38, 0.85]) #[left, bottom, width, height]

slow_grd.plot(ax=ax1,cmap=cm.lipari_r,add_colorbar=False,vmin=0,vmax=region_slow[5]/plot_amp_factor,mouseover=True)
slow_grd.plot.contour(ax=ax1,cmap='Greys_r',linewidths=.65,add_colorbar=False,levels=np.linspace(region_slow[5]/8, region_slow[5]/plot_amp_factor, 4))

baz_grd.plot(ax=ax2,cmap=cm.lipari_r,add_colorbar=False,vmin=0,vmax=region_slow[5]/plot_amp_factor,mouseover=True)
baz_grd.plot.contour(ax=ax2,cmap='Greys_r',linewidths=.65,add_colorbar=False,levels=np.linspace(region_slow[5]/8, region_slow[5]/plot_amp_factor, 4))

ax1.grid(which='minor',axis='x',color='dimgrey', linestyle='--',linewidth=.65,alpha=.75)
ax1.grid(which='major',axis='both',color='dimgrey', linestyle='--',linewidth=.75,alpha=.75)

ax2.grid(which='minor',axis='x',color='dimgrey', linestyle='--',linewidth=.65,alpha=.75)
ax2.grid(which='major',axis='both',color='dimGrey', linestyle='--',linewidth=.75,alpha=.75)

ax1.xaxis.set_minor_locator(MultipleLocator(10))
ax1.xaxis.set_major_locator(MultipleLocator(30))
ax1.yaxis.set_minor_locator(MultipleLocator(.25))
ax1.yaxis.set_major_locator(MultipleLocator(1))
ax2.xaxis.set_minor_locator(MultipleLocator(10))
ax2.xaxis.set_major_locator(MultipleLocator(30))
ax2.yaxis.set_minor_locator(MultipleLocator(5))
ax2.yaxis.set_major_locator(MultipleLocator(10))
### add PP
ax1.scatter(deets['PredPP'][0],deets['PredPP'][1],marker='o',c='violet',s=55,edgecolors='white',zorder=10)
ax2.scatter(deets['PredPP'][0],0,marker='o',c='violet',s=50,edgecolors='white',zorder=10)
ax2.scatter(arr_pP.time,0,marker='o',c='violet',s=50,edgecolors='white',zorder=10)
ax2.scatter(arr_sP.time,0,marker='o',c='violet',s=50,edgecolors='white',zorder=10)

# sys.exit()
ax1.set_ylabel('Slowness (s/deg)')
ax1.set_xlabel('Time (s)')
ax2.set_ylabel('Bazi (deg)')
ax2.set_xlabel('Time (s)')
# plt.show()
# sys.exit()

#####
# fig.text(0.42, -0.05, 'Time (s)',fontsize=12, ha='center', va='center')
ax3=fig.add_axes([0.4, .95, .3, 0.035]) #[left, bottom, width, height]
sm = plt.cm.ScalarMappable(norm=plt.Normalize(vmin=0, vmax=region_slow[5]/plot_amp_factor),cmap=cm.lipari_r )
sm.set_array(np.arange(0,1.33))
cbar = plt.colorbar(sm,cax=ax3,orientation='horizontal')
#####

zoom_factory(ax1)
ph = panhandler(fig, button=2)
klicker = clicker(
   ax1,
   ["Q1", "Q2"],
   markers=["+", "x"], markersize=10,
   colors=['skyblue','thistle'] #thistle, crimson
)

sys.exit()

#####SECOND BIT TO COPY and PASTE
# zoom_factory(ax2)
slow_click=klicker.get_positions()
print(klicker.get_positions())
klicker = clicker(
   ax2,
   ["Q1", "Q2"],
   markers=["+", "x"], markersize=10,
   colors=['skyblue','thistle'])

##############
sys.exit()
##############Third BIT TO COPY and PASTE
baz_click=klicker.get_positions()
print(baz_click)
max_slow_global=slow_grd.max().item()
slow_xf_pick=slow_grd.where((slow_grd.x > slow_click['Q1'][0][0]) & (slow_grd.x < slow_click['Q1'][1][0]), drop=True)
max_index_s = slow_xf_pick.argmax().item()
max_coords_slow = np.unravel_index(max_index_s, slow_xf_pick.shape) ##gets x,y position of Z max in slow_xf_pick
xf_pick_slow=[slow_xf_pick.x[max_coords_slow[1]].item(), slow_xf_pick.y[max_coords_slow[0]].item(), slow_xf_pick.max().item()]
print('xf slow time, slow, xf val=',xf_pick_slow)
###same for baz
baz_xf_pick=baz_grd.where((baz_grd.x > baz_click['Q1'][0][0]) & (baz_grd.x < baz_click['Q1'][1][0]), drop=True)
max_index_b = baz_xf_pick.argmax().item()
max_coords_baz = np.unravel_index(max_index_b, baz_xf_pick.shape) ##gets x,y position of Z max in baz_xf_pick
xf_pick_baz=[baz_xf_pick.x[max_coords_baz[1]].item(), baz_xf_pick.y[max_coords_baz[0]].item(), baz_xf_pick.max().item()]
print('xf baz time, baz, xf val=',xf_pick_baz)
#######
##check if the xf_pick_slow time and xf_pick_baz are within 2 sec of each other
if abs(xf_pick_slow[0] - xf_pick_baz[0]) > 2:
    raise ValueError(f"slow_pick and baz_pick are not within 2 sec of each other.")
else:
    print(f"slow_pick and baz_pick are within 2 sec of each other.")

# sys.exit()
utc_dt=''.join(str(int(x)) for x in deets['Origin'])

fig_name=pick_folder+'picks_gridnum_{}_{}_{}.jpg'.format(grid_number,utc_dt,'AK')
plt.savefig(fig_name,dpi=300,bbox_inches='tight', pad_inches=0.1)

# Write the extracted values (deets) to a new file in the specified format
outfile=pick_folder+'grid_num_{}_{}_{}_PICKS_amp_f_{}.dat'.format(grid_number,utc_dt,'AK',plot_amp_factor)
with open(outfile, 'w') as file:
    if 'ArrCen' in deets:
        file.write(f"arrayLLN {deets['ArrCen'][0]:.4f} {deets['ArrCen'][1]:.4f} AK\n")
    if "Event" in deets:
        file.write(f"eventLLD {deets['Event'][0]:.4f} {deets['Event'][1]:.4f} {deets['Event'][2]}\n")
    if "TrcesSNR" in deets:
        file.write(f"MAXFampPPR {deets['TrcesSNR'][3]:.2f} MAXFampP 0.00\n")
    if "Baz" in deets and "Dist" in deets:
        file.write(f"gcp {deets['Baz'][0]:.1f} dist {deets['Dist'][0]:.1f}\n")
    ##
    #C1-'SRC_LAT' C2-'SRC_LON' C3-'SRC_DEP' C4-'REC_LAT' C5-'REC_LON' C6-'DIST' C7-'BAZ' C8-'SCAT_TIME' C9-'SCAT_SLOW' C10-'SCAT_BAZ' C11-'ABS_BAZ' C12-'SNR_BEAM' C13-'PNR_BEAM' C14-'SIG_MAX' C15-'NOI_MAX' C16-'PREC_MAX'# 12xf_av_SNR, 13xf_max_SNR, 15filt_av_SNR, 16filt_max_SNR#
    #-8.2600 124.9300 10.0 -19.94 134.34 14.8 -39.4 2053.50 1.7 -61 80 1.27906 0.0923374 55.061 3.14051 0.356627
    baz_abs=xf_pick_baz[1]+deets['Baz'][0]
    file.write(f"{deets['Event'][0]:.4f} {deets['Event'][1]:.4f} {deets['Event'][2]} {deets['ArrCen'][0]:.4f} {deets['ArrCen'][1]:.4f} {deets['Dist'][0]:.1f} {deets['Baz'][0]:.1f} {xf_pick_slow[0]:.2f} {xf_pick_slow[1]:.1f} {xf_pick_baz[1]:.1f} {baz_abs:.1f} {deets['TrcesSNR'][3]:.2f} {deets['TrcesSNR'][3]:.2f} {deets['TrcesSNR'][3]:.2f} {deets['TrcesSNR'][3]:.2f} {deets['TrcesSNR'][3]:.2f}\n")

file.close()
plt.close()
#
