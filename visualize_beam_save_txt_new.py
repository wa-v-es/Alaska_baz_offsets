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
import matplotlib.colors as mcolors
from cmcrameri import cm
# import ScientificColourMaps8 as SCM8
import matplotlib.cm as cmm
from cmaptools import readcpt, joincmap, DynamicColormap
import glob as glob
import re
from obspy.core import UTCDateTime as utc
from IPython import get_ipython
from obspy.taup.taup_geo import calc_dist,calc_dist_azi
from obspy.taup import TauPyModel
import subprocess
from matplotlib.widgets import CheckButtons
import xarray as xr
import scipy.signal as sci
import warnings
import scipy.stats as stats
from scipy.stats import norm, skewnorm, kurtosis
from matplotlib.colors import ListedColormap
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
                # print(file_name)
                return file_name
###
def extract_region_from_grid(grid_number,grid_folder,type,beam_type):
    file_pattern = os.path.join(grid_folder, '{}_{}_grid_{}*.grd'.format(beam_type,type,grid_number))

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
    arr_PP=arr_pP=arr_sP=arr_pPP=float('nan')

    try:
        arr_P = model.get_travel_times(source_depth_in_km=eq_depth,distance_in_degree=dist,phase_list=["P"])[0]
    except:
        arr_P = model.get_travel_times(source_depth_in_km=eq_depth,distance_in_degree=dist,phase_list=["Pdiff"])[0]
    try:
        arr_PP = model.get_travel_times(source_depth_in_km=eq_depth,distance_in_degree=dist,phase_list=["PP"])[0]
    except:
        print('PP phase didnt arrive')
    try:
        arr_pP=model.get_travel_times(source_depth_in_km=eq_depth,distance_in_degree=dist,phase_list=["pP"])[0]
    except:
        arr_pP=model.get_travel_times(source_depth_in_km=eq_depth,distance_in_degree=dist,phase_list=["pPdiff"])[0]
    try:
        arr_sP=model.get_travel_times(source_depth_in_km=eq_depth,distance_in_degree=dist,phase_list=["sP"])[0]
    except:
        arr_sP=model.get_travel_times(source_depth_in_km=eq_depth,distance_in_degree=dist,phase_list=["sPdiff"])[0]

    arr_pPP=model.get_travel_times(source_depth_in_km=eq_depth,distance_in_degree=dist,phase_list=["pPP"])[0]

    return(arr_P,arr_PP,arr_pP,arr_sP,arr_pPP)

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

##
def set_locators(ax, axis_type='default'):
    if axis_type == 'slow':
        ax.xaxis.set_minor_locator(MultipleLocator(10))
        ax.xaxis.set_major_locator(MultipleLocator(30))
        ax.yaxis.set_minor_locator(MultipleLocator(0.25))
        ax.yaxis.set_major_locator(MultipleLocator(1))
    elif axis_type == 'baz':
        ax.xaxis.set_minor_locator(MultipleLocator(10))
        ax.xaxis.set_major_locator(MultipleLocator(30))
        ax.yaxis.set_minor_locator(MultipleLocator(5))
        ax.yaxis.set_major_locator(MultipleLocator(10))
    else:
        raise ValueError("Unsupported axis type: {}. Use 'default' or 'alt'.".format(axis_type))

def extract_grid_nums(main_folder):
    pattern = os.path.join(main_folder, '*.jpg')

    # List to store the gridnum values
    gridnum_list = []
    for file_path in glob.glob(pattern):
    # Extract the filename from the file path
        filename = os.path.basename(file_path)

        # Use regex to find the number after 'gridnum'
        match = re.search(r'gridnum(\d+)', filename)
        if match:
            gridnum_list.append(int(match.group(1)))

    return(gridnum_list)

def get_moving_avgs(xarray, time_step=5,overlap=2.5):

    z=baz_grd
    window_size = time_step
    step_size = overlap

    num_x_slices = len(z.x) #total x slices
    max_values = []
    y_values=[]
    midpoints = []
    i = 0
    while i < num_x_slices - window_size + 1:
        # Extract the current window slice
        window_slice = z.isel(x=slice(int(i), int(i + window_size)))

        # Find the max value in the current window
        max_value = window_slice.max().item()
        max_values.append(max_value)

        midpoint = z.x[int(i + window_size // 2)].item()
        midpoints.append(midpoint)

        i += step_size # Move the window by the step size

    max_values = np.array(max_values)

    # If you want to align it with the coordinates of the original x dimension
    # midpoints = z.x.isel(x=slice((window_size//2), int(num_x_slices - window_size + 1), (step_size)))

    # Create a new xarray.DataArray with the max values
    result = xr.DataArray(max_values, coords=[midpoints], dims=['x'])

    return(result)

def get_max_Z(xarray, time_step=5,overlap=2.5):

    z=xarray
    window_size = time_step
    step_size = overlap
    num_x_slices = len(z.x)

    midpoints = []
    max_values = []
    y_values = []

    for start_idx in np.arange(0, num_x_slices - window_size + 1, step_size):
        slice_x = z.isel(x=slice(int(start_idx), int(start_idx + window_size)))
        max_val = slice_x.max(dim=['y', 'x'])
        max_idx_y = slice_x.where(slice_x == max_val, drop=True).y
        if len(max_idx_y) >1: #sometimes, the index returned more than one position..
            max_idx_y=max_idx_y[0]

        # Store the midpoint x, max z value, and corresponding y value
        midpoints.append(slice_x.coords['x'].mean().item())
        max_values.append(max_val.values.max())
        y_values.append(max_idx_y.item())

    # Convert lists to numpy arrays for easier manipulation if needed
    midpoints = np.array(midpoints)
    max_values = np.array(max_values)
    y_values = np.array(y_values)

    return(midpoints,y_values,max_values)

#
def areEqual(arr1, arr2):
    N = len(arr1)
    M = len(arr2)

    # If lengths of array are not
    # equal means array are not equal
    if (N != M):
        return False

    # Sort both arrays
    arr1.sort()
    arr2.sort()

    # Linearly compare elements
    for i in range(0, N):
        if (arr1[i] != arr2[i]):
            return False

        # If all elements were same.
        return True

###
def extract_max_coher_clicks(grd,time_st,time_end):
    #for a grd array, extracts the max coherence for the manually
    # clicked crosses for ray tracing.
    slow_xf_pick=grd.where((grd.x > time_st) & (grd.x < time_end), drop=True)
    max_index_s = slow_xf_pick.argmax().item()
    max_coords_slow = np.unravel_index(max_index_s, slow_xf_pick.shape) ##gets x,y position of Z max in slow_xf_pick
    xf_pick_slow=[slow_xf_pick.x[max_coords_slow[1]].item(), slow_xf_pick.y[max_coords_slow[0]].item(), slow_xf_pick.max().item()]
    return xf_pick_slow
###

def get_contour_around_max(grd,x_max,window_size,percent):
    #percent in (0,1) #window size in sec
    x_min = x_max - window_size
    x_max = x_max + window_size

    window_data = grd.sel(x=slice(x_min, x_max))

    # Flatten the data array within the window
    flattened = window_data.values.flatten()
    # Sort the values in descending order
    sorted_values = np.sort(flattened)[::-1]
    # Calculate the index for the top 5% of the highest values
    top_5_percent_index = int(len(sorted_values) * percent)
    # Get the threshold value for the top 5% of the highest values
    threshold_value = sorted_values[top_5_percent_index]
    # Mask the original array to keep only values above the threshold within the window
    masked_array = window_data.where(window_data >= threshold_value)
    # Extract the y values where z values are within the top 5%
    y_values = masked_array['y'].values[masked_array.notnull().any(dim='x')]

    return masked_array,y_values
####
def get_peaks_grd(grd):

    midpoints,y_values,max_values=get_max_Z(grd,time_step=5,overlap=2.5)
    # midpoints_slow,y_values_slow,max_values_slow=get_max_Z(slow_grd_curtail,time_step=5,overlap=2.5)
    midpoints = np.array(midpoints)
    max_values = np.array(max_values)
    y_values = np.array(y_values)

    # get peaks in max_values (baz and slow grds)
    indexes, dict = sci.find_peaks(np.array(max_values),height=.25*np.max(max_values),prominence=.1*np.max(max_values))
    # indexes_slow, dict_slow = sci.find_peaks(np.array(max_values_slow),height=.25*np.max(max_values_slow),prominence=.1*np.max(max_values_slow))
    return midpoints,max_values,y_values,indexes
# %reset -f

###
cptfile='/Users/keyser/Documents/cmaptools/Andy_GIlmore_2.cpt'#
# cptfile='/Users/keyser/Documents/cmaptools/blue-tan-d14.cpt'#purple-orange-d09.cpt
cptfile_='/Users/keyser/Documents/cmaptools/blue-yellow.cpt'#
cptfile='/Users/keyser/Documents/cmaptools/green-purple-d09.cpt'#

# cptfile='/Users/keyser/Documents/cmaptools/voxpop.cpt'
cmap_try= readcpt(cptfile)
cmap_slow= readcpt(cptfile_)

##
## comment this if running from terminal
# get_ipython().magic('reset -sf')


folder_pattern = "sac_files_with_P/*_inc2_r2.5"
# folder_pattern = "sac_files/*_inc2_r2.5"
folder_pattern = "sac_files/*_*A_wrong"

folder_pattern = "/Users/keyser/Research/TA_arrays/sac_files/*_*A"



matching_folders = glob.glob(folder_pattern)

##
max_mean_gl=[]
# matching_folders=['120101_052755_PA_inc2_r2.5','120428_100807_PA_inc2_r2.5']
# matching_folders=['sac_bf_noise']
# matching_folders=['sac_files_with_P/230402_180411_PA_inc2_r2.5']
# matching_folders=['sac_files/200717_025022_PA_inc2_r2.5']


# sys.exit()
just_P=False

for folder in matching_folders:
    # main_folder='/Users/keyser/Research/AK_all_stations/sac_files/'+folder+'/'
    # main_folder='/Users/keyser/Research/AK_all_stations/'+folder+'/'
    # main_folder='/Users/keyser/Research/axisem/moho_3d/moho_dip_prllN_10s_dir_no_smooth/simu3D/output/stations/AK_81/'+folder+'/'
    # main_folder='/Users/keyser/Research/axisem/loyalty_isl/output_10sec_expl_source/stations/AK_81/'+folder+'/'
    main_folder=folder+'/' # for TA

    folder_datapack=main_folder+'data_pack/'
    grid_folder=main_folder+'grid_folder'
    pick_folder=main_folder+'py_picks/'
    py_figs=main_folder+'py_figs/'

    print('Main folder:',main_folder)
    gridnum_list=extract_grid_nums(main_folder)
    gridnum_list.sort()

    grid_baz_offset=[]
    grid_baz_offset_low_slow=[]
    grid_baz_offset_high_slow=[]

    for grid_number in gridnum_list:
    # for grid_number in [23]:

        plot_amp_factor=3
        # plot_amp_factor=10

        print('plot_amp_factor=',plot_amp_factor)
        print('grid_number=',grid_number)

        grid_list=extract_grid_list(grid_folder)

        #################
        beam_deets=folder_datapack+extract_datapackfile(grid_number,folder_datapack)
        print(beam_deets,'\n')
        #beam_deets '/Users/keyser/Research/sub_array_alaska/sac_files/230702_102743_PA/Datapack_20230702_1027_.05_.5Hz_60samps_Zcomp_WRHbase_gridnum2_num8_PP_Y_N_0.0_Y_-1.txt'
        patterns = {
            "Origin": re.compile(r"Origin: (\d+) (\d+) (\d+) (\d+):(\d+)"),
            "ArrCen": re.compile(r"ArrCen la/lo/elv: (\d+\.\d+) (-?\d+\.\d+) (\d+) Nst:(\d+)"),
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
        print('-------------\n')
        print(deets["ArrCen"])

        #############
        #------------------------
        interp_slow=[0.1,0.05]
        # interp_baz="0.1/0.5"

        slow_grd,region_slow=extract_region_from_grid(grid_number,grid_folder,'slow','xf')
        # print(region_slow)
        #####
        baz_grd,region_baz=extract_region_from_grid(grid_number,grid_folder,'baz','xf')
        # print(region_baz)
        ####
        slow_grd_bm,region_slow_bm=extract_region_from_grid(grid_number,grid_folder,'slow','beam')
        # print(region_slow_bm)
        #####
        baz_grd_bm,region_baz_bm=extract_region_from_grid(grid_number,grid_folder,'baz','beam')
        ##
        if np.size(slow_grd) != np.size(baz_grd):
            warnings.warn('Baz and slow mid-points array unequal :/; skipping this grid')
            continue

        for line in open(beam_deets,'r'):
            line=line.split()
            if line[6]=='Pred':
                PP_t_s=[float(line[10]),float(line[11])]

        ###
        max_position = baz_grd.argmax(dim=['y', 'x'])
        max_position_slow = slow_grd.argmax(dim=['y', 'x'])

        # Extract the indices for 'y' and 'x'
        y_max = baz_grd['y'][max_position['y']].item()
        x_max = baz_grd['x'][max_position['x']].item()

        x_max_slow = slow_grd['x'][max_position_slow['x']].item()
        y_max_slow = slow_grd['y'][max_position_slow['y']].item()


        print("----------------------\n")
        print(f"Max baz_grd for grid {grid_number} is at time: {x_max:.2f}s, baz: {y_max}")
        print("----------------------\n")

        try:
            arr_P,arr_PP,arr_pP,arr_sP,arr_pPP=calc_tt(deets['Event'][0],deets['Event'][1],deets['ArrCen'][0],deets['ArrCen'][1],deets['Event'][2])
        except:
            print('one or more phases didnt arrive')
        # sys.exit()
        # plt.ion()

        fig = plt.figure(figsize=(14, 9))

        # sys.exit()
        ax1 = fig.add_axes([0.05, 0.05, 0.38, 0.52]) #[left, bottom, width, height]
        ax2 = fig.add_axes([0.53, 0.05, 0.38, 0.52]) #[left, bottom, width, height]
        ax4 = fig.add_axes([0.05, 0.65, 0.38, 0.32])
        ax5 = fig.add_axes([0.53, 0.65, 0.38, 0.32])

        # ###check box stuff
        # check_ax = fig.add_axes([0.45, 0.4, 0.05, 0.05])
        # labels = ['slow', 'baz']
        # visibility = [False, False]  # Default visibility
        #
        # check = CheckButtons(check_ax, labels, visibility)
        # for label in check.labels:
        #     label.set_fontsize(11)  # Increase font size
        #     label.set_x(0.4)
        # for rect in check.rectangles:
        #     rect.set_width(0.2)  # Increase width of the tick box
        #     rect.set_height(0.2)
        # ###check box stuff
        ####
        cmap_lip=cm.lipari_r
        colA = cmap_lip(np.arange(cmap_lip.N))
        colA[:,-1] = np.linspace(0.7, 1, cmap_lip.N) #replaces the last column (alpha channel) in colA with values from 0.75 to 1, creating a gradient in transparency across the color map.

        sm_alpha = ListedColormap(colA)

        # the following adds white at the start of a color map.
        # cmapp=plt.get_cmap(sm_alpha)
        # colors = cmapp(np.linspace(0, 1, 256))
        # colors[0:2] = [1, 1, 1, 1]  # White color
        # new_lipari = plt.cm.colors.ListedColormap(colors)

        slow_grd.plot(ax=ax1,cmap=sm_alpha,add_colorbar=False,vmin=0,vmax=region_baz[5]/plot_amp_factor,mouseover=True)
        slow_grd.plot.contour(ax=ax1,cmap='Greys_r',linewidths=.65,add_colorbar=False,levels=np.linspace(region_baz[5]/8, region_baz[5]/plot_amp_factor, 4))

        baz_grd.plot(ax=ax2,cmap=sm_alpha,add_colorbar=False,vmin=0,vmax=region_baz[5]/plot_amp_factor,mouseover=True)
        baz_grd.plot.contour(ax=ax2,cmap='Greys_r',linewidths=.65,add_colorbar=False,levels=np.linspace(region_baz[5]/8, region_baz[5]/plot_amp_factor, 4))
        #cm.lipari_r
        slow_grd_bm.plot(ax=ax4,cmap=cm.broc,add_colorbar=False,mouseover=True)
        # slow_grd_bm.plot.contour(ax=ax4,cmap='Greys_r',linewidths=.65,add_colorbar=False,levels=np.linspace(region_slow[5]/8, region_slow[5]/plot_amp_factor, 4))
        baz_grd_bm.plot(ax=ax5,cmap=cm.broc,add_colorbar=False,mouseover=True)

        # gets max 5% around max in slow/baz grids!!!
        # grd,x_max,window_size,percent=slow_grd,x_max_slow,5,.05

        grd_5_slow,slow_5_vals=get_contour_around_max(slow_grd,x_max_slow,5,.05)
        grd_5_baz,baz_5_vals=get_contour_around_max(baz_grd,x_max,5,.05)


        ax1.scatter([x_max_slow,x_max_slow],[slow_5_vals.min(),slow_5_vals.max()],marker='_',s=100,c='white',zorder=10)
        ax2.scatter([x_max,x_max],[baz_5_vals.min(),baz_5_vals.max()],marker='_',s=100,c='white',zorder=10)

        # ax1.scatter(x_max_slow,slow_5_vals.max(),marker='_',s=80,c='white',zorder=10)

        # TEXT FOR MIN/MAX baz slowness!
        # ax1.text(region_baz[0]+5, 10.1, 'max/min/diff slow (5%) around max coh. ($\pm$5 sec): {:.1f}/{:.1f}/{:.1f}'.format(slow_5_vals.max(),slow_5_vals.min(),(slow_5_vals.max()-slow_5_vals.min())),c='navy',size=11,weight='roman')
        # ax2.text(region_baz[0]+5, 51, 'max/min/diff baz (5%) around max coh. ($\pm$5 sec): {:.1f}/{:.1f}/{:.1f}'.format(baz_5_vals.max(),baz_5_vals.min(),(baz_5_vals.max()-baz_5_vals.min())),c='navy',size=11,weight='roman')

        ######
        ax1.grid(which='minor',axis='x',color='dimgrey', linestyle='--',linewidth=.65,alpha=.75)
        ax1.grid(which='major',axis='both',color='dimgrey', linestyle='--',linewidth=.75,alpha=.75)

        ax2.grid(which='minor',axis='x',color='dimgrey', linestyle='--',linewidth=.65,alpha=.75)
        ax2.grid(which='major',axis='both',color='dimGrey', linestyle='--',linewidth=.75,alpha=.75)

        ### plot horiz line at baz max_baz
        ax2.axhline(y=0, color='darkred', linestyle='--')
        ax2.scatter(x_max,y_max,marker='d',c='darkred',s=55,edgecolors='white',zorder=10)
        ax2.text(region_baz[0]+10, 40, 'max ({}) at {}$^\circ$ bazi'.format(int(baz_grd.max().item()),y_max),c='darkred',size=12,weight='roman',bbox={'facecolor': 'white', 'alpha': 0.85, 'pad': 1.5})
        ##
        set_locators(ax1, 'slow')
        set_locators(ax4, 'slow')
        set_locators(ax2, 'baz')
        set_locators(ax5, 'baz')

        # for phase in [arr_P,arr_pP,arr_sP,arr_PP]:
        for phase in [arr_pP,arr_sP,arr_PP]:

            if 'diff' in phase.name:
                ax1.scatter(phase.time,phase.ray_param*0.0174533,marker='o',c='CORNFLOWERBLUE',s=50,edgecolors='white',zorder=10)
                ax1.text(phase.time, 1.5+phase.ray_param * 0.0174533, phase.name, bbox={'facecolor': 'white', 'alpha': 0.85, 'pad': 1.5},fontsize=9.5,c='CORNFLOWERBLUE', rotation='vertical',ha='center')
            else:
                ax1.scatter(phase.time,phase.ray_param*0.0174533,marker='o',c='violet',s=50,edgecolors='white',zorder=10)
                ax1.text(phase.time, 1.5+phase.ray_param * 0.0174533, phase.name, bbox={'facecolor': 'white', 'alpha': 0.85, 'pad': 1.5},fontsize=9.5,c='violet', rotation='vertical',ha='center')

        plt.rcParams['axes.labelsize'] = 15
        ax1.set_ylabel('Slowness (s/$^\circ$)')
        ax1.set_xlabel('Time (s)')
        # ax2.set_xticklabels([])
        ax2.set_ylabel('Bazi ($^\circ$)')
        ax2.set_xlabel('Time (s)')
        ax4.set_ylabel('Slowness (s/$^\circ$)')
        ax4.set_xlabel('')
        ax5.set_ylabel('Bazi ($^\circ$)')
        ax5.set_xlabel('')
        #####
        # fig.text(0.42, -0.05, 'Time (s)',fontsize=12, ha='center', va='center')
        ax3=fig.add_axes([0.94, .27, .01, 0.3]) #[left, bottom, width, height]
        sm = plt.cm.ScalarMappable(norm=plt.Normalize(vmin=0, vmax=region_baz[5]/plot_amp_factor),cmap=cm.lipari_r )
        sm = plt.cm.ScalarMappable(norm=plt.Normalize(vmin=0, vmax=region_baz[5]/plot_amp_factor),cmap=sm_alpha )

        sm.set_array(np.arange(0,1.33))
        cbar = plt.colorbar(sm,cax=ax3,orientation='vertical')
        cbar.set_label('Coherence', rotation=90, labelpad=2,loc='center')
        #####
        # start bit for click boxes

        # zoom_factory(ax1)
        # ph = panhandler(fig, button=2)
        # klicker = clicker(
        #    ax1,
        #    ["Q1", "Q2"],
        #    markers=["+", "x"], markersize=10,
        #    colors=['skyblue','thistle'] #thistle, crimson
        # )

        # end bit for click boxes
        utc_dt=''.join(str(int(x)) for x in deets['Origin'])

        ############# getting max Z and corresponding baz for each time slice
        slow_grd_curtail=slow_grd.where((slow_grd.x > arr_pP.time - 10) & (slow_grd.x < arr_pPP.time + 20), drop=True)
        baz_grd_curtail=baz_grd.where((baz_grd.x > arr_pP.time - 10) & (baz_grd.x < arr_pPP.time + 20), drop=True)

        slow_grd_justP=slow_grd.where((slow_grd.x < arr_pP.time - 10) , drop=True)
        baz_grd_justP=baz_grd.where((baz_grd.x < arr_pP.time - 10), drop=True)

        #
        midpoints,max_values,y_values,indexes=get_peaks_grd(baz_grd_curtail)
        midpoints_slow,max_values_slow,y_values_slow,indexes_slow=get_peaks_grd(slow_grd_curtail)

        if just_P:
            P_midpoints,P_max_values,P_y_values,P_indexes=get_peaks_grd(baz_grd_justP)
            P_midpoints_slow,P_max_values_slow,P_y_values_slow,P_indexes_slow=get_peaks_grd(slow_grd_justP)

        ##### finding if the energy is everywhere or not

        avg_coherence = baz_grd_curtail.mean(dim=['x', 'y'])
        std_coherence= baz_grd_curtail.std(dim=['x', 'y'])
        z_values_coh = baz_grd_curtail.values.flatten()
        # skew_coherence = stats.skew(z_values_coh)
        # kurtosis_coherence = stats.kurtosis(z_values_coh)
        max_mean= baz_grd_curtail.max(dim=['x', 'y'])/baz_grd_curtail.mean(dim=['x', 'y'])

        max_mean_gl.append(round(max_mean.item(),2))

        ## dividing peaks with slow > 6 and less than 6
        # slowness=6 was chosen coz mostly PP comes at slow > 7 and based on trends
        # 6 seemed like an appropriate value to distinguis precursory phases and PP.
        indexes_less6=[]
        indexes_gr6=[]
        # y_values_slow contains peaks in slowness
        # indexes contains index of selected peaks. Thus, y_values_slow[i] gets slowness vals of peaks in baz offsets.
        if len(indexes) ==0:
            print('-----------------\n')
            warnings.warn('Baz and slow mid-points array unequal :/; skipping this grid')
            print('-----------------\n')
            continue
        for i in indexes:
            if y_values_slow[i] > 5.99:
                indexes_gr6.append(i)
            else:
                indexes_less6.append(i)

        # Create discrete colormap
        num_bins = 8 # for discrete colour bar

        norm = mcolors.BoundaryNorm(np.linspace(-4,4, num_bins + 1), cmap_try.N)
        norm_slow = mcolors.BoundaryNorm(np.linspace(1,9, num_bins + 1), cmap_slow.N)


        ax9=fig.add_axes([0.91, 0.65, 0.07, 0.32])  #[left, bottom, width, height]
        # ax9.plot(y_skewed,x, linestyle='-',color='orchid')#label=f"Skewed (skew={skew_coherence:.1f})",
        ax9.hist(z_values_coh, density=True, bins='auto', histtype='stepfilled', alpha=0.65,color='orchid')
        ax9.axvline(avg_coherence, color='slateblue', linestyle='--', lw=1.25)#label=f'{avg_coherence.item():.1f}',
        ax9.set_yticklabels([])
        ax9.set_ylim((0, None))
        ax9.set_xlim(0,baz_grd_curtail.max().item()/5 )
        # plt.ylabel('coherence')
        plt.ylabel('PDF')
        fig.text(0.95, .95, f'mean={avg_coherence.item():.1f}',fontsize=10, ha='center', va='center',color='slateblue')
        fig.text(0.95, .925, f'std={std_coherence.item():.1f}',fontsize=10, ha='center', va='center')
        fig.text(0.95, .90, f'm/m={max_mean.item():.1f}',fontsize=10, ha='center', va='center')

        # plt.legend()
        ax9.set_yticks([])
        ax9.set_xticks([])

        plt.grid(True)
        # plt.show()


        ### Start of plotting peaks ax7 and ax8
        # ax7 is baz
        ax7=fig.add_axes([0.53, 0.05, 0.475, 0.095],sharex=ax2) #[left, bottom, width, height]
        ax7.set_facecolor("whitesmoke")
        ax7.grid(True,alpha=.25)
        ax7.plot(midpoints, max_values, '-',lw=.25,c='black',alpha=.65)
        ax7.axvline(arr_pP.time - 10, color='hotpink', linestyle='--', lw=1.85)#label=f'{avg_coherence.item():.1f}',
        #ax8 is slow
        ax8=fig.add_axes([0.05, 0.05, 0.475, 0.095],sharex=ax1) #[left, bottom, width, height]
        ax8.set_facecolor("whitesmoke")
        ax8.grid(True,alpha=.25)
        ax8.plot(midpoints_slow, max_values_slow, '-',lw=.25,c='black',alpha=.65)
        ax8.axvline(arr_pP.time - 10, color='hotpink', linestyle='--', lw=1.85)#label=f'{avg_coherence.item():.1f}',

        #plot all mid points and color by baz
        scatter_7 = ax7.scatter(midpoints, max_values, c=y_values, cmap=cmap_try, norm=norm, edgecolor='white', s=15,alpha=.88,linewidth=.15)
        scatter_8 = ax8.scatter(midpoints_slow, max_values_slow, c=y_values_slow, cmap=cmap_slow, norm=norm_slow, edgecolor='white', s=15,alpha=.88,linewidth=.15)
        if just_P:
            ax7.plot(P_midpoints, P_max_values, '-',lw=.25,c='black',alpha=.65)
            ax8.plot(P_midpoints_slow, P_max_values_slow, '-',lw=.25,c='black',alpha=.65)
            scatter_7_P = ax7.scatter(P_midpoints, P_max_values, c=P_y_values, cmap=cmap_try, norm=norm, edgecolor='white', s=15,alpha=.88,linewidth=.15)
            scatter_8_P = ax8.scatter(P_midpoints_slow, P_max_values_slow, c=P_y_values_slow, cmap=cmap_slow, norm=norm_slow, edgecolor='white', s=15,alpha=.88,linewidth=.15)

        ### first if fro low slowness...thats' why P shenanigans
        if len(indexes_less6) ==0:
            print('-----------------\n')
            warnings.warn('nothing in slow < 6 :/; not plotting')
            print('-----------------\n')
        else:
            if max_mean.item() > 20:
                if just_P:
                    low_slow_combined=np.concatenate([y_values[indexes_less6], P_y_values[P_indexes]])
                else:
                    low_slow_combined=np.concatenate([y_values[indexes_less6]])

                grid_baz_offset_low_slow.append((grid_number,np.max(low_slow_combined),np.mean(low_slow_combined),np.std(low_slow_combined),deets["ArrCen"][0],deets["ArrCen"][1],deets["ArrCen"][2],deets["Event"][0],deets["Event"][1],deets["Event"][2],deets["Dist"][0],deets["Baz"][0],deets["ArrCen"][3],(baz_5_vals.max()-baz_5_vals.min()),(slow_5_vals.max()-slow_5_vals.min())))
            ax7.scatter(midpoints[indexes_less6], max_values[indexes_less6], marker='+',color='darkred',s=40,linewidth=1.25)
            ax8.scatter(midpoints_slow[indexes_less6], max_values_slow[indexes_less6], marker='+',color='darkred',s=40,linewidth=1.25)
            # for P..
            if just_P:
                ax7.scatter(P_midpoints[P_indexes], P_max_values[P_indexes], marker='+',color='darkred',s=40,linewidth=1.25)
                ax8.scatter(P_midpoints_slow[P_indexes_slow], P_max_values_slow[P_indexes_slow], marker='+',color='darkred',s=40,linewidth=1.25)

                ax7.text(region_slow[0]+5,1.16*max(max(max_values), max(P_max_values)),f'{np.mean(low_slow_combined):.2f};{np.std(low_slow_combined):.2f}',c='darkred',size=11,weight='roman',bbox={'facecolor': 'white', 'alpha': 0.95, 'pad': 1.75})
            else:
                ax7.text(region_slow[0]+5,1.16*max(max_values),f'{np.mean(low_slow_combined):.2f};{np.std(low_slow_combined):.2f}',c='darkred',size=11,weight='roman',bbox={'facecolor': 'white', 'alpha': 0.95, 'pad': 1.75})

        # second if for high slowness
        if len(indexes_gr6) ==0:
            print('-----------------\n')
            warnings.warn('nothing in slow > 6 :/; not plotting')
            print('-----------------\n')
        else:
            if max_mean.item() > 20:
                grid_baz_offset_high_slow.append((grid_number,np.max(y_values[indexes_gr6]),np.mean(y_values[indexes_gr6]),np.std(y_values[indexes_gr6]),deets["ArrCen"][0],deets["ArrCen"][1],deets["ArrCen"][2],deets["Event"][0],deets["Event"][1],deets["Event"][2],deets["Dist"][0],deets["Baz"][0],deets["ArrCen"][3],(baz_5_vals.max()-baz_5_vals.min()),(slow_5_vals.max()-slow_5_vals.min())))
            ax7.scatter(midpoints[indexes_gr6], max_values[indexes_gr6], marker='+',color='black',s=40,linewidth=1.25)
            ax8.scatter(midpoints[indexes_gr6], max_values[indexes_gr6], marker='+',color='black',s=40,linewidth=1.25)
            if just_P:
                ax7.text(region_baz[1]-45, 1.16*max(max(max_values), max(P_max_values)),f'{np.mean(y_values[indexes_gr6]):.2f};{np.std(y_values[indexes_gr6]):.2f}',c='black',size=11,weight='roman',bbox={'facecolor': 'white', 'alpha': 0.95, 'pad': 1.75})
            else:
                ax7.text(region_baz[1]-45, 1.16*max(max_values),f'{np.mean(y_values[indexes_gr6]):.2f};{np.std(y_values[indexes_gr6]):.2f}',c='black',size=11,weight='roman',bbox={'facecolor': 'white', 'alpha': 0.95, 'pad': 1.75})


        if max_mean.item() > 20:
            if just_P:
                all_slow_combined=np.concatenate([y_values[indexes], P_y_values[P_indexes]])
            else:
                all_slow_combined=np.concatenate([y_values[indexes]])
            grid_baz_offset.append((grid_number,y_max,np.mean(all_slow_combined),np.std(all_slow_combined),deets["ArrCen"][0],deets["ArrCen"][1],deets["ArrCen"][2],deets["Event"][0],deets["Event"][1],deets["Event"][2],deets["Dist"][0],deets["Baz"][0],deets["ArrCen"][3],(baz_5_vals.max()-baz_5_vals.min()),(slow_5_vals.max()-slow_5_vals.min())))
        # grid_slow_offset.append((grid_number,y_max_slow,np.mean(y_values_slow[indexes_slow]),np.std(y_values_slow[indexes_slow]),deets["ArrCen"][0],deets["ArrCen"][1],deets["ArrCen"][2]))

        cbar = plt.colorbar(scatter_7)
        cbar.set_label('Baz$^\circ$', fontsize=12)
        cbar.set_ticks([-4,-2,0,2,4])

        cbar_slow = plt.colorbar(scatter_8)
        cbar_slow.set_label('slow(s/$^\circ$)', fontsize=12)
        cbar_slow.set_ticks([1,3,5,7,9])
        # ax7.set_xticklabels([])
        ax7.set_yticklabels([])
        ax8.set_yticklabels([])
        ax7.set_yticks([])
        ax8.set_yticks([])

        ##
        # ax7.set_xlabel('Time (s)'), ax7.set_ylabel('Coherence')

        fig.text(0.15, .98, 'Grid #{} for eq: {}'.format(grid_number,utc_dt),fontsize=12, ha='center', va='center')
        fig_name=py_figs+'picks_gridnum_{}_{}_{}.jpg'.format(grid_number,utc_dt,'AK')

        plt.savefig(fig_name,dpi=300,bbox_inches='tight', pad_inches=0.1)
        plt.close('all')
        # sys.exit()


    plt.close()
    # print(grid_baz_offset)
    # sys.exit()

    ## contains grid number; baz at max coherence va; mean of peaks;std of peaks; array cen lat long elevation
    ## originally, it was folder[:13] as this script was within the sac folder.
    # np.savetxt('max_vals_coherence/{}_maxVals.txt'.format(folder[17:33]), grid_baz_offset, fmt='%i %.1f %.2f %.2f %.2f %.2f %.1f %.2f %.2f %.1f %.1f %.1f %i %.1f %.1f')
    # if len(grid_baz_offset_low_slow) !=0:
    #     np.savetxt('max_vals_coherence/{}_maxVals_low_slow.txt'.format(folder[17:33]), grid_baz_offset_low_slow, fmt='%i %.1f %.2f %.2f %.2f %.2f %.1f %.2f %.2f %.1f %.1f %.1f %i %.1f %.1f')
    # if len(grid_baz_offset_high_slow) !=0:
    #     np.savetxt('max_vals_coherence/{}_maxVals_high_slow.txt'.format(folder[17:33]), grid_baz_offset_high_slow, fmt='%i %.1f %.2f %.2f %.2f %.2f %.1f %.2f %.2f %.1f %.1f %.1f %i %.1f %.1f')

    # for TA
    folder_='/Users/keyser/Research/TA_arrays/sac_files/TA_max_vals_coherence'

    np.savetxt('{}/{}_maxVals.txt'.format(folder_,folder[-16:]), grid_baz_offset, fmt='%i %.1f %.2f %.2f %.2f %.2f %.1f %.2f %.2f %.1f %.1f %.1f %i %.1f %.1f')
    if len(grid_baz_offset_low_slow) !=0:
        np.savetxt('{}/{}_maxVals_low_slow.txt'.format(folder_,folder[-16:]), grid_baz_offset_low_slow, fmt='%i %.1f %.2f %.2f %.2f %.2f %.1f %.2f %.2f %.1f %.1f %.1f %i %.1f %.1f')
    if len(grid_baz_offset_high_slow) !=0:
        np.savetxt('{}/{}_maxVals_high_slow.txt'.format(folder_,folder[-16:]), grid_baz_offset_high_slow, fmt='%i %.1f %.2f %.2f %.2f %.2f %.1f %.2f %.2f %.1f %.1f %.1f %i %.1f %.1f')

    print('--------------------------------\n')
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
### extract times of max coherence for picked clicks
xf_pick_slow=[]
for i in range(0, len(slow_click['Q1']), 2):
    # Extract two rows at a time
    rows = slow_click['Q1'][i:i+2]
    xf_pick_slow.append(extract_max_coher_clicks(slow_grd,rows[0][0],rows[1][0]))
print(xf_pick_slow)
# print('xf slow time, slow, xf val=',xf_pick_slow)
###same for baz
xf_pick_baz=[]

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
    for picks in xf_pick_slow:
        #C1-'SRC_LAT' C2-'SRC_LON' C3-'SRC_DEP' C4-'REC_LAT' C5-'REC_LON' C6-'DIST' C7-'BAZ' C8-'SCAT_TIME' C9-'SCAT_SLOW' C10-'SCAT_BAZ' C11-'ABS_BAZ' C12-'SNR_BEAM' C13-'PNR_BEAM' C14-'SIG_MAX' C15-'NOI_MAX' C16-'PREC_MAX'# 12xf_av_SNR, 13xf_max_SNR, 15filt_av_SNR, 16filt_max_SNR#
        #-8.2600 124.9300 10.0 -19.94 134.34 14.8 -39.4 2053.50 1.7 -61 80 1.27906 0.0923374 55.061 3.14051 0.356627
        # baz_abs=xf_pick_baz[1]+deets['Baz'][0]
        baz_abs=deets['Baz'][0]

        file.write(f"{deets['Event'][0]:.4f} {deets['Event'][1]:.4f} {deets['Event'][2]} {deets['ArrCen'][0]:.4f} {deets['ArrCen'][1]:.4f} {deets['Dist'][0]:.1f} {deets['Baz'][0]:.1f} {picks[0]:.2f} {picks[1]:.2f} {0} {baz_abs:.1f} {deets['TrcesSNR'][3]:.2f} {deets['TrcesSNR'][3]:.2f} {deets['TrcesSNR'][3]:.2f} {deets['TrcesSNR'][3]:.2f} {deets['TrcesSNR'][3]:.2f}\n")

file.close()
plt.close()
#
#### bit to get grid_numbers in a folder from jpg FILES
