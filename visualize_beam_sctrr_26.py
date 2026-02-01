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
from matplotlib.gridspec import GridSpec
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

def calc_tt_lil(eq_lat,eq_long,st_lat,st_long,eq_depth):
    model = TauPyModel(model="ak135")
    dist=calc_dist(eq_lat,eq_long,st_lat,st_long,6400,0)
    lil_p=lil_s=float('nan')


    lil_p = model.get_travel_times(source_depth_in_km=eq_depth,distance_in_degree=dist,phase_list=["p"])[0]
    lil_s = model.get_travel_times(source_depth_in_km=eq_depth,distance_in_degree=dist,phase_list=["s"])[0]

    return(lil_p,lil_s)

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
        ax.yaxis.set_minor_locator(MultipleLocator(0.5))
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

folder_pattern = "sac_noise_latN_Ptime/*_inc2_r2.5"
# folder_pattern = "sac_files/*_inc2_r2.5"

matching_folders = glob.glob(folder_pattern)

##
max_mean_gl=[]
# matching_folders=['120101_052755_PA_inc2_r2.5','120428_100807_PA_inc2_r2.5']
matching_folders=['sac_noise_latN_Ptime']
matching_folders=['sac_files_with_P/220914_110406_PA_inc2_r2.5']
matching_folders=['200717_025022_PA_inc2_r2.5']

# sys.exit()
plt.rcParams.update({'font.size': 15})
for folder in matching_folders:
    main_folder='/Users/keyser/Research/AK_all_stations/sac_files/'+folder+'/'
    # main_folder='/Users/keyser/Research/AK_all_stations/'+folder+'/'
    # main_folder='/Users/keyser/Research/axisem/moho_3d/moho_dip_prllN_10s_dir_no_smooth/simu3D/output/stations/AK_81/'+folder+'/'

    folder_datapack=main_folder+'data_pack/'
    grid_folder=main_folder+'grid_folder'
    pick_folder=main_folder+'py_picks/'
    py_figs=main_folder+'py_figs_new/'

    print('Main folder:',main_folder)
    gridnum_list=extract_grid_nums(main_folder)
    gridnum_list.sort()

    grid_baz_offset=[]
    grid_baz_offset_low_slow=[]
    grid_baz_offset_high_slow=[]

    for grid_number in gridnum_list:
    # for grid_number in [80]:

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

        # -------------------------
        # Colormaps
        # -------------------------
        cmap_lip = cmm.PuBu
        colA = cmap_lip(np.arange(cmap_lip.N))
        sm_alpha = ListedColormap(colA)
        # -------------------------
        # Window grids in time
        # -------------------------
        slow_grd = slow_grd.where(
            (slow_grd.x > arr_sP.time - 10) & (slow_grd.x < arr_PP.time + 10),
            drop=True)
        baz_grd = baz_grd.where(
            (baz_grd.x > arr_sP.time - 10) & (baz_grd.x < arr_PP.time + 10),
            drop=True)
        # -------------------------
        # 5% contours around maxima
        # -------------------------
        grd_5_slow, slow_5_vals = get_contour_around_max(slow_grd, x_max_slow, 5, .05)
        grd_5_baz, baz_5_vals = get_contour_around_max(baz_grd, x_max, 5, .05)
        # -------------------------
        # Curtail grids (between sP and PP)
        # -------------------------
        slow_grd_curtail = slow_grd.where(
            (slow_grd.x > arr_sP.time + 20) & (slow_grd.x < arr_PP.time - 10),drop=True)
        baz_grd_curtail = baz_grd.where(
            (baz_grd.x > arr_sP.time + 20) & (baz_grd.x < arr_PP.time - 10),drop=True)

        V_max_curtail = baz_grd_curtail.max(dim=['x', 'y'])
        plot_amp_factor_curtail = 1
        # -------------------------
        # Peaks
        # -------------------------
        midpoints, max_values, y_values, indexes = get_peaks_grd(baz_grd_curtail)
        midpoints_slow, max_values_slow, y_values_slow, indexes_slow = get_peaks_grd(slow_grd_curtail)
        # -------------------------
        # Coherence statistics
        # -------------------------
        avg_coherence = baz_grd.mean(dim=['x', 'y'])
        std_coherence = baz_grd.std(dim=['x', 'y'])
        z_values_coh = baz_grd.values.flatten()
        max_mean = baz_grd.max(dim=['x', 'y']) / avg_coherence

        max_mean_gl.append(round(max_mean.item(), 2))

        # -------------------------
        # Discrete norms
        # -------------------------
        num_bins = 8
        norm = mcolors.BoundaryNorm(np.linspace(-4, 4, num_bins + 1), cmap_try.N)
        norm_slow = mcolors.BoundaryNorm(np.linspace(1, 9, num_bins + 1), cmap_slow.N)

        fig = plt.figure(figsize=(15, 8))
        # gs = GridSpec(
        #     nrows=9,
        #     ncols=13,
        #     figure=fig,
        #     height_ratios=[0.35, 3, 3, 0.2, 2.6, 2.6, 0.25, 1.6, 1.6],
        #     width_ratios=[1]*11 + [0.15, 0.15])
        #
        # ax1 = fig.add_subplot(gs[1:3, 0:5])   # slow main
        # ax2 = fig.add_subplot(gs[1:3, 6:11])  # baz main
        #
        # ax3 = fig.add_subplot(gs[0, 6:12])     # coherence colorbar
        #
        # ax4 = fig.add_subplot(gs[4:6, 0:5])   # slow curtailed
        # ax5 = fig.add_subplot(gs[4:6, 6:11])  # baz curtailed
        #
        # ax9 = fig.add_subplot(gs[1:3, 11])    # histogram (spans rows)
        #
        # ax8 = fig.add_subplot(gs[7:9, 0:5])   # slow peaks
        # ax7 = fig.add_subplot(gs[7:9, 6:11])  # baz peaks
        # [left, bottom, width, height]
        ax1 = fig.add_axes([0.07, 0.6, 0.38, 0.3])   # slow main
        ax2 = fig.add_axes([0.52, 0.6, 0.38, 0.3])   # baz main
        ax3 = fig.add_axes([0.7, .92, .15, 0.012])    # colorbar (coherence)
        ax4 = fig.add_axes([0.07, 0.27, 0.38, 0.27]) # slow curtailed
        ax5 = fig.add_axes([0.52, 0.27, 0.38, 0.27]) # baz curtailed
        ax7 = fig.add_axes([0.52, 0.1, 0.475, 0.15]) # baz peaks
        ax8 = fig.add_axes([0.07, 0.1, 0.475, 0.15]) # slow peaks
        ax9 = fig.add_axes([0.9, 0.6, 0.07, 0.3])    # histogram

        ## slowness main plot
        slow_grd.plot(
        ax=ax1, cmap=sm_alpha, add_colorbar=False,
        vmin=0, vmax=region_baz[5] / plot_amp_factor, mouseover=True)
        slow_grd.plot.contour(
            ax=ax1, cmap='Greys_r', linewidths=.65, add_colorbar=False,
            levels=np.linspace(region_baz[5] / 8, region_baz[5] / plot_amp_factor, 4))

        ax1.scatter([x_max_slow, x_max_slow],[slow_5_vals.min(), slow_5_vals.max()],
            marker='_', s=100, c='white', zorder=10)

        ax1.axhline(y=6, color='black', linestyle='--', lw=1.85)

        for phase in [arr_sP,arr_PP]:

            if 'diff' in phase.name:
                ax1.scatter(phase.time,phase.ray_param*0.0174533,marker='o',c='CORNFLOWERBLUE',s=50,edgecolors='white',zorder=10)
                ax1.text(phase.time, 1.5+phase.ray_param * 0.0174533, phase.name, bbox={'facecolor': 'white', 'alpha': 0.85, 'pad': 1.5},fontsize=14,c='CORNFLOWERBLUE', rotation='vertical',ha='center')
            else:
                ax1.scatter(phase.time,phase.ray_param*0.0174533,marker='o',c='violet',s=50,edgecolors='white',zorder=10)
                ax1.text(phase.time, 1.5+phase.ray_param * 0.0174533, phase.name, bbox={'facecolor': 'white', 'alpha': 0.85, 'pad': 1.5},fontsize=14,c='violet', rotation='vertical',ha='center')

        #### ax2: baz main

        baz_grd.plot(ax=ax2, cmap=sm_alpha, add_colorbar=False,
            vmin=0, vmax=region_baz[5] / plot_amp_factor, mouseover=True)
        baz_grd.plot.contour(ax=ax2, cmap='Greys_r', linewidths=.65, add_colorbar=False,
            levels=np.linspace(region_baz[5] / 8, region_baz[5] / plot_amp_factor, 4))

        ax2.scatter([x_max, x_max],[baz_5_vals.min(), baz_5_vals.max()],
            marker='_', s=100, c='white', zorder=10)

        ax2.axhline(y=0, color='darkred', linestyle='--')
        ax2.scatter(x_max, y_max, marker='d', c='darkred', s=55,
                    edgecolors='white', zorder=10)

        ax2.text(region_baz[0] + 10, 20,
            f'max ({int(baz_grd.max().item())}) at {y_max}$^\\circ$ Backazimuth',
            c='darkred', size=12,bbox={'facecolor': 'white', 'alpha': 0.85, 'pad': 1.5})

        ### curtailed slow baz ## ax4/5

        slow_grd_curtail.plot(ax=ax4, cmap=sm_alpha, add_colorbar=False,
            vmin=0, vmax=V_max_curtail / plot_amp_factor_curtail)
        slow_grd_curtail.plot.contour(ax=ax4, cmap='Greys_r', linewidths=.65, add_colorbar=False,
            levels=np.linspace(V_max_curtail / 8, V_max_curtail / plot_amp_factor_curtail, 4))
        ###
        baz_grd_curtail.plot(ax=ax5, cmap=sm_alpha, add_colorbar=False,
            vmin=0, vmax=V_max_curtail / plot_amp_factor_curtail)
        baz_grd_curtail.plot.contour(ax=ax5, cmap='Greys_r', linewidths=.65, add_colorbar=False,
            levels=np.linspace(V_max_curtail / 8, V_max_curtail / plot_amp_factor_curtail, 4))

        ###
        ax5.axhline(y=y_max, color='darkred', linestyle='-',lw=1.2)

        ## ax7/8 peaksss

        ax7.plot(midpoints, max_values, '-', lw=.25, c='black', alpha=.65)
        ax8.plot(midpoints_slow, max_values_slow, '-', lw=.25, c='black', alpha=.65)

        scatter_7 = ax7.scatter(midpoints, max_values, c=y_values,
            cmap=cmap_try, norm=norm,edgecolor='white', s=15, alpha=.88, linewidth=.15)
        scatter_8 = ax8.scatter(midpoints_slow, max_values_slow, c=y_values_slow,
            cmap=cmap_slow, norm=norm_slow,edgecolor='white', s=15, alpha=.88, linewidth=.15)

        # saving based on mean max
        if max_mean.item() > 10:
            ax7.scatter(midpoints[indexes], max_values[indexes],
                        marker='+', c='black', s=40, lw=1.25)
            ax8.scatter(midpoints[indexes], max_values[indexes],
                        marker='+', c='black', s=40, lw=1.25)

            grid_baz_offset.append(
                (grid_number, y_max, np.mean(y_values), np.std(y_values),
                 deets["ArrCen"][0], deets["ArrCen"][1], deets["ArrCen"][2],
                 deets["Event"][0], deets["Event"][1], deets["Event"][2],
                 deets["Dist"][0], deets["Baz"][0], deets["ArrCen"][3],
                 baz_5_vals.max() - baz_5_vals.min(),slow_5_vals.max() - slow_5_vals.min()))

        ##### ax9 histo

        ax9.hist(z_values_coh, density=True, bins='auto',
         histtype='stepfilled', alpha=0.65, color='orchid')
        ax9.axvline(avg_coherence, color='slateblue', linestyle='--', lw=1.25)

        ax9.set_xlim(0, baz_grd_curtail.max().item() / 5)
        ax9.set_yticks([])
        ax9.set_xticks([])

        fig.text(0.94, .75, f'mean={avg_coherence.item():.1f}',fontsize=11, ha='center', color='slateblue')
        fig.text(0.94, .725, f'std={std_coherence.item():.1f}',fontsize=11, ha='center')
        fig.text(0.94, .70, f'm/m={max_mean.item():.1f}',fontsize=11, ha='center')

        ### Colorbars, grids, labels, limits

        sm = plt.cm.ScalarMappable(norm=plt.Normalize(vmin=0, vmax=region_baz[5] / plot_amp_factor),
            cmap=sm_alpha)
        sm.set_array([])
        cbar = plt.colorbar(sm, cax=ax3, orientation='horizontal', extend='max')
        cbar.ax.xaxis.set_label_position('top')
        cbar.ax.xaxis.tick_top()
        cbar.ax.set_xlabel('Coherence', labelpad=5, fontsize=14)

        # Peak colorbars
        cbar7 = plt.colorbar(scatter_7)
        cbar7.set_label('Baz. ($^\\circ$)', fontsize=15)
        cbar7.set_ticks([-4, -2, 0, 2, 4])

        cbar8 = plt.colorbar(scatter_8)
        cbar8.set_label('Slow. (s/$^\\circ$)', fontsize=15)
        cbar8.set_ticks([1, 3, 5, 7, 9])

        ###

        for ax in [ax1, ax2, ax4, ax5, ax7, ax8]:
            ax.grid(which='major', linestyle='--', alpha=.75)
            ax.grid(which='minor', axis='x', linestyle='--', alpha=.65)

        set_locators(ax1, 'slow')
        set_locators(ax4, 'slow')
        set_locators(ax2, 'baz')
        set_locators(ax5, 'baz')

        ax1.set_ylim(2, 10)
        ax2.set_ylim(-25, 25)
        ax4.set_ylim(2, 10)
        ax5.set_ylim(-25, 25)

        ax1.set_ylabel('Slowness (s/$^\\circ$)')
        ax2.set_ylabel('Backazimuth ($^\\circ$)')
        # ax4.set_ylabel('Slowness (s/$^\\circ$)')
        # ax5.set_ylabel('Bazi ($^\\circ$)')
        ax8.set_ylabel('Coherence')
        ax7.set_xlabel('Time (s)')
        ax8.set_xlabel('Time (s)')
        for ax in [ax4, ax5]:
            ax.set_xticks([])
            ax.set_ylabel('')
        for ax in [ax1, ax2,ax4, ax5]:
            ax.set_xlabel('')

        ax7.set_yticks([])
        ax7.margins(x=0)
        ax8.margins(x=0)

        utc_dt=''.join(str(int(x)) for x in deets['Origin'])
        time_list=deets['Origin']
        formatted_time = f"Event origin: {int(time_list[0])} {int(time_list[1]):02d} {int(time_list[2]):02d} {int(time_list[3]):02d}:{int(time_list[4]):02d}"


        fig.text(0.2, .95, 'Grid #{}; {}'.format(grid_number,formatted_time),fontsize=16,color='Teal', ha='center', va='center')
        # fig_name='vespa_paper/picks_gridnum_{}_{}_{}_new.jpg'.format(grid_number,utc_dt,'AK')

        fig_name=py_figs+'picks_gridnum_{}_{}_{}.jpg'.format(grid_number,utc_dt,'II')
        # plt.show()
        plt.savefig(fig_name,dpi=400,bbox_inches='tight', pad_inches=0.1)
        plt.close('all')

print('----------DONE------------\n')
sys.exit()

#########


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

# the following adds white at the start of a color map.
# cmapp=plt.get_cmap(sm_alpha)
# colors = cmapp(np.linspace(0, 1, 256))
# colors[0:2] = [1, 1, 1, 1]  # White color
# new_lipari = plt.cm.colors.ListedColormap(colors)

### FOR beams
# slow_grd_bm.plot(ax=ax4,cmap=cm.broc,add_colorbar=False,mouseover=True)
# # slow_grd_bm.plot.contour(ax=ax4,cmap='Greys_r',linewidths=.65,add_colorbar=False,levels=np.linspace(region_slow[5]/8, region_slow[5]/plot_amp_factor, 4))
# baz_grd_bm.plot(ax=ax5,cmap=cm.broc,add_colorbar=False,mouseover=True)

# gets max 5% around max in slow/baz grids!!!
# grd,x_max,window_size,percent=slow_grd,x_max_slow,5,.05

# skew_coherence = stats.skew(z_values_coh)
# kurtosis_coherence = stats.kurtosis(z_values_coh)

## dividing peaks with slow > 6 and less than 6
# slowness=6 was chosen coz mostly PP comes at slow > 7 and based on trends
# 6 seemed like an appropriate value to distinguis precursory phases and PP.
# indexes_less6=[]
# indexes_gr6=[]
# y_values_slow contains peaks in slowness
# indexes contains index of selected peaks. Thus, y_values_slow[i] gets slowness vals of peaks in baz offsets.
# if len(indexes) ==0:
#     print('-----------------\n')
#     warnings.warn('Baz and slow mid-points array unequal :/; skipping this grid')
#     print('-----------------\n')
#     continue
# for i in indexes:
#     if y_values_slow[i] > 5.99:
#         indexes_gr6.append(i)
#     else:
#         indexes_less6.append(i)

# Create discrete colormap

### first if fro low slowness...thats' why P shenanigans
# if len(indexes_less6) ==0:
#     print('-----------------\n')
#     warnings.warn('nothing in slow < 6 :/; not plotting')
#     print('-----------------\n')
# else:
#     low_slow_combined=y_values[indexes_less6]
#     if max_mean.item() > 20:
#         grid_baz_offset_low_slow.append((grid_number,np.max(low_slow_combined),np.mean(low_slow_combined),np.std(low_slow_combined),deets["ArrCen"][0],deets["ArrCen"][1],deets["ArrCen"][2],deets["Event"][0],deets["Event"][1],deets["Event"][2],deets["Dist"][0],deets["Baz"][0],deets["ArrCen"][3],(baz_5_vals.max()-baz_5_vals.min()),(slow_5_vals.max()-slow_5_vals.min())))
#     ax7.scatter(midpoints[indexes_less6], max_values[indexes_less6], marker='+',color='darkred',s=40,linewidth=1.25)
#     ax8.scatter(midpoints_slow[indexes_less6], max_values_slow[indexes_less6], marker='+',color='darkred',s=40,linewidth=1.25)
#
# # second if for high slowness
# if len(indexes_gr6) ==0:
#     print('-----------------\n')
#     warnings.warn('nothing in slow > 6 :/; not plotting')
#     print('-----------------\n')
# else:
#     if max_mean.item() > 20:
#         grid_baz_offset_high_slow.append((grid_number,np.max(y_values[indexes_gr6]),np.mean(y_values[indexes_gr6]),np.std(y_values[indexes_gr6]),deets["ArrCen"][0],deets["ArrCen"][1],deets["ArrCen"][2],deets["Event"][0],deets["Event"][1],deets["Event"][2],deets["Dist"][0],deets["Baz"][0],deets["ArrCen"][3],(baz_5_vals.max()-baz_5_vals.min()),(slow_5_vals.max()-slow_5_vals.min())))
#     ax7.scatter(midpoints[indexes_gr6], max_values[indexes_gr6], marker='+',color='black',s=40,linewidth=1.25)
#     ax8.scatter(midpoints[indexes_gr6], max_values[indexes_gr6], marker='+',color='black',s=40,linewidth=1.25)
