# plots grd output from Vespapack script.
# can plot grd for mutiple earthquakes and arrays created for each earthquake.
# needs:
# - grds (beam slow/baz and XF slow/baz) in grid_folder
# - datapack txt files in folder named 'data_pack' (eg. Datapack_gridnum23_20221109_0938_.05_.5Hz_60samps_Zcomp_BPPCbase_num32_PP_Y_Y_3.0_Y_-1)
# two empty folder named 'pick_folder' and 'py_figs'. Former is for saving chosen values from vespagrams.
#
# tested on Py3.9.
# Shubh Agrawal
# 11th August 2025.
# USC, SC
###
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from matplotlib import ticker
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib.colors as mcolors
# import ScientificColourMaps8 as SCM8
import matplotlib.cm as cmm
from cmaptools import readcpt, joincmap, DynamicColormap
import glob as glob
import re
from obspy.core import UTCDateTime as utc
from obspy.taup.taup_geo import calc_dist,calc_dist_azi
from obspy.taup import TauPyModel
import subprocess
from matplotlib.widgets import CheckButtons
import scipy.signal as sci
import warnings
import scipy.stats as stats
from scipy.stats import norm, skewnorm, kurtosis
from matplotlib.colors import ListedColormap
##############################################################
# you might have to install these.
import xarray as xr
import pygmt
from cmcrameri import cm
# from IPython import get_ipython # uncomment if using Ipython.
####
# to be able to select points of interest on vespagram, you will need to download: https://mpl-point-clicker.readthedocs.io/en/latest/index.html
# from mpl_point_clicker import clicker
# from mpl_interactions import zoom_factory, panhandler
##
# sys.exit()
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

    arr_sS=model.get_travel_times(source_depth_in_km=eq_depth,distance_in_degree=dist,phase_list=["sS"])[0]
    arr_S=model.get_travel_times(source_depth_in_km=eq_depth,distance_in_degree=dist,phase_list=["S"])[0]


    return(arr_P,arr_PP,arr_pP,arr_sP,arr_sS,arr_S)

def extract_grid_list(grid_folder):
    file_pattern = os.path.join(grid_folder, 'xf_slow_grid_*.grd')
    files = glob.glob(file_pattern)

    # List to store the extracted grid numbers
    grid_numbers = []

    #  to extract numbers after 'grid_'
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
        ax.yaxis.set_minor_locator(MultipleLocator(2))
        ax.yaxis.set_major_locator(MultipleLocator(4))
    else:
        raise ValueError("Unsupported axis type: {}. Use 'default' or 'alt'.".format(axis_type))

def extract_grid_nums(main_folder):
    pattern = os.path.join(main_folder, '*.jpg')
    # List to store the gridnum values
    gridnum_list = []
    for file_path in glob.glob(pattern):
    # Extract the filename from the file path
        filename = os.path.basename(file_path)

        #  find the number after 'gridnum'
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

        i += step_size # Move the window by step size

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

    arr1.sort()
    arr2.sort()

    # compare elements
    for i in range(0, N):
        if (arr1[i] != arr2[i]):
            return False

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

##
## comment this if running from terminal
# get_ipython().magic('reset -sf')

# folder_pattern looks for earthquake folders.
# folder_pattern = "sac_files_with_P/*_inc2_r2.5"

# matching_folders = glob.glob(folder_pattern)

max_mean_gl=[]

# if you want to do on just one earthquake, change this to the specific directory.
matching_folders=['40_sac_noise_latN_R']

plt.rcParams.update({'font.size': 12})
plot_amp_factor=2
#
for folder in matching_folders:
    # location of the earthquake folder..change here.
    main_folder='/Users/keyser/Research/axisem_related_projs/plumes/plumes_iaspi91_5sec_new_loc_src_plume_wave/simu1D/output/stations/'+folder+'/'

    folder_datapack=main_folder+'data_pack/'
    grid_folder=main_folder+'grid_folder'
    pick_folder=main_folder+'py_picks/'
    py_figs=main_folder+'py_figs/'

    print('Main folder:',main_folder)
    # beacuse for each earthquake there are multiple arrays, we plot them on grid number basis.
    # so first we extract all grid numbers.
    gridnum_list=extract_grid_nums(main_folder)
    gridnum_list.sort()

    grid_baz_offset=[]
    grid_baz_offset_low_slow=[]
    grid_baz_offset_high_slow=[]

    # Chnage here if you want to just plot a specific grid number.
    for grid_number in gridnum_list:
    # for grid_number in [10]:

        # this values restrains the plot cpt cut-off. Esentially, we are plotting upto max_coherence/plot_amp_factor

        print('plot_amp_factor=',plot_amp_factor)
        print('grid_number=',grid_number)

        grid_list=extract_grid_list(grid_folder)

        #################
        # we read the datapack file containing details about the array.
        beam_deets=folder_datapack+extract_datapackfile(grid_number,folder_datapack)
        print(beam_deets,'\n')
        # one can alternatively directly test using a specific datapack txt file.
        #beam_deets '/Users/keyser/Research/sub_array_alaska/sac_files/230702_102743_PA/Datapack_20230702_1027_.05_.5Hz_60samps_Zcomp_WRHbase_gridnum2_num8_PP_Y_N_0.0_Y_-1.txt'
        patterns = {
            "Origin": re.compile(r"Origin: (\d+) (\d+) (\d+) (\d+):(\d+)"),
            "ArrCen": re.compile(r"ArrCen la/lo/elv: (-?\d+\.\d+) (-?\d+\.\d+) (\d+) Nst:(\d+)"),
            "ArrBaseStn": re.compile(r"ArrBaseStn: (\w+), grid la/lp (\d+), (-?\d+)"),
            "Event": re.compile(r"Event la/lo/dp: (-?\d+\.\d+) (-?\d+\.\d+) (\d+\.\d+)"),
            "Dist": re.compile(r"Dist: (\d+\.\d+)"),
            "Baz": re.compile(r"Baz \(Arr-Evt\): (\d+\.\d+)"),
            "Frequencies": re.compile(r"Frequencies: (\.\d+) - (\.\d+) Hz"),
            "TrcesSNR": re.compile(r"TrcesSNR mn,SD,min,max: (\d+\.\d+) (\d+\.\d+) (\d+\.\d+) (\d+\.\d+)"),
            "PredPP": re.compile(r"Pred PP \(prem\) time/U: (\d+\.\d+) (\d+\.\d+)")
        }

        # this dictionary conatins the extracted values values from datapack.txt
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
        # some shenanigans to read the values.
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

        # continued shenanigans for ArrBaseStn
        deets["ArrBaseStn"] = [float(x) for x in deets["ArrBaseStn"]]
        print('-------------\n')
        print(deets["ArrCen"])

        #############
        #------------------------
        # interp_slow=[0.1,0.05] # dont think I use this.

        # extract all four grd as xarrays.
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
        ## find global max and min in the whole grd.
        # i think there is a better way to do the following way, already used in this script, potentially.
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

        # get travel times for phases of interest to plot on Vespas.
        # change function calc_tt accordingly or just comment out places where these are plotted further.
        try:
            arr_P,arr_PP,arr_pP,arr_sP,arr_sS,arr_S=calc_tt(deets['Event'][0],deets['Event'][1],deets['ArrCen'][0],deets['ArrCen'][1],deets['Event'][2])
        except:
            print('one or more phases didnt arrive')

        #####
        # Plotting bit starts.
        # should move all analysis before this.

        fig = plt.figure(figsize=(15, 9))
        # values are: [left, bottom, width, height]
        ax1 = fig.add_axes([0.05, 0.05, 0.38, 0.52]) # slow XF
        ax2 = fig.add_axes([0.53, 0.05, 0.38, 0.52]) #  baz XF
        ax3=  fig.add_axes([0.94, .27, .01, 0.3]) # color bar
        ax4 = fig.add_axes([0.05, 0.65, 0.38, 0.32]) # slow beam
        ax5 = fig.add_axes([0.53, 0.65, 0.38, 0.32]) # baz beam
        ax7=  fig.add_axes([0.53, 0.05, 0.475, 0.095],sharex=ax2) # baz peaks
        ax8=  fig.add_axes([0.05, 0.05, 0.475, 0.095],sharex=ax1) # slow peaks
        ax9=fig.add_axes([0.91, 0.65, 0.07, 0.32]) # histogram
        ax1.set_ylim(5,20)
        # ###check box stuff.. to enable click based selection.
        #
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
        ####check box stuff
        ####

        #### if you haven't installed cmcrameri, replace cm.lipari_r with a matplotlib colour scheme
        cmap_lip=cm.lipari_r
        # cmap_lip=cmm.BuPu
        colA = cmap_lip(np.arange(cmap_lip.N))
        # #replaces the last column (alpha channel) in colA with values from 0.7 to 1,
        # # creating a gradient in transparency across the color map. For more transparency, reduce the decimal.
        colA[:,-1] = np.linspace(0.7, 1, cmap_lip.N)

        sm_alpha = ListedColormap(colA)

        ### the following adds white at the start of a color map.
        # cmapp=plt.get_cmap(sm_alpha)
        # colors = cmapp(np.linspace(0, 1, 256))
        # colors[0:5] = [1, 1, 1, 1]  # White color
        # new_cmap = plt.cm.colors.ListedColormap(colors)

        # plotting grids
        slow_grd.plot(ax=ax1,cmap=sm_alpha,add_colorbar=False,vmin=0,vmax=region_baz[5]/plot_amp_factor,mouseover=True)
        slow_grd.plot.contour(ax=ax1,cmap='Greys_r',linewidths=.65,add_colorbar=False,levels=np.linspace(region_baz[5]/8, region_baz[5]/plot_amp_factor, 4))

        baz_grd.plot(ax=ax2,cmap=sm_alpha,add_colorbar=False,vmin=0,vmax=region_baz[5]/plot_amp_factor,mouseover=True)
        baz_grd.plot.contour(ax=ax2,cmap='Greys_r',linewidths=.65,add_colorbar=False,levels=np.linspace(region_baz[5]/8, region_baz[5]/plot_amp_factor, 4))

        slow_grd_bm.plot(ax=ax4,cmap=cmm.RdGy,add_colorbar=False,mouseover=True)
        baz_grd_bm.plot(ax=ax5,cmap=cmm.RdGy,add_colorbar=False,mouseover=True)

        # gets max 5% around max in slow/baz grids! A quality measure, if you will.
        # grd,x_max,window_size,percent=slow_grd,x_max_slow,5,.05
        #
        # grd_5_slow,slow_5_vals=get_contour_around_max(slow_grd,x_max_slow,5,.05)
        # grd_5_baz,baz_5_vals=get_contour_around_max(baz_grd,x_max,5,.05)

        # plot the 5% values as white bars.
        # ax1.scatter([x_max_slow,x_max_slow],[slow_5_vals.min(),slow_5_vals.max()],marker='_',s=100,c='white',zorder=10)
        # ax2.scatter([x_max,x_max],[baz_5_vals.min(),baz_5_vals.max()],marker='_',s=100,c='white',zorder=10)

        # ax1.scatter(x_max_slow,slow_5_vals.max(),marker='_',s=80,c='white',zorder=10)

        # TEXT FOR MIN/MAX baz slowness!
        # ax1.text(region_baz[0]+5, 10.1, 'max/min/diff slow (5%) around max coh. ($\pm$5 sec): {:.1f}/{:.1f}/{:.1f}'.format(slow_5_vals.max(),slow_5_vals.min(),(slow_5_vals.max()-slow_5_vals.min())),c='navy',size=11,weight='roman')
        # ax2.text(region_baz[0]+5, 51, 'max/min/diff baz (5%) around max coh. ($\pm$5 sec): {:.1f}/{:.1f}/{:.1f}'.format(baz_5_vals.max(),baz_5_vals.min(),(baz_5_vals.max()-baz_5_vals.min())),c='navy',size=11,weight='roman')

        ######
        ax1.grid(which='minor',axis='x',color='dimgrey', linestyle='--',linewidth=.65,alpha=.75)
        ax1.grid(which='major',axis='both',color='dimgrey', linestyle='--',linewidth=.75,alpha=.75)
        ax2.grid(which='minor',axis='x',color='dimgrey', linestyle='--',linewidth=.65,alpha=.75)
        ax2.grid(which='major',axis='both',color='dimGrey', linestyle='--',linewidth=.75,alpha=.75)

        # plot great circle path on XF baz plot as red dotted line
        ax2.axhline(y=0, color='darkred', linestyle='--')
        ax2.scatter(x_max,y_max,marker='d',c='darkred',s=55,edgecolors='white',zorder=10)
        ax2.text(region_baz[0]+10, 18, 'max ({}) at {}$^\circ$ bazi'.format(int(baz_grd.max().item()),y_max),c='darkred',size=13,weight='roman',bbox={'facecolor': 'white', 'alpha': 0.85, 'pad': 1.5})
        ##
        # change the function set_locators according to your time, slow, baz domain.

        set_locators(ax1, 'slow')
        set_locators(ax4, 'slow')
        set_locators(ax2, 'baz')
        set_locators(ax5, 'baz')

        # change to phases of interest or comment it altogether.
        for phase in [arr_sS,arr_S]:
            # if
        # for phase in [arr_P,arr_pP,arr_sP,arr_PP]:

            if 'diff' in phase.name:
                ax1.scatter(phase.time,phase.ray_param*0.0174533,marker='o',c='CORNFLOWERBLUE',s=50,edgecolors='white',zorder=10)
                ax1.text(phase.time, 1.5+phase.ray_param * 0.0174533, phase.name, bbox={'facecolor': 'white', 'alpha': 0.85, 'pad': 1.5},fontsize=9.5,c='CORNFLOWERBLUE', rotation='vertical',ha='center')
                ax2.scatter(phase.time,0,marker='o',c='CORNFLOWERBLUE',s=50,edgecolors='white',zorder=10)
                ax2.text(phase.time, 2.5, phase.name, bbox={'facecolor': 'white', 'alpha': 0.85, 'pad': 1.5},fontsize=9.5,c='CORNFLOWERBLUE', rotation='vertical',ha='center')

            else:
                ax1.scatter(phase.time,phase.ray_param*0.0174533,marker='o',c='violet',s=50,edgecolors='white',zorder=10)
                ax1.text(phase.time, 1.5+phase.ray_param * 0.0174533, phase.name, bbox={'facecolor': 'white', 'alpha': 0.85, 'pad': 1.5},fontsize=9.5,c='violet', rotation='vertical',ha='center')
                ax2.scatter(phase.time,0,marker='o',c='violet',s=50,edgecolors='white',zorder=10)
                ax2.text(phase.time, 2.5, phase.name, bbox={'facecolor': 'white', 'alpha': 0.85, 'pad': 1.5},fontsize=9.5,c='violet', rotation='vertical',ha='center')


        plt.rcParams['axes.labelsize'] = 14
        ax1.set_ylabel('Slowness (s/$^\circ$)')
        ax1.set_xlabel('Time (s)')
        # ax2.set_xticklabels([])
        ax2.set_ylabel('Bazi ($^\circ$)')
        ax2.set_xlabel('Time (s)')
        # plt.rcParams['axes.labelsize'] = 15
        ax4.set_ylabel('Slowness (s/$^\circ$)')
        ax4.set_xlabel('')
        ax5.set_ylabel('Bazi ($^\circ$)')
        ax5.set_xlabel('')
        #####
        # fig.text(0.42, -0.05, 'Time (s)',fontsize=12, ha='center', va='center')

        # plotting the colour scale
        sm = plt.cm.ScalarMappable(norm=plt.Normalize(vmin=0, vmax=region_baz[5]/plot_amp_factor),cmap=sm_alpha )

        sm.set_array(np.arange(0,1.33))
        cbar = plt.colorbar(sm,cax=ax3,orientation='vertical', extend='max')
        # various ways to print 'coherence' around the scale. Play and choose.
        cbar.set_label('Coherence', rotation=90, labelpad=2,loc='center')
        # cbar.ax.set_xlabel('Coherence',labelpad=2)
        # cbar.ax.xaxis.set_label_position('top')
        # cbar.ax.xaxis.tick_top()
        # cbar.ax.set_xlabel('Coherence',labelpad=5,fontsize=13)

        #############
        ### one can make use of the following lines to choose specific bits of grid based on time.
        ###
        # slow_grd_curtail=slow_grd.where((slow_grd.x > arr_pP.time - 10) & (slow_grd.x < arr_pPP.time + 20), drop=True)
        # baz_grd_curtail=baz_grd.where((baz_grd.x > arr_pP.time - 10) & (baz_grd.x < arr_pPP.time + 20), drop=True)
        #
        # slow_grd_justP=slow_grd.where((slow_grd.x < arr_pP.time - 10) , drop=True)
        # baz_grd_justP=baz_grd.where((baz_grd.x < arr_pP.time - 10), drop=True)

        ### function get_peaks_grd find peaks for a given grd file.
        ### uses sci.find_peaks
        ### parameters that can be changed: time_step=5,overlap=2.5, prominence etc.

        midpoints,max_values,y_values,indexes=get_peaks_grd(baz_grd)
        midpoints_slow,max_values_slow,y_values_slow,indexes_slow=get_peaks_grd(slow_grd)

        ##### finding if the energy is everywhere or not. QUALITY MEASURE.
        #####
        avg_coherence = baz_grd.mean(dim=['x', 'y'])
        std_coherence= baz_grd.std(dim=['x', 'y'])
        z_values_coh = baz_grd.values.flatten()
        # skew_coherence = stats.skew(z_values_coh)
        # kurtosis_coherence = stats.kurtosis(z_values_coh)
        max_mean= baz_grd.max(dim=['x', 'y'])/baz_grd.mean(dim=['x', 'y'])

        max_mean_gl.append(round(max_mean.item(),2))

        ## if for some reason, one has to divide the peaks based on slowness or baz value, the following routine can be used.
        ## otherwise, ignore.
        ## dividing peaks with slow > 6 and less than 6

        # indexes_less6=[]
        # indexes_gr6=[]
        # # y_values_slow contains peaks in slowness
        # # indexes contains index of selected peaks. Thus, y_values_slow[i] gets slowness vals of peaks in baz offsets.
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
        num_bins = 8 # for discrete colour bar
        cmap_try= cmm.PiYG
        cmap_slow= cmm.PiYG
        norm = mcolors.BoundaryNorm(np.linspace(-4,4, num_bins + 1), cmap_try.N)
        norm_slow = mcolors.BoundaryNorm(np.linspace(1,20, num_bins + 1), cmap_slow.N)

        ## ax9 plots the pdf.
        ax9.hist(z_values_coh, density=True, bins='auto', histtype='stepfilled', alpha=0.65,color='orchid')
        ax9.axvline(avg_coherence, color='slateblue', linestyle='--', lw=1.25)#label=f'{avg_coherence.item():.1f}',
        ax9.set_yticklabels([])
        ax9.set_ylim((0, None))
        ax9.set_xlim(0,baz_grd.max().item()/5 )
        # plt.ylabel('coherence')
        # plt.title('Histogram')
        ax9.yaxis.set_label_position("right")
        ax9.yaxis.tick_right()

        # Set the label text
        ax9.set_ylabel("Probability Density",fontsize=11)

        # fig.text(0.94, .75, f'mean={avg_coherence.item():.1f}',fontsize=11, ha='center', va='center',color='slateblue')
        # fig.text(0.94, .725, f'std={std_coherence.item():.1f}',fontsize=11, ha='center', va='center')
        # fig.text(0.94, .70, f'm/m={max_mean.item():.1f}',fontsize=11, ha='center', va='center')

        ax9.set_yticks([])
        ax9.set_xticks([])
        plt.grid(True)

        ### Start of plotting peaks ax7 and ax8
        # ax7 is baz
        # ax7.set_facecolor("whitesmoke")
        ax7.grid(True,alpha=.25)
        ax7.plot(midpoints, max_values, '-',lw=.25,c='black',alpha=.65)

        # ax8.set_facecolor("whitesmoke")
        ax8.grid(True,alpha=.25)
        ax8.plot(midpoints_slow, max_values_slow, '-',lw=.25,c='black',alpha=.65)

        #plot all mid points and color by baz
        scatter_7 = ax7.scatter(midpoints, max_values, c=y_values, cmap=cmap_try, norm=norm, edgecolor='white', s=15,alpha=.88,linewidth=.15)
        scatter_8 = ax8.scatter(midpoints_slow, max_values_slow, c=y_values_slow, cmap=cmap_slow, norm=norm_slow, edgecolor='white', s=15,alpha=.88,linewidth=.15)


        cbar = plt.colorbar(scatter_7)
        cbar.set_label('Bazi$^\circ$', fontsize=12)
        cbar.set_ticks([-4,-2,0,2,4])

        cbar_slow = plt.colorbar(scatter_8)
        cbar_slow.set_label('Slow. (s/$^\circ$)', fontsize=12)
        cbar_slow.set_ticks([5,10,15,20])
        # cbar_slow.set_ticks([5,10,15])

        # ax7.set_xticklabels([])
        ax7.set_yticklabels([])
        # ax8.set_ylabel('Coherence')
        ax8.set_xlabel('Time (s)')
        ax7.set_xlabel('Time (s)')

        ax8.set_yticklabels([])
        ax7.set_yticks([])
        # ax8.set_yticks([])

        ##
        # ax7.set_xlabel('Time (s)'), ax7.set_ylabel('Coherence')
        utc_dt=''.join(str(int(x)) for x in deets['Origin'])
        time_list=deets['Origin']
        formatted_time = f"Event origin: {int(time_list[0])} {int(time_list[1]):02d} {int(time_list[2]):02d} {int(time_list[3]):02d}:{int(time_list[4]):02d}"

        if max_mean.item() > 20:
            # grid_baz_offset.append((grid_number,y_max,np.mean(y_values),np.std(y_values),deets["ArrCen"][0],deets["ArrCen"][1],deets["ArrCen"][2],deets["Event"][0],deets["Event"][1],deets["Event"][2],deets["Dist"][0],deets["Baz"][0],deets["ArrCen"][3],(baz_5_vals.max()-baz_5_vals.min()),(slow_5_vals.max()-slow_5_vals.min())))
            # save the grid_baz_offset as text, if you so desire.
            # np.savetxt('maxVals.txt'.format(), grid_baz_offset, fmt='%i %.1f %.2f %.2f %.2f %.2f %.1f %.2f %.2f %.1f %.1f %.1f %i %.1f %.1f')

            fig.text(0.2, .98, 'Grid #{}; {}'.format(grid_number,formatted_time),fontsize=16,color='Teal', ha='center', va='center')
            fig_name='vespa_paper/picks_gridnum_{}_{}_{}.jpg'.format(grid_number,utc_dt,'AK')

            fig_name=py_figs+'picks_gridnum_{}_{}_{}.jpg'.format(grid_number,utc_dt,'II')

            plt.savefig(fig_name,dpi=400,bbox_inches='tight', pad_inches=0.1)
            # plt.savefig('example_plot_vespapack.jpg',dpi=400,bbox_inches='tight', pad_inches=0.1)




        plt.close('all')
        # sys.exit()
        # plt.show()
        # sys.exit()
    print('--------------------------------\n')
