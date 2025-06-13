# uses saved grids to plot the vespas..
# doesn't actually shift the vespas..it happens in plot_xyx.py..maybe should change the names..

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
from matplotlib.colors import ListedColormap
import math
from mpl_toolkits.mplot3d import Axes3D

#####
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
#
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
####
def normalize(vector):
    return vector / np.linalg.norm(vector)
#
def get_incidence_angle(p,r,v):
    p=p/np.deg2rad(1)
    i=np.arcsin(p*v/r)
    return np.round(np.rad2deg(i),3)
##
def extract_angles(vector):

    vector_normalized = vector / np.linalg.norm(vector)
    x, y, z = vector_normalized
    #  angle of incidence (i)
    i_rad = np.arccos(z)
    i_deg = np.degrees(i_rad)
    # azimuth (phi)
    phi_rad = np.arctan2(x, y) # phi is measured from north
    phi_deg = np.degrees(phi_rad)

    return i_deg, phi_deg
###
def snells_law(p,normal):
    p=p/np.deg2rad(1)
    v=5.8
    r=6371
    i=np.arcsin(p*v/r) # angle of incidene at surface

    # ray going from surafce to Moho. R1=crust.
    x_ray=np.tan(i)
    n1 = 1.38   # R1 is now
    n2 = 1.0
    # convert incidence angle to vector assuming ray coming from az=0=y.
    incident_ray = np.array([x_ray, 0, -1])  # 0.36 is 20 deg incidence angle for parallel plane
    # normal = np.array([0.05, 0, -1])
    # normal = np.array([0.0, 0, -1])

    # Normalize the vectors
    incident_ray = normalize(incident_ray)
    normal = normalize(normal)

    # angle of incidence with the Moho
    cos_theta_i = np.dot(incident_ray, normal)
    theta_i = np.arccos(cos_theta_i)
    sin_theta_i = np.sqrt(1 - cos_theta_i**2)

    # Snell's Law
    sin_theta_r = (n1 / n2) * sin_theta_i
    if sin_theta_r > 1:
        # Total internal reflection
        return None, theta_i, None, None

    cos_theta_r = np.sqrt(1 - sin_theta_r**2)
    theta_r = np.arccos(cos_theta_r)

    #refracted ray
    refracted_ray = (n1 / n2) * incident_ray + (cos_theta_r - (n1 / n2) * cos_theta_i) * normal
    refracted_ray = normalize(refracted_ray)

    theta_i_surf,azimuth=extract_angles(refracted_ray)

    # i=np.deg2rad(theta)
    r_m=6336
    v_m=8
    p_m=r_m*np.sin(theta_i_surf)/v_m
    # return np.round((p_m*np.deg2rad(1)),2)

    return np.round((p_m*np.deg2rad(1)),3),math.degrees(i),incident_ray,refracted_ray, math.degrees(theta_i),  azimuth, theta_i_surf
###
def convert_rayP_surf_to_moho(p):
    p=p/np.deg2rad(1)
    v=5.8
    r=6371
    i=np.arcsin(p*v/r)
    ##
    n1 = 1.0     # Ri of mantle
    n2 = 1.3793     # n2 = v1/v2
    #
    i_m= np.arcsin(n2*np.sin(i)) # snells law
    ##
    # i=np.deg2rad(theta)
    p=6336*np.sin(i_m)/8
    return np.round((p*np.deg2rad(1)),3)
##
###
def plot_plane(normal, point, ax, size=10):
    # Create a grid of points
    d = -point.dot(normal)
    xx, yy = np.meshgrid(range(-size, size), range(-size, size))

    # Calculate corresponding z values for the plane
    zz = (-normal[0] * xx - normal[1] * yy - d) * 1. / normal[2]

    ax.plot_surface(xx, yy, zz, color='y', alpha=0.25)
###

#######
normal = np.array([0.05, 0, -1])

rayP_m,i_surf,incident_ray,refracted_ray,i_m,r_m,az=snells_law(4.73,normal)
print('rayP_m,i_surf,incident_ray,refracted_ray,i_m,r_m,az\n')

print(rayP_m,i_surf,incident_ray,refracted_ray,i_m,r_m,az)
# sys.exit()

plot_vector=False

if plot_vector:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.quiver(-incident_ray[0], -incident_ray[1], -incident_ray[2], incident_ray[0], incident_ray[1], incident_ray[2], color='cadetblue', label='Incident Ray')
    # ax.quiver(0, 0, 0, incident_ray[0], incident_ray[1], incident_ray[2], color='cadetblue', label='Incident Ray')
    ax.quiver(0, 0, 0, normal[0], normal[1], normal[2], color='purple', label='Normal')
    ax.quiver(0, 0, 0, refracted_ray[0], refracted_ray[1], refracted_ray[2], color='indianred', label='Refracted Ray')
    plot_plane(normal, np.array([0, 0, 0]), ax)

    # Set limits
    ax.set_xlim([-1, 1])
    ax.set_ylim([-1, 1])
    ax.set_zlim([-1, 1])
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.legend()
    plt.show()
# sys.exit()




folder_pattern = "sac_files_with_P/*_inc2_r2.5"
folder_pattern = "sac_files/*_inc2_r2.5"

matching_folders = glob.glob(folder_pattern)
matching_folders=['sac_files_with_P/220526_120223_SA_inc2_r2.5']

for folder in matching_folders:
    main_folder='/Users/keyser/Research/AK_all_stations/'+folder+'/'
    # main_folder='/Users/keyser/Research/axisem/moho_3d/moho_dip_prllN_10s_dir_no_smooth/simu3D/output/stations/AK_81/'+folder+'/'
    # main_folder='/Users/keyser/Research/axisem/loyalty_isl/output_10sec_expl_source/stations/AK_81/'+folder+'/'

    folder_datapack=main_folder+'data_pack/'
    grid_folder=main_folder+'grid_folder'
    pick_folder=main_folder+'py_picks/'
    py_figs=main_folder+'py_figs/'
    xyz_folder=main_folder+'xyz_folder/'


    print('Main folder:',main_folder)
    gridnum_list=extract_grid_nums(main_folder)
    gridnum_list.sort()

    plot_amp_factor=2
    # plot_amp_factor=10
    print('plot_amp_factor=',plot_amp_factor)
    # sys.exit()
    # for grid_number in gridnum_list:
    for grid_number in [93]:
        xf_slow_xyz= np.loadtxt(xyz_folder+'xf_slow_xyz_93_220526_.05_.5.xyz')
        xf_baz_xyz= np.loadtxt(xyz_folder+'xf_baz_xyz_93_220526_.05_.5.xyz')

        sys.exit()

        slow_grd,region_slow=extract_region_from_grid(grid_number,grid_folder,'slow','xf')
        print(region_slow,'\n')
        #####
        baz_grd,region_baz=extract_region_from_grid(grid_number,grid_folder,'baz','xf')
        print(region_baz)

        fig = plt.figure(figsize=(12, 5))
        ax1 = fig.add_axes([0.05, 0.05, 0.38, 0.82]) #[left, bottom, width, height]
        ax2 = fig.add_axes([0.53, 0.05, 0.38, 0.82]) #[left, bottom, width, height]

        cmap_lip=cm.oslo_r
        colA = cmap_lip(np.arange(cmap_lip.N))
        colA[:,-1] = np.linspace(0.65, 1, cmap_lip.N) #replaces the last column (alpha channel) in colA with values from 0.75 to 1, creating a gradient in transparency across the color map.

        sm_alpha = ListedColormap(colA)

        slow_grd.plot(ax=ax1,cmap=sm_alpha,add_colorbar=False,vmin=0,vmax=region_baz[5]/plot_amp_factor,mouseover=True)
        slow_grd.plot.contour(ax=ax1,cmap='Greys_r',linewidths=.65,add_colorbar=False,levels=np.linspace(region_baz[5]/8, region_baz[5]/plot_amp_factor, 4))

        baz_grd.plot(ax=ax2,cmap=sm_alpha,add_colorbar=False,vmin=0,vmax=region_baz[5]/plot_amp_factor,mouseover=True)
        baz_grd.plot.contour(ax=ax2,cmap='Greys_r',linewidths=.65,add_colorbar=False,levels=np.linspace(region_baz[5]/8, region_baz[5]/plot_amp_factor, 4))

        max_position = baz_grd.argmax(dim=['y', 'x'])
        max_position_slow = slow_grd.argmax(dim=['y', 'x'])
        y_max = baz_grd['y'][max_position['y']].item()
        x_max = baz_grd['x'][max_position['x']].item()

        # plt.show()
        ax1.grid(which='minor',axis='x',color='dimgrey', linestyle='--',linewidth=.65,alpha=.75)
        ax1.grid(which='major',axis='both',color='dimgrey', linestyle='--',linewidth=.75,alpha=.75)

        ax2.grid(which='minor',axis='x',color='dimgrey', linestyle='--',linewidth=.65,alpha=.75)
        ax2.grid(which='major',axis='both',color='dimGrey', linestyle='--',linewidth=.75,alpha=.75)

        ### plot horiz line at baz max_baz
        ax2.axhline(y=0, color='darkred', linestyle='--')
        ax2.scatter(x_max,y_max,marker='d',c='darkred',s=55,edgecolors='white',zorder=10)
        ax2.text(region_baz[0]+10, 40, 'max ({}) @ {}$^\circ$'.format(int(baz_grd.max().item()),y_max),c='darkred',size=12,weight='roman')
        ##
        set_locators(ax1, 'slow')
        # set_locators(ax4, 'slow')
        set_locators(ax2, 'baz')
        # set_locators(ax5, 'baz')
        ax1.set_ylabel('Slowness (s/$^\circ$)')
        ax1.set_xlabel('Time (s)')
        # ax2.set_xticklabels([])
        ax2.set_ylabel('Bazi ($^\circ$)')
        ax2.set_xlabel('Time (s)')
        # plt.savefig('vespa_{}.png'.format(grid_number),dpi=300,bbox_inches='tight', pad_inches=0.1)
        plt.show()
