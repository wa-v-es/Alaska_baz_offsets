import numpy as np
import matplotlib.pyplot as plt
import sys
from cmcrameri import cm
from scipy.optimize import minimize
import math
import matplotlib
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from obspy.taup import TauPyModel
######
def normalize(vector):
    return vector / np.linalg.norm(vector)
##
def plot_depth_arr(stations_array,moho):
    # time=float(t)
    plt.ion()

    fig,ax = plt.subplots(figsize=(6, 5))
    # ax.set_facecolor(("xkcd:sandy yellow",.08))

    plt.scatter(stations_array[:, 0], stations_array[:, 1], c=moho,marker='^', cmap='cmc.devon_r', s=100,alpha=.98,lw=.5,edgecolors='k')
    # plt.scatter(stations_array[:, 0], stations_array[:, 1], c='grey',marker='^', s=90,alpha=.88,lw=.5,edgecolors='k')
    for i,st in enumerate(stations_array):
        x,y=st
        # plt.text(x+5, y,i,fontsize=9,c='xkcd:dusk blue')

    cbar = plt.colorbar(label='Moho (km)', fraction=0.1, shrink=0.25)
    cbar.ax.set_position([0.65, 0.57, 0.2, 0.29]) # [left, bottom, width, height]
    cbar.set_ticks([25, 35, 45])
    plt.xlabel('X (km)')
    plt.ylabel('Y (km)')
    # plt.title('Max/min Moho ={:.1f}/{:.1f} km'.format(np.max(moho),np.min(moho)))
    plt.show()
###
def plot_time_arr(stations_array,time_arr,az,save=False):
    time=np.max(time_arr)
    fig,ax = plt.subplots(figsize=(6, 5))
    ax.set_facecolor(("xkcd:pale",.08))

    plt.scatter(stations_array[:, 0], stations_array[:, 1], c=time_arr,marker='^', cmap='PiYG', s=90,alpha=.88,lw=.5,edgecolors='k')
    # plt.scatter(stations_array[:, 0], stations_array[:, 1], c='grey',marker='^', s=90,alpha=.88,lw=.5,edgecolors='k')
    for i,st in enumerate(stations_array):
        x,y=st
        plt.text(x+5, y,i,fontsize=9,c='xkcd:dusk blue')
    plt.colorbar(label='Rel. arrival (s)',fraction=.12,shrink=.5)
    plt.xlabel('X (km)')
    plt.ylabel('Y (km)')
    plt.title('Max time ={:.2f} s. phi={}'.format(time,az))
    if save:
        plt.savefig('crust_mantle_time_diff/cr_m_time_diff_{}'.format(az),dpi=300,bbox_inches='tight', pad_inches=0.1)
        plt.close()
    else:
        plt.show()
    plt.show()
###
def calculate_incidence_vector(i_deg, phi_deg):
    i_rad = np.radians(i_deg)
    phi_rad = np.radians(phi_deg)
    # components of the incidence vector
    x_component = np.sin(i_rad) * np.sin(phi_rad) # x is sin component as phi is measured from north..
    y_component = np.sin(i_rad) * np.cos(phi_rad)
    z_component = np.cos(i_rad)
    #
    incidence_vector = np.array([x_component, y_component, z_component])
    return incidence_vector
###
###
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
def snells_law(incident_ray, normal, n1, n2):
    # Normalize the vectors
    incident_ray = normalize(incident_ray)
    normal = normalize(normal)

    # angle of incidence
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

    return refracted_ray, math.degrees(theta_i), math.degrees(theta_r), azimuth, theta_i_surf
###
##
stations_per_arm = 3
arm_length = 279
distance_between_stations = arm_length / stations_per_arm

# Create stations array
stations = [(0, 0)]
for i in range(1, stations_per_arm + 1):
    stations.append((i * distance_between_stations, 0))
    stations.append((-i * distance_between_stations, 0))
    stations.append((0, i * distance_between_stations))
    stations.append((0, -i * distance_between_stations))

stations_array = np.array(stations)
###
# Moho and ray vectors
v_crust = 5.8  # km/s
v_mantle = 8.0
moho_normal = np.array([0.05, 0.0, 1])
# moho_normal = np.array([0.0, 0.1, 1])
n1 = 1.0                             # Ri of mantle
n2 = 1.38
# phi=45
i_deg = 28.164 # this is angle of incidence at MOHO..given i=20.0004 at surface
for phi in range(0,95,5):
    incident_ray = calculate_incidence_vector(i_deg,phi) # i and azimuth.
    moho_base_depth = 35  # km
    refracted_ray, _, _, _,_ = snells_law(incident_ray, moho_normal, n1, n2)
    ##
    # Compute travel times
    travel_times = []
    z_moho_arr=[]
    for x, y in stations_array:
        # Moho depth at station using plane equation
        z_moho = moho_base_depth + moho_normal[0] * x + moho_normal[1] * y
        z_moho_arr.append(z_moho)
        # Distance in crust
        t_crust = z_moho / refracted_ray[2]  # Time to reach surf
        crust_path_length = t_crust * np.linalg.norm(refracted_ray)  # Actual path length in crust

        # Distance in mantle
        mantle_path_length = (100 - z_moho) / incident_ray[2] * np.linalg.norm(incident_ray)  # Mantle travel

        # Compute total travel time
        total_time = crust_path_length / v_crust + mantle_path_length / v_mantle
        travel_times.append(total_time)

    # Convert to numpy array and compute relative arrivals
    travel_times = np.array(travel_times)
    relative_arrivals = travel_times - travel_times[0]  # Subtract central station travel time

    # Print results
    for (x, y), moho, rel_t in zip(stations_array, z_moho_arr, relative_arrivals):
        print(f"Station ({x:.1f}, {y:.1f}) -> Moho: {moho:.1f}km -> Relative Arrival: {rel_t:.2f} s")
    #

    # plot_time_arr(stations_array,relative_arrivals,phi,save=False)
plt.rcParams.update({'font.size': 15})
plot_depth_arr(stations_array,z_moho_arr)
plt.savefig('depth_moho_fig2.png',dpi=300,bbox_inches='tight', pad_inches=0.1)
