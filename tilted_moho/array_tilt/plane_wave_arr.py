# for a plane and titlted moho, this script calculated the travel time difference for an array for changing azi.
# it assumes the angle of incidence is same at all stations.
# which in effect reports same values as the refracted ray for one station by a titlted moho.
# the travel time differene is optimized using equation in  Schweitzer et al. (2012)
import numpy as np
import matplotlib.pyplot as plt
import sys
from cmcrameri import cm
from scipy.optimize import minimize
import math
import matplotlib
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
#########
def normalize(vector):
    return vector / np.linalg.norm(vector)
def plot_time_arr(stations_array,time_arr,t):
    time=float(t)
    plt.scatter(stations_array[:, 0], stations_array[:, 1], c=time_arr,marker='^', cmap='cmc.cork', s=90,alpha=.88)
    plt.colorbar(label='Time of Arrival (s)',fraction=.12,shrink=.5)
    plt.xlabel('X (km)')
    plt.ylabel('Y (km)')
    plt.title('Max time ={:.2f}'.format(time))
    plt.show()
###
def objective(params):
    phi, i = params
    phi_rad = np.radians(phi)
    i_rad = np.radians(i)

    total_error = 0

    for idx, station in enumerate(stations_array):
        x, y = station
        predicted_time = np.sin(i_rad)*(-x * np.sin(phi_rad) - y * np.cos(phi_rad)) / (5.8)
        total_error += (predicted_time - time_arr[idx]) ** 2

    return total_error
###
def get_ray_param(theta,r,v):
    i=np.deg2rad(theta)
    p=r*np.sin(i)/v
    return np.round((p*np.deg2rad(1)),2)
###
def calculate_incidence_vector(i_deg, phi_deg):
    i_rad = np.radians(i_deg)
    phi_rad = np.radians(phi_deg)
    # components of the incidence vector
    x_component = np.sin(i_rad) * np.sin(phi_rad)
    y_component = np.sin(i_rad) * np.cos(phi_rad)
    z_component = np.cos(i_rad)
    #
    incidence_vector = np.array([x_component, y_component, z_component])
    return incidence_vector
###
def extract_angles(vector):

    vector_normalized = vector / np.linalg.norm(vector)
    x, y, z = vector_normalized
    #  angle of incidence (i)
    i_rad = np.arccos(z)
    i_deg = np.degrees(i_rad)
    # azimuth (phi)
    phi_rad = np.arctan2(x, y)
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

    # azimuth angle (angle from y-axis)

    # azimuth = np.arctan2(refracted_ray[1], refracted_ray[0])
    theta_i_surf,azimuth=extract_angles(refracted_ray)
    #angle of incidene at syrface for refracted ray
    # surface  = np.array([0, 0, 1])
    # cos_theta_i_surf = np.dot(refracted_ray, surface)
    # theta_i_surf = np.arccos(cos_theta_i_surf)

    return refracted_ray, math.degrees(theta_i), math.degrees(theta_r), azimuth, theta_i_surf
#
def get_inci_surf(inci_moho):
    n1 = 1.0                      # Ri of mantle
    n2 = 1.38
    sin_theta_i=np.sin(np.deg2rad(inci_moho))
    sin_theta_r = (n1 / n2) * sin_theta_i
    theta_r = np.arcsin(sin_theta_r)
    return np.rad2deg(theta_r)
###
###

stations_per_arm = 3
arm_length = 279
distance_between_stations = arm_length / (stations_per_arm)

stations = []
stations.append((0, 0))
for i in range(1, stations_per_arm + 1):
    stations.append((i * distance_between_stations, 0))
    stations.append((-i * distance_between_stations, 0))
    stations.append((0, i * distance_between_stations))
    stations.append((0, -i * distance_between_stations))

stations_array = np.array(stations)

# print("Stations coordinates (km):")
# print(stations_array)
### calculate travel time difference wrt to centre station.
fig = plt.figure(figsize=(10, 5))
ax = fig.add_subplot(111)
ax1=ax.twinx()
# ax1 = fig.add_subplot(212)
matplotlib.rcParams.update({'font.size': 12})
ax.set_facecolor(("xkcd:light grey blue",.15))

for phi_deg in range(0,95,5):

    # phi_deg = 55
    symbol=['>','*','o']
    for i,i_deg in enumerate([20.928,24.518,28.164]):
    # for i,i_deg in enumerate([24.518,28.164]):
        sym=symbol[i]
    # i_deg = 28.164 # this is angle of incidence at MOHO
        incident_ray = calculate_incidence_vector(i_deg,phi_deg) # i and azimuth.

        # print('incidence angle at Moho=',i_deg,'Backazi=',phi_deg,'\n')
        # print('incident_ray_vector at Moho=',incident_ray,'\n')

        # STEP 1: calculate rap param and all things for flat Moho.
        normal = np.array([0, 0, 1])
        n1 = 1.0                             # Ri of mantle
        n2 = 1.38
        refracted_ray_flat, theta_i_flat, theta_r_flat, azimuth_flat,theta_i_surf_flat = snells_law(incident_ray, normal, n1, n2)
        ray_p_surf_flat=get_ray_param(theta_i_surf_flat,6371,5.8)
        # print('I_inci surf={:.2f}'.format(theta_i_surf_flat))

        print('###################################')
        # print('For a flat Moho:\n')
        # print('I_inci Moho={:.2f}'.format(theta_i_flat))
        # print('I_ref Moho={:.2f}'.format(theta_r_flat))
        # print('azimuth surf={:.2f}'.format(azimuth_flat))
        # print('I_inci surf={:.2f}'.format(theta_i_surf_flat))
        # print('Flat Moho Ray_p_surf={:.2f}'.format(ray_p_surf_flat))
        # print('Refracted ray=',refracted_ray_flat,'\n')
        # print('ray_p_surf,refracted_ray, theta_i, theta_r, azimuth,theta_i_surf')
        # print(f"{ray_p_surf:.2f},{refracted_ray}, {theta_i:.2f}, {theta_r:.2f}, {azimuth:.2f},{theta_i_surf:.2f}",'\n')

        ## STEP 2: calculate rap param and all things for tilted Moho.

        # normal = np.array([0.05, 0.0, 1])         # 0.1 is 50km moho change over 500km i.e., 5.7 degree #
        normal = np.array([0.0, 0.05, 1])         # 0.1 is 50km moho change over 500km i.e., 5.7 degree #
        i_deg_norm, phi_deg_norm = extract_angles(normal)
        # print('titlted Moho i and phi={:.2f}, {:.2f}'.format(i_deg_norm, phi_deg_norm))
        n1 = 1.0                             # Ri of mantle
        n2 = 1.38                           # n2 = v1/v2
        refracted_ray, theta_i, theta_r, azimuth,theta_i_surf = snells_law(incident_ray, normal, n1, n2)
        ray_p_surf=get_ray_param(theta_i_surf,6371,5.8)
        # print('###################################')
        # print('For a tilted Moho normal = {}:\n'.format(normal))
        # print('I_inci Moho={:.2f}'.format(theta_i))
        # print('I_ref Moho={:.2f}'.format(theta_r))
        # print('azimuth surf={:.2f}'.format(azimuth))
        # print('I_inci surf={:.2f}'.format(theta_i_surf))
        # print('Ray_p_surf={:.2f}'.format(ray_p_surf))
        # print('Refracted ray=',refracted_ray,'\n')
        # print('###################################')

        # time_max=np.sin((np.deg2rad(theta_i_surf)))*240/5.8
        # sys.exit()
        i_deg_r, phi_deg_r = extract_angles(refracted_ray)
        phi_rad_r = np.radians(phi_deg_r)
        i_rad_r = np.radians(i_deg_r)

        vp_cr = 5.8
        time_arr = []

        for station in stations_array:
            x, y = station
            # from Schweitzer et al. (2012); eq 9.2. time=(-x*sin(phi) - y*cos(phi))/(Vc*sin(i))
            time = np.sin(i_rad_r)*(-x * np.sin(phi_rad_r) - y * np.cos(phi_rad_r)) / (vp_cr )
            time_arr.append(time)
        ###

        time_arr = np.array(time_arr)
        # plot_time_arr(stations_array,time_arr,np.max(time_arr))
        # sys.exit()
        ###
        #
        # bounds = [(0, 360), (0, 90)] # bounds for phi and i
        bounds = [(None,None), (0, 90)] # bounds for phi and i

        initial_guess = [5, 5] ## phi and i
        result = minimize(objective, initial_guess,method='Nelder-Mead', bounds=bounds)#m
        if np.round(result.fun,2) != 0:
            print('Optimization failed\n')

        baz_opt, i_opt = result.x
        ray_p_opt=get_ray_param(i_opt,6371,5.8)
        baz_opt=np.round(baz_opt,4)
        if baz_opt <0:
            baz_opt+=360

        # print('Estimated slow and baz for array \n')
        # print('Baz={:.2f}'.format(baz_opt),', Incidence={:.2f}, '.format(i_opt), 'RayP={:.2f}'.format(ray_p))
        ###
        delta_az=np.round(azimuth-phi_deg,4) # for single station
        # delta_az=np.round(baz_opt-phi_deg,4) # for array

        print('Delta azi={:.2f}'.format(np.abs(delta_az)))
        print('####################')
        if np.round(delta_az,2) >= 360:
            delta_az=delta_az-360
        ##
        #
        delta_ray=ray_p_opt-ray_p_surf_flat # for array
        delta_ray=ray_p_surf-ray_p_surf_flat # for single station..

        # ray_p_surf
        print('Delta rayP={:.2f}'.format(np.abs(delta_ray)))
        # print('####################')
        # ax1.set_facecolor(("beige",.2))

        ax.scatter(phi_deg,np.abs(delta_az),marker=sym,color='xkcd:metallic blue',alpha=.65)#,label='Plane Moho')
        ax1.scatter(phi_deg,np.abs(delta_ray),marker=sym,color='xkcd:burnt orange',alpha=.65)#,label='Plane Moho')


##
ax.set_xlabel('azimuth of ray ($^\circ$) at Moho')
ax.spines['left'].set_color("xkcd:metallic blue")
ax.yaxis.label.set_color("xkcd:metallic blue")
ax1.yaxis.label.set_color("xkcd:burnt orange")
ax.set_ylabel('delta azi ($^\circ$)')
ax1.set_ylabel('delta rayP (s/$^\circ$)')
ax.xaxis.set_minor_locator(MultipleLocator(5))
ax.xaxis.set_major_locator(MultipleLocator(10))
ax.set_title('P for i=15,17.5,20$^\circ$ at tilted Moho (i={:.2f}, az={:.2f})'.format(i_deg_norm, phi_deg_norm))
ax1.set_ylim([-.015, .30])


ax.grid(alpha=.3)
plt.savefig('Moho_tilt_2.8_az0_3is.png',dpi=300,bbox_inches='tight', pad_inches=0.1)
plt.show()
