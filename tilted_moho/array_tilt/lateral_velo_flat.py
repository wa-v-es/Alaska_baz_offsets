# this tests the result of a lateral velocity heterogeniety for a flat Moho
import numpy as np
import matplotlib.pyplot as plt
import sys
from cmcrameri import cm
from scipy.optimize import minimize
import math
import matplotlib
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from obspy.taup import TauPyModel
import matplotlib.colors as colors
from scipy.optimize import minimize_scalar
from curve_wave_arr import normalize,get_ray_param,calculate_incidence_vector,extract_angles,snells_law
####
def objective(params):
    phi, i = params
    phi_rad = np.radians(phi)
    i_rad = np.radians(i)

    total_error = 0

    for idx, station in enumerate(stations_array):
        x, y = station
        predicted_time = np.sin(i_rad)*(-x * np.sin(phi_rad) - y * np.cos(phi_rad)) / (5.5)
        total_error += (predicted_time - time_arr[idx]) ** 2
    return total_error
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
###
stations_array = np.array(stations)
model = TauPyModel(model="iasp91")

n1 = 1.0                             # Ri of mantle
n2 = 1.38
vp_cr_1 = 5.5
vp_cr_2 = 6

plt.rcParams.update({'font.size': 16})

fig = plt.figure(figsize=(10, 5))
ax = fig.add_subplot(111)
ax1=ax.twinx()
ax.set_xlabel('backazimuth of ray ($^\circ$) at Moho')
ax.set_ylabel('$\delta$ backazimuth ($^\circ$)',fontweight='bold')
ax1.set_ylabel('$\delta$ rayP (s/$^\circ$)',fontweight='bold')

ax.spines['left'].set_color("xkcd:denim blue")
ax.yaxis.label.set_color("xkcd:denim blue")
ax.tick_params(axis='y', colors="xkcd:denim blue")

ax1.spines['right'].set_color("xkcd:burnt orange")
ax1.yaxis.label.set_color("xkcd:burnt orange")
ax1.tick_params(axis='y', colors="xkcd:burnt orange")

ax.xaxis.set_minor_locator(MultipleLocator(5))
ax.xaxis.set_major_locator(MultipleLocator(10))

# ax1 = fig.add_subplot(212)
matplotlib.rcParams.update({'font.size': 15})

for phi_deg in range(0,95,5):
    for i,i_deg in enumerate([18.7]):#[20.928,24.518,28.164]
        # i_deg=20
        incident_ray = calculate_incidence_vector(i_deg,phi_deg) # i and azimuth.

        # print('incidence angle at Moho=',i_deg,'Backazi=',phi_deg,'\n') # print('incident_ray_vector at Moho=',incident_ray,'\n')

        # STEP 1: calculate rap param and all things for flat Moho.
        normal = np.array([0, 0, 1])

        refracted_ray_flat, theta_i_flat, theta_r_flat, azimuth_flat,theta_i_surf_flat = snells_law(incident_ray, normal, n1, n2)
        ray_p_surf_flat=get_ray_param(theta_i_surf_flat,6371,5.8)
        print('###################################')
        print('ray_p_surf_flat,refracted_ray_flat, theta_i_flat, theta_r_flat, azimuth_flat,theta_i_surf_flat')
        print(f"{ray_p_surf_flat:.2f},{refracted_ray_flat}, {theta_i_flat:.2f}, {theta_r_flat:.2f}, {azimuth_flat:.2f},{theta_i_surf_flat:.2f}",'\n')


        i_deg_r, phi_deg_r = extract_angles(refracted_ray_flat)
        phi_rad_r = np.radians(phi_deg_r)
        i_rad_r = np.radians(i_deg_r)

        time_arr = []
        # time_arr_v2 = []



        # loop for time arrival for plane wave
        #calculated at array centre
        for station in stations_array:
            x, y = station
            # from Schweitzer et al. (2012); eq 9.2. time=(-x*sin(phi) - y*cos(phi))/(Vc*sin(i))
            if x>0:
                time = np.sin(i_rad_r)*(-x * np.sin(phi_rad_r) - y * np.cos(phi_rad_r)) / (vp_cr_1)
            else:
                time = np.sin(i_rad_r)*(-x * np.sin(phi_rad_r) - y * np.cos(phi_rad_r)) / (vp_cr_2)

            time_arr.append(time)

        time_arr = np.array(time_arr)

        bounds = [(None,None), (0, 90)] # bounds for phi and i

        initial_guess = [5, 5] ## phi and i
        result = minimize(objective, initial_guess,method='Nelder-Mead', bounds=bounds)#m
        if np.round(result.fun,1) > 1:
            print('Optimization failed\n')

        baz_opt, i_opt = result.x
        ray_p_opt=get_ray_param(i_opt,6371,5.8)
        baz_opt=np.round(baz_opt,4)
        if baz_opt <0:
            baz_opt+=360

        print('Estimated slow and baz for array \n')
        print('Baz={:.2f}'.format(baz_opt),', Incidence={:.2f}, '.format(i_opt), 'RayP={:.2f}'.format(ray_p_opt))
        delta_az=np.round(baz_opt-phi_deg,3)

        print('Delta azi={:.3f}'.format(np.abs(delta_az)))
        if np.round(delta_az,2) >= 360:
            delta_az=delta_az-360

        delta_ray=ray_p_opt-ray_p_surf_flat
        print('Delta rayP={:.3f}'.format(np.abs(delta_ray)))
        print('###################################\n')

        ax.scatter(phi_deg,np.abs(delta_az),color='xkcd:denim blue',alpha=.75)#,label='Plane Moho')
        ###
        ax1.scatter(phi_deg,np.abs(delta_ray),color='xkcd:burnt orange',alpha=.75)#,label='Plane Moho')

ax.grid(alpha=.6)
# plt.show()
plt.savefig('fig2_b_color_spine_lateral_velocity.png',dpi=300,bbox_inches='tight', pad_inches=0.1)

        # sys.exit()


        ###
