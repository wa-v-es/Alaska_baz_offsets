# for a plane and titlted moho, this script calculates the travel time difference for an array for changing azi.
# the travel time is calculated seperately for each station using different angle of incidence at each station. The script is written such that
# the function 'get_inci' seraches for a distance from array centre for distance between 60 and 70.
# as distance is degree needs to be estimated for that incidence angle.
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
from obspy.taup import TauPyModel
import matplotlib.colors as colors
from scipy.optimize import minimize_scalar
#########
def normalize(vector):
    return vector / np.linalg.norm(vector)
def plot_time_arr(stations_array,time_arr,az,plot=False):
    time=np.max(time_arr)
    fig,ax = plt.subplots(figsize=(6, 5))
    ax.set_facecolor(("xkcd:grey teal",.1))

    plt.scatter(stations_array[:, 0], stations_array[:, 1], c=time_arr,marker='^', cmap='PiYG', s=90,alpha=.88,lw=.5,edgecolors='k')
    # plt.scatter(stations_array[:, 0], stations_array[:, 1], c='grey',marker='^', s=90,alpha=.88,lw=.5,edgecolors='k')
    for i,st in enumerate(stations_array):
        x,y=st
        plt.text(x+5, y,i,fontsize=9,c='xkcd:dusk blue')
    plt.colorbar(label='Time of Arrival (s)',fraction=.12,shrink=.5)
    plt.xlabel('X (km)')
    plt.ylabel('Y (km)')
    plt.title('Max time ={:.2f}. phi={}'.format(time,az))
    if plot:
        plt.savefig('time_diff_plane_tilt/time_diff_{}'.format(az),dpi=300,bbox_inches='tight', pad_inches=0.1)
    else:
        plt.show()
###
def plot_time_arr_both(stations_array,time_arr_abs,time_arr_diff,az,save=False):
    time=np.max(time_arr_abs)
    fig= plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(121)

    ax.set_facecolor(("xkcd:grey teal",.1))
    im=ax.scatter(stations_array[:, 0], stations_array[:, 1], c=time_arr_abs,marker='^', cmap='PiYG', s=90,alpha=.88,lw=.5,edgecolors='k')
    # plt.scatter(stations_array[:, 0], stations_array[:, 1], c='grey',marker='^', s=90,alpha=.88,lw=.5,edgecolors='k')
    for i,st in enumerate(stations_array):
        x,y=st
        plt.text(x+5, y,i,fontsize=9,c='xkcd:dusk blue')
    fig.colorbar(im,ax=ax,label='Time titlted Moho (s)',fraction=.12,shrink=.5,location='bottom')
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_title('Max time ={:.2f}. phi={}'.format(time,az))
    ##
    time=np.max(time_arr_diff)
    ax1 = fig.add_subplot(122)
    ax1.set_facecolor(("xkcd:grey teal",.1))
    norm = colors.Normalize(vmin=-.8*time, vmax=.8*time)
    im1=ax1.scatter(stations_array[:, 0], stations_array[:, 1], c=time_arr_diff,marker='^', cmap='PiYG', s=90,alpha=.88,lw=.5,edgecolors='k',norm=norm)
    # plt.scatter(stations_array[:, 0], stations_array[:, 1], c='grey',marker='^', s=90,alpha=.88,lw=.5,edgecolors='k')
    for i,st in enumerate(stations_array):
        x,y=st
        plt.text(x+5, y,i,fontsize=9,c='xkcd:dusk blue')
    # ax1.colorbar(label='dt (tilt-plane) (s)',fraction=.12,shrink=.5)
    fig.colorbar(im1,ax=ax1,label='dt (Tilt - plane) Moho (s)',fraction=.12,shrink=.5,location='bottom')
    ax1.set_xlabel('X (km)')
    ax1.set_ylabel('Y (km)')

    ax1.set_title('Max time ={:.2f}. phi={}'.format(time,az))
    # fig.suptitle("", fontsize=13)
    if save:
        plt.savefig('corrected_time_tilt/time_diff_{}'.format(az),dpi=300,bbox_inches='tight', pad_inches=0.1)
        plt.close()
    else:
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
def objective_st(params):
    phi, i = params
    phi_rad = np.radians(phi)
    i_rad = np.radians(i)

    total_error = 0

    for idx, station in enumerate(stations_array):
        x, y = station
        predicted_time = np.sin(i_rad)*(-x * np.sin(phi_rad) - y * np.cos(phi_rad)) / (5.8)
        total_error += (predicted_time - time_arr_corrected[idx]) ** 2
    return total_error
###
def objective_taup_st(params):
    phi, i = params
    phi_rad = np.radians(phi)
    i_rad = np.radians(i)

    total_error = 0
    arrivals_centre = model.get_ray_paths(source_depth_in_km=100, distance_in_degree=63.929, phase_list=["P"])
    time_P_arr_centre=arrivals_centre[0].time

    for idx, station in enumerate(stations_array):

        x, y = station
        # predicted_time = np.sin(i_rad)*(-x * np.sin(phi_rad) - y * np.cos(phi_rad)) / (5.8)
        # total_error += (predicted_time - time_arr_corrected[idx]) ** 2

        incidence_angle_station,dist,time_st_taup=get_inci_per_station(station,phi_deg_r,i_deg)
        time_rel=time_st_taup-time_P_arr_centre
        print(f"Station {station}: Incidence Angle = {incidence_angle_station:.3f}°, Relative Time = {time_rel:.2f} s")

        time_arr_st_taup.append(time_rel)
    return total_error
###
def get_ray_param(theta,r,v):
    i=np.deg2rad(theta)
    p=r*np.sin(i)/v
    return np.round((p*np.deg2rad(1)),3)
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
def get_inci_surf(inci_moho):
    n1 = 1.0                      # Ri of mantle
    n2 = 1.38
    sin_theta_i=np.sin(np.deg2rad(inci_moho))
    sin_theta_r = (n1 / n2) * sin_theta_i
    theta_r = np.arcsin(sin_theta_r)
    return np.rad2deg(theta_r)
###
# def get_dist_centre(inci_moho): # searches thorugh distances to find distance that matches the requested angle of incidence.
#     inci_surf=get_inci_surf(inci_moho)
#
#     print('inci surf: {:.3f}'.format(inci_surf))
#     model = TauPyModel(model="iasp91")
#     eq_depth = 100
#     #
#     print('Searching for distance to eq....\n')
#     for distance in np.linspace(60, 70, 1500):  # Search
#         arrivals = model.get_ray_paths(source_depth_in_km=eq_depth, distance_in_degree=distance, phase_list=["P"])
#         if arrivals:
#             incidence_angle = arrivals[0].incident_angle
#             # if
#             if np.isclose(incidence_angle, inci_surf, atol=0.005):
#                 print('closest inci angle:{:.3f}'.format(incidence_angle))
#                 break
#     return distance
# #
def find_earthquake_distance(known_incidence_angle, station_depth=100, model="iasp91"):
    """
    Finds the epicentral distance (in degrees) for a given P-wave incidence angle at a station.
    Parameters:
        known_incidence_angle (float): The known incidence angle of the P-wave at the station (degrees).
    Returns:
        float: Epicentral distance in degrees.
    """
    taup_model = TauPyModel(model=model)
    # Function to minimize: absolute difference between known and computed incidence angle
    def incidence_angle_mismatch(distance_deg):
        arrivals = taup_model.get_ray_paths(source_depth_in_km=station_depth, distance_in_degree=distance_deg, phase_list=["P"])
        if not arrivals:
            return np.inf  # Avoid cases where no arrivals are found
        computed_incidence_angle = arrivals[0].incident_angle
        return abs(computed_incidence_angle - known_incidence_angle)

    # Optimize to find the distance where the incidence angle matches
    result = minimize_scalar(incidence_angle_mismatch, bounds=(0.1, 180), method="bounded")

    return result.x if result.success else None

def get_inci_per_station(station,eq_azi,angle_inci_surf):
    """
    Compute the angle of incidence and travel time for a given station
    considering its actual distance from the earthquake.
    Parameters:
        station (tuple): (x, y) coordinates of the station.
        eq_azi (float): Earthquake azimuth in degrees.
        eq_distance_center (float): Distance of the center station from the earthquake in km.
    Returns:
        tuple: (incident angle at station, distance to station in degrees, travel time)
    """
    x,y=station
    model = TauPyModel(model="iasp91")
    eq_distance=111.19*find_earthquake_distance(angle_inci_surf)
    # eq_distance=111.19*63.929 #distance from the center station (in km) for i = 20 at surface!
    eq_depth = 100
    eq_azi_rad = np.radians(eq_azi)
    ###
    distance_to_station_x = x - eq_distance * np.sin(eq_azi_rad)
    distance_to_station_y = y - eq_distance * np.cos(eq_azi_rad)
    distance_to_station = np.sqrt(distance_to_station_x ** 2 + distance_to_station_y ** 2)
    distance_to_station_deg = distance_to_station / 111.19

    arrivals = model.get_ray_paths(source_depth_in_km=eq_depth, distance_in_degree=distance_to_station_deg, phase_list=["P"])

    return arrivals[0].incident_angle,distance_to_station_deg,arrivals[0].time
##
def get_time_diff_crust_mantle(incident_ray,refracted_ray,moho_normal,stations_array):
    # from time_crust.py
    # incident_ray = calculate_incidence_vector(i_deg,phi) # i and azimuth.
    moho_base_depth = 35  # km
    v_crust = 5.8  # km/s
    v_mantle = 8.0
    # refracted_ray, _, _, _,_ = snells_law(incident_ray, moho_normal, n1, n2)
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
    return relative_arrivals

# eq_distance=111.19*get_dist_centre(angle_inci_moho)
# use this to find eq_distance for an angle of incidence..and then put the value in get_inci_per_station function
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

model = TauPyModel(model="iasp91")

# following bit is to test the distance and incidence angle variations across the array for a fixed angle of incidence.
# eq_distance = 111.19*63.929 #distance from the center station (in km) for i =20 at surface
# eq_depth = 100
# eq_azi = 0
# arrivals_centre = model.get_ray_paths(source_depth_in_km=eq_depth, distance_in_degree=63.929, phase_list=["P"])
# time_P_arr_centre=arrivals_centre[0].time
# eq_azi_rad = np.radians(eq_azi)
# incidence_angles = []
# dist=[]
#
# for station in stations_array:
#     incidence_angle,distance_to_station_deg,_=get_inci_per_station(station,eq_azi,20) # here 20 is only used if a line is not commneted in the fucntion.
#     incidence_angles.append(incidence_angle)
#     dist.append(np.round(distance_to_station_deg,3))
#
# # Convert incidence_angles to a numpy array
# incidence_angles_array = np.array(incidence_angles)
# dist = np.array(dist)

# print("Incidence angles of P waves at each station (in degrees):")
# print('Inci angle min/max={:.2f}/{:.2f}'.format(np.min(incidence_angles_array),np.max(incidence_angles_array)))
# print('Distance min/max={:.3f}/{:.3f}'.format(np.min(dist),np.max(dist)))
#
# ray_P_min=get_ray_param(np.min(incidence_angles_array),6371,5.8)
# ray_P_max=get_ray_param(np.max(incidence_angles_array),6371,5.8)
# print('RayP min/max={:.2f}/{:.2f}'.format(ray_P_min,ray_P_max))

# sys.exit()
# print("Stations coordinates (km):")
# print(stations_array)
##
# plt.ion()
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
# ax.set_facecolor(("xkcd:light grey blue",.1))

# sys.exit()
del_az_all=[]
del_slow_all=[]
n1 = 1.0                             # Ri of mantle
n2 = 1.38
vp_cr = 5.8
for phi_deg in range(0,95,5):

    # phi_deg = 5
    # i_deg = 28.164 # this is angle of incidence at MOHO..given i=20.0004 at surface
    # i_deg= 20.928
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

        ## STEP 2: calculate rap param and all things for tilted Moho.
        normal = np.array([0.1, 0.0, 1])         # 0.1 is 50km moho change over 500km i.e., 5.7 degree #
        # normal = np.array([0.0, 0.05, 1])         # 0.1 is 50km moho change over 500km i.e., 5.7 degree #
        i_deg_norm, phi_deg_norm = extract_angles(normal)
        print('titlted Moho i and phi={:.2f}, {:.2f}'.format(i_deg_norm, phi_deg_norm))

        refracted_ray, theta_i, theta_r, azimuth,theta_i_surf = snells_law(incident_ray, normal, n1, n2)
        ray_p_surf=get_ray_param(theta_i_surf,6371,5.8)
        print('###################################')
        print('For a tilted Moho normal = {}:\n'.format(normal))
        print('I_inci Moho={:.2f}'.format(theta_i))
        print('I_ref Moho={:.2f}'.format(theta_r))
        print('azimuth surf={:.2f}'.format(azimuth))
        print('I_inci surf={:.2f}'.format(theta_i_surf))
        print('Ray_p_surf={:.2f}'.format(ray_p_surf))
        print('Refracted ray=',refracted_ray,'\n')
        print('###################################')

        # time_max=np.sin((np.deg2rad(theta_i_surf)))*279/5.8
        # sys.exit()

        i_deg_r, phi_deg_r = extract_angles(refracted_ray_flat)
        phi_rad_r = np.radians(phi_deg_r)
        i_rad_r = np.radians(i_deg_r)

        time_arr = []
        time_arr_plane = []
        time_arr_st_taup = []


        # loop for time arrival for plane wave
        #calculated at array centre
        for station in stations_array:
            x, y = station
            # from Schweitzer et al. (2012); eq 9.2. time=(-x*sin(phi) - y*cos(phi))/(Vc*sin(i))
            time = np.sin(i_rad_r)*(-x * np.sin(phi_rad_r) - y * np.cos(phi_rad_r)) / (vp_cr)
            time_arr_plane.append(time)
            ###    # loop for time arrival for same angle of incidence for the entire array i.e., angle of incidence
                #calculated at array centre
        # loop for time arrival for titlted wave
        #calculated at array centre
        i_deg_r, phi_deg_r = extract_angles(refracted_ray)
        phi_rad_r = np.radians(phi_deg_r)
        i_rad_r = np.radians(i_deg_r)
        for station in stations_array:
            x, y = station
            # from Schweitzer et al. (2012); eq 9.2. time=(-x*sin(phi) - y*cos(phi))/(Vc*sin(i))
            time = np.sin(i_rad_r)*(-x * np.sin(phi_rad_r) - y * np.cos(phi_rad_r)) / (vp_cr)
            time_arr.append(time)

        ###
        ### next bit is to calculate time difference btw stations using  taup###
        # _,_,time_P_arr_centre=get_inci_per_station(stations_array[0],phi_deg_r,i_deg_r)
        # for station in stations_array:# uses distance for a st in array and calculates absolute time. relative time is saved.
        #     incidence_angle_station,dist,time_st_taup=get_inci_per_station(station,phi_deg_r,i_deg_r)
        #     time_rel=time_st_taup-time_P_arr_centre
        #     print(f"Station {station}: Incidence Angle = {incidence_angle_station:.3f}°, Relative Time = {time_rel:.2f} s")
        #
        #     time_arr_st_taup.append(time_rel)

        time_arr = np.array(time_arr)
        # time_arr_st_taup = np.array(time_arr_st_taup)


        # for (x, y), plane, tilt in zip(stations_array, time_arr_plane, time_arr):
        #     print(f"Station ({x:.1f}, {y:.1f}) -> Plane: {plane:.2f}s -> Tilt: {tilt:.2f} s")
        #
        # for (x, y), plane, tilt in zip(stations_array, time_arr_plane, time_arr_st_taup):
        #     print(f"Station ({x:.1f}, {y:.1f}) -> Plane: {plane:.2f}s -> Tilt_taup: {tilt:.2f} s")

        # for i,st in enumerate(stations_array):
        #     print('St# {}, time={:.3f}'.format(i,time_arr_st_taup[i]))

        # sys.exit()
        # plot_time_arr(stations_array,time_arr,phi_deg)
        # plot_time_arr(stations_array,time_arr_st_taup,phi_deg)

        # plot_time_arr_both(stations_array,time_arr,time_arr-time_arr_plane,phi_deg,save=False)
        # plot_time_arr_both(stations_array,time_arr_st_taup,time_arr_st_taup-time_arr_plane,phi_deg,save=False)

        time_crust_mantle = get_time_diff_crust_mantle(incident_ray,refracted_ray,normal,stations_array)

        time_arr_corrected=time_arr+time_crust_mantle
        # plot_time_arr_both(stations_array,time_arr_corrected,time_arr_corrected-time_arr_plane,phi_deg,save=True)

        # sys.exit()
        #
        # bounds = [(0, 360), (0, 90)] # bounds for phi and i
        bounds = [(None,None), (0, 90)] # bounds for phi and i

        initial_guess = [5, 5] ## phi and i
        result = minimize(objective, initial_guess,method='Nelder-Mead', bounds=bounds)#m
        result_st = minimize(objective_st, initial_guess,method='Nelder-Mead', bounds=bounds)#m

        if np.round(result.fun,1) > 1:
            print('Optimization failed\n')

        baz_opt, i_opt = result.x
        ray_p_opt=get_ray_param(i_opt,6371,5.8)
        baz_opt=np.round(baz_opt,4)
        if baz_opt <0:
            baz_opt+=360
        #
        if np.round(result_st.fun,1) > 1:
            print('St Optimization failed\n')

        baz_opt_st, i_opt_st = result_st.x
        ray_p_st=get_ray_param(i_opt_st,6371,5.8)
        baz_opt_st=np.round(baz_opt_st,4)
        if baz_opt_st <0:
            baz_opt_st+=360

        print('Estimated slow and baz for array \n')
        print('Baz={:.2f}'.format(baz_opt),', Incidence={:.2f}, '.format(i_opt), 'RayP={:.2f}'.format(ray_p_opt))
        delta_az=np.round(baz_opt-phi_deg,3)

        print('Delta azi={:.3f}'.format(np.abs(delta_az)))
        if np.round(delta_az,2) >= 360:
            delta_az=delta_az-360

        delta_ray=ray_p_opt-ray_p_surf_flat
        print('Delta rayP={:.3f}'.format(np.abs(delta_ray)))
        print('###################################\n')
        print('Estimated slow and baz for array st wise \n')
        print('Baz={:.2f}'.format(baz_opt_st),', Incidence={:.2f}, '.format(i_opt_st), 'RayP={:.2f}'.format(ray_p_st))
        delta_az_st=np.round(baz_opt_st-phi_deg,3)

        print('Delta azi={:.3f}'.format(np.abs(delta_az_st)))
        if np.round(delta_az_st,2) >= 360:
            delta_az_st=delta_az_st-360
        delta_ray_st=ray_p_st-ray_p_surf_flat
        print('Delta rayP={:.3f}'.format(np.abs(delta_ray_st)))

        # sys.exit()
        ax.scatter(phi_deg,np.abs(delta_az),color='xkcd:denim blue',alpha=.75)#,label='Plane Moho')
        ax1.scatter(phi_deg,np.abs(delta_ray),color='xkcd:burnt orange',alpha=.75)#,label='Plane Moho')
        ###
        ax.scatter(phi_deg,np.abs(delta_az_st),marker='+',s=190,color='xkcd:denim blue',alpha=.75)#,label='Plane Moho')
        ax1.scatter(phi_deg,np.abs(delta_ray_st),marker='+',s=190,color='xkcd:burnt orange',alpha=.75)#,label='Plane Moho')

        del_az_all.append(baz_opt_st-baz_opt)
        del_slow_all.append(ray_p_st-ray_p_opt)
##

ax.set_title('5.7$^\circ$ Moho tilted towards +X'.format(i_deg_norm, phi_deg_norm))
# ax1.set_ylim([-.015, .30])


ax.grid(alpha=.6)
# plt.savefig('fig2_b_color_spine.png',dpi=300,bbox_inches='tight', pad_inches=0.1)
plt.show()
