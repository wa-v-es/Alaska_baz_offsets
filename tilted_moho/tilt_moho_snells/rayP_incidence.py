# using taup, finds the ray parameter difference between consecutive incidence angle (distances) for P and PP.
#and plots them.
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from obspy.taup.taup_geo import calc_dist,calc_dist_azi
from obspy.taup import TauPyModel
import xkcd as xk
import seaborn as sns
###
from obspy.taup import TauPyModel
import logging

logging.getLogger("matplotlib.font_manager").disabled = True

model = TauPyModel(model="iasp91")

dists = np.arange(60, 100, 2.5)

ray_params = []
incidences = []
distances = []
ray_params_PP = []
incidences_PP = []
distances_PP = []
for dist in dists:
    arrivals = model.get_travel_times(source_depth_in_km=100, distance_in_degree=dist, phase_list=["P",'PP'])[0]
    ray_params.append(arrivals.ray_param_sec_degree)
    incidences.append(arrivals.incident_angle)
    distances.append(dist)
    ##
    arrivals_PP = model.get_travel_times(source_depth_in_km=100, distance_in_degree=dist, phase_list=["P",'PP'])[1]
    ray_params_PP.append(arrivals_PP.ray_param_sec_degree)
    incidences_PP.append(arrivals_PP.incident_angle)
    distances_PP.append(dist)

# Calculate the differences in ray parameters for consecutive incidence angles
ray_param_differences = [j - i for i, j in zip(ray_params[:-1], ray_params[1:])]
ray_param_differences_PP = [j - i for i, j in zip(ray_params_PP[:-1], ray_params_PP[1:])]

# Plot incidence vs ray param difference
with plt.xkcd():
    fig, ax1 = plt.subplots(figsize=(15, 6))
    # plt.xkcd()
    # sns.set_theme(context='paper', style='darkgrid', palette='deep', font='sans-serif', font_scale=1.8, color_codes=True)
    # sns.set_style("darkgrid", {"grid.color": ".8", "grid.linestyle": ":"})
    ax1.plot(distances[:-1], ray_param_differences, marker='o', linestyle='-', color='xkcd:pinkish tan',alpha=1,label='P')
    ax1.plot(distances[:-1], ray_param_differences_PP, marker='*', markersize=10,linestyle='-', color='xkcd:poop green',alpha=.9,label='PP')

    plt.xlabel('Distance ($^\circ$)')
    plt.ylabel('Ray Param diff. (s/$^\circ$)')
    # ax1.set_title('Incidence vs Ray Parameter Difference')
    ax1.spines['right'].set_color('none')
    ax1.spines['top'].set_color('none')
    ax1.grid(True,alpha=.75)
    plt.legend()

    # Create a second x-axis on the top for distances
    # ax2 = ax1.twiny()
    # ax2.set_xlim(ax1.get_xlim())
    # # ax2.set_xticks(np.linspace(12, 22, 11))
    # ax2.set_xticks(np.linspace(ax1.get_xlim()[0], ax1.get_xlim()[1], 35))
    # ax2.set_xticklabels([f'{dist:.1f}' for dist in distances[:-1]])
    # ax2.set_xlabel('Distance (degrees)')

# plt.show()
plt.savefig('dist_rayP_xkcd.png',dpi=300,bbox_inches='tight', pad_inches=0.1)
sys.exit()
xk.XKCDify(ax, xaxis_loc=0.0, yaxis_loc=1.0,f1=50, f2=.005, f3=15,
        xaxis_arrow='+-', yaxis_arrow='+-')
plt.show()
