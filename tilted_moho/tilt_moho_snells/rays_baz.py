# script to find interation of incident ray with tilted moho

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys
import math
import matplotlib
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

##
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
def normalize(vector):
    return vector / np.linalg.norm(vector)
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
    #
    # #angle of incidene at syrface for refracted ray
    # surface  = np.array([0, 0, 1])
    # cos_theta_i_surf = np.dot(refracted_ray, surface)
    # theta_i_surf = np.arccos(cos_theta_i_surf)
    theta_i_surf,azimuth=extract_angles(refracted_ray)


    return refracted_ray, math.degrees(theta_i), math.degrees(theta_r), azimuth, theta_i_surf

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
###
def plot_plane(normal, point, ax, size=10):
    # Create a grid of points
    d = -point.dot(normal)
    xx, yy = np.meshgrid(range(-size, size), range(-size, size))

    # Calculate corresponding z values for the plane
    zz = (-normal[0] * xx - normal[1] * yy - d) * 1. / normal[2]

    ax.plot_surface(xx, yy, zz, color='y', alpha=0.25)

def get_incidence_angle(p,r,v):
    p=p/np.deg2rad(1)
    i=np.arcsin(p*v/r)
    return np.round(np.rad2deg(i),2)
##
def get_ray_param(theta,r,v):
    i=np.deg2rad(theta)
    p=r*np.sin(i)/v
    return (p*np.deg2rad(1))

# math.atan(x)
# math.degrees(math.atan(.36))
##
# i=get_incidence_angle(7.58,6371-35,8) # incidence angle at moho (35km). Vp mantle = 8km/sec. Rayp in sec/deg.

## steepest Eq 96 deg & 600 km deep: P (PP) 13.41 (23.29) incidence at surface with rayP 4.445 (7.58). Incidence angle at Moho 18.76 (33.25)

## shallowest eq 65deg & 60 km deep: P (PP) 19.81 (27.11) incidence at surface with rayP 6.5 (8.74). Incidence angle at Moho 28 (39.2)
##
# For P, incident angle at Moho varies from 18.7 to 28 degree.
# for PP, incident angle at Moho varies from 33.25 to 39.2 degree.

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(211)
ax1 = fig.add_subplot(212)
matplotlib.rcParams.update({'font.size': 12})
ax.set_facecolor(("beige",.25))
ax1.set_facecolor(("beige",.25))
plt.ion()

colours=['royalblue','xkcd:slate green','xkcd:dark pink']
for i,az in enumerate([0,45,90]):
    for inci_angle in np.linspace(18.7,28,20): # for P
    # for inci_angle in [.36]:
        incident_ray = calculate_incidence_vector(inci_angle,az) # i and azimuth.
        normal = np.array([0.05, 0.0, 1])         # 0.1 is 50km moho change over 500km i.e., 5.7 degree #
        n1 = 1.0                             # Ri of mantle
        n2 = 1.38                           # n2 = v1/v2

        refracted_ray, theta_i, theta_r, azimuth,theta_i_surf = snells_law(incident_ray, normal, n1, n2)

        ax.scatter(inci_angle,azimuth-az,color=colours[i])
    # ax.scatter(inci_angle,azimuth-az,color='xkcd:slate green',label='P')

    for inci_angle in np.linspace(33.25,39.2,10): # for PP
    # for inci_angle in [.36]:
        incident_ray = calculate_incidence_vector(inci_angle,az) # i and azimuth.
        normal = np.array([0.05, 0.0, 1])         # 0.1 is 50km moho change over 500km i.e., 5.7 degree #
        n1 = 1.0                             # Ri of mantle
        n2 = 1.38                           # n2 = v1/v2

        refracted_ray, theta_i, theta_r, azimuth,theta_i_surf = snells_law(incident_ray, normal, n1, n2)

        ax.scatter(inci_angle,azimuth-az,color=colours[i])
    # ax.scatter(inci_angle,azimuth-az,color='xkcd:faded pink',label='PP')
    ax.set_title('Moho tilted 2.86$^\circ$ towards X-axis')
    ax.text(30,2.5,s='0$^\circ$')
    ax.text(30,1.5,s='45$^\circ$')
    ax.text(30,0,s='90$^\circ$')

##### following is ax2 scatter..

    for inci_angle in np.linspace(18.7,28,20): # for P
    # for inci_angle in [.36]:
        incident_ray = calculate_incidence_vector(inci_angle,az) # i and azimuth.
        normal = np.array([-0.1, 0, 1])         # [-0.1, 0, 1]...0.1 is 50km moho change over 500km i.e., 5.7 degree #
        n1 = 1.0                             # Ri of mantle
        n2 = 1.38                           # n2 = v1/v2

        refracted_ray, theta_i, theta_r, azimuth,theta_i_surf = snells_law(incident_ray, normal, n1, n2)
        ax1.scatter(inci_angle,azimuth-az,color=colours[i])
    # ax1.text(inci_angle,azimuth+.5,s='Ray from az=90. Normal = {}.PP: normal = [0.5, 0, 1]'.format(normal))
    ax1.set_title('Moho tilted -5.7$^\circ$ towards negative X-axis. For PP, Moho has no tilt.')
    temp1,temp2=inci_angle,azimuth-az
    for inci_angle in np.linspace(33.25,39.2,10): # for PP
    # for inci_angle in [.36]:
        incident_ray = calculate_incidence_vector(inci_angle,az) # i and azimuth.
        normal = np.array([-0, 0.0, 1])         # 0.1 is 50km moho change over 500km i.e., 5.7 degree #
        n1 = 1.0                             # Ri of mantle
        n2 = 1.38                           # n2 = v1/v2

        refracted_ray, theta_i, theta_r, azimuth,theta_i_surf = snells_law(incident_ray, normal, n1, n2)
        ax1.scatter(inci_angle,azimuth-az,color=colours[i])
    ax1.text(30,-5,s='0$^\circ$')
    ax1.text(30,-3.8,s='45$^\circ$')
    ax1.text(30,-.1,s='90$^\circ$')

# ax1.scatter(inci_angle,azimuth-az,color='xkcd:faded pink',label='PP')
# ax1.scatter(temp1,temp2,color='xkcd:slate green',label='P')

#########
ax.annotate(    "",
        xy=(28.5, .8),
        xytext=(18.5, .8),
        arrowprops=dict(arrowstyle='|-|', alpha=.5, lw=1.15, mutation_scale=10),)
ax.text(23.5,.8,s='P',fontweight=550,color='white',bbox={'facecolor': 'black', 'pad': 5,'alpha': 0.85},fontsize=14)
##
ax.annotate(    "",
        xy=(39, .8),
        xytext=(33, .8),
        arrowprops=dict(arrowstyle='|-|', alpha=.5, lw=1.15, mutation_scale=10),)
ax.text(35.8,.8,s='PP',fontweight=550,color='white',bbox={'facecolor': 'black', 'pad': 5,'alpha': 0.85},fontsize=14)

# plt.legend(fontsize="12")
ax1.set_xlabel('Incidence angle ($^\circ$) at Moho ')
ax.set_ylabel('$\delta$ backazimuth ($^\circ$)')
ax1.set_ylabel('$\delta$ backazimuth ($^\circ$)')

ax.xaxis.set_minor_locator(MultipleLocator(2.5))
ax1.xaxis.set_minor_locator(MultipleLocator(2.5))
ax.yaxis.set_minor_locator(MultipleLocator(0.5))
ax1.yaxis.set_minor_locator(MultipleLocator(.5))
# ax1.set_title('(2.9 deg) towards from incoming ray')

# ax.grid(alpha=.28)
# ax1.grid(alpha=.28)

plt.show()
# print('incident_ray:',incident_ray,'\n')
# print('normal:',normal,'\n')
# fig.savefig('2.9_degree_moho_azi_new.png', dpi=300,bbox_inches='tight', pad_inches=0.1)

sys.exit()
incident_ray= normalize(incident_ray)
# Visualization
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
# normal=[.2,.2,1]

ax.quiver(-incident_ray[0], -incident_ray[1], -incident_ray[2], incident_ray[0], incident_ray[1], incident_ray[2], color='cadetblue', label='Incident Ray')
ax.quiver(0, 0, 0, normal[0], normal[1], normal[2], color='purple', label='Normal')
if refracted_ray is not None:
    ax.quiver(0, 0, 0, refracted_ray[0], refracted_ray[1], refracted_ray[2], color='indianred', label='Refracted Ray')

    print('refracted_ray:',refracted_ray)
else:
    print("Total internal reflection occurred.")
plot_plane(normal, np.array([0, 0, 0]), ax)

ax.set_xlim([-1.5, 1.5])
ax.set_ylim([-1.5, 1.5])
ax.set_zlim([-1.5, 1.5])
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
# ax.legend()
plt.show()
#
# Print
print(f"Incidence Angle: {math.degrees(theta_i)}")
print(f"Refracted Angle: {math.degrees(theta_r)}")
print(f"Azimuth Angle: {math.degrees(azimuth)}")
print(f"Incidence Angle surface: {math.degrees(theta_i_surf)}")
