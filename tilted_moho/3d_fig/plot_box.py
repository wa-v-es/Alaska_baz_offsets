import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys
import matplotlib.pyplot as plt
import math
# import xkcd as xk

###
def normalize(vector):
    return vector / np.linalg.norm(vector)
def plot_plane(normal, point, ax, size=10):
    # Create a grid of points
    d = -point.dot(normal)
    xx, yy = np.meshgrid(range(-size, size), range(-size, size))

    # Calculate corresponding z values for the plane
    zz = (-normal[0] * xx - normal[1] * yy - d) * 1. / normal[2]

    ax.plot_surface(xx, yy, zz, color='y', alpha=0.25)
######
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
    phi_rad = np.arctan2(x, y) # phi is measured from north
    phi_deg = np.degrees(phi_rad)

    return i_deg, phi_deg
###
def plot_ray(ax,incident_ray,refracted_ray,moho_x, moho_y, moho_z,raycolor):
    incident_length = 5
    incident_ray_scaled = incident_ray * incident_length

    # Incident ray: from below to the point on Moho
    incident_start = np.array([moho_x, moho_y, moho_z]) - incident_ray_scaled
    incident_end = np.array([moho_x, moho_y, moho_z])

    ##Draw incident ray
    ax.plot(
        [incident_start[0], incident_end[0]],
        [incident_start[1], incident_end[1]],
        [incident_start[2], incident_end[2]],
        color='steelblue', linewidth=1.5)
    #
    refracted_length = 5
    refracted_end = incident_end + refracted_ray * refracted_length

    ax.plot(
        [incident_end[0], refracted_end[0]],
        [incident_end[1], refracted_end[1]],
        [incident_end[2], refracted_end[2]],
        color=raycolor, linewidth=1.5,ls='-')
######
def plot_normal(normal,raycolor):
    n1 = 5 * normal
    origin = np.array([0, 0, 0])

    ax.quiver(
    origin[0], origin[1], origin[2],
    n1[0], n1[1], n1[2],
    color=raycolor,
    linewidths=1.5,
    arrow_length_ratio=0.1,
    length=1.0)

incident_ray = calculate_incidence_vector(20,-90) # i and azimuth.
incident_ray_2 = calculate_incidence_vector(20,0) # i and azimuth.

incident_ray = calculate_incidence_vector(0,0) # i and azimuth.


# x_ray=np.tan(np.deg2rad(inci_angle))
# incident_ray = np.array([x_ray, 0, 1])  # 0.36 is 20 deg incidence angle for parallel plane
normal = np.array([-.087, 0, 1])         # 0.1 is 50km moho change over 500km i.e., 5.7 degree #
n1 = 1.0                             # Ri of mantle
n2 = 1.38                           # n2 = v1/v2

refracted_ray, theta_i, theta_r, azimuth,theta_i_surf = snells_law(incident_ray, normal, n1, n2)
refracted_ray_2, _, _, _,_ = snells_law(incident_ray_2, normal, n1, n2)

#Create figure and 3D axis
plt.rcParams.update({'font.size': 14})

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')

# Grid for surface
x = np.linspace(-5, 5, 7)
y = np.linspace(-5, 5, 7)
X, Y = np.meshgrid(x, y)
Z_top = np.full_like(X, 5)  # Top surface
Z_flat = np.full_like(X, 0)  # Top surface


# Draw top surface (transparent)
ax.plot_surface(X, Y, Z_top, alpha=0.1, color='gray', edgecolor='k')
ax.view_init(elev=15, azim=-54)
# Hide grid lines
ax.grid(False)

# Hide axes ticks
ax.set_xticks([])
ax.set_yticks([])
ax.set_zticks([])
# Plot seismic stations in a cross pattern
station_coords = []
for xi in x:
    station_coords.append((xi, 0))
for yi in y:
    station_coords.append((0, yi))

station_coords = list(set(station_coords))  # Remove (0,0) duplicate
for x_s, y_s in station_coords:
    ax.scatter(x_s, y_s, 5.2, color='black', marker='v', s=50)

ax.scatter(0, 0, 5.2, color='darkred', marker='v', s=60)

# Draw Moho interface
Z_moho = -0.087 * X  # 5 degree dip in x direction: tan(5 deg) ~ 0.087
ax.plot_surface(X, Y, Z_moho, alpha=0.25, color='brown', edgecolor='none')
ax.plot_surface(X, Y, Z_flat, alpha=0.25, color='steelblue', edgecolor='none')
### draw normals

plot_normal(np.array([.087, 0, 1]),'brown')
plot_normal(np.array([0, 0, 1]),'steelblue')


######## plot ray

moho_x, moho_y = 0, 0
moho_z = -0.087 * moho_x
moho_z = 0

# plot_ray(ax,incident_ray,incident_ray,moho_x, moho_y, moho_z,'steelblue')# plots normal..
# plot_ray(ax,incident_ray,refracted_ray,moho_x, moho_y, moho_z,'brown')

colours=['royalblue','xkcd:slate green','xkcd:dark pink']

ax.quiver(0, 4, -6, 0, -4, 0, color='royalblue', arrow_length_ratio=0.16)
ax.text(0, 4.2, -6, '0$^\circ$',color='royalblue', fontsize=20)

ax.quiver(2.8, 2.8, -6, -2.8, -2.8, 0, color='mediumseagreen', arrow_length_ratio=0.16)
ax.text(2.9, 2.9, -6, '45$^\circ$',color='mediumseagreen', fontsize=20)

ax.quiver(4, 0, -6, -4, 0, 0, color='mediumvioletred', arrow_length_ratio=0.16)
ax.text(4.2, 0, -6, '90$^\circ$',color='mediumvioletred', fontsize=20)

# Draw coordinate axes
ax.quiver(-5.2, -5.2, 4.9, 3, 0, 0, color='k', arrow_length_ratio=0.12)
ax.text(-2.2, -5.2, 4.9, 'X', fontsize=17)
ax.quiver(-5.2, -5.2, 4.9, 0, 3, 0, color='k', arrow_length_ratio=0.12)
ax.text(-5.2, -2.2, 4.9, 'Y (North)', fontsize=17)
ax.quiver(-5.2, -5.2, 4.9, 0, 0, 3, color='k', arrow_length_ratio=0.12)
ax.text(-5.2, -5.2, 7.9, 'Z', fontsize=17)

# Set labels and view
# ax.set_xlabel('X')
# ax.set_ylabel('Y')
# ax.set_zlabel('Z')
ax.set_xlim([-6, 6])
ax.set_ylim([-6, 6])
# ax.set_zlim([-6.1, 6])

# plt.title('Seismic Array and Moho Interface with Refracted Rays')
plt.tight_layout()

# plt.savefig('3d_box_threecolor.png', dpi=300,bbox_inches='tight', pad_inches=0.01)
# plt.show()
