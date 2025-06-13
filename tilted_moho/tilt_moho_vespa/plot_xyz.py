# using xyz conatining time, slow, XF, baz, abs baz...plotting them as grid.
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from matplotlib import ticker
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib.colors as mcolors
from cmcrameri import cm
import matplotlib.cm as cmm
import xarray as xr
from matplotlib.colors import ListedColormap
import sys

###
def normalize(vector):
    return vector / np.linalg.norm(vector)
#
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
def get_gridded_time_vs(xf_xyz,param):
    if param=='slow':
        time, slowness, zval, baz, abs_baz = xf_xyz.T  # Unpacking columns # time, slow, coherence, baz, abs_baz
    else:
        time, baz, zval, slowness, abs_baz = xf_xyz.T  # Unpacking columns # time, baz, coherence, slow, abs_baz

    # Apply time filter to exclude the first and last 10 seconds..and removing slowness less than 1.
    time_mask = ((time >= time.min() + 10) & (time <= time.max() - 10)) & (slowness >= 1)

    # slow_mask = (slowness >= 1)

    # Filter time, slowness, baz, abs_baz, and zval accordingly
    filtered_time = time[time_mask]
    filtered_slowness = slowness[time_mask]
    filtered_baz = baz[time_mask]
    filtered_abs_baz = abs_baz[time_mask]
    filtered_zval = zval[time_mask]

    # Define grid resolution using filtered time
    time_grid = np.linspace(filtered_time.min(), filtered_time.max(), 500)
    slow_grid = np.linspace(filtered_slowness.min(), filtered_slowness.max(), 500)
    baz_grid = np.linspace(filtered_baz.min(), filtered_baz.max(), 500)

    # Create the meshgrid
    T, S = np.meshgrid(time_grid, slow_grid)
    _, B = np.meshgrid(time_grid, baz_grid)

    # Interpolate Z values onto the grids
    #Z_slow is a 2D array corresponding to the (T, S) grid.
    Z_slow = griddata((time, slowness), zval, (T, S), method='linear')
    Z_baz = griddata((time, baz), zval, (T, B), method='linear')

    return T,S,B,Z_slow,Z_baz,np.max(zval)
###
def create_vector(slow, theta_deg,v,r):
    # takes slow and baz and spits out the vector..p
    p=slow/np.deg2rad(1)
    i_rad=np.arcsin(p*v/r)
    # angles from degrees to radians
    theta_rad = np.radians(theta_deg)
    #
    z = np.cos(i_rad)
    x = np.sin(i_rad) * np.sin(theta_rad)
    y = np.sin(i_rad) * np.cos(theta_rad)
    vector = np.array([x, y, z])

    return vector
#
def extract_angles(vector):
    vector_normalized = vector / np.linalg.norm(vector)
    x, y, z = vector_normalized
    #  angle of incidence (i)
    i_rad = np.arccos(z)
    i_deg = np.degrees(i_rad)
    # azimuth (phi)
    phi_rad = np.arctan2(x, y) # phi is measured from north
    phi_deg = np.degrees(phi_rad)

    return np.absolute(i_deg), phi_deg
###
def snells_law(incident_ray,normal):

    n1 = 1.38   # 1 crust 2 is mantle
    n2 = 1.0

    # Normalize
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

    i=np.deg2rad(theta_i_surf)
    r_m=6336
    v_m=8
    p_m=r_m*np.sin(theta_i_surf)/v_m
    p_m_deg= np.round((p_m/np.rad2deg(1)),2)

    return refracted_ray, p_m_deg, azimuth, theta_r

###
incident_ray=create_vector(4.5,0,5.8,6371)
normal=np.array([0.05, 0, 1])
print('i/phi',extract_angles(incident_ray))
refracted_ray,theta_i_moho,azimuth, theta_r=snells_law(incident_ray,normal)
print('refracted_ray, theta_i_moho, azimuth, theta_r')
print(f"{refracted_ray}, {theta_i_moho:.2f}, {azimuth:.2f},{theta_r:.2f}",'\n')
###########

#following bit takes the baz grid file and converts it into a new Baz for a tilted Moho.
xyz_folder = "/Users/keyser/Research/AK_all_stations/sac_files_with_P/220526_120223_SA_inc2_r2.5/xyz_folder/"
for grid_number in [93]:
    xf_slow_xyz= np.loadtxt(xyz_folder+'xf_slow_xyz_93_220526_.05_.5.xyz') # time, slow, coherence, baz, abs_baz
    xf_baz_xyz= np.loadtxt(xyz_folder+'xf_baz_xyz_93_220526_.05_.5.xyz')

T,S,B,_,Z_baz,max_coher = get_gridded_time_vs(xf_baz_xyz,'baz')
_,_,_,Z_slow,_,_ = get_gridded_time_vs(xf_slow_xyz,'slow')

B_moho = np.zeros_like(B)  # Initialize transformed B grid
S_moho = np.zeros_like(S)  # Initialize transformed B grid


# #Iterate over all grid points

for i in range(B.shape[0]):
    for j in range(B.shape[1]):
        incident_ray = create_vector(S[i, j], B[i, j], 5.8, 6371)  # Create incident ray vector
        refracted_ray, p_m_deg, azimuth, theta_r = snells_law(incident_ray, normal)
        B_moho[i, j] = azimuth  # Store converted azimuth
        S_moho[i, j] = p_m_deg  # Store converted ray_P


# sys.exit()
##
cmap_lip=cm.oslo_r
colA = cmap_lip(np.arange(cmap_lip.N))
colA[:,-1] = np.linspace(0.65, 1, cmap_lip.N) #replaces the last column (alpha channel) in colA with values from 0.75 to 1, creating a gradient in transparency across the color map.

sm_alpha = ListedColormap(colA)

plot_amp_factor=3
# xyz_folder = "/Users/keyser/Research/AK_all_stations/sac_files_with_P/220526_120223_SA_inc2_r2.5/xyz_folder/"
# for grid_number in [93]:
#     xf_slow_xyz= np.loadtxt(xyz_folder+'xf_slow_xyz_93_220526_.05_.5.xyz') # time, slow, coherence, baz, abs_baz
#     xf_baz_xyz= np.loadtxt(xyz_folder+'xf_baz_xyz_93_220526_.05_.5.xyz')
#
# T,S,B,Z_slow,Z_baz,max_coher = get_gridded_time_vs(xf_slow_xyz,'slow')
# Plot
######## Begin figure
fig, axs = plt.subplots(2, 2, figsize=(16, 12), sharex=True)
(ax1, ax2), (ax3, ax4) = axs
# Time vs Slowness
c1 = ax1.contourf(T, S, Z_slow, levels=15, cmap=sm_alpha,vmin=0, vmax=max_coher/plot_amp_factor)
ax1.contour(T, S, Z_slow, levels=np.linspace(max_coher/8, max_coher/plot_amp_factor, 4), cmap='Greys_r', linewidths=0.65)
ax1.set_ylabel("Slowness")
ax1.set_title('Vespa original')
# fig.colorbar(c1, ax=ax1, label="Z value")
# Time vs Baz
c2 = ax3.contourf(T, B, Z_baz, levels=15, cmap=sm_alpha,vmin=0, vmax=max_coher/plot_amp_factor)
ax3.contour(T, B, Z_baz, levels=np.linspace(max_coher/8, max_coher/plot_amp_factor, 4), cmap='Greys_r', linewidths=0.65)
ax3.axhline(y=0, color='darkred', linestyle='--')
ax3.set_xlabel("Time")
ax3.set_ylabel("Back Azimuth")
# fig.colorbar(c2, cax=ax5, label="Z value",location='bottom',orientation='horizontal',fraction=.2,shrink=.1)
ax5=fig.add_axes([0.02, .3, .01, 0.3]) #[left, bottom, width, height]
sm = plt.cm.ScalarMappable(norm=plt.Normalize(vmin=0, vmax=max_coher/plot_amp_factor),cmap=sm_alpha )
sm.set_array(np.arange(0,1.33))
cbar = plt.colorbar(sm,cax=ax5,orientation='vertical',extend='max')
cbar.set_label('Coherence', rotation=90, labelpad=2,loc='center')
###
for ax in [ax1,ax2,ax3,ax4]:
    ax.grid(which='minor',axis='x',color='dimgrey', linestyle='--',linewidth=.65,alpha=.75)
    ax.grid(which='major',axis='both',color='dimgrey', linestyle='--',linewidth=.75,alpha=.75)
##########
########
######## same but for baz file..
# sys.exit()

# Time vs Slowness
c1 = ax2.contourf(T, S, Z_slow, levels=15, cmap=sm_alpha,vmin=0, vmax=max_coher/plot_amp_factor)
ax2.contour(T, S, Z_slow, levels=np.linspace(max_coher/8, max_coher/plot_amp_factor, 4), cmap='Greys_r', linewidths=0.65)
ax2.set_ylabel("Slowness")
ax2.set_title('Vespa rotated for normal:{}'.format(normal))

# fig.colorbar(c1, ax=ax1, label="Z value")
# Time vs Baz
c2 = ax4.contourf(T, B_moho, Z_baz, levels=15, cmap=sm_alpha,vmin=0, vmax=max_coher/plot_amp_factor)
ax4.contour(T, B_moho, Z_baz, levels=np.linspace(max_coher/8, max_coher/plot_amp_factor, 4), cmap='Greys_r', linewidths=0.65)
ax4.axhline(y=0, color='darkred', linestyle='--')
ax4.set_xlabel("Time")
ax4.set_ylabel("Back Azimuth")
# fig.colorbar(c2, ax=ax4, label="Z value",location='bottom',orientation='horizontal')
ax6=fig.add_axes([0.95, .3, .01, 0.3]) #[left, bottom, width, height]
sm = plt.cm.ScalarMappable(norm=plt.Normalize(vmin=0, vmax=max_coher/plot_amp_factor),cmap=sm_alpha )
# sm.set_array(np.arange(0,1.33))
cbar = plt.colorbar(sm,cax=ax6,orientation='vertical',extend='max')
cbar.set_label('Coherence', rotation=90, labelpad=2,loc='center')
set_locators(ax1, 'slow')
set_locators(ax2, 'slow')
set_locators(ax3, 'baz')
set_locators(ax4, 'baz')
ax4.set_ylim([-50, 50])
# plt.savefig('vespa_{}_compare.png'.format(grid_number),dpi=300,bbox_inches='tight', pad_inches=0.1)
# plt.tight_layout()
plt.show()
