#plots baz offsets vs 95% constarined min/max baz and slow values (how weel constarined the vespas are)
import numpy as np
import miller_alaskamoho_srl2018 as alaskamoho
import sys
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import pygmt
import glob as glob
from cmcrameri import cm
from scipy.optimize import curve_fit
from scipy.stats import linregress
###

#change line 42 and add 'posit' in fig name and text name
# tomo_grad_pos_direc='/Users/keyser/Research/TOMOGRAD-main/2.5deg_grad_posit_grds/'
# tomo_grad_direc='/Users/keyser/Research/TOMOGRAD-main/2.5deg_grad_grds/'


# folder_pattern_pa contains all eqs now
folder_pattern_all = "/Users/keyser/Research/AK_all_stations/sac_files_with_P/max_vals_coherence/*_maxVals_low_slow.txt"

folder_pattern_pa = "/Users/keyser/Research/AK_all_stations/sac_files_with_P/max_vals_coherence/*PA_maxVals_low_slow.txt"
folder_pattern_sa = "/Users/keyser/Research/AK_all_stations/sac_files_with_P/max_vals_coherence/*SA_maxVals_low_slow.txt"

matching_files_pa = sorted(glob.glob(folder_pattern_pa))
print('length of files in low slow folder=',len(matching_files_pa),'\n')
matching_files_sa = sorted(glob.glob(folder_pattern_sa))
matching_files_all = sorted(glob.glob(folder_pattern_all))


#contains grid number; baz at max coherence va; mean of peaks;std of peaks; array cen lat long elevation; eq lat long elev; distance; baz; # stations in array# slow maxmin# baz max-min#.
plt.ion()
fig = plt.figure(figsize=(15,10))
sns.set_style("whitegrid")#, {"axes.facecolor": ".9"} {"grid.color": ".6", "grid.linestyle": ":"}
sns.set_style("whitegrid",{"axes.facecolor": ".3", "grid.linestyle": ":"})
# sns.set_context("notebook")
ax1  = plt.subplot(231)
ax2  = plt.subplot(232)
ax3  = plt.subplot(233)
ax4  = plt.subplot(234)
ax5  = plt.subplot(235)
ax6  = plt.subplot(236)

# axes = [ax1, ax4]
# for ax in axes:
for eq_file in matching_files_all:
        # print('Doing eq_file:',eq_file,'\n')
        for line in open(eq_file,'r'):
            line=line.split()
            line = [float(item) for item in line]
            # ax1.errorbar(line[-2], line[2], yerr=line[3],ecolor='brown',marker='.', alpha=.15,markerfacecolor='white', markeredgecolor='white',markersize=.04,linestyle='none')
            ax1.scatter(line[-2], line[2],marker='o',s=25,facecolor='white', edgecolor='white',alpha=.6,linewidth=.5)
            ax4.scatter(line[-1], line[2],marker='o',s=25,facecolor='white', edgecolor='white',alpha=.6,linewidth=.5)

ax1.set_xlabel('95% Baz variation ($^\circ$)')
ax1.set_ylabel('Mean Baz offset ($^\circ$)')
ax4.set_xlabel('95% slow variation (s/$^\circ$)')
ax4.set_ylabel('Mean Baz offset ($^\circ$)')
ax1.set_title('All low slowness eqs')
####
for eq_file in matching_files_pa:
    # print('Doing eq_file:',eq_file,'\n')
    for line in open(eq_file,'r'):
        line=line.split()
        line = [float(item) for item in line]
        # ax2.errorbar(line[-2], line[2], yerr=line[3],ecolor='brown',marker='.', alpha=.15,markerfacecolor='white', markeredgecolor='white',markersize=.04,linestyle='none')
        ax2.scatter(line[-2], line[2],marker='o',s=25,facecolor='white', edgecolor='white',alpha=.6,linewidth=.5)
        ax5.scatter(line[-1], line[2],marker='o',s=25,facecolor='white', edgecolor='white',alpha=.6,linewidth=.5)

ax2.set_title('Pacific eqs')
##
for eq_file in matching_files_sa:
    # print('Doing eq_file:',eq_file,'\n')
    for line in open(eq_file,'r'):
        line=line.split()
        line = [float(item) for item in line]
        # ax3.errorbar(line[-2], line[2], yerr=line[3],ecolor='brown',marker='.', alpha=.15,markerfacecolor='white', markeredgecolor='white',markersize=.04,linestyle='none')
        ax3.scatter(line[-2], line[2],marker='o',s=25,facecolor='white', edgecolor='white',alpha=.6,linewidth=.5)
        ax6.scatter(line[-1], line[2],marker='o',s=25,facecolor='white', edgecolor='white',alpha=.6,linewidth=.5)

ax3.set_title('SA eqs')
for ax in [ax1, ax2, ax3]:
    ax.set_xlim(5,35)
    ax.set_ylim(-15,15)
for ax in [ax4, ax5, ax6]:
    ax.set_xlim(0.4,3.2)
    ax.set_ylim(-15,15)
# ax3.set_ylim(-35,35)
# ax[i].set_xlim(5,35)


sys.exit()

fig.savefig('baz_slow_min_max_offsets.png', dpi=300,bbox_inches='tight', pad_inches=0.1)

sys.exit()
