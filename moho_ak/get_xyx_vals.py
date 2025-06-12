import numpy as np
import miller_alaskamoho_srl2018 as alaskamoho
import sys
import matplotlib.pyplot as plt
import seaborn as sns
###
moho_model = alaskamoho.MohoModel_opt
moho_data_filename="AlaskaMohoOpt"

lons = np.linspace(-165,-135,100)
lats = np.linspace(55, 72, 100)
lons=np.arange(-170,-135,.2)
lats=np.arange(50, 75,.2)


reg_lons, reg_lats = np.meshgrid(lons, lats)

quality  = moho_model.quality_at_lonlat_degrees(reg_lons, reg_lats, order=1)
reg_moho = moho_model.value_at_lonlat_degrees(reg_lons, reg_lats, order=1)
grad_moho = moho_model.slope_at_lonlat_degrees(reg_lons, reg_lats, order=1)


# grad_moho=alaskamoho.surface_model.slope_at_lonlat_degrees(alaskamoho,reg_lons, reg_lats, order=1)

a = np.transpose(np.stack((reg_lons.reshape(-1),
                           reg_lats.reshape(-1),
                           reg_moho.reshape(-1),
                           grad_moho.reshape(-1),
                           quality.reshape(-1))))

# a = a[a[:, 4] != 0] # removes rows with quality =0

# np.savetxt("{}.1-RegGrid.XYZ".format(moho_data_filename), a, fmt='%.2f %.2f %.2f %.2f %.2f',header="Longitude Latitude Depth gradient Quality")
np.savetxt("{}_forGrad_.2.XYZ".format(moho_data_filename), a, fmt='%.2f %.2f %.2f %.2f %.2f',header="Longitude Latitude Depth gradient Quality")

sys.exit()
##
mean_max='mean'
low_high='low'
# filename='MohoSlope_block_{}_{}_slow_PA_new.png'.format(mean_max,low_high)
block_file='/Users/keyser/Research/AK_all_stations/baz_offset_plots/block_mean_vals/block_{}_{}_slow_PA.txt'.format(mean_max,low_high)
combo_file='ALL_low_slow_offset.txt'
# block file conatins long lat mean_bazoffset std min_baz max_baz
matched_rows=[]
# lines=[]
for line in open(block_file,'r'):
    line=line.split()
    # lines.append(line)
    # Find the row in the filtered_array where column 1 (index 0) equals lon and column 2 (index 1) equals lat
    matched_row = a[(a[:, 0] == float(line[0])) & (a[:, 1] == float(line[1]))]
    if matched_row.size > 0:
        line_data = [float(line[2]), float(line[3])]

        # Combine the matched row and line[2], line[3] into one array
        combined_row = np.hstack((matched_row[0], line_data))
    matched_rows.append(combined_row)
##
matched_rows=np.array(matched_rows)
# matched row format: long lat moho grad quality mean_baz_off std_baz_off
grad_pt = matched_rows[:, 3]
moho_pt = matched_rows[:, 2]

mean_baz_pt = matched_rows[:, 5]

# Create the scatter plot:
fig = plt.figure(figsize=(6,6))
sns.set_style("darkgrid", {"axes.facecolor": ".9"})
sns.set_context("notebook")
ax1  = plt.subplot(111)
for i in range(len(grad_pt)):
    if mean_baz_pt[i] <0:
        ax1.scatter(moho_pt[i], mean_baz_pt[i],marker='o', facecolor='rebeccapurple', alpha=.75)
    else:
        ax1.scatter(moho_pt[i], mean_baz_pt[i],marker='o', facecolor='seagreen', alpha=.75)

# plt.rcParams['axes.labelsize'] = 15
ax1.set_ylabel('Block Mean Baz offset (s/$^\circ$)')
# ax1.set_xlabel('Moho grad.')
ax1.set_xlabel('Moho (km)')

# plt.show()
# fig.savefig('grad_v_mean_baz', dpi=300,bbox_inches='tight', pad_inches=0.1)
fig.savefig('moho_v_mean_baz', dpi=300,bbox_inches='tight', pad_inches=0.1)
#####
##
#same thing but plots of indivdual mean offsets..
mean_max='mean'
low_high='low'
# filename='MohoSlope_block_{}_{}_slow_PA_new.png'.format(mean_max,low_high)
# block_file='/Users/keyser/Research/AK_all_stations/baz_offset_plots/block_mean_vals/block_{}_{}_slow_PA.txt'.format(mean_max,low_high)
combo_file='/Users/keyser/Research/AK_all_stations/sac_files/max_vals_coherence/ALL_{}_slow_offset.txt'.format(low_high)
# combo file conatins long lat baz_max_coherence mean of peaks
matched_rows=[]
a[:, 0] = np.round(a[:, 0], 1)
a[:, 1] = np.round(a[:, 1], 1)
lines=[]
for line in open(combo_file,'r'):
    line=line.split()
    lines.append(line)
    # Find the row in the filtered_array where column 1 equals lon and column 2 equals lat
    matched_row = a[(a[:, 0] == round(float(line[0]), 1)) & (a[:, 1] == round(float(line[1]), 1))]
    if matched_row.size > 0:
        line_data = [float(line[2]), float(line[3])]

        # Combine the matched row and line[2], line[3] into one array
        combined_row = np.hstack((matched_row[0], line_data))
    matched_rows.append(combined_row)
##
matched_rows=np.array(matched_rows)
if len(matched_rows) != len(lines):
    print('NOT ALL offsets retuned moho/moho grad')
# matched row format: long lat moho grad quality mean_baz_off std_baz_off
grad_pt = matched_rows[:, 3]
moho_pt = matched_rows[:, 2]

mean_baz_pt = matched_rows[:, 6]
max_baz_pt = matched_rows[:, 5]


# Create the scatter plot:
fig = plt.figure(figsize=(6,6))
sns.set_style("darkgrid", {"axes.facecolor": ".9"})
sns.set_context("notebook")
ax1  = plt.subplot(111)
for i in range(len(grad_pt)):
    if max_baz_pt[i] <0:
        ax1.scatter(grad_pt[i], max_baz_pt[i],marker='o', facecolor='rebeccapurple', alpha=.75)
    else:
        ax1.scatter(grad_pt[i], max_baz_pt[i],marker='o', facecolor='seagreen', alpha=.75)

# plt.rcParams['axes.labelsize'] = 15
ax1.set_ylabel('Indi Max Baz offset (s/$^\circ$)')
ax1.set_xlabel('Moho grad.')
# ax1.set_xlabel('Moho (km)')

# plt.show()
# fig.savefig('grad_v_mean_baz', dpi=300,bbox_inches='tight', pad_inches=0.1)
fig.savefig('mohograd_v_max_baz_indi.png', dpi=300,bbox_inches='tight', pad_inches=0.1)
