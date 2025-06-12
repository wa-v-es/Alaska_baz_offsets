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
tomo_grad_pos_direc='/Users/keyser/Research/TOMOGRAD-main/2.5deg_grad_posit_grds/'
tomo_grad_direc='/Users/keyser/Research/TOMOGRAD-main/2.5deg_grad_grds/'


# folder_pattern_pa contains all eqs now
# folder_pattern_pa = "/Users/keyser/Research/AK_all_stations/sac_files/max_vals_coherence/*maxVals_low_slow.txt"
folder_pattern_pa = "/Users/keyser/Research/AK_all_stations/sac_files_with_P/max_vals_coherence/*_maxVals_low_slow.txt"

# folder_pattern_sa = "/Users/keyser/Research/AK_all_stations/sac_files/max_vals_coherence/SA/*maxVals_low_slow.txt"
matching_files_pa = sorted(glob.glob(folder_pattern_pa))
print('length of files in low slow folder=',len(matching_files_pa),'\n')
# matching_files_sa = sorted(glob.glob(folder_pattern_sa))

# grad_grd=pygmt.load_dataarray(tomo_grad_direc+'grad_posit_mask.60.norm.grd')

# a was long, lat, moho, moho grad, quality

matched_rows=[]
# lines=[]
nan_count=0

for eq_file in matching_files_pa:
    print('Doing eq_file:',eq_file,'\n')
# sys.exit()
#contains grid number; baz at max coherence va; mean of peaks;std of peaks; array cen lat long elevation; eq lat long elev; distance; baz; # stations in array# slow maxmin# baz max-min#.

    for line in open(eq_file,'r'):
        line=line.split()
        line = [float(item) for item in line]
        # print(line)
        bazi=line[11]
        bazi=round(bazi / 10) * 10 # round to nearest 10
        bazi+=180
        # Convert to -180 to 180 range
        bazi_converted = ((bazi + 180) % 360) - 180
        if bazi_converted == -180:
            bazi_converted = 180
        # grad_grd=pygmt.load_dataarray(tomo_grad_direc+'grad_posit_mask.{}.norm.grd'.format(bazi_converted)) # load grd file for the right baz
        grad_grd=pygmt.load_dataarray(tomo_grad_direc+'grad_mask.{}.norm.grd'.format(bazi_converted)) # load grd file for the right baz

        if line[2] > 7.5 or line[2] < - 7.5 or line[3] > 5: # removes large baz offset mean values and large std peaks
            continue
        # lines.append(line)
        # Find the row in the filtered_array where column 1 (index 0) equals lon and column 2 (index 1) equals lat
        # line [5,4] in array centre lon, lat
        # grad_val=-2
        try:
            grad_val = grad_grd.sel(lon=round(line[5], 1), lat=round(line[4], 1),method='nearest').item()
        except:
            print('Array centre does not have a grad value\n')

        # matched_row = a[(a[:, 0] == float(line[0])) & (a[:, 1] == float(line[1]))]
        if ~np.isnan(grad_val): # if not nan..
            # line_data = [grad_val, line]
            grad_val=round(grad_val, 4)
            combined_row=np.hstack((grad_val, line))

            matched_rows.append(combined_row)
        else:
            print('Grad_val=',grad_val,'\n')
            nan_count+=1
    # print(matched_rows)
    ##

print('nan count=',nan_count,'\n')
print('Len of  matched rows =',len(matched_rows),'\n')

matched_rows=np.array(matched_rows)
# matched row format:
grad_pt = matched_rows[:, 0]
# moho_pt = matched_rows[:, 2]
mean_baz_pt = matched_rows[:, 3]
sts_in_array = matched_rows[:, 13]
min_max_baz = matched_rows[:, 14]

# matched_rows format
# contains grad val; grid number; baz at max coherence va; mean of peaks;std of peaks; array cen lat long elevation; eq lat long elev; distance; baz; # stations in array.
matched_rows_neg_off=[]
matched_rows_pos_off=[]
# Create the scatter plot:
plt.ion()
fig = plt.figure(figsize=(8,5))
sns.set_style("darkgrid", {"axes.facecolor": ".9"})
sns.set_context("notebook")
ax1  = plt.subplot(111)
# originally was S=15 and not defined for pos and neg.
# s=.95*min_max_baz[i]
for i in range(len(grad_pt)):
    if -0.5 < mean_baz_pt[i] < 0.5:
        ax1.scatter(grad_pt[i], mean_baz_pt[i],marker='o',s=15,facecolor='black', edgecolor='white',alpha=.45,linewidth=.95)
    if mean_baz_pt[i] <-0.49:
        matched_rows_neg_off.append(matched_rows[i])
        if min_max_baz[i] >20:
            ax1.scatter(grad_pt[i], mean_baz_pt[i],marker='o' , facecolor='seagreen', edgecolor='seagreen',alpha=.51,linewidth=.65)
        else:
            ax1.scatter(grad_pt[i], mean_baz_pt[i],marker='o' , facecolor='none', edgecolor='seagreen',alpha=1,linewidth=.65)

        # ax1.errorbar(grad_pt[i], mean_baz_pt[i], yerr=min_max_baz[i],ecolor='brown',marker='.', alpha=.15,markerfacecolor='white', markeredgecolor='white',markersize=.04,linestyle='none')

    if mean_baz_pt[i] >0.49:
        matched_rows_pos_off.append(matched_rows[i])
        if min_max_baz[i] >20:
            ax1.scatter(grad_pt[i], mean_baz_pt[i],marker='o', facecolor='rebeccapurple', edgecolor='rebeccapurple',alpha=.51,linewidth=.65)
        else:
            ax1.scatter(grad_pt[i], mean_baz_pt[i],marker='o', facecolor='none', edgecolor='rebeccapurple',alpha=1,linewidth=.65)
######
# mean_baz_pt[mean_baz_pt > 0.49]

matched_rows_neg_off=np.array(matched_rows_neg_off)
print('Len of neg grad matched rows =',len(matched_rows_neg_off),'\n')
Grad_mean_neg_off=np.mean(matched_rows_neg_off[:,0])
Grad_median_neg_off=np.median(matched_rows_neg_off[:,0])

matched_rows_pos_off=np.array(matched_rows_pos_off)
print('Len of pos grad matched rows =',len(matched_rows_pos_off),'\n')

ax1.scatter(np.median(matched_rows_neg_off[:,0]), np.median(matched_rows_neg_off[:,3]),marker='d',s=95, facecolor='maroon', edgecolor='white',alpha=.75,linewidth=1.5)
ax1.scatter(np.median(matched_rows_pos_off[:,0]), np.median(matched_rows_pos_off[:,3]),marker='d',s=95, facecolor='maroon', edgecolor='white',alpha=.75,linewidth=1.5,label='Median')

ax1.scatter(np.mean(matched_rows_neg_off[:,0]), np.mean(matched_rows_neg_off[:,3]),marker='d',s=95, facecolor='gold', edgecolor='white',alpha=.65,linewidth=1.5,label='Mean')
ax1.scatter(np.mean(matched_rows_pos_off[:,0]), np.mean(matched_rows_pos_off[:,3]),marker='d',s=95, facecolor='gold', edgecolor='white',alpha=.65,linewidth=1.5)
plt.legend()
# plt.rcParams['axes.labelsize'] = 15
ax1.set_ylabel('Indi Mean Baz offset (s/$^\circ$)')
ax1.set_xlabel('Directional Moho gradient')
plt.show()
# fig.savefig('moho_v_mean_baz_-90deg', dpi=300,bbox_inches='tight', pad_inches=0.1)
# sys.exit()

#####
##
###curve FITTING

# popt, _ = curve_fit(objective, np.concatenate((matched_rows_neg_off[:,0], matched_rows_pos_off[:,0])), np.concatenate((matched_rows_neg_off[:,3], matched_rows_pos_off[:,3])))
y=np.concatenate((matched_rows_neg_off[:,3], matched_rows_pos_off[:,3]))
x=np.concatenate((matched_rows_neg_off[:,0], matched_rows_pos_off[:,0]))
slope, intercept, r_value, p_value, std_err = linregress(x, y)
print(f"Slope: {slope}")
print(f"Intercept: {intercept}")
print(f"R-squared: {r_value**2}")
print(f"P-value: {p_value}")
print(f"Standard Error: {std_err}")
# define new input values
x_new = np.linspace(-.6,1,50)
# unpack optima parameters for the objective function
# a, b, c = popt
# use optimal parameters to calculate new values
y_new = slope * x_new + intercept
str='y= {:.2f} x+ {:.2f}'.format(slope,intercept)
ax1.text(.7,-3.5,str, rotation=-12)
ax1.text(.8,2.3,'R-sq={:.2f}'.format(r_value**2))
ax1.text(.8,3.3,'P-val={:.1e}'.format(p_value))


plt.plot(x_new,y_new,ls='--',c='grey',alpha=.85,linewidth=1.5)
plt.show()
fig.savefig('all_P_grad_v_baz_mean_low_slow_0deg.png', dpi=300,bbox_inches='tight', pad_inches=0.1)

# fig.savefig('grad_v_baz_mean_high_slow_90deg_sts_arr.png', dpi=300,bbox_inches='tight', pad_inches=0.1)

sys.exit()
residuals = y - slope*x - intercept

plt.scatter(x, residuals,alpha=.5,c='cadetblue')
plt.axhline(0, color='firebrick', linestyle='--')
plt.xlabel("Gradient")
plt.ylabel("Residuals")
plt.title("Residual Plot")
plt.show()
rmse = np.sqrt(np.mean(residuals**2))
mae = np.mean(np.abs(residuals))

print(f"RMSE: {rmse}")
print(f"MAE: {mae}")
