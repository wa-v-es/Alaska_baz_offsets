# reads text files and calculates the max topo change for a sub-array
import numpy as np
import os
import sys
import glob as glob
import re


main_folder='/Users/keyser/Research/AK_all_stations/sac_files/230828_195530_PA_inc2_r2.5/st_files/'
pattern = os.path.join(main_folder, '*.txt')
# List to store the gridnum values
topo_diff_max = []
for file_path in glob.glob(pattern):
    match = re.search(r'number(\d+)', file_path)
    num=int(match.group(1))
    alts=[]
    for line in open(file_path,'r'):
        line=line.split()
        alts.append(float(line[8]))

    topo_diff_max.append(max(alts) - min(alts))
        # sys.exit()
print(topo_diff_max)

## avg topo change in arrays = 1560 mt 
## max topo change = 2450 mt
## std topo chnage = 600 mt
