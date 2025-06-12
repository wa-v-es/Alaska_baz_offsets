import os
import glob
import shutil
import sys
import re
###
##
def count_jpg_files_in_folders(parent_folder):
    # Initialize a dictionary to store the count of .jpg files in each folder
    jpg_counts = {}

    # Loop through each item in the parent folder
    for root, dirs, files in os.walk(parent_folder):
        # Loop through each file in the current folder
        for file in files:
            # Check if the file has a .jpg extension (case insensitive)
            if file.lower().endswith('.jpg'):
                # Get the current folder name
                folder_name = os.path.basename(root)
                # Increment the count for the current folder
                if folder_name not in jpg_counts:
                    jpg_counts[folder_name] = 0
                jpg_counts[folder_name] += 1

    return jpg_counts
###
#this copies all py_picks .dat file and creates and save them all in a new_folder
def find_and_copy_files(parent_folder, pattern, new_folder_name):
    # Create the new folder if it doesn't exist
    new_folder_path = os.path.join(parent_folder, new_folder_name)
    os.makedirs(new_folder_path, exist_ok=True)

    # Loop through each folder in the parent folder
    for root, dirs, files in os.walk(parent_folder):
        # Skip the new folder if it exists within the parent folder
        if new_folder_path in root:
            continue

        # Find files matching the pattern
        for file in glob.glob(os.path.join(root, pattern)):
            # Copy each file to the new folder
            shutil.copy(file, new_folder_path)
#####
## used this to find the average number of seismic sattions in each array for an earthquake..
def find_and_average_numbers(parent_folder, pattern):
    # Initialize a dictionary to store the extracted numbers by folder
    folder_numbers = {}

    # Loop through each folder in the parent folder
    for root, dirs, files in os.walk(parent_folder):
        # Find files matching the pattern
        numbers = []
        for file in glob.glob(os.path.join(root, pattern)):
            # Extract the number after 'num' in the file name
            match = re.search(r'_num(\d+)_', file)
            if match:
                numbers.append(int(match.group(1)))

        # Calculate and save the average if numbers were found
        if numbers:
            folder_numbers[root] = sum(numbers) / len(numbers)
            # folder_numbers[root] =  len(numbers)


    return folder_numbers
###

####
# Path to the parent folder
parent_folder = '/Users/keyser/Research/AK_all_stations/sac_files'

# Get the count of .jpg files in each folder
jpg_counts = count_jpg_files_in_folders(parent_folder)

# Print the counts
count_sum=0

for folder, count in jpg_counts.items():
    if folder[0]=='2':
        print(f'{folder}: {count} .jpg files')
        count_sum=count_sum+count
print(f'Total vespas:{count_sum}')

sys.exit()
##########
# second bit to copy all py picks
parent_folder = '/Users/keyser/Research/AK_all_stations/sac_files'

# Pattern to search for
pattern = 'grid_num_*.dat'

# Name of the new folder to copy files into
new_folder_name = 'py_pics_all'

# Call the function to find and copy files
find_and_copy_files(parent_folder, pattern, new_folder_name)
###############

#third bit to find the number of seismic stations in each array
######
parent_folder = '/Users/keyser/Research/AK_all_stations/sac_files'
for grid_num in [76,77,79,102,103,104,106,107,109,125,127,147,150,167,168,169]:
# for grid_num in [104,103]:

    print('------------\n')
    print('Doing grid_number=',grid_num,'\n')
    pattern = 'Vespapack_gridnum{}*_num*_-1.jpg'.format(grid_num)
    # Call the function to find files and calculate averages
    averages = find_and_average_numbers(parent_folder, pattern)
    # print(averages)
    for folder, avg in averages.items():
        print(f'Folder: {folder}, Average: {avg:.2f}')

# Print the average numbers
