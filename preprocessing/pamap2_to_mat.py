import string
import numpy as np
import pandas as pd
from scipy.io import savemat
from pathlib import Path
import re

# change this path to the path to the PAMAP2 dataset on your device
path = "J:/enee439d/datasets/PAMAP2_Dataset"

# do not change these paths, define input and output paths
# mirror WISDM format for output
hand_accel =  "/hand/accel"
hand_gyro = "/hand/gyro"
chest_accel =  "/chest/accel"
chest_gyro = "/chest/gyro"

data_folder = "/Protocol"
opt_folder = "/Optional"
mat_folder = "/mat"

output_paths = [hand_accel, hand_gyro, chest_accel, chest_gyro]


data_path =  Path(path + data_folder)
opt_data_path = Path(path + opt_folder)
mat_paths = [None]*4

for i, output_path in enumerate(output_paths):
    mat_paths[i] = Path(path + mat_folder + output_path)

# ensure necessary folder for output exist
for mat_path in mat_paths:
    if(not mat_path.exists()):
        mat_path.mkdir(parents=True)


# define indices for PAMAP2 columns of interest
timestamp_col = 0
activity_col = 1
imu_cols = np.asarray([1, 2, 3, 7, 8, 9], np.int32)
imu_hand_cols = imu_cols + 3
imu_chest_cols = imu_cols + 20

# define column names for output
output_col_names = ['SubjectID', 'ActivityCode', 'TimeStamp', 'X', 'Y', 'Z']
# define which columns to load into dataframe
col_indices = np.concatenate(([timestamp_col, activity_col], imu_hand_cols, imu_chest_cols))
# define masks for extracting components

components = output_col_names[-3:]
sensors = ['accel', 'gyro']
locations = ['hand', 'chest']
# generate a column name per component per sensor per location
col_names = [f'{location}_{sensor}_{component}' for location in locations for sensor in sensors for component in components]
col_names = ['TimeStamp', 'ActivityCode'] + col_names

# regex for matching column names
# include timestamp and activity in every resulting dataframe
base_pattern = f'({col_names[0]})|({col_names[1]})'
# filter dataframes by sensor and location
col_patterns = [f'{base_pattern}|({location}_{sensor}.*)' for location in locations for sensor in sensors]

# regex for getting subject id from file name
get_subject_id = re.compile('subject(\d+).dat')

# convert every file into a .mat then write files to output
for i, file in enumerate(data_path.glob('*.dat')):
    # subject id the same for the whole file, extract value from filename
    subject_id = get_subject_id.match(file.name).group(1)

    print(f'Loading data for subject: {subject_id}')
    
    raw_df = pd.read_csv(file, sep=' ', usecols=col_indices, names = col_names)

    # check for opt data, add it if it is found
    opt_file = opt_data_path / file.name
    if(opt_file.exists()) :
        print("Appending optional data")
        raw_df = pd.concat((raw_df, pd.read_csv(opt_file, sep=' ', usecols=col_indices, names = col_names)), axis='index')
    for j, col_pattern in enumerate(col_patterns):
        print(f'Processing {locations[j//2]}_{sensors[j%2]} for subject={subject_id} ')
        out_df = raw_df.filter(regex = col_pattern, axis = 'columns')
        full_len = out_df.shape[0]
        out_df = out_df.dropna(axis='index')
        trim_len = out_df.shape[0]
        if(trim_len != full_len):
            print(f'Removed {full_len - trim_len} missing rows in {output_paths[j]}')
        
        # prefixes were needed for full df, remove them in subsets to match the format used for WISDM
        out_df.rename(mapper=lambda x: x[-1] if (x[-2] == '_') else x, axis = 'columns', inplace=True)
        # add subject id column to match WISDM format
        out_df[output_col_names[0]] = subject_id
        # reorder to finish matching format
        out_df = out_df[output_col_names]
        # PAMAP2 defines 0 as a transient activity, remove to match WISDM which does not include transitions
        out_df.drop(out_df.loc[out_df[output_col_names[1]]==0].index, inplace=True)       
        
        # duplicate rows pollute data (overfit certain subjects), remove them
        out_df.drop_duplicates(inplace=True)

        # create dictionary for generating .mat file
        mdict = {'subject': subject_id, 'activity_data': {}}
        
        # split on activity
        split_dfs = out_df.groupby(out_df[output_col_names[1]].values)
        for g, df in split_dfs:
            # transform pd dataframe to dictionary that MATLAB will understand
            # this also forces column vectors (arbitrary)
            mdict['activity_data'][chr(ord('A') - 1 + g)] = {name: np.reshape(col.values, (len(col.values), 1)) for name, col in df.items()}

        # print saving then overwrite once saved, fancy!
        print("saving", end = '\r')
        savemat(mat_paths[j] / (file.stem + '.mat'), mdict)
        print("saved ")