import string
import numpy as np
import pandas as pd
from scipy.io import savemat
from pathlib import Path

# change this path to the path to the WISDM dataset on you device
path = "J:/enee439d/datasets/wisdm-dataset"

# do not change these paths
watch_accel =  "/watch/accel"
watch_gyro = "/watch/gyro"
phone_accel =  "/phone/accel"
phone_gyro = "/phone/gyro"
raw = "/raw"
mat = "/mat"

data_folders = [watch_accel, watch_gyro, phone_accel, phone_gyro]
raw_paths = [None]*4
mat_paths = [None]*4
col_names = ['SubjectID', 'ActivityCode', 'TimeStampNanos', 'X', 'Y', 'Z']

def remove_semi(string):
    if (string[-1] == ';'):
        return string[0:-1]
    else:
        return string

# ensure necessary folders for output exist
for i, data_folder in enumerate(data_folders):
    raw_paths[i] = Path(path + raw + data_folder)
    mat_paths[i] = Path(path + mat + data_folder)
    
for path in mat_paths:
    if(not path.exists()):
        path.mkdir(parents=True)

# iterate over all the raw files
for i, directory in enumerate(raw_paths):
    # convert every file into a .mat then write files to output
    for file in directory.glob('*.txt'):
        raw_df = pd.read_csv(file, names = col_names, converters = {col_names[-1] : remove_semi})
        
        # subject id the same for the whole file, extract the value from the first row
        subject_id = raw_df.head(1)[col_names[0]].values[0]
        mdict = {'subject': subject_id, 'activity_data': {}}
        print(f'Processing {directory.parent.name}_{directory.name} for subject={subject_id} ')
        # split on activity
        split_dfs = raw_df.groupby(raw_df[col_names[1]].values)
        for g, df in split_dfs:
            # transform pd dataframe to dictionary that MATLAB will understand
            mdict['activity_data'][g] = {name: col.values for name, col in df.items()}
        print("saving")
        savemat(mat_paths[i] / (file.stem + '.mat'), mdict)
        print("\r", end="")
        print("saved")