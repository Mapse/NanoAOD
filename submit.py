import os, subprocess, sys
import datetime

################### To edit #####################

# File with the dataset paths
das_path = 'charmonium_dataset_2016.txt'

################### End editing #####################

# Crab template name
crab_template = "crab_config"
# Step
config_template = "nanoanalyzer"

# Read the file and put the paths on a list
with open(das_path) as f:
    datasets = f.read().splitlines()

# Loop on each dataset path to submit to crab.    
for dataset in datasets:
    if dataset.startswith('#'): continue
    # Split the path.
    # Ex: /Charmonium/Run2018A-12Nov2019_UL2018_rsb-v1/AOD -> ['', 'Charmonium', 'Run2018A-12Nov2019_UL2018_rsb-v1', 'AOD']
    splited = dataset.split('/')
    # Stores dataset name
    PRIM_DATASET = splited[1]
    # Stores era name
    ERA = splited[2].split('-')[0]
    # Stores workarea name. Ex: CharmoniumRun2018UL 
    WORKAREA = PRIM_DATASET + ERA[0:-1] + "UL"
    # Stores directory name. Ex: 'data18UL'
    OUT_DIR = "Data" + ERA[-3:-1] + "UL"
    # Stores request name. Ex : 'CharmoniumRun2018A_Run2018A-12Nov2019_UL2018_rsb-v1_AOD'
    REQUEST_NAME = PRIM_DATASET + ERA + "_" + splited[2][9:] + "_" + splited[-1]

    with open(crab_template + ".py", 'r') as f:
        new_file = f.read().replace("WORKAREA", WORKAREA)
        new_file = new_file.replace("OUT_DIR", OUT_DIR)
        new_file = new_file.replace("REQUEST_NAME", REQUEST_NAME)
        new_file = new_file.replace("INPUT_DATASET", dataset)
        new_file = new_file.replace("CONFIG_NAME", config_template + "_" + ERA + "_" + splited[2][9:])

    with open(crab_template + "_" + ERA + "_" + splited[2][9:] + ".py", 'w') as nf:
        nf.write(new_file)

    with open(config_template + ".py", 'r') as f:
        new_file = f.read().replace("REQUEST_NAME", REQUEST_NAME)

    with open(config_template + "_" + ERA + "_" + splited[2][9:] + ".py", 'w') as nf:
        nf.write(new_file)

    os.system("crab submit -c " + crab_template + "_" + ERA + "_" + splited[2][9:] + ".py")


