# NanoAODPlus
This repository contains the codes for conversion from AOD or miniAOD to nanoAOD.

## Running official DPS/SPS

The first thing to do is to organize the text files containing the paths to the AOD datasets on Data Aggregation Service (DAS). The files are orgnized as:

### $DPS-c\overline{c}$

* **dps_ccbar_official_2016_prevfp.txt**
* **dps_ccbar_official_2016.txt**
* **dps_ccbar_official_2017.txt**
* **dps_ccbar_official_2018.txt**

### $DPS-b\overline{b}$

* **dps_bbbar_official_2016_prevfp.txt**
* **dps_bbbar_official_2016.txt**
* **dps_bbbar_official_2017.txt**
* **dps_bbbar_official_2018.txt**

### $SPS-c\overline{c}$

* **sps_ccbar_official_2016_prevfp.txt**
* **sps_ccbar_official_2016.txt**
* **sps_ccbar_official_2017.txt**
* **sps_ccbar_official_2018.txt**

### $SPS-b\overline{b}$

* **sps_bbbar_official_2016_prevfp.txt**
* **sps_bbbar_official_2016.txt**
* **sps_bbbar_official_2017.txt**
* **sps_bbbar_official_2018.txt**

For instance, the path for $DPS-b\overline{b}$ for $p_T$ 50-100 GeV/c is:

```
/D0ToKPi_Jpsi50to100_HardQCD_TuneCP5_13TeV-pythia8-evtgen/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM
```
After this step, you must modify three files before submitting your request:

* **crab_config.py**: edit lines 6 and 7.
* **nanoanalyzer_mc_official.py**: choose the correct global_tag and good_JSON. 
* **submit_mc_official.py**: edit input_file with the wanted .txt file.

Then, to call CRAB, do the following:

```
python submit_mc_official.py
```
## Running Private SPS

If you want to run a private sample you can perform the steps LHE -> ROOT (hadronization) -> GEN-SIM -> DIGI2RAW -> HLT -> AOD on this repository https://github.com/Mapse/SPS_MC/tree/main

To run the NanoAODPlus step for private production you will use three files:

* **crab_config_nano_MC.py**
* **nanoanalyzer_mc.py**
* **submit_nano_monte_carlo.py**

First, you need to get the .txt files with the paths to the AOD files. Go to /NanoAOD/NanoAnalyzer/paths_monte_carlo/.

The first thing to do here is to identify what is the path of your files. Using T2_Caltec_US as an example, the command is,

```
xrdfs k8s-redir.ultralight.org:1094 ls -l path_to_files
```
The <i>path_to_files</i> is the path where the files were produced. To identify it you refer to the variable **out_dir_base** on **crab_config_LHE.py** file (line 7). In this example, it is like this:

```
out_dir_base = '/store/group/uerj/' + getpass.getuser() + '/'
```
where in this case, **getpass.getuser()** is my CERN username, **mabarros**. Therefore, the first part of <i>path_to_files</i> is,

```
/store/group/uerj/mabarros
```
The next part of the path is given on the variable **output_dataset** (line 8)

```
output_dataset = 'CRAB_PrivateMC_RunII_UL_SPS_2017',
```
the next is given by the *.txt* file you created in the previous step (in this example, it is jpsi_ccbar_9to30_VFNS_SPS_2017_13TeV.txt), and the next is the number given by the CRAB job (in this case it is 230720_143214). Therefore, the <i>path_to_files</i> is, 

```
/store/group/uerj/mabarros/CRAB_PrivateMC_RunII_UL_SPS_2017/jpsi_ccbar_9to30_VFNS_SPS_2017_13TeV/230720_143214
```
Now that you have learned the format of the <i>path_to_files</i> we can produce the next **.txt** file that will be the input of the next step. Go to **SPS_MC/CMSSW_10_6_20_patch1/src/GS/path**
/

On the **get_files_xrootd.py** files, you will see four lists (after line 46). In this case, they will be like this:

```
mc = ['CRAB_PrivateMC_RunII_UL_SPS_2017'] -> same as output_dataset.
    
dataset = ['jpsi_ccbar_9to30_VFNS_SPS_2017_13TeV'] -> same as jpsi_ccbar_9to30_VFNS_SPS_2017_13TeV.txt

crab_folder = ['221113_230503'] -> from CRAB  

n_folders = [1] -> number of directories located on /store/group/uerj/mabarros/CRAB_PrivateMC_RunII_UL_SPS_2017/jpsi_ccbar_9to30_VFNS_SPS_2017_13TeV/230720_143214¹

¹The directories we refer to are those named 0000, 0001, 0002, ....

```
Then, you do the following:

```
 python3 get_files_xrootd.py
```

Now, you can edit the files:

* **crab_config_nano_MC.py**: edit lines 7-11.
* **nanoanalyzer_mc.py**: choose the correct global_tag, good_JSON and out_file must be the same as out_file in crab_config_nano_MC-py 
* **submit_nano_monte_carlo.py**: edit input_file with the name you choose on get_files_xrootd.py and put the number of jobs.

* Then, to call CRAB, do the following:

```
python submit_nano_monte_carlo.py
```





