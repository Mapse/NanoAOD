# NanoAODPlus
This repository contains the codes for conversion from AOD or miniAOD to nanoAOD.

## Running Private SPS

If you want to run a private sample you can perform the steps LHE -> ROOT (hadronization) -> GEN-SIM -> DIGI2RAW -> HLT -> AOD on this repository https://github.com/Mapse/SPS_MC/tree/main

To run the NanoAODPlus step for private production you will use three files:

* **crab_config_nano_MC.py**
* **nanoanalyzer_mc.py**
* **submit_nano_monte_carlo.py**

First, you need to get the .txt files with the paths to the AOD files. Go to /NanoAOD/NanoAnalyzer/paths_monte_carlo/.

The first thing to do here is to identify what is the path of your files. Using T2_Caltec_US as an example, the command is,

```
xrdfs xrootd-redir.ultralight.org ls -l path_to_files
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
* **nanoanalyzer_mc.py**: choose the correct global_tal, good_JSON and out_file must be the same as out_file in crab_config_nano_MC-py 
* **submit_nano_monte_carlo.py**: edit input_file with the name you choose on get_files_xrootd.py and put the number of jobs.

* Then, to call CRAB, do the following:

* ```
python submit_nano_monte_carlo.py
  ```





