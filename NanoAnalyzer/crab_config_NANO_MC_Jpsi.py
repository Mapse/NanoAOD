from CRABClient.UserUtilities import config
import getpass

######################### To edit #################################
out_pri_dataset = 'CRAB_PrivateMC_RunII_UL_2017'
pset_name = 'nanoanalyzer_mc.py'
grid_path = '/store/group/uerj/'
out_file = 'Jpsi_Dstar_new_accept_cuts_13TeV.root'
storage_site = 'T2_US_Caltech' # 'T2_BR_UERJ'
###################################################################



## Config
config = config()
config.General.requestName = 'NANO_MC_Jpsi_24-03-2022'
config.General.workArea = 'crab_projects_monte_carlo'
config.General.transferOutputs = True

config.JobType.pluginName = 'Analysis' 
config.JobType.psetName = pset_name
config.JobType.numCores = 1
config.JobType.maxMemoryMB = 2000 #### I changed from 2500!!!!!


config.Data.outputDatasetTag = 'Jpsi'
config.Data.userInputFiles = open('paths_monte_carlo/Jpsi.txt').readlines()
config.Data.inputDBS = 'phys03'
config.Data.publishDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.totalUnits = config.Data.unitsPerJob * 1
config.Data.outLFNDirBase = grid_path + getpass.getuser() + '/'
config.Data.publication = True
config.Data.outputPrimaryDataset = out_pri_dataset  # Comes after /store/user/mabarros  

config.JobType.outputFiles = [out_file] 

config.Site.storageSite = storage_site

