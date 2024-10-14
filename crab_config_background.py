from CRABClient.UserUtilities import config
import getpass

######################### To edit #################################

grid_path = '/store/group/uerj/' # '/store/user/' + user + '/OUT_DIR'
storage_site = 'T2_US_Caltech' # 'T2_BR_UERJ'

######################## End editing ###############################

user = getpass.getuser()

config = config()

config.General.workArea = 'WORKAREA' # or Data

###### Temporary
#config.Site.ignoreGlobalBlacklist = True

config.General.transferOutputs = True
config.JobType.pluginName = 'Analysis'
config.JobType.maxMemoryMB = 2000
config.Data.totalUnits = -1 
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1

config.Data.publication = True
config.Data.inputDBS = 'global'
config.Data.outLFNDirBase = grid_path + getpass.getuser() + '/'
config.Site.storageSite = storage_site
#config.Data.ignoreLocality = False
#config.Site.whitelist = ['T2_US_Caltech']

### ---------------------  AOD Primary Dataset --------------------- ###
config.General.requestName = 'REQUEST_NAME'
config.JobType.psetName = 'CONFIG_NAME.py'
config.Data.inputDataset = 'INPUT_DATASET'
config.Data.outputDatasetTag = 'REQUEST_NAME'
config.JobType.outputFiles = ['REQUEST_NAME.root']
