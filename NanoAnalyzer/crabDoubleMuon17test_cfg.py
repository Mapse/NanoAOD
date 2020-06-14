from CRABClient.UserUtilities import config

##################################################################################
# CRAB configuration file for multiprocessing                                    #
# This configuration is used for Data 2015                                       #
# Instead of arrange by section, we arrange by common and specific               #
# Plus, we add more parameters like max memory, ignorelocality etc               #
# Use to submit job using CRAB                                                   #
# More details about CRAB configuration see:                                     #
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile          #
# Uncomment & comment relevant lines before you use it                           #
##################################################################################

##################################################################################
################################## CRAB COMMON ###################################
##################################################################################
#config.General.workArea = '/home/cms-opendata/CMSSW_4_2_8/src/Demo/DemoAnalyzer/'
# or
config = config()

config.General.workArea = 'Data17ULtest' # or Data

config.General.transferOutputs = True
config.General.transferLogs = True
config.JobType.pluginName = 'Analysis'
config.JobType.maxMemoryMB = 2500
config.Data.totalUnits = 1000000
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1

config.Data.publication = False
config.Data.inputDBS = 'global'
config.Data.outLFNDirBase = '/store/user/mabarros/Data17DoubleMuonSet'
config.Site.storageSite = 'T2_BR_UERJ'
config.Data.ignoreLocality = False

##################################################################################


##################################################################################
################################ CRAB SPECIFIC ###################################
##################################################################################
# For each datasets you want to run

### --------------- Data17DoubleMuonRunE AOD --------------------- ###
config.General.requestName = 'Data17DoubleMuonRunE_AOD'
config.JobType.psetName = 'nanoanalyzer_cfg.py'
config.Data.inputDataset = '/DoubleMuon/Run2017E-09Aug2019_UL2017-v1/AOD'
config.Data.outputDatasetTag = 'Data17DoubleMuonRunE_AOD'
config.JobType.outputFiles = ['Data17DoubleMuonRunE_AOD.root']
#config.Site.whitelist = ['T2_IT_Legnaro']

""" ### --------------- Data17MuOniaRunE MINIAOD --------------------- ###
config.General.requestName = 'Data17DoubleMuon_MINIAOD'
config.JobType.psetName = 'nanoanalyzer_cfg.py'
config.Data.inputDataset = '/MuOnia/Run2017E-09Aug2019_UL2017-v1/MINIAOD'
config.Data.outputDatasetTag = 'Data17MuOniaRunE_MINIAOD'
config.JobType.outputFiles = ['Data17MuOniaRunE_MINIAOD.root'] """