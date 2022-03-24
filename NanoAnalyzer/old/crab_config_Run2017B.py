from CRABClient.UserUtilities import config
import getpass

user = getpass.getuser()

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

config.General.workArea = 'CharmoniumRun2017UL' # or Data

config.General.transferOutputs = True
config.General.transferLogs = True
config.JobType.pluginName = 'Analysis'
config.JobType.maxMemoryMB = 2000
config.Data.totalUnits = -1 
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1

config.Data.publication = False
config.Data.inputDBS = 'global'
config.Data.outLFNDirBase = '/store/user/' + user + '/Data17UL'
#config.Site.storageSite = 'T3_CH_CERNBOX'
config.Site.storageSite = 'T2_BR_UERJ'
config.Data.ignoreLocality = False

### ---------------------  AOD Primary Dataset --------------------- ###
config.General.requestName = 'CharmoniumRun2017B_AOD'
config.JobType.psetName = 'nanoanalyzer_Run2017B.py'
config.Data.inputDataset = '/Charmonium/Run2017B-09Aug2019_UL2017-v1/AOD'
config.Data.outputDatasetTag = 'CharmoniumRun2017B_AOD'
config.JobType.outputFiles = ['CharmoniumRun2017B_AOD.root']
