from CRABClient.UserUtilities import config
import os
import getpass


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
#config.General.requestName = 'NanoAODplus_JpsiToMuMuwithDstar2018_15-08-2021' # For 2018
config.General.requestName = 'NanoAODplus_JpsiToMuMuwithDstar2017_15-08-2021' # For 2017
config.General.workArea = 'crab_projects_monte_carlo' # or Data
config.General.transferOutputs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'nanoanalyzer_Monte_Carlo.py'
config.JobType.numCores = 1
config.JobType.maxMemoryMB = 2500

#config.Data.outputDatasetTag = 'JpsiToMuMuwithDstar' # For 2018
config.Data.outputDatasetTag = 'CRAB_PrivateMC_NewTested_2017' # For 2017
#config.Data.userInputFiles = open('paths_monte_carlo/JpsiToMuMuwithDstar2018.txt').readlines() # For 2018
config.Data.userInputFiles = open('paths_monte_carlo/JpsiToMuMuwithDstar2017.txt').readlines() # For 2017
config.Data.inputDBS = 'phys03'
config.Data.publishDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.totalUnits = config.Data.unitsPerJob * 10000
config.Data.outLFNDirBase = '/store/user/' + getpass.getuser() + '/'
config.Data.publication = False
#config.Data.outputPrimaryDataset = 'CRAB_PrivateMC_new_2018' # For 2018
config.Data.outputPrimaryDataset = 'CRAB_PrivateMC_NewTested_2017' # For 2017

config.JobType.outputFiles = ['JpsiToMuMuwithDstar_MC_2017_NanoAODplus.root'] 
config.Site.storageSite = 'T2_BR_UERJ'
config.Data.ignoreLocality = False

#config.Site.whitelist = ['T2_BR_UERJ']


##################################################################################





