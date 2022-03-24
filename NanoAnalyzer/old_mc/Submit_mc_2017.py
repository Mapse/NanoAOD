from CRABClient.UserUtilities import config
import os
import getpass

import datetime

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

today = datetime.date.today()
today = today.strftime("%d-%m-%Y")

njobs = 10000

# Generation name
name = 'Jpsi_Dstar_new_accept_cut'

config = config()
config.General.requestName = name + '_' + today
config.General.workArea = 'crab_projects_monte_carlo' # or Data
config.General.transferOutputs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'nanoanalyzer_mc_2017.py'
config.JobType.numCores = 1
config.JobType.maxMemoryMB = 2200

config.Data.outputDatasetTag = 'Jpsi'
config.Data.userInputFiles = open('paths_monte_carlo/Jpsi.txt').readlines()
config.Data.inputDBS = 'phys03'
config.Data.publishDBS = 'phys03'
config.Data.splitting = 'FileBased' 
config.Data.unitsPerJob = 1
config.Data.totalUnits = config.Data.unitsPerJob * njobs
config.Data.outLFNDirBase = '/store/group/uerj/' + getpass.getuser() + '/'
config.Data.publication = False
config.Data.outputPrimaryDataset = 'CRAB_PrivateMC_RunII_UL_2017'

config.JobType.outputFiles = ['Jpsi_Dstar_new_accept_cuts_13TeV.root'] 
#config.Site.storageSite = 'T2_BR_UERJ'
config.Site.storageSite = 'T2_US_Caltech'
config.Data.ignoreLocality = False




##################################################################################





