##################################################################################
# Nanoanalyzer configuration file for all years                                  #
# Use for HT Condor and VM                                                       #
# Uncomment & comment relevant lines before you run it                           #
##################################################################################

import sys
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.Types as CfgTypes
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.Utilities.FileUtils as FileUtils
from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
from FWCore.ParameterSet.VarParsing import VarParsing

process = cms.Process("Nano")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
#process.load("Configuration.Geometry.GeometryIdeal_cff") # for 2011
#process.load("Configuration.StandardSequences.Geometry_cff") # for 2010
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load('Configuration.StandardSequences.EndOfProcess_cff')

################################### VARPASSING ###################################
# Use VarParsing to specify your input directory and output file
# Comment Varparsing if you want to submit job using CRAB
# because CRAB does not support input directory as we put input dataset directly
""" options = VarParsing('analysis')
options.register(
    "/afs/cern.ch/work/m/mabarros/public/CMSSW_10_2_14/src/NanoAOD/NanoAnalyzer",
    "",
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "/afs/cern.ch/work/m/mabarros/public/CMSSW_10_2_14/src/NanoAOD/NanoAnalyzer"
    )
options.register(
    "test.root",
    "",
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "teste1.root"
    )
options.parseArguments()

if (options.inputDir == ""):
    sys.exit("Directory to find input file where???")
else:
    # 2020 Mapse Test
    InputDir = "/afs/cern.ch/work/m/mabarros/public/CMSSW_10_2_14/src/NanoAOD/NanoAnalyzer" """
    # 2010 VM
    #InputDir = "/home/cms-opendata/CMSSW_4_2_8/NanoAOD/NanoAnalyzer/" + options.inputDir
    # 2011 NAF/HTC
    # InputDir = "/nfs/dust/cms/user/zulaiha/PhD/CMSSW_5_3_32/src/NanoAOD/NanoAnalyzer/" + options.inputDir
    # 2015 NAF/HTC
    #InputDir = "/nfs/dust/cms/user/zulaiha/PhD/CMSSW_7_6_1/src/NanoAOD/NanoAnalyzer/" + options.inputDir

##################################################################################

###################################### GLOBAL TAG ################################
# Change the global tag accordingly
# 2010 Data_MinB
#process.GlobalTag.globaltag = 'FT_R_42_V10A::All'
# 2010 Data_Mu
#process.GlobalTag.globaltag = 'FT_R_42_V10A::All'

# 2011 Data_MinB
#process.GlobalTag.globaltag = 'FT_53_LV5_AN1::All'
# 2011 Data_DoubleMu
#process.GlobalTag.globaltag = 'FT_53_LV5_AN1::All'
# 2011 MC_QCDpt0to5, QCDpt5to15
#process.GlobalTag.globaltag = 'START53_LV6::All'
# 2015 Data MinB
#process.GlobalTag.globaltag = '75X_dataRun2_Prompt_ppAt5TeV_v0'
# 2015 MC_DstarD0pi
#process.GlobalTag.globaltag = '76X_mcRun2_asymptotic_v12'
# 2015 MC_MinB
#process.GlobalTag.globaltag = '76X_mcRun2_asymptotic_v12'
#process.GlobalTag.globaltag = 'MCRUN2_74_V8'
# 2015 MC_private
#process.GlobalTag.globaltag = '75X_mcRun2_asymptotic_ppAt5TeV_v3'

# 2016 Data_0B
#process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v3'
# 2016 Double Muon 
#process.GlobalTag.globaltag = '80X_dataRun2_2016LegacyRepro_v4'
# 2018 Data
#process.GlobalTag.globaltag = '102X_dataRun2_Prompt_v4'
#process.GlobalTag.globaltag = '102X_dataRun2_Prompt_v11'
#process.GlobalTag.globaltag = '102X_dataRun2_Prompt_v12'

# This is the global tag used for reprocessing of the 2017 datasets. 
# https://twiki.cern.ch/twiki/bin/view/CMS/PdmV2017Analysis#Release

######################################## Global Tag ########################################

# The following GT are used for legacy samples. You can see which GT can be used on pdmV twiki.

process.GlobalTag.globaltag = '106X_dataRun2_v28'

##################################################################################

# Intialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('Nano')
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.options   = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
                                        
# Set the maximum number of events to be processed here
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))

#####################################  JSON FILE #################################
# Change the directory and JSON file accordingly
# Only uncomment if you run in Data
# 2010 JSON
#goodJSON = '/home/cms-opendata/CMSSW_4_2_8/src/NanoAOD/NanoAnalyzer/datasets_2010/Data/Cert_136033-149442_7TeV_Apr21ReReco_Collisions10_JSON_v2.txt'
# 2011 JSON
#goodJSON = '/nfs/dust/cms/user/zulaiha/PhD/CMSSW_5_3_32/src/NanoAOD/NanoAnalyzer/datasets_2011/Data/Cert_160404-180252_7TeV_ReRecoNov08_Collisions11_JSON.txt'
# 2012 JSON
#goodJSON = '/nfs/dust/cms/user/zulaiha/opendata/CMSSW_5_3_32/src/Demo/DemoAnalyzer/datasets_2012/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt'
# 2015 JSON
#goodJSON = '/nfs/dust/cms/user/zulaiha/PhD/CMSSW_7_6_1/src/NanoAOD/NanoAnalyzer/datasets_2015/Cert_262081-262273_5TeV_PromptReco_Collisions15_25ns_JSON_v2.txt'
#goodJSON = '/nfs/dust/cms/user/zulaiha/PhD/CMSSW_7_6_1/src/NanoAOD/NanoAnalyzer/datasets_2015/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_v2.txt'
# 2016 JSON
#goodJSON = '/nfs/dust/cms/user/zulaiha/PhD/CMSSW_8_0_13/src/NanoAOD/NanoAnalyzer/datasets_2016/Data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
# should update for legacy!

# 2018 JSON
#goodJSON = '/nfs/dust/cms/user/zulaiha/PhD/CMSSW_10_2_6/src/NanoAOD/NanoAnalyzer/datasets_2018/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt'
goodJSON = '/afs/cern.ch/work/m/mabarros/public/CMSSW_10_6_12/src/NanoAOD/NanoAnalyzer/jsonfiles/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt'

##################################################################################

# Get the luminosity list from the JSON
# Only uncomment if you run in Data
myLumis = LumiList.LumiList(filename = goodJSON).getCMSSWString().split(',')

#################################### INPUT FILE ##################################
# To test locally or submit batch job using condor, use this:
#fileinPut = FileUtils.loadListFromFile (InputDir)

#process.source = cms.Source("PoolSource",
#                            fileNames = cms.untracked.vstring(*fileinPut)
#)

# To submit batch job using CRAB or test locally (2nd option), use this:
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('')
)
#/store/data/Run2017E/DoubleMuon/AOD/09Aug2019_UL2017-v1/50000/FFCDF27E-DBEA-AD47-85AE-4F214D13EB62.root

##################################################################################

# Process the lumi
# Only uncomment if you run in Data
process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
process.source.lumisToProcess.extend(myLumis)

# *************************************************
# number of events to be skipped (0 by default)   *
# *************************************************
process.source.skipEvents = cms.untracked.uint32(0)

# Process the analyzer
process.nano = cms.EDAnalyzer('NanoAnalyzer',
                              # Change this:
                              # If HTC/VM:
                              #outFile = cms.string(options.outputName),
                              # If interactive:
                              #outFile = cms.string('test.root'), 
                              # If CRAB:
                              #outFile = cms.string('Data15_MinB_mutrig.root'),
                              #outFile = cms.string('Data16_MinB_foroptmu.root'),
                              outFile = cms.string('Data17DoubleMuonRunEAOD.root'),
                              #outFile = cms.string('Data16_DMu_mDD.root'),
                              #outFile = cms.string('Data16_DMu_testnewGeo.root'),
                              #outFile = cms.string('Data16_0B_testmitRunC0BAOD.root'),
                             
                              #outFile = cms.string('MC15_Dstar_foropt.root'),
                              #outFile = cms.string('MC15_MinBnPU_foropt.root'),
                              #outFile = cms.string('MC15_MinBPU_foropt.root'),
                              # make sure the name is same as the crab config

                              # Change this:
                              # If MC:
                              #isData = cms.bool(False)
                              # If Data:
                              isData = cms.bool(True)
)
process.p = cms.Path(process.nano)
