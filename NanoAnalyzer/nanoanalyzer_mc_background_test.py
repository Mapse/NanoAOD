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
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load('Configuration.StandardSequences.EndOfProcess_cff')

#################################################################################################################################
############################################################ To edit ############################################################
#################################################################################################################################

global_tag = '106X_mc2017_realistic_v8' # 2017 ???
#global_tag = '106X_upgrade2018_realistic_v15_L1v1' #2018
#global_tag = '106X_dataRun2_v32' #2016

#### JSON
good_JSON = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt'

# 2016: '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Legacy_2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt'
# 2017: '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt'
# 2018: '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt'

## Input file

#in_file_1 = 'file:///eos/user/m/mabarros/Monte_Carlo/background_mc/AOD/0001E84B-C0AE-EF49-9FE6-6D3C50912969.root'
#in_file_2 = 'file:///eos/user/m/mabarros/Monte_Carlo/background_mc/AOD/002798D1-44C1-7144-812E-B89C96EC0AAB.root'
#in_file_3 = 'file:///eos/user/m/mabarros/Monte_Carlo/background_mc/AOD/009B133F-7CCE-864E-982E-C40E1454C96E.root'
#in_file_4 = 'file:///eos/user/m/mabarros/Monte_Carlo/background_mc/AOD/00BFB6BE-97DC-044E-9C71-A7DB89D21BFE.root'
#in_file_5 = 'file:///eos/user/m/mabarros/Monte_Carlo/background_mc/AOD/00DEAEE3-4FAB-5A40-96F5-4601587ED2D4.root'
#in_file_6 = 'file:///eos/user/m/mabarros/Monte_Carlo/background_mc/AOD/02050794-E23C-7A47-91DD-27FCB5668309.root'
#in_file_7 = 'file:///eos/user/m/mabarros/Monte_Carlo/background_mc/AOD/02C7D6CD-28D8-EF40-80F8-F07687A83607.root'
#in_file_8 = 'file:///eos/user/m/mabarros/Monte_Carlo/background_mc/AOD/030B39FD-C32C-C94C-9345-FA17BF3F09AD.root'
#in_file_9 = 'file:///eos/user/m/mabarros/Monte_Carlo/background_mc/AOD/033C086A-CB3B-2048-B70D-0C6EA47A731A.root'
#in_file_10 = 'file:///eos/user/m/mabarros/Monte_Carlo/background_mc/AOD/03666741-6AD6-2643-A29B-B513FABA772D.root'
#in_file_11 = 'file:///eos/user/m/mabarros/Monte_Carlo/background_mc/AOD/041412F6-4953-3440-AE4F-58D8CBCD392C.root'
#in_file_12 = 'file:///eos/user/m/mabarros/Monte_Carlo/background_mc/AOD/043BB25E-7D0B-AE44-A71E-151FD49F07E7.root'
#in_file_13 = 'file:///eos/user/m/mabarros/Monte_Carlo/background_mc/AOD/04C4E57B-438A-E44A-83B6-B3324D52F3D0.root'


in_file_1 = 'file:///eos/user/m/mabarros/Monte_Carlo/background_mc/AOD/0553A4AA-9373-9149-9519-650CCB9CC022.root'
in_file_2 = 'file:///eos/user/m/mabarros/Monte_Carlo/background_mc/AOD/056EDEC1-6B52-7F44-95EC-DF77142BE6A9.root'
in_file_3 = 'file:///eos/user/m/mabarros/Monte_Carlo/background_mc/AOD/05F46520-4E20-6944-BE36-D7473C115753.root'
in_file_4 = 'file:///eos/user/m/mabarros/Monte_Carlo/background_mc/AOD/06804BD4-24FE-D843-B217-742ADF9F7D47.root'
in_file_5 = 'file:///eos/user/m/mabarros/Monte_Carlo/background_mc/AOD/068417B7-62C3-E744-B6CF-C19BFC68AF12.root'
in_file_6 = 'file:///eos/user/m/mabarros/Monte_Carlo/background_mc/AOD/069C3BD2-8546-7142-9D05-FE84F6DB55F4.root'
in_file_7 = 'file:///eos/user/m/mabarros/Monte_Carlo/background_mc/AOD/06C1093A-237F-274E-B9B3-3F86F54CE3B6.root'
in_file_8 = 'file:///eos/user/m/mabarros/Monte_Carlo/background_mc/AOD/06D8772F-3612-2C46-A9CD-D9BBE18B4FBE.root'
in_file_9 = 'file:///eos/user/m/mabarros/Monte_Carlo/background_mc/AOD/07272228-132A-DF44-BB4D-C0302274D529.root'
in_file_10 = 'file:///eos/user/m/mabarros/Monte_Carlo/background_mc/AOD/074DEBC6-91E1-FA4B-B196-6EE8E74CA123.root'
in_file_11 = 'file:///eos/user/m/mabarros/Monte_Carlo/background_mc/AOD/07580883-1365-0B4D-9B85-1510979FBCF7.root'
in_file_12 = 'file:///eos/user/m/mabarros/Monte_Carlo/background_mc/AOD/075C3C19-4541-1241-9F65-3724F10FD88A.root'
in_file_13 = 'file:///eos/user/m/mabarros/Monte_Carlo/background_mc/AOD/0819789E-4CF6-4E44-BBF4-D220DF78360A.root'


## Output file
path_out = '/eos/user/m/mabarros/Monte_Carlo/background_mc/BcToJPsiMuMu_inclusive_TuneCP5_13TeV'
out_file = path_out + '/' +'BcToJPsiMuMu_inclusive_TuneCP5_13TeV_2016_NanoAODPlus_2.root'

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

###################################### GLOBAL TAG ################################
# Change the global tag accordingly
# UL Data Run II 
#process.GlobalTag.globaltag = '106X_dataRun2_v28'

# MC 
process.GlobalTag.globaltag = global_tag

##################################################################################

# Intialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('Nano')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options   = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
                                        
# Set the maximum number of events to be processed here
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))

#####################################  JSON FILE #################################
# Change the directory and JSON file accordingly

goodJSON = good_JSON

# Get the luminosity list from the JSON
# Only uncomment if you run in Data
#myLumis = LumiList.LumiList(filename = goodJSON).getCMSSWString().split(',')

#################################### INPUT FILE ##################################
# To test locally or submit batch job using condor, use this:
#fileinPut = FileUtils.loadListFromFile (InputDir)

#process.source = cms.Source("PoolSource",
#                            fileNames = cms.untracked.vstring(*fileinPut)
#)

# To submit batch job using CRAB or test locally (2nd option), use this:
#process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(in_file_1, in_file_2, in_file_3, in_file_4, in_file_6, in_file_7, in_file_8, in_file_9, in_file_10))

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(in_file_2, in_file_3, in_file_4, 
                                                                            in_file_6, in_file_7, in_file_8, 
                                                                            in_file_9, in_file_10, in_file_11,
                                                                            in_file_12, in_file_13))
##################################################################################

# Process the lumi
# Only uncomment if you run in Data
#process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
#process.source.lumisToProcess.extend(myLumis)

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
                              #outFile = cms.string('JpsiToMuMuwithDzero_MC_2017_NanoAODplus.root'),
                              outFile = cms.string(out_file),
                              #outFile = cms.string('JpsiToMuMuwithDstar_MC_2017_NanoAODplus.root'),
                              #outFile = cms.string('JpsiToMuMuwithDSPlus_MC_2017_NanoAODplus.root'),
                              #outFile = cms.string('JpsiToMuMuwithLambdac_MC_2017_NanoAODplus.root'),
                              #outFile = cms.string('Data18CharmoniumRunB_AOD.root'),
                              #outFile = cms.string('Data18CharmoniumRunE_AOD.root'),
                              #outFile = cms.string('Data17MuOniaRunE_MINIAOD.root'),
                              #outFile = cms.string('Data16_MinB_foroptmu.root'),
                              #outFile = cms.string('Data16_DMu_mDD.root'),
                              #outFile = cms.string('Data16_DMu_testnewGeo.root'),
                              #outFile = cms.string('Data16_0B_testmitRunC0BAOD.root'),
                              # If data not in CRAB:
                              #outFile = cms.string('JpsiToMuMuwithDZero_13TeV_NanoAOD.root'),
                              #outFile = cms.string('JpsiToMuMuwithDZero_13TeV_NanoAOD_2018.root'),
                              #outFile = cms.string('JpsiToMuMuwithDZero_13TeV_NanoAOD_Data.root'),
                              #outFile = cms.string('MC15_Dstar_foropt.root'),
                              #outFile = cms.string('MC15_MinBnPU_foropt.root'),
                              #outFile = cms.string('MC15_MinBPU_foropt.root'),
                              # make sure the name is same as the crab config
                              # Change this:
                              # If MC:
                              isData = cms.bool(False)
                              # If Data:
                              #isData = cms.bool(True)
)

process.p = cms.Path(process.nano)
