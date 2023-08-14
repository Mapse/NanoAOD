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

#### Global tag
global_tag = '106X_dataRun2_v32' # 2016, 2017, 2018

#### JSON
good_JSON = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Legacy_2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt'

# 2016: '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Legacy_2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt'
# 2017: '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt'
# 2018: '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt'

#################################################################################################################################
#################################################################################################################################
##

# Reads GT
process.GlobalTag.globaltag = global_tag

##################################################################################

# Intialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('Nano')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options   = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Reads JSON
goodJSON = good_JSON

# Get the luminosity list from the JSON
# Only uncomment if you run in Data
myLumis = LumiList.LumiList(filename = goodJSON).getCMSSWString().split(',')

##################################################################################

# Load jet correction services for all jet algoritms
process.load("JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff")
#

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
#MuOnia UL 2017E AOD
#'/store/data/Run2017E/MuOnia/AOD/09Aug2019_UL2017-v1/20000/000FAD22-300F-164E-BADC-5CDDC14A4E53.root',)
# MuOnia UL 2018A AOD
#'/store/data/Run2018A/MuOnia/AOD/12Nov2019_UL2018-v1/00000/00328C0D-B5A2-1448-8895-C2B3D2FA00E0.root',)
# MuOnia UL 2018D AOD
#'/store/data/Run2018C/MuOnia/AOD/12Nov2019_UL2018-v1/00000/01359FE6-ACF8-0649-902C-82BE6BD89BD2.root',)
# Chamonium UL AOD
#          '/store/data/Run2017E/Charmonium/AOD/09Aug2019_UL2017-v1/60000/FFB59A23-0B37-D240-A528-D4C16BC7AF26.root',)
# Chamonium UL miniAOD
#          '/store/data/Run2017E/Charmonium/MINIAOD/09Aug2019_UL2017-v1/60000/6CFDC26E-8628-6446-890E-AC0F2E3A330D.root',)
#                            fileNames = cms.untracked.vstring(''))
			    ))

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
                              #outFile = cms.string('Data17DoubleMuonRunE_AOD.root'),
                              outFile = cms.string('REQUEST_NAME.root'),
                              writeCovariance = cms.bool(False),
                              # Change this:
                              # If MC:
                              #isData = cms.bool(False)
                              # If Data:
                              isData = cms.bool(True)
)
process.p = cms.Path(process.nano)
