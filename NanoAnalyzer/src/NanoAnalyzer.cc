// **********************************************************************//
// Prototype implementation of NanoAnalyzer                              //
// for the creation of Run 1 nanoAOD-like ntuple from AOD or RECO input  //
// and corresp. Run 2 reference/validation ntuples from AOD or miniAOD   //
// **********************************************************************//
//
// Choose appropriate CMSSW flag (only one of the three should be set):
// keyword of H->4lepton electron implementation: // H4lepton

// Main application is Run 1 legacy CMSSW53X (e.g. 2011/2, 5_3_32)
//#define CMSSW53X
// Define here instead in case CMSSW version is CMSSW42X (e.g. 2010, 4_2_8)
//#define CMSSW42X
// Define here instead in case CMSSW version is for Run 2, CMSSW7 or higher
// (e.g. 2015-18), for validation purposes
//#define CMSSW7plus

// afiqaize try to automate the above; CMSSW is tied to particular compiler versions
// 42X taken from 4_2_8, 53X from 5_3_32, 7XX from 7_6_1
#define GCC_VERSION ( 10000 * __GNUC__ + 100 * __GNUC_MINOR__ + __GNUC_PATCHLEVEL__ )
#if GCC_VERSION < 40305
#define CMSSW42X
#elif GCC_VERSION < 40703
#define CMSSW53X
#elif GCC_VERSION > 40902
#define CMSSW7plus
#endif

// activation flag to check and store compatibility with Golden JSON 
// default is on for 2010 data (not MC, automatic), off otherwise
#ifdef CMSSW42X
#define JSONcheck
#endif

// the following are flags which can be activated manually

// activate this to check Golden JSON setting, i.e. abort upon nonJSON event
// activate this only for check when Golden JSON is activated in configuration
// (protection against accidental use of wrong or no JSON in configuration)
//#define JSONcheckabort

// activate this only for data sets for which trigger treatment is already 
// implemented; checks and aborts in case of inconsistency
// should be activated by default if trigger is implemented for dataset
// (protection against inconsistencies in NanoTrigger implementation)
//#define trigcheckabort

// activate this to achieve maximum compatibility with official Run 2 nanoAOD
// (default is off -> better performance, e.g. use of beam spot constraint)
// (mainly for Run 2 validation - not recommended for Run 1)
// only relevant for Run 2 AOD, turn on for better consistency with miniAOD 
//                        (worse performance)
// nanoext flag can also be steered/overruled via configuration file
//#define Compatibility

// activate this when you read from miniAOD (Run 2 only!)
//#define miniAOD

// turn this on only when using old 'Demo' version nanoAnalyzer setup on VM
//#define Demo
// turn this on to activate D meson true info printout on MC
//#define Bingo

// 18.2.2019,  Nur Zulaiha Jomhari, Achim Geiser  

// The following comment lines are out of date and need update
//
// version 0.0:
//
// run, event, luminosityBlock fully implemented
//
// almost all Muon nanoAOD variables implemented, 
// most cross-checked against nanoAODv2
//    + extension to looser muon selection and more variables
//    "Muon_isNano" flag selects official nanoAOD muons only 
//
// most PV nanoAOD variables implemented
// and cross-checked againts nanoAODv2 (w/o BS option)
//
// OtherPV variables implemented but content not yet validated
//    + extension to list of all primary vertices with extra variables
//     
// partial implementation of variables for electron, photon, tau, MET, JET
// electron variables partially validated, others not yet
//
// implementation of GenPart variables, but not yet validated
//    (so far c,b, muons and D mesons only)
//
//    extension to D0 mesons (in SV format) and D* mesons
//    extension to trigger information in nonstandard format
//    
// code for different CMSSW versions is steered by #ifdef flags
// some external non-CMSSW code is also still being used 
//    -> will be removed or integrated 

// keyword: nuha
// system include files
#include <memory>
#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>
#include <stdexcept>
#include <cmath>
#include <string>

#ifdef CMSSW42X
#include <boost/unordered_map.hpp>
using boost::unordered_map;
#else
#include <unordered_map>
using std::unordered_map;
#endif

// general user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// for fillDescriptions
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

// ------ EXTRA HEADER FILES--------------------//
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
//#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Ref.h"
// the following does not seem to be needed (included through dEdx), but also not to hurt
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Provenance/interface/EventAuxiliary.h"
// Afiq, to get release version
#include "FWCore/Version/interface/GetReleaseVersion.h"
// nuha: header for conversion tools
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
// effective area for rho
//https://github.com/cms-sw/cmssw/blob/CMSSW_9_4_X/RecoEgamma/EgammaTools/interface/EffectiveAreas.h
//#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"

// for math and file handling
#include "TMath.h"
#include "TFile.h"
#include "Math/VectorUtil.h"

// for histogramming
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TStyle.h"
#include "TLorentzVector.h"
#include "TTree.h"

// for trigger information
#include "FWCore/Common/interface/TriggerNames.h"
// not needed?
//#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "DataFormats/Common/interface/TriggerResults.h"
// Afiq
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"  //Qun 17-10-19 TriggerObj
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h" //Qun 17-10-19 TriggerObj
#include <cassert>  //Qun

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h" //Qun
#include "HLTrigger/HLTcore/interface/HLTEventAnalyzerAOD.h" //Qun
#include "HLTrigger/HLTcore/interface/TriggerSummaryAnalyzerAOD.h"//Qun

// for tracking information
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
// for dEdx
#include "DataFormats/TrackReco/interface/DeDxData.h" 	

#ifdef miniAOD
// tracks from PATCandidates
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
// there are three track collections in miniAOD (in addition to the muon 
// and electron collections), embedded into particleflow objects:
//   *** update this ***
// packedPFCandidates allows to rebuild tracks (for pt>0.95/0.4 GeV) using
//                    pseudotrack() (object) or besttrack() (pointer)
// PackedPFCandidatesDiscarded presumably contains only discarded duplicate 
//                    muon candidates -> do not use
// lostTracks (high purity only) might contain many of the non-vertex tracks 
//                    needed for the slow pion measurement 
#endif



/// for track parametrization 
//#include "TrackingTools/TrajectoryParametrization/interface/CartesianTrajectoryError.h" 
#include "TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"
//#include "DataFormats/GeometryVector/interface/GlobalVector.h"
//#include "MagneticField/Engine/interface/MagneticField.h"
//#include "FWCore/Utilities/interface/Likely.h"
//#include "DataFormats/Math/interface/AlgebraicROOTObjects.h"
//#include "DataFormats/GeometrySurface/interface/Surface.h"
//#include "TrackingTools/TrajectoryParametrization/interface/LocalTrajectoryParameters.h"
#include "TrackingTools/TrajectoryParametrization/interface/PerigeeTrajectoryParameters.h"
// commented for 10_X
//#include "RecoVZero/VZeroFinding/interface/VZeroFinder.h"

// for vertex information 
// reconstructed primary is typically within 0.02 cm of true primary in z
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include <DataFormats/VertexReco/interface/Vertex.h>

// for vertices refit:
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/AdaptiveVertexFinder/interface/AdaptiveVertexReconstructor.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/TrimmedKalmanVertexFinder/interface/KalmanTrimmedVertexFinder.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/VertexPrimitives/interface/VertexException.h"
#include "TrackingTools/TransientTrack/interface/TrackTransientTrack.h"
//#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexProducerAlgorithm.h"

// for beamspot information (beam pipe radius: 5.8 -> 4.3 cm, 
//                           beam spot at x~0.2, y~0.4, z~0.3 cm,
//                           width x~0.002?, y~0.002?, z~5.5 cm,
// beam-spot constrained vertices: x~0.001 , y~0.001 , z~0.004 cm) 
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

// for error matrix handling
#ifdef Demo  
#include "Demo/DemoAnalyzer/src/inverceM.h"
#endif
#ifndef Demo  
#include "NanoAOD/NanoAnalyzer/src/inverceM.h"
#endif
#include "TMatrixT.h"
#include "Math/SMatrix.h"
#include "Math/StaticCheck.h"
//#include <MatrixRepresentationsStatic.h>
#include "Math/Expression.h"

/* #ifdef Demo
// for 'ZEUS-type' revertexing (needs external code/library)
// #include "Demo/DemoAnalyzer/interface/vxlite/VertexFitter.cc"
#include "Demo/DemoAnalyzer/interface/vxlite/LinearizedTrack.hh"
#include "Demo/DemoAnalyzer/interface/vxlite/VertexFitter.hh"
// #include "Demo/DemoAnalyzer/interface/vxlite/DAFVertexFinder.cc"
#include "Demo/DemoAnalyzer/interface/vxlite/DAFVertexFinder.hh"
#include "Demo/DemoAnalyzer/include/HelixTransClass.h"
#endif */
//#ifndef Demo
// for 'ZEUS-type' revertexing (needs external code/library)
// #include "Demo/DemoAnalyzer/interface/vxlite/VertexFitter.cc"
#include "NanoAOD/NanoAnalyzer/interface/vxlite/LinearizedTrack.hh"
#include "NanoAOD/NanoAnalyzer/interface/vxlite/VertexFitter.hh"
// #include "Demo/DemoAnalyzer/interface/vxlite/DAFVertexFinder.cc"
#include "NanoAOD/NanoAnalyzer/interface/vxlite/DAFVertexFinder.hh"
#include "NanoAOD/NanoAnalyzer/include/HelixTransClass.h"
//  #endif

#ifndef miniAOD
// for muon information
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/MuonIsolation.h"
// CHANGE, this needs to be protected by IFDEF for CMSSW 4.2.8
// or maybe it is not needed? probably it is (commented code)!
//#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"
#endif
#ifdef miniAOD
// for PAT/miniAOD
#include "DataFormats/PatCandidates/interface/Muon.h"
#endif

// for electron informaton
#ifndef miniAOD
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#endif
#ifdef miniAOD
#include "DataFormats/PatCandidates/interface/Electron.h"
#endif

// for photon information
#ifndef miniAOD
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#endif
#ifdef miniAOD
#include "DataFormats/PatCandidates/interface/Photon.h"
#endif

// for tau information
#ifndef miniAOD
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#endif
#ifdef miniAOD
#include "DataFormats/PatCandidates/interface/Tau.h"
#endif

// for MET informaton
#ifndef miniAOD
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETFwd.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#endif
#ifdef miniAOD
//#include "DataFormats/PatCandidates/interface/PFMET.h"
//#include "DataFormats/PatCandidates/interface/CaloMET.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#endif 

// for jet information
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Math/interface/deltaR.h"
#ifndef miniAOD
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/TrackJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/JetExtendedAssociation.h"
#include "DataFormats/JetReco/interface/JetID.h"
#endif
#ifdef miniAOD
#include "DataFormats/PatCandidates/interface/Jet.h"
#endif

// for gen particle information
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

// for particle flow information
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

// using namespace edm;

// class member declaration

// Qun
//#pragma link C++ class std::vector<float> +; 
//#pragma link C++ class std::vector<double> +; 

class NanoAnalyzer : public edm::EDAnalyzer
{
public:
  explicit NanoAnalyzer(const edm::ParameterSet&);
  ~NanoAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);

//      virtual void beginRun(edm::Run const&, edm::EventSetup const&); //Qun

  // this is the place to define global variables parameters 

        // declare global trigger variables

        // The dataset flags indicate the current knowledge about the whole 
        // current dataset. This should become unique after a few events.
        // The `Trig' flags indicate which datasets the current event may 
        // belong to (since an event can have fired triggers from several 
        // datasets, this is often not unique).
        // The other variables indicate thresholds for certain classes of 
        // triggers, which may span several datasets.

        // bit for generic "Good" JSON 
        bool GoodLumisection;

        // bits for generic "Good" triggers 
        bool GoodMinBiasTrigger;
        bool GoodJetTrigger;
        bool GoodMuTrigger;
        bool GoodETrigger;

        // bits for datasets
        bool MCdataset;
        bool ZeroBiasdataset;
        bool MinimumBiasdataset;
        bool Commissioningdataset;
        bool Mudataset;
        bool MuHaddataset;
        bool DoubleMudataset;
        bool MuOniadataset;
        bool Charmoniumdataset;
        bool MuMonitordataset;
        bool Jetdataset;
        bool MultiJetdataset;
        bool JetMETTauMonitordataset;
        bool BTaudataset;
        bool BParkingdataset;
        bool MuEGdataset;
        bool Electrondataset;
        bool DoubleElectrondataset;
        bool Photondataset;
        bool EGMonitordataset;
        bool METFwddataset;
        bool datasetisunique;
        std::string dataset;

        // bits and variables for MinimumBias and Commissioning datasets
        bool ZeroBiasTrig; 
        bool MinimumBiasTrig;
        bool CommissioningTrig;
        bool GoodMinimumBiasTrig;
        bool GoodMuMinimumBiasTrig;
        int MinBiasFlag;
        int MinBiasMult;
        int ZeroBiasFlag;

        // bits and variables for various Muon, Electron and Photon datasets
        bool MuTrig;
        bool MuHadTrig;
        bool DoubleMuTrig;
        bool MuEGTrig;
        bool ElectronTrig;
        bool DoubleElectronTrig;
        bool PhotonTrig;
        bool MuOniaTrig;
        bool CharmoniumTrig;
        bool MuMonitorTrig;
        bool EGMonitorTrig;
        bool GoodMuTrig;
        bool GoodETrig;
        int MuThresh; 
        int MuL1Thresh; 
        int MuL2Thresh; 
        int IsoMuThresh; 
        int DoubleMuThresh; 
        int JpsiThresh;
        int MuHadFlag;
        int MuEGFlag;
        int ElectronThresh;
        int DoubleElectronThresh;
        int PhotonThresh;

        // bits and variables for various Jet and MET datasets 
        bool JetTrig;
        bool MultiJetTrig;
        bool JetMETTauMonitorTrig;
        bool BTauTrig;
        bool BParkingTrig;
        bool METFwdTrig;
        bool GoodJetTrig;
        int JetThresh; 
        int DiJetThresh;
        int TriJetThresh;
        int QuadJetThresh;
        int HTThresh;
        int BThresh;
      
        int METThresh; 

private:

  virtual void beginJob();

  virtual void beginRun(const edm::Run &iRun, const edm::EventSetup &iStp);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void analyzeTrigger(const edm::Event&, const edm::EventSetup&, const std::string& triggerName); //Qun
  virtual void endRun(edm::Run const&, edm::EventSetup const&);  //Qun

  virtual void endJob();
  
  // declare JSON quality method (code included via #include)
  bool providesGoodLumisection(const edm::Event& iEvent);

  // declare trigger methods (code included via #include)
  bool providesGoodMinBiasTrigger(const edm::Event& iEvent);
  bool providesGoodMuTrigger(const edm::Event& iEvent);
  bool providesGoodETrigger(const edm::Event& iEvent);
  bool providesGoodJetTrigger(const edm::Event& iEvent);
  HLTConfigProvider hltConfig_;  //Qun 17-10-19 TriggerObj

  // afiqaize string manipulation methods
  static bool is_non_digit(const char &c) { return !std::isdigit(c); };
  static std::string remove_version(const std::string &path) { return path.substr(0, path.rfind("_v")); };

  // as it says on the tin
  void update_HLT_branch();

  // HLT config for reading the table and its associated process name
  //HLTConfigProvider hlt_cfg;   // superseded above
  std::string hlt_proc;
  unordered_map<std::string, uint8_t> hlt_bit;
  
#ifdef CMSSW7plus
  // for Run 2, enable access via getByToken
#ifndef miniAOD
  // AOD collections
  EDGetTokenT<reco::BeamSpot> beamTkn;
  EDGetTokenT<reco::TrackCollection> trkTkn;
  EDGetTokenT<reco::TrackCollection> gmuTkn;
  EDGetTokenT<reco::VertexCollection> primvtxTkn;
  EDGetTokenT<reco::MuonCollection> muTkn;
  EDGetTokenT<reco::GsfElectronCollection> eTkn;
  EDGetTokenT<reco::PhotonCollection> photTkn;
  EDGetTokenT<reco::PFTauCollection> tauTkn;
  EDGetTokenT<reco::GenParticleCollection> genTkn;
  // why is the next special?
  //EDGetTokenT< std::vector<reco::PFMET> > pfmetTkn;
  EDGetTokenT<reco::PFMETCollection> pfmetTkn;
  // EDGetTokenT< std::vector<reco::CaloMET> > calometTkn;
  EDGetTokenT<reco::CaloMETCollection> calometTkn;
  // not yet needed
  // EDGetTokenT< std::vector<reco::CaloMET> > muCorrmetTkn;
  EDGetTokenT<reco::PFJetCollection> pfjetTkn;
  EDGetTokenT<edm::TriggerResults> trigTkn;
  // for Trigger and Flags (Afiq) *** need to sort out duplication ***
  EDGetTokenT<edm::TriggerResults> trig_tkn;
  EDGetTokenT<edm::TriggerResults> custom_tkn;
  EDGetTokenT<trigger::TriggerEvent> trigEvn; //Qun 
  // get dEdx ValueMaps (from example by M. Soares)
  edm::EDGetTokenT<edm::ValueMap<reco::DeDxData>> dedxMapStripTag_;
  edm::EDGetTokenT<edm::ValueMap<reco::DeDxData>> dedxMapPixelTag_;
  // nuha: conversion collection
  EDGetTokenT<reco::ConversionCollection> hConvTkn;
#endif
#ifdef miniAOD
  // miniAOD collections
  EDGetTokenT<reco::BeamSpot> beamTkn;
  EDGetTokenT<pat::PackedCandidateCollection> trkTkn;
  EDGetTokenT<pat::PackedCandidateCollection> trkTkndisc;
  EDGetTokenT<pat::PackedCandidateCollection> trkTknlost;
  // EDGetTokenT<reco::TrackCollection> gmuTkn;
  EDGetTokenT<reco::VertexCollection> primvtxTkn;
  EDGetTokenT<pat::MuonCollection> muTkn;
  EDGetTokenT<pat::ElectronCollection> eTkn;
  EDGetTokenT<pat::PhotonCollection> photTkn;
  EDGetTokenT<pat::TauCollection> tauTkn;
  EDGetTokenT<reco::GenParticleCollection> genTkn;
  EDGetTokenT<pat::METCollection> pfmetTkn;
  EDGetTokenT<pat::METCollection> calometTkn;
  EDGetTokenT<pat::JetCollection> pfjetTkn;
  EDGetTokenT<edm::TriggerResults> trigTkn;
  EDGetTokenT<trigger::TriggerEvent> trigEvn; //Qun 
#endif
#endif

  std::string   processName_;  //Qun
  std::string   triggerName_;  //Qun
  // already occurs below
  edm::Handle<edm::TriggerResults> triggerResultsHandle_;
  // move below?
  edm::Handle<trigger::TriggerEvent> triggerEventHandle_;  //Qun 17-10-19 TriggerObj

  // member data
  std::string outFile;
  bool isData;
  int CMSSW;
  bool nanoext; // write out nanoAOD extensions
  bool covout; // whether to write or not covariance matrices in nanoAOD extensions
  //bool relative_;
  
  // H4lepton
  int misshits;
  double IP3d_e;
  double ErrIP3d_e;
  double SIP3d_e;  
  double relPFIso_e;

//////////////////////////////////////////////////////////////////////////////
////////////////////////// declare tree, file, hist //////////////////////////
//////////////////////////////////////////////////////////////////////////////
  
  TFile *file;
  TTree *t_event;

  /// define histograms ///
  TH1D *h_trackpt;
  TH1D *h_trackptlow;
  TH1D *h_tracketa;
  
  TH1D *h_d0pt;
  TH1D *h_dstarpt;
  TH1D *h_PS3pt;
  TH1D *h_K1pt;
  TH1D *h_P2pt;
  TH1D *h_K1eta;
  TH1D *h_P2eta;
  TH1D *h_PS3eta;
  TH1D *h_D0masscut;
  TH1D *h_deltaMassD0Dstar;
  TH1D *h_D0masscut_rightDLcut;
  TH1D *h_deltaMassD0Dstar_rightDLcut;

  // H4lepton
  // histograms for some electron vars before pt cut
  // to compare with H-4l example
  TH1D *h_p_e;
  TH1D *h_et_e;
  TH1D *h_pt_e_b4;
  TH1D *h_eta_e_b4;
  TH1D *h_phi_e;
  TH1D *h_sc_eta;
  TH1D *h_sc_rawE;
  TH1D *h_relPFIso_e;
  TH2D *h_relPFIso_pt_e;
  TH1D *h_dxy_e;
  TH1D *h_SIP3d_e_b4;
  TH1D *h_misshite;  
  
//////////////////////////////////////////////////////////////////////////////
///////////////// declare variables you want to put into the tree ////////////
//////////////////////////////////////////////////////////////////////////////

////////////////////////////////// for general ///////////////////////////////

  /// original nanoAOD ///
  UInt_t run;
  ULong64_t event;
  UInt_t luminosityBlock;

  /// nanoAOD extension ///

  // The following is a historical intention and not (yet) used
  Int_t dataIsInclusive;
  Int_t dataIsHighPtJet;
  Int_t dataIsMuonLow;
  Int_t dataIsMuonHigh;
  Int_t dataIsDimuonLow;
  Int_t dataIsDimuonHigh;
  Int_t dataIsElectronLow;
  Int_t dataIsElectronHigh;
  Int_t dataIsDiElectronLow;
  Int_t dataIsDiElectronHigh;
  Int_t dataIsMuELow;  
  Int_t dataIsMuEHigh;
  // useful dataset chains are
  //   inclusive chain:            ZeroBias -> MinimumBias -> HT -> JetMon 
  //                            -> Jet -> JetHT -> MultiJet -> HTMHTParked
  //   high pt jet chain:          Jet -> JetHT -> MultiJet
  //   low pt single muon chain:   ZeroBias -> MinimumBias -> MuMonitor -> Mu 
  //                            -> SingleMu -> MuHad -> MuEG
  //   high pt single muon chain:  Mu -> SingleMu -> MuHad -> MuEG
  //   low pt double muon chain:   ZeroBias -> MinimumBias -> MuMonitor -> Mu 
  //                            -> SingleMu -> DoubleMu(Parked) -> MuHad 
  //                            -> MuOnia(Parked) -> Scouting 
  //   high pt double muon chain:  Mu -> SingleMu -> DoubleMu(Parked) -> MuHad
  //   low pt single electron chain: ZeroBias -> MinimumBias -> EGmonitor 
  //                            -> Electron -> SingleElectron -> ElectronHad 
  //                            -> Photon -> PhotonHad -> MuEG
  //   high pt single isolated electron chain: Electron -> SingleElectron 
  //                            -> ElectronHad -> MuEG
  //   low pt double electron chain: ZeroBias -> MinumumBias -> EGmonitor 
  //                            -> Electron -> SingleElectron -> DoubleElectron
  //                            -> ElectronHad -> Photon -> SinglePhoton 
  //                            -> DoublePhoton -> DoublePhotonHighPt 
  //                            -> PhotonHad 
  //   high pt isolated double electron chain: Electron -> SingleElectron 
  //                            -> DoubleElectron -> ElectronHad  
  //   low/high pt muon-electron chain: "Or" of single muon and single electron
  //                                    chains  (low or high pt) 
  //   Others 2010:      Commissioning, JetMETTaumonitor, BTau, METFwd 
  //   Others 2011:      TauPlusX, BTag, MET, METBTag, Tau 
  //   Others 2012:      Commissioning, MET, BTag, MTMHTParked, HcalNZS, 
  //                     NoBPTX, TauParked, HCalNZS, BJetPlusX, VBF1Parked  

////////////////////////////// for Gen particle //////////////////////////////

  /// nanoAOD ///     
  //                   
  UInt_t nGenPart;
  vector<Float_t> GenPart_pt;
  vector<Float_t> GenPart_eta;
  vector<Float_t> GenPart_phi;
  vector<Float_t> GenPart_mass;
  vector<Int_t> GenPart_pdgId;
  vector<Int_t> GenPart_status;
  vector<Int_t> GenPart_statusFlags;
  vector<Int_t> GenPart_genPartIdxMother; 

  /// nanoAOD extension ///
  // later, store the above for all(?) GenParticles //
  vector<Int_t> GenPart_Id;
  vector<uint8_t> GenPart_isNano;  // satisfies criteria for nanoAOD shortlist 
                                   // (Bool on ntuple)
  vector<Int_t> GenPart_parpdgId;  // PDG id of parent particle
  vector<Int_t> GenPart_sparpdgId; // PDG id of parent stable particle
  vector<Int_t> GenPart_numberOfDaughters; 
                                   // number of daughter particles in decay
  vector<Int_t> GenPart_nstchgdaug; 
                    // number of stable charged daughter particles in decay
  vector<Float_t> GenPart_vx;
  vector<Float_t> GenPart_vy;
  vector<Float_t> GenPart_vz;
  vector<Float_t> GenPart_mvx;
  vector<Float_t> GenPart_mvy;
  vector<Float_t> GenPart_mvz;
  vector<Int_t> GenPart_recIdx;

  // main generated vertex, to be filled from suitably chosen GenPart_vtx
  Float_t GenPV_x;
  Float_t GenPV_y;
  Float_t GenPV_z;
  Int_t GenPV_recIdx;  
  Int_t GenPV_chmult;  

////////////////////////////// for track variables ///////////////////////////
   
  /// nanoAOD ///            *** not yet filled, not yet written out ***
  UInt_t nIsoTrack; // number of isolated tracks (from main prim. vertex only)
  vector<Float_t> IsoTrack_dxy;
  vector<Float_t> IsoTrack_dz;
  vector<Float_t> IsoTrack_ets;
  vector<uint8_t> IsoTrack_isHighPurityTrack;
  vector<uint8_t> IsoTrack_isPFcand;
  vector<Float_t> IsoTrack_miniPFreliso_all;
  vector<Float_t> IsoTrack_miniPFreliso_chg;
  vector<Int_t>   IsoTrack_pdgId;
  vector<Float_t> IsoTrack_PFreliso03_all;
  vector<Float_t> IsoTrack_PFreliso03_chg;
  vector<Float_t> IsoTrack_phi;
  vector<Float_t> IsoTrack_pt;

  /// nanoAOD extension ///  
  UInt_t nTrk;               // number of tracks in event
  // *** from here not yet filled, not yet written out ***
  vector<Int_t> Trk_Id;      // unique track identifier
  vector<Int_t> Trk_pvIdx;   // if nonzero: primary vertex id (PVtx_id) this 
          // track is associated (but not necessarily fitted) to
          // if zero: the track is not associated to any primary vertex 
  vector<Int_t> Trk_pvaltIdx;// if nonzero: alternative primary vertex id 
          // (Pvtx_id) this track might be associated to  
          // i.e. the vertex association might be ambigous
  vector<Int_t> Trk_pvfIdx;  // if nonzero: primary vertex id (PVtx_id) this 
                             // track is fitted to    
  vector<Int_t> Trk_svfIdx;  // if nonzero: secondary vertex id (SVtx_id) this 
                             // track is fitted to
  vector<Int_t> Trk_isIsoTrack; // appears in IsoTrack list yes(1)/no(0)
  vector<Int_t> Trk_charge;  // track charge 
  vector<Float_t> Trk_pt;    // track transverse momentum
  vector<Float_t> Trk_eta;   // track pseudorapidity
  vector<Float_t> Trk_phi;   // track phi
  // Add here some track quality variables //


////////////////////////////// for vertex variables ///////////////////////// 

  /// official nanoAOD ///
  Int_t   PV_npvs;           // number of reconstructed primary vertices
  Int_t   PV_npvsGood;       // number of good reconstructed primary vertices
          // selection: !isFake && ndof>4 && abs(z) <= 24 && position.rho <=2
  Float_t PV_chi2;           // main primary vertex chi2 
  Float_t PV_ndof;           // main primary vertex number of degr. of freedom
  Float_t PV_score;          // main primary vertex score, 
                             // i.e. sum pt2 of clustered objects
  Float_t PV_x;              // main primary vertex position x coordinate
  Float_t PV_y;              // main primary vertex position y coordinate
  Float_t PV_z;              // main primary vertex position z coordinate
  UInt_t nOtherPV;          // number of other primary vertices, 
                             // excluding the main PV
  vector<Float_t> OtherPV_z; // z position of other primary vertices, 
                             // excluding the main PV 

  /// nanoAOD extension ///
  UInt_t nPVtx;                // number of all primary vertex candidates
  vector<Int_t>   PVtx_Id;     // primary vertex identifier
  vector<uint8_t> PVtx_isMain; // corresponds to main prim. in PV yes(1)/no(0)
                               //  (highest score)
  vector<uint8_t> PVtx_isMainSim; // corresponds to main simulated primary 
         // yes(1)/no(0) (=vertex for which generator information is available)
  vector<uint8_t> PVtx_isGood; // vertex counted in PV_npvsGood yes(1)/no(0)
  vector<uint8_t> PVtx_isValid; // vertex fit converged properly? yes(1)/no(0)
  vector<uint8_t> PVtx_isFake; // vertex is made up? yes(1)/no(0)
  vector<uint8_t> PVtx_isTrigUnique; // This unambigously is the vertex of 
         // the event that fired the trigger (e.g. lepton); may contain some 
         // pileup, may not necessarily coincide with main vertex 
  vector<uint8_t> PVtx_isUnbiased; // This is a vertex which did not contribute 
         // to firing the trigger (no trigger pileup). 
         // *** Can be used as next-to-minimum bias. ***
         // In cases where several vertices might have significantly 
         // contributed to the trigger (e.g. multileptons from different 
         // vertices, 
         // pileup contribution to energies, ...) neither of the two previous 
         // flags will be set for these vertices
  vector<Int_t> PVtx_ntrk;   // number of tracks associated to this primary 
         // vertex (including secondaries from this primary) *** to be impl. ***
  vector<Int_t> PVtx_ntrkfit;// number of tracks fitted to this primary vertex
  vector<Float_t> PVtx_chi2;   // primary vertex chi2
  vector<Float_t> PVtx_ndof;   // primary vertex number of degrees of freedom
  vector<Float_t> PVtx_score;  // score (sum pt^2) for determination of PV
  vector<Float_t> PVtx_sumPt;  // hadronic sum Pt for this vertex (tracks only,
                               // excluding leptons)
  vector<Float_t> PVtx_Rho;    // primary vertex Rho (sqrt(x^2+y^2))
  vector<Float_t> PVtx_x;      // primary vertex x
  vector<Float_t> PVtx_y;      // primary vertex y
  vector<Float_t> PVtx_z;      // primary vertex z
  vector<Float_t> PVtx_Covxx;  // primary vertex covariance matrix (6 entries)
  vector<Float_t> PVtx_Covyx; 
  vector<Float_t> PVtx_Covzx; 
  vector<Float_t> PVtx_Covyy; 
  vector<Float_t> PVtx_Covzy; 
  vector<Float_t> PVtx_Covzz; 
  // In order to find tracks associated to a particular vertex PVtx_id, 
  // loop over the track list and check for Trk_pvIdx == PVtx_id   

  Float_t Bsp_x;              // Beam spot x0
  Float_t Bsp_y;              // Beam spot y0
  Float_t Bsp_z;              // Beam spot z0
  Float_t Bsp_sigmaz;         // Beam spot size in z
  Float_t Bsp_dxdz;           // Beam spot slope in xz
  Float_t Bsp_dydz;           // Beam spot slope in yz
  Float_t Bsp_widthx;         // Beam spot width in x
  Float_t Bsp_widthy;         // Beam spot width in y
  // width is of order 15-20 micron in both x and y
  // uncertainties of slope and width?
  
////////////////////////////// for Muon variables ////////////////////////////

  // variables // *   are not yet filled properly
  
  UInt_t b4_nMuon;  // all muon counter;  (b4 = before selection)
  UInt_t Muon_nNano; // counter for official nanoAOD muons

  /// official nanoAOD ///

  UInt_t nMuon;                   // number of stored muons
  vector<Int_t> Muon_charge;      // muon charge
  vector<Int_t> Muon_tightCharge; // *
  // for the following, take the parameters from global track if available, 
  // from tracker or muon track otherwise 
  // *** or take from vertex? no ->revertex? ***
  // seems to be always taken from inner track if available
  vector<Float_t> Muon_pt;        // muon transverse momentum (in GeV)
  vector<Float_t> Muon_ptErr;     // muon transverse momentum error
  vector<Float_t> Muon_eta;       // muon pseudorapidity
  vector<Float_t> Muon_phi;       // muon azimuth angle
  vector<Float_t> Muon_mass;      // muon mass (redundant ...)
  // the following are taken w.r.t. the main primary vertex (CMS default)
  // see also variable extensions w.r.t. associated primary
  vector<Float_t> Muon_dxy;       // distance in xy (in cm)
  vector<Float_t> Muon_dxyErr;    // error in xy
  vector<Float_t> Muon_dz;        // distance in z
  vector<Float_t> Muon_dzErr;     // error in z
  vector<Float_t> Muon_ip3d;      // 3D impact parameter
  vector<Float_t> Muon_sip3d;     // 3D impact parameter significance
  // for the following, take the information from ...
  vector<Float_t> Muon_pfRelIso03_all;   // PF isolation in cone 0.3
  vector<Float_t> Muon_pfRelIso03_chg;   // PF track isolation in cone 0.3
  vector<Float_t> Muon_pfRelIso04_all;   // PF isolation in cone 0.4
  vector<Float_t> Muon_miniPFRelIso_all; // *
  vector<Float_t> Muon_miniPFRelIso_chg; // *
  vector<Int_t> Muon_jetIdx;   // *
  // all nanoAOD bool arrays are declared uint8_t (bool will not work)
  vector<uint8_t> Muon_isGlobal; // muon is global muon; 
                               // track parameters are global parameters;
                               // tracker only track parameters can be found 
                               // from TRK list, when filled   
  vector<uint8_t> Muon_isTracker;// muon is tracker muon;
  vector<uint8_t> Muon_isPFcand; // muon is PF candididate
  vector<uint8_t> Muon_softId;   // muon satisfies soft ID 
  vector<uint8_t> Muon_mediumId; // muon satisfies medium Id 
  vector<uint8_t> Muon_tightId;  // muon satisfies tight Id 
  vector<UChar_t> Muon_highPtId; // muon satisfies high pt Id, not bool!
  vector<Int_t> Muon_nStations;      // number of muon stations hit
  vector<Int_t> Muon_nTrackerLayers; // number of tracker layers hit
  vector<Float_t> Muon_segmentComp;  // muon segnent compatibility
  vector<UChar_t> Muon_cleanmask;    // *
  vector<Float_t> Muon_mvaTTH;       // *
  vector<Int_t> Muon_pdgId;          // +-13, depending on (-1)*charge
  vector<UChar_t> Muon_genPartFlav;  // *
  // clarify overlap with Muon_simIdx
  vector<Int_t> Muon_genPartIdx;     // *

  /// nanoAOD extension ///
  // store all muon candidates
  vector<Int_t> Muon_Id;       // unique muon identifier
  // add more position and impact parameter info
  vector<Float_t> Muon_x;      // muon track x position at dca to beamspot
  vector<Float_t> Muon_y;      // muon track y position at dca to beanspot
  vector<Float_t> Muon_z;      // muon track z position at dca to beamspot 
  vector<Float_t> Muon_dxyBest;   // distance in xy to Muon_vtxIdx vertex (cm)
  vector<Float_t> Muon_dzBest;    // distance in z to Muon_vtxIdx vertex (cm)
  vector<Float_t> Muon_ip3dBest;  // impact parameter to Muon_vtxIdx vertex(cm)
  vector<Float_t> Muon_sip3dBest; // impact parameter significance
  // add also full covariance matrix here?   
  vector<Float_t> Muon_gpt;      // global muon transverse momentum (in GeV)
  vector<Float_t> Muon_geta;     // global muon pseudorapidity
  vector<Float_t> Muon_gphi;     // global muon azimuth angle
  vector<uint8_t> Muon_looseId;  // muon satisfies loose Id
  vector<uint8_t> Muon_softId4;  // muon satisfies "old" (CMSSW4X) soft Id
  vector<uint8_t> Muon_softIdBest; // soft Id for best vertex (instead of main)
  vector<uint8_t> Muon_isNano;   // muon satisfies criteria for off. nanoAOD 
                               // (to restrict to original list)
  vector<uint8_t> Muon_isMini;   // muon satisfies criteria for MiniAOD;
  vector<uint8_t> Muon_isGood;   // muon satisfies TMOneStationTight;
  vector<uint8_t> Muon_isGoodLast; // muon satisfies TMLastStationTight;
  vector<uint8_t> Muon_isGoodAng; // muon satisfies TMLastStationAnyTight;
  vector<uint8_t> Muon_isArbitrated; // muon satisfies TrackerMuonArbitrated;
  vector<uint8_t> Muon_isStandAlone; // muon is standalone muon;
  vector<uint8_t> Muon_isRPCcand; // muon is RPC candididate
  vector<Int_t> Muon_nValid;   // number of valid hits;
  vector<Int_t> Muon_nPix;     // number of pixel hits;
  vector<Float_t> Muon_Chi2;   // tracker muon chi2/ndof;
  vector<Int_t> Muon_gnValid;  // number of valid global muon hits;
  vector<Int_t> Muon_gnPix;    // number of global muon pixel hits;
  vector<Float_t> Muon_gChi2;  // global muon chi2/ndof;
  vector<Int_t> Muon_gnValidMu;// number of valid muon hits in global muon;
  vector<Int_t> Muon_vtxIdx;   // index (PVtx_Id) of vertex in PVtx list,if any
  vector<Int_t> Muon_vtxFlag;  // quality flag for vertex association
  vector<Int_t> Muon_trkIdx;   // index (Trk_Id) of track in TRK list
  // clarify overlap with Muon_genPartIdx
  vector<Int_t> Muon_simIdx;   // index (GenPart_Id) of particle in GenPart
    
  // store dimuon candidates 
  //------------------------------ For Dimu branches ------------------------//
  UInt_t nDimu;  
  // mu 1 from Dimu
  vector<Int_t> Dimut1_muIdx;    // pointer to first muon
  vector<Float_t> Dimut1_dxy;    // dxy of first muon w.r.t. allocated vertex
  vector<Float_t> Dimut1_dz;     // dz of first muon w.r.t. allocated vertex

  // mu 2 from Dimu
  vector<Int_t> Dimut2_muIdx;    // pointer to second muon
  vector<Float_t> Dimut2_dxy;    // dxy of second muon w.r.t. allocated vertex
  vector<Float_t> Dimut2_dz;     // dz of second muon w.r.t. allocated vertex

  // Dimu
  vector<Float_t> Dimu_pt;       // Dimuon pt  after refit
  vector<Float_t> Dimu_eta;      // Dimuon eta after refit
  vector<Float_t> Dimu_phi;      // Dimuon phi after refit
  vector<Float_t> Dimu_rap;      // Dimuon rapidity after refit
  vector<Float_t> Dimu_mass;     // Dimuon mass
  vector<Int_t>   Dimu_charge;   // Dimuon charge
  vector<Int_t>   Dimu_simIdx;   // matched true gamma/Z/meson in genparticle
  vector<Int_t>   Dimu_vtxIdx;   // associated prim. vtx (can differ from 1,2)
  vector<Float_t> Dimu_chi2;     // chi2 of Dimuon vertex
  vector<Float_t> Dimu_dlxy;     // Dimuon decay length in xy
  vector<Float_t> Dimu_dlxyErr;  // Dimuon decay length uncertainty in xy
  vector<Float_t> Dimu_dlxySig;  // Dimuon decay length significance in xy
  vector<Float_t> Dimu_cosphixy; // cosine of angle (momentum, decay length) xy
  vector<Float_t> Dimu_dl;       // Dimuon decay length in 3D (typically > xy)
  vector<Float_t> Dimu_dlErr;    // Dimuon decay length uncertainty in 3D (>xy)
  vector<Float_t> Dimu_dlSig;    // Dimuon decay length significance in 3D 
  vector<Float_t> Dimu_cosphi;   // cosine of angle (momentum, decay length) 3D
  vector<Float_t> Dimu_ptfrac;   // Dimuon_pt/sum pt at vertex *** not final***
  vector<Float_t> Dimu_x;        // Dimuon vertex x
  vector<Float_t> Dimu_y;        // Dimuon vertex y
  vector<Float_t> Dimu_z;        // Dimuon vertex z
  vector<Float_t> Dimu_Covxx;    // Dimuon vertex covariance
  vector<Float_t> Dimu_Covyx;
  vector<Float_t> Dimu_Covzx;
  vector<Float_t> Dimu_Covyy;
  vector<Float_t> Dimu_Covzy;
  vector<Float_t> Dimu_Covzz;


//////////////////////////// for S.Wunsch variables ///////////////

  /// nanoAOD subset ///

// Electrons
  const static int max_el = 128;
  UInt_t value_el_n;
  float value_el_pt[max_el];
  float value_el_eta[max_el];
  float value_el_phi[max_el];
  float value_el_mass[max_el];
  Int_t value_el_charge[max_el];
  // nuha
  Int_t value_el_tightCharge[max_el];
  
  float value_el_pfreliso03all[max_el];
  float value_el_pfreliso03chg[max_el];
  // nuha
  float value_el_dr03TkSumPtOld[max_el];
  float value_el_dr03TkSumPt[max_el];  
  float value_el_dr03EcalRecHitSumEtOld[max_el];
  float value_el_dr03EcalRecHitSumEt[max_el];  

  float value_el_dr03HcalTowerSumEt[max_el];
  // nuha
  float value_el_dr03HcalDepth1TowerSumEtOld[max_el];
  float value_el_dr03HcalDepth1TowerSumEt[max_el];

  bool  value_el_isEB[max_el];
  bool  value_el_isEE[max_el];
  UChar_t value_el_lostHits[max_el];
  float value_el_convDist[max_el];
  float value_el_convDcot[max_el];
  // nuha
  bool  value_el_convVetoOld[max_el];
  bool  value_el_convVeto[max_el];
  
  float value_el_deltaEtaSC[max_el];
  float value_el_deltaPhiSC[max_el];
  float value_el_deltaEtaSCtr[max_el];
  float value_el_deltaPhiSCtr[max_el];
  float value_el_hoe[max_el];
  float value_el_sieieR1[max_el];
  float value_el_sieie[max_el];
  // nuha
  float value_el_eInvMinusPInvOld[max_el];
  float value_el_eInvMinusPInv[max_el];
  
  float value_el_SCeta[max_el];
  Int_t value_el_cutBased[max_el];

  float value_el_dxy[max_el];
  float value_el_dxyErr[max_el];
  float value_el_dz[max_el];
  float value_el_dzErr[max_el];

  float value_el_ip3d[max_el];
  float value_el_sip3d[max_el];
  
  bool value_el_isPFcand[max_el];   
  bool value_el_isNano[max_el];

  UInt_t Electron_nNano; // counter for official nanoAOD electrons

// Photons
  const static int max_ph = 100;
  UInt_t value_ph_n;
  float value_ph_pt[max_ph];
  float value_ph_eta[max_ph];
  float value_ph_phi[max_ph];
  float value_ph_mass[max_ph];
  int value_ph_charge[max_ph];
  float value_ph_pfreliso03all[max_ph];

// Taus
  const static int max_tau = 100;
  UInt_t value_tau_n;
  float value_tau_pt[max_tau];
  float value_tau_eta[max_tau];
  float value_tau_phi[max_tau];
  float value_tau_mass[max_tau];
  int value_tau_charge[max_tau];
  int value_tau_decaymode[max_tau];
  float value_tau_chargediso[max_tau];
  float value_tau_neutraliso[max_tau];

// MET
  float value_met_pt;
  float value_met_phi;
  float value_met_sumEt;
  float value_met_significance;
  float value_met_covxx;
  float value_met_covxy;
  float value_met_covyy;

// CaloMET
  float value_calomet_pt;
  float value_calomet_phi;
  float value_calomet_sumEt;

// Jets
  const static int max_jet = 300;
  UInt_t value_jet_n;
  float value_jet_pt[max_jet];
  float value_jet_eta[max_jet];
  float value_jet_phi[max_jet];
  float value_jet_mass[max_jet];
//  float value_jet_ptD[max_jet];
  float value_jet_area[max_jet];
  int value_jet_nConstituents[max_jet];
  int value_jet_nElectrons[max_jet];
  int value_jet_nMuons[max_jet];
  float value_jet_chEmEF[max_jet];
  float value_jet_chHEF[max_jet];
  float value_jet_neEmEF[max_jet];
  float value_jet_neHEF[max_jet];

// Flags
  edm::InputTag custom_tag;
  std::vector<std::string> custom_flag;
  std::vector<uint8_t> custom_bit;

//////////////////////////// for Dmeson variables ////////////////////////////
  
  /// nanoAOD extension /// 
    
  //------------------------------ For D0 branches --------------------------//
  UInt_t nD0;  
  // track 1 from D0
  vector<Float_t> D0t1_pt;       // track 1 pt  after refit
  vector<Float_t> D0t1_eta;      // track 1 eta after refit
  vector<Float_t> D0t1_phi;      // track 1 phi after refit
  vector<Int_t> D0t1_chg;        // track 1 charge
  vector<Int_t> D0t1_tkIdx;      // *
  vector<Float_t> D0t1_Kprob;    // * temporarily misused for Strip dEdx
  vector<Float_t> D0t1_piprob;   // * temporarily misused for Strip dEdxErr
  vector<Int_t> D0t1_dEdxnmeas;  // # of dEdx measurements
  vector<Int_t> D0t1_dEdxnsat;   // # of saturated dEdx measurements
  vector<Int_t> D0t1_vtxIdx;     // track 1 primary vertex before refit 
  vector<Float_t> D0t1_chindof;  // track 1 chi2/ndof
  vector<Int_t> D0t1_nValid;     // track 1 Valid Tracker hits
  vector<Int_t> D0t1_nPix;       // track 1 Valid Pixel Hits
  vector<uint8_t> D0t1_isHighPurity; // high purity flag
  vector<Float_t> D0t1_dxy;      // track 1 dxy w.r.t. *** what? ***
  vector<Float_t> D0t1_dz;       // track 1 dz w.r.t. *** what? ***
  vector<Int_t> D0t1_pdgId;      // * 

  // track 2 from D0
  vector<Float_t> D0t2_pt;       // track 2 pt  after refit
  vector<Float_t> D0t2_eta;      // track 2 eta after refit
  vector<Float_t> D0t2_phi;      // track 2 phi after refit
  vector<Int_t> D0t2_chg;        // track 2 charge
  vector<Int_t> D0t2_tkIdx;      // *
  vector<Float_t> D0t2_Kprob;    // * temporarily misused for Strip dEdx
  vector<Float_t> D0t2_piprob;   // * temporarily misused for Strip dEdxErr
  vector<Int_t> D0t2_dEdxnmeas;  // # of dEdx measurements
  vector<Int_t> D0t2_dEdxnsat;   // # of saturated dEdx measurements
  vector<Int_t> D0t2_vtxIdx;     // track 2 primary vertex before refit
  vector<Float_t> D0t2_chindof;  // track 2 chi2/ndof
  vector<Int_t> D0t2_nValid;     // track 2 Valid Tracker hits
  vector<Int_t> D0t2_nPix;       // track 2 Valid Pixel hits
  vector<uint8_t> D0t2_isHighPurity; // high purity flag
  vector<Float_t> D0t2_dxy;      // track 2 dxy w.r.t. *** what? ***
  vector<Float_t> D0t2_dz;       // track 2 dz  w.r.t. *** what? ***
  vector<Int_t> D0t2_pdgId;      // *

  // D0
  vector<Float_t> D0_pt;         // D0 pt  after refit
  vector<Float_t> D0_eta;        // D0 eta after refit
  vector<Float_t> D0_phi;        // D0 phi after refit
  vector<Float_t> D0_rap;        // D0 rapidity after refit
  vector<Float_t> D0_mass12;     // D0 mass 1=K
  vector<Float_t> D0_mass21;     // D0 mass 2=K
  vector<Float_t> D0_massKK;     // D0 mass both=K
  vector<Int_t>   D0_simIdx;     // matched true D0 in genparticle list
  vector<Int_t>   D0_DstarIdx;   // equivalent D0 in Dstar list
  vector<uint8_t> D0_ambiPrim;   // flag indicating ambigous primary assignment
  vector<Int_t>   D0_vtxIdx;     // associated prim. vtx (can differ from 1,2)
  vector<uint8_t> D0_hasMuon;    // true if either "K" or "pi" is muon
  vector<Float_t> D0_chi2;       // chi2 of D0 vertex
  vector<Float_t> D0_dlxy;       // D0 decay length in xy
  vector<Float_t> D0_dlxyErr;    // D0 decay length uncertainty in xy
  vector<Float_t> D0_dlxySig;    // D0 decay length significance in xy
  vector<Float_t> D0_cosphixy;   // cosine of angle (momentum, decay length) xy
  vector<Float_t> D0_dl;         // D0 decay length in 3D (typically > xy)
  vector<Float_t> D0_dlErr;      // D0 decay length uncertainty in 3D (>xy)
  vector<Float_t> D0_dlSig;      // D0 decay length significance in 3D 
  vector<Float_t> D0_cosphi;     // cosine of angle (momentum, decay length) 3D
  vector<Float_t> D0_ptfrac;     // D0_pt/sum pt at vertex *** not final ***
  vector<Float_t> D0_ptfrac15;   // D0_pt/sum pt at vertex in cone 1.5
  vector<Float_t> D0_ptfrac10;   // D0_pt/sum pt at vertex in cone 1.0
  vector<Float_t> D0_ptfrac07;   // D0_pt/sum pt at vertex in cone 0.7
  vector<Float_t> D0_ptfrac04;   // D0_pt/sum pt at vertex in cone 0.4
  vector<Float_t> D0_x;          // D0 vertex x
  vector<Float_t> D0_y;          // D0 vertex y
  vector<Float_t> D0_z;          // D0 vertex z
  vector<Float_t> D0_Covxx;      // D0 vertex covariance
  vector<Float_t> D0_Covyx;
  vector<Float_t> D0_Covzx;
  vector<Float_t> D0_Covyy;
  vector<Float_t> D0_Covzy;
  vector<Float_t> D0_Covzz;


  //----------------------------- For Dstar branches ------------------------//
  UInt_t nDstar;

  // Slow Pion from Dstar
  vector<Float_t> Dstarpis_pt;      // Dstar slow pion pt before refit
  vector<Float_t> Dstarpis_eta;     // Dstar slow pion eta before refit
  vector<Float_t> Dstarpis_phi;     // Dstar slow pion phi before refit
  vector<Float_t> Dstarpis_ptr;     // Dstar slow pion pt after refit
  vector<Float_t> Dstarpis_etar;    // Dstar slow pion eta after refit
  vector<Float_t> Dstarpis_phir;    // Dstar slow pion phi after refit
  vector<Int_t> Dstarpis_chg;       // Dstar slow pion charge 
                                    // (=charge of Dstar for right sign)
  vector<Int_t> Dstarpis_tkIdx;     // *
  vector<Float_t> Dstarpis_Kprob;    // * temporarily misused for Strip dEdx
  vector<Float_t> Dstarpis_piprob;   // * temporarily misused for Strip dEdxErr
  vector<Int_t> Dstarpis_dEdxnmeas;  // # of dEdx measurements
  vector<Int_t> Dstarpis_dEdxnsat;   // # of saturated dEdx measurements
  vector<Int_t> Dstarpis_vtxIdx;    // slow pion primary vertex before refit
  vector<Float_t> Dstarpis_chindof; // slow pion chi2/ndof 
  vector<Float_t> Dstarpis_chir;    // slow pion refit chi2
  vector<Int_t> Dstarpis_nValid;    // slow pion Valid Tracker hits
  vector<Int_t> Dstarpis_nPix;      // slow pion Valid Pixel hits
  vector<Float_t> Dstarpis_dxy;     // slow pion dx w.r.t. *** what? ***
  vector<Float_t> Dstarpis_dz;      // slow pion dz w.r.t. *** what? ***


  // D0 from Dstar
  vector<Float_t> DstarD0_pt;       // D0 from D* pt
  vector<Float_t> DstarD0_eta;      // D0 from D* eta
  vector<Float_t> DstarD0_phi;      // D0 from D* phi
  vector<Float_t> DstarD0_mass;     // D0 from D* mass
  vector<Float_t> DstarD0_chi2;     // D0 from D* vertex chi2
  vector<Float_t> DstarD0_dlxy;     // D0 from D* decay length xy
  vector<Float_t> DstarD0_dlxyErr;  // D0 from D* decay length error xy
  vector<Float_t> DstarD0_dlxySig;  // D0 from D* deacy length significance xy
  vector<Float_t> DstarD0_cosphixy; // D0 from D* cosine (momentum, dl) xy
  vector<Float_t> DstarD0_dl;       // D0 from D* decay length 3D
  vector<Float_t> DstarD0_dlErr;    // D0 from D* decay length error 3D
  vector<Float_t> DstarD0_dlSig;    // D0 from D* decay length significance 3D
  vector<Float_t> DstarD0_cosphi;   // D0 from D* cosine (momentum, dl) 3D
  vector<Float_t> DstarD0_ptfrac;   // D0 pt/all had sum pt at vtx (not final)
  vector<Float_t> DstarD0_ptfrac15; // D0 pt/all had sum pt at vtx in cone 1.5
  vector<Float_t> DstarD0_ptfrac10; // D0 pt/all had sum pt at vtx in cone 1.0
  vector<Float_t> DstarD0_ptfrac07; // D0 pt/all had sum pt at vtx in cone 0.7
  vector<Float_t> DstarD0_ptfrac04; // D0 pt/all had sum pt at vtx in cone 0.4
  vector<Float_t> DstarD0_x;        // D0 vertex x
  vector<Float_t> DstarD0_y;        // D0 vertex y
  vector<Float_t> DstarD0_z;        // D0 vertex z
  vector<Int_t>   DstarD0_simIdx;   // matched true D0 in genparticle list
  vector<Int_t>   DstarD0_recIdx;   // equivalent D0 in D0 list
  vector<uint8_t> DstarD0_ambiPrim; // flag indicating ambigous prim. vtx ass. 

  // Kaon from Dstar
  vector<Float_t> DstarK_pt;        // Kaon pt
  vector<Float_t> DstarK_eta;       // Kaon eta
  vector<Float_t> DstarK_phi;       // Kaon phi
  vector<Int_t> DstarK_chg;         // Kaon charge
  vector<Int_t> DstarK_tkIdx;       // *
  vector<Float_t> DstarK_Kprob;     // * temporarily misused for Strip dEdx
  vector<Float_t> DstarK_piprob;    // * temporarily misused for Strip dEdxErr
  vector<Int_t> DstarK_dEdxnmeas;   // # of dEdx measurements
  vector<Int_t> DstarK_dEdxnsat;    // # of saturated dEdx measurements
  vector<Int_t> DstarK_vtxIdx;      // Kaon primary vertex before refit 
  vector<Float_t> DstarK_chindof;   // Kaon chi2/ndof
  vector<Int_t> DstarK_nValid;      // Kaon Valid Tracker hits
  vector<Int_t> DstarK_nPix;        // Kaon Pixel hits
  vector<uint8_t> DstarK_isHighPurity; // high purity flag
  vector<Float_t> DstarK_dxy;       // Kaon dxy w.r.t. *** what? ***
  vector<Float_t> DstarK_dz;        // Kaon dz w.r.t. *** what? ***

  // Pion from Dstar
  vector<Float_t> Dstarpi_pt;       // Pion pt
  vector<Float_t> Dstarpi_eta;      // Pion eta 
  vector<Float_t> Dstarpi_phi;      // Pion phi
  vector<Int_t> Dstarpi_chg;        // Pion charge (wrong sign allowed)
  vector<Int_t> Dstarpi_tkIdx;      // *
  vector<Float_t> Dstarpi_Kprob;    // * temporarily misused for Strip dEdx
  vector<Float_t> Dstarpi_piprob;   // * temporarily misused for Strip dEdxErr
  vector<Int_t> Dstarpi_dEdxnmeas;  // # of dEdx measurements
  vector<Int_t> Dstarpi_dEdxnsat;   // # of saturated dEdx measurements
  vector<Int_t> Dstarpi_vtxIdx;     // Pion primary vertex before refit  
  vector<Float_t> Dstarpi_chindof;  // Pion chi2/ndof
  vector<Int_t> Dstarpi_nValid;     // Pion Valid Tracker hits
  vector<Int_t> Dstarpi_nPix;       // Pion Valid Pixel hits
  vector<uint8_t> Dstarpi_isHighPurity; // high purity flag
  vector<Float_t> Dstarpi_dxy;      // Pion dxy w.r.t *** what? ***
  vector<Float_t> Dstarpi_dz;       // Pion dz w.r.t *** what? ***
  
  // Dstar
  vector<Float_t> Dstar_pt;         // D* pt  after D0 refit, w/o pis refit
  vector<Float_t> Dstar_eta;        // D* eta after D0 refit, w/o pis refit
  vector<Float_t> Dstar_phi;        // D* phi after D0 refit, w/o pis refit
  vector<Float_t> Dstar_rap;        // D* rapidity after D0 refit, w/o pis refit
  vector<Float_t> Dstar_deltam;     // mD*-mD0 after D0 refit, w/o pis refit 
  vector<Float_t> Dstar_deltamr;    // mD*-mD0 after D0 refit, with pis refit 
  vector<Int_t> Dstar_simIdx;       // associated true D* in Genpart, if any
  vector<Int_t> Dstar_vtxIdx;       // ass. prim. vtx (may differ from tracks)
  vector<uint8_t> Dstar_hasMuon;    // true if "K" or "pi" or "pislow" is muon
  vector<Float_t> Dstar_ptfrac;     // pt(D*)/pt(all hadrons at prim. vertex)

  // flag for printout of "good" D0 candidates
  int Bingo;
  float D0vtxx, D0vtxy, D0vtxz;
  float D0motx, D0moty, D0motz;

// *** move this to a different place ***
  // Qun below  Trigger Object
  UInt_t nTrigObj;                   // number of stored Trigger Obj
  vector<Int_t> TrigObj_id;          // ID of the object
  vector<Int_t> TrigObj_filterBits;  // filter bits - *** still to be implemented ***
  vector<Float_t> TrigObj_pt;        // pt
  vector<Float_t> TrigObj_phi;        // phi
  vector<Float_t> TrigObj_eta;        // eta
  //vector<Int_t> TrigObj_id;          // ID of the object
  //vector<Float_t> TrigObj_pt;        // pt
  //vector<Float_t> TrigObj_phi;        // phi
  //vector<Float_t> TrigObj_eta;        // eta
  //Qconst static int max_TrigObj = 300;
  //QInt_t  TrigObj_id[max_TrigObj];          // ID of the object
  //QFloat_t  TrigObj_pt[max_TrigObj];        // pt
  //QFloat_t TrigObj_phi[max_TrigObj];        // phi
  //QFloat_t TrigObj_eta[max_TrigObj];        // eta
  // Qun above  Trigger Object

}; // end of class member

//
// constants, enums and typedefs
//
const float emass = 0.000510999; // [PDG]
const float pi = 3.141593;
const float mumass = 0.105658;
const float Kmass = 0.49367;
const float pimass = 0.13957;
const float mD0Actual = 1.86484;  // [PDG]
const float mDstarActual = 2.0103;  // [PDG]
const float dmDstarActual = 0.145426;  // [PDG]
// reserve array sizes; reconsider when considering new datasets! 
const unsigned nReserve_GenPart = 1024; 
const unsigned nReserve_Track = 2048; 
const unsigned nReserve_IsoTrack = 32; 
// at most 3 "otherPV" vertices are stored on official nanoAOD
// see PVtx structure for full list
const unsigned nReserve_OtherPV = 3; 
//const unsigned nReserve_OtherPV = 128; 
const unsigned nReserve_PVtx = 128; 
const unsigned nReserve_Muon = 128;
const unsigned nReserve_TrigObj = 1024;
const unsigned nReserve_Dimu = 128;
const unsigned nReserve_D0 = 1024;
const unsigned nReserve_Dstar = 1024;
const int maxnmusim = 128;
const int maxnD0sim = 128;
const int maxnDstarsim = 128;
const int maxnmuonlist = 128;

// preset flag for trigger steering
bool skiptrigger = false;
// preset first event flag 
bool firstevent = true;

//
// static data member definitions
//

//
// constructors and destructor
//
// none

//////////////////////////////////////////////////////////////////////////////
//                        set analysis loop parameters                      //
//////////////////////////////////////////////////////////////////////////////

NanoAnalyzer::NanoAnalyzer(const edm::ParameterSet& iConfig)
{
  // configure parameters
  outFile = iConfig.getParameter<std::string>("outFile");
  isData = iConfig.getParameter<bool>("isData");
  hlt_proc = iConfig.getParameter<std::string>("hltProcess");
  nanoext = iConfig.getParameter<bool>("nanoExtension");
  covout = iConfig.getParameter<bool>("writeCovariance");
  custom_flag = iConfig.getParameter<std::vector<std::string> >("customFlag");
  if (!custom_flag.empty())
    custom_tag = iConfig.getParameter<edm::InputTag>("customTag");
  // use hlt_proc above
  //processName_ = iConfig.getParameter<std::string>("processName"); //Qun
  triggerName_ = iConfig.getParameter<std::string>("triggerName"); //Qun

  // Make sure uncertainty in histogram bins is calculated correctly
  TH1::SetDefaultSumw2(true);

  // nuha
  // https://github.com/cms-sw/cmssw/blob/CMSSW_9_4_X/PhysicsTools/NanoAOD/python/electrons_cff.py#L100
  // https://github.com/cms-sw/cmssw/blob/CMSSW_9_4_X/PhysicsTools/NanoAOD/plugins/IsoValueMapProducer.cc#L49
  //relative_ = iConfig.getParameter<bool>("relative");
  
  // afiqaize get actual CMSSW version 
  // however being that the branch above saves int, strip out all non-digits from the output
  // FIXME this ignores the ifdef flag (likely not an issue)
  CMSSW = 0;

  // afiqaize get the actual CMSSW that is used
  // the current approach ignores patch and pre releases
  std::string cmssw_ver = edm::getReleaseVersion();
  std::vector<size_t> usc(4, 0);
  for (uint iU = 0; iU < usc.size(); ++iU)
    usc.at(iU) = cmssw_ver.find("_", (iU == 0) ? iU : usc.at(iU - 1) + 1);
#ifdef CMSSW42X
  CMSSW = 10000 * std::atoi(cmssw_ver.substr(usc.at(0) + 1, usc.at(1) - usc.at(0) - 1).c_str()) + 
    100 * std::atoi(cmssw_ver.substr(usc.at(1) + 1, usc.at(2) - usc.at(1) - 1).c_str()) + 
    std::atoi(cmssw_ver.substr(usc.at(2) + 1, (usc.at(3) == std::string::npos) ? cmssw_ver.size() - usc.at(2) : usc.at(3) - usc.at(2) - 1).c_str());
#else
  CMSSW = 10000 * std::stoi(cmssw_ver.substr(usc.at(0) + 1, usc.at(1) - usc.at(0) - 1)) + 
    100 * std::stoi(cmssw_ver.substr(usc.at(1) + 1, usc.at(2) - usc.at(1) - 1)) + 
    std::stoi(cmssw_ver.substr(usc.at(2) + 1, (usc.at(3) == std::string::npos) ? cmssw_ver.size() - usc.at(2) : usc.at(3) - usc.at(2) - 1));
#endif

  if (CMSSW == 0) {
    cout << "*** ALARM ***: NanoAnalyzer: proper CMSSW ifdef flag not set" << endl;
    exit(1);
  }
  else {
    cout << "CMSSW flag set to CMSSW" << CMSSW << endl;   
  }

  // for debug: get and print some info on map
  //std::cout << "max_size = " << hlt_bit.max_size() << std::endl;
  //std::cout << "max_bucket_count = " << hlt_bit.max_bucket_count() << std::endl;
  //std::cout << "max_load_factor = " << hlt_bit.max_load_factor() << std::endl;

  GenPV_x = -999.;
  GenPV_y = -999.;
  GenPV_z = -999.;
  GenPV_recIdx = -1;
  GenPV_chmult = 0;

  Bingo = 0;
  D0vtxx = 0.; D0vtxy = 0.; D0vtxz = 0.;
  D0motx = 0.; D0moty = 0.; D0motz = 0.;  

// The following code deals with the content of overlapping datasets
// auto-detection of dataset from triggers occurring in it

// Events can occur on several datasets. Therefore set everything to true to start with.
// But if no trigger of a given data set fired, the event is definitely *not* from this dataset.
// If no trigger information available (2010 MC) all flags will remain true 
MCdataset = true;
ZeroBiasdataset = true;
MinimumBiasdataset = true;
Commissioningdataset = true;
Mudataset = true;
MuHaddataset = true;
DoubleMudataset = true;
MuEGdataset = true;
Electrondataset = true;
DoubleElectrondataset = true;
Photondataset = true;
MuOniadataset = true;
Charmoniumdataset = true;
MuMonitordataset = true;
EGMonitordataset = true;
Jetdataset = true;
MultiJetdataset = true;
JetMETTauMonitordataset = true;
BTaudataset = true;
BParkingdataset = true;
METFwddataset = true;
datasetisunique = false;
dataset = "unknown";

#ifdef CMSSW7plus
// consumes initialization for Run 2
#ifndef miniAOD
  // AOD collections
  muTkn = consumes<reco::MuonCollection>(edm::InputTag("muons"));
  gmuTkn = consumes<reco::TrackCollection>(edm::InputTag("globalMuons"));
  trkTkn = consumes<reco::TrackCollection>(edm::InputTag("generalTracks"));
  beamTkn = consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));
  //eTkn = consumes<reco::GsfElectronCollection>(edm::InputTag("gsfElectrons"));
  eTkn = consumes<reco::GsfElectronCollection>(edm::InputTag("gedGsfElectrons"));
  photTkn = consumes<reco::PhotonCollection>(edm::InputTag("photons"));
  tauTkn = consumes<reco::PFTauCollection>(edm::InputTag("hpsPFTauProducer"));
  genTkn = consumes<reco::GenParticleCollection>(edm::InputTag("genParticles"));
  pfmetTkn = consumes<reco::PFMETCollection>(edm::InputTag("pfMet"));
  // pfmetTkn = consumes< std::vector<reco::PFMET> >(edm::InputTag("pfMet"));
  calometTkn = consumes<reco::CaloMETCollection>(edm::InputTag("caloMet"));
  // calometTkn = consumes< std::vector<reco::CaloMET> >(edm::InputTag("caloMet"));
  // muCorrmetTkn = consumes< std::vector<reco::CaloMET> >(edm::InputTag("corMetGlobalMuons"));
  // nuha: input tag for conversion collection
  hConvTkn = consumes<reco::ConversionCollection>(edm::InputTag("allConversions"));
  
#ifdef Compatibility
  // official nanoAOD uses primary vertex without beam spot constraint
  primvtxTkn = consumes<reco::VertexCollection>(edm::InputTag("offlinePrimaryVertices"));
#endif
#ifndef Compatibility
  // to get primary vertices with beam spot constraint 
  primvtxTkn = consumes<reco::VertexCollection>(edm::InputTag("offlinePrimaryVerticesWithBS"));
#endif
  pfjetTkn = consumes<reco::PFJetCollection>(edm::InputTag("ak4PFJets"));
  // dEdx (from example M. Soares)
  //  dedxMapStripTag_ = consumes<edm::ValueMap<reco::DeDxData>> (iConfig.getParameter<edm::InputTag>("dedxHarmonic2"));
  //dedxMapPixelTag_ = consumes<edm::ValueMap<reco::DeDxData>> (iConfig.getParameter<edm::InputTag>("dedxPixelHarmonic2"));
  dedxMapStripTag_ = consumes<edm::ValueMap<reco::DeDxData>>(edm::InputTag("dedxHarmonic2"));
  dedxMapPixelTag_ = consumes<edm::ValueMap<reco::DeDxData>>(edm::InputTag("dedxPixelHarmonic2"));
  // dEdx (from example M. Soares)
  //  dedxMapStripTag_(consumes<edm::ValueMap<reco::DeDxData>> (iConfig.getParameter<edm::InputTag>("dedxHarmonic2")));
  //dedxMapPixelTag_(consumes<edm::ValueMap<reco::DeDxData>> (iConfig.getParameter<edm::InputTag>("dedxPixelHarmonic2")));
#endif

#ifdef miniAOD
  // miniAOD collections
  muTkn = consumes<pat::MuonCollection>(edm::InputTag("slimmedMuons"));
  trkTkn = consumes<pat::PackedCandidateCollection>(edm::InputTag("packedPFCandidates"));
  trkTkndisc = consumes<pat::PackedCandidateCollection>(edm::InputTag("packedPFCandidatesDiscarded"));
  trkTknlost = consumes<pat::PackedCandidateCollection>(edm::InputTag("lostTracks"));
  beamTkn = consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));
  eTkn = consumes<pat::ElectronCollection>(edm::InputTag("slimmedElectrons"));
  photTkn = consumes<pat::PhotonCollection>(edm::InputTag("slimmedPhotons"));
  tauTkn = consumes<pat::TauCollection>(edm::InputTag("slimmedTaus"));
  genTkn = consumes<reco::GenParticleCollection>(edm::InputTag("prunedGenParticles"));
  pfjetTkn = consumes<pat::JetCollection>(edm::InputTag("slimmedJets"));
  pfmetTkn = consumes<pat::METCollection>(edm::InputTag("slimmedMETs"));
  //  calo MET actually also available via slimmedMETs
  calometTkn = consumes<pat::METCollection>(edm::InputTag("slimmedMETs"));
  primvtxTkn = consumes<reco::VertexCollection>(edm::InputTag("offlineSlimmedPrimaryVertices"));
#endif

  // trigger and flags for both AOD and miniAOD  *** fix duplication ***
  trigTkn = consumes< edm::TriggerResults>(edm::InputTag("TriggerResults","","HLT"));
  trigEvn = consumes< trigger::TriggerEvent>(edm::InputTag("hltTriggerSummaryAOD", "", "HLT")); //Qun
  if (!hlt_proc.empty())
    trig_tkn = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults", "", hlt_proc));
  if (!custom_flag.empty())
    custom_tkn = consumes<edm::TriggerResults>(custom_tag);

  // CMSSWplus 
#endif

// *** move to proper place ***
// <Qun>
//  TrigObj_id = new vector<Int_t>;
//  TrigObj_pt = new vector<Float_t>;
//  TrigObj_phi = new vector<Float_t>;
//  TrigObj_eta = new vector<Float_t>;
// </Qun>

  //cout << "init" << endl;

} // end of constructor

NanoAnalyzer::~NanoAnalyzer()
{
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

///////////////////////////////////////////////////////////////////////////////
////////////// main analysis loop: method called for each event ///////////////
///////////////////////////////////////////////////////////////////////////////

void
NanoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // using namespace edm;
  // using namespace reco;
  // using namespace std;
  using namespace muon; // for TMOneStationTight
  using namespace trigger; // Qun
  
//////////////////////////////////////////////////////////////////////////////
///////////////////////// load relevant event information ////////////////////
//////////////////////////////////////////////////////////////////////////////

#ifndef miniAOD
  // AOD collections
  //Handle<reco::GenParticleCollection> genParticles;
  //Handle<reco::BeamSpot> beamSpotHandle;
  //Handle<reco::VertexCollection> Primvertex;
  Handle<reco::TrackCollection> tracks;
  Handle<reco::TrackCollection> gmuons;
  Handle<reco::MuonCollection> muons;
  Handle<reco::GsfElectronCollection> electrons;
  Handle<reco::PhotonCollection> photons;
  Handle<reco::PFTauCollection> taus;
  Handle<reco::CaloMETCollection> calomet;
  Handle<reco::PFMETCollection> met;
  Handle<reco::PFJetCollection> jets;
  // dEdx (from example Giacomo Fedi)
  //Handle<reco::DeDxDataValueMap> energyLossHandle;
  //edm:: Handle<edm::ValueMap<reco::DeDxData>> energyLossHandle;
  Handle<DeDxDataValueMap> energyLossHandle;
  // somehow none of these work
  //DeDxDataValueMap &  eloss  = *energyLossHandle;
  //const DeDxDataValueMap &  eloss  = *energyLossHandle;
  //const DeDxDataValueMap &  eloss  = *energyLossHandle.product();
  //const DeDxDataValueMap eloss  = *energyLossHandle;
  // nuha
  Handle<reco::ConversionCollection> hConversions;
#endif

#ifdef miniAOD
  // miniAOD collections
  //Handle<reco::GenParticleCollection> genParticles;
  //Handle<reco::BeamSpot> beamSpotHandle;
  //Handle<reco::VertexCollection> Primvertex;
  Handle<pat::PackedCandidateCollection> tracks;
  Handle<pat::PackedCandidateCollection> discTracks;
  Handle<pat::PackedCandidateCollection> lostTracks;
  // Handle<reco::TrackCollection> gmuons;
  Handle<pat::MuonCollection> muons;
  Handle<pat::ElectronCollection> electrons;
  Handle<pat::PhotonCollection> photons;
  Handle<pat::TauCollection> taus;
  Handle<pat::METCollection> calomet;
  Handle<pat::METCollection> met;
  Handle<pat::JetCollection> jets;  
#endif

  // common AOD and miniAOD collections
  edm::Handle<reco::GenParticleCollection> genParticles;
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  edm::Handle<reco::VertexCollection> Primvertex;
  // for trigger and flags
  edm::Handle<edm::TriggerResults> trigger_handle;
  edm::Handle<edm::TriggerResults> custom_handle;

  /*
#ifndef CMSSW7plus
  // seems not needed/not to be working for CMSSW7 and higher?
  edm::Handle<edm::TriggerResults> triggerResults;
#endif
  */
  // the following is a duplication which still needs to be treated
  //edm::Handle<edm::TriggerResults> triggerResultsHandle_; //Qun
  //edm::Handle<trigger::TriggerEvent> triggerEventHandle_;  //Qun 17-10-19 TriggerObj

  //cout << "hello get event" << endl; 

#ifndef CMSSW7plus
  // for Run 1, use access via getByLabel
  iEvent.getByLabel("muons", muons);
  iEvent.getByLabel("globalMuons", gmuons);
  iEvent.getByLabel("generalTracks", tracks); 
  iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);
  iEvent.getByLabel("gsfElectrons", electrons);
  iEvent.getByLabel("photons", photons);
  iEvent.getByLabel("hpsPFTauProducer", taus);
  iEvent.getByLabel("pfMet", met);
  //iEvent.getByLabel("caloMet", calomet);
  iEvent.getByLabel("met", calomet);
  iEvent.getByLabel("ak5PFJets", jets);
  //iEvent.getByLabel("ak7PFJets", jets);

  // choose primary vertices with/without beam spot
  // see https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideOfflinePrimaryVertexProduction
#ifndef Compatibility
  // for best performance (assumes beam spot is well simulated by MC)
  iEvent.getByLabel("offlinePrimaryVerticesWithBS",Primvertex);
#endif
#ifdef Compatibility
  // for compatibility with official nanoAOD
  iEvent.getByLabel("offlinePrimaryVertices",Primvertex);
#endif
  // dEdx (from example G. Fedi)
  int dedxexist = 1;
  // exists!
  if (!iEvent.getByLabel("dedxHarmonic2", energyLossHandle)) dedxexist=0;
  // does not exist
  //if (!iEvent.getByLabel("dedxPixelHarmonic2", energyLossHandle)) dedxexist=0;

  // trigger and flags (Afiq)
  if (!skiptrigger)
    iEvent.getByLabel(edm::InputTag("TriggerResults", "", hlt_proc), trigger_handle);
  if (!custom_flag.empty())
    iEvent.getByLabel(custom_tag, custom_handle);
  // not CMSSW7plus
#endif

#ifdef CMSSW7plus 
  // for Run 2, use access via getByToken
#ifndef miniAOD
  // AOD collections
  iEvent.getByToken(muTkn, muons);
  iEvent.getByToken(gmuTkn, gmuons);
  iEvent.getByToken(trkTkn, tracks);
  iEvent.getByToken(eTkn, electrons);
  iEvent.getByToken(photTkn, photons);
  iEvent.getByToken(tauTkn, taus);
  iEvent.getByToken(genTkn, genParticles);
  iEvent.getByToken(beamTkn, beamSpotHandle);
  iEvent.getByToken(pfmetTkn, met);
  iEvent.getByToken(calometTkn, calomet);
  iEvent.getByToken(primvtxTkn, Primvertex);
  // iEvent.getByToken(muCorrmetTkn, muCorrmets);
  iEvent.getByToken(pfjetTkn, jets);
  // trigger
  iEvent.getByToken(trigTkn, triggerResultsHandle_);
  iEvent.getByToken(trigEvn, triggerEventHandle_); //Qun
  // nuha
  iEvent.getByToken(hConvTkn, hConversions);
  
  // dEdx (from example M. Soares)
  int stripmap = 1;
  edm::Handle<edm::ValueMap<reco::DeDxData>> dedxStMap;
  if (!iEvent.getByToken(dedxMapStripTag_,dedxStMap)) stripmap=0;  
  int pixmap = 1;
  edm::Handle<edm::ValueMap<reco::DeDxData>> dedxPixMap;
  if (!iEvent.getByToken(dedxMapPixelTag_,dedxPixMap)) pixmap=0;  
#endif

#ifdef miniAOD
  // miniAOD collections
  iEvent.getByToken(muTkn, muons);
  //  iEvent.getByToken(gmuTkn, gmuons);
  iEvent.getByToken(trkTkn, tracks);
  iEvent.getByToken(trkTkndisc, discTracks);
  iEvent.getByToken(trkTknlost, lostTracks);
  iEvent.getByToken(eTkn, electrons);
  iEvent.getByToken(photTkn, photons);
  iEvent.getByToken(tauTkn, taus);
  iEvent.getByToken(genTkn, genParticles);
  iEvent.getByToken(beamTkn, beamSpotHandle);
  iEvent.getByToken(pfmetTkn, met);
  iEvent.getByToken(calometTkn, calomet);
  iEvent.getByToken(primvtxTkn, Primvertex);
  // iEvent.getByToken(muCorrmetTkn, muCorrmets);
  iEvent.getByToken(pfjetTkn, jets);
  // trigger (*** fix duplication ***)
  iEvent.getByToken(trigTkn, triggerResultsHandle_);  
  iEvent.getByToken(trigEvn, triggerEventHandle_); //Qun
#endif

  // common collections
  if (!hlt_proc.empty())
    iEvent.getByToken(trig_tkn, trigger_handle);
  if (!custom_flag.empty())
    iEvent.getByToken(custom_tkn, custom_handle);

  // CMSSW7plus
#endif

  //cout << "hello get run/event info" << endl; 

// *************************************************************
//------------------ get run/event info ------------------------
// *************************************************************  

  run = (iEvent.id()).run();
  event = (iEvent.id()).event();
  luminosityBlock = (iEvent.id()).luminosityBlock();

// *************************************************************
//------------------ get/set JSON quality flag -----------------
// *************************************************************

  GoodLumisection = true;
#ifdef JSONcheck
  // the following automatically excludes MC
  if (run > 132000 && run < 150000) {
    // 2010 data only for the moment  
    //   (if Golden JSON is selected in configuration
    //    or #define JSON is not selected 
    //    this flag will always remain true)
    GoodLumisection = providesGoodLumisection(iEvent);
  }
#endif

#ifdef CMSSW42X
              // skip trigger for 2010 MC  (not available)
              if (run == 1) {
                if (!skiptrigger)    // on first event only 
                   cout << "Trigger information will be skipped" << endl;
                skiptrigger = true;
              }
#endif
	      // cout << "hello skip trigger " << run << " " << skiptrigger << endl; 
	      if (!skiptrigger) {
// *************************************************************
//------------------ get and check trigger info ----------------
// *************************************************************
//      sets hlt_bit.at pointers, i.e. whether trigger has fired or not

  // afiqaize (still need to make consistent with the following) 
  //cout << "hello get trigger" << endl; 

  if (!hlt_proc.empty()) {
    //cout << "hello hlt_proc.empty" << endl; 
    for (unsigned iP = 0; iP < trigger_handle->size(); ++iP) {
      // get path name with version stripped off
      const std::string path = remove_version( iEvent.triggerNames(*trigger_handle).triggerName(iP) );
      // check whether path occurs in container (0 or 1)
      if (hlt_bit.count(path)) {
        // for debug 
        //cout << "hello accepted trigger path " << iP << " " << path << endl;
        // trigger_handle->accept returns 0 or 1,
        // .at refers to the value of the second element of the pair 
        hlt_bit.at(path) = trigger_handle->accept(iP);
      }
      else if (path.substr(0, 3) == "HLT")
        std::cout << "WARNING: HLT path " << path << " is not present in the bitmap, when it should be!! Check if the logic is "
          "properly implemented!" << std::endl;
      // for debug
      //else 
      //  cout << "hello *rejected* trigger path " << iP << " " << path << endl;
    }
  } // hlt_proc_empty
	      } // skiptrigger

              // initialize all dataset trigger flags to true
              // (for 2010 MC)
              // will be superseeded if trigger menue is present
              ZeroBiasTrig   = true;
              MinimumBiasTrig= true;
              MuTrig         = true;
              MuMonitorTrig  = true;
              MuOniaTrig     = true;
              ElectronTrig   = true;
              EGMonitorTrig  = true;
              BParkingTrig   = true;
              // ... add more!
              // initialize global trigger content flags to none set
              ZeroBiasFlag   = -10;
              MinBiasFlag    = -10;
              MinBiasMult    = -99;
              MuThresh       = -49;
              MuL1Thresh     = -49;
              MuL2Thresh     = -49;
              IsoMuThresh    = -49;
              DoubleMuThresh = -49;
              JpsiThresh     = -49;
              MuHadFlag      = -10;
              MuEGFlag       = -10;
              ElectronThresh = -49;
              DoubleElectronThresh = -49;
              PhotonThresh   = -49;
              JetThresh      = -199;
              DiJetThresh    = -199; 
              TriJetThresh   = -199; 
              QuadJetThresh  = -199; 
              HTThresh       = -199;
              BThresh        = -199;
              METThresh      = -199;

              if (skiptrigger) {
                // declare ZeroBias Trigger to be set (only trigger)
                ZeroBiasTrig = true;
                ZeroBiasFlag = 1;
                MCdataset = true;
                // all other flags remain unchanged ...
              }

	      // Qun below  Trigger Object
	      nTrigObj = 0;                   // number of stored Trigger Obj
	      TrigObj_id.clear();
	      TrigObj_filterBits.clear();
	      TrigObj_pt.clear();
	      TrigObj_eta.clear();
	      TrigObj_phi.clear();
	      //TrigObj_id->clear();
	      //TrigObj_pt->clear();
	      //TrigObj_eta->clear();
	      //TrigObj_phi->clear();
	      // Qun above  Trigger Object

              // if not to be skipped
              if (!skiptrigger) {
#ifndef CMSSW7plus
                // trigger, do not move up! 
                edm::InputTag trigResultsTag("TriggerResults","","HLT");
                iEvent.getByLabel(trigResultsTag,triggerResultsHandle_);
		//below  Qun 17-10-19 TriggerObj 
		//edm::InputTag triggerEventTag_;  
		//edm::InputTag("hltTriggerSummaryAOD", "", "HLT");//
		edm::InputTag triggerEventTag_("hltTriggerSummaryAOD", "", "HLT");  
		iEvent.getByLabel(triggerEventTag_,triggerEventHandle_);
		if (!triggerEventHandle_.isValid()) {
		  cout << "HLTEventAnalyzerAOD::analyze: Error in getting TriggerEvent product from Event!" << endl;
		  return;		  
		}
		if (triggerEventHandle_.isValid()) {
                  /*   reactivate this later!
		  cout << "Used Processname: " << triggerEventHandle_->usedProcessName() << endl;
		  const size_type nC(triggerEventHandle_->sizeCollections());
		  //cout << "Number of packed Collections: " << nC << endl;
		  //cout << "The Collections: #, tag, 1-past-end index" << endl;
		  for (size_type iC=0; iC!=nC; ++iC) {
     		    cout << iC << " "
	    	    << triggerEventHandle_->collectionTag(iC).encode() << " "
	            << triggerEventHandle_->collectionKey(iC) << endl;
     		  }
                  */
		  const size_type nO(triggerEventHandle_->sizeObjects());
		  //cout << "Number of TriggerObjects: " << nO << endl;
		  //cout << "The TriggerObjects: #, id, pt, eta, phi, mass" << endl;
		  const TriggerObjectCollection& TOC(triggerEventHandle_->getObjects());
		  // wrong cout << "Number of TriggerObjects ID: " << TOC.pt.size() << endl;
		  size_type nTO_QQ=0;
		  for (size_type iO=0; iO!=nO; ++iO) {
		    const TriggerObject& TO(TOC[iO]);
		    //cout << iO << " " << TO.id() << " " << TO.pt() << " " << TO.eta() << " " << TO.phi() << " " << TO.mass() << endl;
                    // Ids indicate kind of trigger object, e.g. +-13 = muon
		    if (TO.id()!=0 && fabs(TO.id())<60) {
		      //TrigObj_id[nTO_QQ]=TO.id();
		      //TrigObj_pt[nTO_QQ]=TO.pt();
		      //TrigObj_eta[nTO_QQ]=TO.eta();
		      //TrigObj_phi[nTO_QQ]=TO.phi();
                      if (nTO_QQ < nReserve_TrigObj) {
		        TrigObj_id.push_back(TO.id());
		        TrigObj_filterBits.push_back(-1);
		        TrigObj_pt.push_back(TO.pt());
		        TrigObj_eta.push_back(TO.eta());
		        TrigObj_phi.push_back(TO.phi());
		        //TrigObj_id->push_back(TO.id());
		        //TrigObj_pt->push_back(TO.pt());
		        //TrigObj_eta->push_back(TO.eta());
		        //TrigObj_phi->push_back(TO.phi());
		        nTO_QQ++; 
		        //cout << "The TriggerObjects: #, id, pt, eta, phi, mass of TrigObj " << TrigObj_id[nTO_QQ-1] <<  " " << TrigObj_pt[nTO_QQ-1] <<  " " <<TrigObj_eta[nTO_QQ-1] <<  " " << TrigObj_phi[nTO_QQ-1] << endl;
		      } 
                      else cout << "*** nReserve_TrigObj exceeded" << endl;    
		    }    
		    //cout << "The TriggerObjects: #, id, pt, eta, phi, mass of TrigObj" << TrigObj_id[nTO_QQ-1] <<  " " << TrigObj_pt[nTO_QQ-1] <<  " " <<TrigObj_eta[nTO_QQ-1] <<  " " << TrigObj_phi[nTO_QQ-1] << endl;
		    //cout << "The TriggerObjects: #, id, pt, eta, phi, mass of TrigObj" << TrigObj_id[iO] <<  " " << TrigObj_pt[iO] <<  " " <<TrigObj_eta[iO] <<  " " << TrigObj_phi[iO] << endl;
		  }
		  nTrigObj=nTO_QQ;
		  //cout << "Number of TriggerObjects QQ with ID!=0 : " << nTO_QQ << endl;
                  /*
                  // for debug:
                  if (nTrigObj>0) {
		    const size_type nC(triggerEventHandle_->sizeCollections());
		    //cout << "Number of packed Collections: " << nC << endl;
		    //cout << "The Collections: #, tag, 1-past-end index" << endl;
            	    for (size_type iC=0; iC!=nC; ++iC) {
     		      cout << iC << " "
	    	      << triggerEventHandle_->collectionTag(iC).encode() << " "
	              << triggerEventHandle_->collectionKey(iC) << endl;
     		    }
                  }
                  */
                  // some relevant collection Tags,                   Key
                  // Electrons, 2011 MC:
                  //    hltL1IsoRecoEcalCandidate::HLT                1
                  //    hltL1NonIsoRecoEcalCandidate::HLT             1
                  //    hltRecoEcalSuperClusterActivityCandidate::HLT 0,2 or 4
                  //    hltPixelMatchElectronsL1Iso::HLT              3,5 or 6
                  //    hltPixelMatchElectronsL1NonIso::HLT           3,5 or 6
                  //    hltPixelMatch3HitElectronsL1Iso::HLT          4 or 5
                  //    hltPixelMatch3HitElectronsL1NonIso::HLT       4 or 5     
                  // Muons, 2011 MC:
                  //    hltL2MuonCandidates::HLT                      1 or 2
                  //    hltL3MuonCandidates::HLT                      2 or 3
                  //    hltMuTrackJpsiPixelTrackCands::HLT            2

		  //Qconst size_type nF(triggerEventHandle_->sizeFilters());
		  //Qcout << "Number of TriggerFilters: " << nF << endl;
		  //Qcout << "The Filters: #, tag, #ids/#keys, the id/key pairs" << endl;
		  //Qfor (size_type iF=0; iF!=nF; ++iF) {
		  //Q  const Vids& VIDS (triggerEventHandle_->filterIds(iF));
		  //Q  const Keys& KEYS(triggerEventHandle_->filterKeys(iF));
		  //Q  const size_type nI(VIDS.size());
		  //Q  const size_type nK(KEYS.size());
		  //Q  cout << iF << " " << triggerEventHandle_->filterTag(iF).encode()
		  //Q       << " " << nI << "/" << nK
		  //Q       << " the pairs: ";
		  //Q  const size_type n(max(nI,nK));
		  //Q  for (size_type i=0; i!=n; ++i) {
		  //Q    cout << " " << VIDS[i] << "/" << KEYS[i];
		  //Q  }
		  //Q  cout << endl;
		  //Q  assert (nI==nK);
		  //Q}
		}
	        //const TriggerObjectCollection& TOC(triggerEventHandle_->getObjects());

		//above  Qun
#endif
                const edm::TriggerNames& trigNames = iEvent.triggerNames(*triggerResultsHandle_); 
		//cout << "QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ" << endl; //Qun
                //cout << "test hlt QQQQQQQQQQQQQQQQ " << hltConfig_.size() <<endl; //Qun

                // check whether trigger info was obtained successfully
                if (!triggerResultsHandle_.isValid()) {
                  cout << "Nano::analyze: Error in getting TriggerResults product from Event!" << endl;
                  // should always be available for data, but not necessarily 
                  // for e.g. 2010 MC
                  if (run > 1) {
                    // stop the program
                    exit (1);
                  }
                  else {
                    // do not ask for trigger info from here onwards
                    skiptrigger = true;
                  }
                }

		// below Qun Trigger Object
		using namespace reco;
		using namespace trigger;

		assert(triggerResultsHandle_->size()==hltConfig_.size()); //Qun
		//cout << " triggerResultsHandle_->size()==hltConfig_.size() QQQQ" << endl;
		//const unsigned int n(hltConfig_.size());
		//Get the trigger index for the current trigger
		const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName_));
		//cout << "before QQQQQ n: " << n << " QQQQQ triggerIndex: " << triggerIndex << endl;
		//check that the trigger in the event and in the configuration agree
		assert(triggerIndex==iEvent.triggerNames(*triggerResultsHandle_).triggerIndex(triggerName_));
		//cout << "after  QQQQQ n: " << n << " QQQQQ triggerIndex: " << triggerIndex << endl;
		//Qif (triggerIndex>=n) {
		//Q  cout << "HLTEventAnalyzerAOD::analyzeTrigger: path "
		//Q       << triggerName_ << " - not found!" << endl;
		//Q  return;
		//Q}
		//Qelse {cout << "HLTEventAnalyzerAOD QQQQQ " << triggerName_ << endl;}
		//Qif (triggerName_=="@") {
		//Q  const unsigned int n(hltConfig_.size());
		//Q  cout << "analyzerTrigger testA  " << n << endl;
		//Q  for (unsigned int i=0; i!=n; ++i) {
		//Q    analyzeTrigger(iEvent,iSetup,hltConfig_.triggerName(i));
		//Q  //cout << "analyzerTrigger testB  " << i << endl;
		//Q  }
		//Q} else {
		//Q  analyzeTrigger(iEvent,iSetup,triggerName_);
		//Q  //cout << "analyzerTrigger testC  " << n << endl;
		//Q}
		// Qun above Trgger Object

               // if not to be skipped (duplicate on purpose!)
               if (!skiptrigger) {

//  Here: check for "good" minimum bias, jet, or muon or electron trigger
//  separate calls from subsequent if statement in order not to make them order-dependent
//  always check all triggers on all datasets
                GoodMinBiasTrigger = providesGoodMinBiasTrigger(iEvent);
                GoodJetTrigger     = providesGoodJetTrigger(iEvent);
                GoodMuTrigger      = providesGoodMuTrigger(iEvent);
                GoodETrigger       = providesGoodETrigger(iEvent);

                //// for debug, dump all triggers of menu
                //      for (unsigned i = 0; i<trigNames.size(); i++) {
                //      // dump only accepted triggers 
		//	//                        if (triggerResultsHandle_->accept(i)==1){
                //          std::cout<<" Trigger name "<<trigNames.triggerName(i)<<", Accepted = ";
                //          std::cout<<triggerResultsHandle_->accept(i)<<std::endl;
		//	  //}
                //      }


// datasets from which triggers are being treated are 
//   (* = still to be implemented)
//
//   Commissioning10: (7 TeV and 900 GeV): Zerobias and MinimumBias 
//
//      2010A       2010B         2011A             2012B/C
//      ZeroBias 
//      MinimumBias MinimumBias   MinimumBias       MinimumBias*     <-- next
//   Commissioning* Commissioning                   Commissioning*
//      Mu          Mu          ( SingleMu        ( SingleMu
//                              ( DoubleMu        ( DoubleMuParked
//                              ( MuHad           ( MuHad*
//      next -->                ( MuEG            ( MuEG*
//      MuMonitor*  MuMonitor
//      MuOnia*     MuOnia        MuOnia* <-- next  MuoniaParked*
//      EG*         Electron    ( SingleElectron* ( SingleElectron*
//                              ( DoubleElectron  ( DoubleElectron*  
//		                ( ElectronHad*    ( ElectronHad*
//      EGmonitor*  EGMonitor
//          next--> Photon*     ( Photon*        (( SinglePhoton*
//                                               (( DoublePhoton*
//                                               (( DoublePhotonHighPt*
//                              ( PhotonHad*      ( PhotonHad*
//      JeTMETTau*
//      BTau*       BTau        ( BTag*          (( BTag*
//                                               (( BJetPlusX*
//                              ( Tau*            ( TauParked*
//                              ( TauPlusX*       ( TauPlusX
//      JetMET*     Jet         ( Jet            (( Jet*
//                                               (( JetHT*
//                                               (( JetMon* 
//                              ( HT*             ( HTMHTParked*
//JetMETTauMonitor* JetMETTauMonitor
//                  Multijet      MultiJet*
//                  METFwd*     ( MET*            ( MET*
//                              ( METBTag* 
//                                                  HcalNZS*
//                                                  NoBPTX*
//                                                  VBF1Parked*
//
//     also: 2015E MinimumBias (5 TeV)
//           2016  ZeroBias (13 TeV)
//           2018  B Parking* (13 TeV)
//
                  // sharpen/check dataset info if not yet unique
                  if (!datasetisunique) {
                    // dataset is not yet unique 
		    std::cout<<" dataset is not yet unique, Triggers:"<<std::endl;
                    // dump the trigger info
                      for (unsigned i = 0; i<trigNames.size(); ++i) {
                      // dump only accepted triggers 
                        if (triggerResultsHandle_->accept(i)==1){
                          std::cout<<" Trigger name "<<trigNames.triggerName(i)<<", Accepted = ";
                          std::cout<<triggerResultsHandle_->accept(i)<<std::endl;
	                }
                      }
                    if (isData) MCdataset = false;
                    if (!ZeroBiasTrig) ZeroBiasdataset = false;
                    if (!MinimumBiasTrig) MinimumBiasdataset = false;
                    if (!CommissioningTrig) Commissioningdataset = false;
                    if (!MuTrig) Mudataset = false;
                    if (!MuHadTrig) MuHaddataset = false;
                    if (!DoubleMuTrig) DoubleMudataset = false;
                    if (!MuEGTrig) MuEGdataset = false;
                    if (!ElectronTrig) Electrondataset = false;
                    if (!DoubleElectronTrig) DoubleElectrondataset = false;
                    if (!PhotonTrig) Photondataset = false;
                    if (!MuOniaTrig) MuOniadataset = false;
                    if (!CharmoniumTrig) Charmoniumdataset = false;
		    if (!MuMonitorTrig) MuMonitordataset = false;
                    if (!EGMonitorTrig) EGMonitordataset = false;
                    if (!JetTrig) Jetdataset = false;
                    if (!MultiJetTrig) MultiJetdataset = false;
                    if (!JetMETTauMonitorTrig) JetMETTauMonitordataset = false;
                    if (!BTauTrig) BTaudataset = false;
                    if (!BParkingTrig) BParkingdataset = false;
                    if (!METFwdTrig) METFwddataset = false;
                    int ndataset = 0;
                    if (MCdataset) ++ndataset;
                    if (ZeroBiasdataset) ++ndataset;
                    if (MinimumBiasdataset) ++ndataset;
                    if (Commissioningdataset) ++ndataset;
                    if (Mudataset) ++ndataset;
                    if (MuHaddataset) ++ndataset;
                    if (DoubleMudataset) ++ndataset;
                    if (MuEGdataset) ++ndataset;
                    if (Electrondataset) ++ndataset;
                    if (DoubleElectrondataset) ++ndataset;
                    if (Photondataset) ++ndataset;
                    if (MuOniadataset) ++ndataset;
                    if (Charmoniumdataset) ++ndataset;
                    if (MuMonitordataset) ++ndataset;
                    if (EGMonitordataset) ++ndataset;
                    if (Jetdataset) ++ndataset;
                    if (MultiJetdataset) ++ndataset;
                    if (JetMETTauMonitordataset) ++ndataset;
                    if (BTaudataset) ++ndataset;
                    if (BParkingdataset) ++ndataset;
                    if (METFwddataset) ++ndataset;
		    std::cout<<" dataset not yet unique, "<<ndataset<<" choices."<<std::endl;
		    std::cout<<" MC= "<<MCdataset<<std::endl;
		    std::cout<<" MinimumBias= "<<MinimumBiasdataset<<std::endl;
		    std::cout<<" ZeroBias= "<<ZeroBiasdataset<<std::endl;
		    std::cout<<" Commissioning= "<<Commissioningdataset<<std::endl;
		    std::cout<<" Mu= "<<Mudataset<<std::endl;
		    std::cout<<" MuHad= "<<MuHaddataset<<std::endl;
		    std::cout<<" DoubleMu= "<<DoubleMudataset<<std::endl;
		    std::cout<<" MuEG= "<<MuEGdataset<<std::endl;
		    std::cout<<" Electron= "<<Electrondataset<<std::endl;
		    std::cout<<" DoubleElectron= "<<DoubleElectrondataset<<std::endl;
		    std::cout<<" Photon= "<<Photondataset<<std::endl;
		    std::cout<<" MuOnia= "<<MuOniadataset<<std::endl;
		    std::cout<<" Charmonium= "<<Charmoniumdataset<<std::endl;
		    std::cout<<" MuMonitor= "<<MuMonitordataset<<std::endl;
		    std::cout<<" EGMonitor= "<<EGMonitordataset<<std::endl;
		    std::cout<<" Jet= "<<Jetdataset<<std::endl;
		    std::cout<<" MultiJet= "<<MultiJetdataset<<std::endl;
		    std::cout<<" JetMETTauMonitor= "<<JetMETTauMonitordataset<<std::endl;
		    std::cout<<" BTau= "<<BTaudataset<<std::endl;
		    std::cout<<" BParking= "<<BParkingdataset<<std::endl;
		    std::cout<<" METFwd= "<<METFwddataset<<std::endl;

                    if (MCdataset) {
                      datasetisunique = true;
                      dataset = "MC";
         	      std::cout<<" dataset now unique: "<<dataset<<std::endl;
                    }
                    else if (ndataset==1) {
		      // dataset is now unique!
                      datasetisunique = true;
                      if (ZeroBiasdataset) dataset = "ZeroBias";
                      if (MinimumBiasdataset) dataset = "MinimumBias";
                      if (Commissioningdataset) dataset = "Commissioning";
                      if (Mudataset) dataset = "Mu";
                      if (MuHaddataset) dataset = "MuHad";
                      if (DoubleMudataset) dataset = "DoubleMu";
                      if (MuEGdataset) dataset = "MuEG";
                      if (Electrondataset) dataset = "Electron";
                      if (DoubleElectrondataset) dataset = "DoubleElectron";
                      if (Photondataset) dataset = "Photon";
                      if (MuOniadataset) dataset = "MuOnia";
                      if (Charmoniumdataset) dataset = "Charmonium";
                      if (MuMonitordataset) dataset = "MuMonitor";
                      if (EGMonitordataset) dataset = "EGMonitor";
                      if (Jetdataset) dataset = "Jet";
                      if (MultiJetdataset) dataset = "MultiJet";
                      if (JetMETTauMonitordataset) dataset = "JetMETTauMonitor";
                      if (BTaudataset) dataset = "BTau";
                      if (BParkingdataset) dataset = "BParking";
                      if (METFwddataset) dataset = "METFwd";
         	      std::cout<<" dataset now unique: "<<dataset<<std::endl;
                    }
#ifndef trigcheckabort
                    else if (ndataset==0) {
                      std::cout<<" ******* no candidate dataset recognized for these triggers ******"<<std::endl;
                      // set dataset unique although it is not
                      datasetisunique = true;
                    }
#endif
#ifdef trigcheckabort
                    else if (ndataset==0) {
  	              // trigger info doesn't match any known dataset, dump trigger and abort
                      std::cout<<" ******* no candidate dataset recognized for these triggers ******"<<std::endl;
                      for (unsigned i = 0; i<trigNames.size(); ++i) {
                      // dump only accepted triggers 
                        if (triggerResultsHandle_->accept(i)==1){
                          std::cout<<" Trigger name "<<trigNames.triggerName(i)<<", Accepted = ";
                          std::cout<<triggerResultsHandle_->accept(i)<<std::endl;
	                }
                      }
                      // abort the job. This should never happen if trigger information from respective dataset
                      // is treated properly
                      // (trigger so far implemented for parts of Run 1 only)
                      exit(1);
                    }
#endif
                  } // !datasetisunique
		  else {
                    // data set is already unique. check consistency of trigger and dataset 
                    if ( dataset == "MC"       // no trigger requirement!
                     || (dataset == "ZeroBias" && ZeroBiasTrig)
                     || (dataset == "MinimumBias" && MinimumBiasTrig)
		     || (dataset == "Mu" && MuTrig)
		     || (dataset == "MuHad" && MuHadTrig)
		     || (dataset == "DoubleMu" && DoubleMuTrig)
		     || (dataset == "MuEG" && MuEGTrig)
		     || (dataset == "Electron" && ElectronTrig)
		     || (dataset == "DoubleElectron" && DoubleElectronTrig)
		     || (dataset == "Photon" && PhotonTrig)
		     || (dataset == "MuOnia" && MuOniaTrig)
		     || (dataset == "Charmonium" && CharmoniumTrig)
		     || (dataset == "MuMonitor" && MuMonitorTrig)
		     || (dataset == "EGMonitor" && EGMonitorTrig)
		     || (dataset == "Jet" && JetTrig)
		     || (dataset == "MultiJet" && MultiJetTrig)
		     || (dataset == "JetMETTauMonitor" && JetMETTauMonitorTrig)
		     || (dataset == "BTau" && BTauTrig)
		     || (dataset == "BParking" && BParkingTrig)
		     || (dataset == "METFwd" && METFwdTrig)
	             || (dataset == "Commissioning" && CommissioningTrig)) {} // do nothing
                    else {
  	              // trigger and dataset are not consistent, dump trigger and abort
                      if (firstevent == true) { 
		        std::cout<<" ****** Trigger not consistent with dataset "<<dataset<<" ****** "<<std::endl;
		        std::cout<<ZeroBiasTrig<<MinimumBiasTrig<<MuTrig<<MuHadTrig<<DoubleMuTrig<<MuEGTrig<<ElectronTrig<<DoubleElectronTrig<<PhotonTrig<<MuOniaTrig<<CharmoniumTrig<<MuMonitorTrig<<EGMonitorTrig<<JetTrig<<MultiJetTrig<<JetMETTauMonitorTrig<<BTauTrig<<BParkingTrig<<METFwdTrig<<CommissioningTrig<<std::endl;
                        for (unsigned i = 0; i<trigNames.size(); ++i) {
                        // dump only accepted triggers 
                          if (triggerResultsHandle_->accept(i)==1){
                            std::cout<<" Trigger name "<<trigNames.triggerName(i)<<", Accepted = ";
                            std::cout<<triggerResultsHandle_->accept(i)<<std::endl;
	                  }
                        }
                        firstevent = false;
#ifdef trigcheckabort
                        // abort the job. This should never happen if trigger information from respective dataset
                        // is treated properly
                        // trigger so far implemented for Run 1 only
                        exit(1);
#endif
                      }
                    } // dataset
                  } // datasetisunique  
                 } // skiptrigger
                } // skiptrigger


/////////////////////////////////////////////////////////////////////////////
////////////// and now proceed to analysis of physics content ///////////////
/////////////////////////////////////////////////////////////////////////////

	      //cout << "hello physics" << endl; 

//////////////////////////////////////////////////////////////////////////////
////////////////////////////// Gen Particle Start ////////////////////////////
//////////////////////////////////////////////////////////////////////////////

  GenPart_pt.clear();
  GenPart_eta.clear();
  GenPart_phi.clear();
  GenPart_mass.clear();
  GenPart_pdgId.clear();
  GenPart_status.clear();
  GenPart_statusFlags.clear();
  GenPart_genPartIdxMother.clear();

  GenPart_Id.clear();
  GenPart_isNano.clear();
  GenPart_parpdgId.clear();
  GenPart_sparpdgId.clear();
  GenPart_numberOfDaughters.clear();
  GenPart_nstchgdaug.clear();
  GenPart_vx.clear();
  GenPart_vy.clear();
  GenPart_vz.clear();
  GenPart_mvx.clear();
  GenPart_mvy.clear();
  GenPart_mvz.clear();
  GenPart_recIdx.clear();

  TLorentzVector p4Kt, p4pit, p4D0t;
  p4Kt.SetPtEtaPhiE(0., 0., 0., 0.);
  p4pit.SetPtEtaPhiE(0., 0., 0., 0.);
  p4D0t.SetPtEtaPhiE(0., 0., 0., 0.);

  int nmusim =0;
  int idmusim[maxnmusim];
  float chgmusim[maxnmusim];
  float ptmusim[maxnmusim];
  float etamusim[maxnmusim];
  float phimusim[maxnmusim];

  int nD0sim =0;
  int idD0sim[maxnD0sim];
  float ptD0sim[maxnD0sim];
  float etaD0sim[maxnD0sim];
  float phiD0sim[maxnD0sim];
  int pdgIdD0t1sim[maxnD0sim];
  int chgD0t1sim[maxnD0sim];
  float ptD0t1sim[maxnD0sim];
  float etaD0t1sim[maxnD0sim];
  float phiD0t1sim[maxnD0sim];
  int pdgIdD0t2sim[maxnD0sim];
  int chgD0t2sim[maxnD0sim];
  float ptD0t2sim[maxnD0sim];
  float etaD0t2sim[maxnD0sim];
  float phiD0t2sim[maxnD0sim];

  int nDstarsim =0;
  int idDstarsim[maxnDstarsim];
  float ptDstarsim[maxnDstarsim];
  float etaDstarsim[maxnDstarsim];
  float phiDstarsim[maxnDstarsim];

  // clear flag for printout of "good" D0 candidates for this event
  Bingo = 0;
  // clear GenPV variables for each event; somehow needed for GenPV_chmult?
  GenPV_x = -999.;
  GenPV_y = -999.;
  GenPV_z = -999.;
  GenPV_recIdx = -1;
  GenPV_chmult = 0;

  if (!isData) {

#ifndef CMSSW7plus   
    iEvent.getByLabel("genParticles", genParticles);
#endif

    bool vertexfilled = false;
    for (unsigned int ee = 0; ee < genParticles->size(); ++ee) {
    
      const GenParticle & genp = (*genParticles)[ee];

      // store true vertex = vertex of first nonzero entry
      if (!vertexfilled && genp.vz() != 0){
	GenPV_x = genp.vx();
	GenPV_y = genp.vy();
	GenPV_z = genp.vz();
        GenPV_recIdx = -1;  // will be filled later
        vertexfilled = true;
      }

      // count charged particle multiplicity
      // somehow doesn't work?
      if (genp.status()==1 && genp.charge()!=0) ++GenPV_chmult; 
      
      // add pt, eta, phi and mass distributions for each gen particle
      // store b and c quarks, D0, D*+, and muons
      // *** needs to be changed to official nanoAOD genparticle filter ***
      if ( (std::abs(genp.pdgId()) == 5) || (std::abs(genp.pdgId()) == 4) || (std::abs(genp.pdgId()) == 421) || (std::abs(genp.pdgId()) == 413) || (std::abs(genp.pdgId()) == 13) ){

        // all these particle types are also stored in standard nanoAOD
        bool isNano = true;

        // get pointer to mother (reco::Candidate type!):
        const Candidate * mom = genp.mother();
	// and store type and index 
        int parpdgId = mom->pdgId();

        int sparpdgId = mom->pdgId();
        // if nonstable should add further loops here 
        // until stable parent is found ... not yet implemented

        // get no. of daughters:
        size_t n = genp.numberOfDaughters();
        int nstchgdaug = 0;

        // initialize vertex variables for this event
	D0vtxx=0, D0vtxy=0, D0vtxz=0; 
        // loop over daughters:
        for(size_t j = 0; j < n; ++j) {

          // get pointer d to a daughter  (reco::Candidate type!):
          const Candidate * d = genp.daughter( j );
  
          // fill position of decay vertex (creation vertex of first daughter)
          if (j==0) {
            D0vtxx = d->vx();
            D0vtxy = d->vy();
            D0vtxz = d->vz();
          }

          // check daughter's stability (status 1) and charge, 
          // and increment counter:
          if (d->status() == 1 && d->charge() != 0) ++nstchgdaug; 
        }
        // creation point of D0 is decay point of mother
        D0motx = genp.vx();
        D0moty = genp.vy();
        D0motz = genp.vz();

        // make list of interesting candidates
	// print out interesting candidates

#ifdef Bingo
        // beauty
        if (genp.pdgId()==5) {cout << "Bingo beauty" << endl;}
#endif

        // muons
        if (abs(genp.pdgId())==13 && genp.status()==1 && nmusim < maxnmusim) {
          idmusim[nmusim]=ee; 
          chgmusim[nmusim]=genp.charge();
          ptmusim[nmusim]=genp.pt();
          etamusim[nmusim]=genp.eta();
          phimusim[nmusim]=genp.phi();
	  ++nmusim;
        }
        if (nmusim >= maxnmusim) cout << "!!! maxnmusim exceeded !!!" << endl;

        // D0 mesons decaying to Kpi, KK, or pipi
        if (abs(genp.pdgId())==421 && n==2 && nstchgdaug ==2 && nD0sim < maxnD0sim) {
          idD0sim[nD0sim]=ee; 
          ptD0sim[nD0sim]=genp.pt();
          etaD0sim[nD0sim]=genp.eta();
          phiD0sim[nD0sim]=genp.phi(); 
          pdgIdD0t1sim[nD0sim]=(genp.daughter(0))->pdgId();
          chgD0t1sim[nD0sim]=(genp.daughter(0))->charge();
          ptD0t1sim[nD0sim]=(genp.daughter(0))->pt();
          etaD0t1sim[nD0sim]=(genp.daughter(0))->eta();
          phiD0t1sim[nD0sim]=(genp.daughter(0))->phi(); 
          pdgIdD0t2sim[nD0sim]=(genp.daughter(1))->pdgId();
          chgD0t2sim[nD0sim]=(genp.daughter(1))->charge();
          ptD0t2sim[nD0sim]=(genp.daughter(1))->pt();
          etaD0t2sim[nD0sim]=(genp.daughter(1))->eta();
          phiD0t2sim[nD0sim]=(genp.daughter(1))->phi(); 
	  ++nD0sim;
        } 
        if (nD0sim >= maxnD0sim) cout << "!!! maxnD0sim exceeded !!!" << endl;

        // all charged Dstar mesons
        if (abs(genp.pdgId())==413 && nDstarsim < maxnDstarsim) {
          idDstarsim[nDstarsim]=ee; 
          ptDstarsim[nDstarsim]=genp.pt();
          etaDstarsim[nDstarsim]=genp.eta();
          phiDstarsim[nDstarsim]=genp.phi();
          ++nDstarsim;
        } 
        if (nDstarsim >= maxnDstarsim) cout << "!!! maxnDstarsim exceeded !!!" << endl;

#ifdef Bingo
        if (std::abs(genp.pdgId()) >100 && std::abs(genp.eta())<2.5) {
          cout << "event " << event << " pdgID " << genp.pdgId() << " pt " << genp.pt() << " eta " << genp.eta() << " phi " << genp.phi() << " n " << n << " " << nstchgdaug << " par " << parpdgId << endl;
          if (n == 2 && nstchgdaug == 2) { 
            Bingo = 1; 
            cout << "*** Bingo D0 ***" << endl; 
            if (std::abs(parpdgId) == 413) {
              cout << "*** Bingo D* ***" << endl;
            }
            // print vertex position of D0 and mother
            cout << "D0 vertex " << D0vtxx << " " << D0vtxy << " " << D0vtxz << " mother " << D0motx << " " << D0moty << " " << D0motz << endl;

            // loop over daughters:
            for(size_t j = 0; j < n; ++j) {

              // get pointer d to a daughter  (reco::Candidate type!):
              const Candidate * d = genp.daughter( j );

              // print information
              cout << " pdgID " << d->pdgId() << " pt " << d->pt() << " eta " << d->eta() << " phi " << d->phi() << endl;
              if (abs(d->pdgId()) == 321) p4Kt.SetPtEtaPhiM(d->pt(), d->eta(), d->phi(), Kmass);
              if (abs(d->pdgId()) == 211) p4pit.SetPtEtaPhiM(d->pt(), d->eta(), d->phi(), pimass);

	    } // end for j loop
            p4D0t = p4Kt + p4pit;
            cout << "true D0 mass " << p4D0t.M() << endl;
          } // end if
        } // end interesting
#endif

	if (GenPart_pt.size() < nReserve_GenPart) {
	  GenPart_pt.push_back(genp.pt());
	  GenPart_eta.push_back(genp.eta());
	  GenPart_phi.push_back(genp.phi());
	  GenPart_mass.push_back(genp.mass());
	  GenPart_pdgId.push_back(genp.pdgId());
          // status = 1: final state stable particles
          // status = 2: hadron level, unstable
          // status = 3: parton level
	  GenPart_status.push_back(genp.status());
          GenPart_statusFlags.push_back(9999);
          // The following does not work, do something else
	  //          GenPart_genPartIdxMother.push_back(genp.mother());
          // Rather loop over previous entries in *this* list, find parent
          // (if stored), and store index of that! 
	  GenPart_genPartIdxMother.push_back(-999);
          GenPart_Id.push_back(ee);
	  GenPart_isNano.push_back(isNano);
	  GenPart_parpdgId.push_back(parpdgId);
	  GenPart_sparpdgId.push_back(sparpdgId);
	  GenPart_numberOfDaughters.push_back(genp.numberOfDaughters());
	  GenPart_nstchgdaug.push_back(nstchgdaug);
          GenPart_vx.push_back(D0vtxx);
          GenPart_vy.push_back(D0vtxy);
          GenPart_vz.push_back(D0vtxz);
          GenPart_mvx.push_back(D0motx);
          GenPart_mvy.push_back(D0moty);
          GenPart_mvz.push_back(D0motz);
          GenPart_recIdx.push_back(-1);
	}
	else{cout << "WARNING!!!!! NO. OF GENPART IS MORE THAN YOUR RESERVED NO.!!!!!!" << endl;}
      }
    }
  } // end of event is not data
   nGenPart = GenPart_pt.size();  // truncated to maximum size
    
//////////////////////////////////////////////////////////////////////////////
/////////////////////////////// Gen Particle End /////////////////////////////
//////////////////////////////////////////////////////////////////////////////

   //cout << "hello beam spot" << endl;

/////////////////////////////////////////////////////////////////////////////
///////////////////////////// Beam spot /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
  // to get BeamSpot information:
  // Handle<reco::BeamSpot> recoBeamSpotHandle;
  // iEvent.getByLabel("offlineBeamSpot", recoBeamSpotHandle);
  // reco::BeamSpot vertexBeamSpot= *recoBeamSpotHandle;  
  reco::BeamSpot vertexBeamSpot= *beamSpotHandle;  

  // build dummy beam spot copy for later use
  //reco::BeamSpot dummyBeamSpot = vertexBeamSpot;  

  // store beam spot info 
  Bsp_x = vertexBeamSpot.x0();
  Bsp_y = vertexBeamSpot.y0();
  Bsp_z = vertexBeamSpot.z0();
  Bsp_sigmaz = vertexBeamSpot.sigmaZ();
  Bsp_dxdz = vertexBeamSpot.dxdz();
  Bsp_dydz = vertexBeamSpot.dydz();
  Bsp_widthx = vertexBeamSpot.BeamWidthX();
  Bsp_widthy = vertexBeamSpot.BeamWidthY();

          
  //cout << "hello tracks and vertices" << endl;

//////////////////////////////////////////////////////////////////////////////
///////////////////////////// Tracks and vertices Start //////////////////////
//////////////////////////////////////////////////////////////////////////////


  // reimplement relevant track and vertex reference histograms
  // use includes (can be commented in and out by hand or by IFDEF)

  // nanoAOD

  // fill basic track information (all tracks in event)
  nTrk = tracks->size();
  PV_npvs = Primvertex->size();

  // nanoAOD
  nOtherPV = 0;

  // provide list of isotracks
  nIsoTrack =0;

  // nanoAOD-like extension
  nPVtx=0;

  //  cout << "hello vertices" << endl;

//////////////////////////////////////////////////////////////////////////////
//////////////////////////// Tools for Vertices //////////////////////////////
//////////////////////////////////////////////////////////////////////////////
    
  // load the tools to work with vertices: 
  // declare new track builder for new Transient track collection
  ESHandle<TransientTrackBuilder> theB;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

  // clear the storage containers for the objects in this event
  OtherPV_z.clear();

  PVtx_Id.clear();
  PVtx_isMain.clear();
  PVtx_isMainSim.clear();
  PVtx_isGood.clear();
  PVtx_isValid.clear();
  PVtx_isFake.clear();
  PVtx_isTrigUnique.clear();
  PVtx_isUnbiased.clear();
  PVtx_ntrk.clear();
  PVtx_ntrkfit.clear();
  PVtx_chi2.clear();
  PVtx_ndof.clear();
  PVtx_score.clear();
  PVtx_sumPt.clear();
  PVtx_Rho.clear();
  PVtx_x.clear();
  PVtx_y.clear();
  PVtx_z.clear();
  PVtx_Covxx.clear();
  PVtx_Covyx.clear();
  PVtx_Covzx.clear();
  PVtx_Covyy.clear();
  PVtx_Covzy.clear();
  PVtx_Covzz.clear();

  int nmuonlist = 0;
  int muid[maxnmuonlist];
  int muvtx[maxnmuonlist];
  double mupvx[maxnmuonlist];
  double mupvy[maxnmuonlist];
  double mupvz[maxnmuonlist];

  PV_npvsGood = 0;
  nOtherPV = 0;
  float PVSimMindist = 999.;
  int PVSimMinId = -1;
  //  for (unsigned int t = 0; t<hVtx->size(); ++t) {
  int vtxid =-1;
  int null=0; 

  // set primary vertex reference point pv
  math::XYZPoint pv(Primvertex->begin()->position());
  // to be filled later
  math::XYZPoint bestpv;

  reco::VertexCollection::const_iterator PV_ite;

  for (reco::VertexCollection::const_iterator vite = Primvertex->begin(); 
       vite != Primvertex->end(); ++vite) {
    if (PVtx_z.size() < nReserve_PVtx) {
      // primVertex_tmp = hVtx->at(t);
      // PVtx_Id.push_back(t);
      ++vtxid;
      // fill nanaoAOD vertex structure
      // by convention, first vertex is main primary *** check ***
      if (vtxid == 0) {
        PV_ite = vite;
        PV_chi2 = vite->chi2()/vite->ndof();
        PV_ndof = vite->ndof();
        PV_x = vite->x();
        PV_y = vite->y();
        PV_z = vite->z();
        PVtx_isMain.push_back(1);
      }  
      else if (OtherPV_z.size() < nReserve_OtherPV) {
	PVtx_isMain.push_back(null);
        // should this be sorted according to score first?
        ++nOtherPV;
        OtherPV_z.push_back(vite->z());
      }
      else {        
	PVtx_isMain.push_back(null);
      }
      // fill extended vertex structure
      // just in case not all vertices are stored (otherwise redundant)
      PVtx_Id.push_back(vtxid);
      // PVtx_chi2.push_back(primVertex_temp.chi2());
      // PVtx_ndof.push_back(primVertex_tmp.ndof());
      // PVtx_score.push_back(primVertex_tmp.score());
      // PVtx_x.push_back(primVertex_tmp.x());
      // PVtx_y.push_back(primVertex_tmp.y());
      // PVtx_z.push_back(primVertex_tmp.z());
      // *** find way to fill the following with sumntrk variable later ***
      PVtx_ntrk.push_back(999);
      PVtx_ntrkfit.push_back(vite->nTracks());
      PVtx_chi2.push_back(vite->chi2());
      PVtx_ndof.push_back(vite->ndof());
      float score=0., sumpt=0.;
      // flag to indicate whether vertex has muon
      //   *** should this be changed to dz/dxy requirement? ***
      //       (will also work for miniAOD) 
      bool hasmuon = false;
      float hasmuonpt = 0.;
      // loop over all tracks from this vertex
      // (In miniAOD the vertices have "lost" their tracks -> loop is dummy:
      //  *** should loop over all PackedCandidates and check vertexRef and
      //      muon property *** )
      for (reco::Vertex::trackRef_iterator iTrack = vite->tracks_begin(); iTrack != vite->tracks_end(); ++iTrack) {
	// get track reference to full track structure
        const reco::TrackRef trackRef = iTrack->castTo<reco::TrackRef>();
        // somehow store track-vertex relation here?
        // ...
        // get and sum track pt
        float trackpt = trackRef->pt();
        score += trackpt*trackpt;
        // check whether track is muon candidate *** apply cuts??? ***
        bool ismuon=false;
        int muonid=-1;
#ifndef miniAOD
        for (reco::MuonCollection::const_iterator itMuon = muons->begin(); itMuon != muons->end(); ++itMuon) {
#endif
#ifdef miniAOD
	for (pat::MuonCollection::const_iterator itMuon = muons->begin(); itMuon != muons->end(); ++itMuon) {
#endif
          // indicates position in reco::muon structure 
          ++muonid;
          if((itMuon->track()).isNonnull()){
            if (trackRef == itMuon->track()) {
              ismuon=true;
              hasmuon=true;
              // *** is the muon-vertex association loose enough? ***
              if ((itMuon->track())->pt() > hasmuonpt) hasmuonpt = (itMuon->track())->pt();
              // cout << "has muon! " << (itMuon->track())->pt() << " track " << trackRef->pt() << endl;
              // store muon-vertex association for later use
              if (nmuonlist < maxnmuonlist) {
                muid[nmuonlist]=muonid;
                muvtx[nmuonlist]=vtxid;
                mupvx[nmuonlist]=vite->x();
                mupvy[nmuonlist]=vite->y();
                mupvz[nmuonlist]=vite->z();
                ++nmuonlist;
              }
              else { cout << "Warning !!!! maxnmuonlist too small" << endl;}
              break;
            } // if trackref
          } // if itMuon
        } // for MuonCollection
        // sum track pt excluding muons
        if (!ismuon) sumpt += trackpt;         
	} // for VertexTracks
      // vite->score() does not seem to exist
      if (vtxid==0) PV_score = score;
      PVtx_score.push_back(score);
      PVtx_sumPt.push_back(sumpt);
      PVtx_Rho.push_back(vite->position().Rho());
      PVtx_x.push_back(vite->position().x());
      PVtx_y.push_back(vite->position().y());
      PVtx_z.push_back(vite->position().z());
      // is it possible that vite->z and vite->position.z differ?
      if (vite->z() != vite->position().z())
         cout << "*** vertex Alarm *** " << vite->z() << " " 
              << vite->position().z() << endl;
      PVtx_Covxx.push_back(vite->covariance(0,0));
      PVtx_Covyx.push_back(vite->covariance(1,0));
      PVtx_Covzx.push_back(vite->covariance(2,0));
      PVtx_Covyy.push_back(vite->covariance(1,1));
      PVtx_Covzy.push_back(vite->covariance(2,1));
      PVtx_Covzz.push_back(vite->covariance(2,2));

      // find good vertices
      // for QCD0-5 MC, about 10% of main simulated vertices are not 
      // reconstructed at all (no tracks in acceptance?) and 
      // a further 10% is not GOOD: most with |z|>10, a few with ndof<=4, 
      // all seem to be Valid, so Valid flag not needed?
      // increase "good" criterion from 10->24 cm in z 
      // and add "Rho" cut (according to doc) (radius?)
      // *** should position be absolute or relative to beam spot? ***
      int vtxGood=0;
      if (!vite->isFake() && vite->isValid() &&
          vite->ndof()>4 && fabs(vite->z()-Bsp_z)<24. &&
          vite->position().Rho() < 2.){ 
        vtxGood=1;
        ++PV_npvsGood;
      }
      float rhotest = sqrt(vite->position().x()*vite->position().x() + 
                           vite->position().y()*vite->position().y());
      // to check that the two are the same up to float precision 
      if (rhotest != float(vite->position().Rho()))
        cout << "*** Alarm Rho *** " << rhotest << " " 
             << vite->position().Rho() << endl;
      PVtx_isGood.push_back(vtxGood);
      PVtx_isValid.push_back(vite->isValid());
      PVtx_isFake.push_back(vite->isFake());

      // check vertices against simulated vertex
      //   ... should treat the case that there is more than one ...
      // The following sets the flag to true for the first vertex with
      // distance < 0.1 cm, and then for all subsequent vertices with 
      // consecutively smaller distance
      // GenPV_recIdx will either contain the first or the closest good vertex
      if (fabs(vite->z()-GenPV_z)<0.2 && fabs(vite->z()-GenPV_z)<PVSimMindist) {
        PVSimMindist = abs(vite->z()-GenPV_z);
        // always take the first, supersede if closer good vertex
        if (PVSimMinId==-1 || vtxGood==1) PVSimMinId = vtxid;
        PVtx_isMainSim.push_back(1);
        GenPV_recIdx = PVSimMinId;
        // to get unique simulated vertex, can check offline: 
        // PVtx_isMainSim && PVtx_Id == GenPV_recIdx
      }
      else {
        PVtx_isMainSim.push_back(0);
        // GenPV_recIdx = -1 has been preset elsewhere, don't overwrite here!;
      }

      // check whether this vertex has uniquely triggered a non-minimum-bias 
      // trigger (not yet fully working)
      PVtx_isTrigUnique.push_back(null);
      // For data: 
      // (technically not disabled for MC, which always has a ZeroBias Trigger)
      // Declare a vertex `unbiased' if a genuine MinimumBias trigger has fired
      if (GoodMinimumBiasTrig) PVtx_isUnbiased.push_back(1);
      // or a pure muon trigger has fired, good muons have been found, 
      // and the vertex does not have any muon candidate 
      // (*** should be changed to `does not have a triggered muon', and be 
      //  expanded to Electron datasets ***)
      else if (GoodMuTrig && (((MuThresh>-1 || IsoMuThresh>-1) && nMuon>0) || (DoubleMuThresh >-1 && nMuon>1)) && !hasmuon) PVtx_isUnbiased.push_back(1);
      // or a single muon trigger has fired, there are at least two muons, 
      // and the highest pt vertex muon has less than half the pt of the single 
      // muon trigger threshold
      else if (GoodMuTrig && (MuThresh>2.*hasmuonpt || IsoMuThresh>2.*hasmuonpt) && muons->size()>1) PVtx_isUnbiased.push_back(1);
      // note that MuHad and MuOnia triggers do not yield unbiased vertices 
      else PVtx_isUnbiased.push_back(0);
      // net effect: MB events will always be accepted; 
      // NMB events will be accepted if they do not have muon
    }
    else {cout << "WARNING!!!!! NO. OF vertices IS MORE THAN YOUR RESERVED NO.!!!!!!" << endl;}
  }
  nPVtx = PVtx_z.size();

//////////////////////////////////////////////////////////////////////////////
///////////////////////////// Tracks and vertices End ////////////////////////
//////////////////////////////////////////////////////////////////////////////


  //cout << "hello muon" << endl;

    
//////////////////////////////////////////////////////////////////////////////
///////////////////////////// Muon Collection Start //////////////////////////
//////////////////////////////////////////////////////////////////////////////

  // reintroduce muon reference histograms (pt, eta, quality, mass plots)
  // (reimport 'old' code via include?)

  // define and initialize all relevant temporary muon variables
  Float_t mu_relPFIsoR03 = -9999.;
  Float_t mu_relPFIsoR04 = -9999.;
  // add isolation variables used in CMSSW 4-2-8
  Float_t mu_relIsoR03 = -9999.;
  Float_t mu_relIsoR05 = -9999.;

  Float_t mu_ip3d = -9999.;
  Float_t mu_ip3dBest = -9999.;
  Float_t mu_Errip3d = -9999.;
  Float_t mu_sip3d = -9999.;
  Float_t mu_sip3dBest = -9999.;

  // tlorentzvector
  TLorentzVector p4Mu;  
  p4Mu.SetPtEtaPhiE(0., 0., 0., 0.);

  // clear the storage containers for the objects in this event
  Muon_charge.clear();
  Muon_tightCharge.clear(); 
  Muon_pt.clear();
  Muon_ptErr.clear();
  Muon_eta.clear();
  Muon_phi.clear();  
  Muon_mass.clear();
  Muon_dxy.clear();
  Muon_dxyErr.clear();
  Muon_dz.clear();
  Muon_dzErr.clear();
  Muon_ip3d.clear();
  Muon_sip3d.clear();
  Muon_pfRelIso03_all.clear();
  Muon_pfRelIso03_chg.clear();
  Muon_pfRelIso04_all.clear();
  Muon_miniPFRelIso_all.clear(); 
  Muon_miniPFRelIso_chg.clear(); 
  Muon_jetIdx.clear(); 
  Muon_isGlobal.clear();
  Muon_isTracker.clear();
  Muon_isPFcand.clear();
  Muon_softId.clear();
  Muon_mediumId.clear(); 
  Muon_tightId.clear(); 
  Muon_highPtId.clear(); 
  Muon_nStations.clear();
  Muon_nTrackerLayers.clear(); 
  Muon_segmentComp.clear(); 
  Muon_cleanmask.clear(); 
  Muon_mvaTTH.clear(); 
  Muon_pdgId.clear(); 
  Muon_genPartFlav.clear(); 
  Muon_genPartIdx.clear(); 

  Muon_Id.clear(); 
  Muon_x.clear(); 
  Muon_y.clear(); 
  Muon_z.clear(); 
  Muon_dxyBest.clear();
  Muon_dzBest.clear();
  Muon_ip3dBest.clear();
  Muon_sip3dBest.clear();
  Muon_gpt.clear();
  Muon_geta.clear();
  Muon_gphi.clear();  
  Muon_looseId.clear();
  Muon_softId4.clear();
  Muon_softIdBest.clear();
  Muon_isNano.clear();
  Muon_isMini.clear();
  Muon_isGood.clear();
  Muon_isGoodLast.clear();
  Muon_isGoodAng.clear();
  Muon_isArbitrated.clear();
  Muon_isStandAlone.clear();
  Muon_isRPCcand.clear();
  Muon_nValid.clear(); 
  Muon_nPix.clear(); 
  Muon_Chi2.clear(); 
  Muon_gnValid.clear(); 
  Muon_gnPix.clear(); 
  Muon_gChi2.clear(); 
  Muon_gnValidMu.clear(); 
  Muon_vtxIdx.clear(); 
  Muon_vtxFlag.clear(); 
  Muon_trkIdx.clear();
  Muon_simIdx.clear(); 
  
  b4_nMuon = muons->size();
  Muon_nNano = 0;
  
  int muonid = -1;

  // set cuts *** to be tuned ***
  float mindzvtxcut = 1.;

#ifndef miniAOD
  for (reco::MuonCollection::const_iterator recoMuon = muons->begin(); recoMuon != muons->end(); ++recoMuon) {
#endif
#ifdef miniAOD
  for (pat::MuonCollection::const_iterator recoMuon = muons->begin(); recoMuon != muons->end(); ++recoMuon) {
#endif
    ++muonid;
    
    // preselection: make "or" of 
    // - Run 2 miniAOD preselection (CMSSW 7X and above, works for 53X):
    //   pT>5 or (pt>3 and (PF or global or tracker or standal. or RPC)) or PF
    //   -> isMini
    // - Run 1 Loose selection (CMSSW 52X and above):
    //      isPFMuon and (isGlobalMuon or is TrackerMuon)
    //   -> isLoose 
    // - Run 1 Soft ID (CMSSW 44X or below):
    //      OneStation Tight, >5 hits, >0 pixel hits, high purity,
    //      dxy <0.3, dz <20.
    //      also works for CMSSW 53X and Run 2
    //   -> isSoft
    // - 
    // 

    // prepare variables:
    bool mu_isGlobal = recoMuon->isGlobalMuon();
    bool mu_isTracker = recoMuon->isTrackerMuon();
    bool mu_isStandAlone = recoMuon->isStandAloneMuon();
    // n/a for 4_2_8:
#ifndef CMSSW42X
    bool mu_isPFMuon = recoMuon->isPFMuon();
    bool mu_isRPCMuon = recoMuon->isRPCMuon();
#endif
#ifdef CMSSW42X
    // set true or false?
    bool mu_isPFMuon = true;
    bool mu_isRPCMuon = false;
#endif

    // Loose ID, see
    // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonID
    // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonIDRun2
    bool mu_looseId = (mu_isGlobal || mu_isTracker) && mu_isPFMuon;

    // prepare variables:
    bool mu_isGood = muon::isGoodMuon(*recoMuon, TMOneStationTight);
    bool mu_isGoodLast = muon::isGoodMuon(*recoMuon, TMLastStationTight);
    bool mu_isGoodAng = muon::isGoodMuon(*recoMuon, TMLastStationAngTight);
    bool mu_isArbitrated = muon::isGoodMuon(*recoMuon, TrackerMuonArbitrated);
    int mu_nStations = recoMuon->numberOfMatchedStations(); 
    // will the following work for Standalone muons?
    int mu_nTrackerLayers = 0;
    int mu_nValidHits = 0;
    int mu_gnValidHits = 0;
    int mu_nPixelLayers = 0; 
    bool mu_highPurity = false;
    float mu_dxy = 0.;
    float mu_dxyBest = 0.;
    float mu_dxyError = -1.;
    float mu_dz = 0.;
    float mu_dzBest = 0.;
    float mu_dzError = -1.;
    int mu_vtxId = -1.;

    // set quality variables
    int  mu_bestflag = -1;
    if (mu_isGlobal || mu_isTracker) {
      mu_nTrackerLayers = recoMuon->innerTrack()->hitPattern().trackerLayersWithMeasurement();
      mu_nValidHits = recoMuon->innerTrack()->hitPattern().numberOfValidTrackerHits();
      mu_nPixelLayers = recoMuon->innerTrack()->hitPattern().pixelLayersWithMeasurement();
      mu_highPurity = recoMuon->innerTrack()->quality(reco::TrackBase::highPurity);
      if (recoMuon->globalTrack().isNonnull()) {
        mu_gnValidHits = recoMuon->globalTrack()->hitPattern().numberOfValidTrackerHits();
      }
      
      // w.r.t. main primary vertex
      // somehow often gives wrong sign 

      // global does not seem to give best result -> commented
      // if (mu_isGlobal) {
      //  mu_dxy = recoMuon->globalTrack()->dxy(pv);
      //  // take this from Stefan Wunsch setup 
      //  mu_dxyError = recoMuon->globalTrack()->d0Error();
      //  mu_dz = recoMuon->globalTrack()->dz(pv);
      //  mu_dzError = recoMuon->globalTrack()->dzError();
      // }
      //else {

#ifdef CMSSW42X
      // innertrack seems to work better than globaltrack
      // does it always exist??
      mu_dxy = recoMuon->innerTrack()->dxy(pv);
      // the following does not give better results
      // mu_dxy = recoMuon->innerTrack()->dxy(pv) * recoMuon->charge();
      // if (recoMuon->phi()<0) mu_dxy = -mu_dxy; 
      // dxy or d0??  *** investigate!
      //mu_dxyError = recoMuon->innerTrack()->dxyError();
      //   d0Error seems to work better than dxyError, 
      //   neither innertrack nor globaltrack work well,
      //   slightly better for innertrack 
      mu_dxyError = recoMuon->innerTrack()->d0Error();
      //   both globaltrack and inner track are reasonable, 
      //   slightly better for innertrack? 
      mu_dz = recoMuon->innerTrack()->dz(pv);
      // neither is good, but global track is somewhat better??
      mu_dzError = recoMuon->innerTrack()->dzError();
#else
      // use bestTrack when available in CMSSW
      mu_dxy = recoMuon->bestTrack()->dxy(pv);
      mu_dxyError = recoMuon->bestTrack()->d0Error();
      mu_dz = recoMuon->bestTrack()->dz(pv);
      mu_dzError = recoMuon->bestTrack()->dzError();
#endif
      // } // else

      // w.r.t. associated primary vertex, preset to beamspot
      //   beware, this will not work for muons close to, 
      //   but not associated to, vertex!
      // preset w.r.t. beamspot
      mu_dxyBest = recoMuon->innerTrack()->dxy(vertexBeamSpot);
      // does not work for z -> keep 0.
      // mu_dzBest = recoMuon->innerTrack()->dz(vertexBeamSpot);
      for (int mm = 0; mm<nmuonlist; ++mm) {
        if (muonid == muid[mm]) {
          // set primary vertex reference point bestpv
          bestpv.SetXYZ(mupvx[mm],mupvy[mm],mupvz[mm]);
          // fill distance
          mu_dxyBest = recoMuon->innerTrack()->dxy(bestpv);
          mu_dzBest = recoMuon->innerTrack()->dz(bestpv);
          mu_vtxId = muvtx[mm];
          mu_bestflag = 0;
          break;    // candidate found -> stop loop
        } //if  muonid
      } // for mm
    } // if global or tracker

    else if (mu_isStandAlone) {
      // standalone only muons
      mu_dxy = recoMuon->outerTrack()->dxy(pv);
      //mu_dxy = recoMuon->outerTrack()->dxy(pv) * recoMuon->charge();
      //if (recoMuon->phi()<0) mu_dxy = -mu_dxy; 
      mu_dxyBest = recoMuon->outerTrack()->dxy(vertexBeamSpot);
      mu_dxyError = recoMuon->outerTrack()->dxyError();
      mu_dz = recoMuon->outerTrack()->dz(pv);
      // somehow beamspot does not work for z (and does not make much sense)
      mu_dzBest = recoMuon->outerTrack()->dz();
      mu_dzError = recoMuon->outerTrack()->dzError();
      mu_bestflag = -2;
    } // else if standalone    

    // treat nonvertex-fitted muons 
    // if not directly associated to any vertex, find closest PV
    if (mu_bestflag<0) {
      float mindzvtx = 99.;
      int idmindzvtx = -1;
      // loop over PVs
      for (uint pvtx=0; pvtx<nPVtx; ++pvtx) {
        float dzvtx = recoMuon->vz()-PVtx_z[pvtx];
        if (abs(dzvtx)<abs(mindzvtx)) {
          mindzvtx = dzvtx;
          idmindzvtx = pvtx;
        } // if
      } // for
      // and store it if closer than some cut value
      if (abs(mindzvtx) < mindzvtxcut) {
        mu_bestflag = 1;
        mu_dzBest = mindzvtx;
        if (abs(mu_dzBest)>0.1) mu_bestflag=2;
        mu_vtxId = idmindzvtx;
      } // if
    } // if bestflag


    // Soft ID, see
    // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonID
    // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonIDRun2
    // for 53X and Run 2:
    bool mu_softId = mu_isGood && mu_nTrackerLayers > 5 && mu_nPixelLayers >0 && mu_highPurity && fabs(mu_dxy) < 0.3 && fabs(mu_dz)<20.;
    // for 42X:  (both defined for both)
    bool mu_softId4 = mu_isGood && mu_nValidHits > 10 && mu_nPixelLayers >1 && fabs(mu_dxy) < 3. && fabs(mu_dz)<30.;
    // relative to best (not main) pv; probably makes little difference
    bool mu_softIdBest = mu_isGood && mu_nValidHits > 10 && mu_nPixelLayers >1 && fabs(mu_dxyBest) < 3. && fabs(mu_dzBest)<30.;

    // MiniAOD preselection
    // to recover MiniAOD content, where available, use isMini
    bool mu_isMini = false;
#ifndef CMSSW42X
    // for 53X and higher
    if (recoMuon->pt() > 5. or (recoMuon->pt() > 3. and (recoMuon->isPFMuon() or recoMuon->isGlobalMuon() or recoMuon->isTrackerMuon() or recoMuon->isStandAloneMuon() or recoMuon->isRPCMuon())) or recoMuon->isPFMuon()) {
      mu_isMini = true;
    } // Mini
#endif
#ifdef CMSSW42X
    // for 42X: isPFMuon and isRPCMuon do not exist yet
    // *** need to add something to accept low pt muons? ***
    if (recoMuon->pt() > 5. or (recoMuon->pt() > 3. and (recoMuon->isGlobalMuon() or recoMuon->isTrackerMuon() or recoMuon->isStandAloneMuon())) ) {
      mu_isMini = true;
    } // Mini
#endif

    // ************************************************************
    // make preselection for ntuple storage: Mini or Loose or Soft:
    // ************************************************************
    // actually, removes almost nothing -> do not make cut, i.e. keep all AOD muons  
    //if (!mu_isMini && !mu_looseId && !mu_softId && !mu_softId4 && !mu_softIdBest) continue;

    // to recover nanoAOD content, where available, use isNano
    bool mu_isNano = false;
    if (mu_isMini && mu_looseId && recoMuon->pt() > 3.) {
      mu_isNano = true;
      ++Muon_nNano;
    } // Nano

    // medium ID
    // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonIDRun2
    float mu_validFraction = 0;
    if (recoMuon->innerTrack().isNonnull()) {
      mu_validFraction = recoMuon->innerTrack()->validFraction();
    }
    float mu_Chi2 = 0;
    if (recoMuon->innerTrack().isNonnull()) {
      mu_Chi2 = recoMuon->innerTrack()->normalizedChi2();
    }
    float mu_globalChi2 = 0;
    if (recoMuon->globalTrack().isNonnull()) {
      mu_globalChi2 = recoMuon->globalTrack()->normalizedChi2();
    }
    float mu_standPosMatch = recoMuon->combinedQuality().chi2LocalPosition;
    float mu_kink = recoMuon->combinedQuality().trkKink;  
    float mu_segmentComp = muon::segmentCompatibility(*recoMuon);
    bool mu_mediumId = false;
    if (mu_looseId && mu_validFraction>0.8 && ((mu_isGlobal && mu_globalChi2 < 3 && mu_standPosMatch<12 && mu_kink<20 && mu_segmentComp > 0.303) || mu_segmentComp >0.451)) {
      mu_mediumId = true;
    } // medium id

    // tight ID (with respect to main primary?), see
    // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonID
    // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonIDRun2
    int mu_gnValidMuHits = 0;
    if (recoMuon->globalTrack().isNonnull()) {
      mu_gnValidMuHits = recoMuon->globalTrack()->hitPattern().numberOfValidMuonHits();
    }

    // probably all of the following is obsolete
#ifndef CMSSW42X
    // ?? why not for 42X ??
    // tight id always uses main primary vertex //
    //float mu_dxyprim = recoMuon->muonBestTrack()->dxy(PV_ite->position());
    //float mu_dB = ...; // dB on PAT::Muon -> check example
    //float mu_dzprim = recoMuon->muonBestTrack()->dz(PV_ite->position());
#endif
    // actually not needed
    //#ifdef CMSSW42X
    // // tight id always uses main primary vertex and global muon //
    // float mu_dxyprim = 999.;
    // float mu_dzprim = 999.;
    // if (recoMuon->globalTrack().isNonnull()) {
    //   mu_dxyprim = recoMuon->globalTrack()->dxy(PV_ite->position());
    //  //float mu_dB = ...; // dB on PAT::Muon -> check example
    //  mu_dzprim = recoMuon->globalTrack()->dz(PV_ite->position());
    //}
    //#endif
    //float mu_dxyprim = mu_dxy;
    //float mu_dzprim = mu_dz;

    // *** fill this properly ***
    float mu_dB = mu_dxy;

    int mu_nValidPixel=0;
    if (recoMuon->innerTrack().isNonnull()) {
      mu_nValidPixel = recoMuon->innerTrack()->hitPattern().numberOfValidPixelHits();
    } 
    int mu_gnValidPixel=0;
    if (recoMuon->globalTrack().isNonnull()) {
      mu_gnValidPixel = recoMuon->globalTrack()->hitPattern().numberOfValidPixelHits();
    } 
#ifndef CMSSW42X
    // 53X or higher:
    bool mu_tightId = mu_isGlobal && mu_isPFMuon && mu_globalChi2<10 && mu_gnValidMuHits>0 && mu_nStations>1 && (fabs(mu_dxy)<0.2 || fabs(mu_dB)<0.2)
      && fabs(mu_dz)<0.5 && mu_nValidPixel>0 && mu_nTrackerLayers>5;     
#endif
#ifdef CMSSW42X
    // 42X:
    bool mu_tightId = mu_isGlobal && mu_globalChi2<10 && mu_gnValidMuHits>0 && mu_nStations>1 && (fabs(mu_dxy)<0.2 || fabs(mu_dB)<0.2)
    && mu_nValidPixel > 0 && mu_nTrackerLayers>8;
    // the latter is roughly equivalent to mu_nValidHits>10     
#endif
    // compare to muon::isTightMuon(*recomuon,*vertices->begin());
    // compare to muon::isSoftMuon(*recomuon,*vertices->begin());

    // high pt ID
    // choose one of the following:
    // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonID
    // use main vertex vertex
    //bool mu_highPtId = muon::isHighPtMuon(*recoMuon,*PV_ite, reco::TunePType = muon::improvedTuneP);
#ifdef CMSSW7plus
    UChar_t mu_highPtId = muon::isHighPtMuon(*recoMuon,*PV_ite);
#endif
#ifdef CMSSW53X
    UChar_t mu_highPtId = muon::isHighPtMuon(*recoMuon,*PV_ite, improvedTuneP);
    //bool mu_highPtId = muon::isHighPtMuon(*recoMuon,*PV_ite, muon::improvedTuneP);
#endif
#ifdef CMSSW42X
    // possibly for 42X only:
    //bool mu_highPtId = muon::isHighPtMuon(*recoMuon,*PV_ite);
    UChar_t mu_highPtId = false;
    //bool mu_highPtId = muon::isHighPtMuon(*recoMuon,*PV_ite, muon::TunePType= muon::defaultTuneP);
#endif

    // calculate relative isolation from PF
    // n/a in 4-2-8, *** need to IFDEF! ***  see commented include!
    // (search for CHANGE)
    // need to add determination from PAT::
#ifndef CMSSW42X
    // should we check recoMuon->isPFisolationValid?
    mu_relPFIsoR03 = ((recoMuon->pfIsolationR03()).sumChargedHadronPt +
    		      (recoMuon->pfIsolationR03()).sumNeutralHadronEt + 
    		      (recoMuon->pfIsolationR03()).sumPhotonEt) / recoMuon->pt();	
    mu_relPFIsoR04 = ((recoMuon->pfIsolationR04()).sumChargedHadronPt +
    		      (recoMuon->pfIsolationR04()).sumNeutralHadronEt + 
    		      (recoMuon->pfIsolationR04()).sumPhotonEt) / recoMuon->pt();
#endif

    // add isolation without PF instead 
    mu_relIsoR03 = ((recoMuon->isolationR03()).sumPt +
    		    (recoMuon->isolationR03()).hadEt + 
    		    (recoMuon->isolationR03()).emEt) / recoMuon->pt();	

    mu_relIsoR05 = ((recoMuon->isolationR05()).sumPt +
    		    (recoMuon->isolationR05()).hadEt + 
    		    (recoMuon->isolationR05()).emEt) / recoMuon->pt();
    
    //cout << "Muon pt " << recoMuon->pt() << " best " << recoMuon->bestTrack()->pt() << endl;
    //if (mu_isGlobal) { cout << "global " << recoMuon->globalTrack()->pt() << endl;} 
    //if (mu_isGlobal || mu_isTracker) { cout << "inner " << recoMuon->innerTrack()->pt() << endl;} 

    // calculate impact parameter and its significance
    mu_sip3d = 0;
    if (recoMuon->isGlobalMuon()) {  
      mu_ip3d = sqrt((recoMuon->globalTrack()->dxy(pv) * recoMuon->globalTrack()->dxy(pv)) + (recoMuon->globalTrack()->dz(pv) * recoMuon->globalTrack()->dz(pv)));
      if (mu_bestflag>-1) mu_ip3dBest = sqrt((recoMuon->globalTrack()->dxy(bestpv) * recoMuon->globalTrack()->dxy(bestpv)) + (recoMuon->globalTrack()->dz(bestpv) * recoMuon->globalTrack()->dz(bestpv))); 
      mu_Errip3d = sqrt((recoMuon->globalTrack()->d0Error() * recoMuon->globalTrack()->d0Error()) + (recoMuon->globalTrack()->dzError() * recoMuon->globalTrack()->dzError()));
      mu_sip3d = mu_ip3d / mu_Errip3d;
      if (mu_bestflag>-1) mu_sip3dBest = mu_ip3dBest / mu_Errip3d;
    } // global muon
    else if (recoMuon->isTrackerMuon()) {
      mu_ip3d = sqrt((recoMuon->innerTrack()->dxy(pv) * recoMuon->innerTrack()->dxy(pv)) + (recoMuon->innerTrack()->dz(pv) * recoMuon->innerTrack()->dz(pv)));
      if (mu_bestflag>-1) mu_ip3dBest = sqrt((recoMuon->innerTrack()->dxy(bestpv) * recoMuon->innerTrack()->dxy(bestpv)) + (recoMuon->innerTrack()->dz(bestpv) * recoMuon->innerTrack()->dz(bestpv)));
      mu_Errip3d = sqrt((recoMuon->innerTrack()->d0Error() * recoMuon->innerTrack()->d0Error()) + (recoMuon->innerTrack()->dzError() * recoMuon->innerTrack()->dzError()));
      mu_sip3d = mu_ip3d / mu_Errip3d;
      if (mu_bestflag>-1) mu_sip3dBest = mu_ip3dBest / mu_Errip3d;
    } // tracker muon
 
      // now store the muon candidate 
      if (Muon_pt.size() < nReserve_Muon) {

       // store nanoAOD extensions only if flag set
       if (nanoext || mu_isNano) {
         
        // official nanoAOD variables

        // kinematic quantities, default is tracker muon quantities, 
        // when available, best other otherwise 
	Muon_charge.push_back(recoMuon->charge());
	Muon_pt.push_back(recoMuon->pt());
        float mu_ptErr = -1.;
#ifndef CMSSW42X
	mu_ptErr = recoMuon->muonBestTrack()->ptError();
#endif
#ifdef CMSSW42X	
        if (recoMuon->innerTrack().isNonnull()) {
	  mu_ptErr = recoMuon->innerTrack()->ptError();
        }
        else if (recoMuon->globalTrack().isNonnull()) {
	  mu_ptErr = recoMuon->globalTrack()->ptError();
        }
#endif	
	Muon_ptErr.push_back(mu_ptErr);
        // ?(muonBestTrack().ptError()/muonBestTrack().pt() < 0.2)?2:0
        if (mu_ptErr/recoMuon->pt()<0.2) Muon_tightCharge.push_back(2);
        else Muon_tightCharge.push_back(0);
	Muon_eta.push_back(recoMuon->eta());
	Muon_phi.push_back(recoMuon->phi());
	Muon_mass.push_back(mumass);
        //cout << "Muon pt, eta, phi " << recoMuon->charge()*recoMuon->pt() << " " << recoMuon->eta() << " " << recoMuon->phi() << endl; 

	// dxy, dz, ip
	Muon_dxy.push_back(mu_dxy);
	Muon_dxyBest.push_back(mu_dxyBest);
	Muon_dxyErr.push_back(mu_dxyError);
	Muon_dz.push_back(mu_dz);
	Muon_dzBest.push_back(mu_dzBest);
	Muon_dzErr.push_back(mu_dzError);	
	Muon_ip3d.push_back(mu_ip3d);
	Muon_sip3d.push_back(mu_sip3d);
	Muon_ip3dBest.push_back(mu_ip3dBest);
	Muon_sip3dBest.push_back(mu_sip3dBest);

        // isolation
        if (mu_relPFIsoR04 != -9999.) {
	  Muon_pfRelIso03_all.push_back(mu_relPFIsoR03);
	  Muon_pfRelIso04_all.push_back(mu_relPFIsoR04);
        } // if
	else { // for CMSSW 4-2-8
	  Muon_pfRelIso03_all.push_back(mu_relIsoR03);
	  Muon_pfRelIso04_all.push_back(mu_relIsoR05);
        } // else
#ifndef CMSSW42X     
	Muon_pfRelIso03_chg.push_back((recoMuon->pfIsolationR03()).sumChargedHadronPt / recoMuon->pt());
#endif
#ifdef CMSSW42X
        // not really appropriate (hadEt is presumably HCAL component) 
	// Muon_pfRelIso03_chg.push_back((recoMuon->isolationR03()).hadEt / recoMuon->pt());
	Muon_pfRelIso03_chg.push_back((recoMuon->isolationR03()).sumPt / recoMuon->pt());
#endif
	Muon_miniPFRelIso_all.push_back(999.); // *
	Muon_miniPFRelIso_chg.push_back(999.); // *

	Muon_jetIdx.push_back(999);

        Muon_isGlobal.push_back(mu_isGlobal); 
        Muon_isTracker.push_back(mu_isTracker); 
	Muon_isPFcand.push_back(mu_isPFMuon);

	// IDs (not directly available in CMSSW53X)
	// Muon_mediumId.push_back(recoMuon->passed(recoMuon->CutBasedIdMedium));
	// Muon_softId.push_back(recoMuon->passed(recoMuon->SoftCutBasedId));
	// Muon_tightId.push_back(recoMuon->passed(recoMuon->CutBasedIdTight));
        // calculated by hand above //
	Muon_softId.push_back(mu_softId);
	Muon_mediumId.push_back(mu_mediumId); 
	Muon_tightId.push_back(mu_tightId);
	Muon_highPtId.push_back(mu_highPtId);

	Muon_nStations.push_back(mu_nStations);

	// unsigned should be int & can't be (-)
	// Fill these variables with something nonsense atm
        // not yet being sorted!
	Muon_cleanmask.push_back('9'); // *
	Muon_genPartFlav.push_back('9'); // *
	Muon_genPartIdx.push_back(999); // *
	Muon_mvaTTH.push_back(999.); // *
	Muon_nTrackerLayers.push_back(mu_nTrackerLayers);
	// this needs to be changed to id after truth matching 
	Muon_pdgId.push_back(-13*recoMuon->charge());
	Muon_segmentComp.push_back(mu_segmentComp);

        ///////////////////////////////
        // nanoAOD-like extensions
        ///////////////////////////////

        // this is relevant when not all muon candidates are stored, 
        // redundant otherwise
        Muon_Id.push_back(muonid);

        Muon_x.push_back(recoMuon->vx());
        Muon_y.push_back(recoMuon->vy());
        Muon_z.push_back(recoMuon->vz());

        if (recoMuon->globalTrack().isNonnull()) {
	  Muon_gpt.push_back(recoMuon->globalTrack()->pt());
          Muon_geta.push_back(recoMuon->globalTrack()->eta());
          Muon_gphi.push_back(recoMuon->globalTrack()->phi());
        }

        Muon_looseId.push_back(mu_looseId);
	Muon_softId4.push_back(mu_softId4);
	Muon_softIdBest.push_back(mu_softIdBest);
        Muon_isMini.push_back(mu_isMini); 
        Muon_isNano.push_back(mu_isNano); 
        Muon_isGood.push_back(mu_isGood); 
        Muon_isGoodLast.push_back(mu_isGoodLast); 
        Muon_isGoodAng.push_back(mu_isGoodAng); 
        Muon_isArbitrated.push_back(mu_isArbitrated); 
        Muon_isStandAlone.push_back(mu_isStandAlone); 
	Muon_isRPCcand.push_back(mu_isRPCMuon);
	Muon_nValid.push_back(mu_nValidHits);
	Muon_nPix.push_back(mu_nValidPixel);
	Muon_Chi2.push_back(mu_Chi2);
	Muon_gnValid.push_back(mu_gnValidHits);
	Muon_gnPix.push_back(mu_gnValidPixel);
	Muon_gChi2.push_back(mu_globalChi2);
        Muon_gnValidMu.push_back(mu_gnValidMuHits);

        Muon_vtxIdx.push_back(mu_vtxId);
        Muon_vtxFlag.push_back(mu_bestflag);

#ifndef miniAOD
	// find muon in generaltracks list
        int itcount = -1;
        int mu_trkIdx = -1;
        if ((recoMuon->innerTrack()).isNonnull()) {
          for (reco::TrackCollection::const_iterator itrack = tracks->begin(); itrack != tracks->end(); ++itrack) {
            ++itcount;
            //if (*(itrack) == *(recoMuon->track())) {
	    //if (itrack->track() == recoMuon->track()) {
	    //if (itrack->get<TrackRef>() == recoMuon->track()) {
	    //if (*(itrack)->get<TrackRef>() == recoMuon->track()) {
	    //if (*(itrack).get<TrackRef>() == recoMuon->track()) {
	    //if (itrack == (recoMuon->track())) {
	    //if (itrack == (recoMuon->track())->get()) {
	    //if (itrack == (recoMuon->track()).get()) {
	    //if (itrack == recoMuon->track().get()) {
	    //if (itrack->ref() == (recoMuon->track())) {
	    //if (*(itrack) == (recoMuon->track())) {          
	    //if (itrack->pt() == (recoMuon->track())->pt()) {
            // capitulate: (is this safe?)
	    if (itrack->pt() == recoMuon->innerTrack()->pt()) {
              // muon found
              mu_trkIdx = itcount;
              break; 
	    } // if
          } // for
        } // nonnull
        Muon_trkIdx.push_back(mu_trkIdx);
#endif 
#ifdef miniAOD
        // not yet treated
        Muon_trkIdx.push_back(-1); // *
#endif

        // check for associated simulated muon from prestored list
        bool musimfound = false;
        for (int im=0; im<nmusim; ++im) {
          if (recoMuon->charge()==chgmusim[im] && abs(recoMuon->pt()-ptmusim[im])/ptmusim[im] < 0.1 && abs(recoMuon->eta()-etamusim[im])<0.1 && fmod(abs(recoMuon->phi()-phimusim[im]),2.*pi)<0.1) {
            musimfound = true;
            Muon_simIdx.push_back(idmusim[im]);
            break;
          } // if
	} // im
        if (!musimfound) {
          Muon_simIdx.push_back(-1);
        } // !musimfound    

       } // nanoext
      } // muon size
      else {
        cout << "WARNING!!!!! NO. OF Mu IS MORE THAN YOUR RESERVED NO.!!!!!!" << endl;
      } // else

  } // end of loop recoMuon

  // set actual size
  nMuon = Muon_pt.size();

  // Resort muons: isNano, ordered according to pt, first! 
  // then the others
  for (int i1=0; i1<int(nMuon)-1; ++i1) {
    for (int i2= i1+1; i2<int(nMuon); ++i2) {
      if (Muon_isNano[i2] && (!Muon_isNano[i1] || 
			      Muon_pt[i2] > Muon_pt[i1])) {

        Int_t tempi = Muon_charge[i1];
        Muon_charge[i1] = Muon_charge[i2];
        Muon_charge[i2] = tempi;    
        tempi = Muon_tightCharge[i1];
        Muon_tightCharge[i1] = Muon_tightCharge[i2];
        Muon_tightCharge[i2] = tempi;    
        Float_t temp = Muon_pt[i1];
        Muon_pt[i1] = Muon_pt[i2];
        Muon_pt[i2] = temp;    
        temp = Muon_ptErr[i1];
        Muon_ptErr[i1] = Muon_ptErr[i2];
        Muon_ptErr[i2] = temp;    
        temp = Muon_eta[i1];
        Muon_eta[i1] = Muon_eta[i2];
        Muon_eta[i2] = temp;    
        temp = Muon_phi[i1];
        Muon_phi[i1] = Muon_phi[i2];
        Muon_phi[i2] = temp;    
        temp = Muon_mass[i1];
        Muon_mass[i1] = Muon_mass[i2];
        Muon_mass[i2] = temp;    

        temp = Muon_dxy[i1];
        Muon_dxy[i1] = Muon_dxy[i2];
        Muon_dxy[i2] = temp;    
        temp = Muon_dxyBest[i1];
        Muon_dxyBest[i1] = Muon_dxyBest[i2];
        Muon_dxyBest[i2] = temp;    
        temp = Muon_dxyErr[i1];
        Muon_dxyErr[i1] = Muon_dxyErr[i2];
        Muon_dxyErr[i2] = temp;    
        temp = Muon_dz[i1];
        Muon_dz[i1] = Muon_dz[i2];
        Muon_dz[i2] = temp;    
        temp = Muon_dzBest[i1];
        Muon_dzBest[i1] = Muon_dzBest[i2];
        Muon_dzBest[i2] = temp;    
        temp = Muon_dzErr[i1];
        Muon_dzErr[i1] = Muon_dzErr[i2];
        Muon_dzErr[i2] = temp;    
        temp = Muon_ip3d[i1];
        Muon_ip3d[i1] = Muon_ip3d[i2];
        Muon_ip3d[i2] = temp;    
        temp = Muon_sip3d[i1];
        Muon_sip3d[i1] = Muon_sip3d[i2];
        Muon_sip3d[i2] = temp;    
        temp = Muon_ip3dBest[i1];
        Muon_ip3dBest[i1] = Muon_ip3dBest[i2];
        Muon_ip3dBest[i2] = temp;    
        temp = Muon_sip3dBest[i1];
        Muon_sip3dBest[i1] = Muon_sip3dBest[i2];
        Muon_sip3dBest[i2] = temp;    

        temp = Muon_pfRelIso03_all[i1];
        Muon_pfRelIso03_all[i1] = Muon_pfRelIso03_all[i2];
        Muon_pfRelIso03_all[i2] = temp;
        temp = Muon_pfRelIso03_chg[i1];
        Muon_pfRelIso03_chg[i1] = Muon_pfRelIso03_chg[i2];
        Muon_pfRelIso03_chg[i2] = temp;
        temp = Muon_pfRelIso04_all[i1];
        Muon_pfRelIso04_all[i1] = Muon_pfRelIso04_all[i2];
        Muon_pfRelIso04_all[i2] = temp;
        temp = Muon_miniPFRelIso_all[i1];
        Muon_miniPFRelIso_all[i1] = Muon_miniPFRelIso_all[i2];
        Muon_miniPFRelIso_all[i2] = temp;
        temp = Muon_miniPFRelIso_chg[i1];
        Muon_miniPFRelIso_chg[i1] = Muon_miniPFRelIso_chg[i2];
        Muon_miniPFRelIso_chg[i2] = temp;

        tempi = Muon_jetIdx[i1];
        Muon_jetIdx[i1] = Muon_jetIdx[i2];
        Muon_jetIdx[i2] = tempi;    

        uint8_t tempi8 = Muon_isGlobal[i1];
        Muon_isGlobal[i1] = Muon_isGlobal[i2];
        Muon_isGlobal[i2] = tempi8;    
        tempi8 = Muon_isTracker[i1];
        Muon_isTracker[i1] = Muon_isTracker[i2];
        Muon_isTracker[i2] = tempi8;    
        tempi8 = Muon_isPFcand[i1];
        Muon_isPFcand[i1] = Muon_isPFcand[i2];
        Muon_isPFcand[i2] = tempi8;    

        tempi8 = Muon_softId[i1];
        Muon_softId[i1] = Muon_softId[i2];
        Muon_softId[i2] = tempi8;    
        tempi8 = Muon_mediumId[i1];
        Muon_mediumId[i1] = Muon_mediumId[i2];
        Muon_mediumId[i2] = tempi8;    
        tempi8 = Muon_tightId[i1];
        Muon_tightId[i1] = Muon_tightId[i2];
        Muon_tightId[i2] = tempi8;    
        UChar_t tempU = Muon_highPtId[i1];
        Muon_highPtId[i1] = Muon_highPtId[i2];
        Muon_highPtId[i2] = tempU;    

        tempi = Muon_nStations[i1];
        Muon_nStations[i1] = Muon_nStations[i2];
        Muon_nStations[i2] = tempi;    
        tempi = Muon_nTrackerLayers[i1];
        Muon_nTrackerLayers[i1] = Muon_nTrackerLayers[i2];
        Muon_nTrackerLayers[i2] = tempi;    
        temp = Muon_segmentComp[i1];
        Muon_segmentComp[i1] = Muon_segmentComp[i2];
        Muon_segmentComp[i2] = temp;    

        tempi = Muon_pdgId[i1];
        Muon_pdgId[i1] = Muon_pdgId[i2];
        Muon_pdgId[i2] = tempi;    

        tempU = Muon_cleanmask[i1];
        Muon_cleanmask[i1] = Muon_cleanmask[i2];
        Muon_cleanmask[i2] = tempU;    
        temp = Muon_mvaTTH[i1];
        Muon_mvaTTH[i1] = Muon_mvaTTH[i2];
        Muon_mvaTTH[i2] = temp;    
        tempU = Muon_genPartFlav[i1];
        Muon_genPartFlav[i1] = Muon_genPartFlav[i2];
        Muon_genPartFlav[i2] = tempU;    
        tempi = Muon_genPartIdx[i1];
        Muon_genPartIdx[i1] = Muon_genPartIdx[i2];
        Muon_genPartIdx[i2] = tempi;    

        tempi = Muon_Id[i1];
        Muon_Id[i1] = Muon_Id[i2];
        Muon_Id[i2] = tempi;    

        temp = Muon_x[i1];
        Muon_x[i1] = Muon_x[i2];
        Muon_x[i2] = temp;    
        temp = Muon_y[i1];
        Muon_y[i1] = Muon_y[i2];
        Muon_y[i2] = temp;    
        temp = Muon_z[i1];
        Muon_z[i1] = Muon_z[i2];
        Muon_z[i2] = temp;    

        temp = Muon_gpt[i1];
        Muon_gpt[i1] = Muon_gpt[i2];
        Muon_gpt[i2] = temp;    
        temp = Muon_geta[i1];
        Muon_geta[i1] = Muon_geta[i2];
        Muon_geta[i2] = temp;    
        temp = Muon_gphi[i1];
        Muon_gphi[i1] = Muon_gphi[i2];
        Muon_gphi[i2] = temp;    

        tempi8 = Muon_looseId[i1];
        Muon_looseId[i1] = Muon_looseId[i2];
        Muon_looseId[i2] = tempi8;    
        tempi8 = Muon_softId4[i1];
        Muon_softId4[i1] = Muon_softId4[i2];
        Muon_softId4[i2] = tempi8;    
        tempi8 = Muon_softIdBest[i1];
        Muon_softIdBest[i1] = Muon_softIdBest[i2];
        Muon_softIdBest[i2] = tempi8;    
        tempi8 = Muon_isMini[i1];
        Muon_isMini[i1] = Muon_isMini[i2];
        Muon_isMini[i2] = tempi8;    
        tempi8 = Muon_isNano[i1];
        Muon_isNano[i1] = Muon_isNano[i2];
        Muon_isNano[i2] = tempi8;    
        tempi8 = Muon_isGood[i1];
        Muon_isGood[i1] = Muon_isGood[i2];
        Muon_isGood[i2] = tempi8;    
        tempi8 = Muon_isGoodLast[i1];
        Muon_isGoodLast[i1] = Muon_isGoodLast[i2];
        Muon_isGoodLast[i2] = tempi8;    
        tempi8 = Muon_isGoodAng[i1];
        Muon_isGoodAng[i1] = Muon_isGoodAng[i2];
        Muon_isGoodAng[i2] = tempi8;    
        tempi8 = Muon_isArbitrated[i1];
        Muon_isArbitrated[i1] = Muon_isArbitrated[i2];
        Muon_isArbitrated[i2] = tempi8;    
        tempi8 = Muon_isStandAlone[i1];
        Muon_isStandAlone[i1] = Muon_isStandAlone[i2];
        Muon_isStandAlone[i2] = tempi8;    
        tempi8 = Muon_isRPCcand[i1];
        Muon_isRPCcand[i1] = Muon_isRPCcand[i2];
        Muon_isRPCcand[i2] = tempi8;    
        tempi = Muon_nValid[i1];
        Muon_nValid[i1] = Muon_nValid[i2];
        Muon_nValid[i2] = tempi;    
        tempi = Muon_nPix[i1];
        Muon_nPix[i1] = Muon_nPix[i2];
        Muon_nPix[i2] = tempi;    
        temp = Muon_Chi2[i1];
        Muon_Chi2[i1] = Muon_Chi2[i2];
        Muon_Chi2[i2] = temp;    
        tempi = Muon_gnValid[i1];
        Muon_gnValid[i1] = Muon_gnValid[i2];
        Muon_gnValid[i2] = tempi;    
        tempi = Muon_gnPix[i1];
        Muon_gnPix[i1] = Muon_gnPix[i2];
        Muon_gnPix[i2] = tempi;    
        temp = Muon_gChi2[i1];
        Muon_gChi2[i1] = Muon_gChi2[i2];
        Muon_gChi2[i2] = temp;    
        tempi = Muon_gnValidMu[i1];
        Muon_gnValidMu[i1] = Muon_gnValidMu[i2];
        Muon_gnValidMu[i2] = tempi;    

        tempi = Muon_vtxIdx[i1];
        Muon_vtxIdx[i1] = Muon_vtxIdx[i2];
        Muon_vtxIdx[i2] = tempi;    
        tempi = Muon_vtxFlag[i1];
        Muon_vtxFlag[i1] = Muon_vtxFlag[i2];
        Muon_vtxFlag[i2] = tempi;    
        tempi = Muon_trkIdx[i1];
        Muon_trkIdx[i1] = Muon_trkIdx[i2];
        Muon_trkIdx[i2] = tempi;    
        tempi = Muon_simIdx[i1];
        Muon_simIdx[i1] = Muon_simIdx[i2];
        Muon_simIdx[i2] = tempi;    
      } // isNano
    } // i2
  } // i1

  // cout << "IS MUON LOOP OK? " << nMuon << endl;
  
  //cout << "hello dimuon" << endl;

  //------------------------------- For Dimu branches -----------------------//
  nDimu = 0;

  // 1st muon
  Dimut1_muIdx.clear();
  Dimut1_dxy.clear();
  Dimut1_dz.clear();

  // 2nd muon
  Dimut2_muIdx.clear();
  Dimut2_dxy.clear();
  Dimut2_dz.clear();

  // Dimu
  Dimu_pt.clear();
  Dimu_eta.clear();
  Dimu_phi.clear();
  Dimu_rap.clear();
  Dimu_mass.clear();
  Dimu_charge.clear();
  Dimu_simIdx.clear();
  Dimu_vtxIdx.clear();
  Dimu_chi2.clear();
  Dimu_dlxy.clear();
  Dimu_dlxyErr.clear();
  Dimu_dlxySig.clear();
  Dimu_cosphixy.clear();
  Dimu_dl.clear();
  Dimu_dlErr.clear();
  Dimu_dlSig.clear();
  Dimu_cosphi.clear();
  Dimu_ptfrac.clear();
  Dimu_x.clear();
  Dimu_y.clear();
  Dimu_z.clear();
  Dimu_Covxx.clear();
  Dimu_Covyx.clear();
  Dimu_Covzx.clear();
  Dimu_Covyy.clear();
  Dimu_Covzy.clear();
  Dimu_Covzz.clear();


  // declare and initialize variables that you want to use for Dimuon only 
  float vcmu[3] = {0.};
  float vcDimu[3] = {0.};
  float vcdimuvtx[3] = {0.};

  int simidDimu = -1; 
  int vtxidDimu = -1; 
 
  float sumptdimu = 0;
  float zDimu = 0.;

  float rapDimu = 999.;
  float chi2Dimu = 999.;
  float dlxyDimu = 0.;
  float dlxyerrDimu = 0.;
  float dlxysigDimu = 0.;
  float cosphixyDimu = 0.;
  float dlDimu = 0.;
  float dlerrDimu = 0;
  float dlsigDimu = 0; 
  float cosphiDimu = 0;
  float mDimu = 0.;

  vector<TransientTrack> mytracksforDimu;
  TransientVertex myDimuVertex;

  int imuvtx =0;
  float mudistzmin = 999.;
  int imuvtxmin = -1;

  // tlorentzvector
  TLorentzVector p4mu1;  
  p4mu1.SetPtEtaPhiE(0., 0., 0., 0.);
  TLorentzVector p4mu2;  
  p4mu2.SetPtEtaPhiE(0., 0., 0., 0.);
  TLorentzVector p4Dimu;  
  p4Dimu.SetPtEtaPhiE(0., 0., 0., 0.);

  reco::VertexCollection::const_iterator iteForDimu;
  reco::VertexCollection::const_iterator itemumin;


////////////////////////////////////////////////////////////////////////////
// List of cuts or selections applied to Dimuon candidates                //
////////////////////////////////////////////////////////////////////////////
  UInt_t muminDimu = 2;     // minimum number of tracks for Dimuon
  float detamax = 3.; 
  // WorkBook: For matching, vertices & track origins should be separated by 
  // at most 1 mm and at most 3 sigma; 
  // here we use 1 mm only throughout, from a previous empirical study,
  // and since we want to associate also close secondaries
  // *** this does not seem to work for secondary vertices! ***
  // -> increase to 5 mm! 
  Float_t vdimu_xymax = 0.5; // maximum track distance in xy for mu1 and mu2 
  Float_t vdimu_zmax = 0.5; // maximum track distance in z for mu1 and mu2
  Float_t vdimutrksum_xymax = 0.5; // maximum vertex distance in xy for ptfrac
  Float_t vdimutrksum_zmax = 0.5;  // maximum vertex distance in z for ptfrac
  Float_t muvdistzmax = 0.5; // maximum distance from primary vertex
  Float_t mucosphimin = -1.; // minimum cosphi of dimuon (no cut)

  // *** need to define parameter and protect! ***
  uint imupoint[1000];

  // check for muon collection with at least 2 entries 
  if (nMuon >= muminDimu) {
    
    // create dimuon pair from each stored (unlike sign?) muon combination
    // with delta eta < 3 (exclude forward/backward configurations)
    //
    // create a list of prestored pointers

    // preset with -1 
    for (uint imu = 0; imu<b4_nMuon; ++imu) {
      imupoint[imu] = -1;
    }  
    // fill with proper pointers (MuId contains place in original list)
    for (uint imu = 0; imu<nMuon; ++imu) {
      imupoint[Muon_Id[imu]] = imu;
    }  

    // cout << "1st MUON LOOP " << nMuon << endl;

    //---------- 1st loop over muons -----------//
  int muonid1 = -1;
#ifndef miniAOD
  for (reco::MuonCollection::const_iterator recoMuon1 = muons->begin(); recoMuon1 != muons->end(); ++recoMuon1) {
#endif
#ifdef miniAOD
  for (pat::MuonCollection::const_iterator recoMuon1 = muons->begin(); recoMuon1 != muons->end(); ++recoMuon1) {
#endif
    ++muonid1;
    int muindex1 = imupoint[muonid1];
    if (muindex1 < 0) continue;
    float mu1dxy = 0;
    float mu1dz = 0;

    // put cuts on first muon here (none for the time being)

    // cout << "2nd MUON LOOP " << muonid1 << endl;

    //------ 2nd loop over muons  ------------//
    // loop over each pair only once;
    int muonid2 = muonid1;
#ifndef miniAOD
    for (reco::MuonCollection::const_iterator recoMuon2 = recoMuon1; recoMuon2 != muons->end(); ++recoMuon2) {
#endif
#ifdef miniAOD
    for (pat::MuonCollection::const_iterator recoMuon2 = recoMuon1; recoMuon2 != muons->end(); ++recoMuon2) {
#endif
      // don't fit track with itself
      if (recoMuon1 == recoMuon2) continue;
      ++muonid2;
      int muindex2 = imupoint[muonid2];
      if (muindex2 < 0) continue;
      float mu2dxy = 0;
      float mu2dz = 0;
          
      // apply cuts
      //cout << "check FW " << muonid2 << endl; 
      // don't accept extreme forward backward pairs
      if (abs(Muon_eta[muindex1]-Muon_eta[muindex2]) > detamax) continue;
      //cout << "check sign" << endl; 
      // accept only unlike sign pairs -> now accept also like sign
      // if (Muon_charge[muindex1]*Muon_charge[muindex2]>0) continue;
      int chargeDimu = Muon_charge[muindex1] + Muon_charge[muindex2]; 

      //cout << "check origin" << endl; 
      // calculate track origin separations in x, y and z
      vcmu[0] = abs(recoMuon1->vx() - recoMuon2->vx());
      vcmu[1] = abs(recoMuon1->vy() - recoMuon2->vy());
      vcmu[2] = abs(recoMuon1->vz() - recoMuon2->vz());

      // check they are separated at most by a radius of 1mm
      // in x and y and a distance of 0.1 in z
      // *** this seems to be too tight! ***
      if (sqrt ((vcmu[0] * vcmu[0]) + (vcmu[1] * vcmu[1])) > vdimu_xymax || vcmu[2] > vdimu_zmax) continue;
		  
      // build Dimu track collection for revertexing
      mytracksforDimu.clear();
      // build transient tracks from relevant (mini)AOD tracks
      TransientTrack  transientTrack1 = theB->build( recoMuon1->bestTrack() );
      mytracksforDimu.push_back(transientTrack1);
      TransientTrack  transientTrack2 = theB->build( recoMuon2->bestTrack() );
      mytracksforDimu.push_back(transientTrack2); 

      //cout << "dimuon fit " << muonid1 << " " << muonid2 << endl;

//////////////////////////////////////////////////////////////////////////////
///////////////// Do the dimuon vertex refit /////////////////////////////////
      // allow track parameter reevaluation 
      // (Twiki: currently only the KalmanVertex can do this)
      KalmanVertexFitter  theFitter2m(true);
      // do the fit
      myDimuVertex = theFitter2m.vertex(mytracksforDimu);  
      if (!myDimuVertex.isValid()) continue;

      //cout << "after dimuon fit " << muonid1 << " " << muonid2 << endl;

      // keep this combination
      // store pointers ...

      // vertex parameters
      float ndof = myDimuVertex.degreesOfFreedom();
      // for Adaptive Vertex Fitter this sum of weights*2-3
      // for Kalman Vertex Fitter this is 1 for two tracks?
      chi2Dimu = myDimuVertex.totalChiSquared();
      // confidence level
      float CL=999.;
      if (ndof >= 1) {
        CL = TMath::Prob(chi2Dimu,(int)ndof); 
      }
      // cut on chisquared (Adaptive) or CL (Kalman) here?
      if (CL<0.01) continue;

      double Dimux = 0, Dimuy=0, Dimuz=0;
      double cov_sv[3][3], cov_pv[3][3], deriv[3];
      double dlerr = 0, dlxyerr=0;

      // good vertex, get position
      Dimux = myDimuVertex.position().x();
      Dimuy = myDimuVertex.position().y();
      Dimuz = myDimuVertex.position().z();

      // track parameters
      if (myDimuVertex.hasRefittedTracks()){
        // refitted tracks are TransientTracks
        // get updated momenta of refitted tracks
        // and recalculate muon and dimuon quantities 
        int itcount = 0;
        vector<TransientTrack> tracks = myDimuVertex.refittedTracks();
        for (vector<TransientTrack>::const_iterator trackIt = tracks.begin(); trackIt != tracks.end(); ++trackIt) {
          ++itcount;
          const Track & track = trackIt->track();
          if (itcount==1) {
            p4mu1.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), mumass);
          }
	  else if (itcount==2) { 
            p4mu2.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), mumass);
          }		              
        } //for

	// get four-vector, and from this dimuon pt, eta, phi, and mass 
	p4Dimu = p4mu1 + p4mu2;
        // investigate cut on momentum asymmetry (costheta*)?
        // (instead of Delta eta above)
	mDimu = p4Dimu.M();

        //cout << "dimuon mass: " << mDimu << endl;

        //cout << "Track pt: " << track.pt() << endl;
        //cout << recmuon1->pt() << " " << recmuon2->pt() << endl;
      } //if 

      // roughly locate Dimu 'vertex' using average coordinate of tracks
      vcDimu[0] = 0.5 * (recoMuon1->vx() + recoMuon2->vx());
      vcDimu[1] = 0.5 * (recoMuon1->vy() + recoMuon2->vy());
      vcDimu[2] = 0.5 * (recoMuon1->vz() + recoMuon2->vz());

////////////////////// find best (closest) primary vertex /////////////////////

      // loop over primary vertices
      // imuvtx sets vertex index = unique id
      imuvtx=-1;
      mudistzmin=999.;
      imuvtxmin=-1;
      vtxidDimu=-1;

      for (reco::VertexCollection::const_iterator ite = Primvertex->begin(); ite != Primvertex->end(); ++ite) {
        ++imuvtx;

        // get vertex distance in z
        float distz = fabs(ite->z()-Dimuz);
        // skip if too far away
        if (distz > muvdistzmax) continue;
        // should we apply quality cuts?
        if (ite->isFake() ||ite->ndof()<2 || fabs(ite->z()-Bsp_z)>20.) continue;
        //if (ite->isFake() || !ite->isValid() ||
        //    ite->ndof()<2 || fabs(ite->z()-Bsp_z)>20.) continue;
        // isFake means beam spot, isValid means fit converged
        // and result is within tracker boundaries (always true?), 
        // the last two are essentially no cut;
        // (ndof>4 for PVtx_isGood means >2 tracks)
        // use PVtx variables to cut harder offline if wished

        // save for later if closer than closest good vertex so far
        if (distz < mudistzmin) {
          mudistzmin = distz;
          imuvtxmin = imuvtx;
          itemumin = ite;
	  //cout << "imuvtxmin " << imuvtxmin << " distz " << distz << endl;
        }

        // find out how to check track-vertex association ...
      }

      // no vertex was directly associated to either track,
      // use closest if found
      if (imuvtxmin >-1) {
        iteForDimu = itemumin;
        vtxidDimu = imuvtxmin;
      }
      else {
        // no vertex associated, skip candidate
        // (beware of exclusive J/psis, do they always make a vertex?)
        continue;
      }
      // for 'direct' candudates, should the dimuon vertex be refitted 
      // with vertex and/or beamspot constraint? 


      // found 'best' primary for this vertex
      // get position
      double xprim = iteForDimu->position().x(); 
      double yprim = iteForDimu->position().y(); 
      double zprim = iteForDimu->position().z();
      // calculate distance of tracks to primary
      // set primary vertex reference point pvdimu
      math::XYZPoint pvdimu(iteForDimu->position());
      mu1dxy = recoMuon1->bestTrack()->dxy(pvdimu); 
      mu1dz = recoMuon1->bestTrack()->dz(pvdimu); 
      mu2dxy = recoMuon2->bestTrack()->dxy(pvdimu); 
      mu2dz = recoMuon2->bestTrack()->dz(pvdimu); 
      // *** rather call it1->dxy(iteForD0)? (see mu_dxy) ***
      // t1dxy = sqrt((it1->vx()-xprim)*(it1->vx()-xprim)
      //	 + (it1->vy()-yprim)*(it1->vy()-yprim));
      // t1dz = it1->vz()-zprim; 
      // t2dxy = sqrt((it2->vx()-xprim)*(it2->vx()-xprim)
      //	 + (it2->vy()-yprim)*(it2->vy()-yprim));
      // t2dz = it2->vz()-zprim; 
      // calculate distance of Dimu to primary
      double vdx = Dimux - xprim;
      double vdy = Dimuy - yprim;
      double vdz = Dimuz - zprim; 
      double dist = sqrt(pow(vdx,2) + pow(vdy,2) + pow(vdz,2));
      double distxy = sqrt(pow(vdx,2) + pow(vdy,2));
      // get direction vector 
      double Dpx = p4Dimu.Px(); 
      double Dpy = p4Dimu.Py();
      double Dpz = p4Dimu.Pz();
      double Dp  = p4Dimu.P();
      double Dpt  = p4Dimu.Pt();
      // calculate 3D and 2D decay length
      dlDimu = (Dpx*vdx + Dpy*vdy + Dpz*vdz)/Dp;
      dlxyDimu = (Dpx*vdx + Dpy*vdy)/Dpt;
      // calculate cosphi
      cosphiDimu = dlDimu/dist;
      cosphixyDimu = dlxyDimu/distxy;
      // cut on cosphi (0.99)
      dlerrDimu=0;
      dlsigDimu=0;
      dlxyerrDimu=0;
      dlxysigDimu=0;
      if (cosphiDimu > mucosphimin) {
        // calculate significance in space
        // transient Dimu vertex
	cov_sv[0][0] = myDimuVertex.positionError().cxx();
	cov_sv[1][0] = myDimuVertex.positionError().cyx();
	cov_sv[2][0] = myDimuVertex.positionError().czx();
        cov_sv[0][1] = cov_sv[1][0];
	cov_sv[1][1] = myDimuVertex.positionError().cyy();
	cov_sv[2][1] = myDimuVertex.positionError().czy();
        cov_sv[0][2] = cov_sv[2][0];
        cov_sv[1][2] = cov_sv[2][1];
	cov_sv[2][2] = myDimuVertex.positionError().czz();
        // primary vertex
	cov_pv[0][0] = iteForDimu->covariance(0,0);
	cov_pv[1][0] = iteForDimu->covariance(1,0);
 	cov_pv[2][0] = iteForDimu->covariance(2,0);
        cov_pv[0][1] = cov_pv[1][0];
	cov_pv[1][1] = iteForDimu->covariance(1,1);
	cov_pv[2][1] = iteForDimu->covariance(2,1);
        cov_pv[0][2] = cov_pv[2][0];
        cov_pv[1][2] = cov_pv[2][1];
	cov_pv[2][2] = iteForDimu->covariance(2,2);
        // distance direction unit vector
        deriv[0] = vdx/dist;
        deriv[1] = vdy/dist;
        deriv[2] = vdz/dist;
        // decay length error
        dlerr=0;
        for (int m=0; m<3; ++m){
          for (int n=0; n<3; ++n){
            dlerr += deriv[m]*deriv[n]*(cov_pv[m][n]+cov_sv[m][n]);
	  } // end for
	} // end for
	dlerrDimu = sqrt(dlerr);
	// decay length significance
        dlsigDimu = dlDimu/dlerrDimu;

        // decay length error in xy
        dlxyerr=0;
        for (int m=0; m<2; ++m){
          for (int n=0; n<2; ++n){
            dlxyerr += deriv[m]*deriv[n]*(cov_pv[m][n]+cov_sv[m][n]);
	  } // end for
	} // end for
	dlxyerrDimu = sqrt(dlxyerr);
	// decay length significance
        dlxysigDimu = dlxyDimu/dlxyerrDimu;
     
      } // end of cosphi cut

      // calculate dimuon rapidity
      //rapDimu = log((sqrt(p4Dimu.M()*p4Dimu.M()+p4Dimu.P()*p4Dimu.P())+p4Dimu.Pz())/(sqrt(mD0Actual*mD0Actual+p4D012.P()*p4D012.P())-p4D012.Pz()))/2.;
      rapDimu = log((p4Dimu.E()+p4Dimu.Pz())/(p4Dimu.E()-p4Dimu.Pz()))/2.;

      // calculate Dimuon momentum fraction
      sumptdimu = 0;
      // loop over all tracks
#ifndef miniAOD
      // loop over AOD track collection
      for (reco::TrackCollection::const_iterator itSum1 = tracks->begin(); itSum1 != tracks->end(); ++itSum1) {
#endif
#ifdef miniAOD
      // loop over miniAOD PackedCandidates with full track info
      // (some low pt tracks will be missed)
      for (pat::PackedCandidateCollection::const_iterator icSum1 = tracks->begin(); icSum1 != tracks->end(); ++icSum1) {
	if (!(icSum1->hasTrackDetails())) continue;
        auto itSum1 = icSum1->bestTrack();
        if (itSum1 == nullptr) continue; 
#endif
	// check track origins
	vcdimuvtx[0] = abs(vcDimu[0] - itSum1->vx());
	vcdimuvtx[1] = abs(vcDimu[1] - itSum1->vy());
	vcdimuvtx[2] = abs(vcDimu[2] - itSum1->vz());
		
	if ((sqrt(vcdimuvtx[0] * vcdimuvtx[0]) + (vcdimuvtx[1] * vcdimuvtx[1])) < vdimutrksum_xymax && vcdimuvtx[2] < vdimutrksum_zmax) {
							
	  sumptdimu += abs(itSum1->pt()); // sum pt for all tracks
	} // end of Sumpt vertex check
      } // end of itSum loop

      // will be slightly different for AOD and miniAOD
      // use bestTrack pt since innerTrack pt does not always exist
      // so this variable can technically be > 1
      //zDimu = (recoMuon1->innerTrack()->pt() + recoMuon2->innerTrack()->pt()) / sumptdimu;
      zDimu = (recoMuon1->bestTrack()->pt() + recoMuon2->bestTrack()->pt()) / sumptdimu;


      // fill ntuple variables
      if (Dimu_pt.size() < nReserve_Dimu) {

        Dimut1_muIdx.push_back(muindex1);
        Dimut1_dxy.push_back(mu1dxy);
        Dimut1_dz.push_back(mu1dz);
        Dimut2_muIdx.push_back(muindex2);
        Dimut2_dxy.push_back(mu2dxy);
        Dimut2_dz.push_back(mu2dz);

        // Dimu
        Dimu_pt.push_back(p4Dimu.Pt());
        Dimu_eta.push_back(p4Dimu.Eta());
        Dimu_phi.push_back(p4Dimu.Phi());
        Dimu_rap.push_back(rapDimu);
        Dimu_mass.push_back(mDimu);
        Dimu_charge.push_back(chargeDimu);
        Dimu_simIdx.push_back(simidDimu);
        Dimu_vtxIdx.push_back(vtxidDimu);
        Dimu_chi2.push_back(chi2Dimu);
        Dimu_dlxy.push_back(dlxyDimu);
        Dimu_dlxyErr.push_back(dlxyerrDimu);
        Dimu_dlxySig.push_back(dlxysigDimu);
        Dimu_cosphixy.push_back(cosphixyDimu);
        Dimu_dl.push_back(dlDimu);
        Dimu_dlErr.push_back(dlerrDimu);
        Dimu_dlSig.push_back(dlsigDimu);
	Dimu_cosphi.push_back(cosphiDimu);
	Dimu_ptfrac.push_back(zDimu);
	Dimu_x.push_back(Dimux);
	Dimu_y.push_back(Dimuy);
	Dimu_z.push_back(Dimuz);
	Dimu_Covxx.push_back(cov_sv[0][0]);
	Dimu_Covyx.push_back(cov_sv[1][0]);
	Dimu_Covzx.push_back(cov_sv[2][0]);
	Dimu_Covyy.push_back(cov_sv[1][1]);
	Dimu_Covzy.push_back(cov_sv[2][1]);
	Dimu_Covzz.push_back(cov_sv[2][2]);
      }
      else {cout << "WARNING!!!!! NO. OF Dimu IS MORE THAN YOUR RESERVED NO.!!!!!!" << endl;}
      nDimu = Dimu_pt.size();

    } // for muon 2
  } // for muon 1 
} // if nMuon



//////////////////////////////////////////////////////////////////////////////
////////////////////////////// Muon Collection End ///////////////////////////
//////////////////////////////////////////////////////////////////////////////

    //cout << "hello electron" << endl; 

//////////////////////////////////////////////////////////////////////////////
//////////// Electron, Photon, Tau, MET and Jet collections //////////////////
//////////////////////////////////////////////////////////////////////////////

// Electrons
//  Handle<GsfElectronCollection> electrons;
//  iEvent.getByLabel(InputTag("gsfElectrons"), electrons);

// number of electrons already Ok, but order seemingly not -> reorder
  value_el_n = 0;
  Electron_nNano = 0;
  // H4lepton
  misshits = -9999;
  IP3d_e = -9999.;
  ErrIP3d_e = -9999.;
  SIP3d_e = -9999.;
  relPFIso_e = -9999.;
  
  //"official" cut
  float el_min_pt = 5;
  // NanoAODplus cut
  if (nanoext) el_min_pt = 0;

  // for (auto it = electrons->begin(); it != electrons->end(); ++it) {
#ifndef miniAOD
  for (reco::GsfElectronCollection::const_iterator it = electrons->begin(); it != electrons->end(); ++it) {
#endif
#ifdef miniAOD
  for (pat::ElectronCollection::const_iterator it = electrons->begin(); it != electrons->end(); ++it) {
#endif

    // H4lepton
    //**************************** Part for H->4l start ******************************//
    // This section is dedicated to plot histograms of some vars before cuts
    // might be deleted after H->4l example comparison finished
    // (beware, some of the variables defined here are also used elsewere)

#ifndef CMSSW42X 
    // satisfies isPFcand   
    if (it->passingPflowPreselection())
      {
#endif
#ifndef CMSSW7plus
	// I need to define some vars before cut
        // see https://twiki.cern.ch/twiki/bin/view/CMSPublic/EgammaPublicData
	misshits = ((it->gsfTrack())->trackerExpectedHitsInner()).numberOfHits();
#endif
	IP3d_e = sqrt ( (it->gsfTrack()->dxy(pv) * it->gsfTrack()->dxy(pv)) + (it->gsfTrack()->dz(pv) * it->gsfTrack()->dz(pv)) );
	ErrIP3d_e = sqrt ( (it->gsfTrack()->d0Error() * it->gsfTrack()->d0Error()) + (it->gsfTrack()->dzError() * it->gsfTrack()->dzError()) );
	SIP3d_e = IP3d_e / ErrIP3d_e;
    
#ifdef CMSSW53X    
	relPFIso_e = ((it->pfIsolationVariables()).chargedHadronIso +
		      (it->pfIsolationVariables()).neutralHadronIso +
		      (it->pfIsolationVariables()).photonIso) /it->pt();
#endif
    
	// save hist for all vars b4 cut in HiggsDemoAnalyzer
	h_SIP3d_e_b4->Fill(SIP3d_e);
	h_p_e->Fill(it->p());
	h_et_e->Fill(it->et());
	h_pt_e_b4->Fill(it->pt());
	h_eta_e_b4->Fill(it->eta());
	h_phi_e->Fill(it->phi());
	h_sc_eta->Fill((it->superCluster())->eta());
	h_sc_rawE->Fill(std::abs((it->superCluster())->rawEnergy()));
	h_misshite->Fill(misshits);   
	h_relPFIso_e->Fill(relPFIso_e);
	h_relPFIso_pt_e->Fill(relPFIso_e, it->pt());
	h_dxy_e->Fill((it->gsfTrack())->dxy(pv));
#ifndef CMSSW42X    
      }
#endif
    //**************************** Part for H->4l end ********************************//     
    
    value_el_isNano[value_el_n] = false;
    if (it->pt() > el_min_pt) {
      if (it->pt() > 5) {
        value_el_isNano[value_el_n] = true;
        ++Electron_nNano;
      }
      value_el_pt[value_el_n] = it->pt();
      value_el_eta[value_el_n] = it->eta();
      value_el_phi[value_el_n] = it->phi();
      value_el_charge[value_el_n] = it->charge();
      // nuha
      value_el_tightCharge[value_el_n] = it->isGsfCtfScPixChargeConsistent() + it->isGsfScPixChargeConsistent();
      
      // the following was validated by comparison to the Run 2 nanoAOD
      value_el_mass[value_el_n] = it->mass();
      //value_el_mass[value_el_n] = emass;

      // need to define this for other CMSSW versions
#ifdef CMSSW42X
      // non pf based approximation, from 
      // https://twiki.cern.ch/twiki/bin/view/CMSPublic/EgammaPublicData
      // improve?
      if (it->isEB()) {
        // better to cut explicitly on eta of supercluster?
        value_el_pfreliso03all[value_el_n] =
	  (it->dr03TkSumPt() + max(0.,it->dr03EcalRecHitSumEt() -1.)
           + it->dr03HcalTowerSumEt() )/it->pt();
      }
      else {
        value_el_pfreliso03all[value_el_n] =
	  (it->dr03TkSumPt() + it->dr03EcalRecHitSumEt()
           + it->dr03HcalTowerSumEt() )/it->pt();
      }
      value_el_pfreliso03chg[value_el_n] = it->dr03TkSumPt()/it->pt();
#endif

#ifdef CMSSW53X
      // see https://twiki.cern.ch/twiki/bin/view/CMS/EgammaPFBasedIsolation
      // verify whether the isPFcand condition is needed here (in Run 2 it is not)
      if (it->passingPflowPreselection()) {
        auto iso03 = it->pfIsolationVariables();
        value_el_pfreliso03all[value_el_n] =
            (iso03.chargedHadronIso + iso03.neutralHadronIso + iso03.photonIso)/it->pt();
        value_el_pfreliso03chg[value_el_n] = iso03.chargedHadronIso/it->pt();
      } 
      else {
        value_el_pfreliso03all[value_el_n] = -999.;
        value_el_pfreliso03chg[value_el_n] = -999.;
      }
#endif

#ifdef CMSSW7plus
      // *** to be fixed!!! ***
      // from https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2017
      // value_el_pfreliso03all[value_el_n] = (it->ecalPFClusterIso()+it->hcalPFClusterIso())/it->pt();
      // do *not* require isPFcand at this stage (validated with official nanoAOD)
      auto iso03 = it->pfIsolationVariables();
      value_el_pfreliso03all[value_el_n] = (iso03.sumChargedHadronPt + iso03.sumNeutralHadronEt + iso03.sumPhotonEt) / it->pt();
      value_el_pfreliso03chg[value_el_n] = iso03.sumChargedHadronPt / it->pt();
#endif

      // https://twiki.cern.ch/twiki/bin/view/CMS/EgammaCutBasedIdentification
      // nuha
      value_el_dr03TkSumPtOld[value_el_n] = it->dr03TkSumPt();
      value_el_dr03TkSumPt[value_el_n] = (it->pt() > 35.) ? it->dr03TkSumPt() : 0.;      
      value_el_dr03EcalRecHitSumEtOld[value_el_n] = it->dr03EcalRecHitSumEt();
      value_el_dr03EcalRecHitSumEt[value_el_n] = (it->pt() > 35.) ? it->dr03EcalRecHitSumEt() : 0.;
      
      // combine the next two variables?
      value_el_dr03HcalTowerSumEt[value_el_n] = it->dr03HcalTowerSumEt();
      // for the moment fill with same input (4_2_8)
     // nuha  *** A.G., the next might need to be reactivated ***
      //value_el_dr03HcalDepth1TowerSumEt[value_el_n] = it->dr03HcalTowerSumEt();
      //#ifdef CMSSW7plus     
      value_el_dr03HcalDepth1TowerSumEtOld[value_el_n] = it->dr03HcalDepth1TowerSumEt();
      value_el_dr03HcalDepth1TowerSumEt[value_el_n] = (it->pt() > 35.) ? it->dr03HcalDepth1TowerSumEt() : 0;      
      //#endif
            
      // do the next three contain redundancy?
      value_el_isEB[value_el_n] = it->isEB();
      value_el_isEE[value_el_n] = it->isEE();
      value_el_SCeta[value_el_n] = (it->superCluster())->eta();

      value_el_isPFcand[value_el_n] = 0;
#ifndef CMSSW42X
      value_el_isPFcand[value_el_n] = it->passingPflowPreselection();
#endif    
      // *** need to fix for Run 2 ***
#ifndef CMSSW7plus
      // https://twiki.cern.ch/twiki/bin/view/CMS/EgammaCutBasedIdentification
      value_el_lostHits[value_el_n] = ((it->gsfTrack())->trackerExpectedHitsInner()).numberOfHits();
      // nuha  (seems to yield same result, although not documented)
      //value_el_lostHits[value_el_n] = ((it->gsfTrack())->trackerExpectedHitsInner()).numberOfLostHits();
#endif
#ifdef CMSSW7plus
      value_el_lostHits[value_el_n] = ((it->gsfTrack())->hitPattern()).numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);
#endif
#ifdef CMSSW7plus
      //this is how it was before
      //value_el_lostHits[value_el_n] = 0;

      // only available from 94X onwards (do we need further ifdef?)
      // https://github.com/cms-sw/cmssw/blob/CMSSW_9_4_X/RecoEgamma/ElectronIdentification/plugins/cuts/GsfEleMissingHitsCut.cc
      value_el_lostHits[value_el_n] = ((it->gsfTrack())->hitPattern()).numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);
      // for 80X, they use this:
      // value_el_lostHits[value_el_n] = ((it->gsfTrack())->hitPattern()).numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
      // https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_X/RecoEgamma/ElectronIdentification/plugins/cuts/GsfEleMissingHitsCut.cc
      // not sure which one to use, need to discuss/check, see also 
      // https://github.com/cms-sw/cmssw/blob/CMSSW_7_6_X/RecoEgamma/ElectronIdentification/plugins/cuts/GsfEleMissingHitsCut.cc
      // https://github.com/cms-sw/cmssw/blob/CMSSW_10_2_X/RecoEgamma/ElectronIdentification/plugins/cuts/GsfEleMissingHitsCut.cc
#endif

      // are the next two redundant with convVeto?
      value_el_convDist[value_el_n] = it->convDist();
      value_el_convDcot[value_el_n] = it->convDcot();
      // flag to *pass* veto (i.e. accept if 1)
      // deal with -10000 entries (not treated so far)

      // nuha
      value_el_convVetoOld[value_el_n] = it->convDist()<0.02 && it->convDcot()<0.02;
      // conversion veto definition for Run 2 and Run1
      // Copy paste the def in: https://github.com/cms-sw/cmssw/blob/CMSSW_9_4_X/PhysicsTools/PatAlgos/plugins/PATElectronProducer.cc
      // set conversion veto selection
      bool passconversionveto = false;
      if (hConversions.isValid()) {
	// this is recommended method
	passconversionveto = !ConversionTools::hasMatchedConversion( *it, *hConversions, beamSpotHandle->position());
      } else {
	// use missing hits without vertex fit method
#ifndef CMSSW7plus
	passconversionveto = it->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() < 1;
#endif
#ifdef CMSSW7plus
	passconversionveto = it->gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) < 1;
#endif	
      }
      value_el_convVeto[value_el_n] = passconversionveto;
	
      // *** recheck the following ***
      // see https://github.com/cms-sw/cmssw/blob/82e3ed4ea0b8756c412d1f56758b7698d717104e/PhysicsTools/NanoAOD/python/electrons_cff.py
      value_el_deltaEtaSC[value_el_n] = (it->superCluster()->eta()) - it->eta();
      value_el_deltaPhiSC[value_el_n] = (it->superCluster()->phi()) - it->phi();
      if (value_el_deltaPhiSC[value_el_n] > 3.1415) value_el_deltaPhiSC[value_el_n] = value_el_deltaPhiSC[value_el_n]-2.*3.1415;  
      if (value_el_deltaPhiSC[value_el_n] < -3.1415) value_el_deltaPhiSC[value_el_n] = value_el_deltaPhiSC[value_el_n]+2.*3.1415;  
      // nanoAODplus variables for Run 1
      value_el_deltaEtaSCtr[value_el_n] = it->deltaEtaSuperClusterTrackAtVtx();
      value_el_deltaPhiSCtr[value_el_n] = it->deltaPhiSuperClusterTrackAtVtx();
      // what about it->hadronicOverEM() ? see
      // https://twiki.cern.ch/twiki/bin/view/CMS/EgammaCutBasedIdentification
      value_el_hoe[value_el_n] = it->hcalOverEcal();
#ifdef CMSSW7plus
      // https://github.com/cms-sw/cmssw/blob/CMSSW_9_4_X/PhysicsTools/NanoAOD/python/electrons_cff.py#L268
      value_el_sieie[value_el_n] = it->full5x5_sigmaIetaIeta();
#endif
#ifndef CMSSW7plus
      value_el_sieie[value_el_n] = -9999.;
#endif
      // nanoAODplus variable for Run 1 (exists also for Run 2)
      value_el_sieieR1[value_el_n] = it->sigmaIetaIeta();
      // what about it->ecalEnergy()  and it->trackMomentumAtVtx().p() ? see
      // https://twiki.cern.ch/twiki/bin/view/CMS/EgammaCutBasedIdentification
      // what about it->ecalEnergy()  and it->trackMomentumAtVtx().p() ? see
      // https://twiki.cern.ch/twiki/bin/view/CMS/EgammaCutBasedIdentification
      // nuha
      value_el_eInvMinusPInvOld[value_el_n] = 1/it->p4().E() - 1/it->p4().P();
      value_el_eInvMinusPInv[value_el_n] = (1-(it->eSuperClusterOverP()) ) / it->ecalEnergy();

      // auto trk = it->gsfTrack();
      // see https://twiki.cern.ch/twiki/bin/view/CMSPublic/EgammaPublicData
      value_el_dxy[value_el_n] = it->gsfTrack()->dxy(pv);
      // value_el_dxy[value_el_n] = trk->dxy(pv) * it->charge();
      // if (it->phi()<0) value_el_dxy[value_el_n] = -value_el_dxy[value_el_n];
      if (it->phi()<0) value_el_dxy[value_el_n] = -value_el_dxy[value_el_n];
      value_el_dz[value_el_n] = it->gsfTrack()->dz(pv);
      value_el_dxyErr[value_el_n] = it->gsfTrack()->d0Error();
      value_el_dzErr[value_el_n] = it->gsfTrack()->dzError();
  
      value_el_ip3d[value_el_n] = IP3d_e;
      value_el_sip3d[value_el_n] = SIP3d_e;
      
      // for cut based selection cuts (2011 data and higher), see
      // https://twiki.cern.ch/twiki/bin/view/CMS/EgammaCutBasedIdentification
      // for 2010 data, see
      // https://twiki.cern.ch/twiki/bin/view/CMSPublic/EgammaPublicData
      bool eisVeto = true;
      bool eisLoose = true;
      bool eisMedium = true;
      bool eisTight = true;
      if (fabs(value_el_SCeta[value_el_n])<=1.479) {
      // barrel (are SCeta and isEB/isEE redundant?)
      // verify which of the two variable variants (with or without tr) 
      // should be used!
        // dEtaIn   *** recheck whether EtaSc or EtaSCtr! ***
        if (fabs(value_el_deltaEtaSC[value_el_n])>0.004) {
	  eisMedium = false;  
	  eisTight = false;  
        }
        if (fabs(value_el_deltaEtaSC[value_el_n])>0.007) {
	  eisVeto = false;  
	  eisLoose = false;  
        }
        // dPhiIn  *** recheck whether PhiSc or PhiSCtr! ***
        if (fabs(value_el_deltaPhiSC[value_el_n])>0.03) {
	  eisTight = false;
        } 
        if (fabs(value_el_deltaPhiSC[value_el_n])>0.06) {
	  eisMedium = false;
        } 
        if (fabs(value_el_deltaPhiSC[value_el_n])>0.15) {
	eisLoose = false;
        } 
        if (fabs(value_el_deltaPhiSC[value_el_n])>0.8) {
	  eisVeto = false;
        } 
        // sigmaIEtaIEta  *** recheck whether sieie or sieieR1! ***
        if (value_el_sieie[value_el_n]>0.01) {
	  eisTight = false;
	  eisMedium = false;
	  eisLoose = false;
	  eisVeto = false;
        }
        // H/E
        if (value_el_hoe[value_el_n]>0.12) {
	  eisTight = false;
	  eisMedium = false;
	  eisLoose = false;
        }
        if (value_el_hoe[value_el_n]>0.15) {
	  eisVeto = false;
        }
        // d0(vtx)
        if (fabs(value_el_dxy[value_el_n])>0.02) {
	  eisTight = false;
	  eisMedium = false;
	  eisLoose = false;
        }
        if (fabs(value_el_dxy[value_el_n])>0.04) {
	  eisVeto = false;
        }
        // dz(vtx)
        if (fabs(value_el_dz[value_el_n])>0.1) {
	  eisTight = false;
	  eisMedium = false;
        }
        if (fabs(value_el_dz[value_el_n])>0.2) {
	  eisLoose = false;
	  eisVeto = false;
        }
        // 1/E-1/p
        if (fabs(value_el_eInvMinusPInvOld[value_el_n])>0.05) { // nuha
	  eisTight = false;
	  eisMedium = false;
	  eisLoose = false;
        }
        // pf isolation/pt 
        if (value_el_pfreliso03all[value_el_n]>0.10) {
	  eisTight = false;
        }
        if (value_el_pfreliso03all[value_el_n]>0.15) {
	  eisMedium = false;
	  eisLoose = false;
	  eisVeto = false;
        }
        // conversion rejection
        // implement properly!
        if (!value_el_convVetoOld[value_el_n]) {
	  eisTight = false;
	  eisMedium = false;
	  eisLoose = false;
        }
        // missing hits
        if (value_el_lostHits[value_el_n]>0) {
	  eisTight = false;
        }
        if (value_el_lostHits[value_el_n]>1) {
	  eisMedium = false;
	  eisLoose = false;       
        }
      }
      else if (fabs(value_el_SCeta[value_el_n])<=2.5) {
      // end cap (are SCeta and isEB/isEE redundant?)
        // dEtaIn
        if (fabs(value_el_deltaEtaSC[value_el_n])>0.005) {
	  eisTight = false;  
        }
        if (fabs(value_el_deltaEtaSC[value_el_n])>0.007) {
	  eisMedium = false;  
        }
        if (fabs(value_el_deltaEtaSC[value_el_n])>0.009) {
	  eisLoose = false;  
        }
        if (fabs(value_el_deltaEtaSC[value_el_n])>0.010) {
	  eisVeto = false;  
        }
        // dPhiIn
        if (fabs(value_el_deltaPhiSC[value_el_n])>0.02) {
	  eisTight = false;
        } 
        if (fabs(value_el_deltaPhiSC[value_el_n])>0.03) {
	  eisMedium = false;
        } 
        if (fabs(value_el_deltaPhiSC[value_el_n])>0.10) {
	eisLoose = false;
        } 
        if (fabs(value_el_deltaPhiSC[value_el_n])>0.7) {
	  eisVeto = false;
        } 
        // sigmaIEtaIEta
        if (value_el_sieie[value_el_n]>0.03) {
	  eisTight = false;
	  eisMedium = false;
	  eisLoose = false;
	  eisVeto = false;
        }
        // H/E
        if (value_el_hoe[value_el_n]>0.10) {
	  eisTight = false;
	  eisMedium = false;
	  eisLoose = false;
        }
        // d0(vtx)
        if (fabs(value_el_dxy[value_el_n])>0.02) {
	  eisTight = false;
	  eisMedium = false;
	  eisLoose = false;
        }
        if (fabs(value_el_dxy[value_el_n])>0.04) {
	  eisVeto = false;
        }
        // dz(vtx)
        if (fabs(value_el_dz[value_el_n])>0.1) {
	  eisTight = false;
	  eisMedium = false;
        }
        if (fabs(value_el_dz[value_el_n])>0.2) {
	  eisLoose = false;
	  eisVeto = false;
        }
        // 1/E-1/p
        if (fabs(value_el_eInvMinusPInvOld[value_el_n])>0.05) { // nuha
	  eisTight = false;
	  eisMedium = false;
	  eisLoose = false;
        }
        // pf isolation/pt 
        if (value_el_pt[value_el_n]>20.) {
          if (value_el_pfreliso03all[value_el_n]>0.10) {
	    eisTight = false;
          }
          if (value_el_pfreliso03all[value_el_n]>0.15) {
	    eisMedium = false;
	    eisLoose = false;
	    eisVeto = false;
          }
        }
	else {
          if (value_el_pfreliso03all[value_el_n]>0.07) {
	    eisTight = false;
          }
          if (value_el_pfreliso03all[value_el_n]>0.10) {
	    eisMedium = false;
	    eisLoose = false;
          }
          if (value_el_pfreliso03all[value_el_n]>0.15) {
	    eisVeto = false;
          }
        }
        // conversion rejection
        // implement properly!
        if (value_el_convVetoOld[value_el_n]) {
	  eisTight = false;
	  eisMedium = false;
	  eisLoose = false;
        }
        // missing hits
        if (value_el_lostHits[value_el_n]>0) {
	  eisTight = false;
        }
        if (value_el_lostHits[value_el_n]>1) {
	  eisMedium = false;
	  eisLoose = false;       
        }
      }
      else {
        // this should not happen
        cout << "NanoAnalyzer: Electron cluster outside bound" << endl; 
	eisTight = false;
	eisMedium = false;
	eisLoose = false;
	eisVeto = false;
      }
      // set variable
      if (eisTight) value_el_cutBased[value_el_n] = 4;
      else if (eisMedium) value_el_cutBased[value_el_n] = 3;
      else if (eisLoose) value_el_cutBased[value_el_n] = 2;
      else if (eisVeto) value_el_cutBased[value_el_n] = 1;
      else value_el_cutBased[value_el_n] = 0;

      ++value_el_n;
      if (int(value_el_n) >= max_el) {
        cout << "NanoAnalyzer: max_el exceeded" << endl;
        continue;
      }

    } // pt 
  } // for 

  //  for 2010 data and nanoAOD, implement
  // Electron_dr03TkSumPt
  // Electron_dr03EcalRecHitSumEt
  // Electron_dr03HcalTowerSumEt  (2010)
  // Electron_dr03HcalDepth1TowerSumEt (nanoAOD)
  // Electron_isEB (2010)
  // Electron_isEE (2010)
  // Electron_lostHits (=misshits)
  // Electron_convDist (2010)
  // Electron_convDcot (2010)
  // Electron_convVeto 
  // Electron_cutBased (nanoAOD)
  // Electron_deltaEtaSC (it->superCluster()->eta()) - it->eta() (nanoAOD)
  // Electron_deltaPhiSC (it->superCluster()->phi()) - it->phi() (nanoAOD)
  // Electron_deltaEtaSCtr (it->deltaEtaSuperClusterTrackAtVtx())  (2010)
  // Electron_deltaPhiSCtr (it->deltaPhiSuperClusterTrackAtVtx())  (2010)
  // Electron_hoe (it->hcalOverEcal())
  // Electron_sieie (it->sigmaIetaIeta())
  // Electron_SCeta
  // Electron_eInvMinusPInv

  //  for 2011
  // Electron_pfRelIso03_chg

  //  for Run 2 
  // Electron_pfRelIso03_all
  // Electron_pfRelIso03_chg

  //  Sort electrons according to pt
  for (int i1=0; i1<int(value_el_n)-1; ++i1) {
    for (int i2= i1+1; i2<int(value_el_n); ++i2) {
      if (value_el_pt[i1]<value_el_pt[i2]) {

        Float_t temp = value_el_pt[i1];
        value_el_pt[i1] = value_el_pt[i2];
        value_el_pt[i2] = temp;    
        temp = value_el_eta[i1];
        value_el_eta[i1] = value_el_eta[i2];
        value_el_eta[i2] = temp;    
        temp = value_el_phi[i1];
        value_el_phi[i1] = value_el_phi[i2];
        value_el_phi[i2] = temp;    
        Int_t tempi = value_el_charge[i1];
        value_el_charge[i1] = value_el_charge[i2];
        value_el_charge[i2] = tempi;    
        temp = value_el_mass[i1];
        value_el_mass[i1] = value_el_mass[i2];
        value_el_mass[i2] = temp;
	// nuha
        tempi = value_el_tightCharge[i1];
        value_el_tightCharge[i1] = value_el_tightCharge[i2];
        value_el_tightCharge[i2] = tempi;

        temp = value_el_pfreliso03all[i1];
        value_el_pfreliso03all[i1] = value_el_pfreliso03all[i2];
        value_el_pfreliso03all[i2] = temp;    
        temp = value_el_pfreliso03chg[i1];
        value_el_pfreliso03chg[i1] = value_el_pfreliso03chg[i2];
        value_el_pfreliso03chg[i2] = temp;

	// nuha
        temp = value_el_dr03TkSumPtOld[i1];
        value_el_dr03TkSumPtOld[i1] = value_el_dr03TkSumPtOld[i2];
        value_el_dr03TkSumPtOld[i2] = temp;
        temp = value_el_dr03TkSumPt[i1];
        value_el_dr03TkSumPt[i1] = value_el_dr03TkSumPt[i2];
        value_el_dr03TkSumPt[i2] = temp;	
        temp = value_el_dr03EcalRecHitSumEtOld[i1];
        value_el_dr03EcalRecHitSumEtOld[i1] = value_el_dr03EcalRecHitSumEtOld[i2];
        value_el_dr03EcalRecHitSumEtOld[i2] = temp;
        temp = value_el_dr03EcalRecHitSumEt[i1];
        value_el_dr03EcalRecHitSumEt[i1] = value_el_dr03EcalRecHitSumEt[i2];
        value_el_dr03EcalRecHitSumEt[i2] = temp;
		
        temp = value_el_dr03HcalTowerSumEt[i1];
        value_el_dr03HcalTowerSumEt[i1] = value_el_dr03HcalTowerSumEt[i2];
        value_el_dr03HcalTowerSumEt[i2] = temp;
	// nuha
        temp = value_el_dr03HcalDepth1TowerSumEtOld[i1];
        value_el_dr03HcalDepth1TowerSumEtOld[i1] = value_el_dr03HcalDepth1TowerSumEtOld[i2];
        value_el_dr03HcalDepth1TowerSumEtOld[i2] = temp;    
        temp = value_el_dr03HcalDepth1TowerSumEt[i1];
        value_el_dr03HcalDepth1TowerSumEt[i1] = value_el_dr03HcalDepth1TowerSumEt[i2];
        value_el_dr03HcalDepth1TowerSumEt[i2] = temp;
	
        bool tempb = value_el_isEB[i1];
        value_el_isEB[i1] = value_el_isEB[i2];
        value_el_isEB[i2] = tempb;    
        tempb = value_el_isEE[i1];
        value_el_isEE[i1] = value_el_isEE[i2];
        value_el_isEE[i2] = tempb;
	tempb = value_el_isPFcand[i1];
        value_el_isPFcand[i1] = value_el_isPFcand[i2];
        value_el_isPFcand[i2] = tempb;
	tempb = value_el_isNano[i1];
        value_el_isNano[i1] = value_el_isNano[i2];
        value_el_isNano[i2] = tempb;
        UChar_t tempU = value_el_lostHits[i1];
        value_el_lostHits[i1] = value_el_lostHits[i2];
        value_el_lostHits[i2] = tempU;    
        temp = value_el_convDist[i1];
        value_el_convDist[i1] = value_el_convDist[i2];
        value_el_convDist[i2] = temp;    
        temp = value_el_convDcot[i1];
        value_el_convDcot[i1] = value_el_convDcot[i2];
        value_el_convDcot[i2] = temp;
	// nuha
        tempb = value_el_convVetoOld[i1];
        value_el_convVetoOld[i1] = value_el_convVetoOld[i2];
        value_el_convVetoOld[i2] = tempb;
        tempb = value_el_convVeto[i1];
        value_el_convVeto[i1] = value_el_convVeto[i2];
        value_el_convVeto[i2] = tempb;
	
        temp = value_el_deltaEtaSC[i1];
        value_el_deltaEtaSC[i1] = value_el_deltaEtaSC[i2];
        value_el_deltaEtaSC[i2] = temp;    
        temp = value_el_deltaPhiSC[i1];
        value_el_deltaPhiSC[i1] = value_el_deltaPhiSC[i2];
        value_el_deltaPhiSC[i2] = temp;    
        temp = value_el_deltaEtaSCtr[i1];
        value_el_deltaEtaSCtr[i1] = value_el_deltaEtaSCtr[i2];
        value_el_deltaEtaSCtr[i2] = temp;    
        temp = value_el_deltaPhiSCtr[i1];
        value_el_deltaPhiSCtr[i1] = value_el_deltaPhiSCtr[i2];
        value_el_deltaPhiSCtr[i2] = temp;    
        temp = value_el_hoe[i1];
        value_el_hoe[i1] = value_el_hoe[i2];
        value_el_hoe[i2] = temp;    
        temp = value_el_sieie[i1];
        value_el_sieie[i1] = value_el_sieie[i2];
        value_el_sieie[i2] = temp; 
        temp = value_el_sieieR1[i1];
        value_el_sieieR1[i1] = value_el_sieieR1[i2];
        value_el_sieieR1[i2] = temp;       
        temp = value_el_SCeta[i1];
        value_el_SCeta[i1] = value_el_SCeta[i2];
        value_el_SCeta[i2] = temp;
	// nuha
        temp = value_el_eInvMinusPInvOld[i1];
        value_el_eInvMinusPInvOld[i1] = value_el_eInvMinusPInvOld[i2];
        value_el_eInvMinusPInvOld[i2] = temp;
        temp = value_el_eInvMinusPInv[i1];
        value_el_eInvMinusPInv[i1] = value_el_eInvMinusPInv[i2];
        value_el_eInvMinusPInv[i2] = temp;
	
        tempi = value_el_cutBased[i1];
        value_el_cutBased[i1] = value_el_cutBased[i2];
        value_el_cutBased[i2] = tempi;    
        //cout << value_el_cutBased[i1] << " " << value_el_cutBased[i2] << endl;

        temp = value_el_dxy[i1];
        value_el_dxy[i1] = value_el_dxy[i2];
        value_el_dxy[i2] = temp;    
        temp = value_el_dxyErr[i1];
        value_el_dxyErr[i1] = value_el_dxyErr[i2];
        value_el_dxyErr[i2] = temp;    
        temp = value_el_dz[i1];
        value_el_dz[i1] = value_el_dz[i2];
        value_el_dz[i2] = temp;    
        temp = value_el_dzErr[i1];
        value_el_dzErr[i1] = value_el_dzErr[i2];
        value_el_dzErr[i2] = temp;    

        temp = value_el_ip3d[i1];
        value_el_ip3d[i1] = value_el_ip3d[i2];
        value_el_ip3d[i2] = temp;    
        temp = value_el_sip3d[i1];
        value_el_sip3d[i1] = value_el_sip3d[i2];
        value_el_sip3d[i2] = temp;
      }
    }
  }

  //cout << "hello tau" << endl; 

  // Taus

  // to be fully implemented as indicated in 
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePFTauID#4XX
  //                                                               53X

  //Handle<PFTauCollection> taus;
  //iEvent.getByLabel(InputTag("hpsPFTauProducer"), taus);

  const float tau_min_pt = 15;
  value_tau_n = 0;
  //  for (auto it = taus->begin(); it != taus->end(); ++it) {
#ifndef miniAOD
  for (reco::PFTauCollection::const_iterator it = taus->begin(); it != taus->end(); ++it) {
#endif
#ifdef miniAOD
  for (pat::TauCollection::const_iterator it = taus->begin(); it != taus->end(); ++it) {
#endif
    if (it->pt() > tau_min_pt) {
      value_tau_pt[value_tau_n] = it->pt();
      value_tau_eta[value_tau_n] = it->eta();
      value_tau_phi[value_tau_n] = it->phi();
      value_tau_charge[value_tau_n] = it->charge();
      value_tau_mass[value_tau_n] = it->mass();
      value_tau_decaymode[value_tau_n] = it->decayMode();
#ifndef miniAOD
      value_tau_chargediso[value_tau_n] = it->isolationPFChargedHadrCandsPtSum();
      value_tau_neutraliso[value_tau_n] = it->isolationPFGammaCandsEtSum();
#endif
      ++value_tau_n;
      if (int(value_tau_n) >= max_tau) {
        cout << "NanoAnalyzer: max_tau exceeded" << endl;
        continue;
      }
    }
  }

  //cout << "hello photon" << endl; 

  // Photons
  //Handle<PhotonCollection> photons;
  //iEvent.getByLabel(InputTag("photons"), photons);

  value_ph_n = 0;
  const float ph_min_pt = 5;
  // for (auto it = photons->begin(); it != photons->end(); ++it) {
#ifndef miniAOD
  for (reco::PhotonCollection::const_iterator it = photons->begin(); it != photons->end(); ++it) {
#endif
#ifdef miniAOD
  for (pat::PhotonCollection::const_iterator it = photons->begin(); it != photons->end(); ++it) {
#endif
    if (it->pt() > ph_min_pt) {
      value_ph_pt[value_ph_n] = it->pt();
      value_ph_eta[value_ph_n] = it->eta();
      value_ph_phi[value_ph_n] = it->phi();
      value_ph_charge[value_ph_n] = it->charge();
      value_ph_mass[value_ph_n] = it->mass();
      value_ph_pfreliso03all[value_ph_n] = it->ecalRecHitSumEtConeDR03();
      ++value_ph_n;
      if (int(value_ph_n) >= max_ph) {
        cout << "NanoAnalyzer: max_ph exceeded" << endl;
        continue;
      }
    }
  }

  //cout << "hello MET" << endl; 

  // MET
  //Handle<PFMETCollection> met;
  //iEvent.getByLabel(InputTag("pfMet"), met);
  value_met_pt = met->begin()->pt();
  value_met_phi = met->begin()->phi();
  value_met_sumEt = met->begin()->sumEt();
  // avoid occasional crash on 2010 data due to singular matrix 
  if (value_met_sumEt == 0 || value_met_sumEt > 0.001)
    value_met_significance = met->begin()->significance();
  else {
    cout << "missing pt/ET: " << value_met_pt << " " << value_met_sumEt << endl;
    cout << "nanoAnalyzer: MET significance set to 0 to avoid singular matrix" << endl; 
    value_met_significance = 0.;
  }
  //cout << "hello MET 3" << endl; 
#ifdef CMSSW53X
  auto cov = met->begin()->getSignificanceMatrix();
  value_met_covxx = cov[0][0];
  value_met_covxy = cov[0][1];
  value_met_covyy = cov[1][1];
#else
  value_met_covxx = 0;
  value_met_covxy = 0;
  value_met_covyy = 0;
#endif

  //cout << "hello caloMET" << endl; 

  // caloMET
#ifndef miniAOD
  value_calomet_pt = calomet->begin()->pt();
  value_calomet_phi = calomet->begin()->phi();
  value_calomet_sumEt = calomet->begin()->sumEt();
#endif
#ifdef miniAOD
  value_calomet_pt = calomet->front().caloMETPt();
  value_calomet_phi = calomet->front().caloMETPhi();
  value_calomet_sumEt = calomet->front().caloMETSumEt();
#endif

  //cout << "hello Jet" << endl; 

  // Jets
  //Handle<PFJetCollection> jets;
  //iEvent.getByLabel(InputTag("ak5PFJets"), jets);

  const float jet_min_pt = 15;
  value_jet_n = 0;
  // for (auto it = jets->begin(); it != jets->end(); ++it) {
#ifndef miniAOD
  for (reco::PFJetCollection::const_iterator it = jets->begin(); it != jets->end(); ++it) {
#endif
#ifdef miniAOD
    for (pat::JetCollection::const_iterator it = jets->begin(); it != jets->end(); ++it) {
#endif 
    if (it->pt() > jet_min_pt) {
      value_jet_pt[value_jet_n] = it->pt();
      value_jet_eta[value_jet_n] = it->eta();
      value_jet_phi[value_jet_n] = it->phi();
      value_jet_mass[value_jet_n] = it->mass();

      int pfConstituents = it->getPFConstituents().size();
      // quark-gluon separator (not yet working)
/*      double ptD = 0, pfqt2 = 0, pfqt = 0;
      for (int ii = 0; ii < pfConstituents; ii ++) {


        double apt = 0;
        const reco::PFCandidatePtr myptr =    (it->getPFConstituent(ii));
        apt = myptr->pt();
        pfqt2 += pow(apt, 2);
        pfqt += apt;
      }
      ptD = sqrt( pfqt2/pow(pfqt, 2) );*/

/*cout<<"ptD: "<<ptD<<endl;
cout<<"area: "<<it->jetArea()<<endl;
cout<<"constituents: "<<pfConstituents<<endl;
cout<<"el mult: "<<it->electronMultiplicity()<<endl;
cout<<"mu mult: "<<it->muonMultiplicity()<<endl;

cout<<"chEmEF: "<<it->chargedEmEnergyFraction()<<endl;
cout<<"chHEF: "<<it->chargedHadronEnergyFraction()<<endl;
cout<<"neEmEF: "<<it->neutralEmEnergyFraction()<<endl;
cout<<"neHEF: "<<it->neutralHadronEnergyFraction()<<endl;

cout<<"muEF: "<<it->muonEnergyFraction()<<endl;
cout<<"---------------------------------"<<endl;*/

//
//      value_jet_ptD[value_jet_n] = ptD;
      value_jet_area[value_jet_n] = it->jetArea();
      value_jet_nConstituents[value_jet_n] = pfConstituents;
      value_jet_nElectrons[value_jet_n] = it->electronMultiplicity();
      value_jet_nMuons[value_jet_n] = it->muonMultiplicity();
      value_jet_chEmEF[value_jet_n] = it->chargedEmEnergyFraction();
      value_jet_chHEF[value_jet_n] = it->chargedHadronEnergyFraction();
      value_jet_neEmEF[value_jet_n] = it->neutralEmEnergyFraction();
      value_jet_neHEF[value_jet_n] = it->neutralHadronEnergyFraction();

      ++value_jet_n;
      if (int(value_jet_n) >= max_jet) {
        cout << "NanoAnalyzer: max_jet exceeded" << endl;
        continue;
      }
    }
  }

  // afiqaize
  // *************************************************************
  // ----------------- get custom filter -------------------------
  // ************************************************************* 

  if (!custom_flag.empty()) {
    const edm::TriggerResults &custom_result = *custom_handle;
    const edm::TriggerNames &custom_list = iEvent.triggerNames(custom_result);

    for (size_t iFlag = 0; iFlag < custom_flag.size(); ++iFlag) {
      size_t iPath = custom_list.triggerIndex(custom_flag.at(iFlag));

      if (iPath < custom_result.size() and custom_result.accept(iPath))
        custom_bit.at(iFlag) = 1;
      else 
        custom_bit.at(iFlag) = 0;
    }
  }

  //cout << "hello D meson" << endl;
    if (nanoext) {
      
//////////////////////////////////////////////////////////////////////////////
//////////////////////////// D Meson Analysis Start //////////////////////////
//////////////////////////////////////////////////////////////////////////////

  // nanoAOD extension //

  // D0 meson object building strategy: (not yet fully implemented)
  // handles: track pt, eta, quality
  //          primary vertex/beam spot selection and quality
  //          secondary vertex track preselection (kinematics, high purity?, uncertainty?)
  //          secondary vertex fit + quality (chi2, prob/CL) (Adaptive, Kalman, DAF)
  //          distance, pointing and pT fraction between primary and secondary
  //          kinematic cuts on D0 and Dstar

  // Preparation: 
  // Loop over all primary vertices, then over all tracks associated to each 
  // vertex and build sum pt for each. 
  // Store a vertex counter and id assignment on the way. 
  // Classify one of the primary vertices as 'muon trigger vertex' (if any),
  // and tag it. If the situation is ambigous, select a main primary vertex 
  // anyway, and mark all other vertices as not to be used (possible trigger 
  // bias). For the main primary vertex only, recalculate sum pt excluding 
  // the good muon(s).
  // Generally exclude (flag) vertices which are not the trigger vertex and 
  // which are within 1 cm of the trigger vertex. 
  // Store list of vertices in ntuple.

  // On the way, fill reference histograms for tracks and vertices.

  // Build transient track list from Reco track list.
  // First loop over all possible pairs of Reco tracks with loose 
  // 'distance cut' and build two invariant masses assuming K/pi or pi/K.
  // If either of them is within a loose D0 mass window (cut), loop over all 
  // primary vertices and find the closest one in z. 
  // Calculate zD0 from track pair and vertex sum pt, apply cut (if any).

  // If unlike sign: 
  // Create a list of the two corresponding transient tracks and apply the
  // Adaptive Vertex Fitter or the Kalman Vertex Fitter.
  // Both perform about equally.
  // If the vertex fit quality is good, recalculate all kinematic D0 quantities 
  // from the refitted tracks. 
  // Reapply the relevant cuts more tightly. 
  // Calculate decay length, uncertainty, significance, and cosphi w.r.t. the 
  // chosen vertex. Check for ambigous close by alternative vertices.
  // Apply cuts. 

  // Apply final cuts to store the D0 candidate in the ntuple, without stopping
  // loop if not satisfied. Store pointer to vertex.
  // Consider to store only one D0 candidate within mass window per vertex.
  // Deal with the Kpi/piK ambiguity.  
  
  // D* object building strategy: (not yet fully implemented)

  // Apply separate set of preselection cuts for D0 from D*.

  // If like sign:
  // First make some rough check whether this can lead to a Dstar candidate,
  // i.e. wheter there are suitable slow pion candidates.
  // If there are, proceed with revertexing as above.
  // If there aren't, go to next pair.

  // Loop over all primary tracks from relevant vertex for slow pion, and 
  // select appropriate D0 combination (Kpi or piK) according to charge.
  // Check whether also some non-primary tracks should be considered.
  
  // Calculate D* parameters, cut, and store the result for both "right" 
  // and "wrong" sign candidates.
  // Allow more than one entry per vertex to avoid bias of wrong sign vs. 
  // right sign (but should be rare). 


  // 'D+' branch can be forked from D0 branch after initial selection:
  // Loop over all possible third track candidates and revertex + cut
  // (not yet implemented).

  // For all meson types, cuts may differ depending on whether 'trigger
  // vertex' is used or not. 


  // declare and initialize variables that u want to use for Dmeson only 
  float vc[3] = {0.};
  float vcD0[3] = {0.};
  float vc2[3] = {0.};
  int sumntrk = 0;
  float sumpt1 = 0.;
  float sumpt2 = 0.;
  float sumptdR15 = 0.;
  float sumptdR10 = 0.;
  float sumptdR07 = 0.;
  float sumptdR04 = 0.;
  
  float vc1[3] = {0.};
  float vc3[3] = {0.};

  int simidD0 = -1; 
  int t1pdgid = 0; 
  int t2pdgid = 0; 
  int vtxidD0 = -1; 
 
  float zD0 = 0.;
  float zD015 = 0.;
  float zD010 = 0.;
  float zD007 = 0.;
  float zD004 = 0.;
  float zDstar = 0.;

  float rapD0 = 999.;
  float chi2D0 = 999.;
  float dlxyD0 = 0.;
  float dlxyerrD0 = 0.;
  float dlxysigD0 = 0.;
  float cosphixyD0 = 0.;
  float dlD0 = 0.;
  float dlerrD0 = 0;
  float dlsigD0 = 0; 
  float cosphiD0 = 0;
  
  float deltam, deltam2;
  float deltamr, deltamr2;
  float chi2Pis = 999.;
  float rapDstar = 999.;
  int   simidDstar = -1;
  float mD012 = 0.;
  float mD021 = 0.;
  float mD0KK = 0.;

  // for both AdaptiveVertexFitter and KalmanVertexFitter options  
  vector<TransientTrack> mytracksforD0;
  vector<TransientTrack> mytracksforpislow;
  vector<TransientTrack> mytracksforDstar;
  TransientVertex myVertex;
  TransientVertex myPisVertex;

  bool VtxFound = false;
  int VtxCounter = 0;
  bool VtxFound2 = false;
  int VtxCounter2 = 0;
  int ivtx =0;
  float distzmin = 999.;
  int ivtxmin = -1;

  // define four-vectors
  TLorentzVector p4K1;
  TLorentzVector p4K2;
  TLorentzVector p4pi1;
  TLorentzVector p4pi2;
  TLorentzVector p4D012;
  TLorentzVector p4D021;
  TLorentzVector p4D0KK;
  TLorentzVector p4pis, p4pisr;
  TLorentzVector p4Dstar, p4Dstar2, p4Dstarr, p4Dstarr2;
  TLorentzVector p4DstarD0;
  
  p4K1.SetPtEtaPhiM(0., 0., 0., 0.);
  p4K2.SetPtEtaPhiM(0., 0., 0., 0.);
  p4pi1.SetPtEtaPhiM(0., 0., 0., 0.);
  p4pi2.SetPtEtaPhiM(0., 0., 0., 0.);
  p4D012.SetPtEtaPhiM(0., 0., 0., 0.);
  p4D021.SetPtEtaPhiM(0., 0., 0., 0.);
  p4pis.SetPtEtaPhiM(0., 0., 0., 0.);
  p4Dstar.SetPtEtaPhiM(0., 0., 0., 0.);
  p4DstarD0.SetPtEtaPhiM(0., 0., 0., 0.);
      
  reco::VertexCollection::const_iterator iteForD0;
  reco::VertexCollection::const_iterator itemin;
  //  vector<TransientTrack>::iterator d1;
  //  vector<TransientTrack>::iterator d2;


////////////////////////////////////////////////////////////////////////////
// List of cuts or selections applied to Dmeson candidates                //
////////////////////////////////////////////////////////////////////////////
  UInt_t trkminD0 = 2;     // minimum number of tracks for D0 
  UInt_t trkminDstar = 3;  // minimum number of tracks for Dstar
  //Float_t ptKpimin = 0.5;  // minimum pT for kaon and pion candidate
  Float_t ptKpimin = 0.3;  // minimum pT for kaon and pion candidate
                           // (before type is decided, before revertexing)
  //Float_t ptKmin = 0.4;    // minimum pT for Kaon candidate (not yet used)
  //Float_t ptpimin = 0.3;   // minimum pT for pion candidate (not yet used)
  // (the softer cuts might be relevant for low pt D0's)
  // (for these: use also dEdx information?)

  // WorkBook: For matching, vertices & track origins should be separated by 
  // at most 1 mm and at most 3 sigma; 
  // here we use 1 mm only throughout, from a previous empirical study,
  // and since we want to associate also close secondaries 
  Float_t vKpi_xymax = 0.1; // maximum vertex distance in xy for K and pi 
  //Float_t vKpi_xymax = 0.15; // maximum vertex distance in xy for K and pi 
  Float_t vKpi_zmax = 0.1;  // maximum vertex distance in z for K and pi
  Float_t vD0pis_xymax = 0.1; // maximum vertex distance in xy for pis
  Float_t vD0pis_zmax = 0.1;  // maximum vertex distance in z for pis 
  Float_t vD0trksum_xymax = 0.1; // maximum vertex distance in xy for D0
  Float_t vD0trksum_zmax = 0.1;  // maximum vertex distance in z for D0
  Float_t vdistzmax = 1.; // maximum z distance of prim. vtx to be be considered
  // Float_t vKpi_xymax = 0.2; // maximum vertex distance in xy for K and pi 
  // Float_t vKpi_zmax = 0.2;  // maximum vertex distance in z for K and pi
  // Float_t vD0pis_xymax = 0.2; // maximum vertex distance in xy for pis
  // Float_t vD0pis_zmax = 0.2;  // maximum vertex distance in z for pis 
  // Float_t vD0trksum_xymax = 0.2; // maximum vertex distance in xy for D0
  // Float_t vD0trksum_zmax = 0.2;  // maximum vertex distance in z for D0
  Float_t mD0min1 = 1.5;  // loose minimum D0 mass for D*
  Float_t mD0max1 = 2.3;  // loose maximum D0 mass for D* 
  Float_t D0ptmin = 3.5;  // minimum D0 pt to accept low z (also for D*)
  //Float_t zD0min  = 0.1;  // z cut on low pt D0's (also for D*)
  Float_t zD0min  = 0.0;  // z cut on low pt D0's (also for D*)
  Float_t D0ptmin1= 1.5;  // minimum D0 pt for medium z (also for D*)
  // Float_t zD0min1 = 0.2;  // z cut on very low pt D0's (also for D*)
  Float_t zD0min1 = 0.15;  // z cut on very low pt D0's (also for D*)
  Float_t DstarD0ptmin = 1.4; // minimum D0 pt for Dstar only
                              // reflects implicit slow pion cutoff
  // Float_t dlxyD0min = -99; // minimum decay length to store D0 (not for Dstar)
  Float_t dlxyD0min = 0.02; // minimum decay length to store D0 (not for Dstar)
  Float_t dlsigD0min = 1.5; // minimum significance to store D0 (not for Dstar)
  //  Float_t zDstarmin = 0.05; // final z cut on all Dstar candidates
  Float_t zDstarmin = 0.00; // final z cut on all Dstar candidates
                            // (beware for multi-parton interactions, 
                            //  associated production, 
                            //  and/or fragmentation function measurements) 
  // Float_t mD0min = 1.75;  // minimum D0 mass for D0
  // Float_t mD0max = 1.98;  // maximum D0 mass for D0
  Float_t mD0min = 1.68;  // minimum D0 mass for D0
  Float_t mD0max = 2.05;  // maximum D0 mass for D0
  Float_t ptpismin = 0.;  // minimum slow pi pt
  Float_t dmDstarmax = 0.17; // standard maximum D* mass
  //Float_t cosphimin = 0.99; // minimum cosphi for D0 w.r.t. primary vertex
  Float_t cosphimin = -1.; // minimum cosphi for D0 w.r.t. primary vertex

  // tight cuts for D* candidate 'cross' storage
  // Float_t mD0tmin = mD0Actual-0.025;  // tight minimum D0 mass
  // Float_t mD0tmax = mD0Actual+0.025;  // tight maximum D0 mass
  Float_t mD0tmin = mD0Actual-0.04;  // tight minimum D0 mass
  Float_t mD0tmax = mD0Actual+0.04;  // tight maximum D0 mass
  // Float_t dmDstartmin = dmDstarActual-0.001; // tight mimimum D* mass
  // Float_t dmDstartmax = dmDstarActual+0.001; // tight maximum D* mass
  Float_t dmDstartmin = dmDstarActual-0.002; // tight mimimum D* mass
  Float_t dmDstartmax = dmDstarActual+0.002; // tight maximum D* mass
///////////////////////////
/// end of list of cuts ///
///////////////////////////

    // clear the storage containers for this objects in this event
  //------------------------------- For D0 branches -------------------------//
  nD0 = 0;

  // Kaon from D0
  D0t1_pt.clear();
  D0t1_eta.clear();
  D0t1_phi.clear();
  D0t1_chg.clear();
  D0t1_tkIdx.clear();
  D0t1_Kprob.clear();
  D0t1_piprob.clear();
  D0t1_dEdxnmeas.clear();
  D0t1_dEdxnsat.clear();
  D0t1_vtxIdx.clear();
  D0t1_chindof.clear();
  D0t1_nValid.clear();
  D0t1_nPix.clear();
  D0t1_isHighPurity.clear();
  D0t1_dxy.clear();
  D0t1_dz.clear();
  D0t1_pdgId.clear();

  // Pion from D0
  D0t2_pt.clear();
  D0t2_eta.clear();
  D0t2_phi.clear();
  D0t2_chg.clear();
  D0t2_tkIdx.clear();
  D0t2_Kprob.clear();
  D0t2_piprob.clear();
  D0t2_dEdxnmeas.clear();
  D0t2_dEdxnsat.clear();
  D0t2_vtxIdx.clear();
  D0t2_chindof.clear();
  D0t2_nValid.clear();
  D0t2_nPix.clear();
  D0t2_isHighPurity.clear();
  D0t2_dxy.clear();
  D0t2_dz.clear();
  D0t2_pdgId.clear();

  // D0
  D0_pt.clear();
  D0_eta.clear();
  D0_phi.clear();
  D0_rap.clear();
  D0_mass12.clear();
  D0_mass21.clear();
  D0_massKK.clear();
  D0_simIdx.clear();
  D0_DstarIdx.clear();
  D0_ambiPrim.clear();
  D0_vtxIdx.clear();
  D0_hasMuon.clear();
  D0_chi2.clear();
  D0_dlxy.clear();
  D0_dlxyErr.clear();
  D0_dlxySig.clear();
  D0_cosphixy.clear();
  D0_dl.clear();
  D0_dlErr.clear();
  D0_dlSig.clear();
  D0_cosphi.clear();
  D0_ptfrac.clear();
  D0_ptfrac15.clear();
  D0_ptfrac10.clear();
  D0_ptfrac07.clear();
  D0_ptfrac04.clear();
  D0_x.clear();
  D0_y.clear();
  D0_z.clear();
  D0_Covxx.clear();
  D0_Covyx.clear();
  D0_Covzx.clear();
  D0_Covyy.clear();
  D0_Covzy.clear();
  D0_Covzz.clear();

  //----------------------------- For Dstar branches ------------------------//
  nDstar = 0;

  // Slow Pion from Dstar
  Dstarpis_pt.clear();
  Dstarpis_eta.clear();
  Dstarpis_phi.clear();
  Dstarpis_ptr.clear();
  Dstarpis_etar.clear();
  Dstarpis_phir.clear();
  Dstarpis_chg.clear();
  Dstarpis_tkIdx.clear();
  Dstarpis_Kprob.clear();
  Dstarpis_piprob.clear();  
  Dstarpis_dEdxnmeas.clear();
  Dstarpis_dEdxnsat.clear();
  Dstarpis_vtxIdx.clear();
  Dstarpis_chindof.clear();
  Dstarpis_chir.clear();
  Dstarpis_nValid.clear();
  Dstarpis_nPix.clear();
  Dstarpis_dxy.clear();
  Dstarpis_dz.clear();

  // D0 from Dstar
  DstarD0_pt.clear();
  DstarD0_eta.clear();
  DstarD0_phi.clear();
  DstarD0_mass.clear();
  DstarD0_chi2.clear();
  DstarD0_dlxy.clear();
  DstarD0_dlxyErr.clear();
  DstarD0_dlxySig.clear();
  DstarD0_cosphixy.clear();
  DstarD0_dl.clear();
  DstarD0_dlErr.clear();
  DstarD0_dlSig.clear();
  DstarD0_cosphi.clear();
  DstarD0_ptfrac.clear();
  DstarD0_ptfrac15.clear();
  DstarD0_ptfrac10.clear();
  DstarD0_ptfrac07.clear();
  DstarD0_ptfrac04.clear();
  DstarD0_x.clear();
  DstarD0_y.clear();
  DstarD0_z.clear();
  DstarD0_simIdx.clear();
  DstarD0_recIdx.clear();
  DstarD0_ambiPrim.clear();

  // Kaon from Dstar  
  DstarK_pt.clear();
  DstarK_eta.clear();
  DstarK_phi.clear();
  DstarK_chg.clear();
  DstarK_tkIdx.clear();
  DstarK_Kprob.clear();
  DstarK_piprob.clear();
  DstarK_dEdxnmeas.clear();
  DstarK_dEdxnsat.clear();
  DstarK_vtxIdx.clear();
  DstarK_chindof.clear();
  DstarK_nValid.clear();
  DstarK_nPix.clear();
  DstarK_isHighPurity.clear();
  DstarK_dxy.clear();
  DstarK_dz.clear();

  // Pion from Dstar
  Dstarpi_pt.clear();
  Dstarpi_eta.clear();
  Dstarpi_phi.clear();
  Dstarpi_chg.clear();
  Dstarpi_tkIdx.clear();
  Dstarpi_Kprob.clear();
  Dstarpi_piprob.clear();
  Dstarpi_dEdxnmeas.clear();
  Dstarpi_dEdxnsat.clear();
  Dstarpi_vtxIdx.clear();
  Dstarpi_chindof.clear();
  Dstarpi_nValid.clear();
  Dstarpi_nPix.clear();
  Dstarpi_isHighPurity.clear();
  Dstarpi_dxy.clear();
  Dstarpi_dz.clear();
  
  // Dstar
  Dstar_pt.clear();
  Dstar_eta.clear();
  Dstar_phi.clear();
  Dstar_rap.clear();
  Dstar_deltam.clear();
  Dstar_deltamr.clear();
  Dstar_simIdx.clear();
  Dstar_vtxIdx.clear();
  Dstar_hasMuon.clear();
  Dstar_ptfrac.clear();

  // check for track collection with at least 2 tracks (for D0) 
  if (tracks->size() >= trkminD0) {
    
#ifndef miniAOD
    // build transient track list from all tracks for later revertexing
    vector<reco::TransientTrack> genralTracks_forD = (*theB).build(tracks);
#endif

    //---------- 1st loop: over tracks (K or pi)  -----------//
    int it1count = -1;
#ifndef miniAOD
    // loop over AOD track collection
    for (reco::TrackCollection::const_iterator it1 = tracks->begin(); it1 != tracks->end(); ++it1) {
      ++it1count;
#endif
#ifdef miniAOD
    // loop over miniAOD PackedCandidates and find tracks with track details
    // (pt > 0.4 GeV)
    for (pat::PackedCandidateCollection::const_iterator ic1 = tracks->begin(); ic1 != tracks->end(); ++ic1) {
      ++it1count;
      if (!(ic1->hasTrackDetails())) continue; 
      auto tk1 = ic1->pseudoTrack();  // needed for transient track build
      auto it1 = ic1->bestTrack(); 
      if (it1 == nullptr) continue; 
#endif
      // it1 is now track iterator for AOD tracks, and track pointer for 
      // 'full info' miniAOD tracks

      // print info about all tracks within 1 cm of true D0 //
      if (Bingo==1 && abs(it1->vx()-D0vtxx)< 1. && abs(it1->vy()-D0vtxy)<1. && abs(it1->vz()-D0vtxz)<1.) {
      cout << "track pt " << it1->pt() << " eta " << it1->eta() << " phi " << it1->phi() << " vtx " << it1->vx() << " " << it1->vy() << " " << it1->vz() << endl; 
      }

      // fill inclusive pt and eta histograms
      h_trackpt->Fill(it1->pt());
      h_trackptlow->Fill(it1->pt());
      h_tracketa->Fill(it1->eta());

      // pt cut on first track
      if (it1->pt() > ptKpimin) {

	//cout << "KELUAR X 0?" << endl;

        // store track quality variables //
        // need to deal with fact that first track might be used with 
        // several 2nd tracks ...
        int t1vtxid = -1;
#ifndef miniAOD
        float t1chindof = it1->chi2()/it1->ndof();
#endif
#ifdef miniAOD
        float t1chindof = it1->normalizedChi2();
#endif
        int t1nvalid = it1->hitPattern().numberOfValidHits();
        int t1npix = it1->hitPattern().numberOfValidPixelHits();
        //bool t1highp = it1->trackHighPurity();
        bool t1highp = it1->quality(Track::highPurity);
        float t1dxy = 0;
        float t1dz = 0;

        float t1piprob = -9999.;
        float t1Kprob = -9999.;
        int dedxtrknmeas = -9999;
        int dedxtrknsat = -9999;
#ifndef CMSSW7plus
        // for dEdx (from example by G. Fedi)
        // beware of potential slowdown!
        //cout << "Hello dEdx" << endl;
	reco::TrackRef trackRef = reco::TrackRef(tracks,it1count);   
        //cout << "Hello dEdx " << dedxexist << endl;
        if (dedxexist>0) {
          // the following is OK
          //cout << "Hello dEdx " << trackRef->pt() << endl;
          //cout << "Hello dEdx " << energyLossHandle->size() << " " << tracks->size() << endl;
          //cout << "Hello dEdx " << (*energyLossHandle)[trackRef].dEdx() << endl;
          // crash on the following line (eloss):
          //float dedxtrk = eloss[trackRef7].dEdx();
          //float dedxtrk = (*energyLossHandle)[trackRef].dEdx();
          //float dedxtrkerr = (*energyLossHandle)[trackRef].dEdxError();
          double dedxtrk = (*energyLossHandle)[trackRef].dEdx();
          double dedxtrkerr = (*energyLossHandle)[trackRef].dEdxError();
          dedxtrknsat = (*energyLossHandle)[trackRef].numberOfSaturatedMeasurements();
          dedxtrknmeas = (*energyLossHandle)[trackRef].numberOfMeasurements();
          // temporarily "misuse" the piprob and Kprob variables for direct dEdx:
          t1piprob = dedxtrkerr;    
          t1Kprob = dedxtrk;
        }
#endif    
#ifdef CMSSW7plus
        // for dEdx (from example by M. Soares)
        // beware of potential slowdown!
	reco::TrackRef trackRef = reco::TrackRef(tracks,it1count);   
        double stripMapdedx = ( stripmap >0 ? (*dedxStMap)[trackRef].dEdx() : 0.);
        double pixMapdedx = ( pixmap >0 ? (*dedxPixMap)[trackRef].dEdx() : 0.);
        // temporarily "misuse" the piprob and Kprob variables for direct dEdx:
        t1piprob = pixMapdedx;    
        t1Kprob = stripMapdedx;    
        dedxtrknsat = ( stripmap >0 ? (*dedxStMap)[trackRef].numberOfSaturatedMeasurements() : 0.);
	dedxtrknmeas = ( stripmap >0 ? (*dedxStMap)[trackRef].numberOfMeasurements() : 0.);
#endif

        // this would be the place to cut on these variables, if wished! //
        // skip unnecessary low pt tracks
        if (it1->pt() < 0.5 && t1Kprob < 4) continue; 
        //cout << "Hello dEdx survived" << endl;

	//------ 2nd loop: over tracks (pi or K) ------------//
        // loop over each pair only once;
        // which is K and which is pi is to be decided later
        int it2count = it1count -1;
#ifndef miniAOD
        // loop over AOD track collection
	for (reco::TrackCollection::const_iterator it2 = it1; it2 != tracks->end(); ++it2) {
          ++it2count;
#endif
#ifdef miniAOD
        // loop over miniAOD PackedCandidate collection
	for (pat::PackedCandidateCollection::const_iterator ic2 = ic1; ic2 != tracks->end(); ++ic2) {
          ++it2count;
          if (!(ic2->hasTrackDetails())) continue; 
          auto tk2 = ic2->pseudoTrack();  // needed for transient track build
          auto it2 = ic2->bestTrack(); 
          if (it2 == nullptr) continue; 
#endif	  
          // it2 is now iteractor for AOD tracks or pointer for miniAOD 
          // tracks with full information

          // both unlike sign (for D0 and D*) and like sign (for wrong sign D*)
          // will be accepted

	  // check track 2 is not track 1 
          // (do this more elegantly directly in iterator?)
	  if (it2 == it1) continue;
          // apply pt cut (no explicit eta cut!)
          if (it2->pt() > ptKpimin) {

            // store track quality variables //
            // logic of resetting t1 should be reconsidered
            t1vtxid = -1;
            int t2vtxid = -1;
            bool ambiprimary = false;
#ifndef miniAOD
            float t2chindof = it2->chi2()/it2->ndof();
#endif
#ifdef miniAOD
            float t2chindof = it2->normalizedChi2();
#endif
            int t2nvalid = it2->hitPattern().numberOfValidHits();
            int t2npix = it2->hitPattern().numberOfValidPixelHits();
            //bool t2highp = it2->trackHighPurity();
            bool t2highp = it2->quality(Track::highPurity);
            float t2dxy = 0;
            float t2dz = 0;

            float t2piprob = -9999.;
            float t2Kprob = -9999.;
            int dedxtrk2nmeas = -9999;
            int dedxtrk2nsat = -9999;
#ifndef CMSSW7plus
            // for dEdx (from example by G. Fedi)
            // beware of potential slowdown!
	    reco::TrackRef trackRef = reco::TrackRef(tracks,it2count);   
            double dedxtrk = (*energyLossHandle)[trackRef].dEdx();
            double dedxtrkerr = (*energyLossHandle)[trackRef].dEdxError();
            dedxtrk2nsat = (*energyLossHandle)[trackRef].numberOfSaturatedMeasurements();
            dedxtrk2nmeas = (*energyLossHandle)[trackRef].numberOfMeasurements();
            // temporarily "misuse" the piprob and Kprob variables for direct dEdx:
            t2piprob = dedxtrkerr;    
            t2Kprob = dedxtrk;
#endif    
#ifdef CMSSW7plus
            // for dEdx (from example by M. Soares)
            // beware of potential slowdown!
	    reco::TrackRef trackRef = reco::TrackRef(tracks,it2count);   
            double stripMapdedx = ( stripmap >0 ? (*dedxStMap)[trackRef].dEdx() : 0.);
            double pixMapdedx = ( pixmap >0 ? (*dedxPixMap)[trackRef].dEdx() : 0.);
            // temporarily "misuse" the piprob and Kprob variables for direct dEdx:
            t2piprob = pixMapdedx;    
            t2Kprob = stripMapdedx;
            dedxtrk2nsat = ( stripmap >0 ? (*dedxStMap)[trackRef].numberOfSaturatedMeasurements() : 0.);
            dedxtrk2nmeas = ( stripmap >0 ? (*dedxStMap)[trackRef].numberOfMeasurements() : 0.);
#endif

           // this would be the place to cut on these variables, if wished! //
           // skip unnecessary low pt tracks
	   if (it2->pt() < 0.5 && t2Kprob < 4) continue; 
           //cout << "Hello dEdx 2 survived" << endl;

	    // calculate track origin separations in x, y and z
	    vc[0] = abs(it1->vx() - it2->vx());
	    vc[1] = abs(it1->vy() - it2->vy());
	    vc[2] = abs(it1->vz() - it2->vz());

	    // check they are separated at most by a radius of 1mm
	    // in x and y and a distance of 0.1 in z
	    if (sqrt ((vc[0] * vc[0]) + (vc[1] * vc[1])) < vKpi_xymax && vc[2] < vKpi_zmax) {
		  
	      // Calculate 4 momentum for kaon and pion using TLorentzVector
              // treat the two possible combinations (Kpi or piK)
	      p4K1.SetPtEtaPhiM(it1->pt(), it1->eta(), it1->phi(), Kmass);
	      p4K2.SetPtEtaPhiM(it2->pt(), it2->eta(), it2->phi(), Kmass);
	      p4pi1.SetPtEtaPhiM(it1->pt(), it1->eta(), it1->phi(), pimass);
	      p4pi2.SetPtEtaPhiM(it2->pt(), it2->eta(), it2->phi(), pimass);
		              
	      // get four-vector, and from this pt, eta, phi, and mass 
              // (2 possible combinations)
	      p4D012 = p4K1 + p4pi2;
	      p4D021 = p4K2 + p4pi1;
	      p4D0KK = p4K1 + p4K2;
              // the two masses can kinematically differ by up to a factor 
              // mK/mpi for very asymmetric momenta.
              // for completely symmetric momenta, they will be identical.
              // investigate cut on momentum asymmetry (costheta*)?
	      mD012 = p4D012.M();
	      mD021 = p4D021.M();
	      mD0KK = p4D0KK.M();

	      if (Bingo) cout << "2nd track pt " << it2->pt() << " dist "<< vc[0] << " " << vc[1] << " " << vc[2] << " mass " << mD012 << " " << mD021 << endl;

              // apply loose mass cut
              if (   (mD012 < mD0min1 || mD012 > mD0max1)
		  && (mD021 < mD0min1 || mD021 > mD0max1) ) continue;

              //cout << "mD0 " << mD012 << " " << mD021 << " " << mD0KK << endl;
	      
              // determine D0 px, py, pz
	      // paxisD0[0] = it1->px() + it2->px();
	      // paxisD0[1] = it1->py() + it2->py();
	      // paxisD0[2] = it1->pz() + it2->pz();
			      
	      // roughly locate D0 'vertex' using average coordinate of tracks
	      vcD0[0] = 0.5 * (it1->vx() + it2->vx());
	      vcD0[1] = 0.5 * (it1->vy() + it2->vy());
	      vcD0[2] = 0.5 * (it1->vz() + it2->vz());

	      // move this outside loop after having cleaned up vertex structure
              // *** even better, directly use vertex variable ***
              // ... but this does not allow partial sums, 
              //     so keep here after all?
              sumntrk = 0;
	      sumpt1 = 0;
	      sumptdR15 = 0;
	      sumptdR10 = 0;
	      sumptdR07 = 0;
	      sumptdR04 = 0;
	      // loop over all tracks
#ifndef miniAOD
              // loop over AOD track collection
	      for (reco::TrackCollection::const_iterator itSum1 = tracks->begin(); itSum1 != tracks->end(); ++itSum1) {
#endif
#ifdef miniAOD
              // loop over miniAOD PackedCandidates with full track info
              // (some low pt tracks will be missed)
	      for (pat::PackedCandidateCollection::const_iterator icSum1 = tracks->begin(); icSum1 != tracks->end(); ++icSum1) {
		if (!(icSum1->hasTrackDetails())) continue;
                auto itSum1 = icSum1->bestTrack();
                if (itSum1 == nullptr) continue; 
#endif
		// check track origins
		vc1[0] = abs(vcD0[0] - itSum1->vx());
		vc1[1] = abs(vcD0[1] - itSum1->vy());
		vc1[2] = abs(vcD0[2] - itSum1->vz());
		
	        if ((sqrt(vc1[0] * vc1[0]) + (vc1[1] * vc1[1])) < vD0trksum_xymax && vc1[2] < vD0trksum_zmax) {
                  ++sumntrk;  // number of tracks associated to this vertex
		  sumpt1 += abs(itSum1->pt()); // sum pt for all tracks
                  float deltaeta1 = p4D012.Eta()-itSum1->eta();
                  float deltaphi1 = p4D012.Phi()-itSum1->phi();
                  if (deltaphi1>3.1415) deltaphi1 = deltaphi1-2.*3.1415;
		  if (deltaphi1<-3.1415) deltaphi1 = deltaphi1+2.*3.1415;
		  float deltaR1 = sqrt(deltaeta1*deltaeta1+deltaphi1*deltaphi1);
                  if (deltaR1<1.5) {
                    sumptdR15 += abs(itSum1->pt());
                    if (deltaR1<1.0) { 
                      sumptdR10 += abs(itSum1->pt());
                      if (deltaR1<0.7) {
                        sumptdR07 += abs(itSum1->pt());
                        if (deltaR1<0.4) { 
                          sumptdR04 += abs(itSum1->pt());
                        }
                      }
                    }
                  }
		} // end of Sumpt vertex check
	      } // end of itSum loop
	      //cout << "KELUAR X 1?" << endl;

              // will be slightly different for AOD and miniAOD
	      // old: zD0 = (p4D0.Pt()) / sumpt1;
	      zD0   = (it1->pt() + it2->pt()) / sumpt1;
	      zD015 = (it1->pt() + it2->pt()) / sumptdR15;
	      zD010 = (it1->pt() + it2->pt()) / sumptdR10;
	      zD007 = (it1->pt() + it2->pt()) / sumptdR07;
	      zD004 = (it1->pt() + it2->pt()) / sumptdR04;

              // reject if both low pt and low z
              // (to limit revertexing, this will also propagate to D*)
              // this is the same for the 12 and 21 combinations
              if (p4D012.Pt() < D0ptmin  && zD0 < zD0min)  continue;
              if (p4D012.Pt() < D0ptmin1 && zD0 < zD0min1) continue;

	      // build D0 track collection for revertexing
		mytracksforD0.clear();
#ifndef miniAOD
                // loop over prebuilt list of AOD transient tracks
		for (vector<TransientTrack>::iterator gt_trans = genralTracks_forD.begin(); gt_trans != genralTracks_forD.end(); ++gt_trans) {
					
		  const reco::TrackRef trackRef1 = (gt_trans->trackBaseRef()).castTo<reco::TrackRef>();
		  if (&*it1 == trackRef1.get() || &*it2 == trackRef1.get()) {
					    
		    TransientTrack  transientTrack1 = theB->build(trackRef1);
		    mytracksforD0.push_back(transientTrack1);
		  }
		}
#endif
#ifdef miniAOD 
                // build transient tracks from relevant miniAOD tracks
	        TransientTrack  transientTrack1 = theB->build( tk1 );
		mytracksforD0.push_back(transientTrack1);
	        TransientTrack  transientTrack2 = theB->build( tk2 );
		mytracksforD0.push_back(transientTrack2);
#endif              

		//cout << "KELUAR X 2?" << endl;
		  
                double D0x = 0, D0y=0, D0z=0;
                double cov_sv[3][3], cov_pv[3][3], deriv[3];
                double dlerr = 0, dlxyerr=0;

                // important to initialize already here
                ivtxmin=-1;

		// turn if into continue? //
		if (mytracksforD0.size() > 1/*&& !VtxFound2*/) {
  					
                  // initialize vertex finding flags
                  // VtxFound indicates match to particular vertex
                  // VtxFound2 indicates that any match has been found
                  // Does this logic still make sense after structural change?
		  VtxFound = false;
  		  VtxFound2 = false;

//////////////////////////////////////////////////////////////////////////////
///////////////// Do the D0 vertex refit /////////////////////////////////////
// 		  AdaptiveVertexFitter  theFitter1;
//		  // if you don't/do want the beam constraint
//                // (since we look for secondary vertices we don't) 
//		  myVertex = theFitter1.vertex(mytracksforD0/*, vertexBeamSpot*/);  
//                // can the latter option also be used to add a track to an
//                // already existing vertex (e.g. a muon to a DO)?
//
//                // track weights can be accessed through 
//                // TransientVertex::trackWeight(trk)

//////////////////////////////////////////////////////////////////////////////
///////////////// Do the D0 vertex refit /////////////////////////////////////
                  // allow track parameter reevaluation 
                  // (Twiki: currently only the KalmanVertex can do this)
 		  KalmanVertexFitter  theFitter2(true);
                  // do the fit
		  myVertex = theFitter2.vertex(mytracksforD0);  
                  if (!myVertex.isValid()) continue;

		  //cout << "haha refitted" << endl;

///////////////// retrieve information ////////////////////////////////////

                  // vertex parameters
		  float ndof = myVertex.degreesOfFreedom();
                  // for Adaptive Vertex Fitter this sum of weights*2-3
                  // for Kalman Vertex Fitter this is 1 for two tracks?
                  chi2D0 = myVertex.totalChiSquared();
                  // confidence level
                  float CL=999.;
                  if (ndof >= 1) {
                    CL = TMath::Prob(chi2D0,(int)ndof); 
                  }
                  // cut on chisquared (Adaptive) or CL (Kalman) here?
                  if (CL<0.01) continue;

                  // good vertex, get position
                  D0x = myVertex.position().x();
                  D0y = myVertex.position().y();
                  D0z = myVertex.position().z();

                  // track parameters
                  if (myVertex.hasRefittedTracks()){
                    // refitted tracks are TransientTracks
                    // get updated momenta of refitted tracks
                    // and recalculate K, pi and D0 quantities 
                    int itcount = 0;
                    vector<TransientTrack> tracks = myVertex.refittedTracks();
                    for (vector<TransientTrack>::const_iterator trackIt = tracks.begin(); trackIt != tracks.end(); ++trackIt) {
                      ++itcount;
                      const Track & track = trackIt->track();
                      if (itcount==1) {
                        p4K1.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), Kmass);
                        p4pi1.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), pimass);
                      }
		      else if (itcount==2) { 
                        p4K2.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), Kmass);
                        p4pi2.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), pimass);
                      }		              
	              // get four-vector, and from this pt, eta, phi, and mass 
                      // (2 possible combinations)
	              p4D012 = p4K1 + p4pi2;
	              p4D021 = p4K2 + p4pi1;
	              p4D0KK = p4K1 + p4K2;
                      // the two masses can kinematically differ by up to a factor 
                      // mK/mpi for very asymmetric momenta.
                      // for completely symmetric momenta, they will be identical.
                      // investigate cut on momentum asymmetry (costheta*)?
	              mD012 = p4D012.M();
	              mD021 = p4D021.M();
	              mD0KK = p4D0KK.M();

                      //cout << "Track pt: " << track.pt() << endl;
                      //cout << it1->pt() << " " << it2->pt() << endl;
                    }
                  }  

                  // get point of closest approach to beam spot 
                  // *** to be calculated ***
                  // for the moment, directly use D0 vertex position

/////////////////////////// find best primary vertex /////////////////////////

                  // loop over primary vertices
                  // assume each track is only associated to one vertex
                  //                             true? 
                  // TrackRefs are CPU-expensive -> minimize

                  // VtxCounter2 counts number of associated vertices found
                  // 0 means no vertex association found
                  // 1 means association found to one vertex (ideal)
                  // 2 means the tracks are associated to different vertices 
  		  VtxCounter2 = 0;
                  // ivtx sets vertex index = unique id
                  ivtx=-1;
                  distzmin=999.;
                  ivtxmin=-1;
                  vtxidD0=-1;

		  for (reco::VertexCollection::const_iterator ite = Primvertex->begin(); ite != Primvertex->end(); ++ite) {
                    ++ivtx;
			  // if (VtxFound2 && (d1 == mytracksforD0.begin() || d2 == mytracksforD0.begin()) ) {continue;}

                    // get vertex distance in z *** update ***
                    float distz = fabs(ite->z()-D0z);
                    // skip if too far away
                    if (distz > vdistzmax) continue;
                    // should we apply quality cuts?
                    if (ite->isFake() ||
                        ite->ndof()<2 || fabs(ite->z()-Bsp_z)>20.) continue;
                    //if (ite->isFake() || !ite->isValid() ||
                    //    ite->ndof()<2 || fabs(ite->z()-Bsp_z)>20.) continue;
                    // isFake means beam spot, isValid means fit converged
                    // and result is within tracker boundaries (always true?), 
                    // the last two are essentially no cut;
                    // (ndof>4 for PVtx_isGood means >2 tracks)
                    // use PVtx variables to cut harder offline if wished

                    // save for later if closer than closest good vertex so far
                    if (distz < distzmin) {
                      distzmin = distz;
                      ivtxmin = ivtx;
                      itemin = ite;
		      //cout << "ivtxmin " << ivtxmin << " distz " << distz << endl;
                    }

                    // loop over tracks from this vertex
                    // (*** this will not work for miniAOD, 
                    //  loop over Packed Candidates and use vertexRef! ***)
                    // VtxCounter counts how many matches to D0 tracks have 
                    // been found for this vertex. Should be 0, 1, or 2  
		    VtxCounter = 0;
#ifndef miniAOD
		    for (reco::Vertex::trackRef_iterator iTrack = ite->tracks_begin(); iTrack != ite->tracks_end(); ++iTrack) {
		      // get track reference	    
		      const reco::TrackRef trackRef = iTrack->castTo<reco::TrackRef>();
		      // loop over D0 tracks	    
                      // VtxFound indicates that a D0 track is at this vertex
                      VtxFound = false;
                      int itcount = 0;
		      for (vector<TransientTrack>::iterator gt1 = mytracksforD0.begin(); gt1 != mytracksforD0.end(); ++gt1) {
                        ++itcount; 
		        //cout<<"in  "<<(gt1->track()).pt()<<endl;
                        // do not check if vertex already found
                        if (itcount == 1 && t1vtxid >-1) continue;
                        if (itcount == 2 && t2vtxid >-1) continue;
                        // get track reference
			const reco::TrackRef trackRef2 = (gt1->trackBaseRef()).castTo<reco::TrackRef>();
			if (trackRef.get() == trackRef2.get()) {
                          VtxFound = true;
                          if (itcount == 1 && t1vtxid >-1) {
                            // should never happen if track at one vertex only
                            cout << "*** Alarm: track at more than one vertex?? ***" << endl;
                            exit(1);
                          }
                          if (itcount == 1) t1vtxid = ivtx;
                          else if (itcount == 2) t2vtxid = ivtx;
                          else {
                            cout << "*** nanoanalyzer: number of tracks inconsistent ***" << endl;
                          } // else  
                        } // trackref
		      } // gt1
		      if (VtxFound) {++VtxCounter;}
                      // stop track loop if both tracks have been found at 
                      // this vertex
                      if (VtxCounter == 2) break;
                      // or at two different vertices
                      if (t1vtxid >-1 && t2vtxid >-1) break;
		    } // itrack
#endif
#ifdef miniAOD
                    // need to find out how to compare vertexref to iterator ...
                    // *** the follwing needs cleanup ***
                    //const reco::VertexRef iteRef = ite->castTo<reco::VertexRef>(); 
                    //const reco::VertexRef iteRef = ite->Vertex();
                    //const reco::VertexRef iteRef = (ite -> key);
                    //const reco::VertexRef iteRef = ite->VertexRef(); 
                    //if (ic1->vertexRef() == ite) {
                    //if (ic1->vertexRef() == iteRef) {
                    //if (ic1->vertexRef() == ite->first()) {
                    //if (ic1->vertexRef() == (ite->key)) {
                    // use z position as poor mans's solution; will equality always work? (rounding)
                    if (ic1->vertexRef()->z() == ite->z()) {
                      // make compiler happy ...
                      if (VtxFound) cout << "blabla" << endl;
                      ++VtxCounter;
		      t1vtxid = ivtx; 
                    }               
		    if (ic2->vertexRef()->z() == ite->z()) {
                      ++VtxCounter;
		      t2vtxid = ivtx; 
                    }               
                    // debug printout
                    //cout << "vertex z " << ite->z() << " " << ic1->vertexRef()->z() << " " << ic2->vertexRef()->z() << endl;   
                    //cout << "vertex for first/2nd D0 track " << t1vtxid << " " << t2vtxid << endl;
#endif 
	            if (VtxCounter>0) VtxCounter2 = VtxCounter2+VtxCounter;
                    // VtxFound2 indicates that match has been found
                    // save vertex pointer 
                    if (!VtxFound2 && VtxCounter > 0) {
                      // this is the first vertex found, store it
                      iteForD0 = ite; 		    
                      vtxidD0 = ivtx;
                      VtxFound2 = true;
                    } // VtxCounter
                    else if (VtxFound2 && VtxCounter > 0) {
                      // this is the second vertex found
                      // cout << "*** nanoanalyzer: vertex association inconsistent ***" << endl;
                      ambiprimary = true; 
                      // find out which is closer
                      if (vtxidD0 == ivtxmin) {
			// previous is minimal, keep previous, change nothing
                        iteForD0 = iteForD0;
                      }
		      else if (ivtx == ivtxmin) {
                        // current is closest, supersede
                        iteForD0 = ite; 		    
                        vtxidD0 = ivtx;
                      }
		      else {
                        // if there is a minimum vertex
                        // (should always be, since two vertices found)
                        if (ivtxmin > -1) {
                          // use current closest if in between iteForD0 and ite,
                          // (use >= and <= to avoid problems with rare coincides)
                          if ((itemin->z() >= iteForD0->z() && 
                               itemin->z() <= ite->z()) ||  
                              (itemin->z() <= iteForD0->z() && 
                               itemin->z() >= ite->z())) {
                            iteForD0 = itemin;
                            vtxidD0 = ivtxmin;
			  }
                          // otherwise use the one closer to closest
                          else if ((itemin->z() >= iteForD0->z() &&
                                    iteForD0->z() >= ite->z()) ||
                                   (itemin->z() <= iteForD0->z() && 
                                    iteForD0->z() <= ite->z())) {
                            // iteForD0 is closest to itemin
                            iteForD0 = iteForD0;
                          }
                          else if ((itemin->z() >= ite->z() &&
                                    ite->z() >= iteForD0->z()) ||
                                   (itemin->z() <= ite->z() && 
                                    ite->z() <= iteForD0->z())) {
                            // ite is closest to itemin
                            iteForD0 = ite;
                            vtxidD0 = ivtx;
                          }
			  else {
                            // should not happen
                            cout << itemin->z() << " " << ite->z() << " " << iteForD0-> z() << endl;
                            cout << "*** nanoAnalyzer: Alarm in vertexing logic 1 ***" << endl;
                            exit(1);
                          }
                        } // ivtxmin
			else {
                          // should not happen
                          cout << "*** nanoAnalyzer: Alarm in vertexing logic 2 ***" << endl;
                          exit(1);
                        } // ivtxmin  
		      } // vtxidD0 
		    } // vtxCounter
                    // stop vertex loop if both tracks have been associated
                    if (VtxCounter2 == 2) break;
		  } // vertex loop
	        } // mytracksforD0size

                // loop over primary vertices and store sumntrk
                for (uint pvtx = 0; pvtx<nPVtx; ++pvtx) {
                  if (PVtx_Id[pvtx] == vtxidD0) PVtx_ntrk[pvtx] = sumntrk;
                }  

                // algorithm to calculate ambiprimary might be messed up before
                // recalculate here:
                if (t1vtxid==t2vtxid && t1vtxid > -1) ambiprimary = false;
                else ambiprimary = true;

                if (!VtxFound2) {
                  // this should not any more always be the case for miniAOD
                  //cout << "*** nanoanalyzer: no associated vertex found ***" << endl;
		  // no vertex was associated to either track,
                  // use closest if found
                  ambiprimary = true;
                  if (ivtxmin >-1) {
                    iteForD0 = itemin;
                    vtxidD0 = ivtxmin;
                    VtxFound2 = true;
                    //cout << "ivtxmin " << ivtxmin << "D0 " << vtxidD0 << endl;
		  } // ivtxmin
                  // otherwise leave VtxFound2 false, i.e. skip D0 candidate
                } // VtxFound2

		//cout << "KELUAR X 5? " << VtxFound2 << endl;

		//		 // this part will reduce number of the loops
		//		 // in "find out closest PV to this D0 tracks" part:
		//		 if (VtxFound2) {		
		//		   d1 = mytracksforD0.begin();
		//		   d2 = ++d1;
		//		 }
		//		 // cout << "ATAU CRASH SINI 5?" << endl;

/////////////////////////////////////////////////////////////////////////////
                 // still need to treat case where no primary vertex was found
                 // ignore? (as done now) or use beam spot? 
                 if (VtxFound2) {
                    // found 'best' primary for this vertex
		    //cout << "hello vtxfound2 " << endl;
                    // get position
                    double xprim = iteForD0->position().x(); 
                    double yprim = iteForD0->position().y(); 
                    double zprim = iteForD0->position().z();
                    //cout << "position" << endl;
                    // calculate distance of tracks to primary
                    // *** rather call it1->dxy(iteForD0)? (see mu_dxy) ***
                    t1dxy = sqrt((it1->vx()-xprim)*(it1->vx()-xprim)
		      + (it1->vy()-yprim)*(it1->vy()-yprim));
                    t1dz = it1->vz()-zprim; 
                    t2dxy = sqrt((it2->vx()-xprim)*(it2->vx()-xprim)
		      + (it2->vy()-yprim)*(it2->vy()-yprim));
                    t2dz = it2->vz()-zprim; 
                    // calculate distance of D0 to primary
                    double vdx = D0x - xprim;
                    double vdy = D0y - yprim;
                    double vdz = D0z - zprim; 
		    double dist = sqrt(pow(vdx,2) + pow(vdy,2) + pow(vdz,2));
		    double distxy = sqrt(pow(vdx,2) + pow(vdy,2));
		    // get direction vector 
		    double Dpx = p4D012.Px(); 
		    double Dpy = p4D012.Py();
                    double Dpz = p4D012.Pz();
		    double Dp  = p4D012.P();
		    double Dpt  = p4D012.Pt();
		    // calculate 3D and 2D decay length
                    dlD0 = (Dpx*vdx + Dpy*vdy + Dpz*vdz)/Dp;
                    dlxyD0 = (Dpx*vdx + Dpy*vdy)/Dpt;
                    //cout << "dlD0" << endl;
		    // calculate cosphi
		    cosphiD0 = dlD0/dist;
		    cosphixyD0 = dlxyD0/distxy;
		    // cut on cosphi (0.99)
                    dlerrD0=0;
                    dlsigD0=0;
                    dlxyerrD0=0;
                    dlxysigD0=0;
                    //cout << "before position error" << endl;
		    if (cosphiD0 > cosphimin) {
                      // calculate significance in space
                      // transient D0 vertex
		      cov_sv[0][0] = myVertex.positionError().cxx();
		      cov_sv[1][0] = myVertex.positionError().cyx();
		      cov_sv[2][0] = myVertex.positionError().czx();
                      cov_sv[0][1] = cov_sv[1][0];
		      cov_sv[1][1] = myVertex.positionError().cyy();
		      cov_sv[2][1] = myVertex.positionError().czy();
                      cov_sv[0][2] = cov_sv[2][0];
                      cov_sv[1][2] = cov_sv[2][1];
		      cov_sv[2][2] = myVertex.positionError().czz();
                      // primary vertex
		      cov_pv[0][0] = iteForD0->covariance(0,0);
		      cov_pv[1][0] = iteForD0->covariance(1,0);
 		      cov_pv[2][0] = iteForD0->covariance(2,0);
                      cov_pv[0][1] = cov_pv[1][0];
		      cov_pv[1][1] = iteForD0->covariance(1,1);
		      cov_pv[2][1] = iteForD0->covariance(2,1);
                      cov_pv[0][2] = cov_pv[2][0];
                      cov_pv[1][2] = cov_pv[2][1];
		      cov_pv[2][2] = iteForD0->covariance(2,2);
                      // distance direction unit vector
                      deriv[0] = vdx/dist;
                      deriv[1] = vdy/dist;
                      deriv[2] = vdz/dist;
                      // decay length error
                      dlerr=0;
                      for (int m=0; m<3; ++m){
                        for (int n=0; n<3; ++n){
                          dlerr += deriv[m]*deriv[n]*(cov_pv[m][n]+cov_sv[m][n]);
		        } // end for
		      } // end for
		      dlerrD0 = sqrt(dlerr);
		      // decay length significance
                      dlsigD0 = dlD0/dlerrD0;

                      // decay length error in xy
                      dlxyerr=0;
                      for (int m=0; m<2; ++m){
                        for (int n=0; n<2; ++n){
                          dlxyerr += deriv[m]*deriv[n]*(cov_pv[m][n]+cov_sv[m][n]);
		        } // end for
		      } // end for
		      dlxyerrD0 = sqrt(dlxyerr);
		      // decay length significance
                      dlxysigD0 = dlxyD0/dlxyerrD0;
     
                    } // end of cosphi cut
                 } // end of if VtxFound2

                 // calculate the decay length and get the vertex chi2
   //  		 DLxD0 = 0; DLyD0 = 0; dlxyD0 = 0;
   //            chi2D0 = 999.;
                 //cout << "before valid" << endl; 
   		 
                 if (myVertex.isValid() && VtxFound2) {
   //			
   //		   DLxD0 = (myVertex.position()).x() - iteForD0->x();
   //		   DLyD0 = (myVertex.position()).y() - iteForD0->y();
   //		   dlxyD0 = ((paxisD0[0] * DLxD0) + (DLyD0 * paxisD0[1])) / (p4D012.Pt());
			// DLmodD0 = sqrt(DLxD0*DLxD0 + DLyD0*DLyD0);
/*  From workbook:
    Although TransientVertex::originalTracks().size() will tell you the number of tracks in the vertex, many of these may have very small weights, meaning that they are unlikely to belong to it. Therefore you should check the individual tracks weights using TransientVertex::trackWeight(trk).
    TransientVertex::degreesOfFreedom() is defined as twice the sum of the tracks weights minus 3. If many of the track weights are small, this can be negative ! Similarly, the vertex chi2 is calculated including the weights. It can therefore be very different from that obtained using a KalmanVertexFit. There is no good reason why it should have a chi2 distribution.
    Vertices will only be kept if at least two tracks have weights exceeding parameter weightthreshold. Note that some of the found vertices may therefore have many/all tracks with very small weights. These should be treated with some distrust.

 The vertex returned may not be valid in some cases. The user had to check the validity of the vertex with the method isValid(). In each case, an error message is put into the log:

    The maximum number of iterations is exceeded
    The fitted position is out of the tracker bounds
    Too many tracks have been downweighted, and fewer than two significant tracks remain.   
*/
		   // preliminary, does not look very useful
		   // chi2DstarD0 = myVertex.degreesOfFreedom();
                   chi2D0 = myVertex.totalChiSquared();

		 }
		 // good!!! 

                 //cout << "after valid" << endl; 
	     
              // go to the next pair in loop here if no useful vertex was found
              // does not seem to reject anything?  
	      if (!(myVertex.isValid() && VtxFound2)) continue;

              // calculate rapidity
              rapD0 = log((sqrt(mD0Actual*mD0Actual+p4D012.P()*p4D012.P())+p4D012.Pz())/(sqrt(mD0Actual*mD0Actual+p4D012.P()*p4D012.P())-p4D012.Pz()))/2.;

              // check whether D0 candidate uses muons
              bool hasmuD0 = false;
              for (uint mm = 0; mm<nMuon; ++mm) {
                if (it1count == Muon_trkIdx[mm] || it2count == Muon_trkIdx[mm]) hasmuD0 = true;
              }

              // get dEdx info (Giacomo Fedi, mail CERN 12.6.19)
	      //Is the dedx collection available in AOD? this is the code we used in Run1:
	      //  Handle<DeDxDataValueMap> energyLossHandle;
	      //  iEvent.getByLabel("dedxHarmonic2", energyLossHandle);
	      //  const DeDxDataValueMap &  eloss  = *energyLossHandle;
	      //  double dedxTrk = eloss[trk1Ref].dEdx();
	      //  double errdedxTrk = eloss[trk1Ref].dEdxError();
	      //  int NumdedxTrk = eloss[trk1Ref].numberOfMeasurements();
              // and this is something more recent
              //  https://cmssdt.cern.ch/lxr/source/DQM/TrackingMonitor/src/dEdxAnalyzer.cc#0172
	      //   if (doDeDxPlots_ || doAllPlots_) {
	      //0162     edm::Handle<reco::TrackCollection> trackCollectionHandle;
	      //0163     iEvent.getByToken(trackToken_, trackCollectionHandle);
	      //0164     if (!trackCollectionHandle.isValid())
	      //0165       return;
	      //0166 
	      //0167     for (unsigned int i = 0; i < dEdxInputList_.size(); i++) {
	      //0168       edm::Handle<reco::DeDxDataValueMap> dEdxObjectHandle;
	      //0169       iEvent.getByToken(dEdxTokenList_[i], dEdxObjectHandle);
	      //0170       if (!dEdxObjectHandle.isValid())
	      //0171         continue;
	      //0172       const edm::ValueMap<reco::DeDxData> dEdxColl = *dEdxObjectHandle.product();
	      //0173 
	      //0174       for (unsigned int t = 0; t < trackCollectionHandle->size(); t++) {
	      //0175         reco::TrackRef track = reco::TrackRef(trackCollectionHandle, t);
	      //0176 
	      //0177         if (track->quality(reco::TrackBase::highPurity)) {
	      //0178           //MIPs
	      //0179           if (track->pt() >= 5.0 && track->numberOfValidHits() > TrackHitMin) {
	      //0180             dEdxMEsVector[i].ME_MipDeDx->Fill(dEdxColl[track].dEdx());
	      //0181             dEdxMEsVector[i].ME_MipDeDxNHits->Fill(dEdxColl[track].numberOfMeasurements());
	      //0182             if (dEdxColl[track].numberOfMeasurements() != 0)
	      //0183               dEdxMEsVector[i].ME_MipDeDxNSatHits->Fill((1.0 * dEdxColl[track].numberOfSaturatedMeasurements()) 
	      //0184                                                         dEdxColl[track].numberOfMeasurements());
	      //0185             dEdxMEsVector[i].ME_MipDeDxMass->Fill(mass(track->p(), dEdxColl[track].dEdx()));
	      //0186 
	      //0187             if (track->pt() >= HighPtThreshold) {
	      //0188               dEdxMEsVector[i].ME_MipHighPtDeDx->Fill(dEdxColl[track].dEdx());
	      //0189               dEdxMEsVector[i].ME_MipHighPtDeDxNHits->Fill(dEdxColl[track].numberOfMeasurements());
	      //0190             }
	      //0191 
	      //0192             //HighlyIonizing particles
	      //0193           } else if (track->pt() < 2 && dEdxColl[track].dEdx() > HIPdEdxMin) {
	      //0194             dEdxMEsVector[i].ME_HipDeDxMass->Fill(mass(track->p(), dEdxColl[track].dEdx()));
	      //0195           }
	      //0196         }
	      //0197       }
	      //0198     }
	      //0199   }
	      //and the name of the collection is: "dedxHarm2"
              
              // A.G.: Have seen dEdx info from miniAOD somewhere?

              // print out interesting candidates
              if (Bingo == 1 && mD012 > mD0tmin && mD012 < mD0tmax) {
                cout << "event " << event << " D0 " << mD012 << " pt " << p4D012.Pt() << " eta " << p4D012.Eta() << " phi " << p4D012.Phi() << " z " << zD0 << " dl " << dlD0 << " signif " << dlsigD0 << " cosphi " << cosphiD0 << endl;
              }
              if (Bingo == 1 && mD021 > mD0tmin && mD021 < mD0tmax) {
                cout << "event " << event << " D0 " << mD021 << " pt " << p4D021.Pt() << " eta " << p4D021.Eta() << " phi " << p4D021.Phi() << " z " << zD0 << " dl " << dlD0 << " signif " << dlsigD0 << " cosphi " << cosphiD0 << endl;
              }

              // check for associated simulated D0 from prestored list
              simidD0 = -1;
              t1pdgid = 0;
              t2pdgid = 0;
              if (Bingo) {
                cout << "nD0sim " << nD0sim << endl;
                if (nD0sim>0) {
                  cout << "pts " << ptD0sim[0] << " eta " << etaD0sim[0] << " phi " << phiD0sim[0] << endl;
                  cout << "pt " << p4D012.Pt() << " eta " << p4D012.Eta() << " phi " << p4D012.Phi() << endl;
                  cout << "pts1 " << ptD0t1sim[0] << " eta " << etaD0t1sim[0] << " phi " << phiD0t1sim[0] << " charge " << chgD0t1sim[0] << endl;
                  cout << "pts2 " << ptD0t2sim[0] << " eta " << etaD0t2sim[0] << " phi " << phiD0t2sim[0] << " charge " << chgD0t2sim[0] << endl;
                  cout << "pt1 " << it1->pt() << " eta " << it1->eta() << " phi " << it1->phi() << " charge " << it1->charge() << endl;
                  cout << "pt2 " << it2->pt() << " eta " << it2->eta() << " phi " << it2->phi() << " charge " << it2->charge() << endl;
                }
              }
              for (int iD0=0; iD0<nD0sim; ++iD0) {
                if (abs(p4D012.Pt()-ptD0sim[iD0])/ptD0sim[iD0] < 0.2 && abs(p4D012.Eta()-etaD0sim[iD0])<0.3 && fmod(abs(p4D012.Phi()-phiD0sim[iD0]),2.*pi)<0.3) {
#ifdef Bingo
                  cout << "D0 found" << endl;
#endif
                  if (it1->charge()==chgD0t1sim[iD0] && it2->charge()==chgD0t2sim[iD0]) {
                    // first charge alternative
                    if (abs(it1->pt()-ptD0t1sim[iD0])/ptD0t1sim[iD0] < 0.2 && abs(it1->eta()-etaD0t1sim[iD0])<0.3 && fmod(abs(it1->phi()-phiD0t1sim[iD0]),2.*pi)<0.3 && abs(it2->pt()-ptD0t2sim[iD0])/ptD0t2sim[iD0] < 0.2 && abs(it2->eta()-etaD0t2sim[iD0])<0.3 && fmod(abs(it2->phi()-phiD0t2sim[iD0]),2.*pi)<0.3) {
                      simidD0 = idD0sim[iD0];
                      t1pdgid = pdgIdD0t1sim[iD0];
                      t2pdgid = pdgIdD0t2sim[iD0];
#ifdef Bingo
                      cout << "D0 identified " << mD012 << " " << mD021 << " " << mD0KK << " " << pdgIdD0t1sim[iD0] << " " << pdgIdD0t2sim[iD0] << endl;
#endif
                      break;
		    }
		  }
                  else if (it2->charge()==chgD0t1sim[iD0] && it1->charge()==chgD0t2sim[iD0]) {
                    // second charge alternative
                    if (abs(it2->pt()-ptD0t1sim[iD0])/ptD0t1sim[iD0] < 0.2 && abs(it2->eta()-etaD0t1sim[iD0])<0.3 && fmod(abs(it2->phi()-phiD0t1sim[iD0]),2.*pi)<0.3 && abs(it1->pt()-ptD0t2sim[iD0])/ptD0t2sim[iD0] < 0.2 && abs(it1->eta()-etaD0t2sim[iD0])<0.3 && fmod(abs(it1->phi()-phiD0t2sim[iD0]),2.*pi)<0.3) {
                      simidD0 = idD0sim[iD0];
                      t1pdgid = pdgIdD0t1sim[iD0];
                      t2pdgid = pdgIdD0t2sim[iD0];
#ifdef Bingo
                      cout << "D0 identified " << mD012 << " " << mD021 << " " << mD0KK << " " << pdgIdD0t1sim[iD0] << " " << pdgIdD0t2sim[iD0] << endl;
#endif
                      break;
		    }
		  }
                  else if (it1->charge()==it2->charge() && it1->charge()==chgD0t2sim[iD0]) {
                    // allow also wrong charge to check fake rate (should we use both combinations?)
                    if (abs(it2->pt()-ptD0t1sim[iD0])/ptD0t1sim[iD0] < 0.2 && abs(it2->eta()-etaD0t1sim[iD0])<0.3 && fmod(abs(it2->phi()-phiD0t1sim[iD0]),2.*pi)<0.3 && abs(it1->pt()-ptD0t2sim[iD0])/ptD0t2sim[iD0] < 0.2 && abs(it1->eta()-etaD0t2sim[iD0])<0.3 && fmod(abs(it1->phi()-phiD0t2sim[iD0]),2.*pi)<0.3) {
                      simidD0 = idD0sim[iD0];
                      t1pdgid = pdgIdD0t1sim[iD0];
                      t2pdgid = pdgIdD0t2sim[iD0];
#ifdef Bingo
                      cout << "fake D0 identified " << mD012 << " " << mD021 << " " << mD0KK << " " << pdgIdD0t1sim[iD0] << " " << pdgIdD0t2sim[iD0] << endl;
#endif
                      break;
		    }
		  }
		}
              }

              bool storeD0 = false;  // flag whether D0 candidate is stored
	      if (D0_pt.size() < nReserve_D0) {

		// fill unlike sign (only), if extra cuts satisfied
                // store both positive and negative dl to allow mirroring
		if (((mD012 > mD0min && mD012 < mD0max) ||
		     (mD021 > mD0min && mD021 < mD0max)) &&
		    it1->charge() != it2->charge() && 
                    abs(dlxyD0)>dlxyD0min && abs(dlsigD0)>dlsigD0min) {
		storeD0 = true;
		D0t1_pt.push_back(it1->pt());
		D0t1_eta.push_back(it1->eta());
		D0t1_phi.push_back(it1->phi());
		D0t1_chg.push_back(it1->charge());
		//D0t1_tkIdx.push_back(-9999);
		D0t1_tkIdx.push_back(it1count);
		D0t1_Kprob.push_back(t1Kprob);    // temporarily strip
		D0t1_piprob.push_back(t1piprob);  // temporarily error or pix
                D0t1_dEdxnmeas.push_back(dedxtrknmeas);
                D0t1_dEdxnsat.push_back(dedxtrknsat);
		D0t1_vtxIdx.push_back(t1vtxid);
		D0t1_chindof.push_back(t1chindof);
		D0t1_nValid.push_back(t1nvalid);
		D0t1_nPix.push_back(t1npix);
		D0t1_isHighPurity.push_back(t1highp);
		D0t1_dxy.push_back(t1dxy);
		D0t1_dz.push_back(t1dz);
                D0t1_pdgId.push_back(t1pdgid);

		D0t2_pt.push_back(it2->pt());
		D0t2_eta.push_back(it2->eta());
		D0t2_phi.push_back(it2->phi());
		D0t2_chg.push_back(it2->charge());
		//D0t2_tkIdx.push_back(-9999);
		D0t2_tkIdx.push_back(it2count);
		D0t2_Kprob.push_back(t2Kprob);
		D0t2_piprob.push_back(t2piprob);
                D0t2_dEdxnmeas.push_back(dedxtrk2nmeas);
                D0t2_dEdxnsat.push_back(dedxtrk2nsat);
		D0t2_vtxIdx.push_back(t2vtxid);
		D0t2_chindof.push_back(t2chindof);
		D0t2_nValid.push_back(t2nvalid);
		D0t2_nPix.push_back(t2npix);
		D0t2_isHighPurity.push_back(t2highp);
		D0t2_dxy.push_back(t2dxy);
		D0t2_dz.push_back(t2dz);
                D0t2_pdgId.push_back(t2pdgid);

		// D0
		D0_pt.push_back(p4D012.Pt());
		D0_eta.push_back(p4D012.Eta());
		D0_phi.push_back(p4D012.Phi());
		D0_rap.push_back(rapD0);
		D0_mass12.push_back(mD012);
		D0_mass21.push_back(mD021);
		D0_massKK.push_back(mD0KK);
		D0_simIdx.push_back(simidD0);
                D0_DstarIdx.push_back(-1); // might be superseded later
                // will contain pointer to last reconstructed Dstar for D0 
                D0_ambiPrim.push_back(ambiprimary);
		D0_vtxIdx.push_back(vtxidD0);
		D0_hasMuon.push_back(hasmuD0);
                D0_chi2.push_back(chi2D0);
		D0_dlxy.push_back(dlxyD0);
		D0_dlxyErr.push_back(dlxyerrD0);
		D0_dlxySig.push_back(dlxysigD0);
		D0_cosphixy.push_back(cosphixyD0);
		D0_dl.push_back(dlD0);
		D0_dlErr.push_back(dlerrD0);
		D0_dlSig.push_back(dlsigD0);
		D0_cosphi.push_back(cosphiD0);
		D0_ptfrac.push_back(zD0);
		D0_ptfrac15.push_back(zD015);
		D0_ptfrac10.push_back(zD010);
		D0_ptfrac07.push_back(zD007);
		D0_ptfrac04.push_back(zD004);
		D0_x.push_back(D0x);
		D0_y.push_back(D0y);
		D0_z.push_back(D0z);
		D0_Covxx.push_back(cov_sv[0][0]);
		D0_Covyx.push_back(cov_sv[1][0]);
		D0_Covzx.push_back(cov_sv[2][0]);
		D0_Covyy.push_back(cov_sv[1][1]);
		D0_Covzy.push_back(cov_sv[2][1]);
		D0_Covzz.push_back(cov_sv[2][2]);
               }
	      }
	      else {cout << "WARNING!!!!! NO. OF D0 IS MORE THAN YOUR RESERVED NO.!!!!!!" << endl;}
              nD0 = D0_pt.size();

	      ///////////////////////
              ////     now D*    ////
              ///////////////////////


	      // check for track collection with at least 3 tracks (for Dstar)
              // reconstructing D*'s with D0 pt < 1 is hopeless (slow pion) 
	      if (tracks->size() >= trkminDstar && p4D012.Pt()>DstarD0ptmin) {

		//--- 3rd loop:over tracks assuming they are slow pions ---//
                // no vertex requirement since slow piosn are often not 
                // associated to primary, and actually should not be 
                // if D*s from B are to be included 
                int it3count = -1;

#ifndef miniAOD
                // loop over AOD track collection
		for (reco::TrackCollection::const_iterator itPS3 = tracks->begin(); itPS3 != tracks->end() ; ++itPS3) {
                  ++it3count;
#endif

#ifdef miniAOD
                // loop over both PackedCandidates 
                //            and lost track miniAOD collections
		// assume that standard collection is always larger than "lost" collection
		// cout << "tracks " << tracks->size() << " lost " << lostTracks->size() << endl;
		// alternatingly provide "good" and "lost" tracks, as long as there are any
                // initialize "lost" tracks
                bool lostended = false;
                pat::PackedCandidateCollection::const_iterator iclost = lostTracks->begin();
                if (iclost == lostTracks->end()) lostended = true;
                // loop over "good" tracks
		for (pat::PackedCandidateCollection::const_iterator icPS3 = tracks->begin(); icPS3 != tracks->end() ; ++icPS3) {
                  ++it3count;
                 // actually, the following shows that low pt tracks are 
                 // stored if associated to a vertex, but only few details 
                 // available. Find vertex and compare to D0 vertex?
		 //if (icPS3->charge()!=0) {
                 //  cout << icPS3->pt() << endl;
                 //  cout << icPS3->vertexRef()->z() << endl;
                 //  cout << icPS3->pvAssociationQuality() << endl;
		 // pvAssociationQuality() method returns the quality of PV-candidate association, in particular (please note that in the following the PV means the associated PV returned by vertexRef()):

		 // UsedInFitTight = 7: the track is used in the PV fit and the weight is above 0.5 -> OK
		 // UsedInFitLoose = 6: the track is used in the PV fit and the weight is below 0.5 -> OK
		 // CompatibilityDz = 5: the track is not used in fit but is very close in dZ to the PV (significance of dZ < 5 and dZ < 1mm). This is used only if the track dz uncertainty is below 500um. -> OK
		 // CompatibilityBTag =4: the track is not compatible with the PV but it is close to the nearest jet axis starting from the PV (distance to jet axis < 700um). This is used only if the track dz uncertainty is below 500um. -> OK but will work only in case of jet
		 // NotReconstructedPrimary = 0: the track is not associated to any PV and is compatible with the BeamSpot hence it is likely to be originating from an interaction for which we did not reconstruct the PV (beamspot compatiblity: dxy sig < 2sigma and dxy < 200um) -> not useful?
                 // OtherDeltaZ = 1: none of the above criteria is satisfied, hence the closest in dZ vertex is associated -> not useful?
                 // require icPS3->pvAssociationQuality() > 1?
                 //     no -> will catch only primary vertex ... 
		 //}  // ICPS3
                 // icPS3->dz() gives dz to main PV, not to actual PV!
                 //
                 // loop over "good" and "lost" track alternatives 
		 for (int ialt = 1; ialt<3; ++ialt) {
		   //if (!(icPS3->hasTrackDetails())) continue;
		   //auto tkPS3 = icPS3->pseudoTrack();
		  auto itPS3 = icPS3->bestTrack(); 
		  if (ialt == 1) {    
                    // use "good" track unless null
                    //if (itPS3 != nullptr)
         	    //   cout << ialt << " pt pis " << itPS3->pt() << endl; 
                    // actually, use all PF candidates!
                  }
		  else if (!lostended) {
                    // use "lost" track, if any remain 
                    // this seems to have a pt cut of 0.5 GeV, and higher for track info, do we need these at all?  
         	    // if (iclost->charge()!=0) cout << ialt << " pt pis c " << iclost->pt() << endl; 
         	    //if (!(iclost->hasTrackDetails())) continue;
                    //tkPS3 = iclost->pseudoTrack();
		    itPS3 = iclost->bestTrack(); 
                    ++iclost;
                    if (iclost == lostTracks->end()) lostended = true; 
                    //if (itPS3 != nullptr)
         	    //  cout << ialt << " pt pis " << itPS3->pt() << endl; 
                    // again, use all candidates
                  }
                  else continue;
#endif

                  // replace the variables used below by variables precomputed 
                  // above, separately for the AOD and miniAOD case?
                  // three categories:
                  //  - tracks with details (covariance)
                  //  - tracks without details (hits and position)
                  //  - four-vector and vertex ref + flag only
    
                  float t3pt, t3ptr;
#ifndef miniAOD
                  // track iterator always defined on AOD
                  t3pt = itPS3->pt();
#endif
#ifdef miniAOD
                  // on miniAOD, check whether track pointer exists and how
                  if (itPS3 != nullptr) {
                    // track defined
		    // check that it is not t1 or t2		    
		    if (itPS3 == it1 || itPS3 == it2)  continue;
                    //cout << "hello it1it2" << endl; 
                    t3pt = itPS3->pt(); 
                  }
                  else {
                    // no track defined, PF candidate, require charged
                    if (icPS3->charge() == 0) continue;
                    //cout << "hello notrack" << endl;
                    // check that the surviving tracks are all 'Pions' (+-211) 
                    // cout << icPS3->pdgId() << endl;  // they are!
                    t3pt = icPS3->pt(); 
                  }
#endif

		    // cut on absolute pis pt
		    if (t3pt < ptpismin) continue;
                    //cout << "hello ptpismin " << t3pt << endl; 

                    // also cut on pis/D0 pt ratio
                    // always lies between .032 and 0.18
                    // (mpi/mD0 = 0.075)
                    if (t3pt/p4D012.Pt() < .03) continue; 
                    if (t3pt/p4D012.Pt() > .20) continue;
                    // for B0 -> D0pi, ask for pt>1 and pi/D0 pt > 0.33 ?

                    // store track quality variables //
                    int t3vtxid = -1;
                    float t3chindof;
                    int t3nvalid;
                    int t3npix;
                    float t3x;
                    float t3y;
                    float t3z;
                    float t3charge;
                    float t3eta, t3etar;
                    float t3phi, t3phir;
#ifndef miniAOD
                    t3chindof = itPS3->chi2()/itPS3->ndof();
                    t3nvalid = itPS3->hitPattern().numberOfValidHits();
                    t3npix = itPS3->hitPattern().numberOfValidPixelHits();
                    t3x = itPS3->vx();  
                    t3y = itPS3->vy();  
                    t3z = itPS3->vz();
                    t3charge = itPS3->charge();
                    t3eta = itPS3->eta();
                    t3phi = itPS3->phi();
#endif
#ifdef miniAOD
		    if (itPS3 != nullptr) {
                      // (reduced) track info exists (pt>0.4)
                      t3chindof = itPS3->normalizedChi2();
                      t3nvalid = itPS3->hitPattern().numberOfValidHits();
                      t3npix = itPS3->hitPattern().numberOfValidPixelHits();
                      t3x = itPS3->vx();  
                      t3y = itPS3->vy();  
                      t3z = itPS3->vz();
                      t3charge = itPS3->charge();
                      t3eta = itPS3->eta();
                      t3phi = itPS3->phi();
                    }
		    else {
                      // no (reduced) track info for this PF candidate (pt<0.4)
                      t3chindof = - icPS3->pvAssociationQuality();
                      t3nvalid = 0;
                      t3npix = 0;
                      // position of associated vertex, not track, 
                      // even if track is far away from it 
                      // -> increased background!
                      t3x = icPS3->vertexRef()->x(); 
                      t3y = icPS3->vertexRef()->y();  
                      t3z = icPS3->vertexRef()->z();
                      t3charge = icPS3->charge();
                      t3eta = icPS3->eta();
                      t3phi = icPS3->phi();
                    }    
#endif

                    float t3piprob = -9999.;
                    float t3Kprob = -9999.;
                    int dedxtrk3nmeas = -9999;
                    int dedxtrk3nsat = -9999;

#ifndef CMSSW7plus
                    // for dEdx (from example by G. Fedi)
                    // beware of potential slowdown!
	            reco::TrackRef trackRef = reco::TrackRef(tracks,it3count);   
                    double dedxtrk = (*energyLossHandle)[trackRef].dEdx();
                    double dedxtrkerr = (*energyLossHandle)[trackRef].dEdxError();
                    dedxtrk3nsat = (*energyLossHandle)[trackRef].numberOfSaturatedMeasurements();
                    dedxtrk3nmeas = (*energyLossHandle)[trackRef].numberOfMeasurements();
                    // temporarily "misuse" the piprob and Kprob variables for direct dEdx:
                    t3piprob = dedxtrkerr;    
                    t3Kprob = dedxtrk;
#endif    
#ifdef CMSSW7plus
                    // for dEdx (from example by M. Soares)
                    // beware of potential slowdown!
	            reco::TrackRef trackRef = reco::TrackRef(tracks,it3count);  
                    double stripMapdedx = ( stripmap >0 ? (*dedxStMap)[trackRef].dEdx() : 0.);
                    double pixMapdedx = ( pixmap >0 ? (*dedxPixMap)[trackRef].dEdx() : 0.);
                    // temporarily "misuse" the piprob and Kprob variables for direct dEdx:
                    t3piprob = pixMapdedx;    
                    t3Kprob = stripMapdedx;    
                    dedxtrk3nsat = ( stripmap >0 ? (*dedxStMap)[trackRef].numberOfSaturatedMeasurements() : 0.);
	            dedxtrk3nmeas = ( stripmap >0 ? (*dedxStMap)[trackRef].numberOfMeasurements() : 0.);
#endif


                    // this would be the place to cut on these variables, if wished! //
                    
		    vc2[0] = abs(vcD0[0] - t3x);
		    vc2[1] = abs(vcD0[1] - t3y);
                    float t3dxy = sqrt((vc2[0] * vc2[0]) + (vc2[1] * vc2[1]));
		    vc2[2] = abs(vcD0[2] - t3z);
                    float t3dz = vcD0[2] - t3z;
 		    //cout << "ATAU CRASH SINI 2.5?" << endl;

		    // check the slow pion originates from the same region as the D0
		    if (t3dxy < vD0pis_xymax && fabs(t3dz) < vD0pis_zmax) {
			
		      // Calc 4 momentum for slow pion using TLorentzVector
		      p4pis.SetPtEtaPhiM(t3pt, t3eta, t3phi, pimass);

                      // determine which is Kaon and which is pion 
                      // depending on pis charge, and recut on mass
                      // keep K-pi+pis+, K+pi-pis- (right sing)
                      // and  K-pi-pis+, K+pi+pis- (wrong sign)
                      // i.e. pis always has opposite sign to Kaon
                      // in wrong sign case, use both possible combinations
                      int Kis = 0;
                      if (t3charge != it1->charge()) {
                        if (t3charge == it2->charge()) { 
                        // first track is Kaon candidate (right sign) 
                          Kis = 1;
                          p4DstarD0 = p4D012;
                          if (mD012 < mD0min1 || mD012 > mD0max1) continue;
                        }
                        else {
			// wrong sign, allow both K pi combinations if satisfied
                          if (mD012 > mD0min1 && mD012 < mD0max1) {
                            // first satisfies mass cuts
                            Kis = 1;
                            p4DstarD0 = p4D012;
                            if (mD021 > mD0min1 && mD021 < mD0max1) {
			      // both satisfied, default is first, need special treatment for second
			      Kis =12;
			    }
                          }
			  else if (mD021 > mD0min1 && mD021 < mD0max1) {
                            // 2nd satisfies mass cuts
                            Kis = 2;
                            p4DstarD0 = p4D021;
                          }
                          else {
                            // none satisfies mass cuts
                            continue; 
                          }
			}
                      }
		      else if (t3charge != it2->charge()) { 
                        // second track is kaon candidate (right sign)
                        Kis = 2;
                        p4DstarD0 = p4D021;
                        if (mD021 < mD0min1 || mD021 > mD0max1) continue;
                      }
		      else {
                        // all three have same sign, reject
                        continue;
		      } 
		      p4Dstar = p4DstarD0 + p4pis;
                      // the 2nd just in case;	  
                      if (Kis==12) p4Dstar2 = p4D021 + p4pis;
                      // cout << " Kis " << Kis << " event " << event << " " << (event/2)*2 << endl;

                      // calculate rapidy
                      rapDstar = log((sqrt(mDstarActual*mDstarActual+p4Dstar.P()*p4Dstar.P())+p4Dstar.Pz())/(sqrt(mDstarActual*mDstarActual+p4Dstar.P()*p4Dstar.P())-p4Dstar.Pz()))/2.; 

		      // calculate D*-D0 mass difference;
		      deltam = (p4Dstar.M()) - p4DstarD0.M();
                      deltam2 = 9999.;
                      // the 2nd just in case;	  
		      if (Kis==12) deltam2 = (p4Dstar2.M()) - p4D021.M();	  

                      // apply cut
		      if (deltam > dmDstarmax && deltam2 > dmDstarmax) continue;

		      // cout << "KELUAR X 3?" << endl;

                      // should be moved outside loop when vertex available
                      // shouldn't sumpt1 and sumpt2 be the same?
                      // i.e., is it necessary to recalculate?
		      sumpt2=0;
#ifndef miniAOD
		      for (reco::TrackCollection::const_iterator itSum2 = tracks->begin(); itSum2 != tracks->end(); ++itSum2) {
#endif
#ifdef miniAOD
		      for (pat::PackedCandidateCollection::const_iterator icSum2 = tracks->begin(); icSum2 != tracks->end(); ++icSum2) {
         		if (!(icSum2->hasTrackDetails())) continue;
                        auto itSum2 = icSum2->bestTrack(); 
                        if (itSum2 == nullptr) continue;
#endif
			// check track origins
			vc3[0] = abs(vcD0[0] - itSum2->vx());
			vc3[1] = abs(vcD0[1] - itSum2->vy());
			vc3[2] = abs(vcD0[2] - itSum2->vz());
		
			if ((sqrt(vc3[0] * vc3[0]) + (vc3[1] * vc3[1])) < vD0trksum_xymax && vc3[2] < vD0trksum_zmax) {
			  
			  sumpt2 += abs(itSum2->pt()); // sum pt for all tracks
			} // end of Sumpt vertex check
		      } // end of itSum loop

		      // variable z for Dstar and D0 from Dstar
		      // old: zDstarD0 = (p4D0.Pt()) / sumpt2;
	              zD0 = (it1->pt() + it2->pt()) / sumpt2;
		      // old: zDstar = (p4Dstar.Pt()) / sumpt2;
	              zDstar = (it1->pt() + it2->pt() + t3pt) / sumpt2;

		      // cut on z
                      if (zDstar < zDstarmin) continue; 

		      // cout << "KELUAR X 4?" << endl;

		    // initialize
                    t3ptr = -1.;
                    t3etar = 9999.;
                    t3phir = 9999.;
                    deltamr = -1.; 
                    deltamr2 = -1.; 

#ifdef miniAOD
                  if (icPS3->hasTrackDetails()) {
#endif
		    // make the pislow more precise by fitting it to the 
                    // primary vertex (appropriate for primary D*)
		    // (CPU intensive! -> after all cuts)

                    // currently this uses a fudged variant of the beam 
                    // constraint
                    // could/should use SingleTrackVertexConstraint instead
//     https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideVertexFitTrackRefit

		    //cout << "pislow before: " << t3pt << " " << t3eta << " " << t3phi << endl; 

         	    // build pislow track collection for revertexing
		    mytracksforpislow.clear();
#ifndef miniAOD
                    // loop over prebuilt list of AOD transient tracks
         	    for (vector<TransientTrack>::iterator gt_trans = genralTracks_forD.begin(); gt_trans != genralTracks_forD.end(); ++gt_trans) {
					
		      const reco::TrackRef trackRef1 = (gt_trans->trackBaseRef()).castTo<reco::TrackRef>();
		      if (&*itPS3 == trackRef1.get()) {			    
		        TransientTrack  transientTrack1 = theB->build(trackRef1);
		        mytracksforpislow.push_back(transientTrack1);
		      }
		    }
#endif
#ifdef miniAOD 
                    // build transient track from relevant miniAOD pseudotrack
	            TransientTrack  transientTrack1 = theB->build(icPS3->pseudoTrack());
		    mytracksforpislow.push_back(transientTrack1);
#endif              

                    // want to fit track to primary vertex, but Fitter does not 
                    // seem to allow to pass a vertex constraint
                    // -> "misuse" bem spot constraint instead 

                    // build dummy beam spot from primary vertex
                    reco::BeamSpot dummyBeamSpot = vertexBeamSpot;  

                    // set primary vertex reference point pvXYZ
                    math::XYZPoint pvXYZ(iteForD0->position());

                    // get beam spot covariance matrix (to define variable)
		    reco::BeamSpot::CovarianceMatrix CovMatBeam = vertexBeamSpot.covariance(); 
                    // supersede 3x3 part of 7x7 covariance matrix by vertex info
                    // and set the remainder to 0 (should some info be kept?)
                    for (int icov = 0; icov < 7; ++icov) {
                      for (int jcov = 0; jcov < 7; ++jcov) {
                        if (icov < 3 && jcov < 3) CovMatBeam[icov][jcov] = iteForD0->covariance(icov,jcov);
                        else CovMatBeam[icov][jcov]=0;
                      }
                    }
                
                      // refill dummybeamspot with vertex info
                      dummyBeamSpot = BeamSpot(pvXYZ, sqrt(iteForD0->covariance(2,2)), vertexBeamSpot.dxdz(), vertexBeamSpot.dydz(), sqrt(iteForD0->covariance(0,0)), CovMatBeam, reco::BeamSpot::Unknown);

		      // do the revertexing of the slow pion 
		      KalmanVertexFitter theFitter3(true);
                      myPisVertex = theFitter3.vertex(mytracksforpislow,dummyBeamSpot);
		    // proceed only if fit was successful
		    if (myPisVertex.isValid()) {

		    // and if confidence level is reasonable
                  // vertex parameters
		  float ndofPis = myPisVertex.degreesOfFreedom();
                  // for Adaptive Vertex Fitter this sum of weights*2-3
                  // for Kalman Vertex Fitter this is 1 for two tracks?
                  chi2Pis = myVertex.totalChiSquared();
                  // confidence level
                  float CLPis=999.;
                  if (ndofPis >= 1) {
                    CLPis = TMath::Prob(chi2Pis,(int)ndofPis); 
                  }
                  // cut on chisquared (Adaptive) or CL (Kalman) here?
                  if (CLPis>0.01) {

                      // *** should cut on the chi2 here, or at least store it ... ***
                
                      // get updated track parameters
                      vector<TransientTrack> trackspis = myPisVertex.refittedTracks();
                      vector<TransientTrack>::const_iterator trackpsIt = trackspis.begin();
                      const Track & trackpis = trackpsIt->track();
                      t3ptr = trackpis.pt();
                      t3etar = trackpis.eta();
                      t3phir = trackpis.phi();

		      //cout << "pislow after: " << t3ptr << " " << t3etar << " " << t3phir << endl; 

		      // Calculate 4 momentum for slow pion using TLorentzVector
		      p4pisr.SetPtEtaPhiM(t3ptr, t3etar, t3phir, pimass);

                      // Calculate D* four-momentum 
		      p4Dstarr = p4DstarD0 + p4pisr;
                      // the 2nd just in case;	  
                      if (Kis==12) p4Dstarr2 = p4D021 + p4pisr;

		      // calculate D*-D0 mass difference;
		      deltamr = (p4Dstarr.M()) - p4DstarD0.M();
                      deltamr2 = 9999.;
                      // the 2nd just in case;	  
		      if (Kis==12) deltamr2 = (p4Dstarr2.M()) - p4D021.M();
		     } // CLPis
		    } // myPisVertex valid  
#ifdef miniAOD
		  } // hasTrackDetails
#endif

		      // cout << "ATAU CRASH SINI 4?" << endl;

  
//      ********  fill Dstar histograms *****************

		      // if ( (p4Dstar.Pt()) > 3.5 && itK1->pt() > 1 && itP2->pt() > 1 && itPS3->pt() > 0.25 && abs(deltam - 0.1454) < 0.001 && itP2->charge() != itK1->charge()) {h_Ncand_D0mass_alSpectrum->Fill(1);}
		      
		      // ATLAS cuts
		      // if (DmesonsDL && (p4Dstar.Pt()) > 3.5 && itK1->pt() > 1 && itP2->pt() > 1 && itPS3->pt() > 0.25 && abs(deltam - 0.1454) < 0.001) {
		      // if ((p4Dstar.Pt()) > 3.5 && itK1->pt() > 1 && itP2->pt() > 1 && itPS3->pt() > 0.25 && myVertex.isValid() && VtxFound2) { // cutmacro

                      // this is a dummy if to keep the structure ...
		      if (myVertex.isValid() && VtxFound2) {

        	        // Fill histograms for right sign only 
                        // (don't need to deal with wrong sign ambiguity)
			if (it1->charge() != it2->charge()) { // cutmacro
                          // fill data for all vertices, 
                          // fill MC only for main simulated vertex
		          if (run!=1 || (vtxidD0>-1 && PVtx_isMainSim[vtxidD0])) {
		            h_d0pt->Fill(p4DstarD0.Pt());
		            h_dstarpt->Fill(p4Dstar.Pt());
		            h_PS3pt->Fill(t3pt);
                            if (Kis == 1) {
                            // first track is Kaon candidate (K-pi+pi+ or K+pi-pi-) 
           		    h_K1pt->Fill(it1->pt());
		            h_P2pt->Fill(it2->pt());
		            h_K1eta->Fill(it1->eta());
		            h_P2eta->Fill(it2->eta());
                            }
                            else if (Kis ==2) {                
                              // second track is Kaon candidate 
                              h_K1pt->Fill(it2->pt());
		              h_P2pt->Fill(it1->pt());
		              h_K1eta->Fill(it2->eta());
		              h_P2eta->Fill(it1->eta());              
                            }
                            else {
                              cout << "*** ALARM !!! Kaon not assigned ***" << endl;
                            }
		            h_PS3eta->Fill(t3eta);
			  
         	            // h_Vertex_Multiplicity_D0mass_allSpectrum->Fill(1);
			    if ((p4Dstar.Pt()) > 3.5 && it1->pt() > 1 && it2->pt() > 1 && t3pt > 0.25) {

			      if (abs(deltam - dmDstarActual) < 0.001) {
			        h_D0masscut->Fill(p4DstarD0.M());

			        // right side decay length D0 from Dstar
			        if (dlxyD0 > 0.02) {
				  h_D0masscut_rightDLcut->Fill(p4DstarD0.M());
			        } // this cut need to put in macro
			      } // this cut need to put in macro

			      if (abs(p4DstarD0.M() - mD0Actual) < 0.025) {
			        h_deltaMassD0Dstar->Fill(deltam);
			    
			        if (dlxyD0 > 0.02) {
				  h_deltaMassD0Dstar_rightDLcut->Fill(deltam);
			        } // this cut need to put in macro
			      } // this cut need to put in macro (the mD0 must be mD0 from Dstar)
                            } // end of tight cuts
			  } // end of cut for histo
			  // cout << "ATAU CRASH SINI 6?" << endl;

			  // print out interesting candidates
                          if (Bingo == 1 && deltam > dmDstartmin && deltam < dmDstartmax) {
                            cout << "event " << event << " DstarD0 " << " pt " << p4DstarD0.Pt() << " eta " << p4DstarD0.Eta() << " phi " << p4DstarD0.Phi() <<  " z " << zD0 << " dl " << dlxyD0 << endl;
		          } // Bingo
		        } // end of unlike sign cut 

//  *** more selections/corrections *** 

                        // check for associated simulated Dstar from prestored list
                        simidDstar = -1;
                        if (Bingo) {
                          cout << "nDstarsim " << nD0sim << endl;
                          if (nDstarsim>0) {
                            cout << "pts " << ptDstarsim[0] << " eta " << etaDstarsim[0] << " phi " << phiDstarsim[0] << endl;
                            cout << "pt " << p4Dstar.Pt() << " eta " << p4Dstar.Eta() << " phi " << p4Dstar.Phi() << endl;
                          } //Dstarsim
                        } // Bingo
                        for (int iDstar=0; iDstar<nDstarsim; ++iDstar) {
                          if (abs(p4Dstar.Pt()-ptDstarsim[iDstar])/ptDstarsim[iDstar] < 0.2 && abs(p4Dstar.Eta()-etaDstarsim[iDstar])<0.3 && fmod(abs(p4Dstar.Phi()-phiDstarsim[iDstar]),2.*pi)<0.3) {
                            simidDstar = idDstarsim[iDstar];
#ifdef Bingo
                            cout << "Dstar identified " << p4Dstar.M() << endl;
#endif
                            break;
                          } // if
                        }  // for

                        // find out whether pis is from D0 primary 
                        // CPU-expensive! do only after all cuts
                        // loop over all tracks from D0 primary vertex 
                        // (does not work for miniAOD: empty loop)
			// *** should use primary track instead of general track? ***
                        for (reco::Vertex::trackRef_iterator iTrack = iteForD0->tracks_begin(); iTrack != iteForD0->tracks_end(); ++iTrack) {
                          // get track reference	    
		          const reco::TrackRef trackRef = iTrack->castTo<reco::TrackRef>();
                          if (trackRef.get() == &*itPS3) {
                            t3vtxid = ivtx;

			    //  cout << "pislow vertex: " << trackRef->pt() << " " << trackRef->eta() << " " << trackRef->phi() << endl;
                            // effective the same as the generaltrack
 
                          } // trackref
                        } // itrack

                        // check whether Dstar candidate uses muons
                        bool hasmuDstar = hasmuD0;
                        for (uint mm = 0; mm<nMuon; ++mm) {
                          if (it3count == Muon_trkIdx[mm]) hasmuDstar = true; 
                        }


//   ******************** fill Dstar ntuple info *******************
 
                        if (Dstar_pt.size() < nReserve_Dstar) {
				      
                          // Save the D0 from Dstar, kaon and pion after revertexing
                          if (deltam < dmDstarmax && 
                           ( (p4DstarD0.M() > mD0tmin && p4DstarD0.M() < mD0tmax) || 
                             (deltam > dmDstartmin && deltam < dmDstartmax) ) ) {
#ifdef miniAOD
			// some check for slow pion
                        //if (icPS3->pt()<0.25 && vtxidD0==0) {
                        //  cout << "pis pt " << icPS3->pt() << " dz " << icPS3->dz() << " vtx " << vtxidD0 << endl;
                        //}
#endif 
                            // D0 parameters 
                            DstarD0_pt.push_back(p4DstarD0.Pt());
                            DstarD0_eta.push_back(p4DstarD0.Eta());
                            DstarD0_phi.push_back(p4DstarD0.Phi());
                            DstarD0_mass.push_back(p4DstarD0.M());
                            DstarD0_chi2.push_back(chi2D0);
                            DstarD0_dlxy.push_back(dlxyD0);
                            DstarD0_dlxyErr.push_back(dlxyerrD0);
                            DstarD0_dlxySig.push_back(dlxysigD0);
                            DstarD0_cosphixy.push_back(cosphixyD0);
                            DstarD0_dl.push_back(dlD0);
                            DstarD0_dlErr.push_back(dlerrD0);
                            DstarD0_dlSig.push_back(dlsigD0);
                            DstarD0_cosphi.push_back(cosphiD0);
                            DstarD0_ptfrac.push_back(zD0);
                            DstarD0_ptfrac15.push_back(zD015);
                            DstarD0_ptfrac10.push_back(zD010);
                            DstarD0_ptfrac07.push_back(zD007);
                            DstarD0_ptfrac04.push_back(zD004);
                            DstarD0_x.push_back(D0x);
                            DstarD0_y.push_back(D0y);
                            DstarD0_z.push_back(D0z);
                            DstarD0_simIdx.push_back(simidD0);
                            if (storeD0) {
                              DstarD0_recIdx.push_back(nD0-1);
                              D0_DstarIdx[nD0-1]=Dstar_pt.size(); // not yet incremented
                            }
			    else {
                              DstarD0_recIdx.push_back(-1);
                            }
                            DstarD0_ambiPrim.push_back(ambiprimary);

                            // K and pi parameters
                            if (Kis == 1 || Kis == 12) {
                              // first track is Kaon candidate
                              DstarK_pt.push_back(it1->pt());
                              DstarK_eta.push_back(it1->eta());
                              DstarK_phi.push_back(it1->phi());
                              DstarK_chg.push_back(it1->charge());
                              DstarK_tkIdx.push_back(it1count);
                              DstarK_Kprob.push_back(t1Kprob);
                              DstarK_piprob.push_back(t1piprob);
                              DstarK_dEdxnmeas.push_back(dedxtrknmeas);
                              DstarK_dEdxnsat.push_back(dedxtrknsat);
                              DstarK_vtxIdx.push_back(t1vtxid);
                              DstarK_chindof.push_back(t1chindof); 
                              DstarK_nValid.push_back(t1nvalid);
                              DstarK_nPix.push_back(t1npix);
		              DstarK_isHighPurity.push_back(t1highp);
                              DstarK_dxy.push_back(t1dxy);
                              DstarK_dz.push_back(t1dz);

                              Dstarpi_pt.push_back(it2->pt());
                              Dstarpi_eta.push_back(it2->eta());
                              Dstarpi_phi.push_back(it2->phi());
                              Dstarpi_chg.push_back(it2->charge());
                              Dstarpi_tkIdx.push_back(it2count);
                              Dstarpi_Kprob.push_back(t2Kprob);
                              Dstarpi_piprob.push_back(t2piprob);
                              Dstarpi_dEdxnmeas.push_back(dedxtrk2nmeas);
                              Dstarpi_dEdxnsat.push_back(dedxtrk2nsat);
                              Dstarpi_vtxIdx.push_back(t2vtxid);
                              Dstarpi_chindof.push_back(t2chindof);
                              Dstarpi_nValid.push_back(t2nvalid);
                              Dstarpi_nPix.push_back(t2npix);
		              Dstarpi_isHighPurity.push_back(t2highp);
                              Dstarpi_dxy.push_back(t2dxy);
                              Dstarpi_dz.push_back(t2dz);
                            } // Kis
                            else if (Kis == 2) {
                              // second track is Kaon candidate
                              DstarK_pt.push_back(it2->pt());
                              DstarK_eta.push_back(it2->eta());
                              DstarK_phi.push_back(it2->phi());
                              DstarK_chg.push_back(it2->charge());
                              DstarK_tkIdx.push_back(it2count);
                              DstarK_Kprob.push_back(t2Kprob);
                              DstarK_piprob.push_back(t2piprob);
                              DstarK_dEdxnmeas.push_back(dedxtrk2nmeas);
                              DstarK_dEdxnsat.push_back(dedxtrk2nsat);
                              DstarK_vtxIdx.push_back(t2vtxid);
                              DstarK_chindof.push_back(t2chindof);
                              DstarK_nValid.push_back(t2nvalid);
                              DstarK_nPix.push_back(t2npix);
		              DstarK_isHighPurity.push_back(t2highp);
                              DstarK_dxy.push_back(t2dxy);
                              DstarK_dz.push_back(t2dz);

                              Dstarpi_pt.push_back(it1->pt());
                              Dstarpi_eta.push_back(it1->eta());
                              Dstarpi_phi.push_back(it1->phi());
                              Dstarpi_chg.push_back(it1->charge());
                              Dstarpi_tkIdx.push_back(it1count);
                              Dstarpi_Kprob.push_back(t1Kprob);
                              Dstarpi_piprob.push_back(t1piprob);
                              Dstarpi_dEdxnmeas.push_back(dedxtrknmeas);
                              Dstarpi_dEdxnsat.push_back(dedxtrknsat);
                              Dstarpi_vtxIdx.push_back(t1vtxid);
                              Dstarpi_chindof.push_back(t1chindof);
                              Dstarpi_nValid.push_back(t1nvalid);
                              Dstarpi_nPix.push_back(t1npix);
		              Dstarpi_isHighPurity.push_back(t1highp);
                              Dstarpi_dxy.push_back(t1dxy);
                              Dstarpi_dz.push_back(t1dz);
	                    } // Kis
                            else {
                              cout << "*** ALARM 2 !!! Kaon not assigned ***" << endl;
                            } // else Kis

                            // slow pion parameters
                            Dstarpis_pt.push_back(t3pt);
                            Dstarpis_eta.push_back(t3eta);
                            Dstarpis_phi.push_back(t3phi);
                            Dstarpis_ptr.push_back(t3ptr);
                            Dstarpis_etar.push_back(t3etar);
                            Dstarpis_phir.push_back(t3phir);
                            Dstarpis_chg.push_back(t3charge);
                            Dstarpis_tkIdx.push_back(it3count);
                            Dstarpis_Kprob.push_back(t3Kprob);
                            Dstarpis_piprob.push_back(t3piprob);
                            Dstarpis_dEdxnmeas.push_back(dedxtrk3nmeas);
                            Dstarpis_dEdxnsat.push_back(dedxtrk3nsat);
                            Dstarpis_vtxIdx.push_back(t3vtxid);
                            Dstarpis_chindof.push_back(t3chindof);
                            Dstarpis_chir.push_back(chi2Pis);
                            Dstarpis_nValid.push_back(t3nvalid);
                            Dstarpis_nPix.push_back(t3npix);
                            Dstarpis_dxy.push_back(t3dxy);
                            Dstarpis_dz.push_back(t3dz);
	
                            // Dstar parameters
                            Dstar_pt.push_back(p4Dstar.Pt());
                            Dstar_eta.push_back(p4Dstar.Eta());
                            Dstar_phi.push_back(p4Dstar.Phi());
                            Dstar_rap.push_back(rapDstar);
                            Dstar_deltam.push_back(deltam);
                            Dstar_deltamr.push_back(deltamr);
                            Dstar_simIdx.push_back(simidDstar);
                            Dstar_vtxIdx.push_back(vtxidD0);
                            Dstar_hasMuon.push_back(hasmuDstar);
                            Dstar_ptfrac.push_back(zDstar);
                          } //end of mass and deltam cut
                        } // end of reserve Dstar 
	                else {
                          cout << "WARNING!!!!! NO. OF D* IS MORE THAN YOUR RESERVED NO.!!!!!!" << endl;
                        } // else reserve Dstar
	                if (Kis==12) { 
                          // Save other candidate for wrong sign combinations
                          if (Dstar_pt.size() < nReserve_Dstar) {
                            p4DstarD0 = p4D021; 

                            // Save the D0 from Dstar, kaon and pion after revertexing
                            if (deltam2 < dmDstarmax &&
                               ( (p4DstarD0.M() > mD0tmin && p4DstarD0.M() < mD0tmax) || 
                               ( deltam2 > dmDstartmin && deltam2 < dmDstartmax) ) ) {

                              // D0 parameters 
                              DstarD0_pt.push_back(p4DstarD0.Pt());
                              DstarD0_eta.push_back(p4DstarD0.Eta());
                              DstarD0_phi.push_back(p4DstarD0.Phi());
                              DstarD0_mass.push_back(p4DstarD0.M());
                              DstarD0_chi2.push_back(chi2D0);
                              DstarD0_dlxy.push_back(dlxyD0);
                              DstarD0_dlxyErr.push_back(dlxyerrD0);
                              DstarD0_dlxySig.push_back(dlxysigD0);
                              DstarD0_cosphixy.push_back(cosphixyD0);
                              DstarD0_dl.push_back(dlD0);
                              DstarD0_dlErr.push_back(dlerrD0);
                              DstarD0_dlSig.push_back(dlsigD0);
                              DstarD0_cosphi.push_back(cosphiD0);
                              DstarD0_ptfrac.push_back(zD0);
                              DstarD0_ptfrac15.push_back(zD015);
                              DstarD0_ptfrac10.push_back(zD010);
                              DstarD0_ptfrac07.push_back(zD007);
                              DstarD0_ptfrac04.push_back(zD004);
                              DstarD0_x.push_back(D0x);
                              DstarD0_y.push_back(D0y);
                              DstarD0_z.push_back(D0z);
                              DstarD0_simIdx.push_back(simidD0);
                              if (storeD0) {
                                // this should actually never happen?
                                DstarD0_recIdx.push_back(nD0-1);
                                D0_DstarIdx[nD0-1]=Dstar_pt.size(); // not yet incremented 
                              }
	              	      else {
                                DstarD0_recIdx.push_back(-1);
                              }
                              DstarD0_ambiPrim.push_back(ambiprimary);

                              // second track is Kaon candidate
                              DstarK_pt.push_back(it2->pt());
                              DstarK_eta.push_back(it2->eta());
                              DstarK_phi.push_back(it2->phi());
                              DstarK_chg.push_back(it2->charge());
                              DstarK_tkIdx.push_back(it2count);
                              DstarK_Kprob.push_back(t2Kprob);
                              DstarK_piprob.push_back(t2piprob);
                              DstarK_dEdxnmeas.push_back(dedxtrk2nmeas);
                              DstarK_dEdxnsat.push_back(dedxtrk2nsat);
                              DstarK_vtxIdx.push_back(t2vtxid);
                              DstarK_chindof.push_back(t2chindof);
                              DstarK_nValid.push_back(t2nvalid);
                              DstarK_nPix.push_back(t2npix);
		              DstarK_isHighPurity.push_back(t2highp);
                              DstarK_dxy.push_back(t2dxy);
                              DstarK_dz.push_back(t2dz);

                              Dstarpi_pt.push_back(it1->pt());
                              Dstarpi_eta.push_back(it1->eta());
                              Dstarpi_phi.push_back(it1->phi());
                              Dstarpi_chg.push_back(it1->charge());
                              Dstarpi_tkIdx.push_back(it1count);
                              Dstarpi_Kprob.push_back(t1Kprob);
                              Dstarpi_piprob.push_back(t1piprob);
                              Dstarpi_dEdxnmeas.push_back(dedxtrknmeas);
                              Dstarpi_dEdxnsat.push_back(dedxtrknsat);
                              Dstarpi_vtxIdx.push_back(t1vtxid);
                              Dstarpi_chindof.push_back(t1chindof);
                              Dstarpi_nValid.push_back(t1nvalid);
                              Dstarpi_nPix.push_back(t1npix);
		              Dstarpi_isHighPurity.push_back(t1highp);
                              Dstarpi_dxy.push_back(t1dxy);
                              Dstarpi_dz.push_back(t1dz);

                              // slow pion parameters
                              Dstarpis_pt.push_back(t3pt);
                              Dstarpis_eta.push_back(t3eta);
                              Dstarpis_phi.push_back(t3phi);
                              Dstarpis_ptr.push_back(t3ptr);
                              Dstarpis_etar.push_back(t3etar);
                              Dstarpis_phir.push_back(t3phir);
                              Dstarpis_chg.push_back(t3charge);
                              Dstarpis_tkIdx.push_back(it3count);
                              Dstarpis_Kprob.push_back(t3Kprob);
                              Dstarpis_piprob.push_back(t3piprob);	
                              Dstarpis_dEdxnmeas.push_back(dedxtrk3nmeas);
                              Dstarpis_dEdxnsat.push_back(dedxtrk3nsat);
                              Dstarpis_vtxIdx.push_back(t3vtxid);
                              Dstarpis_chindof.push_back(t3chindof);  
                              Dstarpis_chir.push_back(chi2Pis);  
                              Dstarpis_nValid.push_back(t3nvalid);
                              Dstarpis_nPix.push_back(t3npix);
                              Dstarpis_dxy.push_back(t3dxy);
                              Dstarpis_dz.push_back(t3dz);

                              // Dstar parameters (2nd candidate)
                              Dstar_pt.push_back(p4Dstar2.Pt());
                              Dstar_eta.push_back(p4Dstar2.Eta());
                              Dstar_phi.push_back(p4Dstar2.Phi());
                              Dstar_rap.push_back(rapDstar);
                              Dstar_deltam.push_back(deltam2);
                              Dstar_deltamr.push_back(deltamr2);
                              Dstar_simIdx.push_back(simidDstar);
                              Dstar_vtxIdx.push_back(vtxidD0);
                              Dstar_hasMuon.push_back(hasmuDstar);
                              Dstar_ptfrac.push_back(zDstar);
                            } //end of mass and deltam cut
                          } // end of reserve Dstar 
                          else {
                            cout << "WARNING 2!!!!! NO. OF D* IS MORE THAN YOUR RESERVED NO.!!!!!!" << endl;
                          } // else reserve Dstar
			} // end of Kis == 12
		      } // end of dummy/ATLAS cuts check
		    } // if (sqrt(vc2[0]*vc2[0] +......
    //		  } // end of not t1 or t2 check
#ifdef miniAOD
         	 } // end of ialt 
#endif
		} // end of PS3 loop
	      } // end of track size >=3
	    } // if ( sqrt (vc[0]*vc[0] + vc[1]*vc[1])<0.1 && vc[2] < 0.1)
	  } // if ( itP2 !=itK1 && itP2->pt() > 0.5 )
	} // end of t2 loop			    
      } //end of t1 Pt cut
    } // end of track collection
  } // end of track size >=2
  nDstar = Dstar_pt.size();

  // cout << "IS DMESON LOOP OK?" << endl;
  
//////////////////////////////////////////////////////////////////////////////
//////////////////////////// D* Meson Analysis End ///////////////////////////
//////////////////////////////////////////////////////////////////////////////
	} // nanoext

  //cout << "fill tree" << endl;

  // Fill the tree for muon and dmeson
  t_event->Fill();

  //cout << "hello analysis end" << endl;

}//NanoAnalyzer::analyze ends

//**************************************************
//---------------------------Actual trigger analysis-------------  Qun below
//**************************************************
void NanoAnalyzer::analyzeTrigger(const edm::Event& iEvent, const edm::EventSetup& iSetup, const std::string& triggerName)
{
  using namespace std;
  using namespace edm;
  using namespace reco;
  using namespace trigger;

  cout<<"Currently analyzing trigger "<<triggerName<<endl;
  //Check the current configuration to see how many total triggers there are
  const unsigned int n(hltConfig_.size());
  //Get the trigger index for the current trigger
  const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName));
  //check that the trigger in the event and in the configuration agree
  assert(triggerIndex==iEvent.triggerNames(*triggerResultsHandle_).triggerIndex(triggerName));
  // abort on invalid trigger name
  if (triggerIndex>=n) {
    cout << "HLTEventAnalyzerAOD::analyzeTrigger: path "
         << triggerName << " - not found!" << endl;
    return;
  }
  //else {cout << "HLTEventAnalyzerAOD inside loop QQQQQ " << triggerName_ << endl;}

  //const std::pair<int,int> prescales(hltConfig_.prescaleValues(iEvent,iSetup,triggerName));
  //cout << "HLTEventAnalyzerAOD::analyzeTrigger: path "
  //    << triggerName << " [" << triggerIndex << "] "
  //    << "prescales L1T,HLT: " << prescales.first << "," << prescales.second
  //    << endl;

  //  One could find the list of modules in a given trigger from the HLTConfigProvider as explained somewhere else in this code.
  //  Get index (slot position) of module giving the decision of the path as described in "DataFormats/Common/interface/HLTGlobalStatus.h"
  const unsigned int m(hltConfig_.size(triggerIndex));
  const vector<string>& moduleLabels(hltConfig_.moduleLabels(triggerIndex));
  const unsigned int moduleIndex(triggerResultsHandle_->index(triggerIndex));
  //cout << " Last active module - label/type: "
  //     << moduleLabels[moduleIndex] << "/" << hltConfig_.moduleType(moduleLabels[moduleIndex])
  //     << " [" << moduleIndex << " out of 0-" << (m-1) << " on this path]"
  //     << endl;

  assert (moduleIndex<m);


  for (unsigned int j=0; j<=moduleIndex; ++j) {
    const string& moduleLabel(moduleLabels[j]);
    const string  moduleType(hltConfig_.moduleType(moduleLabel));

    const unsigned int filterIndex(triggerEventHandle_->filterIndex(InputTag(moduleLabel,"",processName_)));
    if (filterIndex<triggerEventHandle_->sizeFilters()) {
//      cout << " 'L3' filter in slot " << j << " - label/type " << moduleLabel << "/" << moduleType << endl;
      const Vids& VIDS (triggerEventHandle_->filterIds(filterIndex));
      const Keys& KEYS(triggerEventHandle_->filterKeys(filterIndex));
      const size_type nI(VIDS.size());
      const size_type nK(KEYS.size());
      cout << "QQQ before  " << " accepted 'L3' objects found: " << " nI " << nI << " nK " << nK << endl;
      assert(nI==nK);
      const size_type n(max(nI,nK));
      cout << "QQQ   " << n  << " accepted 'L3' objects found: " << " nI " << nI << " nK " << nK << endl;
      const TriggerObjectCollection& TOC(triggerEventHandle_->getObjects());
      for (size_type i=0; i!=n; ++i) {
        const TriggerObject& TO(TOC[KEYS[i]]);
        cout << "QQQ Obj   " << i << " " << VIDS[i] << "/" << KEYS[i] << ": id, pt, eta, phi, mass"
             << TO.id() << " " << TO.pt() << " " << TO.eta() << " "
             << TO.phi() << " " << TO.mass()
             << endl;
      }



    } //end of filterIndex
  } //end of moduleIndex
}// NanoAnalyzer::analyzeTrigger ends
//above Qun


//**************************************************
//************* additional methods *****************
//**************************************************

// ------------ method called once each job just before starting event loop  ------------
void 
//NanoAnalyzer::beginJob(const edm::Event& iEvent)
NanoAnalyzer::beginJob()
{
  // use this if your c++ < 14 (depend on your CMSSW)
  // file = unique_ptr<TFile>(new TFile("nAODish.root", "recreate")); // check
  
  // CHANGE
  /* 
  
  file = unique_ptr<TFile>(new TFile(outFile.c_str(), "recreate")); // check
  
  // use this if your c++ => 14
  // make_unique<TFile>("nAODish.root", "recreate");
  // use this if you didn't use unique pointer
  // file = new TFile("nAODish.root", "recreate");
  
  t_event = unique_ptr<TTree>(new TTree("Events", "Events"));
  
  // h_ = uniques_ptr<TH1D>(new TH1D("h_D0masscut", "", 60, 1.6, 2.2));
  h_d0pt = unique_ptr<TH1D>(new TH1D("h_d0pt","",100, 0 ,15));
  h_dstarpt = unique_ptr<TH1D>(new TH1D("h_dstarpt","",100, 0 ,15));
  h_PS3pt = unique_ptr<TH1D>(new TH1D("h_PS3pt","",100, 0 ,2));
  h_K1pt = unique_ptr<TH1D>(new TH1D("h_K1pt","",100, 0 ,10));
  h_P2pt = unique_ptr<TH1D>(new TH1D("h_P2pt","",100, 0 ,10));
  h_PS3eta = unique_ptr<TH1D>(new TH1D("h_PS3eta","",50, -3 ,3));
  h_P2eta = unique_ptr<TH1D>(new TH1D("h_P2eta","",50, -3 ,3));
  h_K1eta = unique_ptr<TH1D>(new TH1D("h_K1eta","",50, -3 ,3));
  h_D0masscut = unique_ptr<TH1D>(new TH1D("h_D0masscut", "", 60, 1.6, 2.2));
  h_deltaMassD0Dstar = unique_ptr<TH1D>(new TH1D("h_deltaMassD0Dstar", "", 64, 0.138, 0.17));
  h_D0masscut_rightDLcut = unique_ptr<TH1D>(new TH1D("h_D0masscut_rightDLcut", "", 60, 1.6, 2.2));
  h_deltaMassD0Dstar_rightDLcut = unique_ptr<TH1D>(new TH1D("h_deltaMassD0Dstar_rightDLcut", "", 64, 0.138, 0.17));

  */
  
  file = new TFile(outFile.c_str(), "recreate"); // check
  t_event = new TTree("Events", "Events");
  //t_event->SetAutoSave(-500000000);
  t_event->SetAutoSave(0);

#if ROOT_VERSION_CODE > ROOT_VERSION(6, 6, 0)
  t_event->SetImplicitMT(false);
#endif

  h_trackpt = new TH1D("h_trackpt","",500, 0.,50.);
  h_trackptlow = new TH1D("h_trackptlow","",100, 0.,1.);
  h_tracketa = new TH1D("h_tracketa","",160, -4.,4.);

  h_d0pt = new TH1D("h_d0pt","",100, 0. ,15.);
  h_dstarpt = new TH1D("h_dstarpt","",100, 0. ,15.);
  h_PS3pt = new TH1D("h_PS3pt","",100, 0. ,2.);
  h_K1pt = new TH1D("h_K1pt","",100, 0. ,10.);
  h_P2pt = new TH1D("h_P2pt","",100, 0. ,10.);
  h_PS3eta = new TH1D("h_PS3eta","",50, -3. ,3.);
  h_P2eta = new TH1D("h_P2eta","",50, -3. ,3.);
  h_K1eta = new TH1D("h_K1eta","",50, -3. ,3.);
  h_D0masscut = new TH1D("h_D0masscut", "", 60, 1.6, 2.2);
  h_deltaMassD0Dstar = new TH1D("h_deltaMassD0Dstar", "", 64, 0.138, 0.17);
  h_D0masscut_rightDLcut = new TH1D("h_D0masscut_rightDLcut", "", 60, 1.6, 2.2);
  h_deltaMassD0Dstar_rightDLcut = new TH1D("h_deltaMassD0Dstar_rightDLcut", "", 64, 0.138, 0.17);

  // H4lepton
  h_p_e = new TH1D("e_momentum", "Electron momentum", 200, 0., 200.);
  h_et_e = new TH1D("e_eT", "Electron eT", 200, 0., 200.);
  h_pt_e_b4 = new TH1D("b4_e_pT", "Electron pT", 200, 0, 200.);
  h_eta_e_b4 = new TH1D("b4_e_eta", "Electron eta", 140, -3.5, 3.5);
  h_phi_e = new TH1D("e_phi", "Electron phi", 314, -3.17, 3.17);
  h_sc_eta = new TH1D("e_SC_eta", "Electron SC eta", 140, -3.5, 3.5);
  h_sc_rawE = new TH1D("e_SC_rawE", "Electron SC rawE", 200, 0., 200.);
  h_relPFIso_e = new TH1D("e_RelPFIso", "R.PFIso", 100, 0., 5.);
  double Iso[12] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.5, 10.};
  h_relPFIso_pt_e = new TH2D("e_RelPFIso_pT","R.PFIso 2D",11,Iso,100,0., 50.);
  h_dxy_e = new TH1D("e_dxy", "Electron dxy", 100, 0., 1.);
  h_SIP3d_e_b4 = new TH1D("SIP3d_e", "SIP_3D for Electron", 100, 0., 10.);
  h_misshite = new TH1D("e_misshit", "e track missing hits", 5, 0., 5.);
  
  // h_Vertex_Multiplicity_D0mass_allSpectrum = unique_ptr<TH1D>(new TH1D("h_Vertex_Multiplicity_D0mass_allSpectrum","", 10, 0, 10));
  // h_Ncand_D0mass_allSpectrum = unique_ptr<TH1D>(new TH1D("h_Ncand_D0mass_allSpectrum","", 10, 0, 10));  
  // t_Mu = new TTree("Events", "Events"); // same conversion as official nAOD  
  // t_Dmeson = new TTree("Dmeson", "Dmeson");

  //---------------------------- Gen Particle reserve --------------------------//
  
  GenPart_pt.reserve(nReserve_GenPart);
  GenPart_eta.reserve(nReserve_GenPart);
  GenPart_phi.reserve(nReserve_GenPart);
  GenPart_mass.reserve(nReserve_GenPart);
  GenPart_pdgId.reserve(nReserve_GenPart);
  GenPart_status.reserve(nReserve_GenPart);
  GenPart_statusFlags.reserve(nReserve_GenPart);
  GenPart_genPartIdxMother.reserve(nReserve_GenPart);

  GenPart_Id.reserve(nReserve_GenPart);
  GenPart_isNano.reserve(nReserve_GenPart);
  GenPart_parpdgId.reserve(nReserve_GenPart);
  GenPart_sparpdgId.reserve(nReserve_GenPart);
  GenPart_numberOfDaughters.reserve(nReserve_GenPart);
  GenPart_nstchgdaug.reserve(nReserve_GenPart);
  GenPart_vx.reserve(nReserve_GenPart);
  GenPart_vy.reserve(nReserve_GenPart);
  GenPart_vz.reserve(nReserve_GenPart);
  GenPart_mvx.reserve(nReserve_GenPart);
  GenPart_mvy.reserve(nReserve_GenPart);
  GenPart_mvz.reserve(nReserve_GenPart);
  GenPart_recIdx.reserve(nReserve_GenPart);

  //--------------------------------- IsoTrack reserve -----------------------------//
  
  IsoTrack_dxy.reserve(nReserve_IsoTrack);
  IsoTrack_dz.reserve(nReserve_IsoTrack);
  IsoTrack_ets.reserve(nReserve_IsoTrack);
  IsoTrack_isHighPurityTrack.reserve(nReserve_IsoTrack);
  IsoTrack_isPFcand.reserve(nReserve_IsoTrack);
  IsoTrack_miniPFreliso_all.reserve(nReserve_IsoTrack);
  IsoTrack_miniPFreliso_chg.reserve(nReserve_IsoTrack);
  IsoTrack_pdgId.reserve(nReserve_IsoTrack);
  IsoTrack_PFreliso03_all.reserve(nReserve_IsoTrack);
  IsoTrack_PFreliso03_chg.reserve(nReserve_IsoTrack);
  IsoTrack_phi.reserve(nReserve_IsoTrack);
  IsoTrack_pt.reserve(nReserve_IsoTrack);

  //--------------------------------- OtherPV reserve -----------------------------//
  
  OtherPV_z.reserve(nReserve_OtherPV);

  //--------------------------------- PVtx reserve -----------------------------//
  
  PVtx_Id.reserve(nReserve_PVtx);
  PVtx_isMain.reserve(nReserve_PVtx);
  PVtx_isMainSim.reserve(nReserve_PVtx);
  PVtx_isGood.reserve(nReserve_PVtx);
  PVtx_isValid.reserve(nReserve_PVtx);
  PVtx_isFake.reserve(nReserve_PVtx);
  PVtx_isTrigUnique.reserve(nReserve_PVtx);
  PVtx_isUnbiased.reserve(nReserve_PVtx);
  PVtx_ntrk.reserve(nReserve_PVtx);
  PVtx_ntrkfit.reserve(nReserve_PVtx);
  PVtx_chi2.reserve(nReserve_PVtx);
  PVtx_ndof.reserve(nReserve_PVtx);
  PVtx_score.reserve(nReserve_PVtx);
  PVtx_sumPt.reserve(nReserve_PVtx);
  PVtx_Rho.reserve(nReserve_PVtx);
  PVtx_x.reserve(nReserve_PVtx);
  PVtx_y.reserve(nReserve_PVtx);
  PVtx_z.reserve(nReserve_PVtx);
  PVtx_Covxx.reserve(nReserve_PVtx);
  PVtx_Covyx.reserve(nReserve_PVtx);
  PVtx_Covzx.reserve(nReserve_PVtx);
  PVtx_Covyy.reserve(nReserve_PVtx);
  PVtx_Covzy.reserve(nReserve_PVtx);
  PVtx_Covzz.reserve(nReserve_PVtx);

  //--------------------------------- Muon reserve -----------------------------//
  
  Muon_charge.reserve(nReserve_Muon);
  Muon_tightCharge.reserve(nReserve_Muon); 
  Muon_pt.reserve(nReserve_Muon);
  Muon_ptErr.reserve(nReserve_Muon);
  Muon_eta.reserve(nReserve_Muon);
  Muon_phi.reserve(nReserve_Muon);  
  Muon_mass.reserve(nReserve_Muon);
  Muon_dxy.reserve(nReserve_Muon);
  Muon_dxyBest.reserve(nReserve_Muon);
  Muon_dxyErr.reserve(nReserve_Muon);
  Muon_dz.reserve(nReserve_Muon);
  Muon_dzBest.reserve(nReserve_Muon);
  Muon_dzErr.reserve(nReserve_Muon);
  Muon_ip3d.reserve(nReserve_Muon);
  Muon_sip3d.reserve(nReserve_Muon);
  Muon_ip3dBest.reserve(nReserve_Muon);
  Muon_sip3dBest.reserve(nReserve_Muon);
  Muon_pfRelIso03_all.reserve(nReserve_Muon);
  Muon_pfRelIso03_chg.reserve(nReserve_Muon);
  Muon_pfRelIso04_all.reserve(nReserve_Muon);
  Muon_miniPFRelIso_all.reserve(nReserve_Muon); 
  Muon_miniPFRelIso_chg.reserve(nReserve_Muon); 
  Muon_jetIdx.reserve(nReserve_Muon); 
  Muon_isGlobal.reserve(nReserve_Muon);
  Muon_isTracker.reserve(nReserve_Muon);
  Muon_isPFcand.reserve(nReserve_Muon);
  Muon_softId.reserve(nReserve_Muon);
  Muon_mediumId.reserve(nReserve_Muon); 
  Muon_tightId.reserve(nReserve_Muon);
  Muon_highPtId.reserve(nReserve_Muon); 
  Muon_nStations.reserve(nReserve_Muon);
  Muon_nTrackerLayers.reserve(nReserve_Muon); 
  Muon_segmentComp.reserve(nReserve_Muon); 
  Muon_cleanmask.reserve(nReserve_Muon); 
  Muon_mvaTTH.reserve(nReserve_Muon); 
  Muon_pdgId.reserve(nReserve_Muon); 
  Muon_genPartFlav.reserve(nReserve_Muon); 
  Muon_genPartIdx.reserve(nReserve_Muon); 

  Muon_Id.reserve(nReserve_Muon); 
  Muon_x.reserve(nReserve_Muon); 
  Muon_y.reserve(nReserve_Muon); 
  Muon_z.reserve(nReserve_Muon); 
  Muon_gpt.reserve(nReserve_Muon);
  Muon_geta.reserve(nReserve_Muon);
  Muon_gphi.reserve(nReserve_Muon);  
  Muon_looseId.reserve(nReserve_Muon); 
  Muon_softId4.reserve(nReserve_Muon);
  Muon_softIdBest.reserve(nReserve_Muon);
  Muon_isNano.reserve(nReserve_Muon);
  Muon_isMini.reserve(nReserve_Muon);
  Muon_isGood.reserve(nReserve_Muon);
  Muon_isGoodLast.reserve(nReserve_Muon);
  Muon_isGoodAng.reserve(nReserve_Muon);
  Muon_isArbitrated.reserve(nReserve_Muon);
  Muon_isStandAlone.reserve(nReserve_Muon);
  Muon_isRPCcand.reserve(nReserve_Muon);
  Muon_nValid.reserve(nReserve_Muon); 
  Muon_nPix.reserve(nReserve_Muon); 
  Muon_Chi2.reserve(nReserve_Muon); 
  Muon_gnValid.reserve(nReserve_Muon); 
  Muon_gnPix.reserve(nReserve_Muon); 
  Muon_gChi2.reserve(nReserve_Muon); 
  Muon_gnValidMu.reserve(nReserve_Muon); 
  Muon_vtxIdx.reserve(nReserve_Muon); 
  Muon_vtxFlag.reserve(nReserve_Muon); 
  Muon_trkIdx.reserve(nReserve_Muon); 
  Muon_simIdx.reserve(nReserve_Muon); 

  //------------------------------- Dimuon reserve --------------------------//

  Dimut1_muIdx.reserve(nReserve_Dimu);
  Dimut1_dxy.reserve(nReserve_Dimu);
  Dimut1_dz.reserve(nReserve_Dimu);
  Dimut2_muIdx.reserve(nReserve_Dimu);
  Dimut2_dxy.reserve(nReserve_Dimu);
  Dimut2_dz.reserve(nReserve_Dimu);
  Dimu_pt.reserve(nReserve_Dimu);
  Dimu_eta.reserve(nReserve_Dimu);
  Dimu_phi.reserve(nReserve_Dimu);
  Dimu_rap.reserve(nReserve_Dimu);
  Dimu_mass.reserve(nReserve_Dimu);
  Dimu_charge.reserve(nReserve_Dimu);
  Dimu_simIdx.reserve(nReserve_Dimu);
  Dimu_vtxIdx.reserve(nReserve_Dimu);
  Dimu_chi2.reserve(nReserve_Dimu);
  Dimu_dlxy.reserve(nReserve_Dimu);
  Dimu_dlxyErr.reserve(nReserve_Dimu);
  Dimu_dlxySig.reserve(nReserve_Dimu);
  Dimu_cosphixy.reserve(nReserve_Dimu);
  Dimu_dl.reserve(nReserve_Dimu);
  Dimu_dlErr.reserve(nReserve_Dimu);
  Dimu_dlSig.reserve(nReserve_Dimu);
  Dimu_cosphi.reserve(nReserve_Dimu);
  Dimu_ptfrac.reserve(nReserve_Dimu);
  Dimu_x.reserve(nReserve_Dimu);
  Dimu_y.reserve(nReserve_Dimu);
  Dimu_z.reserve(nReserve_Dimu);
  Dimu_Covxx.reserve(nReserve_Dimu);
  Dimu_Covyx.reserve(nReserve_Dimu);
  Dimu_Covzx.reserve(nReserve_Dimu);
  Dimu_Covyy.reserve(nReserve_Dimu);
  Dimu_Covzy.reserve(nReserve_Dimu);
  Dimu_Covzz.reserve(nReserve_Dimu);


  //------------------------------- Dmeson reserve --------------------------//
  
  D0t1_pt.reserve(nReserve_D0);
  D0t1_eta.reserve(nReserve_D0);
  D0t1_phi.reserve(nReserve_D0);
  D0t1_chg.reserve(nReserve_D0);
  D0t1_tkIdx.reserve(nReserve_D0);
  D0t1_Kprob.reserve(nReserve_D0);
  D0t1_piprob.reserve(nReserve_D0);
  D0t1_dEdxnmeas.reserve(nReserve_D0);
  D0t1_dEdxnsat.reserve(nReserve_D0);
  D0t1_vtxIdx.reserve(nReserve_D0);
  D0t1_chindof.reserve(nReserve_D0);
  D0t1_nValid.reserve(nReserve_D0);
  D0t1_nPix.reserve(nReserve_D0);
  D0t1_isHighPurity.reserve(nReserve_D0);
  D0t1_pdgId.reserve(nReserve_D0);
  D0t1_dxy.reserve(nReserve_D0);
  D0t1_dz.reserve(nReserve_D0);
  D0t2_pt.reserve(nReserve_D0);
  D0t2_eta.reserve(nReserve_D0);
  D0t2_phi.reserve(nReserve_D0);
  D0t2_chg.reserve(nReserve_D0);
  D0t2_tkIdx.reserve(nReserve_D0);
  D0t2_Kprob.reserve(nReserve_D0);
  D0t2_piprob.reserve(nReserve_D0);
  D0t2_dEdxnmeas.reserve(nReserve_D0);
  D0t2_dEdxnsat.reserve(nReserve_D0);
  D0t2_vtxIdx.reserve(nReserve_D0);
  D0t2_chindof.reserve(nReserve_D0);
  D0t2_nValid.reserve(nReserve_D0);
  D0t2_nPix.reserve(nReserve_D0);
  D0t2_isHighPurity.reserve(nReserve_D0);
  D0t2_pdgId.reserve(nReserve_D0);
  D0t2_dxy.reserve(nReserve_D0);
  D0t2_dz.reserve(nReserve_D0);
  D0t2_pdgId.reserve(nReserve_D0);
  D0_pt.reserve(nReserve_D0);
  D0_eta.reserve(nReserve_D0);
  D0_phi.reserve(nReserve_D0);
  D0_rap.reserve(nReserve_D0);
  D0_mass12.reserve(nReserve_D0);
  D0_mass21.reserve(nReserve_D0);
  D0_massKK.reserve(nReserve_D0);
  D0_simIdx.reserve(nReserve_D0);
  D0_DstarIdx.reserve(nReserve_D0);
  D0_ambiPrim.reserve(nReserve_D0);
  D0_vtxIdx.reserve(nReserve_D0);
  D0_hasMuon.reserve(nReserve_D0);
  D0_chi2.reserve(nReserve_D0);
  D0_dlxy.reserve(nReserve_D0);
  D0_dlxyErr.reserve(nReserve_D0);
  D0_dlxySig.reserve(nReserve_D0);
  D0_cosphixy.reserve(nReserve_D0);
  D0_dl.reserve(nReserve_D0);
  D0_dlErr.reserve(nReserve_D0);
  D0_dlSig.reserve(nReserve_D0);
  D0_cosphi.reserve(nReserve_D0);
  D0_ptfrac.reserve(nReserve_D0);
  D0_ptfrac15.reserve(nReserve_D0);
  D0_ptfrac10.reserve(nReserve_D0);
  D0_ptfrac07.reserve(nReserve_D0);
  D0_ptfrac04.reserve(nReserve_D0);
  D0_x.reserve(nReserve_D0);
  D0_y.reserve(nReserve_D0);
  D0_z.reserve(nReserve_D0);
  D0_Covxx.reserve(nReserve_D0);
  D0_Covyx.reserve(nReserve_D0);
  D0_Covzx.reserve(nReserve_D0);
  D0_Covyy.reserve(nReserve_D0);
  D0_Covzy.reserve(nReserve_D0);
  D0_Covzz.reserve(nReserve_D0);
  
  Dstarpis_pt.reserve(nReserve_Dstar);
  Dstarpis_eta.reserve(nReserve_Dstar);
  Dstarpis_phi.reserve(nReserve_Dstar);
  Dstarpis_ptr.reserve(nReserve_Dstar);
  Dstarpis_etar.reserve(nReserve_Dstar);
  Dstarpis_phir.reserve(nReserve_Dstar);
  Dstarpis_chg.reserve(nReserve_Dstar);
  Dstarpis_tkIdx.reserve(nReserve_Dstar);
  Dstarpis_Kprob.reserve(nReserve_Dstar);
  Dstarpis_piprob.reserve(nReserve_Dstar);
  Dstarpis_dEdxnmeas.reserve(nReserve_Dstar);
  Dstarpis_dEdxnsat.reserve(nReserve_Dstar);
  Dstarpis_vtxIdx.reserve(nReserve_Dstar);
  Dstarpis_chindof.reserve(nReserve_Dstar);
  Dstarpis_chir.reserve(nReserve_Dstar);
  Dstarpis_nValid.reserve(nReserve_Dstar);
  Dstarpis_nPix.reserve(nReserve_Dstar);
  Dstarpis_dxy.reserve(nReserve_Dstar);
  Dstarpis_dz.reserve(nReserve_Dstar);

  DstarD0_pt.reserve(nReserve_Dstar);
  DstarD0_eta.reserve(nReserve_Dstar);
  DstarD0_phi.reserve(nReserve_Dstar);
  DstarD0_mass.reserve(nReserve_Dstar);
  DstarD0_chi2.reserve(nReserve_Dstar);
  DstarD0_dlxy.reserve(nReserve_Dstar);
  DstarD0_dlxyErr.reserve(nReserve_Dstar);
  DstarD0_dlxySig.reserve(nReserve_Dstar);
  DstarD0_cosphixy.reserve(nReserve_Dstar);
  DstarD0_dl.reserve(nReserve_Dstar);
  DstarD0_dlErr.reserve(nReserve_Dstar);
  DstarD0_dlSig.reserve(nReserve_Dstar);
  DstarD0_cosphi.reserve(nReserve_Dstar);
  DstarD0_ptfrac.reserve(nReserve_Dstar); 
  DstarD0_ptfrac15.reserve(nReserve_Dstar); 
  DstarD0_ptfrac10.reserve(nReserve_Dstar); 
  DstarD0_ptfrac07.reserve(nReserve_Dstar); 
  DstarD0_ptfrac04.reserve(nReserve_Dstar); 
  DstarD0_x.reserve(nReserve_Dstar); 
  DstarD0_y.reserve(nReserve_Dstar); 
  DstarD0_z.reserve(nReserve_Dstar); 
  DstarD0_simIdx.reserve(nReserve_Dstar); 
  DstarD0_recIdx.reserve(nReserve_Dstar); 
  DstarD0_ambiPrim.reserve(nReserve_Dstar);

  DstarK_pt.reserve(nReserve_Dstar);
  DstarK_eta.reserve(nReserve_Dstar);
  DstarK_phi.reserve(nReserve_Dstar);
  DstarK_chg.reserve(nReserve_Dstar);
  DstarK_tkIdx.reserve(nReserve_Dstar);
  DstarK_Kprob.reserve(nReserve_Dstar);
  DstarK_piprob.reserve(nReserve_Dstar);
  DstarK_dEdxnmeas.reserve(nReserve_Dstar);
  DstarK_dEdxnsat.reserve(nReserve_Dstar);
  DstarK_vtxIdx.reserve(nReserve_Dstar);
  DstarK_chindof.reserve(nReserve_Dstar);
  DstarK_nValid.reserve(nReserve_Dstar);
  DstarK_nPix.reserve(nReserve_Dstar);
  DstarK_isHighPurity.reserve(nReserve_Dstar);
  DstarK_dxy.reserve(nReserve_Dstar);
  DstarK_dz.reserve(nReserve_Dstar);

  Dstarpi_pt.reserve(nReserve_Dstar);
  Dstarpi_eta.reserve(nReserve_Dstar);
  Dstarpi_phi.reserve(nReserve_Dstar);
  Dstarpi_chg.reserve(nReserve_Dstar);
  Dstarpi_tkIdx.reserve(nReserve_Dstar);
  Dstarpi_Kprob.reserve(nReserve_Dstar);
  Dstarpi_piprob.reserve(nReserve_Dstar);  
  Dstarpi_dEdxnmeas.reserve(nReserve_Dstar);
  Dstarpi_dEdxnsat.reserve(nReserve_Dstar);
  Dstarpi_vtxIdx.reserve(nReserve_Dstar);
  Dstarpi_chindof.reserve(nReserve_Dstar);
  Dstarpi_nValid.reserve(nReserve_Dstar);
  Dstarpi_nPix.reserve(nReserve_Dstar);
  Dstarpi_isHighPurity.reserve(nReserve_Dstar);
  Dstarpi_dxy.reserve(nReserve_Dstar);
  Dstarpi_dz.reserve(nReserve_Dstar);
  
  Dstar_pt.reserve(nReserve_Dstar);
  Dstar_eta.reserve(nReserve_Dstar);
  Dstar_phi.reserve(nReserve_Dstar);
  Dstar_rap.reserve(nReserve_Dstar);
  Dstar_deltam.reserve(nReserve_Dstar);
  Dstar_deltamr.reserve(nReserve_Dstar);
  Dstar_simIdx.reserve(nReserve_Dstar);
  Dstar_vtxIdx.reserve(nReserve_Dstar);
  Dstar_hasMuon.reserve(nReserve_Dstar);
  Dstar_ptfrac.reserve(nReserve_Dstar);

  TrigObj_id.reserve(nReserve_TrigObj);
  TrigObj_filterBits.reserve(nReserve_TrigObj);
  TrigObj_pt.reserve(nReserve_TrigObj);
  TrigObj_eta.reserve(nReserve_TrigObj);
  TrigObj_phi.reserve(nReserve_TrigObj);

  //---------------------------- Create branch for tree ------------------------//
  
  // nanoAOD run/event structure 
  t_event->Branch("run", &run, "run/I");
  t_event->Branch("event", &event, "event/l");
  t_event->Branch("luminosityBlock", &luminosityBlock, "luminosityBlock/I");

  t_event->Branch("CMSSW", &CMSSW, "CMSSW/I");

    // nanoAOD extension
    // store JSON info
    t_event->Branch("GoodLumisection", &GoodLumisection, "GoodLumisection/O");
    // store dataset info
    t_event->Branch("MCdataset", &MCdataset, "MCdataset/O");
    t_event->Branch("ZeroBiasdataset", &ZeroBiasdataset, "ZeroBiasdataset/O");
    t_event->Branch("MinimumBiasdataset", &MinimumBiasdataset, "MinimumBiasdataset/O");
    t_event->Branch("Jetdataset", &Jetdataset, "Jetdataset/O");
    t_event->Branch("MultiJetdataset", &MultiJetdataset, "MultiJetdataset/O");
    t_event->Branch("Mudataset", &Mudataset, "Mudataset/O");
    t_event->Branch("MuMonitordataset", &MuMonitordataset, "MuMonitordataset/O");
    t_event->Branch("DoubleMudataset", &DoubleMudataset, "DoubleMudataset/O");
    t_event->Branch("MuHaddataset", &MuHaddataset, "MuHaddataset/O");
    t_event->Branch("MuOniadataset", &MuOniadataset, "MuOniadataset/O");
    t_event->Branch("Charmoniumdataset", &Charmoniumdataset, "Charmoniumdataset/O");
    t_event->Branch("BParkingdataset", &BParkingdataset, "BParkingdataset/O");
    t_event->Branch("BTaudataset", &BTaudataset, "BTaudataset/O");
    t_event->Branch("Electrondataset", &Electrondataset, "Electrondataset/O");
    t_event->Branch("DoubleElectrondataset", &DoubleElectrondataset, "DoubleElectrondataset/O");
    t_event->Branch("Photondataset", &Photondataset, "Photondataset/O");
    t_event->Branch("EGMonitordataset", &EGMonitordataset, "EGMonitordataset/O");
    t_event->Branch("MuEGdataset", &MuEGdataset, "MuEGdataset/O");
    t_event->Branch("Commissioningdataset", &Commissioningdataset, "Commissioningdataset/O");
    // store trigger bits //
    t_event->Branch("GoodMinBiasTrigger", &GoodMinBiasTrigger, "GoodMinBiasTrigger/O");
    t_event->Branch("GoodJetTrigger", &GoodJetTrigger, "GoodJetTrigger/O");
    t_event->Branch("GoodMuTrigger", &GoodMuTrigger, "GoodMuTrigger/O");
    t_event->Branch("GoodETrigger", &GoodETrigger, "GoodETrigger/O");
    // *** check whether the "Trig" and "Trigger" variables are the same 
    //     (they should be)
    //t_event->Branch("GoodMinimumBiasTrig", &GoodMinimumBiasTrig, "GoodMinimumBiasTrig/O");
    //t_event->Branch("GoodJetTrig", &GoodJetTrig, "GoodJetTrig/O");
    //t_event->Branch("GoodMuTrig", &GoodMuTrig, "GoodMuTrig/O");
    //t_event->Branch("GoodETrig", &GoodETrig, "GoodETrig/O");
    // generic bits
    t_event->Branch("ZeroBiasTrig", &ZeroBiasTrig, "ZeroBiasTrig/O");
    t_event->Branch("MinimumBiasTrig", &MinimumBiasTrig, "MinimumBiasTrig/O");
    t_event->Branch("JetTrig", &JetTrig, "JetTrig/O");
    t_event->Branch("MultiJetTrig", &MultiJetTrig, "MultiJetTrig/O");
    t_event->Branch("JetMETTauMonitorTrig", &JetMETTauMonitorTrig, "JetMETTauMonitorTrig/O");
    t_event->Branch("MuTrig", &MuTrig, "MuTrig/O");
    t_event->Branch("MuHadTrig", &MuHadTrig, "MuHadTrig/O");
    t_event->Branch("DoubleMuTrig", &DoubleMuTrig, "DoubleMuTrig/O");
    t_event->Branch("MuEGTrig", &MuEGTrig, "MuEGTrig/O");
    t_event->Branch("ElectronTrig", &ElectronTrig, "ElectronTrig/O");
    t_event->Branch("DoubleElectronTrig", &DoubleElectronTrig, "DoubleElectronTrig/O");
    t_event->Branch("PhotonTrig", &PhotonTrig, "PhotonTrig/O");
    t_event->Branch("MuMonitorTrig", &MuMonitorTrig, "MuMonitorTrig/O");
    t_event->Branch("EGMonitorTrig", &EGMonitorTrig, "EGMonitorTrig/O");
    t_event->Branch("MuOniaTrig", &MuOniaTrig, "MuOniaTrig/O");
    t_event->Branch("CharmoniumTrig", &CharmoniumTrig, "CharmoniumTrig/O");
    t_event->Branch("BTauTrig", &BTauTrig, "BTauTrig/O");
    t_event->Branch("BParkingTrig", &BParkingTrig, "BParkingTrig/O");
    t_event->Branch("METFwdTrig", &METFwdTrig, "METFwdTrig/O");
    t_event->Branch("CommissioningTrig", &CommissioningTrig, "CommissioningTrig/O");

    t_event->Branch("ZeroBiasFlag", &ZeroBiasFlag, "ZeroBiasFlag/I");
    t_event->Branch("MinBiasFlag", &MinBiasFlag, "MinBiasFlag/I");
    t_event->Branch("MinBiasMult", &MinBiasMult, "MinBiasMult/I");
    t_event->Branch("MuThresh", &MuThresh, "MuThresh/I");
    t_event->Branch("MuL1Thresh", &MuL1Thresh, "MuL1Thresh/I");
    t_event->Branch("MuL2Thresh", &MuL2Thresh, "MuL2Thresh/I");
    t_event->Branch("IsoMuThresh", &IsoMuThresh, "IsoMuThresh/I");
    t_event->Branch("DoubleMuThresh", &DoubleMuThresh, "DoubleMuThresh/I");
    t_event->Branch("JpsiThresh", &JpsiThresh, "JpsiThresh/I");
    t_event->Branch("MuHadFlag", &MuHadFlag, "MuHadFlag/I");
    t_event->Branch("MuEGFlag", &MuEGFlag, "MuEGFlag/I");
    t_event->Branch("ElectronThresh", &ElectronThresh, "ElectronThresh/I");
    t_event->Branch("DoubleElectronThresh", &DoubleElectronThresh, "DoubleElectronThresh/I");
    t_event->Branch("PhotonThresh", &PhotonThresh, "PhotonThresh/I");
    t_event->Branch("JetThresh", &JetThresh, "JetThresh/I");
    t_event->Branch("DiJetThresh", &DiJetThresh, "DiJetThresh/I");
    t_event->Branch("TriJetThresh", &TriJetThresh, "TriJetThresh/I");
    t_event->Branch("QuadJetThresh", &QuadJetThresh, "QuadJetThresh/I");
    t_event->Branch("HTThresh", &HTThresh, "HTThresh/I");
    t_event->Branch("BThresh", &BThresh, "BThresh/I");
    t_event->Branch("METThresh", &METThresh, "METThresh/I");

  //}

  if (!isData) {

    cout << "This is MC" << endl;
    //---------------------- Create branch of GenPart's tree -------------------//
    
    // official nanoAOD structure
    t_event->Branch("nGenPart", &nGenPart, "nGenPart/I");
    t_event->Branch("GenPart_pt", GenPart_pt.data(), "GenPart_pt[nGenPart]/F");
    t_event->Branch("GenPart_eta", GenPart_eta.data(), "GenPart_eta[nGenPart]/F");
    t_event->Branch("GenPart_phi", GenPart_phi.data(), "GenPart_phi[nGenPart]/F");
    t_event->Branch("GenPart_mass", GenPart_mass.data(), "GenPart_mass[nGenPart]/F");
    t_event->Branch("GenPart_pdgId", GenPart_pdgId.data(), "GenPart_pdgId[nGenPart]/I");
    t_event->Branch("GenPart_status", GenPart_status.data(), "GenPart_status[nGenPart]/I");
    t_event->Branch("GenPart_statusFlags", GenPart_statusFlags.data(), "GenPart_statusFlags[nGenPart]/I");
    t_event->Branch("GenPart_genPartIdxMother", GenPart_genPartIdxMother.data(), "GenPart_genPartIdxMother[nGenPart]/I");

    if (nanoext) {
    // GenPart extension
      t_event->Branch("GenPart_Id", GenPart_Id.data(), "GenPart_Id[nGenPart]/I");
      t_event->Branch("GenPart_isNano", GenPart_isNano.data(), "GenPart_isNano[nGenPart]/O");
      t_event->Branch("GenPart_parpdgId", GenPart_parpdgId.data(), "GenPart_parpdgId[nGenPart]/I");
      t_event->Branch("GenPart_sparpdgId", GenPart_sparpdgId.data(), "GenPart_sparpdgId[nGenPart]/I");
      t_event->Branch("GenPart_numberOfDaughters", GenPart_numberOfDaughters.data(), "GenPart_numberOfDaughters[nGenPart]/I");
      t_event->Branch("GenPart_nstchgdaug", GenPart_nstchgdaug.data(), "GenPart_nstchgdaug[nGenPart]/I");
      t_event->Branch("GenPart_vx", GenPart_vx.data(), "GenPart_vx[nGenPart]/F");
      t_event->Branch("GenPart_vy", GenPart_vy.data(), "GenPart_vy[nGenPart]/F");
      t_event->Branch("GenPart_vz", GenPart_vz.data(), "GenPart_vz[nGenPart]/F");
      t_event->Branch("GenPart_mvx", GenPart_mvx.data(), "GenPart_mvx[nGenPart]/F");
      t_event->Branch("GenPart_mvy", GenPart_mvy.data(), "GenPart_mvy[nGenPart]/F");
      t_event->Branch("GenPart_mvz", GenPart_mvz.data(), "GenPart_mvz[nGenPart]/F");
      t_event->Branch("GenPart_recIdx", GenPart_recIdx.data(), "GenPart_recIdx[nGenPart]/I");

      // GenPV extension
      t_event->Branch("GenPV_x", &GenPV_x, "GenPV_x/F");
      t_event->Branch("GenPV_y", &GenPV_y, "GenPV_y/F");
      t_event->Branch("GenPV_z", &GenPV_z, "GenPV_z/F");
      t_event->Branch("GenPV_recIdx", &GenPV_recIdx, "GenPV_recIdx/I");
      t_event->Branch("GenPV_chmult", &GenPV_chmult, "GenPV_chmult/I");
    } // end of nanoext
    
  } // end of event !isData


  // nanoAOD track structure
  t_event->Branch("nTrk", &nTrk, "nTrk/I");

  // nanoAOD IsoTrack structure
  t_event->Branch("nIsoTrack", &nIsoTrack, "nIsoTrack/I");

  // nanoAOD PV structure
  t_event->Branch("PV_npvs", &PV_npvs, "PV_npvs/I");
  t_event->Branch("PV_npvsGood", &PV_npvsGood, "PV_npvsGood/I");
  t_event->Branch("PV_chi2", &PV_chi2, "PV_chi2/F");
  t_event->Branch("PV_ndof", &PV_ndof, "PV_ndof/F");
  t_event->Branch("PV_score", &PV_score, "PV_score/F");
  t_event->Branch("PV_x", &PV_x, "PV_x/F");
  t_event->Branch("PV_y", &PV_y, "PV_y/F");
  t_event->Branch("PV_z", &PV_z, "PV_z/F");

  // nanoAOD OtherPV structure 
  t_event->Branch("nOtherPV", &nOtherPV, "nOtherPV/I");
  t_event->Branch("OtherPV_z", OtherPV_z.data(), "OtherPV_z[nOtherPV]/F");

  if (nanoext) {

    // PVtx extension
    t_event->Branch("nPVtx", &nPVtx, "nPVtx/I");
    t_event->Branch("PVtx_Id", PVtx_Id.data(), "PVtx_Id[nPVtx]/I");
    t_event->Branch("PVtx_isMain", PVtx_isMain.data(), "PVtx_isMain[nPVtx]/O");
    t_event->Branch("PVtx_isMainSim", PVtx_isMainSim.data(), "PVtx_isMainSim[nPVtx]/O");
    t_event->Branch("PVtx_isGood", PVtx_isGood.data(), "PVtx_isGood[nPVtx]/O");
    t_event->Branch("PVtx_isValid", PVtx_isValid.data(), "PVtx_isValid[nPVtx]/O");
    t_event->Branch("PVtx_isFake", PVtx_isFake.data(), "PVtx_isFake[nPVtx]/O");
    t_event->Branch("PVtx_isTrigUnique", PVtx_isTrigUnique.data(), "PVtx_isTrigUnique[nPVtx]/O");
    t_event->Branch("PVtx_isUnbiased", PVtx_isUnbiased.data(), "PVtx_isUnbiased[nPVtx]/O");
    t_event->Branch("PVtx_ntrk", PVtx_ntrk.data(), "PVtx_ntrk[nPVtx]/I");
    t_event->Branch("PVtx_ntrkfit", PVtx_ntrkfit.data(), "PVtx_ntrkfit[nPVtx]/I");
    t_event->Branch("PVtx_chi2", PVtx_chi2.data(), "PVtx_chi2[nPVtx]/F");
    t_event->Branch("PVtx_ndof", PVtx_ndof.data(), "PVtx_ndof[nPVtx]/F");
    t_event->Branch("PVtx_score", PVtx_score.data(), "PVtx_score[nPVtx]/F");
    t_event->Branch("PVtx_sumPt", PVtx_sumPt.data(), "PVtx_sumPt[nPVtx]/F");
    t_event->Branch("PVtx_Rho", PVtx_Rho.data(), "PVtx_Rho[nPVtx]/F");
    t_event->Branch("PVtx_x", PVtx_x.data(), "PVtx_x[nPVtx]/F");
    t_event->Branch("PVtx_y", PVtx_y.data(), "PVtx_y[nPVtx]/F");
    t_event->Branch("PVtx_z", PVtx_z.data(), "PVtx_z[nPVtx]/F");
    if (covout) {
      t_event->Branch("PVtx_Covxx", PVtx_Covxx.data(), "PVtx_Covxx[nPVtx]/F");
      t_event->Branch("PVtx_Covyx", PVtx_Covyx.data(), "PVtx_Covyx[nPVtx]/F");
      t_event->Branch("PVtx_Covzx", PVtx_Covzx.data(), "PVtx_Covzx[nPVtx]/F");
      t_event->Branch("PVtx_Covyy", PVtx_Covyy.data(), "PVtx_Covyy[nPVtx]/F");
      t_event->Branch("PVtx_Covzy", PVtx_Covzy.data(), "PVtx_Covzy[nPVtx]/F");
      t_event->Branch("PVtx_Covzz", PVtx_Covzz.data(), "PVtx_Covzz[nPVtx]/F");
    }

    // beamspot extension
    t_event->Branch("Bsp_x",&Bsp_x, "Bsp_x/F");
    t_event->Branch("Bsp_y",&Bsp_y, "Bsp_y/F");
    t_event->Branch("Bsp_z",&Bsp_z, "Bsp_z/F");
    t_event->Branch("Bsp_sigmaz",&Bsp_sigmaz, "Bsp_sigmaz/F");
    t_event->Branch("Bsp_dxdz",&Bsp_dxdz, "Bsp_dxdz/F");
    t_event->Branch("Bsp_dydz",&Bsp_dydz, "Bsp_dydz/F");
    t_event->Branch("Bsp_widthx",&Bsp_widthx, "Bsp_widthx/F");
    t_event->Branch("Bsp_widthy",&Bsp_widthy, "Bsp_widthy/F");
  } // end of nanoext
  
  //-------------------------- Create branch of Muons's tree -------------------//
  // official nanoAOD
  t_event->Branch("nMuon", &nMuon, "nMuon/i");
  t_event->Branch("Muon_charge", Muon_charge.data(), "Muon_charge[nMuon]/I");
  t_event->Branch("Muon_tightCharge", Muon_tightCharge.data(), "Muon_tightCharge[nMuon]/I");
  t_event->Branch("Muon_pt", Muon_pt.data(), "Muon_pt[nMuon]/F");
  t_event->Branch("Muon_ptErr", Muon_ptErr.data(), "Muon_ptErr[nMuon]/F");
  t_event->Branch("Muon_eta", Muon_eta.data(), "Muon_eta[nMuon]/F"); // betul
  t_event->Branch("Muon_phi", Muon_phi.data(), "Muon_phi[nMuon]/F");
  t_event->Branch("Muon_mass", Muon_mass.data(), "Muon_mass[nMuon]/F");
  t_event->Branch("Muon_dxy", Muon_dxy.data(), "Muon_dxy[nMuon]/F");
  t_event->Branch("Muon_dxyBest", Muon_dxyBest.data(), "Muon_dxyBest[nMuon]/F");
  t_event->Branch("Muon_dxyErr", Muon_dxyErr.data(), "Muon_dxyErr[nMuon]/F");
  t_event->Branch("Muon_dz", Muon_dz.data(), "Muon_dz[nMuon]/F");
  t_event->Branch("Muon_dzBest", Muon_dzBest.data(), "Muon_dzBest[nMuon]/F");
  t_event->Branch("Muon_dzErr", Muon_dzErr.data(), "Muon_dzErr[nMuon]/F");
  t_event->Branch("Muon_ip3d", Muon_ip3d.data(), "Muon_ip3d[nMuon]/F");
  t_event->Branch("Muon_sip3d", Muon_sip3d.data(), "Muon_sip3d[nMuon]/F");
  t_event->Branch("Muon_ip3dBest", Muon_ip3dBest.data(), "Muon_ip3dBest[nMuon]/F");
  t_event->Branch("Muon_sip3dBest", Muon_sip3dBest.data(), "Muon_sip3dBest[nMuon]/F");
  t_event->Branch("Muon_pfRelIso03_all", Muon_pfRelIso03_all.data(), "Muon_pfRelIso03_all[nMuon]/F");
  t_event->Branch("Muon_pfRelIso03_chg", Muon_pfRelIso03_chg.data(), "Muon_pfRelIso03_chg[nMuon]/F");  
  t_event->Branch("Muon_pfRelIso04_all", Muon_pfRelIso04_all.data(), "Muon_pfRelIso04_all[nMuon]/F");
  t_event->Branch("Muon_miniPFRelIso_all", Muon_miniPFRelIso_all.data(), "Muon_miniPFRelIso_all[nMuon]/F");
  t_event->Branch("Muon_miniPFRelIso_chg", Muon_miniPFRelIso_chg.data(), "Muon_miniPFRelIso_chg[nMuon]/F");
  t_event->Branch("Muon_jetIdx", Muon_jetIdx.data(), "Muon_jetIdx[nMuon]/I");
  t_event->Branch("Muon_isPFcand", Muon_isPFcand.data(), "Muon_isPFcand[nMuon]/O");
  t_event->Branch("Muon_softId", Muon_softId.data(), "Muon_softId[nMuon]/O");
  t_event->Branch("Muon_mediumId", Muon_mediumId.data(), "Muon_mediumId[nMuon]/O");
  t_event->Branch("Muon_tightId", Muon_tightId.data(), "Muon_tightId[nMuon]/O");
  t_event->Branch("Muon_highPtId", Muon_highPtId.data(), "Muon_highPtId[nMuon]/b");
  t_event->Branch("Muon_nStations", Muon_nStations.data(), "Muon_nStations[nMuon]/I");
  t_event->Branch("Muon_nTrackerLayers", Muon_nTrackerLayers.data(), "Muon_nTrackerLayers[nMuon]/I");
  t_event->Branch("Muon_segmentComp", Muon_segmentComp.data(), "Muon_segmentComp[nMuon]/F");
  t_event->Branch("Muon_cleanmask", Muon_cleanmask.data(), "Muon_cleanmask[nMuon]/b");
  t_event->Branch("Muon_mvaTTH", Muon_mvaTTH.data(), "Muon_mvaTTH[nMuon]/F");
  t_event->Branch("Muon_pdgId", Muon_pdgId.data(), "Muon_pdgId[nMuon]/I");
  t_event->Branch("Muon_genPartFlav", Muon_genPartFlav.data(), "Muon_genPartFlav[nMuon]/b");
  t_event->Branch("Muon_genPartIdx", Muon_genPartIdx.data(), "Muon_genPartIdx[nMuon]/I");
  t_event->Branch("Muon_isGlobal", Muon_isGlobal.data(), "Muon_isGlobal[nMuon]/O");
  t_event->Branch("Muon_isTracker", Muon_isTracker.data(), "Muon_isTracker[nMuon]/O");

  if (nanoext) {
    // nanoAOD extension
    t_event->Branch("Muon_Id", Muon_Id.data(), "Muon_Id[nMuon]/I");
    t_event->Branch("Muon_x", Muon_x.data(), "Muon_x[nMuon]/F");
    t_event->Branch("Muon_y", Muon_y.data(), "Muon_y[nMuon]/F");
    t_event->Branch("Muon_z", Muon_z.data(), "Muon_z[nMuon]/F");
    t_event->Branch("Muon_gpt", Muon_gpt.data(), "Muon_gpt[nMuon]/F");
    t_event->Branch("Muon_geta", Muon_geta.data(), "Muon_geta[nMuon]/F");
    t_event->Branch("Muon_gphi", Muon_gphi.data(), "Muon_gphi[nMuon]/F");
    t_event->Branch("Muon_looseId", Muon_looseId.data(), "Muon_looseId[nMuon]/O");
    t_event->Branch("Muon_softId4", Muon_softId4.data(), "Muon_softId4[nMuon]/O");
    t_event->Branch("Muon_softIdBest", Muon_softIdBest.data(), "Muon_softIdBest[nMuon]/O");
    t_event->Branch("Muon_isNano", Muon_isNano.data(), "Muon_isNano[nMuon]/O");
    t_event->Branch("Muon_isMini", Muon_isMini.data(), "Muon_isMini[nMuon]/O");
    t_event->Branch("Muon_isGood", Muon_isGood.data(), "Muon_isGood[nMuon]/O");
    t_event->Branch("Muon_isGoodLast", Muon_isGoodLast.data(), "Muon_isGoodLast[nMuon]/O");
    t_event->Branch("Muon_isGoodAng", Muon_isGoodAng.data(), "Muon_isGoodAng[nMuon]/O");
    t_event->Branch("Muon_isArbitrated", Muon_isArbitrated.data(), "Muon_isArbitrated[nMuon]/O");
    t_event->Branch("Muon_isStandAlone", Muon_isStandAlone.data(), "Muon_isStandAlone[nMuon]/O");
    t_event->Branch("Muon_isRPCcand", Muon_isRPCcand.data(), "Muon_isRPCcand[nMuon]/O");
    t_event->Branch("Muon_nValid", Muon_nValid.data(), "Muon_nValid[nMuon]/I");
    t_event->Branch("Muon_nPix", Muon_nPix.data(), "Muon_nPix[nMuon]/I");
    t_event->Branch("Muon_Chi2", Muon_Chi2.data(), "Muon_Chi2[nMuon]/F");
    t_event->Branch("Muon_gnValid", Muon_gnValid.data(), "Muon_gnValid[nMuon]/I");
    t_event->Branch("Muon_gnPix", Muon_gnPix.data(), "Muon_gnPix[nMuon]/I");
    t_event->Branch("Muon_gChi2", Muon_gChi2.data(), "Muon_gChi2[nMuon]/F");
    t_event->Branch("Muon_gnValidMu", Muon_gnValidMu.data(), "Muon_gnValidMu[nMuon]/I");
    t_event->Branch("Muon_vtxIdx", Muon_vtxIdx.data(), "Muon_vtxIdx[nMuon]/I");
    t_event->Branch("Muon_vtxFlag", Muon_vtxFlag.data(), "Muon_vtxFlag[nMuon]/I");
    t_event->Branch("Muon_trkIdx", Muon_trkIdx.data(), "Muon_trkIdx[nMuon]/I");
    t_event->Branch("Muon_simIdx", Muon_simIdx.data(), "Muon_simIdx[nMuon]/I");
    // temporary
    t_event->Branch("b4_nMuon", &b4_nMuon, "b4_nMuon/i");  
    t_event->Branch("Muon_nNano", &Muon_nNano, "Muon_nNano/i");  

    //
    // Dimuon branches 
    //
    t_event->Branch("nDimu", &nDimu, "nDimu/I");  
    t_event->Branch("Dimut1_muIdx", Dimut1_muIdx.data(), "Dimut1_muIdx[nDimu]/I");
    t_event->Branch("Dimut1_dxy", Dimut1_dxy.data(), "Dimut1_dxy[nDimu]/F");
    t_event->Branch("Dimut1_dz", Dimut1_dz.data(), "Dimut1_dz[nDimu]/F");
    t_event->Branch("Dimut2_muIdx", Dimut2_muIdx.data(), "Dimut2_muIdx[nDimu]/I");
    t_event->Branch("Dimut2_dxy", Dimut2_dxy.data(), "Dimut2_dxy[nDimu]/F");
    t_event->Branch("Dimut2_dz", Dimut2_dz.data(), "Dimut2_dz[nDimu]/F");

    t_event->Branch("Dimu_pt", Dimu_pt.data(), "Dimu_pt[nDimu]/F");
    t_event->Branch("Dimu_eta", Dimu_eta.data(), "Dimu_eta[nDimu]/F");
    t_event->Branch("Dimu_phi", Dimu_phi.data(), "Dimu_phi[nDimu]/F");
    t_event->Branch("Dimu_rap", Dimu_rap.data(), "Dimu_rap[nDimu]/F");
    t_event->Branch("Dimu_mass", Dimu_mass.data(), "Dimu_mass[nDimu]/F");
    t_event->Branch("Dimu_charge", Dimu_charge.data(), "Dimu_charge[nDimu]/I");
    t_event->Branch("Dimu_simIdx", Dimu_simIdx.data(), "Dimu_simIdx[nDimu]/I");
    t_event->Branch("Dimu_vtxIdx", Dimu_vtxIdx.data(), "Dimu_vtxIdx[nDimu]/I");
    t_event->Branch("Dimu_chi2", Dimu_chi2.data(), "Dimu_chi2[nDimu]/F");
    t_event->Branch("Dimu_dlxy", Dimu_dlxy.data(), "Dimu_dlxy[nDimu]/F");
    t_event->Branch("Dimu_dlxyErr", Dimu_dlxyErr.data(), "Dimu_dlxyErr[nDimu]/F");
    t_event->Branch("Dimu_dlxySig", Dimu_dlxySig.data(), "Dimu_dlxySig[nDimu]/F");
    t_event->Branch("Dimu_cosphixy", Dimu_cosphixy.data(), "Dimu_cosphixy[nDimu]/F");
    t_event->Branch("Dimu_dl", Dimu_dl.data(), "Dimu_dl[nDimu]/F");
    t_event->Branch("Dimu_dlErr", Dimu_dlErr.data(), "Dimu_dlErr[nDimu]/F");
    t_event->Branch("Dimu_dlSig", Dimu_dlSig.data(), "Dimu_dlSig[nDimu]/F");
    t_event->Branch("Dimu_cosphi", Dimu_cosphi.data(), "Dimu_cosphi[nDimu]/F");
    t_event->Branch("Dimu_ptfrac", Dimu_ptfrac.data(), "Dimu_ptfrac[nDimu]/F");
    t_event->Branch("Dimu_x", Dimu_x.data(), "Dimu_x[nDimu]/F");
    t_event->Branch("Dimu_y", Dimu_y.data(), "Dimu_y[nDimu]/F");
    t_event->Branch("Dimu_z", Dimu_z.data(), "Dimu_z[nDimu]/F");
    if (covout) {
      t_event->Branch("Dimu_Covxx", Dimu_Covxx.data(), "Dimu_Covxx[nDimu]/F");
      t_event->Branch("Dimu_Covyx", Dimu_Covyx.data(), "Dimu_Covyx[nDimu]/F");
      t_event->Branch("Dimu_Covzx", Dimu_Covzx.data(), "Dimu_Covzx[nDimu]/F");
      t_event->Branch("Dimu_Covyy", Dimu_Covyy.data(), "Dimu_Covyy[nDimu]/F");
      t_event->Branch("Dimu_Covzy", Dimu_Covzy.data(), "Dimu_Covzy[nDimu]/F");
      t_event->Branch("Dimu_Covzz", Dimu_Covzz.data(), "Dimu_Covzz[nDimu]/F");
    } // covout

  } // nanoext

  //----------------------- Stefan Wunsch variables + more -------------------//
  // official nanoAOD subset //

  // Electrons
  t_event->Branch("nElectron", &value_el_n, "nElectron/i");
  t_event->Branch("Electron_pt", value_el_pt, "Electron_pt[nElectron]/F");
  t_event->Branch("Electron_eta", value_el_eta, "Electron_eta[nElectron]/F");
  t_event->Branch("Electron_phi", value_el_phi, "Electron_phi[nElectron]/F");
  t_event->Branch("Electron_mass", value_el_mass, "Electron_mass[nElectron]/F");
  t_event->Branch("Electron_charge", value_el_charge, "Electron_charge[nElectron]/I");
  // nuha
  t_event->Branch("Electron_tightCharge", value_el_tightCharge, "Electron_tightCharge[nElectron]/I");
 
  t_event->Branch("Electron_pfRelIso03_all", value_el_pfreliso03all, "Electron_pfRelIso03_all[nElectron]/F");
  t_event->Branch("Electron_pfRelIso03_chg", value_el_pfreliso03chg, "Electron_pfRelIso03_chg[nElectron]/F");
  // nuha
  t_event->Branch("Electron_dr03TkSumPtOld", value_el_dr03TkSumPtOld, "Electron_dr03TkSumPtOld[nElectron]/F");
  t_event->Branch("Electron_dr03TkSumPt", value_el_dr03TkSumPt, "Electron_dr03TkSumPt[nElectron]/F");  
  t_event->Branch("Electron_dr03EcalRecHitSumEtOld", value_el_dr03EcalRecHitSumEtOld, "Electron_dr03EcalRecHitSumEtOld[nElectron]/F");
  t_event->Branch("Electron_dr03EcalRecHitSumEt", value_el_dr03EcalRecHitSumEt, "Electron_dr03EcalRecHitSumEt[nElectron]/F");
  
  // nanoAOD extension
  t_event->Branch("Electron_dr03HcalTowerSumEt", value_el_dr03HcalTowerSumEt, "Electron_dr03HcalTowerSumEt[nElectron]/F");
  // nuha
  t_event->Branch("Electron_dr03HcalDepth1TowerSumEtOld", value_el_dr03HcalDepth1TowerSumEtOld, "Electron_dr03HcalDepth1TowerSumEtOld[nElectron]/F");
  t_event->Branch("Electron_dr03HcalDepth1TowerSumEt", value_el_dr03HcalDepth1TowerSumEt, "Electron_dr03HcalDepth1TowerSumEt[nElectron]/F");

  // nanoAOD extension (next two)
  t_event->Branch("Electron_isEB", value_el_isEB, "Electron_isEB[nElectron]/O");
  t_event->Branch("Electron_isEE", value_el_isEE, "Electron_isEE[nElectron]/O");
  t_event->Branch("Electron_lostHits", value_el_lostHits, "Electron_lostHits[nElectron]/b");
  t_event->Branch("Electron_isPFcand", value_el_isPFcand, "Electron_isPFcand[nElectron]/O");
  t_event->Branch("Electron_isNano", value_el_isNano, "Electron_isNano[nElectron]/O");
    
  // nanoAOD extension (next two)
  t_event->Branch("Electron_convDist", value_el_convDist, "Electron_convDist[nElectron]/F");
  t_event->Branch("Electron_convDcot", value_el_convDcot, "Electron_convDcot[nElectron]/F");
  // nuha
  t_event->Branch("Electron_convVetoOld", value_el_convVetoOld, "Electron_convVetoOld[nElectron]/O");
  t_event->Branch("Electron_convVeto", value_el_convVeto, "Electron_convVeto[nElectron]/O");
  
  t_event->Branch("Electron_deltaEtaSC", value_el_deltaEtaSC, "Electron_deltaEtaSC[nElectron]/F");
  t_event->Branch("Electron_deltaPhiSC", value_el_deltaPhiSC, "Electron_deltaPhiSC[nElectron]/F");
  t_event->Branch("Electron_deltaEtaSCtr", value_el_deltaEtaSCtr, "Electron_deltaEtaSCtr[nElectron]/F");
  t_event->Branch("Electron_deltaPhiSCtr", value_el_deltaPhiSCtr, "Electron_deltaPhiSCtr[nElectron]/F");
  t_event->Branch("Electron_hoe", value_el_hoe, "Electron_hoe[nElectron]/F");
  t_event->Branch("Electron_sieie", value_el_sieie, "Electron_sieie[nElectron]/F");
  t_event->Branch("Electron_sieieR1", value_el_sieieR1, "Electron_sieieR1[nElectron]/F");
  // nuha
  t_event->Branch("Electron_eInvMinusPInvOld", value_el_eInvMinusPInvOld, "Electron_eInvMinusPInvOld[nElectron]/F");
  t_event->Branch("Electron_eInvMinusPInv", value_el_eInvMinusPInv, "Electron_eInvMinusPInv[nElectron]/F");

  // nanoAOD extension
  t_event->Branch("Electron_SCeta", value_el_SCeta, "Electron_SCeta[nElectron]/F");
  t_event->Branch("Electron_cutBased", value_el_cutBased, "Electron_cutBased[nElectron]/I");

  t_event->Branch("Electron_dxy", value_el_dxy, "Electron_dxy[nElectron]/F");
  t_event->Branch("Electron_dxyErr", value_el_dxyErr, "Electron_dxyErr[nElectron]/F");
  t_event->Branch("Electron_dz", value_el_dz, "Electron_dz[nElectron]/F");
  t_event->Branch("Electron_dzErr", value_el_dzErr, "Electron_dzErr[nElectron]/F");
  t_event->Branch("Electron_ip3d", value_el_ip3d, "Electron_ip3d[nElectron]/F");
  t_event->Branch("Electron_sip3d", value_el_sip3d, "Electron_sip3d[nElectron]/F");
  t_event->Branch("Electron_nNano", &Electron_nNano, "Electron_nNano/i");  

  // Taus
  t_event->Branch("nTau", &value_tau_n, "nTau/i");
  t_event->Branch("Tau_pt", value_tau_pt, "Tau_pt[nTau]/F");
  t_event->Branch("Tau_eta", value_tau_eta, "Tau_eta[nTau]/F");
  t_event->Branch("Tau_phi", value_tau_phi, "Tau_phi[nTau]/F");
  t_event->Branch("Tau_mass", value_tau_mass, "Tau_mass[nTau]/F");
  t_event->Branch("Tau_charge", value_tau_charge, "Tau_charge[nTau]/I");
  t_event->Branch("Tau_decayMode", value_tau_decaymode, "Tau_decayMode[nTau]/I");
  t_event->Branch("Tau_chargedIso", value_tau_chargediso, "Tau_chargedIso[nTau]/F");
  t_event->Branch("Tau_neutralIso", value_tau_neutraliso, "Tau_neutralIso[nTau]/F");

  // Photons
  t_event->Branch("nPhoton", &value_ph_n, "nPhoton/i");
  t_event->Branch("Photon_pt", value_ph_pt, "Photon_pt[nPhoton]/F");
  t_event->Branch("Photon_eta", value_ph_eta, "Photon_eta[nPhoton]/F");
  t_event->Branch("Photon_phi", value_ph_phi, "Photon_phi[nPhoton]/F");
  t_event->Branch("Photon_mass", value_ph_mass, "Photon_mass[nPhoton]/F");
  t_event->Branch("Photon_charge", value_ph_charge, "Photon_charge[nPhoton]/I");
  t_event->Branch("Photon_pfRelIso03_all", value_ph_pfreliso03all, "Photon_pfRelIso03_all[nPhoton]/F");

  // MET
  t_event->Branch("MET_pt", &value_met_pt, "MET_pt/F");
  t_event->Branch("MET_phi", &value_met_phi, "MET_phi/F");
  t_event->Branch("MET_sumEt", &value_met_sumEt, "MET_sumEt/F");
  t_event->Branch("MET_significance", &value_met_significance, "MET_significance/F");
  t_event->Branch("MET_covXX", &value_met_covxx, "MET_covXX/F");
  t_event->Branch("MET_covXY", &value_met_covxy, "MET_covXY/F");
  t_event->Branch("MET_covYY", &value_met_covyy, "MET_covYY/F");

  // caloMET
  t_event->Branch("CaloMET_pt", &value_calomet_pt, "CaloMET_pt/F");
  t_event->Branch("CaloMET_phi", &value_calomet_phi, "CaloMET_phi/F");
  t_event->Branch("CaloMET_sumEt", &value_calomet_sumEt, "CaloMET_sumEt/F");

  // Jets
  t_event->Branch("nJet", &value_jet_n, "nJet/i");
  t_event->Branch("Jet_pt", value_jet_pt, "Jet_pt[nJet]/F");
  t_event->Branch("Jet_eta", value_jet_eta, "Jet_eta[nJet]/F");
  t_event->Branch("Jet_phi", value_jet_phi, "Jet_phi[nJet]/F");
  t_event->Branch("Jet_mass", value_jet_mass, "Jet_mass[nJet]/F");
//  t_event->Branch("Jet_ptD", value_jet_ptD, "Jet_ptD[nJet]/F");
  t_event->Branch("Jet_area", value_jet_area, "Jet_area[nJet]/F");
  t_event->Branch("Jet_nConstituents", value_jet_nConstituents, "Jet_nConstituents[nJet]/i");
  t_event->Branch("Jet_nElectrons", value_jet_nElectrons, "Jet_nElectrons[nJet]/i");
  t_event->Branch("Jet_nMuons", value_jet_nMuons, "Jet_nMuons[nJet]/i");
  t_event->Branch("Jet_chEmEF", value_jet_chEmEF, "Jet_chEmEF[nJet]/F");
  t_event->Branch("Jet_chHEF", value_jet_chHEF, "Jet_chHEF[nJet]/F");
  t_event->Branch("Jet_neEmEF", value_jet_neEmEF, "Jet_neEmEF[nJet]/F");
  t_event->Branch("Jet_neHEF", value_jet_neHEF, "Jet_neHEF[nJet]/F");

  // Flags
  if (!custom_flag.empty())
    custom_bit.reserve(custom_flag.size());
  for (size_t iFlag = 0; iFlag < custom_flag.size(); ++iFlag) {
    custom_bit.push_back(0);
    t_event->Branch(("Flag_" + custom_flag.at(iFlag)).c_str(), &custom_bit.at(iFlag), ("Flag_" + custom_flag.at(iFlag) + "/O").c_str());
  }

  //----------------------- Create branch of Dmeson's tree ---------------------//
  
  if (nanoext) {
    // * = only booked branch, filled it with nonsense atm 
    t_event->Branch("nD0", &nD0, "nD0/I");  
    t_event->Branch("D0t1_pt", D0t1_pt.data(), "D0t1_pt[nD0]/F");
    t_event->Branch("D0t1_eta", D0t1_eta.data(), "D0t1_eta[nD0]/F");
    t_event->Branch("D0t1_phi", D0t1_phi.data(), "D0t1_phi[nD0]/F");
    t_event->Branch("D0t1_chg", D0t1_chg.data(), "D0t1_chg[nD0]/I");
    t_event->Branch("D0t1_tkIdx", D0t1_tkIdx.data(), "D0t1_tkIdx[nD0]/I"); // *
    t_event->Branch("D0t1_Kprob", D0t1_Kprob.data(), "D0t1_Kprob[nD0]/F"); // *
    t_event->Branch("D0t1_piprob", D0t1_piprob.data(), "D0t1_piprob[nD0]/F"); // *
    t_event->Branch("D0t1_dEdxnmeas", D0t1_dEdxnmeas.data(), "D0t1_dEdxnmeas[nD0]/I"); // *
    t_event->Branch("D0t1_dEdxnsat", D0t1_dEdxnsat.data(), "D0t1_dEdxnsat[nD0]/I"); // *
    t_event->Branch("D0t1_vtxIdx", D0t1_vtxIdx.data(), "D0t1_vtxIdx[nD0]/I");
    t_event->Branch("D0t1_chindof", D0t1_chindof.data(), "D0t1_chindof[nD0]/F");
    t_event->Branch("D0t1_nValid", D0t1_nValid.data(), "D0t1_nValid[nD0]/I");
    t_event->Branch("D0t1_nPix", D0t1_nPix.data(), "D0t1_nPix[nD0]/I");
    t_event->Branch("D0t1_isHighPurity", D0t1_isHighPurity.data(), "D0t1_isHighPurity[nD0]/O");
    t_event->Branch("D0t1_dxy", D0t1_dxy.data(), "D0t1_dxy[nD0]/F");
    t_event->Branch("D0t1_dz", D0t1_dz.data(), "D0t1_dz[nD0]/F");
    t_event->Branch("D0t1_pdgId", D0t1_pdgId.data(), "D0t1_pdgId[nD0]/I");
    t_event->Branch("D0t2_pt", D0t2_pt.data(), "D0t2_pt[nD0]/F");
    t_event->Branch("D0t2_eta", D0t2_eta.data(), "D0t2_eta[nD0]/F");
    t_event->Branch("D0t2_phi", D0t2_phi.data(), "D0t2_phi[nD0]/F");
    t_event->Branch("D0t2_chg", D0t2_chg.data(), "D0t2_chg[nD0]/I");
    t_event->Branch("D0t2_tkIdx", D0t2_tkIdx.data(), "D0t2_tkIdx[nD0]/I"); // *
    t_event->Branch("D0t2_Kprob", D0t2_Kprob.data(), "D0t2_Kprob[nD0]/F"); // *
    t_event->Branch("D0t2_piprob", D0t2_piprob.data(), "D0t2_piprob[nD0]/F"); // *
    t_event->Branch("D0t2_dEdxnmeas", D0t2_dEdxnmeas.data(), "D0t2_dEdxnmeas[nD0]/I"); // *
    t_event->Branch("D0t2_dEdxnsat", D0t2_dEdxnsat.data(), "D0t2_dEdxnsat[nD0]/I"); // *
    t_event->Branch("D0t2_vtxIdx", D0t2_vtxIdx.data(), "D0t2_vtxIdx[nD0]/I");
    t_event->Branch("D0t2_chindof", D0t2_chindof.data(), "D0t2_chindof[nD0]/F");
    t_event->Branch("D0t2_nValid", D0t2_nValid.data(), "D0t2_nValid[nD0]/I");
    t_event->Branch("D0t2_nPix", D0t2_nPix.data(), "D0t2_nPix[nD0]/I");
    t_event->Branch("D0t2_isHighPurity", D0t2_isHighPurity.data(), "D0t2_isHighPurity[nD0]/O");
    t_event->Branch("D0t2_dxy", D0t2_dxy.data(), "D0t2_dxy[nD0]/F");
    t_event->Branch("D0t2_dz", D0t2_dz.data(), "D0t2_dz[nD0]/F");
    t_event->Branch("D0t2_pdgId", D0t2_pdgId.data(), "D0t2_pdgId[nD0]/I");
    t_event->Branch("D0_pt", D0_pt.data(), "D0_pt[nD0]/F");
    t_event->Branch("D0_eta", D0_eta.data(), "D0_eta[nD0]/F");
    t_event->Branch("D0_phi", D0_phi.data(), "D0_phi[nD0]/F");
    t_event->Branch("D0_rap", D0_rap.data(), "D0_rap[nD0]/F");
    t_event->Branch("D0_mass12", D0_mass12.data(), "D0_mass12[nD0]/F");
    t_event->Branch("D0_mass21", D0_mass21.data(), "D0_mass21[nD0]/F");
    t_event->Branch("D0_massKK", D0_massKK.data(), "D0_massKK[nD0]/F");
    t_event->Branch("D0_simIdx", D0_simIdx.data(), "D0_simIdx[nD0]/I");
    t_event->Branch("D0_DstarIdx", D0_DstarIdx.data(), "D0_DstarIdx[nD0]/I");
    t_event->Branch("D0_ambiPrim", D0_ambiPrim.data(), "D0_ambiPrim[nD0]/O");
    t_event->Branch("D0_vtxIdx", D0_vtxIdx.data(), "D0_vtxIdx[nD0]/I");
    t_event->Branch("D0_hasMuon", D0_hasMuon.data(), "D0_hasMuon[nD0]/O");
    t_event->Branch("D0_chi2", D0_chi2.data(), "D0_chi2[nD0]/F");
    t_event->Branch("D0_dlxy", D0_dlxy.data(), "D0_dlxy[nD0]/F");
    t_event->Branch("D0_dlxyErr", D0_dlxyErr.data(), "D0_dlxyErr[nD0]/F");
    t_event->Branch("D0_dlxySig", D0_dlxySig.data(), "D0_dlxySig[nD0]/F");
    t_event->Branch("D0_cosphixy", D0_cosphixy.data(), "D0_cosphixy[nD0]/F");
    t_event->Branch("D0_dl", D0_dl.data(), "D0_dl[nD0]/F");
    t_event->Branch("D0_dlErr", D0_dlErr.data(), "D0_dlErr[nD0]/F");
    t_event->Branch("D0_dlSig", D0_dlSig.data(), "D0_dlSig[nD0]/F");
    t_event->Branch("D0_cosphi", D0_cosphi.data(), "D0_cosphi[nD0]/F");
    t_event->Branch("D0_ptfrac", D0_ptfrac.data(), "D0_ptfrac[nD0]/F");
    t_event->Branch("D0_ptfrac15", D0_ptfrac15.data(), "D0_ptfrac15[nD0]/F");
    t_event->Branch("D0_ptfrac10", D0_ptfrac10.data(), "D0_ptfrac10[nD0]/F");
    t_event->Branch("D0_ptfrac07", D0_ptfrac07.data(), "D0_ptfrac07[nD0]/F");
    t_event->Branch("D0_ptfrac04", D0_ptfrac04.data(), "D0_ptfrac04[nD0]/F");
    t_event->Branch("D0_x", D0_x.data(), "D0_x[nD0]/F");
    t_event->Branch("D0_y", D0_y.data(), "D0_y[nD0]/F");
    t_event->Branch("D0_z", D0_z.data(), "D0_z[nD0]/F");
    if (covout) {
      t_event->Branch("D0_Covxx", D0_Covxx.data(), "D0_Covxx[nD0]/F");
      t_event->Branch("D0_Covyx", D0_Covyx.data(), "D0_Covyx[nD0]/F");
      t_event->Branch("D0_Covzx", D0_Covzx.data(), "D0_Covzx[nD0]/F");
      t_event->Branch("D0_Covyy", D0_Covyy.data(), "D0_Covyy[nD0]/F");
      t_event->Branch("D0_Covzy", D0_Covzy.data(), "D0_Covzy[nD0]/F");
      t_event->Branch("D0_Covzz", D0_Covzz.data(), "D0_Covzz[nD0]/F");
    }

    t_event->Branch("nDstar", &nDstar, "nDstar/I");    
    t_event->Branch("Dstarpis_pt", Dstarpis_pt.data(), "Dstarpis_pt[nDstar]/F");
    t_event->Branch("Dstarpis_eta", Dstarpis_eta.data(), "Dstarpis_eta[nDstar]/F");
    t_event->Branch("Dstarpis_phi", Dstarpis_phi.data(), "Dstarpis_phi[nDstar]/F");
    t_event->Branch("Dstarpis_ptr", Dstarpis_ptr.data(), "Dstarpis_ptr[nDstar]/F");
    t_event->Branch("Dstarpis_etar", Dstarpis_etar.data(), "Dstarpis_etar[nDstar]/F");
    t_event->Branch("Dstarpis_phir", Dstarpis_phir.data(), "Dstarpis_phir[nDstar]/F");
    t_event->Branch("Dstarpis_chg", Dstarpis_chg.data(), "Dstarpis_chg[nDstar]/I");
    t_event->Branch("Dstarpis_tkIdx", Dstarpis_tkIdx.data(), "Dstarpis_tkIdx[nDstar]/I");
    t_event->Branch("Dstarpis_Kprob", Dstarpis_Kprob.data(), "Dstarpis_Kprob[nDstar]/F"); // *
    t_event->Branch("Dstarpis_piprob", Dstarpis_piprob.data(), "Dstarpis_piprob[nDstar]/F"); // *
    t_event->Branch("Dstarpis_dEdxnmeas", Dstarpis_dEdxnmeas.data(), "Dstarpis_dEdxnmeas[nDstar]/I"); // *
    t_event->Branch("Dstarpis_dEdxnsat", Dstarpis_dEdxnsat.data(), "Dstarpis_dEdxnsat[nDstar]/I"); // *
    t_event->Branch("Dstarpis_vtxIdx", Dstarpis_vtxIdx.data(), "Dstarpis_vtxIdx[nDstar]/I");
    t_event->Branch("Dstarpis_chindof", Dstarpis_chindof.data(), "Dstarpis_chindof[nDstar]/F");
    t_event->Branch("Dstarpis_chir", Dstarpis_chir.data(), "Dstarpis_chir[nDstar]/F");
    t_event->Branch("Dstarpis_nValid", Dstarpis_nValid.data(), "Dstarpis_nValid[nDstar]/I");
    t_event->Branch("Dstarpis_nPix", Dstarpis_nPix.data(), "Dstarpis_nPix[nDstar]/I");
    t_event->Branch("Dstarpis_dxy", Dstarpis_dxy.data(), "Dstarpis_dxy[nDstar]/F");
    t_event->Branch("Dstarpis_dz", Dstarpis_dz.data(), "Dstarpis_dz[nDstar]/F");
  
    t_event->Branch("DstarD0_pt", DstarD0_pt.data(), "DstarD0_pt[nDstar]/F");
    t_event->Branch("DstarD0_eta", DstarD0_eta.data(), "DstarD0_eta[nDstar]/F");
    t_event->Branch("DstarD0_phi", DstarD0_phi.data(), "DstarD0_phi[nDstar]/F");
    t_event->Branch("DstarD0_mass", DstarD0_mass.data(), "DstarD0_mass[nDstar]/F");
    t_event->Branch("DstarD0_chi2", DstarD0_chi2.data(), "DstarD0_chi2[nDstar]/F");
    t_event->Branch("DstarD0_dlxy", DstarD0_dlxy.data(), "DstarD0_dlxy[nDstar]/F");
    t_event->Branch("DstarD0_dlxyErr", DstarD0_dlxyErr.data(), "DstarD0_dlxyErr[nDstar]/F");
    t_event->Branch("DstarD0_dlxySig", DstarD0_dlxySig.data(), "DstarD0_dlxySig[nDstar]/F");
    t_event->Branch("DstarD0_cosphixy", DstarD0_cosphixy.data(), "DstarD0_cosphixy[nDstar]/F");
    t_event->Branch("DstarD0_dl", DstarD0_dl.data(), "DstarD0_dl[nDstar]/F");
    t_event->Branch("DstarD0_dlErr", DstarD0_dlErr.data(), "DstarD0_dlErr[nDstar]/F");
    t_event->Branch("DstarD0_dlSig", DstarD0_dlSig.data(), "DstarD0_dlSig[nDstar]/F");
    t_event->Branch("DstarD0_cosphi", DstarD0_cosphi.data(), "DstarD0_cosphi[nDstar]/F");
    t_event->Branch("DstarD0_ptfrac", DstarD0_ptfrac.data(), "DstarD0_ptfrac[nDstar]/F");
    t_event->Branch("DstarD0_ptfrac15", DstarD0_ptfrac15.data(), "DstarD0_ptfrac15[nDstar]/F");
    t_event->Branch("DstarD0_ptfrac10", DstarD0_ptfrac10.data(), "DstarD0_ptfrac10[nDstar]/F");
    t_event->Branch("DstarD0_ptfrac07", DstarD0_ptfrac07.data(), "DstarD0_ptfrac07[nDstar]/F");
    t_event->Branch("DstarD0_ptfrac04", DstarD0_ptfrac04.data(), "DstarD0_ptfrac04[nDstar]/F");
    t_event->Branch("DstarD0_x", DstarD0_x.data(), "DstarD0_x[nDstar]/F");
    t_event->Branch("DstarD0_y", DstarD0_y.data(), "DstarD0_y[nDstar]/F");
    t_event->Branch("DstarD0_z", DstarD0_z.data(), "DstarD0_z[nDstar]/F");
    t_event->Branch("DstarD0_simIdx", DstarD0_simIdx.data(), "DstarD0_simIdx[nDstar]/I");
    t_event->Branch("DstarD0_recIdx", DstarD0_recIdx.data(), "DstarD0_recIdx[nDstar]/I");
    t_event->Branch("DstarD0_ambiPrim", DstarD0_ambiPrim.data(), "DstarD0_ambiPrim[nDstar]/O");

    t_event->Branch("DstarK_pt", DstarK_pt.data(), "DstarK_pt[nDstar]/F");
    t_event->Branch("DstarK_eta", DstarK_eta.data(), "DstarK_eta[nDstar]/F");
    t_event->Branch("DstarK_phi", DstarK_phi.data(), "DstarK_phi[nDstar]/F");
    t_event->Branch("DstarK_chg", DstarK_chg.data(), "DstarK_chg[nDstar]/I");
    t_event->Branch("DstarK_tkIdx", DstarK_tkIdx.data(), "DstarK_tkIdx[nDstar]/I");
    t_event->Branch("DstarK_Kprob", DstarK_Kprob.data(), "DstarK_Kprob[nDstar]/F"); // *
    t_event->Branch("DstarK_piprob", DstarK_piprob.data(), "DstarK_piprob[nDstar]/F"); // *
    t_event->Branch("DstarK_dEdxnmeas", DstarK_dEdxnmeas.data(), "DstarK_dEdxnmeas[nDstar]/I"); // *
    t_event->Branch("DstarK_dEdxnsat", DstarK_dEdxnsat.data(), "DstarK_dEdxnsat[nDstar]/I"); // *
    t_event->Branch("DstarK_vtxIdx", DstarK_vtxIdx.data(), "DstarK_vtxIdx[nDstar]/I");
    t_event->Branch("DstarK_chindof", DstarK_chindof.data(), "DstarK_chindof[nDstar]/F");
    t_event->Branch("DstarK_nValid", DstarK_nValid.data(), "DstarK_nValid[nDstar]/I");
    t_event->Branch("DstarK_nPix", DstarK_nPix.data(), "DstarK_nPix[nDstar]/I");
    t_event->Branch("DstarK_isHighPurity", DstarK_isHighPurity.data(), "DstarK_isHighPurity[nDstar]/O");
    t_event->Branch("DstarK_dxy", DstarK_dxy.data(), "DstarK_dxy[nDstar]/F");
    t_event->Branch("DstarK_dz", DstarK_dz.data(), "DstarK_dz[nDstar]/F");

    t_event->Branch("Dstarpi_pt", Dstarpi_pt.data(), "Dstarpi_pt[nDstar]/F");
    t_event->Branch("Dstarpi_eta", Dstarpi_eta.data(), "Dstarpi_eta[nDstar]/F");
    t_event->Branch("Dstarpi_phi", Dstarpi_phi.data(), "Dstarpi_phi[nDstar]/F");
    t_event->Branch("Dstarpi_chg", Dstarpi_chg.data(), "Dstarpi_chg[nDstar]/I");
    t_event->Branch("Dstarpi_tkIdx", Dstarpi_tkIdx.data(), "Dstarpi_tkIdx[nDstar]/I");
    t_event->Branch("Dstarpi_Kprob", Dstarpi_Kprob.data(), "Dstarpi_Kprob[nDstar]/F"); // *
    t_event->Branch("Dstarpi_piprob", Dstarpi_piprob.data(), "Dstarpi_piprob[nDstar]/F"); // *
    t_event->Branch("Dstarpi_dEdxnmeas", Dstarpi_dEdxnmeas.data(), "Dstarpi_dEdxnmeas[nDstar]/I"); // *
    t_event->Branch("Dstarpi_dEdxnsat", Dstarpi_dEdxnsat.data(), "Dstarpi_dEdxnsat[nDstar]/I"); // *
    t_event->Branch("Dstarpi_vtxIdx", Dstarpi_vtxIdx.data(), "Dstarpi_vtxIdx[nDstar]/I");
    t_event->Branch("Dstarpi_chindof", Dstarpi_chindof.data(), "Dstarpi_chindof[nDstar]/F");
    t_event->Branch("Dstarpi_nValid", Dstarpi_nValid.data(), "Dstarpi_nValid[nDstar]/I");
    t_event->Branch("Dstarpi_nPix", Dstarpi_nPix.data(), "Dstarpi_nPix[nDstar]/I");
    t_event->Branch("Dstarpi_isHighPurity", Dstarpi_isHighPurity.data(), "Dstarpi_isHighPurity[nDstar]/O");
    t_event->Branch("Dstarpi_dxy", Dstarpi_dxy.data(), "Dstarpi_dxy[nDstar]/F");
    t_event->Branch("Dstarpi_dz", Dstarpi_dz.data(), "Dstarpi_dz[nDstar]/F");
  
    t_event->Branch("Dstar_pt", Dstar_pt.data(), "Dstar_pt[nDstar]/F");
    t_event->Branch("Dstar_eta", Dstar_eta.data(), "Dstar_eta[nDstar]/F");
    t_event->Branch("Dstar_phi", Dstar_phi.data(), "Dstar_phi[nDstar]/F");
    t_event->Branch("Dstar_rap", Dstar_rap.data(), "Dstar_rap[nDstar]/F");
    t_event->Branch("Dstar_deltam", Dstar_deltam.data(), "Dstar_deltam[nDstar]/F");
    t_event->Branch("Dstar_deltamr", Dstar_deltamr.data(), "Dstar_deltamr[nDstar]/F");
    t_event->Branch("Dstar_simIdx", Dstar_simIdx.data(), "Dstar_simIdx[nDstar]/I");
    t_event->Branch("Dstar_vtxIdx", Dstar_vtxIdx.data(), "Dstar_vtxIdx[nDstar]/I");
    t_event->Branch("Dstar_hasMuon", Dstar_hasMuon.data(), "Dstar_hasMuon[nDstar]/O");
    t_event->Branch("Dstar_ptfrac", Dstar_ptfrac.data(), "Dstar_ptfrac[nDstar]/F");
  } // end of nanoext

// *** move to different place? ***
  // Qun below  Trigger Object
  t_event->Branch("nTrigObj", &nTrigObj, "nTrigObj/I");    
  t_event->Branch("TrigObj_id", TrigObj_id.data(), "TrigObj_id[nTrigObj]/I");
  t_event->Branch("TrigObj_filterBits", TrigObj_filterBits.data(), "TrigObj_FilterBits[nTrigObj]/I");
  t_event->Branch("TrigObj_pt", TrigObj_pt.data(), "TrigObj_pt[nTrigObj]/F");
  t_event->Branch("TrigObj_phi", TrigObj_phi.data(), "TrigObj_phi[nTrigObj]/F");
  t_event->Branch("TrigObj_eta", TrigObj_eta.data(), "TrigObj_eta[nTrigObj]/F");
  //t_event->Branch("TrigObj_id", &TrigObj_id);
  //t_event->Branch("TrigObj_pt", &TrigObj_pt);
  //t_event->Branch("TrigObj_phi", &TrigObj_phi);
  //t_event->Branch("TrigObj_eta", &TrigObj_eta);
  //t_event->Branch("TrigObj_id", TrigObj_id, "TrigObj_id[nTrigObj]/I");
  //t_event->Branch("TrigObj_pt", TrigObj_pt, "TrigObj_pt[nTrigObj]/F");
  //t_event->Branch("TrigObj_phi", TrigObj_phi, "TrigObj_phi[nTrigObj]/F");
  //t_event->Branch("TrigObj_eta", TrigObj_eta, "TrigObj_eta[nTrigObj]/F");
  //Qt_event->Branch("TrigObj_id", TrigObj_id.data(), "TrigObj_id[nTrigObj]/I");
  //Qt_event->Branch("TrigObj_pt", TrigObj_pt.data(), "TrigObj_pt[nTrigObj]/F");
  //Qt_event->Branch("TrigObj_phi", TrigObj_phi.data(), "TrigObj_phi[nTrigObj]/F");
  //Qt_event->Branch("TrigObj_eta", TrigObj_eta.data(), "TrigObj_eta[nTrigObj]/F");
  // Qun above  Trigger Object   

} // end of beginJob

// ------------ method called once each job just after ending the event loop  ------------
void 
NanoAnalyzer::endJob()
{
  file->cd();
  t_event->Write();
  
  if (!nanoext) return;

  h_trackpt->Write();
  h_trackptlow->Write();
  h_tracketa->Write();
  
  if (!nanoext) return;

  h_d0pt->Write();
  h_dstarpt->Write();
  h_PS3pt->Write();
  h_K1pt->Write();
  h_P2pt->Write();
  h_K1eta->Write();
  h_P2eta->Write();
  h_PS3eta->Write();
  h_D0masscut->Write();
  h_deltaMassD0Dstar->Write();
  h_D0masscut_rightDLcut->Write();
  h_deltaMassD0Dstar_rightDLcut->Write();

  // H4lepton
  h_p_e->Write();
  h_et_e->Write();
  h_pt_e_b4->Write();
  h_eta_e_b4->Write();
  h_phi_e->Write();
  h_sc_eta->Write();
  h_sc_rawE->Write();
  h_relPFIso_e->Write();
  h_relPFIso_pt_e->Write();
  h_dxy_e->Write();
  h_SIP3d_e_b4->Write();
  h_misshite->Write();  
}

/*  commented!
//    ------------ method called when starting to processes a run  ------------  Qun   below 
void NanoAnalyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
//--------------------------------------------------------------------------
{
    using namespace std;
    using namespace edm;
  //If the hltConfig can be initialized, then the below is an example of how to extract the config information for the trigger from the so-called provenance.
  //The trigger configuration can change from run to run (during the run is the same), so it needs to be called here.
  //"init" return value indicates whether intitialisation has succeeded
  //"changed" parameter indicates whether the config has actually changed
    bool changed(true);
 // *** this might be duplicate, investigate ***

    if (hltConfig_.init(iRun,iSetup,processName_,changed)) {
      if (changed) {
	if (triggerName_!="@") { // "@" means: analyze all triggers in config
	  const unsigned int n(hltConfig_.size());
	  const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName_));
	  if (triggerIndex>=n) {
            cout << "HLTEventAnalyzerAOD::analyze:"
               << " TriggerName " << triggerName_
               << " not available in (new) config!" << endl;
            cout << "Available TriggerNames are: " << endl;
            hltConfig_.dump("Triggers");
          }
            cout << "Available TriggerNames are:QQQQQQ " << endl;
	}
      }
    } else { 
  cout << "NanoAnalyzer::beginRun: config extraction failure with process name " << processName_<< endl;
    }

} //NanoAnalyzer::beginRun ends
// above Qun
// end of comment 
*/


// ------------ method called once each job just before starting run  ------------
void NanoAnalyzer::beginRun(const edm::Run &iRun, const edm::EventSetup &iStp)
{
  //If the hltConfig can be initialized, then the below is an example of how to extract the config information for the trigger from the so-called provenance.
  //The trigger configuration can change from run to run (during the run is the same), so it needs to be called here.
  //"init" return value indicates whether intitialisation has succeeded
  //"changed" parameter indicates whether the config has actually changed
  bool changed = false;
  cout << "beginRun" << endl;
  //if (!hlt_cfg.init(iRun, iStp, hlt_proc, changed)) {
  if (!hltConfig_.init(iRun, iStp, hlt_proc, changed)) {
    std::cout << "Initialization of HLTConfigProvider failed!!" << std::endl;
    std::cout << "This is normal on 2010 MC only" << endl;
    return;
  }
  if (changed) {
    // std::cout << "HLT menu used in run " << iRun.run() << " is " << hlt_cfg.tableName() << std::endl << std::endl;
    std::cout << "HLT menu used in run " << iRun.run() << " is " << hltConfig_.tableName() << std::endl << std::endl;
    /*
    std::cout << "The trigger paths in this menu are: " << std::endl;
    for (auto &p : hlt_cfg.triggerNames())
      std::cout << remove_version(p) << std::endl;

    std::cout << "There are " << hlt_cfg.triggerNames().size() << " paths in total." << std::endl << std::endl;
    */
    update_HLT_branch();
  }

// Qun, for trigger objects
	if (triggerName_!="@") { // "@" means: analyze all triggers in config
	  const unsigned int n(hltConfig_.size());
	  const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName_));
	  if (triggerIndex>=n) {
            std::cout << "HLTEventAnalyzerAOD::analyze:"
               << " TriggerName " << triggerName_
               << " not available in (new) config!" << endl;
            std::cout << "Available TriggerNames are: " << endl;
            hltConfig_.dump("Triggers");
          }
            cout << "Available TriggerNames are:QQQQQQ " << endl;
	}
// end Qun

  //cout << "hello beginRun end" << endl; 
}

// ------------ method for creating and updating HLT branches as necessary run by run ------------
void NanoAnalyzer::update_HLT_branch()
{
  // zero out all existing bits between runs, so that the triggers that stop existing is set to 0
  // instead of whatever value it has in the last event
  for (unordered_map<std::string, uint8_t>::iterator it = hlt_bit.begin(); it != hlt_bit.end(); ++it)
    it->second = 0;

  // for (uint iP = 0; iP < hlt_cfg.triggerNames().size(); ++iP) {
  //  std::string path = remove_version(hlt_cfg.triggerNames().at(iP));
  for (uint iP = 0; iP < hltConfig_.triggerNames().size(); ++iP) {
    std::string path = remove_version(hltConfig_.triggerNames().at(iP));

    // check if the path already exists
    // if not, create a new branch and pad it with 0 for previous events
    // only allow paths that start with HLT (presumably AlCa etc are not needed)
    if (!hlt_bit.count(path) and path.substr(0, 3) == "HLT") {
      hlt_bit.insert( std::make_pair(path, 0) );
      TBranch *hlt_path = t_event->Branch(path.c_str(), &hlt_bit.at(path), (path + "/O").c_str());  
      for (uint iE = 0; iE < t_event->GetEntries(); ++iE)
        hlt_path->Fill(); 
    }
  }
}

void NanoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions & descriptions) 
{
  edm::ParameterSetDescription desc;
  desc.add<std::string>("outFile", "test.root");
  desc.add<bool>("isData", true);
  desc.add<std::string>("hltProcess", "HLT");
  desc.add<std::string>("hlt_proc", "HLT");
#ifdef Compatibility
  desc.add<bool>("nanoExtension", false);
#else 
  desc.add<bool>("nanoExtension", true);
#endif
  desc.add<bool>("writeCovariance", false);
  desc.add<std::vector<std::string> >("customFlag", std::vector<std::string>());
  desc.add<edm::InputTag>("customTag", edm::InputTag());
  desc.add<std::string>("triggerName_", "@");
  desc.add<std::string>("triggerName", "@");
  descriptions.add("nanoAnalyzer", desc);
}

// ------------ method called when ending the processing of a run  ------------  Qun below
void NanoAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
// above Qun

//include JSON quality check
#include "NanoJSON.cc.forinclude"

//include trigger methods
#include "NanoTrigger.cc.forinclude"

//define this as a plug-in
DEFINE_FWK_MODULE(NanoAnalyzer);
