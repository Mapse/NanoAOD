// to be save of double includes;
#pragma once

// to use transient tracks:
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TrackTransientTrack.h"

// for tracking information
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"

#include "TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"
#include "TrackingTools/TrajectoryParametrization/interface/PerigeeTrajectoryParameters.h"

// to work with this: GlobalVector bvec = FieldB->inTesla(PointRef);
// does not exist in CMSSW_10_X
//#include "RecoVZero/VZeroFinding/interface/VZeroFinder.h"

// header files needed when u comment RecoVZero
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

// to work with VXlite: 
#include "NanoAOD/NanoAnalyzer/interface/vxlite/DAFVertexFinder.hh"
#include "NanoAOD/NanoAnalyzer/interface/vxlite/LinearizedTrack.hh"
#include "NanoAOD/NanoAnalyzer/interface/vxlite/VertexFitter.hh"



#include "NanoAOD/NanoAnalyzer/interface/LorentzVector.h"

/// CMS Fitter: 
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

// general:
#include <iostream>
#include <map>

   using namespace edm;
   using namespace reco;
   using namespace std;
   using namespace VxLite;
   using CLHEP::HepLorentzVector;

class HelixTrans
{

   public:
      	~HelixTrans();
	HelixTrans();
	HelixTrans(bool );

	double 	getRho(double,double,double);
	void 	PrintHelixParsZeus();
	void 	PrintCovarianceMatrixZEUS();
	void 	PrintFitResults();
	double  PrintdRtoVertex(vector <double>, double, double, double);
	void 	ResetFinder();

	double *getHelixParsCovTransform(vector<TransientTrack>::iterator);
	double *getHelixParsCov(vector<TransientTrack>::iterator);
	
	double *getHelixParsTransform(vector<TransientTrack>::iterator);
	double *getHelixPars(vector<TransientTrack>::iterator);

	double *FitVertexXYZ();
	void 	SetTrackCollectionToFit(vector<TransientTrack>);	
	void 	SetBuilder(ESHandle<TransientTrackBuilder>);
	void 	SetMagneticField(ESHandle<MagneticField>);
	double 	GetBz(const reco::TrackRef);



	/// + algorithm: 
	void 				SetGeneralTrackCollection(vector<TransientTrack>);		/// set all track collection
	void 				SetTrackCollectionForAlgorithm(vector<TransientTrack>);		/// set tracks from PV only	
	vector<TransientTrack> 		FindSortedHighPtTracks();					/// find sorted HighPt track collection; 
	vector<TransientTrack> 		FindSortedHighPtTracksSeparated();
	double 				getdRTracks(vector<TransientTrack>::iterator , vector<TransientTrack>::iterator );
	void 				FindTracksCluster();
	void 				CreateTrackClusters();
	vector <double>			*MakeAFitwithVxlite(reco::VertexCollection::const_iterator);
	
	vector<double> GetDLCov(reco::VertexCollection::const_iterator);
	
	void SetPVcollection(Handle<reco::VertexCollection>);
	vector<TransientTrack> FindTracksFromPV();
	vector<TransientTrack> GetTracksFromPV(reco::VertexCollection::const_iterator);
		
	///
	double GetSimilarity(double *, double *, ROOT::Math::SMatrix<double,3> );
	void SetBeamSpotInformation(Handle<reco::BeamSpot>);
// 	void 				ExtractDL(int, reco::VertexCollection::const_iterator);
	vector <double>			ExtractDLVxlite(vector<TransientTrack> , reco::VertexCollection::const_iterator , const VxLite::VertexFitter* , bool);
	void 				ExtractDLwithCov(int, reco::VertexCollection::const_iterator);
	vector <double>			ExtractDLCMSFitter(vector<TransientTrack> , reco::VertexCollection::const_iterator, TransientVertex, bool);
	//bool SortPair(const pair <double, double> &, const pair <double, double> &);


	/// CMS Fitter
	vector <double> 		*FitWithFitterCMS(reco::VertexCollection::const_iterator);


	/// 2 tracks vertex id finder
	void FindPVfor2Tracks(reco::VertexCollection::const_iterator, reco::VertexCollection::const_iterator,  vector<TransientTrack> );
	int  FindPVfor2TracksGetNvertexTracks();
	int  FindPVfor2TracksGetNvertex();	
	VertexCollection::const_iterator FindPVfor2TracksGetVertexIte();
	vector <double> SimpleFitWithFitterCMS(vector<TransientTrack> , reco::VertexCollection::const_iterator, bool );
	vector <double> SimpleFitWithVxlite(vector<TransientTrack> , reco::VertexCollection::const_iterator, bool );

	/// adds : 
	vector <int> GetIdTracks(vector<TransientTrack>);
	void SetNhighPtTracks(int);
	vector <TransientTrack> * GetClusterTracks();
	int GetNclusters();	

   private:
	/// general vars: 
	ESHandle<TransientTrackBuilder> theTrkBuilder;
	vector<reco::TransientTrack> tracksToFit;
	DAFVertexFinder finder;
	double HelixPars[5];
	double HelixParsZeus[5];	
	double VertexPos[4];
	int track_counter;
	bool debuger;

	/// For rho:
	double rho;

	/// For Bz: 
	const MagneticField *B;
	GlobalVector bvec;

	/// For cov matrix; 
	ROOT::Math::SMatrix<double,5> helix_cov_M;
	ROOT::Math::SMatrix<double,5> helix_Jacoby;
	ROOT::Math::SMatrix<double,5> helix_cov_M_zeus;
	double helix_par_cov[15];
	
	/// For vxlite: 
// 	float chi2_cut = 200;
// 	DAFVertexFinder finder(chi2_cut);
	const VxLite::VertexFitter* fitter;
	CLHEP::HepSymMatrix cov_mat;
	float ZEUS_Pars[500][5];
	float ZEUS_Pars_Cov[500][15];

	/// for algorithm: 
	TransientTrack  transientTr;
	
	vector<reco::TransientTrack> GeneralTrackCol;
	vector<reco::TransientTrack> TrackColForAlgorithm;
	vector<TransientTrack> TracksCollection[10];
	vector<TransientTrack> TracksCollectionClean[10];
	vector<TransientTrack>::iterator  it_tmp;
	int track1;
	int track2;
	double dR;	
	double dRtmp;
	bool IsoCheck;

	///
	Handle<reco::VertexCollection> PVertices;
	vector<TransientTrack> tracksVertex;
	double *VertPos;
	reco::BeamSpot beamSpot;
	double numerator;
	double denominator;	
	double DLCov[3];
	vector <double> DLvector[4];

	/// CMS fitter
	AdaptiveVertexFitter  theCMSFitter;
	TransientVertex myVertexCMSFit;


	/// 2 tracks vertex id finder
	int N2foundedVtxtracks;
	int N2foundedVtx;

	reco::VertexCollection::const_iterator iteFor2tracks;
	reco::VertexCollection::const_iterator iteFor2trackstmp;


	/// adds: 
	bool GeneralTrackColcheck;
	bool theTrkBuildercheck;
	bool magfieldCheck;
	int  NhighPtTracks;
	int  NclusterCounter;



};
