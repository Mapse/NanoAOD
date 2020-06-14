#include "NanoAOD/NanoAnalyzer/include/HelixTransClass.h"

/// Simple constructor:
HelixTrans::HelixTrans(bool debugger )
{
    debuger = debugger;	
    GeneralTrackColcheck = false;
    magfieldCheck = false;
    magfieldCheck = false;	
    NhighPtTracks = 10;	
}

/// Destructor (clean now)
HelixTrans::~HelixTrans()
{
}

/// Set track collection: 
void HelixTrans::SetTrackCollectionToFit(vector<TransientTrack> TrTrs)	{
	tracksToFit = TrTrs;
}


/// Set Track Builder by User;
void HelixTrans::SetBuilder(ESHandle<TransientTrackBuilder> TrkBuilder)	{
	theTrkBuilder = TrkBuilder;
	theTrkBuildercheck = true;
		if (debuger) {cout<<"Builder was set here. "<<endl;}
}


/// Set Magnetic fiend by User;
void HelixTrans::SetMagneticField(ESHandle<MagneticField> magfield)	{
	B = &*magfield;	
	magfieldCheck = true;
		if (debuger) {cout<<"Magnetic Field was set here. "<<endl;}
}


/// Extract Bz for each trackRef;
double HelixTrans::GetBz(const reco::TrackRef trackRef)	{

/*
	/// to change reference point on the helix;
	GlobalPoint CScenter(0., 0., 0.);
	transientTr = theTrkBuilder->build(trackRef);
        TrajectoryStateClosestToPoint  traj = transientTr.trajectoryStateClosestToPoint(CScenter);
	
	/// Set a Point to Extract magnetic field value:	
 	GlobalPoint PointRef(traj.position().x(),traj.position().y(),traj.position().z());
*/
	GlobalPoint PointRef(trackRef->vx(), trackRef->vy(), trackRef->vz());
	bvec = B->inTesla(PointRef);
// 	cout<<"Bz = "<<bvec.z()<<endl;
	return bvec.z();
}


/// Set Number of highPt Tracks
void HelixTrans::SetNhighPtTracks(int NhighPtTracks_)	{
	NhighPtTracks = NhighPtTracks_;
}


/// Rho [cm] calculation (q*R), q - charge, R- helix radius: 
double HelixTrans::getRho(double par0, double par1, double Bz)	{
	rho = cos(par1)/(0.0029979*Bz*par0);
	return rho;
}


/// return Helix parameters for reference point on the helix;
double *HelixTrans::getHelixPars(vector<TransientTrack>::iterator TrTrksIter){


	const reco::TrackRef trackRef = (TrTrksIter->trackBaseRef()).castTo<reco::TrackRef>();
	TransientTrack  transientTr = theTrkBuilder->build(trackRef);

	enum index { i_transverseCurvature = 0 , i_theta, i_phi0, i_d0, i_dz };	

        HelixPars[0] = trackRef->parameter(i_transverseCurvature);
        HelixPars[1] = trackRef->parameter(i_theta);
        HelixPars[2] = trackRef->parameter(i_phi0);
        HelixPars[3] = -1*trackRef->parameter(i_d0);
        HelixPars[4] = trackRef->parameter(i_dz);

	return HelixTrans::HelixPars;
}


/// transform Helix CMS format(parameters) into the Helix ZEUS format: 
double *HelixTrans::getHelixParsTransform(vector<TransientTrack>::iterator TrTrksIter)	{

	const reco::TrackRef trackRef = (TrTrksIter->trackBaseRef()).castTo<reco::TrackRef>();
	rho = getRho(HelixPars[0],HelixPars[1],GetBz(trackRef));

	HelixParsZeus[0] = HelixPars[2];								/// phi
	if (HelixParsZeus[0] < 0)	{HelixParsZeus[0] = 2*M_PI-fabs(HelixParsZeus[0]);}
	HelixParsZeus[1] = 1./rho; 									/// q/rho
	HelixParsZeus[2] = HelixPars[3]*HelixPars[0]/fabs(HelixPars[0]);				/// q*Dh	
	HelixParsZeus[3] = HelixPars[4]/cos(HelixPars[1]);						/// zh
	HelixParsZeus[4] = tan(HelixPars[1]);								/// cotan(theta)	

	return HelixTrans::HelixParsZeus;
}


/// transform Helix Covariance matrix from CMS format into Helix ZEUS format: 
double *HelixTrans::getHelixParsCovTransform(vector<TransientTrack>::iterator TrTrksIter){
	
	reco::TrackBase::CovarianceMatrix covma = (TrTrksIter->track()).covariance();
	for (int iM = 0; iM< 5; iM++)	{
		for (int jM = 0; jM<5; jM++)	{
			helix_cov_M[iM][jM] = covma(iM,jM);
			helix_Jacoby[iM][jM] = 0.;
						}
					}
	/// Set elements which are not == 0;
	helix_Jacoby[0][2] = 1.;
	helix_Jacoby[1][0] = HelixParsZeus[1]/HelixPars[0]; 
	helix_Jacoby[1][1]=  HelixParsZeus[1]*HelixParsZeus[4];
	helix_Jacoby[2][3] = HelixPars[0]/fabs(HelixPars[0]);
	helix_Jacoby[3][1] = HelixParsZeus[3]*HelixParsZeus[4];
	helix_Jacoby[3][4] = 1./cos(HelixPars[1]);
	helix_Jacoby[4][1] = 1./(cos(HelixPars[1])*cos(HelixPars[1]));

	/// Transformed M_Zeus = J*M_CMS*Jt
	helix_cov_M_zeus = helix_Jacoby* helix_cov_M * ROOT::Math::Transpose(helix_Jacoby);

	/// Making proper format for VXlite input: (15 elements in the array)
        helix_par_cov[0]  = helix_cov_M_zeus[0][0]; helix_par_cov[1]  = helix_cov_M_zeus[0][1]; helix_par_cov[2]  = helix_cov_M_zeus[0][2];
        helix_par_cov[3]  = helix_cov_M_zeus[0][3]; helix_par_cov[4]  = helix_cov_M_zeus[0][4]; helix_par_cov[5]  = helix_cov_M_zeus[1][1];
        helix_par_cov[6]  = helix_cov_M_zeus[1][2]; helix_par_cov[7]  = helix_cov_M_zeus[1][3]; helix_par_cov[8]  = helix_cov_M_zeus[1][4];
        helix_par_cov[9]  = helix_cov_M_zeus[2][2]; helix_par_cov[10] = helix_cov_M_zeus[2][3]; helix_par_cov[11] = helix_cov_M_zeus[2][4];
        helix_par_cov[12] = helix_cov_M_zeus[3][3]; helix_par_cov[13] = helix_cov_M_zeus[3][4]; helix_par_cov[14] = helix_cov_M_zeus[4][4];


	return helix_par_cov;
}


/// Print 5 helix parameters in ZEUS format: 
void HelixTrans::PrintHelixParsZeus()	{
	cout<<"HelixParsZeus : "<<endl;
// 	int trkcount = 0;
// 	for( vector<TransientTrack>::iterator  TrTrksIter = tracksToFit.begin(); TrTrksIter !=tracksToFit.end(); TrTrksIter++)
//         {	
// 	cout<<"Track #"<<trkcount<<"	";
// 	getHelixParsTransform(TrTrksIter);
	cout<<HelixParsZeus[0]<<"	"<<HelixParsZeus[1]<<"	"<<HelixParsZeus[2]<<"	"<<HelixParsZeus[3]<<"	"<<HelixParsZeus[4]<<endl;
// 	trkcount++;
// 	}
}


/// Print covariance matrix from ZEUS:
void HelixTrans::PrintCovarianceMatrixZEUS()	{
	
	cout<<"Covariance matrix in ZEUS format. Obtained from CMS covariance matrix with formula: M_Zeus = J*M_CMS*Jt"<<endl;
	for (int i=0; i<5; i++) {

		for (int j=0; j<5; j++)	{
			cout<<helix_cov_M_zeus[i][j]<<"	";
					} cout<<endl;	
				}
}


/// Print ZEUS Fit results: 
void HelixTrans::PrintFitResults()	{

		cout<<"INFO: Refitted vertex, xyz coordinates: "<<endl;
		cout<<"X: "<<VertexPos[0]<<" +- "<<sqrt(cov_mat[0][0])<<endl;
		cout<<"Y: "<<VertexPos[1]<<" +- "<<sqrt(cov_mat[1][1])<<endl;
		cout<<"Z: "<<VertexPos[2]<<" +- "<<sqrt(cov_mat[2][2])<<endl;
		cout<<"Ntracks in vertex = "<<tracksToFit.size()<<endl;
		cout<<"Chi2/ndof = "<<VertexPos[3]<<endl;
		cout<<"Chi2 = "<<fitter->chi2()<<endl;
		cout<<"Ndof = "<<fitter->ndof()<<endl;

}


/// Distanse to the defined vertex from the Refitted one;
double HelixTrans::PrintdRtoVertex(vector <double> VtxPars, double x, double y, double z)	{
	double dR= 0.;
	dR = dR + (VtxPars[3]-x)*(VtxPars[3]-x);
	dR = dR + (VtxPars[4]-y)*(VtxPars[4]-y);
	dR = dR + (VtxPars[5]-z)*(VtxPars[5]-z);
	// cout<<"Distance to chosen vertex with coordinates: "<<x<<"	"<<y<<"	"<<z<<"  = "<<sqrt(dR)<<endl;
	return sqrt( dR );
}


/// Fit Vertex
double *HelixTrans::FitVertexXYZ()	{
	
// 	finder(0.);

	
	track_counter = 0;
	for( vector<TransientTrack>::iterator  TrTrksIter = tracksToFit.begin(); TrTrksIter !=tracksToFit.end(); TrTrksIter++)
        {	
		HelixTrans::getHelixPars(TrTrksIter);
// 			cout<<"	Helix parameters from CMS data obtained;"<<endl;
		HelixTrans::getHelixParsTransform(TrTrksIter);
// 			cout<<"	Helix parameters from CMS transformed into ZEUS format ;"<<endl;
		HelixTrans::getHelixParsCovTransform(TrTrksIter);
// 			cout<<"	Helix covariance matrix from CMS transformed into ZEUS format ;"<<endl;

// 		cout<<HelixParsZeus[0]<<"	"<<HelixParsZeus[1]<<"	"<<HelixParsZeus[2]<<"	"<<HelixParsZeus[3]<<"	"<<HelixParsZeus[4]<<endl;
// 		PrintCovarianceMatrixZEUS();

		for (int k=0;k<5;k++)	{ZEUS_Pars[track_counter][k] = HelixParsZeus[k];}
		for (int k=0;k<15;k++)	{ZEUS_Pars_Cov[track_counter][k] = helix_par_cov[k];}

		finder.addTrack(new LinearizedTrack(track_counter+1, ZEUS_Pars[track_counter], ZEUS_Pars_Cov[track_counter]));

		track_counter++;
	
        } 

	if (track_counter > 1)	{
		if (debuger) {cout<<"	Fitter starts here;"<<endl;}
	finder.findVertex();
	fitter = finder.fitter();
	cov_mat = fitter->covariance();

		if (debuger) {cout<<"	Vertex parameters obtained here; "<<endl;}
	VertexPos[0] = fitter->vertex()[0];
	VertexPos[1] = fitter->vertex()[1];
	VertexPos[2] = fitter->vertex()[2];
	VertexPos[3] = fitter->chi2()/fitter->ndof();
	
//	cout<<fitter->vertex()[0]<<"	"<<fitter->vertex()[1]<<"	"<<fitter->vertex()[2]<<endl;
// 	PrintFitResults();
// 	finder.reset();
// 	cout<<"	N tracks = "<<track_counter<<endl;
	
				}	else { cout<<" Size of track collection < 2 "<<endl; exit(0);}
	return  HelixTrans::VertexPos;	
}


/// Reset finder;
void HelixTrans::ResetFinder()	{
	finder.reset();
}


/// test method fit with standart fitter from CMS: 
vector <double> *HelixTrans::FitWithFitterCMS(reco::VertexCollection::const_iterator ite)	{

	CreateTrackClusters();
	vector<TransientTrack> tracksToFitCMSVtx;


	
//	double DLXYCMSFit; 
// 	double dDLXYCMSFit;
//	double DLx; double dDLx;
//	double DLy; double dDLy;
	ROOT::Math::SMatrix<double,3> DataVertexM;

	for( int i=0; i<10;i++)
        {	

		cout<<"Tracks in Cluster "<<i<<" = "<<TracksCollectionClean[i].size()<<endl;
//		DLx = 0.; DLy = 0.; DLXYCMSFit = 0.; dDLx = 0.; dDLy = 0.;


		if (TracksCollectionClean[i].size() < 2) {continue;}

	if (debuger)	{
		cout<<"*****************************"<<endl;
		cout<<"N tracks in the new track collection: "<<TracksCollectionClean[i].size()<<endl;
			}	
		tracksToFitCMSVtx = TracksCollectionClean[i];

		myVertexCMSFit = theCMSFitter.vertex(tracksToFitCMSVtx/*, beamSpot*/);
		if (myVertexCMSFit.isValid())	{
		  
	if (debuger)	{		  
		cout<<"Ntrack = "<<tracksToFitCMSVtx.size()<<endl;
		cout<<"X = "<<myVertexCMSFit.position().x()<<" +- "<<sqrt(myVertexCMSFit.positionError().cxx())<<endl;
		cout<<"Y = "<<myVertexCMSFit.position().y()<<" +- "<<sqrt(myVertexCMSFit.positionError().cyy())<<endl;
		cout<<"Z = "<<myVertexCMSFit.position().z()<<" +- "<<sqrt(myVertexCMSFit.positionError().czz())<<endl;
		cout<<"Chi2 = "<<myVertexCMSFit.normalisedChiSquared()*myVertexCMSFit.degreesOfFreedom()<<endl;
	}

	if (debuger)	{	

		vector<TransientTrack>::iterator trk_tst = TracksCollectionClean[i].begin();
		cout<<"Leading track : Px, Py:  "<<(trk_tst->track()).px()<<"	"<<(trk_tst->track()).py()<<endl;

		/// My test for debug here:
		cout<<"Cluster tracks: "<<endl<<endl;
		int ccount = 0;
		for( vector<TransientTrack>::iterator  ittm1 = GeneralTrackCol.begin(); ittm1 !=GeneralTrackCol.end(); ittm1++)
		{  
		   const reco::TrackRef ref1 = (ittm1->trackBaseRef()).castTo<reco::TrackRef>();
		   for (vector<TransientTrack>::iterator inpa = TracksCollectionClean[i].begin(); inpa != TracksCollectionClean[i].end(); inpa++)	
		   {
			const reco::TrackRef ref2 = (inpa->trackBaseRef()).castTo<reco::TrackRef>();
			if (ref1.get() != ref2.get())	{continue;}
			
			double Drr = ((ittm1->track()).eta() -(trk_tst->track()).eta())*((ittm1->track()).eta()-(trk_tst->track()).eta());
			Drr += ((ittm1->track()).phi()-(trk_tst->track()).phi())*((ittm1->track()).phi()-(trk_tst->track()).phi());
			Drr = sqrt(Drr);
			
// 			cout<<"Phi = "<<(inpa->track()).phi()<<endl;
			cout<<"N = "<<ccount<<"  Pt = "<<(inpa->track()).pt()<<"	Phi = "<<(inpa->track()).phi()<<"	Eta = "<<(inpa->track()).eta()<<"	dR = "<<Drr<<endl;
		
		   }
		ccount++;

		}
			}

// 		double Rerrx = sqrt(myVertexCMSFit.positionError().cxx())/myVertexCMSFit.position().x();
// 		double Rerry = sqrt(myVertexCMSFit.positionError().cyy())/myVertexCMSFit.position().y();
// 		double Rerrz = sqrt(myVertexCMSFit.positionError().czz())/myVertexCMSFit.position().z();


// 		if (Rerrx > 0.3 || Rerry > 0.3 || Rerrz > 0.3)	{continue;}

		/// Extract DL and dDL;
		/// ExtractDLCMSFitter(TracksCollectionClean[i], ite);		
		
		/// check for valid Error matrix, in composition with DL in could extract sqrt from a negative values;
// 		if ( isnan(DLCov[0]) || isnan(DLCov[1])) {continue;}		
		
		/// this method will return:   
		DLvector[0].push_back(DLCov[0]);	/// DL
		DLvector[1].push_back(DLCov[1]);	/// DL error
// 		DLvector[2].push_back(DLCov[2]);	/// mass
		DLvector[3].push_back(myVertexCMSFit.normalisedChiSquared());  
		
		} /// check if vertex is valid

	}

	


	return  DLvector;	

}


vector <double> HelixTrans::SimpleFitWithFitterCMS(vector<TransientTrack> tracksVector, reco::VertexCollection::const_iterator ite, bool withbeam)       {

                vector <double> OutPutVector;
                myVertexCMSFit = theCMSFitter.vertex(tracksVector/*, beamSpot*/);
                if (myVertexCMSFit.isValid())   {

        if (debuger)    {
                cout<<"Ntrack = "<<tracksVector.size()<<endl;
                cout<<"X = "<<myVertexCMSFit.position().x()<<" +- "<<sqrt(myVertexCMSFit.positionError().cxx())<<endl;
                cout<<"Y = "<<myVertexCMSFit.position().y()<<" +- "<<sqrt(myVertexCMSFit.positionError().cyy())<<endl;
                cout<<"Z = "<<myVertexCMSFit.position().z()<<" +- "<<sqrt(myVertexCMSFit.positionError().czz())<<endl;
                cout<<"Chi2 = "<<myVertexCMSFit.normalisedChiSquared()*myVertexCMSFit.degreesOfFreedom()<<endl;
        		}

		OutPutVector = ExtractDLCMSFitter(tracksVector, ite, myVertexCMSFit, withbeam);		/// 0 - DL, 1 - DLsigma, 2 - DL sign;
		OutPutVector.push_back(myVertexCMSFit.position().x());					/// 3 - vtx pos x	
		OutPutVector.push_back(myVertexCMSFit.position().y());					/// 4 - vtx pos y
		OutPutVector.push_back(myVertexCMSFit.position().z());					/// 5 - vtx pos z	
		OutPutVector.push_back(sqrt(myVertexCMSFit.positionError().cxx()));			/// 6 - vtx pos err x
                OutPutVector.push_back(sqrt(myVertexCMSFit.positionError().cyy()));			/// 7 - vtx pos err y
                OutPutVector.push_back(sqrt(myVertexCMSFit.positionError().czz()));			/// 8 - vtx pos err z	
		OutPutVector.push_back(myVertexCMSFit.normalisedChiSquared());				/// 9 - vtx normalized Chi2
		OutPutVector.push_back(myVertexCMSFit.degreesOfFreedom());				/// 10 - Ndof	


		/// mass reconstruction:
		/*
		/// ******************************************************************************
		vector<TransientTrack>::iterator  ittmp = tracksVector.begin();

		HepLorentzVector p(0,0,0,0);
    		for(unsigned i = 0 ; i < tracksVector.size(); i++) 
		{
        	const LinearizedTrack* t = myVertexCMSFit.track(i);
	
		//cout<<(ittmp->track()).p()<<endl;
        	HepLorentzVector tp = t->momentum( (ittmp->track()).p() );
        	p +=  t->weight() * tp;
		ittmp++;
    		}

 		double mass = p.m();
		cout<<"mass = "<<mass<<endl;	
		//
		

		typedef std::map<reco::TransientTrack, float> TransientTrackToFloatMap;
		TransientTrackToFloatMap weightMAP = myVertexCMSFit.weightMap();
		

		for( vector<TransientTrack>::iterator  trIt = tracksVector.begin(); trIt != tracksVector.end(); trIt++)
		{
			TransientTrack  transientTMPtr;
			const reco::TrackRef trItref = (trIt->trackBaseRef()).castTo<reco::TrackRef>();
			transientTMPtr = theTrkBuilder->build(trItref);
			float trackWeightRef = myVertexCMSFit.trackWeight(transientTMPtr);
			cout<<"Weight = "<<trackWeightRef<<endl;
			
// 			myVertexCMSFit.vertex().p4();

// 			momentum()
			
// 			cout<<weightMAP[transientTMPtr]<<endl;
			
	
// 			cout<<"refited ? "<<myVertexCMSFit.refittedTrack(transientTMPtr)<<endl;
		}



		

// 		cout<<myVertexCMSFit.vertexState().weightInMixture()<<endl;
//  		cout<<"Weights? "<<myVertexCMSFit.hasTrackWeight()<<endl;
		*/
		/// ******************************************************************************
		//OutPutVector.push_back(mass);								/// 11 - vertex mass here
		OutPutVector.push_back(1);								/// 11 - vertex valid 1/0
						}
		else {OutPutVector.push_back(0);}
		

		return OutPutVector;

}



vector<double> HelixTrans::ExtractDLCMSFitter(vector<TransientTrack> tracksVector, reco::VertexCollection::const_iterator ite, TransientVertex CMSFittedVtx, bool withbeam)	{

  	vector <double> ExtractDLCMSFitterOutput;  ExtractDLCMSFitterOutput.clear();
	double Lx = 0.;
	double Ly = 0.;
	double Lz = CMSFittedVtx.position().z() - ite->z();

	/// extract position of the PV from data
	double PVposX = ite->x();
	double PVposY = ite->y();

	/// make a sum of all track momentum in track collection: 
	double PxSum = 0.;
	double PySum = 0.;
	double PtSum = 0.;

	for (vector<TransientTrack>::iterator itSumP = tracksVector.begin(); itSumP != tracksVector.end(); itSumP++)	
	{
		PxSum = PxSum + (itSumP->track()).px();
		PySum = PySum + (itSumP->track()).py();
// 		cout<<"Px"<<PxSum<<endl;
	}
	
	PtSum = sqrt(PxSum*PxSum + PySum*PySum);

	if (!withbeam)
	{
	/// calculate DL projection components on XY plane;  
	Lx = CMSFittedVtx.position().x() - PVposX;
	Ly = CMSFittedVtx.position().y() - PVposY;
	}	else 
	{
 	Lx = (CMSFittedVtx.position().x() - (beamSpot.x0() + beamSpot.dxdz()*Lz) );
	Ly = (CMSFittedVtx.position().y() - (beamSpot.y0() + beamSpot.dydz()*Lz) );	
										
	}

	double ProjectedDecayLength  =0;
	double ProjectedDecayLengthError = 0.;	
	
	/// calculate DL xy : 
	ProjectedDecayLength = (Lx*PxSum + Ly*PySum)/PtSum;

	ProjectedDecayLength = ProjectedDecayLength /*+ 0.001*/;
	// cout<<"DL xy proj = "<<ProjectedDecayLength<<endl;

	typedef ROOT::Math::SMatrix<double,3>       SMatrix33;

	SMatrix33 DataVertexM;	
	
	if (!withbeam)
	{
	DataVertexM = ite->covariance();
	} else
	{
	DataVertexM(0,0) = beamSpot.covariance(0, 0);
	DataVertexM(0,1) = beamSpot.covariance(0, 1);	
	DataVertexM(0,2) = 0;
	DataVertexM(1,0) = beamSpot.covariance(1, 0);
	DataVertexM(1,1) = beamSpot.covariance(1, 1);
	DataVertexM(1,2) = 0;
	DataVertexM(2,0) = 0;
	DataVertexM(2,1) = 0;
	DataVertexM(2,2) = 0;
	}



/// error DL calculation in matrices:
/// *********************************  *********************************
	SMatrix33 SumVertexCovM;

	double MyVertexError[3][3];

	MyVertexError[0][0] = CMSFittedVtx.positionError().cxx();
	MyVertexError[0][1] = CMSFittedVtx.positionError().cyx();	
	MyVertexError[0][2] = 0;
	MyVertexError[1][0] = CMSFittedVtx.positionError().cyx();
	MyVertexError[1][1] = CMSFittedVtx.positionError().cyy();
	MyVertexError[1][2] = 0;
	MyVertexError[2][0] = 0;
	MyVertexError[2][1] = 0;
	MyVertexError[2][2] = 0;


	for (int i =0; i<3; i++)	{
		for (int j=0; j<3; j++)		{	

		SumVertexCovM[i][j] = MyVertexError[i][j] + DataVertexM[i][j];
						}	
					}


	/// Extract DLxy sigma;
	double DLSigma = 0.;
	double PClusterV[3] = {PxSum, PySum, 0.};					/// vector components of the DL 3D

	DLSigma = GetSimilarity(PClusterV, PClusterV, SumVertexCovM);
	DLSigma = DLSigma/(PtSum*PtSum);
	DLSigma = sqrt(DLSigma);
	ProjectedDecayLengthError = DLSigma;

/// *********************************  *********************************

	ExtractDLCMSFitterOutput.push_back(ProjectedDecayLength);
	ExtractDLCMSFitterOutput.push_back(ProjectedDecayLengthError);
	ExtractDLCMSFitterOutput.push_back(ProjectedDecayLength/ProjectedDecayLengthError);
		
	/// couts: 
	if (debuger)	{
// 	cout<<"Ptima Vtx  = "<<ite->x()<<"	VtxY = "<<ite->y()<<"	VtxZ = "<<ite->z()<<endl;
	cout<<"Refit Vtx = "<<CMSFittedVtx.position().x()<<"	VtxY = "<<CMSFittedVtx.position().y()<<"	VtxZ = "<<CMSFittedVtx.position().z()<<endl;
	cout<<"DL = "<<ProjectedDecayLength<<" DL error = "<<ProjectedDecayLengthError<<endl;

			}

	return ExtractDLCMSFitterOutput;

}


/// Methods for track sorting algorithm 
/// **********************************************************************************************************

/// Set general track collection:
/// result : 		set Global var GeneralTrackCol of Transient Track Collection;
/// input: 		external Track collection from data sample;
/// use global vars : 	GeneralTrackCol; 
void HelixTrans::SetGeneralTrackCollection(vector<TransientTrack> GenTrkCol)	{
	GeneralTrackCol = GenTrkCol;
	GeneralTrackColcheck = true;
}


/// Set track collection from PV: 
/// result : 		set Global var "TrackColForAlgorithm" of Transient Track Collection;
/// input: 		external Track collection from data sample for one Primary Vertex;
/// use global vars : 	TrackColForAlgorithm; 
void HelixTrans::SetTrackCollectionForAlgorithm(vector<TransientTrack> trackCol)	{
	TrackColForAlgorithm = trackCol;
}


/// this method return sorted vector of Transient Tracks by its greates Pt:
/// result : 		return Transient vector of HighPt tracks;
/// input: 		general var "TrackColForAlgorithm";
/// use global vars : 	NO; 
vector<TransientTrack> HelixTrans::FindSortedHighPtTracks()	{

	vector<TransientTrack> mytracks;
	vector<TransientTrack> tracks = TrackColForAlgorithm;
	
	//vector<pair<double, const reco::TrackRef> > VectPtTransient;
	/// create a map which sort instantly content by greater values; 
	typedef std::multimap<double, TransientTrack, greater<double> > mimap;
	typedef std::pair<double,TransientTrack> ip;
	mimap VectPtTransient;

	//cout<<"Ntracks = "<<tracksToFit.size()<<endl;

	/// tmp var
	TransientTrack  transientTr;	
	
	/// loop over sorted vector and fill new vector <Transient>	
	if (debuger) {cout<<"Vector of pairs <double, Transient> fill here: "<<endl;  }
	for( vector<TransientTrack>::iterator  TrTrksIter = tracks.begin(); TrTrksIter !=tracks.end(); TrTrksIter++)
        {
		const reco::TrackRef trackRef = (TrTrksIter->trackBaseRef()).castTo<reco::TrackRef>();
		transientTr = theTrkBuilder->build(trackRef);
		
		//cout<<"	Not sorted Pt: "<<(TrTrksIter->track()).pt()<<endl;
		VectPtTransient.insert(ip( (TrTrksIter->track()).pt(), transientTr ) );
				
	}


	if (debuger) {cout<<"Vector sort start here: "<<endl;}

	std::multimap<double, TransientTrack, greater<double> >::iterator it;
	for ( it=VectPtTransient.begin(); it!=VectPtTransient.end(); ++it)
    	{
	// cout << (*it).first<<endl;
	mytracks.push_back((*it).second);	
    	}

	/// example how to use new sorten Transient vector (as usual Transient vector...):
	/*
	for( vector<TransientTrack>::iterator  TrTrksIter2 = mytracks.begin(); TrTrksIter2 !=mytracks.end(); TrTrksIter2++)
        {
		const reco::TrackRef trackRef2 = (TrTrksIter2->trackBaseRef()).castTo<reco::TrackRef>();
		cout<<"	Sorted Pt: "<<(TrTrksIter2->track()).pt()<<endl;
	}
	*/

	/// clean vector here;
	VectPtTransient.clear();

	return mytracks;

}


/// this method return dR = sqrt(dphi*dphi+deta*deta) for two tracks use its iterator from the external loop
/// result : 		return dR of two tracks use it iterators;
/// input: 		two iterators for two tracks;
/// use global vars : 	NO; 
double HelixTrans::getdRTracks(vector<TransientTrack>::iterator it1, vector<TransientTrack>::iterator it2 )	{

	const reco::TrackRef trackRef_it1 = (it1->trackBaseRef()).castTo<reco::TrackRef>();
      	const reco::TrackRef trackRef_it2 = (it2->trackBaseRef()).castTo<reco::TrackRef>();
		
	double dR = 999.;

	if (trackRef_it1 != trackRef_it2)	{
		dR =  ( ((it1->track()).phi() - (it2->track()).phi())*((it1->track()).phi() - (it2->track()).phi())  );
		dR += ( ((it1->track()).eta() - (it2->track()).eta())*((it1->track()).eta() - (it2->track()).eta())  );

	} else dR = 0.0001;

	return sqrt(dR);
}


/// Make vector of Transient tracks which are separated from each other :
/// Uses FindSortedHighPtTracks() method to get input track collection;
/// result : 		return Transient vector of HighPt separated tracks;
/// input: 		tracks from the method FindSortedHighPtTracks();
/// use global vars : 	NO; 
vector<TransientTrack> HelixTrans::FindSortedHighPtTracksSeparated()	{


	vector<TransientTrack> mytracksSorted =  FindSortedHighPtTracks() ;
	vector<TransientTrack> mytracksSortedIsolated;
	
	track1 = 0; track2 = 0; dR = 0.;
	IsoCheck = true;


	for( vector<TransientTrack>::iterator  TrTrksIter_it1 = mytracksSorted.begin(); TrTrksIter_it1 !=mytracksSorted.end(); TrTrksIter_it1++)
        {
		if ( (int) mytracksSortedIsolated.size() >= NhighPtTracks)	{continue;}


			const reco::TrackRef trackRef_it1 = (TrTrksIter_it1->trackBaseRef()).castTo<reco::TrackRef>();
			transientTr = theTrkBuilder->build(trackRef_it1);

		if (track1 == 0) { mytracksSortedIsolated.push_back(transientTr); }

		vector<TransientTrack> mytracksSortedIsolatedTmp = mytracksSortedIsolated;
// 		cout<<"loop 1, track "<<track1<<endl;
		track2 = 0;	

		IsoCheck = true;

		if ((TrTrksIter_it1->track()).pt() < 0.5 ) {continue;}

		for( vector<TransientTrack>::iterator  TrTrksIter_it2 = mytracksSortedIsolatedTmp.begin(); TrTrksIter_it2 !=mytracksSortedIsolatedTmp.end(); TrTrksIter_it2++)
        	{

				
// 				cout<<"loop 2, track "<<track2<<endl;
				dR = getdRTracks(TrTrksIter_it1, TrTrksIter_it2);
// 				cout<<(TrTrksIter_it1->track()).pt()<<"	"<<(TrTrksIter_it2->track()).pt()<<"	"<<dR<<endl;
				if (dR < 1 )	{ IsoCheck = false;}	
				track2++;	
		}
		
		if (IsoCheck) {mytracksSortedIsolated.push_back(transientTr);}
// 		cout<<"mytracks size here: "<<mytracksSortedIsolated.size()<<endl;


		track1++ ;
	}

	return mytracksSortedIsolated;

}


/// find tracks which belongs to the isolated high Pt tracks;
/// result : 		fill the global var TracksCollection[10][vector Transient tracks];
/// input: 		tracks from the method FindSortedHighPtTracksSeparated();
/// use global vars : 	fill TracksCollection[][]; 
void HelixTrans::FindTracksCluster(){

	int track11 = 0;
	int track22 = 0;
	dR = 0.;	
	dRtmp = 1000.;
	int jetCounter = 0;

	vector<TransientTrack> fTracks =  FindSortedHighPtTracksSeparated() ;
// 	TransientTrack transientTrTmp;

	vector<TransientTrack> TracksCollectionTmp;
	
	int zz  = 0;
	
		for( vector<TransientTrack>::iterator  iterator1 = fTracks.begin(); iterator1 !=fTracks.end(); iterator1++)
        	{
			const reco::TrackRef trackRefIT = (iterator1->trackBaseRef()).castTo<reco::TrackRef>();
// 			cout<<(iterator1->track()).pt()<<endl;
			TransientTrack transientTrTmp = theTrkBuilder->build(trackRefIT);
			TracksCollection[zz].push_back(transientTrTmp);
			zz++;
		}

// 	for(int i=0;i<10;i++)	{cout<<TracksCollection[i].size()<<endl;}

	for( vector<TransientTrack>::iterator  it1 = TrackColForAlgorithm.begin(); it1 !=TrackColForAlgorithm.end(); it1++)
        {


		track22 = 0; dRtmp = 1000; jetCounter = -999;
		if (fTracks.size() == 0)	{continue;}
		if ((it1->track()).pt() < 0.5 ) {continue;}

		for( vector<TransientTrack>::iterator  it2 = fTracks.begin(); it2 !=fTracks.end(); it2++)
        	{
			dR = getdRTracks(it1, it2);
// 			cout<<"dR = "<<dR<<endl;
			if (dR < 1. && dR < dRtmp && dR != 0.01)	{dRtmp = dR; jetCounter = track22;}
// 			cout<<"track "<<track11<<"	track "<<track22<<"	dR = "<<dR<<"	loop 2 track Pt =  "<<(it2->track()).pt()<<" JetCounter = "<<jetCounter<<endl;
			track22++;
		}
		if (dRtmp<1 /*|| track1 == 0*/ && jetCounter >= 0)	{

			const reco::TrackRef trackRef = (it1->trackBaseRef()).castTo<reco::TrackRef>();
			transientTr = theTrkBuilder->build(trackRef);
			TracksCollection[jetCounter].push_back(transientTr); 
					
						}	
		track11++;
	}
	
}


/// This method fill a new array "TracksCollectionClean" same to the "TracksCollection";
/// But now new array will have good order of the values without 0 entries (only at the end)
/// result : 		fill the global var TracksCollectionClean[10][vector Transient tracks];
/// input: 		tracks from the method FindTracksCluster();
/// use global vars : 	fill TracksCollectionClean[][]; 
void HelixTrans::CreateTrackClusters()	{
	
	FindTracksCluster();

	NclusterCounter = 0;	
	for( int i=0; i<10;i++)
        {
// 		cout<<" Size = "<<TracksCollection[i].size()<<endl;
		if ( TracksCollection[i].size() == 0 ) continue; 
			else 	{
				TracksCollectionClean[NclusterCounter] =  TracksCollection[i];
				NclusterCounter++;
				}
// 		cout<<" Size = "<<TracksCollectionClean[i].size()<<endl;
					
	}

}

vector <TransientTrack> * HelixTrans::GetClusterTracks()	{
	return TracksCollectionClean;
}

int HelixTrans::GetNclusters()	{
	return NclusterCounter;
}


/// GetIdTracks: 
vector <int> HelixTrans::GetIdTracks(vector<TransientTrack> inputVectorTracks)	{

	int TrackId = 0;
	vector <int> TrackIds;	

	if (GeneralTrackColcheck )	{

		for( vector<TransientTrack>::iterator  iter = GeneralTrackCol.begin(); iter != GeneralTrackCol.end(); iter++)
		{
			
		   const reco::TrackRef ref1 = (iter->trackBaseRef()).castTo<reco::TrackRef>();
		   for (vector<TransientTrack>::iterator iter2 = inputVectorTracks.begin(); iter2 != inputVectorTracks.end(); iter2++)	
		   {
			const reco::TrackRef ref2 = (iter2->trackBaseRef()).castTo<reco::TrackRef>();
			if (ref1.get() == ref2.get())	{TrackIds.push_back(TrackId); } else {continue;}
		
		   }
		   TrackId++;	
		}
					} else {cout<<"Warning: General Track Collection was not set."<<endl;}

	return TrackIds;
}


vector <double> HelixTrans::SimpleFitWithVxlite(vector<TransientTrack> tracksVector, reco::VertexCollection::const_iterator ite, bool withbeam)       {

                vector <double> OutPutVectorZEUS;
// 		SetTrackCollectionToFit(tracksVector);

// 	cout<<"Start fit procedure: "<<endl;
	track_counter = 0;
	for( vector<TransientTrack>::iterator  TrTrksIter = tracksVector.begin(); TrTrksIter !=tracksVector.end(); TrTrksIter++)
        {	
		/// important! cut on Ntracks in the PV;
		if (track_counter > 450)	{continue;}	
	
		HelixTrans::getHelixPars(TrTrksIter);
// 			cout<<"	Helix parameters from CMS data obtained;"<<endl;
		HelixTrans::getHelixParsTransform(TrTrksIter);
// 			cout<<"	Helix parameters from CMS transformed into ZEUS format ;"<<endl;
		HelixTrans::getHelixParsCovTransform(TrTrksIter);
// 			cout<<"	Helix covariance matrix from CMS transformed into ZEUS format ;"<<endl;

// 		cout<<HelixParsZeus[0]<<"	"<<HelixParsZeus[1]<<"	"<<HelixParsZeus[2]<<"	"<<HelixParsZeus[3]<<"	"<<HelixParsZeus[4]<<endl;
// 		PrintCovarianceMatrixZEUS();

		for (int k=0;k<5;k++)	{ZEUS_Pars[track_counter][k] = HelixParsZeus[k];}
		for (int k=0;k<15;k++)	{ZEUS_Pars_Cov[track_counter][k] = helix_par_cov[k];}

// 		cout<<"Px trans = "<<(TrTrksIter->track()).px()<<"	"<<(TrTrksIter->track()).p()<<endl;

		finder.addTrack(new LinearizedTrack(track_counter+1, ZEUS_Pars[track_counter], ZEUS_Pars_Cov[track_counter]));

		track_counter++;
	
        } 

	if (track_counter > 1)	{
		if (debuger) {cout<<"	Fitter starts here;"<<endl;}
	finder.findVertex();
	fitter = finder.fitter();
	cov_mat = fitter->covariance();
	
//                 if (myVertexCMSFit.isValid())   {

        if (debuger)    {
                cout<<"Ntrack = "<<tracksVector.size()<<endl;
		cout<<"X: "<<fitter->vertex()[0]<<" +- "<<sqrt(cov_mat[0][0])<<endl;
		cout<<"Y: "<<fitter->vertex()[1]<<" +- "<<sqrt(cov_mat[1][1])<<endl;
		cout<<"Z: "<<fitter->vertex()[2]<<" +- "<<sqrt(cov_mat[2][2])<<endl;
                cout<<"Chi2 = "<<fitter->chi2()<<endl;
		cout<<"Ndof = "<<fitter->ndof()<<endl;
        		}


	/// extract mass of the vertex here:  ??? need explanations here! 
	vector<TransientTrack>::iterator  ittmp = tracksVector.begin();

	HepLorentzVector p(0,0,0,0);
    	for(unsigned i = 0 ; i < tracksVector.size(); i++) 
	{
        const LinearizedTrack* t = finder.track(i);
// 	cout<<"weight = "<<t->weight()<<endl;
        HepLorentzVector tp = t->momentum( (ittmp->track()).p() );
// 	cout<<"Before Lin P = "<<(ittmp->track()).p()<<"  after, P track = "<<tp.vect().mag()<<"	"<<sqrt(tp.x()*tp.x()+tp.y()*tp.y()+tp.z()*tp.z())<<"	"<<tp.t()<<endl;
        p = p +  t->weight() * tp;
	ittmp++;
    	}
    	double mass = p.m();

// 	cout<<"P = "<<sqrt(p.px()*p.px()+p.py()*p.py()+p.pz()*p.pz())<<endl;

// 	cout<<"Mass = "<<mass<<endl;
	if (debuger)	{ cout<<"Mass = "<<mass<<endl; }

// 		ExtractDL(i, ite);
		OutPutVectorZEUS = ExtractDLVxlite(tracksVector, ite, fitter, withbeam);			/// 0 - DL, 1 - DLsigma, 2 - DL sign;
		OutPutVectorZEUS.push_back(fitter->vertex()[0]);						/// 3 - vtx pos x	
		OutPutVectorZEUS.push_back(fitter->vertex()[1]);						/// 4 - vtx pos y
		OutPutVectorZEUS.push_back(fitter->vertex()[2]);						/// 5 - vtx pos z	
		OutPutVectorZEUS.push_back(sqrt(cov_mat[0][0]));						/// 6 - vtx pos err x
                OutPutVectorZEUS.push_back(sqrt(cov_mat[1][1]));						/// 7 - vtx pos err y
                OutPutVectorZEUS.push_back(sqrt(cov_mat[2][2]));						/// 8 - vtx pos err z	
		OutPutVectorZEUS.push_back(fitter->chi2()/fitter->ndof());					/// 9 - vtx normalized Chi2
		OutPutVectorZEUS.push_back(fitter->chi2());							/// 10 - Ndof	
		OutPutVectorZEUS.push_back(mass);								/// 11 - vertex mass
		OutPutVectorZEUS.push_back(1);									/// 12 - vertex valid 1/0

						}
		else {OutPutVectorZEUS.push_back(0);}
		
		ResetFinder();

		return OutPutVectorZEUS;



}


vector <double> *HelixTrans::MakeAFitwithVxlite(reco::VertexCollection::const_iterator ite)	{

	/// find tracks clusters from 	"TrackColForAlgorithm" track collection; (fill TracksCollectionClean[][])
	CreateTrackClusters();

	for( int i=0; i<10;i++)
        {	
		if (TracksCollectionClean[i].size() < 2) {continue;}

	if (debuger)	{	
		cout<<"*****************************"<<endl;
		cout<<"N tracks in the new track collection: "<<TracksCollectionClean[i].size()<<endl;
			}
			
		SetTrackCollectionToFit(TracksCollectionClean[i]);
		VertPos = FitVertexXYZ();
	
		vector<TransientTrack>::iterator trk_tst = TracksCollectionClean[i].begin();
		/// cout<<"Pt = "<<(trk_tst->track()).pt()<<"	Phi = "<<(trk_tst->track()).phi()<<"	Eta = "<<(trk_tst->track()).eta()<<endl;
		
		/// Extract DL and dDL;
// 		ExtractDL(i, ite);

		
		/// check for valid Error matrix, in composition with DL in could extract sqrt from a negative values;
		
		
	if (debuger)	{	
		/// My test for debug here:
		cout<<"Cluster tracks: "<<endl<<endl;
		int ccount = 0;
		for( vector<TransientTrack>::iterator  ittm1 = GeneralTrackCol.begin(); ittm1 !=GeneralTrackCol.end(); ittm1++)
		{  
		   const reco::TrackRef ref1 = (ittm1->trackBaseRef()).castTo<reco::TrackRef>();
		   for (vector<TransientTrack>::iterator inpa = TracksCollectionClean[i].begin(); inpa != TracksCollectionClean[i].end(); inpa++)	
		   {
			const reco::TrackRef ref2 = (inpa->trackBaseRef()).castTo<reco::TrackRef>();
			if (ref1.get() != ref2.get())	{continue;}
			
			double Drr = ((ittm1->track()).eta()-(trk_tst->track()).eta())*((ittm1->track()).eta()-(trk_tst->track()).eta());
			Drr += ((ittm1->track()).phi()-(trk_tst->track()).phi())*((ittm1->track()).phi()-(trk_tst->track()).phi());
			Drr = sqrt(Drr);
			
// 			cout<<"Phi = "<<(inpa->track()).phi()<<endl;
			cout<<"N = "<<ccount<<"  Pt = "<<(inpa->track()).pt()<<"	Phi = "<<(inpa->track()).phi()<<"	Eta = "<<(inpa->track()).eta()<<"	dR = "<<Drr<<endl;
		
		   }
		ccount++;

		}
			}	

	if (debuger)	{
		PrintFitResults();
			}

		ResetFinder();

		// next line is commented only for SL6 because of problems with isnan function: 	
		///if ( isnan(DLCov[0]) || isnan(DLCov[1])) {continue;}

		/// this method will return:   
		DLvector[0].push_back(DLCov[0]);	/// DL
		DLvector[1].push_back(DLCov[1]);	/// DL error
		DLvector[2].push_back(DLCov[2]);	/// mass	
		DLvector[3].push_back(VertPos[3]);	/// chi2/ndof




	}

	return DLvector;
}


/// extract 		tracks from general track collection which are the same for one chosen PV.
/// result: 		return Transient vector of all tracks from the one Primary vertex
/// input: 		ite is PV iterator; 
/// use global vars:  	"PVertices" = SetPVcollection(); 
vector<TransientTrack> HelixTrans::GetTracksFromPV(reco::VertexCollection::const_iterator ite)	{


	bool checkIDentity = false;
	int countss = 0;
	
	
	for( vector<TransientTrack>::iterator  it = GeneralTrackCol.begin(); it !=GeneralTrackCol.end(); it++)
        {       
        checkIDentity = false;
        const reco::TrackRef trackRef_tracks = (it->trackBaseRef()).castTo<reco::TrackRef>();

//      for(reco::VertexCollection::const_iterator ite = PVertices->begin(); ite !=PVertices->end(); ite++)     {

		
//              cout<<"PV size = "<<PVertices->size()<<endl;
                for(reco::Vertex::trackRef_iterator iTrack  =ite->tracks_begin(); iTrack != ite->tracks_end();++iTrack)         {

	                        const reco::TrackRef trackRef = iTrack->castTo<reco::TrackRef>();
                                if (trackRef.get() == trackRef_tracks.get()) {checkIDentity = true;}
//                              cout<<"Track is here: "<<"      "<<(it->track()).pt()<<"        "<<(*iTrack)->pt()<<endl;}      

                                                                                                                                }

//                                                                                                              }

                if (checkIDentity)      {    
// 		cout<<"Track Nr "<<countss<<endl;   
                transientTr = theTrkBuilder->build(trackRef_tracks);
                tracksVertex.push_back(transientTr);
                                        }

	countss++;
        }
                if (debuger)    { cout<<"N tracks from PV: "<<tracksVertex.size()<<endl; }

        return tracksVertex;

}


/// Set BeamSpot information here: 
/// result : 		global var "beamSpot" will contain BeamSpot object from the data file;
/// input : 		beamSpotHandle - BeamSpot object from data file;
/// use global vars : 	NO; 
void HelixTrans::SetBeamSpotInformation(Handle<reco::BeamSpot> beamSpotHandle)	{

 	if ( beamSpotHandle.isValid() )	{beamSpot = *beamSpotHandle;} else { cout<< "No beam spot available from EventSetup "<<endl;}

}


/// GetSimilarity: 	(vector1*M^-1*vector2)
/// result : 		return scalar value of calculation v1*M^{-1}*v2^{T}
/// input: 		v1, v2, inversed matrix M;
/// use global vars : 	NO;
double HelixTrans::GetSimilarity(double *v1, double *v2, ROOT::Math::SMatrix<double,3> M)	{


	double results = 0.;
	double v_tmp[3] = {0.};

	for (int i =0 ; i<3; i++)	{
		for (int j=0;j<3;j++)		{
			v_tmp[i] += M[i][j]*v2[j];
						}
					}

	for (int z=0;z<3; z++)	{	
		results += 	v1[z]*v_tmp[z];
				}

	return results;
}


/// Calculation of the DL for simple procedure of its calculation: 
/// result: 		return scalar DL value calcualted with simple projection method; 
/// input: 		ite is PV iterator; 
/// use global vars : 	DLCov;

vector <double> HelixTrans::ExtractDLVxlite(vector<TransientTrack> tracksVector, reco::VertexCollection::const_iterator ite, const VxLite::VertexFitter* fitterintrin, bool withbeam)	{

	vector <double> DLparameters; 

  	/// Define 3 covariance matrices
	typedef ROOT::Math::SMatrix<double,3>       SMatrix33;

	SMatrix33 DataVertexM;
	SMatrix33 SumVertexCovM;
	CLHEP::HepSymMatrix MyVertexM = fitterintrin->covariance();
// 	DataVertexM = ite->covariance();

	if (!withbeam)
	{
	DataVertexM = ite->covariance();
	} else
	{
	DataVertexM(0,0) = beamSpot.covariance(0, 0);
	DataVertexM(0,1) = beamSpot.covariance(0, 1);	
	DataVertexM(0,2) = 0;
	DataVertexM(1,0) = beamSpot.covariance(1, 0);
	DataVertexM(1,1) = beamSpot.covariance(1, 1);
	DataVertexM(1,2) = 0;
	DataVertexM(2,0) = 0;
	DataVertexM(2,1) = 0;
	DataVertexM(2,2) = 0;
	}


  
	/// set covariance matrices with sum of two matrix from PV and refited one;  
		for (int i =0; i<3; i++)		{
		for (int j=0; j<3; j++)		{	
	// cout<<MyVertexM[i][j]<<"	";
		SumVertexCovM[i][j] = MyVertexM[i][j] + DataVertexM[i][j];
						}	
	// cout<<endl;
							}

	/// **** ///
	/// test part
	double PxSum = 0.;
	double PySum = 0.;
	double PzSum = 0.;
	double PtSum = 0.;	
	
	for (vector<TransientTrack>::iterator itSumP = tracksVector.begin(); itSumP != tracksVector.end(); itSumP++)	
	{
		PxSum = PxSum + (itSumP->track()).px();
		PySum = PySum + (itSumP->track()).py();
		PzSum = PzSum + (itSumP->track()).pz();
	}

	PtSum = sqrt(PxSum*PxSum + PySum*PySum);

// 	cout<<"P standart = "<<sqrt(PxSum*PxSum + PySum*PySum + PzSum*PzSum)<<endl;  

	/// extract position of the PV from data
	double PVposX = ite->x();
	double PVposY = ite->y();
	double Lx = 0.;
	double Ly = 0.;
	double Lz = fitterintrin->vertex()[2]-ite->z();	
	
	if (!withbeam)
	{
	/// calculate DL projection components on XY plane;  
	Lx = fitterintrin->vertex()[0] - PVposX ;
	Ly = fitterintrin->vertex()[1] - PVposY;
	}	else 
	{
 	Lx = (fitterintrin->vertex()[0] - (beamSpot.x0() + beamSpot.dxdz()*Lz) );
	Ly = (fitterintrin->vertex()[1] - (beamSpot.y0() + beamSpot.dydz()*Lz) );	
										
	}

	/// calculate DL here (as a projection DL xy onto projection of the highest Pt track in XY plane)
	Float_t ProjectedDecayLength  =0;
	Float_t ProjectedDecayLengthError = 0.;	

	ProjectedDecayLength = (Lx*PxSum + Ly*PySum)/PtSum;

	ProjectedDecayLength = ProjectedDecayLength /*+ 0.001*/;

	/// Extract DLxy sigma;
	double DLSigma = 0.;
	double PClusterV[3] = {PxSum, PySum, 0.};					/// vector components of the DL 3D

	DLSigma = GetSimilarity(PClusterV, PClusterV, SumVertexCovM);
	DLSigma = DLSigma/(PtSum*PtSum);
	DLSigma = sqrt(DLSigma);
	ProjectedDecayLengthError = DLSigma;
	
	DLparameters.push_back(ProjectedDecayLength);
	DLparameters.push_back(ProjectedDecayLengthError);
	DLparameters.push_back(ProjectedDecayLength/ProjectedDecayLengthError);
	

	/// couts: 
	if (debuger)	{
	cout<<"Prima Vtx  = "<<ite->x()<<"	VtxY = "<<ite->y()<<"	VtxZ = "<<ite->z()<<endl;
// 	cout<<"Refit Vtx = "<<VertPos[0]<<"	VtxY = "<<VertPos[1]<<"	VtxZ = "<<VertPos[2]<<endl;
	cout<<"DLxyProj = "<<ProjectedDecayLength<<" DLsign = "<<ProjectedDecayLength/ProjectedDecayLengthError<<"  DLxy = "<<sqrt(Lx*Lx+Ly*Ly)<<endl;
			}

	return DLparameters;

}

/// test it for 1 track findings
void HelixTrans::FindPVfor2Tracks(reco::VertexCollection::const_iterator PVbegin, reco::VertexCollection::const_iterator PVend,  vector<TransientTrack> TracksForCheck)	{

        bool JPsiVtxFound = false;
	int GenCount = 0;

                int JPsiVtxCounter2 = 0;
                for(reco::VertexCollection::const_iterator ite = PVbegin; ite != PVend; ite++)
                {

                        int JPsiVtxCounter = 0;
                        for(reco::Vertex::trackRef_iterator iTrack  =ite->tracks_begin(); iTrack != ite->tracks_end();++iTrack)
                        {
                                const reco::TrackRef trackRef = iTrack->castTo<reco::TrackRef>();

                                        JPsiVtxFound = false;
                                        for( vector<TransientTrack>::iterator gt1 = TracksForCheck.begin(); gt1 !=TracksForCheck.end(); gt1++)
                                        {

                                                const reco::TrackRef trackRef2 = (gt1->trackBaseRef()).castTo<reco::TrackRef>();
                                                if (trackRef.get() == trackRef2.get())  {JPsiVtxFound = true; }
                                        }
                                        if (JPsiVtxFound) {JPsiVtxCounter++; GenCount++;}

                        }

                        JPsiVtxCounter2++;
			if (GenCount == 1 && (GenCount - JPsiVtxCounter) == 0 )	
				{
					iteFor2tracks = ite; 
					N2foundedVtxtracks = JPsiVtxCounter;
					N2foundedVtx = JPsiVtxCounter2-1;
				} else if (JPsiVtxCounter == 1) { iteFor2tracks = iteFor2trackstmp; N2foundedVtxtracks = 0; N2foundedVtx = -999;}

			if (JPsiVtxCounter > 1)		{ iteFor2tracks = ite; N2foundedVtxtracks = JPsiVtxCounter; N2foundedVtx = JPsiVtxCounter2-1;}
                }

/*
		if (iteFor2tracks != iteFor2trackstmp) 	{
			cout<<"FF : pos2 vtx = "<<iteFor2tracks->x()<<"	"<<iteFor2tracks->y()<<"	"<<iteFor2tracks->z()<<endl;
			cout<<"FF : Ntrack = "<<N2foundedVtxtracks<<endl;
			cout<<"FF : Vertex Nr# "<<N2foundedVtx<<endl;
							}
*/


}


int HelixTrans::FindPVfor2TracksGetNvertexTracks()	{

	if (iteFor2tracks != iteFor2trackstmp)  { return N2foundedVtxtracks; } else return 0; 
}

int HelixTrans::FindPVfor2TracksGetNvertex()	{

	if (iteFor2tracks != iteFor2trackstmp)  { return N2foundedVtx; } else return -999; 
}


VertexCollection::const_iterator HelixTrans::FindPVfor2TracksGetVertexIte()	{

	return iteFor2tracks; 
}

/*
/// Calculation of DL for procedure of its calculation with Covariance matrices: 
/// result : 		return scalar DL value calcualted with Covariance matrices; 
/// input: 		ite is PV iterator; 
/// use global vars : 	NO;
void HelixTrans::ExtractDLwithCov(int iterator, reco::VertexCollection::const_iterator ite)	{

	typedef ROOT::Math::SMatrix<double,3>       SMatrix33;

	SMatrix33 DataVertexM;
	SMatrix33 SumVertexCovM;
	CLHEP::HepSymMatrix MyVertexM = fitter->covariance();

	numerator = 0.;
	denominator = 0.;
	DLCov[0] = 0.;
	DLCov[1] = 0.;
	DLCov[2] = 0.;

	DataVertexM = ite->covariance();

// 	cout<<endl<<"Cov matrix from Data:"<<endl;
// 	DataVertexM.Print(std::cout);


// 	cout<<endl<<"Cov matrix from my refit: "<<endl;
	for (int i =0; i<3; i++)	{
		for (int j=0; j<3; j++)		{	
// 		cout<<MyVertexM[i][j]<<"	";
		SumVertexCovM[i][j] = MyVertexM[i][j] + DataVertexM[i][j];
						}	
// 						cout<<endl;
					}

// 	cout<<endl<<"Cov matrix from Data + myRefit:"<<endl;
// 	SumVertexCovM.Print(std::cout); cout<<endl;

	double PVposX = ite->x();
	double PVposY = ite->y();
	
// 	double PVposX = beamSpot.x0() + beamSpot.dxdz()* ite->z();
// 	double PVposY = beamSpot.y0() + beamSpot.dydz()* ite->z();
	

// 	double PVposZ = ite->z();

	vector<TransientTrack>::iterator  it = TracksCollectionClean[iterator].begin();
	

// 	double PxAx = (it->track()).px();
// 	double PyAx = (it->track()).py();
// // 	double PzAx = (it->track()).pz();
// 	double Ptx  = (it->track()).pt();

	double Lx = VertPos[0] - PVposX;
	double Ly = VertPos[1] - PVposY;
// 	double Lz = VertPos[2] - PVposZ;

// 	cout<<VertPos[0]<<"	"<<VertPos[1]<<"	"<<VertPos[2]<<endl;

	double DLSigma = 0.;
	double DL3D[3] = {Lx, Ly, 0.};					/// vector components of the DL 3D
//	double Axvector[3] = {PxAx/Ptx, PyAx/Ptx, 0.};			/// unit vector of the HighPt track axis 

	DLSigma = GetSimilarity(DL3D, DL3D, SumVertexCovM);


	double phiHPtrack = (it->track()).phi() ;
	if (phiHPtrack < 0)	{phiHPtrack =  2*M_PI - fabs(phiHPtrack);}
	cout<<" phiHPtrack = "<<phiHPtrack<<endl;
/// ************************
    Float_t ProjectedDecayLength=0;

    Float_t AxisSinPhi = sin(phiHPtrack);
    Float_t AxisCosPhi = cos(phiHPtrack);
    ProjectedDecayLength = Lx*AxisCosPhi+Ly*AxisSinPhi;


    // ProjectedDecayLengthError squared consists from three terms
    Float_t ProjectedDecayLengthError = 0.;	
    Float_t d1=cos(phiHPtrack);
    Float_t d2=sin(phiHPtrack);

    Double_t a = d1*d1*(SumVertexCovM[0][0]);
    Double_t b = d2*d2*(SumVertexCovM[1][1]);
    Double_t c = 2*d1*d2*SumVertexCovM[1][2];
    ProjectedDecayLengthError = sqrt(a + b + c);


	cout<<"ProjectedDecayLengthError here : "<<ProjectedDecayLengthError<<endl;
/// ************************

// 	bool mInverse = SumVertexCovM.Invert();
	
// 	cout<<"Inverse sum matrix here:"<<endl;
// 	SumVertexCovM.Print(std::cout);

	vector<TransientTrack>::iterator  ittmp = TracksCollectionClean[iterator].begin();

	HepLorentzVector p(0,0,0,0);
    	for(unsigned i = 0 ; i < TracksCollectionClean[iterator].size()-1; ++i) {
        const LinearizedTrack* t = finder.track(i);
	
	
// 	cout<<(ittmp->track()).p()<<endl;
        HepLorentzVector tp = t->momentum( (ittmp->track()).p() );
        p +=  t->weight() * tp;
	ittmp++;
    	}
    	double mass = p.m();
	if (debuger)	{cout<<"Mass = "<<mass<<endl;}


// 	if (mInverse)	{
// 
// 		numerator = GetSimilarity(Axvector, DL3D, SumVertexCovM);
// 		denominator = GetSimilarity(Axvector, Axvector, SumVertexCovM);
// 
// 	}
// 
// 	DLCov[0] = numerator/denominator;
// 	DLCov[1] =  DLSigma;



	double DLXY;
	
	DLXY = ((it->track()).px()*Lx + (it->track()).py()*Ly  )/( (it->track()).pt() );

// 	DLCov[0] = sqrt(Lx*Lx+Ly*Ly);
// 	DLCov[1] = (Lx*sqrt(SumVertexCovM[0][0]) + Ly*sqrt(SumVertexCovM[1][1]))/DLCov[0];

// 	DLCov[0] = DLXY;
// 	DLCov[1] = fabs(DLXY)*sqrt( (it->track()).px()*(it->track()).px()*SumVertexCovM[0][0]  + (it->track()).py()*(it->track()).py()*SumVertexCovM[1][1]) / ( (it->track()).pt() );
// 	DLCov[1] = fabs(DLCov[0])*sqrt(SumVertexCovM[1][1]+SumVertexCovM[0][0])/sqrt(Lx*Lx + Ly*Ly);


	DLCov[0] = ProjectedDecayLength;
	DLCov[1] = ProjectedDecayLengthError;

	DLCov[2] = mass;
if (debuger)	{
	cout<<"PVx = "<<ite->x()<<"	PVy = "<<ite->y()<<"	PVz = "<<ite->z()<<endl;
	cout<<"VtxX = "<<VertPos[0]<<"	VtxY = "<<VertPos[1]<<"	VtxZ = "<<VertPos[2]<<endl;
	cout<<"DLxy = "<<sqrt(Lx*Lx + Ly*Ly)<<endl;
// 	cout<<"Px "<<(it->track()).px()<<"	Py = "<<(it->track()).py()<<endl;
// 	cout<<"Pt = "<<(it->track()).pt()<<endl;
// 	cout<<"DLx = "<<Lx<<"	DLy = "<<Ly<<endl;
// 	cout<<"dDLx = "<<sqrt(SumVertexCovM[0][0])<<"	dDLx = "<<sqrt(SumVertexCovM[1][1])<<endl;
// 	cout<<"cos a 1= "<<(it->track()).px()*Lx + (it->track()).py()*Ly<<endl;
// 	cout<<"cos a 2= "<<((sqrt(Lx*Lx + Ly*Ly))*((it->track()).pt()))<<endl;

	cout<<"DL = "<<DLCov[0]<<"   SigmaDL = "<<DLCov[1]<<endl;
	cout<<"DL2 = "<<ProjectedDecayLength<<" DL2 error = "<<ProjectedDecayLengthError<<endl;
// 	cout<<"Nan check = "<<std::isnan(ProjectedDecayLengthError)<<endl;

	if (std::isnan(ProjectedDecayLengthError))	{
		cout<<"Matrix oo and 11 : "<<SumVertexCovM[0][0]<<"	"<<SumVertexCovM[1][1]<<"	"<<SumVertexCovM[1][2]<<endl;
		cout<<"Vertex fitted pos: "<<VertPos[0]<<"	"<<VertPos[1]<<"	"<<VertPos[2]<<endl;
		cout<<"cos and sin: "<<d1<<"	"<<d2<<endl;
		cout<<"abc "<<a<<"	"<<b<<"	"<<c<<"	"<<a + b + c<<endl;
	}
}
// 	return DLCov;

}
*/



