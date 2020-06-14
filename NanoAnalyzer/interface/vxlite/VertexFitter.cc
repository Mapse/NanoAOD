//*****************************************************************************
// File:      VertexFitter.cc
// ----------------------------------------------------------------------------
//=============================================================================
// RCS Current Revision Record
//-----------------------------------------------------------------------------
// $Source: /afs/desy.de/user/s/stadie/zeus/cvsroot/tlite/src/vxlite/VertexFitter.cc,v $
// $Revision: 1.17 $
// $Date: 2008/04/07 11:57:20 $
// $Author: stadie $
// $State: Exp $
// $Locker:  $
//*****************************************************************************
//=============================================================================
// Declarations and Definitions
//=============================================================================
#include "VertexFitter.hh"
#include "LinearizedTrack.hh"

#include "../Matrix.h"

#include <cassert>
#include <iostream> 

using namespace VxLite;
using CLHEP::HepVector;
using CLHEP::HepMatrix;
using CLHEP::HepSymMatrix;
using CLHEP::HepLorentzVector;

#define VXLITE_MAXTRACKS 100

extern "C" {
  void vxlitefit_(int& ntracks,float track_par[VXLITE_MAXTRACKS][5],
                  float track_cov[VXLITE_MAXTRACKS][15],float start_point[3],
                  float vertex[3],float vertex_cov[6],float& chi2,
                  float constrained_track_par[VXLITE_MAXTRACKS][3] )
  {
    static VertexFitter fitter;
    static std::vector<LinearizedTrack*> tracks;
    tracks.reserve(ntracks);
    for(int i = 0 ; i < ntracks ; ++i) {
      tracks.push_back(new LinearizedTrack(i+1,track_par[i],track_cov[i]));
    }
    HepVector start(3);
    start[0] =  start_point[0];
    start[1] =  start_point[1];
    start[2] =  start_point[2];
    fitter.setMaxIterations(1);
    fitter.setCalcMomentumCov(false);
    chi2 = fitter.fit(tracks.begin(),tracks.end(),start);
    //do not return fit result when a track was dropped
    //this might happen if one track has a bad covariance
    if( fitter.nTracks() != ntracks ) {
      chi2 = 9999999;
    }
    if(chi2 > 9999998) {
      for(int i = 0 ; i < 3 ; ++i ) vertex[i] = 0;
      for(int i = 0 ; i < 6 ; ++i ) vertex_cov[i] = 0;
      for(int i = 0 ; i < ntracks ; ++i) {
	constrained_track_par[i][0] = 0;
	constrained_track_par[i][1] = 0;
	constrained_track_par[i][2] = 0;
	delete tracks[i];
      }
      tracks.clear();
      return ;
    }
    vertex[0] = fitter.vertex()[0];
    vertex[1] = fitter.vertex()[1];
    vertex[2] = fitter.vertex()[2];
    vertex_cov[0] = fitter.covariance()[0][0];
    vertex_cov[1] = fitter.covariance()[0][1];
    vertex_cov[2] = fitter.covariance()[1][1];
    vertex_cov[3] = fitter.covariance()[0][2];
    vertex_cov[4] = fitter.covariance()[1][2];
    vertex_cov[5] = fitter.covariance()[2][2];
    for(int i = 0 ; i < ntracks ; ++i) {
      constrained_track_par[i][0] = M_PI_2 - atan(tracks[i]->q()[2]);//theta
      constrained_track_par[i][1] = (tracks[i]->q()[0] < 0) ? 2 * M_PI + tracks[i]->q()[0] : tracks[i]->q()[0];//phi
      constrained_track_par[i][2] = tracks[i]->q()[1];//Q/R
      delete tracks[i];
    }
    tracks.clear();
  }   

  void vxlitefitp_(int& ntracks,float track_par[VXLITE_MAXTRACKS][5],
		   float track_cov[VXLITE_MAXTRACKS][15], float track_p[VXLITE_MAXTRACKS],
		   float start_point[3],float vertex[3],float vertex_cov[6],float& chi2,
		   float constrained_track_par[VXLITE_MAXTRACKS][3], 
                                float constrained_track_p[VXLITE_MAXTRACKS], float& mass )
  {
    static VertexFitter fitter;
    static std::vector<LinearizedTrack*> tracks;
    tracks.reserve(ntracks);
    for(int i = 0 ; i < ntracks ; ++i) {
      tracks.push_back(new LinearizedTrack(i+1,track_par[i],track_cov[i]));
    }
    HepVector start(3);
    start[0] =  start_point[0];
    start[1] =  start_point[1];
    start[2] =  start_point[2];
    fitter.setMaxIterations(1);
    fitter.setCalcMomentumCov(false);
    chi2 = fitter.fit(tracks.begin(),tracks.end(),start);
    //do not return fit result when a track was dropped
    //this might happen if one track has a bad covariance
    if( fitter.nTracks() != ntracks ) {
      chi2 = 9999999;
    }
    if(chi2 > 9999998) {
      for(int i = 0 ; i < 3 ; ++i ) vertex[i] = 0;
      for(int i = 0 ; i < 6 ; ++i ) vertex_cov[i] = 0;
      for(int i = 0 ; i < ntracks ; ++i) {
	constrained_track_par[i][0] = 0;
	constrained_track_par[i][1] = 0;
	constrained_track_par[i][2] = 0;
	constrained_track_p[i] = 0;
	mass = 0;
	delete tracks[i];
      }
      tracks.clear();
      return ;
    }
    vertex[0] = fitter.vertex()[0];
    vertex[1] = fitter.vertex()[1];
    vertex[2] = fitter.vertex()[2];
    vertex_cov[0] = fitter.covariance()[0][0];
    vertex_cov[1] = fitter.covariance()[0][1];
    vertex_cov[2] = fitter.covariance()[1][1];
    vertex_cov[3] = fitter.covariance()[0][2];
    vertex_cov[4] = fitter.covariance()[1][2];
    vertex_cov[5] = fitter.covariance()[2][2];
    HepLorentzVector p(0,0,0,0);
    for(int i = 0 ; i < ntracks ; ++i) {
      constrained_track_par[i][0] = M_PI_2 - atan(tracks[i]->q()[2]);//theta
      constrained_track_par[i][1] = (tracks[i]->q()[0] < 0) ? 2 * M_PI + tracks[i]->q()[0] : tracks[i]->q()[0];//phi
      constrained_track_par[i][2] = tracks[i]->q()[1];//Q/R
      HepLorentzVector tp = tracks[i]->momentum(track_p[i]);
      constrained_track_p[i] = tp.v().mag();
      p += tp;
      delete tracks[i];
    }
    tracks.clear();
    mass = p.m();
  }   
    
  void vxlitefit2d_(int& ntracks,float track_par[VXLITE_MAXTRACKS][5],
		    float track_cov[VXLITE_MAXTRACKS][15], float track_p[VXLITE_MAXTRACKS],
		    float start_point[3], float vertex[3],float vertex_cov[6],float& chi2,
		    float constrained_track_par[VXLITE_MAXTRACKS][3], float constrained_track_p[VXLITE_MAXTRACKS], float& mass )
  {
    static VertexFitter fitter;
    static std::vector<LinearizedTrack*> tracks;
    tracks.reserve(ntracks);
    for(int i = 0 ; i < ntracks ; ++i) {
      tracks.push_back(new LinearizedTrack(i+1,track_par[i],track_cov[i],LinearizedTrack::PI_PLUS_MASS,
                                true));
    }
    HepVector start(3);
    start[0] =  start_point[0];
    start[1] =  start_point[1];
    start[2] =  start_point[2];
    fitter.setMaxIterations(1);
    fitter.setCalcMomentumCov(false);
    chi2 = fitter.fit(tracks.begin(),tracks.end(),start);
    //do not return fit result when a track was dropped
    //this might happen if one track has a bad covariance
    if( fitter.nTracks() != ntracks ) {
      chi2 = 9999999;
    }
    if(chi2 > 9999998) {
      for(int i = 0 ; i < 3 ; ++i ) vertex[i] = 0;
      for(int i = 0 ; i < 6 ; ++i ) vertex_cov[i] = 0;
      for(int i = 0 ; i < ntracks ; ++i) {
	constrained_track_par[i][0] = 0;
	constrained_track_par[i][1] = 0;
	constrained_track_par[i][2] = 0;
	delete tracks[i];
      }
      tracks.clear();
      return ;
    }
    vertex[0] = fitter.vertex()[0];
    vertex[1] = fitter.vertex()[1];
    vertex[2] = fitter.vertex()[2];
    vertex_cov[0] = fitter.covariance()[0][0];
    vertex_cov[1] = fitter.covariance()[0][1];
    vertex_cov[2] = fitter.covariance()[1][1];
    vertex_cov[3] = fitter.covariance()[0][2];
    vertex_cov[4] = fitter.covariance()[1][2];
    vertex_cov[5] = fitter.covariance()[2][2];
    HepLorentzVector p(0,0,0,0);
    for(int i = 0 ; i < ntracks ; ++i) {
      constrained_track_par[i][0] = M_PI_2 - atan(track_par[i][4]);//theta
      constrained_track_par[i][1] = (tracks[i]->q()[0] < 0) ? 2 * M_PI + tracks[i]->q()[0] : tracks[i]->q()[0];//phi
      constrained_track_par[i][2] = tracks[i]->q()[1];//Q/R
      double oldpt = track_p[i]  * std::abs(sin(atan(1/track_par[i][4])));
      double pt = oldpt * track_par[i][1] / tracks[i]->q()[1];
      double E = sqrt(LinearizedTrack::PI_PLUS_MASS * LinearizedTrack::PI_PLUS_MASS + 
		      pt * pt * (1 + track_par[i][4] * track_par[i][4]));
      HepLorentzVector tp(cos(tracks[i]->q()[0]) * pt,sin(tracks[i]->q()[0]) * pt,tracks[i]->q()[2] * pt,E);
      p += tp;
      constrained_track_p[i] = tp.v().mag();
      delete tracks[i];
    }
    tracks.clear();
    mass = p.m();
  }   

  void vxlitefitcov_(int& ntracks,float track_par[VXLITE_MAXTRACKS][5],
		     float track_cov[VXLITE_MAXTRACKS][15],
		     float start_point[3],float vertex[3],float vertex_cov[6],
		     float& chi2,
		     float constrained_track_par[VXLITE_MAXTRACKS][3],
		     float constrained_track_cov[VXLITE_MAXTRACKS][6])
  {
    static VertexFitter fitter;
    static std::vector<LinearizedTrack*> tracks;
    tracks.reserve(ntracks);
    for(int i = 0 ; i < ntracks ; ++i) {
      tracks.push_back(new LinearizedTrack(i+1,track_par[i],track_cov[i]));
    }
    HepVector start(3);
    start[0] =  start_point[0];
    start[1] =  start_point[1];
    start[2] =  start_point[2];
    fitter.setMaxIterations(1);
    fitter.setCalcMomentumCov(false);
    chi2 = fitter.fit(tracks.begin(),tracks.end(),start);   
    if( fitter.nTracks() != ntracks ) {
      chi2 = 9999999;
    }
    if(chi2 > 9999998) {     
      for(int i = 0 ; i < 3 ; ++i ) vertex[i] = 0;
      for(int i = 0 ; i < 6 ; ++i ) vertex_cov[i] = 0;
      for(int i = 0 ; i < ntracks ; ++i) {
	constrained_track_par[i][0] = 0;
        constrained_track_par[i][1] = 0;
	constrained_track_par[i][2] = 0;
	delete tracks[i];
      }
      tracks.clear();
      return ;
    }
    vertex[0] = fitter.vertex()[0];
    vertex[1] = fitter.vertex()[1];
    vertex[2] = fitter.vertex()[2];
    vertex_cov[0] = fitter.covariance()[0][0];
    vertex_cov[1] = fitter.covariance()[0][1];
    vertex_cov[2] = fitter.covariance()[1][1];
    vertex_cov[3] = fitter.covariance()[0][2];
    vertex_cov[4] = fitter.covariance()[1][2];
    vertex_cov[5] = fitter.covariance()[2][2];      
    HepMatrix dZeusdq(3,3,0);
    dZeusdq[1][0] = 1;
    dZeusdq[2][1] = 1;
    for(int i = 0 ; i < ntracks ; ++i) {
      LinearizedTrack* t = tracks[i];
      constrained_track_par[i][0] = M_PI_2 - atan(t->q()[2]);//theta
      constrained_track_par[i][1] = (t->q()[0] < 0) ? 2 * M_PI + t->q()[0] : t->q()[0];//phi
      constrained_track_par[i][2] = t->q()[1];//Q/R
      
      HepSymMatrix covii = t->E() + 
	fitter.covariance().similarity(t->EBG() * t->A());
      dZeusdq[0][2] = -1 /(1 + t->q()[2] * t->q()[2]);
      covii = covii.similarity(dZeusdq);
      constrained_track_cov[i][0] =  covii[0][0];
      constrained_track_cov[i][1] =  covii[0][1];
      constrained_track_cov[i][2] =  covii[1][1];
      constrained_track_cov[i][3] =  covii[2][0];
      constrained_track_cov[i][4] =  covii[2][1];
      constrained_track_cov[i][5] =  covii[2][2];
      delete t;
    }
    //std::cout << "PCov:" << fitter.momentum_covariance();
    //std::cout << "mass:" << fitter.mass() << "\n";;  
    tracks.clear();
  } 

  void vxlitefitit_(int& ntracks,float track_par[VXLITE_MAXTRACKS][5],
		    float track_cov[VXLITE_MAXTRACKS][15],float start_point[3],
		    int& iterations,float vertex[3],float vertex_cov[6],
		    float& chi2,
		    float constrained_track_par[VXLITE_MAXTRACKS][3])
  {
    static VertexFitter fitter;
    static std::vector<LinearizedTrack*> tracks;
    tracks.reserve(ntracks);
    for(int i = 0 ; i < ntracks ; ++i) {
      tracks.push_back(new LinearizedTrack(i+1,track_par[i],track_cov[i]));
    }
    HepVector start(3);
    start[0] =  start_point[0];
    start[1] =  start_point[1];
    start[2] =  start_point[2];
    if(iterations > 0) {
      fitter.setMaxIterations(iterations);
    } else {
      fitter.setMaxIterations(1);
    }
    fitter.setCalcMomentumCov(false);
    chi2 = fitter.fit(tracks.begin(),tracks.end(),start);   
    if( fitter.nTracks() != ntracks ) {
      chi2 = 9999999;
    }
    if(chi2 > 9999998) {     
      for(int i = 0 ; i < 3 ; ++i ) vertex[i] = 0;
      for(int i = 0 ; i < 6 ; ++i ) vertex_cov[i] = 0;
      for(int i = 0 ; i < ntracks ; ++i) {
	constrained_track_par[i][0] = 0;
        constrained_track_par[i][1] = 0;
	constrained_track_par[i][2] = 0;
	delete tracks[i];
      }
      tracks.clear();
      return ;
    }
    vertex[0] = fitter.vertex()[0];
    vertex[1] = fitter.vertex()[1];
    vertex[2] = fitter.vertex()[2];
    vertex_cov[0] = fitter.covariance()[0][0];
    vertex_cov[1] = fitter.covariance()[0][1];
    vertex_cov[2] = fitter.covariance()[1][1];
    vertex_cov[3] = fitter.covariance()[0][2];
    vertex_cov[4] = fitter.covariance()[1][2];
    vertex_cov[5] = fitter.covariance()[2][2];
    for(int i = 0 ; i < ntracks ; ++i) {
      constrained_track_par[i][0] = M_PI_2 - atan(tracks[i]->q()[2]);//theta
      constrained_track_par[i][1] = (tracks[i]->q()[0] < 0) ? 2 * M_PI + tracks[i]->q()[0] : tracks[i]->q()[0];//phi
      constrained_track_par[i][2] = tracks[i]->q()[1];//Q/R
      delete tracks[i];
    }
    //std::cout << "PCov:" << fitter.momentum_covariance();
    //std::cout << "mass:" << fitter.mass() << "\n";;  
    tracks.clear();
  }  

  void vxlitefititcstr_(int& ntracks,float track_par[VXLITE_MAXTRACKS][5],
			float track_cov[VXLITE_MAXTRACKS][15],int& iterations,
			float vertex[3],float vertex_cov[6],float& chi2,
			float constrained_track_par[VXLITE_MAXTRACKS][3])
  {
    static VertexFitter fitter;
    static std::vector<LinearizedTrack*> tracks;
    tracks.reserve(ntracks);
    for(int i = 0 ; i < ntracks ; ++i) {
      tracks.push_back(new LinearizedTrack(i+1,track_par[i],track_cov[i]));
    }
    HepVector start(3);
    start[0] =  vertex[0];
    start[1] =  vertex[1];
    start[2] =  vertex[2];
    HepSymMatrix startcov(3);
    startcov[0][0] = vertex_cov[0];
    startcov[0][1] = vertex_cov[1];
    startcov[1][1] = vertex_cov[2];
    startcov[0][2] = vertex_cov[3];
    startcov[1][2] = vertex_cov[4];
    startcov[2][2] = vertex_cov[5];
    if(iterations > 0) {
      fitter.setMaxIterations(iterations);
    } else {
      fitter.setMaxIterations(1);
    }
    fitter.setCalcMomentumCov(false);
    chi2 = fitter.constrained_fit(tracks.begin(),tracks.end(),start,start,
				  startcov);   
    if( fitter.nTracks() != ntracks ) {
      chi2 = 9999999;
    }
    if(chi2 > 9999998) {     
      for(int i = 0 ; i < 3 ; ++i ) vertex[i] = 0;
      for(int i = 0 ; i < 6 ; ++i ) vertex_cov[i] = 0;
      for(int i = 0 ; i < ntracks ; ++i) {
	constrained_track_par[i][0] = 0;
        constrained_track_par[i][1] = 0;
	constrained_track_par[i][2] = 0;
	delete tracks[i];
      }
      tracks.clear();
      return ;
    }
    vertex[0] = fitter.vertex()[0];
    vertex[1] = fitter.vertex()[1];
    vertex[2] = fitter.vertex()[2];
    vertex_cov[0] = fitter.covariance()[0][0];
    vertex_cov[1] = fitter.covariance()[0][1];
    vertex_cov[2] = fitter.covariance()[1][1];
    vertex_cov[3] = fitter.covariance()[0][2];
    vertex_cov[4] = fitter.covariance()[1][2];
    vertex_cov[5] = fitter.covariance()[2][2];
    for(int i = 0 ; i < ntracks ; ++i) {
      constrained_track_par[i][0] = M_PI_2 - atan(tracks[i]->q()[2]);//theta
      constrained_track_par[i][1] = (tracks[i]->q()[0] < 0) ? 2 * M_PI + tracks[i]->q()[0] : tracks[i]->q()[0];//phi
      constrained_track_par[i][2] = tracks[i]->q()[1];//Q/R
      delete tracks[i];
    }
    //std::cout << "PCov:" << fitter.momentum_covariance();
    //std::cout << "mass:" << fitter.mass() << "\n";;  
    tracks.clear();
  }
}

VertexFitter::VertexFitter()
  : m_x(3,0),m_covx(3,0),m_p(0,0,0,0),m_covp(4,0),m_chi2(0),m_iteration(0), 
    m_charge(0),m_sumofweights(0),m_ndof(0),m_max_iter(1),m_calc_covp(false),
    m_c00inv(3,0),m_xsum(3,0)
{
}

double VertexFitter::fit(const std::vector<LinearizedTrack*>::iterator& begin, 
			 const std::vector<LinearizedTrack*>::iterator& end,
			 const HepVector& point)
{
  const int ntracks = distance(begin,end);
  if((ntracks < 2) ||
      (point[0] * point[0] + point[1] * point[1] > 40000)) {
    m_iteration = 0; 
    m_x = HepVector(3,0);
    m_covx = HepSymMatrix(3,0);
    m_c00inv = HepSymMatrix(3,0);
    m_p.set(0,0,0,0);
    m_charge = 0;
    m_sumofweights = 0;
    m_ndof = 0;
    m_chi2 = 9999999;    
    return m_chi2;
  }
  m_x = point;
  const double chiold =  ( m_iteration == 0 ) ? 9999999 : m_chi2;
  calcXSum(begin,end,point);
  if(m_sumofweights < 1.5) {
    m_iteration = 0; 
    m_x = HepVector(3,0);
    m_covx = HepSymMatrix(3,0);
    m_c00inv = HepSymMatrix(3,0);
    m_p.set(0,0,0,0);
    m_charge = 0;
    m_sumofweights = 0;
    m_ndof = 0;
    m_chi2 = 9999999;    
    return m_chi2;
  }
  //std::cout << "Xsum:" << m_xsum[0] << ", " << m_xsum[1] << ", " << m_xsum[2] << "\n";
  //std::cout << "c00inv:" << c00inv;
  m_covx = m_c00inv;
  int errorcode;
  m_covx.invert(errorcode);
  assert(errorcode == 0);
  m_x =  m_covx * m_xsum;
  calcChi2(begin,end);
  ++m_iteration;
  double deltachi = chiold - m_chi2;
  if((m_iteration < m_max_iter) && 
      (std::abs(deltachi) > 0.01 * (2 * m_sumofweights - 1))) {
    //std::cout << "repeating fit:" << m_x[0] << ", " << m_x[1] << ", "
    //	      << m_x[2] << " chi2:" << m_chi2 << "\n";
    if (deltachi<0) { //Fit may be oscillating
      const HepVector newpoint = (m_x + point)/2;
      for(std::vector<LinearizedTrack*>::iterator t = begin ; t < end ; ++t) {
	(*t)->expandAround(newpoint);
      }
      return fit(begin,end,newpoint);
    }
    return fit(begin,end,m_x);
  }
  m_iteration = 0; 
  calcP(begin,end);
  if(m_calc_covp) {
    calcPCov(begin,end);
  }
  m_ndof = m_sumofweights * 2 - 3;
  return m_chi2;
}

double VertexFitter::constrained_fit(const std::vector<LinearizedTrack*>::iterator& begin, const std::vector<LinearizedTrack*>::iterator& end,const HepVector& point, const HepVector& xin, const HepSymMatrix& covxin)
{
  int ntracks = distance(begin,end);
  if((ntracks < 1) ||
      (point[0] * point[0] + point[1] * point[1] > 40000)) { 
    m_iteration = 0; 
    m_x = HepVector(3,0);
    m_covx = HepSymMatrix(3,0);
    m_c00inv = HepSymMatrix(3,0);
    m_p.set(0,0,0,0);
    m_charge = 0;
    m_sumofweights = 0;
    m_ndof = 0;
    m_chi2 = 9999999;
    return m_chi2;
  }
  m_x = point;
  double chiold =  ( m_iteration == 0 ) ? 9999999 : m_chi2;
  calcXSum(begin,end,point);  
  if(m_sumofweights == 0) {
    m_iteration = 0; 
    m_x = HepVector(3,0);
    m_covx = HepSymMatrix(3,0);
    m_c00inv = HepSymMatrix(3,0);
    m_p.set(0,0,0,0);
    m_charge = 0;
    m_sumofweights = 0;
    m_ndof = 0;
    m_chi2 = 9999999;
    return m_chi2;
  }
  //std::cout << "Xsum:" << xsum[0] << ", " << xsum[1] << ", " << xsum[2] << "\n";
  //std::cout << "c00inv:" << c00inv;
 
  //add contribution from vertex constraint
  HepSymMatrix covxin_inv = covxin;
  int errorcode;
  covxin_inv.invert(errorcode);
  assert(errorcode == 0);
  m_c00inv += covxin_inv;
  m_xsum += covxin_inv * xin;

  m_covx = m_c00inv;
  m_covx.invert(errorcode);
  assert(errorcode == 0);
  m_x =  m_covx * m_xsum;
  calcChi2(begin,end);
  m_chi2 += covxin_inv.similarity(m_x - xin);
  ++m_iteration;
  double deltachi = chiold - m_chi2;
  if((m_iteration < m_max_iter) && 
      (std::abs(deltachi) > 0.01 * (2 * m_sumofweights - 1))) {
    //std::cout << "repeating fit:" << m_x[0] << ", " << m_x[1] << ", "
    //	      << m_x[2] << " chi2:" << m_chi2 << "\n";
    if (deltachi<0) { //Fit may be oscillating
      HepVector newpoint = (m_x + point)/2;
      for(std::vector<LinearizedTrack*>::iterator t = begin ; t < end ; ++t) {
	(*t)->expandAround(newpoint);
      }
      return constrained_fit(begin,end,newpoint,xin,covxin);
    }
    return constrained_fit(begin,end,m_x,xin,covxin);
  }
  m_iteration = 0;
  calcP(begin,end);
  if(m_calc_covp) {
    calcPCov(begin,end);
  } 
  m_ndof = m_sumofweights * 2; 

  return m_chi2;
}

void VertexFitter::calcXSum(const std::vector<LinearizedTrack*>::iterator &begin,
			    const std::vector<LinearizedTrack*>::iterator &end,
			    const HepVector &point)
{
  m_c00inv = HepSymMatrix(3,0);
  m_xsum = HepVector(3,0);
  m_sumofweights = 0;
  for(std::vector<LinearizedTrack*>::iterator t = begin ; t < end ; ++t) {
    LinearizedTrack* track = *t;
    if(m_iteration == 0) track->expandAround(point);
    if(! track->expanded()) { 
      assert(track->weight() == 0);
      continue;
    }
    //track->show(); 
    
    HepVector p =  track->alpha() - track->h();
    //std::cout << "p: " << p[0] << ", " << p[1] << ", " << p[2] << ", "
    //	      << p[3] << ", " << p[4] << "\n";
    HepMatrix AW(track->A().T() * track->W());
    HepSymMatrix AWA(track->W().similarityT(track->A()));
    m_c00inv += track->weight() * track->weight() * AWA;
    m_xsum += track->weight() *  track->weight() * AW * p;
    m_sumofweights += track->weight();
  }
}

void VertexFitter::calcChi2(const std::vector<LinearizedTrack*>::iterator &begin,
			     const std::vector<LinearizedTrack*>::iterator &end)
{
  m_chi2 = 0;
  for(std::vector<LinearizedTrack*>::iterator t = begin ; t < end ; ++t) {
    LinearizedTrack* track = *t;  
    if(! track->expanded()) {
      assert(track->weight() == 0);
      continue;
    }
    HepVector p(track->alpha());
    p -= track->h(); 
    HepVector q(track->EBG() * (p - track->A() * m_x));
    p -= track->AxBq(m_x,q);
    track->setChi2(track->G().similarity(p));
    m_chi2 += track->weight() * track->weight() * track->chi2();
    if(m_iteration < m_max_iter) {
      track->expandAround(m_x,q);
    }
  }
}

void VertexFitter::calcP(const std::vector<LinearizedTrack*>::iterator &begin,
			 const std::vector<LinearizedTrack*>::iterator &end)
{
  m_p.set(0,0,0,0);
  m_charge = 0; 
  for(std::vector<LinearizedTrack*>::iterator t = begin ; t < end ; ++t) {
    LinearizedTrack* track = *t;  
    if(! track->expanded()) {
      assert(track->weight() == 0);
      continue;
    }
    HepVector q = track->EBG() * (track->alpha() - track->h() - 
				  track->A() * m_x);
    m_p += track->weight() * track->momentum(q);
    m_charge += track->weight() *  track->charge();
  }
}

void VertexFitter::calcPCov(const std::vector<LinearizedTrack*>::iterator &begin,
			    const std::vector<LinearizedTrack*>::iterator &end)
{
  //calculate momentum cov
  const int ntracks = distance(begin,end);
  HepMatrix** CQ = new HepMatrix*[ntracks];
  for(int i = 0 ; i < ntracks ; ++i) {
    CQ[i] = new HepMatrix[ntracks];
  }
  int i = 0;
  for(std::vector<LinearizedTrack*>::iterator t1 = begin ; t1 < end ; ++t1) {
    const LinearizedTrack* track1 = *t1;    
    if(! track1->expanded()) {
      assert(track1->weight() == 0);
      continue;
    }
    HepMatrix d1 = track1->EBG() * track1->A();
    int j = i;
    for(std::vector<LinearizedTrack*>::iterator t2 = t1 ; t2 < end ; ++t2) {
      const LinearizedTrack* track2 = *t2;   
      if(! track2->expanded()) {
	assert(track2->weight() == 0);
	continue;
      }
      if(i == j) {
	CQ[i][j] = track1->E();
	CQ[i][j] +=  m_covx.similarity(d1);
      } else {
	HepMatrix d2 = track2->EBG() * track2->A();
	CQ[i][j] = d1 * m_covx * d2.T();
	CQ[j][i] = CQ[i][j].T();
      }
      ++j;
    }
    ++i;
  }
  HepMatrix covp = HepMatrix(4,4,0); 
  i = 0;
  for(std::vector<LinearizedTrack*>::iterator t1 = begin ; t1 < end ; ++t1) {
    const LinearizedTrack* track1 = *t1;   
    if(! track1->expanded()) {
      assert(track1->weight() == 0);
      continue;
    }
    int j = 0;
    for(std::vector<LinearizedTrack*>::iterator t2 = begin ; t2 < end ; ++t2) {
      const LinearizedTrack* track2 = *t2;  
      if(! track2->expanded()) {
	assert(track2->weight() == 0);
	continue;
      }
      covp += track1->dPdq() * CQ[i][j] * track2->dPdq().T();
      ++j;
    }
    ++i;
  }
  for(int j = 0 ; j < 4 ; ++j) {
    for(int k = j ; k < 4 ; ++k) {
      m_covp[j][k] = covp[j][k];
    }
  }
  for(int j = 0 ; j < ntracks ; ++j) {
    delete[]  CQ[j];
  }
  delete[] CQ;
}

void VertexFitter::show() const 
{
  std::cout << "Vertex : (" << m_x[0] << ", " << m_x[1] << ", " << m_x[2]
	    << ") Chi2: " << m_chi2 << " ndof: " << ndof() << "\n";
  std::cout << "Cov:" << m_covx;
  std::cout << "nTracks: " << nTracks() << " charge: " << charge() 
	    << " mass: " << mass() << "\n";
}
 
const double VertexFitter::chi2ForTrack(const LinearizedTrack* const track, 
					bool add)
{
  if(! track->expanded()) {
    assert(track->weight() == 0);
    return 9999999;
  }
  if(!add) {
    const HepVector p =  track->alpha() - track->h();
    const HepVector q = track->EBG() * (p - track->A() * m_x);
    return track->weight() * track->G().similarity(p - track->AxBq(m_x,q));
  }
  HepSymMatrix c00(m_c00inv);
  int errorcode;
  HepVector xsum(c00 * m_x);
  const HepVector p =  track->alpha() - track->h();
  const HepMatrix AW(track->A().T() * track->W());
  const HepSymMatrix AWA(track->W().similarityT(track->A()));
  c00 += track->weight() * track->weight() * AWA;
  xsum += track->weight() * track->weight() * AW * p;
  c00.invert(errorcode);
  assert(errorcode == 0);
  const HepVector x = c00 * xsum;
  const HepVector q = track->EBG() * (p - track->A() * x);
  return track->weight() * track->weight() * track->G().similarity(p - track->AxBq(x,q));
}
 
double VertexFitter::mass_error() const
{
  const double m = m_p.m();
  if(m == 0) return 0;
  HepVector dmdp(4);
  dmdp[0] = - m_p[0]/m;
  dmdp[1] = - m_p[1]/m;
  dmdp[2] = - m_p[2]/m;
  dmdp[3] =   m_p[3]/m;
  return sqrt(m_covp.similarity(dmdp));
}
