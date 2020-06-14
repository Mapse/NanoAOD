//*****************************************************************************
// File:      DAFVertexFinder.cc
// ----------------------------------------------------------------------------
//=============================================================================
// RCS Current Revision Record
//-----------------------------------------------------------------------------
// $Source: /afs/desy.de/user/s/stadie/zeus/cvsroot/tlite/src/vxlite/DAFVertexFinder.cc,v $
// $Revision: 1.8 $
// $Date: 2008/04/07 11:57:20 $
// $Author: stadie $
// $State: Exp $
// $Locker:  $
//*****************************************************************************
//=============================================================================
// Declarations and Definitions
//=============================================================================
#include "DAFVertexFinder.hh"
#include "VertexFitter.hh"
#include "LinearizedTrack.hh"

#include <cassert>
#include <iostream> 

  /*
    Implementation of a deterministic annealing filter
    The filter performs the annealing at four steps with fixed
    temperatures. Tracks that have a chi square above chi_cut are
    downweighted.
  */

using namespace VxLite;
using CLHEP::HepVector;
using CLHEP::HepMatrix;
using CLHEP::HepSymMatrix;
using CLHEP::HepLorentzVector;

#define VXLITE_MAXTRACKS 100

extern "C" {
  void vxlitedaffit_(int& ntracks,float track_par[VXLITE_MAXTRACKS][5],
		     float track_cov[VXLITE_MAXTRACKS][15],float vertex[3],
		     float vertex_cov[6],float& chi2_cut,int& mode,float& chi2,
		     float& ndof,float weights[VXLITE_MAXTRACKS],
		     float constrained_track_par[VXLITE_MAXTRACKS][3])
  {
    static DAFVertexFinder finder;

    finder.setChi2Cut(chi2_cut);
    for(int i = 0 ; i < ntracks ; ++i) {
      finder.addTrack(new LinearizedTrack(i+1,track_par[i],track_cov[i]));
    }    
    
    double sumofweights = 0;
    
    if(mode == 1) {
      //with vertex constraint
      HepVector xin(3);
      xin[0] =  vertex[0];
      xin[1] =  vertex[1];
      xin[2] =  vertex[2];
      HepSymMatrix covxin(3);
      covxin[0][0] = vertex_cov[0];
      covxin[0][1] = vertex_cov[1];
      covxin[1][1] = vertex_cov[2];
      covxin[0][2] = vertex_cov[3];
      covxin[1][2] = vertex_cov[4];
      covxin[2][2] = vertex_cov[5];
      
      sumofweights =  finder.findVertex(xin,covxin);
    } else {
      sumofweights =  finder.findVertex();
    }
    const VertexFitter* fitter = finder.fitter();
    if(sumofweights == 0) {     
      for(int i = 0 ; i < 3 ; ++i ) vertex[i] = 0;
      for(int i = 0 ; i < 6 ; ++i ) vertex_cov[i] = 0;
      for(int i = 0 ; i < ntracks ; ++i) {
	constrained_track_par[i][0] = 0;
	constrained_track_par[i][1] = 0;
	constrained_track_par[i][2] = 0;
	weights[i] = 0;
      }
      ndof = 0;
      chi2 = 9999999;
      finder.reset();
      return ;
    }
    vertex[0] = fitter->vertex()[0];
    vertex[1] = fitter->vertex()[1];
    vertex[2] = fitter->vertex()[2];
    vertex_cov[0] = fitter->covariance()[0][0];
    vertex_cov[1] = fitter->covariance()[0][1];
    vertex_cov[2] = fitter->covariance()[1][1];
    vertex_cov[3] = fitter->covariance()[0][2];
    vertex_cov[4] = fitter->covariance()[1][2];
    vertex_cov[5] = fitter->covariance()[2][2];
    chi2 = fitter->chi2();
    ndof = fitter->ndof();
    for(int i = 0 ; i < ntracks ; ++i) {
      const LinearizedTrack* t = finder.track(i);
      if(! t->expanded()) {
	assert(t->weight() == 0);
	constrained_track_par[i][0] = 0;//theta
	constrained_track_par[i][1] = 0;//phi
	constrained_track_par[i][2] = t->charge() * 999999;//Q/R
	weights[i] = 0;
	continue;
      }
      constrained_track_par[i][0] = M_PI_2 - atan(t->q()[2]);//theta
      constrained_track_par[i][1] = (t->q()[0] < 0) ? 2 * M_PI + t->q()[0] : t->q()[0];//phi
      constrained_track_par[i][2] = t->q()[1];//Q/R
      weights[i] = t->weight();
    }
    finder.reset();
  }

  void vxlitedaffitp_(int& ntracks,float track_par[VXLITE_MAXTRACKS][5],
		      float track_cov[VXLITE_MAXTRACKS][15],float track_p[VXLITE_MAXTRACKS],
		      float vertex[3],float vertex_cov[6],float& chi2_cut,int& mode,float& chi2,
		      float& ndof,float weights[VXLITE_MAXTRACKS],
		      float constrained_track_par[VXLITE_MAXTRACKS][3],float constrained_track_p[VXLITE_MAXTRACKS], float& mass )
  { 
    static DAFVertexFinder finder;
    
    finder.setChi2Cut(chi2_cut);
    for(int i = 0 ; i < ntracks ; ++i) {
      finder.addTrack(new LinearizedTrack(i+1,track_par[i],track_cov[i]));
    }    
    
    double sumofweights = 0;
    
    if(mode == 1) {
      //with vertex constraint
      HepVector xin(3);
      xin[0] =  vertex[0];
      xin[1] =  vertex[1];
      xin[2] =  vertex[2];
      HepSymMatrix covxin(3);
      covxin[0][0] = vertex_cov[0];
      covxin[0][1] = vertex_cov[1];
      covxin[1][1] = vertex_cov[2];
      covxin[0][2] = vertex_cov[3];
      covxin[1][2] = vertex_cov[4];
      covxin[2][2] = vertex_cov[5];
      
      sumofweights =  finder.findVertex(xin,covxin);
    } else {
      sumofweights =  finder.findVertex();
    }
    const VertexFitter* fitter = finder.fitter();
    if(sumofweights == 0) {     
      for(int i = 0 ; i < 3 ; ++i ) vertex[i] = 0;
      for(int i = 0 ; i < 6 ; ++i ) vertex_cov[i] = 0;
      for(int i = 0 ; i < ntracks ; ++i) {
	constrained_track_par[i][0] = 0;
	constrained_track_par[i][1] = 0;
	constrained_track_par[i][2] = 0;
	constrained_track_p[i] = 0;
	weights[i] = 0;
      }
      ndof = 0;
      chi2 = 9999999;
      mass = 0;
      finder.reset();
      return ;
    }
    vertex[0] = fitter->vertex()[0];
    vertex[1] = fitter->vertex()[1];
    vertex[2] = fitter->vertex()[2];
    vertex_cov[0] = fitter->covariance()[0][0];
    vertex_cov[1] = fitter->covariance()[0][1];
    vertex_cov[2] = fitter->covariance()[1][1];
    vertex_cov[3] = fitter->covariance()[0][2];
    vertex_cov[4] = fitter->covariance()[1][2];
    vertex_cov[5] = fitter->covariance()[2][2];
    chi2 = fitter->chi2();
    ndof = fitter->ndof();
    HepLorentzVector p(0,0,0,0);
    for(int i = 0 ; i < ntracks ; ++i) {
      const LinearizedTrack* t = finder.track(i);
      if(! t->expanded()) {
	assert(t->weight() == 0);
	constrained_track_par[i][0] = 0;//theta
	constrained_track_par[i][1] = 0;//phi
	constrained_track_par[i][2] = t->charge() * 999999;//Q/R
	weights[i] = 0;
	continue;
      }
      constrained_track_par[i][0] = M_PI_2 - atan(t->q()[2]);//theta
      constrained_track_par[i][1] = (t->q()[0] < 0) ? 2 * M_PI + t->q()[0] : t->q()[0];//phi
      constrained_track_par[i][2] = t->q()[1];//Q/R
      weights[i] = t->weight();   
      HepLorentzVector tp = t->momentum(track_p[i]);
      constrained_track_p[i] = tp.v().mag();
      p +=  t->weight() * tp;
    }    
    mass = p.m();
    finder.reset();
  }

  void vxlitedaffit2d_(int& ntracks,float track_par[VXLITE_MAXTRACKS][5],
		       float track_cov[VXLITE_MAXTRACKS][15],float track_p[VXLITE_MAXTRACKS],
		       float vertex[3],float vertex_cov[6],float& chi2_cut,int& mode,float& chi2,
		       float& ndof,float weights[VXLITE_MAXTRACKS],
		       float constrained_track_par[VXLITE_MAXTRACKS][3], float constrained_track_p[VXLITE_MAXTRACKS], float& mass )
  {
     static DAFVertexFinder finder;
    
    finder.setChi2Cut(chi2_cut);
    for(int i = 0 ; i < ntracks ; ++i) {
      finder.addTrack(new LinearizedTrack(i+1,track_par[i],track_cov[i],LinearizedTrack::PI_PLUS_MASS,
				  true));
    }    
    
    double sumofweights = 0;
    
    if(mode == 1) {
      //with vertex constraint
      HepVector xin(3);
      xin[0] =  vertex[0];
      xin[1] =  vertex[1];
      xin[2] =  vertex[2];
      HepSymMatrix covxin(3);
      covxin[0][0] = vertex_cov[0];
      covxin[0][1] = vertex_cov[1];
      covxin[1][1] = vertex_cov[2];
      covxin[0][2] = vertex_cov[3];
      covxin[1][2] = vertex_cov[4];
      covxin[2][2] = vertex_cov[5];
      
      sumofweights =  finder.findVertex(xin,covxin);
    } else {
      sumofweights =  finder.findVertex();
    }
    const VertexFitter* fitter = finder.fitter();
    if(sumofweights == 0) {     
      for(int i = 0 ; i < 3 ; ++i ) vertex[i] = 0;
      for(int i = 0 ; i < 6 ; ++i ) vertex_cov[i] = 0;
      for(int i = 0 ; i < ntracks ; ++i) {
	constrained_track_par[i][0] = 0;
	constrained_track_par[i][1] = 0;
	constrained_track_par[i][2] = 0;
	constrained_track_p[i] = 0;
	weights[i] = 0;
      }
      ndof = 0;
      chi2 = 9999999;
      mass = 0;
      finder.reset();
      return ;
    }
    vertex[0] = fitter->vertex()[0];
    vertex[1] = fitter->vertex()[1];
    vertex[2] = fitter->vertex()[2];
    vertex_cov[0] = fitter->covariance()[0][0];
    vertex_cov[1] = fitter->covariance()[0][1];
    vertex_cov[2] = fitter->covariance()[1][1];
    vertex_cov[3] = fitter->covariance()[0][2];
    vertex_cov[4] = fitter->covariance()[1][2];
    vertex_cov[5] = fitter->covariance()[2][2];
    chi2 = fitter->chi2();
    ndof = sumofweights - 2;
    if(mode == 1) ndof += 2;
    if(ndof <= 0) {
      for(int i = 0 ; i < 3 ; ++i ) vertex[i] = 0;
      for(int i = 0 ; i < 6 ; ++i ) vertex_cov[i] = 0;
      for(int i = 0 ; i < ntracks ; ++i) {
	constrained_track_par[i][0] = 0;
	constrained_track_par[i][1] = 0;
	constrained_track_par[i][2] = 0;
	constrained_track_p[i] = 0;
	weights[i] = 0;
      }
      ndof = 0;
      chi2 = 9999999;
      mass = 0;
      finder.reset();
      return ;
    }
    HepLorentzVector p(0,0,0,0);
    for(int i = 0 ; i < ntracks ; ++i) {
      const LinearizedTrack* t = finder.track(i);
      if(! t->expanded()) {
	assert(t->weight() == 0);
	constrained_track_par[i][0] = 0;//theta
	constrained_track_par[i][1] = 0;//phi
	constrained_track_par[i][2] = t->charge() * 999999;//Q/R
	weights[i] = 0;
	continue;
      }
      constrained_track_par[i][0] = M_PI_2 - atan(track_par[i][4]);//theta
      constrained_track_par[i][1] = (t->q()[0] < 0) ? 2 * M_PI + t->q()[0] : t->q()[0];//phi
      constrained_track_par[i][2] = t->q()[1];//Q/R
      weights[i] = t->weight();  
      double oldpt = track_p[i]  * std::abs(sin(atan(1/track_par[i][4])));
      double pt = oldpt * track_par[i][1] / t->q()[1];
      double E = sqrt(LinearizedTrack::PI_PLUS_MASS * LinearizedTrack::PI_PLUS_MASS + 
		      pt * pt * (1 + track_par[i][4] * track_par[i][4]));
      HepLorentzVector tp(cos(t->q()[0]) * pt,sin(t->q()[0]) * pt,t->q()[2] * pt,E);  
      constrained_track_p[i] = tp.v().mag();
      p +=  t->weight() * tp;
    }    
    mass = p.m();
    finder.reset();
  }
}

DAFVertexFinder::DAFVertexFinder(double chi_cut)
  : m_chi2_cut(chi_cut)
{
  m_fitter = new VertexFitter(); 
  m_annealing_temps.push_back(100);
  m_annealing_temps.push_back(25);
  m_annealing_temps.push_back(9);
  m_annealing_temps.push_back(4);
  //m_annealing_temps.push_back(2);
  //m_annealing_temps.push_back(1);
  //m_annealing_temps.push_back(0.25); 
}

DAFVertexFinder::DAFVertexFinder(const std::vector<double>& temps,
				 double chi_cut)
  : m_chi2_cut(chi_cut),m_annealing_temps(temps.begin(),temps.end())
{
  m_fitter = new VertexFitter();
}


DAFVertexFinder::~DAFVertexFinder()
{
  reset();
  delete m_fitter;
}
  
void DAFVertexFinder::reset()
{ 
  for(unsigned int k = 0 ; k < m_tracks.size() ; ++k) {
    delete m_tracks[k];
  }
  m_tracks.clear();
}

void DAFVertexFinder::show() const
{
  if(m_tracks.empty()) return;

  m_fitter->show();
  std::cout << "Tracks: \n";
  for(unsigned int i = 0 ; i < m_tracks.size() ; ++i) {
    std::cout << "Id: " << m_tracks[i]->id() << " weight:" 
	      <<  m_tracks[i]->weight() << "\n";
  }
}

double DAFVertexFinder::findVertex(const HepVector& xin)
{
  if(m_tracks.size() < 2) {
    return 0;
  }
  HepVector vertex = xin;

  m_fitter->setMaxIterations(2);
  m_fitter->fit(m_tracks.begin(),m_tracks.end(),vertex);
  if(m_tracks.size() < 3) {
    return m_fitter->nTracks();
  }
  m_fitter->setMaxIterations(1);
  
  double* w = new double[m_tracks.size()];
  for(unsigned int i = 0 ; (i < m_annealing_temps.size()) && (m_fitter->nTracks() > 0) ; ++i) {
    //m_fitter->show();
    double weight_sum = 0;
    double beta = m_annealing_temps[i];
    for(unsigned int j = 0 ; j < m_tracks.size() ; ++j) {
      double chi2;
      if(m_tracks[j]->weight() > 0) {
	chi2 = m_tracks[j]->chi2();
      } else {
	m_tracks[j]->setWeight(1);
	chi2 = m_fitter->chi2ForTrack(m_tracks[j],true);
      }
      w[j] = 1/(1 + exp((chi2 - m_chi2_cut) *0.5/beta));
      //std::cout << "chi2:" << chi2 << "  weight:" << w[j] << "\n";
      weight_sum += w[j];
    }
    if(weight_sum < 1.5) break;
    for(unsigned int j = 0 ; j < m_tracks.size() ; ++j) {
       m_tracks[j]->setWeight(w[j]);
    }
    if(i == m_annealing_temps.size() - 1) {
      m_fitter->setMaxIterations(2);
    }
    vertex = m_fitter->vertex();
    m_fitter->fit(m_tracks.begin(),m_tracks.end(),vertex);
   }
  delete[] w;
  return m_fitter->nTracks();
}

double DAFVertexFinder::findVertex(const HepVector& vtx,
				   const HepSymMatrix& cov)
{
  if(m_tracks.size() < 1) return 0;

  HepVector vertex = vtx;
  m_fitter->setMaxIterations(2);
  m_fitter->constrained_fit(m_tracks.begin(),m_tracks.end(),vertex,vtx,cov);
  m_fitter->setMaxIterations(1);
  for(unsigned int i = 0 ; (i < m_annealing_temps.size()) && (m_fitter->nTracks() > 0) ; ++i) {
    double weight_sum = 0;
    double beta = m_annealing_temps[i];
    for(unsigned int j = 0 ; j < m_tracks.size() ; ++j) {
      double chi2;
      if(m_tracks[j]->weight() > 0) {
	chi2 = m_tracks[j]->chi2();
      } else {
	m_tracks[j]->setWeight(1);
	chi2 = m_fitter->chi2ForTrack(m_tracks[j],true);
      }
      double w = 1/(1 + exp((chi2 - m_chi2_cut) *0.5/beta));
      m_tracks[j]->setWeight(w);
      weight_sum += w;
    }
    if(i == m_annealing_temps.size() - 1) {
      m_fitter->setMaxIterations(2);
    }
    vertex = m_fitter->vertex();
    m_fitter->constrained_fit(m_tracks.begin(),m_tracks.end(),vertex,vtx,cov);
    if(weight_sum < 1.0) break;
  }
  return m_fitter->nTracks();
}

