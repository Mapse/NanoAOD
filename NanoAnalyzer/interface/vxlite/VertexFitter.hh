#ifndef VERTEXFITTER_HH
#define VERTEXFITTER_HH  

//$Id: VertexFitter.hh,v 1.6 2007/05/02 13:44:12 stadie Exp $)

#include "../Vector.h"
#include "../SymMatrix.h"
#include "../LorentzVector.h"
#include <vector>

namespace VxLite {

  class LinearizedTrack;
  
  class VertexFitter
  {
  public:
    VertexFitter();
    virtual ~VertexFitter() {}
    
    void setMaxIterations(int max_iter) { m_max_iter = max_iter;}
    void setCalcMomentumCov(bool calc_cov) { m_calc_covp = calc_cov;}
    // returns the chi2 of the fit
    double fit(const std::vector<LinearizedTrack*>::iterator& begin,
	       const std::vector<LinearizedTrack*>::iterator& end,
	       const CLHEP::HepVector& point = CLHEP::HepVector(3,0)); 
    double constrained_fit(const std::vector<LinearizedTrack*>::iterator& begin, const std::vector<LinearizedTrack*>::iterator& end,const CLHEP::HepVector& point,const CLHEP::HepVector& xin, const CLHEP::HepSymMatrix& covxin);
    double fit(std::vector<LinearizedTrack*>& tracks,
	       const CLHEP::HepVector& point = CLHEP::HepVector(3,0)) { 
      return fit(tracks.begin(),tracks.end(),point); 
    } 
    const double chi2ForTrack(const LinearizedTrack* const track,
			      bool add = true);
    const CLHEP::HepVector &vertex() const { return m_x;}
    const CLHEP::HepSymMatrix &covariance() const { return m_covx;}  
    const CLHEP::HepLorentzVector &momentum() const { return m_p;}
    const CLHEP::HepSymMatrix &momentum_covariance() const { return m_covp;}  
    double charge() const { return m_charge;}
    double nTracks() const { return m_sumofweights;}
    double ndof() const { return m_ndof;}
    double mass() const { return m_p.m();}
    double mass_error() const;
    double chi2() const { return m_chi2;}
    void show() const;
  private:
    VertexFitter(const VertexFitter &fitter) {}//should never be called!
    CLHEP::HepVector m_x;
    CLHEP::HepSymMatrix m_covx;
    CLHEP::HepLorentzVector m_p;
    CLHEP::HepSymMatrix m_covp;
    double m_chi2;
    int m_iteration;
    double m_charge;
    double m_sumofweights;
    double m_ndof;
    int m_max_iter;
    bool m_calc_covp; 
    CLHEP::HepSymMatrix m_c00inv;
    CLHEP::HepVector m_xsum;

    void calcXSum(const std::vector<LinearizedTrack*>::iterator &begin,
		  const std::vector<LinearizedTrack*>::iterator &end,
		  const CLHEP::HepVector &point);
    void calcChi2(const std::vector<LinearizedTrack*>::iterator &begin,
		  const std::vector<LinearizedTrack*>::iterator &end);    
    void calcP(const std::vector<LinearizedTrack*>::iterator &begin,
	       const std::vector<LinearizedTrack*>::iterator &end);
    void calcPCov(const std::vector<LinearizedTrack*>::iterator &begin,
		  const std::vector<LinearizedTrack*>::iterator &end);
  };
}    
#endif
