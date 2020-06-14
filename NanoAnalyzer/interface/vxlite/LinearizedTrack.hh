#ifndef LINEARIZEDTRACK_HH
#define LINEARIZEDTRACK_HH  

//($Id: LinearizedTrack.hh,v 1.8 2008/04/07 11:57:20 stadie Exp $)

#include "../Matrix.h"
#include "../SymMatrix.h"
#include "../Vector.h"
#include "../LorentzVector.h"

namespace VxLite {
  class LinearizedTrack 
  {
  public:
    static const double CURVATURE_CONSTANT;
    static const double BFIELD_TESLA;
    static const double PI_PLUS_MASS;
    
    LinearizedTrack(int track_id, const CLHEP::HepVector &parameters,
		    const CLHEP::HepSymMatrix &covariance,
		    double part_mass = PI_PLUS_MASS);
    LinearizedTrack(int track_id,const float* const parameters,
		    const float* const covariance,
		    double part_mass = PI_PLUS_MASS, bool twoD = false);
    virtual ~LinearizedTrack() {}
    
    // set x or q. This invalidates the matrices
    void expandAround(const CLHEP::HepVector& point);
    // set x or q. This invalidates the matrices
    void expandAround(const CLHEP::HepVector& x,const CLHEP::HepVector& q);

    // allows to scale the weight matrix
    // keep zero weight for corrupted tracks
    void setWeight(double w) { 
      m_weight = ((m_weight == 0) && (m_G[0][0] == 0)) ? 0 : w;
    }
    // sets chi2 with respect to a vertex
    void setChi2(double chi2) { m_chi2 = chi2;}
    // set the momentum scale
    void setMomentumScale(double scale) { m_scale = scale;}

    //the full 5-by-5 weight matrix of track parameters  
    const CLHEP::HepSymMatrix &G() const {return m_G;}
    //expansion point
    const CLHEP::HepVector &x() const { return m_x;}
    //momentum at the expansion point
    const CLHEP::HepVector &q() const { return m_q;}
    //alpha(x,q) - Ax - Bq
    const CLHEP::HepVector &h() const { return m_h;}
    //partial derivative of parameters with respect to expansion point coordinates
    const CLHEP::HepMatrix &A() const { return m_A;}
    //partial derivative of parameters with respect to momentum coordinates
    const CLHEP::HepMatrix &B() const { return m_B;}
    CLHEP::HepVector AxBq(const CLHEP::HepVector& u, const CLHEP::HepVector& v) const;
    CLHEP::HepVector AxBq() const {return AxBq(m_x,m_q);}
    const CLHEP::HepVector &alpha() const { return m_alpha;}
    const CLHEP::HepSymMatrix &covalpha() const { return m_covalpha;}
    const CLHEP::HepSymMatrix &E() const { return m_E;}
    const CLHEP::HepSymMatrix &W() const { return m_W;}
    const CLHEP::HepMatrix &EBG() const { return m_EBG;}
    CLHEP::HepLorentzVector momentum(const CLHEP::HepVector& vq) const;
    CLHEP::HepLorentzVector momentum(const CLHEP::HepVector& vq, double oldmom) const;
    CLHEP::HepLorentzVector momentum() const { return momentum(m_q);}
    CLHEP::HepLorentzVector momentum(double oldmom) const { return momentum(m_q,oldmom);}
    CLHEP::HepMatrix dPdq(const CLHEP::HepVector& vq) const;
    CLHEP::HepMatrix dPdq() const  { return dPdq(m_q);}
    int charge() const { return m_charge;}
    double mass() const { return m_mass;}
    double weight() const { return m_weight;}
    double chi2() const { return m_chi2;}
    int id() const { return m_id;}
    void show() const;
    bool expanded() const { return m_expanded;}
  private:
    LinearizedTrack(const LinearizedTrack &track) {}//should never be called!
    
    int m_id;
    CLHEP::HepVector m_alpha;
    CLHEP::HepSymMatrix m_covalpha;
    double m_mass;//assumed mass of the particle
    int m_charge;
    double m_weight;
    double m_phi;
    CLHEP::HepVector m_x;//EXPANSION POINT COORDINATES
    CLHEP::HepVector m_q;//MOMENTUM AT THE EXPANSION POINT
    CLHEP::HepVector m_h;
    CLHEP::HepSymMatrix m_G;//the full 5-by-5 weight matrix of track parameters
    CLHEP::HepMatrix m_A;//partial derivative of parameters with respect to expansion point coordinates
    CLHEP::HepMatrix m_B;//partial derivative of parameters with respect to momentum coordinates
    CLHEP::HepSymMatrix m_E;//(BGB)^-1
    CLHEP::HepMatrix m_EBG;//EBG
    CLHEP::HepSymMatrix m_W;//G - GBEBG
    bool m_expanded,m_posupdated;
    double m_chi2,m_scale;

    bool computeAandB();
    double findPhi2D(const CLHEP::HepVector& point); 
    double findPhi3D(const CLHEP::HepVector& point);
  };
}

#endif
