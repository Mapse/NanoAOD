#ifndef DAFVERTEXFINDER_HH
#define DAFVERTEXFINDER_HH  

//($Id: DAFVertexFinder.hh,v 1.4 2007/05/02 13:44:12 stadie Exp $)

#include <vector>
#include "../Vector.h"
#include "../SymMatrix.h"
namespace VxLite {
  class LinearizedTrack;
  class VertexFitter;

  class DAFVertexFinder
  {
  public:
    DAFVertexFinder(double chi_cut = 16.0);
    DAFVertexFinder(const std::vector<double>& temps, double chi_cut = 16.0);
    virtual ~DAFVertexFinder();
    
    void setChi2Cut(double chi_cut) { m_chi2_cut = chi_cut;}
    // find vertex. returns the sum of weights of all tracks
    double findVertex(const CLHEP::HepVector& vtx = CLHEP::HepVector(3,0));
    // find vertex using vertex constraint. 
    //returns the sum of weights of all tracks
    double findVertex(const CLHEP::HepVector& vtx,
		      const CLHEP::HepSymMatrix& cov);
    // add a track
//    void DAFVertexFinder::addTrack(LinearizedTrack* t) { m_tracks.push_back(t);}
    void addTrack(LinearizedTrack* t) { m_tracks.push_back(t);}
    unsigned int number_of_tracks() const { return m_tracks.size();}
    double chi2Cut() const { return m_chi2_cut;};
    // return track from track list
    const VxLite::LinearizedTrack* const track(unsigned int i) const { return m_tracks[i];}
    //reset the track list
    void reset();
    const VertexFitter* const fitter() const { return m_fitter;}
    void show() const;
    
  private:
    DAFVertexFinder(const DAFVertexFinder &finder) {}//should never be called!
    
    VertexFitter *m_fitter;
    double m_chi2_cut;
    double m_T;
    std::vector<LinearizedTrack*> m_tracks;
    std::vector<double> m_annealing_temps;
  };
}
#endif
