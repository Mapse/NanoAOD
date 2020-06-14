//*****************************************************************************
// File:      DAFVertexFinder.cc
// ----------------------------------------------------------------------------
//=============================================================================
// RCS Current Revision Record
//-----------------------------------------------------------------------------
// $Source: /afs/desy.de/user/s/stadie/zeus/cvsroot/tlite/src/vxlite/vxlitebs.cc,v $
// $Revision: 1.2 $
// $Date: 2006/12/07 14:41:54 $
// $Author: stadie $
// $State: Exp $
// $Locker:  $
//*****************************************************************************
//=============================================================================
// Declarations and Definitions
//=============================================================================

#include <cassert>
#include <iostream> 


  // function to read the beam line information using bstilt and constructing
  // a primary vertex and its covariance matrix based on that information.
  // This vertex can be used for a beam-constrained vertex fit
extern "C" {
  /*
              ierr        :  error flag (0 = OK,
                                 1 = GAF unavailable
                                 2 = run not in GAF)
  */
  void bstilt_(float pos[3],float sig[3],float err[3],float tilt[2],
	       float tiltError[2], float offset[2], float offsetError[2],
	       int& numEntries,int& ierr);

  void vxlitebs_(float beam_width[2],float vertex[3],float vertex_cov[6],
		 int& ierr)
  {
    float sig[3],err[3],tilt[2],tiltError[2],offset[2],offsetError[2];
    int numEntries;
    //get the beam line
    bstilt_(vertex,sig,err,tilt,tiltError,offset,offsetError,numEntries,ierr);
    if(ierr != 0) {
      return;
    }
    //std::cout << "mean:" << vertex[0] << "," << vertex[1] << "," << vertex[2]
    //	      << " tilt:" << tilt[0] << "," << tilt[1] << "\n";
    //std::cout << "offset:" << offset[0] << "," << offset[1] << "\n";
    // override beam widths from GAF
    if((beam_width[0] > 0) && (beam_width[1] > 0)) {
      sig[0] = beam_width[0];
      sig[1] = beam_width[1];
    }

    //vertex[0] = offset[0] + tilt[0] * vertex[2];
    //vertex[1] = offset[1] + tilt[1] * vertex[2];
    //compute covariance
    vertex_cov[5] = err[2] * err[2] + sig[2] * sig[2];
    vertex_cov[0] = err[0] * err[0] + sig[0] * sig[0] + 
      tilt[0] * tilt[0] * vertex_cov[5] + 
      vertex[2] * vertex[2] * tiltError[0] * tiltError[0];
    vertex_cov[2] = err[1] * err[1] + sig[1] * sig[1] + 
      tilt[1] * tilt[1] * vertex_cov[5] + 
      vertex[2] * vertex[2] * tiltError[1] * tiltError[1];
    vertex_cov[1] = vertex_cov[5] * tilt[0] * tilt[1];
    vertex_cov[3] = vertex_cov[5] * tilt[0];
    vertex_cov[4] = vertex_cov[5] * tilt[1];
  }
}
