//*****************************************************************************
// File:      LinearizedTrack.cc
// ----------------------------------------------------------------------------
//=============================================================================
// RCS Current Revision Record
//-----------------------------------------------------------------------------
// $Source: /afs/desy.de/user/s/stadie/zeus/cvsroot/tlite/src/vxlite/LinearizedTrack.cc,v $
// $Revision: 1.10 $
// $Date: 2008/04/07 11:57:20 $
// $Author: stadie $
// $State: Exp $
// $Locker:  $
//*****************************************************************************
//=============================================================================
// Declarations and Definitions
//=============================================================================
#include "LinearizedTrack.hh"

#include <cassert>
#include <iostream> 

using namespace VxLite;
using CLHEP::HepVector;
using CLHEP::HepMatrix;
using CLHEP::HepSymMatrix;
using CLHEP::HepLorentzVector;

const double LinearizedTrack::CURVATURE_CONSTANT=0.0029979;
const double LinearizedTrack::BFIELD_TESLA=1.4324;
const double LinearizedTrack::PI_PLUS_MASS=0.13957018;

LinearizedTrack::LinearizedTrack(int track_id, const HepVector &parameters,
		  const HepSymMatrix &covariance,double part_mass) 
  : m_id(track_id), m_alpha(parameters), m_covalpha(covariance),m_mass(part_mass),
    m_weight(1.0),m_phi(parameters[0]),m_x(3,0),m_q(3,0),m_h(5),m_A(5,3),
    m_B(5,3),m_expanded(false),m_chi2(0),m_scale(1.0)
{
  m_charge = (m_alpha[1] > 0) ? 1 : -1 ;  
  int inversionerror;
  m_G.assign(m_covalpha);
  m_G.invertCholesky5(inversionerror);   
  if(inversionerror) {
    std::cout << "VxLite::LinearizedTrack::LinearizedTrack:  ignored track id:"
	      << track_id << ". covariance matrix not positive definite.\n";
    m_weight = 0;
    m_G = HepSymMatrix(5,0);
  } 
}

LinearizedTrack::LinearizedTrack(int track_id,const float* const parameters,
				 const float* const covariance,
				 double part_mass, bool twoD) 
  : m_id(track_id),m_alpha(5),m_covalpha(5),m_mass(part_mass),
    m_weight(1.0),m_x(3,0),m_q(3,0),m_h(5),m_A(5,3),m_B(5,3),m_expanded(false),
    m_chi2(0),m_scale(1.0)
{
  for(int i = 0, k = 0 ; i < 5 ; ++i) {
    m_alpha[i] = parameters[i];
    for(int j = i ; j < 5 ; ++j) {
      m_covalpha[i][j] = covariance[k]; 
      m_covalpha[j][i] = covariance[k];
      ++k;
    }
  }
  if(twoD) {    
    m_alpha[3] = 0;
    m_alpha[4] = 0;
    for(int i = 0 ; i < 5 ; ++i) {
      for(int j = i +1 ; j < 5 ; ++j) {
	if((i == 3) ||(j == 3) || (i==4)||(j==4)) {
	  m_covalpha[i][j] = 0; 
	  m_covalpha[j][i] = 0;
	}
      }
    }
  }
  if(m_alpha[0] > M_PI) m_alpha[0] -= 2 * M_PI;
  m_phi =  m_alpha[0];
  m_charge = (m_alpha[1] > 0) ? 1 : -1 ;
  int inversionerror;
  m_G.assign(m_covalpha);
  m_G.invertCholesky5(inversionerror);  
  if(inversionerror) {    
    std::cout << "VxLite::LinearizedTrack::LinearizedTrack:  ignored track id:"
	      << track_id << ". covariance matrix not positive definite.\n";
    m_weight = 0;
    m_G = HepSymMatrix(5,0);
  } 
}
  
void LinearizedTrack::show() const 
{
  std::cout << "Track id=" << m_id << " p=(" << m_alpha[0] << ", "
	    <<  m_alpha[1] << ", " <<  m_alpha[2] << ", " <<  m_alpha[3] 
	    << ", " <<  m_alpha[4] << ")  phi = " << m_phi << "\n";
  std::cout << "expansion point x=(" << m_x[0] << ", " << m_x[1] << ", " 
	    << m_x[2] << ")  q=(" << m_q[0] << ", " << m_q[1] << ", " 
	    << m_q[2] << ")\n";
  if(m_expanded) {
    std::cout << "A: " << m_A;
    std::cout << "B: " << m_B;
    std::cout << "AxBq:" << AxBq();
  }
}

void LinearizedTrack::expandAround(const HepVector& point)
{ 
  if(m_expanded && (point == m_x)) {
    //std::cout << "no new expansion for point " << point << "\n";
    return;
  }
  m_phi = findPhi2D(point);
  m_q[0] = m_phi;
  m_q[1] = m_alpha[1];
  m_q[2] = m_alpha[4];
  m_x[0] = point[0];
  m_x[1] = point[1];
  m_x[2] = point[2];
  
  //  m_x[0] = m_alpha[2] * sin(m_alpha[0]) + 
  // 1/m_alpha[1] * (-sin(m_phi) + sin(m_alpha[0]));
  //m_x[1] = -m_alpha[2] * cos(m_alpha[0]) + 
  //  1/m_alpha[1] * (cos(m_phi) - cos(m_alpha[0]));
  //m_x[2] = m_alpha[3] - m_alpha[4]/m_alpha[1] *(m_phi - m_alpha[0]);
  //std::cout << "expansion point: " << point << "\n";
  //std::cout << "point on helix: " << m_x << "\n";
  //std::cout << "new q:" << m_q[0] << ", " << m_q[1] << ", " << m_q[2] << "\n";
  m_expanded = computeAandB();
}

void LinearizedTrack::expandAround(const HepVector& x,
				   const HepVector& q)
{ 
  if(m_expanded && (x == m_x) && (q == m_q)) {
    //std::cout << "no new expansion for point " << point << "\n";
    return;
  }
  m_x = x;
  m_q = q;
  m_expanded = computeAandB();
}

bool LinearizedTrack::computeAandB()
{
  const double cos_phi = cos(m_q[0]);
  const double sin_phi = sin(m_q[0]);
  const double C = m_q[1];
  const double OneOverC = 1/C;
  const double a = 1/(1+ C * C * (m_x[0] * m_x[0] + m_x[1] * m_x[1]) 
		      + 2 * C * (m_x[0]*sin_phi - m_x[1]*cos_phi));
  const double b = sqrt(a);
  //phi1 - phi0= ArcTan[C(xCos[phi1] + ySin[phi1])/(-1+CyCos[phi1]-CxSin[phi1])
  const double gamma = atan((m_x[0] * cos_phi + m_x[1] * sin_phi)/
			    (-OneOverC + m_x[1] * cos_phi - m_x[0] * sin_phi));

  //dphi0 / dx
  // C (-C y + Cos[phi1]) /
  // (1 + C^2(x^2 + y^2) - 2Cy Cos[phi1] + 2Cx Sin[phi1])
  m_A[0][0] = C * (-C * m_x[1] + cos_phi) * a;
  //dphi0 / dy
  // C(C x + Sin[phi1]) /
  // (1 + C^2(x^2 + y^2) - 2Cy Cos[phi1] + 2Cx Sin[phi1])
  m_A[0][1] =  C * (C * m_x[0] + sin_phi) * a;
  // dphi0 / dz
  m_A[0][2] = 0;
  //dQRinv / dx
  m_A[1][0] = 0;
  //dQRinv/dy
  m_A[1][1] = 0;
  //dQRinv / dz
  m_A[1][2] = 0;
  //dQDH / dx
  //(C x + Sin[phi1])/
  //@(1 + C^2(x^2 + y^2) - 2Cy Cos[phi1] + 2Cx Sin[phi1])
  m_A[2][0] = (C * m_x[0] + sin_phi) * b;
  //dQDH / dy 
  //(C y - Cos[phi1])/
  //@(1 + C^2(x^2 + y^2) - 2Cy Cos[phi1] + 2Cx Sin[phi1])
  m_A[2][1] = (C * m_x[1] - cos_phi) * b;
  //dQDH / dz
  m_A[2][2] = 0;
  //dZ / dx
  // cottheta(Cy - Cos[phi1])/(1+C^2(x^2+y^2)-2CyCos[phi1]+2CxSin[phi1])
  m_A[3][0] = m_q[2] * ( C * m_x[1] - cos_phi) * a;
  //dZ / dy
  // -cottheta(Cx + Sin[phi1])/(1+C^2(x^2+y^2)-2CyCos[phi1]+2CxSin[phi1])
  m_A[3][1] = -m_q[2] * (C * m_x[0] + sin_phi) * a;
  //dZH / dz
  // 1
  m_A[3][2] = 1;
  //dcottheta / dx
  m_A[4][0] = 0; 
  //dcottheta / dy
  m_A[4][1] = 0; 
  //dcottheta / dz
  m_A[4][2] = 0;


  //dphi0 / dphi1
  //(1 - C y Cos[phi1] + C x Sin[phi1])/
  //(1 + C^2(x^2 + y^2)- 2Cy Cos[phi1] + 2Cx Sin[phi1])
  m_B[0][0] = 1 - C*(C*(m_x[0]*m_x[0]+m_x[1]*m_x[1])-m_x[1]*cos_phi +
		     m_x[0]*sin_phi)*a;
  // dphi0 /dCurvature
  //x Cos[phi1] + y Sin[phi1]/
  //(1 + C^2(x^2 + y^2)- 2Cy Cos[phi1] + 2Cx Sin[phi1])
  m_B[0][1] = (m_x[0] * cos_phi + m_x[1] * sin_phi) * a;
  //dphi1 / dcottheta
  m_B[0][2] = 0;
  //dQRinv / dphi1
  m_B[1][0] = 0;
  //dQRinv / dC
  m_B[1][1] = 1;
  //dQRinv / dcottheta
  m_B[1][2] = 0;
  //dQDH / dphi1
  // (x Cos[phi1] + y Sin[phi1]) /
  //@(1 + C^2(x^2 + y^2) - 2Cy Cos[phi1] + 2Cx Sin[phi1])
  m_B[2][0] = (m_x[0] * cos_phi + m_x[1] * sin_phi) * b;
  //dQDH / dC
  // OneOverC^2 - (1 - C y Cos[phi1] + C x Sin[phi1])/
  //@(1 + C^2(x^2 + y^2) - 2Cy Cos[phi1] + 2Cx Sin[phi1])/C^2
  m_B[2][1] = (1 - (1 + C * (m_x[0] * sin_phi - m_x[1] * cos_phi))*b) * OneOverC * OneOverC;
  //dQDH / dcottheta
  m_B[2][2] = 0;
  //dZH / dphi1 = dgamma/dphi1 cottheta/C = (1- dphi0/dphi1) cottheta/C
  m_B[3][0] = (1-m_B[0][0]) * m_q[2] * OneOverC;
  //dZH / dC
  //cottheta(ArcTan[(C(xCos[phi1]+ySin[phi1])/(1 - CyCos[phi1] + CxSin[phi1])]
  // - C(xCos[phi1] + ySin[phi1])/
  //(1 + C^2(x^2+y^2)- 2CyCos[phi1] + 2CxSin[phi1]))/C^2)
  m_B[3][1] = m_q[2]*(-gamma - C*m_B[0][1]) * OneOverC * OneOverC;
  //dZH / dcottheta
  //ArcTan[(C(XCos[phi0] + YSin[phi0])/(1 + Cd0 + CYCos[phi0] - CXSin[phi0])/C
  m_B[3][2] = gamma * OneOverC;
  //dcottheta / dphi1
  m_B[4][0] = 0;
  //dcottheta / dC
  m_B[4][1] = 0;
  //dcottheta / dcottheta
  m_B[4][2] = 1;


  m_h[0] = m_q[0] - gamma - m_A[0][0] * m_x[0] - m_A[0][1] * m_x[1] - 
    m_B[0][0] * m_q[0] - m_B[0][1] * m_q[1];
  m_h[1] = 0;
  m_h[2] = (-1 + 1/b) * OneOverC - m_A[2][0] * m_x[0] - m_A[2][1] * m_x[1] - 
    m_B[2][0] * m_q[0] - m_B[2][1] * m_q[1];
  m_h[3] = m_q[2] * gamma / C - m_A[3][0] * m_x[0] - m_A[3][1] * m_x[1] -
    m_B[3][0] * m_q[0] - m_B[3][1] * m_q[1] - m_B[3][2] * m_q[2];
  m_h[4] = 0;

  //std::cout << "z0 - z0(x,q):" << m_alpha[3] - m_x[2] - m_q[2] * gamma / C
  //	    << "   gamma:" << gamma << "\n";
  m_E.assign(m_G.similarityT(m_B));
  int errorcode;
  m_E.invert(errorcode);
  if(errorcode) {
    if(m_weight) {
      std::cout << "VxLite::LinearizedTrack::computeAandB(): cannot invert BGB.\n";
      assert(errorcode == 0);
    }
    return false;
  }
  m_EBG = m_E * m_B.T() * m_G;
  m_W.assign(m_G);
  m_W -= m_E.similarity(m_G * m_B);
  return true;
}

double LinearizedTrack::findPhi3D(const HepVector& point)
{
  //use the method by Newton
  double phi = m_alpha[0];
  /*
    \!\(\((cottheta\^2\ errD\^2\ phi - cottheta\^2\ errD\^2\ phiH + 
    cottheta\ errD\^2\ QRinv\ z - cottheta\ errD\^2\ QRinv\ ZH + 
    errZ\^2\ QRinv\ x\ Cos[phi] + errZ\^2\ QRinv\ y\ Sin[phi] + 
    errZ\^2\ Sin[phi - phiH] + 
    errZ\^2\ QDH\ QRinv\ Sin[phi - phiH])\)/\((cottheta\^2\ errD\^2 + 
    errZ\^2\ QRinv\ y\ Cos[phi] + 
    errZ\^2\ \((1 + QDH\ QRinv)\)\ Cos[phi - phiH] - 
    errZ\^2\ QRinv\ x\ Sin[phi])\)\)
  */

  //some constants
  const double cottheta2 = m_alpha[4] * m_alpha[4];
  const double const_numerator = m_alpha[4]* m_alpha[1] * (point[2] - m_alpha[3]);
  const double one_plus_QDH_QRinv = 1 +  m_alpha[2] * m_alpha[1];
  const double errZ2 = m_covalpha[3][3];
  const double errD2 = m_covalpha[2][2];
  
  for(int i = 0 ; i < 5 ; ++i) {
    const double sin_phi = sin(phi);
    const double cos_phi = cos(phi);
    const double delta_phi = phi - m_alpha[0];
    const double numerator =   errD2 * (cottheta2 * delta_phi + const_numerator )
      + errZ2 * ( m_alpha[1] * ( point[0] * cos_phi + point[1] * sin_phi) 
		  + sin(delta_phi) * one_plus_QDH_QRinv);
    const double denominator = errD2 * cottheta2 
      + errZ2 * (  cos(delta_phi) * one_plus_QDH_QRinv
		   + m_alpha[1] * ( point[1] * cos_phi - point[0] * sin_phi));
    
    const double dphi = - numerator/denominator;
    phi += dphi;  
    if(std::abs(dphi) < 0.0001) break;
  }
  while (phi < -2 * M_PI) phi += 2 * M_PI;
  while (phi >  2 * M_PI) phi -= 2 * M_PI;
//   if((phi > 2 * M_PI)||(phi < - 2 * M_PI)) {
//     int n = phi / (2 * M_PI);
//     phi -= n * 2 * M_PI;
//     std::cout << " n:" << n << "  phi:" << phi;
//   }
//   std::cout << "\n";
//   if(std::abs(phi - m_alpha[0]) > std::abs(phi - m_alpha[0] + 2 * M_PI)) {
//     phi += 2 * M_PI;
//   } 
//   if(std::abs(phi - m_alpha[0]) > std::abs(phi - m_alpha[0] - 2 * M_PI)) {
//     phi -= 2 * M_PI;
//   } 
  return phi;
}


double LinearizedTrack::findPhi2D(const HepVector& point)
{  
  //phi1 - phi0= -ArcTan[C(Cos[phi0] x + y Sin[phi0])/
  //                     (1 + C d0+CyCos[phi0]-Cx Sin[phi0])]
  
  const double sin_phi = sin(m_alpha[0]);
  const double cos_phi = cos(m_alpha[0]);
  const double k = 1 + m_alpha[2] * m_alpha[1];
  
  const double gamma = -atan2(m_alpha[1] * (point[0] * cos_phi + point[1] * sin_phi),
			k + m_alpha[1] * (point[1] * cos_phi - point[0] * sin_phi));

  return m_alpha[0] + gamma;
}

HepVector LinearizedTrack::AxBq(const HepVector& u, const HepVector& v) const
{
  return m_A * u + m_B * v;
}
 
HepLorentzVector LinearizedTrack::momentum(const HepVector& vq) const
{
  const double pt = std::abs(m_scale * CURVATURE_CONSTANT * BFIELD_TESLA/vq[1]);
  const double E = sqrt(m_mass * m_mass + pt * pt * (1 + vq[2] * vq[2]));
  return HepLorentzVector(cos(vq[0]) * pt,sin(vq[0]) * pt,vq[2] * pt,E);
}
 
HepLorentzVector LinearizedTrack::momentum(const CLHEP::HepVector& vq, double oldmom) const
{
  double oldpt =  oldmom * std::abs(sin(atan(1/m_alpha[4])));
  double pt = oldpt * m_alpha[1] / vq[1];
  const double E = sqrt(m_mass * m_mass + pt * pt * (1 + vq[2] * vq[2]));
  return HepLorentzVector(cos(vq[0]) * pt,sin(vq[0]) * pt,vq[2] * pt,E);
}

HepMatrix LinearizedTrack::dPdq(const HepVector& vq) const
{
  const double pt = std::abs(m_scale * CURVATURE_CONSTANT * BFIELD_TESLA/vq[1]);
  const double a = pt / vq[1];
  const double b = pt * pt / sqrt(m_mass * m_mass + pt * pt * (1 + vq[2] * vq[2]));
  const double cosphi = cos(vq[0]);
  const double sinphi = sin(vq[0]);
  HepMatrix dpdq(4,3);
  //dP/dphi
  dpdq[0][0] = - sinphi * pt;
  dpdq[1][0] =   cosphi * pt;
  dpdq[2][0] = 0;
  dpdq[3][0] = 0;
  //dP/dC
  dpdq[0][1] = - cosphi * a;
  dpdq[1][1] = - sinphi * a;
  dpdq[2][1] = - vq[2] * a;
  dpdq[3][1] = - b * ( 1 + vq[2] * vq[2]) / vq[1];
  //dP/dcot theta
  dpdq[0][2] = 0;
  dpdq[1][2] = 0;
  dpdq[2][2] = pt;
  dpdq[3][2] = vq[2] * b;
  return dpdq;
}
