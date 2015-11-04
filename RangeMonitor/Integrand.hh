//
//  Function Integrand:  Helper function for SenseMonitor.  Reads an 
//  FSpectrum "psd" and copies [f^{7/3}*psd(f)]^{-1} into a second 
//  FSpectrum "out".
//     author:  Kevin C. Schlaufman (kcs149@psu.edu)
//              Patrick J. Sutton (psutton@ligo.caltech.edu)
//
//  $Id: Integrand.hh 4509 2006-08-16 22:46:45Z patrick.sutton $
//
///////////////////////////////////////////////////////////////////////////
#ifndef Integrand_HH
#define Integrand_HH


#ifndef __CINT__
#include "FSpectrum.hh"
#endif  //------ !def(__CINT__)

 
// This segment of code assumes that the FSpectrum object has already been
// PSDed and calibrated.  This function inverts the FSpectrum and divides
// by f^(7/3).

/** Helper function for SenseMonitor.  
  * Assigns to out the FSpectrum for $(f^{7/3} \times psd)^{-1}$.
  * @memo Computes integrand for range estimate integral.
  * @param out Output power spectrum.
  * @param psd Input power spectrum.
  */ 
void integrand(FSpectrum& out, FSpectrum& psd);

#endif     //  Integrand_HH
