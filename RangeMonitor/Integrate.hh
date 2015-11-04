//
//  Function Integrate:  Helper function for SenseMonitor.  
//  Integrates an FSpectrum over all frequencies.
//     author:  Kevin C. Schlaufman (kcs149@psu.edu)
//              Patrick J. Sutton (psutton@ligo.caltech.edu)
//
//  $Id: Integrate.hh 4509 2006-08-16 22:46:45Z patrick.sutton $
//
///////////////////////////////////////////////////////////////////////////


#ifndef Integrate_HH
#define Integrate_HH


#ifndef __CINT__
#include "FSpectrum.hh"
#endif  //------ !def(__CINT__)
 

/** Frequency-integration helper function for SenseMonitor.  
  * Returns $\int_{f\_low}^{f\_high} df X(f)$ calculated using 
  * Simpson's Rule (if the number of frequency steps in X(f) 
  * is even) or the Trapezoidal Rule (otherwise).
  * @memo Frequency integration function.
  * @param X Function to be integrated.  Assumed to be real.
  * @param f_low Lower limit of integration (default 0).
  * @param f_high Upper limit of integration (default 10000).
  */ 
double integrate(FSpectrum const& X, double f_low = 0.0, double f_high = 10000.0);


//  Returns $\int_{f\_low}^{f\_high} df X(f)$ calculated using 
//  Simpson's Rule (if the number of frequency steps in X(f) 
//  is even) or the Trapezoidal Rule (otherwise).
double integrate(FSpectrum const& X, double f_low, double f_high)
{
	//----- Extract the subsection of the function to integrate over.
        FSpectrum temp = X.extract(f_low, f_high - f_low);

	//----- Find number of steps and frequency step in temp FSpectrum.
	int nsteps = temp.getNStep();
	float *data = new float[nsteps + 1];
	double coeffecient = temp.getFStep();
	double total = 0.0;

	//----- Copy data from PSD to array data.
	temp.getData(nsteps + 1, data);

        //----- Compute integral.
	if((nsteps % 2) == 0) {
		//----- Apply Simpson's Rule
		total += data[0];
		for(int i = 1; i < nsteps; i++) {
			if((i % 2) == 1)
				total += 4 * data[i];
			else
				total += 2 * data[i];
		}
		total += data[nsteps];
		total *= (coeffecient / 3.0);
	}
	else {
		//----- Apply Trapezoidal Rule
		total += data[0];
		for (int i = 1; i < nsteps; i++)
			total += 2 * data[i];
		total += data[nsteps];
		total *= (coeffecient / 2.0);
	}

	//----- Return sum
        delete[] data;
	return(total);
}

#endif     //  Integrate_HH
