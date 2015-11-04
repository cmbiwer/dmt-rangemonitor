//
//  Function range():  Helper function for SenseMonitor.
//  Scales the integral of [f^{7/3} * calibrated_PSD]^{-1} 
//  (calculated by Integrate.hh and SenseMonitor.cc) 
//  by the various dimensionful constants needed to convert 
//  the range to physical units of megaparsecs.
//     author:  Kevin C. Schlaufman (kcs149@psu.edu)
//
//  $Id: Range.hh 5982 2009-09-29 20:11:56Z john.zweizig $
//
///////////////////////////////////////////////////////////////////////////


#ifndef Range_HH
#define Range_HH


/** Helper function for SenseMonitor.  This function takes the integral 
  * $\int_{f\_low}^{f\_high} df [f^{7/3} * calibrated_PSD]^{-1}$ 
  * computed by the integrate function and scales it by the necessary 
  * dimensionful constants to convert it to a range in megaparsecs.
  * @memo Compute range.
  * @param channel Channel AS_Q data taken from (eg: L1:LSC-AS_Q).
  * @param f_7_3 Factor returned by the integrate function for $\int df Cal\_PSD(f)$.
  */
double range(double f_7_3, const char *channel);


//  Helper function for SenseMonitor.  This function takes the integral 
//  $\int_{f\_low}^{f\_high} df [f^{7/3} * calibrated_PSD]^{-1}$ 
//  computed by the integrate function and scales it by the necessary 
//  dimensionful constants to convert it to a range in megaparsecs.
    /* Range formula in megaparsecs is (Finn & Chernoff 1993)  
         range = (3 * 1.84)^{1/3} * [
                   (5 M_chirp^{5/3} c^{1/3} f_{7/3})/(96 \pi^{4/3}\rho_0^2)
                 ]^{1/2} * ARM_LENGTH(nm) / METERS_PER_MEGAPARSEC
       where f_{7/3} is given by Calibrate, Integrand, Integrate in 
       (nm)^{-2}*Hz^{-1/3} so ARM_LENGTH is measured in nm, and where 
       the chirp mass for two 1.4 solar mass objects is 
       M_chirp = 1.219*SOLARMASS*NEWTONS_G/C^2
    */
double range(double f_7_3, const char *channel)
{
	//----- Set physical contants and parameters
	const double SOLARMASS = 1.989E30;
	const double NEWTONS_G = 6.67E-11; 
	const double C = 299792458;
	const double METER_PER_MEGAPARSEC = 3.086E22;
	const double RHO_o = 8;
	const double PI_CONST = 3.141592654;
	double ARM_LENGTH;	//Inteferometer arm length in nanometres!

	//----- Assign proper arm length (in nm) based on channel
        if (channel[0] == 'H' && channel[1] == '2')
                ARM_LENGTH = 2E12;
        else if(channel[0] == 'V' && channel[1] == '1')
                ARM_LENGTH = 3E12;
        else
                ARM_LENGTH = 4E12;

  	double M = (1.219 * SOLARMASS * NEWTONS_G) / pow(C,2);
	double r_o = pow(M,5.0/6.0);
	r_o *= pow(5,0.5) * pow(96.0*pow(PI_CONST,4.0/3.0)*pow(RHO_o,2),-0.5);
	r_o *= pow(C,1.0/6.0) * pow(f_7_3,0.5);
	r_o *= pow((3*1.84),1.0/3.0);
        r_o *= 1.0/METER_PER_MEGAPARSEC*ARM_LENGTH;
	
	return(r_o);  
}

#endif     //  Range_HH
