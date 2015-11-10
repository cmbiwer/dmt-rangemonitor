//
//  Function Integrand:  Helper function for SenseMonitor.  Reads an
//  FSpectrum "psd" and copies [f^{7/3}*psd(f)]^{-1} into a second
//  FSpectrum "out".
//     author:  Kevin C. Schlaufman (kcs149@psu.edu)
//
//  $Id: Integrand.cc 4509 2006-08-16 22:46:45Z patrick.sutton $
//
///////////////////////////////////////////////////////////////////////////

#include "Integrand.hh"
#include <math.h>

#include <iostream>
using namespace std;

void integrand(FSpectrum& out, FSpectrum& psd)
{
	//----- Set FSpectrum out equal to FSpectrum psd.
	out = psd;

	//----- Initialize items necessary to prepare the integrand.
	int nsteps = out.getNStep();
	double low_freq = out.getLowFreq();
	double f_step = out.getFStep();
	float *data = new float[nsteps + 1];

	//----- Copy data from out object to data array.
	out.getData(nsteps + 1, data);

	//----- Loop to invert PSDed data.  Also divide by f^(7/3).
	for (int i = 0; i < (nsteps + 1); i++)
	{
		data[i] = pow(data[i],-1);
		data[i] *= (pow(low_freq + (i * f_step),(-7.0/3.0)));
	}

	//----- Send data back to FSpectrum object.
	out.setData(nsteps + 1, data);
	
        delete[] data;
}

void integrand_waveform(FSpectrum& out, FSpectrum& psd, FSpectrum& waveform)
{
        //----- Set FSpectrum out equal to FSpectrum psd.
        out = psd;

        //----- Initialize items necessary to prepare the integrand.
        int nsteps = out.getNStep();
        double low_freq = out.getLowFreq();
        double f_step = out.getFStep();
        float *data = new float[nsteps + 1];

        //----- Copy data from out object to data array.
        out.getData(nsteps + 1, data);

        //----- Get waveform values.
        float* waveform_values = new float[waveform.getNStep() + 1];
        waveform.getData(waveform.getNStep()+1, waveform_values);

        //----- Loop to invert PSDed data.  Also divide by f^(7/3).
//        for (int i = 0; i < (nsteps + 1); i++)
//        {
//                cerr << data[i] << " ";
//                data[i] = pow(data[i],-1);
//                data[i] *= waveform_values[i];
//                cerr << data[i] << " " << waveform_values[i] << endl;
//        }

        //----- Send data back to FSpectrum object.
        out.setData(nsteps + 1, data);

        delete[] data;
}
