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
#include <fstream>
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

        double dyn_range_fac = 5.9029581035870565e+20;

        ofstream myfile;
        myfile.open("waveform_rangemonitor.txt");
        ofstream psdfile;
        psdfile.open("psd_rangemonitor.txt");
        ofstream integrandfile;
        integrandfile.open("integrand_rangemonitor.txt");

        //----- Loop to invert PSDed data.  Also divide by f^(7/3).
        for (int i = 0; i < (nsteps + 1); i++)
        {
//                data[i] = data[i] * pow(1e-10 / 4000, 2);
                myfile << i*0.25 << " " << waveform_values[i] << endl;
                data[i] *= pow(dyn_range_fac, 2);
                psdfile << i*0.25 << " " << data[i] << endl;
                data[i] = pow(data[i],-1);
//                cerr << data[i] << " ";
                data[i] *= waveform_values[i];
//                cerr << waveform_values[i] << " " << data[i] << endl;
                integrandfile << i*0.25 << " " << data[i] << endl;
        }

        myfile.close();
        psdfile.close();
        integrandfile.close();

        //----- Send data back to FSpectrum object.
        out.setData(nsteps + 1, data);

        delete[] data;
}
