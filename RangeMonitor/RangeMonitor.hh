/* -*- mode: c++; c-basic-offset: 4; -*- */
//
//   Class SenseMonitor: Base class for Binary Inspiral Sensitivity Monitor 
//   version: 5.0 (2004.09.24)
//   authors: Kevin C. Schlaufman (kcs149@psu.edu)
//            Patrick J. Sutton (psutton@ligo.caltech.edu)
//
//   $Id: SenseMonitor.hh 7395 2015-06-19 01:47:14Z john.zweizig@LIGO.ORG $
//
///////////////////////////////////////////////////////////////////////////


#ifndef SenseMonitor_HH
#define SenseMonitor_HH

//----- Include standard library stuff.
#include <string>
#include <fstream>
#include <iostream>
#include <stdio.h> 
#include <signal.h>

//----- Include DMT-specific headers.
#include "DatEnv.hh"
#include "FSeries.hh"
#include "FSpectrum.hh"
#include "DecimateBy2.hh"
#include "MonServer.hh"
#include "OperStateCondList.hh"
#include "Time.hh"
#include "Trend.hh"
#include "TSeries.hh"

#include "TrigRslt.hh"
#include "TrigClient.hh"
#include "Segment.hh"

//----- Include monitor-specific headers.
#include "FDEasyCalibrate.hh"
#include "PSD.hh"
#include "RangeSummary.hh"

/** Helper struct for SenseMonitor class. 
  * Holds various parameters for the current 
  * invocation of SenseMonitor.
  * @memo Holds parameters for the current invocation of SenseMonitor.
  */
struct Parameter
{

    /// If true then generate new alpha, beta calibration parameters on-the-fly.
    bool genCalParam;

    /** Number of intervals to divide time-series data into when calculating 
      * the noise PSD.  Must divide the length 'T' of the time 
      * series for the PSD class to be invoked. 
      */
    int num_ints;
    
    /// Length (sec) of data to use for each range estimate.
    double T; 

    /// Highest frequency used in range estimate.
    double h_freq;  

    /// Lowest frequency used in range estimate.
    double l_freq;  

    /// Name of uber xml file holding reference calibration data.
    std::string xml_file_name; 

    /// Name of configuration file holding monitor parameters.
    std::string config_file_name;  

    ///  Interferometer (one of H1, H2, or L1).
    std::string IFO;            

    ///  Monitor channel prefix (default: "SNSM")
    std::string monitor_prefix;            

    /// Channel name.  Holds <IFO>:LSC-AS_Q, where <IFO> is one of H1, H2, or L1.
    std::string mChannel;  

    /** Name of configuration file holding operating state condition (OSC) 
      * definitions for the IFO.
      */
    std::string osc_file_name; 

    /** Operating state condition that must be satisfied for range estimates
      * to be computed.  Must be defined in osc_file_name.
      */
    std::string osc_cond;

    /// Name of window type to use in calculating the noise PSD.
    std::string window_name; 

    /// Monitor name sent to MonServer.
    std::string MonitorName;

    /// Enable segment writing
    bool write_segs;

};


/** Helper struct for SenseMonitor class. 
  * Holds all information on where and what output files will be written.
  * @memo Holds information on output files to be written.
  */
struct Output_Files
{
    //  If true then the monitor will write error messages to file 
    //  <IFO>_SenseMonitor_Error.html (default true).
    bool error_log; 

    //  If true then dump output to local directory
    //  instead of to directory specified by 
    //  environment variable DMTHTMLOUT. 
    bool local; 

    //  If true, output is echoed to the screen (default false).
    bool screen_output;  

    //  If true, SenseMonitor writes to a trend file (default false).
    bool trend;

    //  If true then summary file <IFO>_SenseMonitor_Log.html
    //  is written (default true).
    bool write_log; 

    //  Output stream to file holding calibrated noise spectrum.
    std::ofstream calpsd_file;  

    //  Name of file holding calibrated noise spectrum.
    //  Holds name of calibrated pSD file name.
    std::string calpsd_file_name;     

    /// Name (including path) for DMTViewer history file.
    std::string dmtviewer_file_name; 

    //----- Notes: 
    //      1. By default cout is redirected to the log file,
    //         and cerr to the error file. 
    //      2. dmtviewer_file is used as both input stream 
    //         (in the constructor) and output stream (in the 
    //         destructor) so the stream objects are declared 
    //         locally.
    //      3. Naming conventions:  
    //          file streams typically end with _file. 
    //          file names typically end with _file_name (with/without path). 
    //          file names ending with _file_link do not include a path. 

    /// Directory specified in DMTHTMLOUT, where output files go.
    char *html_dir; 
    
    /// Name (including path) for current log file.
    std::string log_file_name;       

    /// Name (excluding path) for current log file.
    std::string log_file_link;  

    /// Output stream for CumLog.txt file.
    std::ofstream cumlog_file;

    /// Name (excluding path) for cumulative log file.
    std::string cumlog_file_link;  

    /// Name (including path) for cumulative log file.
    std::string cumlog_file_name;  

    /// Name (including path) for error file.
    std::string error_file_name;     

    /// Name (excluding path) for error file
    std::string error_file_link;

    //  Name (including path) for html revolver file.
    std::string revolver_file_name;	

    //  Name (including path) for html summary file.
    std::string summary_file_name;   
};


/** Helper struct for SenseMonitor class. 
  * Holds data related to recent range estimates. 
  * @memo Holds data related to recent range estimates. 
  */
struct Range_Data
{
    /// Current value of calibration parameter alpha.
    double alpha;

    /// Current value of calibration parameter beta.
    double beta; 

    /// Current value of calibration line amplitude.
    double ampl;

    /** Current value of effective detector range,
      * adjusted by calibration-line amplitude.
      */ 
    double range; 

    /// Contains last 5 range estimates. (CAL or UNCAL?)
    double T5_array[5];

    /// Contains last 15 range estimates. (CAL or UNCAL?)
    double T15_array[15];   

    /// Average range over last 5 time steps.
    double T5_distance;

    /// Average range over last 15 time steps.
    double T15_distance;
};


/** The SenseMonitor class is the base class for a DMT monitor of the 
  * same name which produces estimates of the range at which an interferometer  
  * is sensitive with SNR > 8 to the inspiral of a pair of 1.4 solar-mass 
  * neutron stars.
  * @memo SenseMonitor class for DMT range monitor.
  * @author Kevin C. Schlaufman and Patrick J. Sutton (psutton@ligo.caltech.edu)
  * @version 5.0: Last modified 2004.09.24.  
  */ 
class SenseMonitor 
  : public DatEnv 
  , public MonServer 
    // , private TrigClient 
{
  public:

    /** Construct a SenseMonitor object.
      * @memo Constructor.
      */ 
     SenseMonitor(int argc, const char *argv[]); 

    /** Destroy a SenseMonitor object.
      * @memo Destructor.
      */ 
    ~SenseMonitor(void); 

    /** Process data statement.  The monitor reads a specified amount 
      * of time-series data for the LSC-AS_Q channel and computes the average 
      * range to which the interferometer could detect the inspiral of a 
      * binary neutron star system.
      * @memo Process data statement.  
      */
    void ProcessData(void);

    /** Write a brief html page summarizing the monitor and interferometer 
      * status.  The file is named <IFO>_SenseMonitor_Summary.revolver.html, 
      * where <IFO> is one of L1, H1, or H2.  FILE IS WRITTEN TO WHAT DIRECTORY?
      * @memo Write an abbreviated html page summarizing the monitor and interferometer status.
      * @param t Start time of data for which summary is written. 
      * @param Run_Par Struct Parameter holding SenseMonitor parameters.
      * @param R_Dat Struct Range_Data holding data on recent range and calibration.
      */ 
    virtual void RevolverSummary(const Time& t, Parameter Run_Par,
				 Range_Data R_Dat);

    /** Write an html page summarizing the monitor and interferometer 
      * status.  The file is named index.html.
      * FILE IS WRITTEN TO WHAT DIRECTORY?
      * @memo Write an html page summarizing the monitor and interferometer status.
      * @param t Start time of data for which summary is written. 
      * @param Run_Par Struct Parameter holding SenseMonitor parameters.
      * @param R_Dat Struct Range_Data holding data on recent range and calibration.
      */ 
    virtual void Summary(const Time& t, Parameter Run_Par, Range_Data R_Dat);

    /** Read run parameters from specified configuration file.
      * See documentation for syntax.
      * @memo Read run parameters from specified configuration file.
      * @param filename Name of configuration file. 
      */ 
    virtual void Configure(std::string& filename);  

    /** Read run parameters from command-line arguments.  
      * See documentation for syntax.
      * @memo Read run parameters from command-line arguments.
      */ 
    virtual void Configure(int argc, const char *argv[]);

    /** Build up a dmt viewer channel name
     */
    std::string DMTView_Channel(const std::string& txt) const;

    /** Build up a dmt viewer channel name
     */
    std::string Trend_Channel(const std::string& txt) const;

    /** Perform various checks on monitor parameters and set  
      * bool fatal_error to TRUE if any of the following hold:  
      * stride is not an integer multiple of number of averages;
      * low frequency is not less than high frequency;
      * calibration file(s) cannot be opened;
      * operating state condition (OSC) configuration file cannot be opened. 
      * Also prints warnings under various conditions.
      * @memo Perform various checks on monitor parameters.
      * @param Run_Par Struct Parameter holding SenseMonitor parameters.  
      * @param fatal_error Set to TRUE if parameters in Run_Par are not valid.
      */
    virtual void VerifyParameters(Parameter& Run_Par, bool& fatal_error); 

    /** Dump run parameters to specified output stream. 
      * Dumps all members of specified Parameter except dyncal and 
      * config_file_name, plus MaxStride, 
      * frequency of the calibration line tracked (if any), 
      * frequency and amplitude of the simulated calibration line (if any), 
      * plus open-loop gain and sensing function data if dyncal is TRUE.  
      * Output is human-readable, and each line begins with  
      * the # character (except for calbration data).
      * @memo Dump run parameters to specified output stream.
      * @param file Ouput stream.   
      */ 
    virtual void DumpParameters(std::ostream& file); 

    /** Write column headers for range data to specified output stream.  The precise 
      * header information depends on whether the calibration line is tracked.
      * @memo Write column headers for range data to specified output stream.
      * @param out Output stream to which to dump headers.
      */ 
    virtual void DumpReportHeader(std::ostream& out);

    /** Write column headers for range data to specified output stream.  The precise 
      * header information depends on whether the calibration line is tracked.
      * @memo Write column headers for range data to specified output stream.
      * @param out Output stream to which to dump headers.
      * @param frequency Frequency of calibration line to be tracked for calibration adjustments (use zero if no such line).
      */ 
    virtual void DumpReportHeader(std::ostream& out, double frequency);

    /** Write current range and calibration information to the specified 
      * output stream.
      * @memo Write current range and calibration information to the specified output stream.
      * @param out Output stream to which to write information. 
      */
    virtual void ReportResults(std::ostream& out, const Time& t, Parameter& Run_Par);

    /** Act on SIGUSR1 signal (reread configuration file).
      * @memo Act on SIGUSR1 signal.
      */
    //virtual void HandleUsr1();     

    /** Handle signals to monitor.
      * @memo Handle signals to monitor.
      */
    virtual void Attention(void) {
        //----- This call serves data to the DMTViewer. 
    //--- Replaced SIGPOLL with SIGIO Keith Thorne June, 2011
    // on Linux, Solaris SIGPOLL same as SIGIO, Mac OS X only has SIGIO 
	if (testAttnSignal(SIGIO)) MonServer::Attention();
    }

  private:

  //    TrigClient mTrigC;

    //  If true then quit cleanly ASAP.
    bool fatal_error; 

    //  Holds last 12 hrs of calibrated range values.
    double *calvalues;

    //  Holds last 12 hrs of cal-line-ampl values.
    double *amplvalues;	

    //  Holds last 12 hrs of alpha values.
    double *alphavalues;

    //  Holds last 12 hrs of alpha*beta values.
    double *alphabetavalues;	

    //  Holds end time of last stride of data read.  Used to check for gaps in data.
    double prev_endtime;

    //  1000-sec average range reported to 2-week summary file.
    float running_average;

    //  Number of strides contributing to running_avg
    float running_count; 

    //  Total number of strides to be processed.
    int MaxStride; 

    //  Debug flag
    int Debug; 

    //  Number of strides of data missing between last and current stride.
    int missing_strides;

    //  Number of strides already processed.
    int NStride;

    //  Stuff for DQ flags 
    int NDQflag;
    char DQflagname[100][256];
    double DQflagthresh[100];
    int DQactivity[100];


    //  Time since last update of 2-week summary.
    int time_since_update;

    //  Number of time steps in DMTViewer history.
    int tserieslength;

    //  "Easy" calibration object for handling ASQ->strain(nm) conversion.
    FDEasyCalibrate *mCalibrate;

    //  FSpectrum object that holds the integrand of the range integral,  
    //  (f^(7_3) * PSD)^{-1}
    //  See Finn & Chernoff, PRD, 1993.
    FSpectrum f_7_3; 

    // Holds the PSD of the AS_Q data.  It is filled by PSD::generate().
    FSpectrum PSDed_data;  

    // Holds the calibrated amplitude spectrum of the detector noise (PSD^0.5).
    FSpectrum strain_noise_ampl;

    /// GPS time (sec) at which the monitor is started.
    long int starttime; 
    
    //  Monitors lock state of IFO.
    OperStateCondList *mOsclist;

    //  Struct holds metadata on output files.
    Output_Files mOutput; 

    //  Struct (defined above) that keeps all the run parameters in one place.
    Parameter Run_Par;		

    //  Object calculates the power-spectral density
    //  (stored in an FSpectrum object) of the
    //  TSeries data fed to it.  See PSD.hh for 
    //  further explanation.
    PSD *psd;                   

    //  Struct (defined above) that transports all
    //  of the calculated range data in one piece
    //  from the process data method to the summary
    //  method, where the webpage is created.
    Range_Data R_Dat;

    //  Range summary data for Szabi's status web page
    RangeSummary mRangeSumm;	

    //  Struct holds frequency and amplitude of a simulated calibration line.
    //Simulated_Calibration_Line Fake_Line; 

    /** Stores monitor arguments as read from command-line 
      * or configuration file; used in writing summary html files.
      */
    std::string command_line;

    /// Directory where monitor was started.
    std::string pwd_dir;  

    //  Trend object to which trend information is sent.
    Trend mTrend;

    //  History of alpha parameter sent to DMTViewer.
    TSeries *history_alpha; 

    //  History of alpha*beta parameter sent to DMTViewer.
    TSeries *history_alphabeta; 

    //  History of calibration-line amplitude sent to DMTViewer.
    TSeries *history_ampl;  

    //  History of calibration-adjusted range estimates sent to DMTViewer.
    TSeries *history_range;

    //  Time series asq data.
    TSeries ts_asq;

    // Trigger Client
    TrigClient mTC;

    //decimator
    DecimateBy2 mDecimate;
  
};

#endif     //----------  SenseMonitor_HH
