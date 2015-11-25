/* -*- mode: c++; c-basic-offset: 4; -*- */
//
//   Binary Inspiral Sensitivity Monitor 
//   authors: Christopher M. Biwer (cmbiwer@syr.edu)
//            Kevin C. Schlaufman (kcs149@psu.edu)
//            Patrick J. Sutton (psutton@ligo.caltech.edu)
//
///////////////////////////////////////////////////////////////////////////


//----- Defined constants:
//---------- Time (sec) in lock between each dump of a calibrated AS_Q PSD.
#define CalPSD_Dump_Period 1000  
//---------- Time (sec) over which to display range histories on DMTViewer.
#define DMTViewer_Period 43200  

//----- Include this first
#include "RangeMonitor.hh"

//----- Include standard headers
#include <stdlib.h>
#include <iomanip>
#include <sstream>
#include <vector>

//----- Include LIGO-specific headers
#include "Dacc.hh"
#include "html/writer.hh"
#include "html/Attrib.hh"
#include "html/align.hh"
#include "html/color.hh"
#include "html/document.hh"
#include "html/font.hh"
#include "html/hline.hh"
#include "html/image.hh"
#include "html/label.hh"
#include "html/link.hh"
#include "html/size.hh"
#include "html/style.hh"
#include "html/table.hh"
#include "html/text.hh"
#include "TrigRslt.hh"
//----- Windows
#include "Hamming.hh"
#include "Hanning.hh"
#include "Blackman.hh"
#include "FlatTop.hh"

//----- Include monitor-specific headers
#include "Integrand.hh"
#include "Integrate.hh"
#include "Range.hh"

//-->  The next three lines are needed if you are going to generate triggers.
//     The descriptive title in PIDTITLE should the monitor function.
#define PIDCVSHDR "$Header: https://redoubt.ligo-wa.caltech.edu/svn/gds/trunk/Monitors/SenseMonitor/SenseMonitor.cc 7395 2015-06-19 01:47:14Z john.zweizig@LIGO.ORG $"
#define PIDTITLE  "Inspiral range monitor"
#include "ProcIdent.hh"

#ifndef __CINT__
//======================================  Generate the main routine.
EXECDAT(SenseMonitor)
#endif               // !def(__CINT__)

using namespace std;

const char *USAGE_INFO =
"\n"
"SenseMonitor: Binary Inspiral Sensitivity Monitor version 5.0 (2004.09.24)"
"\n\n"
"Usage:\n\n"
"./SenseMonitor [optional arguments] <IFO>\n\n"
"where <IFO> is one of H1, H2, or L1 and the optional arguments are:\n"
"                \n"
"-fmax #         Specify maximum frequency to include in range estimate.\n"  
"                Default 1400Hz.\n"
"                \n"
"-fmin #         Specify minimum frequency to include in range estimate.\n"
"                Default  20Hz.\n"
"                \n"
"-h, --help      Print this usage information.\n"
"                \n"
"-local          Local mode:  All output files are dumped in local directory\n"
"                instead of the default directory specified by $DMTHTMLOUT;\n"
"                plain-text log file renamed <GPS_Start_Time>.log.\n"
"                \n"
"-logfile <name> Name plain-text version of log file <name>.  Only works if\n"
"                -local option is also selected.  Default filename is\n"
"                <GPS_start_time>.log (-local mode) or\n"
"                <IFO>_SenseMonitor_CumLog.txt (otherwise).\n"
"                \n"
"-max #          Specify maximum number of strides to process before program\n"
"                exits.\n"
"                \n"
"-n #            NOT FUNCTIONAL.  Currently set to 1.\n"
"                Specify number of sections to divide data into for \n" 
"                calculating PSD.  THIS NUMBER MUST DIVIDE 'stride' EVENLY!\n"
"                Default 15.\n"
"                \n"
"-OSCfile <file> REQUIRED.  Specify the OSC configuration file defining the\n"
"                various operating conditions for the interferometer. \n"
"                (If not specified, SenseMonitor will attempt to open\n"
"                $SENSEMONITOR_OSCCONF/SenseMonitor_LockLoss.conf,\n"
"                then ./SenseMonitor_LockLoss.conf.)  SenseMonitor will exit\n"
"                if no OSC configuration file can be found.\n"
"                \n"
"-OSCcond <cond> Specify the OSC (operating state condition) for which \n"
"                SenseMonitor will compute range estimates.  Default is \n"
"                SV_IFO_UP.  This condition must be defined in the file \n"
"                specified by -OSCfile.\n"
"                \n"
"-procname <name>  Set the online process (service) name\n"
"                \n" 
"-prefix <name>  Set the monitor prefix (SNSM) (currently used for dmtviewer\n"
"                channel names only\n" 
"                \n" 
"-readalphabeta  Use stored calibrations from reference calibration file\n"
"                instead of generating new alpha, beta on-the-fly.      \n"
"-screen         Update and error messages are sent to the screen instead of\n"
"                to the default files <IFO>_SenseMonitor_Log.txt and \n"
"                <IFO>_SenseMonitor_Errors.txt.\n"
"                \n"
"-stride #       Specify length in seconds of data to use for each range \n"
"                estimate.  Default 60.\n"
"                \n"
"-trend          Write range and calibration data to trend files.\n"
"                \n"
"-window <type>  Specify data window type, where <type> is: \n"
"                  hanning for Hanning window,\n"
"                  blackman for Blackman window,\n"
"                  flattop for FlatTop window (NOT RECOMMENDED),\n"
"                  hamming for Hamming window (NOT RECOMMENDED),\n"
"                  square for no window (NOT RECOMMENDED).\n"
"                Default Hanning.\n"
"                \n"
"-writesegs <tf> Enable/disable segment writing \n"
"                \n"
"-xmlcal <file>  REQUIRED. Use dynamical calibration based on a reference \n"
"                calibration file <file>. \n"
"                \n";


//----- Html Text
const char* const kHtmlTitle = "SenseMonitor (version 5.0): Binary Inspiral "
			       "Sensitivity Monitor";
const char* const kHtmlAuthor = "Kevin Schlaufman (kcs149@psu.edu) and "
			        "Patrick Sutton (psutton@ligo.caltech.edu)";
const char* const kHtmlDescription =
   "This monitor calculates the approximate range at which the interferometer "
   "can detect the inspiral of a pair of 1.4 solar-mass neutron stars with "
   "SNR > 8.";
const char* const kHtmlLIGO =
   "Laser Interferometer Gravitational-wave Observatory ";
const char* const kHtmlLIGOLink = "http://www.ligo.caltech.edu";


//======================================  SenseMonitor constructor.
SenseMonitor::SenseMonitor(int argc, const char *argv[]) 
: DatEnv(argc, argv), fatal_error(false), prev_endtime(0.0), 
  running_average(0.0), running_count(0.0), MaxStride(99999), Debug(0), 
  missing_strides(0), NStride(0), NDQflag(0), time_since_update(0), 
  history_alpha(0), history_alphabeta(0), history_ampl(0), history_range(0)
{

  cout << "Entering constructor" << endl;
  fflush(stdout);

  //-------------------  Defined data quality flag names and thresholds
  //-------------------  on standard deviation in inspiral range
  NDQflag = 0;
  sprintf(DQflagname[NDQflag],"DMT-INSPIRAL_RANGE_STDEV_GT_0P50_MPC");
  DQactivity[NDQflag] = 0;
  DQflagthresh[NDQflag++] = 0.50;
  sprintf(DQflagname[NDQflag],"DMT-INSPIRAL_RANGE_STDEV_GT_0P75_MPC");
  DQactivity[NDQflag] = 0;
  DQflagthresh[NDQflag++] = 0.75;
  sprintf(DQflagname[NDQflag],"DMT-INSPIRAL_RANGE_STDEV_GT_1_MPC");
  DQactivity[NDQflag] = 0;
  DQflagthresh[NDQflag++] = 1.00;
  sprintf(DQflagname[NDQflag],"DMT-INSPIRAL_RANGE_STDEV_GT_2_MPC");
  DQactivity[NDQflag] = 0;
  DQflagthresh[NDQflag++] = 2.00;
  sprintf(DQflagname[NDQflag],"DMT-INSPIRAL_RANGE_STDEV_GT_3_MPC");
  DQactivity[NDQflag] = 0;
  DQflagthresh[NDQflag++] = 3.00;
  sprintf(DQflagname[NDQflag],"DMT-INSPIRAL_RANGE_STDEV_GT_4_MPC");
  DQactivity[NDQflag] = 0;
  DQflagthresh[NDQflag++] = 4.00;
  sprintf(DQflagname[NDQflag],"DMT-INSPIRAL_RANGE_STDEV_GT_5_MPC");
  DQactivity[NDQflag] = 0;
  DQflagthresh[NDQflag++] = 5.00;

    //---------- Initialize structures.
    mOutput.error_log = true;
    mOutput.local = false;
    mOutput.screen_output = false;
    mOutput.trend = false;
    mOutput.write_log = true; 
    Run_Par.genCalParam = true;
    Run_Par.T = 60.0;
    Run_Par.l_freq = 20.0;
    Run_Par.h_freq = 1400.0;
    Run_Par.num_ints = 15;
    Run_Par.monitor_prefix = "SNSM";
    Run_Par.write_segs = true;
    R_Dat.alpha = 0.0;
    R_Dat.beta = 0.0;
    R_Dat.ampl = 0.0;
    R_Dat.range = 0.0;
    R_Dat.T5_distance = 0.0;
    R_Dat.T15_distance = 0.0;
    for(int i=0; i<5; i++)  R_Dat.T5_array[i] = 0.0;
    for(int i=0; i<15; i++)  R_Dat.T15_array[i] = 0.0;

    //---------- Look for command-line arguments.
    //           First check if configuration file is being used.
    //           In this case the command line should read 
    //             SenseMonitor -config <filename>
    //           If file <filename> does not exist or cannot be parsed 
    //           then exit with error message.
    //           If not using config file, make sure we have some arguments and that
    //           the last one does not start with a hyphen.
    if ((argc >= 3) && (!strcmp("-config", argv[1]))) {
        Run_Par.config_file_name = argv[2];
        Configure(Run_Par.config_file_name);
    } else {
        Configure(argc, argv); 
    }

    //---------- Send monitor name to MonServer.
    if (Run_Par.MonitorName.empty()) {
        Run_Par.MonitorName = "SenseMonitor_";
	Run_Par.MonitorName += Run_Par.IFO;
	if (Debug >= 1) {
	  Run_Par.MonitorName += "_Debug";
	}
	cout << "Setting monitor name to " << Run_Par.MonitorName << endl;
    }
    MonServer::setServerName(Run_Par.MonitorName.c_str());

    if (Debug >= 2) { 
      cout << "Returned from setServerName" << endl;
    }

    //---------- Do some elementary checks on the parameters.
    VerifyParameters(Run_Par,fatal_error);

    if (Debug > 1) { 
      cout << "Returned from VerifyParameters with fatal_error = " 
	   << fatal_error << endl;
    }

    //---------- Finish setup of monitor.
    if (!fatal_error) {
        //-------- Construct the names of the various output files.  These
        //         are stored in the mOutput struct.
        //-------- First construct directory path portion of the names. 
        //         If running locally, this is trivial.  If not, it is 
        //         $DMTHTMLOUT.  Note that the 
	//         plain-text log file is named in ProcessData() method.

      if (Debug >= 2) { 
	cout << "Inside finish block with mOutput.local =" << mOutput.local << endl;
      }

        if (!mOutput.local) {
	    mOutput.html_dir = getenv("DMTHTMLOUT");

	    if (Debug >= 1) { 
	      cout << "mOutput.html_dir =" << mOutput.html_dir << endl;
	    }
            if (mOutput.html_dir) {
                mOutput.dmtviewer_file_name = string(mOutput.html_dir) + "/";
		mOutput.error_file_name = string(mOutput.html_dir) + "/";
		mOutput.log_file_name = string(mOutput.html_dir) + "/";
		mOutput.revolver_file_name = string(mOutput.html_dir) + "/";
		mOutput.summary_file_name = string(mOutput.html_dir) + "/";

		if (Debug >= 2) { 
		  cout << "Defined basic mOutput stuff" << endl;
		}
	    } else { 
		std::cerr << "SenseMonitor WARNING: Environment variable "
			  << "DMTHTMLOUT not set; switching to\n";
		std::cerr << "                      -local mode (all output "
			  << "files go to the local directory).\n";
		mOutput.local = true;
	    }
	}

	if (Debug >= 2) { 
	  cout << "More mOutput stuff to set up..." << endl;
	}

        //-------- Now add on actual file name.
        mOutput.dmtviewer_file_name += Run_Par.IFO + "_DMTViewer_Data.txt";
        mOutput.error_file_name += Run_Par.IFO + "_SenseMonitor_Errors.txt";
        mOutput.log_file_name += Run_Par.IFO + "_SenseMonitor_Log.txt";
        mOutput.revolver_file_name += Run_Par.IFO;
        mOutput.revolver_file_name += "_SenseMonitor_Summary.revolver.html";
	mOutput.summary_file_name += "index.html";
        //-------- Links to log files.
	mOutput.error_file_link = Run_Par.IFO + "_SenseMonitor_Errors.txt";
	mOutput.log_file_link = Run_Par.IFO + "_SenseMonitor_Log.txt";  
        if (!mOutput.local) mOutput.cumlog_file_link = Run_Par.IFO + 
			      "_SenseMonitor_CumLog.txt"; 

	if (Debug >= 1) { 
	  cout << "Done with setting up mOutput stuff with mOutput.error_log ="
	       << mOutput.error_log << endl;
          cout << "mOutput.error_file_name = " << mOutput.error_file_name 
	       << endl;
          cout << "mOutput.log_file_name = " << mOutput.log_file_name << endl;
	}

        //-------- Divert log and error messages to files, if desired (default).
        if (mOutput.error_log) {
	    if (!freopen(mOutput.error_file_name.c_str(), "w", stderr))
		perror("Error reopening stderr");
	}

        if (mOutput.write_log) {
	    if (!freopen(mOutput.log_file_name.c_str(), "w", stdout))
		perror("Error reopening stdout");
	}

	if (Debug >= 2) { 
	  cout << "About to set up history Tseries" << endl;
	}

        //---------- Set up array that will hold data for DMTViewer history 
        //           TSeries.  Used to insure that fixed length of data is 
        //           displayed.  Initialize data to -1 to indicate SenseMonitor 
        //           was not running previously.
        tserieslength = int( DMTViewer_Period / Run_Par.T);
        calvalues = new double[tserieslength];
        amplvalues = new double[tserieslength];
        alphavalues = new double[tserieslength];
        alphabetavalues = new double[tserieslength];
        //---------- Initialize TSeries using data from previous run, if possible.
        ifstream dmtviewer_data;
        dmtviewer_data.open(mOutput.dmtviewer_file_name.c_str());
        if (dmtviewer_data.fail()) {

	  if (Debug >= 2) { 
	    cout << "In dmtviewer_data.fail true block" << endl;
	  }

            for (int i=0; i<tserieslength; i++) {
                calvalues[i] = -1.0;
                amplvalues[i] = -1.0;
                alphavalues[i] = -1.0;
                alphabetavalues[i] = -1.0;
            }
        } else {

	  if (Debug >= 2) { 
	    cout << "In dmtviewer_data.fail false block" << endl;
	  }
            //:KLUDGE:  This could go wrong if there is a long time gap between runs of program.
            int i=0;
            while (i<tserieslength && dmtviewer_data >> calvalues[i]) { 
                dmtviewer_data >> amplvalues[i];
                dmtviewer_data >> alphavalues[i];
                dmtviewer_data >> alphabetavalues[i];
                ++i;
            }
            //----- Reset latest values to -2 to indicate that histories have been read from a file.
            calvalues[i-1] = -2.0;
            amplvalues[i-1] = -2.0;
            alphavalues[i-1] = -2.0;
            alphabetavalues[i-1] = -2.0;
            //----- Check that our new histories are of the desired length.
            if (i != tserieslength) {
                cerr << "SenseMonitor WARNING:  Stored DMTViewer history from previous run not of ";
                cerr << "required length.  Re-initializing DMTViewer histories.\n";
                for (int i=0; i<tserieslength; i++) {
                    calvalues[i] = -1.0;
                    amplvalues[i] = -1.0;
                    alphavalues[i] = -1.0;
                    alphabetavalues[i] = -1.0;
                }
            }
        }
        dmtviewer_data.close();

	if (Debug >= 2) { 
	  cout << "Done with dmtviewer_data_close()" << endl;
	}

        //-------- Set window type and intialize PSD object.
        if (Run_Par.window_name == "hamming") {
            Hamming ham;
            psd = new PSD(&ham, Run_Par.num_ints);
        } else if (Run_Par.window_name == "hanning") {
            Hanning han;
            psd = new PSD(&han, Run_Par.num_ints);
        } else if (Run_Par.window_name == "flattop") {
            FlatTop flat;
            psd = new PSD(&flat, Run_Par.num_ints);
        } else if (Run_Par.window_name == "blackman") {
            Blackman black;
            psd = new PSD(&black, Run_Par.num_ints);
        } else if (Run_Par.window_name == "square") {
            psd = new PSD("square", Run_Par.num_ints);
        } else {
            Run_Par.window_name = "hanning";
            Hanning han;
            psd = new PSD(&han, Run_Par.num_ints);
        }

        //---------- Set the time stride.
        getDacc().setStride(Interval(Run_Par.T));
	setStrideAlignment(60, 0.0);

	if (Debug >= 1) { 
	  cout << "Time stride is set..." << endl;
	}

	//---------- Get AS_Q channel; decimate to maximum frequency 4096Hz.
	//           (Keep frequencies above 2kHz for the calibrated PSD dumps.)
	mDecimate.setDecimation(1,4);
        getDacc().addChannel(Run_Par.mChannel.c_str());

        //---------- Initialize FDEasyCalibrate object.  
        //           This should be in monitor constructor (ie, before any data 
        //           is requested) because FDEasyCalibrate adds new channels to
        //           the Dacc. 
        //           KLUDGE: Hard-wire requested frequencies to [0,4096]Hz.
        mCalibrate = new FDEasyCalibrate(
                            &getDacc(),
                            Run_Par.xml_file_name.c_str(),
                            Run_Par.genCalParam,
                            0.0,
                            Run_Par.num_ints/Run_Par.T,
			    // Dubious code: What is supposed to be inted?
                            int(4096.0-0.0)*Run_Par.T/Run_Par.num_ints + 1 
                     );

        //----- Check that calibration file applies to channel we are analysing.
        if ( Run_Par.mChannel != mCalibrate->getChannel() ) {
            cerr << "SenseMonitor ERROR: "
		 << "Reference calibration file is not for desired channel." 
		 << endl;
            cerr << "Reference calibration is for: " 
		 << mCalibrate->getChannel() << endl;
            cerr << "SenseMonitor is analysing: " << Run_Par.mChannel << endl;
            cerr << "Exiting.\n";
            cout << "SenseMonitor ERROR: "
		 << "Reference calibration file is not for desired channel." 
		 << endl;
            cout << "Reference calibration is for: " 
		 << mCalibrate->getChannel() << endl;
            cout << "SenseMonitor is analysing: " << Run_Par.mChannel << endl;
            cout << "Exiting.\n";
            exit(707);
        }

	//---------- Initialize OSC & read config file.
        mOsclist = new OperStateCondList(getDacc());
        mOsclist->readConfig(Run_Par.osc_file_name.c_str());
        mOsclist->ignoreAllExcept(Run_Par.IFO);

	if (Debug >= 1) { 
	  cout << "Done with mOsclist stuff" << endl;
	}

	if (mOutput.trend)
	{

	  if (Debug >= 2) { 
	    cout << "Starting mTrend stuff" << endl;
	  }
	    //---------- Intialize Trend frame.
	    mTrend.setIFO(Run_Par.IFO.c_str());
	    mTrend.setName(Run_Par.MonitorName.c_str());
            mTrend.setType(Trend::kMinute);

	    //---------- Add trend channels.  Force calibration-line-related 
	    //           trends even if not tracking a line (we don't want to 
	    //           change the trends that are being written during a science
	    //           run, as this causes the trend writer to crash). 
	    string trendName;
	    trendName = Trend_Channel("EFFECTIVE_RANGE_MPC");
	    mTrend.addChannel(trendName.c_str());
	    trendName = Trend_Channel("CAL_LINE_AMPL_ASQ");
	    mTrend.addChannel(trendName.c_str());
	    //---------- Trends to carry alpha, alpha*beta 
	    //           (names requested by Inspiral group).
	    trendName = Trend_Channel("CAL_CAV_FAC");
	    mTrend.addChannel(trendName.c_str());
	    trendName = Trend_Channel("CAL_OLOOP_FAC");
	    mTrend.addChannel(trendName.c_str());
            //---------- Tell manager to make list of trended channels.
            mTrend.writeIndex();
	}

    }   //---------- !fatal_error

    //---------- Flush any lingering data

    if (Debug >= 1) { 
      cout << "About to flush output from initialization" << endl;
    }

    fflush(stderr);
    fflush(stdout);
}


//======================================  SenseMonitor destructor.
SenseMonitor::~SenseMonitor() 
{
    //---------- No pointers are assigned if a fatal error is registered.
    if (!fatal_error) {
        //---------- Export DMTViewer data to a file to be read if monitor is restarted.
        //           Name of file set elsewhere; uses mOutput.html_dir if defined.
        ofstream dmtviewer_data;
        dmtviewer_data.open(mOutput.dmtviewer_file_name.c_str());
        for (int i=0; i<tserieslength; i++) {
            dmtviewer_data << calvalues[i] << "  ";
            dmtviewer_data << amplvalues[i] << "  ";
            dmtviewer_data << alphavalues[i] << "  ";
            dmtviewer_data << alphabetavalues[i] << "  ";
            dmtviewer_data << endl; 
        } 
        dmtviewer_data.close();

        //----- Pointers assigned in SenseMonitor::SenseMonitor().
        //---------- These pointers are always assigned.
        delete[] calvalues;
        delete[] amplvalues;
        delete[] alphavalues;
        delete[] alphabetavalues;
        delete psd;
        delete mCalibrate;
        delete mOsclist; 

        //----- Pointers assigned in SenseMonitor::ProcessData().
	delete history_alpha;
	delete history_alphabeta; 
	delete history_range;
        delete history_ampl; 

        //----- Close files.
	if (mOutput.trend) mTrend.close();
        mOutput.cumlog_file.close(); 
    }
    //----- Send "finished happily" message.
    cout << "SenseMonitor MESSAGE: SenseMonitor is finished after "
	      << NStride << " strides.\n";
}

//======================================  Build a dmtviewer channel name.
std::string 
SenseMonitor::DMTView_Channel(const std::string& txt) const {
    ostringstream ostr;
    ostr << Run_Par.IFO << Run_Par.monitor_prefix << " " << txt;
    return ostr.str();
};

//======================================  Build a trend channel name.
std::string 
SenseMonitor::Trend_Channel(const std::string& txt) const {
    ostringstream ostr;
    ostr << Run_Par.IFO << ":DMT-" << Run_Par.monitor_prefix << "_" << txt;
    return ostr.str();
};

//======================================  Frame processing function.
void
SenseMonitor::ProcessData(void) {

        if (Debug >= 1) {
	  cout << endl << "Entered ProcessData with NStride = " << NStride << endl;
	}

	//---------- Get current GPS time and check that it is an integer.
        Time time = getDacc().getFillTime();
        if (time.getS() != time.totalS()) {
            cerr << "SenseMonitor WARNING: GPS start time not an integer: ";
            cerr << time.totalS() << endl;
        }

	//---------- Check for gap in data.
        if ( (time.totalS() != prev_endtime) && (NStride != 0) ) { 
          char msg[150];
          sprintf(
              msg,"Missing data between GPS times (%0.3f, %0.3f).\n",
              prev_endtime,time.totalS()
          );
	  cerr << "SenseMonitor ERROR: " << msg;
          missing_strides = (int) 1+(time.totalS()-prev_endtime)/Run_Par.T;
        } else { 
	  missing_strides = 0;
	}
        prev_endtime = time.totalS();
        prev_endtime += Run_Par.T; 

	//---------- Calculate start time of histories to DMTViewer.
	//           Note:  TSeries start time = "time" - "DMTViewer_Period",
        //           where DMTViewer_Period is a defined quantity copied to 
        //           tserieslength. 
	Time history_start = time;
	history_start -= (tserieslength-1)*(Run_Par.T);
        
	//---------- Set cumlog_file_name.  Choices in order of preference are 
	//           (1) "$DMTHTMLOUT/<IFO>_SenseMonitor_CumLog.txt" if 
        //               -local option NOT selected (default).
	//           (2) User's choice if -local and -logfile options selected.
	//           (3) Use <GPS_start_time>.log if only -local selected. 
        if (!mOutput.local) {  
	    mOutput.cumlog_file_name = string(mOutput.html_dir) + "/" 
		                     + Run_Par.IFO + "_SenseMonitor_CumLog.txt";
	}
        if ( (mOutput.cumlog_file_name).empty() && !NStride ) {
            char temp_file_name[50]; 
            sprintf(temp_file_name,"%ld",time.getS());
            mOutput.cumlog_file_name = temp_file_name;
	    mOutput.cumlog_file_name += ".log";
        }

        //---------- Get AS_Q data to TSeries.
        TSeries raw_asq = *getDacc().refData(Run_Par.mChannel);
	if (!mDecimate.isDataValid(raw_asq)) {
	    cerr << "SensMonitor: decimation filter reset after gap" << endl;
	    mDecimate.reset();
	}
	ts_asq = mDecimate(raw_asq);

        FSpectrum waveform;

        //--------- Read waveform file that contains the frequency and
        //          amplitude series to use in the integrand.
        if (not Run_Par.waveform_file_name.empty())
        {
//            cerr << "Reading " << Run_Par.waveform_file_name << endl;
            waveform = ReadWaveform(Run_Par.waveform_file_name, PSDed_data.getLowFreq(), PSDed_data.getFStep());

            double low_freq = waveform.getLowFreq();
            double df = waveform.getFStep();

            //---------- Get values of waveform.
//            float* waveform_values = new float[waveform.getNStep() + 1];
//            waveform.getData(waveform.getNStep()+1, waveform_values);
//            for(size_t i = 0; i < waveform.getNStep()+1; ++i)
//            {
//                cerr << "LINE" << i << " " << low_freq + (i * df) << " " << waveform_values[i] << endl;
//            }

        }

	//---------- Things to do once only.
    	if (NStride == 0) {

	    //---------- Calcuate PSD.
	    psd->generate(PSDed_data, &ts_asq);

	    //---------- Output run parameters.  
	    //----- Set starttime
 	    starttime = time.getS();

	    //----- Open logfile
	    mOutput.cumlog_file.open(mOutput.cumlog_file_name.c_str(), ios::app);

	    //----- Dump parameters to html log file and cout. 
            DumpParameters(cout);
            DumpParameters(mOutput.cumlog_file);

            DumpReportHeader(cout);
            DumpReportHeader(mOutput.cumlog_file);

            if (!mOutput.local) {  
                //----- Set up two-week summary file
	        string summary_name = string(mOutput.html_dir) + "/" + Run_Par.IFO +
		    		           "_summary.txt";
	        mRangeSumm.init(summary_name, time);
	    }

            //---------- Write histories to DMTViewer.
            //----- Initialize and publish history of
	    //      range with calibration line.
            history_range = new TSeries(history_start, Run_Par.T, 
					tserieslength, calvalues);
            string DMTViewer_name = DMTView_Channel("EFFECTIVE RANGE (MPC)");

//            serveData(DMTViewer_name.c_str(), history_range);

            //----- Initialize and publish calibration-line amplitude history.
            history_ampl = new TSeries(history_start, Run_Par.T, 
				       tserieslength, amplvalues);
            DMTViewer_name = DMTView_Channel("CAL LINE AMPL (ASQ)");
//            serveData(DMTViewer_name.c_str(), history_ampl);

            //----- Initialize and publish alpha parameter history.
            history_alpha = new TSeries(history_start, Run_Par.T,
					tserieslength, alphavalues);
            DMTViewer_name = DMTView_Channel("CAL CAV FAC (ALPHA)");
//            serveData(DMTViewer_name.c_str(), history_alpha);

            //----- Initialize and publish beta parameter history.
            history_alphabeta = new TSeries(history_start, Run_Par.T, 
					    tserieslength, alphabetavalues);
            DMTViewer_name = DMTView_Channel("CAL OLOOP FAC (ALPHAxBETA)");
//            serveData(DMTViewer_name.c_str(), history_alphabeta);

	    //---------- Make fake strain_noise_ampl to initialize DMTViewer.
	    strain_noise_ampl = PSDed_data;
	    DMTViewer_name = DMTView_Channel("CALIBRATED NOISE SPECTRUM");
//	    serveData(DMTViewer_name.c_str(), &strain_noise_ampl);

	    //---------- Generate fake f_7_3 for DMTViewer initialization. 
            if ( Run_Par.waveform_file_name.empty() )
            {
    	        integrand(f_7_3, PSDed_data);
            }
            else
            {
                integrand_waveform(f_7_3, PSDed_data, waveform);
            }
	    DMTViewer_name = DMTView_Channel("RANGE INTEGRAND (ARB UNITS)");
//	    serveData(DMTViewer_name.c_str(), &f_7_3);
        }

	//---------- Check that IFO is locked before proceeding.
        if (mOsclist->satisfied(Run_Par.osc_cond.c_str())) {

	    //---------- Update and apply calibration.  Extracts subsection 
	    //           of PSD for which we have calibration info.
            if (Run_Par.genCalParam) {
                //----- Use alpha, beta computed on-the-fly.
                mCalibrate->UpdateResponseFunction();
            } else {
                //----- Use alpha, beta stored in reference calibration file.
                mCalibrate->UpdateResponseFunction(time);
            }

	    //---------- Calcuate PSD... Perform calibration before.
	    psd->generate(PSDed_data, &ts_asq, *mCalibrate);
            //PSDed_data.Dump(cerr);

            //---------- Also retrieve calibration data to report.
            //           MUST BE DONE AFTER "UpdateResponseFunction"!
	    R_Dat.alpha = mCalibrate->GetAlpha();
	    R_Dat.beta = mCalibrate->GetBeta();
	    R_Dat.ampl = mCalibrate->GetLineAmplitudeASQ();

	    //---------- Compute integrand of range integral ([PSD*f^7/3]^-1).
            if ( Run_Par.waveform_file_name.empty() )
            {
                integrand(f_7_3, PSDed_data);
            }
            else
            {
                integrand_waveform(f_7_3, PSDed_data, waveform);
            }

	    //---------- Calculate distance at which detector is sensitive
            R_Dat.range = integrate(f_7_3, Run_Par.l_freq, Run_Par.h_freq);
            if (Run_Par.waveform_file_name.empty())
            {
	        R_Dat.range = range(R_Dat.range, (Run_Par.mChannel).c_str());
            }
            else
            {
                R_Dat.range = range_waveform(R_Dat.range, (Run_Par.mChannel).c_str());
            }

            //---------- Get values of f_7_3.
//            float* f_7_3_values = new float[f_7_3.getNStep() + 1];
//            f_7_3.getData(f_7_3.getNStep()+1, f_7_3_values);
//            for(size_t i = 0; i < f_7_3.getNStep()+1; ++i)
//            {
//                cerr << "LINE" << i << " " << 0 + (i * 0.25) << " " << f_7_3_values[i] << endl;
//            }

	    //---------- Compute strain_noise_ampl from PSDed_data, for DMTViewer.
            //           This code effectively duplicates the Calibrated PSD dump code.
            //           I should remove the latter.

	    //----- Set up array to hold values of calibrated AS_Q PSD.
            float* strainvalues = new float[PSDed_data.getNStep() + 1];

	    //----- Get values of calibrated PSD.
            PSDed_data.getData(PSDed_data.getNStep() + 1, strainvalues);

            //----- Specify arm length in nm for conversion of PSD units to strain:
            float arm_length = 1.0;
	    if (Run_Par.IFO == "L1" || Run_Par.IFO == "H1") { 
                arm_length = 4.0e12;
            } else if (Run_Par.IFO == "V1") {
                arm_length = 3.0e12;
            } else if (Run_Par.IFO == "H2") {
                arm_length = 2.0e12;
            }

            //----- Compute sqrt of power spectrum, converted to units 
	    //      of strain/Hz^1/2.  
            for (unsigned int j=0; j<PSDed_data.getNStep()+1; j++) {
	      strainvalues[j] = sqrt(double(strainvalues[j])) / arm_length;
            }

            //----- Copy data to strain_noise_ampl and clean up.
	    strain_noise_ampl = PSDed_data;
	    strain_noise_ampl.setData(PSDed_data.getNStep() + 1, strainvalues);
	    delete[] strainvalues;

	} else { 

	    //---------- IFO not Locked; set range to zero, don't increment
	    //		 4-volume.
	    R_Dat.alpha = 0.0;
	    R_Dat.beta = 0.0;
      	    R_Dat.range = 0.0;
	    R_Dat.ampl = 0.0;	

	    //----- Set to 1 PSDs to be send to DMTViewer.
	    int nsteps = PSDed_data.getNStep();
	    float *data = new float[nsteps + 1];
	    for (int i=0; i<(nsteps + 1); i++) {
		data[i] = 1;
	    }
	    strain_noise_ampl.setData(nsteps + 1, data);
	    f_7_3.setData(nsteps + 1, data);
	    delete[] data;
      	}

	//---------- Shift DMTViewer TSeries data one time step.
        int nshift = 1+missing_strides;
	if (nshift > tserieslength) nshift = tserieslength;

        for (int ii=0; ii<(tserieslength-nshift); ii++) {
            calvalues[ii] = calvalues[ii+nshift];
            amplvalues[ii] = amplvalues[ii+nshift];
            alphavalues[ii] = alphavalues[ii+nshift];
            alphabetavalues[ii] = alphabetavalues[ii+nshift];
        }

	//---------- Fill any missing strides with -1.0.
	if (missing_strides) {
            for (int ii=(tserieslength-nshift); ii<(tserieslength-1); ii++) {
                calvalues[ii]   = -1.0;
                amplvalues[ii]  = -1.0;
                alphavalues[ii] = -1.0;
                alphabetavalues[ii]  = -1.0;
            } 
	}

	//---------- Store most recent values.
        calvalues[tserieslength-1] = R_Dat.range;
	amplvalues[tserieslength-1] = R_Dat.ampl;
	alphavalues[tserieslength-1] = R_Dat.alpha;
	alphabetavalues[tserieslength-1] = (R_Dat.alpha)*(R_Dat.beta);

	//---------- Export data to DMTViewer. 
	history_range->setData(history_start, Run_Par.T, calvalues, tserieslength);            
	history_ampl->setData(history_start, Run_Par.T, amplvalues, tserieslength);
	history_alpha->setData(history_start, Run_Par.T, alphavalues, tserieslength);
	history_alphabeta->setData(history_start, Run_Par.T, alphabetavalues, tserieslength);

	//---------- Write the range and 4-volume scanned to the log file
	//	     and cout.
	ReportResults(mOutput.cumlog_file, time, Run_Par);
        if (mOutput.screen_output || mOutput.write_log) {
	    ReportResults(cout, time, Run_Par);
        }

        //---------- Compute standard deviation of range over last 10 strides
	int sumcount = 0;
        double sum = 0.;
        double sumsq = 0.;
	double thisvalue;
        for (int istride=1; istride<=10; istride++) {
	  thisvalue = max(calvalues[tserieslength-istride],0.);
          if (thisvalue > 0.) {
	    sumcount++;
	    sum += thisvalue;
	    sumsq += thisvalue*thisvalue;
	  }
	}
	double sumdevsq;
	double stddev_range = 0.;
	if (sumcount>=2) {
	  sumdevsq = sumsq/(sumcount-1) - sum*sum/((sumcount-1)*sumcount);
	  if (sumdevsq > 0.) {
	    stddev_range = sqrt(sumdevsq);
	  }
	}

	if (Debug >= 1) {
	  cout << "Standard deviation of inspiral range = " << stddev_range << " Mpc based on last 10 measurements (most recent first):" << endl;
	  for (int istride=1; istride<=10; istride++) {
	    cout << istride << ": " << calvalues[tserieslength-istride] << "  ; " ;
	  }
	  cout << endl;
	}

	//---------- Set data quality flags if needed
	if (Run_Par.write_segs) {
	    if (Debug >= 2) cout << "About to check the " << NDQflag 
				 << " DQ flags" << endl;
	    int version = 1;
	    Time tBegin = time;
	    Time tEnd = time + (Interval)60;
	    for (int DQflag=0; DQflag < NDQflag; DQflag++) {
		if (Debug >=2) cout << "About to create segment s" << endl;
		trig::Segment s(DQflagname[DQflag],version,tBegin,tEnd);
		if (Debug >=2) cout << "Created segment s" << endl;
		ostringstream comment;
		comment << "Standard deviation of inspiral range exceeds "
			<< DQflagthresh[DQflag]
			<< " Mpc over last 10 minutes.";
		s.setIfos(Run_Par.IFO.c_str());
		if (Debug >=2) cout << "Called s.setIfos" << endl;
		s.setComment(comment.str());
		if (Debug >=2) cout << "Called s.setComment" << endl;
		if (stddev_range>=DQflagthresh[DQflag]) {
		    DQactivity[DQflag] = 1;
		} else {
		    DQactivity[DQflag] = 0;
		}
		s.setActivity(DQactivity[DQflag]);
		if (Debug >=2) cout << "Called s.setActivity" << endl;

	  
		lmsg::error_type rc = mTC.sendSegment(s);
		if (Debug >=2) cout << "Just tried sending a segment..." << endl;
		if (rc) {
		    cerr << DQflagname[DQflag] 
			 << ": Error sending segment: error code "
			 << rc << endl;
		} else if ( Debug >= 1) {
		    cout << DQflagname[DQflag] 
			 << ": Segment successfully sent with activity = " 
			 << DQactivity[DQflag] << endl;
		}
	    }

	}

	//---------- Update Trends.
	if (mOutput.trend) {
	    string trendName;
	    trendName = Trend_Channel("EFFECTIVE_RANGE_MPC");
	    mTrend.trendData(trendName.c_str(), time, R_Dat.range);
	    trendName = Trend_Channel("CAL_LINE_AMPL_ASQ");
	    mTrend.trendData(trendName.c_str(), time, R_Dat.ampl);
	    trendName = Trend_Channel("CAL_CAV_FAC");
	    mTrend.trendData(trendName.c_str(), time, R_Dat.alpha);
	    trendName = Trend_Channel("CAL_OLOOP_FAC");
	    mTrend.trendData(trendName.c_str(), time, (R_Dat.alpha)*(R_Dat.beta));
	}

        //---------- Do periodic range and calibrated PSD dumps if running in 
        //           background (!local) mode:    
	//           -> two-week summary file  
	//           -> calibrated psd dumped to ascii file
	//           NOTE: Since stride "Run_Par.T" may not exactly divide
        //           the specified update period "CalPSD_Dump_Period",
	//	     not all strides will contribute with equal rate
	//	     to a given average.  Eg: If CalPSD_Dump_Period is 1000 
        //           and Run_Par.T = 300, the first through 
	//	     fourth strides contribute to first 1000-sec average with
	//	     weights 1, 1, 1, 1/3.  Fourth stride also contributes to
	//	     second 1000-sec average with weight 2/3.
	//           Note also that report time (GPS truncated to last multiple
        //           of CalPSD_Dump_Period) is 
	//           guaranteed to lie in the 1000-sec interval reported on.
        //           Times when the IFO is not locked contribute zero to the average.
	if (!mOutput.local) {
	    time_since_update += (int) Run_Par.T;
	    if (time_since_update >= CalPSD_Dump_Period) {

	        //---------- 2-week summary file update:

	        //----- Add current range to running average with
	        //      weight = fraction of stride coming before 
	        //     (CalPSD_Dump_Period)-sec mark.
	        running_average += (Run_Par.T - time_since_update 
				    + CalPSD_Dump_Period) /
		   		    Run_Par.T * R_Dat.range;
                running_count += (Run_Par.T - time_since_update 
				  + CalPSD_Dump_Period)/Run_Par.T;
	        running_average /= running_count;

	        //----- Dump running average to two-week summary file.
                //      Report time to last multiple of CalPSD_Dump_Period 
                //      (eg, to last 000 sec for CalPSD_Dump_Period = 1000).
                //cout <<  (time.getS() % CalPSD_Dump_Period) << endl;
                RangeDataSummary tmpstatus = {time.getS() - (time.getS() % CalPSD_Dump_Period), running_average};
                mRangeSumm.append(tmpstatus);
                mRangeSumm.dumpList();

	        //----- Initialize next running average with whatever is left
	        //      of current stride.
	        time_since_update -= CalPSD_Dump_Period;
		running_average = (time_since_update) / Run_Par.T * R_Dat.range;
	        running_count = (time_since_update)/Run_Par.T;

		//----- calibrated psd file dump:  
                //if (mOsclist->satisfied(osc_cond1.c_str())  && 
                //    mOsclist->satisfied(osc_cond2.c_str())  ) {
                if (mOsclist->satisfied(Run_Par.osc_cond.c_str())) {

	            //---------- Set up array to hold calibrated AS_Q PSD.
                    float* CalPSDValues = new float[PSDed_data.getNStep() + 1];
                    char*  time_holder  = new char[128];

		    //---------- Get values of calibrated PSD.
            	    PSDed_data.getData(PSDed_data.getNStep()+1, CalPSDValues);

	    	    //---------- Print values to text file
                    //----- Specify arm length in nm for conversion of PSD units to strain:
                    float arm_length = 1.0;
                    if ((Run_Par.IFO == "L1") || (Run_Par.IFO == "H1")) {
			arm_length = 4.0e12;
                    } else if (Run_Par.IFO == "V1") {
			arm_length = 3.0e12;
                    } else if (Run_Par.IFO == "H2") {
			arm_length = 2.0e12;
                    } else {
			cerr << "SenseMonitor ERROR: IFO " << Run_Par.IFO 
			     << " not recognized.  Calibrated PSD files will"
			     << endl;
			cerr << "                     not have units of strain."
			     << endl;
                    }
//          	    TimeStr(ts_asq.getStartTime(), time_holder, "%Y%M%d.%H:%N:%S.txt");
            	    TimeStr(ts_asq.getStartTime(), time_holder, "%s.txt");
//          	    TimeStr( operator-(ts_asq.getStartTime(),Interval(time_since_update,0)), time_holder, "%s.txt");
       		    mOutput.calpsd_file_name = time_holder;
            	    mOutput.calpsd_file_name = string(mOutput.html_dir) + "/" + Run_Par.IFO + "_CalPSD_GPS_" + mOutput.calpsd_file_name;
            	    mOutput.calpsd_file.open(mOutput.calpsd_file_name.c_str());
            	    mOutput.calpsd_file << "#Calibrated PSD of AS_Q Channel, produced by SenseMonitor version 5.0." << endl;
            	    mOutput.calpsd_file << "#Units: (strain)/Hz^{-0.5}" << endl;
                    mOutput.calpsd_file << "#Number of Averages for PSD: " << Run_Par.num_ints << endl;
                    mOutput.calpsd_file << "#Window: " << Run_Par.window_name << endl;
            	    mOutput.calpsd_file << "#GPS start time of this data: " << time.getS() << endl;
                    int npoints = PSDed_data.getNStep() + 1; 
                    mOutput.calpsd_file << "#Number of points: " << npoints << endl;
            	    mOutput.calpsd_file << "#Frequency step (Hz): " << PSDed_data.getFStep() << endl;
            	    mOutput.calpsd_file << "#Min Frequency (Hz): " << PSDed_data.getLowFreq() << endl;
            	    mOutput.calpsd_file << "#Max Frequency (Hz): " << PSDed_data.getLowFreq() + (npoints-1)*(PSDed_data.getFStep())<< endl;
            	    mOutput.calpsd_file << "#Range (Mpc): " << R_Dat.range << endl;
            	    mOutput.calpsd_file << "#Calibration file comment: " << mCalibrate->getComment() << endl;
            	    mOutput.calpsd_file << "#alpha parameter: " << mCalibrate->GetAlpha() << endl;
            	    mOutput.calpsd_file << "#beta parameter: " << mCalibrate->GetBeta() << endl;

            	    //----- Dump sqrt of power spectrum, converted to units 
		    //      of strain/Hz^1/2.
            	    for (unsigned int j=0; j<PSDed_data.getNStep()+1; j++) {
		      mOutput.calpsd_file << sqrt(double(CalPSDValues[j]))/arm_length << endl;
                    }
            	    mOutput.calpsd_file.close();

                    //---------- Clean up.
                    delete[] CalPSDValues;
                    delete[] time_holder;
	        }
            } else { 

	        //----- Add current range to running average with weight 1.
    	        //running_average += R_Dat.range_uncal;
    	        running_average += R_Dat.range;
    	        running_count   += 1.0;
	    }
	}

	//----------  Calculate average of last five time segments.
	R_Dat.T5_array[NStride % 5] = R_Dat.range;
	for(int i = 0; i < 5; i++)
	    R_Dat.T5_distance += R_Dat.T5_array[i];
	R_Dat.T5_distance /= 5.0;

	//----------  Calculate average of last fifteen time segments.
	R_Dat.T15_array[NStride % 15] = R_Dat.range;
	for(int i = 0; i < 15; i++)
            R_Dat.T15_distance += R_Dat.T15_array[i];
        R_Dat.T15_distance /= 15.0;

	//---------- Print summary webpages.
	RevolverSummary(time, Run_Par, R_Dat);
	Summary(time, Run_Par, R_Dat);

	//---------- Clear old data.
	R_Dat.T15_distance = 0.0;
	R_Dat.T5_distance = 0.0;

	//---------- Check if should finish.
	if (++NStride >= MaxStride) finish();

	//----------  Flush any lingering data.
        fflush(stderr);
	fflush(stdout);
}


//====================================== Write revolver.html page
void
SenseMonitor::RevolverSummary(const Time& t, Parameter Run_Par,
			      Range_Data R_Dat)
{
    //---------- Html document
    html::document doc (kHtmlTitle);
    char buf[128];

    //---------- Title
    html::block titleblk ("center");
    titleblk.lineBreak();
    html::text title (kHtmlTitle);
    html::font fonttitle ("Helvetica");
    title.setFont (fonttitle);
    title.setColor (html::color (0, 0, 128));
    title.setSize (html::size (+4));
    titleblk.addObject (title);
    titleblk.lineBreak();
    titleblk.addObject (html::text ("&nbsp;"));
    doc.addObject (titleblk);
    doc.addObject (html::hline());

    //---------- Current status section.
    html::block status_blk ("center");
    status_blk.lineBreak();
    html::text status_title ("Current Range");
    status_title.setFont (fonttitle);
    status_title.setSize (html::size (+2));
    status_blk.addObject (status_title);
    status_blk.lineBreak();

    //---------- Result table.
    html::table results;
    results.addColumn ("GPS Time");
    results.addColumn ("15-Stride Average Range (Mpc)");
    results.addColumn ("Next Update Time (UTC)");
    results.addRow();
    results.insertData(0, 0, html::text (TimeStr(t, buf, "%s")));
    results.insertData(0, 1, html::text (R_Dat.T15_distance));
    results.insertData(0, 2, html::text (t + Run_Par.T));
    for (int j=0; j<3; ++j) {
        results.refCell(0, j).setAlign ("right");
    }
    results.setBorder(true);
    status_blk.addObject(results);
    doc.addObject(status_blk);

/*
    //---------- Calibration line info.
    if (!mOutput.screen_output) {  
        if ( (R_Dat.range!=0.0) && (Run_Par.refcal) && (!(calibrate->IsCalibrationGood())) ) {
            //---------- Insert error message link
            html::block errorinfo ("center");
	    errorinfo.lineBreak();
	    html::text GAIN_WARNING ("WARNING:  ");
	    GAIN_WARNING << "Computed calibration parameter alpha is non-physical.  "
                         << "Calibration lines may have died.  Using alpha=1.";
            GAIN_WARNING.setSize(html::size (+1));
            GAIN_WARNING.setColor(html::color ("red"));
	    errorinfo.addObject(GAIN_WARNING);
	    errorinfo.lineBreak();
	    errorinfo.lineBreak();
            errorinfo.lineBreak();
            doc.addObject (errorinfo);
        } else if ( (R_Dat.range != 0.0) && ((R_Dat.gain < 0.5) || (R_Dat.gain > 1.5))  &&  (Cal_Line.freq != 0.0)) {
            //---------- Insert error message link
            html::block errorinfo ("center");
	    errorinfo.lineBreak();
	    html::text GAIN_WARNING ("WARNING:  ");
	    if (R_Dat.gain < 0.5) { 
	        GAIN_WARNING << "Calibration-line amplitude is more than 50% below nominal; "
                             << "calibration lines may have died.";
	    } else if (R_Dat.gain > 1.5) { 
	        GAIN_WARNING << "Calibration-line amplitude is more than 50% above nominal. "
	                     << "SenseMonitor may have been given the wrong nominal amplitude.";
	    } 
            GAIN_WARNING.setSize(html::size (+1));
            GAIN_WARNING.setColor(html::color ("red"));
	    errorinfo.addObject(GAIN_WARNING);
	    errorinfo.lineBreak();
	    errorinfo.lineBreak();
            errorinfo.lineBreak();
            doc.addObject (errorinfo);
	}
    }
*/

    //---------- Last update time.
    doc.addObject (html::hline());
    html::block updateblk ("center");
    html::text lasttime ("This page was last updated: ");
    lasttime << TimeStr (t, buf, "%M %d, %Y at %H:%N:%S");
    updateblk.addObject (lasttime);
    updateblk.lineBreak();
    doc.addObject (updateblk);
    doc.addObject (html::hline());

    //---------- Write web page to file.
    ofstream out (mOutput.revolver_file_name.c_str());
    if (out) {
        html::writer htmlout (out);
        doc.write (htmlout);
    }
}


//====================================== Write summary html page.
void
SenseMonitor::Summary(const Time& t, Parameter Run_Par, Range_Data R_Dat)
{   
    //---------- Html document
    html::document doc (kHtmlTitle);
    char buf[128];

    //---------- Title
    html::block titleblk ("center");
    titleblk.lineBreak();
    html::text title (kHtmlTitle);
    html::font fonttitle ("Helvetica");
    title.setFont (fonttitle);
    title.setColor (html::color (0, 0, 128));
    title.setSize (html::size (+4));
    titleblk.addObject (title);
    titleblk.lineBreak();
    titleblk.addObject (html::text ("&nbsp;"));
    doc.addObject (titleblk);
    doc.addObject (html::hline());
   
    //---------- Short description
    doc.addObject (html::text (kHtmlDescription));
    doc.addObject (html::hline());

    //---------- Last update time.
    html::block firstupdateblk ("center");
    html::text lasttime1 ("This page was last updated: "); 
    lasttime1 << TimeStr (t, buf, "%M %d, %Y at %H:%N:%S %Z");
    firstupdateblk.addObject (lasttime1);
    firstupdateblk.lineBreak();
    html::text lasttime2 ("Next scheduled update time: ");
    lasttime2 << TimeStr (t+Run_Par.T, buf, "%M %d, %Y at %H:%N:%S %Z");
    firstupdateblk.addObject (lasttime2);
    firstupdateblk.lineBreak();
    doc.addObject (firstupdateblk);
    doc.addObject (html::hline());

    //---------- Run parameters
    html::block run_param_blk ("center");
    run_param_blk.lineBreak();
    html::text param_title ("Run Parameters");
    param_title.setFont (fonttitle);
    param_title.setSize (html::size (+2));
    run_param_blk.addObject (param_title);
    run_param_blk.lineBreak();
    //---------- Parameters table.
    html::table params; 
    //----- WARNING:  Must make columns first, then rows.
    params.addColumn("");
    params.addColumn("");
    for (int j=0; j<11; ++j) {
	params.addRow();
    }
    int i=-1;
    params.insertData(++i, 0, html::text ("Channel: "));
    params.insertData(  i, 1, html::text (Run_Par.mChannel));
    params.insertData(++i, 0, html::text ("Stride: "));
    params.insertData(  i, 1, html::text (Run_Par.T));
    params.insertData(  i, 1, html::text (" sec"));
    params.insertData(++i, 0, html::text ("Number of Strides Already "
					  "Processed: "));
    params.insertData(  i, 1, html::text (NStride+1));
    params.insertData(++i, 0, html::text ("Maximum Number of Strides to "
					  "Process: "));
    params.insertData(  i, 1, html::text (MaxStride));
    params.insertData(++i, 0, html::text ("Number of Averages for PSD: "));
    params.insertData(  i, 1, html::text (Run_Par.num_ints));
    params.insertData(++i, 0, html::text ("Window: "));
    params.insertData(  i, 1, html::text (Run_Par.window_name));
    params.insertData(++i, 0, html::text ("Low Frequency : "));
    params.insertData(  i, 1, html::text (Run_Par.l_freq));
    params.insertData(  i, 1, html::text (" Hz"));
    params.insertData(++i, 0, html::text ("High Frequency : "));
    params.insertData(  i, 1, html::text (Run_Par.h_freq));
    params.insertData(  i, 1, html::text (" Hz"));
    params.insertData(++i, 0, html::text ("Monitor Start Time (GPS): "));
    params.insertData(  i, 1, html::text (TimeStr(Time(starttime), buf,
					  "%s")));
    params.insertData(++i, 0, html::text ("(UTC): "));
    params.insertData(  i, 1, html::text (TimeStr(Time(starttime), buf,
					  "%M %d, %Y at %H:%N:%S")));
/*
    //---------- Calibration-Line info.
    if (Cal_Line.freq != 0.0) {
	params.addRow();
	params.addRow();
	params.addRow();
	params.addRow();
	params.insertData(++i, 0, html::text ("Calibration Line Monitored at: "
					      ));
	params.insertData(  i, 1, html::text (Cal_Line.freq));
	params.insertData(  i, 1, html::text (" Hz"));
	params.insertData(++i, 0, html::text ("Reference Amplitude: "));
	params.insertData(  i, 1, html::text (Cal_Line.ampl));
	params.insertData(  i, 1, html::text (" ASQ counts"));
	params.insertData(++i, 0, html::text ("Measured Amplitude: "));
	params.insertData(  i, 1, html::text (R_Dat.ampl));
	params.insertData(  i, 1, html::text (" ASQ counts"));
	params.insertData(++i, 0, html::text ("Current Gain: "));
	params.insertData(  i, 1, html::text (R_Dat.gain));
    }
    //---------- Simulated Calibration Line info.
    if (Fake_Line.freq != 0.0) {
	params.addRow();
	params.addRow();
	html::text warning1 ("Simulated Calibration Line Inserted at: ");
	html::text warning2 ("Amplitude: ");
	warning1.setColor(html::color ("red"));
	warning2.setColor(html::color ("red"));
	params.insertData(++i, 0, warning1);
	params.insertData(  i, 1, html::text (Fake_Line.freq));
	params.insertData(  i, 1, html::text (" Hz"));
	params.insertData(++i, 0, warning2);
	params.insertData(  i, 1, html::text (Fake_Line.ampl));
	params.insertData(  i, 1, html::text (" ASQ counts"));
    }
*/
    //----- WARNING:  It appears that alignments can only be set AFTER the cell
    //	    is filled.
    //      WARNING:  The (html::table).getNRow() method returns 1+(# of rows).
    //	    Beware!
    for (int j=0; j<params.getNRow()-1; ++j) {
        params.refCell(j, 0).setAlign("right");
        params.refCell(j, 1).setAlign ("left");
    }
    run_param_blk.addObject(params);
    doc.addObject (run_param_blk);
   
    //---------- Current status section.
    html::block status_blk ("center");
    status_blk.lineBreak();
    html::text status_title ("Current Range");
    status_title.setFont (fonttitle);
    status_title.setSize (html::size (+2));
    status_blk.addObject (status_title);
    status_blk.lineBreak();
   
    //---------- Result table.
    html::table results;
    results.addColumn ("GPS Time");
    results.addColumn ("1-Stride Range (Mpc)");
    results.addColumn ("5-Stride Average Range (Mpc)");
    results.addColumn ("15-Stride Average Range (Mpc)");
    results.addColumn ("Next Update Time (UTC)");
    results.addRow();
    results.insertData(0, 0, html::text (TimeStr(t, buf, "%s")));
    //results.insertData(0, 1, html::text (R_Dat.range_uncal));
    results.insertData(0, 1, html::text (R_Dat.range));
    results.insertData(0, 2, html::text (R_Dat.T5_distance));
    results.insertData(0, 3, html::text (R_Dat.T15_distance));
    results.insertData(0, 4, html::text (t + Run_Par.T));
    for (int j=0; j<5; ++j) {
        results.refCell(0, j).setAlign ("right");
    }
    results.setBorder(true);
    status_blk.addObject(results);
    doc.addObject(status_blk);
   
    //---------- DQ flag status section.
    html::block dq_blk ("center");
    dq_blk.lineBreak();
    html::text dq_title ("Data Quality Flags this Stride");
    dq_title.setFont (fonttitle);
    dq_title.setSize (html::size (+2));
    dq_blk.addObject (dq_title);
    dq_blk.lineBreak();
   
    //---------- DQ flag table.
    html::table flags;
    flags.addColumn ("Data quality flag");
    flags.addColumn ("Status for this time stride (0-OFF, 1-ON)");
    for (int DQflag=0; DQflag<NDQflag; DQflag++) {
      flags.addRow();
      flags.insertData(DQflag, 0, html::text (DQflagname[DQflag]));
      flags.insertData(DQflag, 1, html::text (DQactivity[DQflag]));
      flags.refCell(DQflag,0).setAlign("left");
      flags.refCell(DQflag,1).setAlign("center");
    }
    flags.setBorder(true);
    dq_blk.addObject(flags);
    doc.addObject(dq_blk);

    if (!mOutput.screen_output) {  
        //---------- Insert error message link
        html::block errorinfo ("center");

/*
        //---------- Calibration line info.
        if ( (R_Dat.range!=0.0) && (Run_Par.refcal) && (!(calibrate->IsCalibrationGood())) ) {
            //---------- Insert error message link
            html::block errorinfo ("center");
	    errorinfo.lineBreak();
	    html::text GAIN_WARNING ("WARNING:  ");
	    GAIN_WARNING << "Computed calibration parameter alpha is non-physical.  "
                         << "Calibration lines may have died.  Using alpha=1.";
            GAIN_WARNING.setSize(html::size (+1));
            GAIN_WARNING.setColor(html::color ("red"));
	    errorinfo.addObject(GAIN_WARNING);
	    errorinfo.lineBreak();
	    errorinfo.lineBreak();
            errorinfo.lineBreak();
            doc.addObject (errorinfo);
        } else if ( (R_Dat.range != 0.0) && ((R_Dat.gain < 0.5) || (R_Dat.gain > 1.5))  &&  (Cal_Line.freq != 0.0)) {
	    errorinfo.lineBreak();
	    html::text GAIN_WARNING ("WARNING:  ");
	    if (R_Dat.gain < 0.5) { 
	        GAIN_WARNING << "Calibration-line amplitude is more than 50% below nominal; "
                             << "calibration lines may have died.";
	    } else if (R_Dat.gain > 1.5) { 
	        GAIN_WARNING << "Calibration-line amplitude is more than 50% above nominal. "
	                     << "SenseMonitor may have been given the wrong nominal amplitude.";
	    } 
            GAIN_WARNING.setSize(html::size (+1));
            GAIN_WARNING.setColor(html::color ("red"));
	    errorinfo.addObject(GAIN_WARNING);
	    errorinfo.lineBreak();
	    errorinfo.lineBreak();
	}
*/

        //---------- Link to error-message log.
	html::link errormessage ("Error Messages", mOutput.error_file_link.c_str() );
	errorinfo.addObject (errormessage);
        errorinfo.lineBreak();
        doc.addObject (errorinfo);

        //---------- Insert links to current (this run of monitor) and 
        //           cumulative (all runs) log files.
        html::block output ("center");
	html::link outputlink1 ("Monitor Output: more details", (mOutput.log_file_link).c_str() );
        output.addObject (outputlink1);
        output.lineBreak();
        if (!mOutput.local) {
	    //html::link outputlink2 ("Cumulative Monitor Output from all runs", (cumlog_file_link).c_str() );
	    html::link outputlink2 ("Cumulative Monitor Output from all runs", (mOutput.cumlog_file_link).c_str() );
            output.addObject (outputlink2);
            output.lineBreak();
        }
        output.lineBreak();
        doc.addObject (output);
    }

    //---------- Last update time.
    doc.addObject (html::hline());
    html::block updateblk ("center");
    html::text lasttime ("This page was last updated: ");
    lasttime << TimeStr (t, buf, "%M %d, %Y at %H:%N:%S %Z");
    updateblk.addObject (lasttime);
    updateblk.lineBreak();

    //---------- Author info.
    updateblk.addObject( html::text ("This monitor was written by ") );
    updateblk.addObject( html::link ("Patrick J. Sutton",
				     "mailto:psutton@ligo.caltech.edu") );
    updateblk.addObject( html::text (" and ") );
    updateblk.addObject( html::link ("Kevin C. Schlaufman",
			             "mailto:kcs149@psu.edu") );
    doc.addObject (updateblk);
    doc.addObject (html::hline());

    //---------- LIGO logo.
    html::block logoblk ("div");
    logoblk.addAttr ("align", html::align ("right"));
    logoblk.addObject (html::link (kHtmlLIGO, kHtmlLIGOLink));
    html::image logo;
    logo.setSource ("ligo_logo.gif");
    logo.setWidth ("80");
    logoblk.addObject (logo);
    doc.addObject (logoblk);

    //---------- Write web page to file.
    ofstream out (mOutput.summary_file_name.c_str());
    if (out) {
        html::writer htmlout (out);
        doc.write (htmlout);
    }
}


//====================================== Dump run parameters to output stream.
void
SenseMonitor::DumpParameters(ostream& out)
{   
    //----- Dump run parameters.
    out << "# Command Line: " << command_line << endl;
    out << "# Launch Directory: " << pwd_dir <<  endl;
    out << "# Log File: " << mOutput.cumlog_file_name << endl;
    out << "# OSC Configuration File: " << Run_Par.osc_file_name << endl; 
    out << "# Reference Calibration File Name: " << Run_Par.xml_file_name << endl; 
    out << "# Reference Calibration File Comment: " << mCalibrate->getComment() << endl;
    out << "# Minimum physically allowed alpha*beta: " << mCalibrate->GetMinAlphaBeta() << endl;
    out << "# Maximum physically allowed alpha*beta: " << mCalibrate->GetMaxAlphaBeta() << endl;
    out << "# Channel: " << Run_Par.mChannel << endl;
    out << "# Stride: " << Run_Par.T << " sec" << endl;
    out << "# Max Number of Strides: " << MaxStride << " sec" << endl;
    out << "# Debug level: " << Debug << endl;
    out << "# Number of Averages for PSD: " << Run_Par.num_ints << endl;
    out << "# Low Freq.  : " << Run_Par.l_freq << " Hz" << endl;
    out << "# High Freq. : " << Run_Par.h_freq << " Hz" << endl;
    out << "# Window: " << Run_Par.window_name << endl;
    out << "# Monitor Start Time: " << starttime << endl;
    out << "# " <<endl;
}


//====================================== Dump column headers to output.
void
SenseMonitor::DumpReportHeader(ostream& out) 
{
    //----- Dump headers for range data that will follow.
    out << "# 1) Data value of 0 means IFO not locked." << endl;
    out << "# 2) Data values of -1, -2 mean monitor not functioning." << endl;
    out << "# " << endl;
    out << "#                                             CAL CAV FAC   CAL OLOOP FAC" << endl;
    out << "#   Start                      High-Freq.     Calibration     Calibration" << endl;
    out << "# Time of           Range      Line Ampl.       Parameter       Parameter" << endl;
    out << "# Segment           (Mpc)    (ASQ Counts)           alpha      alpha*beta" << endl;
//  out << "#       -               -               -               -               -" << endl;
}


//====================================== Dump column headers to output.
void
SenseMonitor::DumpReportHeader(ostream& out, double freq) 
{

    //----- Dump headers for range data that will follow.
    out << "# Notes:" << endl;
    out << "# 1) Range of 0 means IFO not locked." << endl;
    out << "# 2) Ranges quoted are adjusted by the" << endl;
    out << "#    calibration-line amplitude only if the -dyncal option is used." << endl;
    out << "# " << endl;
}


//====================================== Read configuration file and pass contents 
//                                       to command-line parser.
void
SenseMonitor::Configure(string& file_name) 
{
    int NEWargc = 1;
    char **NEWargv;
    //cout << "Trying to read config file " << file_name << ".\n";
    //---------- Verify that configuration file can be opened.
    ifstream input;
    input.open(file_name.c_str());
    if (input.fail()) {
        cerr << "SenseMonitor ERROR: Cannot open configuration file ";
        cerr << file_name  << ".  Exiting.\n";
        fatal_error = true;
        finish();
    } else {
        cout << "SenseMonitor MESSAGE: Opened configuration file ";
        cout << file_name << endl;
        //---------- Get number of words in configuration file.
        string word;
        while (input >> word) {
            ++NEWargc;
            //cout << word << endl;
        }
        //cout << "NEWargc = " << NEWargc << endl;
    }
    input.close();
    NEWargv = new char*[NEWargc];
    NEWargv[0] = const_cast<char*>("Configuration File: ");
    //cout << NEWargv[0] << endl;
    //---------- Open file and copy contents to char arrays.
    ifstream inp2;
    inp2.open(file_name.c_str());
    if (inp2.fail()) {
        cout << "Yikes!!!!! Can't open " << file_name << endl;
	return;
    }
    string word;
    for (int i=1; i<NEWargc; i++) {
        inp2 >> word;
        NEWargv[i] = new char[80];
        strcpy(NEWargv[i],word.c_str());
        //cout << NEWargv[i] << endl;
    }
    inp2.close();
    Configure(NEWargc,(const char **) NEWargv);

    for (int i=1; i<NEWargc; i++) delete [] NEWargv[i];
    delete[] NEWargv;
    cout << "SenseMonitor MESSAGE: Done \n";
}


//====================================== Parse contents of configuration file or command line.
void
SenseMonitor::Configure(int argc, const char *argv[]) 
{
/*
    cout << "SenseMonitor::Configure(int argc, const char *argv[]) called.\n";
    cout << "argc, argv[] are : " << argc << "  ";
    for (int i=0; i<argc; i++) {
        cout << argv[i] << "  ";
    }
    cout << endl;
*/

    //---------- Get command line, for the record.
    for (int i=0; i < argc; i++) {
        string arg_i = argv[i];
        command_line += arg_i;        
        command_line += "  ";         
    }
    
    //---------- Get local directory, for the record.
    char *tempdir = getenv("PWD");
    if (tempdir) {
	pwd_dir = string(tempdir) + "/";
    } 
    
    //---------- Get Default channel name and IFO.
    Run_Par.IFO = argv[argc - 1];
    Run_Par.mChannel  = Run_Par.IFO; 
    Run_Par.mChannel += ":LSC-AS_Q";

    //----- Set default parameters.
    Run_Par.osc_cond = Run_Par.IFO + ":SV_IFO_UP";

    //----- Parse command-line arguments, and compare to list of known 
    //      values.  Ignore unrecognized options.
    if ( (argc <= 1) || (argv[argc-1][0] == '-') ) {
        cout << USAGE_INFO << endl;
        fatal_error = true;
        finish();
    } else {
        for (int i=1 ; i<argc-1 ; i++) {
            string argi = argv[i];
            if (argi == "-max") {
	        MaxStride = strtol(argv[++i], 0, 0);
            } else if (argi == "-debug") {
	        Debug = strtol(argv[++i], 0, 0);
 	    } else if (argi == "-stride") {
	        Run_Par.T = strtod(argv[++i], 0);
	    } else if (argi == "-n") {
	        Run_Par.num_ints = (int)strtod(argv[++i], 0);
	    } else if (argi == "-fmin") {
	        Run_Par.l_freq = strtod(argv[++i], 0);
	    } else if (argi == "-fmax") {
	        Run_Par.h_freq = strtod(argv[++i], 0);
 	    } else if (argi == "-xmlfile") { 
                Run_Par.xml_file_name = argv[++i];
	    } else if (argi == "-logfile") {
	        mOutput.cumlog_file_name = argv[++i];
 	    } else if (argi == "-OSCcond") { 
                Run_Par.osc_cond = Run_Par.IFO + ":" + argv[++i];
 	    } else if (argi == "-OSCfile") { 
                Run_Par.osc_file_name = argv[++i];
 	    } else if (argi == "-procname") { 
                Run_Par.MonitorName = argv[++i];
 	    } else if (argi == "-prefix") { 
                Run_Par.monitor_prefix = argv[++i];
 	    } else if (argi == "-channel") { 
                Run_Par.mChannel = argv[++i];
	    } else if(argi == "-readalphabeta") {
                Run_Par.genCalParam = false; 
	    } else if(argi == "-screen") {
                mOutput.screen_output = true;
     	        mOutput.write_log = false;
 	        mOutput.error_log = false;
	    } else if(argi == "-trend") {
	        mOutput.trend = true;
            } else if (argi == "-waveformfile") {
                Run_Par.waveform_file_name = argv[++i];
  	    } else if(argi == "-window") {
                Run_Par.window_name = argv[++i];
  	        for(unsigned int j=0; j < (Run_Par.window_name).length(); j++) {
                    Run_Par.window_name[j] = tolower(Run_Par.window_name[j]);
	        }
  	    } else if(argi == "-writesegs") {
		if (strtol(argv[++i], 0, 0) == 0) Run_Par.write_segs = false;
		else                              Run_Par.write_segs = true;
	    } else if ( (argi == "-h") || (argi == "--help") ) {
	        cout << USAGE_INFO << endl;
                fatal_error = true;
	        finish();
	    } else if(argi == "-local") {
                mOutput.local = true;
            } else if (isDatEnvArg(argv[i])) {
                //cout << "isDatEnvArg == true:  " << argv[i] << endl;
	        i++;
            } else {
	        cerr << "SenseMonitor WARNING:  Argument: " << argi;
	        cerr << " not recognized." << endl;
                cout << USAGE_INFO << endl;
                //fatal_error = true;
	        //finish();
	    }
        }
    }
}


//====================================== Run list of checks on parameters.
//:KLUDGE: Since this is a method of the SenseMonitor class, shouldn't need to 
//hand it data members of the same class (Run_Par, fatal_error).
void
SenseMonitor::VerifyParameters(Parameter& Run_Par, bool& fatal_error)
{
    if ( ((int)(Run_Par.T) % Run_Par.num_ints) != 0 ) {
	cout << "SenseMonitor ERROR: The number of averages 'n' must divide the length 'stride'.\n"; 
	cout << "                    Exiting.\n";
	cout << USAGE_INFO << endl;
	fatal_error = true;
	finish();
    } 
    if ( Run_Par.l_freq >= Run_Par.h_freq ) {
        cout << "SenseMonitor ERROR: Low frequency must be less than high ";
	cout << "frequency.  Exiting.\n";
        cout << USAGE_INFO << endl;
        fatal_error = true;
        finish();
    }
/*
    if (((int)(86400.0 / Run_Par.T) != (86400.0 / Run_Par.T)) && mOutput.trend) {
	cout << "SenseMonitor ERROR: 'stride' must divide 86400 (number of seconds in 1 day).\n";
	cout << "                    Exiting.\n";
	cout << USAGE_INFO << endl;
	fatal_error = true;
	finish();
    }
*/
    if (Run_Par.window_name == "square") {
	cout << "SenseMonitor WARNING: Square window (ie, no windowing) selected for\n"; 
	cout << "                      power-spectrum estimation.  This is a VERY BAD\n";
	cout << "                      IDEA; you should select anything else.\n";
    }
    if (!fatal_error) {
        //if ( !mOutput.local && !((Run_Par.cumlog_file_name).empty()) ) {
        if ( !mOutput.local && !((mOutput.cumlog_file_name).empty()) ) {
	    cout << "SenseMonitor WARNING: -local option has not been selected; specification\n";
            cout << "                      of log-file name will be ignored.\n";
        }
//----- :KLUDGE: Put in corresponding check that xml calibration file can be read.
/*
        //---------- Get name of calibration file.
	if (Run_Par.statcal) {
	    if ( (Run_Par.cal_file_name).empty() ) {
  	        cout << "SenseMonitor WARNING: No calibration file specified; will look for default\n";
                cout << "                      calibration file.\n";
                //---------- Look for calibration file in directory specified
                //           in $SENSEMONITOR_CALIBRATION instead of in local directory.
                const char* cal_dir = getenv("SENSEMONITOR_CALIBRATION");
                if (cal_dir) {
                    Run_Par.cal_file_name = string(cal_dir) + "/";
                } else { 
	            cout << "SenseMonitor WARNING: Environment variable $SENSEMONITOR_CALIBRATION not set;\n"; 
                    cout << "                      looking in local directory for default calibration file.\n";
                }
                Run_Par.cal_file_name += "SenseMonitor_Calibration_";
                Run_Par.cal_file_name += Run_Par.IFO;
                Run_Par.cal_file_name += ".txt";
            }
            //---------- Verify that calibration file can be opened.
            ifstream input;
            input.open((Run_Par.cal_file_name).c_str());
            if (input.fail()) {
                cerr << "SenseMonitor ERROR: Cannot open calibration file ";
                cerr << Run_Par.cal_file_name  << ".  Exiting.\n";
                fatal_error = true;
                finish();
            } else {
	        cout << "SenseMonitor MESSAGE: Found calibration file ";
	        cout << Run_Par.cal_file_name << endl;
	    }
            input.close();
	} else {
            //---------- Verify that open-loop-gain and sensing-function files can be opened.
            ifstream input;
            input.open((Run_Par.olg_file_name).c_str());
            if (input.fail()) {
                cerr << "SenseMonitor ERROR: Cannot open open-loop-gain file ";
                cerr << Run_Par.olg_file_name  << ".  Exiting.\n";
                fatal_error = true;
                finish();
            } else {
	        cout << "SenseMonitor MESSAGE: Found open-loop-gain file ";
	        cout << Run_Par.olg_file_name << endl;
	    }
            input.close();
            input.open((Run_Par.sf_file_name).c_str());
            if (input.fail()) {
                cerr << "SenseMonitor ERROR: Cannot open sensing-function file ";
                cerr << Run_Par.sf_file_name  << ".  Exiting.\n";
                fatal_error = true;
                finish();
            } else {
	        cout << "SenseMonitor MESSAGE: Found sensing-function file ";
	        cout << Run_Par.olg_file_name << endl;
	    }
            input.close();
	}
*/
        //---------- Get name of OSC configuration files.
        if ( (Run_Par.osc_file_name).empty() ) {
  	    cout << "SenseMonitor WARNING: No OSC configuration file specified; will look for\n";
            cout << "                      default OSC configuration file.\n";
            //---------- Look for OSC config file in directory specified
            //           in $SENSEMONITOR_OSCCONF instead of in local directory.
            const char* osc_dir = getenv("SENSEMONITOR_OSCCONF");
            if (osc_dir) {
                Run_Par.osc_file_name = string(osc_dir) + "/";
            } else { 
	        cout << "SenseMonitor WARNING: Environment variable $SENSEMONITOR_OSCCONF not set;\n"; 
                cout << "                      looking in local directory for default OSC file.\n";
            }
            //Run_Par.osc_file_name += "LockLoss_";
            //Run_Par.osc_file_name += Run_Par.IFO;
            Run_Par.osc_file_name += "SenseMonitor_LockLoss.conf";
        }
//----- :KLUDGE: Put in corresponding check that xml calibration file can be read.
/*
        //---------- Verify that calibration file can be opened.
        ifstream input;  //----- previous declaration inside if () statement therefore gone.
        input.open((Run_Par.osc_file_name).c_str());
        if (input.fail()) {
            cerr << "SenseMonitor ERROR: Cannot open OSC configuration file ";
            cerr << Run_Par.osc_file_name  << ".  Exiting.\n";
            fatal_error = true;
            finish();
        } else {
	    cout << "SenseMonitor MESSAGE: Found OSC configuration file ";
	    cout << Run_Par.osc_file_name << endl;
	}
        input.close();
*/
    }
}

//====================================== Read waveform file.
//std::string
FSpectrum
SenseMonitor::ReadWaveform(const std::string& filename, double low_freq, double df) const {

    //---------- Initializations.
    FSpectrum waveform;
    vector <float> freq;
    vector <float> val;

    //---------- Read columns.
    ifstream in(filename);
    while (!in.eof())
    {
        float tmp_freq, tmp_val;
        in >> tmp_freq >> tmp_val;
        freq.push_back(tmp_freq);
        val.push_back(tmp_val);
    }

    //---------- Set low-frequency cutoff and frequency step.
//    double low_freq = freq[0];
//    double df = freq[1] - freq[0];
    waveform.clear(low_freq, df);

    //---------- Set the waveform values.
    waveform.setData(val.size(), val.data());

    return waveform;
};

//====================================== Write current range, etc to specified
//                                       output stream.
void
SenseMonitor::ReportResults(ostream& out, const Time& time, Parameter& Run_Par)
{
	//---------- Write the range, alpha, etc. to the specified 
	//	     output stream.
	out << time.getS();
	out << ' ' << setw(15) << R_Dat.range;
	out << ' ' << setw(15) << R_Dat.ampl;
	out << ' ' << setw(15) << R_Dat.alpha;
	out << ' ' << setw(15) << (R_Dat.alpha)*(R_Dat.beta);
        out << endl;
}

