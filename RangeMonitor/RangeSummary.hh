//
//  Class Declaration: RangeSummary.hh 
//
//  $Id: RangeSummary.hh 7094 2014-06-10 08:16:49Z john.zweizig@LIGO.ORG $
//
//////////////////////////////////////////////////////////////////////


#ifndef RangeSummary_HH
#define RangeSummary_HH

#include <fstream>
#include <list>

/** Helper struct for class RangeSummary.
  * @memo Helper struct for class RangeSummary.
  */
struct RangeDataSummary
{
        unsigned long gpstime;
        double range;
};


/** Helper class for SenseMonitor.
  * RangeSummary produces two-week summary data files containing 
  * the range estimates and other monitor output averaged over 
  * 1000 second intervals.
  * @author 
  * @memo Produces two-week summary data files.
  */
class RangeSummary
{
public:
    /** RangeSummary constructor.
      * @memo Constructor.
      */
    RangeSummary(void) {};

    /** RangeSummary destructor.
      * @memo Destructor.
      */
    ~RangeSummary(void);

    /// What the hell does this do?
    void init(const std::string &htmldir,
              const Time &now);

    /// Append data point to end of list.
    void append(RangeDataSummary state);

    /// Dump data to file
    void dumpList(void);

private:
    //  length of list
    static const unsigned long mLength = 1210;

    //  Summary "sampling" interval
    static const unsigned long mDT = 1000;

    // list of range data for last 2 weeks
    typedef std::list<RangeDataSummary> Range_List;
    typedef Range_List::const_iterator list_iter;
    Range_List mRangeList;

    //  status file name
    std::string        mFilename;

    //  status file stream
    std::fstream       mFile;

    //  time this object was created
    Time               mCreationTime;
};


//=============================================================================
//Class Member Definitions
//=============================================================================


void RangeSummary::init(const std::string &summary_name,
                        const Time &now)
{
//  mFilename = htmldir + "/summary.txt";
    mFilename = summary_name;

    // set the current time back to the last '000,
    // ditto the start time
    unsigned long nowRound = now.getS() - (now.getS() % mDT);
    unsigned long startTime = nowRound - ( mLength - 1 ) * mDT;

    std::ifstream File(mFilename.c_str(), std::ios::in);
    if ( !File.is_open( ) ) {
        // the summary file does not exist; fill in the list with -1's
        // and create the file for output

        RangeDataSummary tmpstate;
        for (unsigned long i=0; i < mLength; ++i) {
            tmpstate.gpstime  = startTime + i*mDT;
            tmpstate.range = -1.;

            mRangeList.push_back(tmpstate);
        }

        dumpList();
    } else {
        // the summary file exists; need to read it and restore
        // history;
        // there are two possibilities:
        //   1. the start of the summary file is what the start time
        //      should be
        //   2. the start of the summary file has "dropped off the
        //      screen", i.e. it is further back in the past than
        //      we need.

        // read file
        RangeDataSummary tmpstate;
        while ( File.is_open( ) ) {
            File >> tmpstate.gpstime >> tmpstate.range;
            if (File.good() && !File.eof()) {
                mRangeList.push_back(tmpstate);
	    } else {
	        // done reading. close file
                File.close( );
	    }
        }

        // unsigned long fileStartTime = mRangeList.front().gpstime;

        // pop off data that is too old
        while (!mRangeList.empty() && mRangeList.front().gpstime < startTime) {
            mRangeList.pop_front();
        }

        // if the summary list is too short, pad the missing data,
        // and write out the file
        // WARNING:  This assumes any gaps in the data are at the 
	// most recent end, and not, say, in the middle.  Those gaps 
	// will not be filled.
	if (mRangeList.size() < mLength) {

            RangeDataSummary tmpstate;
            while (mRangeList.back().gpstime < nowRound) {
                tmpstate.gpstime = mRangeList.back().gpstime + mDT;
                tmpstate.range = -1.;
                mRangeList.push_back(tmpstate);
            }

            // write out list
            dumpList();
        }
    }
}

//
// RangeSummary::~RangeSummary()
//
RangeSummary::~RangeSummary(void)
{
  //    if (mFile.is_open()) {
  //      mFile.close();
  //  }
}

//
// RangeSummary::append()
//
void RangeSummary::append(RangeDataSummary state)
{
    mRangeList.pop_front();
    mRangeList.push_back(state);
}


//
// RangeSummary::dumpList()
//
void RangeSummary::dumpList(void)
{
    // make sure output file is open
    std::ofstream File(mFilename.c_str(), std::ios::out);
    
    // set output format
    File.setf(std::ios::scientific);
    File.precision(4);
    if (!File.good() || !File.is_open()) {
        std::cerr << "RangeSummary::dumpList(): summary status file is "
                  << " not open: opening" << std::endl;
	return;
    }

    for (list_iter iter = mRangeList.begin(); iter != mRangeList.end(); ++iter)
        File << (*iter).gpstime << "\t" << (*iter).range << std::endl;

    File.close();
}

#endif     //  RangeSummary_HH
