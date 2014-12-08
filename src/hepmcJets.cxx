//std things 
#include <map>
#include <string>
#include <sstream>
#include <iostream>
#include <time.h>
//#include <unistd.h>

//things josh installed
#include "boost/program_options.hpp"

#include "HepMC/IO_HEPEVT.h"
#include "HepMC/IO_GenEvent.h"

//physics things josh wrote 
#include "JetAnalysis.h"

#define HAS_NS_TIME
#define _USE_MATH_DEFINES

using std::cout;
using std::endl;
using std::string;
using std::map;
namespace po = boost::program_options;

int getSeed(int optSeed){
    if (optSeed > -1) return optSeed;

    int timeSeed = time(NULL);

#ifdef HAS_NS_TIME
    {
        timespec tp;
        if ( 0 == clock_gettime(CLOCK_REALTIME, &tp)) timeSeed = tp.tv_nsec;
    }
#endif

    //http://arxiv.org/abs/1005.4117 supposed to keep collisions down to
    //My seed is strong...
    return abs(((timeSeed*181)*((getpid()-83)*359))%104729);
}


int main(int argc, char* argv[]) {
    cout << "Called as: ";
    for(int ii=0; ii<argc; ++ii){
        cout << argv[ii] << " ";
    }
    cout << endl;

    string optSource = "HepMC";
    char optME = 'z';
    int  optLepQ = 1;
    int  optSeed = -1;
    int  optCellV= 1;
    int  optSaveMaxNJets = -1;
    int  optNEvents = 100;
    float optMEMinPt = 50.0;
    float optMEMaxPt = 20000;
    int   optTrimJet = 1;
    float optJetMinPt = 50.0;
    float optJetMaxPt = 99999.0;
    float optJetMinM  = 0;
    float optJetMaxM  = 20*1000;
    float optJetRadius = 0.4;
    float optMu = 0.0;
    int optNB = -1;
    string LHEFile = "";
    po::options_description desc("Allowed options");
    desc.add_options()
      ("help", "produce help message")
      ("source", po::value<string>(&optSource)->default_value("HepMC"),
       "Generator for HepMC input")
      ("metype",   po::value<char>(&optME)->default_value('g')
       , "g,q,v,t,W")
      ("lepcharge", po::value<int>  (&optLepQ)       ->default_value(1)
       , "Charge of lepton desired in final state, only for ME ={v, t, W}.  1=l^+, -1=l^-, 0= All lepton flavors, any other value = no lepton (all hadronic),")
      ("seed", po::value<int>(&optSeed)               ->default_value(-1)
       , "Random seed.  The default, -1, means the HepMC source seed is unknown")
      ("savemaxnjets", po::value<int>(&optSaveMaxNJets) ->default_value(-1)
       , "Max number of jets to save.  Default, -1, saves all jets (which satisfy other selectons)")
      ("nevents", po::value<int>  (&optNEvents)       ->default_value(1)
       , "# Events to generate")
      ("cellv", po::value<int>(&optCellV)               ->default_value(1)
       , "Control how many and which cells to output.  0 shows no cells, 1 (default) shows only those that survive trimming and 2 shows all cells in a jet.")
      ("meminpt"   , po::value<float>(&optMEMinPt)    ->default_value(-1.0)
       , "PTHatMin used in HepMC File(GeV)")
      ("memaxpt"   , po::value<float>(&optMEMaxPt)    ->default_value(-1.0)
       , "PTHatMax used in HepMC File(GeV)")
      ("trimjet", po::value<int>  (&optTrimJet)       ->default_value(1)
       , "To run trimming or not.  Any value other than 1 (default) disables trimming.")
      ("jetminpt"  , po::value<float>(&optJetMinPt)   ->default_value(25.0)
       , "Jet PTMin (GeV)")
      ("jetmaxpt"  , po::value<float>(&optJetMaxPt)   ->default_value(99999.0)
       , "Jet PTMax (GeV)")
      ("jetminm"  , po::value<float>(&optJetMinM)     ->default_value(0)
       , "Jet Min Mass (GeV)")
      ("jetmaxm"  , po::value<float>(&optJetMaxM)     ->default_value(20000)
       , "Jet Max Mass (GeV)")
      ("rad"    , po::value<float>(&optJetRadius)     ->default_value(0.4)
       , "Jet Finder Radius")
      ("mu"     , po::value<float>(&optMu)            ->default_value(0.0)
       , "Average Int./Cross")
      ("smear"  , "Turn on radial energy smearing")
      ("wtagger", "Turn on MVA Tagger output")
      ("input-file", po::value< vector<string> >(), "input file")
      ("n_B"     , po::value<int>(&optNB)            ->default_value(-1)
       , "Exactl number of b-hadrons per-jet (-1 for no cut)")
      ;

    po::positional_options_description p;
    p.add("input-file", -1);
    
    po::variables_map vm;
    //po::store(po::parse_command_line(argc, argv, desc), vm);
    po::store(po::command_line_parser(argc, argv).
              options(desc).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help")>0){
        cout << desc << "\n";
        return 1;
    }

    const char MATRIX_ELEMENT = tolower(optME);
    const int LEP_CHARGE = optLepQ;
    const int N_EVENTS = optNEvents;
    const int R_SEED = getSeed(optSeed);
    const int ME_PT_MIN = optMEMinPt;
    const int ME_PT_MAX = optMEMaxPt;
    const int DO_TRIMMING = optTrimJet;
    const float MU = optMu;

    std::stringstream ss;

    JetAnalysis jet_analysis(optSaveMaxNJets, 
                             DO_TRIMMING, 
                             optJetMinPt, 
                              optJetMaxPt, 
                             optJetMinM, 
                             optJetMaxM, 
                             optJetRadius, 
                             optNB,
                             MU, 
                             R_SEED,
                             (vm.count("smear")>0 ? true : false),
                             (vm.count("wtagger")>0 ? true : false),
                             optCellV
                             );

    map<string, string> props;
    props["source"] = optSource;

    ss.clear(); ss.str("");
    ss << MATRIX_ELEMENT;
    props["ME"]   = ss.str();

    ss.clear(); ss.str("");
    ss << R_SEED;
    props["seed"]   = ss.str();

    HepMC::GenEvent *gevt = new HepMC::GenEvent();

    if (vm.count("input-file")) {
      cout << "Input files are: " << endl;
      vector<string>::const_iterator iFile = vm["input-file"].as< vector<string> >().begin();
      vector<string>::const_iterator iFileEnd = vm["input-file"].as< vector<string> >().end();
      for (; iFile != iFileEnd; ++iFile) {
        cout << *iFile << "\n";

        HepMC::IO_GenEvent ascii_in(iFile->c_str(),std::ios::in);

        gevt = ascii_in.read_next_event();
        int iEvent = 0;
        while ( (N_EVENTS <= 0 || iEvent < N_EVENTS) && gevt) {
          iEvent++;
          if ( iEvent%50==1 ) std::cout << "Processing Event Number " << iEvent
                                        << " its # " << gevt->event_number() 
                                        << std::endl;
          

          ss.clear(); ss.str("");
          ss << iEvent;
          props["eventNumber"] = ss.str();

          jet_analysis.analyze(gevt, props);

          delete gevt;
          gevt = ascii_in.read_next_event();

        } // while events

        delete gevt;
      } // for input files
      
    } // if input-file


//   for (int iEvent = 0; iEvent < N_EVENTS; ++iEvent) {

//     ss.clear(); ss.str("");
//     ss << iEvent;
//     props["eventNumber"] = ss.str();

//     delete gevt;
//     gevt = new HepMC::GenEvent();

//     jet_analysis.analyze(gevt, props);
//     continue;

//   }
  return 0;
}
