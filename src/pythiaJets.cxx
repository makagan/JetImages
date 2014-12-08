//std things 
#include <map>
#include <string>
#include <sstream>
#include <iostream>
#include <time.h>
//#include <unistd.h>

//things josh installed
#include "boost/program_options.hpp"

//physics things josh installed
#include "Pythia.h"
#include "JetMatching.h"
#include "CombineMatchingInput.h"
#include "HepMCInterface.h"
//#include "CellBooster.h"

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


void setLeptonDecays(const char MATRIX_ELEMENT, const int LEP_CHARGE, Pythia8::Pythia &pythia){

  if( !(MATRIX_ELEMENT=='w' || MATRIX_ELEMENT=='t' || MATRIX_ELEMENT=='v') )
    return;


  if(LEP_CHARGE==0)
    {
      pythia.readString("24:onIfAny = 13"); //dilepton
      
      if(MATRIX_ELEMENT=='w') cout << "Turning on dilepton W^+ W^- --> l^+ + nu + l^- + nu  production" << endl;
      if(MATRIX_ELEMENT=='t') cout << "Turning on dilepton t tbar  --> l^+ + nu + l^- + nu b + b production" << endl;
      if(MATRIX_ELEMENT=='v') cout << "Turning on W(->l nu)+jet production" << endl;
    }
  else if(LEP_CHARGE==1)
    {
      pythia.readString("24:onPosIfAny = 13"); //W^+ decays leptonically to electron or muon
      if(MATRIX_ELEMENT=='t' || MATRIX_ELEMENT=='w') 
	pythia.readString("24:onNegIfAny = 1 2 3 4 5");// W^- decays hadronically
      
      if(MATRIX_ELEMENT=='w') cout << "Turning on semi-leptonic W^+ W^- -> l^+ + nu + j+j production" << endl;
      if(MATRIX_ELEMENT=='t') cout << "Turning on semi-leptonic t tbar -> l^+ + nu + b + b + j+j  production" << endl;
      if(MATRIX_ELEMENT=='v') cout << "Turning on W^-(->l^+ nu)+jet production" << endl;  
    }
  else if(LEP_CHARGE==-1)
    {
      pythia.readString("24:onNegIfAny = 13"); //W^- decays leptonically to electron or muon
      if(MATRIX_ELEMENT=='t' || MATRIX_ELEMENT=='w') 
	pythia.readString("24:onPosIfAny = 1 2 3 4 5");// W^+  decays hadronically

      
      if(MATRIX_ELEMENT=='w') cout << "Turning on semi-leptonic W^+ W^- -> l^- + nu + j+j production" << endl;
      if(MATRIX_ELEMENT=='t') cout << "Turning on semi-leptonic t tbar -> l^- + nu + b + b + j+j  production" << endl;
      if(MATRIX_ELEMENT=='v') cout << "Turning on W^-(->l^- nu)+jet production" << endl;     
    }
  else
    {
      pythia.readString("24:onIfAny = 1 2 3 4 5"); //all hadronic
      
      if(MATRIX_ELEMENT=='w') cout << "Turning on W^+ W^- -> all hadronic production" << endl;
      if(MATRIX_ELEMENT=='t') cout << "Turning on t tbar -> all hadronic  production" << endl;
      if(MATRIX_ELEMENT=='v') cout << "Turning on W(->hadronic)+jet production" << endl;            
    }

  return;
}

int main(int argc, char* argv[]) {
    cout << "Called as: ";
    for(int ii=0; ii<argc; ++ii){
        cout << argv[ii] << " ";
    }
    cout << endl;


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
      ("LHEFile", po::value<string>(&LHEFile)->default_value("")
       , "Location of an LHE file, over-rides all other Pythia event generation options")
      ("matching",  "Turn on jet matching (e.g. for Madgraph LHE input)")
      ("metype",   po::value<char>(&optME)->default_value('g')
       , "g,q,l,b,c,v,t,W")
      ("lepcharge", po::value<int>  (&optLepQ)       ->default_value(1)
       , "Charge of lepton desired in final state, only for ME ={v, t, W}.  1=l^+, -1=l^-, 0= All lepton flavors, any other value = no lepton (all hadronic),")
      ("seed", po::value<int>(&optSeed)               ->default_value(-1)
       , "Random seed.  The default, -1, causes the application to "
       "choose one from the ns time as pid")
      ("savemaxnjets", po::value<int>(&optSaveMaxNJets) ->default_value(-1)
       , "Max number of jets to save.  Default, -1, saves all jets (which satisfy other selectons)")
      ("nevents", po::value<int>  (&optNEvents)       ->default_value(1)
       , "# Events to generate")
      ("cellv", po::value<int>(&optCellV)               ->default_value(1)
       , "Control how many and which cells to output.  0 shows no cells, 1 (default) shows only those that survive trimming and 2 shows all cells in a jet.")
      ("meminpt"   , po::value<float>(&optMEMinPt)    ->default_value(50.0)
       , "Pythia PTHatMin (GeV)")
      ("memaxpt"   , po::value<float>(&optMEMaxPt)    ->default_value(20000)
       , "Pythia PTHatMax (GeV)")
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
      ("disablempi", "Dislable Multiple Parton Interactions")
      ("n_B"     , po::value<int>(&optNB)            ->default_value(-1)
       , "Exactl number of b-hadrons per-jet (-1 for no cut)")
      ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
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

//    {
//    CellBooster cb(R_SEED);
//    cb.split_cell(100, 10, .1, 0.2, 0, 0.1);
//    }
//    return 1;

    std::stringstream ss;
    Pythia8::Pythia pythia;

    pythia.readString("Random:setSeed = on");
    ss.clear(); ss.str("");
    ss << "Random:seed = " << R_SEED;
    cout << ss.str() << endl;
    pythia.readString(ss.str());  

    if (LHEFile == "") {
      optME='?';

      pythia.readString("HardQCD:all = off");

      if (MATRIX_ELEMENT == 'g'){
        pythia.readString("HardQCD:gg2gg = on");    
        pythia.readString("HardQCD:qqbar2gg = on");    
        cout << "Turning on gg production" << endl;
      } else if (MATRIX_ELEMENT == 'q'){
        pythia.readString("HardQCD:gg2qqbar = on");    
        pythia.readString("HardQCD:qq2qq = on");    
        cout << "Turning on qq production" << endl;
      } else if (MATRIX_ELEMENT == 'l') {
	pythia.readString("HardQCD:gg2gg = on");
	pythia.readString("HardQCD:gg2qqbar = on");
	pythia.readString("HardQCD:qg2qg = on");
	pythia.readString("HardQCD:qq2qq = on");
	pythia.readString("HardQCD:qqbar2gg = on");
	pythia.readString("HardQCD:qqbar2gg = on");
	pythia.readString("HardQCD:qqbar2qqbarNew = on");
      } else if (MATRIX_ELEMENT == 'c') {
	pythia.readString("HardQCD:gg2ccbar = on"); 
	pythia.readString("HardQCD:qqbar2ccbar = on"); 
      } else if (MATRIX_ELEMENT == 'b') {
	pythia.readString("HardQCD:gg2bbbar = on");   
	pythia.readString("HardQCD:qqbar2bbbar = on"); 
      } else if (MATRIX_ELEMENT == 'v'){
        pythia.readString("WeakBosonAndParton:qqbar2Wg = on");    
        pythia.readString("WeakBosonAndParton:qg2Wq = on"); 
        pythia.readString("24:onMode = off");
        cout << "When ME = " << MATRIX_ELEMENT << " lepton charge forced to 0" << endl;
        setLeptonDecays(MATRIX_ELEMENT, 0, pythia);

      } else if (MATRIX_ELEMENT == 't'){
        pythia.readString("Top:gg2ttbar = on");    
        pythia.readString("Top:qqbar2ttbar = on");  
        //pythia.readString("Top:ffbar2ttbar = on");  
        //pythia.readString("Top:gmgm2ttbar = on");
        pythia.readString("24:onMode = off");
        setLeptonDecays(MATRIX_ELEMENT, LEP_CHARGE, pythia);

      } else if (MATRIX_ELEMENT == 'w'){ 
        pythia.readString("WeakDoubleBoson:ffbar2WW = on");
        pythia.readString("24:onMode = off");
        setLeptonDecays(MATRIX_ELEMENT, LEP_CHARGE, pythia);

      } else {
        cout << "Unknown Matrix Element type " << MATRIX_ELEMENT << endl;
        cout << "Try q or g or v or t or W" << endl;
        return -1;
      }

      //Only applicable to 2->2 processes, ie use WW
      ss.clear(); ss.str("");
      ss << "PhaseSpace:pTHatMin = " << ME_PT_MIN;
      cout << ss.str() << endl;
      pythia.readString(ss.str());  

      ss.clear(); ss.str("");
      ss << "PhaseSpace:pTHatMax = " << ME_PT_MAX;
      cout << ss.str() << endl;
      pythia.readString(ss.str());  

      if (vm.count("disablempi")>0) pythia.readString("PartonLevel:MI = off");
      pythia.init(2212, 2212, 8000);

    } else {
      if (vm.count("matching") > 0) {
        cout << "Set matching parameters" << endl;
        pythia.readString("JetMatching:merge = on");
        pythia.readString("JetMatching:setMad = off:");
        pythia.readString("JetMatching:scheme = 1");
        pythia.readString("JetMatching:qCut = 30");
        pythia.readString("JetMatching:nQmatch = 5");
        pythia.readString("JetMatching:clFact = 1");

        CombineMatchingInput combined;
        UserHooks* matching = combined.getHook(pythia);
        if (!matching) return 1;
        pythia.setUserHooksPtr(matching);
      }
      pythia.init(LHEFile.c_str());
    }

    pythia.settings.listAll();
    

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
    props["source"] = "Pythia8";
    ss.clear(); ss.str("");
    ss << MATRIX_ELEMENT;
    props["ME"]     = ss.str();
    ss.clear(); ss.str("");
    ss << R_SEED;
    props["seed"]   = ss.str();

    HepMC::GenEvent *gevt = new HepMC::GenEvent();


  for (int iEvent = 0; iEvent < N_EVENTS; ++iEvent) {
    if (!pythia.next()) continue;
    //pythia.event.list();

    ss.clear(); ss.str("");
    ss << iEvent;
    props["eventNumber"] = ss.str();

    delete gevt;
    gevt = new HepMC::GenEvent();

    HepMC::I_Pythia8().fill_next_event(pythia, gevt); 
    jet_analysis.analyze(gevt, props);
    continue;

  }
  pythia.stat();
  return 0;
}
