#include <sstream>
#include "Pythia.h"

#include "PileUp.h"
#include "PythiaPileUp.h"
#include "Calorimeter.h"


void
PythiaPileUpGenerator::init(){

    cout << "INIT from PPU"<< endl;
    pythia.readString("HardQCD:all = off");
    pythia.readString("SoftQCD:minBias = on");
    pythia.readString("Random:setSeed = on");
    pythia.readString("PhaseSpace:pTHatMin  = .1");

    std::stringstream ss;
    ss << "Random:seed = " << m_rSeed;
    pythia.readString(ss.str());

    pythia.init(2212, 2212, 8000);
}

PythiaPileUpGenerator::PythiaPileUpGenerator(
                    double averageIntPerCrossing, 
                    Calorimeter& refCalo,
                    int randomSeed)
: PileUpGenerator(averageIntPerCrossing, refCalo, randomSeed)
{
    init();//child class init needs explicit calling
}


Calorimeter 
PythiaPileUpGenerator::next(int exactIntPerCrossing=-1){
    const int usrInts = exactIntPerCrossing;
    const unsigned int nInts = usrInts > -1 ? usrInts : poisson_sample(m_mu);
    Pythia8::Particle pp;

    Calorimeter calo_fill = Calorimeter(m_refCalo);
    calo_fill.set_nInteractions(nInts);
    vector<PseudoJet> pjs;
    PseudoJet pj;
    for(unsigned int i_pu = 0; i_pu < nInts; ){
        if (!pythia.next()) continue;

        for (int iPart = 0; iPart < pythia.event.size(); ++iPart) {
            pp = pythia.event[iPart]; 
            if (pp.isFinal() && pp.isVisible()) {
              pj.reset_momentum(pp.px(), pp.py(), pp.pz(), pp.e());
              pj.set_user_index(pp.id());
              pjs.push_back(pj);
            }
        }
        //don't count events pythia crapped out on
        ++i_pu;
    }

    calo_fill.push(pjs);

    return calo_fill;
}

vector<Calorimeter> 
PythiaPileUpGenerator::next_vector(int exactIntPerCrossing=-1){
    const int usrInts = exactIntPerCrossing;
    const unsigned int nInts = usrInts > -1 ? usrInts : poisson_sample(m_mu);
    Pythia8::Particle pp;

    vector<Calorimeter> ret_vec;
    vector<PseudoJet> pjs;
    PseudoJet pj;

    for(unsigned int i_pu = 0; i_pu < nInts; ){
        if (!pythia.next()) continue;

        pjs.clear();

        for (int iPart = 0; iPart < pythia.event.size(); ++iPart) {
            pp = pythia.event[iPart]; 
            if (pp.isFinal() && pp.isVisible()) {
              pj.reset_momentum(pp.px(), pp.py(), pp.pz(), pp.e());
              pj.set_user_index(pp.id());
              pjs.push_back(pj);
            }
        }

        ret_vec.push_back(Calorimeter(m_refCalo));
        ret_vec[i_pu].set_nInteractions(1);
        ret_vec[i_pu].push(pjs);

        //don't count events pythia crapped out on
        ++i_pu;
    }


    return ret_vec;
}
