#include <assert.h>
#include <vector>
#include <algorithm>

//physics things josh found 
#include "HepMC/GenEvent.h"

//physics things josh installed
#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/Filter.hh"
#include "Nsubjettiness/Nsubjettiness.hh"

//physics things josh wrote
#include "JetAnalysis.h"
#include "Calorimeter.h"
#include "PythiaPileUp.h"
#include "ExtendableJet.h"

//Harvard MVA Tagger
#include "wtag.h"

void Cuts::passed(std::string cut_name) {
  if (cut_name != "" && m_counter >= m_name.size()) {
    m_name.push_back(cut_name);
    m_passed.push_back(0);
  }
  if (cut_name == "" || m_name[m_counter] == cut_name) {
    m_passed[m_counter++] += 1;
  }
}

void Cuts::print() {
  cout << "Cut Summary: " << endl;
  for (int i = 0, iEnd = m_passed.size(); i < iEnd; ++i) {
    cout << "  " << m_name[i] << ": " << m_passed[i] << endl;
  }
}

void 
JetAnalysis::init(int saveMaxNJets, int trimJets, float jetMinPt, float jetMaxPt,
                  float jetMinM, float jetMaxM, 
                  float jetRadius, int nB, int mu, int puSeed, bool doSmear, bool doWTagger, int cellV){

    m_saveMaxNJets = saveMaxNJets; 
    m_trimJets = trimJets; 
    m_jetMinPt = jetMinPt; 
    m_jetMaxPt = jetMaxPt; 
    m_jetMinM = jetMinM; 
    m_jetMaxM = jetMaxM; 
    m_jetRadius = jetRadius; 
    m_nB = nB;
    m_mu = mu;
    m_doSmear = doSmear;
    m_doWTagger = doWTagger;
    m_cellV = cellV;

    m_minDR  = std::min(2.5*m_jetRadius, 1.5);
    m_trimRadius  = 0.3;

    //recall this is needed for BDRS MassDropTagger
    m_primaryDef = fastjet::JetDefinition(fastjet::cambridge_algorithm, m_jetRadius);
    m_trimDef  = fastjet::JetDefinition(fastjet::kt_algorithm, m_trimRadius);
    m_selPtMin = fastjet::SelectorPtMin(m_jetMinPt);
    m_selPtMax = fastjet::SelectorPtMax(m_jetMaxPt);
    m_basePtMin= fastjet::SelectorPtMin(25);
    m_selDR    = fastjet::operator!(fastjet::SelectorDoughnut(.00001, m_minDR));
    m_selEta   = fastjet::SelectorAbsRapRange(0, 2.5);

    std::stringstream ss;
    ss << m_jetRadius << "_" << m_jetMinPt << "_" << m_jetMaxPt << "_" << m_minDR ;
    m_desc = ss.str();

    m_hs_calo = Calorimeter(50, 63, false, 2.5);
    m_pug = new PythiaPileUpGenerator(m_mu, m_hs_calo, puSeed);
    m_hs_calo.setSmear(m_doSmear);

    m_wtagEffs.clear();
    m_wtagEffs.push_back(10);
    m_wtagEffs.push_back(15);
    m_wtagEffs.push_back(20);
    m_wtagEffs.push_back(25);
    m_wtagEffs.push_back(30);
    m_wtagEffs.push_back(35);
    m_wtagEffs.push_back(40);
    m_wtagEffs.push_back(45);
    m_wtagEffs.push_back(50);
    m_wtagEffs.push_back(55);
    m_wtagEffs.push_back(60);
    m_wtagEffs.push_back(65);
    m_wtagEffs.push_back(70);
    m_wtagEffs.push_back(75);
    m_wtagEffs.push_back(80);
    m_wtagEffs.push_back(85);
    m_wtagEffs.push_back(90);
    m_wtagEffs.push_back(95);

    for (int ii=1; ii<=5; ++ii){
        m_nSubCalcs[ii] = new fastjet::contrib::Nsubjettiness::Nsubjettiness(
                                ii, 
                                fastjet::contrib::Njettiness::onepass_kt_axes, 
                                1, 
                                m_jetRadius, 
                                10000.0
                                );
    }
}

JetAnalysis::~JetAnalysis() {
    delete m_pug; 
    m_pug = 0;

    std::map<int, fastjet::contrib::Nsubjettiness *>::iterator it;
    for (it=m_nSubCalcs.begin(); it!=m_nSubCalcs.end(); ++it) delete it->second;

    m_jetCuts.print();
}

bool
JetAnalysis::is_qhadron(int qid, int pid){
    //working from http://pdg.lbl.gov/2002/montecarlorpp.pdf
    //if meson either 10 or 100's digit need to be 5
    //but no top quark bound state, so b always heaviest 
    //--> hundrends digit
    //if baryon then 5 in thousands digit!
    qid = abs(qid);
    pid = abs(pid);
    assert(qid > 0 && qid < 7);
    return (((pid/100)%10)==qid || ((pid/1000)%10)==qid);
}


bool
JetAnalysis::is_final_qhadron(int qid, HepMC::GenParticle *gp){
    //is b_hadron with no b's in children

    if (! is_qhadron(qid, gp->pdg_id())) return false;
    //cout << "#BPasses bhadron with id: " << abs(gp->pdg_id()) << endl;

    HepMC::GenVertex *end = gp->end_vertex();
    HepMC::GenVertex::particles_out_const_iterator it;
    for(it =end->particles_out_const_begin(); 
        it!=end->particles_out_const_end  (); ++it){
        
        if (is_qhadron(qid, (*it)->pdg_id())) return false;

    }
    return true;
}

//I'm gonna smack the snot out of a hepmc developer.
//Who creates a particle class and disallows the ability to store charge
//Do you *really* think that one last int...SHORT event is gonna break
//the bank?
//But WHATEVER i just redid the math myself thank you.
// only works for longer lived particles, not checked on
// susy/technicolor/etc
bool
JetAnalysis::is_charged(int pdg_id){
    int id = abs(pdg_id);
    if (
        id == 12 ||
        id == 14 ||
        id == 16 ||
        id == 22 
        ) return false;

    int q1 = id/10   % 10;
    int q2 = id/100  % 10;
    int q3 = id/1000 % 10;

    int charge  = (q1 % 2 == 0? 2: -1);
    if (q3 == 0) {
        charge -= (q2 % 2 == 0? 2: -1);
    } else {
        charge += (q2 % 2 == 0? 2: -1);
        charge += (q3 % 2 == 0? 2: -1);
    }

    return charge != 0;
}


void
JetAnalysis::analyze(HepMC::GenEvent* event, map<string, string> props){

    std::vector<PseudoJet> particlesForJets;
    std::vector<PseudoJet> finalJets  ;

    m_hs_calo.reset();
    m_hs_calo.set_nInteractions(1);

    Calorimeter pu_calo (m_hs_calo);
    Calorimeter all_calo(m_hs_calo);
    all_calo.reset();
    pu_calo .reset();
    ExtendableJet ext;


    vector<PseudoJet> pjs;
    vector<PseudoJet> bhadrons;
    vector<PseudoJet> chadrons;
    vector<PseudoJet> tracks;

    HepMC::GenEvent::particle_const_iterator it;
    HepMC::FourVector mom;
    //cout << "#B~~~~~~~~~~~~~~~~~~~~" << endl;
    for(it=event->particles_begin(); it != event->particles_end(); it++){ 
        int id = abs((*it)->pdg_id());
        mom = (*it)->momentum();
        PseudoJet pj(mom.px(), mom.py(), mom.pz(), mom.e());

        //if its heavy flavor its ID will be -1*pdgid
        pj.set_user_index(-1*id);

        if (pj.perp() > 10 && is_final_qhadron(5, (*it))) {
            bhadrons.push_back(pj);
        }
        if (pj.perp() > 10 && is_final_qhadron(4, (*it))) {
            chadrons.push_back(pj);
        }

        //final state?
        if (!( (*it)->end_vertex() == NULL && (*it)->status()== 1)) continue;

        if (pj.perp() > 0.4 && is_charged(id)){
            pj.set_user_index(-1);
            tracks.push_back(pj);
        }
        
        //neutrino or muon?
        if (id == 12 || id == 13 || id == 14 || id == 16) continue;

        //Drop electrons from a W decay
        //probably dont' want to do this for a b decay
//        if (id == 11
//            && (*it)->production_vertex() 
//            && (*it)->production_vertex()->particles_in_size() == 1
//            && abs((*(*it)->production_vertex()->particles_in_const_begin())->pdg_id()) == 24) continue;
//
        //if its not a heavy flavor its ID will become the cellindex
        pj.set_user_index(id);
        pjs.push_back(pj);
    }

    m_hs_calo.push(pjs);
    //cout << "#BThere are " << chadrons.size() << " c hadrons" << endl;

    //pu_calo = m_pug->next();
    vector<Calorimeter> pu_calos = m_pug->next_vector();
    for(unsigned int ii=0; ii<pu_calos.size(); ++ii) pu_calo += pu_calos[ii];
    all_calo = m_hs_calo + pu_calo;
    particlesForJets = all_calo.getNonZeroCells();

    for (unsigned int ii=0; ii<bhadrons.size(); ++ii){
        PseudoJet temp(bhadrons[ii]);
        temp.reset_momentum_PtYPhiM(1e-10,temp.rap(), temp.phi(), 1e-10);
        particlesForJets.push_back(temp);
    }
    for (unsigned int ii=0; ii<chadrons.size(); ++ii){
        PseudoJet temp(chadrons[ii]);
        temp.reset_momentum_PtYPhiM(1e-10,temp.rap(), temp.phi(), 1e-10);
        particlesForJets.push_back(temp);
    }
    for (unsigned int ii=0; ii<tracks.size(); ++ii){
        PseudoJet temp(tracks[ii]);
        temp.reset_momentum_PtYPhiM(1e-10,temp.rap(), temp.phi(), 1e-10);
        particlesForJets.push_back(temp);
    }

    //declaring these variables here (in addition to init) is VERY
    //important to allow jets to see their internal constitutents
    //I have no idea why.
    fastjet::ClusterSequence csForJets(particlesForJets, m_primaryDef); 
    finalJets = selectGoodJets(csForJets);
    fastjet::Filter filter(m_trimDef, fastjet::SelectorPtFractionMin(0.05));

    int nsaved=0;

    for(unsigned int iJet=0; iJet<finalJets.size(); ++iJet){
        m_jetCuts.begin();
        m_jetCuts.passed("Initial");
        if(m_saveMaxNJets>=0 && nsaved >= m_saveMaxNJets ){
            //cout << "too many jets in this event, skipping" << endl;
            break;
        }
        m_jetCuts.passed("NJets");
        PseudoJet tj;
        double originalMass = finalJets[iJet].m();

        if (m_trimJets == 1) tj = filter(finalJets[iJet]);
        else tj = finalJets[iJet];

        if (! m_selPtMin.pass(tj)) {continue;}
        m_jetCuts.passed("PtMin");

        if (! m_selPtMax.pass(tj)) {continue;}
        m_jetCuts.passed("PtMax");

        if (! inMassRange(tj.m())) {continue;}
        m_jetCuts.passed("MassRange");

        HepMC::GenParticle *bestParton;
        bestParton = findBestParton(tj, *event);
        if (! bestParton) {
            cout << "error: failed to find a bestparton" << endl;
            continue;
        }
        m_jetCuts.passed("BestParton");

        unsigned int n_ga_bs = 0;
        unsigned int n_ga_cs = 0;
        unsigned int n_ga_tracks = 0;
        vector<PseudoJet> cons = tj.constituents();
        for (unsigned int ii=0; ii<cons.size(); ++ii){
            int id = cons[ii].user_index();
            if (id >= 0) continue;
            if(id == -1)               ++n_ga_tracks;
            if(is_qhadron(4, -1 * id)) ++n_ga_cs;
            if(is_qhadron(5, -1 * id)) ++n_ga_bs;
        }

        nsaved+=1;

        ext = ExtendableJet(tj, all_calo, pu_calos, m_jetRadius, 2, m_trimRadius);
        ext.setProperty(props);
        ext.setProperty("pdgIDHardParton", bestParton->pdg_id());
        ext.setProperty("preTrimMass", originalMass);
        ext.setProperty("n_b", n_ga_bs);
        ext.setProperty("n_c", n_ga_cs);
        ext.setProperty("n_tracks", n_ga_tracks);
        //ext.setProperty("bdrs_mass", get_bdrs_mass(finalJets[iJet], false));
        ext.setProperty("bdrs_mass", get_bdrs_mass(finalJets[iJet], true));

        std::map<int, fastjet::contrib::Nsubjettiness *>::iterator it;
        std::stringstream ss;
        for (it=m_nSubCalcs.begin(); it!=m_nSubCalcs.end(); ++it){
            ss.str("");
            ss << "tau_" << it->first;
            ext.setProperty(ss.str(), it->second->result(tj));
        }
        //int id;
        //double dr;
        //if      (id ==  24) dr = eventInfo["wp_dr_prod"];
        //else if (id == -24) dr = eventInfo["wm_dr_prod"];
        //else                dr = -999;
        //ext.setProperty("w_jet_dr",    dr);

        if (m_doWTagger) {
          for (vector<int>::const_iterator iEff = m_wtagEffs.begin(), iEffEnd = m_wtagEffs.end(); iEff != iEffEnd; ++iEff) {
            bool wtagged = wtag(csForJets, finalJets[iJet], *iEff / 100.);
            ss.str("");
            ss << "wtag_" << *iEff;
            ext.setProperty(ss.str(), wtagged);
          }
          ext.setProperty("wtag_opt", wtag(csForJets, finalJets[iJet], -1));
        }

        cout << "#G" << ext.toString(m_cellV) << endl;
        m_jetCuts.passed("Final");
    }
  return;
}

std::vector<PseudoJet> 
JetAnalysis::selectGoodJets(fastjet::ClusterSequence &cs){
    std::vector<PseudoJet> temp;

    std::vector<PseudoJet> everShrinking = fastjet::sorted_by_pt(cs.inclusive_jets());
    everShrinking = m_basePtMin(everShrinking);
    temp = everShrinking;

//    fastjet::Selector dr = m_selDR;
//    for(unsigned int iJet=0; iJet<everShrinking.size(); ++iJet){
//        dr.set_reference(everShrinking[iJet]);
//        temp = dr(temp);
//    }
    everShrinking = temp;
    everShrinking = m_selEta(everShrinking);
    return everShrinking;
    //nevermore.
}

HepMC::GenParticle*
JetAnalysis::findBestParton(PseudoJet &jet, HepMC::GenEvent &event) const {
    PseudoJet pj;
    HepMC::FourVector mom;
    double maxMetric = 0.0;
    int indexOfMax = -1;
    
    HepMC::GenParticle *pp, *best;
    HepMC::GenEvent::particle_const_iterator it;
    for(it=event.particles_begin(); it != event.particles_end(); it++){
        pp = *it;
        if (pp->status() > 80) continue;

        mom = pp->momentum();
        pj = PseudoJet(mom.px(), mom.py(), mom.pz(), mom.e());
        if (pj.delta_R(jet) > m_jetRadius) continue;

        if (best == NULL or pj.pt() > maxMetric){
            best = pp;
            maxMetric  = pj.pt();
        }
    }

    return best;
}

double
JetAnalysis::get_bdrs_mass(const PseudoJet &pj, bool require_two_bs) const{
    unsigned int n_ga_bs;// = 0;
    PseudoJet mdt_jet = m_mdt.result(pj);

    //failed mass drop
    if (mdt_jet.pt() == 0) return -1;

    if (require_two_bs){
        for(unsigned int i_jet=0; i_jet < mdt_jet.pieces().size(); ++i_jet){
            n_ga_bs = 0;
            vector<PseudoJet> cons = mdt_jet.pieces()[i_jet].constituents();
            for (unsigned int i_con=0; i_con<cons.size(); ++i_con){
                int id = cons[i_con].user_index();
                if (id >= 0) continue;
                if(is_qhadron(5, -1 * id)) ++n_ga_bs;
            }
            //failed to have at least one b per subjet
            if (n_ga_bs < 1) return -2;
        }
    }
    double r_filt = min((double)0.3, (double)m_jetRadius/2);
    fastjet::Filter ff(fastjet::JetDefinition(fastjet::cambridge_algorithm,r_filt), 
                       fastjet::SelectorNHardest(3));
    double mm = ff(mdt_jet).m();
    return mm;
}

inline bool 
JetAnalysis::inMassRange(double jetMass) const {
    return m_jetMinM <= jetMass && jetMass <= m_jetMaxM;
}
