
#include <assert.h>
#include <map>
#include <set>
#include <string>
#include <sstream>
#include <algorithm>

#include "ExtendableJet.h"

using namespace std;
using fastjet::PseudoJet;

void 
ExtendableJet::init(PseudoJet &pJet, 
                Calorimeter &filledCalo, 
                vector<Calorimeter> &pu_calos, 
                double radius,
                unsigned int nSubjets,
                double subjetRadius){
    jet  = pJet;
    //evt  = event;
    calo = filledCalo;
    m_pu_calos = pu_calos;
    m_nSubjets = nSubjets;
    m_subjetRadius = subjetRadius;
    m_radius = radius;
    //http://www.fastjet.fr/repo/doxygen-3.0.3/10-subjets_8cc_source.html
    m_dcut = pow(subjetRadius/radius, 2); //don't know what dcut means

    setProperty("radius", radius);
    setProperty("subjetRadius", subjetRadius);
    setProperty("nParts", jet.constituents().size());
    setProperty("pt" , jet.pt ());
    setProperty("eta", jet.eta());
    setProperty("phi", jet.phi_std());
    setProperty("rap", jet.rap());
    setProperty("e"  , jet.e  ());
    setProperty("m"  , jet.m ());

    //setProperty("nTracks", getNTracks());
    setProperty("npv", calo.get_nInteractions());


    unsigned int _idArray[] = {211, 321, 2212, 11, 13};
    chargedIds = set<unsigned int>(_idArray, _idArray + 5);
}

////this is broken at the moment, since
//// user_id was hijacked for calo index
//unsigned int 
//ExtendableJet::getNTracks() const {
//    return 0;
//    unsigned int pdgId   = 0;
//    unsigned int nTracks = 0;
//    int          index   = 0;
//
//    std::vector<PseudoJet> parts = this->jet.constituents();
//    std::vector<PseudoJet>::iterator it;
//    for (it=parts.begin(); it!=parts.end(); it++){
//
//        index = (*it).user_index();
//        if(index < 0) continue;
//
//        pdgId = abs(this->evt[index].id());
//        if (chargedIds.count(pdgId) == 1) ++nTracks;
//    }
//
//    return nTracks;
//}
//

double
ExtendableJet::get_cos_theta(const vector<PseudoJet> &subjets) const
{
    assert(subjets.size()>1);

    PseudoJet boost = subjets[0] + subjets[1];
    PseudoJet unb_0 = PseudoJet(subjets[0]).unboost(boost);
    PseudoJet unb_1 = PseudoJet(subjets[1]).unboost(boost);

    double ct_0 = 
        unb_0.px() * boost.px() + 
        unb_0.py() * boost.py() + 
        unb_0.pz() * boost.pz();
    double ct_1 = 
        unb_1.px() * boost.px() + 
        unb_1.py() * boost.py() + 
        unb_1.pz() * boost.pz();

    ct_0 /= sqrt(
        unb_0.px() * unb_0.px() + 
        unb_0.py() * unb_0.py() + 
        unb_0.pz() * unb_0.pz());
    ct_1 /= sqrt(
        unb_1.px() * unb_1.px() + 
        unb_1.py() * unb_1.py() + 
        unb_1.pz() * unb_1.pz());

    double boost_mag = sqrt(
        boost.px() * boost.px() + 
        boost.py() * boost.py() + 
        boost.pz() * boost.pz());

    ct_0 /= boost_mag;
    ct_1 /= boost_mag;

    // ct_0 > ct_1 98% of the time in W-200GeV r=1.2
    //cout << "#B" << ct_0 << ", " << ct_1 << endl;

    return ct_0;
}

//level==0 no cells
//level==1 only whitelisted cells
//level==2 all cells in jet regardless of trimming
string 
ExtendableJet::toString(int verbosity_level) const {
    map<string, string>::const_iterator it;

    std::stringstream ss;
    int pdgId       = 0;
    int index       = 0;

    double origin_eta = calo.get_eta_jet_offset(jet.eta(), m_radius);
    double origin_phi = calo.get_phi_jet_offset(jet.phi(), m_radius);
    double rel_phi;

    ss << "{";
    for (it=properties.begin() ; it != properties.end(); it++)
    {
        ss << "'" << (*it).first << "': " ;
        ss << "'" << (*it).second << "', ";
    }

    if (verbosity_level>0){
        std::set<int> whitelist;
        for(unsigned int iCon = 0; iCon < jet.constituents().size(); ++iCon)
        {
            whitelist.insert(jet.constituents()[iCon].user_index());
        }

        ss << calo.jetToString(jet.eta(),
                               jet.phi_std(),
                               m_radius,
                               whitelist,
                               verbosity_level==1
                               );
    }

    //vector<PseudoJet> subjets = fastjet::sorted_by_pt(
    //                            jet.exclusive_subjets(m_dcut));
    vector<PseudoJet> subjets = jet.pieces();

    subjets = fastjet::sorted_by_pt(jet.pieces());

    ss << "'cos_theta': '";
    if (subjets.size() > 1) ss << get_cos_theta(subjets);
    else                    ss << -2;
    ss << "', ";

    ss << "'subjet_dr': '";
    if (subjets.size() < 2) ss << -1;
    else ss << subjets[0].delta_R(subjets[1]);
    ss << "', ";

    ss << "'y_balance': '";
    if (subjets.size() < 2) ss << -1;
    else ss << subjets[1].pt() * subjets[0].delta_R(subjets[1])/jet.m();
    ss << "', ";

    ss << "'log_balance': '";
    if (subjets.size() < 2) ss << -1;
    else ss << -1 * log(1-subjets[0].pt()/jet.pt());
    ss << "', ";

    ss << "'subjets': [";
    for (unsigned int ii=0; ii < min(m_nSubjets, subjets.size()); ++ii){
        PseudoJet sj = subjets[ii];
        ss << "("  << sj.pt();
        ss << ", " << sj.eta();
        ss << ", " << sj.phi_std();
        ss << "),";
    }
    ss << "],";


    ss << "'rel_subjets': [";
    for (unsigned int ii=0; ii < min(m_nSubjets, subjets.size()); ++ii){
        PseudoJet sj = subjets[ii];

        rel_phi = sj.phi_std() - origin_phi;
        rel_phi += (rel_phi < 0 ? 2*M_PI : 0);

        ss << "("  << sj.pt();
        ss << ", " << sj.eta()-origin_eta;
        ss << ", " << rel_phi;
        ss << "),";
    }

    vector<double> pu_energies;
    double pu_e = 0;
    for(unsigned int i_calo=0; i_calo < m_pu_calos.size(); ++i_calo){
        pu_energies.push_back(m_pu_calos[i_calo].get_jet_energy(jet));
        pu_e += pu_energies[i_calo];
    }

    ss << "], 'pu_energies': [";
    for(unsigned int i_e=0; i_e<pu_energies.size(); ++i_e) 
        ss << pu_energies[i_e] << ", ";
    ss << "]";
    ss << ", 'pu_energy': " << pu_e;
    ss << "}";
    return ss.str();
}
