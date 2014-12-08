
#define _USE_MATH_DEFINES

#include <math.h>
#include <vector>
#include <string>
#include <sstream>
#include <set>

#include "fastjet/ClusterSequence.hh"

#include "Calorimeter.h"
//#include "CellBooster.h"

using namespace std;
using fastjet::PseudoJet;

DetectorModel::RadialSmearing Calorimeter::m_smearer;

void
Calorimeter::init(unsigned int nEtaCells, unsigned int nPhiCells,
          bool phiGoesTo2Pi, double etaMax){

    m_nEtaCells     = nEtaCells;
    m_nPhiCells     = nPhiCells;
    m_phiGoesTo2Pi  = phiGoesTo2Pi;
    m_etaMax        = etaMax;

    m_dEta          = 2*etaMax/nEtaCells;
    m_dPhi          = 2*M_PI/nPhiCells;

    m_energy = vector<double>(nEtaCells*nPhiCells, 0.0);
    m_nInteractions = 0;

    m_smear         = false;
}

double 
Calorimeter::rectify_phi(double phi, bool goes_to_2PI){
    double ans = phi;
    double phiMin = goes_to_2PI ? 0 : -1 * M_PI;
    double phiMax = goes_to_2PI ? 2*M_PI : M_PI;

    ans += (ans < phiMin ? 2 * M_PI : 0);
    ans -= (ans > phiMax ? 2 * M_PI : 0);
    return ans;
}

double 
Calorimeter::etaCloseCell(double eta) const {
    //double deta = eta - -1*m_etaMax;
    double nCells = (eta - -1*m_etaMax)/m_dEta;
    int purportedCellIndex = floor(nCells);

    return -1 * m_etaMax + m_dEta * (0.5 + purportedCellIndex);
}

double 
Calorimeter::phiCloseCell(double phi) const {
    double phiMin = m_phiGoesTo2Pi ? 0 : -1 * M_PI;
    phi = rectify_phi(phi, m_phiGoesTo2Pi);

    double nCells = (phi - phiMin)/m_dPhi;
    int purportedCellIndex = floor(nCells);

    return phiMin + m_dPhi * (0.5 + purportedCellIndex);
}

double
Calorimeter::etaFromIndex(unsigned int index) const {
    unsigned int ii = index / m_nPhiCells;
    double val = (ii+0.5)*2.0*m_etaMax/m_nEtaCells - m_etaMax;

    return val;
}

double
Calorimeter::phiFromIndex(unsigned int index) const {
    unsigned int ii = index % m_nPhiCells;
    double val      = (ii + 0.5) * 2.0*M_PI/m_nPhiCells;

    if (! m_phiGoesTo2Pi) val -= M_PI;
    return val;
}

int
Calorimeter::makeIndex(double eta, double phi) const{
    if (fabs(eta) > m_etaMax) return -1;

    double pcheck = rectify_phi(phi, m_phiGoesTo2Pi);
    if (fabs(pcheck - phi) > 0.1) return -1;

    return makeEtaIndex(eta)*m_nPhiCells + makePhiIndex(phi);
}

unsigned int
Calorimeter::makeEtaIndex(double eta) const {
    unsigned int index = makeIndex(eta, -1*m_etaMax, m_etaMax, m_nEtaCells);
    return index;
}

unsigned int 
Calorimeter::makePhiIndex(double phi) const {
    unsigned int index;
    if(m_phiGoesTo2Pi)
        index = makeIndex(phi, 0, 2*M_PI, m_nPhiCells);
    else
        index = makeIndex(phi, -1*M_PI, 1*M_PI, m_nPhiCells);

    return index;
}

unsigned int 
Calorimeter::makeIndex(double val, double low, double high, double n) const {
    val = min(val, high);
    val = max(val, low);

    unsigned int index = (n*(val-low)/(high-low));
    if (index == n) --index;
    return index;
} 

void Calorimeter::push(const std::vector<PseudoJet> & jets) {
  std::vector<PseudoJet>::const_iterator jet = jets.begin();
  std::vector<PseudoJet>::const_iterator jet_end = jets.end();
  
  if (m_smear) {
    const std::vector<PseudoJet> & input = m_smearer.smear(jets);
    jet = input.begin();
    jet_end = input.end();
  }

  for (; jet != jet_end; ++jet) {
    if (m_phiGoesTo2Pi) 
      addEnergy(jet->e(), jet->pseudorapidity(), jet->phi());
    else
      addEnergy(jet->e(), jet->pseudorapidity(), jet->phi_std());
  }
}

void Calorimeter::addEnergy(double e, double eta, double phi){
  int ii = makeIndex(eta, phi);
  if (ii > -1) m_energy[ii] += e;
}

const Calorimeter 
Calorimeter::operator+(const Calorimeter &other) const {
    return Calorimeter(*this) += other;
}

Calorimeter&
Calorimeter::operator+=(const Calorimeter &other){
    //check all the indices are equal

    for(unsigned int ii=0; ii<m_energy.size(); ++ii){
        (this->m_energy)[ii] += other.m_energy[ii];
    }
    this->m_nInteractions += other.m_nInteractions;

    return *this;
}


vector<PseudoJet> 
Calorimeter::getNonZeroCells(){
    vector<PseudoJet> cells;
    float e, eta, phi;
    for (unsigned int ii=0; ii < m_nEtaCells * m_nPhiCells; ++ii){

            e   = m_energy[ii];
            if (e<0.01) continue;

            eta = etaFromIndex(ii);
            phi = phiFromIndex(ii);
            PseudoJet pj(e*cos(phi)/cosh(eta),
                         e*sin(phi)/cosh(eta),
                         e*tanh(eta),
                         e);    
            pj.set_user_index(ii);
            cells.push_back(pj);
    }

    return cells;
}

double
Calorimeter::get_cone_energy(double eta, double phi, double radius) const {
    double e_tot = 0;
    double maxEtaStep = ceil(radius/m_dEta - 0.1);
    double maxPhiStep = ceil(radius/m_dPhi - 0.1);
    double t_eta, t_phi;
    int cellIndex;

    //int howmany=0;
    for (int phiStep = -1*maxPhiStep; phiStep <= maxPhiStep; ++phiStep){
        t_phi = phi + phiStep * m_dPhi;
        t_phi = rectify_phi(t_phi, m_phiGoesTo2Pi);

        for (int etaStep = -1*maxEtaStep; etaStep <= maxEtaStep; ++etaStep){
            t_eta = eta + etaStep * m_dEta;

            cellIndex = makeIndex(t_eta, t_phi);
            e_tot += (cellIndex < 0 ? 0 : m_energy[cellIndex]);
            //if (cellIndex >=0) howmany++;
        }
    }

    //cout << "get_cone_energy took: " << howmany << endl;
    return e_tot;
}

string
Calorimeter::toString() const {

    stringstream ss;
    ss << "[";

    for (int ii=0; ii<m_nEtaCells * m_nPhiCells; ++ii)
        ss << m_energy[ii] << ", ";

    ss << "]";
    return ss.str();
}

string
Calorimeter::jetToString(double eta, double phi, double jet_radius, std::set<int> whitelist, bool use_whitelist) const {
    double maxEtaStep = ceil(jet_radius/m_dEta - 0.1);
    double maxPhiStep = ceil(jet_radius/m_dPhi - 0.1);
    double t_eta, t_phi, ee;
    double rel_eta;
    double rel_phi;
    double origin_eta = get_eta_jet_offset(eta, jet_radius);
    double origin_phi = get_phi_jet_offset(phi, jet_radius);

    int cellIndex;
    unsigned int nCells = 0;


    stringstream  ss, sFirst;
    ss << "'cells': [";


    for (int phiStep = -1*maxPhiStep; phiStep <= maxPhiStep; ++phiStep){
        for (int etaStep = -1*maxEtaStep; etaStep <= maxEtaStep; ++etaStep){
            t_eta = eta + etaStep * m_dEta;
            t_phi = phi + phiStep * m_dPhi;
            t_phi = rectify_phi(t_phi, m_phiGoesTo2Pi);
            
            cellIndex = makeIndex(t_eta, t_phi);
            if (cellIndex < 0 || (use_whitelist && whitelist.count(cellIndex) == 0)) ee = 0;
            else ee = m_energy[cellIndex];

            rel_eta = etaCloseCell(t_eta) - origin_eta;
            rel_phi = phiCloseCell(t_phi) - origin_phi;
            rel_phi += (rel_phi < -1 * m_dPhi ? 2*M_PI : 0);

            ss << "("  << ee/cosh(etaCloseCell(t_eta));
            ss << ", " << rel_eta;
            ss << ", " << rel_phi;
            ss << "),";
            ++nCells;
        }
    }

    ss << "], 'n_cells': '" << nCells << "',";
    ss << "'jet_offset': (" << origin_eta << ", " << origin_phi;
    ss << "), ";
    return ss.str();
}
