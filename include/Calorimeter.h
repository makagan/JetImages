#ifndef CALORIMETER_H
#define CALORIMETER_H

#include <set>
#include <vector>
#include <math.h>
#include <string>

#include "fastjet/PseudoJet.hh"

//Peter Loch's radial energy smearing
#include "DetectorModel.h"

using namespace std;
using fastjet::PseudoJet;

class Calorimeter {
    private:
        //double *m_energy;
        vector<double> m_energy;
        unsigned int m_nEtaCells;
        unsigned int m_nPhiCells;
        bool m_phiGoesTo2Pi;
        double m_etaMax;
        double m_dEta;
        double m_dPhi;

        int m_nInteractions;

        bool m_smear;
        static DetectorModel::RadialSmearing m_smearer;
    

        void
        init(unsigned int nEtaCells, unsigned int nPhiCells,
                  bool phiGoesTo2Pi, double etaMax);

        static double
        rectify_phi(double phi, bool goes_to_2PI);

        double etaCloseCell(double eta) const ;

        double phiCloseCell(double phi) const ;

        double
        etaFromIndex(unsigned int index) const ;

        double
        phiFromIndex(unsigned int index) const ;

        int
        makeIndex(double eta, double phi) const;

        unsigned int
        makeEtaIndex(double eta) const ;

        unsigned int 
        makePhiIndex(double phi) const ;

        unsigned int 
        makeIndex(double val, double low, double high, double n) const ;



    public:
        Calorimeter(){
            init(50, 63, false, 2.5);
        }
          
        Calorimeter(unsigned int nEtaCells, unsigned int nPhiCells, 
                    bool phiGoesTo2Pi, double etaMax){
            init(nEtaCells, nPhiCells, phiGoesTo2Pi, etaMax);
        }

        void setSmear(bool value) { 
          m_smear = value; 
          if (m_smear) {
            DetectorModel::Grid gridEM(m_nEtaCells, -m_etaMax, m_etaMax, m_nPhiCells);
            DetectorModel::Grid gridHAD(m_nEtaCells, -m_etaMax, m_etaMax, m_nPhiCells);
            m_smearer.setEmGrid(gridEM);
            m_smearer.setHadGrid(gridHAD);
          }
        }

        const Calorimeter operator+(const Calorimeter &other) const;
        Calorimeter & operator+=(const Calorimeter &other);


        void addEnergy(double e, double eta, double phi);

        void push(const std::vector<PseudoJet> & jets);

        void push(const PseudoJet & jet) {
          std::vector<PseudoJet> jets;
          jets.push_back(jet);
          push(jets);
        }

        void push(double px, double py, double pz, double e, int pdgid) {
          PseudoJet pj(px, py, pz, e);
          pj.set_user_index(pdgid);
          push(pj);
        }

        void
        set_nInteractions(const int & nInts) {m_nInteractions = nInts;}

        int 
        get_nInteractions() const {return m_nInteractions;}

        std::vector<PseudoJet> getNonZeroCells();

        void reset(){
            m_energy = vector<double>(m_nEtaCells * m_nPhiCells, 0.0);
            set_nInteractions(0);
        }

        double get_eta_jet_offset(double eta, double radius) const {
            static const double maxEtaStep = ceil(radius/m_dEta - 0.1);
            return etaCloseCell(eta - maxEtaStep * m_dEta);
        }

        double get_phi_jet_offset(double phi, double radius) const {
            static const double maxPhiStep = ceil(radius/m_dPhi - 0.1);
            return phiCloseCell(phi - maxPhiStep * m_dPhi);
        }

        double
        get_cone_energy(double eta, double phi, double radius) const;

        double
        get_jet_energy(PseudoJet pj) const{
            double total = 0;
            vector<PseudoJet> pjs = pj.constituents();
            for(unsigned int i_pj; i_pj<pjs.size(); ++i_pj){
                total += m_energy[pjs[i_pj].user_index()]; 
            }
            return total;
        }

        string
        toString() const;

        string
        jetToString(double eta, double phi, double jet_radius, std::set<int> whitelist, bool use_whitelist) const;
};

#endif
