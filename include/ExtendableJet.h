#ifndef EXTENDABLEJET_H
#define EXTENDABLEJET_H

#include <map>
#include <string>
#include <sstream>

#include "fastjet/ClusterSequence.hh"

#include "Calorimeter.h"

using namespace std;
using fastjet::PseudoJet;

class ExtendableJet{
    private:
        map<string, string> properties;
        Calorimeter calo;
        vector<Calorimeter> m_pu_calos;
        double m_subjetRadius;
        double m_radius;
        unsigned int m_nSubjets;
        double m_dcut;

        std::set<unsigned int> chargedIds;

        void init(PseudoJet &pJet, 
                        Calorimeter &filledCalo, 
                        vector<Calorimeter> &pu_calos, 
                        double radius,
                        unsigned int nSubjets,
                        double subjetRadius);

    public:
        PseudoJet jet;

        ExtendableJet(){}

//        ExtendableJet(PseudoJet &pJet, 
//                        Calorimeter &filledCalo, double radius){
//            init(pJet, filledCalo, radius, 2, 0.2);
//        }

        ExtendableJet(PseudoJet &pJet, 
                        Calorimeter &filledCalo, 
                        vector<Calorimeter> &pu_calos, 
                        double radius,
                        unsigned int nSubjets,
                        double subjetRadius){
            init(pJet, filledCalo, pu_calos, radius, nSubjets, subjetRadius);
        }

        double radius() const { return m_radius;}

        static bool
        compare(const PseudoJet& aa, const PseudoJet& bb){ return aa.pt() < bb.pt();}

        unsigned int getNTracks() const ;

        template <class T>
        void setProperty(string key, T value){
            stringstream  ss;
            ss << value;
            properties[key] = ss.str();
        }

        void setProperty(map<string, string> user_props){
            map<string, string>::iterator it;
            for(it=user_props.begin(); it != user_props.end(); ++it)
                properties[it->first] = it->second;
        }

        inline string toString() const {toString(1);}
        string toString(int level) const ;
        void operator<<(ostream& out) {out << toString();}
        double get_cos_theta(const vector<PseudoJet> &subjets) const;
};

#endif
