#ifndef JETANALYSIS_H
#define JETANALYSIS_H

#include <string>

//physics things josh found 
#include "HepMC/GenEvent.h"

//physics things josh installed
#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "Nsubjettiness/Nsubjettiness.hh"

//physics things josh wrote
#include "JetAnalysis.h"
#include "Calorimeter.h"
#include "PileUp.h"


class Cuts{
 public:
  Cuts() { clear(); }
  
  void add_cut(std::string name) {
    m_passed.push_back(0);
    m_name.push_back(name);
  }
  
  void passed(std::string cut_name = "");

  void clear() {
    m_counter = 0;
    m_passed.clear();
    m_name.clear();
  }

  void begin() { m_counter = 0; }
  
  void print();

 private:
  std::vector<int> m_passed;
  std::vector<std::string> m_name;
  int m_counter;
};

class JetAnalysis{

    private:
        std::vector<PseudoJet> selectGoodJets(fastjet::ClusterSequence &cs);
        inline bool inMassRange(double jetMass) const;


    private:
        int m_saveMaxNJets; 
        int m_cellV; 
        int m_trimJets; 
        float m_jetMinPt; 
        float m_jetMaxPt; 
        float m_jetMinM; 
        float m_jetMaxM; 
        float m_jetRadius; 
        int m_mu;
        int m_nB;
        bool m_doSmear;
        bool m_doWTagger;

        float m_minDR;
        float m_trimRadius;

        fastjet::MassDropTagger m_mdt;
        fastjet::JetDefinition m_primaryDef;
        fastjet::JetDefinition m_trimDef;
        fastjet::Selector m_selPtMin;
        fastjet::Selector m_selPtMax;
        fastjet::Selector m_basePtMin;
        fastjet::Selector m_selDR;
        fastjet::Selector m_selEta;
        std::map<int, fastjet::contrib::Nsubjettiness::Nsubjettiness*> m_nSubCalcs;

        std::string m_desc;

        Calorimeter m_hs_calo;
        PileUpGenerator *m_pug;

        std::vector<int> m_wtagEffs;

        void init(int saveMaxNJets, int trimJets, float jetMinPt, float jetMaxPt, float jetMinM, float jetMaxM, float jetRadius, int nB, int mu, int puSeed, bool doSmear, bool doWTagger, int cellV);
        double get_bdrs_mass(const PseudoJet &pj, bool require_two_bs) const;

        Cuts m_jetCuts;
        Cuts m_pclCuts;

    public:
    // probably also want algorithm
        JetAnalysis(int saveMaxNJets, int trimJets, float jetMinPt, float jetMaxPt, float jetMinM, float jetMaxM, float jetRadius, int nB, int mu, int puSeed, bool doSmear, bool doWTagger, int cellV) { 
          init(saveMaxNJets, trimJets, jetMinPt, jetMaxPt, jetMinM, jetMaxM, jetRadius, nB, mu, puSeed, doSmear, doWTagger, cellV);
        }

        ~JetAnalysis();

        inline static bool is_qhadron(int qid, int pid);
        static bool is_final_qhadron(int qid, HepMC::GenParticle *gp);
        static bool is_charged(int pdg_id);
        void analyze(HepMC::GenEvent* event, map<string, string> props);
        HepMC::GenParticle* findBestParton(PseudoJet &jet, HepMC::GenEvent &event) const;
};
#endif
