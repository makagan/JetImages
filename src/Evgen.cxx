// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Evgen class.
//

#include <iostream>
#include <sstream>
#include "Evgen.h"
#include "ThePEG/Config/HepMCHelper.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

//#include "JetAnalysis.h"



using namespace HerwigWrapper;

Evgen::Evgen() {

    //JetAnalysis(int saveMaxNJets, int trimJets, float jetMinPt, float jetMinM, float jetMaxM, float jetRadius, int mu, int puSeed, bool doSmear, bool doWTagger, int cellV) { 
    //jet_analysis = new JetAnalysis(1, 1, 200, 0, 20000, 1.2, 0, 0, false, false, 0);

    props["source"] = "Herwig++";
    props["ME"]     = "NoME";
    props["seed"]   = "NoSeed";
}

Evgen::~Evgen() {
    //delete jet_analysis; jet_analysis = 0;
}



#ifndef LWH_AIAnalysisFactory_H
#ifndef LWH 
#define LWH ThePEGLWH
#endif
#include "ThePEG/Analysis/LWH/AnalysisFactory.h"
#endif

void Evgen::analyze(tEventPtr event, long ieve, int loop, int state) {
    AnalysisHandler::analyze(event, ieve, loop, state);
    std::cout << endl;
    std::cout << ieve << ", " << loop << ", " << state << std::endl;
    HepMC::GenEvent * hepmc 
        = HepMCConverter<HepMC::GenEvent>::convert(*event, false, GeV, millimeter);

    hepmc->print();

    delete hepmc;

    std::stringstream ss;
    ss << ieve;
    props["eventNumber"] = ss.str();
    //jet_analysis->analyze(hepmc, props);

//  HepMC::GenEvent::particle_const_iterator it;
//  for(it=hepmc->particles_begin(); it!=hepmc->particles_end(); ++it){
//        cout << (*it)->pdg_id() << endl;
//  }
//  hepmc->signal_process_id();
  // Rotate to CMS, extract final state particles and call analyze(particles).
}

LorentzRotation Evgen::transform(tcEventPtr event) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void Evgen::analyze(const tPVector & particles, double weight) {
  AnalysisHandler::analyze(particles);
  // Calls analyze() for each particle.
}

void Evgen::analyze(tPPtr, double weight) {}

void Evgen::dofinish() {
  AnalysisHandler::dofinish();
  // *** ATTENTION *** Normalize and post-process histograms here.
}

void Evgen::doinitrun() {
  AnalysisHandler::doinitrun();

  // *** ATTENTION *** histogramFactory().registerClient(this); // Initialize histograms.
  // *** ATTENTION *** histogramFactory().mkdirs("/SomeDir"); // Put histograms in specal directory.
}


IBPtr Evgen::clone() const {
  return new_ptr(*this);
}

IBPtr Evgen::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).



// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeNoPIOClass<Evgen,AnalysisHandler>
  describeHerwigWrapperEvgen("HerwigWrapper::Evgen", "Evgen.so");


void Evgen::Init() {

  static ClassDocumentation<Evgen> documentation
    ("There is no documentation for the Evgen class");

}

