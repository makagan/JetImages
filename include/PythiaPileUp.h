#ifndef PYTHIAPILEUP_H
#define PYTHIAPILEUP_H

#include <sstream>
#include "Pythia.h"
#include "PileUp.h"
#include "Calorimeter.h"


class PythiaPileUpGenerator : public PileUpGenerator{

    private:
        void init();

    public:
        Pythia8::Pythia pythia;

        PythiaPileUpGenerator(double averageIntPerCrossing, 
                             Calorimeter& refCalo,
                             int randomSeed);
        virtual Calorimeter next(int exactIntPerCrossing);
        virtual vector<Calorimeter> next_vector(int exactIntPerCrossing);
};
#endif
