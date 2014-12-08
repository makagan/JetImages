#ifndef PILEUP_H
#define PILEUP_H

#define _USE_MATH_DEFINES

#include <math.h>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <vector>

#include "Calorimeter.h"

class PileUpGenerator {
    
    protected:
        double m_mu;
        double m_max_for_poisson;
        int m_rSeed;
        int m_poissonSeed;
        Calorimeter m_refCalo;

        virtual void init();
        int rand_int(int low, int high);
        double rand_range(double low, double high);

        //Will begin to break for mean > 10^9
        //Safe all the way down to and including mean =0
        unsigned int poisson_sample(double mean);
        double poisson(unsigned int nn, double mean);

    public:
        PileUpGenerator(double averageIntPerCrossing, 
                        Calorimeter& refCalo,
                        int randomSeed);
        ~PileUpGenerator(){}

        double mu() const {return m_mu;}
        virtual Calorimeter next(int exactIntPerCrossing=-1) = 0;
        virtual vector<Calorimeter> next_vector(int exactIntPerCrossing=-1) = 0;
};

#endif
