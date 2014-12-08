#define _USE_MATH_DEFINES

#include <math.h>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <vector>

#include "Calorimeter.h"
#include "PileUp.h"

void 
PileUpGenerator::init(){
    m_refCalo.reset();

    m_max_for_poisson  = poisson(m_mu, (unsigned int) m_mu);
    m_poissonSeed = m_rSeed;
    srand( m_poissonSeed );
}

int 
PileUpGenerator::rand_int(int low, int high){
    int rr = rand();
    srand(rr);

    return (int)((rr/(float)RAND_MAX) * abs(high - low)) + min(high, low);
}

double 
PileUpGenerator::rand_range(double low, double high){
    int rr = rand();
    srand(rr);

    return (rr/(float)RAND_MAX) * abs(high - low) + min(high, low);
}

//Will begin to break for mean > 10^9
//Safe all the way down to and including mean =0
unsigned int 
PileUpGenerator::poisson_sample(double mean){
    const unsigned int max_n = max((int)(6*mean), 2);

    unsigned int nn;
    double yy;

    //disallow possibility of two or more
    if (mean < 1e-3) {
        if (rand_range(0,1) < mean) return 1;    
        else                        return 0;
    }

    //Average iteations through loop:
    //Means >~ 10: 2.4*sqrt(mean) 
    //Means << 1 : 0.8/sqrt(mean)
    for(unsigned int ii=0; ii<=1000*1000; ++ii){
        nn = rand_int(0, max_n);//(unsigned int)(rand() * max_n);
        yy = rand_range(0, m_max_for_poisson);//rand() * m_max_for_poisson;

        if (yy > poisson(nn, mean)) continue;
        
        return nn;
    }

    return 0;
}

double 
PileUpGenerator::poisson(unsigned int nn, double mean){
    double acc = 1;
    for(unsigned int ii=1; ii<=nn; ++ii) acc *= (mean/ii);

    return exp(-1*mean)*acc;
}

PileUpGenerator::PileUpGenerator(double averageIntPerCrossing, 
                        Calorimeter& refCalo,
                        int randomSeed){
    m_mu = averageIntPerCrossing;
    m_refCalo = refCalo;
    m_rSeed = randomSeed;
    init();
}
