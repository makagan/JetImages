#!/bin/bash

echo "You must source this script.  Do not ./thisScript"
export PYTHIA8LOCATION=/afs/slac/g/atlas/c/sw/lcg/external/MCGenerators/pythia8/150/i686-slc5-gcc43-opt/
export PYTHIA8DATA=${PYTHIA8LOCATION}xmldoc/
export BOOSTINCDIR=/usr/include
export BOOSTLIBLOCATION=/usr/lib

export PYTHONPATH=/u/eb/joshgc/mynfs/scikit-learn/install/lib/python2.7/site-packages:$PYTHONPATH
export PYTHONHOME=/afs/slac/g/glast/ground/GLAST_EXT/redhat5-i686-32bit-gcc41/python/2.7.2
export PATH=$PYTHONHOME/bin:$PATH
