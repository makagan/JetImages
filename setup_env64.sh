#!/bin/bash

setup_g64python() {
    export PYTHONHOME=/afs/slac/g/glast/ground/GLAST_EXT/redhat5-x86_64-64bit-gcc41/python/2.7.2
    export PATH=$PYTHONHOME/bin:$PATH
}

setup_64sklearn() {
    export PYTHONPATH=/u/eb/joshgc/mynfs/scikit-learn/install_64/lib/python2.7/site-packages:$PYTHONPATH
    setup_g64python
}

setup_csjets() {
    setup_64sklearn
}

setup_csjets
