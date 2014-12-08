#!/bin/bash

setup_gpython() {
    export PYTHONHOME=/afs/slac/g/glast/ground/GLAST_EXT/redhat5-i686-32bit-gcc41/python/2.7.2
    export PATH=$PYTHONHOME/bin:$PATH
}
setup_g64python() {
    export PYTHONHOME=/afs/slac/g/glast/ground/GLAST_EXT/redhat5-x86_64-64bit-gcc41/python/2.7.2
    export PATH=$PYTHONHOME/bin:$PATH
}

setup_pythia() {
    #export PYTHIA8LOCATION=/afs/slac/g/atlas/c/sw/lcg/external/MCGenerators/pythia8/150/i686-slc5-gcc43-opt/
    export PYTHIA8LOCATION=/u/at/estrauss/afs/WTagger/Pythia8176/
    export PYTHIA8DATA=${PYTHIA8LOCATION}xmldoc/
    export LD_LIBRARY_PATH=${PYTHIA8LOCATION}lib/:$LD_LIBRARY_PATH
}

setup_hepmc() {
    export LD_LIBRARY_PATH=/u/eb/joshgc/mynfs/HepMC/install/lib:$LD_LIBRARY_PATH
    export PATH=${PWD}:${PATH}
}

setup_gcc32() {
    export PATH=/afs/slac/g/atlas/d/gcc-alt-435/i686-slc5-gcc43-opt/bin:$PATH
}

setup_boost() {
    export BOOSTINCDIR=/usr/include
    export BOOSTLIBLOCATION=/usr/lib
}

setup_sherpa() {
    export SHERPA=/afs/slac/g/atlas/c/sw/lcg/external/MCGenerators/sherpa/1.3.0.2/i686-slc5-gcc43-opt
    export SHERPA_INCLUDE_PATH=$SHERPA/include/SHERPA-MC
    export SHERPA_SHARE_PATH=$SHERPA/share/SHERPA-MC
    export SHERPA_LIBRARY_PATH=$SHERPA/lib/SHERPA-MC
    export LD_LIBRARY_PATH=$SHERPA_LIBRARY_PATH:$LD_LIBRARY_PATH

    setup_gcc32
}

setup_madgraph() {
    export MG5=/u/at/estrauss/afs/WTagger/MadGraph5_v1_5_11
}

setup_sklearn() {
    export PYTHONPATH=/u/eb/joshgc/mynfs/scikit-learn/install/lib/python2.7/site-packages:$PYTHONPATH
    setup_gpython
}

setup_photos() {
    setup_gcc32
    setup_pythia
}

setup_feynhiggs() {
    fpath=/u/eb/joshgc/mynfs/FeynHiggs/FeynHiggs-2.9.1
    export MANPATH=$fpath/man/:$MANPATH
    export PATH=$fpath/i686-Linux/bin:$PATH
}

setup_root() {
    source /afs/slac/g/glast/ground/GLAST_EXT/redhat5-i686-32bit-gcc41/ROOT/v5.34.03/bin/thisroot.sh
}

setup_fastjet() {
    export LIBPATH=/afs/slac.stanford.edu/package/vdt/vol5/globus/lib:/usr/lib:/lib:$LIBPATH
    export LD_LIBRARY_PATH=/afs/slac.stanford.edu/package/vdt/vol5/atlasosgcompat/keepexpat/lib:/afs/slac.stanford.edu/package/vdt/vol5/openldap/lib:/afs/slac.stanford.edu/package/vdt/vol5/lcg/lib:/afs/slac.stanford.edu/package/vdt/vol5/curl/lib:/afs/slac.stanford.edu/package/vdt/vol5/glite/lib64:/afs/slac.stanford.edu/package/vdt/vol5/glite/lib:/afs/slac.stanford.edu/package/vdt/vol5/globus/lib:/afs/slac.stanford.edu/package/vdt/vol5/berkeley-db/lib:/afs/slac.stanford.edu/package/vdt/vol5/expat/lib:/afs/slac/g/atlas/d/gcc-alt-435/x86_64-slc5-gcc43-opt/lib64:/afs/slac/g/atlas/d/gcc-alt-435/x86_64-slc5-gcc43-opt/lib:/afs/slac/g/atlas/d/gcc-alt-432/x86_64-slc5-gcc43-opt/lib64:/afs/slac/g/atlas/d/gcc-alt-432/x86_64-slc5-gcc43-opt/lib:/afs/slac.stanford.edu/package/vdt/vol5/lcg/lib64:$LD_LIBRARY_PATH
}

setup_herwig() {
    export PATH=~joshgc/mynfs/Herwig/install-dir/bin/:${PATH}
    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/afs/slac/g/atlas/c/sw/lcg/external/MCGenerators/lhapdf/5.8.5/i686-slc5-gcc43-opt/lib
    export LHAPATH=/afs/slac/g/atlas/c/sw/lcg/external/MCGenerators/lhapdf/5.8.5/share/PDFsets/
}

setup_csjets() {
    setup_boost
    setup_pythia
    setup_madgraph
    setup_gcc32
    setup_sklearn
    setup_hepmc
    setup_root
    setup_fastjet
    setup_herwig

}

setup_csjets
