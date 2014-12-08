#!/bin/bash

[ "$USER" == "joshgc" ] && WorkDir=/afs/slac.stanford.edu/u/eb/joshgc/trunk/CSJets && LogDir=/u/eb/joshgc/mynfs/CSJets/logs
[ "$USER" == "estrauss" ] && WorkDir=${PWD} && LogDir=/u/at/estrauss/afs/WTagger/logs_herwig_v4

SubFileLoc=`pwd`/_batchSingleSubHepMC.sh
DateSuffix=`date +%s`

echo '#!/bin/bash
set -e
echo WORKDIR to $1
echo CMD is $2

ID=`uuidgen`
mkdir -p /scratch/${USER}/${ID}

cd $1
source setup_env.sh
cd /scratch/${USER}/${ID}
cmd=$1/$2

shift
shift
echo Calling $cmd $*
$cmd $*

cd /scratch
rm -rf /scratch/${USER}/${ID}
' > $SubFileLoc
chmod u+x $SubFileLoc

rad=1.2
for metype in v w; do
for mu in 0 ; do
    Queue=long
    nevents=50000
    njobs=20
    [ "$mu" == "30" ] && Queue=long && nevents=5000  && njobs=1000
    [ "$mu" == "0"  ] && Queue=medium && nevents=10000  && njobs=40
    [ "$metype" == "v" ] && njobs=$((njobs*4))
for pt_start in 200 500; do
LogPrefix=${LogDir}/bsub_${metype}_${pt_start}_${mu}_${rad}_${DateSuffix}_
mkdir -p `dirname $LogPrefix`
echo
echo "Submitting $njobs jobs each with $nevents events to $Queue"
echo $LogPrefix
for (( ii=1; ii<=$njobs; ii++ )) ;  do
    echo $ii
    bsub -R rhel50 -q ${Queue} -o $LogPrefix${ii}.log $SubFileLoc           \
        ${WorkDir} \
        ./runHerwig.py    --nevents ${nevents}            \
                          --seed ${ii}                    \
                          --savemaxnjets 1                \
                          --metype $metype                \
                          --mu $mu                        \
                          --rad $rad                      \
                          --trimjet 1                     \
                          --wtagger                       \
                          --lepcharge 1                   \
                          --meminpt $((pt_start-50))             \
                          --memaxpt $((pt_start+100))      \
                          --jetminpt $pt_start            \
                          --jetmaxpt $((pt_start+50))     \
                          --jetminm 60                    \
                          --jetmaxm 100

done
done
done
done

#files=(~/mynfs/CSJets/logs/bsub_*.log)
#
#for in_file in ${files[*]} ; do
#    out_file=${in_file%.*}
#
#    echo $in_file $out_file
#    bsub -q long -o ${out_file}_rotate.log $SubFileLoc \
#        /afs/slac.stanford.edu/u/eb/joshgc/trunk/CSJets \
#        python rotate_jets.py -1 $in_file $out_file.npy
#
#done
