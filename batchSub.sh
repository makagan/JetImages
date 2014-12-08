#!/bin/bash

echo "Begin at `date`"

[ "$USER" == "joshgc" ] && WorkDir=/afs/slac.stanford.edu/u/eb/joshgc/trunk/CSJets && LogDir=/u/eb/joshgc/mynfs/CSJets/logs
[ "$USER" == "estrauss" ] && WorkDir=${PWD} && LogDir=/u/at/estrauss/afs/WTagger/logs_17052013/

SubFileLoc=`pwd`/_batchSingleSub.sh
DateSuffix=`date +%s`

echo '#!/bin/bash
echo CD to $1
echo CMD is $2

cd $1
source setup_env.sh
cmd=$2

shift
shift
echo Calling $cmd $*
$cmd $*' > $SubFileLoc
chmod u+x $SubFileLoc

rad=1.2
gen="p"

total_nevents=50000

for metype in w v g; do
for trim in 0 1; do
for mu in 0 30; do
    Queue=long
    nevents=$((total_nevents * 10))
    njobs=100
    #Correct number of events for approx efficiency of a 10 GeV mass window
    [[ "$metype" == "g"  ]] && nevents=$((total_nevents * 50))
    [[ "$metype" == "v"  ]] && nevents=$((total_nevents * 100))
    [[ "$metype" == "w"  ]] && nevents=$((total_nevents * 20))
    #mu = 30 jobs run longer
    [[ "$mu" == "30" ]] && njobs=$((njobs*5))
    nevents=$((nevents/njobs + 1))
for pt_start in 200 500 ; do
LogPrefix=${LogDir}/bsub_${DateSuffix}_${gen}_${metype}_${rad}_${pt_start}_${mu}_T${trim}_ID
mkdir -p `dirname $LogPrefix`
echo
echo "Submitting $njobs jobs each with $nevents events to $Queue with prefix $LogPrefix"
echo $LogPrefix
bsub -R rhel50 -q ${Queue} -o "${LogPrefix}%I.out" -J "${LogPrefix}[1-${njobs}]" $SubFileLoc\
        ${WorkDir} \
        ./pythiaJets.exe    --nevents ${nevents}            \
                            --savemaxnjets 1                \
                            --metype $metype                \
                            --mu $mu                        \
                            --rad $rad                      \
                            --trimjet ${trim}               \
                            --wtagger                       \
                            --lepcharge 1                   \
                            --meminpt $((pt_start-50))      \
                            --memaxpt $((pt_start+100))     \
                            --jetminpt $pt_start            \
                            --jetmaxpt $((pt_start+50))     \
                            --jetminm 55                    \
                            --jetmaxm 105

done
done
done
done

echo "End at `date`"

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
