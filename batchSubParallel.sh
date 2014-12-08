#!/bin/bash

echo "Begin at `date`"

[ "$USER" == "joshgc" ] && WorkDir=/afs/slac.stanford.edu/u/eb/joshgc/trunk/CSJets && LogDir=/u/eb/joshgc/mynfs/CSJets/logs
[ "$USER" == "estrauss" ] && WorkDir=${PWD} && LogDir=/u/at/estrauss/afs/WTagger/logs_subjobs/

SubFileLoc=${WorkDir}/_parallelMaster.sh
DateSuffix=`date +%s`

rad=1.2
gen="p"
[[ "$gen" == "p" ]] && prog="./pythiaJets.exe"
[[ "$gen" == "h" ]] && prog="./runHerwig.py"


for metype in w v g; do
for trim in 0 1; do
for mu in 0 30; do
    Queue=long
    nevents=50000
    nsubjobs=10
    njobs=15
    [[ "$metype" == "g"  ]] && nevents=$(($nevents * 50))
    [[ "$metype" == "v"  ]] && nevents=$(($nevents * 100))
    [[ "$metype" == "w"  ]] && nevents=$(($nevents * 20))
    [[ "$mu" == "30" ]] && njobs=$(($njobs*4))
    nevents=$(($nevents / ($njobs*$nsubjobs) ))
for pt_start in 200 500; do
LogPrefix=${LogDir}/bsub_${DateSuffix}_${gen}_${metype}_${rad}_${pt_start}_${mu}_T${trim}
mkdir -p `dirname $LogPrefix`
echo
echo "Submitting $njobs jobs each with $nsubjobs subjobs running $nevents events to $Queue with prefix $LogPrefix"
echo $LogPrefix
bsub -J "${LogPrefix}[1-${njobs}]" -R "rhel50" -n ${nsubjobs} -q ${Queue} $SubFileLoc                  \
    ${LogPrefix}                                        \
    ${WorkDir}                                          \
    ${nsubjobs}                                         \
    ${prog}             --nevents ${nevents}            \
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

