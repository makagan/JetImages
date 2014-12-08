#!/bin/bash

echo "Begin at `date`"

[ "$USER" == "estrauss" ] && WorkDir=${PWD} && InputDir=/u/at/estrauss/afs/WTagger/logs_bhadron && OutputDir=/u/at/estrauss/afs/WTagger/npy_madgraph_bhadron

SubFileLoc=`pwd`/_batchRotate.sh

echo '
#!/bin/bash

INDEX=`echo $LSB_BATCH_JID | sed -e "s/^[0-9]*//; s/\[//; s/\]//"`
BASENAME=$2
OUTPUT=$3/`basename ${BASENAME}`

echo "Index: $INDEX"
echo "Basename: $BASENAME"
echo "Output: $OUTPUT"

cd $1

if [ `ls ${BASENAME}_${INDEX}.* | wc -l` == 0]; then 
  echo "No such file ${BASENAME}_${INDEX}.*"
  exit 0
fi

. setup_env.sh
T=$(($RANDOM%600))
echo "Sleep for $T seconds"
sleep $T
echo python rotate_jets.py -n -1 -i ${BASENAME}_${INDEX}.* -o ${OUTPUT}_${INDEX}.npy
python rotate_jets.py -n -1 -i ${BASENAME}_${INDEX}.* -o ${OUTPUT}_${INDEX}.npy
' > $SubFileLoc
chmod u+x $SubFileLoc

mkdir -p ${OutputDir}

SAMPLES=`ls $InputDir/* | sed -e 's/_[0-9]*.\(log\|tgz\)//' | sort -u`

for SAMPLE in $SAMPLES; do
    echo `bjobs | wc -l`
    while [[ `bjobs|wc -l` -ge 120 ]]; do
        echo "More than 120 jobs on the queues, waiting for some to clear"
        sleep 300
    done
    echo "Submit $SAMPLE"
    N=`ls ${SAMPLE}* | wc -l`
    bsub -q medium -J `basename ${SAMPLE}`[1-$N] -R rhel50 $SubFileLoc $WorkDir $SAMPLE $OutputDir
    sleep 2
done

echo "Done at `date`"
