#!/bin/bash

echo "Begin at `date`"

[ "$USER" == "estrauss" ] && HomeDir=${PWD} && InputDir=/u/at/estrauss/atlint02/logs_madgraph/ && OutputDir=/u/at/estrauss/atlint02/logs_madgraph_v2_fixed/

SubFileLoc=`pwd`/_batchRerunLHE.sh

echo '
#!/bin/bash

INDEX=`echo $LSB_BATCH_JID | sed -e "s/^[0-9]*//; s/\[//; s/\]//"`

BASENAME=$2
OUTPUT=$3/`basename ${BASENAME}`

ID=`uuidgen`
[[ -e /tmp ]] && WORKDIR=/tmp/${USER}/${ID}
[[ -e /scratch ]] && WORKDIR=/scratch/${USER}/${ID}

echo "Index: $INDEX"
echo "Basename: $BASENAME"
echo "Output: $OUTPUT"
echo "WorkDir: $WORKDIR"

cd $1

if [ `ls ${BASENAME}_${INDEX}.* | wc -l` == 0 ]; then 
  echo "No such file ${BASENAME}_${INDEX}.*"
  exit 0
fi

. setup_env.sh
T=$(($RANDOM%600))
echo "Sleep for $T seconds"
sleep $T

mkdir -p ${WORKDIR}
if [[ "$?" != "0" ]]; then
   echo "Could not create WorkDir"
   exit 1
fi

cd ${WORKDIR}

LOGNAME=`tar -xzvf ${BASENAME}_${INDEX}.*`
awk "/Start LHE File/{x=\"F\"++i;next}{print  > \"${LOGNAME}_split\"x;}" ${LOGNAME}

for FILE in `ls ${LOGNAME}_splitF*`; do
  echo "Process ${FILE}"
  echo "@@ Start LHE File" >> regen.log
  awk "//{f=1;}f&&/End LHE File/{exit}f" ${FILE} >> regen.log
  echo "@@ END LHE File" >> regen.log
  `grep "Called as" ${FILE} | sed -e "s/Called as: //; s/LHEFile=.*lhe/LHEFile=${FILE}/"` >> regen.log
done

tar -czvf ${OUTPUT}_${INDEX}.tgz regen.log
cd $1
rm -rf ${WORKDIR}

' > $SubFileLoc
chmod u+x $SubFileLoc

mkdir -p ${OutputDir}

SAMPLES=`ls $InputDir/* | sed -e 's/_[0-9]*.\(log\|tgz\)//' | sort -u`

for SAMPLE in $SAMPLES; do
    while [[ `bjobs|wc -l` -ge 100 ]]; do
        echo "More than 100 jobs on the queues, waiting for some to clear"
        sleep 300
    done
    echo "Submit $SAMPLE"
    N=`ls ${SAMPLE}* | wc -l`
    bsub -q long -J `basename ${SAMPLE}`[1-$N] -R rhel50 $SubFileLoc $HomeDir $SAMPLE $OutputDir
    sleep 30
done

echo "Done at `date`"
