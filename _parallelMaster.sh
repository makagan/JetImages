#!/bin/bash

setup_job() {
    echo "Setting up job"
    cd ${SourceDir}
    . setup_env.sh
    
    ID=`uuidgen`
    [[ -e /tmp ]] && WorkDir=/tmp/${USER}/${ID}
    [[ -e /scratch ]] && WorkDir=/scratch/${USER}/${ID}

    echo "WorkDir = ${WorkDir}"

    mkdir -p ${WorkDir}
    cd ${WorkDir}
}

cleanup_job() {
    echo "Cleaning up job"
    cd ${WorkDir}
    echo "================"
    echo "ls $PWD"
    ls
    echo "================"
    find . -name '*.gen.log' -print0 | xargs -0 -I % sh -c 'echo == % ==; cat %' > ${LocalLogFile}.log
    tar -czvf ${LocalLogFile}.tgz ${LocalLogFile}.log
    echo mv ${LocalLogFile}.tgz ${LogFile}.tgz
    mv ${LocalLogFile}.tgz ${LogFile}.tgz
    cd ..
    rm -rf ${WorkDir}
}

run_subjobs() {
    cd
    PIDS=""
    for HOST in $LSB_HOSTS; do
        ID=`uuidgen`
        (lsrun -m "$HOST" ${SourceDir}/_parallelSlave.sh ${SourceDir} ${Cmd} ${Args}) &> ${WorkDir}/${ID}.gen.log &
        PIDS="$PIDS $!"
    done
    cd -
}

mill_around() {
    echo "Waiting for subjobs to end with PID: ${PIDS}"
    Running=1
    while [[ ${Running} -eq 1 ]]; do
        sleep $(($RANDOM%300))
        Running=0
        for PID in $PIDS; do
            kill -0 ${PID} &> /dev/null && Running=1
        done
    done
}


LogFile=${1}
shift
SourceDir=${1}
shift
NSubJobs=${1}
shift
Cmd=${SourceDir}/${1}
shift
Args=$*

LogFile="${LogFile}_`echo $LSB_BATCH_JID | sed -e 's/\[/_/; s/\]//'`"
LocalLogFile=`basename ${LogFile}`

echo "Begin at `date`"
echo SourceDir = $SourceDir
echo LogFile = ${LogFile}
echo Cmd = $Cmd
echo Args = $Args

echo "Subjob hosts: $LSB_HOSTS"
setup_job
run_subjobs ${NSubJobs}

trap cleanup_job SIGHUP SIGINT SIGTERM

mill_around
cleanup_job

echo "Let's blow this popsicle stand"
echo "End at `date`"
