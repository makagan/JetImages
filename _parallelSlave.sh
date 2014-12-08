#!/bin/bash

uname -a

export HomeDir=${1}
shift
export Cmd=${1}
shift
export Args=$*

echo "HomeDir: $HomeDir"
echo "Cmd    : $Cmd"
echo "Args   : $Args"

ID=`uuidgen`
echo $HOSTNAME} ${PWD}

[[ -e /tmp ]] && WorkDir=/tmp/${USER}/${ID}
[[ -e /scratch ]] && WorkDir=/scratch/${USER}/${ID}

cd ${HomeDir}
. setup_env.sh

mkdir -p ${WorkDir}
cd ${WorkDir}
${Cmd} ${Args}
cd -
rm -rf ${WorkDir}

echo "${ID} Done"
