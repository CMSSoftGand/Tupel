#!/bin/bash

#determine CMSSW config
SCRIPT=$(readlink -f $0)
SCRIPTPATH=`dirname $SCRIPT`
ARCH=${SCRIPTPATH##/*/}
WORKDIR=${SCRIPTPATH}/../

cmssw=`$VO_CMS_SW_DIR/cmsset_default.sh`
#configure environment
cd $WORKDIR
#export SCRAM_ARCH=$ARCH
export SCRAM_ARCH=$ARCH
#eval `scram r -sh`
source $cmssw
eval `scram runtime -sh`

#run with the arguments passed
$*
