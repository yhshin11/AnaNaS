#!/bin/bash
#
# variables from arguments string in jdl
#

RUN_DIR=$1
SAMPLE_NAME=$2

#
# header 
#
cd $RUN_DIR
source /home/yhshin/.bashrc
eval `scramv1 runtime -sh`
source ./setup
analysis -a DM -f ~/hadoop/DMNtuples/$SAMPLE_NAME.root -o workdir/root/DM/$SAMPLE_NAME/
