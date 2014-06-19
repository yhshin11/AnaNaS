#!/bin/bash
#
# variables from arguments string in jdl
#

RUN_DIR=$1
SAMPLE_NAME=$2
TARGET_DIR=$3

#
# header 
#
source ~/.bashrc
cd $RUN_DIR
eval `scramv1 runtime -sh`
source ./setup
analysis -a DM -f ~/hadoop/DMNtuples/processed/MC/${SAMPLE_NAME}.root -o ${TARGET_DIR}
