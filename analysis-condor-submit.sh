#!/bin/bash

TIME=`date +'%Y-%m-%d-%H%M'`
# DATA="DoubleEle_A_190456_208686 DoubleEle_B_190456_208686 DoubleEle_C_190456_208686 DoubleEle_D_190456_208686 DoubleMu_A_190456_208686 DoubleMu_B_190456_208686 DoubleMu_C_190456_208686 DoubleMu_D_190456_208686 EleMu_A_190456_208686 EleMu_B_190456_208686 EleMu_C_190456_208686 EleMu_D_190456_208686"
# MC="Tbar_schan Tbar_tchan Tbar_tWchan T_schan TTbarJets T_tchan T_tWchan WGamma WW WZ ZGamma ZJets_2l1 ZJets_2l2 ZJets_2l3 ZJets_2l4 ZZ"
# Find files in directory and print their path-stripped file names with spaces in between.
MC=`find /hadoop/store/user/yoshin/DMNtuples/processed/MC -type f -printf '%f '`
DATA=`find /hadoop/store/user/yoshin/DMNtuples/processed/Data -type f -printf '%f '`
WZ=`find /hadoop/store/user/yoshin/DMNtuples/processed/MC -type f -name 'WZ*' -printf '%f ' | grep WZ`
DIR=AnaNaS

mkdir -p /data/users/yhshin/condor/${DIR}/${TIME}
mkdir -p ./condor

for samplefile in $WZ
do
	sample=`basename $samplefile .root`
	TARGET_DIR=/data/users/yhshin/DMSkims/${TIME}/${sample}/
	mkdir -p $TARGET_DIR
	cat ./analysis-condor.jdl \
	| sed -e s~PREFIX_NAME~/data/users/yhshin/condor/${DIR}/${TIME}~ \
	| sed -e s~RUN_DIR~/home/yhshin/CMSSW_5_3_9_patch3/src/AnaNaS~ \
	| sed -e s~SAMPLE_NAME~${sample}~ \
	| sed -e s~TARGET_DIR~${TARGET_DIR}~ \
	> condor/analysis-${sample}.jdl
	echo "condor_submit condor/analysis-${sample}.jdl"
	condor_submit condor/analysis-${sample}.jdl
done
