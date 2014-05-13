#!/bin/bash

DATA="DoubleEle_A_190456_208686 DoubleEle_B_190456_208686 DoubleEle_C_190456_208686 DoubleEle_D_190456_208686 DoubleMu_A_190456_208686 DoubleMu_B_190456_208686 DoubleMu_C_190456_208686 DoubleMu_D_190456_208686 EleMu_A_190456_208686 EleMu_B_190456_208686 EleMu_C_190456_208686 EleMu_D_190456_208686"
MC="Tbar_schan Tbar_tchan Tbar_tWchan T_schan TTbarJets T_tchan T_tWchan WGamma WW WZ ZGamma ZJets_2l1 ZJets_2l2 ZJets_2l3 ZJets_2l4 ZZ"
DIR=AnaNaS

mkdir -p /data/users/yhshin/condor/${DIR}
mkdir -p ./condor

for sample in $MC
do
	cat ./analysis-condor.jdl \
	| sed -e s~PREFIX_NAME~/data/users/yhshin/condor/${DIR}~ \
	| sed -e s~RUN_DIR~/home/yhshin/CMSSW_5_3_9_patch3/src/AnaNaS~ \
	| sed -e s~SAMPLE_NAME~${sample}~ \
	> condor/analysis-${sample}.jdl
	echo "condor_submit condor/analysis-${sample}.jdl"
done
