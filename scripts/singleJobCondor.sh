#!/bin/bash

path=${PWD}
#export STAGE_SVCCLASS=t1transfer
cd $6
source ./setup
cd $path


#stager_get -M /castor/cern.ch/user/m/mmarionn/$1/$2/$3
#rfcp /castor/cern.ch/user/m/mmarionn/$1/$2/$3 .


if [ ! -d $5/$2 ]; then
mkdir $5/$2
fi

### tmp!!! copy in EOS in the same time, permits to make the transfert once
#EOSPATH=/eos/cms/store/user/mmarionn/AnaNaSNtuples
# tmp=`/afs/cern.ch/project/eos/installation/0.2.5/bin/eos.select ls $EOSPATH/$1`
#if [[ -z $tmp ]]; then
#    xrd eoscms mkdir $EOSPATH/$1
#fi
# tmp=`/afs/cern.ch/project/eos/installation/0.2.5/bin/eos.select ls $EOSPATH/$1/$2`
#if [[ -z $tmp ]]; then
#    xrd eoscms mkdir $EOSPATH/$1/$2
#fi
# tmp=`/afs/cern.ch/project/eos/installation/0.2.5/bin/eos.select ls $EOSPATH/$1/$2/$3`
#if [[ -z $tmp ]]; then
#    xrdcp $3 root://eoscms/$EOSPATH/$1/$2/$3
#fi
### ==============================================================================


#cd /afs/cern.ch/user/m/mmarionn/private/AnaNaS/ #commande nécessaire en dépendance de l'initialisation ou pas
#ln -sf $3 /afs/cern.ch/user/m/mmarionn/private/AnaNaS/workdir/data/$2/$3

### this line
analysis -a $4 -f $1/$2/$3 -o $5/$2/AnaTuple_$3

#analysis -a $4 -f $path/$3 -o $path/ #commande nécessaire en dépendance de l'initialisation ou pas

#N=`ls -d $5/$2 2>&1 | grep directory`


#ls -ltr $path/
#ls -ltr $path/$4/

#cp $path/$4/test.root $5/$2/AnaTuple_$3
#rm /afs/cern.ch/user/m/mmarionn/private/AnaNaS/workdir/data/$2/$3
