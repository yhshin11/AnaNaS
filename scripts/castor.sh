#!/bin/bash

path=${PWD}

cd /afs/cern.ch/user/m/mmarionn/private/AnaNaS/
source rootscript
cd $path


stager_get -M /castor/cern.ch/user/m/mmarionn/$1/$2/$3
rfcp /castor/cern.ch/user/m/mmarionn/$1/$2/$3 .
#cd /afs/cern.ch/user/m/mmarionn/private/AnaNaS/ #commande nécessaire en dépendance de l'initialisation ou pas
#ln -sf $3 /afs/cern.ch/user/m/mmarionn/private/AnaNaS/workdir/data/$2/$3
analysis -a $4 -f $3 -o $path/
#analysis -a $4 -f $path/$3 -o $path/ #commande nécessaire en dépendance de l'initialisation ou pas

N=`rfdir /castor/cern.ch/user/m/mmarionn/Ntuples/v9_2011/$2 2>&1 | grep directory`

if [ -z "$N" ]; then
rfmkdir /castor/cern.ch/user/m/mmarionn/Ntuples/v9_2011/$2
fi


rfcp $path/$4/$3.root /castor/cern.ch/user/m/mmarionn/Ntuples/$5/$2/AnaTuple_$3
#rm /afs/cern.ch/user/m/mmarionn/private/AnaNaS/workdir/data/$2/$3
