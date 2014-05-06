#!/bin/bash

BASE=/hadoop/store/user/mmarionn/MET
VERSION=Prod8TeV_v1  #Prod8TeV_v4   #Prod8TeV_53X_v1
ListDS=DatasetFile

#ls $ListDS

ANTYPE=MET #MET

TMPOPATH=/data/users/mmarionn/MET
OVERSION=ANSNtuplesv1

OPATH=$TMPOPATH/$OVERSION

ANPATH=$ANANAS
N=0
nohup ./CleanLSF.sh &

if [ ! -d $OPATH ]; then
    mkdir -p $OPATH
fi


while read line; do


echo $line
  #  REP=` rfdir /castor/cern.ch/user/m/mmarionn/Ntuples/$OPATH/$DATASET`
    
    DATASET=$line
    
#FIXME
 #   if [ "$DATASET" == "TTbarJets" ]; then
#	VERSION=Prod8TeV_v2
#    fi


    if [ "${DATASET:0:2}" == "//" ]; then
	echo "skipping $DATASET"
	continue
    fi


    if [ -z "$OPATH/$line" ]; then
	mkdir $OPATH/$line
    fi

    echo $BASE/$VERSION/$DATASET/

      for j in `ls $BASE/$VERSION/$DATASET/`; do

    	FILE=$j
#	MONTH=`ls $BASE/$VERSION/$DATASET/$j | awk '{print $6}'`
	#DAY=`ls $BASE/$VERSION/$DATASET/$j }'`
	
	OFILE=$OPATH/$DATASET/AnaTuple_$FILE

	#echo $OFILE

	#protection against existing files
	if [ -e $OFILE ]; then
	  echo "file existing : $OFILE" 
	  continue
	fi

#	echo $FILE
	
	#echo $VERSION $DATASET $FILE $ANTYPE $OPATH $ANPATH
	#bsub -R "pool>10000" -q 8nh source singleJob.sh $BASE/$VERSION $DATASET $FILE $ANTYPE $OPATH $ANPATH
	
# Condor submit file                                                                                                                                            
cat > tmpFiles/subfile_$N <<EOF

universe = vanilla
Executable = /home/mmarionn/AnaNaS/scripts/singleJobCondor.sh
getenv = True
Should_Transfer_Files = NO
Output = /home/mmarionn/AnaNaS/scripts/logsProd/logOut_$N.log
Error = /home/mmarionn/AnaNaS/scripts/logsProd/logErr_$N.log
Log = /home/mmarionn/AnaNaS/scripts/logsProd/log_$N.log
notify_user = matthieu.marionneau@cern.ch
Arguments = $BASE/$VERSION $DATASET $FILE $ANTYPE $OPATH $ANPATH
Queue 1

EOF

# submit the job
condor_submit tmpFiles/subfile_$N

	N=`echo $N + 1 | bc`

	if [[ $N == 80 ]]; then
	    N=0
	    sleep 10m
	fi

	break;

    done


done < $ListDS

