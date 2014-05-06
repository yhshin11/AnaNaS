#!/bin/bash

VERSION=v18_2011_PhiMod


for i in `rfdir /castor/cern.ch/user/m/mmarionn/Ntuples/$VERSION | awk '{print $9}'`; do

#    if [ "${i:0:8}" != "Z_2m_MCX" ]; then
#	continue
#    fi
	
mkdir /tmp/mmarionn/$i


for j in `rfdir /castor/cern.ch/user/m/mmarionn/Ntuples/$VERSION/$i | awk '{print $9}'`; do

echo $i/$j
if [ ! -e /tmp/mmarionn/$i/$j ]; then
rfcp /castor/cern.ch/user/m/mmarionn/Ntuples/$VERSION/$i/$j /tmp/mmarionn/$i/. &
pid=$!
sleep 3 && tmp=`ps -p $pid | grep $pid`; echo $tmp; if [ -n "$tmp" ]; then  kill $pid; fi &&  echo -n " fail copie " $i/$j
fi

done

done
