#!/bin/bash

path=/castor/cern.ch/user/m/mmarionn/ReprocessMET/DoubleMu_160404_162000

for i in `rfdir $path | awk '{print $9}'`; do


stager_qry -M $path/$i

done
