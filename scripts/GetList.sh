#!/bin/bash


for i in `ls /tmp/mmarionn/`; do

    ls "/tmp/mmarionn/$i/" > log_$i  

done

ls /tmp/mmarionn/ > log_dataset
