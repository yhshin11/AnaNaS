#!/bin/bash

n=0

while [ ! $n -eq 100 ] ; do


    rm -r LSF*
    sleep 3m

    n=`echo $n + 1 | bc` 
done
