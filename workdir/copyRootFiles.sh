#!/bin/sh
DIR=root/DM
cd $DIR
for i in *; do cp $i/DM/test.root ~/data/DMSkims/$i.root >& copy_log.txt; done
cd ../..

