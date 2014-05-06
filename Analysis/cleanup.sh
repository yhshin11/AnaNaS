#!/bin/bash
make clean
rm -f *~
rm -f bin/*
rm -f lib/*
list="core utils tools selectors src exe fit"
for i in $list; do
    rm -f $i/*~
    rm -f $i/*_dict*
    rm -f .o/$i/*
    rm -f .d/$i/*
done
echo '> Done.'

