#!/bin/sh
#echo "setup in bash shell"
export Analysis=${PWD}
export LD_LIBRARY_PATH=${ROOTSYS}/lib:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${Analysis}/lib:${LD_LIBRARY_PATH}
export PATH=${Analysis}/bin:${PATH}


