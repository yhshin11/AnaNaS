#!/bin/sh
#echo "setup in bash shell"
export MECore=${PWD}
export LD_LIBRARY_PATH=${ROOTSYS}/lib:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${MECore}/lib:${LD_LIBRARY_PATH}
export PATH=${MECore}/bin:${PATH}


