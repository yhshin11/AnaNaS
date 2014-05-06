#!/bin/sh
#echo "setup in bash shell"
#ghmif [ -z $MUSECAL ]; then
#ghm    echo "Warning MUSECAL is not yet defined"
#ghm    echo "set MUSECAL to default"
#ghm    MUSECAL=../MusEcal
#ghmelse    
#ghm    echo "MUSECAL="$MUSECAL
#ghmfi
export Display=${PWD}
export LD_LIBRARY_PATH=${ROOTSYS}/lib:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${Display}/lib:${LD_LIBRARY_PATH}
export PATH=${Display}/bin:${PATH}


