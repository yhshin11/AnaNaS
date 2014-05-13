#!/usr/bin/python

import os,sys,string

jdlfile = """
universe = vanilla
Executable = ./analysis-condor.sh
should_transfer_files = NO
notify_user = 
Output = PREFIX_NAME/$(cluster)_$(process).stdout
Error  = PREFIX_NAME/$(cluster)_$(process).stderr
Log    = PREFIX_NAME/$(cluster)_$(process).condor
Arguments = RUN_DIR SAMPLE_NAME
Queue 1
"""

jdlfile = jdlfile.replace()

    print cmd
    #Execute command
    os.system(cmd)

