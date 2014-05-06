# Source files to compile 
FILES := ParticleData Config Constants RooUtils RooFitUtils ObjectStore KineUtils StatUtils CfgParser FileUtils StringUtils

# Header files to use for dictionary generation
DICTFILES := $(FILES) LinkDef

# Executable files
PROGRAMS := sacroot test_TTree test_Tuple test_Scan test_Cuts test_Dataset

