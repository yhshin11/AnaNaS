# Source files to compile 
FILES += SelectorBase ConeSelector PtSelector
FILES += ElectronSelector MuonSelector
FILES += LeptonSelector 
FILES += IsolationSelector
FILES += MCMatchingSelector
FILES += RandomConeSelector
FILES += VBTFElectronSelector
FILES += OfficialIsolationSelector

# Header files to use for dictionary generation
DICTFILES := $(FILES) LinkDef

