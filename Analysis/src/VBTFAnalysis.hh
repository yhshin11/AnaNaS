#ifndef VBTFAnalysis_h
#define VBTFAnalysis_h

#include <vector>
#include <map>
#include <string>
#include <TNtuple.h>

#include "Analysis/core/SampleAnalysis.hh"
#include "Analysis/core/PseudoExp.hh"

class VBTFAnalysis : public SampleAnalysis, PseudoExp
{
public:

  // Lepton modes (aka imode)
  enum LMode { kElectron=0, kMuon, kNLModes };
  static std::string lmodeName[kNLModes];
  static std::string analysisName[kNLModes];

  // Vector boson modes
  enum Mode { kAll=0,
	      kW_en, kZ_2e,   // imode =0
	      kW_mn, kZ_2m,   // imode =1
	      kNModes };
  static std::string modeName[kNModes];

  // level of selection
  enum Level { kNosel=0, kHLT, kPresel, kSel, kNLevels };
  static std::string levelName[kNLevels];

  VBTFAnalysis( Sample& sample, unsigned imode, const std::string & collectionFileName );
  virtual ~VBTFAnalysis();
 
private:

  virtual void bookHistograms();
  virtual bool analyzeEvent();
  virtual void writeHistograms();

  // build the Z->ll solutions
  //  imode=0 (Z->ee) or 1 (Z->mumu)
  void   buildVBTFLeptonLists();
  size_t buildZllLists();

  const char* leptonListName( const char* level="Presel" ) const;
  const char* ZListName() const;

  bool _pseudo;

  // the lepton mode
  unsigned _imode;

  size_t nZHypos;

  // occurence of modes at each level of selection
  typedef std::map< Mode, int > ModeStat;
  typedef ModeStat::const_iterator ModeStatIterator;
  ModeStat _nMode[kNLevels]; 

  std::map< std::pair< std::string, std::string >, int > _sel;
  std::map< std::string, int > _all;
  std::map< std::string , int > _in;

  // declare selection of a Hypo for a Mode at a given Level
  void select( Mode, Level );

  // analyses
  //  -- return true if the event hypothesis pointed at by the 
  //                 iterator is selected
  bool    Z_2l__analysis();
  bool    W_ln__analysis();

  // stat histograms
  void fillStatHistograms();

  // print-outs
  void print( ostream& o ) const;
  void debugPrintouts( ostream& o );

  // ntuple
  std::vector< std::string > _indexVar;
  std::map< std::string, float > _tupleVar;
  TNtuple* _theTuple;

  ClassDef( VBTFAnalysis, 0 )  
};

#endif

     
