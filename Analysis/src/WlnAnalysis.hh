#ifndef WlnAnalysis_h
#define WlnAnalysis_h

#include <vector>
#include <map>
#include <string>

class TH1;
class TH2;

#include "Analysis/core/SampleAnalysis.hh"
#include "Analysis/core/PseudoExp.hh"

class WlnAnalysis : public SampleAnalysis, PseudoExp
{
public:

  // Diboson modes
  enum Mode { kAll=0,
              kW_en,      kW_mn,
	      kNModes };
  static std::string modeName[kNModes];

  // level of selection
  enum Level { kNosel=0, kPresel, kSel, kNLevels };
  static std::string levelName[kNLevels];

  WlnAnalysis( Sample& sample, const std::string & collectionFileName );
  virtual ~WlnAnalysis();
 
private:

  virtual void bookHistograms();
  virtual bool analyzeEvent();
  virtual void writeHistograms();

  std::map< std::string, TH1* > h_;
  std::map< std::string, TH2* > h2_;

  // occurence of modes at each level of selection
  typedef std::map< Mode, int > ModeStat;
  typedef ModeStat::const_iterator ModeStatIterator;
  ModeStat _nMode[kNLevels]; 

  std::map< std::pair< std::string, std::string >, int > _sel;
  std::map< std::string, int > _all;
  std::map< std::string , int > _in;

  // Corrected Met (temporary)
  float mT( const Candidate* cand, const Candidate* met );

  // analyses
  //  -- return true if the event hypothesis pointed at by the 
  //                 iterator is selected
  bool    W_en__analysis();
  bool    W_mn__analysis();

  // stat histograms
  void fillStatHistograms();

  // print-outs
  void print( ostream& o ) const;
  void debugPrintouts( ostream& o );

  ClassDef( WlnAnalysis, 0 )  
};

#endif

     
