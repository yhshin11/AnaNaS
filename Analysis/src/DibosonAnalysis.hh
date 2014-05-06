#ifndef DibosonAnalysis_h
#define DibosonAnalysis_h

#include <vector>
#include <map>
#include <string>

class TH1;
class TH2;

#include "Analysis/core/SampleAnalysis.hh"

class DibosonAnalysis : public SampleAnalysis
{
public:

  // Diboson modes
  enum Mode { kAll=0,
	      kZZ_2e2n,   kZZ_2m2n,   
	      kZZ_4e,     kZZ_2e2m,     kZZ_4m,
	      kWZ_3en,    kWZ_2emn,   kWZ_e2mn,  kWZ_3mn,
	      kWW_2e2n,   kWW_em2n,   kWW_2m2n, 
	      kZ_2e,      kZ_2m,
              kW_en,      kW_mn,
	      kNModes };
  static std::string modeName[kNModes];

  // level of selection
  enum Level { kNosel=0, kPresel, kSel, kNLevels };
  static std::string levelName[kNLevels];

  //  DibosonAnalysis( Sample& sample, const std::string & collectionFileName );
  DibosonAnalysis( Sample& sample, EventManager& manager );
  virtual ~DibosonAnalysis();
 
public:

  virtual void bookHistograms();
  virtual bool analyzeEvent();
  virtual void writeHistograms();

  std::map< std::string, TH1* > h_;
  std::map< std::string, TH2* > h2_;

  // HLT decision
  virtual bool passHlt();

  // build the Z->ll solutions
  //  imode=0 (Z->ee),  1 (Z->mumu) or 2 (Z-> tautau -> emu)
  enum Zmode { kZeeMode=0, kZmmMode, kZemMode, kNZmode };
  size_t buildZllLists( int imode );
  size_t nZHypos[kNZmode];
  static std::string ZmodeName(          int imode );
  static std::string rejectedZListName(  int imode );
  static std::string ZListName(          int imode, size_t ihyp );
  static std::string leptonName(         int ilepton );
  static std::string leptonListName(     int ilepton );
  static std::string leptonListName(     int ilepton, size_t ihyp );

  // build the list of Zem candidates
  size_t buildZemList();

  // An event hypothesis is a pair of Z->ee and Z->mm solutions
  typedef std::pair< size_t, size_t > Hypo;

  // multimap of event hypos in association with modes
  typedef std::multimap< Hypo, Mode > HypoList;
  typedef HypoList::const_iterator HypoListIterator;

  // HypoList at each level of selection
  HypoList _hypoList[kNLevels]; 

  // occurence of modes at each level of selection
  typedef std::map< Mode, int > ModeStat;
  typedef ModeStat::const_iterator ModeStatIterator;
  ModeStat _nMode[kNLevels]; 

  // occurence of hypotheses at each level of selection
  typedef std::map< Hypo, int > HypoStat;
  typedef HypoStat::const_iterator HypoStatIterator;
  HypoStat _nHypo[kNLevels]; 

  std::map< std::pair< std::string, std::string >, int > _sel;
  std::map< std::string, int > _all;
  std::map< std::string , int > _in;

  // declare selection of a Hypo for a Mode at a given Level
  void select( Hypo, Mode, Level );

  //MM Photon IDIso
  bool isPhotonIDIso(Candidate* ph);

  // Corrected Met (temporary)
  Candidate* correctedMet( HypoListIterator );
  float mT( const Candidate* cand, const Candidate* met );

  // analyses
  //  -- return true if the event hypothesis pointed at by the 
  //                 iterator is selected
  bool   ZZ_4l__analysis( HypoListIterator );
  bool ZZ_2l2n__analysis( HypoListIterator );
  //  bool WZ_3ln__analysis( HypoListIterator );
  bool    Z_2l__analysis( HypoListIterator );
  bool    W_ln__analysis( HypoListIterator );
  bool WW_em2n__analysis( HypoListIterator );

  // stat histograms
  void fillStatHistograms();

  // print-outs
  void print( ostream& o ) const;
  void print( ostream& o, HypoStatIterator, Level ) const;
  void debugPrintouts( ostream& o );

//   //Weighting;
//   float SumWPt;
//   float SumWPtOrig;
//   int Nacc;
//   int Nall;
//   float SumWPtEnd;
//   float SumWPtOrigEnd;
//   int NaccEnd;
//   int NallEnd;
//   float SumWPtEndZZ;
//   float SumWPtOrigEndZZ;
//   int NaccEndZZ;
//   int NallEndZZ;

//   float SumWPtEndE;
//   float SumWPtOrigEndE;
//   int NaccEndE;
//   int NallEndE;
//   float SumWPtEndZZE;
//   float SumWPtOrigEndZZE;
//   int NaccEndZZE;
//   int NallEndZZE;

//   float SumWPtEndM;
//   float SumWPtOrigEndM;
//   int NaccEndM;
//   int NallEndM;
//   float SumWPtEndZZM;
//   float SumWPtOrigEndZZM;
//   int NaccEndZZM;
//   int NallEndZZM;



//   //  void LoadDataBases();
//   vector<vector<float> > DBWeightsZZ;
//   vector<vector<float> > DBWeightsWZ;

  //  float findZWeight(float Zpt, string categ);


  ClassDef( DibosonAnalysis, 0 )  
};

#endif

     
