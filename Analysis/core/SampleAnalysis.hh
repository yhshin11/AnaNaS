#ifndef SampleAnalysis_h
#define SampleAnalysis_h

#include <string>
#include <vector>
#include <map>

#include "Analysis/utils/HistoManager.hh"
#include "Analysis/utils/TupleManager.hh"
#include "Analysis/core/Sample.hh"
#include "Analysis/core/EventManager.hh"

// Notes on protected methods available in SampleAnalysis
// 
// BuildElectronLists: uses the VBTFElectron selector
//    creates lists of electrons
//      with ET>20GeV within ECAL fiducial [0,1.442]U[1.56,2.5]
//          VeryLoose = WP95 
//          Loose     = WP90 
//          Tight     = WP80
//          VeryTight = WP70
//    all the electrons info blocks are updated with:
//      isEB/B, isEE/B,
//      isEt32/B, isEt30/B, isEt27/B, isEt25/B, isEt20/B,
//      eidXX/B, seieiXX/B, dphiXX/B, detaXX/B, hoeXX/B,  
//      eidAnti/B,  (reversed= dPhiIn and dEtaIn)
//      sumTrk/F, sumEcal/F, sumHcal/F
//      isoXX/B, trkXX/B, ecalXX/B, hcalXX/B,
//      isoAnti/B, (looser isolation)
//      selXX/B, 
//      selAnti/B
//    with XX=60,70,80,85,90,95,   
//    and calls fillElectron histograms for "ref", "VL", "L", "T", "VT"
//
//


//
// event loop on a data sample
//
class SampleAnalysis      
{
public:

  SampleAnalysis( const std::string & name, Sample& sample, 
		  const std::string & collectionFileName ); // deprecated
  SampleAnalysis( const std::string & name, Sample& sample, EventManager& manager );

  virtual ~SampleAnalysis();
 
  void analyze( int nevt=-1, int nskip=-1, bool random=false );

  // virtual functions
  virtual void bookHistograms();
  virtual bool analyzeEvent();
  virtual void writeHistograms();

  virtual void print( ostream& o ) const {}

  // for interactive analysis
  bool nextEvent();
  bool sameEvent();

  void setDebug( int nDebug ) { _nDebug=nDebug; }

  // Z mass windows
  static const float Z0_fullDY[2];
  static const float Z0_largeWindow[2];
  static const float Z0_window[2];
  static const float Z0_veto[2];
  static const float Z0_tightVeto[2];

  // comparators
  struct CandCompare
  {
    bool operator()( Candidate* , Candidate*  ) const;
  };

  struct LeptonCompare
  {
    bool operator()( Candidate* , Candidate*  ) const;
  };

  struct ZllCompare
  {
    bool operator()( Candidate* , Candidate*  ) const;    
  };

protected:

  TString   _name;
  Sample& _sample;
  EventManager& _e;

  bool _verbose;
  int _nDebug;

  // HLT trigger lines
  bool _hlt; // if true, filter on Hlt decision
  std::vector< std::string > _hltLines; // list of Hlt lines
  virtual bool passHlt(); // true if >=1 HLT line is fired

  // get the matched MC candidate for lepton, it exists
  const Candidate* mcCand( const Candidate& cand ) const;

  HistoManager hm;
  void defineTemplate( const std::string & templ, 
		       int nbin, float min, float max );
  void defineTemplate( const std::string & templ, 
		       int nbin1, float min1, float max1,
		       int nbin2, float min2, float max2 );
  void fill( const std::string & templ, const std::string & id, float, const std::string & prefix="" );
  void fill( const std::string & templ, const std::string & id, float, float, const std::string & prefix="" );

  // New GHM (15/05/11)
  TupleManager tm;

  // higer level functions
  void buildLeptonLists();

  void buildElectronLists();
  void fillElectronHistograms( const CandList&, const std::string & );
  void fillElectronHistograms( const Candidate& cand, const std::string &, const std::string & );

  void buildMuonLists();
  void fillMuonHistograms( const CandList&, const std::string & );
  void fillMuonHistograms( const Candidate& cand, const std::string &, const std::string & );

  static bool selectZCand( const Candidate&, float mmin, float mmax, 
			   unsigned nTV=0, unsigned nT=0, unsigned nL=0, unsigned nVL=0 );

  void fillMultHists( const std::string & prefix, const std::string & suffix );

  // event counter
  int  _ievt;

  // output file for histograms
  TFile* _outputFile;

public:
  
  void bookHistograms_();
  void initEvent_();
  bool analyzeEvent_();
  void updateCounters_();
  void writeHistograms_();
  void printStats_();

  bool _trigHlt;
  bool _isSelected;
  std::map< std::string, int > _nAnalyzed;
  std::map< std::string, int > _nHlt;
  std::map< std::string, int > _nSelected;

  // get an histogram from the manager
  TH1* hist( const std::string & templ, const std::string & id, const std::string & dir="", const std::string & prefix="" );
  
  ClassDef( SampleAnalysis, 0 )  
};


#endif

     
