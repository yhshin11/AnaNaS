#ifndef DMAnalysis_h
#define DMAnalysis_h

#include <vector>
#include <map>
#include <string>

class TH1;
class TH2;

#include "Analysis/core/SampleAnalysis.hh"

#include "Analysis/selectors/MuonSelector.hh"
#include "Analysis/selectors/ElectronSelector.hh"

class DMAnalysis : public SampleAnalysis
{
public:

  //  DMAnalysis( Sample& sample, const std::string & collectionFileName );
  DMAnalysis( Sample& sample, EventManager& manager );
  virtual ~DMAnalysis();
 
public:

  virtual void bookHistograms();
  virtual bool analyzeEvent();
  virtual void writeHistograms();

private:
  // insert private members here
  //XXX TH1* hist_mll;
  //XXX double mll_tree;
  float getRapidity(const Candidate* zCandidate);

  MuonSelector MuonSel;
  ElectronSelector ElectronSel;

  bool HLTFiredMultiChan();

  int _zChan; // 0 for Z->mumu; 1 for Z->ee; 2 for Z->mue

  ClassDef( DMAnalysis, 0 )  
};

#endif

     
