#ifndef JETAnalysis_h
#define JETAnalysis_h

#include <vector>
#include <map>
#include <string>

#include "Analysis/core/SampleAnalysis.hh"

#include "Analysis/selectors/MuonSelector.hh"

class JETAnalysis : public SampleAnalysis
{
public:

  JETAnalysis( Sample& sample, EventManager& manager );
  virtual ~JETAnalysis();

  virtual void bookHistograms();
  virtual bool analyzeEvent();
  virtual void writeHistograms();

private:

  Candidate* BuildBestZCand();

  Vertex* JetVtxAssoc(const Candidate* jet, std::string opt);

  Vertex* mainVtx;
  Candidate* _ZCand;

  MuonSelector MuonSel;

  bool isLepton(const Candidate* jet);
  bool isJetMcMatched(const Candidate* jet);
  bool isMcMatched(const Candidate* pfc);

  void FillJets();

  int Nevt_;
  std::string SampleName;


  ClassDef( JETAnalysis, 0 )  


};



#endif
