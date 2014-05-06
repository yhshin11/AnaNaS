#ifndef METAnalysis_h
#define METAnalysis_h

#include <vector>
#include <map>
#include <string>

#include <TRandom3.h>

#include "Analysis/core/SampleAnalysis.hh"

#include "Analysis/selectors/MuonSelector.hh"
#include "Analysis/selectors/ElectronSelector.hh"

class METAnalysis : public SampleAnalysis
{
public:

  METAnalysis( Sample& sample, EventManager& manager );
  virtual ~METAnalysis();

  virtual void bookHistograms();
  virtual bool analyzeEvent();
  virtual void writeHistograms();

private:

  Candidate* BuildBestZCand();
  Candidate* BuildBestZCandMultiChan();
  Candidate* BuildGSFZCand();

  Candidate* MVAJetMET(std::string type,float fact, std::string idLvl);
  Candidate* AtlasMET(float jpt, float fact);
  Candidate* AtlasMETAssoc(std::string type, float fact);
  Candidate* AssocMET(std::string type,float fact);

  Candidate* AtlasMETPt(string type, float jpt, float fact);
  
  Candidate* coreMET(std::string lvlId, float& sumEt);
  Candidate* puMET(std::string lvlId, float& sumEt);


  Vertex* JetVtxAssoc(const Candidate* jet, std::string opt);

  float AtlasFactor();

  void FillMET(string METtype, string redLvl, string factType, float pTVal=20, std::string idLvl="");

  Vertex* mainVtx;
  Candidate* _ZCand;

  MuonSelector MuonSel;
  ElectronSelector ElectronSel;

  bool isLepton(const Candidate* jet);
  bool isJetMcMatched(const Candidate* jet);
  bool isMcMatched(const Candidate* pfc);

  bool HLTFired();
  bool HLTFiredMultiChan();

  void missEDensity();

  //bool ElectronSel(const Candidate* c);

  int Nevt_;
  int _zChan;

  std::string SampleName;


  //Smearing
  void LoadJetSmearingDB();
  Candidate* oldToNewsmearing(const Candidate* smJet);
  Candidate* smearing(const Candidate* jet);
  const Candidate* genJetMatch(const Candidate* jet);
  const Candidate* rawJetMatch(const Candidate* jet);
  const Candidate* smearedJetMatch(const Candidate* jet);
  const Candidate* corJetMatch(const Candidate* jet);
  float smearBy_;
  float shiftBy_;
  TH2D* lut_;
  float skipRawJetPtThreshold_;
  float skipCorrJetPtThreshold_;
  TRandom3 rnd_;


  Candidate* metJetSmearing(Candidate* met);


  Candidate* pfMatch(Candidate* el);
  TLorentzVector pfConeMatch(Candidate* el);
  Candidate* lepTrkMap(Candidate* lep);


  bool _useMu;

  ClassDef( METAnalysis, 0 )  

};



#endif
