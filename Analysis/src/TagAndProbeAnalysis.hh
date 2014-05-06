#ifndef TAPAnalysis_h
#define TAPAnalysis_h

#include <vector>
#include <map>
#include <string>

#include <TString.h>

#include "Analysis/core/SampleAnalysis.hh"
#include "Analysis/selectors/IsolationSelector.hh"
#include "Analysis/selectors/TagAndProbeElectronSelector.hh"
#include "Analysis/utils/KineUtils.hh"
#include "Analysis/tools/CandUtil.hh"

//typedef std::vector< std::string, std::vector<std::string> > CutMap;
//typedef CutMap::const_iterator CutMapIterator;

typedef std::vector<std::vector<int> > NumTP;

class TagAndProbeAnalysis : public SampleAnalysis
{
public:

  TagAndProbeAnalysis( Sample& sample, const std::string & collectionFileName );
  virtual ~TagAndProbeAnalysis();

 
private:
  
  virtual void bookHistograms();
  virtual bool analyzeEvent();
  virtual void writeHistograms();

  //Elements
  //CandList 

  void IsoEfficiency();
  void EndIsolationEfficiency();
  NumTP Ntt_Iso;
  NumTP Ntf_Iso;
  NumTP Ntp_Iso;

  void IDEfficiency();
  void EndIdentificationEfficiency();
  NumTP Ntt_ID;
  NumTP Ntf_ID;
  NumTP Ntp_ID;

  void RecoEfficiency();
  void EndReconstructionEfficiency();
  NumTP Ntt_Reco;
  NumTP Ntf_Reco;
  NumTP Ntp_Reco;


  void AddCondition(CutMap& map, std::string name, std::vector<std::string> values);

  void EfficiencyMeasurement( const CandList& tagList, const CandList& probeList,
			      std::vector<int> _TagList,
			      std::vector<int> _ProbeList, std::string typeCut,
			      std::vector<std::string> valCut,
			      NumTP& ntt, NumTP& ntf, NumTP& ntp);

  
  float ptStep;
  float etaStep;
  int npt;
  int neta;

  int _NVertex;
  int _NJet;

  //Tree variables
  TTree* outTreeID;
  TTree* outTreeIso;
  TTree* outTreeReco;


  float PptID;
  float PetaID;
  float ZMassID;
  int IDPass;
  int IDFail;
  int IDTT; 
  int IDNvtx;
  int IDNJet;

  float PptIso;
  float PetaIso;
  float ZMassIso;
  int IsoPass;
  int IsoFail;
  int IsoTT;
  int IsoNvtx;
  int IsoNJet;

  float PptReco;
  float PetaReco;
  float ZMassReco;
  int RecoPass;
  int RecoFail;
  int RecoTT; 
  int RecoNvtx;
  int RecoNJet;

  void initRootTree();
  void fillTree(std::string type, float mass, bool tf, bool tp, bool tt, float pt, float eta);
  
  void GetNJet();
  void GetNvertex();

  std::string SampleName;

  bool HLTFiring();

  Candidate * BestZCand( const CandList& zllList );
  
  ClassDef( TagAndProbeAnalysis, 0 )  
};

#endif

     
