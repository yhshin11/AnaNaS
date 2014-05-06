#ifndef WAnalysis_h
#define WAnalysis_h

#include <vector>
#include <map>
#include <string>

#include <TString.h>

#include "Analysis/core/SampleAnalysis.hh"
#include "Analysis/selectors/IsolationSelector.hh"
#include "Analysis/utils/KineUtils.hh"
#include "Analysis/tools/CandUtil.hh"

#include <fstream>

class WAnalysis : public SampleAnalysis
{
public:

  WAnalysis( Sample& sample, const std::string & collectionFileName );
  virtual ~WAnalysis();

  static int imode;
 
private:

  //Outputfile
  ofstream file;
  TFile RFile_;
  TTree outTree_;  
  float dPhiJetLep;
  float IsoVariable;
  float MET[2];
  float leptonPt[3];
  void initRootTree();
  void fillTree(float iso ); 

  //Compteurs 
  int Nevt_;
  int Nacceptance_;
  int NHLT_;
  int NPresel_[3];
  int NIsol_[3];
  int NID_[3];
  int NSecondCutE_[3];
  int NConvRej_[3];
  int NGamJetCut_[3];

  //Iso Compteurs
  float NTrkIso_[3];
  float NEcalIso_[3];
  float NHcalIso_[3];

  //Best Electron
  int unsigned BestLepton_;

  //Categories
  std::string EP_;
  int EcalP_;
  std::string IsoCateg_;

  //Elements
  const Candidate* met_;
  const Candidate* theJet_;
  Candidate* theLepton_;

  //Base Matrix
  float LCBM_[2][2];

  //Isolation
  IsolationSelector IsoElec;

  //Herited functions
  virtual void bookHistograms();
  virtual bool analyzeEvent();
  virtual void writeHistograms();


  //Functions
  bool HLTFiring();
  bool EventPreselection();
  bool LeptonIsolation();
  bool LeptonID();
  bool SecondEcut();
  void FillingHistos();
  float PFBaseIso(const Candidate *cand, int n, const std::string & ECAL, float& cIso, float& nIso, float& gIso);
  bool GamJetRejection();
  std::string Wmode;

  // Underlying event
  void FillUnderlyingJets(Candidate* WCand );
  void FillUnderlyingPhoton(Candidate* WCand );

  //MET Commissionning
  void FillSeveralMET(const std::string &);
  TVector2 ComputeRecoil(const Candidate* met);
  void FillRecoil(TVector2 Recoil, Candidate* W, const std::string & type);
  void FillMETProjection(const Candidate* met, Candidate* W, const std::string & type);
  void FillMETVariables(const Candidate* met, const std::string & type);

  TVector2 ComputeRecoilWithCalo(const Candidate* met);
  TVector2 ComputeEtVectorFromRecHit(ecalRecHit rh);
  TVector2 ComputeRecoilWithCaloRecHit(const Candidate* met);
  TVector2 ComputeRecoilWithPF(const Candidate* met);

  //ID
  bool IDByHands();
  //void fillIDObs(int IEp);

  //MC Matching
  void MCMatching(const std::string & P);

  float TrackIsolation();
  float MuonTrackIsolation();
  void Test();
  int FindBestdRCandidate(std::map<int, std::pair<float,float> > AssocMap );
  int NOpposite_;

  ClassDef( WAnalysis, 0 )  
};

#endif

     
