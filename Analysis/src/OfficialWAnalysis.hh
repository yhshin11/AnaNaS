#ifndef OfficialWAnalysis_h
#define OfficialWAnalysis_h

#include <vector>
#include <map>
#include <string>

#include <TString.h>

#include "Analysis/core/SampleAnalysis.hh"
#include "Analysis/selectors/OfficialIsolationSelector.hh"
#include "Analysis/selectors/VBTFElectronSelector.hh"
#include "Analysis/utils/KineUtils.hh"
#include "Analysis/tools/CandUtil.hh"

#include <stdlib.h>
#include <fstream>

class OfficialWAnalysis : public SampleAnalysis
{
public:

  OfficialWAnalysis( Sample& sample, const std::string & collectionFileName );
  virtual ~OfficialWAnalysis();

  static int imode;
 
private:

  //Protection double comptage
  std::vector<std::pair<int,int> > Events;

  //Compteurs 
  int Nevt_;
  int Nacceptance_;
  int NHLT_;
  int NPresel_[3];
  int NID_[3];
  int NIso_[3];
  int NSecondCutE_[3];
  int NConvRej_[3];

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
  OfficialIsolationSelector IsoElec;
  VBTFElectronSelector vbtfSel;

  //Herited functions
  virtual void bookHistograms();
  virtual bool analyzeEvent();
  virtual void writeHistograms();

  //Name of sample (for deta corrections
  std::string SampleName;

  //Functions
  bool HLTFiring();
  bool EventPreselection();
  bool SpikeCleaning();
  bool LeptonID();
  bool LeptonIso();
  bool SecondEcut();
  bool ConversionRejection();
  void FillingHistos(const std::string & type);
  float PFBaseIso(const Candidate *cand, int n, const std::string & ECAL, float& cIso, float& nIso, float& gIso);
  bool GamJetRejection();
  std::string Wmode;

  void VBTFSelection();

  //outputs
  void InitASCIIS();
  void FillASCII(const std::string & type);
  std::map<std::string, std::ofstream* > ASCIIMap;
  TFile RFile_;
  TTree outTree_; 
  float dPhiJetLep; //ok
  float dPhiMETLep[6]; //ok
  float dPhiMETJet[6]; //ok
  float IsoVar[3]; //ok
  float IDVar[6]; //ok
  //Met and recoil
  float Mets[6][2]; //ok
  float MetProj[6][2]; //ok
  float Recoils[6][2]; //ok
  float RecoilProj[6][2]; //ok
  float dPhiRecoilLep[6]; //ok
  float BasicRecoils[6][2]; //ok
  float BasicRecoilProj[6][2]; //ok
  float BasicdPhiRecoilLep[6]; //ok

  float CTRecoil[6][2];

  float MT[6]; //ok
  float SumET[6];
  float WPt[6];
  
  //Leptons
  float PFLepton[3]; //ok
  float Lepton[3]; 
  float SCLepton[3]; //ok
  float EtSC;

  int charge;
  
  //Selection
  bool VetoE;
  bool ConvRej;
  //bool ConvRej[3];
  bool inHit;
  bool expInHit;

 
  //Underlying event
  bool VetoJet;
  int Nvertex;
  float dZ;
  float z2V;
  int JetMultiplicity;
  int PhotonMultiplicity;
  float LJetPt;
  float LPhotonPt;
  float dPhiJetW[6];
  float dPhiPhotonW[6];
  

  //Event
  int Run;
  int Event;
  int AbsEvent;
  char sample[21];
  char fileName[21];

  float GenRecoil[2][2];
  void ComputeGenRecoil();
  float GenWPt[2];
  float GenEPt[2];
  float GenCharge;
  float GenEnergy;
  
  void initRootTree();
  void PrepareTree(float, float);
  void TreeFillMETProjection(const Candidate* met, int metType);
  void TreeFillRecoil(TVector2 Recoil, int metType);
  void TreeFillBasicRecoil(TVector2 Recoil, int metType);
  void TreeFillCTRecoil(TVector2 Recoil, int metType);
  void TreeFillUnderlyingJetsLepton(const Candidate* met, int metType );
  void TreeFillUnderlyingPhoton(const Candidate* met, int metType );
  void TreeFillLepton();
  bool TreeFillVertices();
  void TreeFillUnderlyingEvent();
  void fillTree(float Iso[3], float ID[6]); 
  
  
  // Underlying event
  void FillUnderlyingJets(Candidate* WCand );
  void FillUnderlyingPhoton(Candidate* WCand );
  bool VetoJetFunction();

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
  TVector2 ComputeRecoilWithTowers(const Candidate* met);

  ClassDef( OfficialWAnalysis, 0 )  
};

#endif

     
