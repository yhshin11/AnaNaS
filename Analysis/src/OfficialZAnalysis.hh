#ifndef _OfficialZAnalysis_h
#define _OfficialZAnalysis_h

#include <vector>
#include <map>
#include <string>
#include <stdlib.h>

#include <TString.h>
#include <TProfile.h>
#include <TGraph.h>

#include "Analysis/core/SampleAnalysis.hh"
#include "Analysis/core/EventManager.hh"
#include "Analysis/core/CandInfo.hh"

#include "Analysis/utils/KineUtils.hh"
#include "Analysis/tools/CandUtil.hh"
#include "Analysis/utils/Constants.hh"
//#include "Analysis/utils/Config.hh"
#include "Analysis/selectors/IsolationSelector.hh"
//#include "Analysis/selectors/ConeSelector.hh"
#include "Analysis/utils/StatUtils.hh"




class OfficialZAnalysis : public SampleAnalysis
{
public:

  OfficialZAnalysis( Sample& sample, const std::string & collectionFileName  );
  virtual ~OfficialZAnalysis();

  static int imode;

private:

  //Protection double comptage
  std::vector<std::pair<int,int> > Events;
  
  //Compteurs 
  int Nevt_;
  int NAcceptance_;
  int NevtCand_;
  int NHLT_;
  int NPresel_;
  int NIsol_;
  int NID_;
  int NMassCut_;

  int NZ_;
 
  std::string SampleName;

  //Profiltes
  TProfile* PFUparaQpara;
  TProfile* TcUparaQpara;
  TProfile* CaloUparaQpara;
  TProfile* RawCaloUparaQpara;

  std::vector<TGraph*> dPhiRecZ;
  std::vector<TGraph*> RecVsPt;
  std::vector<TGraph*> RecpVsPt;
  std::vector<TGraph*> UQTvsQT;
  std::vector<TGraph*> UpQTvsQT;
  std::vector<TGraph*> RecpCorVsPt;
  std::vector<TGraph*> RecpeCorVsPt;
  std::vector<TGraph*> RecpCorVsGenPt;
  std::vector<TGraph*> RecpeCorVsGenPt;
  //Elements
  //  const CandList& electrons_;
  const Candidate* met_;
  const Candidate* theJet_;
   
  CandList dau;

  //Categories
  std::string CatEta;

  //Isolation
  //  OfficialIsolationSelector IsoElec;
  //  VBTFElectronSelector vbtfSel;

  //Herited functions
  virtual void bookHistograms();
  virtual bool analyzeEvent();
  virtual void writeHistograms();

  TVector2 genRecoil;

  //Functions
  CandList MakeZCandidates();
  CandList MakeZCandidatesWithCharge();
  Candidate* BestZCand(const CandList& zllList );
  bool HLTFiring();
  bool EventPreselection(Candidate* ZCand);
  bool ConversionRejection();
  bool LeptonIso();
  bool LeptonID();
  bool MassSelection(Candidate* ZCand);

  bool SpikeCleaning();

  void FillingHistos(Candidate* ZCand, const std::string &);
  void FillUnderlyingPhoton(Candidate* ZCand );
  void FillUnderlyingJets(Candidate* ZCand );
  void FillUnderlyingMET(Candidate* ZCand );

  std::string Zmode;

  //Commissioning MET
  float ZCBM_[2][2];
  TVector2 ComputeRecoil(const Candidate* met, Candidate* ZCand,double cor);
  void FillRecoil(TVector2 Recoil, Candidate* ZCand, const std::string & type);
  void FillMETProjection(const Candidate* met, Candidate* ZCand, const std::string & type);
  void FillMETVariables(const Candidate* met, Candidate* ZCand, const std::string & type);
  void FillSeveralMET( Candidate* ZCand, const std::string & ALevel);
 
  TVector2 ComputeRecoilWithCalo(const Candidate* met, Candidate* ZCand,double cor);
  TVector2 ComputeRecoilWithPF(const Candidate* met, Candidate* ZCand,double cor);
  TVector2 ComputeEtVectorFromRecHit(ecalRecHit rh);
  TVector2 ComputeRecoilWithCaloRecHit(const Candidate* met, Candidate* ZCand);
  TVector2 ComputeRecoilWithTowers(const Candidate* met, Candidate* ZCand,double cor);
  TVector2 ComputeGenRecoil(const Candidate* met);

  void FillMETWithCorrections(const Candidate* met,int metType, Candidate* ZCand);

  std::vector<Candidate*> RemoveLaserCor(Candidate* ZCand);

  double electron_br1EB(double brLinear, double e);
  double electron_br1EE(double brLinear, double e);
  double electron_br1_completeEB(double brLinear, double e);
  double electron_br1_completeEE(double brLinear, double e);
  double f5x5(double iEta);

//outputs
  void InitASCIIS();
  void FillASCII(const std::string & type);
  std::map<std::string, std::ofstream* > ASCIIMap;
  //TFile RFile_;
  TTree* outTree_; 

  float dPhiMETZ[6]; //ok
  float dPhiMETJet[6]; //ok
  float IsoVar[2][3]; //ok
  float IDVar[2][6]; //ok
  //Met and recoil
  float Mets[6][2]; //ok
  float MetProj[6][2]; //ok
  float Recoils[6][2]; //ok
  float RecoilProj[6][2]; //ok
  float dPhiRecoilZ[6]; //ok
  float BasicRecoils[6][2]; //ok
  float BasicRecoilProj[6][2]; //ok
  float BisectorRecoilProj[6][2]; //ok
  float BasicdPhiRecoilZ[6]; //ok

  float Corrections[6];

  float CTRecoil[6][2];

  float ZPt;
  float ZMass;

  //Spec laser Mass

  float ZMassSC;
  float ZPtSC;
  float ZMassSCUnCor;
  float ZPtSCUnCor;
  float SCCor[2][3];
  float UnCorSC[2][3];
  
  //Leptons
  float PFLepton[2][3]; //ok
  float Lepton[2][3]; 
  float SCLepton[2][3]; //ok
  float EtSC[2];

  int charge[2];
  
  //Selection
  bool VetoE;
  bool inHit1;
  bool expInHit1;
  bool ConvRej1;
  //bool ConvRej[3];
  bool inHit2;
  bool expInHit2;
  bool ConvRej2;
 
  //Underlying event
  bool VetoJet;
  int Nvertex;
  float dZ;
  int JetMultiplicity;
  int PhotonMultiplicity;
  float LJetPt;
  float LPhotonPt;
  float dPhiJetZ;
  float dPhiPhotonZ;
  
  int NPhot1;
  int NPhot2;
  
  //MET corrected
  float MetCor[6][2];
  float MetCorProj[6][2];
  float RecoilCor[6][2];
  float RecoilCorProj[6][2];
  float SumET[6];
  
  float cEmFrac;
  float nEmFrac;
  float cHadFrac;
  float nHadFrac;
  float cMuFrac;

  //Decomposed MET
  float CaloSumETEB;
  float CaloSumETEE;
  float CaloSumETHB;
  float CaloSumETHE;
  float CaloSumETHF;
  
  float CaloMETEB[2];
  float CaloMETEE[2];
  float CaloMETHB[2];
  float CaloMETHE[2];
  float CaloMETHF[2];
  
  float PfSumETEB;
  float PfSumETEE;
  float PfSumETHB;
  float PfSumETHE;
  float PfSumETHF;
  
  float PfMETEB[2];
  float PfMETEE[2];
  float PfMETEEx[2];
  float PfMETHB[2];
  float PfMETHE[2];
  float PfMETHF[2];

  float PfEcalMET[2];
  float PfHcalMET[2];
  float CaloEcalMET[2];
  float CaloHcalMET[2];

  //Event
  int Run;
  int Event;
  int AbsEvent;
  char sampleName[21];
  char fileName[21];
  int _time;

  
  float GenRecoil[2][2];
  void ComputeGenRecoil(Candidate* ZCand);
  float GenZPt;
  float dPhiGenZ;

  void initRootTree();
  void PrepareTree( Candidate* ZCand);
  void TreeFillMETProjection(const Candidate* met, int metType, Candidate* ZCand);
  void TreeFillRecoil(TVector2 Recoil, int metType, Candidate* ZCand);
  void TreeFillBasicRecoil(TVector2 Recoil, int metType, Candidate* ZCand);
  void TreeFillBisectorRecoil(TVector2 Recoil, int metType, Candidate* ZCand);
  void TreeFillCTRecoil(TVector2 Recoil, int metType, Candidate* ZCand);
  void TreeFillUnderlyingJetsZ(const Candidate* met, int metType, Candidate* ZCand );
  void TreeFillUnderlyingPhoton(const Candidate* met, int metType, Candidate* ZCand );
  void TreeFillLepton();
  bool TreeFillVertices();
  void TreeFillUnderlyingEvent();
  void fillTree(float Iso[2][3], float ID[2][6], Candidate* ZCand); 
  void SeparateMETElements();
  void SpecPfMETVariables();
  
  // Underlying event
  bool VetoJetFunction();
  
  ClassDef( OfficialZAnalysis, 0 )  
};

#endif

     
