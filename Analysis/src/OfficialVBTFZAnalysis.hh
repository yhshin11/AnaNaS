#ifndef OfficialVBTFZAnalysis_h
#define OfficialVBTFZAnalysis_h

#include <vector>
#include <map>
#include <string>
#include <fstream>

#include <TString.h>
#include <TProfile.h>

#include "Analysis/core/SampleAnalysis.hh"
#include "Analysis/selectors/OfficialIsolationSelector.hh"
#include "Analysis/selectors/VBTFElectronSelector.hh"
#include "Analysis/utils/KineUtils.hh"
#include "Analysis/tools/CandUtil.hh"

#include "TGraph.h"

class OfficialVBTFZAnalysis : public SampleAnalysis
{
public:

  OfficialVBTFZAnalysis( Sample& sample, const std::string & collectionFileName  );
  virtual ~OfficialVBTFZAnalysis();

  static int imode;
 
private:

  //Map Leptons
  void LoadElecDB();
  std::map< unsigned int , std::vector<float> > mapLepton;
  std::map< unsigned int , std::vector<float> >::const_iterator itermapLep; 

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

  int Nz;
 
  std::string SampleName;

  //Elements
  //  const CandList& electrons_;
  const Candidate* met_;
  const Candidate* theJet_;
   
  CandList dau;

  //Categories
  std::string CatEta;

  //Isolation
  OfficialIsolationSelector IsoElec;
  VBTFElectronSelector vbtfSel;

  //Herited functions
  virtual void bookHistograms();
  virtual bool analyzeEvent();
  virtual void writeHistograms();

  TVector2 genRecoil;

  //Functions
  CandList MakeZCandidates();
  CandList MakeZCandidatesWithCharge();
  CandList MakeZCandidatesUsingLeptonDB();
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

//outputs
  void InitASCII();
  void FillASCII();
  std::ofstream* ASCIIFile;
  TFile RFile_;
  TTree outTree_; 
  
  float dPhiMETZ[5]; //ok
  float dPhiMETJet[5]; //ok
  float IsoVar[2][3]; //ok
  float IDVar[2][6]; //ok
  //Met and recoil
  float Mets[5][2]; //ok
  float MetProj[5][2]; //ok
  float Recoils[5][2]; //ok
  float RecoilProj[5][2]; //ok
  float dPhiRecoilZ[5]; //ok
  float BasicRecoils[5][2]; //ok
  float BasicRecoilProj[5][2]; //ok
  float BasicdPhiRecoilZ[5]; //ok

  float Corrections[5];

  float CTRecoil[5][2];

  float ZPt;
  float ZMass;
  
  //Leptons
  float PFLepton[2][3]; //ok
  float Lepton[2][3]; 
  float SCLepton[2][4]; //ok
  float EtSC[2];
  int iSeed[2][3];
  int charge[2];
  
  //Quality variables
  float FBrem[2];
  float EOP[2];
  float ESeedOPout[2];
  float pIn[2];

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
  float MetCor[5][2];
  float MetCorProj[5][2];
  float RecoilCor[5][2];
  float RecoilCorProj[5][2];
  float SumET[5];
  
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
  int time;
  char sample[21];
  char fileName[21];


  float GenRecoil[2][2];
  void ComputeGenRecoil();
  float GenZPt;
  float GenZMass;
  float GenLepton[2][3];

  void initRootTree();
  void PrepareTree(float[2], float[2], Candidate* ZCand);
  void TreeFillMETProjection(const Candidate* met, int metType, Candidate* ZCand);
  void TreeFillRecoil(TVector2 Recoil, int metType, Candidate* ZCand);
  void TreeFillBasicRecoil(TVector2 Recoil, int metType, Candidate* ZCand);
  void TreeFillCTRecoil(TVector2 Recoil, int metType, Candidate* ZCand);
  void TreeFillUnderlyingJetsZ(const Candidate* met, int metType, Candidate* ZCand );
  void TreeFillUnderlyingPhoton(const Candidate* met, int metType, Candidate* ZCand );
  void TreeFillLepton();
  void TreeFillVertices();
  void TreeFillUnderlyingEvent();
  void fillTree(float Iso[2][3], float ID[2][6], Candidate* ZCand); 
  void SeparateMETElements();
  void SpecPfMETVariables();
  
  // Underlying event
  bool VetoJetFunction();
  
  ClassDef( OfficialVBTFZAnalysis, 0 )  
};

#endif

     
