#ifndef WWCrossCheck_h
#define WWCrossCheck_h

#include <vector>
#include <map>
#include <string>

#include <TString.h>
#include <TObjString.h>
#include <TProfile.h>

#include "Analysis/core/SampleAnalysis.hh"
#include "Analysis/utils/Constants.hh"
#include "Analysis/utils/KineUtils.hh"
#include "Analysis/tools/CandUtil.hh"

#include "TGraph.h"

#include "METNN.hh"

class WWCrossCheck : public SampleAnalysis
{
public:

  WWCrossCheck( Sample& sample, const std::string & collectionFileName  );
  virtual ~WWCrossCheck();

  static int imode;
 
private:

  //For matching
  typedef std::map< Candidate* , Candidate* > Matching;
  typedef Matching::const_iterator MatchingIterator;

  //Protection double comptage
  std::vector<std::pair<int,int> > Events;
  
  //Compteurs 
  int Nevt_;
  int NAcceptance_[2];
  int NevtCand_;
  int NHLT_;
  int NPresel_;
  int NIsol_;
  int NID_;
  int NMassCut_;

  int Nz;
 
  std::string SampleName;

  CandList dau;

  //Categories
  std::string CatEta;

  //Herited functions
  virtual void bookHistograms();
  virtual bool analyzeEvent();
  virtual void writeHistograms();

  TVector2 genRecoil;

  //PFT1 met number
  int _vmet;

  //Functions
  CandList MakeZCandidates();
  CandList MakeZCandidatesWithCharge();
  Candidate* BestZCand(const CandList& zllList );
  bool HLTFiring();
  bool EventPreselection(Candidate* ZCand);
  bool ConversionRejection();
  bool LeptonIso();
  bool LeptonID();
  bool MuonSpecID(Candidate* mu);
  bool MassSelection(Candidate* ZCand);

  bool SpikeCleaning();

  std::string Zmode;

  
  float ZCBM_[2][2];
  float L1CBM_[2][2];
  float L2CBM_[2][2];
  TVector2 ComputeRecoil(const Candidate* met, Candidate* ZCand);
  void ComputeCBMs( Candidate* ZCand );

  std::vector<float> GetRValues(float ZPt, float Recoil[2], int Nv);
  std::vector<float> GetRValuesFromData(float ZPt, float Recoil[2], int Nv);
  float GetResponseForPFT1(float ZPt );

  //MCTruth
  Matching MCToGen;
  Matching GenToMC;
  bool MCok;
  Candidate* FindBestdRCandidate(std::map<Candidate*, std::pair<float,float> > ); 
  int OriginOfMCCandidate(Candidate* cand);
  Candidate* MotherOfMCCandidate(Candidate* cand);
  bool MCTruthChain();
  
  //Neural Net
  IClassifierReader* metNN_;
  void InitializeElNNMET();
  void GetElNNMETValue();

  //Primary vertex
 Vertex* Primary; 

  //outputs
  TTree outTree_; 
  
  // Leptons
  float IDVar[2][4];
  float IsoVar[2][3];
  float Lepton[2][4];
  float SCLepton[2][3];
  int charge[2];
  int pdgId[2];
  float EOP[2];
  bool inHit1;
  bool expInHit1;
  bool ConvRej1;
  bool inHit2;
  bool expInHit2;
  bool ConvRej2;

  int ElecMult;
  int MuonMult;
  float LeadAddElec[3];
  float LeadAddMuon[3];
  bool IsIsoMu;
  bool IsIDMu;
  bool IsIsoElec;
  bool IsIDElec;
  
  float ElecMT;
  float MuMT;

  //MCTruth
  float McLepton[2][5];
  float McNeutrino[2][5];
  float McZn[4];
  float McZl[4];

  //Z
  float Z[4];
  float ZMass;
  float ZMT;
  int ZCateg;
  int ZMultiplicity;
  float PfZ[4];
  float PfZMass;

  //MET
  float MET[2];
  float Recoil[2];

  float SpecMET;
  float SumEt;

  float TrkMET[2];

  float METPrime[2];
  float METPrimeRecoil[2];
  float METPrimeCor[2];
 
  float METPrimeMETM[2];
  float METPrimeJetM[2];
  float METPrimeJetP[2];
  float METPrimeUnclusM[2];
  float METPrimeUnclusP[2];


  float ProjMET;
  float CorProjMET;

  //MET NeuralNet
  float outMETNN;

  //Resolution tails
  float RValMC[2];
  float RValData[2];

  //Jet
  int JetMultiplicity;
  int MultiJetMult[10];
  int TcJetMultiplicity;
  float LJet[4];
  float bTag[2];

  //Photon
  int PhotonMultiplicity;
  float LPhoton[5];

  //Lepton related variables
  float dPhiLepton;
  float dPhiLeptonZ;
  float dRLepton;
  float dEtaLepton;
  float dPtLepton;
  float dPLepton;

  float dPhiLeptonS;
  float dRLeptonS;
  float dEtaLeptonS;
  float dPtLeptonS;
  float dPLeptonS;
  float CosThetaM;
  float CosThetaP;

  float CosThetaCSM;
  float CosThetaCSP;

  //MET related variables
  float ZMETMass;
  float ZMassMET;
  float ZPtMET;
  float ZPMET;
  float dPhiMETZ;
  float dPhiMETL1;
  float dPhiMETL2;
  float METProjZ[2];
  float METProjL1[2];
  float METProjL2[2];
  float ZRecoil;
  float RecoilProjZ[2];
  float RecoilProjL1[2];
  float RecoilProjL2[2];
  float dPhiRecoilZ;
  float dPhiRecoilL1;
  float dPhiRecoilL2;

  float HT[2];

  //Jet related variables
  float dPhiJetMET;
  float dPhiJetZ;
  float dEtaJetZ;
  float ZPtJet;
  float ZJetMass;
  float dPhiJetL1;
  float dPhiJetL2;
  float dEtaJetL1;
  float dEtaJetL2;
  float dPhiPhotonMET;
  float dPhiPhotonZ;
  float dEtaPhotonZ;
  float ZPtPhoton;
  float ZPhotonMass;
  float dPhiPhotonL1;
  float dPhiPhotonL2;
  float dEtaPhotonL1;
  float dEtaPhotonL2;
  
  //Other variables
  int Nvertex;
  
  //Miscellaneous
  int Run ;
  int Event;
  int AbsEvent;
  char sample[50];
  char fileName[50];
  char EvtCateg[50];
  
  void initRootTree();
  void fillTree(Candidate* ZCand); 
  
  void fillLeptonTree(Candidate* ZCand); 
  void fillZTree(Candidate* ZCand); 
  void fillMETTree(Candidate* ZCand); 
  void fillJetTree(Candidate* ZCand); 
  void fillPhotonTree(Candidate* ZCand); 
  void fillEventTree(Candidate* ZCand); 
  void ZFrameVariables(Candidate* ZCand);
  void FillGenEvent();
  void fillMETPrime(Candidate* ZCand);
  void fillProjectedMET();
  bool fillVertices();

  float ComputeResolution(float E);

  ClassDef( WWCrossCheck, 0 )  
};

#endif

     
