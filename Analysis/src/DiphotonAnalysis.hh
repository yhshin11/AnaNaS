#ifndef DiphotonAnalysis_h
#define DiphotonAnalysis_h

#include "Analysis/core/SampleAnalysis.hh"
#include "Analysis/core/CandInfo.hh"
#include "Analysis/core/EventManager.hh"
#include "Analysis/utils/KineUtils.hh"
#include "Analysis/utils/Constants.hh"
#include "Analysis/utils/StatUtils.hh"
#include "Analysis/utils/KineUtils.hh"
#include "Analysis/selectors/ConeSelector.hh"
#include "Analysis/selectors/MCMatchingSelector.hh"
#include "Analysis/selectors/RandomConeSelector.hh"
#include "Analysis/selectors/ConeSelector.hh"
#include "Analysis/tools/CandUtil.hh"

#include <vector>
#include <map>
#include <string>
#include <TString.h>
#include <TH1I.h>
#include <TH1F.h>
#include <TH2I.h>
#include <TH2F.h>
#include <TRandom3.h>
#include <stdlib.h>
#include <fstream>
#include <time.h>
#include <Math/VectorUtil.h>
#include <Math/Vector4D.h>

class DiphotonAnalysis : public SampleAnalysis
{
public:

  DiphotonAnalysis( Sample& sample, const std::string & collectionFileName );
  virtual ~DiphotonAnalysis();

  static int imode;
 
private:

  //Protection double comptage
  std::vector<std::pair<int,int> > Events;
  bool debug;
  bool doDiphotons;
  bool doShuffle;
  bool doFSel;
  int addSigCut; // 0 doNothing 1 cut<.011 2 flip>.011
  int noPixSeed; // 0 doNothing 1 vetoPix 2 requirePix 

  //Compteurs 
  int Nevt_;
  int Nacceptance_;
  int NHLT_;
  int Npresel_;
  
  //Preselection
  float pho1_minEt;
  float pho2_minEt;
  float pho_maxHoe;

  float barrel_maxEta;
  float endcap_minEta;
  float endcap_maxEta;

  float pho_maxHcal;
  float pho_maxTrack;
  float pho_maxSietaEB;
  float pho_maxSietaEE;
  float pho_maxSieta;

  //Elements
  Candidate* tmpPhoton1_; Candidate* tmpPhoton2_;
  Candidate* thePhoton1_; Candidate* thePhoton2_;
  bool hasPhoton1;
  bool hasPhoton2;
  Candidate* theMCpho1_;  Candidate* theMCpho2_;
  Candidate* theJet1_;    Candidate* theJet2_;
  
  const Candidate* hardestJet_; // highest pt jet
  
  //Inerited functions
  virtual bool analyzeEvent();
  virtual void writeHistograms();
  virtual void bookHistograms();

  //Name of sample (for deta corrections)
  std::string SampleName;

  //Functions
  bool  HLTFiring();
  bool  EventPreselection();
  bool  GoodVertex();
  bool  SpikeCleaning(Candidate* cand);
  void  initRootTree();
  void  TreeFillJet();
  void  TreeFillPhoton();
  void  TreeFillVertices();
  void  TreeFillUnderlyingEvent();
  void  fillTree(); 
  void  MCtruthMatching(); 
  void  JetMatching(); 
  float getGenIso(Candidate*,bool);
  bool  inCrack(int, int, int);
  void  makeRHvector();
  bool  jetInCone(Candidate*);
  float cosThetaStar(TLorentzVector, TLorentzVector);


  // to loop only once, keep a vector of the indices of the interesting rechits
  vector<int>   rhVector;
  vector<float> keepEta, keepPhi;

  //outputs
  void InitASCIIS();
  void FillASCII(std::string type);
  std::map<std::string, std::ofstream* > ASCIIMap;
  TTree outTree_;
  float dPhiJetLep;
  float IsoVar[3];
  float IDVar[6];

  // random
  TRandom3 random;

  //Selection
  bool VetoE;
  bool ConvRej;
  bool inHit;
  bool expInHit;
 
  // Event
  int Run;
  int Event;
  int LumiSec;
  int nPhotons;
  int nVertex;
  int nGoodVertex;
  int nJets;
  int nTracks;
  int nClusters;
  int nPixSeeds;
  int nSC;
  int nRecHits;
  float met;
  float met_phi;
  float mtrans;
  float ptHat;
  bool  isShuffled;
    
  // trigger bits
  bool HLT_Photon10_Cleaned_L1R;
  bool HLT_Photon17_Isol_SC17HE_L1R_v1;
  bool HLT_Photon17_SC17HE_L1R_v1;
  bool HLT_Photon20_Isol_Cleaned_L1R_v1;
  bool HLT_Photon25_Cleaned_L1R;
  bool HLT_Photon35_Isol_Cleaned_L1R_v1;
  bool HLT_Photon40_Isol_Cleaned_L1R_v1;
  bool HLT_Photon40_CaloId_Cleaned_L1R_v1;
  bool HLT_Photon15_Cleaned_L1R;
  bool HLT_Photon15_L1R;
  bool HLT_Photon22_SC22HE_L1R;
  bool HLT_Photon20_Cleaned_L1R;
  bool HLT_Photon30_Cleaned_L1R;
  bool HLT_Photon50_Cleaned_L1R;
  bool HLT_Photon70_Cleaned_L1R;
  bool HLT_DoublePhoton22_L1R_v1;
  bool HLT_DoublePhoton20_L1R;
  bool HLT_DoublePhoton15_L1R;
  bool HLT_DoublePhoton5_L1R;
  bool HLT_DoublePhoton17_L1R;
  bool HLT_DoublePhoton17_SingleIsol_L1R_v1;
  bool HLT_DoublePhoton10_L1R;

  // Primary vertex
  float vX;
  float vY;
  float vZ;
  int   vNtracks;
  int   vNout;
  bool  vIsGood;
  
  //Photon1&2
  float p1;        float p2;
  float px1;       float px2;
  float py1;       float py2;
  float pz1;       float pz2;
  float et1;       float et2;
  float energy1;   float energy2;
  float eta1;      float eta2;
  float phi1;      float phi2;
  bool isEBGap1;   bool isEBGap2;
  bool isEEGap1;   bool isEEGap2;
  bool isEBEEGap1; bool isEBEEGap2;
  bool isEB1;      bool isEB2;
  bool isEE1;      bool isEE2;
  float sigmaEtaEta1;             float sigmaEtaEta2;
  float sigmaIetaIeta1;           float sigmaIetaIeta2;
  float e1x51;                    float e1x52;
  float e2x51;                    float e2x52;
  float ecalRecHitSumEtConeDR031; float ecalRecHitSumEtConeDR032;
  float ecalRecHitSumEtConeDR041; float ecalRecHitSumEtConeDR042;
  float hadronicDepth1OverEm1;    float hadronicDepth1OverEm2;
  float hadronicDepth2OverEm1;    float hadronicDepth2OverEm2;
  float hadronicOverEm1;          float hadronicOverEm2;
  float hcalDepth1TowerSumEtConeDR031; float hcalDepth1TowerSumEtConeDR032;
  float hcalDepth1TowerSumEtConeDR041; float hcalDepth1TowerSumEtConeDR042;
  float hcalDepth2TowerSumEtConeDR031; float hcalDepth2TowerSumEtConeDR032;
  float hcalDepth2TowerSumEtConeDR041; float hcalDepth2TowerSumEtConeDR042;
  float hcalTowerSumEtConeDR031; float hcalTowerSumEtConeDR032;
  float hcalTowerSumEtConeDR041; float hcalTowerSumEtConeDR042;
  float nTrkHollowConeDR031;     float nTrkHollowConeDR032;
  float nTrkHollowConeDR041;     float nTrkHollowConeDR042;
  float nTrkSolidConeDR031;      float nTrkSolidConeDR032;
  float nTrkSolidConeDR041;      float nTrkSolidConeDR042;
  float trkSumPtHollowConeDR031; float trkSumPtHollowConeDR032;
  float trkSumPtHollowConeDR041; float trkSumPtHollowConeDR042;
  float trkSumPtSolidConeDR031;  float trkSumPtSolidConeDR032;
  float trkSumPtSolidConeDR041;  float trkSumPtSolidConeDR042;
  float patSumIso1;  float patSumIso2;
  float hcalIso1;    float hcalIso2;
  int   hcalIso1_n,  hcalIso2_n;  
  float hcalIsoH1,   hcalIsoH2;
  int   hcalIsoH1_n, hcalIsoH2_n;
  float trackIso1;   float trackIso2;
  float trackIsoRd1; float trackIsoRd2;
  float trackIsoH1;  float trackIsoH2;
  int trackIso1_n;
  int trackIsoH1_n;
  float sumIso1;    float sumIso2;

  float scPhi1;      float scPhi2;
  float scEta1;      float scEta2;
  float scRawE1;     float scRawE2;
  int   nXtalsSC1;   int   nXtalsSC2;
  float scPhiWidth1; float scPhiWidth2;
  float scEtaWidth1; float scEtaWidth2;
  int   scIxSeed1;   int   scIxSeed2;
  int   scIySeed1;   int   scIySeed2;
  int   scIzSeed1;   int   scIzSeed2;
  float r91;         float r92;
  float fBrem1;      float fBrem2;
  float r41;         float r42;
  
  bool isConv1;     bool isConv2;
  bool hasPixSeed1; bool hasPixSeed2;
  bool inCrack1;    bool inCrack2;
  float nXtlSC1; float nXtlSC2;
  float preshE1; float preshE2;
  float nESclu1; float nESclu2;
    
  bool  hasGenMatch1; bool  hasGenMatch2;
  float genDR1;       float genDR2; 
  float genEnergy1;   float genEnergy2; 
  float genEta1;      float genEta2;
  float genPhi1;      float genPhi2;
  float genEt1;       float genEt2;
  int genMother1;     int genMother2;
  float genIso1;      float genIso2;
  float genIsoAll1;   float genIsoAll2;
  bool isFrag1;       bool isFrag2;
  
  // Complementary Cone
  float cc1_eta;
  float cc1_phi;
  float cc1_hcalIso; int cc1_Nhcal;
  float cc1_trackIso;int cc1_Ntrack;
  float cc1_sumIso;  
  bool cc1_hasJet;
  
  float cc2_eta;
  float cc2_phi;
  float cc2_ecalIso; int cc2_Necal;
  float cc2_hcalIso; int cc2_Nhcal;
  float cc2_trackIso;int cc2_Ntrack;
  float cc2_sumIso;  
  bool  cc2_hasJet;
  
  float cc3_eta;
  float cc3_phi;
  float cc3_ecalIso; int cc3_Necal;
  float cc3_hcalIso; int cc3_Nhcal;
  float cc3_trackIso;int cc3_Ntrack;
  float cc3_sumIso;  
  vector<int> cc3_Vieta, cc3_Viphi;
  vector<float> cc3_Vet, cc3_Ve;
  bool cc3_hasJet;
  
  float cc4_eta;
  float cc4_phi;
  float cc4_ecalIso; int cc4_Necal;
  float cc4_hcalIso; int cc4_Nhcal;
  float cc4_trackIso;int cc4_Ntrack;
  float cc4_sumIso;  
  vector<int> cc4_Vieta, cc4_Viphi;
  vector<float> cc4_Vet, cc4_Ve;
  bool cc4_hasJet;
  
  // EcalIso of the photons  
  float ecalIso1, ecalIso2;
  float ecalIso1_280, ecalIso2_280;
  float ecalIso1_320, ecalIso2_320;
  float ecalIso1_350, ecalIso2_350;
  float ecalIso1_400, ecalIso2_400;
  float ecalIsoCMS1;
  int ecalIso1_n, ecalIso2_n;
  int ecalIso1_ntk, ecalIso2_ntk;
  vector<int>   ecalIso1_Vieta, ecalIso2_Vieta;
  vector<int>   ecalIso1_Viphi, ecalIso2_Viphi;
  vector<float> ecalIso1_Veta, ecalIso2_Veta;
  vector<float> ecalIso1_Vphi, ecalIso2_Vphi;
  vector<float> ecalIso1_Vet, ecalIso2_Vet;
  vector<float> ecalIso1_Ve, ecalIso2_Ve;
  vector<int>   ecalIso1_Vieta_nv, ecalIso2_Vieta_nv;
  vector<int>   ecalIso1_Viphi_nv, ecalIso2_Viphi_nv;
  vector<float> ecalIso1_Veta_nv, ecalIso2_Veta_nv;
  vector<float> ecalIso1_Vphi_nv, ecalIso2_Vphi_nv;
  vector<float> ecalIso1_Vet_nv, ecalIso2_Vet_nv;
  vector<float> ecalIso1_Ve_nv, ecalIso2_Ve_nv;
  vector<float> ecalIso1_Vtkpt, ecalIso2_Vtkpt;
  vector<bool>  ecalIso1_VtkCnt, ecalIso2_VtkCnt; // count in trackIso ?
  vector<float> ecalIso1_Vtkd0, ecalIso2_Vtkd0;
  vector<float> ecalIso1_Vtkz0, ecalIso2_Vtkz0;
  vector<float> ecalIso1_VtkEtaH, ecalIso2_VtkEtaH;
  vector<float> ecalIso1_VtkPhiH, ecalIso2_VtkPhiH;
  vector<float> ecalIso1_VtkDEta, ecalIso2_VtkDEta;
  vector<float> ecalIso1_VtkDPhi, ecalIso2_VtkDPhi;
  vector<float> ecalIso1_VtkDR, ecalIso2_VtkDR;
  vector<int>   ecalIso1_VtkNh, ecalIso2_VtkNh;
  vector<int>   ecalIso1_VtkIh, ecalIso2_VtkIh; // expInnHits
  vector<float> ecalIso1_VtkChi2, ecalIso2_VtkChi2;
  vector<int>   ecalIso1_VtkNdof, ecalIso2_VtkNdof;
  vector<int>   ecalIso1_VtkQ, ecalIso2_VtkQ;
  vector<bool>  ecalIso1_VtkFP, ecalIso2_VtkFP;
  vector<bool>  ecalIso1_VtkSV, ecalIso2_VtkSV;
  vector<bool>  ecalIso1_VtkInV, ecalIso2_VtkInV;
  vector<int>   ecalIso1_VtkVid, ecalIso2_VtkVid;
  vector<float> ecalIso1_VtkDepE, ecalIso2_VtkDepE;
  vector<float> ecalIso1_VtkDepEt, ecalIso2_VtkDepEt;
  int ecalIso1_nTkGinV, ecalIso2_nTkGinV;
  int ecalIso1_nTkGnoZ, ecalIso2_nTkGnoZ;
  int ecalIso1_nTkG, ecalIso2_nTkG;
  int ecalIso1_nTkG_28, ecalIso2_nTkG_28;
  int ecalIso1_nTkG_32, ecalIso2_nTkG_32;
  int ecalIso1_nTkG_40, ecalIso2_nTkG_40;
  float eIso1_RECO, eIso2_RECO;

  // EcalIso of the randomCone cc1
  float cc1_ecalIso; 
  int cc1_ecalIso_n;
  int cc1_ecalIso_ntk;
  int cc1_ecalIso_nTkG;
  vector<int>   cc1_ecalIso_Vieta;
  vector<int>   cc1_ecalIso_Viphi;
  vector<float> cc1_ecalIso_Veta;
  vector<float> cc1_ecalIso_Vphi;
  vector<float> cc1_ecalIso_Vet;
  vector<float> cc1_ecalIso_Ve;
  vector<int>   cc1_ecalIso_Vieta_nv;
  vector<int>   cc1_ecalIso_Viphi_nv;
  vector<float> cc1_ecalIso_Veta_nv;
  vector<float> cc1_ecalIso_Vphi_nv;
  vector<float> cc1_ecalIso_Vet_nv;
  vector<float> cc1_ecalIso_Ve_nv;
  vector<float> cc1_ecalIso_Vtkpt;
  vector<float> cc1_ecalIso_Vtkd0;
  vector<float> cc1_ecalIso_Vtkz0;
  vector<int>   cc1_ecalIso_VtkNh;
  vector<int>   cc1_ecalIso_VtkIh;
  vector<float> cc1_ecalIso_VtkChi2;
  vector<int>   cc1_ecalIso_VtkNdof;
  vector<int>   cc1_ecalIso_VtkQ;
  vector<bool>  cc1_ecalIso_VtkFP;
  vector<bool>  cc1_ecalIso_VtkSV;
  vector<float> cc1_ecalIso_VtkEtaH;
  vector<float> cc1_ecalIso_VtkPhiH;
  vector<bool>  cc1_ecalIso_VtkInV;
  vector<int>   cc1_ecalIso_VtkVid;
  vector<float> cc1_ecalIso_VtkDepE;
  vector<float> cc1_ecalIso_VtkDepEt;
  
  // Diphoton
  float diphoMass;
  float diphoQt;
  float diphoEta;
  float deltaY;
  float diphoY;
  float diphoPhi;
  float deltaEta;
  float deltaPhi; float protoDphi;
  float deltaR;
  float cosTheta;
  float cosThetaH;
  float cosThetaS;
  float cosThetaP;
     
  // Jets
  float hardestJetPt;
  float secondJetPt;
  bool VetoJet;
  
  bool  hasJetMatch1; bool  hasJetMatch2;
  float jetDR1;       float jetDR2; 
  float jetEnergy1;   float jetEnergy2; 
  float jetEta1;      float jetEta2;
  float jetPhi1;      float jetPhi2;
  float jetEt1;       float jetEt2;
  
  ClassDef( DiphotonAnalysis, 0 )  
};

#endif

     
