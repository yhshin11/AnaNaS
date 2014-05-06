
#ifndef ZZ_2l2n_Analysis_hh
#define ZZ_2l2n_Analysis_hh


#include "Analysis/utils/CutUtils.hh"
#include "Analysis/utils/Config.hh"

#include "Analysis/ghmzz/TupleAnalysis.hh"

class ZZ_2l2n_Analysis : public TupleAnalysis
{
public:

  ZZ_2l2n_Analysis( int im=0, string ip="2011", bool ZJets_cs = false ); 

  virtual ~ZZ_2l2n_Analysis();

  void setATGC( float f4Z, float f5Z ) 
  {
    if( ew==0 ) 
      {
	string aTgcDB = Config::confPath + "aTGCDatabase.txt";
	ew = new ATGCWeight( aTgcDB );
      }
    ew->setATGC( f4Z, f5Z );
  }
          
  void analyze_( float lumi );

private:

  ATGCWeight* ew;

  bool _ZJets_cs;

  int run;
  int event;
  
  int mc_event;
  int mc_run;
  float mc_mVV;  
  float mc_mV1;
  float mc_ptV1;
  float mc_pzV1;
  float mc_EV1;
  float mc_yV1;
  float mc_mV2;
  float mc_ptV2;
  float mc_pzV2;
  float mc_EV2;
  float mc_yV2;
  vector<int>* mc_pdgL;
  vectorFloat* mc_ptL;
  vectorFloat* mc_etaL;
  vectorFloat* mc_phiL;
  
  int nZCand;        // number of Z candidates
  int nVertex;       // number of vertices
  float mll;           // di-lepton mass
  float qTll;          // di-lepton transverse momentum
  float phill;         // di-lepton phi
  float yll;           // di-lepton rapidity
  float hll;           // di-lepton helicity
  float MET;           // PF-MET type-I
  float METPhi;        // PF-MET type-I - phi
  float sumEt;         // sum of transverse energies
  float mEtSig;        // sigma of MET
  float mTZZ;
  float projMET;       // projected MET
  float corProjMET;    // PU-vertex corrected MET
  float sigMET;        // corrected MET significance
  float ghmMET;        // all kind of "corrected" METs
  float cpuMETSum;
  float cpuMETSqrt;
  float cpuMETVtx;
  float cpuMETSumMin;
  float cpuMETSqrtMin;
  float cpuMETVtxMin;
  float trueMET;
  float unclMET;
  float ghmSig;
  float dPhiMin;       // minimum angle of MET to a lepton
  float balZMET;       // "balance" Z/MET
  bool isTight;   // is the Z candidate tight?
  bool isLLTight;   // additional criteria on electrons (NEW)
  int nLoose;        // number of additional loose leptons
  int jMult40;       // number of jets above 40 GeV-Et
  int jMult35;       // number of jets above 35 GeV-Et
  int jMult30;       // number of jets above 30 GeV-Et
  int jMult25;       // number of jets above 25 GeV-Et
  int jMult20;       // number of jets above 25 GeV-Et
  int jMult15;       // number of jets above 15 GeV-Et
  int jMult10;       // number of jets above 10 GeV-Et
  float JZB15;
  float JZB20;
  float JZB25;
  float dPhiJMET;      // angle between leading jet and MET
  float jEt_1;
  float jEt_2;
  float jEt_3;
  float jEt_4;
  float dPhiSumJMET;
  float dPhiZMET;
  float thirdL_MT;
  float jBtagTkC_max;
  float jBtagSoftM_max;
  float bal;

  vectorFloat* jEt;
  vectorFloat* jEta;
  vectorFloat* jPhi;
  vectorFloat* jDR1;
  vectorFloat* jDR2;
  vectorFloat* jDPhiMET;
  vectorFloat* jBtagTkC;
  vectorFloat* jBtagSoftM;
  
  vectorFloat* LPt;
  vector<int>* LPdg;
  vectorFloat* LEta;
  vectorFloat* LPhi;
  vectorFloat* LMT;
  vector<int>* LID;
  vectorFloat* LTrkIso;
  vectorFloat* LEcalIso;
  vectorFloat* LHcalIso;

  // event weights
  float eventWeight;
  float aTGCfactor;
  float kfactor;
  float kfactor_jv30;

  //
  bool isMuMu;
  bool isEE;
  bool isLL;
  bool isDYSimulation;
  bool isDibosonSimulation;
  bool isData;
  bool isOtherSimulation;
  bool analyzeSignalMC();

  // is it in acceptance?
  bool _signalOnly;
  bool isZZ;
  bool isWZ;
  bool isWW;
  bool isQt30;
  bool isInAcceptance;

  

  
  ClassDef( ZZ_2l2n_Analysis, 0 )  
};

#endif
