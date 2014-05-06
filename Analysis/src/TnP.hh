#ifndef TnP_h
#define TnP_h

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

enum eop { kE=0, kP=1 };

class TnP : public SampleAnalysis
{
public:

  TnP( Sample& sample, const std::string & collectionFileName );
  virtual ~TnP();

  static int imode;
 
private:

  //Gerenal variables
  std::vector<std::pair<int,int> > Events;
  bool debug;
  bool isData;
  TRandom3 rg;

  //Compteurs 
  int Nevt_;
  int Nnosc_;
  int NHLT_;
  int Npresel_;
  
  //Elements
  vector<Candidate*> tags;
  vector<Candidate*> probes;
  vector<Candidate*> eprobes;
  Candidate* tag_;
  Candidate* probe_;
  
  //Inerited functions
  virtual bool analyzeEvent();
  virtual void writeHistograms();
  virtual void bookHistograms();
  virtual bool preselection();

  //Name of sample (for deta corrections)
  std::string SampleName;

  //Functions
  bool  EventPreselection();
  void  initRootTree();
  void  fillTree(); 
  bool  isSpike(Candidate*);
  bool  isTag(Candidate*);

  // constants
  float _minPt;       // cut on tag and probe
  float _minM, _maxM; // Z mass window
  
  // outputs
  TTree outTree_;

  // general stuff
  int run;
  int event;
  int nProbes;

  // vertex
  int nGoodVtx;
  float vX, vY, vZ;
    
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

  // tag variables
  

  // probe variables
  float p_eta;
  float p_phi;
  float p_pt;
  float p_sigIeta;
  float p_hcalIso;
  float p_ecalIso;
  float p_trackIso;
  float p_r9;
  float p_hoe;

  // common variables
  float mass;


  ClassDef( TnP, 0 )  
};

#endif

     
