#include <cassert>
#include <string>
#include <typeinfo>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <cassert>
#include <vector>
#include <map>
#include <algorithm>

using namespace std;

#include <TVector2.h>
#include <TVector3.h>
#include <TFile.h>
#include <TBranch.h>
#include <TKey.h>
#include <TTree.h>
#include <TString.h>
#include <TBits.h>
#include <Math/VectorUtil.h>
#include <TMath.h>
#include <TVector3.h>
#include <TLorentzVector.h>

#include "Analysis/utils/Constants.hh"

using namespace std;

class CrossSection
{
  string _name;
  float _val;
  float _err;
  float _err1;
  float _err2;
  bool _asymErr;
public:
  CrossSection( string name, float val, float err1, float err2=-1 )
    : _name(name), _val(val), _err1(err1), _err2(err2), _asymErr(_err2>0)
  {
    _err = _err1;
    if( _asymErr ) _err = 0.5*(_err1+_err2);
  }

  virtual ~CrossSection();
  void print()
  {
    printf( "XS[%-10s]= %6.2f +/- %-5.2f pb", _name.c_str(), _val, _err );
  }
};

int main()
{
  map< string, float >   xsect;
  map< string, float >   BR;
  map< string, float>    kf;

  //  map< string, Sample* > samples;

  /*
    Cross sections (pb) from
    https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSections
    and Matthieu's twiki
    https://twiki.cern.ch/twiki/bin/view/Main/SaclayEWKAnalysis
  */
  xsect.clear();

  xsect["W+_ln_NNLO"] = 18456;  // NNLO (FEWZ)
  xsect["W-_ln_NNLO"] = 12858;  // NNLO (FEWZ)
  xsect["W_ln_NNLO"]  = 31314;  // NNLO (FEWZ)
  xsect["W_Ln_NNLO"]  = 31314*2./3.;  // NNLO (FEWZ)

  xsect["W_ln"]  = 24640;        // W+jets (Madgraph LO)  kfacteur = 1.27   
  xsect["W_Ln"]  = 24640*2./3.;  // W+jets (Madgraph LO)
  
  //  xsect["Z_2l"]  = 4998; // > 20 GeV
  //  xsect["Z_2l"]  = 3048; // > 50 GeV
  xsect["Z_2l_[60-120]"]  = 2916; // in 60-120 window
  xsect["Z_2L_[60-120]"]  = 2916*2./3.; // in 60-120 window
  xsect["Z_2l"]  = 1614*3; // from Matthieu's twiki (POWHEG NLO)
  xsect["Z_2L"]  = 1614*2; // from Matthieu's twiki (POWHEG NLO)
  xsect["ZJets"] = 2321*3.837; // from Matthieu's twiki * ratio GHM!!!
  
  //  xsect["ttbar"] = 157.5; // inclusive  NLO (MCFM)
  xsect["ttbar"] = 173.0; // CMS measurement 2010
  
  xsect["t(s-channel)"] =  4.6;
  xsect["t(t-channel)"] = 64.6;
  xsect["tW"]           = 10.6;
  
  // gluon-gluon contributions to ZZ and WW (applied to LO)
  float ggWW = 1.20; // a la louche
  float ggZZ = 1.12; // comme Higgs en ZZ

  xsect["WW"]      = 27.8*ggWW;
  xsect["WW_2l2n"] =  2.9*ggWW; // 2l means ee(1), mm(1), tt(1), em(2), et(2), mt(2)
  xsect["WW_NLO"]  = 43;  // inclusive  MCFM (NLO)
  
  xsect["WZ"]     = 10.4;  // Pythia (LO)
  xsect["WZ_3ln"] = 0.34;  // Pythia (LO) 
  xsect["WZ_[60-120]"]     = xsect["WZ"]     * 0.9664; // GHM from mcTruth
  xsect["WZ_3ln_[60-120]"] = xsect["WZ_3ln"] * 0.9642; // ---
  xsect["W+Z_NLO_40"] = 11.8;  // mll>40 inclusive MCFM (NLO)
  xsect["W-Z_NLO_40"] =  6.4;  // mll>40 inclusive MCFM (NLO)
  xsect["WZ_NLO_40"]  = 18.2;  // mll>40 inclusive MCFM (NLO)
  xsect["WZ_LO_40"]   = xsect["WZ"] * 0.9838;  // GHM from mcTruth
  
  xsect["ZZ"]       =  4.3*ggZZ;     // Pythia (LO) mll>~12.5 GeV
  xsect["ZZ_2l2n"]  =  0.1913*ggZZ;  // Pythia (LO) mll>~12.5 GeV
  xsect["ZZ_[60-120]"]      = xsect["ZZ"]     *0.9340*ggZZ; // GHM from mcTruth !!!
  xsect["ZZ_2l2n_[60-120]"] = xsect["ZZ_2l2n"]*0.8979*ggZZ; // ---              !!!
  xsect["ZZ_NLO_40"]  =  5.9;  // mll>40 inclusive MCFM (NLO)
  xsect["ZZ_LO_40"]  =  xsect["ZZ"]*0.9546*ggZZ;  // mll>40 inclusive MCFM (NLO) !!!

  xsect["VGam"] = 173; // Madgraph, PhotonVJets
  
  /*
    branching ratios
  */
  BR.clear();
  
  BR["W_en"] = Constants::W_br_enu;
  BR["W_mn"] = Constants::W_br_munu;
  BR["W_tn"] = Constants::W_br_taunu;
  BR["W_ln"] = BR["W_en"]+BR["W_mn"]+BR["W_tn"];
  
  BR["Z_2e"] = Constants::Z0_br_ee;
  BR["Z_2m"] = Constants::Z0_br_mumu;
  BR["Z_2t"] = Constants::Z0_br_tautau;
  BR["Z_2l"] = BR["Z_2e"]+BR["Z_2m"]+BR["Z_2t"]; 
  
  BR["WW_2l2n"] =    
    Constants::W_br_enu   * Constants::W_br_enu 
    +     Constants::W_br_munu  * Constants::W_br_munu
    +     Constants::W_br_taunu * Constants::W_br_taunu
    + 2 * Constants::W_br_enu   * Constants::W_br_munu 
    + 2 * Constants::W_br_enu   * Constants::W_br_taunu 
    + 2 * Constants::W_br_munu  * Constants::W_br_taunu; 
  
  BR["WZ_3ln"] = Constants::Z0_br_ll * Constants::W_br_lnu;
  
  BR["ZZ_4e"]   = Constants::Z0_br_ee       * Constants::Z0_br_ee;
  BR["ZZ_4m"]   = Constants::Z0_br_mumu       * Constants::Z0_br_mumu;
  BR["ZZ_2e2m"] = 2 * Constants::Z0_br_ee       * Constants::Z0_br_mumu;
  BR["ZZ_4L"] = BR["ZZ_4e"] + BR["ZZ_4m"] + BR["ZZ_2e2m"]; 
  BR["ZZ_4l"] = 
    BR["ZZ_4L"]
    +     Constants::Z0_br_tautau   * Constants::Z0_br_tautau
    + 2 * Constants::Z0_br_ee       * Constants::Z0_br_tautau 
    + 2 * Constants::Z0_br_mumu     * Constants::Z0_br_tautau;
  
  BR["ZZ_2l2n"] = 2*Constants::Z0_br_nunu*Constants::Z0_br_ll;
  BR["ZZ_2e2n"] = 2*Constants::Z0_br_nunu*Constants::Z0_br_ee;
  BR["ZZ_2m2n"] = 2*Constants::Z0_br_nunu*Constants::Z0_br_mumu;
  BR["ZZ_2L2n"] = BR["ZZ_2e2n"] + BR["ZZ_2m2n"];
  
  kf.clear();
  kf[ "Z_2t" ]    = 1.032*0.95; // FEZW->r mll_20->60 : 0.6098->931/984=0.95 (NNLOkF=1.03)
  kf[ "Z_2e" ]    = 1.032*0.95;
  kf[ "Z_2m" ]    = 1.032*0.95;

  //  if( VVDynamicKFactors )
  //    {
  //      kf[ "ZZ_2l2n" ] = 1.;
  //      kf[ "WW_2l2n" ] = 1.;
  //      kf[ "WZ_3ln" ]  = 1.;
  //    }
  //  else
    {
      kf[ "ZZ_2l2n" ] = 1.37;
      kf[ "WW_2l2n" ] = 1.48; //2.9 /   4.3 pb tot  -> 1.48
      kf[ "WZ_3ln" ]  = 1.73;
    }
  kf[ "VGam" ]    = 1.11; //0.95 
  kf[ "WJets" ]   = 1.27;
  //  kf[ "ttbar" ]   = 173./165.; //1.3
  kf[ "ttbar" ]   = 1.; // this is a measurement!
  //  kf[ "T_tchan" ] = 64.6/21; //1.75
  //  kf[ "T_schan" ] = 4.6/0.99; //1.41
  //  kf[ "tW" ]      = 10.6/10.56; //1.75    //all top : 153.55 / 173  -> 1.13
  kf[ "T_tchan" ] = 1;
  kf[ "T_schan" ] = 1;
  kf[ "tW" ]      = 1;


  cout << "Hello World " << endl;
  //  cout << "This is a simple test for lumi=" << lumi << endl;

  cout << "xsect[ZZ_2l2n]=" << xsect["ZZ_2l2n_[60-120]"] << " pb, Mll in [60-120] GeV" << endl;
  cout << "xsect[ZZ_2l2n]=" << xsect["ZZ_2l2n"] << " pb" << endl;
  cout << "xsect[ZZ]=" << xsect["ZZ_[60-120]"] << " pb, Mll in [60-120] GeV" << endl;
  cout << "xsect[ZZ]*BR[2l2n]=" << xsect["ZZ_[60-120]"] * BR["ZZ_2l2n"] << " pb, Mll in [60-120] GeV" << endl;
  cout << "xsect[ZZ]*BR[4l]=" << xsect["ZZ_[60-120]"] * BR["ZZ_4l"] << " pb, Mll in [60-120] GeV" << endl;
  cout << "xsect[WZ]*BR[3ln]=" << xsect["WZ_[60-120]"] * BR["WZ_3ln"] << " pb, Mll in [60-120] GeV" << endl;
  cout << "xsect[WW]*BR[2l2n]=" << xsect["WW"] * BR["WW_2l2n"] << " pb" << endl;


  cout << "xsect[Z_2l]=" <<  xsect["Z_2l_[60-120]"] << " pb, l=e,mu,tau, Mll in [60-120] GeV" << endl; 
  cout << "xsect[Z_2l]=" <<  xsect["Z_2l"] << " pb, l=e,mu,tau, Mll > 20 GeV" << endl; 
  cout << "xsect[Z_2L]=" <<  xsect["Z_2L"] << " pb, l=e,mu      Mll > 20 GeV" << endl; 

  cout << "Z_2l/ZZ_2l2n=" <<  xsect["Z_2l_[60-120]"]/(xsect["ZZ_2l2n_[60-120]"]*kf[ "ZZ_2l2n" ]) 
       << ", l=e,mu,tau, Mll in [60-120] GeV with kfact" <<  endl;
  cout << "Z_2l/ZZ_2l2n=" <<  xsect["Z_2l"]/xsect["ZZ_2l2n"] << ", l=e,mu,tau, Pythia x-sect without kfact" <<  endl;

  cout << "W_ln/ZZ_2l2n=" <<  (xsect["W_ln"]*kf["WJets"])/(xsect["ZZ_2l2n"]*kf["ZZ_2l2n"]) << ", l=e,mu,tau, with kfact" <<  endl;
  cout << "W_ln/ZZ_2l2n=" <<  (xsect["W_Ln"]*kf["WJets"])/(xsect["ZZ_2l2n"]*BR["ZZ_2L2n"]/BR["ZZ_2l2n"]*kf["ZZ_2l2n"]) << ", l=e,mu, with kfact" <<  endl;

  cout << "ttbar/ZZ_2l2n=" <<  (xsect["ttbar"]*kf["ttbar"])/(xsect["ZZ_2l2n"]*kf["ZZ_2l2n"]) << ", l=e,mu,tau, with kfact" <<  endl;
  cout << "WW_2l2n/ZZ_2l2n=" <<  (xsect["WW_2l2n"]*kf["WW_2l2n"])/(xsect["ZZ_2l2n"]*kf["ZZ_2l2n"]) << ", l=e,mu,tau, with kfact" <<  endl;
  cout << "WZ_3ln/ZZ_2l2n=" <<  (xsect["WZ_3ln"]*kf["WZ_3ln"])/(xsect["ZZ_2l2n"]*kf["ZZ_2l2n"]) << ", l=e,mu,tau, with kfact " <<  endl;


  cout << "BR[ZZ_2e2n]=" << BR["ZZ_2e2n"] << endl;
  cout << "BR[ZZ_2m2n]=" << BR["ZZ_2m2n"] << endl;
  cout << "BR[ZZ_2e2m]=" << BR["ZZ_2e2m"] << endl;
  cout << "BR[ZZ_4e]=" << BR["ZZ_4e"] << endl;
  cout << "BR[ZZ_4m]=" << BR["ZZ_4m"] << endl;
  cout << "BR[ZZ_4L]=" << BR["ZZ_4L"] << endl;
  cout << "BR[ZZ_2l2n]=" << BR["ZZ_2l2n"] << endl;
  cout << "BR[ZZ_2L2n]=" << BR["ZZ_2L2n"] << endl;
  cout << "BR[ZZ_2L2n]/BR[ZZ_4L]=" << BR["ZZ_2L2n"]/BR["ZZ_4L"] << endl;

  float xsect_ZZ_2l2n =  xsect["ZZ"] * BR["ZZ_2l2n"]; 
  cout << "Z_2l/ZZ_2l2n  = " << xsect["Z_2l"] / xsect_ZZ_2l2n << endl;
  cout << "W_ln/ZZ_2l2n  = " << xsect["W_ln"] / xsect_ZZ_2l2n << endl; 
  cout << "ttbar/ZZ_2ln = " << xsect["ttbar"] / xsect_ZZ_2l2n << endl; 

  return 0;
}
