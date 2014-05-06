#include <unistd.h> 
#include <iostream>
#include <sstream>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <typeinfo>
#include <string>
#include <cmath>
#include <algorithm>
#include <functional>
#include <vector>
#include <map>

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

#include "Analysis/ghmzz/TupleAnalysis.hh"

#include "Analysis/core/Sample.hh"
#include "Analysis/utils/Constants.hh"
#include "Analysis/utils/Config.hh"
#include "Analysis/utils/KineUtils.hh"
#include "Analysis/utils/RooUtils.hh"
#include "Analysis/utils/HistoManager.hh"
#include "Analysis/utils/TupleManager.hh"
#include "Analysis/utils/DatasetManager.hh"
#include "Analysis/utils/CutUtils.hh"

ClassImp( TupleAnalysis )

void
TupleAnalysis::defineTemplate( string templ, 
			       int nbin, float min, float max )
{
  cout << "1D-Template[" << templ << "]-->(nbin=" << nbin << "/min=" << min << "/max=" << max << ")" << endl; 

  TH1F* h_ = new TH1F( templ.c_str(), templ.c_str(), 
		       nbin, min, max );
  h_->Sumw2();
  hm.addTemplate<TH1F>( templ, h_ );
}

void
TupleAnalysis::defineTemplate( string templ, 
			       int nbin, float* xbins )
{
  cout << "1D-Template[" << templ << "]-->(nbin=" << nbin << "variable size" << endl; 
  TH1F* h_ = new TH1F( templ.c_str(), templ.c_str(), 
		       nbin, xbins );
  h_->Sumw2();
  hm.addTemplate<TH1F>( templ, h_ );
}

TH1* 
TupleAnalysis::hist( string templ, string cut, string dir, string prefix )
{
  TH1* h_ = (TH1*)hm.h<TH1F>( templ, cut, dir, prefix );
  assert( h_!=0 );
  return h_;
}

void
TupleAnalysis::fill( string templ, string cut, string dir, string prefix,
		     float val, float weight, bool overflow, bool norm )
{
  TH1* h_ = hist( templ, cut, dir, prefix );
  
  TAxis* ax_ = h_->GetXaxis();
  int ibin_ = ax_->FindBin(val);
  int N_ = ax_->GetNbins();
  if( ibin_>N_ )  // overflow
    {
      if( overflow ) 
	{
	  h_->SetBinContent( N_, h_->GetBinContent(N_) 
			     + weight * ax_->GetBinWidth(1)/ax_->GetBinWidth(N_));
	  return;
	}
    }  
  if( ax_->IsVariableBinSize() )
    {
      float w_ = 1;
      if( norm )
	{
	  if( ibin_>0 && ibin_<=N_ ) 
	    w_  = ax_->GetBinWidth(1)/ax_->GetBinWidth(ibin_);
	}
      h_-> Fill( val, weight * w_ );
    }
  else
    {
      h_->Fill( val, weight );
    }
}

void 
TupleAnalysis::setOutputFile()
{
  string histfile = Config::histPath + analysis;
  histfile += ".root";
  
  // output file 
  TFile* outputfile = TFile::Open( histfile.c_str(), "RECREATE" );
  
  // histo manager
  hm.setFile( outputfile );

  // dataset manager
  dm.setFile( outputfile, DatasetManager::kWrite );

  // ntuple manager
  tmsel.setFile( outputfile, TupleManager::kWrite );

}

void
TupleAnalysis::defineSamples()
{
  cout << "Entering defineSamples " << endl;

  VVDynamicKFactors = true;

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

  if( VVDynamicKFactors )
    {
      kf[ "ZZ_2l2n" ] = 1.;
      kf[ "WW_2l2n" ] = 1.;
      kf[ "WZ_3ln" ]  = 1.;
    }
  else
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

  // define the samples
  samples.clear();

  //
  // the number of events are extracted from the log file
  // ---> eventually, from the root file itself
  // ---> create a "stat" ntuple
  //

  if( period=="2010" )
    {
      samples["Electron_1"] = 
	new Sample("EG_Dec22_132440_141000/data", "data" );
      samples["Electron_2"] = 
	new Sample("EG_Dec22_141001_143000/data", "data" );
      samples["Electron_3"] =
	new Sample("EG_Dec22_143001_144114/data", "", "", 0, 0, 0, 0,
		   4743800 );
      samples["Electron_4"] = 
	new Sample("Electron_Dec22_146240_147000/data", "", "", 0, 0, 0, 0,
		   2424964 );  
      samples["Electron_5"] = 
	new Sample("Electron_Dec22_147001_148000/data", "", "", 0, 0, 0, 0,
		   3545069 );  
      samples["Electron_6"] = 
	new Sample("Electron_Dec22_148001_149000/data", "", "", 0, 0, 0, 0,
		   2652499 );      
      samples["Electron_7"] = 
	new Sample("Electron_Dec22_149001_149442/data", "", "", 0, 0, 0, 0,
		   2955861 );
    
      samples["Muon_1"] = 
	new Sample("Mu_Dec22_132440_141000/data", 0, 0, 0, 0,
		   100668 );
      samples["Muon_2"] = 
	new Sample("Mu_Dec22_141001_143000/data", 0, 0, 0, 0,
		   179706 );
      samples["Muon_3"] = 
	new Sample("Mu_Dec22_143001_144114/data", 0, 0, 0, 0,
		   330878 );
      samples["Muon_4"] = 
	new Sample("Mu_Dec22_146240_147000/data", 0, 0, 0, 0,
		   446775 );
      samples["Muon_5"] = 
	new Sample("Mu_Dec22_147001_148000/data", 0, 0, 0, 0,
		   683728 );
      samples["Muon_6"] = 
	new Sample("Mu_Dec22_148001_149000/data", 0, 0, 0, 0,
		   623134 );
      samples["Muon_7"] = 
	new Sample("Mu_Dec22_149001_149442/data", 0, 0, 0, 0,
		   758453 );
    }
  else if( period=="2011A" )
    {
      samples["SingleEle_1"] = 
	new Sample("SingleEle_160404_162000/data", 0, 0, 0, 0,
		   0 );
      samples["SingleEle_2"] = 
	new Sample("SingleEle_162001_163000/data", 0, 0, 0, 0,
		   0 );
      samples["SingleEle_3"] = 
	new Sample("SingleEle_163001_163869/data", 0, 0, 0, 0,
		   0 );

      // double electron
      samples["DoubleEle_1"] = 
	new Sample("DoubleEle_160404_162000/data", 0 );
      samples["DoubleEle_2"] = 
	new Sample("DoubleEle_162001_163000/data", 0 );
      samples["DoubleEle_3"] = 
	new Sample("DoubleEle_163001_163869/data", 0 );
      samples["DoubleEle_4"] = 
	new Sample("DoubleEle_165071_165400/data", 0 );
      samples["DoubleEle_5"] = 
	new Sample("DoubleEle_165401_165700/data", 0 );
      samples["DoubleEle_6"] = 
	new Sample("DoubleEle_165543_165900/data", 0 );
      samples["DoubleEle_7"] = 
	new Sample("DoubleEle_165901_166200/data", 0 );
      samples["DoubleEle_8"] = 
	new Sample("DoubleEle_166201_166502/data", 0 );
      samples["DoubleEle_9"] = 
	new Sample("DoubleEle_166503_166600/data", 0 );
      samples["DoubleEle_10"] = 
	new Sample("DoubleEle_166601_166700/data", 0 );
      samples["DoubleEle_11"] = 
	new Sample("DoubleEle_166701_166861/data", 0 );
      samples["DoubleEle_12"] = 
	new Sample("DoubleEle_166862_167000/data", 0 );
      samples["DoubleEle_13"] = 
	new Sample("DoubleEle_167001_167151/data", 0 );
      samples["DoubleEle_14"] = 
	new Sample("DoubleEle_167152_167450/data", 0 );
      samples["DoubleEle_15"] = 
	new Sample("DoubleEle_167451_167784/data", 0 );
      samples["DoubleEle_16"] = 
	new Sample("DoubleEle_167785_167913/data", 0 );
      // no data for runs 167914-170500
      samples["DoubleEle_17"] = 
	new Sample("DoubleEle_170501_171000/data", 0 );
      samples["DoubleEle_18"] = 
	new Sample("DoubleEle_171001_171500/data", 0 );
      samples["DoubleEle_19"] = 
	new Sample("DoubleEle_171501_172000/data", 0 );
      samples["DoubleEle_20"] = 
	new Sample("DoubleEle_172001_172619/data", 0 );
      samples["DoubleEle_21"] = 
	new Sample("DoubleEle_172620_172900/data", 0 );
      samples["DoubleEle_22"] = 
	new Sample("DoubleEle_172901_173244/data", 0 );
      samples["DoubleEle_23"] = 
	new Sample("DoubleEle_173245_173450/data", 0 );
      samples["DoubleEle_24"] = 
	new Sample("DoubleEle_173451_173692/data", 0 );
      samples["DoubleEle_25"] = 
	new Sample("DoubleEle_175701_176023/data", 0 );
      samples["DoubleEle_26"] = 
	new Sample("DoubleEle_176023_176309/data", 0 );

      samples["SingleMu_1"] = 
	new Sample("SingleMu_160404_162000/data", 0, 0, 0, 0,
		   49455 );
      samples["SingleMu_2"] = 
	new Sample("SingleMu_162001_163000/data", 0, 0, 0, 0,
		   48526 );
      samples["SingleMu_3"] = 
	new Sample("SingleMu_163001_163869/data", 0, 0, 0, 0,
		   1461300 );

      samples["DoubleMu_1"] = 
	new Sample("DoubleMu_160404_162000/data", 0 );
      samples["DoubleMu_2"] = 
	new Sample("DoubleMu_162001_163000/data", 0 );
      samples["DoubleMu_3"] = 
	new Sample("DoubleMu_163001_163869/data", 0 );
      samples["DoubleMu_4"] = 
	new Sample("DoubleMu_165071_165400/data", 0 );
      samples["DoubleMu_5"] = 
	new Sample("DoubleMu_165401_165700/data", 0 );
      samples["DoubleMu_6"] = 
	new Sample("DoubleMu_165543_165900/data", 0 );
      samples["DoubleMu_7"] = 
	new Sample("DoubleMu_165901_166200/data", 0 );
      samples["DoubleMu_8"] = 
	new Sample("DoubleMu_166201_166502/data", 0 );
      samples["DoubleMu_9"] = 
	new Sample("DoubleMu_166503_166600/data", 0 );
      samples["DoubleMu_10"] = 
	new Sample("DoubleMu_166601_166700/data", 0 );
      samples["DoubleMu_11"] = 
	new Sample("DoubleMu_166701_166861/data", 0 );
      samples["DoubleMu_12"] = 
	new Sample("DoubleMu_166862_167000/data", 0 );
      samples["DoubleMu_13"] = 
	new Sample("DoubleMu_167001_167151/data", 0 );
      samples["DoubleMu_14"] = 
	new Sample("DoubleMu_167152_167450/data", 0 );
      samples["DoubleMu_15"] = 
	new Sample("DoubleMu_167451_167784/data", 0 );
      samples["DoubleMu_16"] = 
	new Sample("DoubleMu_167785_167913/data", 0 );
      // no data for runs 167914-170500
      samples["DoubleMu_17"] = 
	new Sample("DoubleMu_170501_171000/data", 0 );
      samples["DoubleMu_18"] = 
	new Sample("DoubleMu_171001_171500/data", 0 );
      samples["DoubleMu_19"] = 
	new Sample("DoubleMu_171501_172000/data", 0 );
      samples["DoubleMu_20"] = 
	new Sample("DoubleMu_172001_172619/data", 0 );
      samples["DoubleMu_21"] = 
	new Sample("DoubleMu_172620_172900/data", 0 );
      samples["DoubleMu_22"] = 
	new Sample("DoubleMu_172901_173244/data", 0 );
      samples["DoubleMu_23"] = 
	new Sample("DoubleMu_173245_173450/data", 0 );
      samples["DoubleMu_24"] = 
	new Sample("DoubleMu_173451_173692/data", 0 );
      samples["DoubleMu_25"] = 
	new Sample("DoubleMu_175701_176023/data", 0 );
      samples["DoubleMu_26"] = 
	new Sample("DoubleMu_176023_176309/data", 0 );
    }

  samples["WW_2l2n"] =
    new Sample( "WW_2l2n/all", "", "", xsect["WW_2l2n"], 1, 
		1.,
		//1.,  // use kf=1 because correct versus pT
		kf["WW_2l2n"],  // !!!GHM temporary!!! 
		110000 );
  
  samples["WZ_3ln"] =
    new Sample( "WZ_3ln/all", "", "", xsect["WZ_3ln"], 1, 
		1.,
		//		1.,  // use kf=1 because correct versus pT
		kf["WZ_3ln"],  // !!!GHM temporary!!! 
		110000 );
  
  samples["ZZ_2l2n"] =
    new Sample( "ZZ_2l2n/all", "", "", xsect["ZZ_2l2n"], 1, 
		1.,
		// 1.,  // use kf=1 because correct versus pT
		kf["ZZ_2l2n"], // !!!GHM temporary!!! 
		140884 );
  
  samples["ZZ_2e2n"] =
    new Sample( "ZZ_2l2n/ZZ_2e2n", "", "", xsect["ZZ_2l2n"], 1, 
		BR["ZZ_2e2n"]/BR["ZZ_2l2n"], 
		//		1.,
		kf["ZZ_2l2n"],  
		// 46844 // this is because of a bug 
		47054
		);
  
  samples["ZZ_2m2n"] =
    new Sample( "ZZ_2l2n/ZZ_2m2n", "", "", xsect["ZZ_2l2n"], 1, 
		BR["ZZ_2m2n"]/BR["ZZ_2l2n"], 
		//1.,
		kf["ZZ_2l2n"], 
		// 47045,
		47244 
		);
  
  samples["Z_2e"] =
    new Sample( "Z_2e/all", "", "", xsect["Z_2l"], 1, 
		BR["Z_2e"]/BR["Z_2l"], kf["Z_2e"], 1992276 );
  
  samples["Z_2m"] =
    new Sample( "Z_2m/all", "", "", xsect["Z_2l"], 1, 
		BR["Z_2m"]/BR["Z_2l"], kf["Z_2m"], 1914154 );

  samples["Z_2t"] =
    new Sample( "Z_2t/all", "", "", xsect["Z_2l"], 1, 
		BR["Z_2t"]/BR["Z_2l"], kf["Z_2t"], 1388101 );
  
  samples["ZJets_2e"] =
    new Sample( "ZJets/all", "", "", xsect["ZJets"], 1, 
		BR["Z_2e"]/BR["Z_2l"], kf["Z_2e"], 2000000 );
  
  samples["ZJets_2m"] =
    new Sample( "ZJets/all", "", "", xsect["ZJets"], 1, 
		BR["Z_2m"]/BR["Z_2l"], kf["Z_2m"], 2000000 );

  samples["VGam"] =
    new Sample( "VGam/all", "", "", xsect["VGam"], 1, 
		1, kf["VGam"], 1102193 );
  
  samples["ttbar"] =
    new Sample( "ttbar/all", "", "", xsect["ttbar"], 1, 
		1, kf["ttbar"],   1165716  );

  // TChan
  // 31.59% - Sel/Hlt/Tot[W_en        ]=   131/107414/152921 --> (  0.09+/-0.01 )%
  // 31.78% - Sel/Hlt/Tot[W_mn        ]=  1975/125006/153855 --> (  1.28+/-0.03 )%
  // 31.65% - Sel/Hlt/Tot[W_tn        ]=   124/ 30142/153209 --> (  0.08+/-0.01 )%
  //100.00% - Sel/Hlt/Tot[all         ]=  2339/275832/484060 --> (  0.48+/-0.01 )%
  //  4.97% - Sel/Hlt/Tot[bkg         ]=   109/ 13270/ 24075 --> (  0.45+/-0.04 )%
  samples["t_tchan"] =
    new Sample( "T_tchan/all", "", "", xsect["t(t-channel)"], 1, 
		1, kf["T_tchan"],   484060  );

  // SChan
  // 32.11% - Sel/Hlt/Tot[W_en        ]=   122/110368/158938 --> (  0.08+/-0.01 )%
  //  0.00% - Sel/Hlt/Tot[W_en2b      ]=     0/     4/    13 --> (  0.00+/-0.00 )%
  // 32.09% - Sel/Hlt/Tot[W_mn        ]=  2377/126930/158818 --> (  1.50+/-0.03 )%
  //  0.00% - Sel/Hlt/Tot[W_mn2b      ]=     0/     0/     3 --> (  0.00+/-0.00 )%
  // 30.82% - Sel/Hlt/Tot[W_tn        ]=   231/ 33777/152551 --> (  0.15+/-0.01 )%
  //  0.00% - Sel/Hlt/Tot[W_tn2b      ]=     0/     0/     8 --> (  0.00+/-0.00 )%
  //100.00% - Sel/Hlt/Tot[all         ]=  2845/284766/494967 --> (  0.57+/-0.01 )%
  //  4.98% - Sel/Hlt/Tot[bkg         ]=   115/ 13687/ 24636 --> (  0.47+/-0.04 )%
  samples["t_schan"] =
    new Sample( "T_schan/all", "", "", xsect["t(s-channel)"], 1, 
		1, kf["T_schan"],   494967  );

//   1.14% - Sel/Hlt/Tot[WW_2e2n     ]=  2553/  5163/  5578 --> ( 45.77+/-0.67 )%
//   1.14% - Sel/Hlt/Tot[WW_2m2n     ]=  3264/  5453/  5602 --> ( 58.26+/-0.66 )%
//   0.68% - Sel/Hlt/Tot[WW_2t2n     ]=    17/  1094/  3328 --> (  0.51+/-0.12 )%
//  40.00% - Sel/Hlt/Tot[WW_4j       ]=    23/ 20503/195747 --> (  0.01+/-0.00 )%
//   2.26% - Sel/Hlt/Tot[WW_em2n     ]=   185/ 10618/ 11044 --> (  1.68+/-0.12 )%
//  13.53% - Sel/Hlt/Tot[WW_en2j     ]=    52/ 49516/ 66232 --> (  0.08+/-0.01 )%
//   2.26% - Sel/Hlt/Tot[WW_et2n     ]=   424/  8778/ 11063 --> (  3.83+/-0.18 )%
//  13.48% - Sel/Hlt/Tot[WW_mn2j     ]=   783/ 56547/ 65957 --> (  1.19+/-0.04 )%
//   2.25% - Sel/Hlt/Tot[WW_mt2n     ]=   491/  9649/ 10990 --> (  4.47+/-0.20 )%
//  13.32% - Sel/Hlt/Tot[WW_tn2j     ]=    59/ 15174/ 65179 --> (  0.09+/-0.01 )%
//   6.49% - Sel/Hlt/Tot[W_2j        ]=    64/  8763/ 31784 --> (  0.20+/-0.03 )%
//   1.07% - Sel/Hlt/Tot[W_en        ]=   299/  4191/  5248 --> (  5.70+/-0.32 )%
//   1.07% - Sel/Hlt/Tot[W_mn        ]=   429/  4637/  5261 --> (  8.15+/-0.38 )%
//   1.03% - Sel/Hlt/Tot[W_tn        ]=    59/  1954/  5061 --> (  1.17+/-0.15 )%
// 100.00% - Sel/Hlt/Tot[all         ]=  8725/202592/489417 --> (  1.78+/-0.02 )%
//   0.27% - Sel/Hlt/Tot[bkg         ]=    23/   552/  1343 --> (  1.71+/-0.35 )%  
  samples["tW"] =
    new Sample( "TW/all", "", "", xsect["tW"], 1, 
		1, kf["tW"],   489417  );

}

void
TupleAnalysis::setVVKFactors()
{
  cout << "\nK-factors for ZZ" << endl;
  string kfDB_ZZ_jv30 = Config::confPath + "kFactorDatabaseZZ_jv30.txt";
  kFactors["ZZ_jv30"] = new KFactor(kfDB_ZZ_jv30);
  string kfDB_ZZ = Config::confPath + "kFactorDatabaseZZ.txt";
  kFactors["ZZ"] = new KFactor(kfDB_ZZ);
  
  cout << "\nK-factors for WZ" << endl;
  string kfDB_WZ_jv30 = Config::confPath + "kFactorDatabaseWZ_jv30.txt";
  kFactors["WZ_jv30"] = new KFactor(kfDB_WZ_jv30);
  string kfDB_WZ = Config::confPath + "kFactorDatabaseWZ.txt";
  kFactors["WZ"] = new KFactor(kfDB_WZ);
  
  cout << "\nK-factors for WW" << endl;
  string kfDB_WW_jv30 = Config::confPath + "kFactorDatabaseWW_jv30.txt";
  kFactors["WW_jv30"] = new KFactor(kfDB_WW_jv30);
  string kfDB_WW = Config::confPath + "kFactorDatabaseWW.txt";
  kFactors["WW"] = new KFactor(kfDB_WW);
}

void
TupleAnalysis::setMomentumWeights()
{
  momWeights["data"] = new MomentumWeight("data");
  momWeights["Z_2l"] = new MomentumWeight("Z_2l");
  momWeights["Z_2e"] = new MomentumWeight("Z_2e");
  momWeights["Z_2m"] = new MomentumWeight("Z_2m");
  momWeights["Z_2t"] = new MomentumWeight("Z_2t");
  momWeights["VGam"] = new MomentumWeight("VGam");
  momWeights["tt"] = new MomentumWeight("tt");
  momWeights["t"] = new MomentumWeight("t");
  momWeights["WW_2l2n"] = new MomentumWeight("WW_2l2n");
  momWeights["WZ_3ln"] = new MomentumWeight("WZ_3ln");
  momWeights["ZZ_2l2n"] = new MomentumWeight("ZZ_2l2n");
  momWeights["ZZ_2e2n"] = new MomentumWeight("ZZ_2l2n");
  momWeights["ZZ_2m2n"] = new MomentumWeight("ZZ_2l2n");
}
void
TupleAnalysis::init_()
// 
{
  {
    // variable bin met
    float eps_ = 0.0001;
    vector<float> bins_;
    float xbin_(0);
    //    while( xbin_<=400 )
    //      {
    //	bins_.push_back(xbin_);
    //	if(      xbin_< 50-eps_ ) xbin_+=2.5; 
    //	else if( xbin_<100-eps_ ) xbin_+=5; 
    //	else if( xbin_<200-eps_ ) xbin_+=10; 
    //	else if( xbin_<300-eps_ ) xbin_+=50; 
    //	else if( xbin_<400-eps_ ) xbin_+=100; 
    //	else break;
    //      }
    while( xbin_<=150 )
      {
	bins_.push_back(xbin_);
	if(      xbin_< 50-eps_ ) xbin_+=2.5; 
	else if( xbin_<100-eps_ ) xbin_+=5; 
	else if( xbin_<150-eps_ ) xbin_+=10; 
	//	else if( xbin_<300-eps_ ) xbin_+=50; 
	//	else if( xbin_<400-eps_ ) xbin_+=100; 
	else break;
      }
    bins_.push_back(175.0);
    bins_.push_back(225.0);
    bins_.push_back(300.0);
    bins_.push_back(397.5);
    bins_.push_back(400.0);
    float bins[1000];
    size_t nbin = bins_.size();
    for( size_t ibin=0; ibin<nbin; ibin++ )
      {
	cout << "xbin[" << ibin << "]=" << bins_[ibin] << endl;
	bins[ibin] = bins_[ibin];
      }
    defineTemplate( "vmet", nbin-1, bins );
  }
  {
    // variable bin met
    vector<float> bins_;
    float xbin_(0);
    while( xbin_<=20 )
      {
	bins_.push_back(xbin_);
	if(      xbin_<4.99   ) xbin_+=0.25; 
	else if( xbin_<9.99   ) xbin_+=0.5; 
	else if( xbin_<14.99  ) xbin_+=1.0; 
	else break;
      }
    bins_.push_back(15);
    bins_.push_back(19.75);
    bins_.push_back(20);
    float bins[1000];

    size_t nbin = bins_.size();
    for( size_t ibin=0; ibin<nbin; ibin++ )
      {
	cout << "xbin[" << ibin << "]=" << bins_[ibin] << endl;
	bins[ibin] = bins_[ibin];
      }
    defineTemplate( "vsig", nbin-1, bins );
  }
  {
    //    float mmin  = 5.;
    //    int nbin    = 50;
    //    float bin   = 0.5;
    float mmin  = 5.;
    int nbin    = 63;
    float bin   = 0.5;
    float lmmin = log(mmin);
    float lbin  = log(mmin+bin) - lmmin; 

    float lm = lmmin;
    vector<float> bins_;
    bins_.push_back(mmin);
    for( int ibin=0; ibin<nbin; ibin++ )
      {
	lm += lbin;
	bins_.push_back( exp(lm) );
      }
    float bins[1000];
    int nbin_ = bins_.size();
    for( int ibin=0; ibin<nbin_; ibin++ )
      {
	cout << "xbin[" << ibin << "]=" << bins_[ibin] << endl;
	bins[ibin] = bins_[ibin];
      }
    defineTemplate( "Vmet", nbin-1, bins );
  }
  {
    //    float mmin  = 5.;
    //    int nbin    = 50;
    //    float bin   = 0.5;
    float s0   = 4.2;
    float bin  = 0.5;
    int   nbin = 60;
    float lbin = log( s0+bin ) - log( s0 );
    float lmmin = log( s0 ) - (nbin/2)*lbin;
    float mmin = exp( lmmin );
    
    //      float mmin  = 0.00001;
    //      float bin   = 0.000005;
    //      float lmmin = log(mmin);
    //      float lbin  = log(mmin+bin) - lmmin; 
    
    float lm = lmmin;
    vector<float> bins_;
    bins_.push_back(mmin);
    for( int ibin=0; ibin<nbin; ibin++ )
      {
	lm += lbin;
	bins_.push_back( exp(lm) );
      }
    float bins[1000];
    int nbin_ = bins_.size();
    for( int ibin=0; ibin<nbin_; ibin++ )
      {
	cout << "xbin[" << ibin << "]=" << bins_[ibin] << endl;
	bins[ibin] = bins_[ibin];
      }
    defineTemplate( "Vsig", nbin-1, bins );
  }    
  
  defineTemplate( "met", 400, 0., 400.   );
  defineTemplate( "mt",  1000, 0., 1000.   );
  defineTemplate( "corsig", 1000, -10., 90.   );
  defineTemplate( "cormet", 100, 0., 100.   );
  defineTemplate( "set", 700, 0., 7000. );
  defineTemplate( "sig",  80, 0., 20.    );
  defineTemplate( "mll",  60, 60, 120.   );
  defineTemplate( "pt",  1000, 0., 1000.   );
  defineTemplate( "y",   120, -3., 3.    );
  defineTemplate( "phi", 180, -180, 180. );
  defineTemplate( "bal", 60, 0., 6.      );
  defineTemplate( "n10",   11, -0.5, 10.5   );
  defineTemplate( "n20",   21, -0.5, 20.5   );
  defineTemplate( "JZB", 400, -100, 300  );
  defineTemplate( "btagmu", 100, 0, 1  );
  defineTemplate( "btagtk", 120, -10, 50  );

  // Novosibisrk fits with cpuMETVtx
  bool fixedTail = true;
  if( fixedTail )
    {
      // jet veto at 30 GeV
      {
	// data
		    
	// electrons
	tail["data_0_j30veto"] = 0.425;
	alpha_0["data_0_j30veto"] = 4.997;
	alpha_PU["data_0_j30veto"] = 0.234;
	sigma_0["data_0_j30veto"] = 3.962;
	sigma_PU["data_0_j30veto"] = 1.069;
		
	// muons
	tail["data_1_j30veto"] = 0.425;
	alpha_0["data_1_j30veto"] = 4.900;
	alpha_PU["data_1_j30veto"] = 0.245;
	sigma_0["data_1_j30veto"] = 3.773;
	sigma_PU["data_1_j30veto"] = 1.170;

	// leptons
	tail["data_2_j30veto"] = 0.425;
	alpha_0["data_2_j30veto"] = 4.946;
	alpha_PU["data_2_j30veto"] = 0.240;
	sigma_0["data_2_j30veto"] = 3.896;
	sigma_PU["data_2_j30veto"] = 1.103;

	// DY simulation
		    
	// electrons
	tail["Z_2e_0_j30veto"] = 0.425;
	alpha_0["Z_2e_0_j30veto"] = 3.807;
	alpha_PU["Z_2e_0_j30veto"] = 0.305;
	sigma_0["Z_2e_0_j30veto"] = 2.924;
	sigma_PU["Z_2e_0_j30veto"] = 1.267;
		    
	// muons
	tail["Z_2m_1_j30veto"] = 0.425;
	alpha_0["Z_2m_1_j30veto"] = 3.776;
	alpha_PU["Z_2m_1_j30veto"] = 0.315;
	sigma_0["Z_2m_1_j30veto"] = 2.908;
	sigma_PU["Z_2m_1_j30veto"] = 1.262;
		    
	// leptons
	tail["Z_2l_2_j30veto"] = 0.425;
	alpha_0["Z_2l_2_j30veto"] = 3.797;
	alpha_PU["Z_2l_2_j30veto"] = 0.310;
	sigma_0["Z_2l_2_j30veto"] = 2.905;
	sigma_PU["Z_2l_2_j30veto"] = 1.275;
		    
      }		  
      {
	// fits with fixed tail

	// data

	// electrons
	tail["data_0_j25veto"] = 0.425;
	alpha_0["data_0_j25veto"] = 4.939;
	alpha_PU["data_0_j25veto"] = 0.222;
	sigma_0["data_0_j25veto"] = 3.908;
	sigma_PU["data_0_j25veto"] = 1.014;

	//muons
	tail["data_1_j25veto"] = 0.425;
	alpha_0["data_1_j25veto"] = 4.753;
	alpha_PU["data_1_j25veto"] = 0.248;
	sigma_0["data_1_j25veto"] = 3.743;
	sigma_PU["data_1_j25veto"] = 1.085;

	// leptons
	tail["data_2_j25veto"] = 0.425;
	alpha_0["data_2_j25veto"] = 4.904;
	alpha_PU["data_2_j25veto"] = 0.225;
	sigma_0["data_2_j25veto"] = 3.862;
	sigma_PU["data_2_j25veto"] = 1.027;

	// DY simulation

	// electrons
	tail["Z_2e_0_j25veto"] = 0.425;
	alpha_0["Z_2e_0_j25veto"] = 3.823;
	alpha_PU["Z_2e_0_j25veto"] = 0.283;
	sigma_0["Z_2e_0_j25veto"] = 2.943;
	sigma_PU["Z_2e_0_j25veto"] = 1.177;

	// muons
	tail["Z_2m_1_j25veto"] = 0.425;
	alpha_0["Z_2m_1_j25veto"] = 3.691;
	alpha_PU["Z_2m_1_j25veto"] = 0.313;
	sigma_0["Z_2m_1_j25veto"] = 2.946;
	sigma_PU["Z_2m_1_j25veto"] = 1.179;

	// leptons
	tail["Z_2l_2_j25veto"] = 0.425;
	alpha_0["Z_2l_2_j25veto"] = 3.753;
	alpha_PU["Z_2l_2_j25veto"] = 0.297;
	sigma_0["Z_2l_2_j25veto"] = 2.922;
	sigma_PU["Z_2l_2_j25veto"] = 1.182;
      }
      {
	// electrons, fixed tails
	tail["data_0_qT30"] = 0.425;
	alpha_0["data_0_qT30"] = 6.824;
	alpha_PU["data_0_qT30"] = 0.236;
	sigma_0["data_0_qT30"] = 5.382;
	sigma_PU["data_0_qT30"] = 1.174;

	// muons, fixed tails
	tail["data_1_qT30"] = 0.425;
	alpha_0["data_1_qT30"] = 7.014;
	alpha_PU["data_1_qT30"] = 0.236;
	sigma_0["data_1_qT30"] = 5.770;
	sigma_PU["data_1_qT30"] = 1.046;

	// leptons, fixed tails
	tail["data_2_qT30"] = 0.425;
	alpha_0["data_2_qT30"] = 6.965;
	alpha_PU["data_2_qT30"] = 0.234;
	sigma_0["data_2_qT30"] = 5.644;
	sigma_PU["data_2_qT30"] = 1.109;


	// electrons, fixed tails
	tail["Z_2e_0_qT30"] = 0.425;
	alpha_0["Z_2e_0_qT30"] = 5.221;
	alpha_PU["Z_2e_0_qT30"] = 0.350;
	sigma_0["Z_2e_0_qT30"] = 4.225;
	sigma_PU["Z_2e_0_qT30"] = 1.360;

	// muons, fixed tails
	tail["Z_2m_1_qT30"] = 0.425;
	alpha_0["Z_2m_1_qT30"] = 5.310;
	alpha_PU["Z_2m_1_qT30"] = 0.357;
	sigma_0["Z_2m_1_qT30"] = 4.292;
	sigma_PU["Z_2m_1_qT30"] = 1.388;

	// leptons, fixed tails
	tail["Z_2l_2_qT30"] = 0.425;
	alpha_0["Z_2l_2_qT30"] = 5.248;
	alpha_PU["Z_2l_2_qT30"] = 0.353;
	sigma_0["Z_2l_2_qT30"] = 4.268;
	sigma_PU["Z_2l_2_qT30"] = 1.367;



	// electrons, fixed tails
	tail["Z_2e_0_qT30"] = 0.534;
	alpha_0["Z_2e_0_qT30"] = 4.722;
	alpha_PU["Z_2e_0_qT30"] = 0.301;
	sigma_0["Z_2e_0_qT30"] = 3.895;
	sigma_PU["Z_2e_0_qT30"] = 1.271;

	// muons, fixed tails
	tail["Z_2m_1_qT30"] = 0.534;
	alpha_0["Z_2m_1_qT30"] = 4.763;
	alpha_PU["Z_2m_1_qT30"] = 0.320;
	sigma_0["Z_2m_1_qT30"] = 3.960;
	sigma_PU["Z_2m_1_qT30"] = 1.331;

	// leptons, fixed tails
	tail["Z_2l_2_qT30"] = 0.534;
	alpha_0["Z_2l_2_qT30"] = 4.743;
	alpha_PU["Z_2l_2_qT30"] = 0.304;
	sigma_0["Z_2l_2_qT30"] = 3.946;
	sigma_PU["Z_2l_2_qT30"] = 1.271;

      }
      {

	// gamma+jets qT>30 no jet veto, fixed tails
	tail["data_4_qT30"] = 0.460;
	alpha_0["data_4_qT30"] = 6.817;
	alpha_PU["data_4_qT30"] = 0.283;
	sigma_0["data_4_qT30"] = 5.420;
	sigma_PU["data_4_qT30"] = 1.354;



	// gamma+jets qT>30 with jet veto (25 GeV), fixed tails
	tail["data_3_qT30"] = 0.400;
	alpha_0["data_3_qT30"] = 7.040;
	alpha_PU["data_3_qT30"] = 0.431;
	sigma_0["data_3_qT30"] = 5.259;
	sigma_PU["data_3_qT30"] = 1.608;

      }
    }
  else
    {
      {
	//
	// jet veto 30 GeV
	//

	// fits with floating tail


	// data

	// electrons
	tail["data_0_j30veto"] = 0.418;
	alpha_0["data_0_j30veto"] = 4.996;
	alpha_PU["data_0_j30veto"] = 0.201;
	sigma_0["data_0_j30veto"] = 4.191;
	sigma_PU["data_0_j30veto"] = 0.475;

	// muons
	tail["data_1_j30veto"] = 0.426;
	alpha_0["data_1_j30veto"] = 4.870;
	alpha_PU["data_1_j30veto"] = 0.256;
	sigma_0["data_1_j30veto"] = 3.803;
	sigma_PU["data_1_j30veto"] = 1.173;

	// leptons
	tail["data_2_j30veto"] = 0.429;
	alpha_0["data_2_j30veto"] = 4.842;
	alpha_PU["data_2_j30veto"] = 0.268;
	sigma_0["data_2_j30veto"] = 3.748;
	sigma_PU["data_2_j30veto"] = 1.223;

	// DY MC

	// electrons
	tail["Z_2e_0_j30veto"] = 0.419;
	alpha_0["Z_2e_0_j30veto"] = 3.834;
	alpha_PU["Z_2e_0_j30veto"] = 0.312;
	sigma_0["Z_2e_0_j30veto"] = 2.934;
	sigma_PU["Z_2e_0_j30veto"] = 1.274;

	// muons
	tail["Z_2m_1_j30veto"] = 0.423;
	alpha_0["Z_2m_1_j30veto"] = 3.749;
	alpha_PU["Z_2m_1_j30veto"] = 0.313;
	sigma_0["Z_2m_1_j30veto"] = 2.984;
	sigma_PU["Z_2m_1_j30veto"] = 1.248;

	// leptons
	tail["Z_2l_2_j30veto"] = 0.442;
	alpha_0["Z_2l_2_j30veto"] = 3.762;
	alpha_PU["Z_2l_2_j30veto"] = 0.329;
	sigma_0["Z_2l_2_j30veto"] = 2.984;
	sigma_PU["Z_2l_2_j30veto"] = 1.253;

      }
      {
	//
	// jet veto 25 GeV
	//

	// fits with floating tail

	// data

	// electrons
	tail["data_0_j25veto"] = 0.413;
	alpha_0["data_0_j25veto"] = 4.962;
	alpha_PU["data_0_j25veto"] = 0.234;
	sigma_0["data_0_j25veto"] = 3.847;
	sigma_PU["data_0_j25veto"] = 1.066;

	// muons
	tail["data_1_j25veto"] = 0.417;
	alpha_0["data_1_j25veto"] = 4.914;
	alpha_PU["data_1_j25veto"] = 0.231;
	sigma_0["data_1_j25veto"] = 3.849;
	sigma_PU["data_1_j25veto"] = 1.069;

	// leptons
	tail["data_2_j25veto"] = 0.425;
	alpha_0["data_2_j25veto"] = 4.870;
	alpha_PU["data_2_j25veto"] = 0.248;
	sigma_0["data_2_j25veto"] = 3.794;
	sigma_PU["data_2_j25veto"] = 1.125;

	// DY simulation

	// electrons
	tail["Z_2e_0_j25veto"] = 0.408;
	alpha_0["Z_2e_0_j25veto"] = 3.921;
	alpha_PU["Z_2e_0_j25veto"] = 0.290;
	sigma_0["Z_2e_0_j25veto"] = 2.952;
	sigma_PU["Z_2e_0_j25veto"] = 1.207;

	// muons
	tail["Z_2m_1_j25veto"] = 0.431;
	alpha_0["Z_2m_1_j25veto"] = 3.848;
	alpha_PU["Z_2m_1_j25veto"] = 0.290;
	sigma_0["Z_2m_1_j25veto"] = 2.927;
	sigma_PU["Z_2m_1_j25veto"] = 1.216;

	// leptons
	tail["Z_2l_2_j25veto"] = 0.422;
	alpha_0["Z_2l_2_j25veto"] = 3.738;
	alpha_PU["Z_2l_2_j25veto"] = 0.319;
	sigma_0["Z_2l_2_j25veto"] = 2.907;
	sigma_PU["Z_2l_2_j25veto"] = 1.222;
      }


      {
	// qT>30 and jet veto 25 GeV

	// data

	// electrons
	tail["data_0_qT30_j25veto"] = 0.443;
	alpha_0["data_0_qT30_j25veto"] = 6.380;
	alpha_PU["data_0_qT30_j25veto"] = 0.114;
	sigma_0["data_0_qT30_j25veto"] = 4.660;
	sigma_PU["data_0_qT30_j25veto"] = 1.208;

	// muons
	tail["data_1_qT30_j25veto"] = 0.418;
	alpha_0["data_1_qT30_j25veto"] = 6.151;
	alpha_PU["data_1_qT30_j25veto"] = 0.327;
	sigma_0["data_1_qT30_j25veto"] = 4.530;
	sigma_PU["data_1_qT30_j25veto"] = 1.499;

	// leptons
	tail["data_2_qT30_j25veto"] = 0.424;
	alpha_0["data_2_qT30_j25veto"] = 6.188;
	alpha_PU["data_2_qT30_j25veto"] = 0.317;
	sigma_0["data_2_qT30_j25veto"] = 4.554;
	sigma_PU["data_2_qT30_j25veto"] = 1.492;

	// DY simulation

	// electrons
	tail["Z_2e_0_qT30"] = 0.466;
	alpha_0["Z_2e_0_qT30"] = 5.053;
	alpha_PU["Z_2e_0_qT30"] = 0.368;
	sigma_0["Z_2e_0_qT30"] = 4.098;
	sigma_PU["Z_2e_0_qT30"] = 1.394;

	// muons
	tail["Z_2m_1_qT30"] = 0.494;
	alpha_0["Z_2m_1_qT30"] = 4.828;
	alpha_PU["Z_2m_1_qT30"] = 0.422;
	sigma_0["Z_2m_1_qT30"] = 4.032;
	sigma_PU["Z_2m_1_qT30"] = 1.489;

	// leptons
	tail["Z_2l_2_qT30"] = 0.489;
	alpha_0["Z_2l_2_qT30"] = 4.921;
	alpha_PU["Z_2l_2_qT30"] = 0.403;
	sigma_0["Z_2l_2_qT30"] = 3.967;
	sigma_PU["Z_2l_2_qT30"] = 1.513;

      }

      {
	// qT>30 and no jet veto

	// data

	// electrons, floating tail
	tail["data_0_qT30"] = 0.543;
	alpha_0["data_0_qT30"] = 5.841;
	alpha_PU["data_0_qT30"] = 0.382;
	sigma_0["data_0_qT30"] = 5.083;
	sigma_PU["data_0_qT30"] = 1.295;

	// muons, floating tails
	tail["data_1_qT30"] = 0.535;
	alpha_0["data_1_qT30"] = 6.270;
	alpha_PU["data_1_qT30"] = 0.268;
	sigma_0["data_1_qT30"] = 5.374;
	sigma_PU["data_1_qT30"] = 1.197;

	// lepton, floating tails
	tail["data_2_qT30"] = 0.534;
	alpha_0["data_2_qT30"] = 6.203;
	alpha_PU["data_2_qT30"] = 0.297;
	sigma_0["data_2_qT30"] = 5.261;
	sigma_PU["data_2_qT30"] = 1.257;


	// DY simulation

	// electrons, floating tails
	tail["Z_2e_0_qT30"] = 0.466;
	alpha_0["Z_2e_0_qT30"] = 5.053;
	alpha_PU["Z_2e_0_qT30"] = 0.368;
	sigma_0["Z_2e_0_qT30"] = 4.098;
	sigma_PU["Z_2e_0_qT30"] = 1.394;

	// muons, floating tail
	tail["Z_2m_1_qT30"] = 0.507;
	alpha_0["Z_2m_1_qT30"] = 4.801;
	alpha_PU["Z_2m_1_qT30"] = 0.418;
	sigma_0["Z_2m_1_qT30"] = 3.949;
	sigma_PU["Z_2m_1_qT30"] = 1.489;

	// leptons, floating tails
	tail["Z_2l_2_qT30"] = 0.492;
	alpha_0["Z_2l_2_qT30"] = 4.932;
	alpha_PU["Z_2l_2_qT30"] = 0.409;
	sigma_0["Z_2l_2_qT30"] = 4.058;
	sigma_PU["Z_2l_2_qT30"] = 1.467;


      }

      {
	// gamma+jets qT>30 no jet veto, floating tails
	tail["data_4_qT30"] = 0.463;
	alpha_0["data_4_qT30"] = 6.774;
	alpha_PU["data_4_qT30"] = 0.344;
	sigma_0["data_4_qT30"] = 5.375;
	sigma_PU["data_4_qT30"] = 1.489;

	// gamma+jets qT>30 with jet veto (25 GeV), floating tails
	tail["data_3_qT30"] = 0.250;
	alpha_0["data_3_qT30"] = 7.980;
	alpha_PU["data_3_qT30"] = 0.243;
	sigma_0["data_3_qT30"] = 5.576;
	sigma_PU["data_3_qT30"] = 1.450;
      }

    }

  // Dynamical K factors
  setVVKFactors();

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++
  //+++++++++   end of initialisation    ++++++++++++++++
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++

}

