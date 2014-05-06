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

#include "Analysis/ghmzz/ZZ_2l2n_Analysis.hh"
#include "Analysis/utils/CutUtils.hh"
#include "Analysis/utils/Config.hh"

#include "Analysis/utils/Constants.hh"
#include "Analysis/utils/Config.hh"
#include "Analysis/utils/KineUtils.hh"

ClassImp( ZZ_2l2n_Analysis )

bool verbose = true;

void 
ZZ_2l2n_Analysis::analyze_( float lumi )
{    
  //
  // loop on the samples
  //
  if( _ZJets_cs ) setMomentumWeights();

  string sampleAnalysis("Diboson");
  string tname("ZZtuple");
  string mc_tname("VV_mcTruth");

  // list of samples samples
  //  Sample::SampleList sampleList;

  typedef multimap< string, Sample* > SampleMap;
  typedef multimap< string, Sample* >::iterator SampleIt;
  typedef pair< string, Sample* > SamplePair;
  SampleMap sampleList;
  vector<string> ZZanal;

  //   if( dummy )
  //     {
  //       // the signal
  //       sampleList.push_back(samples["Z_2e"]);
  //       sampleList.push_back(samples["ZJets_2e"]);
  //       //      sampleList.push_back(samples["WZ_3ln"]);
  //       //      sampleList.push_back(samples["WW_2l2n"]);
  //       verbose = true;
  //     }
  //   else 

  if( _signalOnly )
    {
      if( imode==0 || imode==2 )
	// electron analysis
	{
	  ZZanal.push_back("ZZ_2e2n");
	}
      if( imode==1 || imode==2 )
	// electron analysis
	{
	  ZZanal.push_back("ZZ_2m2n");
	}
        
      // ZZ
      sampleList.insert(SamplePair("ZZ_2l2n",samples["ZZ_2l2n"]));

      // WW
      sampleList.insert(SamplePair("WW_2l2n",samples["WW_2l2n"]));
  
      // WZ
      sampleList.insert(SamplePair("WZ_3ln",samples["WZ_3ln"]));      
    }
  else
    {
      if( imode==0 || imode==2 )
	// electron analysis
	{
	  ZZanal.push_back("ZZ_2e2n");
	  if( period=="2010" )
	    {
	      sampleList.insert(SamplePair("data",samples["Electron_1"]));
	      sampleList.insert(SamplePair("data",samples["Electron_2"]));
	      sampleList.insert(SamplePair("data",samples["Electron_3"]));
	      sampleList.insert(SamplePair("data",samples["Electron_4"]));
	      sampleList.insert(SamplePair("data",samples["Electron_5"]));
	      sampleList.insert(SamplePair("data",samples["Electron_6"]));
	      sampleList.insert(SamplePair("data",samples["Electron_7"]));
	    }
	  else if( period=="2011A" )
	    {
	      sampleList.insert(SamplePair("data",samples["DoubleEle_1"]));	  
	      sampleList.insert(SamplePair("data",samples["DoubleEle_2"]));	  
	      sampleList.insert(SamplePair("data",samples["DoubleEle_3"]));
	      sampleList.insert(SamplePair("data",samples["DoubleEle_4"]));
	      sampleList.insert(SamplePair("data",samples["DoubleEle_5"]));
	      sampleList.insert(SamplePair("data",samples["DoubleEle_6"]));
	      sampleList.insert(SamplePair("data",samples["DoubleEle_7"]));
	      sampleList.insert(SamplePair("data",samples["DoubleEle_8"]));
	      sampleList.insert(SamplePair("data",samples["DoubleEle_9"]));
	      sampleList.insert(SamplePair("data",samples["DoubleEle_10"]));
	      sampleList.insert(SamplePair("data",samples["DoubleEle_11"]));
	      sampleList.insert(SamplePair("data",samples["DoubleEle_12"]));
	      sampleList.insert(SamplePair("data",samples["DoubleEle_13"]));
	      sampleList.insert(SamplePair("data",samples["DoubleEle_14"]));
	      sampleList.insert(SamplePair("data",samples["DoubleEle_15"]));
	      sampleList.insert(SamplePair("data",samples["DoubleEle_16"]));
	    }

	  // DY

	}
      if( imode==1 || imode==2 )
	// muon analysis
	{
	  ZZanal.push_back("ZZ_2m2n");
	  if( period=="2010" )
	    {
	      sampleList.insert(SamplePair("data",samples["Muon_1"]));
	      sampleList.insert(SamplePair("data",samples["Muon_2"]));
	      sampleList.insert(SamplePair("data",samples["Muon_3"]));
	      sampleList.insert(SamplePair("data",samples["Muon_4"]));
	      sampleList.insert(SamplePair("data",samples["Muon_5"]));
	      sampleList.insert(SamplePair("data",samples["Muon_6"]));
	      sampleList.insert(SamplePair("data",samples["Muon_7"]));
	    }
	  else if( period=="2011A" )
	    {
	      sampleList.insert(SamplePair("data",samples["DoubleMu_1"]));
	      sampleList.insert(SamplePair("data",samples["DoubleMu_2"]));
	      sampleList.insert(SamplePair("data",samples["DoubleMu_3"]));
	      sampleList.insert(SamplePair("data",samples["DoubleMu_4"]));
	      sampleList.insert(SamplePair("data",samples["DoubleMu_5"]));
	      sampleList.insert(SamplePair("data",samples["DoubleMu_6"]));
	      sampleList.insert(SamplePair("data",samples["DoubleMu_7"]));
	      sampleList.insert(SamplePair("data",samples["DoubleMu_8"]));
	      sampleList.insert(SamplePair("data",samples["DoubleMu_9"]));
	      sampleList.insert(SamplePair("data",samples["DoubleMu_10"]));
	      sampleList.insert(SamplePair("data",samples["DoubleMu_11"]));
	      sampleList.insert(SamplePair("data",samples["DoubleMu_12"]));
	      sampleList.insert(SamplePair("data",samples["DoubleMu_13"]));
	      sampleList.insert(SamplePair("data",samples["DoubleMu_14"]));
	      sampleList.insert(SamplePair("data",samples["DoubleMu_15"]));
	      sampleList.insert(SamplePair("data",samples["DoubleMu_16"]));
	    }

	  // DY
	}
      if( imode==0 )
	{
	  sampleList.insert(SamplePair("Z_2e",samples["Z_2e"]));
	}
      else if( imode==1 )
	{
	  sampleList.insert(SamplePair("Z_2m",samples["Z_2m"]));
	}
      else if( imode==2 )
	{
	  sampleList.insert(SamplePair("Z_2l",samples["Z_2e"]));
	  sampleList.insert(SamplePair("Z_2l",samples["Z_2m"]));
	}



      // Z to tau-tau
      sampleList.insert(SamplePair("Z_2t",samples["Z_2t"]));

      // top
      sampleList.insert(SamplePair("tt",samples["ttbar"]));
      sampleList.insert(SamplePair("t",samples["t_schan"]));
      sampleList.insert(SamplePair("t",samples["t_tchan"]));
      sampleList.insert(SamplePair("t",samples["tW"]));
  
      // ZZ
      sampleList.insert(SamplePair("ZZ_2l2n",samples["ZZ_2l2n"]));

      // WW
      sampleList.insert(SamplePair("WW_2l2n",samples["WW_2l2n"]));
  
      // WZ
      sampleList.insert(SamplePair("WZ_3ln",samples["WZ_3ln"]));      
  
      // VGam
      sampleList.insert(SamplePair("VGam",samples["VGam"]));      
  
    }

  cout << "Start loop on samples " << endl;
  
  ofstream ifile;
  string fsel(analysis);
  fsel += ".sel";
  ifile.open(fsel.c_str());

  //
  // Cuts
  //
  bool presel=true;
  OneSidedCut<bool> c_presel   ( "presel",      presel, true, "==" );

  // Z quality
  OneSidedCut<int>   c_nZCand  ( "nZCand",      nZCand,    1,   "==" );
  TwoSidedCut<float> c_ZWindow ( "ZWindow",     mll, 60., 120., "[]" );
  OneSidedCut<bool>  c_ZLQual  ( "ZLQual",     isTight, true,   "==" );
  OneSidedCut<bool>  c_ZLLTight ( "ZLLTight",   isLLTight, true,   "==" );
  CompositeCut        c_ZQual   ( "ZQual", "&&" );
  c_ZQual.add( c_nZCand );
  c_ZQual.add( c_ZWindow );
  c_ZQual.add( c_ZLQual );
  c_ZQual.add( c_ZLLTight );

  // Z qT
  OneSidedCut<float> c_qT30 ("qT30", qTll, 30., ">=" );

  // jet veto
  OneSidedCut<int>   c_j30Veto ( "j30Veto",    jMult30, 0,      "==" );
  OneSidedCut<int>   c_j30One  ( "j30One",     jMult30, 1,      "==" );

  // MET cut
  OneSidedCut<float> c_ghmMET50 ( "ghmMET50",  ghmMET,  50.,   ">"  );
  OneSidedCut<float> c_ghmMET60 ( "ghmMET60",  ghmMET,  60.,   ">"  );

  // balance
  TwoSidedCut<float> c_bal  ( "bal",    bal, 0.4, 1.8, "[]" ); //!!! GHM

  // jet/MET quality
  TwoSidedCut<float> c_dPhiJMET( "dPhiJMET", dPhiJMET, -20.,20, "![]");

  // Z/MET quality
  TwoSidedCut<float> c_dPhiZMET ( "dPhiZMET", dPhiZMET, -60, 60, "[]" ); //!!! GHM

  // B-tagging
  OneSidedCut<float> c_BtagTkCVeto   ( "BtagTkCVeto",   jBtagTkC_max, 2.5, "<=" );
  OneSidedCut<float> c_BtagSoftMVeto ( "BtagSoftMVeto", jBtagSoftM_max, 0.30, "<=" );
  CompositeCut c_BtagVeto( "BTagVeto", "&&" );
  c_BtagVeto.add( c_BtagTkCVeto    );
  c_BtagVeto.add( c_BtagSoftMVeto  );	  

  // diboson reduction
  OneSidedCut<int>   c_LVeto   ( "LVeto",       nLoose, 0,      "==" );
  CompositeCut c_dibosonVeto( "dib", "&&" );
  c_dibosonVeto.add( c_LVeto );

  TwoSidedCut<float> c_massCut  ( "massCut", mll, 80., 100, "[]");

  // final selection
  CompositeCut c_finalSelection( "final", "&&" );
  c_finalSelection.add( c_ZQual       );
  c_finalSelection.add( c_bal         );
  c_finalSelection.add( c_j30Veto     );
  c_finalSelection.add( c_ghmMET50    ); 

  // show content:
  for ( SampleIt it=sampleList.begin() ; it != sampleList.end(); it++ )
    {
      string sampleDir = it->first;  // defines the histogram directory
      Sample* sample_ = it->second;

      cout << "Histo Directory = " << sampleDir << endl;

      //  for( size_t isample=0; isample<sampleList.size(); isample++ )
      //    {

      sample_->print( cout );
      cout << "n_tot(" << lumi << ")=" << sample_->n_tot( lumi ) << endl; 
      cout << "weight(" << lumi << ")=" << sample_->w( lumi ) << endl; 

      float sampleWeight = sample_->w( lumi );

      string sampleName = sample_->name();
      size_t pos = sampleName.find('/');  
      string subsample = sampleName.substr(pos+1);
      string sampleRoot = sampleName.substr(0,pos);
      
      isZZ = sampleRoot.find("ZZ")!=string::npos;
      isWZ = sampleRoot.find("WZ")!=string::npos;
      isWW = sampleRoot.find("WW")!=string::npos;

      // create the sample weight
      VertexWeight vw( sampleRoot, period );      

      string filename = Config::rootPath;
      filename += sampleAnalysis;
      filename += "/";
      filename += sampleRoot;
      filename += ".root";
            
      if( subsample=="data" )
	ifile << sampleName << endl;

      //      string dir(subsample);
      //  string fname("ZZ_2l2n.root");
      //  string dir("ZZ_2e2n");
      TFile* inputFile  = new TFile( filename.c_str(), "READ" );

      // tuple manager
      TupleManager tm;
      tm.setFile( inputFile, TupleManager::kRead );

      int nsel=0;
      for( size_t ianal=0; ianal<ZZanal.size(); ianal++ )
	{
	  isMuMu = ( ZZanal[ianal]=="ZZ_2m2n" );
	  isEE   = ( ZZanal[ianal]=="ZZ_2e2n" );
	  assert( isMuMu || isEE );

	  //	  if( sampleRoot!="ttbar" ) continue;
	  isDYSimulation = 
	    (!_ZJets_cs) &&
	    (sampleRoot=="Z_2l" || sampleRoot=="Z_2e" || sampleRoot=="Z_2m");
	  isDibosonSimulation = 
	    (sampleRoot=="WW" || sampleRoot=="WZ" || sampleRoot=="ZZ" ||
	     sampleRoot=="WW_2l2n" || sampleRoot=="WZ_3ln" || 
	     sampleRoot=="ZZ_2l2n" );	  
	  isData = (subsample=="data");
	  isOtherSimulation = !(isData || isDYSimulation || isDibosonSimulation );
	 
	  // tree
	  TTree* tree = tm.getTree( tname, ZZanal[ianal] );
	  if( tree==0 ) continue;

	  TBits   c_bits;
	  TBits*  c_bits_ptr = &c_bits;
	  string  c_str;
	  string* c_str_ptr  = &c_str;
	  string  categ;
	  string* categ_ptr = & categ;

 	  tm.add< int >(    "run",   &run       );
 	  tm.add< int >(    "event", &event     );
  	  tm.add< string*>( "categ", &categ_ptr );
  	  tm.add< string*>( "c_str", &c_str_ptr  );      
  	  tm.add< TBits* >( "c_bits", &c_bits_ptr  );      
  	  tm.add< int >(  "nZCand",   &nZCand        );
  	  tm.add< int >(  "nVertex",   &nVertex      );
  	  tm.add< bool  >(  "isTight",&isTight    );
  	  tm.add< float >(  "mll",   &mll        );
  	  tm.add< float >(  "qTll",   &qTll        );
  	  tm.add< float >(  "phill",  &phill        );
  	  tm.add< float >(  "yll",   &yll        );
  	  tm.add< float >(  "hll",   &hll        );
  	  tm.add< int >(  "nLoose",&nLoose    );

	  tm.add< vectorFloat* >(  "LPt",   &LPt   );
	  tm.add< vector<int>* >(  "LPdg",   &LPdg   );
	  tm.add< vectorFloat* >(  "LEta",  &LEta  );
	  tm.add< vectorFloat* >(  "LPhi",  &LPhi  );
	  tm.add< vector<int>* >(  "LID",   &LID   );
	  tm.add< vectorFloat* >(  "LMT",   &LMT   );
	  tm.add< vectorFloat* >( "LTrkIso", &LTrkIso );
	  tm.add< vectorFloat* >( "LEcalIso", &LEcalIso );
	  tm.add< vectorFloat* >( "LHcalIso", &LHcalIso );

 	  tm.add< float >(  "MET",   &MET        );
 	  tm.add< float >(  "METPhi",&METPhi    );
 	  tm.add< float >(  "sumEt", &sumEt );
 	  tm.add< float >(  "mEtSig", &mEtSig );
 	  tm.add< float >(  "mTZZ", &mTZZ );
 	  tm.add< float >(  "projMET",&projMET    );
 	  tm.add< float >(  "corProjMET",&corProjMET    );
 	  tm.add< float >(  "sigMET",&sigMET    );
 	  tm.add< float >(  "cpuMETSum",&cpuMETSum    );
 	  tm.add< float >(  "cpuMETSqrt",&cpuMETSqrt    );
 	  tm.add< float >(  "cpuMETVtx",&cpuMETVtx    );
 	  tm.add< float >(  "dPhiMin",&dPhiMin    );
 	  tm.add< float >(  "balZMET",&balZMET    );


	  tm.add< vectorFloat* > ( "jEt", &jEt );
	  tm.add< vectorFloat* > ( "jEta", &jEta );
	  tm.add< vectorFloat* > ( "jPhi", &jPhi );
	  tm.add< vectorFloat* > ( "jDR1", &jDR1 );
	  tm.add< vectorFloat* > ( "jDR2", &jDR2 );
	  tm.add< vectorFloat* > ( "jDPhiMET", &jDPhiMET );
	  tm.add< vectorFloat* > ( "jBtagTkC", &jBtagTkC );
	  tm.add< vectorFloat* > ( "jBtagSoftM", &jBtagSoftM );

 	  tm.add< int >(  "jMult40",&jMult40  );
 	  tm.add< int >(  "jMult35",&jMult35  );
 	  tm.add< int >(  "jMult30",&jMult30  );
 	  tm.add< int >(  "jMult25",&jMult25  );
 	  tm.add< int >(  "jMult20",&jMult20  );
 	  tm.add< int >(  "jMult15",&jMult15  );
 	  tm.add< int >(  "jMult10",&jMult10  );
 	  tm.add< float >(  "JZB25",&JZB25  );
 	  tm.add< float >(  "JZB20",&JZB20  );
 	  tm.add< float >(  "JZB15",&JZB15  );
 	  tm.add< float >(  "dPhiJMET",&dPhiJMET );	  

	  int n = tree->GetEntriesFast();

	  int nVtx_max(20);
	  bool vtxEq[20];

	  TTree* mc_tree = tm.getTree( mc_tname, ZZanal[ianal] );
	  int n_mc = 0;

	  if( mc_tree )
	    {
	      n_mc = mc_tree->GetEntriesFast();

	      tm.add< int   >(    "run",   &mc_run   );
	      tm.add< int   >(    "event", &mc_event );
	      tm.add< string*>( "categ", &categ_ptr );
	      tm.add< float >(    "mVV",   &mc_mVV   );
	      tm.add< float >(    "mV1",   &mc_mV1   );
	      tm.add< float >(    "ptV1",  &mc_ptV1  );
	      tm.add< float >(    "pzV1",  &mc_pzV1  );
	      tm.add< float >(    "EV1",   &mc_EV1   );
	      tm.add< float >(    "yV1",   &mc_yV1   );
	      tm.add< float >(    "mV2",   &mc_mV2   );
	      tm.add< float >(    "ptV2",  &mc_ptV2  );
	      tm.add< float >(    "pzV2",  &mc_pzV2  );
	      tm.add< float >(    "EV2",   &mc_EV2   );
	      tm.add< float >(    "yV2",   &mc_yV2   );
	      tm.add< vector<int>* > ( "pdgL", &mc_pdgL );
	      tm.add< vectorFloat* > ( "ptL",  &mc_ptL  );
	      tm.add< vectorFloat* > ( "etaL", &mc_etaL );
	      tm.add< vectorFloat* > ( "phiL", &mc_phiL );	      
	    }

	  // maps 
	  //      map< string, bool >  sel;  // cut results
	  //      map< AbsCut*, bool > cutres;  // cut results
	  vector< HBool > hbool;

	  map< string, bool > sel;
	  map< string, bool > sel_ntu;

	  int ievent(0);
	  int ievent_0(0);
	  int currun(-1);
	  int curevt(-1);
	  for( int ii=0; ii<n; ii++ )
	    {
	      tree->GetEntry( ii );

	      eventWeight = vw( nVertex );		  		  
	      aTGCfactor   = 1;
	      kfactor      = 1;
	      kfactor_jv30 = 1;

	      if( isDibosonSimulation && mc_tree!=0 )
		{
		  mc_tree->GetEntry( ii );
		  assert( mc_run==run && mc_event==event);

		  analyzeSignalMC();

		  if( _ZJets_cs && !isInAcceptance ) continue; 
		  
		}
	      
	      // electron channel: determine if both leptons are well identified
	      isLLTight = true;
	      if( isEE )
		{
		  for( size_t ij=0; ij<2; ij++ )
		    {
		      float LPt_   = (*LPt)[ij];
		      float LEta_  = (*LEta)[ij];
		      float LPhi_  = (*LPhi)[ij];
		      float LPdg_  = (*LPdg)[ij];
		      int   LID_   = (*LID)[ij];
		      TVector2 L2D_( LPt_*cos(LPhi_), LPt_*sin(LPhi_) );
		      if( LID_<7 ) isLLTight=false;
		    }		  
		}

	      if( !c_ZQual() ) continue; 
	      nsel++;

	      eventWeight *= kfactor_jv30;
	      //eventWeight *= kfactor;
	      eventWeight *= aTGCfactor;

	 
	      ievent_0++;
	      if( run!=currun || event!=curevt )
		{
		  ievent++;
		  currun = run;
		  curevt = event;
		}

	      // fill the vertex booleans
	      bool addCut = c_j30Veto();  // additional cut for vertex studies
	      addCut = addCut && c_qT30();  
	      for( int iv=0; iv<nVtx_max; iv++ )
		{		  
		  vtxEq[iv]=false;
		  if( addCut && nVertex==iv ) vtxEq[iv]=true;
		}
	      if( nVertex>nVtx_max-1 ) vtxEq[nVtx_max-1]=true;



	      string sigStr_;
	      string refStr_;
	      isLL = false; //!!!
	      if( isLL )
		{
		  refStr_ = "data_2_j30veto";
		}
	      else if( isEE )
		{
		  refStr_ = "data_0_j30veto";
		}
	      else if( isMuMu )
		{
		  refStr_ = "data_1_j30veto";
		}
	      //	      if( isDYSimulation )
	      if( !isData )
		{
		  if( isLL )
		    {
		      sigStr_ = "Z_2l_2_j30veto";
		    }
		  else if( isEE )
		    {
		      sigStr_ = "Z_2e_0_j30veto";
		    }
		  else if( isMuMu )
		    {
		      sigStr_ = "Z_2m_1_j30veto";
		    }
		}
	      // if( isData || isDibosonSimulation || isOtherSimulation )
	      else
		{
		  if( isLL )
		    {
		      sigStr_ = "data_2_j30veto";		      
		    }
		  else if( isEE )
		    {
		      sigStr_ = "data_0_j30veto";
		    }
		  else if( isMuMu )
		    {
		      sigStr_ = "data_1_j30veto";
		    }
		}
// 	      else
// 		{
// 		  assert(0);
// 		}

	      float alpha_0_   = alpha_0 [sigStr_];
	      float alpha_PU_  = alpha_PU[sigStr_];
	      float sigma_0_   = sigma_0 [sigStr_];
	      float sigma_PU_  = sigma_PU[sigStr_];
	      float alpha_ref_ = alpha_0 [refStr_];
	      float sigma_ref_ = sigma_0 [refStr_];
	      
	      TVector2 sumJet2D;
	      TVector2 sumJet2D_25;

	      float METPhi_ = METPhi; 
	      TVector2 MET2D( MET*cos( METPhi_), MET*sin( METPhi_ ) );
	      TVector2 METVtx2D( cpuMETVtx*cos( METPhi_), cpuMETVtx*sin( METPhi_ ) ); // !!!  should be another phi !!!
	      TVector2 ll2D( qTll*cos(phill), qTll*sin(phill) );	     
	      
	      TVector2 llMET2D = ll2D + MET2D;

	      jEt_1 = -1;
	      jEt_2 = -1;
	      jEt_3 = -1;
	      jEt_4 = -1;

	      jBtagTkC_max   = -100.;
	      jBtagSoftM_max = -100.;

	      bool ok_ = true;
	      if( _ZJets_cs && !isDibosonSimulation ) 
		{
		  if( !c_j30One() ) continue;
		  ok_ = false;
		}
	      int ji=-1;
	      for( size_t ij=0; ij<jEt->size(); ij++ )
		{
		  ji++;
		  float jPhi_ = (*jPhi)[ij];
		  float jEt_  = (*jEt)[ij];
		  TVector2 j2D_( jEt_*cos(jPhi_), jEt_*sin(jPhi_) );
		  if( _ZJets_cs && !isDibosonSimulation )
		    //if( _ZJets_cs )
		    {		      
		      if( ij==0  )
			{
			  assert( jEt_>=30. );
			  //			  if( jEt_<30. )	       
			  //			    {
			  //			      break;
			  //			    }
			  //			  else
			    {
			      ok_ = true;
			      MET2D += j2D_;
			      METVtx2D += j2D_;
			      MET = MET2D.Mod();
			      //			  cpuMETSum = MET;  //!!!
			      cpuMETVtx = METVtx2D.Mod();
			      //			  cpuMETSqrt = MET;
			      jMult30--;
			      jMult25--;
			      jMult20--;
			      jMult15--;
			      jMult10--;
			      ji--;
			      
			      eventWeight *= (*momWeights[sampleDir])(qTll);
			      
			      continue;
			    }
			}
		    }

		  if( ji==0 ) jEt_1 = jEt_;
		  if( ji==1 ) jEt_2 = jEt_;
		  if( ji==2 ) jEt_3 = jEt_;
		  if( ji==3 ) jEt_4 = jEt_;
		  
		  if( (*jBtagTkC)[ij]>jBtagTkC_max ) 
		    jBtagTkC_max = (*jBtagTkC)[ij];
		  if( (*jBtagSoftM)[ij]>jBtagSoftM_max ) 
		    jBtagSoftM_max = (*jBtagSoftM)[ij];

		  sumJet2D += j2D_;
		  if( jEt_>=25 ) 
		    sumJet2D_25 += j2D_;		    
		}
	      if( !ok_ ) continue;
	      
	 	  
	      float cpuMET = cpuMETVtx;
	      float cpuMETMin = MET;
	      if( cpuMET<MET ) cpuMETMin = cpuMET;
	      float theMET = cpuMETMin;	      
	      float sigma_T  = 
		sqrt( pow(sigma_0_,2) + (nVertex-1)*pow(sigma_PU_,2) );
	      float peak_T   = alpha_0_ + (nVertex-1)*alpha_PU_;
	      ghmSig = ( theMET - peak_T ) / sigma_T;
	      ghmMET = sigma_ref_ * ghmSig + alpha_ref_;

	      // angle between the MET and the sum of all the jets with ET>10
	      TVector2 minusMET2D = MET2D;
	      minusMET2D *= -1;
	      dPhiSumJMET = KineUtils::dPhi( sumJet2D.Phi(), minusMET2D.Phi() )
		* Constants::radToDeg;	      
	      dPhiZMET    = KineUtils::dPhi( ll2D.Phi(), minusMET2D.Phi() )
		* Constants::radToDeg;	      

	      //	      TVector2 trueMET2D = ll2D + sumJet2D_25;
	      TVector2 trueMET2D = ll2D + sumJet2D;
	      trueMET2D *= -1;

	      TVector2 unclMET2D = MET2D - trueMET2D;


	      trueMET = trueMET2D.Mod();
	      unclMET = unclMET2D.Mod();
	      
	      // new definition of the balance
	      bal = MET / qTll;

	      // new definition of the transverse mass
	      float m_Z2 = pow( Constants::Z0_m, 2 );
	      float qT2  = pow( qTll, 2 );
	      float MET2 = pow( MET, 2 );
	      float mTZZ2 = 
		pow( sqrt( m_Z2 + qT2 ) + sqrt( m_Z2 + MET2 ), 2 )
		- 
		pow( llMET2D.Mod(), 2 );
	      mTZZ = 0;
	      if( mTZZ2>0 ) mTZZ = sqrt( mTZZ2 );
		
	      thirdL_MT = 0;
	      if( LPt->size()>2 )
		{
		  thirdL_MT = (*LMT)[2];
		}

	      //!!!
	      sel["ptV30"] = mc_ptV2>=30;

	      CompositeCut& cut = c_finalSelection;
	      sel[cut.name()] = cut();

	      bool b_q30 = c_qT30();
	      bool b_jet = c_j30Veto();  //!!! NEW     
	      bool b_m60 = c_ghmMET60();  //!!! NEW
	      bool b_m50 = c_ghmMET50();  //!!! NEW
	      bool b_met = b_m50; 
	      bool b_bal = c_bal();
	      bool b_phj = c_dPhiJMET();
	      bool b_phz = c_dPhiZMET();
	      bool b_btg = c_BtagVeto();
	      bool b_dib = c_dibosonVeto();
	      bool b_mas = c_massCut();
	      bool b_mtl = !b_dib && thirdL_MT>50;

	      sel["presel"] = true;
	      sel["jet"] = b_jet;

	      // special successive cuts
	      sel["q30"] 
		= b_q30;
	      sel["q30_jet"] 
		= b_q30 && b_jet;
	      sel["q30_jet_met"] 
		= b_q30 && b_jet && b_met;
	      sel["q30_jet_met_bal"] 
		= b_q30 && b_jet && b_met && b_bal;
	      sel["q30_jet_met_bal"] 
		= b_q30 && b_jet && b_met && b_bal;
	      sel["q30_jet_met_bal_phj"] 
		= b_q30 && b_jet && b_met && b_bal && b_phj;
	      sel["q30_jet_met_bal_phj_phz"] 
		= b_q30 && b_jet && b_met && b_bal && b_phj && b_phz;
	      sel["q30_jet_met_bal_phj_phz_btg"] 
		= b_q30 && b_jet && b_met && b_bal && b_phj && b_phz && b_btg;
	      sel["q30_jet_met_bal_phj_phz_btg_dib"] 
		= b_q30 && b_jet && b_met && b_bal && b_phj && b_phz && b_btg && b_dib;
	      sel["q30_jet_met_bal_phj_phz_btg_dib_mas"] 
		= b_q30 && b_jet && b_met && b_bal && b_phj && b_phz && b_btg && b_dib && b_mas;

	      sel["q30_dib"] 
		= b_q30 && b_dib;
	      sel["q30_dib_jet"] 
		= b_q30 && b_dib && b_jet;
	      sel["q30_dib_jet_met"] 
		= b_q30 && b_dib && b_jet && b_met;
	      sel["q30_dib_jet_met_bal"] 
		= b_q30 && b_dib && b_jet && b_met && b_bal;
	      sel["q30_dib_jet_met_bal"] 
		= b_q30 && b_dib && b_jet && b_met && b_bal;
	      sel["q30_dib_jet_met_bal_phj"] 
		= b_q30 && b_dib && b_jet && b_met && b_bal && b_phj;
	      sel["q30_dib_jet_met_bal_phj_phz"] 
		= b_q30 && b_dib && b_jet && b_met && b_bal && b_phj && b_phz;
	      sel["q30_dib_jet_met_bal_phj_phz_btg"] 
		= b_q30 && b_dib && b_jet && b_met && b_bal && b_phj && b_phz && b_btg;
	      sel["q30_dib_jet_met_bal_phj_phz_btg_mas"] 
		= b_q30 && b_dib && b_jet && b_met && b_bal && b_phj && b_phz && b_btg && b_mas;

	      // N-1 
	      sel["jet_met_bal_phj_phz_btg_dib_mas"] 
		=          b_jet && b_met && b_bal && b_phj && b_phz && b_btg && b_dib && b_mas;
	      sel["q30_met_bal_phj_phz_btg_dib_mas"] 
		= b_q30 &&          b_met && b_bal && b_phj && b_phz && b_btg && b_dib && b_mas;
	      sel["q30_jet_bal_phj_phz_btg_dib_mas"] 
		= b_q30 && b_jet          && b_bal && b_phj && b_phz && b_btg && b_dib && b_mas;
	      sel["q30_jet_met_phj_phz_btg_dib_mas"] 
		= b_q30 && b_jet && b_met          && b_phj && b_phz && b_btg && b_dib && b_mas;
	      sel["q30_jet_met_bal_phz_btg_dib_mas"] 
		= b_q30 && b_jet && b_met && b_bal          && b_phz && b_btg && b_dib && b_mas;
	      sel["q30_jet_met_bal_phj_btg_dib_mas"] 
		= b_q30 && b_jet && b_met && b_bal && b_phj          && b_btg && b_dib && b_mas;
	      sel["q30_jet_met_bal_phj_phz_dib_mas"] 
		= b_q30 && b_jet && b_met && b_bal && b_phj && b_phz          && b_dib && b_mas;
	      sel["q30_jet_met_bal_phj_phz_btg_mas"] 
		= b_q30 && b_jet && b_met && b_bal && b_phj && b_phz && b_btg          && b_mas;

	      // more 
	      sel["q30_met_bal_phj_phz_mas"] 
		= b_q30 && b_met && b_bal && b_phj && b_phz && b_mas;
	      sel["q30_met_bal_phj_phz_sam"] 
		= b_q30 && b_met && b_bal && b_phj && b_phz && !b_mas;
	      sel["q30_met_btg_mas"] 
		= b_q30 && b_met && b_btg && b_mas;
	      sel["q30_met_gtb_sam"] 
		= b_q30 && b_met && !b_btg && !b_mas;

	      // other combinations of cuts
	      sel["met"] =  b_met;
	      sel["bal"] =  b_bal;
	      sel["dib"] =  b_dib;

	      sel["q30_met"] =  b_q30 &&  b_met;

	      // anticuts
	      sel["tej"]                 =                             !b_jet;
	      sel["q30_tej"]             =  b_q30                   && !b_jet;
	      sel["met_bal_tej"]         =                    b_met && b_bal && !b_jet;
	      sel["q30_met_bal_tej"]     =  b_q30          && b_met && b_bal && !b_jet;
	      sel["bid"]                 =                                      !b_dib;
	      sel["bid"]                 =                                      b_mtl;
	      sel["jet_bid"]             =                             b_jet && !b_dib;
	      sel["jet_mtl"]             =                             b_jet && b_mtl;
	      sel["q30_jet_bid"]         =  b_q30 && b_jet          && !b_dib;
	      sel["jet_met_bid"]         =           b_jet && b_met          && !b_dib;
	      sel["q30_jet_met_bid"]     =  b_q30 && b_jet && b_met          && !b_dib;
	      sel["q30_jet_met_bal_bid"] =  b_q30 && b_jet && b_met && b_bal && !b_dib;
	      sel["q30_jet_met_bal_phj_phz_btg_bid"] 
		= b_q30 && b_jet && b_met && b_bal && b_phj && b_phz && b_btg          && !b_dib;
	      sel["q30_jet_phj_phz_btg_dib_mas"] 
		= b_q30 && b_jet          && b_phj && b_phz && b_btg && b_dib && b_mas;
	      sel["q30_jet_phz_btg_dib_mas"] 
		= b_q30 && b_jet         && b_phz && b_btg && b_dib && b_mas;
	      sel["q30_jet_dib_mas"] 
		= b_q30 && b_jet          && b_dib && b_mas;

	      if( hbool.size()==0 )
		{
		  string dir_;
	      
		  // fill the vertex booleans
		  for( int iv=0; iv<nVtx_max; iv++ )
		    {
		      //		  dir_ =  sampleName;
		      //		  if( subsample=="data" ) dir_ = "data";
		      dir_ = sampleDir;
		      dir_ += string("/");
		      string str_ = "vertex/vtx_";
		      ostringstream o;
		      o << str_ << iv;
		      dir_ += o.str();

		      hbool.push_back( 
				      HBool(
					    hist(
						 "Vmet", "MET", dir_, "" 
						 ),
					    &MET,
					    &(vtxEq[iv]), false, false
					    )  
				       );	  
		      hbool.push_back( 
				      HBool(
					    hist(
						 "Vmet", "projMET", dir_, "" 
						 ),
					    &projMET,
					    &(vtxEq[iv]), false, false
					    )  
				       );	  
		      hbool.push_back( 
				      HBool(
					    hist(
						 "Vmet", "ghmMET", dir_, "" 
						 ),
					    &ghmMET,
					    &(vtxEq[iv]), false, false
					    )  
				       );	  
		      hbool.push_back( 
				      HBool(
					    hist(
						 "n20", "jMult30", dir_, "" 
						 ),
					    &jMult30,
					    &(vtxEq[iv]), false, false
					    )  
				       );	  
		      hbool.push_back( 
				      HBool(
					    hist(
						 "n20", "jMult25", dir_, "" 
						 ),
					    &jMult25,
					    &(vtxEq[iv]), false, false
					    )  
				       );	  
		      hbool.push_back( 
				      HBool(
					    hist(
						 "n20", "jMult20", dir_, "" 
						 ),
					    &jMult20,
					    &(vtxEq[iv]), false, false
					    )  
				       );	  
		      hbool.push_back( 
				      HBool(
					    hist(
						 "n20", "jMult15", dir_, "" 
						 ),
					    &jMult15,
					    &(vtxEq[iv]), false, false
					    )  
				       );	  
		      hbool.push_back( 
				      HBool(
					    hist(
						 "n20", "jMult10", dir_, "" 
						 ),
					    &jMult10,
					    &(vtxEq[iv]), false, false
					    )  
				       );	  
		  
		      // 		  hbool.push_back( 
		      // 				  HBool(
		      // 					hist(
		      // 					     "Vmet", "MET_ge", dir_, "" 
		      // 					     ),
		      // 					&MET,
		      // 					&(vtxGe[iv]), false, false
		      // 					)  
		      // 				   );	  
		      // 		  hbool.push_back( 
		      // 				  HBool(
		      // 					hist(
		      // 					     "Vmet", "MET_lt", dir_, "" 
		      // 					     ),
		      // 					&MET,
		      // 					&vtxLt[iv], false, false
		      // 					)  
		      // 				   );	  		  		  
		      hbool.push_back( 
				      HBool(
					    hist(
						 "met", "MET", dir_, "" 
						 ),
					    &MET,
					    &vtxEq[iv], false, true
					    )  
				       );	  
		      hbool.push_back( 
				      HBool(
					    hist(
						 "met", "trueMET", dir_, "" 
						 ),
					    &trueMET,
					    &vtxEq[iv], false, true
					    )  
				       );	  
		      hbool.push_back( 
				      HBool(
					    hist(
						 "met", "unclMET", dir_, "" 
						 ),
					    &unclMET,
					    &vtxEq[iv], false, true
					    )  
				       );	  
		      hbool.push_back( 
				      HBool(
					    hist(
						 "met", "projMET", dir_, "" 
						 ),
					    &corProjMET,
					    &vtxEq[iv], false, true
					    )  
				       );	  
		      hbool.push_back( 
				      HBool(
					    hist(
						 "met", "ghmMET", dir_, "" 
						 ),
					    &ghmMET,
					    &vtxEq[iv], false, true
					    )  
				       );	  
		      hbool.push_back( 
				      HBool(
					    hist(
						 "met", "cpuMETSum", dir_, "" 
						 ),
					    &cpuMETSum,
					    &vtxEq[iv], false, true
					    )  
				       );	  
		      hbool.push_back( 
				      HBool(
					    hist(
						 "met", "cpuMETSqrt", dir_, "" 
						 ),
					    &cpuMETSqrt,
					    &vtxEq[iv], false, true
					    )  
				       );	  
		      hbool.push_back( 
				      HBool(
					    hist(
						 "met", "cpuMETVtx", dir_, "" 
						 ),
					    &cpuMETVtx,
					    &vtxEq[iv], false, true
					    )  
				       );	  
		      hbool.push_back( 
				      HBool(
					    hist(
						 "met", "cpuMETSumMin", dir_, "" 
						 ),
					    &cpuMETSumMin,
					    &vtxEq[iv], false, true
					    )  
				       );	  
		      hbool.push_back( 
				      HBool(
					    hist(
						 "met", "cpuMETSqrtMin", dir_, "" 
						 ),
					    &cpuMETSqrtMin,
					    &vtxEq[iv], false, true
					    )  
				       );	  
		      hbool.push_back( 
				      HBool(
					    hist(
						 "met", "cpuMETVtxMin", dir_, "" 
						 ),
					    &cpuMETMin,
					    &vtxEq[iv], false, true
					    )  
				       );	  
		      hbool.push_back( 
				      HBool(
					    hist(
						 "corsig", "ghmSig", dir_, "" 
						 ),
					    &ghmSig,
					    &vtxEq[iv], false, true
					    )  
				       );	  
		      hbool.push_back( 
				      HBool(
					    hist(
						 "met", "mEtSig", dir_, "" 
						 ),
					    &mEtSig,
					    &vtxEq[iv], false, true
					    )  
				       );	  
		      hbool.push_back( 
				      HBool(
					    hist(
						 "set", "sumEt", dir_, "" 
						 ),
					    &sumEt,
					    &vtxEq[iv], false, true
					    )  
				       );	  
		      hbool.push_back( 
				      HBool(
					    hist(
						 "bal", "balZMET", dir_, "" 
						 ),
					    &bal,
					    &vtxEq[iv], false, true 
					    )  
				       );
		  
		      // 		  hbool.push_back( 
		      // 				  HBool(
		      // 					hist(
		      // 					     "met", "MET_ge", dir_, "" 
		      // 					     ),
		      // 					&MET,
		      // 					&vtxGe[iv], false, true
		      // 					)  
		      // 				   );	  
		      // 		  hbool.push_back( 
		      // 				  HBool(
		      // 					hist(
		      // 					     "met", "MET_lt", dir_, "" 
		      // 					     ),
		      // 					&MET,
		      // 					&vtxLt[iv], false, true
		      // 					)  
		      // 				   );	  		  		  
		    }
	      


		  for( map<string,bool>::iterator it=sel.begin(); it!=sel.end(); it++ )
		    {
		      string cutName = it->first;

		      //		  dir_ =  sampleName;
		      //		  if( subsample=="data" ) dir_ = "data";
		      dir_ = sampleDir;

		      dir_ += string("/") + cutName;

		      hbool.push_back( 
				      HBool(
					    hist(
						 "met", "MET", dir_, "" 
						 ),
					    &MET,
					    &(it->second), false, true
					    )  
				       );
		      hbool.push_back( 
				      HBool(
					    hist(
						 "met", "ghmMET", dir_, "" 
						 ),
					    &ghmMET,
					    &(it->second), false, true
					    )  
				       );
		      hbool.push_back( 
				      HBool(
					    hist(
						 "met", "cpuMETSum", dir_, "" 
						 ),
					    &cpuMETSum,
					    &(it->second), false, true
					    )  
				       );
		      hbool.push_back( 
				      HBool(
					    hist(
						 "met", "cpuMETSqrt", dir_, "" 
						 ),
					    &cpuMETSqrt,
					    &(it->second), false, true
					    )  
				       );
		      hbool.push_back( 
				      HBool(
					    hist(
						 "met", "cpuMETVtx", dir_, "" 
						 ),
					    &cpuMETVtx,
					    &(it->second), false, true
					    )  
				       );
		      hbool.push_back( 
				      HBool(
					    hist(
						 "met", "cpuMETSumMin", dir_, "" 
						 ),
					    &cpuMETSumMin,
					    &(it->second), false, true
					    )  
				       );
		      hbool.push_back( 
				      HBool(
					    hist(
						 "met", "cpuMETSqrtMin", dir_, "" 
						 ),
					    &cpuMETSqrtMin,
					    &(it->second), false, true
					    )  
				       );
		      hbool.push_back( 
				      HBool(
					    hist(
						 "met", "cpuMETVtxMin", dir_, "" 
						 ),
					    &cpuMETMin,
					    &(it->second), false, true
					    )  
				       );
		      hbool.push_back( 
				      HBool(
					    hist(
						 "met", "trueMET", dir_, "" 
						 ),
					    &trueMET,
					    &(it->second), false, true
					    )  
				       );
		      hbool.push_back( 
				      HBool(
					    hist(
						 "met", "unclMET", dir_, "" 
						 ),
					    &unclMET,
					    &(it->second), false, true
					    )  
				       );
		      hbool.push_back( 
				      HBool(
					    hist(
						 "Vsig", "sigMET", dir_, "" 
						 ),
					    &sigMET,
					    &(it->second), false, false
					    )  
				       );
		      hbool.push_back( 
				      HBool(
					    hist(
						 "Vsig", "ghmSig", dir_, "" 
						 ),
					    &ghmSig,
					    &(it->second), false, false
					    )  
				       );
		      hbool.push_back( 
				      HBool(
					    hist(
						 "corsig", "ghmSig", dir_, "" 
						 ),
					    &ghmSig,
					    &(it->second), false, true
					    )  
				       );	  
		      hbool.push_back( 
				      HBool(
					    hist(
						 "mll", "mll", dir_, "" 
						 ),
					    &mll,
					    &(it->second), false, false  
					    )  
				       );
		      hbool.push_back( 
				      HBool(
					    hist(
						 "pt", "qTll", dir_, "" 
						 ),
					    &qTll,
					    &(it->second), false, true  
					    )  
				       );
		      hbool.push_back( 
				      HBool(
					    hist(
						 "y", "yll", dir_, "" 
						 ),
					    &yll,
					    &(it->second), false, false 
					    )  
				       );
		      hbool.push_back( 
				      HBool(
					    hist(
						 "bal", "balZMET", dir_, "" 
						 ),
					    &bal,
					    &(it->second), false, true 
					    )  
				       );
		      hbool.push_back( 
				      HBool(
					    hist(
						 "n20", "nVertex", dir_, "" 
						 ),
					    &nVertex,
					    &(it->second), false, true 
					    )  
				       );
		      hbool.push_back( 
				      HBool(
					    hist(
						 "cormet", "ghmMET", dir_, "" 
						 ),
					    &ghmMET,
					    &(it->second), false, true
					    )  
				       );	  
		      hbool.push_back( 
				      HBool(
					    hist(
						 "Vmet", "MET", dir_, "" 
						 ),
					    &MET,
					    &(it->second), false, false
					    )  
				       );	  
		      hbool.push_back( 
				      HBool(
					    hist(
						 "Vmet", "ghmMET", dir_, "" 
						 ),
					    &ghmMET,
					    &(it->second), false, false
					    )  
				       );	  
		      hbool.push_back( 
				      HBool(
					    hist(
						 "Vmet", "trueMET", dir_, "" 
						 ),
					    &trueMET,
					    &(it->second), false, false
					    )  
				       );	  
		      hbool.push_back( 
				      HBool(
					    hist(
						 "Vmet", "unclMET", dir_, "" 
						 ),
					    &unclMET,
					    &(it->second), false, false
					    )  
				       );	  
		      hbool.push_back( 
				      HBool(
					    hist(
						 "Vmet", "qTll", dir_, "" 
						 ),
					    &qTll,
					    &(it->second), false, false
					    )  
				       );	  
		      hbool.push_back( 
				      HBool(
					    hist(
						 "mt", "MT", dir_, "" 
						 ),
					    &thirdL_MT,
					    &(it->second), false, true
					    )  
				       );	  
		      hbool.push_back( 
				      HBool(
					    hist(
						 "mt", "mTZZ", dir_, "" 
						 ),
					    &mTZZ,
					    &(it->second), false, true
					    )  
				       );	  
		      hbool.push_back( 
				      HBool(
					    hist(
						 "n20", "jMult30", dir_, "" 
						 ),
					    &jMult30,
					    &(it->second), false, true
					    )  
				       );	  
		      hbool.push_back( 
				      HBool(
					    hist(
						 "n20", "jMult25", dir_, "" 
						 ),
					    &jMult25,
					    &(it->second), false, true
					    )  
				       );	  
		      hbool.push_back( 
				      HBool(
					    hist(
						 "n20", "jMult20", dir_, "" 
						 ),
					    &jMult20,
					    &(it->second), false, true
					    )  
				       );	  
		      hbool.push_back( 
				      HBool(
					    hist(
						 "n20", "jMult15", dir_, "" 
						 ),
					    &jMult15,
					    &(it->second), false, true
					    )  
				       );	  
		      hbool.push_back( 
				      HBool(
					    hist(
						 "n20", "jMult10", dir_, "" 
						 ),
					    &jMult10,
					    &(it->second), false, true
					    )  
				       );	  

		      hbool.push_back( 
				      HBool(
					    hist(
						 "pt", "Etj1", dir_, "" 
						 ),
					    &jEt_1,
					    &(it->second), false, true  
					    )  
				       );
		      hbool.push_back( 
				      HBool(
					    hist(
						 "pt", "Etj2", dir_, "" 
						 ),
					    &jEt_2,
					    &(it->second), false, true  
					    )  
				       );
		      hbool.push_back( 
				      HBool(
					    hist(
						 "pt", "Etj3", dir_, "" 
						 ),
					    &jEt_3,
					    &(it->second), false, true  
					    )  
				       );
		      hbool.push_back( 
				      HBool(
					    hist(
						 "pt", "Etj4", dir_, "" 
						 ),
					    &jEt_4,
					    &(it->second), false, true  
					    )  
				       );
		      hbool.push_back( 
				      HBool(
					    hist(
						 "phi", "dPhiJMET", dir_, "" 
						 ),
					    &dPhiJMET,
					    &(it->second), false, false  
					    )  
				       );
		      hbool.push_back( 
				      HBool(
					    hist(
						 "phi", "dPhiZMET", dir_, "" 
						 ),
					    &dPhiZMET,
					    &(it->second), false, false  
					    )  
				       );
		      hbool.push_back( 
				      HBool(
					    hist(
						 "phi", "dPhiSumJMET", dir_, "" 
						 ),
					    &dPhiSumJMET,
					    &(it->second), false, false  
					    )  
				       );
		      hbool.push_back( 
				      HBool(
					    hist(
						 "btagtk", "jBtagTkC", dir_, "" 
						 ),
					    &jBtagTkC_max,
					    &(it->second), false, true  
					    )  
				       );
		      hbool.push_back( 
				      HBool(
					    hist(
						 "btagmu", "jBtagSoftM", dir_, "" 
						 ),
					    &jBtagSoftM_max,
					    &(it->second), false, true  
					    )  
				       );
		    }
		}
	      float w_ = sampleWeight * eventWeight;
	      for( size_t ih=0; ih<hbool.size(); ih++ )
		{
		  hbool[ih].fill( w_ );
		}

	      sel_ntu["cut"] = cut();
	      sel_ntu["jet_met_bal"] = b_jet && b_met && b_bal;
	      sel_ntu["jet_bid"] = b_jet && !b_dib;
	      sel_ntu["q30_jet_met_bal_phj_phz_btg_dib"] 
		= b_q30 && b_jet && b_met && b_bal && b_phj && b_phz && b_btg && b_dib;
	      sel_ntu["q30_jet_met_bal_phj_phz_btg_bid"] 
		= b_q30 && b_jet && b_met && b_bal && b_phj && b_phz && b_btg && !b_dib;
	      
	      if( subsample=="data" && 
		  ( sel_ntu["met_bal_jet"] || sel_ntu["jet_bid"] ) )
		{
		  ifile << run << " " << event;
		  ifile << endl;
		}

	      
	      for( map<string,bool>::iterator it=sel_ntu.begin(); it!=sel_ntu.end(); it++ )
		{	      		      
		  if( it->second )
		    {
		      string cutName = it->first;

		      //		  dir_ =  sampleName;
		      //		  if( subsample=="data" ) dir_ = "data";
		      string dir_ = sampleDir;

		      dir_ += string("/") + cutName;
  
		      string tuplename_ = "ZZtuple";
		      tmsel.setTree( tuplename_, dir_ );
		      tmsel.add< int >(    "run",   &run       );
		      tmsel.add< int >(    "event", &event     );
		      tmsel.add< string*>( "categ", &categ_ptr );
		      tmsel.add< string*>( "c_str", &c_str_ptr  );      
		      tmsel.add< TBits* >( "c_bits", &c_bits_ptr  );      
		      tmsel.add< int >(  "nZCand",   &nZCand        );
		      tmsel.add< int >(  "nVertex",   &nVertex      );
		      tmsel.add< bool  >(  "isTight",&isTight    );
		      tmsel.add< float >(  "mll",   &mll        );
		      tmsel.add< float >(  "qTll",   &qTll        );
		      tmsel.add< float >(  "phill",  &phill        );
		      tmsel.add< float >(  "yll",   &yll        );
		      tmsel.add< float >(  "hll",   &hll        );
		      tmsel.add< int >(  "nLoose",&nLoose    );
		      tmsel.add< vectorFloat* >(  "LPt",   &LPt   );
		      tmsel.add< vector<int>* >(  "LPdg",   &LPdg   );
		      tmsel.add< vectorFloat* >(  "LEta",  &LEta  );
		      tmsel.add< vectorFloat* >(  "LPhi",  &LPhi  );
		      tmsel.add< vector<int>* >(  "LID",   &LID   );
		      tmsel.add< vectorFloat* >(  "LMT",   &LMT   );
		      tmsel.add< vectorFloat* >( "LTrkIso", &LTrkIso );
		      tmsel.add< vectorFloat* >( "LEcalIso", &LEcalIso );
		      tmsel.add< vectorFloat* >( "LHcalIso", &LHcalIso );
		      tmsel.add< float >(  "MET",   &MET        );
		      tmsel.add< float >(  "METPhi",&METPhi    );
		      tmsel.add< float >(  "sumEt", &sumEt );
		      tmsel.add< float >(  "mEtSig", &mEtSig );
		      tmsel.add< float >(  "mTZZ", &mTZZ );
		      tmsel.add< float >(  "projMET",&projMET    );
		      tmsel.add< float >(  "corProjMET",&corProjMET    );
		      tmsel.add< float >(  "sigMET",&sigMET    );
		      tmsel.add< float >(  "dPhiMin",&dPhiMin    );
		      tmsel.add< float >(  "balZMET",&bal    );  //!!!
		      tmsel.add< vectorFloat* > ( "jEt", &jEt );
		      tmsel.add< vectorFloat* > ( "jEta", &jEta );
		      tmsel.add< vectorFloat* > ( "jPhi", &jPhi );
		      tmsel.add< vectorFloat* > ( "jDR1", &jDR1 );
		      tmsel.add< vectorFloat* > ( "jDR2", &jDR2 );
		      tmsel.add< vectorFloat* > ( "jDPhiMET", &jDPhiMET );
		      tmsel.add< vectorFloat* > ( "jBtagTkC", &jBtagTkC );
		      tmsel.add< vectorFloat* > ( "jBtagSoftM", &jBtagSoftM );
		      tmsel.add< int >(  "jMult25",&jMult25  );
		      tmsel.add< int >(  "jMult20",&jMult20  );
		      tmsel.add< int >(  "jMult15",&jMult15  );
		      tmsel.add< int >(  "jMult10",&jMult10  );
		      tmsel.add< float >(  "JZB25",&JZB25  );
		      tmsel.add< float >(  "JZB20",&JZB20  );
		      tmsel.add< float >(  "JZB15",&JZB15  );
		      tmsel.add< float >(  "dPhiJMET",&dPhiJMET );	  
		      tmsel.flush();
		      
		      // 		  string dir_ = sampleDir;
	      
		      dm.add( "mll", mll );
		      dm.flush( dir_ );		      
		    }
		}		  
	    }
	  cout << "Number of entries           " << ievent_0 << endl;
	  cout << "Number of different events  " << ievent << endl;

    
	}
      cout << "nsel=" << nsel << endl;
      inputFile->Close();
    }

  ifile.close();

  hm.save();
  dm.save();
  tmsel.save();

}

bool
ZZ_2l2n_Analysis::analyzeSignalMC()
{
  // is the MC event in acceptance?

  // the four leptons
  float ptL_0 = (*mc_ptL)[0];
  float ptL_1 = (*mc_ptL)[1];
  float ptL_2 = (*mc_ptL)[2];
  float ptL_3 = (*mc_ptL)[3];
  float phiL_0 = (*mc_phiL)[0];
  float phiL_1 = (*mc_phiL)[1];
  float phiL_2 = (*mc_phiL)[2];
  float phiL_3 = (*mc_phiL)[3];
  float etaL_0 = (*mc_etaL)[0];
  //  float etaL_1 = (*mc_etaL)[1];
  //  float etaL_2 = (*mc_etaL)[2];
  //  float etaL_3 = (*mc_etaL)[3];
  int pdgL_0 = (*mc_pdgL)[0];
  //  int pdgL_1 = (*mc_pdgL)[1];
  //  int pdgL_2 = (*mc_pdgL)[2];
  //  int pdgL_3 = (*mc_pdgL)[3];
  TVector2 ptL_0_2D( ptL_0*cos( phiL_0 ), ptL_0*sin( phiL_0 ) );
  TVector2 ptL_1_2D( ptL_1*cos( phiL_1 ), ptL_1*sin( phiL_1 ) );
  TVector2 ptL_2_2D( ptL_2*cos( phiL_2 ), ptL_2*sin( phiL_2 ) );
  TVector2 ptL_3_2D( ptL_3*cos( phiL_3 ), ptL_3*sin( phiL_3 ) );

  isInAcceptance = true;
  isQt30 = true;
  TVector2 ptV_0_2D;
  TVector2 ptV_1_2D;
  if( isZZ )
    {
      if( ptL_0<20 || ptL_1<20 ) isInAcceptance = false;
      ptV_0_2D = ptL_0_2D + ptL_1_2D;
      ptV_1_2D = ptL_2_2D + ptL_3_2D;
    }
  else if( isWZ )
    {
      if( ptL_2<20 || ptL_3<20 ) isInAcceptance = false;
      ptV_0_2D = ptL_2_2D + ptL_3_2D;
      ptV_1_2D = ptL_1_2D;
      if( ( abs( pdgL_0 )==11 && fabs(etaL_0)<2.4 && ptL_0>20 )
	  ||
	  ( abs( pdgL_0 )==13 && fabs(etaL_0)<2.1 && ptL_0>20 ) )
	ptV_1_2D += ptL_0_2D;
    }
  else if( isWW )
    {
      if( ptL_1<20 || ptL_3<20 ) isInAcceptance = false;
      ptV_0_2D = ptL_0_2D + ptL_2_2D;
      ptV_1_2D = ptL_1_2D + ptL_3_2D;
    }
  if( ptV_0_2D.Mod()<30 ) isQt30 = false;
  if( ptV_1_2D.Mod()<30 ) isInAcceptance = false;      

  float ptV = 0;
  KFactor* vvkf(0);
  KFactor* vvkf_jv30(0);
  if( isZZ )
    {
      if( ew!=0 )  aTGCfactor = (*ew)( mc_mVV );
      vvkf      = kFactors["ZZ"];
      vvkf_jv30 = kFactors["ZZ_jv30"];
      ptV = mc_ptV1;
      
    }
  else if( isWZ )
    {
      vvkf      = kFactors["WZ"];
      vvkf_jv30 = kFactors["WZ_jv30"];
      ptV = mc_ptV2;
    }
  else if( isWW )
    {
      vvkf      = kFactors["WW"];
      vvkf_jv30 = kFactors["WW_jv30"];
      ptV = mc_ptV1;
    }

  if( vvkf ) 
    {
      kfactor      *= (*vvkf)( ptV );		  
      kfactor_jv30 *= (*vvkf_jv30)( ptV );
    }

  return true;
}

ZZ_2l2n_Analysis::ZZ_2l2n_Analysis( int im, string ip, bool ZJets_cs ) 
    : TupleAnalysis(im,ip), ew(0), _ZJets_cs(ZJets_cs)
{  
  if( imode==0 )
    analysis = "ZZ_2e2n";
  else if( imode==1 )
    analysis = "ZZ_2m2n";
  else if( imode==2 )
      analysis = "ZZ_2l2n";
  
  _signalOnly = false;

  if( _ZJets_cs ) analysis += "_cs";
  if( _signalOnly ) analysis += "_sig";

  // initialisations
  {

    run=0;
    event=0;  
    mc_event=0;
    mc_run=0;
    mc_mVV=0;  
    mc_mV1=0;
    mc_ptV1=0;
    mc_pzV1=0;
    mc_EV1=0;
    mc_yV1=0;
    mc_mV2=0;
    mc_ptV2=0;
    mc_pzV2=0;
    mc_EV2=0;
    mc_yV2=0;
    nZCand=0;      
    nVertex=0;     
    mll=0;         
    qTll=0;        
    phill=0;       
    yll=0;         
    hll=0;         
    MET=0;         
    METPhi=0;      
    sumEt=0;       
    mEtSig=0;      
    mTZZ=0;
    projMET=0;     
    corProjMET=0;  
    sigMET=0;      
    ghmMET=0;      
    cpuMETSum=0;
    cpuMETSqrt=0;
    cpuMETVtx=0;
    cpuMETSumMin=0;
    cpuMETSqrtMin=0;
    cpuMETVtxMin=0;
    trueMET=0;
    unclMET=0;
    ghmSig=0;
    dPhiMin=0;     
    balZMET=0;     
    isTight=false;  
    isLLTight=false;
    nLoose=0;   
    jMult40=0;  
    jMult35=0;  
    jMult30=0;  
    jMult25=0;  
    jMult20=0;  
    jMult15=0;  
    jMult10=0;  
    JZB15=0;
    JZB20=0;
    JZB25=0;
    dPhiJMET=0; 
    jEt_1=0;
    jEt_2=0;
    jEt_3=0;
    jEt_4=0;
    dPhiSumJMET=0;
    dPhiZMET=0;
    thirdL_MT=0;
    jBtagTkC_max=0;
    jBtagSoftM_max=0;
    bal=0;

    // vectors
    jEt      = new vectorFloat;
    jEta     = new vectorFloat;
    jPhi     = new vectorFloat;
    jDR1     = new vectorFloat;
    jDR2     = new vectorFloat;
    jDPhiMET = new vectorFloat;
    jBtagTkC   = new vectorFloat;
    jBtagSoftM = new vectorFloat;
    LPt  = new vectorFloat;
    LPdg = new vector<int>;
    LEta = new vectorFloat;
    LPhi = new vectorFloat;
    LMT   = new vectorFloat;
    LID  = new vector<int>;
    LTrkIso  = new vectorFloat;
    LEcalIso = new vectorFloat;
    LHcalIso = new vectorFloat;
    mc_pdgL = new vector<int>;
    mc_ptL  = new vectorFloat;
    mc_etaL = new vectorFloat;
    mc_phiL = new vectorFloat;

    eventWeight=1;
    aTGCfactor=1;
    kfactor=1;
    kfactor_jv30=1;
    
    // is it in acceptance?
    isInAcceptance=false;
    
  }
}

ZZ_2l2n_Analysis::~ZZ_2l2n_Analysis() 
{ 
  delete ew; 

  // vectors
  delete jEt;
  delete jEta;
  delete jPhi;
  delete jDR1;
  delete jDR2;
  delete jDPhiMET;
  delete jBtagTkC;
  delete jBtagSoftM;
  delete LPt;
  delete LPdg;
  delete LEta;
  delete LPhi;
  delete LMT;
  delete LID;
  delete LTrkIso;
  delete LEcalIso;
  delete LHcalIso;
  delete mc_pdgL;
  delete mc_ptL;
  delete mc_etaL;
  delete mc_phiL;
}
