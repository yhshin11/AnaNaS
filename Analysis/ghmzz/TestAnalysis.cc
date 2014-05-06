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

#include "Analysis/utils/CutUtils.hh"
#include "Analysis/utils/Config.hh"
#include "Analysis/ghmzz/TestAnalysis.hh"
		
ClassImp( TestAnalysis )
	     
void 
TestAnalysis::analyze_( float lumi )
{
  cout << "Hello World " << endl;
  cout << "This is a simple test for lumi=" << lumi << endl;


  cout << "xsect[ZZ_2l2n]=" << xsect["ZZ_2l2n_[60-120]"] << " pb, Mll in [60-120] GeV" << endl;
  cout << "xsect[ZZ_2l2n]=" << xsect["ZZ_2l2n"] << " pb" << endl;
  cout << "xsect[ZZ]=" << xsect["ZZ_[60-120]"] << " pb, Mll in [60-120] GeV" << endl;
  cout << "xsect[ZZ]*BR[2l2n]=" << xsect["ZZ_[60-120]"] * BR["ZZ_2l2n"] << " pb, Mll in [60-120] GeV" << endl;


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

  // list of samples samples
  typedef multimap< string, Sample* > SampleMap;
  typedef multimap< string, Sample* >::iterator SampleIt;
  typedef pair< string, Sample* > SamplePair;
  SampleMap sampleList;

  vector<string> ZZanal;

  if( imode==0 )
    // electron analysis
    {
      ZZanal.push_back("ZZ_2e2n");

      // DY
      //      sampleList.insert(SamplePair("Z_2e",samples["Z_2e"]));
	    
      // the signal
      sampleList.insert(SamplePair("ZZ_2l2n",samples["ZZ_2l2n"]));
      sampleList.insert(SamplePair("WZ_3ln",samples["WZ_3ln"]));
    }
  else if( imode==1 )
    {
      ZZanal.push_back("ZZ_2m2n");

      // DY
      //      sampleList.insert(SamplePair("Z_2m",samples["Z_2m"]));
      
      // the signal
      sampleList.insert(SamplePair("ZZ_2l2n",samples["ZZ_2l2n"]));
      sampleList.insert(SamplePair("WZ_3ln",samples["WZ_3ln"]));
    }
  else if( imode==2 )
    {
      ZZanal.push_back("ZZ_2e2n");
      ZZanal.push_back("ZZ_2m2n");

      // DY
      //      sampleList.insert(SamplePair("Z_2l",samples["Z_2e"]));
      //      sampleList.insert(SamplePair("Z_2l",samples["Z_2m"]));
	
      sampleList.insert(SamplePair("ZZ_2l2n",samples["ZZ_2l2n"]));
      sampleList.insert(SamplePair("WZ_3ln",samples["WZ_3ln"]));
    }

  string sampleAnalysis("Diboson");
  string tname("ZZtuple");
  string mc_tname("VV_mcTruth");

  for ( SampleIt it=sampleList.begin() ; it != sampleList.end(); it++ )
    {
      string sampleDir = it->first;  // defines the histogram directory
	
      cout << "Histo Directory = " << sampleDir << endl;
	
      //  for( size_t isample=0; isample<sampleList.size(); isample++ )
      //    {
      Sample* sample_ = it->second;
	
      string sampleName = sample_->name();
      size_t pos = sampleName.find('/');  
      string subsample = sampleName.substr(pos+1);
      string sampleRoot = sampleName.substr(0,pos);

      bool isZZ =  sampleRoot.find("ZZ")!=string::npos;
      bool isWZ =  sampleRoot.find("WZ")!=string::npos;
      bool isWW =  sampleRoot.find("WW")!=string::npos;
	
      //    for( size_t isample=0; isample<sampleList.size(); isample++ )
      //	string sampleRoot = "ZZ_2l2n";
      //	Sample* sample_ = samples["ZZ_2e2n"];
      sample_->print( cout );
      cout << "n_tot(" << lumi << ")=" << sample_->n_tot( lumi ) << endl; 
      cout << "weight(" << lumi << ")=" << sample_->w( lumi ) << endl;

      TString filename = Config::rootPath;
      filename += sampleAnalysis;
      filename += "/";
      filename += sampleRoot;
      filename += ".root";
            
      TFile* inputFile  = new TFile( filename, "READ" );
	
      // tuple manager
      TupleManager tm;
      tm.setFile( inputFile, TupleManager::kRead );
    
      float sampleWeight = sample_->w( lumi );
      sampleWeight = 1;
      // create the sample weight
      VertexWeight vw( sampleRoot, period );      

      //*************

      int nsel(0);
      for( size_t ianal=0; ianal<ZZanal.size(); ianal++ )
	{
	  bool isMuMu = ( ZZanal[ianal]=="ZZ_2m2n" );
	  bool isEE   = ( ZZanal[ianal]=="ZZ_2e2n" );
	  assert( isMuMu || isEE );

	  // tree
	  TTree* tree = tm.getTree( tname, ZZanal[ianal] );
	  if( tree==0 ) continue;
	  int n = tree->GetEntriesFast();

	  // declare the main analysis variables
	
	  int run;
	  int event;
	  string  categ;
	  string* categ_ptr = & categ;
	  int   nZCand(0);        // number of Z candidates
	  int nVertex(0);
	  bool  isTight(false);   // is the Z candidate tight?
	  bool  isLLTight(false);   // is the Z candidate tight?
	  float mll(0);           // di-lepton mass
	  float qTll(0);          // di-lepton transverse momentum
	  float phill(0);         // di-lepton phi
	  float yll(0);           // di-lepton rapidity
	  float hll(0);           // di-lepton helicity
	  float mTZZ(0);          // transverse mass of Z pair
	  float MET(0);           // PF-MET type-I
	  int   jMult30(0);       // number of jets above 30 GeV-Et
	  int   jMult25(0);       // number of jets above 25 GeV-Et
	  vector<int>* LID = new vector<int>;

	  tm.add< int >(    "run",   &run       );
	  tm.add< int >(    "event", &event     );
	  tm.add< int >(    "nZCand", &nZCand     );
	  tm.add< int >(  "nVertex",   &nVertex      );
	  tm.add< bool  >(  "isTight",&isTight    );
	  tm.add< float >(  "mll",   &mll        );
	  tm.add< float >(  "qTll",   &qTll        );
	  tm.add< float >(  "phill",  &phill        );
	  tm.add< float >(  "yll",   &yll        );
	  tm.add< float >(  "hll",   &hll        );
	  tm.add< float >(  "MET",   &MET        );
	  tm.add< float >(  "mTZZ", &mTZZ );
	  tm.add< int >(  "jMult30",&jMult30  );
	  tm.add< int >(  "jMult25",&jMult25  );
	  tm.add< vector<int>* >(  "LID",   &LID   );

	  //
	  // Cuts
	  //
	  bool presel=true;
	  OneSidedCut<bool>  c_presel    ( "presel",      presel, true,     "==" );	

	  CompositeCut        c_ZQual     ( "ZQual", "&&" );
	  OneSidedCut<int>   c_nZCand    ( "nZCand",      nZCand,    1,   "==" );
	  TwoSidedCut<float> c_ZWindow   ( "ZWindow",     mll, 60., 120., "[]" );
	  OneSidedCut<bool>  c_ZLQual    ( "ZLQual",     isTight, true,   "==" );
	  OneSidedCut<bool>  c_ZLLTight ( "ZLLTight",   isLLTight, true,   "==" );
	  c_ZQual.add( c_nZCand );
	  c_ZQual.add( c_ZWindow );
	  c_ZQual.add( c_ZLQual );
	  c_ZQual.add( c_ZLLTight );

	  OneSidedCut<float> c_qT30      ( "qT30",          qTll, 30.,      ">=" );
	  OneSidedCut<int>   c_j30Veto   ( "j30Veto",    jMult30, 0,        "==" );
	
	  TTree* mc_tree = tm.getTree( mc_tname, ZZanal[ianal] );

	  int mc_event(0);
	  int mc_run(0);
	  float mc_mVV(0);
	
	  float  mc_mV1(0);
	  float mc_ptV1(0);
	  float mc_pzV1(0);
	  float  mc_EV1(0);
	  float  mc_yV1(0);
	
	  float  mc_mV2(0);
	  float mc_ptV2(0);
	  float mc_pzV2(0);
	  float  mc_EV2(0);
	  float  mc_yV2(0);

	  //new
	  vectorFloat* mc_pdgL = new vectorFloat;
	  vectorFloat* mc_ptL  = new vectorFloat;
	  vectorFloat* mc_etaL = new vectorFloat;
	  vectorFloat* mc_phiL = new vectorFloat;
	  //new
	
	  int n_mc(0);
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
	      
	      //new
	      tm.add< vectorFloat* > ( "pdgL", &mc_pdgL );
	      tm.add< vectorFloat* > ( "ptL",  &mc_ptL  );
	      tm.add< vectorFloat* > ( "etaL", &mc_etaL );
	      tm.add< vectorFloat* > ( "phiL", &mc_phiL );	      
	      //new
	    }		


	  cout << "Number of events in tree    = " << n << endl;
	  cout << "Number of events in MC tree = " << n_mc << endl;

	  float ptV(0);
	    
	  for( int ii=0; ii<n; ii++ )
	    {
	      tree->GetEntry( ii );

	      float eventWeight = vw( nVertex );
	      float eventWeight_jv30 = eventWeight;

	      // acceptance
	      bool leptonEta(true);
	      bool leptonPt(true);
	      bool dileptonM(true); 
	      bool bosonQt(true);
	      bool leptonOK(true);
	      bool dileptonOK(true);
	      bool inAcceptance(true);
	      if( mc_tree )
		{
		  mc_tree->GetEntry( ii );
		  assert( mc_run==run && mc_event==event);


		  //
		  // acceptance cuts
		  //

		  // acceptance
		  leptonEta   = true;
		  leptonPt    = true;
		  int il0 = 0;
		  if( isWZ ) il0=2;
		  if( isWZ )
		    {
		      if( categ == "WZ_3en" || categ == "WZ_2emn" ||
			  categ == "WZ_e2mn" || categ== "WZ_3mn"     )
			{
			  // OK
			}
		      else
			continue;
		    }
		  for( int il=il0; il<il0+2; il++ )
		    {
		      if( fabs((*mc_ptL)[il])<20 ) 
			{
			  leptonPt = false;
			}
		      if( abs((*mc_pdgL)[il]) == 11 )
			{
			  if( fabs((*mc_etaL)[il])>2.5 ) 
			    {
			      leptonEta = false;
			    }	      
			}	 
		      else if( abs((*mc_pdgL)[il]) == 13 )
			{
			  if( fabs((*mc_etaL)[il])>2.4 ) 
			    {
			      leptonEta = false;
			    }	      
			}
		      else
			leptonEta = false;
		    }
		  dileptonM   = mc_mV1>=60 && mc_mV1<=120;
		  bosonQt     = mc_ptV1>30.;
		  leptonOK     = leptonEta  && leptonPt;
		  dileptonOK   = leptonOK   && dileptonM;
		  inAcceptance = dileptonOK && bosonQt;
		    
		  KFactor* vvkf(0);
		  KFactor* vvkf_jv30(0);
		  if( sampleRoot.find("ZZ")!=string::npos )
		    {
		      //		if( ew!=0 )
		      //		  eventWeight *= (*ew)( mc_mVV );
		      vvkf      = kFactors["ZZ"];
		      vvkf_jv30 = kFactors["ZZ_jv30"];
		      ptV = mc_ptV1;
		    }
		  else if( sampleRoot.find("WZ")!=string::npos )
		    {
		      vvkf      = kFactors["WZ"];
		      vvkf_jv30 = kFactors["WZ_jv30"];
		      ptV = mc_ptV2;
		    }
		  else if( sampleRoot.find("WW")!=string::npos )
		    {
		      vvkf      = kFactors["WW"];
		      vvkf_jv30 = kFactors["WW_jv30"];
		      ptV = mc_ptV1;
		    }
		  if( vvkf ) 
		    {
		      eventWeight      *= (*vvkf)( ptV );		  
		      eventWeight_jv30 *= (*vvkf_jv30)( ptV );		  
			
		      if( ii<100 )
			{
			  cout << "run/event " << run  << "/" << event 
			       << " MVV=" << mc_mVV << " pTV=" << ptV 
			       << " pdg0=" << (*mc_pdgL)[0]
			       << " pdg1=" << (*mc_pdgL)[1]
			       << " weight " << eventWeight
			       << " weight(JV30) " << eventWeight_jv30 
			       << endl;
			}
		    }
		}

	      if( !dileptonOK ) continue;
	      nsel++;

              isLLTight = true;
              if( isEE )
                {
                  for( size_t ij=0; ij<2; ij++ )
                    {
                      int   LID_   = (*LID)[ij];
                      if( LID_<7 ) isLLTight=false;
                    }             
                }

	      bool b_pre = true;
	      bool b_dil = c_ZQual();
	      bool b_q30 = c_qT30();
	      bool b_jet = c_j30Veto();

	      map< string, bool > sel;
	      sel["presel"] = b_pre;
	      sel["ZQual"]  = b_pre && b_dil;
	      sel["q30"]    = b_pre && b_dil && b_q30;
	      sel["jet"]    = b_pre && b_dil && b_q30 && b_jet;
	       	    	    
	      for( map<string,bool>::iterator it=sel.begin(); 
		   it!=sel.end(); it++ )
		{
		  string cutName = it->first;
		  string dir_ = sampleDir;		    
		  dir_ += string("/") + cutName;
		
		  if( !(it->second) ) continue;
		  hist( "pt", "qTll_noKF",    dir_, "" )->Fill( qTll, sampleWeight );  
		  hist( "pt", "qTll_KF",      dir_, "" )->Fill( qTll, sampleWeight*eventWeight );  
		  hist( "pt", "qTll_KF_jv30", dir_, "" )->Fill( qTll, sampleWeight*eventWeight_jv30 );  
		}
	
	    }

	  //new
	  delete mc_pdgL;
	  delete mc_ptL;
	  delete mc_etaL;
	  delete mc_phiL;      
	  delete LID;
	  //new
	}
      cout << "nsel=" << nsel << endl;
      inputFile->Close();

    }
  hm.save();
}

