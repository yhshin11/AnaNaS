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
#include <TAxis.h>
#include <TH1F.h>
#include <TLorentzVector.h>

#include "Analysis/core/Candidate.hh"
#include "Analysis/core/CandInfo.hh"
#include "Analysis/tools/CandUtil.hh"
#include "Analysis/utils/Config.hh"
#include "Analysis/utils/DatasetManager.hh"
#include "Analysis/utils/TupleManager.hh"
#include "Analysis/utils/HistoManager.hh"
#include "Analysis/utils/Constants.hh"
#include "Analysis/utils/KineUtils.hh"

using namespace std;

typedef vector<double> vectorFloat; 

class DummyTupleAnalysis
{
public:
  TupleManager   tm;
  HistoManager   hm;                                                                            
  void analyze();

  void declareInt(   string var );
  void declareFloat( string var );
  void declareBool(  string var );
  map< string, float > _float; 
  map< string, int   > _int; 
  map< string, bool  > _bool; 

  void defineTemplate( string templ, int nbin, float min, float max );
  TH1* hist( string templ, string cut, string dir, string prefix );
  void fill( string templ, string cut, string dir, string prefix,
	     float val, float weight=1, bool overflow=false, bool norm=true );

};

void
DummyTupleAnalysis::defineTemplate( string templ, 
				    int nbin, float min, float max )
{
  cout << "1D-Template[" << templ << "]-->(nbin=" << nbin << "/min=" << min << "/max=" << max << ")" << endl; 

  TH1F* h_ = new TH1F( templ.c_str(), templ.c_str(), 
		       nbin, min, max );
  h_->Sumw2();
  hm.addTemplate<TH1F>( templ, h_ );
}

TH1* 
DummyTupleAnalysis::hist( string templ, string cut, string dir, string prefix )
{
  TH1* h_ = (TH1*)hm.h<TH1F>( templ, cut, dir, prefix );
  assert( h_!=0 );
  h_->SetMarkerSize(0);
  return h_;
}

void
DummyTupleAnalysis::fill( string templ, string cut, string dir, string prefix,
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
DummyTupleAnalysis::declareInt( string var )
{
  string tname  = tm.tree()->GetName();
  string varKey = tm.varKey( tname, var ); 
  tm.add<int>(    var,    &_int[varKey] );
}

void
DummyTupleAnalysis::declareFloat( string var )
{
  string tname  = tm.tree()->GetName();
  string varKey = tm.varKey( tname, var ); 
  tm.add<float>(    var,    &_float[varKey] );
}

void
DummyTupleAnalysis::declareBool( string var )
{
  string tname  = tm.tree()->GetName();
  string varKey = tm.varKey( tname, var ); 
  tm.add<bool>(    var,    &_bool[varKey] );
}

void 
DummyTupleAnalysis::analyze()
{
  //  string filename = Config::rootPath + "Zll/Fall11/ZZ_2e2m.root";
  string filename = Config::rootPath + "Zll/data2011A/DoubleElectron_173000_173389.root";
  TFile* inputFile  = TFile::Open( filename.c_str(), "READ" );
  inputFile->Print();
  
  tm.reset();
  tm.setFile( inputFile,  TupleManager::kRead );

  string dir_("Zll");

  TTree* eventTuple = tm.getTree( "EventTuple", dir_ );
  {
    declareInt(    "run"     );
    declareInt(    "event"   );
    declareInt(    "leptonVertices"   );
    declareInt(    "Electron"   );
    declareInt(    "Muon"   );
    declareInt(    "Zee"   );
    declareInt(    "Zmm"   );
    declareInt(    "CleanJets"   );
  }

  TTree* vtxTuple = tm.getTree( "leptonVertices", dir_ );
  {
    declareInt(    "run"     );
    declareInt(    "event"   );
    declareInt(    "counter"   );
    declareInt(    "id"   );
    declareFloat( "X" );
    declareFloat( "Y" );
    declareFloat( "Z" );
  }

  TTree* leptonTuple[2];
  string leptonName[2] = { "Electron", "Muon" };
  for( int lmode=0; lmode<2; lmode++ )
    {
      leptonTuple[lmode] = tm.getTree( leptonName[lmode], dir_ );
      {
	declareInt(    "run"     );
	declareInt(    "event"   );
	declareInt(    "counter"   );
	declareInt(    "index"   );
	declareFloat( "pt" );
	declareFloat( "eta" );
	declareFloat( "phi" );
	declareInt(   "vtxId"     );
	declareFloat(    "x0" );
	declareFloat(    "y0" );
	declareFloat(    "z0" );
	declareInt(   "id"       );
      }
      if( lmode==0 )
	{
	  declareFloat( "fBrem"     );
	  declareFloat( "EOP"       );
	  declareFloat( "dPhiIn"    );
	  declareFloat( "dEtaIn"    );
	  declareFloat( "dr03TrkIso"    );
	  declareFloat( "dr04EcalIso"    );
	  declareFloat( "dr04HcalD1Iso"    );
	  declareFloat( "dr04HcalD2Iso"    );
	  declareBool( "isConv"    );
	  declareBool( "isFPHitTrk"    );
	  declareInt( "expInnHits"    );
	  declareInt( "ix"    );
	  declareInt( "iy"    );
	  declareInt( "iz"    );
	}
      else
	{
	}
    }
  
  TTree* ZllTuple[2];
  string ZllName[2] = { "Zee", "Zmm" };
  
  for( int lmode=0; lmode<2; lmode++ )
    {
      //      TTree* ZllTuple_ = tm.getTree( ZllName[lmode], dir_ );
      ZllTuple[lmode] = tm.getTree( ZllName[lmode], dir_ );
      
      declareInt(    "run"     );
      declareInt(    "event"   );
      declareInt(    "counter" );
      declareInt(    "dau0"    );
      declareInt(    "dau1"    );
      
      declareInt(    "mode"    );
      declareInt(    "hyp"     );
      declareInt(    "iZ"      );
      declareBool(  "isTight"  );
      declareFloat(  "mll"     );
      declareFloat(  "qTll"    );
      declareFloat(  "yll"     );
      declareFloat(  "hll"     );

      // vertices
      declareInt(   "ZV_id" );
      declareFloat(  "ZV_X" );
      declareFloat(  "ZV_Y" );
      declareFloat(  "ZV_Z" );      
    }    
  
  int n   = eventTuple->GetEntriesFast();  
  //  int n_max = 10;
  //  if( n>n_max ) n=n_max;
  
  int jVtx=0;
  //  int jJet=0;
  vector<int> jLep(2,0);
  vector<int> jZll(2,0);

  bool _verbose(false);

  for( int ii=0; ii<n; ii++ )
    {
      eventTuple->GetEntry( ii );
      string tupleName = "EventTuple_";
      int run   = _int[tupleName+"run"];
      int event = _int[tupleName+"event"];
      int nVtx  = _int[tupleName+"leptonVertices"];      
      vector<int> nLep(2);
      vector<int> nZll(2);
      for( int lmode=0; lmode<2; lmode++ )
	{
	  nLep[lmode] = _int[tupleName+leptonName[lmode]];
	  nZll[lmode] = _int[tupleName+ZllName[lmode]];
	}
      int nJet  = _int[tupleName+"CleanJets"];

      if( _verbose )
	{
	  cout << run << "/" << event << endl;      
	  cout << "Vertices   --> " <<  nVtx  << endl;
	  cout << "CleanJets  --> " <<  nJet  << endl;
	  cout << "Electrons  --> " <<  nLep[0]  << endl;
	  cout << "Muons      --> " <<  nLep[1]  << endl;
	  cout << "Zee        --> " <<  nZll[0]  << endl;
	  cout << "Zmm        --> " <<  nZll[1]  << endl;
	}

      map<int,Vertex*>     vertices;
      vector< map<int,Candidate*> > leptons(2);

      tupleName = "leptonVertices_";
      for( int iVtx=0; iVtx<nVtx; iVtx++ )
	{
	  vtxTuple->GetEntry( jVtx++ );
	  assert( _int[tupleName+"run"]==run );
	  assert( _int[tupleName+"event"]==event );

	  int index =  _int[tupleName+"id"];
	  assert( vertices.count(index)==0 );
	  
	  Vertex* vtx_ = Vertex::create();
	  vtx_->setXYZ( _float[tupleName+"X"],
			_float[tupleName+"Y"], 
			_float[tupleName+"Z"] );

	  // special track information
	  CandInfo* info_ = CandInfo::create();
	  info_->setInt( "index", index );
	  vtx_->setInfo( info_  );
	  if( index==0 ) vtx_->setAsPrimary();

	  vertices[index] =  vtx_;
	}

      for( int lmode=0; lmode<2; lmode++ )
	{
	  tupleName = leptonName[lmode]+"_";
	  for( int iLep=0; iLep<nLep[lmode]; iLep++ )
	    {
	      leptonTuple[lmode]->GetEntry( jLep[lmode]++ );
	      assert( _int[tupleName+"run"]==run );
	      assert( _int[tupleName+"event"]==event );

	      if( _verbose )
		{
		  cout << tupleName << iLep 
		       << " counter=" << _int[tupleName+"counter"] 
		       << " index=" << _int[tupleName+"index"] 
		       << " pt="    << _float[tupleName+"pt"] 
		       << endl;
		}
		  
	      int pdgId_ = 11;
	      if( _float[tupleName+"pt"]>0 ) pdgId_*=-1; 
	      TVector3 p3_(0,0,0);
	      p3_.SetPtEtaPhi( 
			      fabs( _float[tupleName+"pt"] ), 
			      _float[tupleName+"eta"], 
			      _float[tupleName+"phi"] );
	      int index = _int[tupleName+"index"];
	      int vtxId = _int[tupleName+"vtxId"]; 
	      assert( vertices.count(vtxId)!=0 );
	  
	      Candidate* cand_= Candidate::create( p3_, pdgId_, vertices[vtxId] ); 
		  
	      // special electron information
	      CandInfo* info_ = CandInfo::create();
	      if( lmode==0 )
		{
		  info_->setFloat( "fBrem", _float[tupleName+"fBrem"] );
		  info_->setFloat( "EOP",_float[tupleName+"EOP"] );
		  info_->setFloat( "dPhiIn", _float[tupleName+"dPhiIn"] );
		  info_->setFloat( "dEtaIn", _float[tupleName+"dEtaIn"] );
		  info_->setFloat( "dr03TrkIso", _float[tupleName+"dr03TrkIso"] );
		  info_->setFloat( "dr04EcalIso", _float[tupleName+"dr04EcalIso"] );
		  info_->setFloat( "dr04HcalD1Iso",_float[tupleName+"dr04HcalD1Iso"] );
		  info_->setFloat( "dr04HcalD2Iso",_float[tupleName+"dr04HcalD2Iso"] );
		}
	      else
		{
		}
	      cand_->setInfo( info_ );
	      cand_->lock();
		  
	      leptons[lmode][index] = cand_;

	      //	  cand_->oneLine(cout );
	    }
	}
	  
      for( int lmode=0; lmode<2; lmode++ )
	{
	  tupleName = ZllName[lmode]+"_";
	  for( int iZll=0; iZll<nZll[lmode]; iZll++ )
	    {
	      ZllTuple[lmode]->GetEntry( jZll[lmode]++ );
	      assert( _int[tupleName+"run"]==run );
	      assert( _int[tupleName+"event"]==event );
	      assert( _int[tupleName+"mode"]==lmode );

	      int iDau0 = _int[tupleName+"dau0"];
	      int iDau1 = _int[tupleName+"dau1"];
	      assert( leptons[lmode].count(iDau0)!=0 && 
		      leptons[lmode].count(iDau1)!=0 ); 

	      int ihyp = _int[tupleName+"hyp"];
	      int iZ   = _int[tupleName+"iZ"];

	      Candidate* cand_ = Candidate::create( leptons[lmode][iDau0], 
						    leptons[lmode][iDau1] );
	      if( _verbose )
		{
		  cout << "ihyp/iZ  --> M/qT/eta/phi "
		       << ihyp << "/" << iZ << " --> "
		       << cand_->mass() << "/" 
		       << cand_->pt() << "/" 
		       << cand_->eta() << "/" 
		       << cand_->phi()*Constants::radToDeg << endl; 
		}
 
	      int pdgId_ = 23;
	      cand_->setPdgCode( pdgId_ );
	      cand_->setMass( _float[tupleName+"mll"] );

	      CandInfo* info_ = CandInfo::create();
	      info_->setFloat( "qTll", _float[tupleName+"qTll"] );
	      info_->setFloat( "yll", _float[tupleName+"yll"] );
	      info_->setFloat( "hll", _float[tupleName+"hll"] );
	      info_->setBool( "isTight", _bool[tupleName+"isTight"] );
	      info_->setInt( "hyp", ihyp );
	      info_->setInt( "iZ", iZ );

	      cand_->setInfo( info_  );
	  
	      // OK store Zll candidate
	      cand_ -> lock();

	      if( _verbose )
		{
		  cand_->oneLine( cout );
		  info_->print( cout );
		}
	    }
	}



      // cleaning
      Candidate::reset();
      ObjectStore::clear();
    }
}
  
int main()
{
  DummyTupleAnalysis dum;
  dum.analyze();
  return 0;
}



