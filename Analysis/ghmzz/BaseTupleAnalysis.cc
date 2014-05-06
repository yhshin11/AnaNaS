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

#include "Analysis/ghmzz/BaseTupleAnalysis.hh"

#include "Analysis/core/Sample.hh"
#include "Analysis/core/Candidate.hh"
#include "Analysis/core/CandInfo.hh"
#include "Analysis/tools/CandUtil.hh"
#include "Analysis/tools/PrintTree.hh"
#include "Analysis/utils/Constants.hh"
#include "Analysis/utils/Config.hh"
#include "Analysis/utils/KineUtils.hh"
#include "Analysis/utils/RooUtils.hh"
#include "Analysis/utils/HistoManager.hh"
#include "Analysis/utils/TupleManager.hh"
#include "Analysis/utils/DatasetManager.hh"
#include "Analysis/utils/CutUtils.hh"

string BaseTupleAnalysis::ZllName[2] = { "Zee", "Zmm" };

ClassImp( BaseTupleAnalysis )

BaseTupleAnalysis::BaseTupleAnalysis( const char* sampleSetName, int n ) : 
  _verbose(false), n_max(n), _lumi(0), dirName("Zll"),
  run(0), event(0), categ("data"), categ_long("data"), genCand(0), MetCand(0), 
  nLep(2), leptons(2), nZll(2), Zll(2)
{
  _sampleSetName = "data2011";
  if( sampleSetName!=0 ) _sampleSetName = string(sampleSetName);
  size_t found(0);
  if( _sampleSetName=="data2011" )
    {
      _sampleSet.push_back( "data2011A_DoubleElectron" );
      _sampleSet.push_back( "data2011B_DoubleElectron" );
      _sampleSet.push_back( "data2011A_DoubleMu" );
      _sampleSet.push_back( "data2011B_DoubleMu" );
      _sampleSet.push_back( "data2011A_MuEG" );
      _sampleSet.push_back( "data2011B_MuEG" );
    }
  else if( _sampleSetName.find("data2011_")!=string::npos )
    {
      string str_=_sampleSetName;
      found = str_.find("_");      
      str_.insert(found,"A"); 
      _sampleSet.push_back( str_ );
      str_.replace(found,1,"B");
      _sampleSet.push_back( str_ );
    }
  else
    _sampleSet.push_back( _sampleSetName );

  cout << _sampleSetName << endl;
  for( size_t ii=0; ii<_sampleSet.size(); ii++ )
    {
      cout << "---> " << _sampleSet[ii] << endl;
    }
}

void
BaseTupleAnalysis::go()
{
  init_();
  eventLoop();
  finalize_();
}


void
BaseTupleAnalysis::init_()
{
  setOutputFile();
  defineSamples();
  init();
}

bool 
BaseTupleAnalysis::analyzeEvent()
{
  bool isSelected=false;
  return isSelected;
}

void
BaseTupleAnalysis::finalize_()
{
  finalize();
  for( map<string,int>::iterator it=nAnalyzed.begin();
       it!=nAnalyzed.end(); ++it )
    {
      cout << "nSelected/nAnalyzed[" << it->first << "]=" 
	   << nSelected[it->first] 
	   << "/" 
	   << nAnalyzed[it->first] 
	   << endl;
    }

  hm_out.save();
  tm_out.save();
}

void
BaseTupleAnalysis::init()
{
}

void
BaseTupleAnalysis::finalize()
{
}

void
BaseTupleAnalysis::eventLoop()
{
  string sampleAnalysis = "Zll";

  for( SampleIt it=sampleList.begin(); it!=sampleList.end(); it++ )
    {
      _sampleName = it->first; // this will be the name of the directory
      _sample = it->second;       
      _subSampleName = _sample->name();

      string filename = Config::rootPath;
      // string filename = "/home/gpfs/manip/mnt/cms/data/root2011/";
      filename += sampleAnalysis;
      filename += "/";
      filename += _subSampleName;
      filename += ".root";

      TFile* inputFile  = TFile::Open( filename.c_str(), "READ" );
      if( inputFile==0 ) continue;
      inputFile->Print();

      // tuple manager
      tm.reset();
      tm.setFile( inputFile, TupleManager::kRead );

      // get the HLT line names for this file
      hltLines.clear();
      TTree* hltLineTuple = tm.getTree( "HltLines", dirName );
      int nHltLine   = hltLineTuple->GetEntriesFast();  
      {
	int lineIndex_;
	string lineName_;
	string* lineName_ptr = &lineName_;
	tm.add<    int   > ( "lineIndex", &lineIndex_ );
	tm.add<  string*  > ( "lineName", &lineName_ptr );
	for( int ihlt=0; ihlt<nHltLine; ihlt++ )
	  {
	    hltLineTuple->GetEntry(ihlt);
	    hltLines[lineIndex_] = lineName_;
	  }
      }

      TTree* eventTuple = tm.getTree( "EventTuple", dirName );
      vector< int >* ilines_ = new vector< int >; 
      tm.add< vector<int>* > ( "lines", &ilines_  );
      {
	declareInt(    "run"     );
	declareInt(    "event"   );
	
	// counters
	declareInt(    "leptonVertices"   );
	declareInt(    "Electron"   );
	declareInt(    "Muon"   );
	declareInt(    "Zee"   );
	declareInt(    "Zmm"   );
	declareInt(    "CleanJets"   );
	declareInt(    "MCTruth"   );

	// other stuff
	declareFloat(  "MET" );
	declareFloat(  "METPhi" );
	declareFloat(  "sumEt" );
	declareFloat(  "mEtSig" );
	declareFloat( "rhoFJ" );
	declareInt(      "PV_id" );
	declareFloat(     "PV_X" );
	declareFloat(     "PV_Y" );
	declareFloat(     "PV_Z" );
	declareInt(   "nZHypo_0" );
	declareInt(   "nZHypo_1" );
	declareInt(    "nVertex" );
	declareInt(    "nGoodVertices" );	
      }

      categ      = "data";
      categ_long = "data";
      string* categ_ptr = &categ;
      string* categ_long_ptr = &categ_long;
      TTree* MCEventTuple(0);
      TTree* MCTruthTuple(0);
      if( isMC() )
	{
	  MCEventTuple = tm.getTree( "MCEventTuple", dirName );
	  if( MCEventTuple==0 ) break;	  
	  {
	    declareInt(      "run" );
	    declareInt(    "event" );
	    tm.add<  string*  > ( "categ", &categ_ptr );
	    tm.add<  string*  > ( "categ_long", &categ_long_ptr );
	  }

	  MCTruthTuple = tm.getTree( "MCTruth", dirName );
	  if( MCTruthTuple==0 ) break;
	  {
	    declareInt(      "run" );
	    declareInt(    "event" );
	    declareInt(  "counter" );
	    declareInt(  "mother" );
	    declareInt(  "pdgId" );
	    declareFloat(  "pt" );
	    declareFloat(  "eta" );
	    declareFloat(  "phi" );
	    declareFloat(  "mass" );
	    declareFloat(  "E" );
	    declareFloat(  "pz" );
	    declareFloat(  "y" );
	    declareInt(    "vtxId" );
	    declareFloat(  "VX" );
	    declareFloat(  "XY" );
	    declareFloat(  "VZ" );
	  }
	}

      TTree* vtxTuple = tm.getTree( "leptonVertices", dirName );
      {
	declareInt(    "run"     );
	declareInt(    "event"   );
	declareInt(    "counter"   );
	declareInt(    "id"   );
	declareFloat( "X" );
	declareFloat( "Y" );
	declareFloat( "Z" );
      }

      TTree* jetTuple = tm.getTree( "CleanJets", dirName );
      {
	declareInt(    "run"     );
	declareInt(    "event"   );
	declareInt(    "counter"   );
	declareFloat(     "Et" );
	declareFloat(    "eta" );
	declareFloat(    "phi" );
	declareInt(    "index" );
	declareInt(    "vtxId" );
	declareFloat( "area" );
	declareFloat( "chargedEmEnergyFraction" );
	declareFloat( "chargedHadronEnergyFraction" );
	declareFloat( "chargedMuEnergyFraction" );
	declareFloat( "chargedMultiplicity" );
	declareFloat( "etaPhiMom" );
	declareFloat( "maxD" );
	declareFloat( "muonMultiplicity" );
	declareFloat( "neutralEmEnergyFraction" );
	declareFloat( "neutralHadronEnergyFraction" );
	declareFloat( "neutralMultiplicity" );
	declareFloat( "phiPhiMom" );
	declareBool( "isCleanJet" );
	declareBool( "CJisCharged" );
	declareBool( "CJisIsolatedLepton" );
	declareFloat( "CJsumPt" );
      }

      TTree* leptonTuple[2];
      string leptonName[2] = { "Electron", "Muon" };
      for( int lmode=0; lmode<2; lmode++ )
	{
	  leptonTuple[lmode] = tm.getTree( leptonName[lmode], dirName );
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
	      declareBool(  "isGlobal" );
	      declareBool(  "muIdTracker" );
	      declareFloat( "normChi2" );
	      declareInt(   "nVHits" );
	      declareInt(   "nMatch" );
	      declareInt(   "nPixHit" );
	      declareInt(   "nMuonHit" );
	      declareFloat( "isoR03strk" );
	      declareFloat( "isoR03sem" );
	      declareFloat( "isoR03shad" );
	      declareFloat( "dPV" );
	    }
	}
  
      TTree* ZllTuple[2];
  
      for( int lmode=0; lmode<2; lmode++ )
	{
	  ZllTuple[lmode] = tm.getTree( ZllName[lmode], dirName );
      
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
      if( n_max!=-1 && n>n_max ) n=n_max; 
 
      int jVtx=0;
      int jJet=0;
      vector<int> jLep(2,0);
      vector<int> jZll(2,0);
      int jMC=0;
      
      //
      // OK -- That is the event list for this sample
      //
      for( int ii=0; ii<n; ii++ )
	{
	  // load the next event
	  eventTuple->GetEntry( ii );
	  
	  // increment event counter
	  nAnalyzed[_sampleName]++;
	  if( nAnalyzed[_sampleName]%50000==0 )
	    {
	      cout << "nEvent[" << _sampleName << "]=" << nAnalyzed[_sampleName] << endl;
	    }

	  // event initializations
	  run = 0;
	  event = 0;
	  genCand=0;
	  nHLT=0;
	  nVtx=0;
	  nJet=0;
	  nMC=0;
	  vertices.clear();
	  leptonAtVertex.clear();
	  caloOnlyJets.clear();
	  jetAtVertex.clear();
	  ZAtVertex.clear();
	  for( int lmode=0; lmode<2; lmode++ )
	    {
	      nLep[lmode] = 0;
	      nZll[lmode] = 0;
	      leptons[lmode].clear();
	      Zll[lmode].clear();
	    }
	  MetCand=0;
	  firedHltLines.clear();
	  // OK
	  
	  string tupleName = "EventTuple_";
	  irun   =  _int[tupleName+"run"];
	  ievent =  _int[tupleName+"event"];
	  run   = (unsigned long) irun;
	  event = (unsigned long) ievent;

	  nVtx  = _int[tupleName+"leptonVertices"];      
	  for( int lmode=0; lmode<2; lmode++ )
	    {
	      nLep[lmode] = _int[tupleName+leptonName[lmode]];
	      nZll[lmode] = _int[tupleName+ZllName[lmode]];
	    }
	  nJet  = _int[tupleName+"CleanJets"];
	  nHLT  = (int) ilines_->size();

	  if( MCEventTuple!=0 )
	    {
	      MCEventTuple->GetEntry(ii);
	      nMC = _int["EventTuple_MCTruth"];
	      tupleName = "MCTruth_";
	      CandList     _mcPrimary;
	      map<int,Candidate*>     _mcCand;
	      map< int, vector<int> > _dauLinks;
	      map< int, TVector3 > _vtxP3;
	      for( int iMC=0; iMC<nMC; iMC++ )
		{
		  MCTruthTuple->GetEntry( jMC++ );
		  assert( ((unsigned long)_int[tupleName+"run"])==run );
		  assert( ((unsigned long)_int[tupleName+"event"])==event );
		 
		  int pdgId_ = _int[tupleName+"pdgId"];
		  TVector3 p3_(0,0,0);
		  float pt_   = _float[tupleName+"pt"];
		  float eta_  = _float[tupleName+"eta"];
		  float phi_  = _float[tupleName+"phi"];
		  float mass_ = _float[tupleName+"mass"];
		  float E_    = _float[tupleName+"E"];
		  float pz_   = _float[tupleName+"pz"];

		  //		  float y_    = _float[tupleName+"y"];
		  if( pt_==0 )
		    {
		      float E2_ = E_*E_;
		      float M2_ = mass_*mass_;
		      float p2_ = E2_;
		      if( M2_<E2_ ) p2_ = E2_-M2_;
		      //		      else
		      //			cout << "Warning : M2/E2 " << M2_ << "/" << E2_ << endl; 
		      //		      float p_ = sqrt( p2_ );
		      float pt2_ = p2_ - pz_*pz_;
		      if( pt2_>0 ) pt_ = sqrt( pt2_ );
		      else
			pt_ = 1.e-05;
		      //			cout << "Warning : p/pz " << p_ << "/" << pz_ << endl;
		    }
		  p3_.SetPtEtaPhi( pt_, eta_, phi_ ); 
		  Candidate* cand_= Candidate::create( p3_, pdgId_ ); 
		  cand_->setMass( mass_ );
	      
		  CandInfo* info_ = CandInfo::create();
		  cand_->setInfo( info_  );
		  
		  int counter = _int[tupleName+"counter"];
		  int mother  = _int[tupleName+"mother"];
		  _mcCand[counter] = cand_;
		  if( mother==0 ) 
		    _mcPrimary.push_back(cand_);
		  else
		    {
		      _dauLinks[mother].push_back(counter);
		    }
		  
		  float VX_   = _float[tupleName+"VX"];
		  float VY_   = _float[tupleName+"XY"]; //!!!!
		  float VZ_   = _float[tupleName+"VZ"]; 
		  TVector3 vtxP3_( VX_, VY_, VZ_ );
		  if( _vtxP3.count(mother)==0 )
		    {
		      _vtxP3[mother] = vtxP3_;
		    }
		} 	      	
	      for( map< int, vector<int> >::iterator it_=_dauLinks.begin();
		   it_!=_dauLinks.end(); it_++ )
		{
		  int mother_ = it_->first;
		  assert( mother_>0 );
		  int vtxId_(-1);
		  Vertex* vtx_(0);
		  if( _vtxP3.count(mother_) )
		    {
		      vtxId_ = mother_;
		      vtx_ = Vertex::create( _vtxP3[mother_] );
		      CandInfo* info_ = CandInfo::create();
		      info_->setInt( "index", mother_ );
		      vtx_->setInfo( info_  );
		    }		  
		  Candidate* theMother_ = _mcCand[mother_];
		  assert( theMother_!=0 );
		  vector<int>& dauCounter_ = it_->second;
		  for( size_t idau=0; idau<dauCounter_.size(); idau++ )
		    {
		      Candidate* theDaughter_ = _mcCand[ dauCounter_[idau] ];
		      assert( theDaughter_!=0 );
		      theMother_->addDaughter( theDaughter_ );
		      theDaughter_->setMother( theMother_ );
		      theDaughter_->setVertex( vtx_ );
		    }
		}
	      assert( _mcPrimary.size()!=0 );
	      genCand = Candidate::create( _mcPrimary );
	      genCand->setName("HardScatter");
	    }
		  
	  for( size_t il=0; il<ilines_->size(); il++ )
	    {
	      int i_ = (*ilines_)[il];
	      firedHltLines.push_back( hltLines[i_] );
	    }


	  // the MET candidate
	  {
	    TVector2 p2(0.,0.);
	    float Et_  = _float["EventTuple_MET"];
	    float phi_ = _float["EventTuple_METPhi"];
	    p2.Set( Et_*cos(phi_), Et_*sin(phi_) );      
	    Candidate* cand_ = Candidate::create( p2 );
	    CandInfo* info_ = CandInfo::create();
	    cand_->setInfo( info_ );
	    
	    info_->setFloat( "sumEt", _float["EventTuple_sumEt"] );
	    info_->setFloat( "mEtSig", _float["EventTuple_mEtSig"] );
	    
	    //
	    MetCand = cand_;
	    MetCand->setName( "MET" );
	  }

	  tupleName = "leptonVertices_";
	  int PV_id = _int["EventTuple_PV_id"];
	  for( int iVtx=0; iVtx<nVtx; iVtx++ )
	    {
	      vtxTuple->GetEntry( jVtx++ );
	      assert( ((unsigned long)_int[tupleName+"run"])==run );
	      assert( ((unsigned long)_int[tupleName+"event"])==event );

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
	      if( index==PV_id ) vtx_->setAsPrimary();

	      vertices[index] =  vtx_;
	    }

	  // the primary vertex is not a lepton vertex...
	  {
	    if( vertices.count(PV_id)==0 )
	      {
		int index = PV_id;
		Vertex* vtx_ = Vertex::create();
		vtx_->setXYZ( _float["EventTuple_PV_X"],
			      _float["EventTuple_PV_Y"], 
			      _float["EventTuple_PV_Z"] );
		
		// special track information
		CandInfo* info_ = CandInfo::create();
		info_->setInt( "index", index );
		vtx_->setInfo( info_  );
		vtx_->setAsPrimary();
		
		vertices[index] =  vtx_;
	      }
	  }

	  tupleName = "CleanJets_";
	  for( int iJet=0; iJet<nJet; iJet++ )
	    {
	      jetTuple->GetEntry( jJet++ );
	      assert( ((unsigned long)_int[tupleName+"run"])==run );
	      assert( ((unsigned long)_int[tupleName+"event"])==event );
		  
	      int pdgId_ = 21;
	      TVector3 p3_(0,0,0);
	      p3_.SetPtEtaPhi( 
			      fabs( _float[tupleName+"Et"] ), 
			      _float[tupleName+"eta"], 
			      _float[tupleName+"phi"] );
	      int index = _int[tupleName+"index"];
	      int vtxId = _int[tupleName+"vtxId"]; 
	      Vertex* vtx(0);
	      if( vertices.count(vtxId)!=0 )
		{
		  vtx = vertices[vtxId];
		}
	      else if( vtxId!=-1 )
 		{
 		  cout << "WARNING vtxid=" << vtxId <<endl;
 		  cout << " -- vectices: " << vertices.size() << endl;
 		  for( map<int,Vertex*>::iterator ij = vertices.begin(); 
 		       ij!=vertices.end(); ij++ )
 		    {
 		      cout << " -- index: " << ij->first << endl;
 		      ij->second->print(cout);
 		    }
		  assert(0);
 		}

	      Candidate* cand_= Candidate::create( p3_, pdgId_, vtx ); 
	      
	      // special jet information
	      CandInfo* info_ = CandInfo::create();
	      cand_->setInfo( info_  );

	      info_->setInt(  "index", index );
	      info_->setFloat( "area", _float[tupleName+"area"] );
	      info_->setFloat( "chargedEmEnergyFraction", _float[tupleName+"chargedEmEnergyFraction"] );
	      info_->setFloat( "chargedHadronEnergyFraction", _float[tupleName+"chargedHadronEnergyFraction"] );
	      info_->setFloat( "chargedMuEnergyFraction", _float[tupleName+"chargedMuEnergyFraction"] );
	      info_->setFloat( "chargedMultiplicity", _float[tupleName+"chargedMultiplicity"] );
	      info_->setFloat( "etaPhiMom", _float[tupleName+"etaPhiMom"] );
	      info_->setFloat( "maxD", _float[tupleName+"maxD"] );
	      info_->setFloat( "muonMultiplicity", _float[tupleName+"muonMultiplicity"] );
	      info_->setFloat( "neutralEmEnergyFraction", _float[tupleName+"neutralEmEnergyFraction"] );
	      info_->setFloat( "neutralHadronEnergyFraction", _float[tupleName+"neutralHadronEnergyFraction"] );
	      info_->setFloat( "neutralMultiplicity", _float[tupleName+"neutralMultiplicity"] );
	      info_->setFloat( "phiPhiMom", _float[tupleName+"phiPhiMom"] );
	      info_->setFloat( "CJsumPt", _float[tupleName+"CJsumPt"] );

	      bool isCleanJet_       = _bool[tupleName+"isCleanJet"];
	      bool isIsolatedLepton_ = _bool[tupleName+"CJisIsolatedLepton"];
	      
	      info_->setBool(  "isCleanJet", isCleanJet_ );
	      //	      info_->setBool(  "CJisCharged", _bool[tupleName+"CJisCharged"] );
	      info_->setBool(  "CJisIsolatedLepton", isIsolatedLepton_ );

	      if( vtxId==-1 )
		{
		  caloOnlyJets.push_back( cand_ );
		}
	      else if( isCleanJet_ && isIsolatedLepton_ )
		{
		  jetAtVertex[vtxId].push_back( cand_ );
		}
	    }

	  for( int lmode=0; lmode<2; lmode++ )
	    {
	      tupleName = leptonName[lmode]+"_";
	      for( int iLep=0; iLep<nLep[lmode]; iLep++ )
		{
		  leptonTuple[lmode]->GetEntry( jLep[lmode]++ );
		  assert( ((unsigned long)_int[tupleName+"run"])==run );
		  assert( ((unsigned long)_int[tupleName+"event"])==event );

		  int pdgId_ = 11;
		  if( lmode==1 ) pdgId_ = 13;
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
		  
		  // special information
		  CandInfo* info_ = CandInfo::create();
		  int id_ = _int[tupleName+"id"];
		  info_->setInt( "id", id_ );
		  bool isVeryLoose = id_>=1;
		  bool isLoose     = id_>=3;
		  bool isTight     = id_>=7;
		  bool isVeryTight = id_==15;
		  info_->setBool("isVeryLoose", isVeryLoose );
		  info_->setBool("isLoose",     isLoose     );
		  info_->setBool("isTight",     isTight     );
		  info_->setBool("isVeryTight", isVeryTight );
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
		      info_->setBool(  "isGlobal", _bool[tupleName+"isGlobal"]  );
		      info_->setBool(  "muIdTracker", _bool[tupleName+"muIdTracker"] );
		      info_->setFloat( "normChi2", _float[tupleName+"normChi2"] );
		      info_->setInt(   "nVHits", _int[tupleName+"nVHits"] );
		      info_->setInt(   "nMatch", _int[tupleName+"nMatch"] );
		      info_->setInt(   "nPixHit", _int[tupleName+"nPixHit"] );
		      info_->setInt(   "nMuonHit", _int[tupleName+"nMuonHit"] );
		      info_->setFloat( "isoR03strk", _float[tupleName+"isoR03strk"] );
		      info_->setFloat( "isoR03sem", _float[tupleName+"isoR03sem"] );
		      info_->setFloat( "isoR03shad", _float[tupleName+"isoR03shad"] );
		      info_->setFloat( "dPV", _float[tupleName+"dPV"] );

		    }
		  cand_->setInfo( info_ );
		  cand_->lock();
		  
		  leptons[lmode][index] = cand_;
		  leptonAtVertex[vtxId].push_back( cand_ );
		  //	  cand_->oneLine(cout );
		}
	    }
	  
	  for( int lmode=0; lmode<2; lmode++ )
	    {
	      tupleName = ZllName[lmode]+"_";
	      for( int iZll=0; iZll<nZll[lmode]; iZll++ )
		{
		  ZllTuple[lmode]->GetEntry( jZll[lmode]++ );
		  assert( ((unsigned long)_int[tupleName+"run"])==run );
		  assert( ((unsigned long)_int[tupleName+"event"])==event );
		  assert( _int[tupleName+"mode"]==lmode );

		  int iDau0 = _int[tupleName+"dau0"];
		  int iDau1 = _int[tupleName+"dau1"];
		  assert( leptons[lmode].count(iDau0)!=0 && 
			  leptons[lmode].count(iDau1)!=0 ); 

		  int ihyp = _int[tupleName+"hyp"];
		  int iZ   = _int[tupleName+"iZ"];

		  Candidate* cand_ = Candidate::create( leptons[lmode][iDau0], 
							leptons[lmode][iDau1] );
		  int pdgId_ = 23;
		  cand_->setPdgCode( pdgId_ );
		  cand_->setMass( _float[tupleName+"mll"] );

		  CandInfo* info_ = CandInfo::create();
		  info_->setFloat( "qTll", _float[tupleName+"qTll"] );
		  info_->setFloat( "yll", _float[tupleName+"yll"] );
		  info_->setFloat( "hll", _float[tupleName+"hll"] );
		  info_->setBool( "isTight", _bool[tupleName+"isTight"] );
		  info_->setInt( "mode", _int[tupleName+"mode"] );
		  info_->setInt( "hyp", ihyp );
		  info_->setInt( "iZ", iZ );

		  cand_->setInfo( info_  );
	  
		  // OK store Zll candidate
		  cand_ -> lock();

		  int vtx_ = cand_->vertexIndex();

		  Zll[lmode].push_back( cand_ );
		  ZAtVertex[vtx_].push_back( cand_ );
		}
	    }

	  bool isSelected = analyzeEvent();
	  if( isSelected ) nSelected[_sampleName]++;
	  
	  // event cleaning
	  Candidate::reset();
	  ObjectStore::clear();
	} // end of event loop

      delete ilines_;
      inputFile->Close();  
    
    }  // end of sample loop
}
// 	  if( hbool.size()==0 )
// 	    {
// 	      // bool histograms
// 	      hbool.push_back(
// 			      HBool(
// 				    hist(
// 					 "mll", "mll", dirName, "" 
// 					 ),            // the histogram
// 				    &_float["Zee_mll"],    // pointer to the variable
// 				    &_bool["Zee_isTight"], // pointer to the boolean
// 				    true,              // underflow bin?
// 				    true               // overflow bin?
// 				    )
// 			      );		
// 	    }
// 	  float w_ = 1;
// 	  for( size_t ih=0; ih<hbool.size(); ih++ )
// 	    {
// 	      hbool[ih].fill( w_ );
// 	    }

void
BaseTupleAnalysis::defineTemplate( string templ, 
				   int nbin, float min, float max )
{
  cout << "1D-Template[" << templ << "]-->(nbin=" << nbin << "/min=" << min << "/max=" << max << ")" << endl; 

  TH1F* h_ = new TH1F( templ.c_str(), templ.c_str(), 
		       nbin, min, max );
  h_->Sumw2();
  h_->SetMarkerSize(0);
  hm_out.addTemplate<TH1F>( templ, h_ );
}

void
BaseTupleAnalysis::defineTemplate( string templ, 
			       int nbin, float* xbins )
{
  cout << "1D-Template[" << templ << "]-->(nbin=" << nbin << "variable size" << endl; 
  TH1F* h_ = new TH1F( templ.c_str(), templ.c_str(), 
		       nbin, xbins );
  h_->Sumw2();
  hm_out.addTemplate<TH1F>( templ, h_ );
}

TH1* 
BaseTupleAnalysis::hist( string templ, string cut, string dir, string prefix )
{
  TH1* h_ = (TH1*)hm_out.h<TH1F>( templ, cut, dir, prefix );
  assert( h_!=0 );
  return h_;
}

void
BaseTupleAnalysis::fill( string templ, string cut, string dir, string prefix,
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
BaseTupleAnalysis::setOutputFile()
{
  string analysis = _sampleSetName;

  string histfile = Config::histPath + analysis;
  histfile += ".root";
  
  // output file 
  TFile* outputfile = TFile::Open( histfile.c_str(), "RECREATE" );
  
  // histo manager
  hm_out.setFile( outputfile );

  tm_out.setFile( outputfile, TupleManager::kWrite );
}

void
BaseTupleAnalysis::declareInt( string var )
{
  string tname  = tm.tree()->GetName();
  string varKey = tm.varKey( tname, var ); 
  tm.add<int>(    var,    &_int[varKey] );
}

void
BaseTupleAnalysis::declareFloat( string var )
{
  string tname  = tm.tree()->GetName();
  string varKey = tm.varKey( tname, var ); 
  tm.add<float>(    var,    &_float[varKey] );
}

void
BaseTupleAnalysis::declareBool( string var )
{
  string tname  = tm.tree()->GetName();
  string varKey = tm.varKey( tname, var ); 
  tm.add<bool>(    var,    &_bool[varKey] );
}

bool 
BaseTupleAnalysis::isData2011A() const 
{
  return _subSampleName.find("data2011A/")!=string::npos;
}

bool 
BaseTupleAnalysis::isData2011B() const
{ 
  return _subSampleName.find("data2011B/")!=string::npos;
}

bool 
BaseTupleAnalysis::isData() const
{      
  return isData2011A() || isData2011B();
}

bool 
BaseTupleAnalysis::isMC() const
{        
  return _subSampleName.find("Fall11/")!=string::npos;
}

bool 
BaseTupleAnalysis::isDoubleElectron() const  
{ 
  return  _subSampleName.find("DoubleElectron_")!=string::npos; 
}

bool BaseTupleAnalysis::isDoubleMu() const
{       
  return  _subSampleName.find("DoubleMu_")!=string::npos;
}

void
BaseTupleAnalysis::defineSamples()
{
  cout << "Entering defineSamples " << endl;

  // define the samples
  samples.clear();

  if( false )
    {
      //      _subSampleName = "data2011A/DoubleElectron_173000_173389";
      //      _subSampleName = "Fall11/Z_2e__1";
      _subSampleName = "Fall11/WZ_3ln";

//       Sample* sample_ = 
// 	new Sample( _subSampleName.c_str(), 
// 		    "Run2011A-16Jan2012-v1/AOD", 
// 		    "data",
// 		    109.44, 0.287, 1., 1., 231369 );
//       if( isData() ) sample_->setAsData();
      Sample* sample_ = 
	new Sample( _subSampleName.c_str(), 
		    "DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia/Fall11-PU_S6_START42_V14B-v1/AODSIM", 
		    "MC",
		    1700., 1.000, 1., 1., 493000 );

      cout 
	<< "n_proc=" << sample_->n_proc()
	<< " sigma=" << sample_->sigma()
	<< " Integrated luminosity = " 	
	<< sample_->integratedLumi()
	<< endl;

      samples[_subSampleName] = sample_; 
            
      //      sampleList.insert( SamplePair( "data2011B", samples[_subSampleName] ) );
      sampleList.insert( SamplePair( "Fall11_Z_2e", samples[_subSampleName] ) );

      return;
    }

  cout << "Reading the samples/txt file " << endl;
  TString sampleFile = Config::confPath + "samples.txt";

  ifstream fin;
  fin.open(sampleFile);
  
  int nline(0);
  char line_[512];
  char c_;

  while( fin.good() )
    {
      c_=fin.peek();
      if( c_=='/' )
	{      
	  fin.getline( line_, 512 );
	  continue;
	} 
      nline++;
      
      string collectionName;
      int   nEvt_(0);
      float sig_tot_(0);
      float filter_eff_(1);
      
      fin >> _subSampleName;
      fin >> nEvt_;
      fin >> sig_tot_;
      fin >> filter_eff_;
      fin >> collectionName;
      fin.getline( line_, 512 );
      string comments(line_);
      
      if( _subSampleName.length()==0 ) continue;
	

      //      if( isMC() ) continue;

      string sampleName_;
      if( isData() )
	{
	  size_t found1_(0), found2_(0);
	  found1_ = _subSampleName.find("_");
	  found2_ = _subSampleName.rfind("_");
	  string firstrun_;
	  string lastrun_;
	  assert( found1_!=string::npos && found2_!=string::npos 
		  && found2_>found1_ ); 
	  sampleName_ = _subSampleName.substr(0,found1_);
	  //	  firstrun_   = _subSampleName.substr(found1_+1,found2_-found1_-1);
	  //	  lastrun_   = _subSampleName.substr(found2_+1);
	  size_t found_ = sampleName_.find("/");
	  if( found_!=string::npos )
	    sampleName_.replace(found_,1,"_");
	  

	  //	  string period = isData2011A() ? "data2001A" : "data2011B";

	  float lumi_ = sig_tot_; 
	  sig_tot_ = nEvt_/lumi_/filter_eff_;

	}
      else if( isMC() )
	{
	  size_t found_(0);
	  found_ = _subSampleName.find("__");
	  sampleName_ = _subSampleName;
	  if( found_!=string::npos ) 
	    sampleName_ = _subSampleName.substr(0,found_);
	  found_ = sampleName_.find("/");
	  if( found_!=string::npos )
	    sampleName_.replace(found_,1,"_");
	}

      int mycount = (int) count( _sampleSet.begin(), _sampleSet.end(), sampleName_ );
      if( mycount==0 ) continue;	
            
      cout << "SubSampleName : " << _subSampleName << endl;
      cout << " ---> sampleName " << sampleName_ << endl;

      Sample* sample_ = 
	new Sample( _subSampleName.c_str(), 
		    collectionName.c_str(), 
		    comments.c_str(),
		    sig_tot_, filter_eff_, 1., 1., nEvt_ );
      if( isData() ) sample_->setAsData();

      if( sample_->n_proc()==0 ) continue;

      cout 
	<< "n_proc=" << sample_->n_proc()
	<< " sigma=" << sample_->sigma()
	<< " Integrated luminosity = " 	
	<< sample_->integratedLumi()
	<< endl;

      if( isData() )
	_lumi += sample_->integratedLumi();

      samples[_subSampleName] = sample_; 
            
      sampleList.insert( SamplePair( sampleName_, 
				     samples[_subSampleName] ) );
    }

  for( SampleIt it=sampleList.begin(); it!=sampleList.end(); it++ )
    {
      string sampleName = it->first; // this will be the name of the directory
      Sample* sample_ = it->second;  
      //      if( !sample_->isData() ) continue;
      integratedLumi[sampleName] += sample_->integratedLumi();
    }  

  for( map< string, float >::iterator it=integratedLumi.begin();
       it!=integratedLumi.end(); ++it )
    {
      cout << "Lumi[" << it->first << "]=" << it->second << endl;
    }  
}

void
BaseTupleAnalysis::print()
{
  cout << "\n+++++++++++++++++++++++++++++++++" << endl;
 
    cout << "run/event: " << run << "/" << event << endl;      
    if( genCand!=0 )
      {
	cout << " ---> " << categ << "/" << categ_long << endl; 
	genCand->oneLine(cout);
	PrintTree prtTree;
	cout << prtTree.Print( *genCand );
      }
    
    cout << "Counters" << endl;
    {
      cout << "\tHLT        --> " <<  nHLT  << endl;
      cout << "\tVertices   --> " <<  nVtx  << endl;
      cout << "\tCleanJets  --> " <<  nJet  << endl;
      cout << "\tElectrons  --> " <<  nLep[0]  << endl;
      cout << "\tMuons      --> " <<  nLep[1]  << endl;
      cout << "\tZee        --> " <<  nZll[0]  << endl;
      cout << "\tZmm        --> " <<  nZll[1]  << endl;     
    }

    if( firedHltLines.size()>0 )
      {
	cout << "Fired HLT Lines --> " << firedHltLines.size() << endl;
	for( size_t il=0; il<firedHltLines.size(); il++ )
	  {
	    cout << firedHltLines[il] << ":";
	  }
	cout << endl;
      }

    cout << "Missing ET" << endl;
    cout << "\t";
    MetCand->oneLine( cout );

    for( int lmode=0; lmode<2; lmode++ )
      {
	for( size_t iZ=0; iZ<Zll[lmode].size(); iZ++ )
	  {
	    Candidate* cand_ = Zll[lmode][iZ];
	    CandInfo* info_ = cand_->info();
	    cand_->oneLine( cout );
	    info_->print( cout );
	  }
      }

    for( map<int,CandList>::iterator it=leptonAtVertex.begin();
	 it!=leptonAtVertex.end(); ++it )
      {
	cout << "Vertex ID = " << it->first << "\t";
	vertices[it->first]->print(cout);
	CandUtil::printList( cout, it->second );
      }

    cout << "=================================\n" << endl;
}



