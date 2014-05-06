#include "Analysis/src/ZllAnalysis.hh"

#include <vector>
#include <map>
#include <list>
#include <algorithm>
#include <cassert>
#include <sstream>
using namespace std;

#include <TH1I.h>
#include <TH1F.h>
#include <TH2I.h>
#include <TH2F.h>

#include "Analysis/utils/Constants.hh"
#include "Analysis/utils/KineUtils.hh"
#include "Analysis/core/Candidate.hh"
#include "Analysis/core/CandInfo.hh"
#include "Analysis/core/EventManager.hh"
#include "Analysis/tools/CandUtil.hh"
#include "Analysis/selectors/IsolationSelector.hh"

ClassImp( ZllAnalysis )

typedef vector<double> vectorFloat; 

ZllAnalysis::ZllAnalysis( Sample& sample, EventManager& manager ) 
  : SampleAnalysis( "Zll", sample, manager )
{
  cout << "\t-------------------------------------------" << endl; 
  cout << "\t---            Zll analysis            ----" << endl; 
  cout << "\t-------------------------------------------" << endl; 

  _nDebug = 0;

}

ZllAnalysis::~ZllAnalysis()
{
}

bool
ZllAnalysis::analyzeEvent()
{
  //
  // Analysis of events with at least one Z boson
  //
  float rhoFJ = _e.getRhoFastJet();
  
  // if isMC, "categ" describes the "true" event
  bool isMC = EventServer::isMC;
  string categ      = _e.categ();
  string categ_long = _e.categ_long();
  string* categ_ptr = &categ;
  string* categ_long_ptr = &categ_long;

  // run and event number
  int run   = _e.run();
  int event = _e.event();

  if( _verbose )
    {
      cout << "run/event " << run << "/" << event << endl;
    }
  
  const VtxList& vertices = _e.vertices();
  int nVertex = vertices.size();
  int nGoodVertices = _e.nGoodVertices();

  // the current primary vertex
  const Vertex& PVertex = _e.primaryVertex();
  int PV_id   = PVertex.index();
  float PV_X  = PVertex.pos().X();
  float PV_Y  = PVertex.pos().Y();
  float PV_Z  = PVertex.pos().Z();
  
  // build the lists of electrons and muons (see core/SampleAnalysis.cc)
  buildLeptonLists();

  // build the combinatorics of Z candidates 
  // (Z-hypothesis = sorted lists of non-overlapping Z candidates and leptons)
  for( int lmode=0; lmode<2; lmode++ )
    {
      nZHypos[lmode] = buildZllLists( lmode );
    }

  // preselection
  //  bool select = true;
  if( nZHypos[0]==0 && nZHypos[1]==0 ) 
    {
      return false;
    }

  vector< CandList > leptonList(2);
  // list of leptons
  for( size_t lmode=0; lmode<2; lmode++ )
    leptonList[lmode] = _e.userCandList( leptonListName(lmode) );
  
  // the missing transverse momentum 
  Candidate* metCand =  _e.met( EventManager::kPfMet )->clone();
  assert( metCand!=0 );
  float MET    = metCand->pt();
  float METPhi = metCand->phi();
  CandInfo* metInfo = metCand->info();
  float sumEt  = metInfo->getFloat( "sumEt" );
  float mEtSig = metInfo->getFloat( "mEtSig" );

  // put all leptons together in a list
  CandList allLeptonList;
  CandUtil::mergeLists( leptonList[0], leptonList[1], allLeptonList );
  //  CandUtil::printList( cout, allLeptonList );
  
  // dispatch the list
  vector< CandList > allLeptonPerVertex;
  CandUtil::dispatchList( allLeptonList, allLeptonPerVertex );

  if( _verbose )
    {
      cout << "All lepton list " << endl; 
      CandUtil::printList( cout, allLeptonList );
      cout << "All lepton list per vertex " << endl; 
      for( size_t ii=0; ii<allLeptonPerVertex.size(); ii++ ) 
	{
	  TString str_("vertex_"); str_+=ii;
	  CandUtil::printList( cout, allLeptonPerVertex[ii], str_ );
	}
    }

  // map between vertices and leptons  
  // list of vertices to consider in the event
  vector< Vertex* > leptonVertices;
  //  vertices.push_back( PVertex );

  CandList cleanJets = _e.userCandList("cleanJets");      
  if( _verbose )
    {
      cout << "Clean jets at PV " << endl; 
      CandUtil::printList( cout, cleanJets );
    }

  for( size_t ilist=0; ilist<allLeptonPerVertex.size(); ilist++ )
    {
      CandList& list_ = allLeptonPerVertex[ilist];

      if( list_.size()==0 ) continue;

      Vertex* vtx = list_[0]->vertex();
      leptonVertices.push_back( vtx );
      
      if( _verbose )
	{
	  cout << "Lepton jets at vertex " << vtx->index() << endl; 
	  CandUtil::printList( cout, _e.userCandList("leptonJets") );
	}
      CandUtil::pruneList( _e.userCandList("leptonJets"), cleanJets, cleanJets );

      if( vtx->index() == PV_id ) continue;
        
      // this is a new vertex
      if( _verbose )
	cout << "New Lepton vertex " << vtx->index() << endl;

      _e.setPrimaryVertex( vtx );
      if( _verbose )
	{
	  cout << "Clean jet list at lepton vertex " << vtx->index() << endl; 
	  CandUtil::printList( cout, _e.userCandList("cleanJets") );
	}
      CandUtil::pruneList( _e.userCandList("cleanJets"), cleanJets, cleanJets );
      CandUtil::pruneList( _e.userCandList("leptonJets"), cleanJets, cleanJets );
    }
  if( _verbose )
    {
      cout << "Current clean jet list after pruning" << endl; 
      CandUtil::printList( cout, cleanJets, "", true );
    }

  string dir_(_name);
  int counter(0); 
  vector<string> tupleName;
  {
    tupleName.push_back("leptonVertices");
    tm.setTree(tupleName.back(),dir_);       
    counter = 0;
    for( size_t ivtx=0; ivtx<leptonVertices.size(); ivtx++ )
      {       
	Vertex* vtx_ = leptonVertices[ivtx];
	int vtx_id   = vtx_->index();
	float vtx_X  = vtx_->pos().X();
	float vtx_Y  = vtx_->pos().Y();
	float vtx_Z  = vtx_->pos().Z();
       
	declareInt(     "run", run        );
	declareInt(   "event", event      );
	declareInt( "counter", ++counter  );
	declareInt(      "id", vtx_id     );
	declareFloat(     "X", vtx_X      );
	declareFloat(     "Y", vtx_Y      );
	declareFloat(     "Z", vtx_Z      );
	tm.fill();
      }
  }
  
  {
    // loop on the list of candidate leptons
    for( size_t lmode=0; lmode<2; lmode++ )
      {
	tupleName.push_back(leptonName(lmode));
	tm.setTree(tupleName.back(),dir_);
	counter = 0;
	for( CandList::iterator it=leptonList[lmode].begin();
	     it!=leptonList[lmode].end(); ++it )
	  {	   
	    Candidate* LCand_ = *it;  
	    int id_(0);
	    CandInfo* info = LCand_->info();	   
	    int index = info->getInt("index");
	    declareInt(     "run", run );
	    declareInt(   "event", event );
	    declareInt( "counter", ++counter );
	    declareInt(   "index", index );
	    declareFloat(    "pt", LCand_->charge()*LCand_->pt() );
	    declareFloat(   "eta", LCand_->eta() );
	    declareFloat(   "phi", LCand_->phi() );
	    declareInt(   "vtxId", LCand_->vertexIndex() );

	    TVector3 v3_ = LCand_->pos();
	    float x0_  = 0;
	    float y0_  = 0;
	    float z0_  = 0;
	    if( !info->getFloat( "vx", x0_ ) ) x0_  = v3_.X();
	    if( !info->getFloat( "vy", y0_ ) ) y0_  = v3_.Y();
	    if( !info->getFloat( "vz", z0_ ) ) z0_  = v3_.Z();
	    declareFloat(    "x0", x0_ );
	    declareFloat(    "y0", y0_ );
	    declareFloat(    "z0", z0_ );

	    if( info->getBool("VeryLoose") ) id_ += 1;
	    if( info->getBool("Loose") )     id_ += 2;
	    if( info->getBool("Tight") )     id_ += 4;
	    if( info->getBool("VeryTight") ) id_ += 8;
	    declareInt( "id", id_ );

	    if( lmode==0 )
	      {
		float fBrem  = info->getFloat("fBrem"   );
		float EOP    = info->getFloat("EOP"     );
		float dPhiIn = info->getFloat("dPhiIn"  );
		float dEtaIn = info->getFloat("dEtaIn"  );
	      
		float dr03TrkIso = info->getFloat("dr03TrkIso");
		float dr04EcalIso = info->getFloat("dr04EcalIso");
		float dr04HcalD1Iso = info->getFloat("dr04HcalD1Iso"); 
		float dr04HcalD2Iso = info->getFloat("dr04HcalD2Iso");	  
	      
		bool isConv = info->getBool("isConv");
		Candidate* gsfTrk_ = _e.getSecond( "electron-GsfTrack", LCand_ );
		CandInfo* gsfInfo = gsfTrk_->info();
		assert( gsfInfo );
		bool isFPHitTrk = gsfInfo->getBool("isFPHitTrk");
		int  expInnHits = gsfInfo->getInt("expInnHits");
	      
		declareFloat( "fBrem",  fBrem  );
		declareFloat( "EOP",    EOP    );
		declareFloat( "dPhiIn", dPhiIn );
		declareFloat( "dEtaIn", dEtaIn );
		declareFloat( "dr03TrkIso", dr03TrkIso );
		declareFloat( "dr04EcalIso", dr04EcalIso );
		declareFloat( "dr04HcalD1Iso", dr04HcalD1Iso );
		declareFloat( "dr04HcalD2Iso", dr04HcalD2Iso );
		declareBool(  "isConv", isConv );
		declareBool(  "isFPHitTrk", isFPHitTrk );
		declareInt(   "expInnHits", expInnHits );

		// the closest cell
		string sc_mapname_ = "electron-SuperCluster";
		string prefix_     = "E-Gamma";
		Candidate* sc_ = _e.getSecond(  sc_mapname_, LCand_ );

		int ix_ = -999;
		int iy_ = -999;
		int iz_ = -999;
		if( sc_!=0 )
		  {
		    const CandInfo* sc_info_ = sc_->info();
		    ix_ = sc_info_->getInt( "ixClosestCell" );
		    iy_ = sc_info_->getInt( "iyClosestCell" );
		    iz_ = sc_info_->getInt( "izClosestCell" );
		  }
		declareInt(   "ix", ix_ );
		declareInt(   "iy", iy_ );
		declareInt(   "iz", iz_ );
	      }
	    else
	      {
		bool isGlobal(true);
		info->getBool("muidGlobalMuonPromptTight",isGlobal);
		bool muIdTracker(true);
		info->getBool("muidTrackerMuonArbitrated",muIdTracker);
		float normChi2 = info->getFloat("normChi2");
		int nVHits = info->getInt("nTrkHits"); 
		int nMatch = info->getInt("nMatch");
		int nPixHit = info->getInt("nPixHits");
		int nMuonHit = info->getInt("nMuonHits");
		float isoR03strk = info->getFloat("isoR03strk");
		float isoR03sem  = info->getFloat("isoR03sem");
		float isoR03shad = info->getFloat("isoR03shad");
		float dPV = info->getFloat("dPV"); //MM not dPV but d0
	      
		declareBool(  "isGlobal",  isGlobal  );
		declareBool(  "muIdTracker",  muIdTracker  );
		declareFloat( "normChi2",    normChi2    );
		declareInt(   "nVHits", nVHits );
		declareInt(   "nMatch", nMatch );
		declareInt(   "nPixHit", nPixHit );
		declareInt(   "nMuonHit", nMuonHit );
		declareFloat( "isoR03strk",    isoR03strk   );
		declareFloat( "isoR03sem",    isoR03sem    );
		declareFloat( "isoR03shad",    isoR03shad    );
		declareFloat( "dPV",    dPV    );
	      }
	    tm.fill();
	  }
      }
  }

  {
    for( size_t lmode=0; lmode<2; lmode++ )
      {
	tupleName.push_back(ZmodeName(lmode));
	tm.setTree(tupleName.back(),dir_);       
	counter=0;

	if( nZHypos[lmode]==0 ) continue;

	for( size_t ihyp=0; ihyp<nZHypos[lmode]; ihyp++ )
	  {
	    // list of Zee candidates
	    CandList& ZList  = _e.userCandList( ZListName(lmode,ihyp) );

	    for ( size_t iZ=0; iZ< ZList.size(); iZ++ )
	      {
		Candidate* ZCand = ZList[iZ];
		assert( ZCand!=0 );

		// get the vertex of the Z and make it the primary vertex
		Vertex* ZVertex = ZCand->vertex();
		assert( ZVertex!=0 );

		int ZV_id = ZVertex->index();
		float ZV_X  = ZVertex->pos().X();
		float ZV_Y  = ZVertex->pos().Y();
		float ZV_Z  = ZVertex->pos().Z();

		// the main Z variables
		float   mll = ZCand->mass();
		float  qTll = ZCand->pt();
		CandInfo* ZInfo = ZCand->info();
		float yll = ZInfo->getFloat("rapidity");
		float hll = ZInfo->getFloat("helicity");

		bool isTight(false);
		{
		  // selectZCand: see core/SampleAnalysis.cc
		  //     full DY window = [20,oo[
		  int nVT=0; // min number of VT leptons
		  int nT =1; // min number of  T leptons
		  int nL =2; // min number of  L leptons
		  int nVL=2; // min number of VL leptons    
		  isTight = selectZCand( *ZCand, 
					 Z0_fullDY[0], Z0_fullDY[1], 
					 nVT, nT, nL, nVL );
		}

		// get the daughter index
		int dauIdx[2];
		for( size_t idau=0; idau<2; idau++ )
		  {
		    const Candidate* dau = ZCand->daughter(idau);
		    CandInfo* dauInfo = dau->info();
		    dauIdx[idau] = dauInfo->getInt("index");
		  }

		declareInt(    "run",    run       );
		declareInt(    "event",  event     );
		declareInt( "counter", ++counter   );
		declareInt(    "mode",  lmode      );
		declareInt(    "hyp",  ihyp        );
		declareInt(    "iZ",  iZ           );
		declareBool(  "isTight", isTight   );
		declareFloat(  "mll",    mll       );
		declareFloat(  "qTll",   qTll      );
		declareFloat(  "yll",    yll       );
		declareFloat(  "hll",    hll       );
		declareInt(  "dau0",    dauIdx[0]  );
		declareInt(  "dau1",    dauIdx[1]  );
	       	       
		// vertices
		declareInt(   "ZV_id", ZV_id   );
		declareFloat(  "ZV_X", ZV_X    );
		declareFloat(  "ZV_Y", ZV_Y    );
		declareFloat(  "ZV_Z", ZV_Z    );
	       
		tm.fill();	       
	      }
	  }
      }
  }

  {
    const CandAssoc& PatJetPfJetAssoc  = _e.candAssoc( "PatJet-PfJet"  ); 
    //   const CandAssoc& PatJetLeptonAssoc = _e.candAssoc( "PatJet-lepton" ); 
    tupleName.push_back("CleanJets");
    tm.setTree(tupleName.back(),dir_);
    counter = 0;
    for( size_t ij=0; ij<cleanJets.size(); ij++ ) 
      { 
	Candidate* jet = cleanJets[ij];
	float jEta_ = jet->eta();
	float jEt_  = jet->Et();
	float jPhi_ = jet->phi();
	int vtxId_ = jet->vertexIndex();

	CandInfo* jinfo = jet->info();
	int index_ = jinfo->getInt("index");

	Candidate* pfJet  = PatJetPfJetAssoc.getSecond( jet );
	assert( pfJet!=0 );
	CandInfo* pfjinfo = pfJet->info();
	float area = pfjinfo->getFloat( "area" );
	float chargedEmEnergyFraction = pfjinfo->getFloat( "chargedEmEnergyFraction" );
	float chargedHadronEnergyFraction = pfjinfo->getFloat( "chargedHadronEnergyFraction");
	float chargedMuEnergyFraction = pfjinfo->getFloat( "chargedMuEnergyFraction" );
	float chargedMultiplicity = pfjinfo->getFloat( "chargedMultiplicity" );
	float etaPhiMom = pfjinfo->getFloat( "etaPhiMom" );
	float maxD = pfjinfo->getFloat( "maxD" );
	float muonMultiplicity = pfjinfo->getFloat( "muonMultiplicity" );
	float neutralEmEnergyFraction = pfjinfo->getFloat( "neutralEmEnergyFraction" );
	float neutralHadronEnergyFraction = pfjinfo->getFloat( "neutralHadronEnergyFraction" );
	float neutralMultiplicity = pfjinfo->getFloat( "neutralMultiplicity" );
	float phiPhiMom = pfjinfo->getFloat( "phiPhiMom" );

	declareInt(      "run",   run    );
	declareInt(    "event", event  );
	declareInt(  "counter", ++counter );
	declareFloat(     "Et", jEt_   ); 
	declareFloat(    "eta", jEta_  );
	declareFloat(    "phi", jPhi_  );
	declareInt(    "index", index_ );
	declareInt(    "vtxId", vtxId_ );
	declareFloat( "area",area);
	declareFloat( "chargedEmEnergyFraction", chargedEmEnergyFraction);
	declareFloat( "chargedHadronEnergyFraction", chargedHadronEnergyFraction);
	declareFloat( "chargedMuEnergyFraction", chargedMuEnergyFraction);
	declareFloat( "chargedMultiplicity", chargedMultiplicity);
	declareFloat( "etaPhiMom", etaPhiMom);
	declareFloat( "maxD", maxD);
	declareFloat( "muonMultiplicity", muonMultiplicity);
	declareFloat( "neutralEmEnergyFraction", neutralEmEnergyFraction);
	declareFloat( "neutralHadronEnergyFraction", neutralHadronEnergyFraction);
	declareFloat( "neutralMultiplicity", neutralMultiplicity);
	declareFloat( "phiPhiMom", phiPhiMom);

	bool isCleanJet(false);
	bool CJisCharged(false);
	bool CJisIsolatedLepton(false);
	float CJsumPt(0);
	//
	jinfo->getBool( "isCleanJet", isCleanJet );
	jinfo->getBool( "CJisCharged", CJisCharged );
	jinfo->getBool( "CJisIsolatedLepton", CJisIsolatedLepton );
	jinfo->getFloat( "CJsumPt", CJsumPt );
	//
	declareBool( "isCleanJet", isCleanJet );
	declareBool( "CJisCharged", CJisCharged );
	declareBool( "CJisIsolatedLepton", CJisIsolatedLepton );
	declareFloat( "CJsumPt", CJsumPt );

	tm.fill();
      }      
  }

  // special ntuples for MC events
  if( isMC )
    {
      Candidate* genCand_ = _e.decayTree();
      assert( genCand_!=0 );

      // initiate with the list of EWK bosons
      CandList VList;
      // ---> the Z bosons
      CandUtil::get(  23, genCand_, VList );
      // ---> the W bosons
      CandUtil::get(  24, genCand_, VList );
      CandUtil::get( -24, genCand_, VList );

      tupleName.push_back("MCTruth");
      tm.setTree(tupleName.back(),dir_);
      counter = 0;      
      // recursive filling of the MC tuple
      fillMCList( VList, counter, 4 );
      
      // the MC event tuple
      tupleName.push_back("MCEventTuple");
      tm.setTree(tupleName.back(),dir_);
      declareInt(      "run",   run    );
      declareInt(    "event", event  );
      tm.declare<  string*  > ( "categ", &categ_ptr );
      tm.declare<  string*  > ( "categ_long", &categ_long_ptr );
      tm.fill();
    }

  // event ntuple: one entry per event
  {
    tm.setTree( "EventTuple",dir_);
    declareInt(        "run", run   );
    declareInt(      "event", event );
    for( size_t ii=0; ii<tupleName.size(); ii++ )
      {
	declareInt( tupleName[ii], _int[tupleName[ii]+"_counter"] );
	_int[tupleName[ii]+"_counter"] = 0;
      }
    declareFloat(  "MET", MET );
    declareFloat(  "METPhi", METPhi );
    declareFloat(  "sumEt", sumEt );
    declareFloat(  "mEtSig", mEtSig );     
    declareFloat( "rhoFJ" , rhoFJ );
    declareInt(      "PV_id", PV_id );
    declareFloat(     "PV_X", PV_X  );
    declareFloat(     "PV_Y", PV_Y  );
    declareFloat(     "PV_Z", PV_Z  );
    declareInt(   "nZHypo_0", nZHypos[0] );
    declareInt(   "nZHypo_1", nZHypos[1] );
    declareInt(    "nVertex", nVertex );
    declareInt(    "nGoodVertices", nGoodVertices );

    // trigger ntuple: one per event
    vector< int >* ilines_ = new vector< int >;    
    {
      map< int, string > lines_;
      _e.firedLines( lines_ );
      for( map< int, string >::iterator itHlt = lines_.begin();
	   itHlt != lines_.end(); ++itHlt )
	{
	  string lineName_ =  itHlt->second;
	  if( _verbose )
	    cout << "HLT[" << itHlt->first << "] --> " << lineName_ << endl;
	  int iline_(0);
	  map< string, int >::iterator it_ = _hltLines.find( lineName_ );
	  if( it_ == _hltLines.end() ) _hltLines[ lineName_ ] =  ++_iline; 
	  iline_ = _hltLines[ lineName_ ];
	  assert( iline_>0 );
	  ilines_->push_back(iline_);
	}
      tm.declare< vector<int>* > ( "lines", &ilines_ );
    }

    // fill the event tuple fpr this event
    tm.fill();

    delete ilines_;
  }

  return true;

}

void
ZllAnalysis::bookHistograms()
{
  _iline = 0;
  _hltLines.clear();
}

void
ZllAnalysis::writeHistograms()
{
  // the HLT line ntuple: one entry per HLT line
  string dir_(_name);
  int index_;
  string lineName_;
  string* lineName_ptr = &lineName_;
  for( map< string, int >::iterator it = _hltLines.begin();
       it != _hltLines.end(); ++it )
    {
      lineName_ = it->first;
      index_ = it->second;
      cout << "HLT[" << index_ << "] ---> " << lineName_  << endl;
      tm.setTree( "HltLines",dir_);      
      tm.add<    int   > ( "lineIndex", &index_ );
      tm.add<  string*  > ( "lineName", &lineName_ptr );
      tm.flush();
    }
}

size_t
ZllAnalysis::buildZllLists( int imode )
{
  // get the lepton list
  //  CandList& leptonList = _e.userCandList( leptonListName(imode) );

  const CandList& initialList = _e.compositeCandList( ZmodeName(imode) );

  // the list of selected Z candidates
  CandList list;
  for( size_t iZ=0; iZ<initialList.size(); iZ++ )
    {
      Candidate* cand = initialList[iZ];

      // selectZCand: see core/SampleAnalysis.cc
      //     full DY window = [20,oo[
      int nVT=0; // min number of VT leptons
      int nT =0; // min number of  T leptons
      int nL =0; // min number of  L leptons
      int nVL=2; // min number of VL leptons
      if( selectZCand( *cand, 
		       Z0_fullDY[0], Z0_fullDY[1], 
		       nVT, nT, nL, nVL )          )
	{
	  list.push_back( cand );
	}
    }

  // sort and dispatch the list
  sort( list.begin(), list.end(), ZllCompare() );
  vector< CandList > Zhypos;
  CandUtil::dispatchList( list, Zhypos );
  size_t nZhypos_ = Zhypos.size();

  // create the hypotheses
  for( size_t ihyp=0; ihyp<nZhypos_; ihyp++ )
    {
      CandList& ZList_ = _e.userCandList( ZListName(imode,ihyp) );
      ZList_ = Zhypos[ihyp];
    }
  
  // return the number of hypotheses
  return nZhypos_;
}

string 
ZllAnalysis::ZmodeName(          int imode )
{
  assert( imode>=00 && imode<kNZmode );
  switch ( imode )
    {
    case kZeeMode:
      return "Zee";
    case kZmmMode:
      return "Zmm";
    }    
  return "";
}

string 
ZllAnalysis::ZListName(          int imode, size_t ihyp )
{
  string ZmodeName_( ZmodeName(imode) );
  ZmodeName_ += "_";
  stringstream ss;
  ss << ihyp;
  ZmodeName_ += ss.str();
  return ZmodeName_;
}

string 
ZllAnalysis::leptonName(         int ilepton )
{
  assert( ilepton==0 || ilepton==1 );
  switch ( ilepton )
    {
    case 0: 
      return "Electron";
    case 1:
      return "Muon";
    }
  return 0;
}

string 
ZllAnalysis::leptonListName(     int ilepton )
{
  return leptonName(ilepton)+"VeryLoose";
}

string 
ZllAnalysis::leptonListName(     int ilepton, size_t ihyp )
{
  string leptonListName_( leptonListName( ilepton ) );
  leptonListName_ += "_";
  stringstream ss;//create a stringstream
  ss << ihyp;//add number to the stream
  leptonListName_ += ss.str();

  return leptonListName_;
}

void
ZllAnalysis::declareInt( string var, int value )
{
  string tname  = tm.tree()->GetName();
  string varKey = tm.varKey( tname, var ); 
  if( _int.count(varKey)==0 )  
    tm.declare<int>(    var,    &_int[varKey] );
  _int[varKey] = value;
}

void
ZllAnalysis::declareFloat( string var, float value )
{
  string tname  = tm.tree()->GetName();
  string varKey = tm.varKey( tname, var ); 
  if( _float.count(varKey)==0 )  
    tm.declare<float>(    var,    &_float[varKey] );
  _float[varKey] = value;
}

void
ZllAnalysis::declareBool( string var, bool value )
{
  string tname  = tm.tree()->GetName();
  string varKey = tm.varKey( tname, var ); 
  if( _bool.count(varKey)==0 )  
    tm.declare<bool>(    var,    &_bool[varKey] );
  _bool[varKey] = value;
}

// recursive filling of the MC list...
bool
ZllAnalysis::fillMCList( const CandList& theList, int& counter, int nLevel )
{
  nLevel--;
  if( nLevel==0 ) return false;

  int run   = _e.run();
  int event = _e.event();

  int motherCounter=counter;
  for( size_t iV=0; iV<theList.size(); iV++ ) 
    { 
      Candidate* V_ = theList[iV];
      int pdgId_   = V_->pdgCode();
      float pt_   = V_->pt();
      float eta_  = V_->eta();
      float phi_  = V_->phi();
      float mass_ = V_->mass();
      float E_    = V_->E();
      float pz_   = V_->pz();
      float y_ = KineUtils::y( E_, pz_ );
      
      float xv_(0);
      float yv_(0);
      float zv_(0);
      //      int vtxId_ = -1;
      Vertex* vtx_ = V_->vertex();
      if( vtx_!=0 )
	{
	  //	  vtxId_ = vtx_->index();
	  xv_ = vtx_->pos().X();
	  yv_ = vtx_->pos().Y();
	  zv_ = vtx_->pos().Z();
	}
      
      declareInt(      "run",   run    );
      declareInt(    "event", event  );
      declareInt(  "counter", ++counter );
      declareInt(  "mother", motherCounter );
      declareInt(  "pdgId", pdgId_ );
      declareFloat(  "pt", pt_ );
      declareFloat(  "eta", eta_ );
      declareFloat(  "phi", phi_ );
      declareFloat(  "mass", mass_ );
      declareFloat(  "E", E_ );
      declareFloat(  "pz", pz_ );
      declareFloat(  "y", y_ );
      //      declareInt(    "vtxId", vtxId_ );
      declareFloat(  "VX", xv_ );
      declareFloat(  "VY", yv_ );
      declareFloat(  "VZ", zv_ );
      tm.fill();

      const CandList& dauList = V_->theDaughters();
      fillMCList( dauList, counter, nLevel );

    }
  return true;
}
  // void
// ZllAnalysis::declareFloat( string var )
// {
//   string tname  = tm.tree()->GetName();
//   string varKey = tm.varKey( tname, var ); 
//   tm.declare<float>(    var,    &_float[varKey] );
// }

// void
// ZllAnalysis::declareBool( string var )
// {
//   string tname  = tm.tree()->GetName();
//   string varKey = tm.varKey( tname, var ); 
//   tm.declare<bool>(    var,    &_bool[varKey] );
// }


