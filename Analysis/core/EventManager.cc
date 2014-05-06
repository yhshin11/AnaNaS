#include "Analysis/core/EventManager.hh"

#include <algorithm>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <sstream>
using namespace std;

#include <TH1F.h>

#include "MECore/src/MEEBGeom.hh"
#include "MECore/src/MEEEGeom.hh"
#include "Analysis/utils/Constants.hh"
#include "Analysis/utils/ParticleData.hh"
#include "Analysis/utils/KineUtils.hh"
#include "Analysis/utils/FileUtils.hh"
#include "Analysis/tools/CandUtil.hh"
#include "Analysis/selectors/ConeSelector.hh"
#include "Analysis/core/Sample.hh"

// GHM!!!
#include "Analysis/tools/PrintTree.hh"

EventManager* EventManager::_instance=0;

string EventManager::objectName[EventManager::kNType] = {
  "Track", 
  "GsfTrack", 
  "Photon",
  "Electron", 
  "Muon",
  "Tau",
  "CaloTower", 
  "EcalSC", 
  "EcalBC",
  "EcalPSC",
  "EcalPfSC", 
  "EcalPfBC",
  "PfCand",
  "PfCand_lowPt",
  "McTruthCandidate", 
  "GenEventCandidate"
};

string EventManager::jetName[EventManager::kNJet] = {
  "PatJet", 
  "CaloJet", 
  "TrkJet",
  "PfJet",
  "GenJet"
};

string EventManager::metName[EventManager::kNMet] = {
  "PatMet", 
  "CaloMet",
  "TrkMet"
  "PfMet",
  "GenMet"
};

ClassImp(EventManager)

EventManager::EventManager( const char* filename, const char * collectionFileName ) : 
    _collectionFileName( collectionFileName),
    _a(0)
{ 
  checkFileName( filename );

  // init
  init();
  
  // refresh (just to make sure)
  refresh();
}

EventManager::EventManager( vector< string > filename, const char * collectionFileName ) :  
    _filename(filename), 
    _collectionFileName(collectionFileName),
    _a(0)
{ 

  //MM fixing double couting for datasets with only one file
  // if( filename.size()==1 ) 
  //   {
  //     checkFileName( filename[0].c_str() );
  //   }
 
  // init
  init();

  // refresh (just to make sure)
  refresh();
}


void
EventManager::checkFileName( const char* filename )
{

  TString filename_(filename);
  if( filename_.Contains(".root") )
    {
      // this is a single file
      //MM not all the time! sometimes you run a
      //dataset with only one file, and in this case you
      // currently run twice on the same time
      // fixed in previous function

      cout << "Creating Event Manager for single file " << filename_ << endl;
      _filename.push_back( filename ); // this line was commented??? (2/3/2012)

    }
  else if( filename_.Contains(".") )
    {
      // this is probably a list of events      
      EventList  events;
      FileUtils::getListOfEvents( filename, events );	
      _a = new EventServer( events, _collectionFileName.c_str());
    }
  else
    {
      // this is probably a directory containing root files
      TString files_ = filename_ + "/*.root";                               
      _filename.clear();
      _filename = FileUtils::getFileList( files_.Data() );
      size_t nfiles_ = _filename.size();
      //      assert( nfiles_>0 );
      if(nfiles_ <= 0){
        std::cerr << "No file " << files_ << " found!" << std::endl;
        abort();
      }
      cout << "Creating Event Manager for " << nfiles_ << " files in " << filename_ << endl;
    }

}

EventServer&
EventManager::a()
{ 
  if( _a==0 ) 
    {
      _a = new EventServer( _filename, _collectionFileName.c_str() );
    }
  return *_a; 
}

void
EventManager::init()
{
  if( _instance!=0 )
    {
      cout << "Only one Event Manager can live at a given time" << endl;
      assert(0);
    }
  _instance = this;

  // dimension the vector of cand lists
  _candList.resize( kNType );
  //_jetList.resize(  kNJet );
  _jetList.clear();
  //  _met.resize(      kNMet );
  _met.clear();
  
}

EventManager::~EventManager( ) 
{    
  refresh();
  _instance = 0;
 
}

EventManager*
EventManager::e()
{
  assert( _instance!=0 );
  return _instance;
}

bool
EventManager::nextEvent()
{
  refresh();
  
  bool ok_ = a().nextEvent();
  if( !ok_ ) return false;
  
  loadEvent();  

  return true;
}

bool
EventManager::goToEvent( int ievt )
{
  refresh();

  bool ok_ = a().event( ievt );
  if( !ok_ ) return false;

  loadEvent();  

  return true;
}


bool
EventManager::findEvent( int ievt, int irun )
{
  refresh();

  for(int iev=1;iev<a().getNevtInFile()+1;iev++) {
    
    bool ok_ = a().event( iev );
    if( !ok_ ) return false;
    
    if(event() == ievt && run() == irun ) {
      loadEvent();   break;
    }
  } 
  return true;
}


void
EventManager::refresh()
{
  // clear candidate lists
  for( size_t ii=0; ii<kNType; ii++ )
    { _candList[ii].clear(); }
 //  for( size_t ii=0; ii<kNJet; ii++ )
//     { _jetList[ii].clear(); }
  _jetList.clear();
  //  for( size_t ii=0; ii<kNMet; ii++ )
  //    { _met[ii]=0; }
  _met.clear();
  _ecalRecHits.clear();
  _vertices.clear();
  _compositeCandList.clear();
  _userCandList.clear();

  _uncPfCands.clear();

  _categ        = "";
  _categ_long   = "";
  _inAcceptance = true;
  _neutrinosFromBosonsP4.SetPxPyPzE(0,0,0,0);
  _neutrinosFromTausP4.SetPxPyPzE(0,0,0,0);

  // reset the candidate factory
  Candidate::reset();

  // delete the previous event objects
  ObjectStore::clear();

 
  _primaryVertex = 0;
  _genCand       = 0;

  _caloTowerByIndex.clear();
  _caloJetTowerMap.clear();

  _cleanLeptons.clear();
  _cleanJets.clear();
  _leptonJets.clear();

  for( map<string,CandIdIsoMap>::iterator it=_isoMap.begin(); it!=_isoMap.end(); ++it )
    {
      it->second.clear();
    }
  _isoMap.clear();

//   for( map<string,CandIdEcalMap>::iterator it=_ecalMap.begin(); it!=_ecalMap.end(); ++it )
//     {
//       it->second.clear();
//     }
//   _ecalMap.clear();

  for( map<string,CandIdRecHitMap>::iterator it=_ecalRecHitMap.begin(); it!=_ecalRecHitMap.end(); ++it )
    {
      it->second.clear();
    }
  _ecalRecHitMap.clear();
  
  for( map<string,CandIdPSRecHitMap>::iterator it=_ecalPSRecHitMap.begin(); it!=_ecalPSRecHitMap.end(); ++it )
    {
      it->second.clear();
    }
  _ecalPSRecHitMap.clear();

  for( map<string,CandMap>::iterator it=_candMap.begin(); it!=_candMap.end(); ++it )
    {
      it->second.clear();
    }
  _candMap.clear();
 
  for( map<string,CandAssoc>::iterator it=_candAssoc.begin(); it!=_candAssoc.end(); ++it )
    {
      it->second.clear();
    }
  _candAssoc.clear();
 
  _mcMatchingMap.clear();

  for( map<string,TH1*>::iterator it=_sumOfEt.begin(); 
       it!=_sumOfEt.end(); ++it )
    {
      TH1* h_ = it->second;
      h_->Reset();
    }
}

void
EventManager::loadEvent()
{
  // loop over tracks
  size_t nVertices = a().n("vtx");
  
  for( size_t ii=0; ii<nVertices; ii++ )
    {
      assert( a().load( "vtx", ii ) );
      float vtx_x = a().get_f("vtx","x");
      float vtx_y = a().get_f("vtx","y");
      float vtx_z = a().get_f("vtx","z");
      Vertex* vtx_ = Vertex::create();
      vtx_->setXYZ( vtx_x, vtx_y, vtx_z );
      // special track information
      CandInfo* info_ = CandInfo::create();
      info_->setInt( "index", ii );
      setIntInfo(   info_, "vtx", "nTracks"  );
      setFloatInfo( info_, "vtx", "ndof"     );
      setBoolInfo(  info_, "vtx", "isFake"   );
      setBoolInfo(  info_, "vtx", "isGoodVtx" );
      vtx_->setInfo( info_  );
      _vertices.push_back( vtx_ );
    }

  _iPrimaryVertex = -1;
  _nGoodVertices  =  0;
  if( nVertices==0 )
    _primaryVertex = Vertex::create();
  else
    {
      // by default, take the first one...
      _primaryVertex = _vertices[0];

      // good vertex criteria
      float ndof_min = 4.;
      float z_max    = 24.;
      float perp_max = 2.;

      bool findPV(false);
      for( size_t iv=0; iv<_vertices.size(); iv++ ) 
	{
	  Vertex* Vert = _vertices[iv];
	  CandInfo* infoV = Vert->info();
   
	  float ndof(5.);
	  infoV->getFloat("ndof",ndof);
	  bool isFake(false);
	  infoV->getBool("isFake",isFake);
	  bool goodV=false;
	  if( 
	     fabs((Vert->pos()).Z())<z_max && 
	     (Vert->pos()).Perp()<perp_max && 
	     ndof>ndof_min                 && 
	     !isFake 
	      )
	    goodV=true;

	  infoV->setBool("isGoodVertex",goodV);
	  
	  if( goodV )
	    {
	      if( !findPV )
		{
		  _primaryVertex  = Vert;
		  _iPrimaryVertex = (int)iv;
		  findPV = true;
		} 
	      _nGoodVertices++;
	    }
	}      
    }
  setPrimaryVertex();
}

void
EventManager::setPrimaryVertex( Vertex* vtx )
{
  if( vtx==_primaryVertex ) return;
  if( vtx!=0 ) _primaryVertex = vtx;
  if( _primaryVertex==0 ) return;

  for( size_t iv=0; iv<_vertices.size(); iv++ ) 
    {
      Vertex* vtx_ = _vertices[iv];
      if( vtx_ == _primaryVertex ) vtx_->setAsPrimary();
      else vtx_->unsetAsPrimary();      
    }
  //setListsOfCleanJets(); //FIXME MM :  need to check carefully double counting in MET algorithms -> clean at vertex lvl not useful for MET PU mitigation 
}

void
EventManager::setListsOfCleanJets()
{
  loadJetPfCandMaps();
  const CandList& pfList = candList( kPfCand );
  for( size_t ii=0; ii<pfList.size(); ii++ )
    {
      pfList[ii]->info()->setBool("unClus",true);
    }

  //Disabled, not used MM
  // CandList& cleanLeptons        = userCandList("cleanLeptons");
  CandList& cleanJets           = userCandList("cleanJets");
  CandList& leptonJets          = userCandList("leptonJets");
  CandList& cleanJets_caloOnly  = userCandList("cleanJets_caloOnly");
  CandList& cleanJets_charged   = userCandList("cleanJets_charged");

  //  CandList& cleanJets_central   = userCandList("cleanJets_central");
  //  CandList& cleanJets_forward   = userCandList("cleanJets_forward");

  cleanJets.clear();
  leptonJets.clear();
  //  cleanJets_central.clear();
  //  cleanJets_forward.clear();
  cleanJets_caloOnly.clear(); 
  cleanJets_charged.clear(); 

  // the neutral jets
  if( _cleanJets.count( -1 ) )
    {      
      cleanJets_caloOnly = _cleanJets[-1];
    }

  //  the charged jets at the primary vertex
  int primaryVtxId = _primaryVertex->index();
  if( _cleanJets.count( primaryVtxId ) )
    {
      cleanJets_charged = _cleanJets[primaryVtxId];
    }

  CandUtil::mergeLists( cleanJets_caloOnly, cleanJets_charged, cleanJets );

  // the lepton jets
  if( _leptonJets.count( primaryVtxId ) )
    {
      leptonJets = _leptonJets[primaryVtxId];
    }
  
//   const CandAssoc& PatJetPfJetAssoc  = candAssoc( "PatJet-PfJet"  ); 
//   const CandAssoc& PatJetLeptonAssoc = candAssoc( "PatJet-lepton" ); 
//   const CandAssoc& electronPfAssoc   = candAssoc( "electron-PfCand" ); 
//   const CandAssoc& muonPfAssoc       = candAssoc( "muon-PfCand" ); 
   
//   float pt_min  = 10;
//   float eta_max =  5;
//   const CandList& patJetList = jetList( kPatJet,-1 );
//   for( size_t ii=0; ii<patJetList.size(); ii++ )
//     {
//       Candidate* jet    = patJetList[ii];
      
//       if( jet->pt() < pt_min ) continue;
//       bool isCentral = fabs( jet->eta() )<eta_max;
//       bool isForward = !isCentral;
      
//       Candidate* pfJet  = PatJetPfJetAssoc.getSecond( jet );
//       float dR_pfJet=-1;
//       if( pfJet!=0 )
// 	{
// 	  dR_pfJet =  CandUtil::dR( jet, pfJet );
// 	}
//       else
// 	{
// 	  cout << "Warning, no associated pfJet" << endl;
// 	}
//       Candidate* lepton = PatJetLeptonAssoc.getSecond( jet );
//       Candidate* pfLepton(0);
      
//       // list the constituents
//       const CandMap& jetConst_ = candMap("Jet-PfCand");	  
//       CandMapIterator it_lower_ = jetConst_.lower_bound( jet );
//       CandMapIterator it_upper_ = jetConst_.upper_bound( jet );
//       CandMapIterator it;
      
//       map< int, CandList > chargedAtVertex;
//       CandList neutralConstituents;
//       for( it=it_lower_; it!=it_upper_; it++ )
// 	{ 
// 	  Candidate* cand = it->second;
// 	  bool isNeutral = cand->charge()==0;
// 	  if( isNeutral )
// 	    {
// 	      neutralConstituents.push_back(cand);
// 	    }
// 	  else
// 	    {
// 	      int vtxId = cand->vertexIndex();
// 	      chargedAtVertex[vtxId].push_back( cand );
// 	    }
// 	}
      
//       CandList chargedComponents;
//       Candidate* neutralComponent(0);
//       for( map< int, CandList >::iterator itc=chargedAtVertex.begin();
// 	   itc!=chargedAtVertex.end(); ++itc )
// 	{
// 	  Candidate* sum_ = Candidate::create( itc->second );
// 	  chargedComponents.push_back( sum_ );
// 	}
//       sort( chargedComponents.begin(), chargedComponents.end(), Candidate::SortByPt() );
//       if( neutralConstituents.size()>0 )
// 	{
// 	  neutralComponent = Candidate::create( neutralConstituents );
// 	}
      
//       bool isNeutral = chargedComponents.size()==0;
//       bool isCharged = !isNeutral;
//       bool isAtPrimary = true;
//       int jetVtxId = -1;
//       int primaryVtxId = _primaryVertex->index();
//       if( isCharged ) 
// 	{
// 	  Vertex* jetVtx = chargedComponents[0]->vertex();
// 	  assert( jetVtx!=0 );
// 	  jetVtxId = jetVtx->index();
// 	  if( jetVtxId!=primaryVtxId ) isAtPrimary=false;
// 	  if( jet->vertexIndex()!=jetVtxId )
// 	    {
// 	      jet->setVertex( jetVtx );
// 	    }
// 	} 

//       if( !isAtPrimary ) continue;
      
//       // find the closest lepton candidate
//       float dR_lJet=-1;
//       bool isIsolatedLeptonJet = false;
//       int leptonVtxId = -1;
//       if( lepton!=0 )
// 	{
// 	  dR_lJet =  CandUtil::dR( jet, lepton );
// 	  pfLepton = electronPfAssoc.getSecond( lepton );
// 	  if( pfLepton==0 ) pfLepton = muonPfAssoc.getSecond( lepton );	      
// 	  leptonVtxId = lepton->vertexIndex();
	  
// 	  if( dR_lJet>0 && leptonVtxId==primaryVtxId )
// 	    {
// 	      // this is a lepton at the primary vertex
// 	      if( isCharged && dR_lJet<0.1 && lepton->pt()>0.80*chargedComponents[0]->pt() )
// 		{
// 		  // the lepton take more than 80% of the jet energy
// 		  isIsolatedLeptonJet = true;
// 		  cleanLeptons.push_back( lepton );
// 		}
// 	    }
// 	}
      
//       if( isIsolatedLeptonJet ) continue;

//       cleanJets.push_back(jet);
//       if( isForward ) 
// 	cleanJets_forward.push_back(jet);
//       else
// 	cleanJets_central.push_back(jet);

//       if( isNeutral )
// 	cleanJets_caloOnly.push_back(jet);
//       else
// 	cleanJets_charged.push_back(jet);
	
      
//       if(chargedComponents.size()>0 )
// 	{
// 	  for( size_t jj=0; jj<chargedComponents[0]->nDaughters(); jj++ )
// 	    {
// 	      chargedComponents[0]->daughter(jj)->info()->setBool("unClus",false);
// 	    }
// 	}
//       if( neutralComponent>0 )
// 	{
// 	  for( size_t jj=0; jj<neutralComponent->nDaughters(); jj++ )
// 	    {
// 	      neutralComponent->daughter(jj)->info()->setBool("unClus",false);
// 	    }
// 	}
//     }
}

bool
EventManager::setPrimaryVertex( int vtxid )
{
  bool ok=false;
  for( size_t iv=0; iv<_vertices.size(); iv++ ) 
    {
      Vertex* vtx_ = _vertices[iv];
      if( ok==false && vtx_->index() == vtxid )
	{
	  ok = true;
	  _primaryVertex = vtx_;
	  setPrimaryVertex();
	}
    }
  return ok;
}

void 
EventManager::loadTracks()
{
  
  //***** 
  static int idbg_=0;
  if( idbg_++<1 )
    {
      cout << "Warning -- you are loading the tracks !" << endl;
    }

  CandList& _tracks = _candList[ kTrack ];
  if( _tracks.size()!=0 ) return; // already loaded

  // loop over tracks
  size_t nTracks = a().n("trk");
  for( size_t ii=0; ii<nTracks; ii++ )
    {
      assert( a().load( "trk", ii ) );
      float pt_  = a().get_f("trk","pt");
      float eta_ = a().get_f("trk","eta");
      float phi_ = a().get_f("trk","phi");
      int pdgId_ = 211; // by default, charged pion
      if( pt_<0 ) pdgId_*=-1; 

      TVector3 p3_(0,0,0);
      p3_.SetPtEtaPhi( fabs( pt_), eta_, phi_ );

      float vx_  = a().get_f("trk","vx");
      float vy_  = a().get_f("trk","vy");
      float vz_  = a().get_f("trk","vz");
      Vertex* v3_ = Vertex::create( TVector3( vx_, vy_, vz_ ) );
      Vertex* vOff_ =  vertexMatching(v3_);

      Candidate* cand_ = Candidate::create( p3_, pdgId_, vOff_ ); 

      // special track information
      CandInfo* info_ = CandInfo::create();
      setIntInfo(   info_, "trk", "index"  );
      //      setIntInfo(   info_, "trk", "nTHits" );
      setFloatInfo( info_, "trk", "d0"     );
      setFloatInfo( info_, "trk", "z0"     );
      setFloatInfo( info_, "trk", "Chi2" );
      setFloatInfo( info_, "trk", "vx"     );
      setFloatInfo( info_, "trk", "vy"     );
      setFloatInfo( info_, "trk", "vz" );
      setIntInfo( info_, "trk", "numberOfValidHits"); 
      setIntInfo(info_, "trk", "numberOfLostHits");
      setIntInfo( info_, "trk", "Ndof" );
      setBoolInfo( info_, "trk", "isFPHitTrk" );
      setIntInfo( info_, "trk", "expInnHits" );
      
      setFloatInfo( info_, "trk", "impEta" );
      setFloatInfo( info_, "trk", "impPhi" );
      

      cand_->setInfo( info_  );

      // OK, store track candidate
      cand_ -> lock();
      _tracks.push_back( cand_ );      
    }
}

void
EventManager::loadTrackHits()
{

  if( _isoMap.count( "track-hits" )!=0 ) return;
  //  cout << "You are loading track hits! " << endl;

  //  cout << "Number of tracks " << _tracks.size() << endl;
  _isoMap.insert( make_pair( "track-hits", CandIdIsoMap() ) );

  loadTracks();
  CandList& _tracks = _candList[kTrack];
  for( size_t iTrack=0; iTrack<_tracks.size(); iTrack++ )
    {
      Candidate* cand_ = _tracks[iTrack];
      size_t uid_ = cand_->uid();
      const CandInfo* info_ = cand_->info();
      int index = info_->getInt( "index" );
      if( !a().load( "htrk", index ) ) continue;
      
      size_t nTHits_ = (size_t) a().get_i("htrk","nTHits");

      for( size_t jj=0; jj<nTHits_; jj++ )
	{
	  isoDep iso_;
	  iso_.dr  = a().get_vf("htrk","trho",jj);
	  iso_.val = 0;
	  iso_.eta = a().get_vf("htrk","teta",jj);
	  iso_.phi = a().get_vf("htrk","tphi",jj);
	  _isoMap["track-hits"].insert( make_pair( uid_, iso_ ) );
	}
    }
}

void 
EventManager::loadElectrons()
{
  CandList& _electrons = _candList[ kElectron ];
  if( _electrons.size()!=0 ) return; // already loaded
  
  // loop over electrons
  size_t nElectrons = a().n("el");
  for( size_t ii=0; ii<nElectrons; ii++ )
    {
      assert( a().load( "el", ii ) );
      float pt_  = a().get_f("el","pt");
      float eta_ = a().get_f("el","eta");
      float phi_ = a().get_f("el","phi");
      int pdgId_ = 11;
      if( pt_>0 ) pdgId_*=-1; 
      TVector3 p3_(0,0,0);
      p3_.SetPtEtaPhi( fabs( pt_), eta_, phi_ );

      float vx_  = a().get_f("el","vx");
      float vy_  = a().get_f("el","vy");
      float vz_  = a().get_f("el","vz");
      Vertex* v3_ = Vertex::create( TVector3( vx_, vy_, vz_ ) );
      Vertex* vOff_ =  vertexMatching(v3_);

      Candidate* cand_= Candidate::create( p3_, pdgId_, vOff_ ); 

      // special electron information
      CandInfo* info_ = CandInfo::create();

      setIntInfo( info_, "el", "index" );
      setIntInfo( info_, "el", "class" );

      //PF related variables
     //  setBoolInfo( info_, "el", "ecalDriven" );
//       setBoolInfo( info_, "el", "trkDriven" );

      setBoolInfo( info_, "el", "isConv" );

      setFloatInfo( info_, "el", "vx"          );
      setFloatInfo( info_, "el", "vy"          );
      setFloatInfo( info_, "el", "vz"          );
      setFloatInfo( info_, "el", "pin"         );
      setFloatInfo( info_, "el", "fBrem"       );
      setFloatInfo( info_, "el", "caloEnergy"  );
      setFloatInfo( info_, "el", "EOP"         );
      setFloatInfo( info_, "el", "ESeed"       );
      setFloatInfo( info_, "el", "ESeedOPout"  );
      setFloatInfo( info_, "el", "HOE"         );
      setFloatInfo( info_, "el", "caloRho"     );
      setFloatInfo( info_, "el", "caloEnergy"  );
      setFloatInfo( info_, "el", "caloEta"     );
      setFloatInfo( info_, "el", "caloPhi"     );
      setFloatInfo( info_, "el", "trackRho"    );
      setFloatInfo( info_, "el", "trackEta"    );
      setFloatInfo( info_, "el", "trackPhi"    );
      setFloatInfo( info_, "el", "dPhiIn"      );
      setFloatInfo( info_, "el", "dEtaIn"      );
      setFloatInfo( info_, "el", "sigiEtaiEta" );

      setIntInfo(   info_, "el", "nTrk"     );
      setFloatInfo( info_, "el", "sumOfPt"  );
      setFloatInfo( info_, "el", "sumOfPt2" );

      setIntInfo( info_, "el", "vetoID" );
      setIntInfo( info_, "el", "looseID" );
      setIntInfo( info_, "el", "mediumID" );
      setIntInfo( info_, "el", "tightID" );

      //Electron acceptance flags
      for( size_t ieAcc=0; ieAcc<a().eAcc.size(); ieAcc++ )
	{
	  bool ok_ = a().isEAccBitSet(ieAcc);
	  const string & eAccname_ = a().eAcc[ieAcc];
	  info_->setBool( eAccname_, ok_ );
	}

      //    Electron ID flags
      for( size_t ieid=0; ieid<a().eid.size(); ieid++ )
	{
	  bool ok_ = a().isEidBitSet(ieid);
	  const string & eidname_ = a().eid[ieid];
	  info_->setBool( eidname_, ok_ );
	}

      // electron Isolation
      setFloatInfo( info_, "el", "dr03TrkIso" );
      setFloatInfo( info_, "el", "dr03EcalIso" );
      setFloatInfo( info_, "el", "dr03HcalD1Iso" );
      setFloatInfo( info_, "el", "dr03HcalD2Iso" );
      
      setFloatInfo( info_, "el", "dr04EcalIso" );
      setFloatInfo( info_, "el", "dr04HcalD1Iso" );
      setFloatInfo( info_, "el", "dr04HcalD2Iso" );
      
      setFloatInfo( info_, "el", "chIso" );
      setFloatInfo( info_, "el", "nhIso" );
      setFloatInfo( info_, "el", "phIso" );
      setFloatInfo( info_, "el", "puIso" );

      // 1	-->	eidLoose
      // 2	-->	eidRobustHighEnergy
      // 3	-->	eidRobustLoose
      // 4	-->	eidRobustTight
      // 5	-->	eidTight
      //***** setBoolInfo(  info_, "el", "elID_Loose"  );
      //***** setBoolInfo(  info_, "el", "elID_RLoose" );
      //***** setBoolInfo(  info_, "el", "elID_RHE"    );
      //***** setBoolInfo(  info_, "el", "elID_Rtight" );
      //***** setBoolInfo(  info_, "el", "elID_Tight"  );
      cand_->setInfo( info_  );
      // OK store electron candidate
      cand_->lock();
      _electrons.push_back( cand_ );      

      //isolation is now only Pflow, isodeposit disabled
//       bool loadIsoMaps_(true);

//       if( loadIsoMaps_ )
// 	{
// 	  // now the maps
// 	  size_t uid_ = cand_->uid();
// 	  size_t nClus_ = (size_t) a().get_i("el","nclus");
// 	  _isoMap.insert( make_pair( "electron-EcalCluster", CandIdIsoMap() ) );
// 	  for( size_t jj=0; jj<nClus_; jj++ )
// 	    {
// 	      isoDep iso_;
// 	      iso_.dr  = a().get_vf("el","rhoclus",jj);
// 	      iso_.val = a().get_vf("el","eclus",jj);  // warning, energy 
// 	      iso_.eta = a().get_vf("el","etaclus",jj);
// 	      iso_.phi = a().get_vf("el","phiclus",jj);
// 	      _isoMap["electron-EcalCluster"].insert( make_pair( uid_, iso_ ) );
// 	    }
	  
// 	  size_t nEcalDep_ = (size_t) a().get_i("el","nECALdep");
// 	  _isoMap.insert( make_pair( "electron-EcalIsoDep", CandIdIsoMap() ) );
// 	  for( size_t jj=0; jj<nEcalDep_; jj++ )
// 	    {
// 	      isoDep iso_;
// 	      iso_.dr  = a().get_vf("el","dRECALdep",jj);
// 	      iso_.val = a().get_vf("el","E_ECALdep",jj);
// 	      iso_.eta = a().get_vf("el","etaECALdep",jj);
// 	      iso_.phi = a().get_vf("el","phiECALdep",jj);
// 	      _isoMap["electron-EcalIsoDep"].insert( make_pair( uid_, iso_ ) );
// 	    }
// 	}
    }
}


void 
EventManager::loadMuons()
{
  CandList& _muons = _candList[ kMuon ];
  if( _muons.size()!=0 ) return; // already loaded

  // loop over muons
  size_t nMuons = a().n("mu");
  for( size_t ii=0; ii<nMuons; ii++ )
    {
      assert( a().load( "mu", ii ) );
      float pt_  = a().get_f("mu","pt");
      float eta_ = a().get_f("mu","eta");
      float phi_ = a().get_f("mu","phi");
      int pdgId_ = 13;
      if( pt_>0 ) pdgId_*=-1; 
      TVector3 p3_(0,0,0);
      p3_.SetPtEtaPhi( fabs( pt_), eta_, phi_ );

      float vx_  = a().get_f("mu","vx");
      float vy_  = a().get_f("mu","vy");
      float vz_  = a().get_f("mu","vz");
      Vertex* v3_ = Vertex::create( TVector3( vx_, vy_, vz_ ) );
      // Candidate* cand_= Candidate::create( p3_, pdgId_, v3_ ); 
      
 
      // special muon information
      CandInfo* info_ = CandInfo::create();

      //D0 computation before vertex association //FIXME MM
      float d0 = KineUtils::d0(_primaryVertex->pos() , v3_->pos() ,p3_);
      info_->setFloat("d0", d0 );
      
      Vertex* vOff_ =  vertexMatching(v3_);
      Candidate* cand_= Candidate::create( p3_, pdgId_, vOff_ ); 

      setIntInfo( info_, "mu", "index" );
     
      setIntInfo( info_,   "mu", "nch"        );

      setFloatInfo( info_, "mu", "vx"          );
      setFloatInfo( info_, "mu", "vy"          );
      setFloatInfo( info_, "mu", "vz"          );
      setFloatInfo( info_, "mu", "em"         );
      setFloatInfo( info_, "mu", "had"        );
      setFloatInfo( info_, "mu", "normChi2"   );

      setIntInfo(   info_, "mu", "nTrk"       );
      setFloatInfo( info_, "mu", "sumOfPt"    );
      setFloatInfo( info_, "mu", "sumOfPt2"   );


      //*****setBoolInfo(  info_, "mu", "StandAlone"    );
      //*****setBoolInfo(  info_, "mu", "Global"        );
      //*****setBoolInfo(  info_, "mu", "GlobalPromptT" );
      //*****setBoolInfo(  info_, "mu", "TMLSOLPLoose"  );
      //*****setBoolInfo(  info_, "mu", "TMLSOLPTight"  );
      
      //Muon ID
      for( size_t imid=0; imid<a().muid.size(); imid++ )
	{
	  bool ok_ = a().isMidBitSet(imid);
	  string muidname_ = a().muid[imid];

	  //GHM : trick to remove the "muon:" in CMSSW52X
	  size_t found_ = muidname_.find(":");
	  if( found_!= string::npos )
	    muidname_ = muidname_.substr(found_+1);
	  //GHM
	  // cout << "!!!!GHM " << muidname_ << " ---> " << ok_ << endl; 

	  info_->setBool( muidname_, ok_ );
	}

      setBoolInfo( info_,  "mu", "isPFMuon" );

      //Muon Isolation
      setFloatInfo( info_,  "mu", "isoR03strk" );
      setFloatInfo( info_,  "mu", "isoR03sem" );
      setFloatInfo( info_,  "mu", "isoR03shad" );
      setFloatInfo( info_,  "mu", "isoR03strkVeto" );
      setFloatInfo( info_,  "mu", "isoR03semVeto" );
      setFloatInfo( info_,  "mu", "isoR03shadVeto" );
      setFloatInfo( info_,  "mu", "isoR03sum" );
      setFloatInfo( info_,  "mu", "isoR03pfch" );
      setFloatInfo( info_,  "mu", "isoR03pfcpart" );
      setFloatInfo( info_,  "mu", "isoR03pfnh" );
      setFloatInfo( info_,  "mu", "isoR03pfph" );
      setFloatInfo( info_,  "mu", "isoR03pfnhThres" );
      setFloatInfo( info_,  "mu", "isoR03pfphThres" );
      setFloatInfo( info_,  "mu", "isoR03pfpu" );
      setFloatInfo( info_,  "mu", "isoR04pfch" );
      setFloatInfo( info_,  "mu", "isoR04pfcpart" );
      setFloatInfo( info_,  "mu", "isoR04pfnh" );
      setFloatInfo( info_,  "mu", "isoR04pfph" );
      setFloatInfo( info_,  "mu", "isoR04pfnhThres" );
      setFloatInfo( info_,  "mu", "isoR04pfphThres" );
      setFloatInfo( info_,  "mu", "isoR04pfpu" );
      //End Muon Isolation ====================
      

      setFloatInfo( info_,  "mu", "dxyBS" );
      setFloatInfo( info_,  "mu", "dxyPV" );
      setFloatInfo( info_,  "mu", "dir"   );
      setIntInfo( info_,  "mu", "nValMuHits" );  
      setIntInfo( info_,  "mu", "nValTrackerHits" );  
      
      setIntInfo( info_,  "mu", "nTrkHits" );
      setIntInfo( info_,  "mu", "nPixHits" );
      setIntInfo( info_,  "mu", "nTrkLayerHits" );
      setIntInfo( info_,  "mu", "nPixLayerHits" );
      setIntInfo( info_,  "mu", "nMuonHits" );
      setIntInfo( info_,  "mu", "nMatch" );


      cand_->setInfo( info_  );

      // OK store muon candidate
      cand_->lock();
      _muons.push_back( cand_ );      

      //no more isodeposits, only pflow isoaltion is supported
    //   bool loadIsoMaps_(true);

//       if( loadIsoMaps_ )
// 	{
// 	  // now the maps
// 	  size_t uid_ = cand_->uid();
// 	  size_t nEcalDep_ = (size_t) a().get_i("mu","nECALdep");

// 	  _isoMap.insert( make_pair( "muon-EcalIsoDep", CandIdIsoMap() ) );
// 	  for( size_t jj=0; jj<nEcalDep_; jj++ )
// 	    {
// 	      isoDep iso_;
// 	      iso_.dr  = a().get_vf("mu","dRECALdep",jj);
// 	      iso_.val = a().get_vf("mu","valECALdep",jj);
// 	      iso_.eta = a().get_vf("mu","etaECALdep",jj);
// 	      iso_.phi = a().get_vf("mu","phiECALdep",jj);
// 	      _isoMap["muon-EcalIsoDep"].insert( make_pair( uid_, iso_ ) );
// 	    }
	  
// 	  // track hit map
// 	  _isoMap.insert( make_pair( "mu-trackHits", CandIdIsoMap() ) );
// 	  size_t nTHits_ = (size_t) a().get_i("mu","nTHits");
// 	  for( size_t jj=0; jj<nTHits_; jj++ )
// 	    {
// 	      isoDep iso_;
// 	      iso_.dr  = a().get_vf("mu","trho",jj);
// 	      iso_.val = 0;
// 	      iso_.eta = a().get_vf("mu","teta",jj);
// 	      iso_.phi = a().get_vf("mu","tphi",jj);
// 	      // cout << "mu-trackhit ";
// 	      // cout << iso_.dr << "/" << iso_.eta << "/" << iso_.phi << endl;

// 	      _isoMap["mu-trackHits"].insert( make_pair( uid_, iso_ ) );
// 	    }
// 	}
    }
}


void 
EventManager::loadTaus()
{
  CandList& _taus = _candList[ kTau ];
  if( _taus.size()!=0 ) return; // already loaded

  // loop over taus
  size_t nTaus = a().n("tau");
  for( size_t ii=0; ii<nTaus; ii++ )
    {
      assert( a().load( "tau", ii ) );
      float pt_  = a().get_f("tau","pt");
      float eta_ = a().get_f("tau","eta");
      float phi_ = a().get_f("tau","phi");
      int pdgId_ = 13;
      if( pt_>0 ) pdgId_*=-1; 
      TVector3 p3_(0,0,0);
      p3_.SetPtEtaPhi( fabs( pt_), eta_, phi_ );

      float vx_  = a().get_f("tau","vx");
      float vy_  = a().get_f("tau","vy");
      float vz_  = a().get_f("tau","vz");
      Vertex* v3_ = Vertex::create( TVector3( vx_, vy_, vz_ ) );
      // Candidate* cand_= Candidate::create( p3_, pdgId_, v3_ ); 
      
 
      // special tau information
      CandInfo* info_ = CandInfo::create();

      //D0 computation before vertex association //FIXME MM
      float d0 = KineUtils::d0(_primaryVertex->pos() , v3_->pos() ,p3_);
      info_->setFloat("d0", d0 );
      
      Vertex* vOff_ =  vertexMatching(v3_);
      Candidate* cand_= Candidate::create( p3_, pdgId_, vOff_ ); 

      setIntInfo( info_, "tau", "index" );
     
      setIntInfo( info_,   "tau", "nch"        );

      setFloatInfo( info_, "tau", "vx"          );
      setFloatInfo( info_, "tau", "vy"          );
      setFloatInfo( info_, "tau", "vz"          );
      setFloatInfo( info_,  "tau", "leadPPt"  );
      setFloatInfo( info_,  "tau", "leadCPPt"  );
      setFloatInfo( info_,  "tau", "leadNPPt"  );
      setFloatInfo( info_,  "tau", "trkd0"  );
      setFloatInfo( info_,  "tau", "trkdz"  );
      setIntInfo( info_,  "tau", "nCHSigC"  );
      setIntInfo( info_,  "tau", "nNHSigC"  );
      setIntInfo( info_,  "tau", "nPhSigC"  );
      setIntInfo( info_,  "tau", "nPSigC"  );
      setIntInfo( info_,  "tau", "nCHIsoC"  );
      setIntInfo( info_,  "tau", "nNHIsoC"  );
      setIntInfo( info_,  "tau", "nPhIsoC"  );
      setIntInfo( info_,  "tau", "nPIsoC"  );
      setFloatInfo( info_,  "tau", "sPtCHIsoC"  );
      setFloatInfo( info_,  "tau", "sPtPhIsoC"  );

      
      //Tau ID
      for( size_t itid=0; itid<a().tauid.size(); itid++ )
	{
	  bool ok_ = a().isTidBitSet(itid);
	  string tauidname_ = a().tauid[itid];

	  info_->setBool( tauidname_, ok_ );
	}

      setFloatInfo( info_,  "tau", "pfElMVA"  );
      setFloatInfo( info_,  "tau", "pfJetPt"  );
      setFloatInfo( info_,  "tau", "pfEtaPt"  );
      setFloatInfo( info_,  "tau", "pfPhiPt"  );

      cand_->setInfo( info_  );
      
      // OK store muon candidate
      cand_->lock();
      _taus.push_back( cand_ );      

    }
}




void 
EventManager::loadCaloTowers()
{
  //GHM
  //  abort();
  //***** 
  static int idbg_=0;
  if( idbg_++<1 )
    {
      cout << "Warning -- you are loading the caloTowers !" << endl;
    }

  CandList& _caloTowers = _candList[ kCaloTower ];
  if( _caloTowers.size()!=0 ) return; // already loaded

  TH1* et_h = sumOfEt( "ct_pt"    );
  TH1* em_h = sumOfEt( "ct_em"    );
  TH1* had_h = sumOfEt( "ct_had"   );
  TH1* outer_h = sumOfEt( "ct_outer" );
  
  TH1*  et_eta_h[4];
  TH1*  em_eta_h[4];
  TH1*  had_eta_h[4];
  TH1*  outer_eta_h[4];

  for(int jj=0;jj<4;jj++) {
    et_eta_h[jj] = sumOfEt( "ct_et_eta", jj );
    em_eta_h[jj] = sumOfEt( "ct_em_eta", jj );
    had_eta_h[jj] = sumOfEt( "ct_had_eta", jj );
    outer_eta_h[jj] = sumOfEt( "ct_outer_eta", jj );
  }
  

  // loop over caloTowers
  size_t nTowers = a().n("ct");
  for( size_t ii=0; ii<nTowers; ii++ )
    {
      assert( a().load( "ct", ii ) );
      float et_  = a().get_f("ct","pt"); 
      if( et_<=0 ) 
	{
	  //	  cout << "Warning -- calo tower with energy=0 " << ii << " et=" << et_ << " set to 0.1MeV" << endl;
	  et_ = 0.0001;
	  //	  continue;
	}
      float eta_ = a().get_f("ct","eta");
      float phi_ = a().get_f("ct","phi");
      int pdgId_ = 22; // by default, set as photon
      TVector3 p3_(0,0,0);
      p3_.SetPtEtaPhi( fabs( et_), eta_, phi_ );
      Candidate* cand_= Candidate::create( p3_, pdgId_ ); 
 
      // special caloTower information
      CandInfo* info_ = CandInfo::create();
      setIntInfo(   info_, "ct", "index" );
      setIntInfo(   info_, "ct", "ieta"  );
      setIntInfo(   info_, "ct", "iphi"  );
      setFloatInfo( info_, "ct", "energy"  );
      setFloatInfo( info_, "ct", "emEt"    );
      setFloatInfo( info_, "ct", "hadEt"   );
      setFloatInfo( info_, "ct", "outerEt" );
      cand_->setInfo( info_  );

      // OK, store caloTower candidate
      cand_ -> lock();
      _caloTowerByIndex[a().get_i("ct","index")] = cand_;
      _caloTowers.push_back( cand_ );            

      //fill the sum of Et histograms
      int iphi_    = info_->getInt("iphi");
      int ieta_    = info_->getInt("ieta");
      float em_    = info_->getFloat("emEt");
      float had_   = info_->getFloat("hadEt");
      float outer_ = info_->getFloat("outerEt");

      //cout << "!!!! ieta_ " << ieta_ << "/" << iphi_ << endl;

      //FIXME MM
      int jj(-1);
      if( ieta_> 0 ) 
	{
	  if(      iphi_>0  && iphi_<=36 ) jj=0;
	  else if( iphi_>36 && iphi_<=72 ) jj=1;	  	    
	}
      else if( ieta_<0 )
	{
	  if(      iphi_>0  && iphi_<=36 ) jj=2;
	  else if( iphi_>36 && iphi_<=72 ) jj=3;	  	    
	}
      assert( jj!=-1 );
      //      cout << "jj= " << jj  << endl;
      em_eta_h[jj]    -> AddBinContent( abs(ieta_), em_ );
      had_eta_h[jj]   -> AddBinContent( abs(ieta_), had_ );
      outer_eta_h[jj] -> AddBinContent( abs(ieta_), outer_ );
      et_eta_h[jj]    -> AddBinContent( abs(ieta_), et_ );
      
      if( abs(ieta_)<=20 )
	{
	  em_h    ->AddBinContent( iphi_, em_    ); 
	  had_h   ->AddBinContent( iphi_, had_   ); 
	  outer_h ->AddBinContent( iphi_, outer_ ); 
	  et_h    ->AddBinContent( iphi_, et_    );
	}
      else if( abs(ieta_)<=39 )
	{
	  int dum_ = 2;
	  while( dum_-- )
	    {
	      em_h    ->AddBinContent( iphi_, em_/2.    ); 
	      had_h   ->AddBinContent( iphi_, had_/2.   ); 
	      outer_h ->AddBinContent( iphi_, outer_/2. ); 
	      et_h    ->AddBinContent( iphi_, et_/2.    );
	      iphi_++;
	      if( iphi_>72 ) iphi_-=72;
	    }
	}
      else
	{
	  int dum_ = 4;
	  while( dum_-- )
	    {
	      em_h    ->AddBinContent( iphi_, em_/4.    ); 
	      had_h   ->AddBinContent( iphi_, had_/4.   ); 
	      outer_h ->AddBinContent( iphi_, outer_/4. ); 
	      et_h    ->AddBinContent( iphi_, et_/4.    );
	      iphi_++;
	      if( iphi_>72 ) iphi_-=72;
	    }
	}
      

    }
}

void 
EventManager::loadGsfTracks()
{
  CandList& _gsfTracks = _candList[ kGsfTrack ];
  if( _gsfTracks.size()!=0 ) return; // already loaded
  // loop over gsfTracks
  size_t nGsfTracks = a().n("gsf");
  for( size_t ii=0; ii<nGsfTracks; ii++ )
    {
      assert( a().load( "gsf", ii ) );
      //      int match_ = a().get_i( "gsf", "matchedTrack"  );
      //if( match_==-1 ) continue;  // !!!FIXME!!!
     
      float pt_  = a().get_f("gsf","pt");
      float eta_ = a().get_f("gsf","eta");
      float phi_ = a().get_f("gsf","phi");
      int pdgId_ = 11; // by default, electron
      if( pt_>0 ) pdgId_*=-1; 
      TVector3 p3_(0,0,0);
      p3_.SetPtEtaPhi( fabs( pt_), eta_, phi_ );

      float vx_  = a().get_f("gsf","vx");
      float vy_  = a().get_f("gsf","vy");
      float vz_  = a().get_f("gsf","vz");
      Vertex* v3_ = Vertex::create( TVector3( vx_, vy_, vz_ ) );
      Vertex* vOff_ =  vertexMatching(v3_);
      Candidate* cand_ = Candidate::create( p3_, pdgId_, vOff_ ); 

      // special gsfTrack information
      CandInfo* info_ = CandInfo::create();
      setIntInfo(   info_, "gsf", "matchedTrack"  );
      setIntInfo(   info_, "gsf", "index"       );
      setIntInfo(   info_, "gsf", "nLostHits"   );
      setIntInfo(   info_, "gsf", "nValidHits"  );
      setFloatInfo( info_, "gsf", "d0"     );
      setFloatInfo( info_, "gsf", "z0"     );
      setFloatInfo( info_, "gsf", "innerPoint_x"  );
      setFloatInfo( info_, "gsf", "innerPoint_y"  );
      setFloatInfo( info_, "gsf", "innerPoint_z"  );
      setFloatInfo( info_, "gsf", "outerPoint_x"  );
      setFloatInfo( info_, "gsf", "outerPoint_y"  );
      setFloatInfo( info_, "gsf", "outerPoint_z"  );
      setFloatInfo( info_, "gsf", "normChi2"      );
      setBoolInfo( info_, "gsf", "isFPHitTrk" );
      setIntInfo( info_, "gsf", "expInnHits" );
      
      cand_->setInfo( info_  );
      
      // OK, store gsfTrack candidate
      cand_ -> lock();
      _gsfTracks.push_back( cand_ );      

      bool loadIsoMaps_(true);
      if( loadIsoMaps_ )
	{
	  // now the maps
	  size_t uid_ = cand_->uid();
	  size_t nTHits_ = (size_t) a().get_i("gsf","nTHits");
	  _isoMap.insert( make_pair( "gsfTrack-hits", CandIdIsoMap() ) );

	  //	  cout << "uid=" << uid_ << endl;

	  for( size_t jj=0; jj<nTHits_; jj++ )
	    {
	      isoDep iso_;
	      iso_.dr  = a().get_vf("gsf","trho",jj);
	      iso_.val = 0.; 
	      iso_.eta = a().get_vf("gsf","teta",jj);
	      iso_.phi = a().get_vf("gsf","tphi",jj);
	      //	      cout << iso_.dr << "/" << iso_.eta << "/" << iso_.phi << ":";
	      _isoMap["gsfTrack-hits"].insert( make_pair( uid_, iso_ ) );
	    }	  
	  //	  cout << endl;
	}     
    }
}

void 
EventManager::loadPhotons()
{
  CandList& _photons = _candList[ kPhoton ];
  if( _photons.size()!=0 ) return; // already loaded

  // loop over photons
  size_t nPhotons = a().n("ph");
  for( size_t ii=0; ii<nPhotons; ii++ )
    {
      assert( a().load( "ph", ii ) );
      float et_  = a().get_f("ph","pt");
      float eta_ = a().get_f("ph","eta");
      float phi_ = a().get_f("ph","phi");
      int pdgId_ = 22; // photon

      TVector3 p3_(0,0,0);
      p3_.SetPtEtaPhi( fabs( et_), eta_, phi_ );

      Candidate* cand_ = Candidate::create( p3_, pdgId_, _primaryVertex ); 

      // special photon information
      CandInfo* info_ = CandInfo::create();
      setFloatInfo( info_, "ph", "caloRho" );
      setFloatInfo( info_, "ph", "caloEta" );
      setFloatInfo( info_, "ph", "caloPhi" );

      //FIXME MM pour éviter de se casser les pieds
      info_->setFloat( "caloEnergy", KineUtils::e( et_, eta_ ) );

      /* setFloatInfo( info_, "ph", "r9" );
      setFloatInfo( info_, "ph", "e5x5" );
      setFloatInfo( info_, "ph", "e3x3" );*/
      // setFloatInfo( info_, "ph", "r1x5" );
      setFloatInfo( info_, "ph", "sigmaEtaEta" );
      setFloatInfo( info_, "ph", "sigmaIetaIeta" );
      
      setFloatInfo( info_, "ph", "e1x5" );
      setFloatInfo( info_, "ph", "e2x5" );

      //Isolation
      setFloatInfo( info_, "ph", "ecalRecHitSumEtConeDR03" );
      setFloatInfo( info_, "ph", "ecalRecHitSumEtConeDR04" );
      setFloatInfo( info_, "ph", "hadronicDepth1OverEm" );
      setFloatInfo( info_, "ph", "hadronicDepth2OverEm" );
      setFloatInfo( info_, "ph", "hadronicOverEm" );
      setFloatInfo( info_, "ph", "hcalDepth1TowerSumEtConeDR03" );
      setFloatInfo( info_, "ph", "hcalDepth1TowerSumEtConeDR04" );
      setFloatInfo( info_, "ph", "hcalDepth2TowerSumEtConeDR03" );
      setFloatInfo( info_, "ph", "hcalDepth2TowerSumEtConeDR04" );
      setFloatInfo( info_, "ph", "hcalTowerSumEtConeDR03" );
      setFloatInfo( info_, "ph", "hcalTowerSumEtConeDR04" );
      setFloatInfo( info_, "ph", "nTrkHollowConeDR03" );
      setFloatInfo( info_, "ph", "nTrkHollowConeDR04" );
      setFloatInfo( info_, "ph", "nTrkSolidConeDR03" );
      setFloatInfo( info_, "ph", "nTrkSolidConeDR04" );
      setFloatInfo( info_, "ph", "trkSumPtHollowConeDR03" );
      setFloatInfo( info_, "ph", "trkSumPtHollowConeDR04" );
      setFloatInfo( info_, "ph", "trkSumPtSolidConeDR03" );
      setFloatInfo( info_, "ph", "trkSumPtSolidConeDR04" );
      setFloatInfo( info_, "ph", "ecalIso" );

      setBoolInfo( info_, "ph", "isEB" );
      setBoolInfo( info_, "ph", "isEE" );
      setBoolInfo( info_, "ph", "isEBGap" );
      setBoolInfo( info_, "ph", "isEEGap" );
      setBoolInfo( info_, "ph", "isEBEEGap" );
      
      setBoolInfo( info_, "ph", "hasPixelSeed" );
      setBoolInfo( info_, "ph", "isConverted" );
      
      //Photon IDs
      for( size_t iphid=0; iphid<a().phid.size(); iphid++ )
	{
	  bool ok_ = a().isPhidBitSet(iphid);
	  const string & phidname_ = a().phid[iphid];
	  info_->setBool( phidname_, ok_ );
	}
 
      cand_->setInfo( info_  );

      // OK, store photon candidate
      cand_ -> lock();
      _photons.push_back( cand_ );     

      //Super Cluster Map
      loadEcalSC();
      int scIndex_ = a().get_i("ph", "indexSC" );
      _candAssoc.insert ( make_pair( "photon-SuperCluster", CandAssoc() ) );
      if(scIndex_ != -1) {
	Candidate* sc_ = _candList[ kEcalSC ][ scIndex_ ];
	_candAssoc["photon-SuperCluster"].insert( cand_, sc_  );
      }
    }
}

void 
EventManager::loadEcalSC()
{
  CandList& _ecalSC = _candList[ kEcalSC ];
  if( _ecalSC.size()!=0 ) return; // already loaded

  size_t nEcalSC = a().n("sc"); 
  for ( size_t ii=0; ii<nEcalSC; ++ii )
    {
        assert( a().load( "sc", ii ) );
        float e_ = a().get_f("sc", "energy");
        float eta_ = a().get_f("sc", "eta");
        float theta_ = KineUtils::theta( eta_ );
        float pt_ = e_*sin(theta_);
        float phi_ = a().get_f("sc", "phi");
        int pdgId_ = 22; // photon by construction

        TVector3 p3_(0, 0, 0);
        p3_.SetPtEtaPhi( fabs( pt_), eta_, phi_ );

        Candidate* cand_ = Candidate::create( p3_, pdgId_, _primaryVertex ); 

        // special super cluster information
        CandInfo* info_ = CandInfo::create();
        setFloatInfo( info_, "sc", "rawEnergy" );

	setFloatInfo( info_, "sc", "preshowerEnergy" );
    	
	setFloatInfo( info_, "sc", "R4" );
	setFloatInfo( info_, "sc", "R9" );
	setFloatInfo( info_, "sc", "etaWidth" );
	setFloatInfo( info_, "sc", "phiWidth" );
	setIntInfo(   info_, "sc", "nHits" );
	setIntInfo(   info_, "sc", "index" );
	setIntInfo(   info_, "sc", "ixSeed" );
	setIntInfo(   info_, "sc", "iySeed" );
	setIntInfo(   info_, "sc", "izSeed" );
	setIntInfo(   info_, "sc", "ixClosestCell" );
	setIntInfo(   info_, "sc", "iyClosestCell" );
	setIntInfo(   info_, "sc", "izClosestCell" );

        cand_->setInfo( info_ );

	// OK, store Super Cluster
        cand_->lock();
        _ecalSC.push_back( cand_ );	
    }
}

void 
EventManager::loadEcalPfSC()
{
  CandList& _ecalPfSC = _candList[ kEcalPfSC ];
  if( _ecalPfSC.size()!=0 ) return; // already loaded

  size_t nEcalPfSC = a().n("pfSc");
  for ( size_t ii=0; ii<nEcalPfSC; ++ii )
    {
        assert( a().load( "pfSc", ii ) );
        float e_ = a().get_f("pfSc", "energy");
        float eta_ = a().get_f("pfSc", "eta");
        float theta_ = KineUtils::theta( eta_ );
        float pt_ = e_*sin(theta_);
        float phi_ = a().get_f("pfSc", "phi");
        int pdgId_ = 22; // photon by construction

        TVector3 p3_(0, 0, 0);
        p3_.SetPtEtaPhi( fabs( pt_), eta_, phi_ );

        Candidate* cand_ = Candidate::create( p3_, pdgId_, _primaryVertex ); 

        // special super cluster information
        CandInfo* info_ = CandInfo::create();
        setFloatInfo( info_, "pfSc", "rawEnergy" );

	setFloatInfo( info_, "pfSc", "preshowerEnergy" );
		
	setFloatInfo( info_, "pfSc", "R4" );
	setFloatInfo( info_, "pfSc", "R9" );
	setFloatInfo( info_, "pfSc", "etaWidth" );
	setFloatInfo( info_, "pfSc", "phiWidth" );
	setIntInfo(   info_, "pfSc", "nHits" );
	setIntInfo(   info_, "pfSc", "index" );

        cand_->setInfo( info_ );

	// OK, store Super Cluster
        cand_->lock();
        _ecalPfSC.push_back( cand_ );	
    }
}

void 
EventManager::loadEcalBC()
{
  //GHM
  //  cout << "You are loading the ECAL basic clusters" << endl;
  // abort();

  CandList& _ecalBC = _candList[ kEcalBC ];
  if( _ecalBC.size()!=0 ) return; // already loaded

  size_t nEcalBC = a().n("bc");
  for ( size_t ii=0; ii<nEcalBC; ++ii )
    {
        assert( a().load( "bc", ii ) );
        float e_ = a().get_f("bc", "energy");
        float eta_ = a().get_f("bc", "eta");
        float theta_ = KineUtils::theta( eta_ );
        float pt_ = e_*sin(theta_);
        float phi_ = a().get_f("bc", "phi");
        int pdgId_ = 22; // photon by construction

        TVector3 p3_(0, 0, 0);
        p3_.SetPtEtaPhi( fabs( pt_), eta_, phi_ );

        Candidate* cand_ = Candidate::create( p3_, pdgId_, _primaryVertex ); 

        // special super cluster information
        CandInfo* info_ = CandInfo::create();
	setFloatInfo( info_, "bc", "R4" );
	setFloatInfo( info_, "bc", "R9" );
	setIntInfo(   info_, "bc", "nHits" );
	setIntInfo(   info_, "bc", "index" );

	setFloatInfo( info_, "bc", "covEtaEta" );
	setFloatInfo( info_, "bc", "covEtaPhi" );
	setFloatInfo( info_, "bc", "covPhiPhi" );
	setFloatInfo( info_, "bc", "eBottom" );
	setFloatInfo( info_, "bc", "eLeft" );
	setFloatInfo( info_, "bc", "eTop" );
	setFloatInfo( info_, "bc", "eRight" );
	setFloatInfo( info_, "bc", "eMax" );

        cand_->setInfo( info_ );

	// OK, store Basic Cluster
        cand_->lock();
        _ecalBC.push_back( cand_ );
    }
}

void 
EventManager::loadEcalPSC()
{
  //GHM
  //  cout << "You are loading the ECAL basic clusters" << endl;
  // abort();

  CandList& _ecalPSC = _candList[ kEcalPSC ];
  if( _ecalPSC.size()!=0 ) return; // already loaded

  size_t nEcalPSC = a().n("psc");
  for ( size_t ii=0; ii<nEcalPSC; ++ii )
    {
        assert( a().load( "psc", ii ) );
        float e_ = a().get_f("psc", "energy");
        float eta_ = a().get_f("psc", "eta");
        float theta_ = KineUtils::theta( eta_ );
        float pt_ = e_*sin(theta_);
        float phi_ = a().get_f("psc", "phi");
        int pdgId_ = 22; // photon by construction

        TVector3 p3_(0, 0, 0);
        p3_.SetPtEtaPhi( fabs( pt_), eta_, phi_ );

        Candidate* cand_ = Candidate::create( p3_, pdgId_, _primaryVertex ); 

        // special super cluster information
        CandInfo* info_ = CandInfo::create();
	setIntInfo( info_, "psc", "plane" );
	setIntInfo( info_, "psc", "nHits" );
	setIntInfo( info_, "psc", "index" );

	setIntInfo( info_, "psc", "indexBC" );


        cand_->setInfo( info_ );

	// OK, store Basic Cluster
        cand_->lock();
        _ecalPSC.push_back( cand_ );
    }
}

void 
EventManager::loadEcalPfBC()
{
  //GHM
  // abort();

  CandList& _ecalPfBC = _candList[ kEcalPfBC ];
  if( _ecalPfBC.size()!=0 ) return; // already loaded

  size_t nEcalPfBC = a().n("pfBc");
  for ( size_t ii=0; ii<nEcalPfBC; ++ii )
    {
        assert( a().load( "pfBc", ii ) );
        float e_ = a().get_f("pfBc", "energy");
        float eta_ = a().get_f("pfBc", "eta");
        float theta_ = KineUtils::theta( eta_ );
        float pt_ = e_*sin(theta_);
        float phi_ = a().get_f("pfBc", "phi");
        int pdgId_ = 22; // photon by construction

        TVector3 p3_(0, 0, 0);
        p3_.SetPtEtaPhi( fabs( pt_), eta_, phi_ );

        Candidate* cand_ = Candidate::create( p3_, pdgId_, _primaryVertex ); 

        // special super cluster information
        CandInfo* info_ = CandInfo::create();
	setFloatInfo( info_, "pfBc", "R4" );
	setFloatInfo( info_, "pfBc", "R9" );
	setIntInfo(   info_, "pfBc", "nHits" );
	setIntInfo(   info_, "pfBc", "index" );
        cand_->setInfo( info_ );
	
        // OK, store Basic Cluster
        cand_->lock();
        _ecalPfBC.push_back( cand_ );

    }
}

void
EventManager::loadEcalRecHits()
{
  //GHM 
  // abort();

  //  cout << "Entering loadEcalRecHits" << endl;

  if( _ecalRecHits.size()!=0 ) return; // already loaded
  
  size_t nRecHits = (size_t) a().n("rh");
  for(size_t irh=0;irh<nRecHits;irh++) {
    assert( a().load("rh", irh) );
    
    ecalRecHit recHit_;
    recHit_.e = a().get_f("rh","energy");
    recHit_.ix = a().get_i("rh","ix");
    recHit_.iy = a().get_i("rh","iy");
    recHit_.iz = a().get_i("rh","iz");
    recHit_.index = a().get_i("rh","index");
    recHit_.time = a().get_i("rh","time");

    antiSpike aSpike_;
    aSpike_.rFlag = a().get_i( "rh", "recoFlag" );
    aSpike_.chi2 = a().get_f( "rh", "chi2" );
    aSpike_.ootChi2 = a().get_f( "rh", "outOfTimeChi2" );
    aSpike_.ootE = a().get_f( "rh", "outOfTimeEnergy" );
//     aSpike_.sevLvl = a().get_i( "rh", "severityLevel" );
//     aSpike_.E1oE9 = a().get_f( "rh", "E1OverE9" );
//     aSpike_.swCross = a().get_f( "rh", "swissCross" );
//     aSpike_.dbStatus = a().get_f( "rh", "dbStatus" );

    recHit_.SpikeVar = aSpike_;

    laser lVar_;
    lVar_.laserCorr = a().get_f( "rh", "laserCorr" );
    lVar_.apdpnref = a().get_f( "rh", "apdpnref" );
    lVar_.apdpn_p2 = a().get_f( "rh", "apdpn_p2" );
    lVar_.apdpn_t2 = a().get_f( "rh", "apdpn_t2" );
    lVar_.alpha = a().get_f( "rh", "alpha" );

    recHit_.LaserVar = lVar_;

    //Store rechits
    _ecalRecHits.push_back( recHit_ );
  }

  //  cout << "done" << endl;
}

void
EventManager::loadEcalPreshowerRecHits()
{

  if( _ecalPSRecHits.size()!=0 ) return; // already loaded

  size_t nPSRecHits = (size_t) a().n("rhps");
  for(size_t irh=0;irh<nPSRecHits;irh++) {
    assert( a().load("rhps", irh) );

    ecalPSRecHit psRecHit_;
    psRecHit_.e = a().get_f("rhps","energy");
    psRecHit_.ix = a().get_i("rhps","ix");
    psRecHit_.iy = a().get_i("rhps","iy");
    psRecHit_.iz = a().get_i("rhps","iz");
    psRecHit_.index = a().get_i("rhps","index");
    psRecHit_.plane = a().get_i("rhps","plane");
    psRecHit_.strip = a().get_i("rhps","strip");

    //Store rechits
    _ecalPSRecHits.push_back( psRecHit_ );
  }
}

void
EventManager::refreshPfCandidates() {
  _candList[kPfCand].clear();
}

void
EventManager::reloadPfCandidates( int v ) {
  loadPfCandidates( v );
}

void
EventManager::loadPfCandidates( int v )
{
  //GHM
  //  abort();

  //  cout << "Entering loadPfCandidates " << endl;

  if( v==-1 ) v = EventServer::pfCandidates_default;  
  size_t ipf_ = kPfCand;
  if( v!=0 ) ipf_ = kPfCand_lowPt;
 
  CandList& _pfCandidates = _candList[ ipf_ ];
  if( _pfCandidates.size()!=0 ) return; // already loaded

  size_t nPfCand_ = a().n("pfc",v);

  for ( size_t ii=0; ii<nPfCand_; ++ii )
    {
      assert( a().load( "pfc", ii, v ) );
      float pt_ =  a().get_f("pfc", "pt", v); 
      float eta_ = a().get_f("pfc", "eta", v); 
      float phi_ = a().get_f("pfc", "phi", v);
      int pdgId_ =  a().get_i("pfc", "pdgId", v);
      
      //      cout << "ii=" << ii << " pt=" << pt_ << endl;
      // if(pt_==0) continue; //??? //MM  continue ===> break segfault
      TVector3 p3_(0, 0, 0);
      p3_.SetPtEtaPhi( fabs( pt_), eta_, phi_ );
      
      float vx_  = a().get_f("pfc","vx", v);
      float vy_  = a().get_f("pfc","vy", v);
      float vz_  = a().get_f("pfc","vz", v);
      Vertex* v3_ = Vertex::create( TVector3( vx_, vy_, vz_ ) );

      //Vertex* vOff_ =  vertexMatching(v3_); //MM new "official matching"
      int vtxIndex=-1;
      if(!a().get("pfc","vtxIndex", vtxIndex,  v) )
	vtxIndex = -1;
      //= a().get_i("pfc","vtxIndex", v);
      Vertex* vOff_;


    //   cout<<" new PFCandidate ============== "<<endl;
//       cout<<" basic matching ";
//       vOff_->print( cout );
      if(vtxIndex!=-1) vOff_ = _vertices[ vtxIndex ];
      else vOff_ =  vertexMatching(v3_);
   //    cout<<" new matching "<<vtxIndex<<endl;
//       vOff_->print( cout );
      

      //Vertex* vOff_ =  v3_;
      Candidate* cand_ = Candidate::create( p3_, pdgId_, vOff_ ); 
      
      // special Pf Candidate informations
      CandInfo* info_ = CandInfo::create();
      setIntInfo(   info_, "pfc", "index", v);
      setFloatInfo( info_, "pfc", "vx", v );
      setFloatInfo( info_, "pfc", "vy", v );
      setFloatInfo( info_, "pfc", "vz", v );
      setFloatInfo( info_, "pfc", "ecalE", v);
      setFloatInfo( info_, "pfc", "hcalE", v);
      setFloatInfo( info_, "pfc", "pS1E", v);
      setFloatInfo( info_, "pfc", "pS2E", v);
      setFloatInfo( info_, "pfc", "deltaP", v);
      setFloatInfo( info_, "pfc", "mva_e_mu", v);
      setFloatInfo( info_, "pfc", "mva_e_pi", v);
      setFloatInfo( info_, "pfc", "mva_pi_mu", v);
      setFloatInfo( info_, "pfc", "mva_nothing_gamma", v);
      setFloatInfo( info_, "pfc", "mva_nothing_nh", v);
      setFloatInfo( info_, "pfc", "mva_gamma_nh", v);
      
      setFloatInfo( info_, "pfc", "caloRho", v);
      setFloatInfo( info_, "pfc", "caloEta", v);
      setFloatInfo( info_, "pfc", "caloPhi", v);

      setFloatInfo( info_, "pfc", "vx", v);
      setFloatInfo( info_, "pfc", "vy", v);
      setFloatInfo( info_, "pfc", "vz", v);
      
      //vertex ID
      setIntInfo(   info_, "pfc", "vtxIndex", v);

      //By default, pfCandidates are not matched with a jet, computing the jet-pfcand map 
      //determines in the pfcandidate participates to a jet
      info_->setBool( "unClus", 1 );

      cand_->setInfo( info_ );
	
      // OK, store Pf Candidate
      cand_->lock();
      _pfCandidates.push_back( cand_ );

    }  

}

const CandList&
EventManager::uncPfCandidates() {

  if(_uncPfCands.size() !=0)
    return _uncPfCands;
  
  //Need the jet-pfcand map
  loadJetPfCandMaps();

  //Need pfcandidates high and low pt
  //First high pt
  const CandList& pfHList_ = pfCandidates();
  for( size_t ii=0; ii<pfHList_.size(); ii++ )
    {
      Candidate* pfcand = pfHList_[ii];
      CandInfo* info = pfcand->info();
      if( info->getBool( "unClus") ) {
	_uncPfCands.push_back( 	pfcand );		      
      }
    }
      
  //FIXME low pf cands discarded for the moment 
  //(low processing rate)
  //MM
  
  //Now low pt pfcands
  //refreshPfCandidates(); //to be sure we take the good one
 //  int v=1; //pfCand_lowPt
//   reloadPfCandidates( v ); //protection
//   const CandList& pfLList_ = candList( EventManager::kPfCand_lowPt );
//   for( size_t ii=0; ii<pfLList_.size(); ii++ )
//     {
//       Candidate* pfcand = pfLList_[ii];
//       CandInfo* info = pfcand->info();
//       if( info->getBool( "unClus") ) {
// 	_uncPfCands.push_back( 	pfcand );		      
//       }
//     }

  return _uncPfCands;
}





void 
EventManager::loadJets( int jettype, int v )
{
  if( jettype==kCaloJet )
    loadCaloTowers(); // make sure CaloTowers are loaded

  pair<int,int> ijet(jettype,v) ;
  CandList& jetList_ = _jetList[ ijet ];
  if( jetList_.size()!=0 ) return; // already loaded

  string prefix;
  if(      jettype==kPatJet )  
    {
      prefix="j";
      v = -1;
    }
  else if( jettype==kPfJet )   
    {
      prefix="pfj";
      if( v==-1 ) v = EventServer::pfJets_default;
    }
  else if( jettype==kTrkJet )   
    {
      prefix="trkj";
      if( v==-1 ) v = EventServer::trackJets_default;
    }
  else if( jettype==kCaloJet ) 
    {
      prefix="cj";
      if( v==-1 ) v = EventServer::caloJets_default;
    }
  else if( jettype==kGenJet )  
    {
      prefix="gj";
      if( v==-1 ) v = EventServer::genJets_default;
    }

  // loop over jets
  size_t nJets = a().n(prefix,v);

  for( size_t ii=0; ii<nJets; ii++ )
    {
      assert( a().load( prefix, ii, v ) );
      float et_  = 
	a().get_f(prefix, "pt", v );
      // (jettype==kPatJet||jettype==kPfJet) ? "pt":"Et"); //!!!FIXME!!!
      if( fabs(et_)<=0 ) continue;
      float eta_ = a().get_f(prefix,"eta",v);
      float phi_ = a().get_f(prefix,"phi",v);
      int pdgId_ = 21; // by default, set as gluon
      TVector3 p3_(0,0,0);
      p3_.SetPtEtaPhi( fabs( et_), eta_, phi_ );

      Vertex* vOff_ = findVertex(eta_, phi_);
      //Vertex* vOff_ = 0;
      Candidate* cand_ = Candidate::create( p3_, pdgId_, vOff_ ); 
      
      size_t uid_ = cand_->uid();

      // special caloJet information
      CandInfo* info_ = CandInfo::create();
      cand_->setInfo( info_  );

      if( jettype==kCaloJet )
	{
	  setIntInfo(   info_, prefix,  "index"      , v );
	  setIntInfo(   info_, prefix,  "nCT"        , v );
	  setFloatInfo( info_, prefix,  "emFraction" , v );
	  setFloatInfo( info_, prefix,  "towersArea" , v );
	  
	  // now the maps
	  size_t nCT_ = (size_t) a().get_i(prefix,"nCT", v );

	  for( size_t jj=0; jj<nCT_; jj++  )
	    {
	      int CTIndex = a().get_vi(prefix,"indexCT",jj, v );
	      const Candidate* ctCand_ = _caloTowerByIndex[ CTIndex ];
	      if( ctCand_==0 ) continue;
	      
	      int ctUid_ = ctCand_->uid();
	      _caloJetTowerMap.insert( make_pair( uid_, ctUid_ ) );
	      
	      isoDep iso_;
	      iso_.val = a().get_vf(prefix,"ptCT",jj, v );
	      iso_.eta = a().get_vf(prefix,"etaCT",jj, v );
	      iso_.phi = a().get_vf(prefix,"phiCT",jj, v );	  
	      float deta_ = iso_.eta - cand_->eta();
	      float dphi_ = iso_.phi - cand_->phi();
	      if( dphi_>=Constants::pi ) dphi_ -= Constants::twopi;
	      if( dphi_<-Constants::pi ) dphi_ += Constants::twopi;
	      iso_.dr = sqrt( deta_*deta_ + dphi_*dphi_ );
	      _isoMap["caloJet-caloTower"].insert( make_pair( uid_, iso_ ) );
	    }
	}
      else if( jettype==kPatJet )
	{
	  setIntInfo(   info_, prefix,  "index"      , v );
	  setIntInfo(   info_, prefix,  "partonFlavour" , v );
	  setFloatInfo( info_, prefix,  "charge" , v );
	  setFloatInfo( info_, prefix,  "corrFactor" , v );
	  setFloatInfo( info_, prefix,  "gj_pt" , v );
	  setFloatInfo( info_, prefix,  "gj_eta" , v );
	  setFloatInfo( info_, prefix,  "gj_phi" , v );
	  setFloatInfo( info_, prefix,  "gp_pt" , v );
	  setFloatInfo( info_, prefix,  "gp_eta" , v );
	  setFloatInfo( info_, prefix,  "gp_phi" , v );
	  setFloatInfo( info_, prefix,  "btagTkC" , v );
	  setFloatInfo( info_, prefix,  "btagSoftM" , v );
	  setFloatInfo( info_, prefix,  "btagTkC" , v );
	  setFloatInfo( info_, prefix,  "btagSoftM" , v );
	  setFloatInfo( info_, prefix,  "btagJPB" , v );
	  setFloatInfo( info_, prefix,  "btagSVtx" , v );
	  setFloatInfo( info_, prefix,  "btagCSVBT" , v );
	  setFloatInfo( info_, prefix,  "btagCSVMVA" , v );

	  setBoolInfo(  info_, prefix,  "mvaIdLoose" , v );
	  setBoolInfo(  info_, prefix,  "mvaIdMedium" , v );
	  setBoolInfo(  info_, prefix,  "mvaIdTight" , v );
	  setFloatInfo( info_, prefix,  "mvaJetId" , v );
	  
	}
      else if( jettype==kPfJet )
	{
	  setIntInfo(   info_, prefix,  "index"      , v );
	  setFloatInfo( info_, prefix,  "maxD" , v );
	  setFloatInfo( info_, prefix,  "etaPhiMom" , v );
	  setFloatInfo( info_, prefix,  "phiPhiMom" , v );
	  setFloatInfo( info_, prefix,  "area" , v );
	  setFloatInfo( info_, prefix,  "chargedHadronEnergyFraction" , v );
	  setFloatInfo( info_, prefix,  "neutralHadronEnergyFraction" , v );
	  setFloatInfo( info_, prefix,  "chargedEmEnergyFraction" , v );
	  setFloatInfo( info_, prefix,  "chargedMuEnergyFraction" , v );
	  setFloatInfo( info_, prefix,  "neutralEmEnergyFraction" , v );
	  setFloatInfo( info_, prefix,  "chargedMultiplicity" , v );
	  setFloatInfo( info_, prefix,  "neutralMultiplicity" , v );
	  setFloatInfo( info_, prefix,  "muonMultiplicity" , v );

	  setFloatInfo (info_, prefix, "jetResolution", v );

	  // setIntInfo(   info_, prefix,  "nJC" , v );
	  //!!!FIXME!!!
	}
      else if( jettype==kTrkJet )
	{
	  setIntInfo(   info_, prefix,  "index"      , v );
	  setFloatInfo( info_, prefix,  "maxD" , v );
	}
      else if( jettype==kGenJet )
	{
	}

      // OK, store jet candidate
      //      cand_->lock(); // will be done after vertex is set
      jetList_.push_back( cand_ );            
    }
}

void 
EventManager::loadMet( int mettype, int v )
{
  pair<int,int> imet(mettype,v) ;

  if( _met.count(imet)>0 ) return; // already loaded

  string prefix;
  if(      mettype==kPatMet )  
    {
      prefix="patmet";
      if( v==-1 ) v = EventServer::patMET_default;
    }
  else if( mettype==kPfMet )   
    {
      prefix="pfmet";
      if( v==-1 ) v = EventServer::pfMET_default;
    }
  else if( mettype==kCaloMet ) 
    {
      prefix="cmet";
      if( v==-1 ) v = EventServer::caloMET_default;
    }
  else if( mettype==kGenMet )  
    {
      prefix="gmet";
      if( v==-1 ) v = EventServer::genMET_default;
    }
  else assert(0);
  a().load( prefix, 0, v );
  TVector2 p2(0.,0.);
  float Et_  = a().get_f(prefix,"pt",v);
  float phi_ = a().get_f(prefix,"phi",v);
  p2.Set( Et_*cos(phi_), Et_*sin(phi_) );      
  Candidate* cand_ = Candidate::create( p2 );
  CandInfo* info_ = CandInfo::create();
  cand_->setInfo( info_ );

  setFloatInfo( info_, prefix,  "sumEt" , v );
  setFloatInfo( info_, prefix,  "mEtSig" , v );
  if( mettype==kGenMet )
    {
      setFloatInfo( info_, prefix,  "em" , v );
      setFloatInfo( info_, prefix,  "had" , v );
      setFloatInfo( info_, prefix,  "inv" , v );
    }
  if( mettype==kPfMet || mettype==kPatMet )
    {
      setFloatInfo( info_, prefix,  "cHadEF" , v );
      setFloatInfo( info_, prefix,  "nHadEF" , v );
      setFloatInfo( info_, prefix,  "cEmEF" , v );
      setFloatInfo( info_, prefix,  "nEmEF" , v );
      setFloatInfo( info_, prefix,  "cMuEF" , v );
      
      if( mettype==kPatMet ) {
	setFloatInfo( info_, prefix,  "emEBF"  , v );
	setFloatInfo( info_, prefix,  "emEEF"  , v );
	setFloatInfo( info_, prefix,  "emHFF"  , v );
	setFloatInfo( info_, prefix,  "hadHBF" , v );
	setFloatInfo( info_, prefix,  "hadHEF" , v );
	setFloatInfo( info_, prefix,  "hadHFF" , v );
	setFloatInfo( info_, prefix,  "hadHOF" , v );
	setFloatInfo( info_, prefix,  "hadF"   , v );
      }

    }

  // OK, store the met candidate
  cand_->lock();
  _met[imet] = cand_;
}

void 
EventManager::loadElectronSCMaps()
{
  // check if not already loaded
  if( _candAssoc.count("electron-SuperCluster")!=0 ) return;

  static int idbg_=0;
  if( idbg_++<1 )
    {
      cout << "Warning -- you are loading the electron-SC maps!" << endl;
    }

  // register the maps
  _candAssoc.insert( make_pair( "electron-SuperCluster",       CandAssoc() ) );
  _candAssoc.insert( make_pair( "electron-PfSuperCluster",     CandAssoc() ) );
  _candAssoc.insert( make_pair( "SuperCluster-PfSuperCluster", CandAssoc() ) );
  _candAssoc.insert( make_pair( "electron-GsfTrack",       CandAssoc() ) );

  // dependencies
  loadEcalSC();
  loadEcalPfSC();
  loadGsfTracks();

  const CandList& elList = electrons();
  for( size_t ii=0; ii<elList.size(); ii++ )
    {
      assert( a().load( "el", ii ) );
      Candidate* cand_= elList[ii];

      int scIndex_   = a().get_i("el", "indexSC" );
      int pfScIndex_ = a().get_i("el", "indexPfSC" );
    
      //ECAL SuperCluster map
      //Need Ecal super cluster
      Candidate* egsc_(0);
      if(scIndex_ != -100) 
	{
	  egsc_ = _candList[ kEcalSC ][ scIndex_ ];
	  _candAssoc["electron-SuperCluster"].insert( cand_, egsc_ );
	}
	  
      //PFSuperCluster map
      Candidate* pfsc_(0);
      if( pfScIndex_!=-1 )
	{
	  pfsc_ = _candList[ kEcalPfSC ][ pfScIndex_ ];
	}
      else if( egsc_!=0 )
	{
	  const CandList& pfScList = _candList[ kEcalPfSC ];
	  float drmin(0);
	  pfsc_ = CandUtil::closest( egsc_, pfScList, drmin, 0., true );
	  if( drmin>0.1 ) pfsc_=0;
	}
      if( pfsc_!=0 )
	_candAssoc["electron-PfSuperCluster"].insert( cand_, pfsc_ );
      if( pfsc_!=0 && egsc_!=0 )
	_candAssoc["SuperCluster-PfSuperCluster"].insert( egsc_, pfsc_ );

      //GsfTrack map
      int gsfTrkIndex_ = a().get_i("el", "indexGsfTrack" );
      if( gsfTrkIndex_ != -1 ) 
	{
	  Candidate* gsfTrk_ = _candList[ kGsfTrack ][ gsfTrkIndex_ ];
	  _candAssoc["electron-GsfTrack"].insert( cand_, gsfTrk_ );
	}
    }	  
}

void 
EventManager::loadElectronIsoDepMaps()
{
  static int idbg_=0;
  if( idbg_++<1 )
    {
      cout << "Warning -- you are loading the electron-IsoDep maps!" << endl;
    }
}

void 
EventManager::loadElectronTrackMaps()
{
  // check if not already loaded
  if( _candAssoc.count("electron-Track")!=0 ) return;

  static int idbg_=0;
  if( idbg_++<1 )
    {
      cout << "Warning -- you are loading the electron-Track maps!" << endl;
    }

  // register the maps
  _candAssoc.insert( make_pair( "electron-Track",       CandAssoc() ) );

  // dependencies
  loadTracks();

  const CandList& elList = electrons();
  for( size_t ii=0; ii<elList.size(); ii++ )
    {
      assert( a().load( "el", ii ) );
      Candidate* cand_= elList[ii];

      int trkIndex_ = a().get_i("el", "indexTrack" );
      Candidate* trkCand_(0);
      if( trkIndex_!=-1 )
	{
	  trkCand_ = _candList[ kTrack ][ trkIndex_ ];
	}
      else
	{
	  const CandList& trkList = tracks();
	  float drmin(0);
	  trkCand_ = CandUtil::closest( cand_, trkList, drmin, 0., true ); 
	}
      if( trkCand_!=0 )
	{
	  _candAssoc["electron-Track"].insert( cand_, trkCand_ );
	}      
    }
}

void 
EventManager::loadElectronPfCandMaps()
{
  // check if not already loaded
  if( _candAssoc.count("electron-PfCand")!=0 ) return;

  static int idbg_=0;
  if( idbg_++<1 )
    {
      cout << "Warning -- you are loading the electron-PfCand maps!" << endl;
    }

  // register the maps
  _candAssoc.insert( make_pair( "electron-PfCand",       CandAssoc() ) );

  // dependencies
  loadPfCandidates();

  const CandList& elList = electrons();
  for( size_t ii=0; ii<elList.size(); ii++ )
    {
      assert( a().load( "el", ii ) );
      Candidate* cand_= elList[ii];

      int pfCandIndex_ = a().get_i("el", "indexPfCandidate" );

      Candidate* pfCand_(0);
      if( pfCandIndex_!=-1 )
	{
	  pfCand_ = _candList[ kPfCand ][ pfCandIndex_ ];
	}
      else
	{
	  const CandList& pfList = pfCandidates();
	  float drmin(0);
	  pfCand_ = CandUtil::closest( cand_, pfList, drmin, 0., true ); 
	  if( drmin>0.1 ) pfCand_=0;
	}
      if( pfCand_!=0 )
	{
	  _candAssoc["electron-PfCand"].insert( cand_, pfCand_ );
	}      
    }
}

void 
EventManager::loadEcalSCBCMaps()
{
  // check if not already loaded
  if( _candMap.count("SuperCluster-BasicCluster")!=0 ) return;

  static int idbg_=0;
  if( idbg_++<1 )
    {
      cout << "Warning -- you are loading the SC-BC maps!" << endl;
    }

  // register the maps
  _candMap.insert( make_pair( "SuperCluster-BasicCluster", CandMap() ) );
  _candMap.insert( make_pair( "PfSuperCluster-PfBasicCluster", CandMap() ) );

  // dependencies
  loadEcalSC();
  loadEcalBC();
  loadEcalPfSC();
  loadEcalPfBC();

  const CandList& scList = _candList[kEcalSC];
  for( size_t ii=0; ii<scList.size(); ii++ )
    {
      assert( a().load( "sc", ii ) );
      Candidate* cand_ = scList[ii];
      std::vector<Int_t> bClus_ = a().get_vi("sc", "indexBasicCluster" );
      for( size_t jj=0; jj<bClus_.size(); jj++ ) 
	{	      
	  size_t bcIndex_ = bClus_[jj];
	  Candidate* bc_ = _candList[ kEcalBC ][ bcIndex_ ];
	  _candMap["SuperCluster-BasicCluster"].insert( make_pair( cand_, bc_  ) );
	}      
    }

  const CandList& pfscList = _candList[kEcalPfSC];
  for( size_t ii=0; ii<pfscList.size(); ii++ )
    {
      assert( a().load( "pfSc", ii ) );
      Candidate* cand_ = pfscList[ii];
      std::vector<Int_t> bClus_ = a().get_vi("pfSc", "indexBasicCluster" );
      for( size_t jj=0; jj<bClus_.size(); jj++ ) 
	{	      
	  size_t bcIndex_ = bClus_[jj];
	  Candidate* bc_ = _candList[ kEcalPfBC ][ bcIndex_ ];
	  _candMap["PfSuperCluster-PfBasicCluster"].insert( make_pair( cand_, bc_  ) );
	}      
    }
}

void 
EventManager::loadEcalBCRecHitMaps()
{
  // check if not already loaded
  if( _ecalRecHitMap.count("BasicCluster-EcalRecHit")!=0 ) return;

  static int idbg_=0;
  if( idbg_++<1 )
    {
      cout << "Warning -- you are loading the BC-RecHit maps!" << endl;
    }

  // register the maps
  _ecalRecHitMap.insert( make_pair( "BasicCluster-EcalRecHit", 
				    CandIdRecHitMap() ) );
  _ecalRecHitMap.insert( make_pair( "PfBasicCluster-EcalRecHit", 
				    CandIdRecHitMap() ) );

  // dependencies
  loadEcalRecHits();

  const CandList& bcList = _candList[ kEcalBC ];
  for( size_t ii=0; ii<bcList.size(); ii++ )
    {
      assert( a().load( "bc", ii ) );
      Candidate* cand_ = bcList[ii];

      std::vector<Int_t> eRecHits_ = a().get_vi("bc", "indexHit" );
      size_t uid_ = cand_->uid();
      for( size_t jj=0; jj<eRecHits_.size(); jj++ ) 
	{	  
	  size_t erhIndex_ = eRecHits_[jj];
	  ecalRecHit rh_ = _ecalRecHits[ erhIndex_ ];
	  std::pair< ecalRecHit, float > rechit_ (rh_, 1 );
	  _ecalRecHitMap["BasicCluster-EcalRecHit"].insert( make_pair( uid_, rechit_  ) );
	}
    }

  //GHM -- fixme (if needed)
//   const CandList& pfbcList = _candList[ kEcalPfBC ];
//   for( size_t ii=0; ii<pfbcList.size(); ii++ )
//     {
//       assert( a().load( "pfBc", ii ) );
//       Candidate* cand_ = pfbcList[ii];

//       std::vector<Int_t> eRecHits_ = a().get_vi("pfBc", "indexHit" );
//       std::vector<Float_t> fracHits_ = a().get_vf("pfBc", "fractionHit" );
//       size_t uid_ = cand_->uid();
//       for( size_t jj=0; jj<eRecHits_.size(); jj++ ) 
// 	{	  
// 	  size_t erhIndex_ = eRecHits_[jj];
// 	  ecalRecHit rh_ = _ecalRecHits[ erhIndex_ ];
// 	  float fraction_ =  fracHits_[jj];
// 	  std::pair<ecalRecHit, float> rechit_ (rh_, fraction_ );
// 	  _ecalRecHitMap["PfBasicCluster-EcalRecHit"].insert( make_pair( uid_, rechit_ ) );
// 	}
//     }

}

void 
EventManager::loadEcalBCPSCMaps()
{

  // check if not already loaded
  if( _candAssoc.count("PreshowerCluster-BasicCluster")!=0 ) return;

  static int idbg_=0;
  if( idbg_++<1 )
    {
      cout << "Warning -- you are loading the SC-PSC maps!" << endl;
    }

  // register the maps
  _candAssoc.insert( make_pair( "PreshowerCluster-BasicCluster", CandAssoc() ) );
 
  // dependencies
  loadEcalPSC();
  loadEcalBC();

  const CandList& pscList = _candList[kEcalPSC];
  for( size_t ii=0; ii<pscList.size(); ii++ )
    {
      assert( a().load( "psc", ii ) );
      Candidate* cand_ = pscList[ii];
      int bcIndex_ = a().get_i("psc", "indexBC" );
     
      Candidate* bc_ = _candList[ kEcalBC ][ bcIndex_ ];
      _candAssoc["PreshowerCluster-BasicCluster"].insert( cand_, bc_ );
    }      
}

void 
EventManager::loadEcalPSCPSRecHitMaps()
{
  // check if not already loaded
  if( _ecalPSRecHitMap.count("PreshowerCluster-EcalPSRecHit")!=0 ) return;
  
  static int idbg_=0;
  if( idbg_++<1 )
    {
      cout << "Warning -- you are loading the PSC-PSRecHit maps!" << endl;
    }

  // register the maps
  _ecalPSRecHitMap.insert( make_pair( "PreshowerCluster-EcalPSRecHit", 
				      CandIdPSRecHitMap() ) );
  
  // dependencies
  loadEcalPreshowerRecHits();
  loadEcalPSC();

  const CandList& pscList = _candList[ kEcalPSC ];
  for( size_t ii=0; ii<pscList.size(); ii++ )
    {
      assert( a().load( "psc", ii ) );
      Candidate* cand_ = pscList[ii];

      std::vector<Int_t> ePSRecHits_ = a().get_vi("psc", "indexHit" );
      size_t uid_ = cand_->uid();
    
      for( size_t jj=0; jj<ePSRecHits_.size(); jj++ ) 
	{	  
	  size_t epsrhIndex_ = ePSRecHits_[jj];
	  ecalPSRecHit psrh_ = _ecalPSRecHits[ epsrhIndex_ ];
	  std::pair< ecalPSRecHit, float > psrechit_ (psrh_, 1 );
	  _ecalPSRecHitMap["PreshowerCluster-EcalPSRecHit"].insert( make_pair( uid_, psrechit_  ) );
	}
    }
}

void 
EventManager::loadMuonTrackMaps()
{
  // check if not already loaded
  if( _candAssoc.count("muon-Track")!=0 ) return;

  static int idbg_=0;
  if( idbg_++<1 )
    {
      cout << "Warning -- you are loading the muon-Track maps!" << endl;
    }

  // register the maps
  _candAssoc.insert( make_pair( "muon-Track",       CandAssoc() ) );

  // dependencies
  loadTracks();

  const CandList& muList = muons();
  for( size_t ii=0; ii<muList.size(); ii++ )
    {
      assert( a().load( "mu", ii ) );
      Candidate* cand_= muList[ii];

      int trkIndex_ = a().get_i("mu", "indexTrack" );
      Candidate* trkCand_(0);
      if( trkIndex_!=-1 )
	{
	  trkCand_ = _candList[ kTrack ][ trkIndex_ ];
	}
      else
	{
	  const CandList& trkList = tracks();
	  float drmin(0);
	  trkCand_ = CandUtil::closest( cand_, trkList, drmin, 0., true ); 
	}
      if( trkCand_!=0 )
	{
	  _candAssoc["muon-Track"].insert( cand_, trkCand_ );
	}      
    }
}

void 
EventManager::loadMuonPfCandMaps()
{
  // check if not already loaded
  if( _candAssoc.count("muon-PfCand")!=0 ) return;

  static int idbg_=0;
  if( idbg_++<1 )
    {
      cout << "Warning -- you are loading the muon-PfCand maps!" << endl;
    }

  // register the maps
  _candAssoc.insert( make_pair( "muon-PfCand",       CandAssoc() ) );

  // dependencies
  loadPfCandidates();

  const CandList& muList = muons();
  for( size_t ii=0; ii<muList.size(); ii++ )
    {
      assert( a().load( "mu", ii ) );
      Candidate* cand_= muList[ii];

      int pfCandIndex_ = a().get_i("mu", "indexPfCandidate" );

      Candidate* pfCand_(0);
      if( pfCandIndex_!=-1 )
	{
	  pfCand_ = _candList[ kPfCand ][ pfCandIndex_ ];
	}
      else
	{
	  const CandList& pfList = pfCandidates();
	  float drmin(0);
	  pfCand_ = CandUtil::closest( cand_, pfList, drmin, 0., true ); 
	  if( drmin>0.1 ) pfCand_=0;
	}
      if( pfCand_!=0 )
	{
	  _candAssoc["muon-PfCand"].insert( cand_, pfCand_ );
	}      
    }
}

void
EventManager::loadJetPfCandMaps() 
{
  // dependencies
  // new GHM -- maps PatJets-Pfjets, PatJets-leptons...
  loadJets( EventManager::kPatJet , -1 );
  loadPfCandidates(); //High Pt cands
  loadElectronPfCandMaps();
  loadMuonPfCandMaps();
  loadJetRelatedMaps();

  // check if not already loaded
  if( _candMap.count("Jet-PfCand")!=0 ) return;

  static int idbg_=0;
  if( idbg_++<1 )
    {
      cout << "Warning -- you are loading the Jet-PFC maps!" << endl;
    }

  // register the map
  _candMap.insert( make_pair( "Jet-PfCand", CandMap() ) );
  
 
  const CandList& patJetList = jetList( kPatJet, -1 );
  for( size_t ii=0; ii<patJetList.size(); ii++ )
    {
      assert( a().load( "j", ii ) );
      Candidate* cand_ = patJetList[ii];
      
      //      cout << "+++++" << endl;
      //      cand_->oneLine(cout);

      //High pt pf cands
      std::vector<Int_t> pfH_ = a().get_vi("j", "indexPfConstHPt" );
      
      for( int kk=0; kk<2; kk++ )
	{
	  for( size_t jj=0; jj<pfH_.size(); jj++ ) 
	    {	      
	      size_t pfHIndex_ = pfH_[jj];
	      Candidate* pfhc_ = _candList[ kPfCand ][ pfHIndex_ ];
	      bool isNeutral = pfhc_->charge()==0;
	      if( kk==0 &&  isNeutral ) continue;
	      if( kk==1 && !isNeutral ) continue;
	      //	  pfhc_->oneLine(cout);
	      _candMap["Jet-PfCand"].insert( make_pair( cand_, pfhc_  ) );
	      
	      //Easy solution to know if a pfCandidate is in one jet
	      CandInfo* infoh_ = pfhc_->info();
	   
	      infoh_->setBool( "unClus", 0 );	      
	    }
	}
    }

  // GHM 05/12 -- select the clean jets
  const CandAssoc& PatJetPfJetAssoc  = candAssoc( "PatJet-PfJet"  ); 
  const CandAssoc& PatJetLeptonAssoc = candAssoc( "PatJet-lepton" ); 
  const CandAssoc& electronPfAssoc   = candAssoc( "electron-PfCand" ); 
  const CandAssoc& muonPfAssoc       = candAssoc( "muon-PfCand" ); 
   
  float pt_min  = 10;
  //float eta_max =  5;
  //  const CandList& patJetList = jetList( kPatJet,-1 );
  for( size_t ii=0; ii<patJetList.size(); ii++ )
    {
      Candidate* jet    = patJetList[ii];
      CandInfo* jInfo = jet->info();
      jInfo->setBool("isCleanJet", false );
      
      if( jet->pt() < pt_min ) continue;

      //bool isCentral = fabs( jet->eta() )<eta_max;
      //bool isForward = !isCentral; //FIXME MM not used in the code
      
      Candidate* pfJet  = PatJetPfJetAssoc.getSecond( jet );
      float dR_pfJet=-1;
      if( pfJet!=0 )
	{
	  dR_pfJet =  CandUtil::dR( jet, pfJet );
	}
      else
	{
	  cout << "Warning, no associated pfJet" << endl;
	}
      Candidate* lepton = PatJetLeptonAssoc.getSecond( jet );
      Candidate* pfLepton(0);
      
      // list the constituents
      const CandMap& jetConst_ = candMap("Jet-PfCand");	  
      CandMapIterator it_lower_ = jetConst_.lower_bound( jet );
      CandMapIterator it_upper_ = jetConst_.upper_bound( jet );
      CandMapIterator it;
      
      map< int, CandList > chargedAtVertex;
      CandList neutralConstituents;
      for( it=it_lower_; it!=it_upper_; it++ )
	{ 
	  Candidate* cand = it->second;
	  bool isNeutral = cand->charge()==0;
	  if( isNeutral )
	    {
	      neutralConstituents.push_back(cand);
	    }
	  else
	    {
	      int vtxId = cand->vertexIndex();
	      chargedAtVertex[vtxId].push_back( cand );
	    }
	}
      
      CandList chargedComponents;
      Candidate* neutralComponent(0);
      for( map< int, CandList >::iterator itc=chargedAtVertex.begin();
	   itc!=chargedAtVertex.end(); ++itc )
	{
	  Candidate* sum_ = Candidate::create( itc->second );
	  chargedComponents.push_back( sum_ );
	}
      sort( chargedComponents.begin(), chargedComponents.end(), Candidate::SortByPt() );
      if( neutralConstituents.size()>0 )
	{
	  neutralComponent = Candidate::create( neutralConstituents );
	}
      
      bool isNeutral = chargedComponents.size()==0;
      bool isCharged = !isNeutral;
      //bool isAtPrimary = true;
      float sumPt(0);
      int jetVtxId = -1;      
      Vertex* jetVtx(0);
      if( isCharged ) 
	{
	  sumPt = chargedComponents[0]->pt();
	  jetVtx = chargedComponents[0]->vertex();
	  assert( jetVtx!=0 );
	  jetVtxId = jetVtx->index();
	} 
      jet->setVertex(jetVtx);
      jet->lock();
      
      // find the closest lepton candidate
      float dR_lJet=-1;
      bool isIsolatedLeptonJet = false;
      int leptonVtxId = -1;
      if( lepton!=0 )
	{
	  dR_lJet =  CandUtil::dR( jet, lepton );
	  pfLepton = electronPfAssoc.getSecond( lepton );
	  if( pfLepton==0 ) pfLepton = muonPfAssoc.getSecond( lepton );	      
	  leptonVtxId = lepton->vertexIndex();
	  
	  if( isCharged )
	    {
	      if( dR_lJet>0 && leptonVtxId==jetVtxId )
		{
		  // this is a lepton at the primary vertex
		  if( dR_lJet<0.1 && lepton->pt()>0.80*sumPt )
		    {
		      // the lepton take more than 80% of the jet energy
		      isIsolatedLeptonJet = true;
		      if( _cleanLeptons.count( leptonVtxId ) )
			// clean leptons are lepton not within a jet
			_cleanLeptons[leptonVtxId].push_back( lepton );
		    }
		}
	    }
	}
      
      jInfo->setBool("CJisCharged",isCharged);
      jInfo->setFloat("CJsumPt",sumPt);
      jInfo->setBool("CJisIsolatedLepton",false);
      if( isIsolatedLeptonJet )
	{
	  jInfo->setBool("CJisIsolatedLepton",true);
	  _leptonJets[jetVtxId].push_back( jet );
	}      
      else
	{
	  jInfo->setBool("isCleanJet", true );
	  _cleanJets[jetVtxId].push_back( jet );
	}
    }
}

void
EventManager::loadJetRelatedMaps() 
{
  //GHM -- 04/18/2012

  // check if not already filled
  if( _candAssoc.count("PatJet-PfJet")!=0 ) return;

  static int idbg_=0;
  if( idbg_++<1 )
    {
      cout << "Warning -- you are loading the PatJet-PfJet maps!" << endl;
    }

  // register the map
  _candAssoc.insert( make_pair( "PatJet-PfJet",  CandAssoc() ) );
  _candAssoc.insert( make_pair( "PatJet-lepton", CandAssoc() ) );
  
  // dependencies
  loadElectrons();
  loadMuons();
  loadJets( kPatJet , -1 );
  loadJets( kPfJet  , -1 );
  loadPfCandidates(); //High Pt cands
 
  const CandList& patJetList = jetList( kPatJet );
  const CandList& pfJetList  = jetList( kPfJet  );
  float dR_max=0.03;
  for( size_t ii=0; ii<patJetList.size(); ii++ )
    {
      Candidate* jet = patJetList[ii];

      // find the closest PF jet
      float dR_pfJet=5;
      Candidate* pfJet = CandUtil::closest( jet, pfJetList, dR_pfJet );
      if( pfJet && dR_pfJet<dR_max )
	{
	  _candAssoc["PatJet-PfJet"].insert( jet, pfJet );
	}     

      // find the closest lepton candidate
      float dR_eJet=5; 
      float dR_mJet=5;      
      float dR_lJet=5;      
      Candidate* lJet(0);
      Candidate* eJet = CandUtil::closest( jet, electrons(), dR_eJet );
      Candidate* mJet = CandUtil::closest( jet,     muons(), dR_mJet );
      if( eJet && dR_eJet<dR_mJet )
	{
	  dR_lJet= dR_eJet;
	  lJet = eJet;
	}
      else if( mJet && dR_mJet<0.3 )
	{
	  dR_lJet= dR_mJet;
	  lJet = mJet;	  
	}
      if( lJet!=0 && dR_lJet<0.1 )
	{
	  _candAssoc["PatJet-lepton"].insert( jet, lJet );
	}
    } 
}

void 
EventManager::loadTrackHitMaps()
{
  static int idbg_=0;
  if( idbg_++<1 )
    {
      cout << "Warning -- you are loading the Track-Hit maps!" << endl;
    }
}


void 
EventManager::loadZllCandidates()
{
  for( int imode=0; imode<2; imode++ )
    {
      // need electrons, muons
      string prefix;
      int ilept(0);
      string str_;
      if( imode==EventServer::kZee ) 
	{
	  loadElectrons();
	  prefix = "zee";
	  ilept = kElectron;
	  str_ = "Zee";
	}
      else if( imode==EventServer::kZmm )   
	{
	  loadMuons();
	  prefix = "zmm";
	  ilept = kMuon;
	  str_ = "Zmm";
	}
      else assert(0);
      
      //      const string & str_ = EventServer::Zmode[imode];
      CandList& _Zll = _compositeCandList[ str_ ];
      if( _Zll.size()!=0 ) return; // already loaded
      
      // loop over Z candidates
      size_t nZll = a().n(prefix);
      for( size_t ii=0; ii<nZll; ii++ )
	{
	  assert( a().load( prefix, ii ) );
	  
	  // get the daughters
	  const Candidate* dau_[2];
	  
	  dau_[0] = _candList[ilept][ a().get_i(prefix,"dau_0") ];
	  dau_[1] = _candList[ilept][ a().get_i(prefix,"dau_1") ];
	  Candidate* cand_= Candidate::create( dau_[0], dau_[1] );

	  float mass_  = a().get_f(prefix,"mass");
	  float pt_    = fabs( a().get_f(prefix,"pt") );
	  float theta_ = a().get_f(prefix,"theta");
	  float phi_   = a().get_f(prefix,"phi");
	  float px_    = pt_*cos(phi_);
	  float py_    = pt_*sin(phi_);
	  assert( theta_!=0 );
	  float pz_    = pt_/tan(theta_);
	  int pdgId_   = 23;
	  cand_->setPdgCode( pdgId_ );
	  cand_->setMass( mass_ );
	  cand_->setPxPyPz( px_, py_, pz_ );
	  
	  // special Zll information
	  CandInfo* info_ = CandInfo::create();
	  cand_->setInfo( info_  );

	  setFloatInfo( info_, prefix, "rapidity" );
	  setFloatInfo( info_, prefix, "helicity" );
	  setFloatInfo( info_, prefix, "pt_0"     );
	  setFloatInfo( info_, prefix, "eta_0"    );
	  setFloatInfo( info_, prefix, "phi_0"    );
	  setFloatInfo( info_, prefix, "pt_1"     );
	  setFloatInfo( info_, prefix, "eta_1"    );
	  setFloatInfo( info_, prefix, "phi_1"    );
	  
	  // OK store Zll candidate
	  cand_ -> lock();
	  _Zll.push_back( cand_ );      
	}
    }
}

void 
EventManager::loadMcTruthCandidates()
{

  CandList& _mcTruthCandidates = _candList[ kMcTruthCandidate ];
  if( _mcTruthCandidates.size()!=0 ) return; // already loaded
 
  // loop over mcTruthCandidates
  size_t nMcTruthCandidates = a().n("mc"); //cout<<nMcTruthCandidates<<endl;
  for( size_t ii=0; ii<nMcTruthCandidates; ii++ )
    {
      assert( a().load( "mc", ii ) );
      float pt_  = a().get_f("mc","pt");
      float eta_ = a().get_f("mc","eta");
      float phi_ = a().get_f("mc","phi");
      int pdgId_ = a().get_i("mc","pdgId");

      TVector3 p3_(0,0,0);
      p3_.SetPtEtaPhi( fabs( pt_), eta_, phi_ );

      float vx_  = a().get_f("mc","vx");
      float vy_  = a().get_f("mc","vy");
      float vz_  = a().get_f("mc","vz");
      Vertex* v3_ = Vertex::create( TVector3( vx_, vy_, vz_ ) );
      Candidate* cand_ = Candidate::create( p3_, pdgId_, v3_ ); 
      cand_->setMass( a().get_f("mc","mass") );

      // special mcTruthCandidate information
      CandInfo* info_ = CandInfo::create();
      setIntInfo(   info_, "mc", "nDau"  );
      setIntInfo(   info_, "mc", "motherIndex" );
      setIntInfo(   info_, "mc", "status" );
      setIntInfo(   info_, "mc", "collId" );
      cand_->setInfo( info_  );

      // OK, store mcTruthCandidate candidate
      cand_ -> lock();
      _mcTruthCandidates.push_back( cand_ );      
    }
}

const CandList&
EventManager::candList( int ilist ) const
{
  switch( ilist )
    {
    case kTrack:
      me()->loadTracks();    break;
    case kGsfTrack:
      me()->loadGsfTracks();   break;
    case kElectron:
      me()->loadElectrons();   break;
    case kMuon:
      me()->loadMuons();       break;
    case kTau:
      me()->loadTaus();   break;
    case kCaloTower:
      me()->loadCaloTowers();  break;
    case kPhoton:
      me()->loadPhotons();    break;
    case kEcalPSC:
      me()->loadEcalPSC();     break;
    case kEcalBC:
      me()->loadEcalBC();     break;
    case kEcalSC:
      me()->loadEcalSC();     break;
    case kEcalPfBC:
      me()->loadEcalPfBC();     break;
    case kEcalPfSC:
      me()->loadEcalPfSC();     break;
    case kPfCand:
      me()->loadPfCandidates();     break;
    case kPfCand_lowPt:
      me()->loadPfCandidates( 1 );     break;
    case kMcTruthCandidate:
      me()->loadMcTruthCandidates(); break;
    case kGenEventCandidate:
      break;
    default:
      abort();      
    }
  return _candList[ ilist ];
}

const CandList&
EventManager::jetList( const string& str, int v ) const
{
  int jettype;
  if( str=="" || str=="paf" || str=="Pat" || str=="Pat" ) 
    jettype=kPatJet;
  else if( str=="cj" || str=="Calo" || str=="CaloJet" || str=="caloJet" )    
    jettype = kCaloJet;
  else if( str=="pf" || str=="PF" || str=="ParticleFlow" ) jettype = kPfJet;
  else if( str=="gen" || str=="Gen" ) jettype = kGenJet;
  else { cout << "str=" << str << endl; assert(0); }
  return jetList( jettype, v );
}

const CandList&
EventManager::jetList( int ilist, int v ) const
{
  switch( ilist )
    {
    case kPatJet:
      me()->loadJets( kPatJet, v );    break;
    case kCaloJet:
      me()->loadJets( kCaloJet, v );    break;
    case kPfJet:
      me()->loadJets( kPfJet, v );    break;
    case kTrkJet:
      me()->loadJets( kTrkJet, v );    break;
    case kGenJet:
      me()->loadJets( kGenJet, v );    break;
    default:
      abort();      
    }

  pair<int,int> ijet(ilist,v);
  return me()->_jetList[ ijet ];
}

const CandList&
EventManager::jetList( int ilist, string& str ) const
{
  int v=-1;

  switch( ilist )
    {
    case kPatJet:
      me()->loadJets( kPatJet, v );    break;
    case kCaloJet:
      for(int unsigned i=0;i<a().caloJets.size();i++) {
	if(str ==a().caloJets[i])
	  { v=i; break; }
      }
      me()->loadJets( kCaloJet, v );    break;
    case kPfJet:
      for(int unsigned i=0;i<a().pfJets.size();i++) {
	if(str ==a().pfJets[i])
	  { v=i; break; }
      }
      me()->loadJets( kPfJet, v );    break;
    case kTrkJet:
      for(int unsigned i=0;i<a().trackJets.size();i++) {
	if(str ==a().trackJets[i])
	  { v=i; break; }
      }
      me()->loadJets( kTrkJet, v );    break;
    case kGenJet:
      for(int unsigned i=0;i<a().genJets.size();i++) {
	if(str ==a().trackJets[i])
	  { v=i; break; }
      }
      me()->loadJets( kGenJet, v );    break;
    default:
      abort();      
    }

  pair<int,int> ijet(ilist,v);
  return me()->_jetList[ ijet ];
}


const Candidate* 
EventManager::met(const string& str, const string& name) const {
 int mettype;
 int v=-1;
 
 if( str=="" || str=="paf" || str=="Pat" || str=="pat" || str=="Reco" || str=="reco" ) {
    mettype=kPatMet;
    for(int unsigned i=0;i<a().patMET.size();i++) {
      if(name ==a().patMET[i])
	{ v=i; break; }
    }
 }
 else if( str=="cj" || str=="Calo" || str=="CaloMet" || str=="caloMet" ) {   
   mettype = kCaloMet;
   for(int unsigned i=0;i<a().caloMET.size();i++) {
     if(name ==a().caloMET[i])
       { v=i; break; }
   }
 }
 else if( str=="pf" || str=="PF" || str=="ParticleFlow" ) {
   mettype = kPfMet;
   for(int unsigned i=0;i<a().pfMET.size();i++) {
     if(name ==a().pfMET[i])
       { v=i; break; }
   }
 }
 else if( str=="track" || str=="trk" || str=="Trk" ) {
   mettype = kTrkMet;
   for(int unsigned i=0;i<a().recoMET.size();i++) {
     if(name ==a().recoMET[i])
       { v=i; break; }
   }
 }
 else if( str=="gen" || str=="Gen" ) {
   mettype = kGenMet;
   for(int unsigned i=0;i<a().genMET.size();i++) {
     if(name ==a().genMET[i])
       { v=i; break; }
   }
 }
  else { cout << "str=" << str << endl; assert(0); }
 
 if(v==-1 && name!="") cout<<"Warning, no "<<name<<" MET available, default returned "<<endl; 

  return met( mettype, v );
}

const Candidate*
EventManager::met( const string& str, int v ) const
{
  int mettype;
  if( str=="" || str=="paf" || str=="Pat" || str=="pat" || str=="Reco" || str=="reco" ) 
    mettype=kPatMet;
  else if( str=="cj" || str=="Calo" || str=="CaloMet" || str=="caloMet" )    
    mettype = kCaloMet;
  else if( str=="pf" || str=="PF" || str=="ParticleFlow" ) mettype = kPfMet;
  else if( str=="track" || str=="trk" || str=="Trk" ) mettype = kTrkMet;
  else if( str=="gen" || str=="Gen" ) mettype = kGenMet;
  else { cout << "str=" << str << endl; assert(0); }
  return met( mettype, v );
}

const Candidate*
EventManager::met( int ilist, int v ) const
{
  switch( ilist )
    {
    case kPatMet:
      me()->loadMet( kPatMet, v );    break;
    case kCaloMet:
      me()->loadMet( kCaloMet, v );    break;
    case kTrkMet:
      me()->loadMet( kTrkMet, v );    break;
    case kPfMet:
      me()->loadMet( kPfMet, v );    break;
    case kGenMet:
      me()->loadMet( kGenMet, v );    break;
    default:
      abort();      
    }
  pair<int,int> imet(ilist,v);
  Candidate* cand=me()->_met[imet];
  return cand;
}

const CandIdIsoMap& 
EventManager::isoMap( const std::string& str ) const
{
  if( str=="track-hits" ) me()->loadTrackHits();

  map<string,CandIdIsoMap>::const_iterator it = _isoMap.find( str );
  if( it==_isoMap.end() )
    {
      cout << str << endl;
      assert(0);
    }
  return it->second;
}

// const CandIdEcalMap& 
// EventManager::ecalMap( std::string str ) const
// {
//   map<string,CandIdEcalMap>::const_iterator it = _ecalMap.find( str );
//   if( it==_ecalMap.end() )
//     {
//       cout << str << endl;
//       assert(0);
//     }
//   return it->second;
// }

const CandIdRecHitMap& 
EventManager::ecalRecHitMap( const std::string& str ) const
{

  if(        str=="BasicCluster-EcalRecHit"   ) me()->loadEcalBCRecHitMaps();
  else if (  str=="PfBasicCluster-EcalRecHit" ) me()->loadEcalBCRecHitMaps();

  map<string,CandIdRecHitMap>::const_iterator it = _ecalRecHitMap.find( str );
  if( it==_ecalRecHitMap.end() )
    {
      cout << str << endl;
      assert(0);
    }
  return it->second;
}

const CandIdPSRecHitMap& 
EventManager::ecalPSRecHitMap( const std::string& str ) const
{

  if(  str=="PreshowerCluster-EcalPSRecHit" ){ me()->loadEcalPSCPSRecHitMaps();
   }
  map<string,CandIdPSRecHitMap>::const_iterator it = _ecalPSRecHitMap.find( str );
  if( it==_ecalPSRecHitMap.end() )
    {
      cout << str << endl;
      assert(0);
    }
  return it->second;
}

const CandMap& 
EventManager::candMap( const std::string& str ) const
{
  if(      str=="SuperCluster-BasicCluster" )     me()->loadEcalSCBCMaps();
  else if( str=="PfSuperCluster-PfBasicCluster" ) me()->loadEcalSCBCMaps();
  else if( str=="Jet-PfCand" ) me()->loadJetPfCandMaps();

  map<string,CandMap>::const_iterator it = _candMap.find( str );
  if( it==_candMap.end() )
    {
      cout << str << endl;
      assert(0);
    }
  return it->second;
}

const CandAssoc& 
EventManager::candAssoc( const std::string& str ) const
{
  if(      str=="electron-SuperCluster" )        me()->loadElectronSCMaps();
  else if( str=="electron-PfSuperCluster" )      me()->loadElectronSCMaps();
  else if( str=="electron-GsfTrack" )            me()->loadElectronSCMaps();
  else if( str=="SuperCluster-PfSuperCluster" )  me()->loadElectronSCMaps();
  else if( str=="PreshowerCluster-BasicCluster" ) me()->loadEcalBCPSCMaps();
  else if( str=="PreshowerCluster-EcalPSRecHit" ) me()->loadEcalPSCPSRecHitMaps();
  else if( str=="electron-Track" )               me()->loadElectronTrackMaps();
  else if( str=="electron-PfCand" )              me()->loadElectronPfCandMaps();
  else if( str=="muon-Track" )                   me()->loadMuonTrackMaps();
  else if( str=="muon-PfCand" )                  me()->loadMuonPfCandMaps();
  else if( str=="PatJet-PfJet")                  me()->loadJetRelatedMaps();
  else if( str=="PatJet-lepton")                 me()->loadJetRelatedMaps();

  map<string,CandAssoc>::const_iterator it = _candAssoc.find( str );
  if( it==_candAssoc.end() )
    {
      cout << str << endl;
      assert(0);
    }
  return it->second;
}

TH1*
EventManager::sumOfEt( string str, int jj ) const
{
  if( jj!=-1 )
    {
      assert( jj>=0 && jj<4 );
      TString dumstr("_");
      dumstr += jj; // hack
      str += dumstr.Data();
    }
  if( _sumOfEt.count( str )==0 )
    {
      string hname = "sumOfEt_";
      hname += str;
      if( jj==-1 )
	{
	  //*****
	  //	  cout << "Creating Tower histogram " << hname << endl;
	  TH1* h_;
	  h_ = new TH1F( hname.c_str(), hname.c_str(), 72, 0, 72 );
	  h_->SetDirectory(0);
	  (me()->_sumOfEt)[str] = h_;
	}
      else
	{
	  //	  for( int kk=0; kk<4; kk++ )
	  int kk = jj;
	    {
	      TString hname0( hname.c_str() );
	      hname0 += "_";
	      hname0 += kk;	  
	      //	      cout << hname0 << endl;
	      TH1* h_ = new TH1F( hname0, hname0, 41, 0, 41 );
	      h_->SetDirectory(0);
	      (me()->_sumOfEt)[str] = h_;
	    }
	}
    }
  TH1* h_ = (me()->_sumOfEt)[str];
  if( h_==0 )
    {
      cout << "str=" << str << endl;
      assert(0);
    }
  return h_;
}

const Candidate*
EventManager::decayTree() const
{
  return me()->decayTree();
}

Candidate*
EventManager::decayTreeOld()
{
  bool dbg_ = true;   
  Candidate* firstCand(0);

  CandList& _genEventCandidates = _candList[ kGenEventCandidate ];
  size_t nGenEvent = a().n("mc");
  if( dbg_ ) cout << "nMc=" << nGenEvent << endl;
  if( nGenEvent==0 ) return 0;

  map< int, Vertex* > _vertices;

  // loop over gen Event entries and create vertices
  for( size_t ii=0; ii<nGenEvent; ii++ )
    {
      assert( a().load( "gen", ii ) );

       int status = a().get_i("gen","status");

       //====
       float px_  = a().get_f("gen","px");
       float py_  = a().get_f("gen","py");     
       float pz_  = a().get_f("gen","pz");
       TVector3 tmpp3( px_, py_, pz_ );
       //=====
      
       if( dbg_ )
	 {
	   cout << "Gen " << ii << "-- \t";
	   cout << a().get_i("gen","pdgId") << "\t";
	   cout << a().get_f("gen","e") << "\t";
	   cout << tmpp3.Pt() << "\t";
	   cout << a().get_i("gen","status") << "\t";
	   cout << a().get_i("gen","PV") << "\t";
	   cout << a().get_i("gen","EV") << "\t";
	   cout << endl;
	 }
       if( status!=3 ) continue;

       
       // production vertex
       int ipv_ = a().get_i("gen","PV");
       assert( ipv_!=0 );
       Vertex* v3_(0);
       if( _vertices.count( ipv_ )==0 )
	 {
	   float pvx_ = a().get_f("gen","PV_x");
	   float pvy_ = a().get_f("gen","PV_y");
	   float pvz_ = a().get_f("gen","PV_z");
 	  //	  cout << "New production vertex at " << pvx_ << "/" << pvy_ << "/" << pvz_ << endl;
	   v3_ = Vertex::create( TVector3( pvx_, pvy_, pvz_ ) );
	   _vertices[ipv_] = v3_;
	 }

       // end vertex
       int iev_ = a().get_i("gen","EV");
       if( iev_!=0 &&_vertices.count( iev_ )==0 )
	 {
	   float evx_ = a().get_f("gen","EV_x");
	   float evy_ = a().get_f("gen","EV_y");
	   float evz_ = a().get_f("gen","EV_z");
	   //	  cout << "New end vertex at " << evx_ << "/" << evy_ << "/" << evz_ << endl;
	   v3_ = Vertex::create( TVector3( evx_, evy_, evz_ ) );
	   _vertices[iev_] = v3_;
	 }
    }

  Candidate* lastCand(0);
  
  for( size_t ii=0; ii<nGenEvent; ii++ )
    {
      assert( a().load( "gen", ii ) );
      
      int status = a().get_i("gen","status");
      //      if( status!=3 ) continue;
      
      float px_  = a().get_f("gen","px");
      float py_  = a().get_f("gen","py");     
      float pz_  = a().get_f("gen","pz");
      int pdgId_ = a().get_i("gen","pdgId");

      Vertex* pv3_(0);
      Vertex* ev3_(0);

      int ipv_ = a().get_i("gen","PV");
      assert( ipv_!=0 );
      if( _vertices.count( ipv_ )==0 ) continue;
      pv3_ = _vertices[ipv_];

      int iev_ = a().get_i("gen","EV");
      if( iev_!=0 )
	{
	  if( _vertices.count( iev_ )==1 ) ev3_ = _vertices[iev_];
	}
      
      const TParticlePDG* pdg_ = ParticleData::Table()->GetParticle( pdgId_ );
      assert( pdg_!=0 );

      //      int status_ = a().get_i("gen","status");
      TVector3 p3_( px_, py_, pz_ );
      Candidate* cand_ = Candidate::create( p3_, *pdg_, pv3_ );
      cand_->setMass( a().get_f("gen","m") );
      if( status==3 ) lastCand = cand_;

      if(dbg_) cout<<" LastCand pdgId "<<lastCand->pdgCode()<<endl;

      if( pv3_!=0 ) {
	pv3_->setOutgoingCand( cand_ ); }
      if( ev3_!=0 ) {
	ev3_->setIncomingCand( cand_ ); }
    }
  
  for( map< int, Vertex* >::const_iterator it=_vertices.begin();
       it!=_vertices.end(); ++it )
    {
      if( dbg_) cout << "\n********************\n Vertex " << it->first << endl;
      Vertex* vertex_ = it->second;
      if( dbg_) vertex_->print( cout );
      Candidate* incoming_ = vertex_->incomingCand();
      size_t nout_ = vertex_->nOutgoingCand();
      if( dbg_) cout << "Number of outgoing candidates: " << nout_ << endl;
      for( size_t ii=0; ii<nout_; ii++ )
	{
	  Candidate* outgoing_ = vertex_->outgoingCand( ii );
	  if( dbg_) cout << "Outgoing candidate " << ii << " --> " << outgoing_->name() << endl;
	}
      if( incoming_!=0 )
	{
	  if( dbg_)
	    {
	      cout << "Incoming candidate: ";;
	      cout << incoming_->name();
	      cout << endl;
	    }
	  for( size_t ii=0; ii<nout_; ii++ )
	    {
	      Candidate* outgoing_ = vertex_->outgoingCand( ii );
	      incoming_->addDaughter( outgoing_ );
	      if(dbg_ && outgoing_->theMother()==NULL ) cout<<" No Mother "<<endl;
	      if(dbg_ && outgoing_->theMother()!=NULL ) cout<<" Ya qque chose! "<<endl;
	      if(dbg_) cout << endl;
	    }
	}
    }

  assert( lastCand!=0 );
  //
  // take the last candidate
  //

  if( dbg_ ) 
    lastCand->print(cout);

  while( lastCand->theMother()!=0 )
    {
      firstCand = lastCand->theMother();
      lastCand = firstCand;
    }

  if( dbg_ ) cout << "genEventCandidate = " << _genEventCandidates.size() << endl;
  _genEventCandidates.push_back( firstCand );

  _vertices.clear();

  //  if(1) return firstCand;
  //
  // Now, final particles
  //
  dbg_ = false;
  if(dbg_) cout << "Final particles " << endl;
  for( size_t ii=0; ii<nGenEvent; ii++ )
    {
      assert( a().load( "gen", ii ) );

      int status = a().get_i("gen","status");
      if( status!=1 ) continue;
      int pdgId = a().get_i("gen","pdgId");

      // only create vertices with tau neutrinos
     
      if( abs(pdgId)!=16 ) continue;
     
      if( dbg_ )
	{
	  cout << "Gen " << ii << "-- \t";
	  cout << a().get_i("gen","pdgId") << "\t";
	  cout << a().get_f("gen","e") << "\t";
	  cout << a().get_i("gen","status") << "\t";
	  cout << a().get_i("gen","PV") << "\t";
	  cout << a().get_i("gen","EV") << "\t";
	  cout << endl;
	}

      // production vertex
      int ipv_ = a().get_i("gen","PV");
      assert( ipv_!=0 );
      Vertex* v3_(0);
      if( _vertices.count( ipv_ )==0 )
	{
	  float pvx_ = a().get_f("gen","PV_x");
	  float pvy_ = a().get_f("gen","PV_y");
	  float pvz_ = a().get_f("gen","PV_z");
	  v3_ = Vertex::create( TVector3( pvx_, pvy_, pvz_ ) );
	  _vertices[ipv_] = v3_;
	}

      // end vertex
      int iev_ = a().get_i("gen","EV");
      assert( iev_==0 );
    }

  for( size_t ii=0; ii<nGenEvent; ii++ )
    {
      assert( a().load( "gen", ii ) );

      int status = a().get_i("gen","status");
      if( status!=1 ) continue;

      float px_  = a().get_f("gen","px");
      float py_  = a().get_f("gen","py");     
      float pz_  = a().get_f("gen","pz");
      int pdgId_ = a().get_i("gen","pdgId");

      const TParticlePDG* pdg_ = ParticleData::Table()->GetParticle( pdgId_ );

      Vertex* pv3_(0);

      int ipv_ = a().get_i("gen","PV");
      assert( ipv_!=0 );

      // only keep particles at a tau vertex
      if( _vertices.count( ipv_ )==0 ) continue;

      pv3_ = _vertices[ipv_];

      int iev_ = a().get_i("gen","EV");
      assert( iev_==0 );

      //      int status_ = a().get_i("gen","status");
      TVector3 p3_( px_, py_, pz_ );
      Candidate* cand_ = Candidate::create( p3_, *pdg_, pv3_ );

      assert( pdg_!=0 );

      cand_->setMass( a().get_f("gen","m") );

      pv3_->setOutgoingCand( cand_ );
    }

  for( map< int, Vertex* >::const_iterator it=_vertices.begin();
       it!=_vertices.end(); ++it )
    {
      Vertex* vertex_ = it->second;
      if( vertex_->theMother()!=0 ) 
	{
	  //	  continue;
	}
      if( dbg_) vertex_->print(cout);
      size_t nout_ = vertex_->nOutgoingCand();
      if( dbg_) cout << "Number of outgoing candidates: " << nout_ << endl;
      bool keep = false;
      //      vector< const Candidate* > listOfDau_;
      CandList listOfDau_;
      for( size_t ii=0; ii<nout_; ii++ )
	{
	  Candidate* outgoing_ = vertex_->outgoingCand( ii );
	  int pdgId_ = outgoing_->pdgCode();
	  if( ( pdgId_==11 || pdgId_==-11 ) ||  
	      ( pdgId_==12 || pdgId_==-12 ) ||  
	      ( pdgId_==13 || pdgId_==-13 ) ||  
	      ( pdgId_==14 || pdgId_==-14 ) ||   
	      ( pdgId_==15 || pdgId_==-15 ) ||   
	      ( pdgId_==16 || pdgId_==-16 )
	      ) keep=true;
	  if( dbg_) cout << "Outgoing candidate " << ii << " --> " << outgoing_->name() << endl;
	  listOfDau_.push_back( outgoing_ );
	}
      if( !keep ) continue;

      Candidate* incoming_ = Candidate::create( listOfDau_ );
      //      Candidate* incoming_ = vertex_->incomingCand();
      
      //      if( incoming_==0 ) continue;

      if( dbg_ )
	{
	  cout << "incoming " << incoming_ << endl; 
	  incoming_->oneLine( cout );
	}

      _genEventCandidates.push_back( incoming_ );
    }

  return firstCand;
}

void
EventManager::defineCateg()
{  
  if( _categ!="" ) return;

  int nWp(0), nWm(0), nZ(0);
  int ne(0), nm(0), nt(0), nn(0);
  int nte(0), ntm(0), ntn(0);
  int nj(0);//, nc(0), 
  int nb(0);

  // a priori, an event is within acceptance
  _inAcceptance = true;

  // neutrinos
  _neutrinosFromBosonsP4.SetPxPyPzE(0,0,0,0);
  _neutrinosFromTausP4.SetPxPyPzE(0,0,0,0);

  //
  // the Generator-level decay tree
  //
  _genCand = decayTree();
  const CandList& genEventList_ = _candList[ kGenEventCandidate ];
  if( _genCand==0 ) 
    {
      _categ = "data";
      return;
    }
  _categ = "bkg";

  // 
  // determine the type of event
  //
  ostringstream os;
  CandList ewBosons;
  CandUtil::get(  24, _genCand, ewBosons );
  CandUtil::get( -24, _genCand, ewBosons );
  CandUtil::get(  23, _genCand, ewBosons );
  size_t nbos = ewBosons.size();

  size_t nlep(0);
  size_t ngam(0);
  size_t nlep2(0);
  size_t nqrk(0);
  CandList lepFromBosons;
  CandList qrkFromBosons;
  for( size_t ibos=0; ibos<nbos; ibos++ )
    {
      Candidate* bos_ = ewBosons[ibos];
      if( ibos>0 ) os << ":";
      os << bos_->name();
      int pdgId_ = bos_->pdgCode();

      // counters
      if( pdgId_==24 )       nWp++;
      else if( pdgId_==-24 ) nWm++;
      else if( pdgId_==23 )  nZ++;
      else assert(0);
      
      // search for leptons (11:e; 13:mu; 15:tau)
      for( int ii=11; ii<16; ii+=2 )
	{
	  CandUtil::get(  ii, bos_, lepFromBosons );
	  CandUtil::get( -ii, bos_, lepFromBosons );
	}

      // search for neutrinos (12:nu_e; 14:nu_mu; 16:nu_tau)
      for( int ii=12; ii<17; ii+=2 )
	{
	  CandUtil::get(  ii, bos_, lepFromBosons );
	  CandUtil::get( -ii, bos_, lepFromBosons );
	}

      // search for quarks (1:d, 2:u, 3:s, 4:c, 5:b)
      for( int ii=1; ii<6; ii++ )
	{
	  CandUtil::get(  ii, bos_, qrkFromBosons );
	  CandUtil::get( -ii, bos_, qrkFromBosons );
	}
    }
  nlep = lepFromBosons.size();
  nqrk = qrkFromBosons.size();
  if( nlep!=0 || nqrk!=0 )  os << "_";
  CandList gamFromBosons;
  bool  frstTkn(true);
  for( size_t ilep=0; ilep<nlep; ilep++ )
    {
      Candidate* lep_ = lepFromBosons[ilep];
      if( frstTkn ) frstTkn=false;  else  os << ":";  
      os << lep_->name();
      int pdgId_ = lep_->pdgCode();
      float eta_ = lep_->eta();
      float pt_  = lep_->pt();

     
      // counters
      if( pdgId_==11 || pdgId_==-11 )     
	{
	  ne++;
	  if( pt_<20. || fabs(eta_)>2.5 ) _inAcceptance = false; 
	}
      else if( pdgId_==13 || pdgId_==-13 )     
	{
	  nm++;
	  if( pt_<20. || fabs(eta_)>2.4 ) _inAcceptance = false; 
	}
      else if( pdgId_==15 || pdgId_==-15 )     nt++;
      else if( ( pdgId_==12 || pdgId_==-12 ) || ( pdgId_==14 || pdgId_==-14 ) )
	{	  
	  nn++;
	  // here, criteria on acceptance of MET?
	  _neutrinosFromBosonsP4 += lep_->p4();
	}
      else if(  pdgId_==16 || pdgId_==-16  ) nn++;
      else abort();
      
      // look for photons with pT>ptmin (22:gamma)
      float ptmin = 5.;
      CandUtil::get( 22, lep_, gamFromBosons, ptmin );
    }



  ngam = gamFromBosons.size();
  if( ngam>0 ) 
    {
      if( frstTkn ) frstTkn=false;  else  os << ":";  
      if( ngam>1 ) os << ngam;
      os << "gamma";	      
    }
  for( size_t iqrk=0; iqrk<nqrk; iqrk++ )
    {
      Candidate* qrk_ = qrkFromBosons[iqrk];
      if( frstTkn ) frstTkn=false;  else  os << ":";  
      os << qrk_->name();
      int pdgId_ = qrk_->pdgCode();
	  
      // counters
      if(        pdgId_==5 || pdgId_==-5 )     nb++;
      else if( ( pdgId_==4 || pdgId_==-4 ) ||
	       ( pdgId_==1 || pdgId_==-1 ) ||
	       ( pdgId_==2 || pdgId_==-2 ) ||
	       ( pdgId_==3 || pdgId_==-3 )   ) nj++;
      else abort();
    }

  if( nt>0 )
    {
      CandList lepFromTaus;
      for( size_t igen=1; igen<genEventList_.size(); igen++ )
	{
	  //	  assert( genEventList_[igen]!=0 );
	  //	  genEventList_[igen]->oneLine(cout);
	  Candidate* nt_  = CandUtil::first(  16, genEventList_[igen] ); 
	  Candidate* nbt_ = CandUtil::first( -16, genEventList_[igen] ); 
	  Candidate* ne_  = CandUtil::first(  12, genEventList_[igen] ); 
	  Candidate* nbe_ = CandUtil::first( -12, genEventList_[igen] ); 
	  Candidate* nm_  = CandUtil::first(  14, genEventList_[igen] ); 
	  Candidate* nbm_ = CandUtil::first( -14, genEventList_[igen] ); 
	  if( ( nt_!=0 && nbe_!=0 ) ||
	      ( nt_!=0 && nbm_!=0 ) ||
	      ( nbt_!=0 && ne_!=0 ) ||
	      ( nbt_!=0 && nm_!=0 )   )
	    {
	      CandUtil::get(  11, genEventList_[igen], lepFromTaus );
	      CandUtil::get( -11, genEventList_[igen], lepFromTaus );
	      CandUtil::get(  13, genEventList_[igen], lepFromTaus );
	      CandUtil::get( -13, genEventList_[igen], lepFromTaus );
	      CandUtil::get(  12, genEventList_[igen], lepFromTaus );
	      CandUtil::get( -12, genEventList_[igen], lepFromTaus );
	      CandUtil::get(  14, genEventList_[igen], lepFromTaus );
	      CandUtil::get( -14, genEventList_[igen], lepFromTaus );
	      CandUtil::get(  16, genEventList_[igen], lepFromTaus );
	      CandUtil::get( -16, genEventList_[igen], lepFromTaus );      
	    }
	}
      nlep2 = lepFromTaus.size();
      if( nlep2!=0 )
	{
	  os << "_";
	  for( size_t ilep2=0; ilep2<nlep2; ilep2++ )
	    {
	      Candidate* lep2_ = lepFromTaus[ilep2];
	      int pdgId_ = lep2_->pdgCode();
	      if( ilep2>0 ) os << ":";
	      os << lep2_->name();
	      
	      // counters
	      if(        pdgId_==11 || pdgId_==-11 )     nte++;
	      else if(   pdgId_==13 || pdgId_==-13 )     ntm++;
	      else if( ( pdgId_==12 || pdgId_==-12 ) ||
		       ( pdgId_==14 || pdgId_==-14 ) ||
		       ( pdgId_==16 || pdgId_==-16 )   )
		{
		  ntn++;
		  _neutrinosFromTausP4 += lep2_->p4();
		}
	      else abort();
	    }
	}
    }
  _categ_long = os.str();

  ostringstream os2;
  if( nWp>0 || nWm>0 || nZ>0 )
    {
      int nW = nWp+nWm;
      //      if( nWp==1 && nWm==1 )        os2 << "WW";
      //      else
      //	{
      //	  for( int ii=0; ii<nWp; ii++ ) os2 << "Wp";
      //	  for( int ii=0; ii<nWm; ii++ ) os2 << "Wm";
      //	}
      for( int ii=0; ii<nW; ii++ )  os2 << "W";
      for( int ii=0; ii<nZ;  ii++ ) os2 << "Z";
      if( ne>0 || nm>0 || nt>0 || nn>0 || nb>0 || nj>0 )
	{
	  os2 << "_";
	  if( ne>1 ) os2 << ne; if( ne>0 ) os2 << "e";
	  if( nm>1 ) os2 << nm; if( nm>0 ) os2 << "m";
	  if( nt>1 ) os2 << nt; if( nt>0 ) os2 << "t";
	  if( nn>1 ) os2 << nn; if( nn>0 ) os2 << "n";
	  if( nb>1 ) os2 << nb; if( nb>0 ) os2 << "b";
	  //	  if( nc>1 ) os2 << nc; if( nc>0 ) os2 << "c";
	  if( nj>1 ) os2 << nj; if( nj>0 ) os2 << "j";
	}
      if( nte>0 || ntm>0 )
	{
	  os2 << "_";
	  if( nte>1 ) os2 << nte; if( nte>0 ) os2 << "e";
	  if( ntm>1 ) os2 << ntm; if( ntm>0 ) os2 << "m";
	  if( ntn>1 ) os2 << ntn; if( ntn>0 ) os2 << "n";
	}
    }
  else 
    os2<< "bkg";
  _categ = os2.str();

  // debug special cases
  bool _verbose = false;
  if( _categ.find("WWZZ")!=string::npos ) _verbose=true;
  if( _categ == "ZZ" ) _verbose=true;

//   if( _categ == "ZZ_e2n" || _categ == "ZZ_2en" ) decayTree(true); 
//   if( _categ == "ZZ_m2n" || _categ == "ZZ_2mn" ) decayTree(true); 
//   if( _categ == "ZZ_t2n" || _categ == "ZZ_2tn" ) decayTree(true); 
//   if( _categ == "Z_2n" || _categ == "Z_2m" || _categ == "Z_2t" ) decayTree(true); 
//   if( _categ == "Z_2e2n" ) decayTree(true); 
//   if( _categ == "Z_2m2n" ) decayTree(true); 
//   if( _categ == "Z_2t2n" ) decayTree(true); 

}

void
EventManager::mcMatching()
{
  if( _mcMatchingMap.size()!=0 ) return;
  for( size_t ilep=0; ilep<2; ilep++ )
    {  
      const CandList& leptonList = candList( kElectron+ilep );  
      for( size_t ii=0; ii<leptonList.size(); ii++ )
	{
	  const Candidate& cand = *leptonList[ii];
	  ConeSelector mcsel_( cand, 0.05, cand.pdtEntry() );
	  vector< indexDRPair > p_;
	  mcsel_.getPairs( mcTruthCandidates(), p_ ); 
	  for( size_t jj=0; jj<p_.size(); jj++ )
	    {
	      _mcMatchingMap.insert( make_pair( cand.uid(), p_[jj] ) );
	    }
	}
    }
}

Vertex*
EventManager::vertexMatching(Vertex* vC) {


  Vertex* mVertex(0);

  float dX;
  float dY;
  float dZ;

  float dRtmp=10000;

  for(int unsigned iv=0;iv<_vertices.size();iv++) {

    dX =  vC->pos().X() - _vertices[iv]->pos().X();
    dY =  vC->pos().Y() - _vertices[iv]->pos().Y();
    dZ =  vC->pos().Z() - _vertices[iv]->pos().Z();

    float dR = sqrt(dX*dX + dY*dY + dZ*dZ);

    if(dR<dRtmp) {
      mVertex = _vertices[iv];
      dRtmp = dR;
    }

  }

  
  return mVertex;

}

const Candidate*
EventManager::mcCand( const Candidate& lepton ) const
{
 
  me()->mcMatching();
  CandIdDrMapIterator it_lower = _mcMatchingMap.lower_bound( lepton.uid() );
  CandIdDrMapIterator it_upper = _mcMatchingMap.upper_bound( lepton.uid() );
 
  if( it_lower==it_upper ) return 0;
  size_t uid = it_lower->second.first;
  return Candidate::base( uid );
}

string
EventManager::categ() const
{
  me()->defineCateg();
  return _categ;
}

string
EventManager::categ_long() const
{
  me()->defineCateg();
  return _categ_long;
}

bool 
EventManager::inAcceptance() const
{
  me()->defineCateg();
  return _inAcceptance; 
}

const TLorentzVector&
EventManager::neutrinosFromBosonsP4()
{
  me()->defineCateg();
  return _neutrinosFromBosonsP4; 
}

const TLorentzVector&
EventManager::neutrinosFromTausP4()
{
  me()->defineCateg();
  return _neutrinosFromTausP4; 
}

CandList&
EventManager::userCandList( const string& list )
{
  return _userCandList[ list ];
}

const CandList&
EventManager::compositeCandList( const string& list ) const
{
  me()->loadZllCandidates();
  map< string, CandList >::const_iterator it = _compositeCandList.find(list);
  assert( it!=_compositeCandList.end() );
  return it->second;
}

void
EventManager::setBoolInfo( CandInfo* info_, const string& prefix, const string& var, int v )
{
   info_->setBool( var,  (bool) a().get_b(prefix,var,v) );
}

void
EventManager::setIntInfo( CandInfo* info_, const string& prefix, const string& var, int v )
{
  info_->setInt( var,  a().get_i(prefix,var,v) );
}

void
EventManager::setFloatInfo( CandInfo* info_, const string& prefix, const string& var, int v )
{
  info_->setFloat( var,  a().get_f(prefix,var,v) );
}

bool
EventManager::isSelected()
{
  return true;
}


Candidate* 
EventManager::getFirst( const string& str, Candidate* c ) const
{
  const CandAssoc& candAssoc_ = candAssoc( str ); 
  return candAssoc_.getFirst( c->theBase() );
}

Candidate* 
EventManager::getSecond( const string& str, Candidate* c ) const
{
  const CandAssoc & candAssoc_ = candAssoc( str ); 
  return candAssoc_.getSecond( c->theBase() );
}

void 
CandAssoc::insert( Candidate* c1, Candidate* c2 )
{
  static int ifail1 = 0;
  //  static int ifail2 = 0;
  if( _map1.count(c1)!=0 )
    {
      ifail1++;
      //      if( ifail1++<10 )
      //	{
	  //	  cout << "ifail1=" << ifail1++ << endl;
	  //c1->oneLine( std::cout );
      //	}
      //  abort(); //MM FIXME
    }
  _map1[c1] = c2;
  if( _map2.count(c2)!=0 )
    {
      /* cout << "ifail2=" << ifail2++ << endl;
       std::cout << "WARNING---Multiple associations---" << std::endl;
       std::cout << "--- current candidate" << std::endl;    */  
      //c2->oneLine( std::cout );      
      // std::cout << "--- first candidate" << std::endl;      
      Candidate* c1_=_map2[c2];
      float dr1 = CandUtil::dR( c2, c1_ );
      //  std::cout << "dR=" << dr1 << "...";
      //   c1_->oneLine( std::cout ); 
      //   std::cout << "--- second candidate" << std::endl;      
      Candidate* c2_=c1;
      float dr2 = CandUtil::dR( c2, c2_ );
      //   std::cout << "dR=" << dr2 << "...";
      //   c2_->oneLine( std::cout );
       //              abort();
      if( dr2>dr1 ) return;
    }
  _map2[c2] = c1;
}

Candidate* 
CandAssoc::getSecond( Candidate* c ) const
{
  if( _map1.count(c)!=0 ) return _map1.find(c)->second;
  return 0;
}

Candidate* 
CandAssoc::getFirst( Candidate* c ) const
{
  if( _map2.count(c)!=0 ) return _map2.find(c)->second;
  return 0;
}

void 
CandAssoc::clear()
{
  _map1.clear();
  _map2.clear();
}

Candidate*
EventManager::createMcCand_( int ii )
{
  if( ii==-1 || !(a().load( "mc", ii )) ) return Candidate::create();
      
  int status = a().get_i("mc","status");
  int motherIndex = a().get_i("mc","motherIndex");
  int pdgId  = a().get_i("mc","pdgId");
  int nDau   = a().get_i("mc","nDau");
  int collId = a().get_i("mc","collId");
  //====
  float pt_   = a().get_f("mc","pt");
  float eta_  = a().get_f("mc","eta");     
  float phi_  = a().get_f("mc","phi");
  TVector3 tmpp3;
  tmpp3.SetPtEtaPhi(pt_,eta_,phi_);
  float pvx_ = a().get_f("mc","vx");
  float pvy_ = a().get_f("mc","vy");
  float pvz_ = a().get_f("mc","vz");
  //=====
  const TParticlePDG* pdg_ = ParticleData::Table()->GetParticle( pdgId );
  assert( pdg_!=0 );      
  Candidate* cand_ = Candidate::create( tmpp3, *pdg_ );
  cand_->setMass( a().get_f("mc","mass") );
  CandInfo* info_ = CandInfo::create();
  info_->setInt("status", status );
  info_->setInt("motherIndex", motherIndex );
  info_->setInt("pdgId",  pdgId  );
  info_->setInt("nDau",   nDau   );
  info_->setInt("collId", collId );
  info_->setFloat("pvx", pvx_ );
  info_->setFloat("pvy", pvy_ );
  info_->setFloat("pvz", pvz_ );
  cand_->setInfo( info_  );

  return cand_;
}


Candidate*
EventManager::decayTree( bool dbg_ )
{
  //  dbg_ = true;   

  CandList& _genEventCandidates = _candList[ kGenEventCandidate ];//_candList[ kMcTruthCandidate ];
  if( _genEventCandidates.size()>0 ) return _genEventCandidates[0];

  size_t nMcEvent = a().n("mc");
  if( dbg_ ) cout << "nMc=" << nMcEvent << endl;
  if( nMcEvent==0 ) return 0;

  //  map< int, Vertex* >    _vertices;
  map< int, Candidate* > _candidates;
  map< int, int > _status;
  map< int, int > _pdgId;
  map< int, vector<int> > _mothers;

  // loop over the mc store and create mothers
  _status[-1] = 3;
  _pdgId[-1] = 0;
  for( size_t ii=0; ii<nMcEvent; ii++ )
    {
      assert( a().load( "mc", ii ) );
      int status = a().get_i("mc","status");
      int motherIndex = a().get_i("mc","motherIndex");
      int pdgId       = a().get_i("mc","pdgId");
      if( motherIndex==-1 && status!=3 ) continue;
      _status[ii] = status;
      _pdgId[ii] = pdgId;
      //      if( _status.count(motherIndex) && _status[motherIndex]==3 )
      if( _status.count(motherIndex) )
	if( _status[motherIndex]==3 || abs(_pdgId[motherIndex])==15 )
	  _mothers[motherIndex].push_back(ii);
    }

  for( map< int, vector<int> >::const_iterator itm=_mothers.begin();
       itm!=_mothers.end(); ++itm )
    {
      int im_ =  itm->first;

      // for each mother, create a candidate
      if( _candidates.count(im_)==0 )
	{
	  _candidates[im_] = createMcCand_(im_);
	}
      Candidate*  mother_ = _candidates[im_];
      assert( mother_!=0 );

      CandInfo* info_ = mother_->info();
      if( info_!=0 )
	{
	  int status_ = info_->getInt("status");
	  int pdgId_  = info_->getInt("pdgId");
	  if( status_==3 || abs(pdgId_)==15 ) {}
	  else
	    continue;
	}

      Vertex* vertex_(0);
      for( size_t ii=0; ii<(itm->second).size(); ii++ )
	{
	  int id_ = (itm->second)[ii];

	  if( _candidates.count(id_)==0 )
	    {
	      _candidates[id_] = createMcCand_(id_);
	    }
	  Candidate* daughter_ = _candidates[id_];

	  assert( daughter_!=0 );
	  assert( daughter_->theMother()==0 );
	  if( vertex_==0 )
	    {
	      CandInfo* info_ = daughter_->info();
	      float pvx_ = info_->getFloat( "pvx" );
	      float pvy_ = info_->getFloat( "pvy" );
	      float pvz_ = info_->getFloat( "pvz" );
	      vertex_ = Vertex::create( TVector3( pvx_, pvy_, pvz_ ) );
	      vertex_->setIncomingCand( mother_ );
	    }
	  daughter_->setVertex( vertex_ );
	  vertex_->setOutgoingCand( daughter_ );
	  mother_->addDaughter( daughter_ );
	}
    }

  _genEventCandidates.push_back( _candidates[-1] );
  for( map< int, vector<int> >::const_iterator itm=_mothers.begin();
       itm!=_mothers.end(); ++itm )
    {
      int im_ =  itm->first;
      Candidate* mother_ = _candidates[im_];
      assert( mother_!=0 );

      if( _status[im_]!=3 ) continue;

      if( im_!=-1 )
	_genEventCandidates.push_back( mother_ );

      if( dbg_ )
	{
	  PrintTree prtTree;
	  cout << "Mother: " << im_ << endl;
	  cout << prtTree.Print( *mother_ );
	}
    }

  return _candidates[-1];
}

bool
EventManager::isFired( const std::string & hltLine, bool noVersion ) const  
{   
  return me()->a().isFired( hltLine, noVersion); 

}

void
EventManager::firedLines( map< int, string >& lines ) const
{
  me()->a().firedLines( lines );
}

float
EventManager::getRhoFastJet() {

  return a().getRFJ();

}

float
EventManager::getTrueNInt() {
  return a().getTrueNInt();
}

int
EventManager::getnInteract(string str) {
  return a().getnInteract(str);
}



Vertex* 
EventManager::findVertex( float eta_, float phi_) {
  
  Vertex* mVertex(0);

  const CandList& trkList = tracks();
  
  float dRtmp = 0.5;
  
  // for( size_t ii=0; ii<trkList.size(); ii++ )
  //   {
  //     assert( a().load( "trk", ii ) );
  //     Candidate* trkcand_= trkList[ii];
      
  //     //matching
  //     float dr = KineUtils::dR(eta_, trkcand_->eta(), phi_, trkcand_->phi() );
  //     if(dr < dRtmp ) {
  // 	dRtmp = dr;
  // 	mVertex = trkcand_->vertex();
  //     }
      
  //   }
  
  
  if(fabs(eta_)>3) return mVertex;

  map< Vertex* , float> sumPtAssoc;
  map< Vertex* , float>::const_iterator itmap;
  int cnt=0;
  for( size_t ii=0; ii<trkList.size(); ii++ )
    {
      assert( a().load( "trk", ii ) );
      Candidate* trkcand_= trkList[ii];
      
      //matching
      float dr = KineUtils::dR(eta_, trkcand_->eta(), phi_, trkcand_->phi() );
      if(dr < dRtmp ) {
	cnt++; 
	itmap = sumPtAssoc.find( trkcand_->vertex() );
	if( itmap != sumPtAssoc.end() ) {
	  sumPtAssoc[ trkcand_->vertex() ] += trkcand_->pt();
	  
	}
	else {
	  sumPtAssoc[ trkcand_->vertex() ] = trkcand_->pt();
	}
      }
    }

  float sPttmp=0;
  for(itmap = sumPtAssoc.begin(); itmap != sumPtAssoc.end(); itmap++) {
    if(sPttmp < (*itmap).second ) {
      sPttmp = (*itmap).second;
      mVertex =  (*itmap).first;
    }
  }
  
   return mVertex;

 //  string opt="basic";
//   Vertex* mVertex(0);

//   map< Vertex* , float> sumPtAssoc;
//   map< Vertex* , float>::const_iterator itmap;
//   float sumPt;

//   const CandMap& jetConst_ = candMap("Jet-PfCand");	 

//   CandMapIterator jc_lower_ = jetConst_.lower_bound( const_cast<Candidate*>(jet) );
//   CandMapIterator jc_upper_ = jetConst_.upper_bound( const_cast<Candidate*>(jet) );
//   CandMapIterator itjc;

//   for( itjc=jc_lower_; itjc!=jc_upper_; itjc++ )
//     { 
//        Candidate* cand = itjc->second;
//        if(cand->charge()==0) continue;

//        sumPt += cand->pt();
//        itmap = sumPtAssoc.find( cand->vertex() );
	
//        if( itmap != sumPtAssoc.end() ) {
// 	  sumPtAssoc[ cand->vertex() ] += cand->pt();
//        }
//        else {
// 	 sumPtAssoc[ cand->vertex() ] = cand->pt();
//        }
//     }

//   //now different options
//   if(opt=="basic") {
//     float sPttmp=0;
//     for(itmap = sumPtAssoc.begin(); itmap != sumPtAssoc.end(); itmap++) {
//       if(sPttmp < (*itmap).second ) {
// 	sPttmp = (*itmap).second;
// 	mVertex =  (*itmap).first;
//       }
//     }
//   }
//   if(opt=="atlas") {
    
//     for(itmap = sumPtAssoc.begin(); itmap != sumPtAssoc.end(); itmap++) {
//       if((*itmap).second/sumPt>0.75)
// 	{
// 	  mVertex = (*itmap).first;
// 	  break;
// 	}
//     }
//   }
  
//   return mVertex;

}


const Candidate*
EventManager::cleanPUmet( Vertex* _vtx, Candidate* bosC, float rF ) 
{
  
  TVector2 met(0,0);

  const CandList& jets = jetList( EventManager::kPatJet);
  
  bool findjet=false;
  for(int unsigned i=0;i<jets.size(); i++) {
    Candidate* jet = jets[i];
  
    if( _vtx == jet->vertex() && bosC->pt()>25 &&
       jet->pt()/bosC->pt() > 0.7 &&
       fabs(KineUtils::dPhi( jet->phi(), bosC->phi() ) ) > 2.8 )
      {findjet=true; break;}
  }



  bool assocJet=false;
  for(int unsigned i=0;i<jets.size(); i++) {
    
    Candidate* jet = jets[i];
    assocJet=false;
    
    bool jetbal = !findjet && _vtx != jet->vertex() 
                  && bosC->pt()>25 && jet->pt()/bosC->pt() > 0.7 
                  && fabs(KineUtils::dPhi( jet->phi(), bosC->phi() ) ) > 2.8;

    if( ( _vtx == jet->vertex() || jet->pt()>25 || jetbal ) )  
      {met -=  jet->p2(); assocJet=true;}
    
    if(bosC->isComposite() && !assocJet) {
      if( ( jet->dR( bosC->daughter(0) )<0.1 && jet->pt()/bosC->daughter(0)->pt() >0.8 ) ||
	  ( jet->dR( bosC->daughter(1) )<0.1 && jet->pt()/bosC->daughter(1)->pt() >0.8)  ) 
	{met -=  jet->p2(); assocJet=true;}
    }
    if(!bosC->isComposite() && !assocJet) {
      if( ( jet->dR( bosC) <0.1 && jet->pt()/bosC->pt() >0.8 )  ) 
	{met -=  jet->p2(); assocJet=true;}
    }
    
    if( !assocJet)
      {met -= jet->p2()/rF;}
    
  }
  
  CandList uncPfc_ = uncPfCandidates();
  bool assocpf=false;
  for(int unsigned i=0;i<uncPfc_.size();i++) {
    assocpf=false;
   const Candidate* pfc = uncPfc_[i];
    
    if( ( _vtx == pfc->vertex()  && pfc->charge()!=0 ) )
      { met -=  pfc->p2();assocpf=true;}
    if( !bosC->isComposite() && !assocpf) {
      if( ( pfc->dR( bosC) <0.1 && pfc->pt()/bosC->pt() >0.8 )  ) 
	{met -=  pfc->p2(); assocpf=true;}
    }

    if(!assocpf)
      { met -=  pfc->p2()/rF;}
  }

  Candidate* cmet = Candidate::create(met);

  return cmet;
  
}

Candidate*
EventManager::getCandidate( const string& str, size_t ii )
{
  int kl(-1);
  if( str=="electron" )
    {
      kl = kElectron;
    }
  else if( str=="muon" )
    {
      kl = kMuon;
    }
  else if( str=="photon" )
    {
      kl = kPhoton;
    }
  else if( str=="tracks" )
    {
      kl = kTrack;
    }
  size_t lsize(0);
  Candidate* cand(0);
  if( kl<0 )
    {
      const CandList& list = compositeCandList( str );
      lsize= list.size();
      if( ii<lsize ) cand=list[ii];
    }
  else 
    {
      const CandList& list = candList( kl );
      lsize= list.size();
      if( ii<lsize ) cand=list[ii];
    }
  cout << "Number of candidates in list for " << str << "=" << lsize << endl;
  if( cand!=0 ) cand->oneLine(cout);
  return cand;
}

void 
EventManager::print( string what )
{
  if( what=="all" || what=="MET" )
    {
      cout << "Calo MET" << endl;
      for( size_t ii=0; ii<EventServer::caloMET.size(); ii++ )
	{
	  cout << "caloMET -- " << EventServer::caloMET[ii] << endl;
	  met( kCaloMet, ii )->oneLine(cout);
	}

    //   cout << "Reco MET" << endl;
//       for( size_t ii=0; ii<EventServer::recoMET.size(); ii++ )
// 	{
// 	  cout << "recoMET -- " << EventServer::recoMET[ii] << endl;
// 	  met( kTrkMet,ii )->oneLine(cout);  
// 	}

      cout << "Pf MET" << endl;
      for( size_t ii=0; ii<EventServer::pfMET.size(); ii++ )
	{
	  cout << "pfMET -- " << EventServer::pfMET[ii] << endl;
	  met( kPfMet,ii )->oneLine(cout);
	}

      cout << "Pat MET" << endl;
      for( size_t ii=0; ii<EventServer::patMET.size(); ii++ )
	{
	  cout << "patMET -- " << EventServer::patMET[ii] << endl;
	  met( kPatMet,ii )->oneLine(cout);
	}

      CandList& cleanLeptons = userCandList("cleanLeptons");
      CandList& cleanJets    = userCandList("cleanJets");
      CandList unclusteredCharged;
      const CandList& pfList = pfCandidates();
      for( size_t ii=0; ii<pfList.size(); ii++ )
	{
	  Candidate* pfcand = pfList[ii];
	  if( !pfcand->isAtPrimaryVertex() ) continue;
	  bool isNeutral = pfcand->charge()==0;
	  if( isNeutral ) continue;
	  int pdgId = pfcand->pdgCode();
	  bool isHFHad         =     pdgId==1;
	  bool isHFEM          =     pdgId==2;
	  if( isHFHad || isHFEM )
	    {
	      continue;
	    }
	  else
	    {
	      CandInfo* info = pfcand->info();	  
	      if( info->getBool( "unClus") ) 
		{
		  unclusteredCharged.push_back( pfcand );		      
		}
	    }
	}

      CandList& all = userCandList("summedCandidates");
      Candidate* leptons = Candidate::create( cleanLeptons );
      leptons->setName("sumOfLeptonsAtPV");
      Candidate* jets    = Candidate::create( cleanJets );
      jets->setName("sumOfCleanedJets");
      Candidate* charged = Candidate::create( unclusteredCharged );
      charged->setName("sumOfUnclusteredChargedAtPV");
      all.push_back( leptons );
      all.push_back( jets );
      all.push_back( charged );
      Candidate* sum = Candidate::create( all );
      sum->setName("sumOfAllAtPV");
      all.push_back( sum );

      for( size_t ii=0; ii<all.size(); ii++ )
	{
	  Candidate* cand = all[ii];
	  cout << cand->name() << endl;
	  cand->oneLine(cout);
	}
      
      //      getBool("unClus");

    }
  
  if( what=="all" || what=="jets" )
    {
      CandUtil::printList( cout, electrons(), "electrons" ); 
      CandUtil::printList( cout, muons(),     "muons" ); 
      
      CandList& cleanJets = userCandList("cleanJets");
      const CandMap& jetConst_ = candMap("Jet-PfCand");	  
	  
      cout << "Number of clean jets at Primary vertex --> " << cleanJets.size() << endl;
      for( size_t ii=0; ii<cleanJets.size(); ii++ )
	{
	  Candidate* jet = cleanJets[ii];
	  CandInfo* infoJet = jet->info();

 	  cout << "\nPat jet: " << infoJet->getInt("index") << "\n\t";
	  jet->oneLine(cout);
	  CandMapIterator it_lower_ = jetConst_.lower_bound( jet );
	  CandMapIterator it_upper_ = jetConst_.upper_bound( jet );
	  CandMapIterator it;
	  for( it=it_lower_; it!=it_upper_; it++ )
	    { 
	      Candidate* cand = it->second;
	      cout << "\t\t";
	      cand->oneLine(cout);
	    }

// 	  if( pfJet!=0 )
// 	    {
// 	      cout << "PF-jet at DR=" << dR_pfJet << " \n\t";
// 	      pfJet->oneLine(cout);
// 	    }
// 	  if( lepton!=0 )
// 	    {
// 	      cout << "lepton at DR=" << dR_lJet << " \n\t";	      
// 	      lepton->oneLine(cout);
// 	    }
// 	  if( pfLepton!=0 ) 
// 	    {
// 	      cout << "associated Pf lepton \n\t";	      
// 	      pfLepton->oneLine(cout);		  
// 	    }		
	  	  
// 	  for( size_t ich=0; ich<chargedComponents.size(); ich++ )
// 	    {
// 	      Candidate* sum_= chargedComponents[ich];
// 	      int vtxId = sum_->vertexIndex();
// 	      cout << "Sum of charged constituents for vertex " << vtxId<< endl;
// 	      sum_->oneLine(cout);
// 	    }
// 	  if( neutralComponent!=0 )
// 	    {
// 	      cout << "Sum of neutral constituents " << endl;
// 	      neutralComponent->oneLine(cout);
// 	    }	  
// 	  if( isNeutral ) cout << "Jet is (mostly) neutral" << endl;
// 	  if( isForward ) cout << "Forward jet" << endl;
	}
      

//       CandUtil::printList( cout, jetList( kPfJet ),  "PF jets" ); 
//       cout << "\nPat jets above 20 GeV of Pt (detailed)" << endl;
//       for( size_t ii=0; ii<jetList(kPfJet).size(); ii++ )
//  	{
//  	  Candidate& cand = *jetList(kPfJet)[ii];
//  	  if( cand.pt() < 20. ) continue;
//  	  cand.print(cout);
//  	}

    }

  if( what=="all" || what=="tracks" )
    {
      cout << "\ntracks"  << endl;
      CandUtil::printList( cout, tracks() ); 
    }

  if( what=="all" || what=="PFCandidates" )
    {
      cout << "\npfCandidates"  << endl;
      CandUtil::printList( cout, candList( kPfCand ) ); 
    }

  if( what=="all" || what=="photons" )
    {
      cout << "\nphotons"  << endl;
      CandUtil::printList( cout, photons() ); 
    }

  //   for( size_t ii=0; ii<electrons().size(); ii++ )
  //     {
  //       Candidate* el1_ = electrons()[ii];
  //       for( size_t jj=ii+1; jj<electrons().size(); jj++ )
  // 	{
  // 	  Candidate* el2_ = electrons()[jj];

  // 	  Candidate* diel_ = Candidate::create( el1_, el2_ );
  // 	  diel_->oneLine(cout );
  // 	}
  //     }  

//   for( size_t ii=0; ii<jetList(kPfJet).size(); ii++ )
//     {
//       Candidate* jet1_ = jetList(kPfJet)[ii];
//       for( size_t jj=ii+1; jj<jetList(kPfJet).size(); jj++ )
// 	{
// 	  Candidate* jet2_ = jetList(kPfJet)[jj];

// 	  Candidate* dijet_ = Candidate::create( jet1_, jet2_ );
// 	  dijet_->oneLine(cout );
// 	}
//     }

  //   const CandList& scList = candList( kEcalSC );
  //   for( size_t ii=0; ii<scList.size(); ii++ )
  //     {
  //       Candidate* sc1_ = scList[ii];
  //       for( size_t jj=ii+1; jj<scList.size(); jj++ )
  // 	{
  // 	  Candidate* sc2_ = scList[jj];
	  
  // 	  Candidate* disc_ = Candidate::create( sc1_, sc2_ );
  // 	  disc_->oneLine(cout );
  // 	}
  //     }  

  if( what=="all" || what=="muons" )
    {  
      for( size_t ii=0; ii<muons().size(); ii++ )
	{
	  Candidate* mu_ = muons()[ii];
	  mu_->oneLine(cout);
	  
	  
	  Vertex* vtx_ = mu_->vertex();
	  vtx_->print(cout);

	}

//       if( muons().size()!=0 )
// 	{
// 	  Candidate* mu = muons()[0];
// 	  mu->oneLine(cout);
// 	  const CandList& pfList = pfCandidates();
// 	  cout << "PFCandidates ... " << pfList.size() << endl;
// 	  //      const CandAssoc& pf_assoc_ = candAssoc( "muon-PfCand" );
// 	  Candidate* pf = getSecond( "muon-PfCand", mu );
// 	  float drmin(0);
// 	  float ptmin(5.);
// 	  if( pf==0 )
// 	    {
// 	      cout << "PF muon not found, find the closest muon" << endl;
// 	      pf = CandUtil::closest( mu, pfList, drmin, ptmin, true ); 
// 	    }
// 	  if( pf==0 )
// 	    {
// 	      cout << "PF muon not found, find the closest hadron" << endl;
// 	      pf = CandUtil::closest( mu, pfList, drmin, ptmin, false, true ); 
// 	    }
// 	  if( pf!=0 )
// 	    {
// 	      cout << "at dR=" << drmin << "...";
// 	      pf->oneLine( cout );
// 	    }
// 	  else
// 	    cout << "No hard PF candidate with same charge in vicinity" << endl;
// 	}

//       for( size_t ii=0; ii<muons().size(); ii++ )
// 	{
// 	  Candidate* mu_ = muons()[ii];
// 	  Candidate* trk_ = getSecond( "muon-Track", mu_ );
// 	  if( trk_!=0 )
// 	    {
// 	      cout << "---->";
// 	      trk_->oneLine( cout );
// 	    }
// 	}  

//       cout << "Tracks" << endl;
//       const CandList& trackList = tracks();
//       //  const CandAssoc& trk_assoc_ = candAssoc( "muon-Track" );
//       //  const CandAssoc& gsf_assoc_ = candAssoc( "muon-GsfTrack" );
//       for( size_t ii=0; ii<trackList.size(); ii++ )
// 	{      
// 	  Candidate* trk_ = trackList[ii];
// 	  //      trk_->oneLine( cout );
// 	  Candidate* mu_ = getFirst( "muon-Track", trk_ );
// 	  if( mu_!=0 )
// 	    {
// 	      cout << "---->";
// 	      mu_->oneLine( cout );
// 	    }
// 	}

      cout << "\nmuons"  << endl;
      CandUtil::printList( cout, muons() ); 
    }

  if( what=="all" || what=="electrons" )
    {

      for( size_t ii=0; ii<electrons().size(); ii++ )
	{
	  Candidate* el_ = electrons()[ii];
	  //	  el_->oneLine(cout);
	  el_->print(cout);
	  
	  if( 1 ) continue;
	  
	  for( size_t isc_=0; isc_<2; isc_++ )
	    {
	      string sc_mapname_ = "electron-SuperCluster";
	      string bc_mapname_ = "SuperCluster-BasicCluster";
	      string rh_mapname_ = "BasicCluster-EcalRecHit";
	      string prefix_     = "E-Gamma";
	      if( isc_==1 )
		{
		  sc_mapname_ = "electron-PfSuperCluster";
		  bc_mapname_ = "PfSuperCluster-PfBasicCluster";
		  rh_mapname_ = "PfBasicCluster-EcalRecHit";
		  prefix_     = "Particle-Flow";
		}
	      
	      //	  const CandAssoc& sc_assoc_ = candAssoc( sc_mapname_ );
	      Candidate* sc_ = getSecond(  sc_mapname_, el_ );
	      const CandInfo* sc_info_ = sc_->info();
	      if( sc_==0 )
		{
		  cout << "Warning -- No " << prefix_ 
		       << " Super-Cluster associated: " << endl;
		  continue;
		}
	      cout << prefix_ << " Super-Cluster" << endl;
	      cout << "\t";
	      sc_->oneLine( cout );
	      
	      if( isc_== 0 )
		{
		  int ix_ = sc_info_->getInt( "ixClosestCell" );
		  int iy_ = sc_info_->getInt( "iyClosestCell" );
		  int iz_ = sc_info_->getInt( "izClosestCell" );
		  if( iz_==0 )
		    {
		      int ism_ = MEEBGeom::sm( ix_, iy_ );
		      MEEBGeom::XYCoord local_ = MEEBGeom::localCoord( ix_, iy_ );
		      cout << "Super-Cluster in the Ecal Barrel";
		      cout << "--- SM=" << ism_; 
		      cout << "[" << MEEBGeom::smName(ism_) << "]" ;
		      cout << " TT=" << MEEBGeom::tt( ix_, iy_ );
		      cout << "(" << local_.first << "," << local_.second << ")";
		      cout << endl;
		    }
		  else
		    {
		      // super-crystal
		      int iX_ = (ix_-1)/5+1;
		      int iY_ = (iy_-1)/5+1;
		      int ism_ = MEEEGeom::sm( iX_, iY_, iz_ );
		      int icr_ = MEEEGeom::crystal( ix_, iy_ );
		      
		      cout << "Super-Cluster in the Ecal Endcaps";
		      cout << "--- dee=" << MEEEGeom::dee( iX_, iY_, iz_ ); 
		      cout << "-S=" << MEEEGeom::sector( iX_, iY_ ); 
		      cout << "-SC=" << MEEEGeom::sc( iX_, iY_ ); 
		      cout << "[" << MEEEGeom::smName(ism_) << "]" ;
		      cout << " crystal=" << icr_;
		      cout << endl;
		    }
		}	  
	      cout << "\t... composed of the following basic clusters:" << endl;
	      //GHM const CandIdEcalMap& bc_map_ = ecalMap(bc_mapname_);
	      const CandMap& bc_map_ = candMap(bc_mapname_);
	      CandMapIterator bc_it_lower_ = bc_map_.lower_bound( sc_ );
	      CandMapIterator bc_it_upper_ = bc_map_.upper_bound( sc_ );
	      CandMapIterator bc_it;
	      for( bc_it=bc_it_lower_; bc_it!=bc_it_upper_; bc_it++ )
		{	    
		  Candidate* bc_ = bc_it->second;
		  cout << "\t\t";
		  bc_->oneLine( cout );
		  cout << "\t\t... composed of the following RecHits:" << endl;
		  const CandIdRecHitMap& rh_map_ = ecalRecHitMap(rh_mapname_);
		  CandIdRecHitMapIterator rh_it_lower_ = rh_map_.lower_bound( bc_->uid() );
		  CandIdRecHitMapIterator rh_it_upper_ = rh_map_.upper_bound( bc_->uid() );
		  CandIdRecHitMapIterator rh_it;
		  size_t nrh_(0);
		  for( rh_it=rh_it_lower_; rh_it!=rh_it_upper_; rh_it++ )
		    {	    
		      ecalRecHit rh_   = rh_it->second.first;
		      float      frac_ = rh_it->second.second;
		      float rh_e_ = rh_.e;
		      nrh_++;
		      if( rh_e_>0.5 )
			{
			  cout << "\t\t\t";
			  cout << "\tix/iy/iz=" 
			       << rh_.ix << "/"
			       << rh_.iy << "/"
			       << rh_.iz;
			  cout << "\te=" << rh_.e;
			  cout << "\tfrac=" << frac_;
			  cout << endl;
			}
		    } // rechits
		  if( nrh_>0 ) cout << "\t\t..." << nrh_ << " RecHits " << endl;
		} // basic clusters		  		  
	    } // e/g -> PF	  

	  // tracks
	  for( size_t itrk_=0; itrk_<2; itrk_++ ) {	  
	    string trk_mapname_ = "electron-Track";
	    if( itrk_==1 )
	      {
		trk_mapname_ = "electron-GsfTrack";
		cout << "GSF ";
	      }
	    cout << "track(s)" << endl;
	    //	const CandAssoc& trk_map_ = candAssoc( trk_mapname_ );
	    Candidate* trk_ = getSecond( trk_mapname_, el_ );
	    cout << "\t";
	    //	trk_->oneLine( cout );
	    trk_->print( cout );
	  } // track <-> GSF track      
	}// electrons


//       if( electrons().size()!=0 )
// 	{
// 	  Candidate* el = electrons()[0];
// 	  el->oneLine(cout);
// 	  const CandList& pfList = pfCandidates();
// 	  cout << "PFCandidates ... " << pfList.size() << endl;
// 	  //      const CandAssoc& pf_assoc_ = candAssoc( "electron-PfCand" );
// 	  Candidate* pf = getSecond( "electron-PfCand", el );
// 	  float drmin(0);
// 	  float ptmin(5.);
// 	  if( pf==0 )
// 	    {
// 	      cout << "PF electron not found, find the closest electron" << endl;
// 	      pf = CandUtil::closest( el, pfList, drmin, ptmin, true ); 
// 	    }
// 	  if( pf==0 )
// 	    {
// 	      cout << "PF electron not found, find the closest hadron" << endl;
// 	      pf = CandUtil::closest( el, pfList, drmin, ptmin, false, true ); 
// 	    }
// 	  if( pf!=0 )
// 	    {
// 	      cout << "at dR=" << drmin << "...";
// 	      pf->oneLine( cout );
// 	    }
// 	  else
// 	    cout << "No hard PF candidate with same charge in vicinity" << endl;
// 	}

      cout << "Electron tracks" << endl;
      const CandList& trackList = tracks();
      //  const CandAssoc& trk_assoc_ = candAssoc( "electron-Track" );
      //  const CandAssoc& gsf_assoc_ = candAssoc( "electron-GsfTrack" );
      for( size_t ii=0; ii<trackList.size(); ii++ )
	{      
	  Candidate* trk_ = trackList[ii];
	  //      trk_->oneLine( cout );
	  Candidate* el_ = getFirst( "electron-Track", trk_ );
	  if( el_!=0 )
	    {
	      cout << "---->";
	      el_->oneLine( cout );
	      Candidate* gsf_ = getSecond( "electron-GsfTrack", el_ );
	      assert( gsf_!=0 );
	      cout << "---->";
	      gsf_->oneLine( cout );
	    }
	}
  

      cout << "\nelectrons"  << endl;
      CandUtil::printList( cout, electrons() ); 
    }
  
  if( what=="all" || what=="SuperClusters" )
    {
      cout << "Super Clusters" << endl;
      string sc_mapname_ = "electron-SuperCluster";
      string bc_mapname_ = "SuperCluster-BasicCluster";
      string rh_mapname_ = "BasicCluster-EcalRecHit";
      string prefix_     = "E-Gamma";
      
      const CandList& ecalSCList = candList(kEcalSC);
      //  const CandAssoc& sc_assoc_ = candAssoc( "electron-SuperCluster" );
      for( size_t ii=0; ii<ecalSCList.size(); ii++ )
	{      
	  Candidate* sc_ = ecalSCList[ii];
	  sc_->oneLine( cout );
	  Candidate* el_ = getFirst( "electron-SuperCluster", sc_ );
	  if( el_!=0 )
	    {
	      cout << "---->";
	      el_->oneLine( cout );
	    }
	  const CandInfo* sc_info_ = sc_->info();
	  {
	    int ix_ = sc_info_->getInt( "ixClosestCell" );
	    int iy_ = sc_info_->getInt( "iyClosestCell" );
	    int iz_ = sc_info_->getInt( "izClosestCell" );
	    if( iz_==0 )
	      {
		int ism_ = MEEBGeom::sm( ix_, iy_ );
		MEEBGeom::XYCoord local_ = MEEBGeom::localCoord( ix_, iy_ );
		cout << "Super-Cluster in the Ecal Barrel";
		cout << "--- SM=" << ism_; 
		cout << "[" << MEEBGeom::smName(ism_) << "]" ;
		cout << " TT=" << MEEBGeom::tt( ix_, iy_ );
		cout << "(" << local_.first << "," << local_.second << ")";
		cout << endl;
	      }
	    else
	      {
		// super-crystal
		int iX_ = (ix_-1)/5+1;
		int iY_ = (iy_-1)/5+1;
		int ism_ = MEEEGeom::sm( iX_, iY_, iz_ );
		int icr_ = MEEEGeom::crystal( ix_, iy_ );
		
		cout << "Super-Cluster in the Ecal Endcaps";
		cout << "--- dee=" << MEEEGeom::dee( iX_, iY_, iz_ ); 
		cout << "-S=" << MEEEGeom::sector( iX_, iY_ ); 
		cout << "-SC=" << MEEEGeom::sc( iX_, iY_ ); 
		cout << "[" << MEEEGeom::smName(ism_) << "]" ;
		cout << " crystal=" << icr_;
		cout << endl;
	      }
	  }	  
	  cout << "\t... composed of the following basic clusters:" << endl;
	  //GHM const CandIdEcalMap& bc_map_ = ecalMap(bc_mapname_);
	  const CandMap& bc_map_ = candMap(bc_mapname_);
	  CandMapIterator bc_it_lower_ = bc_map_.lower_bound( sc_ );
	  CandMapIterator bc_it_upper_ = bc_map_.upper_bound( sc_ );
	  CandMapIterator bc_it;
	  for( bc_it=bc_it_lower_; bc_it!=bc_it_upper_; bc_it++ )
	    {	    
	      Candidate* bc_ = bc_it->second;
	      cout << "\t\t";
	      bc_->oneLine( cout );
	      cout << "\t\t... composed of the following RecHits:" << endl;
	      const CandIdRecHitMap& rh_map_ = ecalRecHitMap(rh_mapname_);
	      CandIdRecHitMapIterator rh_it_lower_ = rh_map_.lower_bound( bc_->uid() );
	      CandIdRecHitMapIterator rh_it_upper_ = rh_map_.upper_bound( bc_->uid() );
	      CandIdRecHitMapIterator rh_it;
	      size_t nrh_(0);
	      for( rh_it=rh_it_lower_; rh_it!=rh_it_upper_; rh_it++ )
		{	    
		  ecalRecHit rh_   = rh_it->second.first;
		  float      frac_ = rh_it->second.second;
		  float rh_e_ = rh_.e;
		  nrh_++;
		  if( rh_e_>0.5 )
		  {
		    cout << "\t\t\t";
		    cout << "\tix/iy/iz=" 
			 << rh_.ix << "/"
			 << rh_.iy << "/"
			 << rh_.iz;
		    cout << "\te=" << rh_.e;
		    cout << "\tfrac=" << frac_;
		    cout << endl;
		  }
		} // rechits
	      if( nrh_>0 ) cout << "\t\t..." << nrh_ << " RecHits " << endl;
	    } // basic clusters		  		  
	}
    }
  


  if( what=="all" || what=="Z" )
    {
      
      cout << "\nZee" << endl;
      CandUtil::printList( cout, compositeCandList( "Zee" ) );
      
      cout << "\nZmm" << endl;
      CandUtil::printList( cout, compositeCandList( "Zmm" ) );
      
    }
  

}

