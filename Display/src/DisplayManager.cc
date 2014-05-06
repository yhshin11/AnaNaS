#include "Display/src/DisplayManager.hh"

#include <cassert>
#include <algorithm>
#include <iostream>
#include <sstream>
using namespace std;

#include <TLatex.h>

#include "Analysis/utils/RooUtils.hh"
#include "Analysis/utils/KineUtils.hh"
#include "Analysis/utils/Constants.hh"
#include "Analysis/tools/CandUtil.hh"
#include "Analysis/selectors/ConeSelector.hh"
#include "Analysis/selectors/IsolationSelector.hh"

#include "MECore/src/MEEBGeom.hh"
#include "MECore/src/MEEEGeom.hh"
#include "MECore/src/MEEBDisplay.hh"
#include "MECore/src/MEEBDisplay.hh"


ClassImp(DisplayManager)

float DisplayManager::PtMax       = 100;
float DisplayManager::TrackPtMin   = 0.5;
float DisplayManager::TrackRMin   = 0;
float DisplayManager::TrackRMax   = 120;
float DisplayManager::TrackZMax   = 280;
float DisplayManager::LeptonPtMin  =  5;
float DisplayManager::CaloJetEtMin =  5;
float DisplayManager::PfJetEtMin   =  10;
//float DisplayManager::McPtcut   = 1;
//float DisplayManager::McEtacut  = 5;
//float DisplayManager::Wghcut    = 0.95;
//float DisplayManager::JetEtcut  = 10;
float DisplayManager::Emax      = 20;
float DisplayManager::length    = 25;
float DisplayManager::delta     = 1;
float DisplayManager::BField    = 3.8;
float DisplayManager::r0        = 0.;
float DisplayManager::x0        = 0.;
float DisplayManager::y0        = 0.;
float DisplayManager::z0        = 0.;
float DisplayManager::R         = 200;   // 170.;
float DisplayManager::L         = 350.;
float DisplayManager::AR        = 1.;
float DisplayManager::DE        = 5.3;
float DisplayManager::DPh       = 1.2;
float DisplayManager::eta0      = 0.;
float DisplayManager::phi0      = 1.;
float DisplayManager::CTScale   = 1.;

bool  DisplayManager::showPfCands     = true;
bool  DisplayManager::showEcalRH      = false;
bool  DisplayManager::showEcalSC      = false;
bool  DisplayManager::showCaloMET     = false;
bool  DisplayManager::showTrkMET      = false;
bool  DisplayManager::showCaloTowers  = true;
bool  DisplayManager::showCaloJets    = false;
bool  DisplayManager::showPfJets      = false;
bool  DisplayManager::showPfLeptons   = true;
bool  DisplayManager::showPfMet       = true;
bool  DisplayManager::showVertices    = true;
bool  DisplayManager::showTrackHits   = true;
bool  DisplayManager::showLeptonHits   = true;
bool  DisplayManager::showTracks      = true;
bool  DisplayManager::showElectrons   = true;
bool  DisplayManager::showMuons       = true;
bool  DisplayManager::showPhotons     = true;
bool  DisplayManager::showMcTruth     = false;
bool  DisplayManager::showPhotonIso   = false;
bool  DisplayManager::showElectronIso = true;
bool  DisplayManager::showMuonIso     = true;
bool  DisplayManager::showOnlyPVTracks    = false;
bool  DisplayManager::showOnlyPVTracks_vertex   = false;

float DisplayManager::electronDRin  = 0.3;
float DisplayManager::electronDRout = 0.4;
float DisplayManager::muonDRin      = 0.3;
float DisplayManager::muonDRout     = 0.4;
float DisplayManager::photonDRin    = 0.3;
float DisplayManager::photonDRout   = 0.4;

//int DisplayManager::trackColor = kGray+3;
//int DisplayManager::hitColor = kBlue;
int DisplayManager::lowPtTrackColor = kGray;
int DisplayManager::trackColor = kViolet;
//int DisplayManager::hitColor = kViolet;
int DisplayManager::CTHadColor = kGreen;
int DisplayManager::CTEmColor = kRed;

// colors
int DisplayManager::col_electron      = kRed;    // 
int DisplayManager::col_muon          = kMagenta; // 
//int DisplayManager::col_pvTrack_clu   = kCyan+1;   // kBlue;
//int DisplayManager::col_pvTrack_clu   = kGray;   // kBlue;
//int DisplayManager::col_pvTrack_unclu = kGray+1; // kGreen+1;
int DisplayManager::col_pvTrack_clu   = kAzure+10;   // kBlue;
int DisplayManager::col_pvTrack_unclu = kCyan-8; // kGreen+1;
int DisplayManager::col_puTrack_clu   = kBlue+2;   // kGreen+2;
int DisplayManager::col_puTrack_unclu = kBlue+3; // kGreen+3;


int DisplayManager::firstRun = 1;
int DisplayManager::lastRun = 999999999;
int DisplayManager::firstEvent = -1;
int DisplayManager::lastEvent = 999999999;

bool
DisplayManager::isSelected()
{
  //  int firstRun = 124120;
  //  int lastRun  = 124120;
  //  int firstEvt = 6613074;
  //  int lastEvt  = 6613074;

  if( run()<firstRun )   return false;
  if( run()>lastRun )    return false;
  if( event()<firstEvent ) return false;
  if( event()>lastEvent )  return false;

  //  cout << eventNumber() 
  //       << " --- run " << run() << " --- event " << event() << endl;

  if( 1 ) return 1; //!!! GHM

  cout << "BPTX=" << 
    ( isFired( "HLT_L1_BPTX" ) ? "true":"false" ) << endl;

  const Vertex& vtx_ = primaryVertex();

  cout << "primaryVertex=" << endl;
  vtx_.print(cout);
  cout << "OK " << endl;

  const CandInfo* info = vtx_.info();
  int vtx_nTracks(0);

  if( !info->getInt( "nTracks", vtx_nTracks ) ) return false;
  if( vtx_nTracks==0 ) return false;
  cout << "n[tracks] at vertex=" << vtx_nTracks << endl;

  //  size_t nTracks    = tracks().size();
  
  //   if(0)
  //     {
  //       if( nElectrons==0 ) return false;
  //     }
  
  //   if(0)
  //     {
  //       if( nMuons==0 ) return false;
  //     }
  
  //   if(0)
  //     {
  //       if( nPhotons==0 ) return false;
  //     }
  
  //   if(0)
  //     {
  //       if( nElectrons==0 && nMuons==0 && nPhotons==0 ) return false;
  //     }
  
  return true;
}

DisplayManager::DisplayManager( const char* filename, const char * collectionFileName ) 
  : EventManager( filename, collectionFileName ) 
    //    _view(hXY),
    //    _currentElectron(-1), _currentMuon(-1)
{ 
  setupDisplay();

  // setAtlas();

  setThreshold( 0., kGray );
  //  setThreshold( 1., kViolet );
}

DisplayManager::DisplayManager( std::vector< std::string > filenames, const char * collectionFileName ) 
  : EventManager( filenames, collectionFileName )
    //  , _view(hXY), 
    //    _currentElectron(-1), _currentMuon(-1)
{ 
  setupDisplay();

  // setAtlas();
  setThreshold( 0., kGray );
  //  setThreshold( 1., kViolet );
}

void
DisplayManager::setupDisplay()
{
  _ed = new MultiDisplay();
  _ed->EtaPhi();
  _ed->Pt();
  _ed->RZ("RZ");
  _ed->XY("XY");
}

void
DisplayManager::vertexViews()
{
  assert( _ed!=0 );
  _ed->setVertexLimits( _primaryVertex->pos().X(), 
			_primaryVertex->pos().Y(),
			_primaryVertex->pos().Z()  
			);
  _ed->RZ("VXZ");
  _ed->RZ("VYZ");
  //  _ed->RZ("VRZ");
  _ed->XY("VXY");
  _ed->drawScales();
  view();
}

void
DisplayManager::projectionViews()
{
  assert( _ed!=0 );
  _ed->RZ("XZ");
  _ed->RZ("YZ");
  view();
}


DisplayManager::~DisplayManager() 
{
  delete _ed;
}

bool
DisplayManager::goToEvent( int evt )
{
  EventManager::goToEvent( evt );
  //  _currentElectron=-1;
  //  _currentMuon=-1;
  view();
  return true;
}

bool
DisplayManager::nextEvent()
{
  do
    {
      bool ok_ = EventManager::nextEvent();
      if( !ok_ ) return false;
    }
  while( !isSelected() );

  _ed->refreshAllPanels();
  showOnlyPVTracks = false;
  view();
  return true;
}

bool
DisplayManager::view()
{
  _ed->setEvent( run(), event() );

  if( showCaloTowers )
    {
      drawCaloTowers();
    }

  if( showPfCands )
    {      
      drawPfCandidates();
      drawElectrons();
      drawMuons();
    }
   else 
     {
       if( showTracks )
 	{
 	  drawTracks();
 	}
      
       if( showElectrons )
 	{
 	  drawElectrons();
 	}
       
       if( showMuons )
	 {
 	  drawMuons();
 	}
     }

  drawVertices();

  ptView();

  etaPhiView();

  return true;
}

void
DisplayManager::drawElectrons()
{
  //  const CandAssoc& gsf_assoc_ = candAssoc( "electron-GsfTrack" );
  const CandList& listElectron = candList( kElectron );
  //      cout << "n[electrons]=" << listElectron.size() << endl;
  for( size_t iElectron=0; iElectron<listElectron.size(); iElectron++ )
    {
      Candidate* el_ = listElectron[iElectron];
      //      const Candidate& electron = *el_;
      //	  electron.oneLine(cout );
      //      Candidate* gsf_ = gsf_assoc_.getSecond( el_ );
      Candidate* gsf_ = getSecond( "electron-GsfTrack", el_ );
      assert( gsf_!=0 );
      //	  gsf_->oneLine(cout );
      //      int electronColor = kBlue;
      //      int electronColor = kMagenta;
      int electronColor = kRed;
      int electronWidth = 2;

      //      Vertex* vtx_ = el_->vertex();
      //      bool isAtPrimaryVtx = ( vtx_ == _primaryVertex );
      bool isAtPrimaryVtx = el_->isAtPrimaryVertex();
      if( !isAtPrimaryVtx )
	{
	  if( showOnlyPVTracks ) continue;
	  electronWidth = 1;
	}

      //      if( el_->pt()<20 ) 
      //      	{
      //	  electronColor = kGreen;
      //      	  electronWidth = 1;
      //      	}
      
      drawTrack( *gsf_, electronColor, electronWidth, kSolid );
      
//       if( showLeptonHits )
// 	{
// 	  size_t uid_ = gsf_->uid();  
// 	  CandIdIsoMapIterator it_lower 
// 	    = isoMap("gsfTrack-hits").lower_bound(uid_);
// 	  CandIdIsoMapIterator it_upper 
// 	    = isoMap("gsfTrack-hits").upper_bound(uid_);
// 	  CandIdIsoMapIterator it_;
// 	  map< float, pair<float,float> > map_;
// 	  //	    cout << "in DisplayManager uid=" << uid_ << endl;
// 	  for( it_=it_lower; it_!=it_upper; ++it_ )
// 	    {
// 	      isoDep iso_ = it_->second;
// 	      map_.insert( make_pair( iso_.dr, make_pair( iso_.eta, iso_.phi ) ) );
// 	      //		cout << iso_.dr << "/" << iso_.eta << "/" << iso_.phi << ":";
// 	    }      
// 	  //	    cout << endl;
// 	  float markerSize  = 0.4;
// 	  //GHM	    if( _view==hXY && R>150 ) markerSize=0.2;
// 	  //GHM  if( _view!=hXY && L>400 ) markerSize=0.2;
// 	  //	  int markerColor = kBlue;
// 	  int markerColor = kMagenta;
// 	  for( map< float, pair<float,float> >::const_iterator it_= map_.begin();
// 	       it_!=map_.end(); ++it_ )
// 	    {
// 	      _ed->drawMarkerRhoEtaPhi( it_->first, 
// 					it_->second.first, it_->second.second, 
// 					markerSize, markerColor );
// 	    }
// 	}
      
      
//       // 	  ConeSelector elSel( electron, 0.01 );
//       // 	  CandList tracks;
//       // 	  elSel.getSortedList( candList( kTrack ), tracks );
//       // 	  if( tracks.size()!=0 )
//       // 	    {
//       // 	      // 	  cout << gsfTracks.size() << "gsfTrack found !" << endl;
//       // 	      drawTrack( *tracks[0], kBlue, 1, kDashed );
//       // 	    }

//       if( _view==hEtaPhi )
// 	{
// 	  const CandInfo* info = electron.info();
// 	  float caloEta = info->getFloat("caloEta");
// 	  float caloPhi = info->getFloat("caloPhi");
// 	  float DRin_ecal  = 0.;
// 	  float DRout_ecal = 0.4;
// 	  _ed->isolation( caloEta, caloPhi, DRin_ecal, 0, caloEta, caloPhi, DRout_ecal, kSolid, 2, kBlue ); 
// 	}

//       if( !showEcalSC ) continue;

//       int uid_ = electron.uid();
//       CandIdIsoMapIterator it_lower = 
// 	isoMap("electron-EcalCluster").lower_bound(uid_);
//       CandIdIsoMapIterator it_upper = 
// 	isoMap("electron-EcalCluster").upper_bound(uid_);
//       CandIdIsoMapIterator it_;
//       for( it_=it_lower; it_!=it_upper; ++it_ )
// 	{
// 	  isoDep iso_ = it_->second;
// 	  float cl_rho_ = iso_.dr;
	      
// 	  float cl_e_   = iso_.val;
// 	  float cl_eta_ = iso_.eta;
// 	  float cl_phi_ = iso_.phi;

// 	  cout << "rho/e/eta/phi ";
// 	  cout << cl_rho_ << "/";
// 	  cout << cl_e_ << "/";
// 	  cout << cl_eta_ << "/";
// 	  cout << cl_phi_*Constants::radToDeg << endl;
	      
// 	  float cl_x_, cl_y_, cl_z_, cl_l_, cl_theta_;
// 	  RooUtils::toRTheta( cl_rho_, cl_eta_, cl_l_, cl_theta_ );
// 	  bool ok = RooUtils::toXYZ( cl_rho_, cl_eta_, cl_phi_, 
// 				     cl_x_, cl_y_, cl_z_ );
// 	  assert( ok ); 
	      
// 	  float cl_et_ = cl_e_*sin(cl_theta_);
	      
// 	  int icol = RooUtils::color( cl_et_, Emax );
	      
// 	  int lineWidth = 1;
// 	  int lineColor = icol;
// 	  int fillColor = icol; 
// 	  float dd(0), rin(0), rout(0), phi0(0);
// 	  float sign_ = 1; //!!!GHM!!!

// 	  if( _view==hXY )
// 	    {
// 	      dd = delta;
// 	      rin = cl_rho_;
// 	      rout = cl_rho_+ cl_et_*length/Emax;
// 	      phi0 = cl_phi_*Constants::radToDeg;
// 	    }
// 	  else if( _view==hRZ )
// 	    {
// 	      dd =  0.5;
// 	      phi0 = sign_*cl_theta_*Constants::radToDeg;;
// 	      rin = cl_l_;
// 	      rout = rin + cl_et_*length/Emax;
// 	    }
// 	  else if( _view==hXZ )
// 	    {
// 	      dd =  0.5;
// 	      phi0 = atan2(cl_x_,cl_z_)*Constants::radToDeg;;
// 	      rin = cl_x_/sin(phi0);
// 	      rout = rin + cl_et_*length/Emax;
// 	    }
// 	  else if( _view==hYZ )
// 	    {
// 	      dd   =  0.5;
// 	      phi0 = atan2(cl_y_,cl_z_)*Constants::radToDeg;;
// 	      rin  = cl_y_/sin(phi0);
// 	      rout = rin + cl_et_*length/Emax;
// 	    }
// 	  else
// 	    break;
	      
// 	  _ed->drawWedge( dd, rin, rout, phi0, lineWidth, lineColor, fillColor );    	      
// 	}      
    }
}

void
DisplayManager::drawMuons()
{
  const CandList& listMuon = candList( kMuon );
  //  const CandAssoc& trk_assoc_ = candAssoc( "muon-Track" );
  for( size_t iMuon=0; iMuon<listMuon.size(); iMuon++ )
    {
      Candidate* muon = listMuon[iMuon];
      int muonColor = col_muon;
      int muonWidth = 2;
	  
      //      Vertex* vtx_ = muon->vertex();
      //      bool isAtPrimaryVtx = ( vtx_ == _primaryVertex );
      bool isAtPrimaryVtx = muon->isAtPrimaryVertex();
      if( !isAtPrimaryVtx ) 
	{
	  if( showOnlyPVTracks ) continue;
	  muonWidth = 1;
	}

      //      if( muon->pt()<20 ) 
      //      	{
	  //	  muonColor = kViolet;
      //      	  muonWidth = 1;
      //      	}
      drawTrack( *muon, muonColor, muonWidth, kSolid ); 

      //      Candidate* trk_ = trk_assoc_.getSecond( muon );
      //      assert( trk_!=0 );
      //      drawTrack( *trk_, kRed, 1, kSolid );

//       if( showLeptonHits )
// 	{
// 	  size_t uid_ = muon->uid();  
// 	  CandIdIsoMapIterator it_lower 
// 	    = isoMap("mu-trackHits").lower_bound(uid_);
// 	  CandIdIsoMapIterator it_upper 
// 	    = isoMap("mu-trackHits").upper_bound(uid_);
// 	  CandIdIsoMapIterator it_;
// 	  map< float, pair<float,float> > map_;
// 	  //	    cout << "in DisplayManager uid=" << uid_ << endl;
// 	  for( it_=it_lower; it_!=it_upper; ++it_ )
// 	    {
// 	      isoDep iso_ = it_->second;
// 	      map_.insert( make_pair( iso_.dr, make_pair( iso_.eta, iso_.phi ) ) );
// 	    }      
// 	  float markerSize  = 0.4;
// 	  int markerColor = kRed;
// 	  for( map< float, pair<float,float> >::const_iterator it_= map_.begin();
// 	       it_!=map_.end(); ++it_ )
// 	    {
// 	      _ed->drawMarkerRhoEtaPhi( it_->first, 
// 					it_->second.first, it_->second.second, 
// 					markerSize, markerColor );
// 	    }
// 	}
      
    }
}

void
DisplayManager::drawJet( size_t ii, int color )
{
  const CandList& jList = jetList( kPatJet,-1 );
  if( ii>=jList.size() ) return;
  Candidate* jet    = jList[ii];

  //  cout << "Drawing jet ";
  //  jet->oneLine(cout);

  const CandMap& jetConst_ = candMap("Jet-PfCand");	  
  CandMapIterator it_lower_ = jetConst_.lower_bound( jet );
  CandMapIterator it_upper_ = jetConst_.upper_bound( jet );
  CandMapIterator it;
  for( it=it_lower_; it!=it_upper_; it++ )
    { 
      Candidate* cand = it->second;
      bool isNeutral = cand->charge()==0;
      if(  isNeutral ) continue;
      int pdgId = cand->pdgCode();
      bool isHFHad         =     pdgId==1;
      bool isHFEM          =     pdgId==2;
      if( isHFHad || isHFEM ) continue;
      if( !cand->isAtPrimaryVertex() ) continue;

      //      cand->oneLine(cout);
      drawTrack( *cand, color, 1, kSolid ); 
    }
}


void 
DisplayManager::drawCaloTower( const Candidate& ct )
{
  float et_ = ct.pt();
  const CandInfo* info = ct.info();
  // int ieta_ = info->getInt("ieta");
  // int iphi_ = info->getInt("iphi");

  //  float eta_ = ct.eta();
  //  float phi_ = ct.phi();
  //  float energy_ = info->getFloat("energy");
  //  float emEt_  = info->getFloat("emEt");
  //  float hadEt_ = info->getFloat("hadEt");
  int ieta_ = info->getInt("ieta");
  int iphi_ = info->getInt("iphi");
  //float emEt_  = info->getFloat("emEt");
  //float hadEt_ = info->getFloat("hadEt");
  //  if( hadEt_==0 ) return;  // draw only tower with hadronic energy
  //  float outerEt_ = info.getFloat("outerEt");

  //  int icol    = _ed->color_( et_, 0.1*Emax );
  
  int icol    = RooUtils::color( et_, Emax );
  //  cout << "et/Emax/icol " << et_ << "/" << Emax << "/" << icol << endl;

  //  _ed->drawCT( _view, ieta_, iphi_, 0, icol, icol ); 
  _ed->drawCT( ieta_, iphi_, 1, kBlack, icol ); 
}

void
DisplayManager::drawTrack( const Candidate& track,
			   int trajCol, int trajWidth, int trajStyle )
{
  size_t uid_ = track.uid();
  float pt_   = track.charge()*track.pt();
  //  if( fabs(pt_)<TrackPtMin ) return;
  int pdgCode_ = track.pdgCode();
  float eta_ = track.eta();
  float phi_ = track.phi();

  TVector3 v3_ = track.pos();
  float x0_  = 0;
  float y0_  = 0;
  float z0_  = 0;
  const CandInfo* info = track.info();
  if( !info->getFloat( "vx", x0_ ) ) x0_  = v3_.X();
  if( !info->getFloat( "vy", y0_ ) ) y0_  = v3_.Y();
  if( !info->getFloat( "vz", z0_ ) ) z0_  = v3_.Z();

  float rmin_ = TrackRMin;
  float rmax_ = TrackRMax;
  float zmax_ = TrackZMax;
  if( abs( pdgCode_ )==11 )
    {
      rmax_ = 130;
      zmax_ = 280;
      if( fabs(eta_)>1.5 )
	{
	  rmax_ = 150;
	  zmax_ = 325;
	}
    }
  if( abs( pdgCode_ )==13 )
    {
      rmax_ = 290;
      zmax_ = 820;
      //!!!      rmax_ = 167;
      //!!!      zmax_ = 350;
    }
  
  CandIdIsoMapIterator it_lower = isoMap("track-hits").lower_bound(uid_);
  CandIdIsoMapIterator it_upper = isoMap("track-hits").upper_bound(uid_);
  CandIdIsoMapIterator it_;
  map< float, pair<float,float> > map_;
  for( it_=it_lower; it_!=it_upper; ++it_ )
    {
      isoDep iso_ = it_->second;
      map_.insert( make_pair( iso_.dr, make_pair( iso_.eta, iso_.phi ) ) );
    }      
  if( map_.begin()!=map_.end() )
    {
      rmin_ = map_.begin()->first;
      rmax_ = map_.rbegin()->first;
      rmax_ += 0.05*(rmax_-rmin_);
    }
  
  // OK. plot.
  _ed->drawHelix( x0_, y0_, z0_, pt_, eta_, phi_, 
		  rmax_, zmax_,
		  BField,
		  trajCol, trajWidth, trajStyle );

  //   if( showTrackHits )
  //     {
  //       float markerSize  = 0.35;
  //       if( _view==hXY && R>150 ) markerSize=0.2;
  //       if( _view!=hXY && L>400 ) markerSize=0.2;
  //       int markerColor = kBlue;	
  //       for( map< float, pair<float,float> >::const_iterator it_= map_.begin();
  // 	   it_!=map_.end(); ++it_ )
  // 	{
  // 	  _ed->drawMarkerRhoEtaPhi( _view, 
  // 				    it_->first, 
  // 				    it_->second.first, it_->second.second, 
  // 				    markerSize, markerColor );
  // 	}
  //     }
}

// void
// DisplayManager::drawTrackHits( Candidate* trk_, float markerSize, int markerColor )
// {
//   //  const CandList& listTrack = candList( kTrack );
//   //  for( size_t iTrack=0; iTrack<listTrack.size(); iTrack++ )
//   //    {
//   //      Candidate* trk_ = listTrack[iTrack];
      
//   //      if( trk_->pt()<TrackPtMin ) continue;

//   size_t uid_ = trk_->uid();
  
//   CandIdIsoMapIterator it_lower = isoMap("track-hits").lower_bound(uid_);
//   CandIdIsoMapIterator it_upper = isoMap("track-hits").upper_bound(uid_);
//   CandIdIsoMapIterator it_;
//   map< float, pair<float,float> > map_;
//   for( it_=it_lower; it_!=it_upper; ++it_ )
//     {
//       isoDep iso_ = it_->second;
//       map_.insert( make_pair( iso_.dr, make_pair( iso_.eta, iso_.phi ) ) );
//     }      
//   for( map< float, pair<float,float> >::const_iterator it_= map_.begin();
//        it_!=map_.end(); ++it_ )
//     {
//       _ed->drawMarkerRhoEtaPhi( _view, 
// 				it_->first, 
// 				it_->second.first, it_->second.second, 
// 				markerSize, markerColor );
//     }
//   //    }
// }

// void
// DisplayManager::drawNextElectron( float scale, float scaleE )
// {
//   if( !drawElectron( ++_currentElectron, scale, scaleE ) )
//     {
//       _currentElectron=-1;
//       nextEvent();
//       drawNextElectron( scale, scaleE );
//     }
// }

// bool
// DisplayManager::redrawElectron( float scale, float scaleE )
// {
//   return drawElectron( _currentElectron, scale, scaleE );
// }

// bool
// DisplayManager::drawElectron( size_t ii, float scale, float scaleE )
// {
//   return drawLepton( kElectron, ii, scale, scaleE );
// }

// void
// DisplayManager::drawNextMuon( float scale, float scaleE )
// {
//   if( !drawMuon( ++_currentMuon, scale, scaleE ) )
//     {
//       _currentMuon=-1;
//       nextEvent();
//       drawNextMuon( scale, scaleE );
//     }
// }

// void
// DisplayManager::drawNextLepton( float scale, float scaleE )
// {
//   if( !drawMuon( ++_currentMuon, scale, scaleE ) )
//     {
//       if( !drawElectron( ++_currentElectron, scale, scaleE ) )
// 	{
// 	  _currentElectron=-1;
// 	  _currentMuon=-1;
// 	  nextEvent();
// 	  drawNextLepton( scale, scaleE );
// 	}
//     }
// }

// bool
// DisplayManager::redrawMuon( float scale, float scaleE )
// {
//   return drawMuon( _currentMuon, scale, scaleE );
// }

// bool
// DisplayManager::drawMuon( size_t ii, float scale, float scaleE )
// {
//   return drawLepton( kMuon, ii, scale, scaleE );
// }

// bool
// DisplayManager::drawLepton( int ilepton, size_t ii, float scale, float scaleE )
// {
//   bool isElectron(false);
//   bool isMuon(false);
//   if( ilepton==kElectron )
//     {
//       isElectron = true;
//     }
//   else if( ilepton==kMuon )
//     {
//       isMuon     = true;
//     }
//   else
//     abort();
//   const CandList& leptons_ = candList( ilepton );
//   if( ii>=leptons_.size() ) 
//     {
//       //      cout << "This lepton does not exist " << endl;
//       return false;
//     }
//   Candidate& cand = *( leptons_[ii] ); 
//   if( cand.pt()<LeptonPtMin ) return false;
//   cand.print(cout);

//   // make sure that we're dealing with leptons
//   // and that the settings are correct
//   int pdgCode_ = cand.pdgCode();
//   unsigned settings(0);

//   string mapname_;
//   if( isElectron )
//     {
//       assert( fabs( pdgCode_ )==11 );
//       settings = kElectron;
//       mapname_ = "electron-EcalIsoDep";      
//     }
//   else if( isMuon )
//     {
//       assert( fabs( pdgCode_ )==13 );
//       settings = kMuon;
//       mapname_ = "muon-EcalIsoDep";
//     }
//   else assert(0);

//   float pt  = cand.pt();
//   float chg = cand.charge();
//   assert( pt>0 );
//   float eta = cand.eta();
//   float phi = cand.phi();

//   // NEW!
//   // PF isolation
  
//   float dRmax = 0.4;
//   float dRmin = 0;
//   float pTmin = 0.5;
//   float Jhs   = -1; 
//   int pdg_[3] = { 211, 130, 22 };
//   const CandList& pfCandList = pfCandidates();
//   for( int ii=0; ii<3; ii++ )
//     {
//       ConeSelector coneSel_( cand, dRmax, dRmin, pTmin, Jhs, pdg_[ii] ); 
//       CandList list_;
//       coneSel_.getList( pfCandList, list_ );
//       CandUtil::printList( cout, list_ );
//     }

//   //
//   // to be checked
//   //
//   float r_ecal = 129;
//   float z_ecal = 280;
//   if( fabs(eta)>1.5 )
//     {
//       r_ecal = 150;
//       z_ecal = 325;
//     }
//   float r_hcal = 167;
//   float z_hcal = 350;

//   IsolationSelector isoSel_( settings );
//   isoSel_.computeIso( cand );
  
//   //  cand.oneLine( cout );
//   const CandInfo* info = cand.info();
//   float ESeed  = 5;
//   if( isElectron )
//     {
//       ESeed  = info->getFloat("ESeed");
//     }

//   int view_    = _view;
//   double eta0_ = eta0;
//   double phi0_ = phi0;
//   double DE_   = DE;
//   double DPh_  = DPh;
//   double Emax_ = Emax;

//   _view = hEtaPhi;
//   eta0 = eta;
//   phi0 = phi/Constants::pi;
//   DE   = 1/scale;
//   DPh  = DE*0.3;
//   Emax = ESeed/scaleE;

//   float DRin_trk = isoSel_.dR_trk_in;
//   float DRout_trk = isoSel_.dR_trk_out;
//   float Jhs_trk = isoSel_.strip_trk;

//   float DRin_hcal = isoSel_.dR_hcal_in;
//   float DRout_hcal = isoSel_.dR_hcal_out;
//   float Jhs_hcal = -1;

//   float etaAtEcal(0);
//   float phiAtEcal(0);
//   float etaAtHcal(0);
//   float phiAtHcal(0);
//   TVector3 v3 = cand.pos();
//   float x0  = v3.X();
//   float y0  = v3.Y();
//   float z0  = v3.Z();
//   KineUtils::straightToEnvelop( etaAtEcal, phiAtEcal, eta, phi, x0, y0, z0, r_ecal, z_ecal );
//   KineUtils::adjust( phiAtEcal, phi );
//   KineUtils::straightToEnvelop( etaAtHcal, phiAtHcal, eta, phi, x0, y0, z0, r_hcal, z_hcal );
//   KineUtils::adjust( phiAtHcal, phi );

//   int iecal = isoSel_.iecal;
//   float DRin_ecal  = isoSel_.dR_ecal_in[iecal];
//   float DRout_ecal = isoSel_.dR_ecal_out[iecal];
//   float Jhs_ecal   = isoSel_.strip_ecal[iecal];
  
//   //  view();
//   _ed->showDet = 1;
//   _ed->etaPhiPlane( DE, DPh, eta0, phi0, Emax );
//   //  drawCaloTowers();

//   drawEcalRecHits();
//   //  if( isElectron )
//   //    {
//   //      int index_ = info->getInt("index");
//   //GHM      drawElectronHits( index_ );
//   //    }
//   //  drawEcalSuperClusters();

//   float iso    = 0;
//   float S      = 0;
//   float S_trk  = 0;
//   float S_ecal = 0;
//   float S_hcal = 0;
//   float a_trk  = isoSel_.a_trk;
//   float a_hcal = isoSel_.a_hcal;
//   float a_ecal = isoSel_.a_ecal;
//   int n_trk(0);
//   int n_hcal(0);
//   int n_ecal(0);
  
//   const CandList& listTracks = candList( kTrack );
//   for( size_t ii=0; ii<listTracks.size(); ii++ )
//     {
//       const Candidate* trk = listTracks[ii];
//       float pt_  = trk->pt();
//       float eta_ = trk->eta();
//       float phi_ = trk->phi();
//       float deta_ = eta_-eta;
//       float dR_=KineUtils::dR( eta, eta_, phi, phi_ );
//       bool inCone_=true;
//       if( dR_<DRin_trk )  inCone_=false;
//       if( dR_>DRout_trk ) inCone_=false;
//       if( fabs( deta_ )<Jhs_trk ) inCone_=false;

//       //      cout << "dR/in/out " << dR_ << "/" << DRin_trk << "/" << DRout_trk << endl;
//       TVector3 v3_ = trk->pos();
//       float x0_  = v3_.X();
//       float y0_  = v3_.Y();
//       float z0_  = v3_.Z();
//       float r0_ = sqrt( x0_*x0_ + y0_*y0_ );

//       bool keep_=inCone_;
//       if( pt_<isoSel_.pt0_trk )        keep_=false; 
//       if( r0_>isoSel_.d0_trk )         keep_=false;
//       if( fabs(z0_)>isoSel_.dz0_trk )  keep_=false;

//       if( keep_ ) 
// 	{
// 	  n_trk++;
// 	  S_trk += pt_;
// 	}

//       //
//       // drawing
//       //
//       int icol_;
//       if( inCone_ )
// 	{
// 	  icol_ = ( keep_ ) ? kGreen : kMagenta;
// 	}
//       else 
// 	{
// 	  icol_ = kRed;
// 	}
//       _ed->drawMarker( eta_, phi_/Constants::pi, 5, 2, icol_ );
//     }


//   float eta0_hcal = eta;
//   float phi0_hcal = phi;
//   KineUtils::helixToEnvelop( eta0_hcal, phi0_hcal, chg*pt, eta, phi, x0, y0, z0, r_hcal, z_hcal, BField );

//   const CandList& listCT = candList( kCaloTower );
//   for( size_t iCT=0; iCT<listCT.size(); iCT++ )
//     {
//       const Candidate& ct = *listCT[iCT];
//       const CandInfo* info = ct.info();
//       float eta_ = ct.eta();
//       float phi_ = ct.phi();
//       float e_ = info->getFloat("energy");
//       float emEt_  = info->getFloat("emEt");
//       float hadEt_ = info->getFloat("hadEt");
      
//       _ed->drawCaloTower( eta_, phi_/Constants::pi, e_, emEt_, hadEt_, Emax );

//       if( hadEt_==0 ) continue;

//       float hadE_ = KineUtils::e( hadEt_, eta_ );

//       float dR_out=KineUtils::dR( etaAtHcal, eta_, phiAtHcal, phi_ );
//       float dR_in =KineUtils::dR( eta0_hcal, eta_, phi0_hcal, phi_ );
//       bool inCone_=true;
//       if( dR_out>DRout_hcal ) inCone_=false;
//       if( dR_in <DRin_hcal  ) inCone_=false;

//       bool keep_=inCone_;
//       if( hadEt_<isoSel_.et0_hcal ) keep_=false; 
//       if( hadE_<isoSel_.e0_hcal )   keep_=false; 

//       if( keep_ ) 
// 	{
// 	  n_hcal++;
// 	  S_hcal += hadEt_;
// 	}

//       int icol_;
//       if( inCone_ )
// 	{
// 	  icol_ = ( keep_ ) ? kGreen : kMagenta;
// 	}
//       else 
// 	{
// 	  icol_ = kRed;
// 	}
//       _ed->drawMarker( eta_, phi_/Constants::pi, 3, 4, icol_ );
//     }

//   _ed->showDet = 3;
//   _ed->drawEtaPhiDet();
//   _ed->showDet = 1;


//   float eta0_ecal = etaAtEcal;
//   float phi0_ecal = phiAtEcal;
//   //  cout << eta0_hcal << "/" << phi0_hcal << endl;
			     
//   float caloEta(0);
//   float caloPhi(0);
//   float trackEta(0);
//   float trackPhi(0);
//   if( isElectron )
//     {
//       caloEta = info->getFloat("caloEta");
//       caloPhi = info->getFloat("caloPhi");
//       KineUtils::adjust( caloPhi, phi );
//       trackEta = info->getFloat("trackEta");
//       trackPhi = info->getFloat("trackPhi");      
//       KineUtils::adjust( trackPhi, phi );
//       eta0_ecal = caloEta;
//       phi0_ecal = caloPhi;
//     }
//   else
//     {
//       KineUtils::helixToEnvelop( eta0_ecal, phi0_ecal, chg*pt, eta, phi, x0, y0, z0, r_ecal, z_ecal, BField );
//     }

//   // get the map of iso object for this candidate
//   const CandIdIsoMap& map_ = isoMap(mapname_);
//   CandIdIsoMapIterator it_lower_ = map_.lower_bound( cand.uid() );
//   CandIdIsoMapIterator it_upper_ = map_.upper_bound( cand.uid() );
//   CandIdIsoMapIterator it;
	
//   bool dbg_ = true;
//   if( dbg_ )
//     {
//       cout << "Uid=" << cand.uid() << endl;
//       cout << "Number of associated ECAL deposits: " << map_.count( cand.uid() ) << endl;             
//     }
//   for( it=it_lower_; it!=it_upper_; it++ )
//     {
//       //      float dr_  = (*it).second.dr;
//       float e_   = (*it).second.val; // WARNING!
//       float eta_ = (*it).second.eta;
//       float phi_ = (*it).second.phi;
//       float et_  = KineUtils::et( e_, eta_ );
//       //      float th_ = 2.0*atan(exp(-eta_));
//       //      assert( th_>0 );
//       //      float et_ = e_*sin(th_);

//       float deta_ = eta_-eta0_ecal;
//       float dR_out  = KineUtils::dR( etaAtEcal, eta_, phiAtEcal, phi_ );
//       float dR_in   = KineUtils::dR( eta0_ecal, eta_, phi0_ecal, phi_ );
//       //       //      cout << "deta/dphi " << deta_ << "/" << dphi_ << endl;
//       //       if( dphi_ >  Constants::pi ) dphi_-=2*Constants::pi; 
//       //       if( dphi_ < -Constants::pi ) dphi_+=2*Constants::pi;       
//       //      float dR_ = sqrt( deta_*deta_ + dphi_*dphi_ );
//       //      cout << "dphi/dr " << dphi_ << "/" << dR_ << endl;
      
//       bool inCone_=true;
//       if( dR_out>DRout_ecal ) inCone_=false;
//       if( dR_in <DRin_ecal  ) inCone_=false;
//       if( fabs( deta_ )<Jhs_ecal ) inCone_=false;

//       bool keep_=inCone_;
//       int iecal_ = ( fabs( eta_ )<=1.49 ) ? 0 : 1;
//       if(  e_<isoSel_.e0_ecal[iecal_]  ) keep_=false;
//       if( et_<isoSel_.et0_ecal         ) keep_=false;

//       if( keep_ ) 
// 	{
// 	  n_ecal++;
// 	  S_ecal += et_;
// 	}

//       int icol_;
//       if( inCone_ )
// 	{
// 	  icol_ = ( keep_ ) ? kGreen : kMagenta;
// 	}
//       else 
// 	{
// 	  icol_ = kRed;
// 	}
	
//       _ed->drawMarker( eta_, phi_/Constants::pi, 2, 1, icol_ );

//       //       //	    float phi_ = (*it).second.phi;
//     }
  
//   //  cout << "track : eta/phi " << trackEta << "/" << trackPhi/Constants::pi << endl;

//   if( isElectron )
//     {
//       _ed->drawMarker( trackEta, trackPhi/Constants::pi, 30, 2, 
// 		       //		       kBlue );
// 		       kMagenta );
//       //      cout << "calo : eta/phi " << caloEta << "/" << caloPhi/Constants::pi << endl;
//       _ed->drawMarker( caloEta, caloPhi/Constants::pi, 30, 2, kBlack );
//     }

//   cand.oneLine(cout);
//   S = a_trk*S_trk + a_ecal*S_ecal + a_hcal*S_hcal;
//   iso = pt/( pt + S );
//   cout << "S_trk  = "   << S_trk   << " /" << n_trk << endl;
//   cout << "S_hcal = "   << S_hcal  << " /" << n_hcal << endl;
//   cout << "S_ecal = "   << S_ecal  << " /" << n_ecal << endl;
//   cout << "S      = "   << S       << endl;
//   cout << "iso=" << iso << "/" << isoSel_.iso << endl;
  
//   //  drawTracks( kViolet, 1 );
//   drawTracks();
//   drawTrack( cand, kBlack, 4, kSolid );

//   int linestyle_, linewidth_, linecolor_;

//   linestyle_ = kSolid;
//   linewidth_ = 1;
//   //  linecolor_ = kBlue;
//   linecolor_ = kMagenta;
//   _ed->isolation( eta, phi, DRin_trk, Jhs_trk, eta, phi, DRout_trk, linestyle_, linewidth_, linecolor_ );

//   linestyle_ = kDashed;
//   linewidth_ = 1;
//   linecolor_ = kRed;
//   _ed->isolation( eta0_hcal, phi0_hcal, DRin_hcal, Jhs_hcal, etaAtHcal, phiAtHcal, DRout_hcal, linestyle_, linewidth_, linecolor_ );

//   linestyle_ = kSolid;
//   linewidth_ = 2;
//   linecolor_ = kBlack;
//   //  if( isElectron )
//   //    {
//   _ed->isolation( eta0_ecal, phi0_ecal, DRin_ecal, Jhs_ecal, etaAtEcal, phiAtEcal, DRout_ecal, linestyle_, linewidth_, linecolor_ ); 
//   //    }
//   //  else
//   //    {
//   //      _ed->isolation( etaAtEcal, phiAtEcal, DRin_ecal, DRout_ecal, Jhs_ecal, linestyle_, linewidth_, linecolor_ );
//   //    }
      
//   _view=view_;
//   eta0 = eta0_;
//   phi0 = phi0_;
//   DE = DE_;
//   DPh = DPh_;
//   Emax = Emax_;

//   return true;
// }

void
DisplayManager::drawEcalRecHits( float Emax_, bool applyThresh )
{
  size_t nEcalRH = a().n("rh");
  for( size_t ii=0; ii<nEcalRH; ii++ )
    {
      assert( a().load( "rh", ii ) );
      float e_    = a().get_f("rh","energy");
      int ix = a().get_i("rh","ix");
      int iy = a().get_i("rh","iy");
      int iz = a().get_i("rh","iz");
      
      //      float Emax=30;  // !!!
      float emin=0;
      float etmin=0;
      if( applyThresh )
	{
	  etmin=0.200;
	  if( iz==0 )
	    {
	      emin = 0.120;
	    }
	  else
	    {
	      emin = 0.450;
	    }
	}
      _ed->drawEcalRecHit( ix, iy, iz, e_, Emax_, emin, etmin );
    }      
}

// void
// DisplayManager::drawGlobalEcalHist( string opt )
// {
//   _ed->refreshGlobalEcalHist( Emax );

//   size_t nEcalRH = a().n("rh");
//   for( size_t ii=0; ii<nEcalRH; ii++ )
//     {
//       assert( a().load( "rh", ii ) );
//       float e_    = a().get_f("rh","energy");
//       int ix = a().get_i("rh","ix");
//       int iy = a().get_i("rh","iy");
//       int iz = a().get_i("rh","iz");
      
//       float emin=0;
//       //      float etmin=0.200;
//       if( iz==0 )
// 	{
// 	  emin = 0.120;
// 	}
//       else
// 	{
// 	  emin = 0.450;
// 	}

//       if( e_>emin )
// 	_ed->fillGlobalEcalHist( ix, iy, iz, e_ );
//     }      

//   _ed->drawGlobalEcalHist( opt );
// }

// void
// DisplayManager::drawElectronEcalHist( int hist_size, string opt )
// {
//   size_t ilepton = kElectron; //!!!
//   size_t ii;
//   bool isElectron(false);
//   bool isMuon(false);
//   if( ilepton==kElectron )
//     {
//       isElectron = true;
//       ii = _currentElectron;
//     }
//   else if( ilepton==kMuon )
//     {
//       isMuon     = true;
//       ii = _currentMuon;
//     }
//   else
//     abort();
//   const CandList& leptons_ = candList( ilepton );
//   if( ii>=leptons_.size() ) 
//     {
//       //      cout << "This lepton does not exist " << endl;
//       return;
//     }
//   Candidate& cand = *( leptons_[ii] ); 
//   cand.print(cout);

//   float pt  = cand.pt();
//   assert( pt>0 );
//   float eta = cand.eta();
//   float phi = cand.phi();

//   size_t nEcalRH = a().n("rh");
//   vector<int>  v_ix;
//   vector<int>  v_iy;
//   vector<int>  v_iz;
//   vector<float> v_e;
//   for( size_t ii=0; ii<nEcalRH; ii++ )
//     {
//       assert( a().load( "rh", ii ) );
//       float e_    = a().get_f("rh","energy");
//       int ix = a().get_i("rh","ix");
//       int iy = a().get_i("rh","iy");
//       int iz = a().get_i("rh","iz");
      
//       float emin=-1000;
//       //       if( iz==0 )
//       // 	{
//       // 	  emin = 0.120;
//       // 	}
//       //       else
//       // 	{
//       // 	  emin = 0.450;
//       // 	}

//       if( e_>emin )
// 	{
// 	  v_ix.push_back( ix );
// 	  v_iy.push_back( iy );
// 	  v_iz.push_back( iz );
// 	  v_e.push_back( e_ );
// 	}
//     }      

//   _ed->drawElectronEcalHist( pt, eta, phi, v_ix, v_iy, v_iz, v_e, hist_size, opt );
// }

// void
// DisplayManager::drawEcalSuperClusters()
// {
//   size_t nEcalSC = a().n("sc");
//   for( size_t ii=0; ii<nEcalSC; ii++ )
//     {
//       assert( a().load( "sc", ii ) );
//       //      float esc_  = a().get_f("sc","energy");
//       size_t nHits_ = (size_t) a().get_i("sc","nHits");
//       //      cout << "Ecal Super Cluster: E=" << esc_ << "GeV";
//       //      cout << " nHits=" << nHits_ << endl;	      	      
//       for( size_t jj=0; jj<nHits_; jj++ )
// 	{		  
// 	  int ilab_ = a().get_vi("sc","ilab",jj);
// 	  float e_  = a().get_vf("sc","ehit",jj);
	  
// 	  int iz = ilab_%10-1;
// 	  int ix = ilab_/10000;
// 	  if( iz==0 )
// 	    {
// 	      ix = ilab_/10000-85;
// 	    }
// 	  int iy = (ilab_%10000)/10;

// 	  //	  cout << "ilab/ix/iy/iz/e " << ilab_ << "/" << ix << "/" << iy << "/" << iz << "/" << e_ << endl;
	  
// 	  // 		  int icol;
// 	  // 		  if( e_<=0 )
// 	  // 		    {
// 	  // 		      icol = kBlack;
// 	  // 		    }
// 	  // 		  else
// 	  // 		    {
// 	  // 		      icol = RooUtils::color( e_, Emax );
// 	  // 		    }
	  
// 	  _ed->drawEcalRecHit( ix, iy, iz, e_, Emax , -1000, -1000);
	  
// 	}
//     }
// }

// void
// DisplayManager::drawElectronHits( int ii )
// {
//   string prefix_("el");
  
//   assert( a().load( prefix_.c_str(), ii ) );
//   size_t nHits_ = (size_t) a().get_i(prefix_.c_str(),"nHits");
//   //      cout << "Ecal Super Cluster: E=" << esc_ << "GeV";
//   //      cout << " nHits=" << nHits_ << endl;	      	      
//   for( size_t jj=0; jj<nHits_; jj++ )
//     {		  
//       int ilab_ = a().get_vi(prefix_.c_str(),"ilab",jj);
//       float e_  = a().get_vf(prefix_.c_str(),"ehit",jj);
      
//       int iz = ilab_%10-1;
//       int ix = ilab_/10000;
//       if( iz==0 )
// 	{
// 	  ix = ilab_/10000-85;
// 	}
//       int iy = (ilab_%10000)/10;
      
//       _ed->drawEcalRecHit( ix, iy, iz, e_, Emax , -1000, -1000);
      
//     }
// }

void
DisplayManager::drawPrimaryVertex()
{
  const Vertex& vtx_ = primaryVertex();
  drawVertex( vtx_, kRed, true );
}

void 
DisplayManager::drawVertices()
{
  if( !showOnlyPVTracks_vertex )
    {
      const VtxList& listVtx = vertices();
      for( size_t iVtx=0; iVtx<listVtx.size(); iVtx++ )
	{
	  Vertex* vtx_ = listVtx[iVtx];
	  drawVertex( *vtx_, kGray+1, false );
	}  
    }
  drawPrimaryVertex();
}

void
DisplayManager::drawVertex( const Vertex& vtx_, int markerColor, bool isPrimary )
{
  //  const Vertex& vtx_ = primaryVertex();
  const TVector3&  pos_ = vtx_.pos();
  
  float X_ = pos_.X();
  float Y_ = pos_.Y();
  float Z_ = pos_.Z();

  _ed->drawVertex( X_, Y_, Z_, markerColor, isPrimary );
}

void
DisplayManager::drawTracks()
{ 
  const CandList& listTrack = candList( kTrack );
  //      CandUtil::printList( cout, listTrack );
  //  const CandAssoc& trk_assoc_ = candAssoc( "electron-Track" );

  for( size_t jj=0; jj<2; jj++ )
    {
      if( showOnlyPVTracks && jj==0 ) continue;
      //      for( size_t iTrack=0; iTrack<listTrack.size(); iTrack++ )
      for( int iTrack=listTrack.size()-1; iTrack>=0; iTrack-- )
	{
	  Candidate* trk_ = listTrack[iTrack];
	  float pt_ = trk_->pt();
	  
	  //	  Vertex* vtx_ = trk_->vertex();
	  //	  bool isAtPrimaryVtx = ( vtx_ == _primaryVertex );
	  bool isAtPrimaryVtx = trk_->isAtPrimaryVertex();
	  if(  isAtPrimaryVtx && jj==0 ) continue;
	  if( !isAtPrimaryVtx && jj==1 ) continue;
	  
	  //      if( !isAtPrimaryVtx ) continue;
	  
	  if( pt_<_thres[0].pt ) continue;
	  int color_ = _thres[0].color;
	  int width_ = _thres[0].width;
	  int style_ = _thres[0].style;
	  
	  for( size_t ii=1; ii<_thres.size(); ii++ )
	    {
	      if( pt_<_thres[ii].pt ) break;
	      color_ = _thres[ii].color;
	      width_ = _thres[ii].width;
	      style_ = _thres[ii].style;
	    }
	  
	  if( !isAtPrimaryVtx )
	    {
	      color_ = kGreen+3;
	      width_ = 1;
	    }
	  
	  Candidate* el_ = getFirst( "electron-Track", trk_ );
	  if( el_!=0 ) 
	    {
	      // 	      cout << "found electron track " << endl;
	      // 	      cout << "<--- ";
	      // 	      trk_->oneLine( cout );
	      // 	      cout << "---> ";
	      // 	      el_->oneLine( cout );	  
	      continue;
	    }
	  
	  Candidate* mu_ = getFirst( "muon-Track", trk_ );
	  if( mu_!=0 ) 
	    {
	      //	      cout << "found muon track " << endl;
	      // 	      cout << "<--- ";
	      // 	      trk_->oneLine( cout );
	      // 	      cout << "---> ";
	      // 	      mu_->oneLine( cout );	  
	      continue;
	    }
	  
	  drawTrack( *trk_, color_, width_, style_ );
	  
// 	  if( showTrackHits )
// 	    {
// 	      float markerSize  = 0.2;
// 	      //	      if( _view==hXY && R>150 ) markerSize=0.2;
// 	      //	      if( _view!=hXY && L>400 ) markerSize=0.2;
// 	      int markerColor = color_;
// 	      drawTrackHits( trk_, markerSize, markerColor );
// 	    }
	}
    }
}

void
DisplayManager::drawPfCandidates()
{ 
  loadJetPfCandMaps();

  const CandList& pfList = candList( kPfCand );

  for( int jj=6; jj>=0; jj-- )
    {
      if( jj==2 )
	{
	  CandList& cleanJets = userCandList("cleanJets");
	  size_t nJets = cleanJets.size();
	  for( int ijet=nJets-1; ijet>=0; ijet-- )
	    {
	      int color = col_pvTrack_unclu;
	      if( ijet==3 ) color = kViolet;
	      if( ijet==2 ) color = kGreen;
	      if( ijet==1 ) color = kYellow;
	      if( ijet==0 ) color = col_pvTrack_clu;
	      int jetVtxIdx = cleanJets[ijet]->info()->getInt("index");
	      drawJet(jetVtxIdx,color);
	    }
	  continue;
	}

      for( size_t ii=0; ii<pfList.size(); ii++ )
	{
	  Candidate* cand_ = pfList[ii];
	  int pdgId = cand_->pdgCode();
	  if( cand_->pt()<TrackPtMin ) continue;
	  if( cand_->charge()==0 ) continue;
	  if( pdgId==1 || pdgId==2 ) continue;

	  CandInfo* info_ = cand_->info();
	  bool isUnclus        = info_->getBool( "unClus");
	  bool isAtPrimaryVtx  = cand_->isAtPrimaryVertex();
	  bool isElectron      =     abs(pdgId)==11;
	  bool isMuon          =     abs(pdgId)==13;
	  bool isLepton        = isElectron || isMuon;
	  //bool isPhoton        =     pdgId==22;
	  //bool isNeutralHadron =     pdgId==130;
	  //bool isChargedHadron =     abs(pdgId)==211;
	  bool isHFHad         =     pdgId==1;
	  bool isHFEM          =     pdgId==2;
	  
	  if( isHFHad || isHFEM ) continue;
	  
	  int width = 1;
	  int color = col_puTrack_unclu;
	  if( isElectron ) 
	    {
	      color = col_electron;
	      if( isAtPrimaryVtx ) 
		continue; //!!! draw primary electrons separately
	    }
	  else if( isMuon )
	    {
	      color = col_muon;
	      if( isAtPrimaryVtx ) 
		continue; //!!! draw primary muons separately
	    } 
	  if( !isAtPrimaryVtx )
	    {
	      if( showOnlyPVTracks ) continue;
	      if( !isLepton )
		{
		  if( isUnclus )
		    {
		      if( jj!=6 ) continue;
		      color = col_puTrack_unclu ;
		    }
		  else
		    {
		      if( jj!=5 ) continue;
		      color = col_puTrack_clu;
		    }
		}
	      else 
		{
		  if( jj!=4 ) continue;
		}
	    }
	  else
	    {
	      if( !isLepton )
		{
		  if( isUnclus )
		    {
		      if( jj!=3 ) continue;
		      color = col_pvTrack_unclu;
		    }
		  else
		    {
		      continue;
		      //		      if( jj!=2 ) continue;
		      //		      color = col_pvTrack_clu;		      
		    }
		}
	      else
		{
		  width = 2;
		  if( isElectron )
		    { 
		      if( jj!=1 ) continue;
		    }
		  else if( isMuon )
		    {
		      if( jj!=0 ) continue;
		    }
		  else
		    continue;
		}
	    }
	  
	  drawTrack( *cand_, color, width );

	}
    }

  
  //  loadJetPfCandMaps();

  if(0) 
    {
      const CandMap& jet_cand = candMap("Jet-PFCand");
      
      for( size_t ii=0; ii<jetList(kPatJet).size(); ii++ )
	{
	  Candidate* jet_ = jetList(kPatJet)[ii];
	  if( jet_->pt()<30 ) continue;
	  jet_->oneLine( cout );
	  int nTracks = 0;
	  int nTracksAtPrimaryVertex = 0;
	  CandMapIterator c_it_lower_ = jet_cand.lower_bound( jet_ );
	  CandMapIterator c_it_upper_ = jet_cand.upper_bound( jet_ );
	  CandMapIterator c_it;
	  for( c_it=c_it_lower_; c_it!=c_it_upper_; c_it++ )
	    {	    
	      Candidate* c_ = c_it->second;
	      c_->oneLine(cout);
	      if( fabs(c_->pdgCode())!=211  ) continue;
	      nTracks++;
	      //		  else continue;
	      if( c_->isAtPrimaryVertex() )
		nTracksAtPrimaryVertex++;
	    }
	  if( nTracksAtPrimaryVertex>0 )
	    {
	      //	      cout << "Jet with " << nTracksAtPrimaryVertex << "/" << nTracks << " at primary vertex " << endl;
	      //	      jet_->print(cout);
	    }
	}
    }
}

void
DisplayManager::drawCaloTowers()
{
  loadCaloTowers();
  
  vector< TH1* > _ct_h(10);
  
  size_t ii=0;
  _ct_h[ii++] = sumOfEt("ct_pt");
  _ct_h[ii++] = sumOfEt("ct_em");
  for( int kk=0; kk<4; kk++ )
    {
      _ct_h[ii++] = sumOfEt( "ct_et_eta", kk );
      _ct_h[ii++] = sumOfEt( "ct_em_eta", kk );
    }
  _ed->drawCaloTowerHistos( _ct_h );


  if( _ed->goToPanel("EtaPhi") )
    {
      const CandList& listCT = candList( kCaloTower );
      for( size_t iCT=0; iCT<listCT.size(); iCT++ )
	{
	  const Candidate& ct = *listCT[iCT];
	  drawCaloTower( ct );
	}
    }
}


//    if( _view==hXY )
//     {
//       _ed->drawPhiHisto( ct_pt_h, CTScale, 1, CTHadColor, CTHadColor, false, true );
//       _ed->drawPhiHisto( ct_em_h, CTScale, 1, CTEmColor, CTEmColor, false, true );
//     }
//   else if( _view==hRZ || _view==hXZ || _view==hYZ || _view==hPhiZ )
//     {
// 	{
// 	  TH1* et_h_ = sumOfEt( "ct_et_eta", kk );
// 	  _ed->drawEtaHisto( et_h_, kk, CTScale, 1, CTHadColor, CTHadColor, false, true );
// 	  TH1* em_h_ = sumOfEt( "ct_em_eta", kk );
// 	  _ed->drawEtaHisto( em_h_, kk, CTScale, 1, CTEmColor, CTEmColor, false, true );
// 	}
//     }
//   else
//     {
//       const CandList& listCT = candList( kCaloTower );
//       for( size_t iCT=0; iCT<listCT.size(); iCT++ )
// 	{
// 	  const Candidate& ct = *listCT[iCT];
// 	  drawCaloTower( ct );
// 	}
//     }
//}

// void
// DisplayManager::drawCaloJets()
// {
//   const CandList& list_j = jetList( kCaloJet );
//   for( size_t ij=0; ij<list_j.size(); ij++ )
//     {
//       const Candidate& jCand = *list_j[ij];

//       CandInfo* info_ = jCand.info();
//       float area = info_->getFloat("towersArea");
//       float emf  = info_->getFloat("emFraction");

//       float etj  = jCand.pt();

      
//       if( etj<CaloJetEtMin ) continue;

//       float etaj = jCand.eta();
//       float phij = jCand.phi();
      
//       _ed->drawCaloJet( _view, etj, etaj, phij, area, emf, CTScale, PtMax );
//     }
// }

// void
// DisplayManager::drawCaloMet()
// {
//   const Candidate* met = caloMet();
//   float et = met->pt();
//   float phi = met->phi();
//   _ed->drawMET( et, phi, CTScale, Emax);
// }    

// void
// DisplayManager::drawTrkMet()
// {
//   const Candidate* tmet = met( kPatMet );
//   float et = tmet->pt();
//   float phi = tmet->phi();
//   _ed->drawMET( et, phi, CTScale, Emax);
// }   


// void
// DisplayManager::printGenMet()
// {
 
//   const Candidate* tmet = met( kGenMet );
//   float et = tmet->pt();
//   float phi = tmet->phi();
  
//   cout<<" genMet   "<<" pt="<<et<<"   phi="<<phi<<endl; 
//   //  _ed->drawMET( _view, et, phi, CTScale, Emax);
// }  

void
DisplayManager::ptView()
{  
  //  if( ptmax_>0 ) PtMax = ptmax_;

  TString curPanel_ = _ed->curPanel;
  _ed->goToPanel("Pt");
  
  float scale=1;
  float dphi=0.3;
  int lineWidth=2;
  int fillStyle=3001;
  int color1=kWhite;
  int color2=kWhite;
  bool fillSameColor=true;

  draw2DVectors( scale, dphi, lineWidth, fillStyle, color1, color2, fillSameColor );


  _ed->goToPanel( curPanel_ );

}

void
DisplayManager::etaPhiView()
{
  TString curPanel_ = _ed->curPanel;
  if( !_ed->goToPanel("EtaPhi") ) return;

  float scale_ = 1./200.;
  
  const float etMin_ = PfJetEtMin;
  const CandList& list_j = jetList( kPfJet );
  for( size_t ii=0; ii<list_j.size(); ii++ )
    {
      Candidate& jet_ = *list_j[ii];
      float eta_ = jet_.eta();
      float phi_ = jet_.phi();
      float et_  = jet_.pt();
      if( et_<etMin_ ) continue;
      //      jet_.print(cout);
      const CandInfo* info = jet_.info();
      float emFraction = 
	info->getFloat("chargedEmEnergyFraction")+
	info->getFloat("neutralEmEnergyFraction");
      // int ieta_ = info->getInt("ieta");
      // int iphi_ = info->getInt("iphi");
      
      //  float eta_ = ct.eta();
      //  float phi_ = ct.phi();
      //  float energy_ = info->getFloat("energy");
      //  float emEt_  = info->getFloat("emEt");
      //  float hadEt_ = info->getFloat("hadEt");
      //      int ieta_ = info->getInt("ieta");
      //      int iphi_ = info->getInt("iphi");
      TVector2 jetP2_ = jet_.p2();
      //      int color1 = kOrange;
      //      int color2 = kYellow;
      //      if( fillSameColor ) color2 = color1; 
      _ed->isolation( eta_, phi_, 0., 0., eta_, phi_, scale_*et_*emFraction, kDashed, 1, kOrange );
      _ed->isolation( eta_, phi_, 0., 0., eta_, phi_, scale_*et_, kSolid, 2, kOrange );
    }

  for( size_t ilepton = kElectron; ilepton<=kMuon; ilepton++ )
    {
      const CandList& list_l = candList( ilepton );
      for( size_t ii=0; ii<list_l.size(); ii++ )
	{
	  Candidate& lepton_ = *list_l[ii];
	  //	  lepton_.print(cout);
	  float eta_ = lepton_.eta();
	  float phi_ = lepton_.phi();
	  float pt_  = lepton_.pt();
	  TVector2 leptonP2_ = lepton_.p2();
	  int color_ = col_electron;
	  if( ilepton==kMuon ) color_ = col_muon;
	  //      if( fillSameColor ) color2 = color1; 
	  _ed->isolation( eta_, phi_, 0., 0., eta_, phi_, 0.1*pt_/20., kSolid, 2, color_ );
	}
    }

  _ed->goToPanel( curPanel_ );
}

void
DisplayManager::draw2DVectors( float scale, float dphi, 
			       int lineWidth, int fillStyle, 
			       int color1, int color2, 
			       bool fillSameColor )
{
  const Vertex& vtx_ = primaryVertex();
  x0 = vtx_.pos().X();
  y0 = vtx_.pos().Y();
  z0 = vtx_.pos().Z();


  // ce sont les PAT jets, qui recoivent les corrections d'energie
  //  const CandList& list_j = jetList( kPatJet );
  const CandList& list_j = userCandList("cleanJets");
  // const CandList& list_j = jetList( kPatJet );
  // const CandList& list_j = jetList( kCaloJet );
  
//   vector< int > veto_jUid;
//   for( int ilep=kElectron; ilep<=kMuon; ilep++ )
//     {
//       const CandList& list_l = candList( ilep );
//       for( size_t il=0; il<list_l.size(); il++ )
// 	{
// 	  const Candidate& lCand = *list_l[il];
// 	  TVector2 lCandP2 = lCand.p2();

// 	  // find the closest pf jet
// 	  ConeSelector sel_( lCand, 0.01 );
//           CandList l_jList;
//           sel_.getList( list_j, l_jList );
//           if( l_jList.size()>1 )
//             {
//               cout << "WARNING -- More than one jet associated " 
//                    << "with " << objectName[ilep] << endl;
//             }
//           size_t ij=0;
//           for( ij=0; ij<l_jList.size(); ij++ )
//             {
//               const Candidate& jCand = *l_jList[ij];
// 	      //              float DR = jCand.p3().DeltaR( lCand.p3() );
	      
// 	      int nJC(0);
// 	      const CandInfo* info = jCand.info();
// 	      if( info->getInt( "nJC", nJC ) )
// 		{
// 		  if( nJC!=1 ) continue;
// 		}
// 	      //              veto_jUid.push_back( jCand.uid() );
//             }
// 	}
//     }
  
   if(1)
    {
      const Candidate* metCand = met( kPfMet );
      TVector2 metCandP2 = metCand->p2();
      float R_ = metCandP2.Mod();
      color1 = kBlue;
      if( fillSameColor ) color2 = color1; 
      if( showPfMet )
	{
	  //	  if( _view==hXY || _view==hPt )
	    {
	      _ed->drawRadialArrow( R_/scale, metCandP2.Phi(), 
				    dphi, 0.80, 0.60, 
				    color1, lineWidth, color2, fillStyle );
	    }
// 	  else 
// 	    {
// 	      float theta_ = Constants::pi/2.;
// 	      if( sin( metCandP2.Phi() )<0 ) theta_*=-1;
// 	      _ed->drawRadialArrow( R_/scale, theta_, 
// 				    dphi, 0.80, 0.60, 
// 				    color1, lineWidth, color2, fillStyle,
// 				    z0, 0 );
// 	    }
	}
    }  
  const float ptMin_ = PfJetEtMin;
  for( size_t ij=0; ij<list_j.size(); ij++ )
    {
      const Candidate& jCand = *list_j[ij];
      //int jUid = jCand.uid();

//       vector<int>::const_iterator it_ 
// 	= find( veto_jUid.begin(), veto_jUid.end(), jUid );
//       if( it_ != veto_jUid.end() ) continue;

      float pt_j = jCand.pt();
      if( pt_j<ptMin_ ) continue;

      //      jCand.print( cout );

      //      if( (_view==hXY&&showPfJets) || _view==hPt )
	{
	  TVector2 jCandP2 = jCand.p2();
	  color1 = kOrange;
	  color2 = kYellow;
	  if( fillSameColor ) color2 = color1; 
	  _ed->drawRadialArrow( jCandP2.Mod()/scale, jCandP2.Phi(), 
				dphi, 0.80, 0.60, 
				color2, lineWidth, color2, fillStyle,
				x0, y0 );
	}
//       else if( showPfJets )
// 	{
// 	  TVector2 jCandP2 = jCand.p2();
// 	  float jEta_ = jCand.eta();
// 	  if( fabs( jEta_ )<3. )
// 	    {
// 	      float theta_ = KineUtils::theta( jEta_ );
// 	      float R_ = jCand.pt()/sin(theta_);
// 	      float phi_ = jCand.phi();
// 	      if( sin(phi_)<0 ) theta_ *= -1;
// 	      color1 = kOrange;
// 	      color2 = kYellow;
// 	      if( fillSameColor ) color2 = color1; 
// 	      _ed->drawRadialArrow( R_/scale, theta_, 
// 				    dphi, 0.80, 0.60, 
// 				    color1, lineWidth, color2, fillStyle,
// 				    z0, 0 );
// 	    }
// 	}
    }    

  for( int ilep=kElectron; ilep<=kMuon; ilep++ )
    {
      const CandList& list_l = candList( ilep );
      for( size_t il=0; il<list_l.size(); il++ )
	{
	  const Candidate& lCand = *list_l[il];
	  TVector2 lCandP2 = lCand.p2();

	  if( ilep==kElectron ) color1 = col_electron;
	  else                  color1 = col_muon;
	  if( fillSameColor ) color2 = color1; 

	  //	  if( (_view==hXY&&showPfLeptons) || _view==hPt )
	  if( showPfLeptons )
	    {
	      _ed->drawRadialArrow( lCandP2.Mod()/scale, lCandP2.Phi(), 
				    dphi, 0.80, 0.60, 
				    color1, lineWidth, color2, fillStyle,
				    x0, y0 );
	    }
	  else if( showPfLeptons )
	    {
	      TVector2 lCandP2 = lCand.p2();
	      float theta_ = KineUtils::theta( lCand.eta() );
	      float R_ = lCand.pt()/sin(theta_);
	      float phi_ = lCand.phi();
	      if( sin(phi_)<0 ) theta_ *= -1;
	      _ed->drawRadialArrow( R_/scale, theta_, 
				    dphi, 0.80, 0.60, 
				    color1, lineWidth, color2, fillStyle );
	    }
	}
    }
  return;
}


// void
// DisplayManager::setAtlas( float ptmin )
// {
//   _ed->setAtlas();
//   lowPtTrackColor = kGray;
//   CTHadColor = kYellow-7;
//   showTrackHits = 1;
//   setThreshold( ptmin, kBlue, 1   );
//   //  setThreshold( 1.0, kGray );
//   //  setThreshold( 5.0, kWhite );  

//   view();
// }

void
DisplayManager::setThreshold( float val, int color, int width, int style )
{
  ThresAndColor x_;
  x_.pt     = val;
  x_.color  = color; 
  x_.width  = width;
  x_.style  = style;
  _thres.push_back( x_ );
  sort( _thres.begin(), _thres.end(), SortThresAndColor() );

}

void
DisplayManager::setPtMin( float ptmin )
{
  _thres.clear();
  //  setThreshold( ptmin, kGray );
  //  setAtlas( ptmin );
}

// helper function implemented in the TGV with an headache :(  19/6/11
// void
// DisplayManager::drawTexts()
// {
//   TLatex l_;
//   l_.SetTextFont(42);
//   l_.SetTextColor(kWhite);
//   l_.SetTextAlign(22);
//   float x_(0);
//   float y_(0);

//   const Candidate* metCand = met( kPatMet );
//   {
//     const Candidate* cand = metCand;
//     float  phi_ = cand->phi();
//     float pt_ = cand->pt();
//     if( _view==hPt )
//       {
// 	x_ = pt_*cos(phi_);
// 	y_ = pt_*sin(phi_);
// 	_ed->drawArrow( 0, 0, x_, y_, 
// 			kYellow, 4, kDashed, 0.5 );        
// 	x_ *= 1.2;
// 	y_ *= 1.2;
// 	l_.SetTextSize(0.07);
// 	l_.DrawLatex( x_, y_, "MET" );
//       }
//     else
//       {
// 	if( _view==hXY )
// 	  {
// 	    x_ = R*cos(phi_);
// 	    y_ = R*sin(phi_);	  
// 	    _ed->drawMET( pt_, phi_, CTScale, Emax, kYellow );
// 	    //	  }
// 	    //	else if( _view==hRZ )
// 	    //	  {
// 	    //	    y_ = L;
// 	    //	    if( sin(phi_)<0 ) L*=-1;
// 	    //	    x_ = z0;
// 	    //	  }
// 	    l_.DrawLatex( x_, y_, text_ );
// 	    TString text_;
	    
// 	    ostringstream o;
// 	    o << "MET = ";
// 	    o.precision(1);
// 	    o << fixed << cand->pt();
// 	    o << " GeV";
// 	    text_ = o.str();
// 	    l_.SetTextSize(0.05);
// 	    l_.DrawLatex( x_, y_-0.1*R, text_ );
// 	  }
//       }    
//   } 

//   for( int imode=0; imode<2; imode++ )
//     {
//       const CandList& Zlist = compositeCandList( EventServer::Zmode[imode] );
//       for( size_t ii=0; ii<Zlist.size(); ii++ )
// 	{
// 	  Candidate* cand = Zlist[ii];
// 	  float phi_ = cand->phi();
// 	  float pt_  = cand->pt();

// 	  int icol = kMagenta;
// 	  if( imode==1 ) icol = kRed;

// 	  TString text_("Z");
// 	  if( _view==hXY )
// 	    {
// 	      l_.SetTextSize(0.07);
// 	      _ed->drawMET( _view, pt_, phi_, CTScale, Emax, icol  );
// 	      x_ = R*cos(phi_);
// 	      y_ = R*sin(phi_);	  
// 	    }
// 	  else if( _view==hRZ )
// 	    {
// 	      continue;
// 	    }
// 	  else if( _view==hPt )
// 	    {
// 	      l_.SetTextSize(0.07);
// 	      x_ = pt_*cos(phi_);
// 	      y_ = pt_*sin(phi_);
// 	      _ed->drawArrow( 0, 0, x_, y_,
// 			      icol, 4, kDashed, 0.5 );	      	  
// 	      x_ *= 1.2;
// 	      y_ *= 1.2;
// 	      l_.DrawLatex( x_, y_, text_ );
// 	      continue;
// 	    }
// 	  //      l_.DrawLatex( x_, y_, text_ );
	  
// 	  l_.DrawLatex( x_, y_+0.1*R, text_ );
	  
// 	  ostringstream o;
// 	  o << "q_{T} = ";
// 	  o.precision(1);
// 	  o << fixed << cand->pt();
// 	  o << " GeV";
// 	  text_ = o.str();
// 	  l_.SetTextSize(0.05);
// 	  l_.DrawLatex( x_, y_-0.1*R, text_ );
	  
// 	  //      o << endl;
// 	}
//     }
  
//   const CandList* list_ptr(0);
//   for( int imode=0; imode<2; imode++ )
//     {
//       if( imode==0 ) list_ptr = &electrons();
//       else           list_ptr = &muons();
      
//       const CandList& list = *list_ptr;
//       for( size_t ii=0; ii<list.size(); ii++ )
// 	{
// 	  Candidate* cand = list[ii];
// 	  float phi_ = cand->phi();
// 	  float eta_ = cand->eta();
// 	  float pt_ = cand->pt();
// 	  if( pt_<20 ) continue;
// 	  TString text_;
// 	  // = cand->pdtEntry()->GetName();
// 	  int ipdg = cand->pdgCode();
	  
// 	  if( ipdg==11  ) text_ = "e^{-}";
// 	  if( ipdg==-11 ) text_ = "e^{+}";
// 	  if( ipdg==13  ) text_ = "#mu^{-}";
// 	  if( ipdg==-13 ) text_ = "#mu^{+}";
      
// 	  if( _view==hXY )
// 	    {
// 	      x_ = R*cos(phi_);
// 	      y_ = R*sin(phi_);	  
// 	    }
// 	  else if( _view==hRZ )
// 	    {
// 	      if( fabs(eta_)<2.1 )
// 		{
// 		  x_ = z0 + 0.5*L*cand->tanDip();
// 		  y_ = 0.5*L;
// 		}
// 	      else
// 		{
// 		  int sign_ = 1;
// 		  if( eta_<0 ) sign_ = -1;
// 		  x_ = z0 + 0.5*L*sign_;
// 		  y_ = 0.5*L*cand->tanTheta();
// 		}
// 	      if( sin(phi_)<0 ) y_=-L;
// 	    }
// 	  else if( _view==hPt )
// 	    {
// 	      x_ = 1.2*pt_*cos(phi_);
// 	      y_ = 1.2*pt_*sin(phi_);
	      
// 	      l_.SetTextSize(0.05);
// 	      l_.DrawLatex( x_, y_, text_ );
// 	      continue;
// 	    }
// 	  //      l_.DrawLatex( x_, y_, text_ );
	  
// 	  l_.SetTextSize(0.06);
// 	  l_.DrawLatex( x_, y_+0.1*R, text_ );
	  
// 	  ostringstream o;
// 	  if( _view==hXY )
// 	    {
// 	      o << "p_{T} = ";
// 	      o.precision(1);
// 	      o << fixed << cand->pt();
// 	      o << " GeV";
// 	    }
// 	  else if( _view==hRZ )
// 	    {
// 	      o << "#eta = ";
// 	      o.precision(1);
// 	      o << fixed << cand->eta();
// 	    }
// 	  text_ = o.str();
// 	  l_.SetTextSize(0.04);
// 	  l_.DrawLatex( x_, y_-0.1*R, text_ );
// 	  if( _view==hRZ )
// 	    {
// 	      ostringstream o;
// 	      o << "E = ";
// 	      o.precision(1);
// 	      o << fixed << cand->E();
// 	      o << " GeV";
// 	      o << flush;
// 	      text_ = o.str();
// 	      l_.SetTextSize(0.04);
// 	      l_.DrawLatex( x_, y_-0.2*R, text_ );
// 	    }
	  
// 	  //      o << endl;
// 	}
//     }
// }

void
DisplayManager::displayElectrons( float scale, float scaleE )
{
  const CandList& listElectron = candList( kElectron );
  for( size_t ii=0; ii<listElectron.size(); ii++ )
    {
      displayElectron( ii, scale, scaleE );
    }
  return;
}

bool
DisplayManager::displayElectron( size_t ii, float scale, float scaleE )
{

  bool ok = displayLepton( kElectron, ii, scale, scaleE );
  if(ok) _ed->drawComments();
  return ok;
}

bool
DisplayManager::displayLepton( int ilepton, size_t ii, float zoom, float scaleE )
{
  bool isElectron(false);
  bool isMuon(false);
  if( ilepton==kElectron )
    {
      isElectron = true;
    }
  else if( ilepton==kMuon )
    {
      isMuon     = true;
    }
  else
    abort();
  const CandList& leptons_ = candList( ilepton );
  if( ii>=leptons_.size() ) 
    {
      //      cout << "This lepton does not exist " << endl;
      return false;
    }
  Candidate& cand = *( leptons_[ii] ); 
  //   if( cand.pt()<LeptonPtMin ) return false;
  cand.oneLine(cout);
  ostringstream stream;
  cand.oneLine( stream, true );
  //   // make sure that we're dealing with leptons
  //   // and that the settings are correct
  int pdgCode_ = cand.pdgCode();
  unsigned settings(0);
  string mapname_;
  if( isElectron )
    {
      assert( fabs( pdgCode_ )==11 );
      settings = kElectron;
      mapname_ = "electron-EcalIsoDep";      
    }
  else if( isMuon )
    {
      assert( fabs( pdgCode_ )==13 );
      settings = kMuon;
      mapname_ = "muon-EcalIsoDep";
    }

  float pt  = cand.pt();
  float chg = cand.charge();
  assert( pt>0 );
  float eta = cand.eta();
  float phi = cand.phi();
   
  {
     
    TString panelName = "Hist_";
    if( isElectron )
      {
	panelName += "electron_";
      }
    else if( isMuon )
      {
	panelName += "muon_";
      }
    panelName += ii;
    _ed->setPanel( panelName, 0, 0, 800, 800, _ed->_style[0], false );
    _ed->refreshPanel( panelName );
    _ed->goToPanel( panelName, true );
     
    size_t nEcalRH = a().n("rh");
    vector<int>  v_ix;
    vector<int>  v_iy;
    vector<int>  v_iz;
    vector<float> v_e;
    for( size_t ii=0; ii<nEcalRH; ii++ )
      {   
	assert( a().load( "rh", ii ) );
	float e_    = a().get_f("rh","energy");
	int ix = a().get_i("rh","ix");
	int iy = a().get_i("rh","iy");
	int iz = a().get_i("rh","iz");
	 
	float emin=-1000;
	if( e_>emin )
	  {
	    v_ix.push_back( ix );
	    v_iy.push_back( iy );
	    v_iz.push_back( iz );
	    v_e.push_back( e_ );
	  }
      }      
     
    int hist_size = int( 20./zoom );
    if( hist_size<5 ) hist_size=5;
     
    _ed->drawElectronEcalHist( pt, eta, phi, v_ix, v_iy, v_iz, v_e, hist_size, "LEGO2" );
          
    TPad* fPad_ = goToPad( panelName );
    fPad_->cd();
    TString str_( stream.str().c_str() );
    TLatex ltx_;
    ltx_.SetNDC();
    ltx_.SetTextFont( 42 );
    ltx_.SetTextColor( kBlack );
    ltx_.SetTextSize( 0.03 );
    ltx_.SetTextAlign( 31 );
    ltx_.DrawLatex( 0.97,
		    0.03, str_ );
  }
  {
    TString panelName = "EtaPhi_";
    if( isElectron )
      {
	panelName += "electron_";
      }
    else if( isMuon )
      {
	panelName += "muon_";
      }

    float scale = 10.*zoom;
    float eta0 = eta;
    float phi0 = phi/Constants::pi;
    float DE   = 5.5/scale;
    float DPh  = 1.3/scale;
    _ed->setEtaPhiLimits( eta0, phi0, DE, DPh );

    panelName += ii;
    _ed->setPanel( panelName, 0, 0, 1000, 800, _ed->_style[0], false );
    _ed->refreshPanel( panelName );
    _ed->goToPanel( panelName, true );
     
    // the seed energy...
    const CandInfo* info = cand.info();
    float ESeed  = 5;
    if( isElectron )
      {
	ESeed  = info->getFloat("ESeed");
      }

    // use the cone selector to determine the number of PF candidates
    // of different types in a cone of 0.4 around the lepton
    float dRmax = 0.4;
    float dRmin = 0;
    float pTmin = 0.5;
    float Jhs   = -1; 
    int pdg_[3] = { 211, 130, 22 };
    const CandList& pfCandList = pfCandidates();
    for( int ii=0; ii<3; ii++ )
      {
	ConeSelector coneSel_( cand, dRmax, dRmin, pTmin, Jhs, pdg_[ii] ); 
	CandList list_;
	coneSel_.getList( pfCandList, list_ );
	//	CandUtil::printList( cout, list_ );
      }

    //
    // to be checked
    //
    float r_ecal = 129;
    float z_ecal = 280;
    if( fabs(eta)>1.5 )
      {
	r_ecal = 150;
	z_ecal = 325;
      }
    float r_hcal = 167;
    float z_hcal = 350;

    IsolationSelector isoSel_( settings );
    isoSel_.computeIso( cand );

    float DRin_trk = isoSel_.dR_trk_in;
    float DRout_trk = isoSel_.dR_trk_out;
    float Jhs_trk = isoSel_.strip_trk;

    float DRin_hcal = isoSel_.dR_hcal_in;
    float DRout_hcal = isoSel_.dR_hcal_out;
    float Jhs_hcal = -1;

    float etaAtEcal(0);
    float phiAtEcal(0);
    float etaAtHcal(0);
    float phiAtHcal(0);
    TVector3 v3 = cand.pos();
    float x0  = v3.X();
    float y0  = v3.Y();
    float z0  = v3.Z();
    KineUtils::straightToEnvelop( etaAtEcal, phiAtEcal, eta, phi, x0, y0, z0, r_ecal, z_ecal );
    KineUtils::adjust( phiAtEcal, phi );
    KineUtils::straightToEnvelop( etaAtHcal, phiAtHcal, eta, phi, x0, y0, z0, r_hcal, z_hcal );
    KineUtils::adjust( phiAtHcal, phi );

    int iecal = isoSel_.iecal;
    float DRin_ecal  = isoSel_.dR_ecal_in[iecal];
    float DRout_ecal = isoSel_.dR_ecal_out[iecal];
    float Jhs_ecal   = isoSel_.strip_ecal[iecal];

    drawEcalRecHits( ESeed/scaleE, true );
    _ed->drawEtaPhiEcal();

    float iso    = 0;
    float S      = 0;
    float S_trk  = 0;
    float S_ecal = 0;
    float S_hcal = 0;
    float a_trk  = isoSel_.a_trk;
    float a_hcal = isoSel_.a_hcal;
    float a_ecal = isoSel_.a_ecal;
    int n_trk(0);
    int n_hcal(0);
    int n_ecal(0);
  
    const CandList& listTracks = candList( kTrack );
    for( size_t ii=0; ii<listTracks.size(); ii++ )
      {
	const Candidate* trk = listTracks[ii];
	float pt_  = trk->pt();
	float eta_ = trk->eta();
	float phi_ = trk->phi();
	float deta_ = eta_-eta;
	float dR_=KineUtils::dR( eta, eta_, phi, phi_ );
	bool inCone_=true;
	if( dR_<DRin_trk )  inCone_=false;
	if( dR_>DRout_trk ) inCone_=false;
	if( fabs( deta_ )<Jhs_trk ) inCone_=false;
	//	cout << "dR/in/out " << dR_ << "/" << DRin_trk << "/" << DRout_trk << endl;
	TVector3 v3_ = trk->pos();
	float x0_  = v3_.X();
	float y0_  = v3_.Y();
	float z0_  = v3_.Z();
	float r0_ = sqrt( x0_*x0_ + y0_*y0_ );

	bool keep_=inCone_;
	if( pt_<isoSel_.pt0_trk )        keep_=false; 
	if( r0_>isoSel_.d0_trk )         keep_=false;
	if( fabs(z0_)>isoSel_.dz0_trk )  keep_=false;

	if( keep_ ) 
	  {
	    n_trk++;
	    S_trk += pt_;
	  }

	//
	// drawing
	//
	int icol_;
	if( inCone_ )
	  {
	    icol_ = ( keep_ ) ? kGreen : kMagenta;
	  }
	else 
	  {
	    icol_ = kRed;
	  }
	_ed->drawMarker( eta_, phi_/Constants::pi, 5, 2, icol_ );
      }

    float eta0_hcal = eta;
    float phi0_hcal = phi;
    KineUtils::helixToEnvelop( eta0_hcal, phi0_hcal, chg*pt, eta, phi, x0, y0, z0, r_hcal, z_hcal, BField );

    const CandList& listCT = candList( kCaloTower );
    for( size_t iCT=0; iCT<listCT.size(); iCT++ )
      {
	const Candidate& ct = *listCT[iCT];
	const CandInfo* info = ct.info();
	float eta_ = ct.eta();
	float phi_ = ct.phi();
	//float e_ = info->getFloat("energy");
	//float emEt_  = info->getFloat("emEt");
	float hadEt_ = info->getFloat("hadEt");
      
	//_ed->drawCaloTower( eta_, phi_/Constants::pi, e_, emEt_, hadEt_, Emax );

	if( hadEt_==0 ) continue;

	float hadE_ = KineUtils::e( hadEt_, eta_ );

	float dR_out=KineUtils::dR( etaAtHcal, eta_, phiAtHcal, phi_ );
	float dR_in =KineUtils::dR( eta0_hcal, eta_, phi0_hcal, phi_ );
	bool inCone_=true;
	if( dR_out>DRout_hcal ) inCone_=false;
	if( dR_in <DRin_hcal  ) inCone_=false;

	bool keep_=inCone_;
	if( hadEt_<isoSel_.et0_hcal ) keep_=false; 
	if( hadE_<isoSel_.e0_hcal )   keep_=false; 

	if( keep_ ) 
	  {
	    n_hcal++;
	    S_hcal += hadEt_;
	  }

	int icol_;
	if( inCone_ )
	  {
	    icol_ = ( keep_ ) ? kGreen : kMagenta;
	  }
	else 
	  {
	    icol_ = kRed;
	  }
	_ed->drawMarker( eta_, phi_/Constants::pi, 3, 4, icol_ );
      }


    float eta0_ecal = etaAtEcal;
    float phi0_ecal = phiAtEcal;
    //  cout << eta0_hcal << "/" << phi0_hcal << endl;
			     
    float caloEta(0);
    float caloPhi(0);
    float trackEta(0);
    float trackPhi(0);
    if( isElectron )
      {
	caloEta = info->getFloat("caloEta");
	caloPhi = info->getFloat("caloPhi");
	KineUtils::adjust( caloPhi, phi );
	trackEta = info->getFloat("trackEta");
	trackPhi = info->getFloat("trackPhi");      
	KineUtils::adjust( trackPhi, phi );
	eta0_ecal = caloEta;
	phi0_ecal = caloPhi;
      }
    else
      {
	KineUtils::helixToEnvelop( eta0_ecal, phi0_ecal, chg*pt, eta, phi, x0, y0, z0, r_ecal, z_ecal, BField );
      }

    // get the map of iso object for this candidate
    const CandIdIsoMap& map_ = isoMap(mapname_);
    CandIdIsoMapIterator it_lower_ = map_.lower_bound( cand.uid() );
    CandIdIsoMapIterator it_upper_ = map_.upper_bound( cand.uid() );
    CandIdIsoMapIterator it;
	
    bool dbg_ = false;
    if( dbg_ )
      {
	cout << "Uid=" << cand.uid() << endl;
	cout << "Number of associated ECAL deposits: " << map_.count( cand.uid() ) << endl;             
      }
    for( it=it_lower_; it!=it_upper_; it++ )
      {
	//      float dr_  = (*it).second.dr;
	float e_   = (*it).second.val; // WARNING!
	float eta_ = (*it).second.eta;
	float phi_ = (*it).second.phi;
	float et_  = KineUtils::et( e_, eta_ );
	//      float th_ = 2.0*atan(exp(-eta_));
	//      assert( th_>0 );
	//      float et_ = e_*sin(th_);

	float deta_ = eta_-eta0_ecal;
	float dR_out  = KineUtils::dR( etaAtEcal, eta_, phiAtEcal, phi_ );
	float dR_in   = KineUtils::dR( eta0_ecal, eta_, phi0_ecal, phi_ );
	//       //      cout << "deta/dphi " << deta_ << "/" << dphi_ << endl;
	//       if( dphi_ >  Constants::pi ) dphi_-=2*Constants::pi; 
	//       if( dphi_ < -Constants::pi ) dphi_+=2*Constants::pi;       
	//      float dR_ = sqrt( deta_*deta_ + dphi_*dphi_ );
	//      cout << "dphi/dr " << dphi_ << "/" << dR_ << endl;
      
	bool inCone_=true;
	if( dR_out>DRout_ecal ) inCone_=false;
	if( dR_in <DRin_ecal  ) inCone_=false;
	if( fabs( deta_ )<Jhs_ecal ) inCone_=false;

	bool keep_=inCone_;
	int iecal_ = ( fabs( eta_ )<=1.49 ) ? 0 : 1;
	if(  e_<isoSel_.e0_ecal[iecal_]  ) keep_=false;
	if( et_<isoSel_.et0_ecal         ) keep_=false;

	if( keep_ ) 
	  {
	    n_ecal++;
	    S_ecal += et_;
	  }

	int icol_;
	if( inCone_ )
	  {
	    icol_ = ( keep_ ) ? kGreen : kMagenta;
	  }
	else 
	  {
	    icol_ = kRed;
	  }
	
	_ed->drawMarker( eta_, phi_/Constants::pi, 2, 1, icol_ );

	//       //	    float phi_ = (*it).second.phi;
      }
  
    //  cout << "track : eta/phi " << trackEta << "/" << trackPhi/Constants::pi << endl;

    if( isElectron )
      {
	_ed->drawMarker( trackEta, trackPhi/Constants::pi, 30, 1., 
			 //		       kBlue );
			 kMagenta );
	//      cout << "calo : eta/phi " << caloEta << "/" << caloPhi/Constants::pi << endl;
	_ed->drawMarker( caloEta, caloPhi/Constants::pi, 30, 1., kBlack );
      }

    //    cand.oneLine(cout);
    S = a_trk*S_trk + a_ecal*S_ecal + a_hcal*S_hcal;
    iso = pt/( pt + S );
    cout << "S_trk  = "   << S_trk   << " /" << n_trk << endl;
    cout << "S_hcal = "   << S_hcal  << " /" << n_hcal << endl;
    cout << "S_ecal = "   << S_ecal  << " /" << n_ecal << endl;
    cout << "S      = "   << S       << endl;
    cout << "iso=" << iso << "/" << isoSel_.iso << endl;
  
    //  drawTracks( kViolet, 1 );
    //    cout << "OK-- draw Tracks " << endl;
    //    drawTracks();
    //    drawTrack( cand, kBlack, 4, kSolid );

    int linestyle_, linewidth_, linecolor_;

    //    linestyle_ = kSolid;
    //    linewidth_ = 1;
    //  linecolor_ = kBlue;
    //    linecolor_ = kMagenta;
    //    _ed->isolation( eta, phi, DRin_trk, Jhs_trk, eta, phi, DRout_trk, linestyle_, linewidth_, linecolor_ );

    linestyle_ = kDashed;
    linewidth_ = 1;
    linecolor_ = kRed;
    _ed->isolation( eta0_hcal, phi0_hcal, DRin_hcal, Jhs_hcal, etaAtHcal, phiAtHcal, DRout_hcal, linestyle_, linewidth_, linecolor_ );

    linestyle_ = kSolid;
    linewidth_ = 2;
    linecolor_ = kBlack;
    //  if( isElectron )
    //    {
    _ed->isolation( eta0_ecal, phi0_ecal, DRin_ecal, Jhs_ecal, etaAtEcal, phiAtEcal, DRout_ecal, linestyle_, linewidth_, linecolor_ ); 
    //    }
    //  else
    //    {
    //      _ed->isolation( etaAtEcal, phiAtEcal, DRin_ecal, DRout_ecal, Jhs_ecal, linestyle_, linewidth_, linecolor_ );
    //    }

    //    drawEcalRecHits();



    TPad* fPad_ = goToPad( panelName );
    fPad_->cd();
    TString str_( stream.str().c_str() );
    TLatex ltx_;
    ltx_.SetNDC();
    ltx_.SetTextFont( 42 );
    ltx_.SetTextColor( kBlack );
    ltx_.SetTextSize( 0.03 );
    ltx_.SetTextAlign( 31 );
    ltx_.DrawLatex( 0.97,
		    0.03, str_ );
     
    _ed->setEtaPhiLimits();

  }
  return true;
}


//   return true;
// }

void
DisplayManager::zoomXYZ( float zoom, float x0, float y0, float z0 )
{
  if( zoom<=0 ) return;
  _ed->setXYZLimits(); // resets the defaults
  _ed->setXYZLimits( x0, y0, z0, _ed->_R/zoom );
  _ed->refreshPanel( "XY", true );
  _ed->refreshPanel( "RZ", true );
  _ed->XY();
  _ed->RZ();
  if( goToPanel( "XZ" ) )
    {
      _ed->refreshPanel( "XZ", true );
      _ed->refreshPanel( "YZ", true );
      _ed->RZ("XZ");
      _ed->RZ("YZ");
    }
  view();
}

void
DisplayManager::zoomVertex( float zoom )
{
  if( zoom<=0 ) return;
  _ed->setVertexLimits(); // resets the defaults
  _ed->setVertexLimits( _primaryVertex->pos().X(), 
			_primaryVertex->pos().Y(),
			_primaryVertex->pos().Z(),  
			_ed->_RV/zoom );
  if( goToPanel( "VXY" ) )
    {
      _ed->refreshPanel( "VXY", true );
      _ed->refreshPanel( "VRZ", true );
      _ed->refreshPanel( "VXZ", true );
      _ed->refreshPanel( "VYZ", true );
      _ed->RZ("VXZ");
      _ed->RZ("VYZ");
      _ed->RZ("VRZ");
      _ed->XY("VXY");
      _ed->drawScales();
      view();
    }
}

void
DisplayManager::energyScales( float pt_scale, float E_scale, float CT_scale )
{
  _ed->setEnergyLimits(); // resets the defaults
  _ed->setEnergyLimits( _ed->_ptmax*pt_scale, 
			_ed->_emax*E_scale, 
			_ed->_CTScale*CT_scale ); 
  _ed->refreshPanel("Pt", true);
  _ed->Pt();
  view();
}
