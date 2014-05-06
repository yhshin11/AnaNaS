#include <cassert>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <string>
#include <cmath>
using namespace std;

#include <Math/VectorUtil.h>

#include "Analysis/utils/Config.hh"
#include "Analysis/utils/CfgParser.hh"

#include "utils/TreeManager.hh"
#include "Analysis/core/EventServer.hh"

int    EventServer::verbosity = 0;  // was zero
string EventServer::Zmode[2]  = { "Zee", "Zmumu" };
string EventServer::Wmode[2]  = { "Wen", "Wmn" };
string EventServer::lepton[2] = { "Electron", "Muon" };

bool   EventServer::isMC = true;

std::vector< std::string > EventServer::caloMET;
std::vector< std::string > EventServer::recoMET;
std::vector< std::string > EventServer::pfMET;
std::vector< std::string > EventServer::patMET;
std::vector< std::string > EventServer::genMET;
std::vector< std::string > EventServer::caloJets;
std::vector< std::string > EventServer::pfJets;
std::vector< std::string > EventServer::pfCandidates;
std::vector< std::string > EventServer::trackJets;
std::vector< std::string > EventServer::genJets;

bool EventServer::parsed_ = false;

int EventServer::caloMET_default;
int EventServer::recoMET_default;
int EventServer::pfMET_default;
int EventServer::patMET_default;
int EventServer::genMET_default;
int EventServer::caloJets_default;
int EventServer::pfJets_default;
int EventServer::pfCandidates_default;
int EventServer::trackJets_default;
int EventServer::genJets_default;

ClassImp(EventServer)


EventServer::EventServer( EventList& events, const char * collectionFileName ) : _treeMgr(0)
{
  // show content:
  EventList::iterator it;
  for ( it=events.begin() ; it != events.end(); it++ )
    {
      string filename_ = (*it).first;
      cout << filename_ << endl;

      _filenames.push_back( filename_ );
      int ifile_ = _filenames.size()-1;
      //       es.openFile( filename_.c_str() );
      
       vector<int> vect_ = (*it).second;
       for( size_t ii=0; ii<vect_.size(); ii++ )	    
 	{
 	  cout << "\t" << vect_[ii] << endl;
	  _eventList.insert( pair< int, int >( ifile_, vect_[ii]+1 ) );
 	}
    }
  init( collectionFileName );
}  

EventServer::EventServer( const char* filename, const char * collectionFileName ) : _treeMgr(0)
{
  _filenames.push_back( filename );
  init( collectionFileName );
}  

EventServer::EventServer( vector< string > filenames, const char * collectionFileName ) : _filenames( filenames ), _treeMgr(0)
{
  init( collectionFileName );

}  

void
EventServer::init( const char * collectionFileName )
{
  if ( ! parsed_ ) {
    
    //    string fileName_( Config::confPath );
    //    fileName_ += collectionFileName;

    CfgParser parser( collectionFileName );
    parser.values( "CaloMET", caloMET );
    caloMET_default = parser.defaultValue( "CaloMET" );
    parser.values( "RecoMET", recoMET );
    recoMET_default = parser.defaultValue( "RecoMET" );
    parser.values( "PFMET", pfMET );
    pfMET_default = parser.defaultValue( "PFMET" );
    parser.values( "PatMET", patMET );
    patMET_default = parser.defaultValue( "PatMET" );
    parser.values( "CaloJets", caloJets );
    caloJets_default = parser.defaultValue( "CaloJets" );
    parser.values( "PFCandidates", pfCandidates );
    pfCandidates_default = parser.defaultValue( "PFCandidates" );
    parser.values( "PFJets", pfJets );
    pfJets_default = parser.defaultValue( "PFJets" );
    parser.values( "TrackJets", trackJets );
    trackJets_default = parser.defaultValue( "TrackJets" );
    parsed_ = true;

    isMC = parser.hasKey( "isMC" );
    if( isMC )
      {
	// check if the "mcTruth" tree exists
	//	cout << "This is Monte-Carlo " << endl;
	parser.values( "GenJets", genJets );
	genJets_default = parser.defaultValue( "GenJets" );
	parser.values( "GenMET", genMET );
	genMET_default = parser.defaultValue( "GenMET" );
	//genMET_default = parser.Value( "GenMet" );
      }
    parser.print( std::cout );

  }

  _filesChecked = false;
  assert( checkFiles() );
  curRun      = 0;
  curLumi     = 0;
  curEvent    = 0;
  ifile       = 0;
  ievent      = 0; 
  eventNumber = 0;

  _hasEventList = false;
  if( !_eventList.empty() )
    {
      _hasEventList = true;
      _eventListIt = _eventList.begin();
    }
  if(_filenames.size()==0){
    //    std::cerr << __FILE__ << ":" << __LINE__ << ": no data file! Exits."
    //              << std::endl;
    //    exit(1);
    cout << "No input file " << endl;
    return;
  }
  openFile( _filenames[0] );
}

bool
EventServer::checkFiles()
{
  assert( !_filesChecked );
  _filesChecked = true;

  TFile* _file_(0); 
  TTree* _ttree_(0); 

  for( size_t ii=0; ii<_filenames.size(); ii++ )
    {
      string _filename_ = _filenames[ii];
      if(_filename_.find("eos")!=(size_t)-1) _filename_ = "root://eoscms/"+_filename_;
      _file_ = TFile::Open( _filename_.c_str() );
      if( _file_==0 ) 
	{
	  cout << "file " << _filename_ << " not found " << endl;
	  return false;
	}
      
      // run summary
      string str_ = "eventSummary";
      _ttree_ = (TTree*) _file_->Get( str_.c_str() );
      int n_ = (int)(_ttree_->GetEntriesFast());
      if( ii>0 ) n_+= _nevent[ii-1];
      _nevent.push_back( n_ );             

      // GHM 8/05/11 -- needed in case of skims
      _file_->Close();

    }
  assert( _nevent.size()==_filenames.size() );
  for( size_t ii=0; ii<_filenames.size(); ii++ )
    {
      /*cout << "n[" << _filenames[ii] 
	   << "]=" << _nevent[ii]
	   << endl;*/ // Laurent stop the printing
    }
  cout << "Checking input files....done" << endl;
  return true;
}

void 
EventServer::closeFile() 
{ 
  hlt.clear();
  eid.clear();
  eAcc.clear();
  muid.clear();
  tauid.clear();
  phid.clear();
  evt_n.clear();
  evt_first.clear();

  delete _treeMgr;
  map_i.clear();
  map_f.clear();
  map_s.clear();
  map_b.clear();
  map_vi.clear();
  map_vf.clear();
}


bool
EventServer::openFile( const string & filename )
{
  cout << "Open file " << filename << endl;
 
  closeFile();
  string filename_=filename;
  if(filename_.find("eos")!=(size_t)-1) filename_ = "root://eoscms/"+filename_;
  _treeMgr = new TreeManager( filename_.c_str(), TreeManager::kRead );

  fName =  filename;
  
  vector<string> run_name_;
  
  // event
  register_t( "evt", "eventSummary" );
  register_i( "evt", "run" );
  register_i( "evt", "lumiBlock" );
  register_i( "evt", "event" );
  register_i( "evt", "processID" );
  register_i( "evt", "time" );
  register_f( "evt", "ptHat" );
  register_f( "evt", "rhoFastJet" );

  register_f( "evt", "trueNI" );

  register_i( "evt", "nInteracBXm1" );
  register_i( "evt", "nInteracBX" );
  register_i( "evt", "nInteracBXp1" );

  {
    // hlt
    run_name_.push_back( "hltTriggerNames" );
    register_t( "hlt", "triggerResults"  );
    
    register_b( "hlt", "hltBits" ); 
  }

  if( isMC )
    { 
      // mc truth
      run_name_.push_back( "mcTruth" );
      register_t( "mc", "mcTruth" );
      register_i( "mc", "run"   );
      register_i( "mc", "event" );
      register_i( "mc", "pdgId" );
      register_i( "mc", "nDau"  );
      register_i( "mc", "qx3"   );
      register_f( "mc", "mass"  );
      register_f( "mc", "pt"    );
      register_f( "mc", "eta"   );
      register_f( "mc", "phi"   );
      register_f( "mc", "vx"    );
      register_f( "mc", "vy"    );
      register_f( "mc", "vz"    );
      register_i( "mc", "motherIndex" );
      register_i( "mc", "status" );
      register_i( "mc", "collId" );
    }



  if( _treeMgr->t("ntu_mcTruth")==0 &&  
      _treeMgr->t("ntu_genMET_genMetTrue")==0 && 
      _treeMgr->t("ntu_genJets_ak5GenJets")==0 ) 
    {
      //      cout << "No mcTruth tree -- this is not MC" << endl;
      isMC = false;
    }
  //  cout<<" ====================> "<<isMC<<endl;
  if( isMC )
    {
      // generated event
      run_name_.push_back( "hepMCProduct" );
      register_t( "gen", "genEvent"  );
      register_i( "gen", "run"   );
      register_i( "gen", "event"   );
      register_i( "gen", "pdgId"   );
      register_i( "gen", "status"   );
      register_i( "gen", "PV"   );
      register_i( "gen", "EV"   );
      register_f( "gen", "px"  );
      register_f( "gen", "py"  );
      register_f( "gen", "pz"  );
      register_f( "gen", "e"  );
      register_f( "gen", "m"  );
      register_f( "gen", "PV_x"  );
      register_f( "gen", "PV_y"  );
      register_f( "gen", "PV_z"  );
      register_f( "gen", "EV_x"  );
      register_f( "gen", "EV_y"  );
      register_f( "gen", "EV_z"  );
    }
  if( isMC )
    {
      for( size_t v=0; v<genMET.size(); v++ )
	{
	  // gen met
	  register_t( "gmet", v, "genMET", genMET[v] );      
	  run_name_.push_back( getName( "gmet", v ) );      
	  register_i( "gmet", "run"    , v );
	  register_i( "gmet", "event"  , v );
	  register_f( "gmet", "pt"     , v );
	  register_f( "gmet", "phi"    , v );
	  register_f( "gmet", "sumEt"  , v );
	  register_f( "gmet", "mEtSig" , v );
	  register_f( "gmet", "em"     , v );
	  register_f( "gmet", "had"    , v );
	  register_f( "gmet", "inv"    , v );
	}
    }
 if( isMC ) 
    {
      for( size_t v=0; v<genJets.size(); v++ )
	{
	  // gen jets
	  register_t( "gj", v, "genJets", genJets[v] );      
	  run_name_.push_back( getName( "gj", v ) );      
	  register_i( "gj", "run"  , v );
	  register_i( "gj", "event"  , v );
	  register_f( "gj", "pt"  , v );
	  register_f( "gj", "eta"  , v );
	  register_f( "gj", "phi"  , v );
	  register_f( "gj", "em"  , v );
	  register_f( "gj", "had"  , v );
	  register_f( "gj", "inv"  , v );
	  register_i( "gj", "nGen"  , v );
	  register_vi( "gj", "barcodeGen"  , v );
	}
    }
 
 

  {
    // PF super clusters
    run_name_.push_back( "pfSuperCluster" );
    register_t( "pfSc", "pfSuperClusters" );     
    register_i( "pfSc", "run"  );
    register_i( "pfSc", "event"  );
    register_f( "pfSc", "rawEnergy"  );
    register_f( "pfSc", "preshowerEnergy"  );
    register_f( "pfSc", "energy"  );
    register_f( "pfSc", "eta"  );
    register_f( "pfSc", "phi"  );
    register_f( "pfSc", "etaWidth"  );
    register_f( "pfSc", "phiWidth"  );
    register_i( "pfSc", "nHits"  );
    register_f( "pfSc", "R4" );
    register_f( "pfSc", "R9" );
    register_i( "pfSc", "index" );
    register_vi( "pfSc", "indexBasicCluster" );
  }
  {
    // PF basic Clusters
    run_name_.push_back( "pfBasicCluster" );
    register_t( "pfBc", "pfBasicClusters" );    
    register_i( "pfBc", "run"  );
    register_i( "pfBc", "event"  );
    register_f( "pfBc", "energy"  );
    register_f( "pfBc", "eta"  );
    register_f( "pfBc", "phi"  );
    register_f( "pfBc", "R4" );
    register_f( "pfBc", "R9" );
    register_i( "pfBc", "index" );
    register_vi( "pfBc", "indexHit" );
    register_vf( "pfBc", "fractionHit" );
    register_i( "pfBc", "nHits" );
  }
  {      
    // ECAL super clusters
    run_name_.push_back( "ecalSuperClusterEB" );
    run_name_.push_back( "ecalSuperClusterEE" );
    register_t( "sc", "ecalSuperClusters" );     
    register_i( "sc", "run"  );
    register_i( "sc", "event"  );
    register_f( "sc", "rawEnergy"  );
    register_f( "sc", "preshowerEnergy"  );
    register_f( "sc", "energy"  );
    register_f( "sc", "eta"  );
    register_f( "sc", "phi"  );
    register_f( "sc", "etaWidth"  );
    register_f( "sc", "phiWidth"  );
    register_f( "sc", "R4" );
    register_f( "sc", "R9" );
    register_i( "sc", "nHits"  );
    register_i( "sc", "index" );
    register_vi( "sc", "indexBasicCluster" );
    register_i( "sc", "ixSeed" );
    register_i( "sc", "iySeed" );
    register_i( "sc", "izSeed" );
    register_i( "sc", "ixClosestCell" );
    register_i( "sc", "iyClosestCell" );
    register_i( "sc", "izClosestCell" );
    
  }
  {
    // ECAL basic Clusters
    run_name_.push_back( "ecalBasicClusterEB" );
    run_name_.push_back( "ecalBasicClusterEE" );
    register_t( "bc", "ecalBasicClusters" );    
    register_i( "bc", "run"  );
    register_i( "bc", "event"  );
    register_f( "bc", "energy"  );
    register_f( "bc", "eta"  );
    register_f( "bc", "phi"  );
    register_f( "bc", "R4" );
    register_f( "bc", "R9" );
    register_f( "bc", "covEtaEta" );
    register_f( "bc", "covEtaPhi" );
    register_f( "bc", "covPhiPhi" );
    register_f( "bc", "eBottom" );
    register_f( "bc", "eLeft" );
    register_f( "bc", "eTop" );
    register_f( "bc", "eRight" );
    register_f( "bc", "eMax" );
    register_i( "bc", "index" );
    register_vi( "bc", "indexHit" );
    register_i( "bc", "nHits" );

  }
  {
    // ECAL preshower Clusters
    run_name_.push_back( "ecalPreshowerBasicClusters" );
    register_t( "psc", "ecalPreshowerBasicClusters" );    
    register_i( "psc", "run"  );
    register_i( "psc", "event"  );
    register_f( "psc", "energy"  );
    register_f( "psc", "eta"  );
    register_f( "psc", "phi"  );
    register_i( "psc", "plane" );
    register_i( "psc", "index" );
    register_i( "psc", "indexBC" );
    register_vi( "psc", "indexHit" );
    register_i( "psc", "nHits" );

    }
  {
    // ecal rec hits
    run_name_.push_back( "ecalRecHitEB" );
    run_name_.push_back( "ecalRecHitEE" );
    register_t( "rh", "ecalRecHits" );     
    register_i( "rh", "run" );
    register_i( "rh", "event" );
    register_i( "rh", "ix" );
    register_i( "rh", "iy" );
    register_i( "rh", "iz" );
    register_f( "rh", "energy" );
    register_i( "rh", "index" );
    register_i( "rh", "time" );
    
    //Spikes variables
     register_i( "rh", "recoFlag" );
     register_f( "rh", "chi2" );
     register_f( "rh", "outOfTimeChi2" );
     register_f( "rh", "outOfTimeEnergy" );
     register_i( "rh", "severityLevel" );
     register_f( "rh", "E1OverE9" );
     register_f( "rh", "swissCross" );
     register_f( "rh", "dbStatus" );

     //Laser variables
     register_f( "rh", "laserCorr" );
     register_f( "rh", "apdpnref"  );
     register_f( "rh", "apdpn_p2" );
     register_f( "rh", "apdpn_t2" );
     register_f( "rh", "alpha"     );
     
  }
  {
    // ecal preshower rec hits
    run_name_.push_back( "ecalPreshowerRecHits" );
    register_t( "rhps", "ecalPreshowerRecHits" );     
    register_i( "rhps", "run" );
    register_i( "rhps", "event" );
    register_i( "rhps", "ix" );
    register_i( "rhps", "iy" );
    register_i( "rhps", "iz" );
    register_f( "rhps", "energy" );
    register_i( "rhps", "index" );
    // register_f( "rhps", "time" );
    register_i( "rhps", "plane" );
    register_i( "rhps", "strip" );
  }
  {
    // electrons
    run_name_.push_back( "electron" );
    register_t(  "el"   , "electrons"   );  
    register_i( "el", "run"   );
    register_i( "el", "event" );
    register_i( "el", "index" );
    register_i( "el", "class" );
    register_f( "el", "pt"    );
    register_f( "el", "eta"   );
    register_f( "el", "phi"   );
    register_f( "el", "vx"    );
    register_f( "el", "vy"    );
    register_f( "el", "vz"    );
    register_f( "el", "pin"   );
    register_f( "el", "fBrem" );
    register_f( "el", "caloEnergy" );
    register_f( "el", "EOP" );
    register_f( "el", "ESeed" );
    register_f( "el", "ESeedOPout" );
    register_f( "el", "HOE" );
    register_f( "el", "caloRho" );
    register_f( "el", "caloEta" );
    register_f( "el", "caloPhi" );
    register_f( "el", "trackRho" );
    register_f( "el", "trackEta" );
    register_f( "el", "trackPhi" );
    register_f( "el", "dPhiIn" );
    register_f( "el", "dEtaIn" );
    register_f( "el", "sigiEtaiEta" );
    register_i( "el", "nTrk" );
    register_f( "el", "sumOfPt" );
    register_f( "el", "sumOfPt2" );
    
    register_i( "el", "vetoID" );
    register_i( "el", "looseID" );
    register_i( "el", "mediumID" );
    register_i( "el", "tightID" );

    register_b( "el", "electronIDs" );
    run_name_.push_back( "electronIDs" );


    register_b( "el", "elAccBits" );
    run_name_.push_back( "electronFiduFs" );

    //MM ======= isolation variables
 //    register_i(  "el", "nECALdep"   );
//     register_vf( "el", "dRECALdep"  );
//     register_vf( "el", "E_ECALdep"  );
//     register_vf( "el", "Et_TECALdep" );
//     register_vf( "el", "etaECALdep" );
//     register_vf( "el", "phiECALdep" );

//     register_i(  "el", "nHCALdep"   );
//     register_vf( "el", "dRHCALdep"  );
//     register_vf( "el", "valHCALdep" );
//     register_vf( "el", "etaHCALdep" );
//     register_vf( "el", "phiHCALdep" );
//     register_i(  "el", "ntrkdep"    );
//     register_vf( "el", "dRtrkdep"  );
//     register_vf( "el", "valtrkdep" );
//     register_vf( "el", "etatrkdep" );
//     register_vf( "el", "phitrkdep" );

    // Default isolation variables
    register_f( "el", "dr03TrkIso" );
    register_f( "el", "dr04EcalIso" );
    register_f( "el", "dr04HcalD1Iso" );
    register_f( "el", "dr04HcalD2Iso" );
    register_f( "el", "dr03EcalIso" );
    register_f( "el", "dr03HcalD1Iso" );
    register_f( "el", "dr03HcalD2Iso" );

    register_f( "el", "chIso");
    register_f( "el", "nhIso");
    register_f( "el", "phIso");
    register_f( "el", "puIso");
    //======= 

    //PF related variables
    register_b( "el", "ecalDriven" );
    register_b( "el", "trkDriven" );
    register_i( "el", "indexPfCandidate" );

    //Clusters    
    register_i( "el", "indexSC" );
    register_i( "el", "indexPfSC" );

    register_i(  "el", "nclus"    );
    register_vf( "el", "eclus"   );
    register_vf( "el", "rhoclus" );
    register_vf( "el", "etaclus" );
    register_vf( "el", "phiclus" );
    //  
    //Tracks 
    register_i(  "el", "nTrack"    );
    register_vi( "el", "indexTracks"   );
    register_vf( "el", "trkdr"   );
    
    register_i( "el", "indexGsfTrack" );
    register_i( "el", "indexTrack" );

    register_b( "el", "isConv" );
  }

   {
    // muons
    run_name_.push_back( "muon" );
    register_t(  "mu", "muons" );      
    register_i(  "mu", "run"   );
    register_i(  "mu", "event" );
    register_i(  "mu", "index" );
    register_f(  "mu", "pt"    );
    register_f(  "mu", "eta"   );
    register_f(  "mu", "phi"   );
    register_f(  "mu", "vx"    );
    register_f(  "mu", "vy"    );
    register_f(  "mu", "vz"    );
    register_f(  "mu", "dir"   );
    register_f(  "mu", "dxyBS" );
    register_f(  "mu", "dxyPV" );
    register_i(  "mu", "nch"   );
    register_f(  "mu", "em"    );
    register_f(  "mu", "had"   );
    register_i(  "mu", "nTrk" );
    register_i(  "mu", "nValMuHits" );
    register_i(  "mu", "nValTrackerHits" );
    register_f(  "mu", "sumOfPt" );
    register_f(  "mu", "sumOfPt2" );

    register_b(  "mu", "muonIDs" );
    run_name_.push_back( "muonIDs" );

    register_b(  "mu", "isPFMuon" );

    register_f(   "mu", "normChi2" );
    
    //Quality variables
    register_i(  "mu", "nTrkHits" );
    register_i(  "mu", "nPixHits" );
    register_i(  "mu", "nTrkLayerHits" );
    register_i(  "mu", "nPixLayerHits" );
    register_i(  "mu", "nMuonHits");
    register_i(  "mu", "nMatch");
    
    //Isolation variables ===========
    
    //Default variables
    register_f( "mu", "isoR03strk" );
    register_f( "mu", "isoR03sem" );
    register_f( "mu", "isoR03shad" );
    register_f( "mu", "isoR03strkVeto" );
    register_f( "mu", "isoR03semVeto" );
    register_f( "mu", "isoR03shadVeto" );
    register_f( "mu", "isoR03sum" );
    register_f( "mu", "isoR03pfch" );
    register_f( "mu", "isoR03pfcpart" );
    register_f( "mu", "isoR03pfnh" );
    register_f( "mu", "isoR03pfph" );
    register_f( "mu", "isoR03pfnhThres" );
    register_f( "mu", "isoR03pfphThres" );
    register_f( "mu", "isoR03pfpu" );
    register_f( "mu", "isoR04pfch" );
    register_f( "mu", "isoR04pfcpart" );
    register_f( "mu", "isoR04pfnh" );
    register_f( "mu", "isoR04pfph" );
    register_f( "mu", "isoR04pfnhThres" );
    register_f( "mu", "isoR04pfphThres" );
    register_f( "mu", "isoR04pfpu" );
    //==============================
    
    //Mu track hits
    register_i( "mu", "nTHits" );
    register_vi( "mu", "tid" );
    register_vf( "mu", "trho" );
    register_vf( "mu", "teta" );
    register_vf( "mu", "tphi" );


    register_i(    "mu", "indexTrack" );
    register_i(    "mu", "indexPfCandidate" );

    register_i(    "mu", "nTrack"    );
    register_vi(   "mu", "indexTracks"  );
    register_vf(   "mu", "trkdr"  );
  }

  { //taus
    run_name_.push_back( "tau" );
    register_t( "tau", "taus"     );      
    register_i( "tau", "run"    );
    register_i( "tau", "event"  );
    register_i( "tau", "index" );
    register_f( "tau", "pt"  );
    register_f( "tau", "eta"  );
    register_f( "tau", "phi"  );
    register_f( "tau", "energy"  );
    register_f( "tau", "caloEta"  );
    register_f( "tau", "caloPhi"  );
    register_f( "tau", "vx"  );
    register_f( "tau", "vy"  );
    register_f( "tau", "vz"  );
    register_f( "tau", "leadPPt"  );
    register_f( "tau", "leadCPPt"  );
    register_f( "tau", "leadNPPt"  );
    register_f( "tau", "trkd0"  );
    register_f( "tau", "trkdz"  );
    register_i( "tau", "nCHSigC"  );
    register_i( "tau", "nNHSigC"  );
    register_i( "tau", "nPhSigC"  );
    register_i( "tau", "nPSigC"  );
    register_i( "tau", "nCHIsoC"  );
    register_i( "tau", "nNHIsoC"  );
    register_i( "tau", "nPhIsoC"  );
    register_i( "tau", "nPIsoC"  );
    register_f( "tau", "sPtCHIsoC"  );
    register_f( "tau", "sPtPhIsoC"  );

    register_b( "tau", "tauIDs"  );
    run_name_.push_back( "tauIDs" );

    register_f( "tau", "pfElMVA"  );
    register_f( "tau", "pfJetPt"  );
    register_f( "tau", "pfEtaPt"  );
    register_f( "tau", "pfPhiPt"  );
  }
  
  
  {
    // photons
    run_name_.push_back( "photon" );
    register_t( "ph", "photons"     );      
    register_i( "ph", "run"    );
    register_i( "ph", "event"  );
    register_f( "ph", "pt"  );
    register_f( "ph", "eta"  );
    register_f( "ph", "phi"  );
    register_f( "ph", "caloRho"  );
    register_f( "ph", "caloEta"  );
    register_f( "ph", "caloPhi"  );

    register_b( "ph", "photonIDs"  );
    run_name_.push_back( "photonIDs" );

    register_i(   "ph", "nTrack" );
    register_vi(  "ph", "indexTracks"  );
    register_vf(  "ph", "trkdr" );

    // ID variables ========
    // register_f( "ph", "r9" );
    //  register_f( "ph", "e5x5" );
    //  register_f( "ph", "e3x3" );
    register_f( "ph", "e1x5" );
    register_f( "ph", "e2x5" );
    register_f( "ph", "sigmaEtaEta" );
    register_f( "ph", "sigmaIetaIeta" );
    //======================

    //Isolation variables ====
    register_f( "ph", "ecalRecHitSumEtConeDR03" );
    register_f( "ph", "ecalRecHitSumEtConeDR04" );
    register_f( "ph", "hadronicDepth1OverEm" );
    register_f( "ph", "hadronicDepth2OverEm" );
    register_f( "ph", "hadronicOverEm" );
    register_f( "ph", "hcalDepth1TowerSumEtConeDR03" );
    register_f( "ph", "hcalDepth1TowerSumEtConeDR04" );
    register_f( "ph", "hcalDepth2TowerSumEtConeDR03" );
    register_f( "ph", "hcalDepth2TowerSumEtConeDR04" );
    register_f( "ph", "hcalTowerSumEtConeDR03" );
    register_f( "ph", "hcalTowerSumEtConeDR04" );
    register_f( "ph", "nTrkHollowConeDR03" );
    register_f( "ph", "nTrkHollowConeDR04" );
    register_f( "ph", "nTrkSolidConeDR03" );
    register_f( "ph", "nTrkSolidConeDR04" );
    register_f( "ph", "trkSumPtHollowConeDR03" );
    register_f( "ph", "trkSumPtHollowConeDR04" );
    register_f( "ph", "trkSumPtSolidConeDR03" );
    register_f( "ph", "trkSumPtSolidConeDR04" );
    register_f( "ph", "ecalIso" );
    //========================
    register_b( "ph", "isEB" );
    register_b( "ph", "isEE" );
    register_b( "ph", "isEBGap" );
    register_b( "ph", "isEEGap" );
    register_b( "ph", "isEBEEGap" );
    //    register_b( "ph", "" );
    //========================
    
    register_i( "ph", "indexSC" );
    register_b( "ph", "hasPixelSeed" );
    register_b( "ph", "isConverted" );
  }
  
  for( size_t v=0; v<pfCandidates.size(); v++ )
    {
      //PF Candidates
      register_t(  "pfc", v, "pfCandidates", pfCandidates[v] );
      run_name_.push_back( getName( "pfc", v ) );
      register_i(  "pfc", "run", v );
      register_i(  "pfc", "event", v );
      register_i(  "pfc", "index", v );
      
      register_i(  "pfc", "pdgId", v );
      register_f(  "pfc", "pt", v );
      register_f(  "pfc", "eta", v );
      register_f(  "pfc", "phi", v );
      register_f(  "pfc", "vx", v );
      register_f(  "pfc", "vy", v );
      register_f(  "pfc", "vz", v );
      
      register_i(  "pfc", "vtxIndex", v );

      register_f(  "pfc", "ecalE", v );
      register_f(  "pfc", "hcalE", v );
      register_f(  "pfc", "pS1E", v );
      register_f(  "pfc", "pS2E", v );
      register_f(  "pfc", "deltaP", v );
      register_f(  "pfc", "mva_e_mu", v );
      register_f(  "pfc", "mva_e_pi", v );
      register_f(  "pfc", "mva_pi_mu", v );
      register_f(  "pfc", "mva_nothing_gamma", v );
      register_f(  "pfc", "mva_nothing_nh", v );
      register_f(  "pfc", "mva_gamma_nh", v );
      
      register_f(  "pfc", "caloRho", v );
      register_f(  "pfc", "caloEta", v );
      register_f(  "pfc", "caloPhi", v );
      
    }

  {
    // vertices
    run_name_.push_back( "vertex" );
    register_t(  "vtx", "vertices"    );      
    register_i(  "vtx", "run"   );
    register_i(  "vtx", "event" );
    register_i(  "vtx", "index" );
    register_f(  "vtx", "x"    );
    register_f(  "vtx", "y"   );
    register_f(  "vtx", "z"   );
    register_i(  "vtx", "nTracks"   );
    register_vf( "vtx", "tpt"  );
    register_vf( "vtx", "teta" );
    register_vf( "vtx", "tphi" );
    register_vf( "vtx", "twgh" );

    //quality variables
    register_f(  "vtx", "ndof" );
    register_b(  "vtx", "isFake" );
    register_b(  "vtx", "isGoodVtx" );

  }

  {
      // gsf tracks 
      register_t( "gsf", "gsfTracks"   );     
      register_i( "gsf", "run");
      register_i( "gsf", "event");
      register_i( "gsf", "index");
      register_i( "gsf", "nLostHits");
      register_i( "gsf", "nValidHits");
      register_i( "gsf", "matchedTrack");
      register_f( "gsf", "pt");
      register_f( "gsf", "eta");
      register_f( "gsf", "phi");
      register_f( "gsf", "vx");
      register_f( "gsf", "vy");
      register_f( "gsf", "vz");
      register_f( "gsf", "d0");
      register_f( "gsf", "z0");
      register_f( "gsf", "innerPoint_x");
      register_f( "gsf", "innerPoint_y");
      register_f( "gsf", "innerPoint_z");
      register_f( "gsf", "outerPoint_x");
      register_f( "gsf", "outerPoint_y");
      register_f( "gsf", "outerPoint_z");
      register_f( "gsf", "normChi2");
      register_b( "gsf", "isFPHitTrk" );
      register_i( "gsf", "expInnHits" );

      register_i(  "gsf", "nTHits"    );
      register_vi( "gsf", "tid"   );
      register_vf( "gsf", "trho" );
      register_vf( "gsf", "teta" );
      register_vf( "gsf", "tphi" );
  }  

  {
    // refitted tracks
    register_t(  "trk"  , "tracks"      );     
    register_i(  "trk", "run"   );
    register_i(  "trk", "event" );
    register_i(  "trk", "index" );
    register_f(  "trk", "pt"    );
    register_f(  "trk", "eta"   );
    register_f(  "trk", "phi"   );
    register_f(  "trk", "vx"    );
    register_f(  "trk", "vy"    );
    register_f(  "trk", "vz"    );
    register_f(  "trk", "d0"    );
    register_f(  "trk", "z0"    );
    register_f(  "trk", "Chi2"  );
    register_i(  "trk", "Ndof"  );
    register_i(  "trk", "numberOfValidHits"  );
    register_i(  "trk", "numberOfLostHits"  );
    register_b(  "trk", "isFPHitTrk" );
    register_i(  "trk", "expInnHits" );
    register_f(  "trk", "impEta" );
    register_f(  "trk", "impPhi" );
  }

  {
    register_t(  "htrk"  , "trackRecHits"      );     
    register_i(   "htrk", "nTHits"          );
    register_vi(  "htrk", "tid"   );
    register_vf(  "htrk", "trho"    );
    register_vf(  "htrk", "teta"    );
    register_vf(  "htrk", "tphi"    );
  }

  {
    // calo towers
    run_name_.push_back( "caloTower" );
    register_t( "ct", "caloTowers"  );      
    register_i( "ct", "run");
    register_i( "ct", "event");
    register_i( "ct", "index");
    register_f( "ct", "energy");
    register_i( "ct", "ieta");
    register_i( "ct", "iphi");
    register_f( "ct", "pt");
    register_f( "ct", "eta");
    register_f( "ct", "phi");
    register_f( "ct", "emEt");
    register_f( "ct", "hadEt");
    register_f( "ct", "outerEt");
  }

  {
    // jets
    register_t( "j"    , "jets"        );      
    register_i( "j","run");
    register_i( "j","event");
    register_i( "j","index");
    register_f( "j","pt");
    register_f( "j","eta");
    register_f( "j","phi");
    register_f( "j","charge");
    register_i( "j","partonFlavour");
    register_f( "j","corrFactor");
    register_f( "j","gp_pt");
    register_f( "j","gp_eta");
    register_f( "j","gp_phi");
    register_f( "j","gj_pt");
    register_f( "j","gj_eta");
    register_f( "j","gj_phi");

    register_f( "j","btagTkC");
    register_f( "j","btagSoftM");
    register_f( "j","btagJPB");
    register_f( "j","btagSVtx");
    register_f( "j","btagCSVBT");
    register_f( "j","btagCSVMVA");

    register_b( "j","mvaIdLoose");
    register_b( "j","mvaIdMedium");
    register_b( "j","mvaIdTight");
    register_f( "j","mvaJetId");
    
    register_vi( "j","indexPfConstHPt");
    register_vi( "j","indexPfConstLPt");
    
  }

  for( size_t v=0; v<caloMET.size(); v++ )
    {
      // calo met
      register_t( "cmet", v, "caloMET", caloMET[v] );      
      run_name_.push_back( getName( "cmet", v ) );      
      register_i( "cmet", "run"   , v );
      register_i( "cmet", "event" , v );
      register_f( "cmet", "pt"    , v );
      register_f( "cmet", "phi"   , v );
      register_f( "cmet", "sumEt" , v );
      register_f( "cmet", "mEtSig", v );
    }
  
  for( size_t v=0; v<recoMET.size(); v++ )
    {
      // reco met
      register_t( "rmet", v, "recoMET", recoMET[v] );      
      run_name_.push_back( getName( "rmet", v ) );      
      register_i( "rmet", "run"   , v );
      register_i( "rmet", "event" , v );
      register_f( "rmet", "pt"    , v );
      register_f( "rmet", "phi"   , v );
      register_f( "rmet", "sumEt" , v );
      register_f( "rmet", "mEtSig", v );
    }
  
  for( size_t v=0; v<pfMET.size(); v++ )
    {
      // pf met
      register_t( "pfmet", v, "pfMET", pfMET[v] );      
      run_name_.push_back( getName( "pfmet", v ) );      
      register_i( "pfmet", "run"   , v );
      register_i( "pfmet", "event" , v );
      register_f( "pfmet", "pt"    , v );
      register_f( "pfmet", "phi"   , v );
      register_f( "pfmet", "sumEt" , v );
      register_f( "pfmet", "mEtSig", v );
      register_f( "pfmet", "cEmEF", v );
      register_f( "pfmet", "cHadEF", v );
      register_f( "pfmet", "cMuEF", v );
      register_f( "pfmet", "nEmEF", v );
      register_f( "pfmet", "nHadEF", v );
    }

  for( size_t v=0; v<patMET.size(); v++ )
    {
      // pat met
      register_t( "patmet", v, "patMET", patMET[v] );      
      run_name_.push_back( getName( "patmet", v ) );      
      register_i( "patmet", "run"   , v );
      register_i( "patmet", "event" , v );
      register_f( "patmet", "pt"    , v );
      register_f( "patmet", "phi"   , v );
      register_f( "patmet", "sumEt" , v );
      register_f( "patmet", "mEtSig", v );
      register_f( "patmet", "cEmEF", v );
      register_f( "patmet", "cHadEF", v );
      register_f( "patmet", "cMuEF", v );
      register_f( "patmet", "nEmEF", v );
      register_f( "patmet", "nHadEF", v );
      register_f( "patmet", "emEBF"  , v );
      register_f( "patmet", "emEEF"  , v );
      register_f( "patmet", "emHFF"  , v );
      register_f( "patmet", "hadHBF" , v );
      register_f( "patmet", "hadHEF" , v );
      register_f( "patmet", "hadHFF" , v );
      register_f( "patmet", "hadHOF" , v );
      register_f( "patmet", "hadF"   , v );
    }
  
   for( size_t v=0; v<caloJets.size(); v++ )
    {
      // calo jets
      register_t( "cj", v, "caloJets", caloJets[v] );      
      run_name_.push_back( getName( "cj", v ) );      
      register_i( "cj", "run"  , v );
      register_i( "cj", "event"  , v );
      register_i( "cj", "index"  , v );
      register_f( "cj", "pt"  , v );
      register_f( "cj", "eta"  , v );
      register_f( "cj", "phi"  , v );
      register_f( "cj", "emFraction"  , v );
      register_f( "cj", "towersArea"  , v );
      register_i( "cj", "nCT"  , v );
      register_vi( "cj", "indexCT"  , v );
      register_vf( "cj", "ptCT"  , v );
      register_vf( "cj", "etaCT"  , v );
      register_vf( "cj", "phiCT"  , v );
    }

   for( size_t v=0; v<pfJets.size(); v++ )
    {
      // pf jets
      register_t( "pfj", v, "pfJets", pfJets[v] );      
      run_name_.push_back( getName( "pfj", v ) );      
      register_i( "pfj", "run" , v );
      register_i( "pfj", "event" , v );
      register_i( "pfj", "index" , v );
      register_f( "pfj", "pt" , v );
      register_f( "pfj", "eta" , v );
      register_f( "pfj", "phi" , v );
      register_f( "pfj", "maxD" , v );
      register_f( "pfj", "etaEtaMom" , v );
      register_f( "pfj", "etaPhiMom" , v );
      register_f( "pfj", "phiPhiMom" , v );
      register_f( "pfj", "area" , v );
      register_f( "pfj", "chargedHadronEnergyFraction" , v );
      register_f( "pfj", "neutralHadronEnergyFraction" , v );
      register_f( "pfj", "chargedEmEnergyFraction" , v );
      register_f( "pfj", "chargedMuEnergyFraction" , v );
      register_f( "pfj", "neutralEmEnergyFraction" , v );
      register_f( "pfj", "chargedMultiplicity" , v );
      register_f( "pfj", "neutralMultiplicity" , v );
      register_f( "pfj", "muonMultiplicity" , v );
      register_f( "pfj", "jetResolution" , v );

      /* register_i( "pfj", "nJC"  , v ); //MM discarded for the moment
      register_vi("pfj", "indexJC"  , v );
      register_vf("pfj", "ptJC"   , v );
      register_vf("pfj", "etaJC"  , v );
      register_vf("pfj", "phiJC"  , v );*/
    }

  for( size_t v=0; v<trackJets.size(); v++ )
    {
      // track jets
      register_t( "trkj", v, "trackJets", trackJets[v] );      
      run_name_.push_back( getName( "trkj", v ) );      
      register_i( "trkj", "run" , v );
      register_i( "trkj", "event" , v );
      register_i( "trkj", "index" , v );
      register_f( "trkj", "pt" , v );
      register_f( "trkj", "eta" , v );
      register_f( "trkj", "phi" , v );
      register_f( "trkj", "maxD" , v );
    }

  // Z -> ll
  for( int iZll=0; iZll<2; iZll++ )
    //  int iZll=0; 
    {
      string prefix = "zee";
      if( iZll==1 ) prefix = "zmm";
      register_t(   prefix, Zmode[iZll] );      
      register_i(   prefix, "run" );
      register_i(   prefix, "event" );
      register_f(   prefix, "mass" );
      register_f(   prefix, "pt" );
      register_f(   prefix, "theta" );
      register_f(   prefix, "phi" );
      register_f(   prefix, "rapidity" );
      register_f(   prefix, "helicity" );
      register_i(   prefix, "dau_0" );
      register_f(   prefix, "pt_0" );
      register_f(   prefix, "eta_0" );
      register_f(   prefix, "phi_0" );
      register_i(   prefix, "dau_1" );
      register_f(   prefix, "pt_1" );
      register_f(   prefix, "eta_1" );
      register_f(   prefix, "phi_1" );
    }
 
  string treeName_("eventSummary");
  map< string, pair< string, string > >::const_iterator it;
  for( it=treeNames.begin(); it!=treeNames.end(); it++ )
    {
      string prefix_ = it->first;
      if( prefix_=="run" || prefix_=="evt" ) continue;
      string str1_ = it->second.first; 
      string str2_ = it->second.second;
      string str_ = concatanate( str1_, str2_ );
      evt_n[prefix_]=0; 
      _treeMgr->branchio<int>( treeName_.c_str(), ("n_"+str_).c_str(), 
			       &evt_n[prefix_]  );
      evt_first[prefix_].resize(0);
      evt_first[prefix_].push_back(0);
    }

  // read the eventSummary ntuple to assign first element pointers
  Long64_t  nentries_ = _treeMgr->t(treeName_.c_str())->GetEntriesFast();
  for( unsigned jj=0; jj<nentries_; jj++ )
    {
      _treeMgr->entry( jj, treeName_.c_str() );
      for( map< string, vector<Int_t> >::const_iterator it=evt_first.begin();
	   it!=evt_first.end(); ++it ) 
	{
	  string prefix_ = it->first;
	  Int_t ifirst = evt_first[prefix_][jj];
	  evt_first[prefix_].push_back( ifirst+evt_n[prefix_] );
	}
    }  

  register_t( "run", "runSummary" );
  register_i( "run", "analyzed"   );
  register_i( "run", "selected"   );
  register_i( "run", "nEventInJob" );
  for( size_t ii=0; ii<run_name_.size(); ii++ )
    register_s( "run", run_name_[ii].c_str() );
  register_f( "run","intXSec" );
  register_f( "run","intXSecErr" );
  register_f( "run","extXSecLO" );
  register_f( "run","extXSecLOErr" );
  register_f( "run","extXSecNLO" );
  register_f( "run","extXSecNLOErr" );
  register_f( "run","filterEfficiency" );    
 
  getEntry( "run", 0 );
  // HLT Trigger lines
  hlt.clear();
  bool first__=true;

  {
  
    string triggerName_;
    triggerName_ = get_s("run","hltTriggerNames");
  

    istringstream iss;
    iss.str( triggerName_ );
    size_t ii=0;
    string name_;
    //  while( name_!="HLTriggerFinalPath")
    do
      { 
	iss >> name_;
	if( !iss.good() ) break;
	hlt[name_] = ii++;
	cout << ii << "\t-->\t" << name_ << endl;
      }
    //while( name_!="endjob_step"); // for MC events
    //while( name_!="HLTriggerFinalPath"); // for first collisions
    while( 1 );
  }

  {
    string eAccName_( get_s("run","electronFiduFs") );
    istringstream iss;
    iss.str( eAccName_ );
    //    size_t ii=0;    
    string name_;
    eAcc.clear();
    do     
      { 
	iss >> name_;
	if( !iss.good() ) break;
	eAcc.push_back(name_);
	//	cout << ++ii << "\t-->\t" << name_ << endl;
      }
    while( 1 );
  } 

  {
    string eidName_( get_s("run","electronIDs") );
    istringstream iss;
    iss.str( eidName_ );
    //    size_t ii=0;    
    string name_;
    eid.clear();
    do     
      { 
	iss >> name_;
	if( !iss.good() ) break;
	eid.push_back(name_);
	//	cout << ++ii << "\t-->\t" << name_ << endl;
      }
    while( 1 );
  } 
  
  {
    string muidName_( get_s("run","muonIDs") );
    istringstream iss;
    iss.str( muidName_ );
    string name_;
    muid.clear();
    do     
      { 
	iss >> name_;
	if( !iss.good() ) break;
	muid.push_back(name_);
	//	cout << ++ii << "\t-->\t" << name_ << endl;
      }
    while( 1 );
  }  

  {
    string phidName_( get_s("run","photonIDs") );
    istringstream iss;
    iss.str( phidName_ );
    //    size_t ii=0;
    string name_;
    phid.clear();
    do     
      { 
	iss >> name_;
	if( !iss.good() ) break;
	phid.push_back(name_);
	//	cout << ++ii << "\t-->\t" << name_ << endl;
      }
    while( 1 );
  }  

  {
    string tauidName_( get_s("run","tauIDs") );
    istringstream iss;
    iss.str( tauidName_ );
    string name_;
    tauid.clear();
    do     
      { 
	iss >> name_;
	if( !iss.good() ) break;
	tauid.push_back(name_);
	//	cout << ++ii << "\t-->\t" << name_ << endl;
      }
    while( 1 );
  }  

  if( first__ && verbosity>=2 ) 
    {
      string fmt_("[%-20s]\t=\t%-s\n");
      for( size_t ii=0; ii<run_name_.size(); ii++ )
	{
	  cout << ii << " " << run_name_[ii] << endl;
	  printf( fmt_.c_str(), run_name_[ii].c_str(),
		  get_s("run",run_name_[ii].c_str() ) );
	}

      cout << "Selected/Analyzed=" 
	   << get_i("run","analyzed") << "/" 
	   << get_i("run","selected") << " over "
	   << get_i("run","nEventInJob") <<endl;  
      cout << "intXSec/extXSecLO/extXSecNLO/Filter-Eff " 
	   << get_f("run","intXSec")   <<"+/-"<<get_f("run","intXSecErr")   << "/" 
	   << get_f("run","extXSecLO") <<"+/-"<<get_f("run","extXSecLOErr") << "/" 
	   << get_f("run","extXSecNLO")<<"+/-"<<get_f("run","extXSecNLOErr")<< "/" 
	   << get_f("run","filterEfficiency") << "%" << endl;      
      // MM
      // cout << "Trigger lines:" << endl;
      // for( map< string, size_t >::const_iterator it_=hlt.begin();
      // 	   it_!=hlt.end(); it_++ )
      // 	{
      // 	  string hltLine = it_->first;
      // 	  cout << hltLine << endl;
      // 	}
    
    }
  if( first__ ) {
    //store the number of event processed in the job (includes filter efficiency)
    nEventProc = get_si("run","nEventInJob");
    
    first__=false;
  }

  return true;
}

bool
EventServer::nextEvent()
{
  if( _hasEventList )
    {
      if( _eventListIt==_eventList.end() ) return false;
      bool ok_;
      ok_ = eventInFile( (*_eventListIt).first, (*_eventListIt).second );
      _eventListIt++;
      return ok_;
    }
  int ievt_ = eventNumber+1;
  return event(ievt_);
}

bool
EventServer::previousEvent()
{
  int ievt_ = eventNumber-1;
  return event(ievt_);
}

bool
EventServer::event( int ievt_ )
{

  if( ievt_==0 ) 
    {
      cout << "Must be strictly positive " << endl;      
      return false;
    }

  // MM
  nEventProc = get_si("run","nEventInJob");

  // check if event is in the same file
  //  bool sameFile_=true;
  int fileNum = ifile;
  if( ievt_>_nevent[fileNum] )
    {
      //      sameFile_=false;
      size_t ii=fileNum+1;
      if( ii>=_filenames.size() )
	{
	  ievt_--;
	  cout << "Total number of events reached --> " << ievt_ << endl;
	  return false;
	}
      for( ; ii<_nevent.size(); ii++ )
	{
	  if( ievt_>_nevent[ii] ) continue;
	  break;
	}
      if( ii>=_nevent.size() ) return false;
      fileNum = ii;
     
    }
  else if( fileNum>0 && ievt_<=_nevent[fileNum-1] )
    {
      //      sameFile_=false;
      size_t ii=fileNum-1;
      for( ; ii>=0; ii-- )
	{
	  if( ievt_<=_nevent[ii] ) continue;
	  break;
	}
      if( ii<0 ) return false;
      fileNum = ii;
      
    }
    
  eventNumber = ievt_;
  int evtNum = ievt_;
  if( fileNum>0 )
    evtNum -= _nevent[fileNum-1];

  return eventInFile( fileNum, evtNum );
}

bool
EventServer::eventInFile( int fileNum, int evtNum )
{
  if( fileNum!=ifile )
    {
      ifile = fileNum;
      assert( openFile( _filenames[ifile] ) );
    }
  ievent = evtNum;
  assert( ievent<=_nevent[ifile] );
  eventNumber = ievent;
  if( ifile>0 )  eventNumber += _nevent[ifile-1];

  unsigned jj = ievent-1;

  getEntry("evt",jj);
  curRun      = get_i("evt","run");
  curLumi     = get_i("evt","lumiBlock");
  curEvent    = get_i("evt","event");
  time        = get_i("evt","time");
  hasPrimaryVtx = false;

  ptHat     = get_f("evt","ptHat");
  processID = get_i("evt","processID");
  
  rhoFJ     = get_f("evt","rhoFastJet");

  trueNI    = get_f("evt","trueNI");

  nIBXm1    = get_i("evt","nInteracBXm1");
  nIBX      = get_i("evt","nInteracBX");
  nIBXp1    = get_i("evt","nInteracBXp1");

  // hlt
  getEntry("hlt",jj); 
  
  _hltBits.Clear();

  bool hltBool = get( "hlt", "hltBits", _hltBits );
  assert( hltBool );
 

  size_t nbits_ = _hltBits.GetNbits();
  if( nbits_!=hlt.size() )
    {
      //***  
      //      cout << "hlt_bits/hlt " << nbits_ << "/" << hlt.size() << endl;
      //      assert(0);
    }
 
  if( verbosity>=2 ) 
    {
      cout << "Fired HLT lines" << endl;
      for( map< string, size_t >::const_iterator it_=hlt.begin();
	   it_!=hlt.end(); it_++ )
	{
	  if( isFired( it_->first ) ) 
	    {
	      cout << "==> " << it_->first << endl;
	    }
	}
    }

  _eidBits.Clear();
  _muidBits.Clear();
  _tauidBits.Clear();
  _phidBits.Clear();
  _eAccBits.Clear();
  
  return true;
}

float
EventServer::getRFJ() {
  return rhoFJ;
}


float
EventServer::getTrueNInt() {
  return trueNI;
}

int 
EventServer::getnInteract( string str ) {

  if(str == "BXm1" || str == "m1" ) 
    return nIBXm1;
  else if(str == "BXp1" || str == "p1" ) 
    return nIBXp1;
  else
    return nIBX;

}


bool
EventServer::isFired( const string & hltLine, bool noVersion )
{
  // get the vector index
  if( noVersion ) {
    for(itHlt = hlt.begin();itHlt != hlt.end(); itHlt++) {
      if( (*itHlt).first.find(hltLine) == 0  ) {
	if((*itHlt).first.substr(hltLine.size(),2)=="_v") {
	  return _hltBits[ (*itHlt).second ];
	}
      }
    }
    return false;
  }
  else {
    if( hlt.count( hltLine )==0 ) return false;
    return _hltBits[ hlt[ hltLine ] ];
  }
}

void
EventServer::firedLines( map< int, string >& lines )
{
  lines.clear();

  // get the vector index
  for( itHlt = hlt.begin();itHlt != hlt.end(); itHlt++ ) 
    {
      size_t ibit = itHlt->second;
      if( !_hltBits[ ibit ] ) continue;
      lines[ ibit ] = itHlt->first; 
    }
  return;
}

EventServer::~EventServer()
{
  closeFile();
}

bool
EventServer::load( const string & prefix, int ii, int v )
{
  int n_   = n( prefix, v );

  if( n_==0 ) 
    {
      return false;
    }
  if( ii<0 || ii>=n_ )
    {
      cout << "There are only " << n_ << " " 
	   << getName( prefix, v ) << " in the event " << ii << endl;
      return false;
    }

  Long64_t jentry = n0(prefix,v) + ii;

  getEntry( prefix, jentry, v );

  if( prefix!="hlt" && prefix!="htrk" )
    {
      if( get_i(prefix,"run",v)!=curRun || get_i(prefix,"event",v)!=curEvent )
	return false;
    }
  if( prefix=="vtx" )
    {
      hasPrimaryVtx = true;
    }
  else if( prefix=="el" )
    {
      if( get_i(prefix,"index",v)!=ii ) return false;
      _eidBits.Clear();
      assert( get( "el", "electronIDs", _eidBits ) );
      _eAccBits.Clear();
      assert( get( "el", "elAccBits", _eAccBits ) );
    }
  else if( prefix=="mu" )
    {      
      _muidBits.Clear();
      assert( get( "mu", "muonIDs", _muidBits ) );
    }
  else if( prefix=="ph" )
    {      
      _phidBits.Clear();
      assert( get( "ph", "photonIDs", _phidBits ) );
    }

  return true;
}

void 
EventServer::oneLine( ostream& o )
{ 
  o << "Event " << eventNumber;
  o << " File/Event " << _filenames[ifile] << "/" << ievent;
  o << " " << get_i("evt","run") << "/" << get_i("evt","event");
  o << endl;
} 

//TTree*
//EventServer::getTree( string prefix )
//{
// GHMNEW
//  assert( ntuNames.count(prefix)==1 );
//  string str_; 
//  if( prefix!="run" && prefix!="evt" ) str_ += "ntu_";
//  str_ += ntuNames[prefix];
//  return _treeMgr->t( str_.c_str() );
//}

void 
EventServer::setPrefix( const string & prefix, const string & ntuName )
{
  //   if( ntuNames.count( prefix )!=0 )
  //     {
  //       cout << "Prefix already used " << prefix << " -- abort " << endl;
  //       assert(0);
  //     }
  //   ntuNames.insert( make_pair( prefix , ntuName ) );
}

// GHMNEW
string 
EventServer::concatanate( const string & str1, const string & str2, bool sep )
{
  string str_(str1); 
  if( str2!="" )
    {
      if( sep ) str_ += "_";
      str_ += str2;
    }
  return str_;
}

string
EventServer::getName( const string & prefix, int v, bool sep )
{
  string key_ = getKey( prefix, v );
  if( treeNames.count( key_ )!=1 )
    {
      assert(0);
    }
  return concatanate( treeNames[key_].first, treeNames[key_].second, sep );
}

string
EventServer::getTreeName( const string & prefix, int v )
{
  string str_; 
  if( prefix!="run" && prefix!="evt" ) str_ += "ntu_";
  str_ += getName( prefix, v );
  return str_;
}

void
EventServer::register_t( const string & prefix, const string & str1 )
{
  register_t( prefix, -1, str1, "" );  
}

void
EventServer::register_t( const string & prefix, int v, const string & str1, const string & str2 )
{
  string key_ = getKey( prefix, v );
  //cout<<key_<<endl;
  treeNames.insert( make_pair( key_, make_pair( str1, str2 ) ) );  
}

string
EventServer::getKey( const string & prefix, int v )
{
  if( v<0 ) return prefix;
  ostringstream oss;
  oss << prefix << "_" << v;
  return oss.str();
}

void
EventServer::register_i( const string & prefix, const string & varName, int v )
{
  string treeName = getTreeName(prefix, v );
  string str = getKey( prefix, v );
  str += "_";
  str += varName;
  //cout<<str<<"   "<<map_i.count(str)<<endl;
  assert( map_i.count(str)==0 );
  _treeMgr->branchio<int>( treeName.c_str(), varName.c_str(), &map_i[str] ); 
}

void
EventServer::register_f( const string & prefix, const string & varName, int v )
{
  string treeName = getTreeName(prefix, v );
  string str = getKey( prefix, v );
  str += "_";
  str += varName;
  assert( map_f.count(str)==0 );
  _treeMgr->branchio<float>( treeName.c_str(), varName.c_str(), &map_f[str] ); 
}

void
EventServer::register_s( const string & prefix, const string & varName, int v )
{
  string treeName = getTreeName(prefix, v );
  string str = getKey( prefix, v );
  str += "_";
  str += varName;
  //cout<<str<<"   "<<map_s.count(str)<<endl;
  assert( map_s.count(str)==0 );
  _treeMgr->branchio< string* >( treeName.c_str(), varName.c_str(), &map_s[str] ); 
}

void
EventServer::register_b( const string & prefix, const string & varName, int v )
{
  string treeName = getTreeName(prefix, v );
  string str = getKey( prefix, v );
  str += "_";
  str += varName;
  assert( map_b.count(str)==0 );
  map_b[str] = NULL; //MM FIXME
  _treeMgr->branchio< TBits* >( treeName.c_str(), varName.c_str(), &map_b[str] ); 
}

void
EventServer::register_vi( const string & prefix, const string & varName, int v )
{
  string treeName = getTreeName(prefix, v );
  string str = getKey( prefix, v );
  str += "_";
  str += varName;
  assert( map_vi.count(str)==0 );
  _treeMgr->branchio< vector<int>* >( treeName.c_str(), varName.c_str(), &map_vi[str] ); 
}
  
void
EventServer::register_vf( const string & prefix, const string & varName, int v )
{
  string treeName = getTreeName(prefix, v );
  string str = getKey( prefix, v );
  str += "_";
  str += varName;
  assert( map_vf.count(str)==0 );
  _treeMgr->branchio< vector<float>* >( treeName.c_str(), varName.c_str(), &map_vf[str] ); 
}
  
void
EventServer::getEntry( const string & prefix, unsigned jj, int v )
{
  string treeName = getTreeName(prefix,v);
  _treeMgr->entry( jj, treeName.c_str() );
}


bool
EventServer::get( const string & prefix, const string & varName, float& val, int v )
{  
  string str = getKey( prefix, v );
  str += "_";
  str += varName;  
  if( map_f.count(str)!=1 ) return false;
  val = map_f[str];
  return true;
}

bool
EventServer::get( const string & prefix, const string & varName, int& val, int v )
{  
  string str = getKey( prefix, v );
  str += "_";
  str += varName;    
  if( map_i.count(str)!=1 ) return false;
  val = map_i[str];
  return true;
}

bool
EventServer::get( const string & prefix, const string & varName, string & val, int v )
{  
  string str = getKey( prefix, v );
  str += "_";
  str += varName;    
  if( map_s.count(str)!=1 ) return false;
  val = *map_s[str];
  return true;
}

bool
EventServer::get( const string & prefix, const string & varName, TBits& val, int v )
{  
  string str = getKey( prefix, v );
  str += "_";
  str += varName;    
  if( map_b.count(str)!=1 || map_b[str]==NULL ) return false;  //MM FIXME
  val = *map_b[str];
  return true;
}

int 
EventServer::get_i( const string & prefix, const string & varName, int v )
{  
  string str = getKey( prefix, v );
  str += "_";
  str += varName;
  if( map_i.count(str)!=1 ) { cout << str << endl; assert(0); }
  return map_i[str];
}


int 
EventServer::get_si( const string & prefix, const string & varName, int v )
{  
  string str = getKey( prefix, v );
  str += "_";
  str += varName;
  if( map_i.count(str)!=1 ) { cout << str << endl; assert(0); }
 
  int sum = 0;
  for( size_t ii=0; ii<_treeMgr->entries( "runSummary" ); ii++ ) {
    getEntry( "run", ii );
    sum += map_i[str];
  }
  getEntry( "run", 0 ); //reinit
  return sum;
}


float
EventServer::get_f( const string & prefix, const string & varName, int v )
{  

  // for(map< string, float >::const_iterator it=map_f.begin(); it!=map_f.end();it++)
//     cout<<" --> "<<it->first<<endl;
//   abort();
  string str = getKey( prefix, v );
  str += "_";
  str += varName;  
  if( map_f.count(str)!=1 ) { cout << str << endl; assert(0); }
  return map_f[str];
}

const char*
EventServer::get_s( const string & prefix, const string & varName, int v )
{  
  string str = getKey( prefix, v );
  str += "_";
  str += varName;  
  if( map_s.count(str)!=1 ) { cout << str << endl; assert(0); }
  return map_s[str]->c_str();
}

TBits*
EventServer::get_b( const string & prefix, const string & varName, int v )
{  
  string str = getKey( prefix, v );
  str += "_";
  str += varName;  
  if( map_b.count(str)!=1 ) { cout << str << endl; assert(0); }
  return map_b[str];
}

vector<int>&
EventServer::get_vi( const string & prefix, const string & varName, int v )
{  
  string str = getKey( prefix, v );
  str += "_";
  str += varName;  
  if( map_vi.count(str)!=1 ) { cout << str << endl; assert(0); }
  return *map_vi[str];
}

int
EventServer::get_vi( const string & prefix, const string & varName, size_t ii, int v )
{  
  return get_vi(prefix,varName,v)[ii];
}

vector<float>&
EventServer::get_vf( const string & prefix, const string & varName, int v )
{  
  string str = getKey( prefix, v );
  str += "_";
  str += varName;  
  if( map_vf.count(str)!=1 ) { cout << str << endl; assert(0); }
  return *map_vf[str];
}

float
EventServer::get_vf( const string & prefix, const string & varName, size_t ii, int v )
{  
  return get_vf(prefix,varName,v)[ii];
}

int 
EventServer::n( const string & prefix, int v )
{
  string str = getKey( prefix,  v );
  return evt_n[str]; 
}

int 
EventServer::n0( const string & prefix, int v )
{
  string str = getKey( prefix, v );
  return evt_first[str][ievent-1];
}
