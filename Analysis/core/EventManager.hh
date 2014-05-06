#ifndef EventManager_h
#define EventManager_h

#include <map>
#include <vector>
#include <string>

#include <TString.h>
#include <TH1.h>

#include "Analysis/core/EventServer.hh"
#include "Analysis/core/Candidate.hh"

#include "Analysis/utils/Config.hh"

class CandAssoc
{
public:
  CandAssoc() {}
  void insert( Candidate* c1, Candidate* c2 );
  ~CandAssoc() { clear(); }
  Candidate* getSecond( Candidate* c ) const;
  Candidate* getFirst( Candidate* c ) const;
  void clear();
private:
  std::map< Candidate*, Candidate* > _map1;
  std::map< Candidate*, Candidate* > _map2;
};

class EventManager
{
public:

  // base candidate types
  enum { kTrack=0, 
	 kGsfTrack,
	 kPhoton,
	 kElectron,
	 kMuon,
	 kTau,
	 kCaloTower,
	 kEcalSC,
	 kEcalBC,
	 kEcalPSC,
	 kEcalPfSC,
	 kEcalPfBC,
	 kPfCand,
	 kPfCand_lowPt,
         kMcTruthCandidate, 
         kGenEventCandidate, 
	 kNType };
  static std::string objectName[kNType];

  // jets
  enum { kPatJet=0, 
	 kCaloJet, 
	 kTrkJet, 
	 kPfJet, 
	 kGenJet,
	 kNJet }; 
  static std::string jetName[kNJet];

  // met
  enum { kPatMet=0, 
	 kCaloMet, 
	 kTrkMet, 
	 kPfMet, 
	 kGenMet,
	 kNMet }; 
  static std::string metName[kNMet];

  // access to THE EventManager (created beforehand by SampleAnalysis)
  static EventManager* e();

  // access to the undelying event server
  virtual EventServer& a(); 
  const EventServer& a() const { return *_a; }

  // move to next event and load it
  virtual bool nextEvent();
  virtual bool goToEvent( int ievt ); 

  virtual bool findEvent( int ievt, int irun ); 
  
  virtual bool isSelected();

  // access to global event quantities

  int run()          const { return a().curRun;      } 
  int lumisec()      const { return a().curLumi;     } 
  int timestamp()    const { return a().time;        } 
  int event()        const { return a().curEvent;    } 
  int eventNumber()  const { return a().eventNumber; } 
  int file()         const { return a().ifile;       } 
  int eventInFile()  const { return a().ievent;      } 
  TString fileName() const { return a().fName;       } 
 
  float ptHat()      const { return a().ptHat;       }

  //run and job quantities
  int nEventProc()   const { return a().nEventProc;  }

  // access to lists of base candidates
  const CandList& tracks()        const { return candList( kTrack       );  }
  const CandList& gsfTracks()     const { return candList( kGsfTrack    );  }
  const CandList& electrons()     const { return candList( kElectron    );  }
  const CandList& muons()         const { return candList( kMuon        );  }
  const CandList& taus()          const { return candList( kTau         );  }
  const CandList& photons()       const { return candList( kPhoton      );  }
  const CandList& caloTowers()    const { return candList( kCaloTower   );  }
  const CandList& pfCandidates()  const { return candList( kPfCand      );  }
  const CandList& uncPfCandidates();
  const CandList& candList( int ) const;

  // primary vertex
  void setPrimaryVertex( Vertex* vtx=0 ); 
  bool setPrimaryVertex( int vtxid ); 
  const Vertex* getPrimaryVertex() const { return _primaryVertex; };
  const Vertex& primaryVertex() const { return *_primaryVertex; }
  const VtxList& vertices() const { return _vertices; }
  int nGoodVertices() const { return _nGoodVertices; }
  int iPrimaryVertex() const { return _iPrimaryVertex; } // index of primary vertex

  //Vertex finding
  Vertex* findVertex( float eta, float phi );

  // access to lists of jets 
  const CandList& jetList( int type=kPatJet, int v=-1 ) const;
  const CandList& jetList( int type, std::string& str ) const;
  const CandList& jetList( const std::string &, int v=-1 ) const;
  const CandList& caloJets()      const { return jetList( kCaloJet,-1 );  }

  // access to missing transverse momentum 
  const Candidate* met( int type=kPatMet, int v=-1 ) const;
  const Candidate* met( const std::string &, int v=-1 ) const;
  const Candidate* met(const std::string& str, const std::string& name ) const;
  const Candidate* cleanPUmet( Vertex* _vtx, Candidate* bosC, float rF );
  
  //Refresh pfCandidates
  void refreshPfCandidates();
  void reloadPfCandidates( int v );

  const Candidate* caloMet()      const { return met( kCaloMet, -1 ); }

  // access to lists of pre-computed of composite candidates
  const CandList& compositeCandList( const std::string & ) const;

  // histograms with Et of objects summed in 72 bins in phi
  TH1* sumOfEt( std::string , int jj=-1 ) const;

  // HLT lines 
  bool isFired( const std::string & hltLine, bool noVersion=false ) const;  
  void firedLines( std::map< int, std::string >& lines ) const;

  //Rho Fast jet
  float getRhoFastJet();

  //Ninteraction
  int getnInteract( std::string str="");
  float getTrueNInt();

  // creation and access to user-defined candidate lists
  CandList& userCandList( const std::string & );

  // access to "iso" maps ("iso"=light object with eta, phi, val, dr)
  const CandIdIsoMap& isoMap( const std::string & ) const;

  // access to Ecal maps ( sc->cluster candidates)
  //GHMconst CandIdEcalMap& ecalMap( std::string ) const;

  const CandIdRecHitMap& ecalRecHitMap( const std::string & ) const;
  const CandIdPSRecHitMap& ecalPSRecHitMap( const std::string & ) const;

  // access to candidate -> candidate map (lepton -> track, sc )
  const CandMap& candMap( const std::string & ) const;

  // two way access candidate <-> candidate map (lepton -> sc )
  const CandAssoc& candAssoc( const std::string & ) const;

  // get the partner candidate in a one-to-one association
  Candidate* getFirst ( const std::string &, Candidate* ) const;
  Candidate* getSecond( const std::string &, Candidate* ) const;

  // direct access to candidates in lists
  Candidate* getCandidate( const std::string& , size_t i );

  // get the decay tree from the generated event
  Candidate*       decayTreeOld();
  Candidate*       decayTree( bool dbg=false );
  const Candidate* decayTree() const;

  // get the true event category
  std::string categ()      const;
  std::string categ_long() const;

  // true is all primary leptons (from Z and W) are within acceptance, 
  //  for electrons: pt>10. and |eta|<2.5
  //  for muons:     pt>10. and |eta|<2.4
  bool inAcceptance() const;

  // sum of neutrino four momenta
  const TLorentzVector& neutrinosFromBosonsP4();
  const TLorentzVector& neutrinosFromTausP4();

  // Lepton MC matching
  const CandList& mcTruthCandidates() const 
  { return candList( kMcTruthCandidate   );  }
  const CandList& mcGenEventCandidates() const 
  { return candList( kGenEventCandidate   );  }
  const Candidate* mcCand( const Candidate& lepton ) const;  

  // printing
  void print( std::string what="all" );

protected:
    
  // init, refresh
  void init();
  void refresh();

  // load event
  void loadEvent();

  // load base candidates
  void loadTracks();
  void loadGsfTracks();
  void loadElectrons();

  void loadElectronSCMaps();
  void loadElectronIsoDepMaps();
  void loadElectronTrackMaps();
  void loadElectronPfCandMaps();
  void loadEcalSCBCMaps();
  void loadEcalBCPSCMaps();
  void loadEcalBCRecHitMaps();
  void loadEcalPSCPSRecHitMaps();
  void loadMuonTrackMaps();
  void loadMuonPfCandMaps();
  void loadTrackHitMaps();
  void loadJetPfCandMaps();
  
  // jet cleaning
  void loadJetRelatedMaps();

  void loadMuons();
  void loadTaus();
  void loadPhotons();
  void loadEcalBC();
  void loadEcalSC();
  void loadEcalPSC();
  void loadEcalPfBC();
  void loadEcalPfSC();
  void loadEcalRecHits();
  void loadEcalPreshowerRecHits();
  void loadCaloTowers();
  void loadPfCandidates( int v=-1 );
  void loadJets( int, int v );
  void loadMet( int, int v );
  void loadMcTruthCandidates();
  void loadZllCandidates();

  void loadTrackHits();
  
  // MC matching (for leptons so far)
  void mcMatching();

  //Vertex matchin (for any candidate)
  Vertex* vertexMatching(Vertex* vtx);

  // lists of clean jets at the primary vertex
  void setListsOfCleanJets(); 

  // maps of clean jets and clean leptons, per vertex index
  std::map< int, CandList > _cleanJets;
  std::map< int, CandList > _leptonJets;
  std::map< int, CandList > _cleanLeptons;

  // list of rechits
  ecalRecHitList _ecalRecHits;
  ecalPSRecHitList _ecalPSRecHits;

  // lists of lists of base candidates
  std::vector< CandList > _candList;

  // lists of lists of jet candidates
  //std::vector< CandList > _jetList;
  std::map< std::pair<int, int >, CandList > _jetList;

  // lists of met candidates
  std::map< std::pair< int, int >, Candidate* > _met;

  // map of lists of composite candidates
  std::map< std::string, CandList > _compositeCandList;

  // map of lists of user defined candidates
  std::map< std::string, CandList > _userCandList;

  // list of vertices
  VtxList _vertices;

  // primary vertex 
  Vertex* _primaryVertex;
  int _iPrimaryVertex;
  int _nGoodVertices;

  // Generator level info
  Candidate* _genCand;
  TLorentzVector _neutrinosFromBosonsP4;
  TLorentzVector _neutrinosFromTausP4;

  // event category (from genEvent)
  void        defineCateg();
  std::string _categ_long;
  std::string _categ;

  bool _inAcceptance;

  //unclustered particules
  CandList _uncPfCands;

  // list of input file names
  std::vector< std::string > _filename; 

  // description of available collections
  std::string _collectionFileName;

  // map CaloTower index --> CaloTower candidate
  std::map< int, Candidate* > _caloTowerByIndex;

  // multimap CaloJet --> CaloTowers
  CandIdMap _caloJetTowerMap;

  // "iso" maps 
  std::map< std::string, CandIdIsoMap > _isoMap;

  //Ecal maps
  //GHM  std::map< std::string, CandIdEcalMap > _ecalMap;
  std::map< std::string, CandIdRecHitMap > _ecalRecHitMap;
  std::map< std::string, CandIdPSRecHitMap > _ecalPSRecHitMap;

  //Candidate maps
  std::map< std::string, CandMap > _candMap;

  //Candidate assoc maps
  std::map< std::string, CandAssoc > _candAssoc; 

  // mc association map
  CandIdDrMap _mcMatchingMap;

  // histograms of Et summed in 72 phi bins
  std::map< std::string, TH1* > _sumOfEt;

  
public: // temporary

  friend class SampleAnalysis;
  // constructors, destructor are only for friends
  EventManager( const char* filename, const char * collectionFileName );
  EventManager( std::vector< std::string > filename, const char * collectionFileName );
  virtual ~EventManager( );

protected: // temporary

  // singleton
  static EventManager* _instance;

  // interface to the data
  EventServer* _a;

private:

  // cast away const privately
  EventManager* me() const { return const_cast<EventManager*>(this); }

  void checkFileName( const char* filename );

  // cannot be templated
  void setBoolInfo(  CandInfo*, const std::string &, const std::string &, int v=-1 );
  void setIntInfo(   CandInfo*, const std::string &, const std::string &, int v=-1 );
  void setFloatInfo( CandInfo*, const std::string &, const std::string &, int v=-1 );
 
  // utility function
  Candidate* createMcCand_( int ii);

  // class definition
  ClassDef( EventManager, 0 )    
     
};

#endif
