#ifndef EventServer_hh
#define EventServer_hh

//
// Author  : Gautier Hamel de Monchenault, Saclay
//

//
// Simple CMS analysis
//   based on ROOT input files created by CMSSW module GhmNtupleMaker
//   (source code available in CVS in:
//         UserCode/GautierHdeM/src/GhmAnalysis/GhmNtupleMaker)
//

#include <iostream>

// std headers
#include <string>
#include <vector>
#include <map>
#include <list>

// root headers
#include <TROOT.h>
#include <TChain.h>
#include <TH2I.h>
#include <TH2F.h>
#include <TFile.h>
#include <TString.h>
#include <TBits.h>

#include "Analysis/utils/FileUtils.hh"

class TreeManager;

class EventServer 
{

public:

  static bool isMC;

  enum { kZee=0, kZmm, kZll };
  enum { kWe=0, kWm, kWl };
  static std::string Zmode[kZll];
  static std::string Wmode[kWl]; //MM add W decay mode
  static std::string lepton[kZll];

  static bool parsed_;

  static std::vector< std::string > caloMET;
  static int caloMET_default; 

  static std::vector< std::string > recoMET;
  static int recoMET_default; 

  static std::vector< std::string > pfMET;
  static int pfMET_default; 

  static std::vector< std::string > patMET;
  static int patMET_default; 

  static std::vector< std::string > genMET;
  static int genMET_default; 
 
  static std::vector< std::string > caloJets;
  static int caloJets_default; 

  static std::vector< std::string > pfJets;
  static int pfJets_default; 
 
  static std::vector< std::string > pfCandidates;
  static int pfCandidates_default; 
 
  static std::vector< std::string > trackJets;
  static int trackJets_default; 
 
  static std::vector< std::string > genJets;
  static int genJets_default; 

  // contructor/destructor
  EventServer( EventList& events, const char* collectionFileName );
  EventServer( const char* filename, const char * collectionFileName );
  EventServer( std::vector< std::string >  filenames, const char * collectionFileName );
  virtual ~EventServer();

  // options & parameters
  static int   verbosity;

  // event navigation
  virtual bool event( int ievt );
  virtual bool eventInFile( int fileNum, int evtNum );
  bool nextEvent();
  bool previousEvent();

  int getNevtInFile() { return _nevent[0];};

  // event navigation
  int eventNumber;  // numevent = _nevent[ifile] + ievent;
  int ifile;
  int ievent;
  int curRun;
  int curLumi;
  int curEvent;
  int time;
  TString fName;
  int n(  const std::string &, int v=-1 );

  // get data by prefix and varName (and index for arrays)
  bool        get( const std::string &, const std::string &, float&      , int v=-1 );
  bool        get( const std::string &, const std::string &, int&        , int v=-1 );
  bool        get( const std::string &, const std::string &, std::string&, int v=-1 );
  bool        get( const std::string &, const std::string &, TBits&      , int v=-1 );
  int         get_i( const std::string &, const std::string &, int v=-1 );
  float       get_f( const std::string &, const std::string &, int v=-1 );
  const char* get_s( const std::string &, const std::string &, int v=-1 );
  TBits*      get_b( const std::string &, const std::string &, int v=-1 );
  std::vector<int>&   get_vi( const std::string &, const std::string &, int v=-1 );
  int   get_vi( const std::string &, const std::string &, size_t ii,    int v=-1 );
  std::vector<float>& get_vf( const std::string &, const std::string &, int v=-1 );
  float get_vf( const std::string &, const std::string &, size_t ii,    int v=-1 );


  //special function for nEvntProc summing when files are merged
  int get_si( const std::string &, const std::string &, int v=-1 );


  // object loading
  bool load( const std::string & prefix, int n, int v=-1 ); 

  // access to bit names
  std::map< std::string, size_t >::const_iterator itHlt;
  std::map< std::string, size_t > hlt;
  std::vector< std::string > eid;
  std::vector< std::string > eAcc;
  std::vector< std::string > muid;
  std::vector< std::string > rmuid;
  std::vector< std::string > phid;
  std::vector< std::string > tauid;

  // returns true if an HLT line is fired
  bool isFired( const std::string &, bool noVersion=false );
  void firedLines( std::map< int, std::string >& lines );

  //Get Rho Fast jet
  float getRFJ();
  
  //Get the number of interactions
  int getnInteract( std::string str="" );
  float getTrueNInt();

  //return true if the photon is identified
  bool isPhidBitSet( int ii ) { return _phidBits[ii]; }
  
  // return true if lepton is identified
  bool isEidBitSet( int ii ) { return _eidBits[ii]; }
  bool isMidBitSet( int ii ) { return _muidBits[ii]; }
  bool isTidBitSet( int ii ) { return _tauidBits[ii]; }
  
  //return true if lepton is in acceptance definition
  bool isEAccBitSet( int ii ) { return _eAccBits[ii]; }

  // global event quantities
  bool hasPrimaryVtx;
  float ptHat;
  int processID;
  float rhoFJ;
  
  float trueNI;

  int nIBXm1;
  int nIBX;
  int nIBXp1;

  //run and job quantities
  int nEventProc;

  void oneLine( ostream& o );

  void getEntry( const std::string &, unsigned jj, int v=-1 );

  std::string getName( const std::string & prefix, int v=-1, bool sep=true );
  
public:
  // used to be private
  bool openFile( const std::string & );
  void closeFile(); 
  bool checkFiles();
  void init(const char * collectionFileName);

private:
 
  // data store
  std::map< std::string, int            > map_i;    
  std::map< std::string, float          > map_f;  
  std::map< std::string, std::string*   > map_s;  
  std::map< std::string, TBits*         > map_b;  
  std::map< std::string, std::vector<float>* > map_vf;  
  std::map< std::string, std::vector<int>*   > map_vi;  

  // register variables by "prefix" and "varName"
  void register_i(   const std::string &, const std::string &, int v=-1 );
  void register_f(   const std::string &, const std::string &, int v=-1 );
  void register_s(   const std::string &, const std::string &, int v=-1 );
  void register_b(   const std::string &, const std::string &, int v=-1 );
  void register_vi(  const std::string &, const std::string &, int v=-1 );
  void register_vf(  const std::string &, const std::string &, int v=-1 );

  // bit management
  TBits _hltBits;
  TBits _eidBits;
  TBits _eAccBits;
  TBits _muidBits;
  TBits _phidBits;
  TBits _tauidBits;

  // file management
  std::vector< std::string > _filenames;
  std::vector< int > _nevent;
  bool _filesChecked;
  bool _hasEventList;
  std::multimap< int, int > _eventList;
  std::multimap< int, int >::iterator _eventListIt;

  // event management 
  std::map< std::string, std::vector<Int_t> > evt_first;

  // tree management
  TreeManager* _treeMgr;
  std::string getKey(      const std::string & prefix, int v=-1    );
  std::string getTreeName( const std::string & prefix, int v=-1 );
  //  TTree*      getTree( std::string );
  std::map< std::string, std::string > ntuNames;
  std::map< std::string, std::pair< std::string, std::string> > treeNames;
  void setPrefix( const std::string &, const std::string & );
  void register_t( const std::string & prefix, const std::string & key );
  void register_t( const std::string & prefix, int v, const std::string & key1, const std::string & key2 );
  std::string concatanate( const std::string &, const std::string &, bool sep=true );
  std::map< std::string, Int_t > evt_n;
  int n0( const std::string &, int v=-1 );

ClassDef(EventServer,0) // EventServer 
};

#endif

