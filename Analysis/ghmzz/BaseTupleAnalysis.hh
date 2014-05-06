#ifndef BaseTupleAnalysis_hh
#define BaseTupleAnalysis_hh

#include <string>
#include <vector>
#include <map>

#include "Analysis/core/Sample.hh"
#include "Analysis/core/Candidate.hh"
#include "Analysis/utils/HistoManager.hh"
#include "Analysis/utils/DatasetManager.hh"
#include "Analysis/utils/TupleManager.hh"

using namespace std;

typedef vector<double> vectorFloat; 

typedef multimap< string, Sample* > SampleMap;
typedef multimap< string, Sample* >::iterator SampleIt;
typedef pair< string, Sample* > SamplePair;

class BaseTupleAnalysis
{
public:

  BaseTupleAnalysis( const char* sampleSet=0, int n=-1 );
  virtual ~BaseTupleAnalysis() {}
  void go();
  
protected:
			 
  virtual void init();
  virtual void eventLoop();
  virtual bool analyzeEvent();
  virtual void finalize();
  virtual void print();

  bool _verbose;

  string histfile;
  void setOutputFile();

  static string ZllName[2];

  // input
  TupleManager   tm;

  // output
  HistoManager   hm_out;
  DatasetManager dm_out;
  TupleManager   tm_out;

  // counters
  map<string,int> nAnalyzed;
  map<string,int> nSelected;
  int n_max;

  // samples
  float _lumi;
  map< string, Sample* > samples;
  void defineSamples();
  vector<string> _sampleSet;
  Sample* _sample;
  string _sampleSetName;
  string _sampleName;
  string _subSampleName;
  bool isData() const;
  bool isData2011A() const;
  bool isData2011B() const;
  bool isDoubleElectron() const;
  bool isDoubleMu() const;
  bool isMC() const;
  SampleMap sampleList;
  map< string, float > integratedLumi;
  string dirName;

  // event quantities  
  int irun;
  int ievent;
  unsigned long run;
  unsigned long event;
  int nMC, nHLT, nVtx, nJet;
  map< int, string > hltLines;
  string categ;
  string categ_long;
  Candidate* genCand;
  map<int,Vertex*>     vertices;
  map<int,CandList> leptonAtVertex;
  CandList caloOnlyJets;
  map<int,CandList> ZAtVertex;
  map<int,CandList> jetAtVertex;
  Candidate* MetCand;
  vector< string > firedHltLines;
  vector<int> nLep;
  vector< map<int,Candidate*> > leptons;
  vector<int> nZll;
  vector< CandList > Zll;

  // histograms  
  vector< HBool > hbool;
  void defineTemplate( string templ, 
		       int nbin, float min, float max );
  void defineTemplate( string templ, 
		       int nbin, float* xbins );
  TH1* hist( string templ, string cut, string dir, string prefix );
  void fill( string templ, string cut, string dir, string prefix,
	     float val, float weight=1, bool overflow=false, bool norm=true );

  // input ntuples
  void declareInt(   string var );
  void declareFloat( string var );
  void declareBool(  string var );
  map< string, float > _float; 
  map< string, int   > _int; 
  map< string, bool  > _bool; 

private:
  
  void init_();
  void finalize_();

  ClassDef( BaseTupleAnalysis, 0 )  
};

#endif
