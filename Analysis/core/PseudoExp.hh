#ifndef PseudoExp_h
#define PseudoExp_h

#include <vector>

class RooDataSet;
#include <RooArgSet.h>
#include <RooDataSet.h>
#include <RooRealVar.h>

#include "Analysis/core/Sample.hh"

class PseudoExp
{
public:

  PseudoExp( Sample& sample, float lumi=10. );
  virtual ~PseudoExp();

  void setLumi( float lumi ) { _lumi=lumi; }
 
protected:

  void checkEvent( int ievt );
  void saveEvent();
  void createVar( const char* name, float min=-1000., float max=1000. );
  void setVarValue( const char* name, float val );
  void writeDatasets();

private:

  RooDataSet& dataset();
  void finalisePseudoExperiment();
  void determineLastEvent();

  Sample& _theSample;
  RooArgSet _rooArgSet;
  RooDataSet* _dataset;
  float  _lumi;
  int    _n;
  int    _current;
  int    _first;
  int    _last;
  bool   _random;

  // the variables
  std::map< std::string, RooRealVar* >    _var;

  std::vector< RooDataSet* > _datasetStore;

  ClassDef( PseudoExp, 0 )  
};

#endif

     
