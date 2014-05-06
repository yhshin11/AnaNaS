#ifndef PseudoExperiments_h
#define PseudoExperiments_h

class RooDataSet;
#include <RooArgSet.h>

#include "Analysis/fit/Model.hh"
#include "Analysis/core/SampleAnalysis.hh"

class PseudoExperiments : public SampleAnalysis
{
public:

  PseudoExperiments( Sample& sample, float lumi );
  virtual ~PseudoExperiments();
 
private:

  virtual void bookHistograms();
  virtual bool analyzeEvent();
  virtual void writeHistograms();

  void finalisePseudoExperiment();
  void determineLastEvent();

  Model _model;
  RooArgSet _rooArgSet;
  RooDataSet* _dataset;
  float  _lumi;
  int    _n;
  int    _first;
  int    _last;
  bool   _random;

  ClassDef( PseudoExperiments, 0 )  
};

#endif

     
