#ifndef MCAnalysis_h
#define MCAnalysis_h

#include <vector>
#include <map>
#include <string>

#include <TString.h>
#include <TObjString.h>
#include <TProfile.h>

#include "Analysis/core/SampleAnalysis.hh"
#include "Analysis/utils/Constants.hh"
#include "Analysis/utils/KineUtils.hh"
#include "Analysis/tools/CandUtil.hh"

#include "TGraph.h"

class KFactor;

//#include "METNN.hh"

class MCAnalysis : public SampleAnalysis
{
public:

  MCAnalysis( Sample& sample, const std::string & collectionFileName  );
  virtual ~MCAnalysis();

private:

  //Herited functions
  virtual void bookHistograms();
  virtual bool analyzeEvent();
  virtual void writeHistograms();

  ClassDef( MCAnalysis, 0 )  

  map< string, KFactor* > _kFactors;

  int _N;
  map< string, float > _counters;
};

#endif

     
