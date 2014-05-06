#ifndef AnalysisFactory_h
#define AnalysisFactory_h

#include <string>

#include "Analysis/core/Sample.hh"
#include "Analysis/core/SampleAnalysis.hh"

class AnalysisFactory
{
public:
  //
  // get an analysis object from name
  // 
  //  static SampleAnalysis* get( const std::string &, Sample&, const std::string &);
  static SampleAnalysis* get( const std::string &, Sample&, EventManager&);
  //static SampleAnalysis* get( std::string &, Sample&, std::string &);

};

#endif

     
