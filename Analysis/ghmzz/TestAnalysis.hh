
#ifndef TestAnalysis_hh
#define TestAnalysis_hh

#include "Analysis/ghmzz/TupleAnalysis.hh"

class TestAnalysis : public TupleAnalysis
{
public:

  TestAnalysis( int im=0, string ip="2011" ) : TupleAnalysis(im, ip) 
  {
    analysis = "test";
  }

  virtual ~TestAnalysis() {}
			     
  void analyze_( float lumi );

  ClassDef( TestAnalysis, 0 )
};

#endif
