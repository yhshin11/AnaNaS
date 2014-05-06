#ifndef Commissioning_h
#define Commissioning_h

#include <vector>
#include <map>
#include <string>

class TH1;
class TH2;

#include "Analysis/core/SampleAnalysis.hh"

class Commissioning : public SampleAnalysis
{
public:

  Commissioning( Sample& sample, const std::string & collectionFileName );
  virtual ~Commissioning();
 
private:

  virtual void bookHistograms();
  virtual bool analyzeEvent();
  virtual void writeHistograms();

  std::map< std::string, TH1* > h_;
  std::map< std::string, TH2* > h2_;

  std::map< int, int > _runs;

  ClassDef( Commissioning, 0 )  
};

#endif
