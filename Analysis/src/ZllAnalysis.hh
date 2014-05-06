#ifndef ZllAnalysis_h
#define ZllAnalysis_h

#include "Analysis/core/SampleAnalysis.hh"

class ZllAnalysis : public SampleAnalysis
{
public:

  ZllAnalysis( Sample& sample, EventManager& manager );
  virtual ~ZllAnalysis();

private:

  virtual void bookHistograms();
  virtual bool analyzeEvent();
  virtual void writeHistograms();

private:
  
  // build the Z->ll solutions
  //  imode=0 (Z->ee),  1 (Z->mumu) 
  enum Zmode { kZeeMode=0, kZmmMode, kNZmode };
  size_t buildZllLists( int imode );
  size_t nZHypos[kNZmode];
  static std::string ZmodeName(          int imode );
  static std::string ZListName(          int imode, size_t ihyp );
  static std::string leptonName(         int ilepton );
  static std::string leptonListName(     int ilepton );
  static std::string leptonListName(     int ilepton, size_t ihyp );

  void declareInt(   string var, int   value );
  void declareFloat( string var, float value );
  void declareBool(  string var, bool  value );
  
  map< string, float > _float; 
  map< string, int   > _int; 
  map< string, bool  > _bool; 

  int _iline;
  map< string, int > _hltLines;

  bool fillMCList( const CandList&, int& , int );

  ClassDef( ZllAnalysis, 0 )  
};

#endif

     
