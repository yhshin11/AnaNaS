#ifndef ConeSelector_h
#define ConeSelector_h

#include "Analysis/selectors/SelectorBase.hh"

class ConeSelector : public SelectorBase     
{
  const Candidate& _cand;
  TVector3 _P3;
  float _pt;
  float _eta;
  int     _Uid;

  float _DRmax; 
  float _DRmin;
  float _ptmin;
  float _Jhs;  // Jurassic half-strip in eta
  int   _pdg;
  
public:

  ConeSelector( const Candidate& cand, float DRmax, float DRmin=0., float ptmin=0., float Jhs=-1, int pdg=0 );
  ConeSelector( const Candidate& cand, float DRmax, const TParticlePDG* );
  virtual ~ConeSelector() {}

  float dr(   const Candidate& cand ) const;
  float deta( const Candidate& cand ) const;

  // get the index/DR pair for cand "cand", if within limits
  bool  getPair( const Candidate& cand, indexDRPair& p ) const;

  // given a list, get the sorted list of index/DR pairs with limits
  void  getPairs( const CandList& list, std::vector< indexDRPair >& listOfP ) const;

  bool operator()(  const Candidate* c1, const Candidate* c2 ) const;

  bool sortByDR( const Candidate* c1, const Candidate* c2 ) const
  {
    if( c1==0 || c2==0 ) return false;
    float dr1 = dr( *c1 );
    float dr2 = dr( *c2 );
    return dr1<dr2;      
  }

  // return true if the candidate is selected 
  virtual bool accept( const Candidate& cand ) const;

  // sorting
  virtual void sortList( CandList& list, int sorting ) const;

private:
  
  void init();

  ClassDef( ConeSelector, 0 )  
};


#endif

     
