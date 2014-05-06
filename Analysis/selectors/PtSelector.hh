#ifndef PtSelector_h
#define PtSelector_h

#include <vector>
#include "Analysis/selectors/SelectorBase.hh"

//
// Pt Selector
//
class PtSelector : public SelectorBase     
{
  float _ptmin;
  float _etamax;
  std::vector< float > _etaveto;
  
public:

  PtSelector( float ptmin=0, float etamax=1000 ) : _ptmin( ptmin ), _etamax( etamax ) {}
  virtual ~PtSelector() {}

  void setPtmin(  double ptmin  ) { _ptmin  = ptmin;  }
  void setEtamax( double etamax ) { _etamax = etamax; }
  void vetoRangeInEta( double etamin, double etamax )
  { _etaveto.push_back( etamin ); _etaveto.push_back( etamax ); }

  // return true if the candidate is selected 
  virtual bool accept( const Candidate& cand, bool withCalo="false" ) const;

  ClassDef( PtSelector, 0 )  
};

#endif

     
