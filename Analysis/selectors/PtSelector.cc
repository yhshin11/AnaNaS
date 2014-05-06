#include "Analysis/selectors/PtSelector.hh"

using namespace std;

// Pt selector
ClassImp( PtSelector )

bool 
PtSelector::accept( const Candidate& cand, bool withCalo ) const
{
  float pt  = cand.pt();
  float eta = fabs( cand.eta() );

  CandInfo* info = cand.info();
  if(withCalo)
    {
      pt =  info->getFloat("caloEnergy")/cosh(info->getFloat("caloEta") );
      eta = fabs(info->getFloat("caloEta"));
    }
  
  //cout<<eta<<"   "<<_etamax<<"   "<<_etaveto[0]<<"    "<<_etaveto[1]<<endl;

  if( pt<_ptmin )        return false;

  if( eta>_etamax )      return false;

  for( size_t ii=0; ii<_etaveto.size(); ii+=2 )
    {
      if( eta>_etaveto[ii] && eta<_etaveto[ii+1] ) return false;
    }

  return true;
}

