#include <algorithm>
#include <functional>
#include <list>
#include "Analysis/selectors/ConeSelector.hh"

using namespace std;

ClassImp( ConeSelector )

ConeSelector::ConeSelector( const Candidate& cand, 
			    float DRmax, float DRmin, 
			    float ptmin, float Jhs, int pdg ) : 
_cand(cand), _DRmax(DRmax), _DRmin(DRmin), _ptmin(ptmin), _Jhs(Jhs), _pdg(pdg) 
{
  init();
}

ConeSelector::ConeSelector( const Candidate& cand, 
			    float DRmax, const TParticlePDG* pdtEntry ) : 
_cand(cand), _DRmax(DRmax), _DRmin(0), _ptmin(0), _Jhs(-1), _pdg(0) 
{
  if( pdtEntry!=0 ) _pdg = pdtEntry->PdgCode();
  init();
}

void
ConeSelector::init()
{
  _P3   = _cand.p3();
  _pt   = _cand.pt();
  _eta  = _cand.eta();
  _Uid  = _cand.uid();
}

bool 
ConeSelector::accept( const Candidate& cand ) const
{
  indexDRPair p_;
  return getPair( cand, p_ );
}

bool 
ConeSelector::getPair( const Candidate& cand, indexDRPair& p ) const
{
  int Uid_ = cand.uid();

  if( Uid_==_Uid ) return false;  // same as the reference candidate

  float pt_ = cand.pt();
  if( pt_<=0     ) return false;   // against null energy calo towers
  if( pt_<_ptmin ) return false;   // pt cut
  
  int pdgCode_ = cand.pdgCode();
  if( _pdg!=0 )
    {
      if( fabs( pdgCode_ )!=_pdg ) return false;
    }

  float DR_ = dr( cand );
  if( DR_>_DRmax ) return false;
  if( DR_<_DRmin ) return false;
 
  if( _Jhs>0 )
    {
      float DEta_ = deta( cand );
      if( fabs(DEta_)<_Jhs ) return false;
    }

  p = make_pair( Uid_, DR_ );

  return true;
}

void  
ConeSelector::getPairs( const CandList& list, std::vector< indexDRPair >& listOfP ) const
{
  listOfP.clear();
  for( size_t ii=0; ii<list.size(); ii++ )
    {
      indexDRPair p_;
      if( !getPair( *list[ii], p_ ) ) continue;
      listOfP.push_back( p_ );
    }
  sort( listOfP.begin(), listOfP.end(), SortIndexDRPairs() ); 
}

float
ConeSelector::dr( const Candidate& cand ) const
{
  TVector3 P3_ = cand.p3();
  float DR_ = _P3.DeltaR( P3_ );
  return DR_;
}

float
ConeSelector::deta( const Candidate& cand ) const
{
  float eta_ = cand.eta();
  return eta_-_eta;
}

void
ConeSelector::sortList(  CandList& candlist, int sorting ) const
{
  if( candlist.size()==0 ) return;

  if( sorting!=Candidate::kSortDefault )
    {
      SelectorBase::sortList( candlist, sorting );
      return;
    }
  sort( candlist.begin(), candlist.end(), *this );
}

bool 
ConeSelector::operator()(  const Candidate* c1, const Candidate* c2 ) const
{
  return sortByDR( c1, c2 );
}
