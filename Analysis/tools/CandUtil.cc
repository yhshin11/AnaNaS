#include "Analysis/tools/CandUtil.hh"

#include <algorithm>
#include <cassert>
#include <iostream>

using namespace std;

#include "Analysis/utils/Constants.hh"
#include "Analysis/utils/KineUtils.hh"

ClassImp(CandUtil)

void
CandUtil::getUids( const Candidate& cand, vector<int>& uids )
{
  int candUid_ = cand.uid();

  // check if this Uid is already there 
  vector<int>::const_iterator it_ =  find( uids.begin(), uids.end(), candUid_ );
  if( it_!=uids.end() ) return;

  uids.push_back( candUid_ );

  for( size_t idau=0; idau<cand.nDaughters(); idau++ )
    {
      const Candidate* dau_ = cand.daughter(idau);
      getUids( *dau_, uids );
    } 
//   if( cand.NDaughters()!=0 )
//     {
//       HCandListIterator it_ = cand.DaughterIterator();
//       Candidate* dau_;
//       while( (dau_=it_.Next())!=0 )
// 	{
// 	  getUids( *dau_, uids );
// 	}
//     }
}

void
CandUtil::getUids( const CandList& list, vector<int>& uids )
{
  for( size_t ii=0; ii<list.size(); ii++ )
    {
      const Candidate* cand = list[ii];
      getUids( *cand, uids );
    }
}

void 
CandUtil::dispatchList( const CandList& list, vector< CandList >& Lists )
{
  // fill final list by removing overlapping candidates
  CandList list1;
  CandList list2;

  // vectors of uids
  vector<int> uids;
  int vtx=-1;
  for( size_t ii=0; ii<list.size(); ii++ )
    {
      Candidate* cand = list[ii];
      if( vtx==-1 ) vtx = cand->vertexIndex();
      bool keep_(true);
      if( cand->vertexIndex() != vtx ) keep_ = false;
      for( size_t idau=0; idau<cand->nDaughters(); idau++ )
	{
	  const Candidate* dau = cand->daughter(idau);
	  int dauUid = dau->uid(); 
	  vector<int>::const_iterator it_;  
	  it_ = find( uids.begin(), uids.end(), dauUid );
	  if( it_!=uids.end() ) 
	    {
	      keep_ = false;
	      break;
	    }
	}
      if( !keep_ ) 
	{
	  list2.push_back( cand );
	}
      else
	{
	  getUids( *cand, uids );
	  list1.push_back( cand );
	}
    }
  Lists.push_back( list1 );

  if( list2.size()!=0 )
    {
      dispatchList( list2, Lists );
    }
}

bool
CandUtil::overlap( const Candidate& cand1, const Candidate& cand2 )
{
  vector<int> uids1; getUids( cand1, uids1 );
  vector<int> uids2; getUids( cand2, uids2 );
  for( size_t ii=0; ii<uids1.size(); ii++ )
    {
      vector<int>::const_iterator it_ = find( uids2.begin(), uids2.end(), uids1[ii] );
      if( it_!=uids2.end() )  return true;
    }
  return false;
}

bool
CandUtil::overlaps( const Candidate& cand, const CandList& veto )
{
  for( size_t ii=0; ii<veto.size(); ii++ )
    {
      if( overlap( cand, *veto[ii] ) ) return true;
    }
  return false;
}

void
CandUtil::pruneList( const CandList& list, CandList& pruned, const CandList& veto )
{
  for( size_t ii=0; ii<list.size(); ii++ )
    {
      Candidate* cand = list[ii];
      if( overlaps( *cand, veto ) ) continue;
      pruned.push_back( cand );
    }
}

void
CandUtil::printList( ostream& o, const CandList& list, const char* name, bool full )
{
  o << "Number of " << name << " candidates = " << list.size() << endl;
  for( size_t ii=0; ii<list.size(); ii++ )
    {
      Candidate* cand = list[ii];
      if( full )
	cand->print( o );
      else	
	cand->oneLine( o );
      //      o << endl;
    }
}

void
CandUtil::mergeLists( const CandList& list1, const CandList& list2, CandList& list )
{
  for( size_t ii=0; ii<list1.size(); ii++ )
    list.push_back( list1[ii] );
  for( size_t ii=0; ii<list2.size(); ii++ )
    list.push_back( list2[ii] );
}


void
CandUtil::get( int pdgId, Candidate* cand, CandList& list, float ptmin )
{
  assert( cand!=0 );
  float pt_  = cand->pt();
  int pdgId_ = cand->pdgCode();
  if( pdgId_==pdgId && pt_>ptmin ) 
    {
      list.push_back( cand );
      return;
    }
  for( size_t idau=0; idau<cand->nDaughters(); idau++ )
    {
      Candidate* dau_ = cand->daughter( idau );
      get( pdgId, dau_, list, ptmin );
    }
}

Candidate*
CandUtil::first( int pdgId, Candidate* cand, float ptmin )
{
  assert( cand!=0 );
  float pt_  = cand->pt();
  int pdgId_ = cand->pdgCode();
  if( pdgId_==pdgId && pt_>ptmin ) 
    {
      return cand;
    }
  for( size_t idau=0; idau<cand->nDaughters(); idau++ )
    {
      Candidate* dau_ = cand->daughter( idau );
      Candidate* first_ =  first( pdgId, dau_, ptmin );
      if( first_!=0 ) return first_;
    }
  return 0;
}

float 
CandUtil::dR( Candidate* c1, Candidate* c2 )
{
  return KineUtils::dR( c1->eta(), c2->eta(), c1->phi(), c2->phi() );
}

float 
CandUtil::dPhi( Candidate* c1, Candidate* c2 )
{
  return KineUtils::dPhi( c1->phi(), c2->phi() );
}

Candidate* 
CandUtil::closest( Candidate* cand, const CandList& list, 
		   float& dr, float ptmin, bool checkPdgId, 
		   bool sameCharge, bool checkCharge, int charge )
{
  dr = 5.;
  int pdgId_    = cand->pdgCode();
  float charge_ = cand->charge();
  Candidate* closest_(0);
  for( size_t ii=0; ii<list.size(); ii++ )
    {      
      Candidate* c_ = list[ii];
      if( c_->pt()<ptmin )                         continue;
      if( checkPdgId ) 
	{ if( abs(c_->pdgCode())!= abs(pdgId_) ) continue; }
      else if( sameCharge ) 
	{ if( c_->charge()!= charge_ ) continue; }
      else if( checkCharge )
	{ if( c_->charge()!= charge  ) continue; }
      float dr_ = dR(cand,c_);
      if( dr_<dr )
	{
	  dr=dr_;
	  closest_=c_;
	}
    }
  return closest_;
}

float
CandUtil::mT( const Candidate* cand, const Candidate* met )
{
  //   returns  the transverse mass for a candidate and a missing momentum
  TVector2 metP2  =  met->p2();
  TVector2 candP2 = cand->p2();
  
  float mT2 =  2*( candP2.Mod() * metP2.Mod() -  candP2 * metP2 );

  if( mT2<0 ) mT2=0;
  
  return sqrt( mT2 );
}
