#include <algorithm>
#include <cassert>

#include "Analysis/selectors/SelectorBase.hh"

using namespace std;

ClassImp( SelectorBase )

bool 
SelectorBase::accept( const Candidate& cand ) const
{
  return true;
}

void
SelectorBase::getList( const CandList& inputList,
		       CandList& outputList,
		       const vector<int>* uids ) const
{
  outputList.clear();
  for( size_t ii=0; ii<inputList.size(); ii++ )
    {
      Candidate* cand_  = inputList[ii];
      assert( cand_!=0 );
      
      // check whether the current uid is in the black list
      if( uids!=0 )
	{
	  int candUid_ = cand_->uid();
	  vector<int>::const_iterator it_ 
	    = find( uids->begin(), uids->end(), candUid_ );
	  if ( it_!=uids->end() ) continue;
	}

      // check if the candidate is accepted
      if( !accept( *cand_ ) ) continue;

      // candidate is accepted and not in the black list
      outputList.push_back( cand_ );
    }
}

void
SelectorBase::getSortedList( const CandList& inputList,
			     CandList& outputList,
			     int sorting,
			     const vector<int>* uids ) const
{
  getList( inputList, outputList, uids );
  sortList( outputList, sorting );
}

void
SelectorBase::sortList( CandList& list, int sorting ) const
{
  switch( sorting )
    {
    case Candidate::kSortByPt:
      sort( list.begin(), list.end(), Candidate::SortByPt() );
      return;
    case Candidate::kSortByEnergy:
      sort( list.begin(), list.end(), Candidate::SortByEnergy() );
      return;
    }
  return;
}
