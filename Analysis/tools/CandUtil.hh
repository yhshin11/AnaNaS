#ifndef CandUtil_h
#define CandUtil_h

#include <vector>

#include "Analysis/core/Candidate.hh"

class CandUtil
{
public:

  virtual ~CandUtil() {}

  static bool overlap( const Candidate& cand1, const Candidate& cand2 );

  static bool overlaps( const Candidate& cand, const CandList& veto );

  // get the vector of Uids of candidates overlapping with candidate "cand"
  static void getUids( const Candidate& cand, std::vector<int>& uids );

  // get the vector of Uids of candidates overlapping with candidates in list "list"
  static void getUids( const CandList& list, std::vector<int>& uids );
   
  // dispatch the list "list" to lists of non-overlapping candidates
  static void dispatchList( const CandList& list, std::vector< CandList >& lists );

  // get the list of candidates in list "list" non overlapping with candidates in list "veto"
  static void pruneList( const CandList& list, CandList& pruned, const CandList& veto );

  static void mergeLists( const CandList& list1, const CandList& list2, CandList& list );

  static void printList( ostream& o, const CandList& list, const char* name="", bool full=false );

  static void get( int pdgId, Candidate* cand, CandList& list, float ptmin=0. );
  static Candidate* first( int pdgId, Candidate* cand, float ptmin=0. );

  static float dR( Candidate* c1, Candidate* c2 );
  static float dPhi( Candidate* c1, Candidate* c2 );
  static Candidate* closest( Candidate* cand, const CandList& list, 
			     float& dR, float ptmin=0, 
			     bool checkPdgId=false, 
			     bool sameCharge=false, 
			     bool checkCharge=false, int charge=0 );

  static float mT( const Candidate* cand, const Candidate* met );

  ClassDef( CandUtil, 0 )
};  

#endif
