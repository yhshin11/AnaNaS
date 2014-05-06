#ifndef SelectorBase_h
#define SelectorBase_h

#include "Analysis/core/Candidate.hh"

//
// candidate selector base class
//
class SelectorBase      
{
public:

  SelectorBase() : _verbose(false) {}
  virtual ~SelectorBase() {}

  // return true if the candidate is selected 
  virtual bool accept( const Candidate& cand ) const;

  //
  // get a list of selected candidates 
  //    from a list of candidates,
  void getList( const CandList& inputList, 
		CandList& outputList, 
		const std::vector<int>* uids=0 ) const;

  void getSortedList( const CandList& inputList,
		      CandList& outputList,
		      int sorting=Candidate::kSortByPt,
		      const std::vector<int>* uids=0 ) const;

  void setVerbose() {_verbose = true; } //!!!

  //
  // helpers
  //
  virtual void sortList(  CandList& list, int sorting ) const;

  bool _verbose; //!!!

  ClassDef( SelectorBase, 0 )  

};

#endif

     
