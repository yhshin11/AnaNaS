#ifndef LeptonSelector_h
#define LeptonSelector_h

#include <vector>
#include "Analysis/core/EventManager.hh"
#include "Analysis/selectors/SelectorBase.hh"

//
// Pt Selector
//
class LeptonSelector : public SelectorBase
{
public:

  enum { kVeryLoose=0, kLoose, kMedium, kTight,
	 kVeryLooseID, kLooseID, kMediumID, kTightID,
	 kVeryLooseIso, kLooseIso, kMediumIso, kTightIso, kNSel };
  static std::string selection[kNSel];
  
  LeptonSelector( unsigned settings, unsigned sel );

  virtual ~LeptonSelector() {}

  // return true if the candidate is selected 
  virtual bool accept( const Candidate& cand ) const;

protected:

  unsigned _settings;
  unsigned _sel;

  virtual void setFlags( const Candidate& ) const=0;

private:


  ClassDef( LeptonSelector, 0 )  
};

#endif

     
