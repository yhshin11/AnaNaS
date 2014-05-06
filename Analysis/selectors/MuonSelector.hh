#ifndef MuonSelector_h
#define MuonSelector_h

#include "Analysis/selectors/LeptonSelector.hh"

class MuonSelector : public LeptonSelector     
{  
public:

  enum {  kPt5=0, kPt10,  kPt15,  kPt20,  kPt25,  kPt30, kNPt };

  MuonSelector( unsigned ptcut, unsigned sel );
  virtual ~MuonSelector() {}

private:

  void setFlags( const Candidate& ) const;

  ClassDef( MuonSelector, 0 )  
};

#endif

     
