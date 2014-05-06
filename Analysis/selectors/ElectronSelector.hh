#ifndef ElectronSelector_h
#define ElectronSelector_h

#include "Analysis/selectors/LeptonSelector.hh"

class ElectronSelector : public LeptonSelector     
{

public:

  enum {  kPt5=0, kPt10,  kPt15,  kPt20,  kPt25,  kPt30, kNPt };

  ElectronSelector( unsigned ptcut, unsigned sel );
  virtual ~ElectronSelector() {}

private:
  
  void setFlags( const Candidate& ) const;

  ClassDef( ElectronSelector, 0 )  
};

#endif

     
