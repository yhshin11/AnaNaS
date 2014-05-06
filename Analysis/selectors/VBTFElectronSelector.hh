#ifndef VBTFElectronSelector_h
#define VBTFElectronSelector_h

#include "Analysis/selectors/LeptonSelector.hh"
#include "Analysis/selectors/OfficialIsolationSelector.hh"

class VBTFElectronSelector : public LeptonSelector     
{

public:

  enum {  kEt5=0, kEt10, kEt15, kEt20, kEt25,  kEt27,  kEt30,  kEt32, kNEt };
  enum {  kEid=0,    kEid95, kEid90, kEid85, kEid80, kEid70, kEid60, kEidAnti, kNEid };
  enum {  kIso=0,    kIso95, kIso90, kIso85, kIso80, kIso70, kIso60, kIsoAnti, kNIso };
  enum {  kPresel=0, kSel95, kSel90, kSel85, kSel80, kSel70, kSel60, kSelAnti, kNWP };
  static std::string level[kNWP];

  VBTFElectronSelector( unsigned etcut, unsigned sel );
  virtual ~VBTFElectronSelector() {}

  // return true if the candidate is selected 
  //  virtual bool accept( const Candidate& cand ) const;
  
  bool useStandardVariables;
  bool useOfficialIsoComputer;

  float detaCorrections(float, float); 
  float dphiCorrections(float, float); 

  void computeFlags( const Candidate& cand );
  bool isFlag( const Candidate& cand,  std::string isoid );
  
private:
  
  //  unsigned _etcut;
  //  unsigned _sel;

  VBTFElectronSelector* me() const { return (VBTFElectronSelector*)this; }


  void setFlags( const Candidate& ) const;
  OfficialIsolationSelector IsoElec;
  ClassDef( VBTFElectronSelector, 0 )  
};

#endif

     
