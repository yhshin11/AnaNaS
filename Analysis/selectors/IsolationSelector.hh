#ifndef IsolationSelector_h
#define IsolationSelector_h

#include <vector>
#include "Analysis/core/EventManager.hh"
#include "Analysis/selectors/SelectorBase.hh"
#include "Analysis/utils/Constants.hh"

#include "MECore/src/ME.hh"
#include "MECore/src/MEEBDisplay.hh"
#include "MECore/src/MEEEDisplay.hh"

class IsolationSelector : public SelectorBase
{
public:

  IsolationSelector( unsigned settings_, float iso0_=0. );

  virtual ~IsolationSelector() {}

  // return true if the candidate is selected 
  virtual bool accept( const Candidate& cand ) const;

  // set the isolation limit
  void setIso( float iso0_ ) { iso0 = iso0_; }

  //
  // isolation variables
  //
  void computeIso( const Candidate& );

protected:
  unsigned _settings;

public:

  float r_EB, z_EB, e_EE, r_EE, z_EE, r_HB, z_HE;

  float pt0_trk, dR_trk_in,  dR_trk_out;
  float d0_trk, dz0_trk;
  float strip_trk;
  
  bool  ecalUseIsoDeposits;
  bool isBarrel;
  int iecal;
  float et0_ecal;
  float dR_ecal_in[2], dR_ecal_out[2];
  float e0_ecal[2], strip_ecal[2];

  float dR_hcal_in, dR_hcal_out;
  float et0_hcal;
  float e0_hcal;

  float a_trk;
  float a_ecal;
  float a_hcal;

  float iso;
  float S;
  float S_trk;
  float S_ecal;
  float S_hcal;
  int n_trk;
  int n_hcal;
  int n_ecal;

  float iso0;

  float BField;

private:

  IsolationSelector* me() const { return (IsolationSelector*)this; }

  ClassDef( IsolationSelector, 0 )  
};

#endif

     
