#ifndef OfficialIsolationSelector_h
#define OfficialIsolationSelector_h

#include <vector>
#include "Analysis/core/EventManager.hh"
#include "Analysis/selectors/SelectorBase.hh"
#include "Analysis/utils/Constants.hh"

#include "MECore/src/ME.hh"
#include "MECore/src/MEEBDisplay.hh"
#include "MECore/src/MEEEDisplay.hh"

class OfficialIsolationSelector : public SelectorBase
{
public:

  OfficialIsolationSelector( unsigned settings_, float iso0_=0. );

  virtual ~OfficialIsolationSelector() {}

  // return true if the candidate is selected 
  virtual bool accept( const Candidate& cand ) const;

  // set the isolation limit
  void setIso( float iso0_ ) { iso0 = iso0_; }

  //
  // isolation variables
  //
  void computeIso( const Candidate& );

  //Load Xtal geometry DataBase
  void loadXtalDatabase();
  

protected:
  unsigned _settings;

public:

  float r_EB, z_EB, e_EE, r_EE, z_EE, r_HB, z_HE;

  float pt0_trk, dR_trk_in,  dR_trk_out;
  float d0_trk, dz0_trk;
  float strip_trk;
  
  bool  ecalUseIsoDeposits;
  bool useGeometryDatabase;
  bool isBarrel;
  int iecal;
  float et0_ecal[2];
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

  OfficialIsolationSelector* me() const { return (OfficialIsolationSelector*)this; }

  bool _loadDatabase;
  std::map<int, std::pair<float, float> > MapXtal;
  std::map<int, std::pair<float, float> >::const_iterator itMapXtal;

  ClassDef( OfficialIsolationSelector, 0 )  
};

#endif

     
