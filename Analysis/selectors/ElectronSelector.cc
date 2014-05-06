#include "Analysis/selectors/ElectronSelector.hh"

#include <cassert>

using namespace std;

#include "Analysis/core/EventManager.hh"
#include "Analysis/core/CandInfo.hh"

ClassImp( ElectronSelector )

ElectronSelector::ElectronSelector( unsigned ptcut, unsigned sel ) : LeptonSelector(EventManager::kElectron, sel)
{
}

void
ElectronSelector::setFlags( const Candidate& cand ) const
{
  // only consider electron candidates
  int pdgCode_ = cand.pdgCode();
  assert( fabs( pdgCode_ )==11 );

  float pt  = cand.pt();
  float eta = cand.eta();

  // set all the flags to false
  CandInfo* info = cand.info();
  for( unsigned isel=kVeryLoose; isel<kNSel; isel++ )
    info->setBool( selection[isel], false );

  //fiducial
  bool isEleFiducial= (fabs(eta)<2.5 && !(fabs(eta)>1.444 && fabs(eta)<1.56) );
  info->setBool( "isEleFiducial", isEleFiducial );
 
  if(!isEleFiducial) return;


  //transverse momentum


  bool isPt30 = pt>=30;
  bool isPt25 = pt>=25;
  bool isPt20 = pt>=20;
  bool isPt15 = pt>=15;
  bool isPt10 = pt>=10;
  bool isPt5  = pt>=5;
  info->setBool( "isPt5",  isPt5  );
  info->setBool( "isPt10", isPt10 );
  info->setBool( "isPt15", isPt15 );
  info->setBool( "isPt20", isPt20 );
  info->setBool( "isPt25", isPt25 );
  info->setBool( "isPt30", isPt30 );

 bool isAbovePtThreshold = false;
  if( _settings==kPt5 )
    {
      if( info->getBool( "isPt5" ) ) isAbovePtThreshold = true; 
    }
  else if( _settings==kPt10 )
    {
      if( info->getBool( "isPt10" ) ) isAbovePtThreshold = true; 
    }
  else if( _settings==kPt15 )
    {
      if( info->getBool( "isPt15" ) ) isAbovePtThreshold = true; 
    }
  else if( _settings==kPt20 )
    {
      if( info->getBool( "isPt20" ) ) isAbovePtThreshold = true; 
    }
  else if( _settings==kPt25 )
    {
      if( info->getBool( "isPt25" ) ) isAbovePtThreshold = true; 
    }
  else if( _settings==kPt30 )
    {
      if( info->getBool( "isPt30" ) ) isAbovePtThreshold = true; 
    }
  else
    {
      abort();
    }
  
  if( !isAbovePtThreshold ) return;


  //ID and isolation
 
  int testVL = info->getInt("vetoID");
  int testL = info->getInt("looseID");
  int testM = info->getInt("mediumID");
  int testT = info->getInt("tightID");

  
  //ugly...
  info->setBool( selection[kVeryLoose], testVL==1023 );
  info->setBool( selection[kVeryLooseID], testVL==895 || testVL==1023 );
  info->setBool( selection[kVeryLooseIso], (testVL&128)!=0 );
  
  info->setBool( selection[kLoose], testL==1023 );
  info->setBool( selection[kLooseID], testL==895 || testL==1023 );
  info->setBool( selection[kLooseIso], (testL&128)!=0 );
  
  info->setBool( selection[kMedium], testM==1023 );
  info->setBool( selection[kMediumID], testM==895 || testM==1023 );
  info->setBool( selection[kMediumIso], (testM&128)!=0 );
  
  info->setBool( selection[kTight], testT==1023 );
  info->setBool( selection[kTightID], testT==895 || testT==1023 );
  info->setBool( selection[kTightIso], (testT&128)!=0 );
  
  return;
}
