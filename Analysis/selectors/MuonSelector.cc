#include "Analysis/selectors/MuonSelector.hh"

#include <cassert>

using namespace std;

#include "Analysis/core/CandInfo.hh"
#include "Analysis/core/EventManager.hh"

ClassImp( MuonSelector )

MuonSelector::MuonSelector( unsigned ptcut, unsigned sel ) : LeptonSelector( ptcut, sel)
{
}

void
MuonSelector::setFlags( const Candidate& cand ) const
{
  // only consider muon candidates
  int pdgCode_ = cand.pdgCode();
  assert( fabs( pdgCode_ )==13 );

  float pt  = cand.pt();
  float eta = cand.eta();

  // set all the flags to false
  CandInfo* info = cand.info();
  for( unsigned isel=kVeryLoose; isel<kNSel; isel++ )
    info->setBool( selection[isel] , false );

  bool isMuonFiducial = (fabs(eta)<=2.4);
  info->setBool( "isMuonFiducial", isMuonFiducial );
 
  if(!isMuonFiducial) return;

  EventManager& _e = *(EventManager::e());

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

  bool triggerOK = (fabs(eta)<=2.1);
  info->setBool( "triggerOK", triggerOK );


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
  

  float d0Max        =  0.02;
  float normChi2Max  =  10;
  int   nvHitsMin    =  10;
  int   nMatchMin    =   2;
  int   nPixHitMin   =   1;
  int   nMuonHitMin  =   1;

  float d0 = cand.d0( &(_e.primaryVertex()));
  float dZ = cand.dz( &(_e.primaryVertex()));
  info->setFloat( "dPV", d0 );
  info->setFloat( "dZ", dZ );

  d0 = info->getFloat("d0");

  bool isGlobal(true);
  bool isGood(true);
  bool isPF(true);
  bool  muIdTracker = true;
  float normChi2(normChi2Max);
  int nvHits(nvHitsMin);
  int   nMatch(nMatchMin);
  int   nPixHit(nPixHitMin);
  int   nPixLayerHits(0);
  int   nMuonHit(nMuonHitMin);
  info->getBool("muidGlobalMuonPromptTight", isGlobal);
  info->getBool("muidTrackerMuonArbitrated",muIdTracker);
  info->getBool("isPFMuon", isPF);
  info->getBool("muidTMOneStationTight",isGood);
  info->getFloat("normChi2",normChi2);
  info->getInt("nTrkHits",nvHits); 
  info->getInt("nMatch",nMatch);
  info->getInt("nPixHits",nPixHit);
  info->getInt("nPixLayerHits",nPixLayerHits);
  info->getInt("nMuonHits",nMuonHit);
  
  int nTrkLayerHits =  info->getInt("nTrkLayerHits"); 

  float isoR03ch = info->getFloat("isoR03pfch");
  float isoR03nh = info->getFloat("isoR03pfnh");
  float isoR03ph = info->getFloat("isoR03pfph") ;
  float isoR03pu = info->getFloat("isoR03pfpu");

  float iso = isoR03ch + max(0., isoR03ph + isoR03nh - 0.5*isoR03pu );
  iso /= pt;
  
  if( isGood && nTrkLayerHits > 5 && nPixLayerHits > 1 &&
      normChi2 <1.8 ) {
    
    if(iso < 0.20) {
      info->setBool( selection[kVeryLooseIso], true);
      info->setBool( selection[kLooseIso], true);
    }
    if( d0<3. && dZ<30.) {
      info->setBool( selection[kVeryLooseID], true);
      info->setBool( selection[kLooseID], true);
    }
    if(iso < 0.20 && d0<3. && dZ<30.) {
      info->setBool( selection[kVeryLoose], true);
      info->setBool( selection[kLoose], true);
    }

  }

  if( isGlobal && isPF && normChi2 <10 && nMuonHit>1 &&
      nvHits>0 && nTrkLayerHits > 5 ) {
  
  if(iso < 0.12) {
      info->setBool( selection[kMediumIso], true);
      info->setBool( selection[kTightIso], true);
    }
    if( d0<.2 && dZ<.5) {
      info->setBool( selection[kMediumID], true);
      info->setBool( selection[kTightID], true);
    }
    if(iso < 0.12 && d0<.2 && dZ<.5) {
      info->setBool( selection[kMedium], true);
      info->setBool( selection[kTight], true);
    }

  }

  
  return;
}
