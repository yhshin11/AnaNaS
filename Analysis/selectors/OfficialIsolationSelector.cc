#include "Analysis/selectors/OfficialIsolationSelector.hh"

#include <cassert>

#include "Analysis/utils/KineUtils.hh"
#include "Analysis/selectors/ConeSelector.hh"
#include "Math/GenVector/VectorUtil.h"
#include <fstream>
#include <iostream>
#include <stdlib.h>

using namespace std;

// Pt selector
ClassImp( OfficialIsolationSelector )

OfficialIsolationSelector::OfficialIsolationSelector( unsigned settings_, float iso0_ ) : _settings(settings_), iso0(iso0_)
{  

  //DataBase
  _loadDatabase =false;

  // fixme!!!
  r_EB = 129;
  z_EB = 280;
  e_EE = 1.49;
  r_EE = 150;
  z_EE = 325;
  r_HB = 167;
  z_HE = 350;
  BField = 3.8;

  // settings (see CMS AN-2009/042, p7)
  if( _settings==EventManager::kMuon )
    {
      pt0_trk    = 0.;
      d0_trk     = 1.;
      dz0_trk    = 5.;
      dR_trk_in  = 0.01;
      dR_trk_out = 0.30;
      strip_trk  = -1;
      
      ecalUseIsoDeposits  = false; //MM By default, use RecHits
      et0_ecal[0]    = 0.200;
      et0_ecal[1]    = 0.200;
      
      e0_ecal[0]     = 0.120;
      dR_ecal_in[0]  = 0.07;
      dR_ecal_out[0] = 0.30;
      strip_ecal[0]  = -1; 
      
      e0_ecal[1]     = 0.450; 
      dR_ecal_in[1]  = 0.07;
      dR_ecal_out[1] = 0.30;
      strip_ecal[1]  = -1;
      
      dR_hcal_in  = 0.10;
      dR_hcal_out = 0.30;
      et0_hcal    = 0.500; 
      e0_hcal     = 0.600; 
    }
  else if( _settings==EventManager::kElectron )
    {
      pt0_trk    = 0.7;
      d0_trk  = 9999.;
      dz0_trk = 0.2;
      dR_trk_in  = 0.015;
      dR_trk_out = 0.30;
      strip_trk  = 0.015;
      
      ecalUseIsoDeposits  = false; //MM By default, use RecHits
      useGeometryDatabase = true; //MM By default, use Xtal database

      et0_ecal[0]    = 0.000;
      e0_ecal[0]     = 0.080;
      dR_ecal_in[0]  = 3.0; //# Xtals
      dR_ecal_out[0] = 0.30;
      strip_ecal[0]  = 1.5;  //# Xtals  

      et0_ecal[1]    = 0.100;
      e0_ecal[1]     = 0.000; 
      dR_ecal_in[1]  = 3.0; //# Xtals
      dR_ecal_out[1] = 0.30;
      strip_ecal[1]  = 1.5;  //# Xtals
      
      dR_hcal_in  = 0.15;
      dR_hcal_out = 0.30;
      et0_hcal    = 0.000;
      e0_hcal     = 0.000;
    }
  else
    assert(0);

  a_trk  = 1;
  a_ecal = 1;
  a_hcal = 1;

  iso0 = 0.;

}

bool 
OfficialIsolationSelector::accept( const Candidate& cand ) const
{
  CandInfo* info = cand.info();

  float iso    = 0;
  if( ! info->getFloat( "iso", iso ) )
    {
      me()->computeIso( cand );
      iso = info->getFloat( "iso" );
    }

  return iso>=iso0;
}

void 
OfficialIsolationSelector::computeIso( const Candidate& cand )
{
  bool dbg_ = false;

  CandInfo* info = cand.info();
  
  iso    = 0;
  S      = 0;
  S_trk  = 0;
  S_ecal = 0;
  S_hcal = 0;
  n_trk    = 0;
  n_hcal   = 0;
  n_ecal   = 0;

  float pt  = cand.pt();
  assert( pt>0 );
  float chg = cand.charge();
  float eta = cand.eta();
  float phi = cand.phi();
  TVector3 v3 = cand.pos();
  float x0  = v3.X();
  float y0  = v3.Y();
  float z0  = v3.Z();

  float r_ecal = r_EB;
  float z_ecal = z_EB;
  if( fabs(eta)>e_EE )
    {
      r_ecal = r_EE;
      z_ecal = z_EE;
    }
  float r_hcal = r_HB;
  float z_hcal = z_HE;

  EventManager& _e = *(EventManager::e());

  // make sure that we're dealing with leptons
  // and that the settings are correct
  int pdgCode_ = cand.pdgCode();
  bool isElectron = fabs( pdgCode_ )==11;
  bool isMuon     = fabs( pdgCode_ )==13;
  static string mapname_;

  Candidate* track;

  if( _settings==EventManager::kElectron )
    {
      assert( isElectron );
      mapname_ = "electron-EcalIsoDep";    
      track = _e.getSecond("electron-GsfTrack", (const_cast<Candidate*>(&cand) ) );
    }
  else if( _settings==EventManager::kMuon )
    {
      assert( isMuon );
      mapname_ = "muon-EcalIsoDep";
    }
  
  float etatrk = track->eta();
  TVector3 Tv3_ = track->pos();
  float Tz0_  = Tv3_.Z();
  
  const CandList& listTracks = _e.tracks();
  for( size_t ii=0; ii<listTracks.size(); ii++ )
    {
      


      const Candidate* trk = listTracks[ii];
      float pt_  = trk->pt();
      float eta_ = trk->eta();
      float deta_ = eta_-etatrk;
      
      float dR_= ROOT::Math::VectorUtil::DeltaR(trk->p3(), track->p3() ); 
      bool inCone_=true;
      if( dR_<dR_trk_in  )          inCone_=false;
      if( dR_>dR_trk_out )          inCone_=false;
      if( fabs( deta_ )<strip_trk ) inCone_=false;
      
      TVector3 v3_ = trk->pos();
      float x0_  = v3_.X();
      float y0_  = v3_.Y();
      float z0_  = v3_.Z();
      float r0_ = sqrt( x0_*x0_ + y0_*y0_ );
      
      float dzCut = fabs( z0_ - Tz0_ );

      bool keep_=inCone_;
      if( pt_<pt0_trk )        keep_=false; 
      if( r0_>d0_trk )         keep_=false;
      if( fabs(dzCut)>dz0_trk )  keep_=false;
      
      if( keep_ ) 
	{
	  n_trk++;
	  S_trk += pt_;
	}
    }
  
  float etaAtHcal(0);
  float phiAtHcal(0);
  KineUtils::straightToEnvelop( etaAtHcal, phiAtHcal, eta, phi, x0, y0, z0, r_hcal, z_hcal );
  KineUtils::adjust( phiAtHcal, phi );

  float eta0_hcal = eta;
  float phi0_hcal = phi;
  KineUtils::helixToEnvelop( eta0_hcal, phi0_hcal, chg*pt, eta, phi, x0, y0, z0, r_hcal, z_hcal, BField );
  
  float DRin_hcal  = dR_hcal_in;
  float DRout_hcal = dR_hcal_out;

  const CandList& listCT = _e.caloTowers();
  for( size_t iCT=0; iCT<listCT.size(); iCT++ )
    {
      const Candidate& ct = *listCT[iCT];
      const CandInfo* info_ = ct.info();
      float eta_ = ct.eta();
      float phi_ = ct.phi();
      //      float e_ = info_->getFloat("energy");
      //      float emEt_  = info_->getFloat("emEt");
      float hadEt_ = info_->getFloat("hadEt");
      
      if( hadEt_==0 ) continue;

      float hadE_ = KineUtils::e( hadEt_, eta_ );

      float dR_out=KineUtils::dR( etaAtHcal, eta_, phiAtHcal, phi_ );
      float dR_in =KineUtils::dR( eta0_hcal, eta_, phi0_hcal, phi_ );
      bool inCone_=true;
      if( dR_out>DRout_hcal ) inCone_=false;
      if( dR_in <DRin_hcal  ) inCone_=false;

      bool keep_=inCone_;
      if( hadEt_<et0_hcal ) keep_=false; 
      if( hadE_<e0_hcal )   keep_=false; 

      if( keep_ ) 
	{
	  n_hcal++;
	  S_hcal += hadEt_;
	}
    }
  
  float etaAtEcal(0);
  float phiAtEcal(0);
  KineUtils::straightToEnvelop( etaAtEcal, phiAtEcal, eta, phi, x0, y0, z0, r_ecal, z_ecal );
  KineUtils::adjust( phiAtEcal, phi );

  isBarrel = (fabs(eta)<=e_EE)?true:false;
  iecal = isBarrel?0:1;
  float DRin_ecal  = dR_ecal_in[iecal];
  float DRout_ecal = dR_ecal_out[iecal];
  float Jhs_ecal   = strip_ecal[iecal];
  
  //  if( ecalUseIsoDeposits )
  //    {      
   float eta0_ecal = etaAtEcal;
   float phi0_ecal = phiAtEcal;
  //  float eta0_ecal = eta;
  //float phi0_ecal = phi;
  if( isElectron )
    {
      float caloEta = info->getFloat("caloEta");
      float caloPhi = info->getFloat("caloPhi");
      // KineUtils::adjust( caloPhi, phi );
      eta0_ecal = caloEta;
      phi0_ecal = caloPhi;
    }
  else
    {      
      KineUtils::helixToEnvelop( eta0_ecal, phi0_ecal, chg*pt, eta, phi, x0, y0, z0, r_ecal, z_ecal, BField );
    }


 if( ecalUseIsoDeposits )
   { 
     // get the map of iso object for this candidate
     const CandIdIsoMap& map_ = _e.isoMap( mapname_ );
     CandIdIsoMapIterator it_lower_ = map_.lower_bound( cand.uid() );
     CandIdIsoMapIterator it_upper_ = map_.upper_bound( cand.uid() );
     CandIdIsoMapIterator it;
     
     for( it=it_lower_; it!=it_upper_; it++ )
       {
	 float e_   = (*it).second.val;
	 float eta_ = (*it).second.eta;
	 float phi_ = (*it).second.phi;
	 float et_  = e_;
	 
	 float deta_ = eta_-eta0_ecal;
	 float dR_out  = KineUtils::dR( eta0_ecal, eta_, phi0_ecal, phi_ );
	 float dR_in   = KineUtils::dR( eta0_ecal, eta_, phi0_ecal, phi_ );
	 //       //      cout << "deta/dphi " << deta_ << "/" << dphi_ << endl;
	 //       if( dphi_ >  Constants::pi ) dphi_-=2*Constants::pi; 
	 //       if( dphi_ < -Constants::pi ) dphi_+=2*Constants::pi;       
	 //      float dR_ = sqrt( deta_*deta_ + dphi_*dphi_ );
	 //      cout << "dphi/dr " << dphi_ << "/" << dR_ << endl;

	 bool inCone_=true;
	 if( dR_out>DRout_ecal ) inCone_=false;
	 if( dR_in <DRin_ecal  ) inCone_=false;
	 if( fabs( deta_ )<Jhs_ecal ) inCone_=false;
	 
	 bool keep_=inCone_;
	 int iecal_ = ( fabs( eta_ )<=e_EE ) ? 0 : 1;
	 if(  e_<e0_ecal[iecal_]  ) keep_=false;
	 if( et_<et0_ecal[iecal_] ) keep_=false;
	 
	 if( keep_ ) 
	   {
	   
	     n_ecal++;
	     S_ecal += et_;
	   }
       }
   }
 else { //MM Use RecHits
   
   if( useGeometryDatabase && !_loadDatabase ) {
     loadXtalDatabase();
   }
     


   size_t nRecHits = (size_t) _e.a().n("rh");
   for(size_t irh=0;irh<nRecHits;irh++) {
     assert( _e.a().load("rh", irh) );
     
     float e_ = _e.a().get_f("rh","energy");
     int ieta_ = _e.a().get_i("rh","ix");
     int iphi_ = _e.a().get_i("rh","iy");
     int iz_ = _e.a().get_i("rh","iz");

     float eta_;
     float phi_;
    
     //Spike Rejection
     float swissC = _e.a().get_f("rh","swissCross");				 
     if(swissC >= 0.95 && swissC>0 ) continue; 
     //********************

     if( useGeometryDatabase ) {
       
       int ilab =  10000*(85+ieta_) + 10*iphi_ + iz_ + 1;
       itMapXtal = MapXtal.find( ilab );
       if( itMapXtal == MapXtal.end() )
	 { cout<<" Xtal Missing in database :  "<<ilab<<"\n  Please contact Matthieu if this problem appears "<<endl; abort(); }
       else
	 {
	   eta_ = (*itMapXtal).second.first;
	   phi_ = (*itMapXtal).second.second;
	 }
     }
     else {
       if(iz_==0)
	 {
	   MEEBGeom::EtaPhiPoint point_ = MEEBDisplay::center( ieta_, iphi_ );
	   
	   eta_ = point_.first;
	   phi_ = point_.second * Constants::pi;
	 }
       else //Endcaps
	 {
	   MEEEGeom::EtaPhiPoint point_ = MEEEDisplay::center( ieta_, iphi_, iz_ );
	   
	   eta_ = point_.first;
	   phi_ = point_.second * Constants::pi;
	 }
     }
  
     float deta_ = eta_- eta0_ecal;
     float dR_out  = KineUtils::dR( eta0_ecal, eta_, phi0_ecal, phi_ );
     float dR_in   = KineUtils::dR( eta0_ecal, eta_, phi0_ecal, phi_ );
     
     double phiDiff= KineUtils::dPhi(phi_, phi0_ecal );
     dR_in = sqrt( deta_*deta_ + phiDiff*phiDiff ); 

     float et_ = KineUtils::et( e_, eta_ );
   
     bool inCone_=true;
     if( dR_out>DRout_ecal ) inCone_=false; 
     

     if( fabs(eta0_ecal) < 1.474 ) {  // Barrel num crystals, crystal width = 0.0174
       if ( fabs(deta_) < 0.0174*Jhs_ecal) continue;//inCone_=false;
       if ( dR_in < 0.0174*DRin_ecal) continue;//inCone_=false;
     } else {                       // Endcap num crystals, crystal width = 0.00864*fabs(sinh(eta))
       if ( fabs(deta_) < 0.00864*fabs(sinh(eta_))*Jhs_ecal) continue;//inCone_=false;
       if ( dR_in < 0.00864*fabs(sinh(eta_))*DRin_ecal) continue;//inCone_=false;
     }
     
     bool keep_=inCone_;
     int iecal_ = ( fabs( eta_ )<=e_EE ) ? 0 : 1;
     if(  fabs(e_)<e0_ecal[iecal_]  ) keep_=false;
     if( fabs(et_)<et0_ecal[iecal_] ) keep_=false;
    
     if( keep_ ) 
       {
	 n_ecal++;
	 S_ecal += et_;
       }
     
   } //End loop on Rechits
 }
  // 
  S = a_trk*S_trk + a_ecal*S_ecal + a_hcal*S_hcal;
  iso = pt/( pt + S );

  if( dbg_ )
    {
      cout << "S_trk  = "   << S_trk   << " /" << n_trk << endl;
      cout << "S_hcal = "   << S_hcal  << " /" << n_hcal << endl;
      cout << "S_ecal = "   << S_ecal  << " /" << n_ecal << endl;
      cout << "S      = "   << S       << endl;
      cout << "iso=" << iso << endl;
    }
  
  info->setFloat( "S_trk",  S_trk  );
  info->setFloat( "S_ecal", S_ecal );
  info->setFloat( "S_hcal", S_hcal );
  info->setFloat( "iso",       iso );
  
  return;
}

void
OfficialIsolationSelector::loadXtalDatabase() {

  ifstream database("/home/usr201/mnt/mmarionn/AnaNaS/MECore/XtalDataBase.txt");

  if(database)
    {
      int ilab;
      float eta, phi;
    
      std::pair<float, float> pairCoord;

      while(!database.eof()) {
	database >> ilab >> eta >> phi; 
	pairCoord.first = eta;
	pairCoord.second = phi;
	MapXtal[ ilab ] = pairCoord;
      }
    }
  else {
    cout<<" XtalDataBase not loaded, please check pointing directory in the Official Isolation Selector \n abort "<<endl;
    abort(); //MM
  }

  _loadDatabase = true;
  
}
