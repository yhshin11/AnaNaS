#include "Analysis/selectors/IsolationSelector.hh"

#include <cassert>

#include "Analysis/utils/KineUtils.hh"
#include "Analysis/selectors/ConeSelector.hh"

using namespace std;

// Pt selector
ClassImp( IsolationSelector )

IsolationSelector::IsolationSelector( unsigned settings_, float iso0_ ) : _settings(settings_), iso0(iso0_)
{  
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
      et0_ecal    = 0.200;
      
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
      pt0_trk    = 1.;
      d0_trk  = 1.;
      dz0_trk = 5.;
      dR_trk_in  = 0.015;
      dR_trk_out = 0.30;
      strip_trk  = 0.01;
      
      ecalUseIsoDeposits  = false; //MM By default, use RecHits

      et0_ecal    = 0.200;

      e0_ecal[0]     = 0.120;
      dR_ecal_in[0]  = 0.045;
      dR_ecal_out[0] = 0.40;
      strip_ecal[0]  = 0.02;  

      e0_ecal[1]     = 0.450; 
      dR_ecal_in[1]  = 0.070;
      dR_ecal_out[1] = 0.40;
      strip_ecal[1]  = 0.02;
      
      dR_hcal_in  = 0.10;
      dR_hcal_out = 0.40;
      et0_hcal    = 0.500;
      e0_hcal     = 0.600;
    }
  else
    assert(0);

  a_trk  = 1;
  a_ecal = 1;
  a_hcal = 1;

  iso0 = 0.;

}

bool 
IsolationSelector::accept( const Candidate& cand ) const
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
IsolationSelector::computeIso( const Candidate& cand )
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
  if( _settings==EventManager::kElectron )
    {
      assert( isElectron );
      mapname_ = "electron-EcalIsoDep";      
    }
  else if( _settings==EventManager::kMuon )
    {
      assert( isMuon );
      mapname_ = "muon-EcalIsoDep";
    }
  
  const CandList& listTracks = _e.tracks();
  for( size_t ii=0; ii<listTracks.size(); ii++ )
    {
      const Candidate* trk = listTracks[ii];
      float pt_  = trk->pt();
      float eta_ = trk->eta();
      float phi_ = trk->phi();
      float deta_ = eta_-eta;
      float dR_=KineUtils::dR( eta, eta_, phi, phi_ );
      bool inCone_=true;
      if( dR_<dR_trk_in  )          inCone_=false;
      if( dR_>dR_trk_out )          inCone_=false;
      if( fabs( deta_ )<strip_trk ) inCone_=false;
      
      TVector3 v3_ = trk->pos();
      float x0_  = v3_.X();
      float y0_  = v3_.Y();
      float z0_  = v3_.Z();
      float r0_ = sqrt( x0_*x0_ + y0_*y0_ );
      
      bool keep_=inCone_;
      if( pt_<pt0_trk )        keep_=false; 
      if( r0_>d0_trk )         keep_=false;
      if( fabs(z0_)>dz0_trk )  keep_=false;
      
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
  if( isElectron )
    {
      float caloEta = info->getFloat("caloEta");
      float caloPhi = info->getFloat("caloPhi");
      KineUtils::adjust( caloPhi, phi );
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
	 float dR_out  = KineUtils::dR( etaAtEcal, eta_, phiAtEcal, phi_ );
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
	 if( et_<et0_ecal         ) keep_=false;
	 
	 if( keep_ ) 
	   {
	     n_ecal++;
	     S_ecal += et_;
	   }
       }
   }
 else { //MM Use RecHits
   
   size_t nRecHits = (size_t) _e.a().n("rh");
   for(size_t irh=0;irh<nRecHits;irh++) {
     assert( _e.a().load("rh", irh) );
     
     float e_ = _e.a().get_f("rh","energy");
     int ieta_ = _e.a().get_i("rh","ix");
     int iphi_ = _e.a().get_i("rh","iy");
     int iz_ = _e.a().get_i("rh","iz");
     
     if(iz_==0)
    {
      MEEBGeom::EtaPhiPoint point_ = MEEBDisplay::center( ieta_, iphi_ );
      
      float Eta = point_.first;
      float Phi = point_.second * Constants::pi;

      float et_  = KineUtils::et( e_, point_.first );
     
      float deta_ =  Eta -eta0_ecal;
      float dR_out  = KineUtils::dR( etaAtEcal, Eta, phiAtEcal, Phi );
      float dR_in   = KineUtils::dR( eta0_ecal, Eta, phi0_ecal, Phi );
      bool inCone_=true;
      if( dR_out>DRout_ecal ) inCone_=false;
      if( dR_in <DRin_ecal  ) inCone_=false;
      if( fabs( deta_ )<Jhs_ecal ) inCone_=false;

      bool keep_=inCone_;
      int iecal_ = ( fabs( Eta )<=e_EE ) ? 0 : 1;
       if(  e_<e0_ecal[iecal_]  ) keep_=false;
       if( et_<et0_ecal         ) keep_=false;
      
      if( keep_ ) 
	{
	  n_ecal++;
	  S_ecal += et_;
	}
    }
  else //Endcaps
    {
      MEEEGeom::EtaPhiPoint point_ = MEEEDisplay::center( ieta_, iphi_, iz_ );

      float Eta = point_.first;
      float Phi = point_.second * Constants::pi;

      float et_  = KineUtils::et( e_, point_.first );
      float deta_ =  Eta -eta0_ecal;
      float dR_out  = KineUtils::dR( etaAtEcal, Eta, phiAtEcal, Phi );
      float dR_in   = KineUtils::dR( eta0_ecal, Eta, phi0_ecal, Phi );
      bool inCone_=true;
      if( dR_out>DRout_ecal ) inCone_=false;
      if( dR_in <DRin_ecal  ) inCone_=false;
      if( fabs( deta_ )<Jhs_ecal ) inCone_=false;

      bool keep_=inCone_;
      int iecal_ = ( fabs( Eta )<=e_EE ) ? 0 : 1;
      if(  e_<e0_ecal[iecal_]  ) keep_=false;
      if( et_<et0_ecal         ) keep_=false;
      
      if( keep_ ) 
	{
	  n_ecal++;
	  S_ecal += et_;
	}
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

