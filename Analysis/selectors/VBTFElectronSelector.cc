#include "Analysis/selectors/VBTFElectronSelector.hh"

#include <cassert>

using namespace std;

#include "MECore/src/ME.hh"
#include "MECore/src/MEEBDisplay.hh"
#include "MECore/src/MEEEDisplay.hh"

#include "Analysis/core/EventManager.hh"
#include "Analysis/core/CandInfo.hh"
#include "Analysis/utils/Constants.hh"
#include "Analysis/utils/KineUtils.hh"
#include "Analysis/utils/Config.hh"
#include "Analysis/selectors/ConeSelector.hh"

#include "TF2.h"


ClassImp( VBTFElectronSelector )

float VBTFElectronSelector::detaCorrections(float etaEle,float phiEle)
{
   TF2 
f12_fit2("f12_fit2","[0]+(TMath::TanH(y)/325.)*([1]-([2]*TMath::SinH(y)*TMath::Cos(x-[3])))",-TMath::Pi(),TMath::Pi(),-10.,10.);

   if (etaEle>1.479)
     {
       f12_fit2.SetParameter(0,0.0013);
       f12_fit2.SetParameter(1,-0.06);
       f12_fit2.SetParameter(2,0.52);
       f12_fit2.SetParameter(3,2.17);
     }
   else if (etaEle<-1.479)
     {
       f12_fit2.SetParameter(0,-0.0013);
       f12_fit2.SetParameter(1,-0.32);
       f12_fit2.SetParameter(2,0.45);
       f12_fit2.SetParameter(3,-1.58);
     }
   return f12_fit2.Eval(phiEle,etaEle);
}


float VBTFElectronSelector::dphiCorrections(float etaEle,float phiEle)
{
   TF2 
f12_fit2("f12_fit2","[0]+[1]*(TMath::SinH(y)/325.)*(TMath::Sin([2]-x))",-TMath::Pi(),TMath::Pi(),-10.,10.);

   if (etaEle>1.479)
     {
       f12_fit2.SetParameter(0,0.);
       f12_fit2.SetParameter(1,0.52);
       f12_fit2.SetParameter(2,2.17);
     }
   else if (etaEle<-1.479)
     {
       f12_fit2.FixParameter(0,0.);
       f12_fit2.FixParameter(1,0.45);
       f12_fit2.FixParameter(2,-1.58);
     }
   return f12_fit2.Eval(phiEle,etaEle);
}


string VBTFElectronSelector::level[VBTFElectronSelector::kNWP] = 
{
  "Presel", "Sel95", "Sel90", "Sel85", "Sel80", "Sel70", "Sel65", "SelAnti"
};

VBTFElectronSelector::VBTFElectronSelector( unsigned etcut, unsigned sel ) :
  LeptonSelector(etcut,sel), IsoElec( (EventManager::kElectron) )
{
  useOfficialIsoComputer = false;
  useStandardVariables = true;
  if(useOfficialIsoComputer)
    useStandardVariables =false;
}

void
VBTFElectronSelector::computeFlags( const Candidate& cand) {

  setFlags( cand );
}

bool
VBTFElectronSelector::isFlag( const Candidate& cand,  std::string isoid ){

  CandInfo* info = cand.info();
  
  bool out = false;
  out = info->getBool(isoid);
 
  return out;

}


void
VBTFElectronSelector::setFlags( const Candidate& cand ) const 
{
  //sigiEtaiEta is fixed, 0.01 for EB, 0.03 for EE.
  //Deta is fixed to 0.006 in the EB.

  //EB
  //       70     80     90     95
  //Dphi  0.02   0.02   0.04   0.8 (=no cut)
  //HOE   0.02   0.05   0.05   0.05

  //strk  2.5    3.0    6.0    7.0
  //secal 3.0    4.0    5.0    5.0
  //shcal 5.0    5.0    5.0    5.0 


  //EE
  //      70     80     90     95
  //Deta 0.003  0.006  0.008  0.008
  //Dphi 0.02   0.02   0.025  0.7 (=no cut)
  //HOE  0.0025 0.025  0.025  0.04
  //(the cut at 0.0025 means HOE=0)

  //strk  0.8    1.5    6.0    8.0
  //secal 2.5    2.5    2.5    3.0
  //shcal 0.25   0.7    1.5    2.0 


  // Tracks :
  // Cône externe 0.3
  // Cône interne 0.015
  // threshold pt 1 GeV
  
  // Ecal RecHit :
  // Cône externe 0.4
  // Barrel
  // Cône interne 0.045
  // Veto EtaPhi : eta +- 0.01

    // Endcaps
  // Cône interne 0.070
  // Veto EtaPhi : eta +- 0.02

  // Hcal Towers depth 1 :
  // Cône externe 0.4
  // Cône interne inexistant
  // Pas de seuil
  
  // Hcal Towers depth 2 :
  // Cône externe 0.4
  // Cône interne inexistant
  // Pas de seuil
  
  // Pour le Hcal il n'y a pas de seuil interne
  //  mais il y a une coupure intrinsèque en H/E que je n'ai pas encore trouvée.
  // Le "depth 1Iso" ne contient que les premières parties/couches du HCAL 
  // plus les dernières du "depth 2".
  // Le "depth 2 Iso" est complémentaire et utilise les couches profondes, 
  //exceptées les dernières du depth 2 déja utilisées.
  
  // La variable est la somme non pondérée des dépôts dans tous les cas.

  //*********

  // only consider electron candidates
  int pdgCode_ = cand.pdgCode();
  assert( fabs( pdgCode_ )==11 );

  // set all the flags to false
  CandInfo* info = cand.info();

  float eta   = info->getFloat("caloEta");
  float E     = info->getFloat("caloEnergy");
  float Et    = E/cosh(eta);
  float pt    = /*cand.pt();*/ Et; // !!!
  float phi   = cand.phi();
  float chg   = cand.charge();
  TVector3 v3 = cand.pos();
  float x0    = v3.X();
  float y0    = v3.Y();
  float z0    = v3.Z();

  //  float caloE   = info->getFloat("caloEnergy");
  float caloEta = info->getFloat("caloEta");
  float caloPhi = info->getFloat("caloPhi");
  KineUtils::adjust( caloPhi, phi );

  float dEtaIn      = info->getFloat("dEtaIn");
  float dPhiIn      = info->getFloat("dPhiIn");
  float sigiEtaiEta = info->getFloat("sigiEtaiEta"); 
  float HOE         = info->getFloat("HOE");

  // is SC cluster in fiducial region in eta?
  bool isEB = fabs(caloEta)<=1.442;
  bool isEE = fabs(caloEta)>=1.56 && fabs(caloEta)<=2.5;
  bool isEcalFiducial = isEB || isEE;


  for( unsigned isel=kVeryLoose; isel<kNSel; isel++ )
    info->setBool( selection[isel], false );
  // do not continue is SC is not in fiducial region
  info->setBool( "isEcalFiducial", isEcalFiducial );
  if(!isEcalFiducial) return;

  // OK continue
  info->setBool( "isEB", isEB );
  info->setBool( "isEE", isEE );

  bool isEt32 = Et>=32;
  bool isEt30 = Et>=30;
  bool isEt27 = Et>=27;
  bool isEt25 = Et>=25;
  bool isEt20 = Et>=20;
  bool isEt15 = Et>=15;
  bool isEt10 = Et>=10;
  bool isEt5  = Et>=5;

  // do not continue is energy of super cluster is below 20 GeV
  info->setBool( "isEt05", isEt5  );
  info->setBool( "isEt10", isEt10 );
  info->setBool( "isEt15", isEt15 );
  info->setBool( "isEt20", isEt20 );
  info->setBool( "isEt25", isEt25 );
  info->setBool( "isEt27", isEt27 );
  info->setBool( "isEt30", isEt30 );
  info->setBool( "isEt32", isEt32 );
    
  //
  // identification
  //
  bool eid60(true);
  bool eid70(true);
  bool eid80(true);
  bool eid85(true);
  bool eid90(true);
  bool eid95(true);
  bool eidAnti(true);

  bool sieie60(true);
  bool sieie70(true);
  bool sieie80(true);
  bool sieie85(true);
  bool sieie90(true);
  bool sieie95(true);

  bool dphi60(true);
  bool dphi70(true);
  bool dphi80(true);
  bool dphi85(true);
  bool dphi90(true);
  bool dphi95(true);

  bool deta60(true);
  bool deta70(true);
  bool deta80(true);
  bool deta85(true);
  bool deta90(true);
  bool deta95(true);

  bool hoe60(true);
  bool hoe70(true);
  bool hoe80(true);
  bool hoe85(true);
  bool hoe90(true);
  bool hoe95(true);
  

  { // WP 60%
    if( isEB )
      {
	// identification
	if(       sigiEtaiEta  >0.01  ) { eid60=false; sieie60=false; }
	if( fabs( dPhiIn      )>0.025 ) { eid60=false; dphi60=false;  }
	if( fabs( dEtaIn      )>0.004 ) { eid60=false; deta60=false;  }
	if(       HOE          >0.025 ) { eid60=false; hoe60=false;   }
      }
    else if( isEE )
      {
	// identification
	if(       sigiEtaiEta  >0.03   ) { eid60=false; sieie60=false; }
	if( fabs( dPhiIn      )>0.02   ) { eid60=false; dphi60=false;  }
	if( fabs( dEtaIn      )>0.005  ) { eid60=false; deta60=false;  }
	if(       HOE          >0.15  ) { eid60=false; hoe60=false;   } // 2010/2011
      }
    else
      eid60=false;
  }
  //if( !eid60 )
  {
    if( isEB )
      {
	// identification
	if(       sigiEtaiEta  >0.01  ) { eid70=false; sieie70=false; }
	if( fabs( dPhiIn      )>0.03  ) { eid70=false; dphi70=false;  }
	if( fabs( dEtaIn      )>0.004 ) { eid70=false; deta70=false;  }
	if(       HOE          >0.025 ) { eid70=false; hoe70=false;   }
      }
    else if( isEE )
      {
	// identification
	if(       sigiEtaiEta  >0.03  ) { eid70=false; sieie70=false; }
	if( fabs( dPhiIn      )>0.02 ) { eid70=false; dphi70=false;  }
	if( fabs( dEtaIn      )>0.005 ) { eid70=false; deta70=false;  }
	if(       HOE          >0.15 ) { eid70=false; hoe70=false;   }
      }
    else
      eid70=false;
  }
  // if( !eid70 ) 
  {
    if( isEB )
      {
	// identification
	if(       sigiEtaiEta  >0.01  ) { eid80=false; sieie80=false; }
	if( fabs( dPhiIn      )>0.06  ) { eid80=false; dphi80=false;  }
	if( fabs( dEtaIn      )>0.004 ) { eid80=false; deta80=false;  }
	if(       HOE          >0.04  ) { eid80=false; hoe80=false;   }
      }
    else if( isEE )
      {
	// identification
	if(       sigiEtaiEta  >0.03   ) { eid80=false; sieie80=false; }
	if( fabs( dPhiIn      )>0.03   ) { eid80=false; dphi80=false;  }
	if( fabs( dEtaIn      )>0.007  ) { eid80=false; deta80=false;  }
	if(       HOE          >0.15  ) { eid80=false; hoe80=false;   }
      }
    else
      eid80=false;
  }
// if( !eid80 ) 
  {
    if( isEB )
      {
	// identification
	if(       sigiEtaiEta  >0.01  ) { eid85=false; sieie85=false; }
	if( fabs( dPhiIn      )>0.06  ) { eid85=false; dphi85=false;  }
	if( fabs( dEtaIn      )>0.006 ) { eid85=false; deta85=false;  }
	if(       HOE          >0.04  ) { eid85=false; hoe85=false;   }
      }
    else if( isEE )
      {
	// identification
	if(       sigiEtaiEta  >0.03   ) { eid85=false; sieie85=false; }
	if( fabs( dPhiIn      )>0.04   ) { eid85=false; dphi85=false;  }
	if( fabs( dEtaIn      )>0.007  ) { eid85=false; deta85=false;  }
	if(       HOE          >0.15  ) { eid85=false; hoe85=false;   }
      }
    else
      eid85=false;
  }
  //  if( !eid85 ) 
  {
    if( isEB )
      {
	// identification
	if(       sigiEtaiEta  >0.01  ) { eid90=false; sieie90=false; }
	if( fabs( dPhiIn      )>0.8   ) { eid90=false; dphi90=false;  }
	if( fabs( dEtaIn      )>0.007 ) { eid90=false; deta90=false;  }
	if(       HOE          >0.12  ) { eid90=false; hoe90=false;   }
      }
    else if( isEE )
      {
	// identification
	if(       sigiEtaiEta  >0.03   ) { eid90=false; sieie90=false; }
	if( fabs( dPhiIn      )>0.7    ) { eid90=false; dphi90=false;  }
	if( fabs( dEtaIn      )>0.009  ) { eid90=false; deta90=false;  }
	if(       HOE          >0.15   ) { eid90=false; hoe90=false;   }
      }
    else
      eid90=false;
  }
  // if( !eid90 )
  {
    if( isEB )
      {
	// identification
	if(       sigiEtaiEta  >0.01 ) { eid95=false; sieie95=false; }
	if( fabs( dPhiIn      )>0.8   ) { eid95=false; dphi95=false;  }
	if( fabs( dEtaIn      )>0.007 ) { eid95=false; deta95=false;  }
	if(       HOE          >0.15  ) { eid95=false; hoe95=false;   }
      }
    else if( isEE )
      {
	// identification
	if(       sigiEtaiEta  >0.03  ) { eid95=false; sieie95=false; }
	if( fabs( dPhiIn      )>0.7    ) { eid95=false; dphi95=false;  }
	if( fabs( dEtaIn      )>0.01  ) { eid95=false; deta95=false;  }
	if(       HOE          >0.15   ) { eid95=false; hoe95=false;   }
      }
    else
      eid95=false;
  }
  {
    if( isEB )
      {
	// identification
	if(       sigiEtaiEta  >0.01   ) eidAnti=false;
	
	// reversed
	if( fabs( dPhiIn      )<0.04   ) eidAnti=false;
	if( fabs( dEtaIn      )<0.006  ) eidAnti=false;
	// desrever

	if(       HOE          >0.02   ) eidAnti=false;
      }
    else if( isEE )
      {
	// identification
	if(       sigiEtaiEta  >0.03   ) eidAnti=false;

	// reversed
	if( fabs( dPhiIn      )<0.025  ) eidAnti=false;
	if( fabs( dEtaIn      )<0.008  ) eidAnti=false;
	// desrever

	if(       HOE          >0.0025 ) eidAnti=false;
      }
    else
      eidAnti=false;
  }


  // set the electron-ID flags
  info->setBool( "eid60", eid60 );
  info->setBool( "eid70", eid70 );
  info->setBool( "eid80", eid80 );
  info->setBool( "eid85", eid85 );
  info->setBool( "eid90", eid90 );
  info->setBool( "eid95", eid95 );
  info->setBool( "eidAnti", eidAnti );

  info->setBool( "sieie60", sieie60 );
  info->setBool( "sieie70", sieie70 );
  info->setBool( "sieie80", sieie80 );
  info->setBool( "sieie85", sieie85 );
  info->setBool( "sieie90", sieie90 );
  info->setBool( "sieie95", sieie95 );

  info->setBool( "dphi60", dphi60 );
  info->setBool( "dphi70", dphi70 );
  info->setBool( "dphi80", dphi80 );
  info->setBool( "dphi85", dphi85 );
  info->setBool( "dphi90", dphi90 );
  info->setBool( "dphi95", dphi95 );

  info->setBool( "deta60", deta60 );
  info->setBool( "deta70", deta70 );
  info->setBool( "deta80", deta80 );
  info->setBool( "deta85", deta85 );
  info->setBool( "deta90", deta90 );
  info->setBool( "deta95", deta95 );

  info->setBool( "hoe60", hoe60 );
  info->setBool( "hoe70", hoe70 );
  info->setBool( "hoe80", hoe80 );
  info->setBool( "hoe85", hoe85 );
  info->setBool( "hoe90", hoe90 );
  info->setBool( "hoe95", hoe95 );

  //
  // isolation
  //
 
  //FIXME MM
  bool cIso=true;

  // calculate the sums
  float S_trk  = 0;
  float S_ecal = 0;
  float S_hcal = 0;

  if( useStandardVariables )
    {
      // tracks

      S_trk  += info->getFloat( "dr03TrkIso" ) / Et;
	
	// ecal
      S_ecal += info->getFloat( "dr03EcalIso" ) / Et;
	
	// hcal
      S_hcal += info->getFloat( "dr03HcalD1Iso" ) / Et;
      S_hcal += info->getFloat( "dr03HcalD2Iso" ) / Et;

	
    }
  else if( useOfficialIsoComputer )
    {
      me()->IsoElec.computeIso( cand ); 
      
      S_trk = IsoElec.S_trk / Et;
      S_ecal = IsoElec.S_ecal / Et;
      S_hcal = IsoElec.S_hcal / Et;
    }
  else
    {
      // tracker
      float pt0_trk    = 1.;
      float d0_trk     = 1.; // not in official selection
      float dz0_trk    = 5.; // not in official selection
      float dR_trk_in  = 0.015;
      float dR_trk_out = 0.30;
      float strip_trk  = 0.01;  // not in official selection
      
      // ECAL
      float et0_ecal[2], e0_ecal[2], dR_ecal_in[2], dR_ecal_out[2], strip_ecal[2];   
      // EB
      et0_ecal[0]    = 0.000;
      e0_ecal[0]     = 0.080;  // cut on energy in the barrel
      dR_ecal_in[0]  = 0.045;
      dR_ecal_out[0] = 0.40;
      strip_ecal[0]  = 0.01;   
      // EE
      et0_ecal[1]    = 0.100;  // cut on transverse energy in the endcaps
      e0_ecal[1]     = 0.000; 
      dR_ecal_in[1]  = 0.070;  
      dR_ecal_out[1] = 0.40;
      strip_ecal[1]  = 0.02;
      
      // HCAL
      float dR_hcal_in  = 0.00;  // not internal cone!
      float dR_hcal_out = 0.40;
      float et0_hcal    = 0.000; // no threshold (at this level)!
      float e0_hcal     = 0.000;

      EventManager& _e = *(EventManager::e());

      const CandList& listTracks = _e.tracks();
      int n_trk=0;
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
	      S_trk += pt_ / Et;
	    }
	}

      // ECAL location of the SC
      int iecal = ( isEB ? 0:1 );

      // ?? historical
      float etaAtEcal = caloEta;
      float phiAtEcal = caloPhi;
      float eta0_ecal = etaAtEcal;
      float phi0_ecal = phiAtEcal;
   
      size_t nRecHits = (size_t) _e.a().n("rh");
      int n_ecal=0;
      for( size_t irh=0; irh<nRecHits; irh++ ) 
	{
	  assert( _e.a().load("rh", irh) );
     
	  float e_  = _e.a().get_f("rh","energy");
	  int ieta_ = _e.a().get_i("rh","ix");
	  int iphi_ = _e.a().get_i("rh","iy");
	  int iz_   = _e.a().get_i("rh","iz");

	  int iecal_(0);
	  float eta_(0.), phi_(0.);
      
	  if(iz_==0)
	    {
	      iecal_ = 0;
	      MEEBGeom::EtaPhiPoint point_ = MEEBDisplay::center( ieta_, iphi_ );
	      eta_ = point_.first;
	      phi_ = point_.second * Constants::pi;
	    }
	  else
	    {
	      iecal_ = 1;
	      MEEEGeom::EtaPhiPoint point_ = MEEEDisplay::center( ieta_, iphi_, iz_ );	  
	      eta_ = point_.first;
	      phi_ = point_.second * Constants::pi;
	    }

	  float et_  = KineUtils::et( e_, eta_ );
	  
	  float deta_ =  eta_ -eta0_ecal;
	  float dR_out  = KineUtils::dR( etaAtEcal, eta_, phiAtEcal, phi_ );
	  float dR_in   = KineUtils::dR( eta0_ecal, eta_, phi0_ecal, phi_ );

	  bool inCone_=true;
	  if( dR_out > dR_ecal_out[iecal] )       inCone_=false;
	  if( dR_in  < dR_ecal_in[iecal]  )       inCone_=false;
	  if( fabs( deta_ ) < strip_ecal[iecal] ) inCone_=false;
	  
	  bool keep_=inCone_;
	  
	  // apply thresholds corresponding to the hit location
	  if(  e_<e0_ecal[iecal_]  ) keep_=false;
	  if( et_<et0_ecal[iecal_] ) keep_=false;
	  
	  if( keep_ ) 
	    {
	      n_ecal++;
	      S_ecal += et_ / Et;
	    }
	} //End loop on Rechits

      float r_hcal = 167;
      float z_hcal = 350;
      float BField = 3.8;

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
      int n_hcal = 0;
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
	      S_hcal += hadEt_ / Et;
	    }
	}
    }

  // store the sums
  info->setFloat( "sumTrk",  S_trk  );
  info->setFloat( "sumEcal", S_ecal );
  info->setFloat( "sumHcal", S_hcal );

  //Combined isolation
  EventManager& _e = *(EventManager::e());
  float rhoFJ = _e.getRhoFastJet();
  float combIso = (S_trk + S_ecal + S_hcal) - rhoFJ*0.22/Et ; 
  info->setFloat( "combIso", combIso);


  bool iso60(true);
  bool iso70(true);
  bool iso80(true);
  bool iso85(true);
  bool iso90(true);
  bool iso95(true);
  bool isoAnti(true);

  bool trk60(true);
  bool trk70(true);
  bool trk80(true);
  bool trk85(true);
  bool trk90(true);
  bool trk95(true);

  bool ecal60(true);
  bool ecal70(true);
  bool ecal80(true);
  bool ecal85(true);
  bool ecal90(true);
  bool ecal95(true);
  
  bool hcal60(true);
  bool hcal70(true);
  bool hcal80(true);
  bool hcal85(true);
  bool hcal90(true);
  bool hcal95(true);

  bool comb60(true);
  bool comb70(true);
  bool comb80(true);
  bool comb85(true);
  bool comb90(true);
  bool comb95(true);
  
  {
    if( isEB )
      {
	// isolation
	if( S_trk >0.04 ) { iso60=false; trk60=false; }
	if( S_ecal>0.04 ) { iso60=false; ecal60=false;}
	if( S_hcal>0.03 ) { iso60=false; hcal60=false;}

	if(cIso) iso60=true;
	if( combIso > 0.03 ) { iso60=false; comb60=false;  }

      }
    else if( isEE )
      {
	// isolation
	if( S_trk >0.025  ) { iso60=false; trk60=false; }
	if( S_ecal>0.02  )  { iso60=false; ecal60=false;}
	if( S_hcal>0.02 )   { iso60=false; hcal60=false;}
	
	if(cIso) iso60=true;
	if( combIso > 0.02 ) { iso60=false; comb60=false;  }

      }
  }
  // if( !iso60 )
  {
    if( isEB )
      {
	// isolation
	if( S_trk > 0.05 )  { iso70=false; trk70=false; }
	if( S_ecal > 0.06 ) { iso70=false; ecal70=false;}
	if( S_hcal > 0.03 ) { iso70=false; hcal70=false;}

	if(cIso) iso70=true;
	if( combIso > 0.05 ) { iso70=false; comb70=false;  }

      }
    else if( isEE )
      {
	// isolation
	if( S_trk > 0.025 )  { iso70=false; trk70=false; }
	if( S_ecal > 0.025 ) { iso70=false; ecal70=false;}
	if( S_hcal > 0.02 )  { iso70=false; hcal70=false;}
	
     	if(cIso) iso70=true;
	if( combIso > 0.04 ) { iso70=false; comb70=false;  }
	
      }
  }
  // if( !iso70 ) 
    {
    if( isEB )
      {
	// isolation
	if( S_trk >0.09 ) { iso80=false; trk80=false; }
	if( S_ecal>0.07 ) { iso80=false; ecal80=false;}
	if( S_hcal>0.10 ) { iso80=false; hcal80=false;}

     	if(cIso) iso80=true;
	if( combIso > 0.07 ) { iso80=false; comb80=false;  }

      }
    else if( isEE )
      {
	// isolation
	if( S_trk >0.04  )  { iso80=false; trk80=false; }
	if( S_ecal>0.05  )  { iso80=false; ecal80=false;}
	if( S_hcal>0.025 )  { iso80=false; hcal80=false;}
	
     	if(cIso) iso80=true;
	if( combIso > 0.06 ) { iso80=false; comb80=false;  }

      }
  }
  // if( !iso80 )
    {
    if( isEB )
      {
	// isolation 
	if( S_trk >0.09 ) { iso85=false; trk85=false; }
	if( S_ecal>0.08 ) { iso85=false; ecal85=false;}
	if( S_hcal>0.10 ) { iso85=false; hcal85=false;}
     
     	if(cIso) iso85=true;
	if( combIso > 0.09 ) { iso85=false; comb85=false;  }

      }
    else if( isEE )
      {
	// isolation
	if( S_trk >0.05  ) { iso85=false; trk85=false; }
	if( S_ecal>0.05  ) { iso85=false; ecal85=false;} 
	if( S_hcal>0.025 ) { iso85=false; hcal85=false;}

     	if(cIso) iso85=true;
	if( combIso > 0.07 ) { iso85=false; comb85=false;  }

      }
    }
// if( !iso85 ) 
    {
    if( isEB )
      {
	// isolation
	if( S_trk >0.12 ) { iso90=false; trk90=false; }
	if( S_ecal>0.09 ) { iso90=false; ecal90=false;}
	if( S_hcal>0.10 ) { iso90=false; hcal90=false;}

     	if(cIso) iso90=true;
	if( combIso > 0.15 ) { iso90=false; comb90=false;  }

      }
    else if( isEE )
      {
	// isolation
	if( S_trk >0.05 ) { iso90=false; trk90=false; }
	if( S_ecal>0.06 ) { iso90=false; ecal90=false;}
	if( S_hcal>0.03 ) { iso90=false; hcal90=false;}
	
     	if(cIso) iso90=true;
	if( combIso > 0.1 ) { iso90=false; comb90=false;  }

      }
  }
// if( !iso90 ) 
  {
    if( isEB )
      {
	// isolation
	if( S_trk >0.15 ) { iso95=false; trk95=false; }
	if( S_ecal>2.00 ) { iso95=false; ecal95=false;}
	if( S_hcal>0.12 ) { iso95=false; hcal95=false;}

     	if(cIso) iso95=true;
	if( combIso > 0.2) { iso95=false; comb95=false;  }

      }
    else if( isEE )
      {
	// isolation
	if( S_trk >0.08 ) { iso95=false; trk95=false; }
	if( S_ecal>0.06 ) { iso95=false; ecal95=false;}
	if( S_hcal>0.05 ) { iso95=false; hcal95=false;}
	
	if(cIso) iso95=true;
	if( combIso > 0.15 ) { iso95=false; comb95=false;  }

      }
  }
  {
    if( isEB )
      {
	// isolation
	if( S_trk >2.5 ) isoAnti=false;
	if( S_ecal>3.0 ) isoAnti=false;
	if( S_hcal>5.  ) isoAnti=false;
      }
    else if( isEE )
      {
	// isolation
	if( S_trk >0.8  ) isoAnti=false;
	if( S_ecal>2.5  ) isoAnti=false;
	if( S_hcal>0.25 ) isoAnti=false;
      }
  }
  
  // set the isolation flags
  info->setFloat( "sumTrk",  S_trk  );
  info->setFloat( "sumEcal", S_ecal );
  info->setFloat( "sumHcal", S_hcal );
  info->setFloat( "iso", S_trk+S_ecal+S_hcal ); 
  info->setFloat( "combIso", combIso ); 
  info->setBool( "iso60", iso60 );
  info->setBool( "iso70", iso70 );
  info->setBool( "iso80", iso80 );
  info->setBool( "iso85", iso85 );
  info->setBool( "iso90", iso90 );
  info->setBool( "iso95", iso95 );
  info->setBool( "isoAnti", isoAnti );
  
  info->setBool( "trk60", trk60 );
  info->setBool( "trk70", trk70 );
  info->setBool( "trk80", trk80 );
  info->setBool( "trk85", trk85 );
  info->setBool( "trk90", trk90 );
  info->setBool( "trk95", trk95 );

  info->setBool( "ecal60", ecal60 );
  info->setBool( "ecal70", ecal70 );
  info->setBool( "ecal80", ecal80 );
  info->setBool( "ecal85", ecal85 );
  info->setBool( "ecal90", ecal90 );
  info->setBool( "ecal95", ecal95 );

  info->setBool( "hcal60", hcal60 );
  info->setBool( "hcal70", hcal70 );
  info->setBool( "hcal80", hcal80 );
  info->setBool( "hcal85", hcal85 );
  info->setBool( "hcal90", hcal90 );
  info->setBool( "hcal95", hcal95 );

  info->setBool( "comb60", comb60 );
  info->setBool( "comb70", comb70 );
  info->setBool( "comb80", comb80 );
  info->setBool( "comb85", comb85 );
  info->setBool( "comb90", comb90 );
  info->setBool( "comb95", comb95 );


  // set the selection flags
  bool sel60 = eid60 && iso60;
  bool sel70 = eid70 && iso70;
  bool sel80 = eid80 && iso80;
  bool sel85 = eid85 && iso85;
  bool sel90 = eid90 && iso90;
  bool sel95 = eid95 && iso95;
  bool selAnti = eidAnti && isoAnti;

  info->setBool( "sel60", sel60 );
  info->setBool( "sel70", sel70 );
  info->setBool( "sel80", sel80 );
  info->setBool( "sel85", sel85 );
  info->setBool( "sel90", sel90 );
  info->setBool( "sel95", sel95 );
  info->setBool( "selAnti", selAnti );

  bool isAboveEtThreshold = false;
  if( _settings==kEt5 )
    {
      if( info->getBool( "isEt05" ) ) isAboveEtThreshold = true; 
    }
  else if( _settings==kEt10 )
    {
      if( info->getBool( "isEt10" ) ) isAboveEtThreshold = true; 
    }
  else if( _settings==kEt15 )
    {
      if( info->getBool( "isEt15" ) ) isAboveEtThreshold = true; 
    }
  else if( _settings==kEt20 )
    {
      if( info->getBool( "isEt20" ) ) isAboveEtThreshold = true; 
    }
  else if( _settings==kEt25 )
    {
      if( info->getBool( "isEt25" ) ) isAboveEtThreshold = true; 
    }
  else if( _settings==kEt27 )
    {
      if( info->getBool( "isEt27" ) ) isAboveEtThreshold = true; 
    }
  else if( _settings==kEt30 )
    {
      if( info->getBool( "isEt30" ) ) isAboveEtThreshold = true; 
    }
  else if( _settings==kEt32 )
    {
      if( info->getBool( "isEt32" ) ) isAboveEtThreshold = true; 
    }
  else
    {
      abort();
    }

  bool selectAll = false;
  if( selectAll )
    {
      info->setBool(  selection[kTight], true );
      info->setBool(  selection[kMedium], true );
      info->setBool(  selection[kLoose], true );
      info->setBool(  selection[kVeryLoose], true );
      
      return;
    }

  if( isAboveEtThreshold )
    {
      sel95=true;

      if( sel95 ) 
	{
	  info->setBool(  selection[kVeryLoose], true );
	  if( sel90 )
	    {
	      info->setBool(  selection[kLoose], true );	  
	      if( sel80 )
		{
		  info->setBool(  selection[kMedium], true );	  
		  if( sel70 )
		    {
		      info->setBool(  selection[kTight], true );	  
		    }
		}
	    }
	}
    }

  return;
}
