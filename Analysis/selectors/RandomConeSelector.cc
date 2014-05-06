#include "Analysis/selectors/RandomConeSelector.hh"

#include <cassert>

#include "Analysis/utils/KineUtils.hh"
#include "Analysis/selectors/ConeSelector.hh"
#include "Math/GenVector/VectorUtil.h"
#include <fstream>
#include <iostream>
#include <stdlib.h>

using namespace std;

// Pt selector
ClassImp( RandomConeSelector )

RandomConeSelector::RandomConeSelector( vector<int> interest )
{  
  debug = false;
  
  _dbLoaded = false;
  absCut = false;
  useGeometryDatabase = true;
  doSCcut = true;
  
  RHinterest = interest;
  niceTkCut = 2.;
  
  e_EE = 1.479;
  r_HB = 167;
  z_HE = 350;

  // CMSSW/RecoEgamma/PhotonIdentification/python/isolationCalculator_cfi.py
  
  // track
  pt0_trk    = 0.000; // if 0 all tracks counted
  d0_trk     = 0.100; 
  dz0_trk    = 0.200;
  dR_trk_in  = 0.040;
  dR_trk_out = 0.400;
  strip_trk  = 0.015;
      
  // ecal barrel
  et0_ecal[0]    = 0.280; 
  e0_ecal[0]     = -999.; // ZS 
  dR_ecal_in[0]  = 3.5; // xtal counts
  dR_ecal_out[0] = 0.400;
  strip_ecal[0]  = 2.5; // xtal counts

  // ecal endcap
  et0_ecal[1]    = 0.280;
  e0_ecal[1]     = -999.; 
  dR_ecal_in[1]  = 3.5; // xtal counts
  dR_ecal_out[1] = 0.400;
  strip_ecal[1]  = 2.5; // xtal counts
      
  // hcal
  dR_hcal_in  = 0.150;
  dR_hcal_out = 0.400;
  et0_hcal    = 0.000;
  e0_hcal     = 0.000;

}

void 
RandomConeSelector::computeIso( const Candidate& cand, const Candidate& candSC , bool ptk)
{
  bool dbg_ = false;

  float pt  = cand.pt();
  assert( pt>0 );
  
  // compute the three isolations
  computeTrackIso ( cand );
  computeHcalIso  ( cand );
  computeHcalIsoH ( cand );
  computeEcalIso  ( cand, candSC, ptk );

  
  if( dbg_ ) {
    cout << "S_trk  = "   << S_trk   << " / " << n_trk << endl;
    cout << "S_hcal = "   << S_hcal  << " / " << n_hcal << endl;
    cout << "S_ecal = "   << S_ecal  << " / " << n_ecal << endl;
  }

  S_sum = S_ecal + S_hcal + S_trk;

  return;
}

void
RandomConeSelector::computeTrackIso( const Candidate& cand ) {
  
  S_trk = 0;
  n_trk = 0;
  EventManager& _e = *(EventManager::e());
  const CandList& listTracks = _e.tracks();
  

  for( size_t ii=0; ii<listTracks.size(); ii++ ) {
  
    const Candidate* trk = listTracks[ii];
    float pt_  = trk->pt();
    float deta_ = trk->eta() - cand.eta();
    
    float dR_= ROOT::Math::VectorUtil::DeltaR(trk->p3(), cand.p3() ); 
    if( dR_<dR_trk_in  )        {if(debug)cout << "drin" << endl; continue;} //inCone_=false;
    if( dR_>dR_trk_out )        {if(1==0)cout << "drout " << dR_ << " pt " << pt_ << endl; continue;} //inCone_=false;
    if( fabs(deta_)<strip_trk ) {if(debug)cout << "strip" << endl; continue;} //inCone_=false;
      
    TVector3 v3_ = trk->pos();
    float x0_ = v3_.X()-cand.vertex()->pos().X();
    float y0_ = v3_.Y()-cand.vertex()->pos().Y();
    float z0_ = v3_.Z()-cand.vertex()->pos().Z();
    float r0_ = sqrt( x0_*x0_ + y0_*y0_ );
      
    if( pt_       < pt0_trk ) {if(debug)cout << "pt" << endl; continue;} //keep_=false;
    if( r0_       > d0_trk  ) {if(debug && 1==0)cout << "d0 " << r0_ << endl; continue;} //keep_=false;
    if( fabs(z0_) > dz0_trk ) {if(debug)cout << "dz0" << endl; continue;} //keep_=false;
      
    n_trk++;
    S_trk += pt_;
  }
}

void RandomConeSelector::computeHcalIsoH( const Candidate& cand ) {
  S_hcalH = 0;
  n_hcalH = 0;
  EventManager& _e = *(EventManager::e());
  float eta = cand.eta();
  float phi = cand.phi();
  TVector3 v3 = cand.pos();
  float x0  = v3.X();
  float y0  = v3.Y();
  float z0  = v3.Z();

  float r_hcal = r_HB;
  float z_hcal = z_HE;

  const CandList& listCT = _e.caloTowers();

  // preparing the vector of track impacts to ignore
  vector<float> etaHole, phiHole;
  for( int tk=0; tk<(int)vtkpt.size(); tk++ ) {
    if( vtkpt[tk]<niceTkCut ) continue;
    if( vtkFP[tk]==false )    continue; // check if the track is "good"
    etaHole.push_back(vtkEtaH[tk]);
    phiHole.push_back(vtkPhiH[tk]);
  }
  // if( etaHole.size()>0 ) cout << "etaHole size " << etaHole.size() << " nTow " << listCT.size() << endl;

  float etaAtHcal(0);
  float phiAtHcal(0);
  KineUtils::straightToEnvelop( etaAtHcal, phiAtHcal, eta, phi, x0, y0, z0, r_hcal, z_hcal );
  KineUtils::adjust( phiAtHcal, phi );

  float DRin_hcal  = dR_hcal_in;
  float DRout_hcal = dR_hcal_out;

  for( size_t iCT=0; iCT<listCT.size(); iCT++ ) {
    const Candidate& ct = *listCT[iCT];
    const CandInfo* info_ = ct.info();
    float eta_   = ct.eta();
    float phi_   = ct.phi();
    float hadEt_ = info_->getFloat("hadEt");
    if( hadEt_==0 ) continue;

    float hadE_ = KineUtils::e( hadEt_, eta_ ); // compute energy from et
    float dR_ = KineUtils::dR( etaAtHcal, eta_, phiAtHcal, phi_ );
    if( dR_ > DRout_hcal ) continue; 
    if( dR_ < DRin_hcal  ) continue; 

    // exclude the towers right behind a "good" track impact
    bool exclude = false;
    for( int tk=0; tk<(int)etaHole.size(); tk++) {
      dR_ = KineUtils::dR( etaHole[tk], eta_, phiHole[tk], phi_ );
      //      cout << "dR to impact " << dR_ << "         etaHole " << etaHole[tk] << " etaTow " << eta_ << " etaSC " << eta  << " phiHole " << phiHole[tk] << " phiTow " << phi_ << " phiSC " << phi << endl;
//       if( dR_<4*0.0174 ) { exclude = true; }
      if( dR_<4*0.0174 ) { exclude = true; }
    }
    if( exclude==true )   continue;
    if( hadEt_<et0_hcal ) continue;
    if( hadE_ <e0_hcal  ) continue;

    n_hcalH++;
    S_hcalH += hadEt_;
  }
  //  if( etaHole.size()>0 ) { cout << n_hcalH << endl; cout << endl; }
}

void
RandomConeSelector::computeHcalIso( const Candidate& cand ) {

  S_hcal = 0;
  n_hcal = 0;
  EventManager& _e = *(EventManager::e());
  float eta = cand.eta();
  float phi = cand.phi();
  TVector3 v3 = cand.pos();
  float x0  = v3.X();
  float y0  = v3.Y();
  float z0  = v3.Z();
  
  float r_hcal = r_HB;
  float z_hcal = z_HE;
  
  float etaAtHcal(0);
  float phiAtHcal(0);
  KineUtils::straightToEnvelop( etaAtHcal, phiAtHcal, eta, phi, x0, y0, z0, r_hcal, z_hcal );
  KineUtils::adjust( phiAtHcal, phi );

  float DRin_hcal  = dR_hcal_in;
  float DRout_hcal = dR_hcal_out;

  const CandList& listCT = _e.caloTowers();

  for( size_t iCT=0; iCT<listCT.size(); iCT++ ) {
    const Candidate& ct = *listCT[iCT];
    const CandInfo* info_ = ct.info();
    float eta_   = ct.eta();
    float phi_   = ct.phi();
    float hadEt_ = info_->getFloat("hadEt");
    if( hadEt_==0 ) continue;

    float hadE_ = KineUtils::e( hadEt_, eta_ ); // compute energy from et

    float dR_ = KineUtils::dR( etaAtHcal, eta_, phiAtHcal, phi_ ); 
    bool inCone_=true;
    if( dR_ > DRout_hcal ) continue; //inCone_=false;
    if( dR_ < DRin_hcal  ) continue; //inCone_=false;

    bool keep_=inCone_;
    if( hadEt_<et0_hcal ) continue; //keep_=false; //     et cut
    if( hadE_ <e0_hcal  ) continue; //keep_=false; // energy cut

    if( keep_ ) {
      n_hcal++;
      S_hcal += hadEt_;
    }
  }
  
}

void
RandomConeSelector::computeEcalIso( const Candidate& cand, const Candidate& candSC, bool ptk) {
  
  usedTrkMatch.clear();
  S_ecal = 0;
  n_ecal = 0;
  EventManager& _e = *(EventManager::e());
  ieta.clear();
  iphi.clear();
  veta.clear();
  vphi.clear();
  iz.clear();
  ve.clear();
  vet.clear();
  
  float etaSC_ = candSC.eta(); if( debug )cout << "etaSC " << etaSC_ << " etaPHO " << cand.eta() << endl;
  float phiSC_ = candSC.phi(); if( debug )cout << "phiSC " << phiSC_ << " phiPHO " << cand.phi() << endl;
  
  iecal = (fabs(etaSC_)<=e_EE)? 0 : 1;
  float DRin_ecal  = dR_ecal_in[iecal];
  float DRout_ecal = dR_ecal_out[iecal];
  float Strip_ecal = strip_ecal[iecal];
  
  if( useGeometryDatabase && !_dbLoaded ) loadXtalDatabase();
  
  for(size_t irh=0;irh<RHinterest.size();irh++) {
    if( debug ) cout << "RHinterest loop start." << endl;
    assert( _e.a().load("rh", RHinterest[irh]) );
    
    float e_    = _e.a().get_f("rh","energy");
    int   ieta_ = _e.a().get_i("rh","ix");
    int   iphi_ = _e.a().get_i("rh","iy");
    int   iz_   = _e.a().get_i("rh","iz");

    float etaRH_;
    float phiRH_;
    
    if( useGeometryDatabase ) {
      int ilab =  10000*(85+ieta_) + 10*iphi_ + iz_ + 1;
      itMapXtal = MapXtal.find( ilab );
      if( itMapXtal == MapXtal.end() )  { 
        cout<<" Xtal Missing in database :  " << ilab << endl; 
        abort(); 
      }
      else
      {
        etaRH_ = (*itMapXtal).second.first;
        phiRH_ = (*itMapXtal).second.second;
      }
    }
    else {
      if(iz_==0) { //Barrel
        MEEBGeom::EtaPhiPoint point_ = MEEBDisplay::center( ieta_, iphi_ );
        etaRH_ = point_.first;
        phiRH_ = point_.second * Constants::pi;
      }
      else { //Endcaps
        MEEEGeom::EtaPhiPoint point_ = MEEEDisplay::center( ieta_, iphi_, iz_ );
        etaRH_ = point_.first;
        phiRH_ = point_.second * Constants::pi;
      }
    }
    if( debug ) cout << "Xtal position found." << endl;
    float et_    = KineUtils::et( e_, etaRH_ ); // rechit et
    float deta_ = etaRH_ - etaSC_;
    float dR_   = KineUtils::dR( etaSC_, etaRH_, phiSC_, phiRH_ );

    if( dR_ > DRout_ecal ) { continue; } 


    float xSizeEB = 0.0174;
    float xSizeEE = 0.00864;
    if( fabs(etaSC_) < e_EE ) {  
      if ( fabs(deta_) < Strip_ecal*xSizeEB ) { continue; } // strip veto
      if ( dR_         < DRin_ecal*xSizeEB) { continue; }  // inner cone veto
    }
    else {
      if ( fabs(deta_) < fabs(sinh(etaRH_))*Strip_ecal*xSizeEE ) {continue; } // strip veto
      if ( dR_         < fabs(sinh(etaRH_))*DRin_ecal*xSizeEE) {continue; } // inner cone veto
    }
    
    // apply a cut on e and/or et
    if (absCut) {
      if(  fabs(e_) < e0_ecal [iecal] ) { continue; } // energy cut
      if( fabs(et_) < et0_ecal[iecal] ) { continue; } // et cut
    }
    else {
      if( e_  < e0_ecal[iecal]  ) { continue; } // energy cut
      if( et_ < et0_ecal[iecal] ) { continue; } // et cut
    }
    if( debug ) cout << "Simple eliminations done." << endl;
    
    // do not keep rechits that are part of the supercluster
    bool isInSC = false;
    const CandMap&         Sc_Bc_ = _e.candMap("SuperCluster-BasicCluster");
    const CandIdRecHitMap& Bc_rh_ = _e.ecalRecHitMap("BasicCluster-EcalRecHit");
    Candidate* cSC = const_cast<Candidate*> (&candSC);
    CandMapIterator itbc_lower_ = Sc_Bc_.lower_bound( cSC );
    CandMapIterator itbc_upper_ = Sc_Bc_.upper_bound( cSC );
    CandMapIterator itbc;
    for (itbc=itbc_lower_; itbc!=itbc_upper_; itbc++) {
      CandIdRecHitMapIterator itrh_lower_ = Bc_rh_.lower_bound( (*itbc).second->uid() );
      CandIdRecHitMapIterator itrh_upper_ = Bc_rh_.upper_bound( (*itbc).second->uid() );
      CandIdRecHitMapIterator itrh;
      for (itrh=itrh_lower_; itrh!=itrh_upper_; itrh++) {
        ecalRecHit rh = (*itrh).second.first;
        if (rh.ix==ieta_ && rh.iy==iphi_ && rh.iz==iz_) {
          isInSC = true;
        }
      }
    }
    if (doSCcut && isInSC) continue;
    if( debug ) cout << "SC xtals kicked out." << endl;
    
    // do not keep rechits in other ecal parts - obsolete from 20/01/2011
    //if ( abs(iz_) != iecal ) continue; 
    
    // eventually, add the et to the sum
    n_ecal++;
    S_ecal += et_;
    ieta.push_back(ieta_);
    iphi.push_back(iphi_);
    iz.push_back(iz_);
    veta.push_back(etaRH_);
    vphi.push_back(phiRH_);
    ve.push_back(e_);
    vet.push_back(et_);
  }
  if( debug ) cout << "InterestRH loop done." << endl;
  
  vtkpt.clear();
  vtkCnt.clear();
  vtkd0.clear();
  vtkz0.clear();
  vtkNhits.clear();
  vtkNxIhits.clear();
  vtkNdof.clear();
  vtkQ.clear();
  vtkChi2.clear();
  vtkFP.clear();
  vtkSV.clear();
  vtkEtaH.clear();
  vtkPhiH.clear();
  vtkDEta.clear();
  vtkDPhi.clear();
  vtkDR.clear();
  vtkInVeto.clear();
  vtkVid.clear();
  vtkDepE.clear();
  vtkDepEt.clear();
  n_tkecal = 0;
  
  const CandList& listTracks = _e.tracks();
  for( size_t ii=0; ii<listTracks.size(); ii++ ) {
  
    const Candidate* trk = listTracks[ii];
    float pt_  = trk->pt();
    float hiteta, hitphi;
    
    TVector3 v3_ = trk->pos();
//     float x0_ = cand.vertex()->pos().X(); //v3_.X();
    float x0_ = v3_.X();
    float y0_ = v3_.Y();
    float z0_ = v3_.Z();
    float r0_ = sqrt( x0_*x0_ + y0_*y0_ );
    
    float r_ecal = (fabs(etaSC_)<e_EE)? 129:150;
    float z_ecal = (fabs(etaSC_)<e_EE)? 280:325;
    
    int qtk_ = (int) trk->charge(); 
    KineUtils::helixToEnvelop( hiteta, hitphi, qtk_*pt_, trk->eta(), trk->phi(), x0_, y0_, z0_, r_ecal, z_ecal, 3.8);
    
    float deta_ = hiteta - etaSC_;
    const CandInfo* infot_ = trk->info();

    float pi = TMath::Pi();
    float cmsEta = infot_->getFloat("impEta");
    float cmsPhi = infot_->getFloat("impPhi");
    if( cmsPhi<0 ) cmsPhi += 2*pi;

    //cout << "impactEta     hand " << hiteta << " cmssw " << cmsEta << "     diff " << hiteta-cmsEta << endl;
    //cout << "impactPhi     hand " << hitphi << " cmssw " << cmsPhi << "     diff " << hitphi-cmsPhi << endl;

    isInVeto = 0;
    float xSize = (fabs(etaSC_)<e_EE)? 0.0174:0.00864;

    //float dR_   = KineUtils::dR( hiteta, etaSC_, hitphi, phiSC_ );
    float dR_   = KineUtils::dR( cmsEta, etaSC_, cmsPhi, phiSC_ );

    float dphi_ = TMath::Sqrt(dR_*dR_ - deta_*deta_);
    if( dR_        > DRout_ecal       )      { continue; }
    
    if( dR_<(DRin_ecal)*xSize || fabs(deta_)<(Strip_ecal)*xSize ) { isInVeto = true; }
    if( dR_>(DRin_ecal)*xSize && fabs(deta_)>(Strip_ecal)*xSize ) { isInVeto = false; }

    // impact parameters
    bool samevtx = false ;
    float xip_  = v3_.X()-cand.vertex()->pos().X();
    float yip_  = v3_.Y()-cand.vertex()->pos().Y();
    float lip_ = v3_.Z()-cand.vertex()->pos().Z();
    float tip_ = sqrt( xip_*xip_ + yip_*yip_ );
    if( tip_<d0_trk && fabs(lip_)<dz0_trk ) samevtx = true;

    // save the information

    int   nVhits    = infot_->getInt("numberOfValidHits");
    int   nIhits    = infot_->getInt("expInnHits");
    float chi2    = infot_->getFloat("Chi2");
    int   ndof      =  infot_->getInt("Ndof");
    int   charge    = (int) trk->charge();
    bool  isfp     = infot_->getBool("isFPHitTrk");
    float pttkcut = 2.;
    int   tkId      = ( isfp==true && pt_>pttkcut )? getTrkGenPdgId(trk,ptk) : 0;
    //float tDepE   = ( isfp==true && pt_>pttkcut )? getTkDeposit(hiteta,hitphi,1) : -1.;
    //float tDepEt  = ( isfp==true && pt_>pttkcut )? getTkDeposit(hiteta,hitphi,2) : -1.;
    float tDepE   = ( isfp==true && pt_>pttkcut )? getTkDeposit(cmsEta,cmsPhi,1) : -1.;
    float tDepEt  = ( isfp==true && pt_>pttkcut )? getTkDeposit(cmsEta,cmsPhi,2) : -1.;
    bool  cnt_    = (isfp && pt_>pttkcut && samevtx) ? countsInTkIso( trk, cand ) : false;

    n_tkecal++;
    vtkpt.push_back(pt_);
    vtkCnt.push_back(cnt_);
    vtkd0.push_back(r0_);
    vtkz0.push_back(z0_);
    vtkNhits.push_back(nVhits);
    vtkNxIhits.push_back(nIhits);
    vtkChi2.push_back(chi2);
    vtkNdof.push_back(ndof);
    vtkQ.push_back(charge);
    vtkFP.push_back(isfp);
    vtkSV.push_back(samevtx);
    //vtkEtaH.push_back(hiteta);
    //vtkPhiH.push_back(hitphi);
    vtkEtaH.push_back(cmsEta);
    vtkPhiH.push_back(cmsPhi);
    vtkDEta.push_back(deta_);
    vtkDPhi.push_back(dphi_);
    vtkDR.push_back(dR_);
    vtkInVeto.push_back(isInVeto);
    vtkVid.push_back(tkId);
    vtkDepE.push_back(tDepE);
    vtkDepEt.push_back(tDepEt);
  }
  if( debug ) cout << "Track loop done." << endl;
}

void
RandomConeSelector::loadXtalDatabase() 
{
  char sochn[200];

  const char* mdt = getenv("HOME"); // stdlib.h
  string mydir = mdt;
  if( mydir.substr(0,4)=="/afs" ) mydir += "/scratch0/";

  sprintf(sochn,"%s/AnaNaS/MECore/XtalDataBase.txt",mydir.data());
  ifstream database(sochn);

  //ifstream database("/home/usr201/mnt/lmillisc/AnaNaS/MECore/XtalDataBase.txt");

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
    abort();
  }

  _dbLoaded = true;
  
}

int
RandomConeSelector::getTrkGenPdgId(const Candidate* candMatch, bool p)
{
  float dRopt  = 9999.;
  int   pdgIdOut     = 0;
  float dPtMax       = .4;
  float dRmax        = .4;
  bool  isMatched    = false;
  EventManager& _e = *(EventManager::e());
  const CandList& mcCandList = _e.mcTruthCandidates();
  Candidate* mcTrack;
  bool usedMC = false;
  int imckept = -1;

  for(int imc=0;imc<(int)mcCandList.size();imc++) {
    Candidate* mcCand = mcCandList[imc];
    usedMC = false;

    int idrun = mcCand->pdgCode();
    if( abs(idrun)<7 || idrun==22 || idrun==130 || idrun==21 || idrun==310 ) continue;
    
    for( int u=0; u<(int)usedTrkMatch.size(); u++) { if( imc==usedTrkMatch[u] ) usedMC = true; }
    if( usedMC==true ) continue;

    float ptrun = mcCand->pt();
    if( fabs(ptrun-candMatch->pt())/candMatch->pt() > dPtMax ) continue;

    float dRrun = KineUtils::dR(mcCand->eta(), candMatch->eta(), mcCand->phi(), candMatch->phi() );
    if( dRrun>dRmax )  continue;
    if( dRrun>dRopt ) continue;
    dRopt = dRrun;

    pdgIdOut = mcCand->pdgCode();
    mcTrack = mcCand;
    isMatched = true;
    imckept = imc;
  }

  usedTrkMatch.push_back(imckept);
  char sochn[100];
  if( p==true ) {
    cout << endl;
    sprintf(sochn,"Track   pt %5.2f eta % 1.2f phi % 1.2f",candMatch->pt(),candMatch->eta(),candMatch->phi());
    cout << sochn << endl;
    if( isMatched ) {
      sprintf(sochn,"Matched pt %5.2f eta % 1.2f phi % 1.2f  pdgId % 3i", mcTrack->pt(), mcTrack->eta(), mcTrack->phi(), mcTrack->pdgCode() );
      cout << sochn << endl;
    }
  }
  return pdgIdOut;
}

void
RandomConeSelector::genList(float hiteta, float hitphi)
{
  char sochn[100];
  
  float dEtaCut = 2.;
  float dPhiCut = 2.;
  float xSize = (fabs(hiteta)<1.45)? 0.0174:0.00864;
  float pi = TMath::Pi();
  
  float r_ecal = (fabs(hiteta)<e_EE)? 129:150;
  float z_ecal = (fabs(hiteta)<e_EE)? 280:325;
  
  EventManager& _e = *(EventManager::e());
  const CandList& mcCandList = _e.mcTruthCandidates();  

  for(unsigned int imc=0;imc<mcCandList.size();imc++) {
    Candidate* mcCand = mcCandList[imc];

    float etaMChit, phiMChit;
    
    float x0_ = mcCand->vertex()->pos().X();
    float y0_ = mcCand->vertex()->pos().Y();
    float z0_ = mcCand->vertex()->pos().Z();
    
    int qtk_ = (int) mcCand->charge(); //(mcCand->pdgCode()>0)? 1:-1;
    
    if( qtk_!= 0 )
    KineUtils::helixToEnvelop(etaMChit,phiMChit,qtk_*mcCand->pt(),mcCand->eta(),mcCand->phi(), x0_, y0_, z0_, r_ecal, z_ecal, 3.8);
    else 
    KineUtils::straightToEnvelop(etaMChit,phiMChit,mcCand->eta(),mcCand->phi(),x0_,y0_,z0_,r_ecal,z_ecal);
    
    if( phiMChit<0 ) phiMChit += 2*pi;
    float dphitmp1 = phiMChit - hitphi;
    float dphitmp2 = (dphitmp1>pi)? dphitmp1-2*pi : (dphitmp1<-pi)? dphitmp1+2*pi : dphitmp1;
    float dphi_ = TMath::Abs(dphitmp2);
    float deta_ = fabs(etaMChit-hiteta);
    
    if( dphi_>dPhiCut*xSize ) continue;
    if( deta_>dEtaCut*xSize ) continue;
    if( mcCand->pdgCode()==21 || fabs(mcCand->pdgCode())<7 ) continue;
    
    
    sprintf(sochn,"[list]   eta % 4.3f phi % 4.3f pdg % 4i pt %6.2f    ",etaMChit,phiMChit,mcCand->pdgCode(),mcCand->pt());
    cout << sochn << endl;
  }
  return;
}


float
RandomConeSelector::getTkDeposit(float hiteta, float hitphi, int EorEt)
{
  float dep = 0;
  float deptmp = 0;
  float etmp = 0;
  char sochn[100];
  
  for(size_t irh=0;irh<ieta.size();irh++) {

    int   ieta_ = ieta[irh];
    int   iphi_ = iphi[irh];
    int   iz_   = iz[irh];
    float et_   = vet[irh];
    float e_    = ve[irh];

    float etaRH_;
    float phiRH_;
    
    int ilab =  10000*(85+ieta_) + 10*iphi_ + iz_ + 1;
    itMapXtal = MapXtal.find( ilab );
    if( itMapXtal == MapXtal.end() )  { 
      cout<<" Xtal Missing in database :  " << ilab << endl; 
      abort(); 
    }
    else {
      etaRH_ = (*itMapXtal).second.first;
      phiRH_ = (*itMapXtal).second.second;
    }
    
    float xSizeEB = 0.0174;
    float xSizeEE = 0.00864;
    float xSize = (TMath::Abs(hiteta)<1.45)? xSizeEB : xSizeEE;
    
    float wdw = 10.;
    float dR_ = KineUtils::dR( etaRH_, hiteta, phiRH_, hitphi );
    if( dR_   > wdw*xSize ) continue;
    if( et_   < 0.150     ) continue;
    
    etmp = (EorEt==1)? e_ : et_;
    deptmp += etmp;

  }
  
  dep = deptmp;
  sprintf(sochn,"[getDep] eta % 4.3f phi % 4.3f         %6.3f GeV %s dep",hiteta,hitphi,dep,(EorEt==1)?"E":"Et");
  return dep;
}


bool RandomConeSelector::countsInTkIso( const Candidate* trk, const Candidate& cand ) 
{
  bool res = true;

  float pt_  = trk->pt();
  float deta_ = trk->eta() - cand.eta();
  
  float dR_= ROOT::Math::VectorUtil::DeltaR(trk->p3(), cand.p3() ); 
  if( dR_<dR_trk_in  )        { return false; }
  if( dR_>dR_trk_out )        { return false; }
  if( fabs(deta_)<strip_trk ) { return false; }
  
  TVector3 v3_ = trk->pos();
  float x0_ = v3_.X()-cand.vertex()->pos().X();
  float y0_ = v3_.Y()-cand.vertex()->pos().Y();
  float z0_ = v3_.Z()-cand.vertex()->pos().Z();
  float r0_ = sqrt( x0_*x0_ + y0_*y0_ );
  
  if( pt_       < pt0_trk ) { return false; }
  if( r0_       > d0_trk  ) { return false; }
  if( fabs(z0_) > dz0_trk ) { return false; }
  
  return res;
}



