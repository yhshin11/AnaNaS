#include "Analysis/src/WAnalysis.hh"

#include <algorithm>
#include <cassert>

using namespace std;

#include <TH1I.h>
#include <TH1F.h>
#include <TH2I.h>
#include <TH2F.h>

#include "Analysis/utils/Constants.hh"
#include "Analysis/core/CandInfo.hh"
#include "Analysis/core/EventManager.hh"
#include "Analysis/selectors/IsolationSelector.hh"
#include "Analysis/selectors/ConeSelector.hh"
#include "Analysis/utils/StatUtils.hh"

int WAnalysis::imode = EventServer::kWm;

ClassImp( WAnalysis )
  
WAnalysis::WAnalysis( Sample& sample, const std::string & collectionFileName ) :
  SampleAnalysis( "mm"+EventServer::Wmode[imode] , sample, collectionFileName ),
  //  RFile_( ("/home/gpfs/manip/mnt/cms/mmarionn/NoBinDistri/"+ (string)(sample.name()) ).c_str(),"RECREATE"),
  outTree_("SplotVariables","Variables for MET Splots"),
  IsoElec( (imode==EventServer::kWe)?EventManager::kElectron:EventManager::kMuon)
{
  Wmode  = EventServer::Wmode[imode];
 
  cout << "\t-------------------------------------------" << endl; 
  cout << "\t---             W  analysis            ----" << endl; 
  cout << "\t-------------------------------------------" << endl; 

 
  // _nDebug = 10;
 
 //Init compteurs
  Nevt_=0;
  Nacceptance_=0;
  NHLT_ =0;
  NPresel_[0] =0;
  NPresel_[1] =0;
  NPresel_[2] =0;
  NIsol_[0] =0;
  NIsol_[1] =0;
  NIsol_[2] =0;
  NID_[0] =0;
  NID_[1] =0;
  NID_[2] =0;
  NSecondCutE_[0]=0;
  NSecondCutE_[1]=0;
  NSecondCutE_[2]=0;
  NConvRej_[0]=0;
  NConvRej_[1]=0;
  NConvRej_[2]=0;
  NOpposite_=0;

  //Iso compteurs
  NTrkIso_[0] = 0;
  NEcalIso_[0] = 0;
  NHcalIso_[0] = 0;
  NTrkIso_[1] = 0;
  NEcalIso_[1] = 0;
  NHcalIso_[1] = 0;
  NTrkIso_[2] = 0;
  NEcalIso_[2] = 0;
  NHcalIso_[2] = 0;

  SampleAnalysis::defineTemplate( "Iso", 500,0,5);
  SampleAnalysis::defineTemplate( "IsoVar", 501,0,1.002);
  SampleAnalysis::defineTemplate( "mass", 200, 0,200);
  
  SampleAnalysis::defineTemplate( "dlphi", 500, 0 , 1 ); 
  SampleAnalysis::defineTemplate( "dPhi", 361, -180.5, 180.5);
  SampleAnalysis::defineTemplate( "pdgId", 5001, -2500.5, 2500.5 );

  SampleAnalysis::defineTemplate( "charge", 3, -1.5,1.5);
 
  SampleAnalysis::defineTemplate( "drvspt", 100, 0., 1., 100, 0., 50. );
  SampleAnalysis::defineTemplate( "drvsdpt", 100, 0., 1., 201, -2.01, 2.01 );
  SampleAnalysis::defineTemplate( "dptopt", 301, -3.01,  3.01  );
  SampleAnalysis::defineTemplate( "IsoChi2", 201, -2.01, 2.01, 200,0,10);
  SampleAnalysis::defineTemplate( "zDr", 500, 0, 1);
  SampleAnalysis::defineTemplate( "rvtx", 500, -2.5, 2.5 );
  SampleAnalysis::defineTemplate( "drvsdr", 100, 0., 0.5, 200, -0.5, 0.5 );

  SampleAnalysis::defineTemplate( "d0", 201, -100.5,100.5);
  SampleAnalysis::defineTemplate( "r0",201, -100.5,100.5);
  SampleAnalysis::defineTemplate( "vx",201, -100.5,100.5);
  SampleAnalysis::defineTemplate( "vy",201, -100.5,100.5);
  SampleAnalysis::defineTemplate( "vz",201, -100.5,100.5);
  SampleAnalysis::defineTemplate( "nvsdpt",50,0,50, 201, -2.01, 2.01);

  SampleAnalysis::defineTemplate( "dPhivsMag",100,0,100,361,-180.5, 180.5);
  SampleAnalysis::defineTemplate( "noAbsPt", 401,-200.5,200.5);

  SampleAnalysis::defineTemplate( "sigieie",400, -1, 1 );
  
  SampleAnalysis::defineTemplate( "Iso2D",  501,0,1.002, 500,0,5 );

  // string ext(sample.name());
  //string name = "/home/gpfs/manip/mnt/cms/mmarionn/NoBinDistri/"+ ext;
  // file.open(name.c_str(), ios::out | ios::trunc );
  
  initRootTree();

}

WAnalysis::~WAnalysis()
{
}

bool
WAnalysis::analyzeEvent()
{
  if((_ievt+1)%500 == 0)
    cout<<" Event "<<_ievt+1<<endl;
  //Elements
  met_  = _e.met( EventManager::kPfMet);
  
  Nevt_++;
  if(_e.inAcceptance())
    Nacceptance_++;
  
  float pttmp=10;
  const CandList& Jets = _e.jetList( EventManager::kPfJet);
  for(int unsigned ij=0;ij<Jets.size();ij++) {  // Electrons
    Candidate* jet = Jets[ij];
    if((fabs(jet->eta())<1.4442 ||
	(fabs(jet->eta())<2.5 && fabs(jet->eta())>1.56 ) ) )
      {
	if(pttmp<jet->pt())
	  { theJet_=jet; pttmp = jet->pt(); }
      }
  }
 
  //Fill Before analysis
  fill( "pt", "MET", met_->Et() );
  fill( "phi", "MET", met_->phi() ); 
 
  //Analyze Way
  // if(!HLTFiring() ) return false;
  // else NHLT_++;
  if(!EventPreselection()) return false;
  else { NPresel_[EcalP_]++; NPresel_[2]++;}
 
  //MET Commissioning
  FillSeveralMET("Presel");
 
  if(!SecondEcut() ) return false;
  else { NSecondCutE_[EcalP_] ++;; NSecondCutE_[2]++;}


  {
    Candidate* gsfTrack_ = _e.getSecond( "electron-GsfTrack", theLepton_);
    CandInfo* info_ = gsfTrack_->info();
    fill("n10", "hitgsf", (int)info_->getBool("isFPHitTrk") );
    
    if(!info_->getBool("isFPHitTrk") )
      return false;
    else
      { NConvRej_[EcalP_] ++;; NConvRej_[2]++; }
  }

  // if(!LeptonID() ) return false;
  if(!IDByHands() ) return false;
  else { NID_[EcalP_]++; NID_[2]++;}

  if(!LeptonIsolation() ) return false;
  else { NIsol_[EcalP_]++; NIsol_[2]++;}
   FillingHistos();
  
  

  return true;
  
}

void
WAnalysis::bookHistograms()
{
  //
  // histogram booking
  //
}

void
WAnalysis::writeHistograms()
{
  // file.close();
  outTree_.Write();
  //RFile_.Close();
  
  cout<< " Na " << Nevt_<<endl;
  cout<<" Nacc "<<Nacceptance_<<"  "<<(float)Nevt_*100/Nacceptance_<<endl;
  cout<<" NHLT "<<NHLT_<<"  "<<((float)NHLT_*100)/Nacceptance_<<endl;
 
  
  for(int i=0;i<3;i++) {
    if(i==0)  cout<<"\n"<<" EB "<<endl;
    else if(i==1)  cout<<"\n"<<" EE "<<endl;
    else cout<<"\n"<<" Combined "<<endl;


  cout<<" Npresel "<<NPresel_[i]<<"  "<<((float)NPresel_[i]*100)/NHLT_
       <<"  "<<StatUtils::EBinom(((float)NPresel_[i]*100)/NHLT_, NHLT_)<<endl;
  cout<< " NSec "<<NSecondCutE_[i]<<"  "<<((float)NSecondCutE_[i]*100)/NPresel_[i]
      <<"  "<<StatUtils::EBinom(((float)NSecondCutE_[i]*100)/NPresel_[i], NPresel_[i])<<endl;
  cout<<" NConvRej "<<NConvRej_[i]<<"  "<<((float)NConvRej_[i]*100)/NSecondCutE_[i]
      <<"  "<<StatUtils::EBinom(((float)NConvRej_[i]*100)/NSecondCutE_[i], NSecondCutE_[i])<<endl;

  cout<<" NID "<<NID_[i]<<"  "<<((float)NID_[i]*100)/NConvRej_[i]
       <<"  "<<StatUtils::EBinom(((float)NID_[i]*100)/NConvRej_[i], NConvRej_[i])<<endl;
  cout<<" Niso "<<NIsol_[i]<<"  "<<((float)NIsol_[i]*100)/NID_[i]
       <<"  "<<StatUtils::EBinom(((float)NIsol_[i]*100)/NID_[i], NID_[i])<<endl;

  cout<<" EffAcc "<<NIsol_[i]<<"  "<<((float)NIsol_[i]*100)/NPresel_[i]
	<<"  "<<StatUtils::EBinom(((float)NIsol_[i]*100)/NPresel_[i], NPresel_[i])<<endl;

  cout<<" \n Track Isolation : "<<((float)NTrkIso_[i]*100)/NID_[i]<<" +- "<<StatUtils::EBinom(((float)NTrkIso_[i]*100)/NID_[i], (int)NID_[i])<<endl;
  cout<<" Ecal Isolation : "<<((float)NEcalIso_[i]*100)/NTrkIso_[i]<<" +- "<<StatUtils::EBinom(((float)NHcalIso_[i]*100)/(float)NTrkIso_[i], (int)NTrkIso_[i])<<endl;
  cout<<" Hcal Isolation : "<<((float)NHcalIso_[i]*100)/NEcalIso_[i]<<" +- "<<StatUtils::EBinom(((float)NHcalIso_[i]*100)/NEcalIso_[i], (int)NEcalIso_[i])<<endl;
  cout<<" Bilan : "<<NHcalIso_[i]<<" <=> "<<NIsol_[i]<<endl;

  }

}

bool
WAnalysis::HLTFiring() {

  if(Wmode=="Wen")
    return _e.isFired("HLT_Ele15_SW_L1R"); //HLT_Ele15_SW_L1R //HLT_L1SingleEG8
  else
    return _e.isFired("HLT_Mu9"); //HLT_Mu9 //HLT_L2Mu9
}


bool
WAnalysis::EventPreselection() {

  bool accept=false;
  
  float pttmp=0;
  bool Presel=false;
 
  if(Wmode=="Wen") {
    const CandList& electrons_ = _e.electrons();
    for(int unsigned ie=0;ie<electrons_.size();ie++) {  // Electrons
      Candidate* electron = electrons_[ie];
      CandInfo* info = electron->info();
      float eta = info->getFloat("caloEta");
      float et = info->getFloat("caloEnergy")/cosh(eta);
      if(et > 30 && (fabs(eta)<1.4442 ||
	(fabs(eta)<2.5 && fabs(eta)>1.56 ) ) )
	{Presel=true;
	if(pttmp<et)
	  { BestLepton_=ie; pttmp = et;
	  theLepton_ = electron;
	  }
	}
    }      
  }
  else {
    const CandList& muons_ = _e.muons();
    for(int unsigned im=0;im<muons_.size();im++) {  // Muons
      Candidate* muon = muons_[im];
      fill("pt","initMuon", muon->pt() );
      if(muon->Et() > 30 && (fabs(muon->eta())<2.4 ) )
	{Presel=true;
	if(pttmp<muon->Et())
	  { BestLepton_=im; pttmp = muon->Et();
	  theLepton_ = muon;
	  }
	}
    }
  }
 
  if(Presel) {
    //Control plots
    fill( "pt", "Lepton", theLepton_->Et() );
    fill( "phi", "Lepton", theLepton_->phi() );
    
    Candidate* WCand = Candidate::create(theLepton_, met_ );
    fill( "mass", "WTransverseMass_noSel", WCand->mass() );
    fill( "pt", "Wpt_noSel", WCand->pt() );
    fill( "eta", "Weta_noSel", WCand->eta() );
    fill( "phi", "Wphi_noSel", WCand->phi() );

    //Lepton Base Matrix
    LCBM_[0][0]=-theLepton_->cosPhi();
    LCBM_[0][1]=-theLepton_->sinPhi();
    LCBM_[1][0]=+theLepton_->sinPhi();
    LCBM_[1][1]=-theLepton_->cosPhi();
    
    EP_="EB"; EcalP_=0;
    if( (fabs(theLepton_->eta() )> 1.479 && Wmode=="Wen") || 
	(fabs(theLepton_->eta() )> 1.3 && Wmode=="Wmn") )
      { EP_="EE";EcalP_=1;}
  }
  if(Presel) {
    accept=true;
  }

  return accept;
}

bool
WAnalysis::LeptonIsolation() {

 float etaSep;
  if(Wmode=="Wen")
    etaSep=1.479;
  else
    etaSep=1.3;

  bool accept=false; 
  
  IsoElec.ecalUseIsoDeposits = false;
  IsoElec.computeIso( *theLepton_ );
  float TrkIso;
  if(Wmode=="Wen")
    TrkIso = TrackIsolation();
  else
    TrkIso = MuonTrackIsolation();
  //trk 0.05
  float cute[2][3]={ {0.08, 0.05, 0.02},
		     {0.06, 0.04, 0.02} };
  float cutm[2][3]={ {0.1, 0.05, 0.02},
		     {0.1, 0.04, 0.02} };
  float cuts[2][3];
  for(int i=0;i<2;i++)
    for(int j=0;j<3;j++) {
      if(Wmode=="Wen") cuts[i][j]=cute[i][j];
      else cuts[i][j]=cutm[i][j];
    }
  // file << theLepton_->pt()/(theLepton_->pt() + IsoElec.S_ecal + IsoElec.S_hcal) <<endl;
  if(fabs(theLepton_->eta())<etaSep)
    {
      fill("IsoVar", "GlobalIsoEB",theLepton_->pt()/(theLepton_->pt() + IsoElec.S_ecal + IsoElec.S_hcal + TrkIso) );
      fill( "Iso" , "TrkGhmEB", IsoElec.S_trk/theLepton_->pt() );
	  
      fill("Iso2D", "IsoEB",  IsoElec.S_trk/theLepton_->pt() , theLepton_->pt()/(theLepton_->pt() + IsoElec.S_ecal + IsoElec.S_hcal) );

      if(TrkIso/theLepton_->pt()  < cuts[0][0])
	{
	  NTrkIso_[EcalP_]++; NTrkIso_[2]++;
	  IsoCateg_="Iso";
	  if(TrkIso/theLepton_->pt()<0.015)
	    IsoCateg_="SIso";

	  fill( "Iso" , "ECALGhmEB"+IsoCateg_, IsoElec.S_ecal/theLepton_->pt() );
	  if(IsoElec.S_ecal/theLepton_->pt()  < cuts[0][1]) {

	    NEcalIso_[EcalP_]++;  NEcalIso_[2]++;
	    fill( "Iso" , "HCALGhmEB"+IsoCateg_, IsoElec.S_hcal/theLepton_->pt() );
		
	    if(IsoElec.S_hcal/theLepton_->pt()  < cuts[0][2] ) {
		  
	      NHcalIso_[EcalP_]++;  NHcalIso_[2]++;
	      accept =true;
	    }
	  }
	}
    }
  else
    {
      fill("IsoVar", "GlobalIsoEE",theLepton_->pt()/(theLepton_->pt() + IsoElec.S_ecal + IsoElec.S_hcal + TrkIso) );
      fill( "Iso" , "TrkGhmEE", IsoElec.S_trk/theLepton_->pt() );
	
      fill("Iso2D", "IsoEE",  IsoElec.S_trk/theLepton_->pt() , theLepton_->pt()/(theLepton_->pt() + IsoElec.S_ecal + IsoElec.S_hcal) );

      if(TrkIso/theLepton_->pt()  < cuts[1][0]) {
	    
	NTrkIso_[EcalP_]++; NTrkIso_[2]++;
	IsoCateg_="Iso";
	if(TrkIso/theLepton_->pt()<0.015)
	  IsoCateg_="SIso";

	fill( "Iso" , "ECALGhmEE"+IsoCateg_, IsoElec.S_ecal/theLepton_->pt() );

	if( IsoElec.S_ecal/theLepton_->pt()  < cuts[1][1] ){

	  NEcalIso_[EcalP_]++; NEcalIso_[2]++;
	  fill( "Iso" , "HCALGhmEE"+IsoCateg_, IsoElec.S_hcal/theLepton_->pt() );

	  if(  IsoElec.S_hcal/theLepton_->pt()  < cuts[1][2] ) {
		
	    NHcalIso_[EcalP_]++; NHcalIso_[2]++;
	    accept =true;
	  }
	}
      }
    }
 
  // fillTree( theLepton_->pt()/(theLepton_->pt() + IsoElec.S_ecal + IsoElec.S_hcal + TrkIso ) );

  return accept;

}


bool
WAnalysis::LeptonID() {

  CandInfo* info_ = theLepton_->info();

  if(Wmode=="Wen")
    return info_->getBool("eidRobustTight");
  else
    return info_->getBool("muidGlobalMuonPromptTight");  
}



bool
WAnalysis::SecondEcut() {
  
  bool accept=true;
  
  if(Wmode=="Wen") {
    const CandList& electrons_ = _e.electrons();
    for(int unsigned ie=0;ie<electrons_.size();ie++) {  // Electrons
      Candidate* electron = electrons_[ie];
      CandInfo* info = electron->info();
      float eta = info->getFloat("caloEta");
      float et = info->getFloat("caloEnergy")/cosh(eta);
      bool type = info->getBool("eidRobustLoose");
      if(et > 20 && ie!= BestLepton_ && type)
	{accept =false;break;}
    }
  }
  else {
    const CandList& muons_ = _e.muons();
    for(int unsigned im=0;im<muons_.size();im++) {  // Muons
      Candidate* muon = muons_[im];
      if(muon->Et() > 20 && im!= BestLepton_)
	{accept =false;break;}
    }
  }
  
  return accept;
 
}


void
WAnalysis::FillingHistos() {
  
  //Filling transverse mass
  const Candidate* lepton = theLepton_;
  Candidate* WCand = Candidate::create(lepton, met_ );

  fill( "mass", "WTransverseMass"+EP_, WCand->mass() );
  fill( "pt", "Wpt"+EP_, WCand->pt() );
  fill( "eta", "Weta"+EP_, WCand->eta() );
  fill( "phi", "Wphi"+EP_, WCand->phi() );

  fill( "dPhi", "MetElec"+EP_,  KineUtils::dPhi(met_->phi(),  lepton->phi() )*180/3.1415 );
 
  

  fill( "pt", "MET_selected"+EP_, met_->Et() );
  fill( "phi", "MET_selected"+EP_, met_->phi() );
  fill( "pt", "lepton_selected"+EP_, lepton->pt() );
  fill( "phi", "lepton_selected"+EP_, lepton->phi() );
  fill( "charge", "lepton"+EP_, lepton->charge() );

  //Fill underlying event
  FillUnderlyingPhoton(WCand);
  // FillUnderlyingJets(WCand);


  //MET Commissioning
  FillSeveralMET("End");

 

}



void
WAnalysis::FillUnderlyingPhoton(Candidate* WCand ) {

 const CandList& Photons = _e.photons();
  Candidate* thePhoton_;

  float pttmp=0;
  for(int unsigned ip=0;ip<Photons.size();ip++) { 
    Candidate* photon = Photons[ip];
    {
      float dr_ = KineUtils::dR(theLepton_->eta(), photon->eta(),
				 theLepton_->phi(), photon->phi() ); 
      
      if(pttmp<photon->pt() && dr_>0.1)
	{ thePhoton_=photon; pttmp = photon->pt(); }
    }
    fill( "pt", Wmode+"_AllPhotons", photon->pt() );
  }
  
  fill( "n", Wmode+"_NPhotons" , Photons.size() );
  if( Photons.size()!=0 && pttmp!=0) {
    fill( "phi", Wmode+"_LeadPhoton"+EP_+IsoCateg_, thePhoton_->phi( kDeg ) );
    fill( "pt", Wmode+"_LeadPhoton"+EP_+IsoCateg_ , thePhoton_->pt() );
    fill( "dphi", Wmode+"_LeadPhoton"+EP_+IsoCateg_, KineUtils::dPhi(WCand->phi(), thePhoton_->phi() )*180/3.1415 );
  }

}

void
WAnalysis::FillUnderlyingJets(Candidate* WCand ) {

  const CandList& Jets = _e.jetList( EventManager::kPatJet);
  Candidate* theJet_;

  float pttmp=10;
  for(int unsigned ij=0;ij<Jets.size();ij++) { 
    Candidate* jet = Jets[ij];
    {
      float dr_ = KineUtils::dR(theLepton_->eta(), jet->eta(),
				 theLepton_->phi(), jet->phi() );
      
      if(pttmp<jet->pt() && dr_>0.1  )
	{ theJet_=jet; pttmp = jet->pt(); }
    }
    fill( "pt", Wmode+"_AllJets", jet->pt() );
  }
  fill( "n", Wmode+"_NJets" , Jets.size() );
  if( Jets.size()!=0 && pttmp!=0) {
    fill( "phi", Wmode+"_LeadJet"+EP_, theJet_->phi( kDeg ) );
    fill( "pt", Wmode+"_LeadJet"+EP_ , theJet_->pt() );
    fill( "dphi", Wmode+"_LeadJet"+EP_, KineUtils::dPhi(WCand->phi(), theJet_->phi() )*180/3.1415 );
    fill( "dphi", "MET_LeadJet"+EP_, KineUtils::dPhi(met_->phi(), theJet_->phi() )*180/3.1415 );
    fill( "dPhi","LeptonJet"+EP_, fabs(KineUtils::dPhi(theJet_->phi(), theLepton_->phi())*180/Constants::pi ) );
    //file << fabs(KineUtils::dPhi(theJet_->phi(), theLepton_->phi())*180/Constants::pi )<<"   ";

 }
  
}

void
WAnalysis::FillSeveralMET(const string & ALevel) {

  //string type;

  //With PF MET
  const Candidate* pfmet_ = _e.met( EventManager::kPfMet );
  FillMETVariables(pfmet_, "PF"+ALevel);

  //With other MET : Calo Corrected or Raw Calo =======
  const Candidate* calomet_  = _e.met( EventManager::kCaloMet);
  
  FillMETVariables(calomet_, "Calo"+ALevel);
  //========================================

  //With Gen MET ========
  //const Candidate* genmet_ = _e.met( EventManager::kGenMet);
  //FillMETVariables(genmet_, "Gen"+ALevel);
  //========================================

  //With Track Corr Met
  const Candidate* trkcmet_  = _e.met( EventManager::kPatMet);
  FillMETVariables(trkcmet_, "TrkC"+ALevel);
  //========================================

  //========================================

} 

float 
WAnalysis::PFBaseIso(const Candidate *cand, int n, const string & ECAL, float& cIso, float& nIso, float& gIso) {
  
  float ns=0,cs=0,ps=0;

  float nw = 3* 0.35;
  float cw = 3* 0.25;
  float pw = 3* 0.40;

  assert( _e.a().load( "el", n ) );
  size_t nDeps_ = (size_t) _e.a().get_i("el","nHCALdep");
  float nsum=0;
  for( size_t jj=0; jj<nDeps_; jj++ ) // Deposits
    {		  
      float e_ =  _e.a().get_vf("el","valHCALdep",jj);
      float dR_ =  _e.a().get_vf("el","dRHCALdep",jj);

      if(dR_ < 0.4) { ns+=e_;
	nsum += e_ / cand->pt();
      }
    }
  fill("Iso","Neutral"+ECAL, nsum);
  nIso = nsum;

  size_t neDeps_ = (size_t) _e.a().get_i("el","nECALdep");
  float gsum=0;
  for( size_t jj=0; jj<neDeps_; jj++ ) // Deposits
    {		  
      float e_ =  _e.a().get_vf("el","E_ECALdep",jj);
      float dR_ =  _e.a().get_vf("el","dRECALdep",jj);

      if(dR_ < 0.4) { ps+=e_;
	gsum += e_ / cand->pt() ;}
    }
  fill("Iso","Gamma"+ECAL,gsum);
  gIso = gsum;

  size_t ntDeps_ = (size_t) _e.a().get_i("el","ntrkdep");
  float csum=0;
  for( size_t jj=0; jj<ntDeps_; jj++ ) // Deposits
    {		  
      float e_ =  _e.a().get_vf("el","valtrkdep",jj);
      float dR_ =  _e.a().get_vf("el","dRtrkdep",jj);

      if(dR_ < 0.4) { cs+=e_;
	csum +=e_ / cand->pt() ; }
    }
  fill("Iso","Charged"+ECAL,csum);
  cIso = csum;

  //Now the global iso
  float GIso;
  GIso = cand->pt() / (cand->pt() + (nw)*ns + (pw)*ps + (cw)*cs);
  
  fill("IsoVar", "GlobalIso"+ECAL, GIso);

  return GIso;

}

void
WAnalysis::MCMatching(const string & P) {
  
  std::map< Candidate* , pair<float,float> > MCDrAssocMap;
  const CandList& MCList = _e.mcTruthCandidates();

  for(unsigned int imc=0;imc<MCList.size(); imc++) {
    Candidate* MCCand = MCList[imc];
    
    float dR=KineUtils::dR(MCCand->eta(), theLepton_->eta(), 
			   MCCand->phi(), theLepton_->phi() );
    float dPt = fabs( fabs(MCCand->pt()) - fabs(theLepton_->pt()) )/theLepton_->pt();
	
    
    MCDrAssocMap[ MCCand ].first= dR;
    MCDrAssocMap[ MCCand ].second= dPt;

  }

  std::map<Candidate*, pair<float,float> >::const_iterator it;
  std::map<Candidate*, pair<float,float> >::const_iterator itbegin= MCDrAssocMap.begin();
  std::map<Candidate*, pair<float,float> >::const_iterator itend= MCDrAssocMap.end();
  
  float dRtmp=10;
  Candidate* BestCand;
  for(it = itbegin;it != itend; it++) {
    
    float dR = (*it).second.first;
    if(dR < dRtmp) {
      dRtmp = dR;
      BestCand = (*it).first;
    }
  }
  fill( "eta", "MCMatchingLepton"+P, BestCand->eta() );
  fill( "pdgId", "MCMatchingLepton"+P, BestCand->pdgCode() );


  std::map< Candidate* , pair<float,float> > TrkDrAssocMap;
  const CandList& TracksList = _e.tracks();
  for(unsigned int it=0;it<TracksList.size();it++) {
    Candidate* Track = TracksList[it];

    float dR=KineUtils::dR(Track->eta(), theLepton_->eta(), 
			   Track->phi(), theLepton_->phi() );
    float dPt = fabs( fabs(Track->pt()) - fabs(theLepton_->pt()) )/theLepton_->pt();
    
    if(dR<0.2 && dPt< 0.4) {
      TrkDrAssocMap[ Track ].first= dR;
      TrkDrAssocMap[ Track ].second= dPt;
    }
  }

  std::map<Candidate*, pair<float,float> >::const_iterator itt;
  std::map<Candidate*, pair<float,float> >::const_iterator ittbegin= TrkDrAssocMap.begin();
  std::map<Candidate*, pair<float,float> >::const_iterator ittend= TrkDrAssocMap.end();
  
  float dRttmp=10;
  Candidate* BestTrkCand;
  for(itt = ittbegin;it != ittend; itt++) {
    
    float dR = (*it).second.first;
    
    if(dR < dRttmp) {
      dRttmp = dR;
      BestTrkCand = (*it).first;
    }
  }
  
  CandInfo* inftrk = BestTrkCand->info(); 
  fill("n10", "TRK"+P, TrkDrAssocMap.size() );
  fill( "d0", "Trk"+P, inftrk->getFloat("d0") );
  fill( "r0", "Trk"+P, inftrk->getFloat("r0") );
  fill( "vx", "Trk"+P, inftrk->getFloat("vx") );
  fill( "vy", "Trk"+P, inftrk->getFloat("vy") );
  fill( "vz", "Trk"+P, inftrk->getFloat("vz") );
  
  

}

int WAnalysis::FindBestdRCandidate(std::map<int, std::pair<float,float> > AssocMap ) {

  std::map<int, pair<float,float> >::const_iterator it;
  std::map<int, pair<float,float> >::const_iterator itbegin= AssocMap.begin();
  std::map<int, pair<float,float> >::const_iterator itend= AssocMap.end();
  
  float dRtmp=10;

  int BestCand;
  if(AssocMap.size()==0) return -1;
  for(it = itbegin;it != itend; it++) {
    
    float dR = (*it).second.first;

    if(dR < dRtmp) {
      dRtmp = dR;
      BestCand = (*it).first;
     
    }
  }
  
  return BestCand;
}




float WAnalysis::TrackIsolation() {

  string EP = "EB";
  if(fabs(theLepton_->eta())> 1.479 )
    EP="EE";

  int indexTrk = -1;
  Candidate* elTrk_ = _e.getSecond("electron-Track", theLepton_ );
  if(elTrk_!=NULL) {
    CandInfo* inftrk = elTrk_->info();
    fill("n10", "hittrk", (int)inftrk->getBool("isFPHitTrk") );
    indexTrk = inftrk->getInt("index");
  }
  //=============
  float S_trk=0;
  int n=0;
 

  const CandList& listTracks = _e.tracks();
  for( size_t ii=0; ii<listTracks.size(); ii++ )
    {
      Candidate* trk = listTracks[ii];
      CandInfo* info = trk->info();
      //      bool FPH_ = info->getBool("isFPHitTrk");
      TVector3 vtx = trk->pos();
      float rVtx = vtx.Perp();
      float pt_  = trk->pt();
      float eta_ = trk->eta();
      float phi_ = trk->phi();
      float dR_=KineUtils::dR( theLepton_->eta(), eta_, theLepton_->phi(), phi_ );

      if((int unsigned)indexTrk!=ii && pt_>1. && dR_<0.4 &&  rVtx < 0.45 ) {
	
	fill("dr", "DR"+EP, dR_);
	fill("drvsdpt", "dRVsdPTIso"+EP, dR_, (theLepton_->pt()-pt_)/theLepton_->pt()   );
	
	
	bool inCone_=true;

	if( inCone_ ) 
	  { n++; S_trk+=pt_;
	  fill("IsoChi2", "chi2trk"+EP,
	       (theLepton_->pt()-pt_)/theLepton_->pt(),
	       info->getFloat("Chi2")/info->getInt("Ndof") );
	  fill("drvsdpt", "detaTrk"+EP, fabs(theLepton_->eta()-eta_),(theLepton_->pt()-pt_)/theLepton_->pt() );
	  fill("drvsdpt","dphitrk"+EP, KineUtils::dPhi(theLepton_->phi(), phi_),(theLepton_->pt()-pt_)/theLepton_->pt() );
	  fill("nvsdpt", "NValid"+EP, info->getInt("numberOfValidHits"),(theLepton_->pt()-pt_)/theLepton_->pt() );
	  fill("nvsdpt", "NLost"+EP, info->getInt("numberOfLostHits"),(theLepton_->pt()-pt_)/theLepton_->pt());
	  fill("nvsdpt", "Ndof"+EP, info->getInt("Ndof"),(theLepton_->pt()-pt_)/theLepton_->pt() );
	  fill("drvsdr", "detavsdphi"+EP, fabs(theLepton_->eta()-eta_), KineUtils::dPhi(theLepton_->phi(), phi_) );
	  }   
      }
    }
  fill("Iso", "TestTrk"+EP , S_trk/theLepton_->pt() );
  fill("pt", "TestTrk"+EP, S_trk );
  fill( "dPhi", "MetElecAtIso"+EP_,  KineUtils::dPhi(met_->phi(),  theLepton_->phi() )*180/3.1415 );
  
  Candidate* WCand = Candidate::create(theLepton_, met_ );
  FillUnderlyingJets(WCand);

    
  return S_trk;
  
}


float WAnalysis::MuonTrackIsolation() {

  std::map<int, pair<float,float> > TrackMap;
  
  string EP = "EB";
  if(fabs(theLepton_->eta())> 1.5 )
    EP="EE";

  const CandList& listTracks = _e.tracks();

  for( size_t ii=0; ii<listTracks.size(); ii++ )
    {
      Candidate* trk = listTracks[ii];
      float pt_  = trk->pt();
      float eta_ = trk->eta();
      float phi_ = trk->phi();
      TVector3 vtx = trk->pos();
      float rVtx = vtx.Perp();
    
      float dr=KineUtils::dR(theLepton_->eta(), eta_, theLepton_->phi(), phi_);
      float dpt= fabs(theLepton_->pt()-pt_)/theLepton_->pt() ;
      if( dr<0.3 && rVtx<0.5 && dpt<0.5  ) {
	TrackMap[ ii ].first = dr;
	TrackMap[ ii ].second = dpt;
      }
    }
 
  int besttrk = FindBestdRCandidate(TrackMap);
  Candidate* BestTrkCand;
  if(besttrk != -1) {
    BestTrkCand = listTracks[ besttrk ];
    
    TVector3 vtx = BestTrkCand->pos();
    float rVtx = vtx.Perp();
    fill( "rvtx", "EVtx", rVtx);
    fill( "dptopt", "dPtTrack"+EP, TrackMap[ besttrk ].second );
    fill( "dr", "DRTrack"+EP , TrackMap[besttrk ].first );
    fill( "drvsdpt", "MatchTrack"+EP, TrackMap[ besttrk ].first, TrackMap[ besttrk ].second  );
  }
  
  float S_trk=0;
  int n=0;
 
  for( size_t ii=0; ii<listTracks.size(); ii++ )
    {
      Candidate* trk = listTracks[ii];
      CandInfo* info = trk->info();
      TVector3 vtx = trk->pos();
      float rVtx = vtx.Perp();
      float pt_  = trk->pt();
      float eta_ = trk->eta();
      float phi_ = trk->phi();
      float dR_=KineUtils::dR( theLepton_->eta(), eta_, theLepton_->phi(), phi_ );
      
      if((int unsigned)besttrk!=ii && pt_>1. && dR_<0.4) fill("rvtx", "Vtx", rVtx);
      if((int unsigned)besttrk!=ii && pt_>1. && dR_<0.4 && rVtx < 0.45) {

	
	fill("dr", "DR"+EP, dR_);
	fill("drvsdpt", "dRVsdPTIso"+EP, dR_, (theLepton_->pt()-pt_)/theLepton_->pt()   );
	fill("drvsdpt", "dRVsdPTIsoTrk"+EP, dR_, (BestTrkCand->pt()-pt_)/BestTrkCand->pt()   );
	
	bool inCone_=true;
	
	if( inCone_ ) 
	  { n++; S_trk+=pt_;
	  fill("IsoChi2", "chi2trk"+EP,
	       (theLepton_->pt()-pt_)/theLepton_->pt(),
	       info->getFloat("Chi2")/info->getInt("Ndof") );
	  fill("drvsdpt", "detaTrk"+EP, fabs(theLepton_->eta()-eta_),(theLepton_->pt()-pt_)/theLepton_->pt() );
	  fill("drvsdpt","dphitrk"+EP, KineUtils::dPhi(theLepton_->phi(), phi_),(theLepton_->pt()-pt_)/theLepton_->pt() );
	  fill("nvsdpt", "NValid"+EP, info->getInt("numberOfValidHits"),(theLepton_->pt()-pt_)/theLepton_->pt() );
	  fill("nvsdpt", "NLost"+EP, info->getInt("numberOfLostHits"),(theLepton_->pt()-pt_)/theLepton_->pt());
	  fill("nvsdpt", "Ndof"+EP, info->getInt("Ndof"),(theLepton_->pt()-pt_)/theLepton_->pt() );
	  fill("drvsdr", "detavsdphi"+EP, fabs(theLepton_->eta()-eta_), KineUtils::dPhi(theLepton_->phi(), phi_) );
	  }   
      }
    }
  fill("Iso", "TestTrk"+EP , S_trk/theLepton_->pt() );

  return S_trk;

}


TVector2
WAnalysis::ComputeRecoil(const Candidate* met) {
  
  TVector2 Recoil(0,0);
  Recoil = met->p2() + theLepton_->p2();
  return Recoil;

}

TVector2
WAnalysis::ComputeRecoilWithPF(const Candidate* met) {
  
  TVector2 Recoil(0,0);

  Candidate* pfcand_ = _e.getSecond("electron-PfCand", theLepton_ );
 if(pfcand_!=NULL)
   Recoil = met->p2() + pfcand_->p2();
 return Recoil;

}

TVector2
WAnalysis::ComputeRecoilWithCaloRecHit(const Candidate* met) {
  
  TVector2 Recoil(0,0);

  TVector2 sumRH_(0,0);

  //const CandAssoc& e_Sc_ =  _e.candAssoc("electron-SuperCluster");
  const CandMap& Sc_Bc_ = _e.candMap("SuperCluster-BasicCluster");
  const CandIdRecHitMap& Bc_rh_ = _e.ecalRecHitMap("BasicCluster-EcalRecHit");

  Candidate* superCluster = _e.getSecond("electron-SuperCluster", theLepton_ );
  CandMapIterator itbc_lower_ = Sc_Bc_.lower_bound( superCluster );
  CandMapIterator itbc_upper_ = Sc_Bc_.upper_bound( superCluster );
  CandMapIterator itbc;
  for( itbc=itbc_lower_; itbc!=itbc_upper_; itbc++ )
    { 
      CandIdRecHitMapIterator itrh_lower_ = Bc_rh_.lower_bound( (*itbc).second->uid() );
      CandIdRecHitMapIterator itrh_upper_ = Bc_rh_.upper_bound( (*itbc).second->uid() );
      CandIdRecHitMapIterator itrh;
      for( itrh=itrh_lower_; itrh!=itrh_upper_; itrh++ )
	{ 
	  sumRH_ += ComputeEtVectorFromRecHit( (*itrh).second.first );
	}
    }
 
  Recoil = met->p2() + sumRH_;
  return Recoil;

}

TVector2
WAnalysis::ComputeRecoilWithCalo(const Candidate* met) {
  
  TVector2 Recoil(0,0);
  Candidate* superCluster_ = _e.getSecond("electron-SuperCluster", theLepton_ );
  if(superCluster_!=NULL)
    //  superCluster_ = _e.getSecond("electron-PfSuperCluster", theLepton_ );

    Recoil = met->p2() + superCluster_->p2();
  return Recoil;

}


TVector2
WAnalysis::ComputeEtVectorFromRecHit(ecalRecHit rh) {


  TVector2 etVector(0,0);
  float eta_ =0;
  float phi_ =0;
  float et_ = 0;
  float x_ = 0;
  float y_ = 0;
   
  if(rh.iz==0) //Barrel
    {
      MEEBGeom::EtaPhiPoint point_ = MEEBDisplay::center( rh.ix, rh.iy );
      eta_ = point_.first;
      phi_ = point_.second * Constants::pi;

      et_ = KineUtils::et(rh.e , eta_);
      x_ = cos(phi_)*et_;
      y_ = sin(phi_)*et_;

      etVector.Set( x_, y_);

    }
  else //Endcap
    {
      MEEEGeom::EtaPhiPoint point_ = MEEEDisplay::center( rh.ix, rh.iy, rh.iz );
      eta_ = point_.first;
      phi_ = point_.second * Constants::pi;
      
      et_ = KineUtils::et(rh.e , eta_);
      x_ = cos(phi_)*et_;
      y_ = sin(phi_)*et_;
      
      etVector.Set( x_, y_);
    }
 
  return etVector;

}

void
WAnalysis::FillRecoil(TVector2 Recoil, Candidate* W, const string & type){

  fill( "pt",type+"Recoil",Recoil.Mod() );
  //Recoil angle btw 0,2 pi
  Double_t rphi = TVector2::Phi_0_2pi(Recoil.Phi() );
  fill( "dPhi",type+"Recoil", KineUtils::dPhi(rphi,theLepton_->phi() )*180/3.141592 );
  fill( "dPhivsMag",type+"Recoil", Recoil.Mod() , KineUtils::dPhi(rphi,theLepton_->phi() )*180/3.141592);

  //Recoil Projection
  float Rpara= LCBM_[0][0]*Recoil.Px()+LCBM_[0][1]*Recoil.Py();
  float Rperp= LCBM_[1][0]*Recoil.Px()+LCBM_[1][1]*Recoil.Py();
  
  fill("noAbsPt",type+"Recoil_ParaLepton", Rpara);
  fill("noAbsPt",type+"Recoil_PerpLepton", Rperp);

  fill("noAbsPt",type+"Recoil_X", Recoil.Px() );
  fill("noAbsPt",type+"Recoil_Y", Recoil.Py() );

  //Recoil projection on W axis
  float RWpara = W->cosPhi()*Recoil.Px() + W->sinPhi()*Recoil.Py();
  float RWperp = (- W->sinPhi())*Recoil.Px() + W->cosPhi()*Recoil.Py();

  fill("noAbsPt",type+"Recoil_ParaW", RWpara);
  fill("noAbsPt",type+"Recoil_PerpW", RWperp);

 
  //Recoil projection on cluster direction
  if(type.substr(0,1)=="C") {
    Candidate* superCluster_ = _e.getSecond("electron-SuperCluster", theLepton_ );
    if(superCluster_!=NULL) {
      float RSCpara = superCluster_->cosPhi()*Recoil.Px() + superCluster_->sinPhi()*Recoil.Py();
      float RSCperp = (- superCluster_->sinPhi())*Recoil.Px() + superCluster_->cosPhi()*Recoil.Py();
      //cout<<type<<"  "<< Rpara<<"  "<<RSCpara<<"   "<<theLepton_->phi()<<"  "<<superCluster_->phi()<<endl;
      fill("noAbsPt",type+"Recoil_ParaCalo", RSCpara);
      fill("noAbsPt",type+"Recoil_PerpCalo", RSCperp);
    }
  }
}


void
WAnalysis::FillMETProjection(const Candidate* met, Candidate* W, const string & type) {

  //MET Projection on lepton axis
  float METpara= LCBM_[0][0]*met->px()+LCBM_[0][1]*met->py();
  float METperp= LCBM_[1][0]*met->px()+LCBM_[1][1]*met->py();
  
  fill("noAbsPt",type+"MET_ParaLepton", METpara);
  fill("noAbsPt",type+"MET_PerpLepton", METperp);


  //MET Projections on in frame
  fill("noAbsPt",type+"MET_X", met->px() );
  fill("noAbsPt",type+"MET_Y", met->py() );

  //MET projection on W axis  
  float METWpara = W->cosPhi()*met->px() + W->sinPhi()*met->py();
  float METWperp = (- W->sinPhi())*met->px() + W->cosPhi()*met->py();
  //float METWpara = met->p2()*W->p2()/W->pt();
  //float METWperp = (met->p2() + METWpara*(W->p2().Unit()) ).Mod();
  
  fill("noAbsPt",type+"MET_ParaW", METWpara);
  fill("noAbsPt",type+"MET_PerpW", METWperp);

}



void
WAnalysis::FillMETVariables(const Candidate* met, const string & type) {

 Candidate* W = Candidate::create(theLepton_, met );

  fill( "mass", "W"+type+"_"+EP_, W->mass() );
  fill( "pt", "W"+type+"_"+EP_, W->pt() );
  fill( "pt", type+"MET_"+EP_, met->pt() );
  fill( "eta", "W"+type+"_"+EP_, W->eta() );
  fill( "phi", "W"+type+"_"+EP_, W->phi() );
  fill( "dphi", type+"Dphi_"+EP_, KineUtils::dPhi(met->phi(),theLepton_->phi() )*180/3.141592 );
  
  FillMETProjection(met, W, type); 

  
  TVector2 Recoil = ComputeRecoil(met);
  FillRecoil(Recoil, W, type);

  
  /* TVector2 CaloRecoil = ComputeRecoilWithCalo(met);
  TVector2 PFRecoil = ComputeRecoilWithPF(met);
  TVector2 CalRecHitRecoil = ComputeRecoilWithCaloRecHit(met);
  // if( type.substr(0,5)=="CaloE")
  // cout<<type<<"  "<<CaloRecoil.Px()<<"  "<<Recoil.Px()<<"  "<<CalRecHitRecoil.Px()<<endl;
 
  if( type.substr(0,5)=="CaloE" && CaloRecoil.Px()!=0 && CaloRecoil.Py()!=0 ) {
      FillRecoil(Recoil, W, type);
      FillRecoil(CaloRecoil, W, type+"Calo");
      FillRecoil(CalRecHitRecoil, W, type+"RecHit");
    }
  if( type.substr(0,3)=="PFE" && PFRecoil.Px()!=0 && PFRecoil.Py()!=0 ) {
      FillRecoil(Recoil, W, type);
      FillRecoil(PFRecoil, W, type+"PF");
    }*/
    
}



bool WAnalysis::IDByHands() {

  CandInfo* info_ = theLepton_->info();

  float hoe = info_->getFloat("HOE");
  float dphi = info_->getFloat("dPhiIn");
  float deta = info_->getFloat("dEtaIn");
  float sigieie = info_->getFloat("sigiEtaiEta");

  //HOE, dphi, deta, sigieie
   float dEtaIn[2];
    float dPhiIn[2];
    float sigiEiE[2];
    float HOE[2];

    //Working Point 70 %
    dEtaIn[0]= 0.006; dEtaIn[1]= 0.003; 
    dPhiIn[0]= 0.02; dPhiIn[1]= 0.02; 
    sigiEiE[0]= 0.01; sigiEiE[1]= 0.03;
    HOE[0]= 0.02; HOE[1]= 0.0025;
     //****************

    /* float cuts[2][4]={{0.01,0.025,0.0040,0.0099},
       {0.01,0.020,0.0066,0.028}};*/

  bool acc =false;
  //Inversion deta
  if(EP_=="EB") {
    if(hoe < HOE[0] && fabs(sigieie) < sigiEiE[0] && fabs(dphi) < dPhiIn[0] && fabs(deta) < dEtaIn[0] ) 
      acc =true;
  }
  else {
    if(hoe < HOE[1] && fabs(sigieie) < sigiEiE[1] && fabs(dphi) < dPhiIn[1] && fabs(deta) < dEtaIn[1] )
      acc = true;
  }
 
  return acc;

}



/*
void WAnalysis::fillIDObs(int IEp) {

 CandInfo* info_ = theLepton_->info();

  float hoe = info_->getFloat("HOE");
  float dphi = info_->getFloat("dPhiIn");
  float deta = info_->getFloat("dEtaIn");
  float sigieie = info_->getFloat("sigiEtaiEta");

  string ET;
  if(IEp==0)
    ET="EB";
  else
    ET="EE";

  fill("HOE","elec"+ET,hoe);
  fill("dPhiIn","elec"+ET,dphi);
  fill("dEtaIn","elec"+ET,deta);
  fill("sigieie","elec"+ET,sigieie);

}

*/


void WAnalysis::initRootTree() {
  outTree_.Branch("dPhi",&dPhiJetLep,"dPhiJetLep/F");
  outTree_.Branch("IsoVar",&IsoVariable,"IsoVariable/F");
  outTree_.Branch("MET",MET,"MET[2]/F");
  outTree_.Branch("leptonPt",leptonPt,"leptonPt[3]/F");
}


void WAnalysis::fillTree(float iso) {

  //Compute the dPhi
  const CandList& Jets = _e.jetList( EventManager::kPatJet);
  Candidate* theJet_;
  float pttmp=10;
  for(int unsigned ij=0;ij<Jets.size();ij++) { 
    Candidate* jet = Jets[ij];
    {
      float dr_ = KineUtils::dR(theLepton_->eta(), jet->eta(),
				theLepton_->phi(), jet->phi() );
      if(pttmp<jet->pt() && dr_>0.1  )
	{ theJet_=jet; pttmp = jet->pt(); }
    }
  }
  if( Jets.size()!=0 && pttmp>10) {
    dPhiJetLep = (float)fabs(KineUtils::dPhi(theJet_->phi(), theLepton_->phi())*180/Constants::pi );
  }
  else
    {
      dPhiJetLep = -1;
    }
  
 
  IsoVariable = (float)iso;
  MET[0] = met_->px(); MET[1] = met_->py();
  leptonPt[0] = theLepton_->px(); 
  leptonPt[1] = theLepton_->py();
  leptonPt[2] = theLepton_->pz();

  outTree_.Fill();
}
