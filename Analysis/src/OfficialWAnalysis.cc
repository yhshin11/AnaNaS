#include "Analysis/src/OfficialWAnalysis.hh"

#include <algorithm>
#include <cassert>
#include <sstream>

using namespace std;

#include <TH1I.h>
#include <TH1F.h>
#include <TH2I.h>
#include <TH2F.h>

#include "Analysis/utils/Constants.hh"
#include "Analysis/utils/Config.hh"
#include "Analysis/core/CandInfo.hh"
#include "Analysis/core/EventManager.hh"
#include "Analysis/selectors/IsolationSelector.hh"
#include "Analysis/selectors/ConeSelector.hh"
#include "Analysis/utils/StatUtils.hh"

int OfficialWAnalysis::imode = EventServer::kWm;

ClassImp( OfficialWAnalysis )
  
OfficialWAnalysis::OfficialWAnalysis( Sample& sample, const std::string & collectionFileName ) :
  SampleAnalysis( EventServer::Wmode[imode] , sample, collectionFileName ),
  IsoElec( (imode==EventServer::kWe)?EventManager::kElectron:EventManager::kMuon),
  vbtfSel( VBTFElectronSelector::kEt20, VBTFElectronSelector::kPresel),
  outTree_("METCommVariables","Variables for MET Splots")
{
  Wmode  = EventServer::Wmode[imode];
 
  cout << "\t-------------------------------------------" << endl; 
  cout << "\t---         W  analysis "<<Wmode<<"            ----" << endl; 
  cout << "\t-------------------------------------------" << endl; 

  // InitASCIIS();

  _nDebug = 10;
 
  vbtfSel.useOfficialIsoComputer = true;

  //
  SampleName = sample.name();

 //Init compteurs
  Nevt_=0;
  Nacceptance_=0;
  NHLT_ =0;
  NPresel_[0] =0;
  NPresel_[1] =0;
  NPresel_[2] =0;
  NID_[0] =0;
  NID_[1] =0;
  NID_[2] =0;
  NIso_[0] =0;
  NIso_[1] =0;
  NIso_[2] =0;
  NSecondCutE_[0]=0;
  NSecondCutE_[1]=0;
  NSecondCutE_[2]=0;
  NConvRej_[0]=0;
  NConvRej_[1]=0;
  NConvRej_[2]=0;

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

  SampleAnalysis::defineTemplate( "METvsLepEta", 400, 0., 200. , 240, -3., 3. );
  SampleAnalysis::defineTemplate( "MassvsLepEta", 400, 0., 200. , 240, -3., 3. );
 
  //  SampleAnalysis::defineTemplate( "d0", 200, -50,50);


}

OfficialWAnalysis::~OfficialWAnalysis()
{
}

bool
OfficialWAnalysis::analyzeEvent()
{

  if((_ievt)%1000 == 0)
    cout<<" Event "<<_ievt<<endl;
  //Elements
  met_  = _e.met( EventManager::kPfMet);
  
  Nevt_++;
  if(_e.inAcceptance())
    Nacceptance_++;
  
  float pttmp=0;
  const CandList& Jets = _e.jetList( EventManager::kPfJet);
  for(int unsigned ij=0;ij<Jets.size();ij++) { 
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


  //Vertices
  // if(SampleName.substr(0,5)!="ttbar")
  if(!TreeFillVertices() )return false;

  //Analysis Way
  if(!HLTFiring() ) return false; //HLT Filter 1st step
  else NHLT_++;
  
  if(!EventPreselection()) return false; //Preselection 2nd step
  else { NPresel_[EcalP_]++; NPresel_[2]++;}

  //Ensure no overlapping and spike rejection in electrons  
  if(SampleName.substr(0,8)=="FullLumi" || SampleName.substr(0,2)=="EG" || SampleName.substr(0,5)=="SD_EG" ) {
    
    if(SpikeCleaning()) //Presence of spike 
      return false;

    int run,event;
    for(int unsigned i=0;i<Events.size();i++) {
      run = Events[i].first;
      event = Events[i].second;

      if(run == _e.run() && event==_e.event())
	return false;
    }    
    std::pair<int,int> tmp(_e.run(),_e.event());
    Events.push_back(tmp);
  }
     

  //MET Commissioning
  FillingHistos("Presel");

  if(!SecondEcut() ) return false; //Veto Second electron 3rd step
  else { NSecondCutE_[EcalP_] ++; NSecondCutE_[2] ++; }
 
   if(!ConversionRejection() ) return false; //Conversion Rejection 4th step
  else { NConvRej_[EcalP_] ++;  NConvRej_[2] ++;}
   // ConversionRejection();
  //HLT object Matching missing


  //*** End W Candidate creation ***
  //*** Begin W ploting for several IDs/Working points ***

    if(!LeptonID() ) return false;
  else { NID_[EcalP_]++; NID_[2]++; }
  
  if(!LeptonIso() ) return false;
  else { NIso_[EcalP_]++; NIso_[2]++; }
  

  //if(!LeptonID() ) return false;
  //else{/*fillTree();*/}
  

  PrepareTree(0,0);

  //MM go to Tree
  // VBTFSelection();


  FillingHistos("Et25_Sel80");
  
  

  return true;
  
}

void
OfficialWAnalysis::bookHistograms()
{
  //
  // histogram booking
  //
  //  outTree_.SetDirectory(0);
  outTree_.SetDirectory( _outputFile->GetDirectory(NULL) );
  initRootTree();
}

void
OfficialWAnalysis::writeHistograms()
{

  outTree_.Write();

  cout<< " Na " << Nevt_<<endl;
  cout<<" Nacc "<<Nacceptance_<<"  "<<(float)Nevt_*100/Nacceptance_<<endl;
  cout<<" NHLT "<<NHLT_<<"  "<<((float)NHLT_*100)/Nacceptance_<<endl;
 
  
  for(int i=0;i<3;i++) {
    if(i==0)  cout<<"\n"<<" EB "<<endl;
    else if(i==1)  cout<<"\n"<<" EE "<<endl;
    else  cout<<"\n"<<" Combined "<<endl;

  cout<<" Npresel "<<NPresel_[i]<<"  "<<((float)NPresel_[i]*100)/NHLT_
       <<"  "<<StatUtils::EBinom(((float)NPresel_[i]*100)/NHLT_, NHLT_)<<endl;
  cout<< " NSec "<<NSecondCutE_[i]<<"  "<<((float)NSecondCutE_[i]*100)/NPresel_[i]
      <<"  "<<StatUtils::EBinom(((float)NSecondCutE_[i]*100)/NPresel_[i], NPresel_[i])<<endl;
  cout<<" NConvRej "<<NConvRej_[i]<<"  "<<((float)NConvRej_[i]*100)/NSecondCutE_[i]
      <<"  "<<StatUtils::EBinom(((float)NConvRej_[i]*100)/NSecondCutE_[i], NSecondCutE_[i])<<endl;
  cout<<" NID "<<NID_[i]<<"  "<<((float)NID_[i]*100)/NConvRej_[i]
      <<"  "<<StatUtils::EBinom(((float)NID_[i]*100)/NConvRej_[i], NConvRej_[i])<<endl;
  cout<<" NIso "<<NIso_[i]<<"  "<<((float)NIso_[i]*100)/NID_[i]
      <<"  "<<StatUtils::EBinom(((float)NIso_[i]*100)/NID_[i], NID_[i])<<endl;

  cout<<" \n EffAcc "<<NIso_[i]<<"  "<<((float)NIso_[i]*100)/NPresel_[i]
      <<"  "<<StatUtils::EBinom(((float)NIso_[i]*100)/NPresel_[i], NPresel_[i])<<endl;
  
  }

}

bool
OfficialWAnalysis::HLTFiring() {
  bool HLTfiring=false;

  int run = _e.run();

  if(Wmode=="Wen") {
    if(SampleName.substr(0,1)=="E") {
 
    if( ( (run <= 140401) && _e.isFired("HLT_Ele15_LW_L1R") ) ||
	( (run >= 140402 && run <= 143962 ) && _e.isFired("HLT_Ele15_SW_L1R") ) ||
	( (run >= 143963 && run <= 144114 ) && _e.isFired("HLT_Ele15_SW_CaloEleId_L1R") ) ||
	( (run >= 144115 && run <= 147116 ) && _e.isFired("HLT_Ele17_SW_CaloEleId_L1R") ) ||
	( (run >= 147117 && run <= 148058 ) && _e.isFired("HLT_Ele17_SW_TightEleId_L1R") ) ||
	( (run >= 148059 && run <= 149064 ) && _e.isFired("HLT_Ele17_SW_TighterEleIdIsol_L1R_v2") ) ||
	( (run >= 149065 && run <= 149442 ) && _e.isFired("HLT_Ele17_SW_TighterEleIdIsol_L1R_v3") ) )
      HLTfiring = true;

    }
     else
       // if(_e.isFired("HLT_Ele15_LW_L1R") || _e.isFired("HLT_Ele20_LW_L1R") )
       HLTfiring=true;
    
    //FIXME
  // HLTfiring = true;
    return HLTfiring;//_e.isFired("HLT_Photon15_Cleaned_L1R"); //HLT_Ele15_SW_L1R //HLT_L1SingleEG8

  }
  else
    return true; //FIXME
  //    return _e.isFired("HLT_Mu9"); //HLT_Mu9 //HLT_L2Mu9
}


bool
OfficialWAnalysis::EventPreselection() {

  bool accept=false;
  
  float pttmp=0;
  bool Presel=false;
  float eta=0;
  if(Wmode=="Wen") {
    const CandList& electrons_ = _e.electrons();
    for(int unsigned ie=0;ie<electrons_.size();ie++) {  // Electrons
      Candidate* electron_ = electrons_[ie];
        CandInfo* info = electron_->info();
	eta = info->getFloat("caloEta");
	float et = info->getFloat("caloEnergy")/cosh(eta);
      if(et > 25 && (fabs(eta)<=1.44 ||
	(fabs(eta)<=2.5 && fabs(eta)>=1.56 ) ) )
	{Presel=true;
	if(pttmp<et)
	  { BestLepton_=ie; pttmp = et;
	  theLepton_ = electron_;
	  }
	}
    }      
  }
  else {
    const CandList& muons_ = _e.muons();
    for(int unsigned im=0;im<muons_.size();im++) {  // Muons
      Candidate* muon = muons_[im];
      fill("pt","initMuon", muon->pt() );
      if(muon->pt() > 25 && (fabs(muon->eta())<2.4 ) )
	{Presel=true;
	if(pttmp<muon->pt())
	  { BestLepton_=im; pttmp = muon->pt();
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
    LCBM_[0][0]=theLepton_->cosPhi();
    LCBM_[0][1]=theLepton_->sinPhi();
    LCBM_[1][0]=-theLepton_->sinPhi();
    LCBM_[1][1]=theLepton_->cosPhi();
    
    EP_="EB"; EcalP_=0;
    if( (fabs(eta)> 1.479 && Wmode=="Wen") || 
	(fabs(theLepton_->eta() )> 1.3 && Wmode=="Wmn") )
      { EP_="EE";EcalP_=1;}
  }
  if(Presel) {
    accept=true;
  }
  
  return accept;
}


bool
OfficialWAnalysis::ConversionRejection() {

  //Rejection options
  bool firstPixel =true;
  bool expFirstPixel =true;
  bool EGamConv = true;
  
  bool Conv[3]={true,true,true};


  Candidate* gsfTrack_ = _e.getSecond( "electron-GsfTrack", theLepton_);
  CandInfo* info_ = gsfTrack_->info();
 
  //possibility, use first pixel hit
  if(firstPixel){
    bool FHP_ = info_->getBool("isFPHitTrk");
    if(!FHP_) /*return false; //*/
      Conv[0] = false;
  }

  //other Possibility , using expected hits in pixel layers
  if(expFirstPixel){
    int EIP_ = info_->getInt("expInnHits");
    int maxExpInHit =0;
    if(EIP_ > maxExpInHit)
      Conv[1] = false;
  }
  //other possibility, using EGamma rejection, not available with Summer09 datasets
  if(EGamConv) {
    CandInfo* elInfo = theLepton_->info();
    bool isConv = elInfo->getBool("isConv");
    
    if(isConv)
      Conv[2] = false;
  }
  
  inHit = Conv[0];
  expInHit = Conv[1];
  ConvRej = Conv[2];

  return Conv[1];
}


bool
OfficialWAnalysis::LeptonID() {

  CandInfo* info_ = theLepton_->info();

  float etaSep;

  if(Wmode=="Wen") {
    
    etaSep = 1.479;

    float hoe = info_->getFloat("HOE");
    float dphi = info_->getFloat("dPhiIn");
    float deta = info_->getFloat("dEtaIn");
    float sigieie = info_->getFloat("sigiEtaiEta");

    const CandList& tracks_ = _e.tracks();
    for(int unsigned i=0;i<tracks_.size();i++) {
      Candidate* track = tracks_[i];
      CandInfo* tinfo = track->info();
      fill("d0","tracks",tinfo->getFloat("d0"));
    }

    fill("HOE","BeginID_"+EP_, hoe );
    fill("dPhiIn","BeginID_"+EP_, dphi );
    fill("dEtaIn","BeginID_"+EP_, deta );
    fill("sigiEtaiEta","BeginID_"+EP_, sigieie );
  
    //0 barrel, 1 endcap
    float dEtaIn[2];
    float dPhiIn[2];
    float sigiEiE[2];
    float HOE[2];

    //Working Point 95 %
    /*  dEtaIn[0]= 0.007; dEtaIn[1]= 0.01; 
    dPhiIn[0]= 0.8; dPhiIn[1]= 0.7; 
    sigiEiE[0]= 0.01; sigiEiE[1]= 0.03;
    HOE[0]= 0.15; HOE[1]= 0.07;*/

    dEtaIn[0]= 0.004; dEtaIn[1]= 0.007; 
    dPhiIn[0]= 0.06; dPhiIn[1]= 0.03; 
    sigiEiE[0]= 0.01; sigiEiE[1]= 0.03;
    HOE[0]= 0.04; HOE[1]= 0.025;

  //*****************


    bool acc=false;
    
    if(etaSep > fabs(theLepton_->eta() ) ) { //barrel
    
      if( hoe < HOE[0] && sigieie < sigiEiE[0] && 
	  fabs(deta) < dEtaIn[0] && fabs(dphi) < dPhiIn[0])
	{
	  acc =true;
	}
    }
    else { //Endcap
      
      if( hoe < HOE[1] && sigieie < sigiEiE[1] && 
	  fabs(deta) < dEtaIn[1] && fabs(dphi) < dPhiIn[1] )
	{
	  acc =true;
	}
    }
    
    return acc;
  }
  else // W->mu
    {
      return info_->getBool("muidGlobalMuonPromptTight");  
    }
  
}

bool
OfficialWAnalysis::LeptonIso() {

  CandInfo* info_ = theLepton_->info();

  float etaSep;

  if(Wmode=="Wen") {
    
    etaSep = 1.479;
    
    float TrkIso[2];
    float EcalIso[2];
    float HcalIso[2];

    float trkIso = info_->getFloat("dr03TrkIso")/theLepton_->pt(); //*/info_->getFloat("dr03TrkIso") ;
    float ecalIso = info_->getFloat("dr03EcalIso")/theLepton_->pt(); //*/info_->getFloat("dr04EcalIso");
    float hcalIso = (info_->getFloat("dr03HcalD1Iso") + info_->getFloat("dr03HcalD2Iso"))/theLepton_->pt(); //*/info_->getFloat("dr04HcalD1Iso") + info_->getFloat("dr04HcalD2Iso") ;
  
    //Working point 95%
    /* TrkIso[0]= 0.15; TrkIso[1]= 0.08;
    EcalIso[0]= 2.0; EcalIso[1]= 0.06;
    HcalIso[0]= 0.12; HcalIso[1]= 0.05;*/

    TrkIso[0]= 0.09; TrkIso[1]= 0.04;
    EcalIso[0]= 0.07; EcalIso[1]= 0.05;
    HcalIso[0]= 0.1; HcalIso[1]= 0.025;

    //***********

    fill( "Iso", "NTrkIso"+EP_,trkIso  );
   fill( "Iso", "NEcalIso"+EP_,ecalIso);
   fill( "Iso", "NHcalIso"+EP_, hcalIso);

  bool acc=false;
    
    if(etaSep > fabs(info_->getFloat("caloEta") ) ) { //barrel
      if(trkIso < TrkIso[0] && ecalIso < EcalIso[0] &&
	 hcalIso < HcalIso[0] )
	{
	  acc =true;
	}
    }
    else { //Endcap
      if(trkIso < TrkIso[1] && ecalIso < EcalIso[1] &&
	 hcalIso < HcalIso[1] )
	{
	  acc =true;
	}
    }
    return acc; 
  }      
  else // W->mu
    {
      return true;
    }
}

bool
OfficialWAnalysis::SecondEcut() {
  
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
OfficialWAnalysis::FillingHistos(const string & vbtfSel) {
  
  //Filling transverse mass
  const Candidate* lepton = theLepton_;
  Candidate* WCand = Candidate::create(lepton, met_ );

  fill( "mass", "WTransverseMass"+EP_+"_"+vbtfSel, WCand->mass() );
  fill( "pt", "Wpt"+EP_+"_"+vbtfSel, WCand->pt() );
  fill( "eta", "Weta"+EP_+"_"+vbtfSel, WCand->eta() );
  fill( "phi", "Wphi"+EP_+"_"+vbtfSel, WCand->phi() );

  fill( "dPhi", "MetElec"+EP_+"_"+vbtfSel,  KineUtils::dPhi(met_->phi(),  lepton->phi() )*180/3.1415 );
 
  //fill( "dPhi", "MetJet"+EP_, KineUtils::dPhi(theJet_->phi(),met_->phi() )*180/3.1415 ); 

  fill( "eta", "lepton_selected_"+vbtfSel, lepton->eta() );
  fill( "METvsLepEta", "lepton_selected_"+vbtfSel, met_->pt() , lepton->eta() );
  fill( "MassvsLepEta", "lepton_selected_"+vbtfSel, WCand->mass() , lepton->eta() );
 

  fill( "pt", "MET_selected"+EP_+"_"+vbtfSel, met_->Et() );
  fill( "phi", "MET_selected"+EP_+"_"+vbtfSel, met_->phi() );
  fill( "pt", "lepton_selected"+EP_+"_"+vbtfSel, lepton->pt() );
  fill( "phi", "lepton_selected"+EP_+"_"+vbtfSel, lepton->phi() );
  fill( "charge", "lepton"+EP_+"_"+vbtfSel, lepton->charge() );

  //Fill underlying event
  FillUnderlyingPhoton(WCand);
  FillUnderlyingJets(WCand);

   
  //MET Commissioning
  string type = "End_" +vbtfSel + "_";
  FillSeveralMET( type );

 

}



void
OfficialWAnalysis::FillUnderlyingPhoton(Candidate* WCand ) {

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
OfficialWAnalysis::FillUnderlyingJets(Candidate* WCand ) {

  const CandList& Jets = _e.jetList( EventManager::kPatJet);
  Candidate* theJet_;

  float pttmp=0;
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
    fill( "phi", Wmode+"_LeadJet"+EP_+IsoCateg_, theJet_->phi( kDeg ) );
    fill( "pt", Wmode+"_LeadJet"+EP_+IsoCateg_ , theJet_->pt() );
    fill( "dphi", Wmode+"_LeadJet"+EP_+IsoCateg_, KineUtils::dPhi(WCand->phi(), theJet_->phi() )*180/3.1415 );
  }
  
}

void
OfficialWAnalysis::FillSeveralMET(const string & ALevel) {

  //string type;

  //With PF MET
  const Candidate* pfmet_ = _e.met( EventManager::kPfMet );
  FillMETVariables(pfmet_, "PF"+ALevel);

  /* if(pfmet_->pt() > 60 ) {
    cout<<" ==================================== "<<endl;
    cout<<" *********** Event found ************ "<<endl;
    cout<<" ==================================== "<<endl;
    cout<< _e.run() <<"  "
	<< _e.event() <<"  "
	<< _e.file() <<"  "
	<< _e.eventInFile()<<"  "<<endl;
  }*/

  //With other MET : Calo Corrected or Raw Calo =======
  _e.refreshCaloMET(); //to be sure we take the good one
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

  //With Calo MET
  _e.refreshCaloMET(); //to be sure we take the good one
  int v=-1;
  for(int unsigned i=0;i<_e.a().caloMET.size();i++) {
    if("met" ==_e.a().caloMET[i])
      { v=i; break; }
  }
 
  const Candidate* cmet_ = _e.met( EventManager::kCaloMet, v );
  FillMETVariables(cmet_, "RawCalo"+ALevel);
  //========================================

 //With Calo MET
  _e.refreshCaloMET(); //to be sure we take the good one
  v=-1;
  for(int unsigned i=0;i<_e.a().caloMET.size();i++) {
    if("metMuonJESCorAK5" ==_e.a().caloMET[i])
      { v=i; break; }
  }
 
  const Candidate* cmetT2_ = _e.met( EventManager::kCaloMet, v );
  FillMETVariables(cmetT2_, "CaloTypeII"+ALevel);
  //========================================

} 



TVector2
OfficialWAnalysis::ComputeRecoil(const Candidate* met) {
  
  TVector2 Recoil(0,0);
  Recoil -=  met->p2() + theLepton_->p2();
  return Recoil;

}

TVector2
OfficialWAnalysis::ComputeRecoilWithPF(const Candidate* met) {
  
  TVector2 Recoil(0,0);
  Candidate* pfcand_;
  if(Wmode=="Wen")
    pfcand_ = _e.getSecond("electron-PfCand", theLepton_ );
  if(Wmode=="Wmn")
    pfcand_ = _e.getSecond("muon-PfCand", theLepton_ );
  
if(pfcand_!=NULL)
   Recoil -=  met->p2() + pfcand_->p2();
 return Recoil;

}

TVector2
OfficialWAnalysis::ComputeRecoilWithCaloRecHit(const Candidate* met) {
  
  TVector2 Recoil(0,0);

  TVector2 sumRH_(0,0);
  if(Wmode=="wen") {
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
 
  Recoil -=  met->p2() +  sumRH_;
  }
  return Recoil;

}

TVector2
OfficialWAnalysis::ComputeRecoilWithCalo(const Candidate* met) {
  
  TVector2 Recoil(0,0);
  if(Wmode=="Wen"){
  Candidate* superCluster_ = _e.getSecond("electron-SuperCluster", theLepton_ );
  if(superCluster_!=NULL)
    Recoil -= met->p2() + superCluster_->p2();
  }
  
  return Recoil;

}

TVector2
OfficialWAnalysis::ComputeRecoilWithTowers(const Candidate* met) {

  TVector2 Recoil(0,0);

  const CandList& caloTowers_ = _e.caloTowers();
  float dR;
  for(int unsigned i=0;i<caloTowers_.size();i++) {
    Candidate* caloTower = caloTowers_[i];

    dR = KineUtils::dR(theLepton_->eta(),caloTower->eta(),theLepton_->phi(),caloTower->phi());
    if(dR<0.3)
      Recoil -= caloTower->p2();
  }
  Recoil -= met->p2();
  return Recoil;
  
}



TVector2
OfficialWAnalysis::ComputeEtVectorFromRecHit(ecalRecHit rh) {


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
OfficialWAnalysis::FillRecoil(TVector2 Recoil, Candidate* W, const string & type){

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
      //  cout<<type<<"  "<< Rpara<<"  "<<RSCpara<<"   "<<theLepton_->phi()<<"  "<<superCluster_->phi()<<endl;
      fill("noAbsPt",type+"Recoil_ParaCalo", RSCpara);
      fill("noAbsPt",type+"Recoil_PerpCalo", RSCperp);
    }
  }
}


void
OfficialWAnalysis::FillMETProjection(const Candidate* met, Candidate* W, const string & type) {

  //MET Projection on lepton axis
  float METpara= LCBM_[0][0]*met->px()+LCBM_[0][1]*met->py();
  float METperp= LCBM_[1][0]*met->px()+LCBM_[1][1]*met->py();
  
  fill("noAbsPt",type+"MET_ParaLepton", -METpara);
  fill("noAbsPt",type+"MET_PerpLepton", -METperp);


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
OfficialWAnalysis::FillMETVariables(const Candidate* met, const string & type) {

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

  
  TVector2 CaloRecoil = ComputeRecoilWithCalo(met);
  TVector2 PFRecoil = ComputeRecoilWithPF(met);
  //  TVector2 CalRecHitRecoil = ComputeRecoilWithCaloRecHit(met);
  // if( type.substr(0,5)=="CaloE")
  // cout<<type<<"  "<<CaloRecoil.Px()<<"  "<<Recoil.Px()<<"  "<<CalRecHitRecoil.Px()<<endl;
 
  if( type.substr(0,2)!="PF" && CaloRecoil.Px()!=0 && CaloRecoil.Py()!=0 ) {
    // FillRecoil(Recoil, W, type);
    FillRecoil(CaloRecoil, W, type+"Calo");
    // FillRecoil(CalRecHitRecoil, W, type+"RecHit");
  }
  if( type.substr(0,2)=="PF" && PFRecoil.Px()!=0 && PFRecoil.Py()!=0 ) {
    //FillRecoil(Recoil, W, type);
    FillRecoil(PFRecoil, W, type+"PF");
  }
    
}

void
OfficialWAnalysis::VBTFSelection() 
{
 
  //  string EtT[3]={"Et20_","Et25_","Et30_"};

  /* for(unsigned iet=0;iet<VBTFElectronSelector::kNEt;iet++) {
      
      for( unsigned sel_=1; sel_<VBTFElectronSelector::kNSel; sel_++ )
	{

	  VBTFElectronSelector vbtfSel(iet, sel_);
	  vbtfSel.useOfficialIsoComputer = true;
	  if( vbtfSel.accept( *theLepton_ ) )
	    {
	      string elType = EtT[iet] + VBTFElectronSelector::level[ sel_ ];
	      FillingHistos(elType);

	      // FillASCII(elType);
	    }	  
	}
    }
*/

  /* unsigned etcut_ = VBTFElectronSelector::kEt20;
  unsigned sel_;
  sel_= VBTFElectronSelector::kSel95;		       
  string level_ =  VBTFElectronSelector::level[sel_];
  VBTFElectronSelector vbtfSel(etcut_,sel_); 
  */
  // vbtfSel.useOfficialIsoComputer = true;
  // vbtfSel.setFlags( *theLepton_ );

  bool selection[6]={false, false, false, false, false, false};
  bool eid[6]={false, false, false, false, false, false};
  bool iso[6]={false, false, false, false, false, false};
  bool fidu[2]={false,false};
  bool Et[3]={false, false, false};
  static string EtT[3]={"Et20_","Et25_","Et30_"};
  static string Selec[6]={ "Sel60", "Sel70","Sel80", "Sel85", "Sel90", "Sel95"};
  // string Eid[6]={ "Eid60", "Eid70", "Eid80", "Eid85", "Eid90", "Eid95"};
  // string Iso[6]={ "Iso60", "Iso70", "Iso80", "Eid85", "Iso90", "Iso95"};

  if(vbtfSel.accept( *theLepton_ )) {

    CandInfo* info = theLepton_->info();
    
    selection[0] = info->getBool( "sel60" );
    selection[1] = info->getBool( "sel70" );
    selection[2] = info->getBool( "sel80" );
    selection[3] = info->getBool( "sel85" );
    selection[4] = info->getBool( "sel90" );
    selection[5] = info->getBool( "sel95" );
    //  selection[5] = info->getBool( "selAnti" );

    eid[0] = info->getBool( "eid60" );
    eid[1] = info->getBool( "eid70" );
    eid[2] = info->getBool( "eid80" );
    eid[3] = info->getBool( "eid85" );
    eid[4] = info->getBool( "eid90" );
    eid[5] = info->getBool( "eid95" );
    //  eid[4] = info->getBool( "eidAnti" );

    iso[0] = info->getBool( "iso60" );
    iso[1] = info->getBool( "iso70" );
    iso[2] = info->getBool( "iso80" );
    iso[3] = info->getBool( "iso85" );
    iso[4] = info->getBool( "iso90" );
    iso[5] = info->getBool( "iso95" );
    //   iso[6] = info->getBool( "isoAnti" );
    
    fidu[0] = info->getBool("isEB");
    fidu[1] = info->getBool("isEE");
    
    Et[0] = info->getBool("isEt20");
    Et[1] = info->getBool("isEt25");
    Et[2] = info->getBool("isEt30");
    

    for(int iet=0;iet<3;iet++) {
      
      if(Et[iet]) {
	for(int isel=0;isel<6;isel++) {
	 
	  if(selection[isel]) {
	    
	    static string elType = EtT[iet] + Selec[isel];
	    FillingHistos(elType);
	    
	  }
	}//End Selection
      }
    }// End Et
   
    //Special Selection
    if(Et[0]) {

      float hoe = info->getFloat("HOE");
      float dphi = info->getFloat("dPhiIn");
      float deta = info->getFloat("dEtaIn");
      float sig = info->getFloat("sigiEtaiEta");

      if(iso[2] && met_->pt() > 40) {
	fill("HOE","Iso90_MetCut_"+EP_, hoe );
	fill("dPhiIn","Iso90_MetCut_"+EP_, dphi );
	fill("dEtaIn","Iso90_MetCut_"+EP_, deta );
	fill("sigiEtaiEta","Iso90_MetCut_"+EP_, sig );
      }
      
      if(iso[0]) {//Iso 70%
	
	fill("HOE","Iso70_"+EP_, hoe );
	fill("dPhiIn","Iso70_"+EP_, dphi );
	fill("dEtaIn","Iso70_"+EP_, deta );
	fill("sigiEtaiEta","Iso70_"+EP_, sig );

	static string tt="Iso70";
	FillingHistos(tt); 

	if(eid[2]) {
	  static string elType = "Iso70_eid90";
	  FillingHistos(elType); 
	  if(eid[0]) {
	    elType = "Iso70_eid70";
	    FillingHistos(elType); 
	  }
	}
      }//End Iso 70
      
      if(eid[0]) {//eid 70%

	static string tt="eid70";
	FillingHistos(tt); 

	if(iso[2]) {
	  static string elType = "eid70_Iso90";
	  FillingHistos(elType); 
	  if(iso[0]) {
	    elType = "eid70_Iso70";
	    FillingHistos(elType); 
	  }
	}
      } //end eid 70
    }//End Et
  }
}




void 
OfficialWAnalysis::InitASCIIS() {

  static string EtT[3]={"Et20_","Et25_","Et30_"};

  for(unsigned iet=0;iet<VBTFElectronSelector::kNEt;iet++) {
    
    for( unsigned sel_=1; sel_<VBTFElectronSelector::kNWP; sel_++ )
      {
	string elType = EtT[iet] + VBTFElectronSelector::level[ sel_ ];
	
	ofstream* ofile= new ofstream((elType+".txt").c_str(), ios::out | ios::trunc );
	if(ofile)
	  ASCIIMap[ elType ]  = ofile;
	else
	  {cout <<" error opening ascii file "<<endl; abort(); }
      }
  }
}


void 
OfficialWAnalysis::FillASCII(const string & type) {

  std::map<std::string, std::ofstream* >::iterator file  = ASCIIMap.find(type);
  

  // ostringstream os;
  // os <<  _e.run() + "_" +  _e.event() + "_"

  (* (*file).second) << _e.run() <<"  "
		     << _e.event() <<"  "
		     << _e.file() <<"  "
		     << _e.eventInFile()<<"  "<<endl;
  
}



void 
OfficialWAnalysis::initRootTree() {

  //Event
  outTree_.Branch("Run",&Run,"Run/I");
  outTree_.Branch("Event",&Event,"Event/I");
  outTree_.Branch("AbsEvent",&AbsEvent,"AbsEvent/I");
  outTree_.Branch("sample",&sample,"sample[21]/C"); 
  outTree_.Branch("fileName",&fileName,"fileName[21]/C"); 

  //Selection
  outTree_.Branch("VetoE",&VetoE,"VetoE/B");
  //outTree_.Branch("ConvRej",&ConvRej,"ConvRej/B");
  //outTree_.Branch("ConvRej",con,"con[3]/B");

  outTree_.Branch("InHit",&inHit,"inHit/B");
  outTree_.Branch("ExpInHit",&expInHit,"expInHit/B");
  outTree_.Branch("ConvRej",&ConvRej,"ConvRej/B");

  //W
  outTree_.Branch("WPt",WPt,"WPt[6]/F");
  outTree_.Branch("GenWPt",GenWPt,"GenWPt[2]/F");
  outTree_.Branch("GenEPt",GenEPt,"GenEPt[2]/F");
  outTree_.Branch("GenCharge",&GenCharge,"GenCharge/F");
  outTree_.Branch("GenEnergy",&GenEnergy,"GenEnergy/F");

  //DeltaZ between vertices
  outTree_.Branch("DZvertices",&dZ,"dZ/F");
  outTree_.Branch("zv2",&z2V,"z2V/F");

  //Jet
  outTree_.Branch("VetoJet",&VetoJet,"VetoJet/B");
  outTree_.Branch("JetMultiplicity",&JetMultiplicity,"JetMultiplicity/I");
  outTree_.Branch("LeadingJet",&LJetPt,"LJetPt/F");
  outTree_.Branch("dPhiJetW",dPhiJetW,"dPhiJetW[6]/F");
  
  //Photon
  outTree_.Branch("PhotonMultiplicity",&PhotonMultiplicity,"PhotonMultiplicity/I");
  outTree_.Branch("LeadingPhoton",&LPhotonPt,"LPhotonPt/F");
  outTree_.Branch("dPhiPhotonW",dPhiPhotonW,"dPhiPhotonW[6]/F");
  //  outTree_.Branch("",,"");

  //Angles
  outTree_.Branch("dPhiJetLep",&dPhiJetLep,"dPhiJetLep/F");
  outTree_.Branch("dPhiMETLep",dPhiMETLep,"dPhiMETLep[6]/F");
  outTree_.Branch("dPhiMETJet",dPhiMETJet,"dPhiMETJet[6]/F");

  outTree_.Branch("dPhiRecoilLep",dPhiRecoilLep,"dPhiRecoilLep[6]/F");
  outTree_.Branch("BasicdPhiRecoilLep",BasicdPhiRecoilLep,"BasicdPhiRecoilLep[6]/F");

  //Lepton
  outTree_.Branch("IDVar",IDVar,"IDVar[6]/F");
  outTree_.Branch("IsoVar",IsoVar,"IsoVar[3]/F"); 
  outTree_.Branch("Lepton",Lepton,"Lepton[3]/F");
  outTree_.Branch("PFLepton",PFLepton,"PFLepton[3]/F");
  outTree_.Branch("SCLepton",SCLepton,"SCLepton[3]/F");
  outTree_.Branch("EtSC",&EtSC,"EtSC/F");

  outTree_.Branch("Charge",&charge,"charge/I");
  
  //Init MET
  std::stringstream buffer;
  buffer << "E"
         <<"[" << sizeof(Mets)/sizeof(Mets[0]) << "]"
         << "[" << sizeof(Mets[0])/sizeof(Mets[0][0]) << "]/F";
  outTree_.Branch("Mets", Mets, buffer.str().c_str());
  buffer.str("");
  buffer << "E"
         <<"[" << sizeof(MetProj)/sizeof(MetProj[0]) << "]"
         << "[" << sizeof(MetProj[0])/sizeof(MetProj[0][0]) << "]/F";
  outTree_.Branch("MetProj", MetProj, buffer.str().c_str());
  buffer.str("");
  //Init Recoil
  buffer << "E"
         <<"[" << sizeof(Recoils)/sizeof(Recoils[0]) << "]"
         << "[" << sizeof(Recoils[0])/sizeof(Recoils[0][0]) << "]/F";
  outTree_.Branch("Recoils", Recoils, buffer.str().c_str());
  buffer.str("");
  buffer << "E"
         <<"[" << sizeof(RecoilProj)/sizeof(RecoilProj[0]) << "]"
         << "[" << sizeof(RecoilProj[0])/sizeof(RecoilProj[0][0]) << "]/F";
  outTree_.Branch("RecoilProj", RecoilProj, buffer.str().c_str());
  buffer.str("");

  //Init Basic Recoil
  buffer << "E"
         <<"[" << sizeof(BasicRecoils)/sizeof(BasicRecoils[0]) << "]"
         << "[" << sizeof(BasicRecoils[0])/sizeof(BasicRecoils[0][0]) << "]/F";
  outTree_.Branch("BasicRecoils", BasicRecoils, buffer.str().c_str());
  buffer.str("");
  buffer << "E"
         <<"[" << sizeof(BasicRecoilProj)/sizeof(BasicRecoilProj[0]) << "]"
         << "[" << sizeof(BasicRecoilProj[0])/sizeof(BasicRecoilProj[0][0]) << "]/F";
  outTree_.Branch("BasicRecoilProj", BasicRecoilProj, buffer.str().c_str());
  buffer.str("");
  buffer << "E"
         <<"[" << sizeof(GenRecoil)/sizeof(GenRecoil[0]) << "]"
         << "[" << sizeof(GenRecoil[0])/sizeof(GenRecoil[0][0]) << "]/F";
  outTree_.Branch("GenRecoil", GenRecoil, buffer.str().c_str());
  buffer.str("");
  buffer << "E"
         <<"[" << sizeof(CTRecoil)/sizeof(CTRecoil[0]) << "]"
         << "[" << sizeof(CTRecoil[0])/sizeof(CTRecoil[0][0]) << "]/F";
  outTree_.Branch("CTRecoil", CTRecoil, buffer.str().c_str());
  buffer.str("");
  //  outTree_.Branch("GenRecoil",GenRecoil,"GenRecoil[2]/F");
  
  outTree_.Branch("MT",MT,"MT[6]/F");
  
  outTree_.Branch("SumET",SumET,"SumET[6]/F");

  outTree_.Branch("NVertex",&Nvertex,"Nvertex/I");
  
  
}

void
OfficialWAnalysis::PrepareTree(float detacor, float dphicor ) {
  
  float ID[6];
  float Iso[3];
  
  CandInfo* info = theLepton_->info();
  
  if(Wmode=="Wen") {
 
  
  ID[0] = info->getFloat("sigiEtaiEta");
  ID[1] = info->getFloat("dEtaIn");
  ID[2] = info->getFloat("dPhiIn");
  ID[3] = info->getFloat("HOE");
  ID[4] = detacor;
  ID[5] = dphicor;

  IsoElec.ecalUseIsoDeposits = false;
  IsoElec.computeIso( *theLepton_ );

  Iso[0] = /*IsoElec.S_trk;//*/info->getFloat("dr03TrkIso") ;
  Iso[1] =/* IsoElec.S_ecal;//*/info->getFloat("dr03EcalIso");
  Iso[2] =/* IsoElec.S_hcal;//*/info->getFloat("dr03HcalD1Iso") + info->getFloat("dr03HcalD2Iso") ;

  

  }
  else { //muon

    	ID[0] = info->getFloat("normChi2");
	ID[1] = (float)info->getInt("nValMuHits");
	ID[2] = (float)info->getInt("nValTrackerHits");
	ID[3] = info->getFloat("dxy");
	ID[4] = (float)info->getBool("muidGlobalMuonPromptTight");
	ID[5] = (float)info->getBool("muidTrackerMuonArbitrated");

	Iso[0] = info->getFloat("isoR03strk");
	Iso[1] = info->getFloat("isoR03strkVeto");
	Iso[2] = info->getFloat("isoR03sum");

	inHit    = 0;
	expInHit = 0;
	ConvRej  = 0;

  }

  fillTree(Iso, ID);

}


void
OfficialWAnalysis::fillTree(float Iso[3], float ID[6]) {

  for(int i=0;i<3;i++) {
    IsoVar[i] = Iso[i];
    
  }
  for(int i=0;i<6;i++) {
    IDVar[i] = ID[i];
  }

  TreeFillLepton();
  TreeFillUnderlyingEvent();
  if(SampleName.substr(0,4)=="W_en")
    ComputeGenRecoil();
  else {
    GenRecoil[0][0] =0;
    GenRecoil[0][1] =0;
    GenRecoil[1][0] =0;
    GenRecoil[1][1] =0;
  }
  //Event
  Run = _e.run();
  Event = _e.eventInFile();
  AbsEvent = _e.event();
 
  strcpy (sample , ((string)SampleName).c_str() );
  
  int n = (Config::dataPath).size() + SampleName.size() +1;
  TString tmpName = _e.fileName();
  strcpy (fileName, (((string)tmpName).substr(n, ((string)tmpName).size()-n )).c_str() );
  
  // CurEvent = _e.eventInFile();
  // File = _e.event();

  //VetoJet
  VetoJet = VetoJetFunction();

  //Selection
  VetoE = SecondEcut();  //true c'ets que ça passe la coupure
  bool convrej =false;
  if(Wmode=="Wen")
    convrej= ConversionRejection();
  //ConvRej = ConversionRejection(); //ConvRej = convrej;   //true c'ets que ça passe la coupure
  // cout<<" ConvRej "<<ConvRej[0]<<"  "<<ConvRej[1]<<endl;
    int METType=0;

    //With PF MET
    {
      METType=0;
      const Candidate* pfmet_ = _e.met( EventManager::kPfMet );
      TVector2 PFRecoil = ComputeRecoilWithPF(pfmet_);
      TVector2 Recoil = ComputeRecoil(pfmet_);
      TVector2 CTRecoil = ComputeRecoilWithTowers(pfmet_);

      TreeFillMETProjection(pfmet_,METType);
      TreeFillRecoil(PFRecoil,METType);
      TreeFillUnderlyingJetsLepton(pfmet_,METType);
      TreeFillUnderlyingPhoton(pfmet_,METType);
      TreeFillBasicRecoil(Recoil,METType);
      TreeFillCTRecoil(CTRecoil,METType);
    }
    //With Track Corr Met
    {
      METType=1;
      const Candidate* trkcmet_  = _e.met( EventManager::kPatMet);
      TVector2 Recoil = ComputeRecoil(trkcmet_);
      TVector2 CTRecoil = ComputeRecoilWithTowers(trkcmet_);
      
      TreeFillMETProjection(trkcmet_,METType);
      TreeFillRecoil(Recoil,METType);
      TreeFillUnderlyingJetsLepton(trkcmet_,METType);
      TreeFillUnderlyingPhoton(trkcmet_,METType);
      TreeFillBasicRecoil(Recoil,METType);  
      TreeFillCTRecoil(CTRecoil,METType);

    }
    //With TypeI Calo =======
    {
      METType=3;
      _e.refreshCaloMET(); //to be sure we take the good one
      const Candidate* calomet_  = _e.met( EventManager::kCaloMet);
      TVector2 CaloRecoil = ComputeRecoilWithCalo(calomet_);
      TVector2 Recoil = ComputeRecoil(calomet_);
      TVector2 CTRecoil = ComputeRecoilWithTowers(calomet_);

      TreeFillMETProjection(calomet_,METType);
      TreeFillRecoil(CaloRecoil,METType);
      TreeFillUnderlyingJetsLepton(calomet_,METType);
      TreeFillUnderlyingPhoton(calomet_,METType);
      TreeFillBasicRecoil(Recoil,METType);
      TreeFillCTRecoil(CTRecoil,METType);
      
    }
    //With Raw Calo
    {
      METType=2;
      _e.refreshCaloMET(); //to be sure we take the good one
      int v=-1;
      for(int unsigned i=0;i<_e.a().caloMET.size();i++) {
	if("met" ==_e.a().caloMET[i])
	  { v=i; break; }
      }
      const Candidate* cmet_ = _e.met( EventManager::kCaloMet, v );
      TVector2 CaloRecoil = ComputeRecoilWithCalo(cmet_);
      TVector2 Recoil = ComputeRecoil(cmet_);
      TVector2 CTRecoil = ComputeRecoilWithTowers(cmet_);

      TreeFillMETProjection(cmet_,METType);
      TreeFillRecoil(CaloRecoil,METType);
      TreeFillUnderlyingJetsLepton(cmet_,METType);
      TreeFillUnderlyingPhoton(cmet_,METType);
      TreeFillBasicRecoil(Recoil,METType);
      TreeFillCTRecoil(CTRecoil,METType);
      
    }
     //With TypeII Calo
    {
      METType=4;
      _e.refreshCaloMET(); //to be sure we take the good one
      int v=-1;
      for(int unsigned i=0;i<_e.a().caloMET.size();i++) {
	if("metMuonJESCorAK5" ==_e.a().caloMET[i])
	  { v=i; break; }
      }
      const Candidate* cmetT2_ = _e.met( EventManager::kCaloMet, v );
      TVector2 CaloRecoil = ComputeRecoilWithCalo(cmetT2_);
      TVector2 Recoil = ComputeRecoil(cmetT2_);
      TVector2 CTRecoil = ComputeRecoilWithTowers(cmetT2_);

      TreeFillMETProjection(cmetT2_,METType);
      TreeFillRecoil(CaloRecoil,METType);
      TreeFillUnderlyingJetsLepton(cmetT2_,METType);
      TreeFillUnderlyingPhoton(cmetT2_,METType);
      TreeFillBasicRecoil(Recoil,METType);
      TreeFillCTRecoil(CTRecoil,METType);

    }
    //With TypeI PF
    {
      METType=5;
      _e.refreshRecoMET(); //to be sure we take the good one
      int v=-1;
      for(int unsigned i=0;i<_e.a().recoMET.size();i++) {
	if("metJESCorPFAK5" ==_e.a().recoMET[i])
	  { v=i; break; }
      }
      const Candidate* cmetPFT1_ = _e.met( EventManager::kPatMet, v );
      TVector2 CaloRecoil = ComputeRecoilWithCalo(cmetPFT1_);
      TVector2 Recoil = ComputeRecoil(cmetPFT1_);
      TVector2 CTRecoil = ComputeRecoilWithTowers(cmetPFT1_);

      TreeFillMETProjection(cmetPFT1_,METType);
      TreeFillRecoil(CaloRecoil,METType);
      TreeFillUnderlyingJetsLepton(cmetPFT1_,METType);
      TreeFillUnderlyingPhoton(cmetPFT1_,METType);
      TreeFillBasicRecoil(Recoil,METType);
      TreeFillCTRecoil(CTRecoil,METType);

    }

    // cout<<Mets[0][0]<<"   "<<Mets[1][0]<<"   "<<Mets[2][0]<<"   "
    //	<<Mets[3][0]<<"   "<<Mets[4][0]<<"   "<<Mets[5][0]<<endl;

    outTree_.Fill();
}


void
OfficialWAnalysis::TreeFillMETProjection(const Candidate* met, int metType) {

  //MET Projection on lepton axis
  float METpara= LCBM_[0][0]*met->px()+LCBM_[0][1]*met->py();
  float METperp= LCBM_[1][0]*met->px()+LCBM_[1][1]*met->py();

  MetProj[ metType ][0] = -METpara;
  MetProj[ metType ][1] = -METperp;

  Mets[ metType][0] = met->pt();
  Mets[ metType][1] = met->phi();
  
  CandInfo* info = met->info();
  SumET[ metType ] = info->getFloat("sumEt");

  Candidate* WCand = Candidate::create(met, theLepton_);
  MT[ metType] = WCand->mass();


  //W
  WPt[ metType ]=WCand->pt();


}

void
OfficialWAnalysis::TreeFillRecoil(TVector2 Recoil, int metType){

  if(Recoil.Px()!=0 && Recoil.Py()!=0 ) {

  //Recoil angle btw 0,2 pi
  Double_t rphi = TVector2::Phi_0_2pi(Recoil.Phi() );
 
//Recoil Projection
  float Rpara= LCBM_[0][0]*Recoil.Px()+LCBM_[0][1]*Recoil.Py();
  float Rperp= LCBM_[1][0]*Recoil.Px()+LCBM_[1][1]*Recoil.Py();
  

  Recoils[ metType][0] = Recoil.Mod();
  Recoils[ metType][1] = rphi;

  RecoilProj[ metType ][0] = Rpara;
  RecoilProj[ metType ][1] = Rperp;

  dPhiRecoilLep[ metType ] = KineUtils::dPhi(rphi,theLepton_->phi() )*180/Constants::pi;
 
  }
  else {
    Recoils[ metType][0] = -100000;
    Recoils[ metType][1] = -100000;
    
    RecoilProj[ metType ][0] = -100000;
    RecoilProj[ metType ][1] = -100000;
    
    dPhiRecoilLep[ metType ] = -100000;
  }
  
}

void
OfficialWAnalysis::TreeFillBasicRecoil(TVector2 Recoil, int metType){

  Double_t rphi = TVector2::Phi_0_2pi(Recoil.Phi() );
  float Rpara= LCBM_[0][0]*Recoil.Px()+LCBM_[0][1]*Recoil.Py();
  float Rperp= LCBM_[1][0]*Recoil.Px()+LCBM_[1][1]*Recoil.Py();

  BasicRecoils[ metType][0] = Recoil.Mod();
  BasicRecoils[ metType][1] = rphi;

  BasicRecoilProj[ metType ][0] = Rpara;
  BasicRecoilProj[ metType ][1] = Rperp;

  BasicdPhiRecoilLep[ metType ] = KineUtils::dPhi(rphi,theLepton_->phi() )*180/Constants::pi;

}

void
OfficialWAnalysis::TreeFillCTRecoil(TVector2 Recoil, int metType){

  //Double_t rphi = TVector2::Phi_0_2pi(Recoil.Phi() );
  float Rpara= LCBM_[0][0]*Recoil.Px()+LCBM_[0][1]*Recoil.Py();
  float Rperp= LCBM_[1][0]*Recoil.Px()+LCBM_[1][1]*Recoil.Py();

  CTRecoil[ metType ][0] = Rpara;
  CTRecoil[ metType ][1] = Rperp;
}


void
OfficialWAnalysis::TreeFillUnderlyingJetsLepton(const Candidate* met, int metType ) {

  const CandList& Jets = _e.jetList( EventManager::kPatJet);
  Candidate* theJet_;

  float pttmp=0;
  for(int unsigned ij=0;ij<Jets.size();ij++) { 
    Candidate* jet = Jets[ij];
    {
      float dr_ = KineUtils::dR(theLepton_->eta(), jet->eta(),
				theLepton_->phi(), jet->phi() );
      
      if(pttmp<jet->pt() && dr_>0.2  )
	{ theJet_=jet; pttmp = jet->pt(); }
    }
  }
  if( Jets.size()!=0 && pttmp!=0) {
    LJetPt = theJet_->pt();
    dPhiJetLep =  KineUtils::dPhi(theLepton_->phi(), theJet_->phi() )*180/Constants::pi;
 
    dPhiMETJet[metType] = KineUtils::dPhi(met->phi(), theJet_->phi() )*180/Constants::pi;
  }
  else {
    LJetPt = -1000;
    dPhiJetLep = -1000;
    dPhiMETJet[metType] = -1000;
  }

  //The W
  Candidate* WCand = Candidate::create(met, theLepton_);
  if( Jets.size()!=0 && pttmp!=0) {
  dPhiJetW[ metType ] = KineUtils::dPhi(WCand->phi(), theJet_->phi() );
  }
  else
    dPhiJetW[ metType ] = 0;
  
  //The lepton
  dPhiMETLep[metType] = KineUtils::dPhi(met->phi(), theLepton_->phi() )*180/Constants::pi;
  
}

void
OfficialWAnalysis::TreeFillUnderlyingPhoton(const Candidate* met, int metType ) {

  const CandList& Photons = _e.photons();
  Candidate* thePhoton_;

  float pttmp=0;
  for(int unsigned ij=0;ij<Photons.size();ij++) { 
    Candidate* photon = Photons[ij];
    {
      float dr_ = KineUtils::dR(theLepton_->eta(), photon->eta(),
				theLepton_->phi(), photon->phi() );
      
      if(pttmp<photon->pt() && dr_>0.2  )
	{ thePhoton_=photon; pttmp = photon->pt(); }
    }
  }
  if( Photons.size()!=0 && pttmp!=0) {
    LPhotonPt = thePhoton_->pt();
  }
  else {
    LPhotonPt = -1000;
  }

  //The W
  Candidate* WCand = Candidate::create(met, theLepton_);
  if( Photons.size()!=0 && pttmp!=0) {
  dPhiPhotonW[ metType ] = KineUtils::dPhi(WCand->phi(), thePhoton_->phi() );
  }
  else
    dPhiPhotonW[ metType ] = 0;
  
}

bool OfficialWAnalysis::VetoJetFunction() {

 
  //Jet
  const CandList& Jets = _e.jetList( EventManager::kPatJet);
  //Candidate* theJet_;

  for(int unsigned ij=0;ij<Jets.size();ij++) { 
    Candidate* jet = Jets[ij];
    {
      float dr_ = KineUtils::dR(theLepton_->eta(), jet->eta(),
				theLepton_->phi(), jet->phi() );
      
      if(jet->pt()> 25 && dr_>0.2  )
	return true;
    }
  }
  return false;
}

void OfficialWAnalysis::TreeFillUnderlyingEvent(){

  //Photon

  PhotonMultiplicity=0;

  const CandList& Photons = _e.photons();
  for(int unsigned ij=0;ij<Photons.size();ij++) { 
    Candidate* photon = Photons[ij];
    {
      float dr_ = KineUtils::dR(theLepton_->eta(), photon->eta(),
				theLepton_->phi(), photon->phi() );
      
      if(photon->pt()> 15 && dr_>0.2  )
	{
	  PhotonMultiplicity++;
	}
    }
  }
 
 //Jet

  JetMultiplicity=0;

  const CandList& Jets = _e.jetList( EventManager::kPatJet);
  for(int unsigned ij=0;ij<Jets.size();ij++) { 
    Candidate* jet = Jets[ij];
    {
      float dr_ = KineUtils::dR(theLepton_->eta(), jet->eta(),
				theLepton_->phi(), jet->phi() );
      
      if(jet->pt()> 25 && dr_>0.2  )
	{
	  JetMultiplicity++;
	}
    }
  }
  
}

void OfficialWAnalysis::TreeFillLepton() {

  Lepton[0] = theLepton_->pt();
  Lepton[1] = theLepton_->eta();
  Lepton[2] = theLepton_->phi();

  charge = (int)theLepton_->charge();

  CandInfo* info = theLepton_->info();
  if(Wmode=="Wen")
    EtSC = info->getFloat("caloEnergy")/cosh(info->getFloat("caloEta"));
  else
    EtSC = theLepton_->pt();

  Candidate* pfCand = _e.getSecond("electron-PfCand",theLepton_);
  if(pfCand!=NULL)
    {
      PFLepton[0] = pfCand->pt();
      PFLepton[1] = pfCand->eta();
      PFLepton[2] = pfCand->phi();
    }
  else
    {
      PFLepton[0]=-1; PFLepton[1]=0; PFLepton[2]=0;
    }

  Candidate* SCCand(0);
  if(Wmode=="Wen")
    SCCand= _e.getSecond("electron-SuperCluster",theLepton_);
  if(SCCand!=NULL)
    {
      SCLepton[0] = SCCand->pt();
      SCLepton[1] = SCCand->eta();
      SCLepton[2] = SCCand->phi();
    }
  else {
    SCLepton[0]=-1; SCLepton[1]=0; SCLepton[2]=0;
  }


}



bool OfficialWAnalysis::TreeFillVertices() {

  const VtxList& vertices = _e.vertices();
 
  Vertex* Primary; 
  bool findPV=false;
  //DZ vertices
  if(vertices.size() <2) 
    dZ = -1000;
  
  float dz=1000;
  Nvertex =0;
  for(int unsigned iv=0;iv<vertices.size();iv++) {
      
    Vertex* Vert = vertices[iv];
    CandInfo* infoV = Vert->info();

    bool goodV=false;
    if(fabs((Vert->pos()).Z())<24 && (Vert->pos()).Perp()<2 && infoV->getFloat("ndof")> 4 && !(infoV->getBool("isFake")) )
      goodV=true;

    if(!findPV && goodV )
      {Primary = Vert; Nvertex++; z2V = (Primary->pos()).Z(); findPV=true; continue;}
   
    if(goodV && findPV && vertices.size()>1) {
      Nvertex++;
      if(dz>fabs( (Primary->pos()).Z() - (Vert->pos()).Z()) )
	{ dz = fabs( (Primary->pos()).Z() - (Vert->pos()).Z());
	  dZ = (Primary->pos()).Z() - (Vert->pos()).Z();
	}
    }

  }
 
  if(findPV)
    return true;
  else
    return false;

}




void OfficialWAnalysis::ComputeGenRecoil() {

  Candidate* Wmc;
  Candidate* lepmc;

  bool find[2]={false,false};

  //Find the W
  const CandList& mc_ = _e.mcTruthCandidates();
  for(int unsigned i=0;i<mc_.size();i++) {
    Candidate* mc = mc_[i];

    if(abs(mc->pdgCode())==24)
      { Wmc = mc; find[0]=true; }
  
    float dR= KineUtils::dR(mc->eta(), theLepton_->eta(),mc->phi(), theLepton_->phi() );
    
    if( (abs(mc->pdgCode())==11 || abs(mc->pdgCode())==13 ) && dR<0.1)
      {lepmc = mc; find[1]=true;}
    
    if(find[0] && find[1])
      break;
  }


  if(lepmc!=NULL && Wmc!=NULL) {
     
  float genLCBM_[2][2];

    //Lepton Base Matrix
  genLCBM_[0][0]=lepmc->cosPhi();
    genLCBM_[0][1]=lepmc->sinPhi();
    genLCBM_[1][0]=-lepmc->sinPhi();
    genLCBM_[1][1]=lepmc->cosPhi();

    float Rpara= genLCBM_[0][0]*Wmc->px()+genLCBM_[0][1]*Wmc->py();
    float Rperp= genLCBM_[1][0]*Wmc->px()+genLCBM_[1][1]*Wmc->py();

    GenRecoil[0][0] = -Rpara;
    GenRecoil[0][1] = -Rperp;


    float Rparar= LCBM_[0][0]*Wmc->px()+LCBM_[0][1]*Wmc->py();
    float Rperpr= LCBM_[1][0]*Wmc->px()+LCBM_[1][1]*Wmc->py();
    
    GenRecoil[1][0] = -Rparar;
    GenRecoil[1][1] = -Rperpr;

    //W Gen Pt
    GenWPt[0] = Wmc->pt();
    GenWPt[1] = Wmc->phi();
    GenEnergy = lepmc->E();
     

    //Gen lepton
    GenEPt[0] = lepmc->pt();
    GenEPt[1] = lepmc->phi();
    GenCharge = lepmc->charge();
   
  }
  else {

    GenRecoil[0][0] = -1000;
    GenRecoil[0][1] = -1000;
    GenWPt[0] = -1000;
    GenWPt[1] = -1000;
    GenEnergy = -1000;
    GenEPt[0] = -1000;
    GenEPt[1] = -1000;
    GenCharge = -1000;
  }

}


bool
OfficialWAnalysis::SpikeCleaning() {


  bool spike=false;

   //Load RecHits
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
	  ecalRecHit rh = (*itrh).second.first;
	  if(rh.e > 10 && rh.iz ==0 ) {
	    // GHM fixme	    float swisscross = rh.SpikeVar.swCross;
	    float swisscross = 1.;
	    
	    if(swisscross > 0.95 )
	      { spike=true; return spike; }
	  }
	}
    }

  return spike;

}
