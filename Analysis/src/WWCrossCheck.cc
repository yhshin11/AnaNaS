#include "Analysis/src/WWCrossCheck.hh"

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
#include "Analysis/selectors/OfficialIsolationSelector.hh"
#include "Analysis/selectors/ConeSelector.hh"
#include "Analysis/utils/StatUtils.hh"
#include "Analysis/tools/PrintTree.hh"


int WWCrossCheck::imode = EventServer::kZmm;

ClassImp( WWCrossCheck )
  
  WWCrossCheck::WWCrossCheck( Sample& sample, const std::string & collectionFileName ):
  SampleAnalysis( "WWCCAnalysis" , sample, collectionFileName ),
  metNN_(0), Primary(0),
  outTree_("METCommVariables","Variables for MET Splots")
{
  cout << "\t-------------------------------------------" << endl; 
  cout << "\t---            ZZ  analysis            ----" << endl; 
  cout << "\t-------------------------------------------" << endl; 

  //Init compteurs
  Nevt_=0;
  NAcceptance_[0]=0;
  NAcceptance_[1]=0;
  NevtCand_=0;
  NHLT_ =0;
  NPresel_ =0;
  NIsol_ =0;
  NID_ =0;
  NMassCut_=0;
 
  Nz=0;

  SampleName = sample.name();

  SampleAnalysis::defineTemplate( "Iso", 500,0,5);
  SampleAnalysis::defineTemplate( "IsoVar", 501,0,1.002);
  SampleAnalysis::defineTemplate( "mass", 200, 0,200);
  
  SampleAnalysis::defineTemplate( "pdgId", 5001, -2500.5, 2500.5 );
  //  SampleAnalysis::defineTemplate( "EOP",251 , 0, 5.02);

  //IsoPlots
  SampleAnalysis::defineTemplate( "dlphi", 500, 0 , 1 ); 
  SampleAnalysis::defineTemplate( "dPhi", 361, -180.5, 180.5);
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

  SampleAnalysis::defineTemplate("massvsdpt",500,0,500,401,-2.005,2.005);

  SampleAnalysis::defineTemplate( "dPhivsMag",100,0,100,361,-180.5, 180.5);
  SampleAnalysis::defineTemplate( "dPhivsZMag",100,0,100,361,-180.5, 180.5);
  SampleAnalysis::defineTemplate( "noAbsPt", 400,-200,200);
  SampleAnalysis::defineTemplate( "ptVsPt", 500,0,500,500,0,500 );

  SampleAnalysis::defineTemplate( "UQtVsQt" ,500,0,500,100,0.,10 );

 
  //Get V for PFT1 MET
  _e.refreshRecoMET(); //to be sure we take the good one
  
  for(int unsigned i=0;i<_e.a().recoMET.size();i++) {
    if("metJESCorPFAK5" ==_e.a().recoMET[i])
      { _vmet=i; break; }
  }

  _verbose = true;


  InitializeElNNMET();

}

WWCrossCheck::~WWCrossCheck()
{
}

bool
WWCrossCheck::analyzeEvent()
{
  CatEta="";
  dau.clear();

  if((_ievt+1)%500 == 0)
    cout<<" Event "<<_ievt+1<<endl;
 
  if(!fillVertices()) return false;
  
  //Elements
  Nevt_++;
  
  //Need one Z candidate or more
  const CandList& zllList = MakeZCandidates();
  if( zllList.size()==0) return false;
  else  NevtCand_++;
  Candidate* ZCand = BestZCand(zllList);
  fill( "n", Zmode+"_ZCand", zllList.size() );
  fill( "mass",Zmode+"_Begin", ZCand->mass() );

   if(_e.inAcceptance() )
     if(Zmode=="Zemu")
       NAcceptance_[0]++;
     else
       NAcceptance_[1]++;

   //End preselection events ***********
  // identical to Z VBTF analysis *******

  //Analysis Way
  if(!HLTFiring() ) return false; //1st Step
  else NHLT_++;                   //Disabled for test

  if(!EventPreselection(ZCand)) return false; //2nd Step
  else  NPresel_++;
  
  if(!LeptonID() ) return false; //3rd step 
  else NID_++;
  
  if(!LeptonIso() ) return false; //4th step
  else NIsol_++;
  
  if(!MassSelection( ZCand) ) return false; //Just to give an approximation
  else NMassCut_++;                         // not implemented in Z analysis
  

  fillTree(ZCand);

  return true;
  
}

void
WWCrossCheck::bookHistograms()
{
  //
  // histogram booking
  //
  outTree_.SetDirectory( _outputFile->GetDirectory(NULL) );
  initRootTree();

}

void
WWCrossCheck::writeHistograms()
{
  outTree_.Write();

  cout<<" Number of events analyzed "<<Nevt_<<endl;
  cout<<" Number of events containing one or more ZCand "<<NevtCand_<<endl;
  cout<<" Number of events in Acceptance : elec "<<NAcceptance_[0]<<"  ; muon "<<NAcceptance_[1]<<"  "<<((double)(NAcceptance_[0]+NAcceptance_[1])*100)/Nevt_<<" %"<<endl;
  cout<<" Number of events firing HLT "<<NHLT_<<"  "<<((double)NHLT_*100)/(NAcceptance_[0]+NAcceptance_[1])<<" %;  "<<((double)NHLT_*100)/Nevt_<<endl;
  cout<<" Number of events preselected "<<NPresel_<<"  "<<((double)NPresel_*100)/NHLT_<<" %;  "<<((double)NPresel_*100)/Nevt_<<endl;
  cout<<" Number of events identified "<<NID_<<"  "<<((double)NID_*100)/NPresel_<<" %;  "<<((double)NID_*100)/Nevt_<<" %"<<endl;
  cout<<" Number of events Isolated "<<NIsol_<<"  "<<((double)NIsol_*100)/NID_<<" %;  "<<((double)NIsol_*100)/Nevt_<<" %"<<endl;
  cout<<" Number of events Passing";
  cout<<" trought the Mass window "<<NMassCut_<<"  "<<((double)NMassCut_*100)/NIsol_<<" %; "<<((double)NMassCut_*100)/Nevt_<<" %"<<endl;

}

bool
WWCrossCheck::HLTFiring() {

  bool fired=false;
  
  int run = _e.run();
  
  if(SampleName.substr(0,1)=="M" || SampleName.substr(0,1)=="E") { //FIXME MC...
    
    if( ( (run <= 140401) && _e.isFired("HLT_Ele15_LW_L1R") ) ||
	( (run >= 140402 && run <= 143962 ) && _e.isFired("HLT_Ele15_SW_L1R") ) ||
	( (run >= 143963 && run <= 144114 ) && _e.isFired("HLT_Ele15_SW_CaloEleId_L1R") ) ||
	( (run >= 144115 && run <= 147116 ) && _e.isFired("HLT_Ele17_SW_CaloEleId_L1R") ) ||
	( (run >= 147117 && run <= 148058 ) && _e.isFired("HLT_Ele17_SW_TightEleId_L1R") ) ||
	( (run >= 148059 && run <= 149064 ) && _e.isFired("HLT_Ele17_SW_TighterEleIdIsol_L1R_v2") ) ||
	( (run >= 149065 && run <= 149442 ) && _e.isFired("HLT_Ele17_SW_TighterEleIdIsol_L1R_v3") ) )
      fired =true;
    

    if( ( (run <= 147116) &&  _e.isFired("HLT_Mu9") ) ||
	( (run >= 147117) && (  _e.isFired("HLT_Mu15_v1") || _e.isFired("HLT_Mu15_v2") ) ) )
      fired =true;

    return fired;
  }
  else {
    //FIXME MC
    // return true;

    if(
       _e.isFired("HLT_Ele15_LW_L1R") ||
       _e.isFired("HLT_Ele15_SW_L1R") ||
       _e.isFired("HLT_Ele15_SW_CaloEleId_L1R") ||
       _e.isFired("HLT_Ele17_SW_CaloEleId_L1R") ||
       _e.isFired("HLT_Ele17_SW_TightEleId_L1R") ||
       _e.isFired("HLT_Ele17_SW_TighterEleIdIsol_L1R_v2") ||
       _e.isFired("HLT_Ele17_SW_TighterEleIdIsol_L1R_v3") ||
       _e.isFired("HLT_Mu9") ||
       _e.isFired("HLT_Mu15_v1") ||
       _e.isFired("HLT_Mu15_v2") )      
      fired=true;
  
  }
  
  return fired;

}
 

bool
WWCrossCheck::EventPreselection(Candidate* ZCand) {

  bool accept=false;
  
  bool Presel[2]={false,false};
  string categ_[2];
  
  float pt_0      =  ZCand->daughter(0)->pt();
  float pt_1      =  ZCand->daughter(1)->pt();
  
  if( pt_0>pt_1 )
    {
      dau.push_back( ZCand->daughter(0) );
      dau.push_back( ZCand->daughter(1) );
    }
  else
    {
      dau.push_back( ZCand->daughter(1) );
      dau.push_back( ZCand->daughter(0) );
    }

 if(abs(dau[0]->pdgCode())==11 && abs(dau[1]->pdgCode())==13)
    Zmode = "Zemu";
  else
    Zmode = "Zmue";


 for( int idau=0; idau<2; idau++ )
    {
      Candidate* dau_ = dau[idau];
      
      CandInfo* info = dau_->info();
      
      float eta_;
      float Et_; 
      if((Zmode=="Zemu" && idau==0) ||(Zmode=="Zmue" && idau==1) ) {
	eta_ = info->getFloat("caloEta");
	Et_	= info->getFloat("caloEnergy")/cosh(eta_);
      }
      else {
	eta_ = dau_->eta();
	Et_ = dau_->pt();
      }

      if( (idau==0 && ( (Zmode=="Zemu" && fabs(eta_)<1.479) || (Zmode=="Zmue" && fabs(eta_)<1.3 ) ) ) ||
	(idau==1 && ( (Zmode=="Zmue" && fabs(eta_)<1.479) || (Zmode=="Zemu" && fabs(eta_)<1.3 ) ) ) )
	categ_[idau]="B";
      else if( (idau==0 && ( (Zmode=="Zemu" && fabs(eta_)>=1.479) || (Zmode=="Zmue" && fabs(eta_)>=1.3 ) ) ) ||
	       (idau==1 && ( (Zmode=="Zmue" && fabs(eta_)>=1.479) || (Zmode=="Zemu" && fabs(eta_)>=1.3 ) ) ) )
	categ_[idau]="E";

      if( Et_> 20 && ((Zmode=="Zemu" && idau==0) || (Zmode=="Zmue" && idau==1) ) && ( fabs(eta_)<1.4442 || (fabs(eta_)>1.56 && fabs(eta_)<2.5) ) )
	Presel[idau]=true;
      if( Et_> 20 && ((Zmode=="Zemu" && idau==1) || (Zmode=="Zmue" && idau==0) ) && ( fabs(eta_)<2.1 ) )
	Presel[idau]=true;
    }
 
  CatEta=categ_[0]+categ_[1];
  if(CatEta=="BE")
    CatEta="EB";
  

  if(Presel[0] && Presel[1]) {
    accept=true;
  }
  
  return accept;
}

bool
WWCrossCheck::LeptonIso() {

  float TrkIso[2];
  float EcalIso[2];
  float HcalIso[2];

  float SumIso[2];

  float etaSep;
 
    //Electron Working Point 
  /*  TrkIso[0]= 0.15; TrkIso[1]= 0.05;
    EcalIso[0]= 2.0; EcalIso[1]= 0.06;
    HcalIso[0]= 0.12; HcalIso[1]= 0.05;*/
  TrkIso[0]= 0.09; TrkIso[1]= 0.04;
  EcalIso[0]= 0.07; EcalIso[1]= 0.05;
  HcalIso[0]= 0.10; HcalIso[1]= 0.025;
    //*****************

    //Muon isolation (relative)
    SumIso[0]=0.15; SumIso[1]=0.15;
  

    bool accept=false; 
  bool acc[2]={false,false};
  
  for( int idau=0; idau<2; idau++ )
    {
      Candidate* dau_ = dau[idau];
      CandInfo* info_ = dau_->info();
      
      if( (Zmode=="Zemu" && idau==0) ||
	  (Zmode=="Zmue" && idau==1) ) {
	
	etaSep = 1.479;

	float trkIso=0;
	float ecalIso=0;
	float hcalIso=0;
	
	trkIso = info_->getFloat("dr03TrkIso")/dau_->pt();
	ecalIso = info_->getFloat("dr03EcalIso")/dau_->pt();
	hcalIso = (info_->getFloat("dr03HcalD1Iso") + info_->getFloat("dr03HcalD2Iso") )/dau_->pt() ;
	  
	  
	if(etaSep > fabs(info_->getFloat("caloEta") ) ) { //barrel
	  if(trkIso < TrkIso[0] && ecalIso < EcalIso[0] &&
	     hcalIso < HcalIso[0] )
	    {
	      acc[idau] =true;
	    }
	}
	else { //Endcap
	  if(trkIso < TrkIso[1] && ecalIso < EcalIso[1] &&
	     hcalIso < HcalIso[1] )
	    {
	      acc[idau] =true;
	    }
	}
	
      }//muon  mode
      else {

	//No Iso for easier controla fter
	//FIXME
	//MM
	acc[idau]= true;
	/*
	etaSep = 1.3;

	float sumIso = info_->getFloat("isoR03strk")/dau_->pt(); //em cut only
	
	if(etaSep > fabs(dau_->eta())) { //barrel
	  if(sumIso < SumIso[0] )
	    acc[idau]=true;
	}
	else { //Endcap
	   if(sumIso < SumIso[1] )
	     acc[idau]=true;
	     }*/
      }//End Zmm
    }//End daughter loop

  if(acc[0] && acc[1])
    accept=true;

  return accept;


}


bool
WWCrossCheck::LeptonID() {

 //Cut definition ***
  //0 barrel, 1 endcap
  float dEtaIn[2];
  float dPhiIn[2];
  float sigiEiE[2];
  float HOE[2];

  // Electron Working Point
  dEtaIn[0]= 0.007; dEtaIn[1]= 0.01; //FIXME needed when alignment will be done 
  dPhiIn[0]= 0.8; dPhiIn[1]= 0.7; 
  sigiEiE[0]= 0.01; sigiEiE[1]= 0.03;
  HOE[0]= 0.15; HOE[1]= 0.07;
  //****************

  //Muon Cut ***
  bool muID[2] = {false,false};

  bool accept=false;
  bool acc[2]={false,false};
  for( int idau=0; idau<2; idau++ ) {
    Candidate* dau_ = dau[idau];
    CandInfo* info_ = dau_->info();
    if( (Zmode=="Zemu" && idau==0) || (Zmode=="Zmue" && idau==1) ) {

      float etaSep = 1.479;

         
      float hoe = info_->getFloat("HOE");
      float dphi = info_->getFloat("dPhiIn");
      float deta = info_->getFloat("dEtaIn");
      float sigieie = info_->getFloat("sigiEtaiEta");
      
      if(etaSep > fabs(info_->getFloat("caloEta") ) ) { //barrel
    
	if( hoe < HOE[0] && sigieie < sigiEiE[0] && 
	    fabs(deta) < dEtaIn[0] && fabs(dphi) < dPhiIn[0])
	  {
	    acc[idau] = true;
	  }
      }
      else { //Endcap
	if( hoe < HOE[1] && sigieie < sigiEiE[1] && 
	    fabs(deta) < dEtaIn[1] && fabs(dphi) < dPhiIn[1] )
	  {
	   acc[idau] = true;
	  }
      }
      
    }
    else { // Z -> mumu
      acc[idau] = true;
      float nvHits =  info_->getInt("nTrkHits"); //info_->getInt("nValMuHits") +
      
      if(info_->getBool("muidGlobalMuonPromptTight") == true && (nvHits>=10) )
        acc[idau] = true;
	
      if(MuonSpecID(dau_))
	muID[idau] = true;

    }
  }
  //cout<<Zmode<<"    "<<acc[0]<<"    "<<acc[1]
  //    <<"   "<<muID[0]<<"   "<<muID[1]<<endl;
  if(acc[0] && (acc[1] && muID[1]) && Zmode=="Zemu")
    accept=true;
  if( (acc[0] && muID[0]) && acc[1] &&  Zmode=="Zmue" )
    accept=true;

  return accept;

}

bool
WWCrossCheck::MuonSpecID(Candidate* mu) {

  CandInfo* info = mu->info();

  bool acc=true;

  if( fabs(mu->d0(Primary)) >= 0.2 )     acc=false; // cout<<" fi "<<acc<<"   "<<info->getFloat("dxy")<<endl;
  if(info->getFloat("normChi2") >= 10 ) acc=false;  //cout<<" fn "<<acc<<"   "<<info->getFloat("normChi2")<<endl;
  if(info->getInt("nMatch") <2 )    acc=false;  //cout<<" in "<<acc<<"    "<<info->getInt("nValMuHits")<<endl;
  if(info->getBool("muidTrackerMuonArbitrated")!= true)  acc=false; //  cout<<" f "<<acc<<endl;
  if(info->getInt("nPixHits") <1 ) acc=false;
  if(info->getInt("nMuonHits") <1 ) acc=false;
  
  return acc;

}

bool 
WWCrossCheck::MassSelection(Candidate* ZCand) {

  bool accept=false;

  if(ZCand->mass()> 40 && ZCand->mass()< 200 )
    accept=true;
  
  return accept;

}

Candidate *
WWCrossCheck::BestZCand( const CandList& zllList ) {

  Candidate* BCand;
  vector<float> tmp(2,0);
  vector<vector<float> > Pts;
  for(int unsigned iz=0;iz<zllList.size();iz++) {
   
    Pts.push_back(tmp);

    Candidate* ZCand = zllList[iz];
   
    float pt_0      =  ZCand->daughter(0)->pt();
    float pt_1      =  ZCand->daughter(1)->pt();
    //    cout<<ZCand->daughter(0)->pt()<<"   "<<ZCand->daughter(0)->pdgCode()<<"   <->   "<<ZCand->daughter(1)->pt()<<"   "<<ZCand->daughter(1)->pdgCode()<<endl;
    if( pt_0>pt_1 )
       {
	 Pts[iz][0]=pt_0;
	 Pts[iz][1]=pt_1;	 
       }
     else{
       Pts[iz][1]=pt_0;
       Pts[iz][0]=pt_1;	
     }
  }

  float pttmp=0;
  for(int unsigned ii=0;ii<Pts.size();ii++) {
    if(pttmp<Pts[ii][0]) pttmp=Pts[ii][0];
  }
  float pttmp2=0;  int fpt=0;
  for(int unsigned ii=0;ii<Pts.size();ii++) {
   if(Pts[ii][0]==pttmp) {
     if(pttmp2< Pts[ii][1]) {
       pttmp2=Pts[ii][1];
       fpt=ii;
     }
   }
  }

  BCand = zllList[fpt];


  //get the Type of the event
  int pdgId[2];
  pdgId[0] = abs(BCand->daughter(0)->pdgCode());
  pdgId[1] = abs(BCand->daughter(1)->pdgCode());


  //  cout<<" Selected --->   "<<BCand->daughter(0)->pt()<<"   "<<BCand->daughter(0)->pdgCode()<<"   <->   "<<BCand->daughter(1)->pt()<<"   "<<BCand->daughter(1)->pdgCode()<<endl;

  if(pdgId[0]==11 && pdgId[1]==13)
    Zmode = "Zemu";
  else if(pdgId[0]==13 && pdgId[1]==11)
    Zmode = "Zmue";
  else 
    { cout<<" Strange issue in Z creation, abort "<<endl; abort();}

  return BCand;
}


CandList WWCrossCheck::MakeZCandidates() {


  //Making Z Candidates
  CandList ZCands;

  const CandList electrons = _e.electrons();
  const CandList muons = _e.muons();

  for(int unsigned il1=0;il1<electrons.size();il1++) {
    Candidate* l1 = electrons[il1];
    
    for(int unsigned il2=0;il2<muons.size();il2++) {
      Candidate* l2 = muons[il2];
      Candidate* ZCand = Candidate::create(l1,l2);
      ZCands.push_back(ZCand);
    }
  }

  return ZCands;

}



CandList WWCrossCheck::MakeZCandidatesWithCharge() {

  const CandList* leptons;

  if(Zmode=="Zee")
    leptons = &(_e.electrons());
  else
    leptons = &(_e.muons());

  //Splitting l+ and l-
  CandList lplus;
  CandList lminus;
  for(int unsigned il=0;il<leptons->size();il++) {
    
    Candidate* lepton= (*leptons)[il];
    if(lepton->charge()==1)
      lplus.push_back(lepton);
    if(lepton->charge()==-1)
      lminus.push_back(lepton);
  }

  //Making Z Candidates
  CandList ZCands;
  for(int unsigned ilp=0;ilp<lplus.size();ilp++) {

    Candidate* lp = lplus[ilp];
    
    for(int unsigned ilm=0;ilm<lminus.size();ilm++) {
      Candidate* lm = lminus[ilm];
      Candidate* ZCand = Candidate::create(lp,lm);
      ZCands.push_back(ZCand);
    }
  }
  return ZCands;

}


TVector2
WWCrossCheck::ComputeRecoil(const Candidate* met, Candidate* ZCand) {
  
  TVector2 Recoil(0,0);
 
  Recoil -= met->p2() + ZCand->p2();
  return Recoil;

}



//*************************************=========================
//*************************************=========================

void 
WWCrossCheck::initRootTree() {
  
 //Event
  outTree_.Branch("Run",&Run,"Run/I");
  outTree_.Branch("Event",&Event,"Event/I");
  outTree_.Branch("AbsEvent",&AbsEvent,"AbsEvent/I");
  outTree_.Branch("sample",&sample,"sample[50]/C"); 
  outTree_.Branch("fileName",&fileName,"fileName[50]/C"); 

  //Category of the event
  outTree_.Branch("EvtCateg",EvtCateg,"EvtCateg[50]/C");

  // Leptons 
  outTree_.Branch("IDVar",IDVar,"IDVar[2][4]/F");
  outTree_.Branch("IsoVar",IsoVar,"IsoVar[2][3]/F"); 
  outTree_.Branch("Lepton",Lepton,"Lepton[2][4]/F");
  outTree_.Branch("SCLepton",SCLepton,"SCLepton[2][3]/F");
  outTree_.Branch("Charge",charge,"charge[2]/I");
  outTree_.Branch("PdgId",pdgId,"pdgId[2]/I");

  outTree_.Branch("EOP",EOP,"EOP[2]/F");

  outTree_.Branch("ElecMult",&ElecMult,"ElecMult/I");
  outTree_.Branch("MuonMult",&MuonMult,"MuonMult/I");
  outTree_.Branch("LeadAddElec",LeadAddElec,"LeadAddElec[3]/F");
  outTree_.Branch("LeadAddMuon",LeadAddMuon,"LeadAddMuon[3]/F");

  outTree_.Branch("IsIsoMu",&IsIsoMu,"IsIsoMu/B");
  outTree_.Branch("IsIDMu",&IsIDMu,"IsIDMu/B");
  outTree_.Branch("IsIsoElec",&IsIsoElec,"IsIsoElec/B");
  outTree_.Branch("IsIDElec",&IsIDElec,"IsIDElec/B");
  
  outTree_.Branch("ElecMT",&ElecMT,"ElecMT/F");
  outTree_.Branch("MuMT",&MuMT,"MuMT/F");

  outTree_.Branch("InHit1",&inHit1,"inHit1/B");
  outTree_.Branch("ExpInHit1",&expInHit1,"expInHit1/B");
  outTree_.Branch("ConvRej1",&ConvRej1,"ConvRej1/B");
  outTree_.Branch("InHit2",&inHit2,"inHit2/B");
  outTree_.Branch("ExpInHit2",&expInHit2,"expInHit2/B");
  outTree_.Branch("ConvRej2",&ConvRej2,"ConvRej2/B");

  //Z
  outTree_.Branch("Z",Z,"Z[4]/F");
  outTree_.Branch("ZMass",&ZMass,"ZMass/F");
  outTree_.Branch("ZMT",&ZMT,"ZMT/F");
  outTree_.Branch("ZCateg",&ZCateg,"ZCateg/I");
  outTree_.Branch("ZMultiplicity",&ZMultiplicity,"ZMultiplicity/I");
  outTree_.Branch("PfZ",PfZ,"PfZ[4]/F");
  outTree_.Branch("PfZMass",&PfZMass,"PfZMass/F");

  //MET
  outTree_.Branch("MET", MET,"MET[2]/F");
  outTree_.Branch("Recoil", Recoil,"Recoil[2]/F");

  outTree_.Branch("RValMC", RValMC,"RValMC[2]/F");
  outTree_.Branch("RValData", RValData,"RValData[2]/F");

  outTree_.Branch("SumEt", &SumEt,"SumEt/F");

  outTree_.Branch("TrkMET", TrkMET,"TrkMET[2]/F");
  outTree_.Branch("METPrime", METPrime,"METPrime[2]/F");
  outTree_.Branch("METPrimeRecoil", METPrimeRecoil,"METPrimeRecoil[2]/F");
  outTree_.Branch("METPrimeCor", METPrimeCor,"METPrimeCor[2]/F");
 
  outTree_.Branch("METPrimeMETM",METPrimeMETM,"METPrimeMETM[2]/F");
  outTree_.Branch("METPrimeJetM",METPrimeJetM,"METPrimeJetM[2]/F");
  outTree_.Branch("METPrimeJetP",METPrimeJetP,"METPrimeJetP[2]/F");
  outTree_.Branch("METPrimeUnclusM",METPrimeUnclusM,"METPrimeUnclusM[2]/F");
  outTree_.Branch("METPrimeUnclusP",METPrimeUnclusP,"METPrimeUnclusP[2]/F");

  
  outTree_.Branch("ProjMET",&ProjMET,"ProjMET/F");
  outTree_.Branch("CorProjMET",&CorProjMET,"CorProjMET/F");

  outTree_.Branch("outMETNN",&outMETNN,"outMETNN/F");

  //SpecMET var
  outTree_.Branch("SpecMET", &SpecMET, "SpecMET/F");

  //Jet
  outTree_.Branch("JetMultiplicity",&JetMultiplicity,"JetMultiplicity/I");
  outTree_.Branch("MultiJetMultiplicity",&MultiJetMult,"MultiJetMult[10]/I");
  outTree_.Branch("TcJetMultiplicity",&TcJetMultiplicity,"TcJetMultiplicity/I");
  outTree_.Branch("LeadingJet",LJet,"LJet[4]/F");
  outTree_.Branch("bTag",bTag,"bTag[2]/F");

  //Photon
  outTree_.Branch("PhotonMultiplicity",&PhotonMultiplicity,"PhotonMultiplicity/I");
  outTree_.Branch("LeadingPhoton",LPhoton,"LPhoton[5]/F");
  
  //Lepton related variables
  outTree_.Branch("dPhiLepton",&dPhiLepton,"dPhiLepton/F");
  outTree_.Branch("dPhiLeptonZ",&dPhiLeptonZ,"dPhiLeptonZ/F");
  outTree_.Branch("dEtaLepton",&dEtaLepton,"dEtaLepton/F");
  outTree_.Branch("dPtLepton",&dPtLepton,"dPtLepton/F");
  outTree_.Branch("dPLepton",&dPLepton,"dPLepton/F");
  outTree_.Branch("dRLepton",&dRLepton,"dRLepton/F");

  outTree_.Branch("dPhiLeptonS",&dPhiLeptonS,"dPhiLeptonS/F");
  outTree_.Branch("dEtaLeptonS",&dEtaLeptonS,"dEtaLeptonS/F");
  outTree_.Branch("dPtLeptonS",&dPtLeptonS,"dPtLeptonS/F");
  outTree_.Branch("dPLeptonS",&dPLeptonS,"dPLeptonS/F");
  outTree_.Branch("dRLeptonS",&dRLeptonS,"dRLeptonS/F");
  outTree_.Branch("CosThetaM",&CosThetaM,"CosThetaM/F");
  outTree_.Branch("CosThetaP",&CosThetaP,"CosThetaP/F");

  outTree_.Branch("CosThetaCSM",&CosThetaCSM,"CosThetaCSM/F");
  outTree_.Branch("CosThetaCSP",&CosThetaCSP,"CosThetaCSP/F");

  //MET related variables
  outTree_.Branch("ZMETMass",&ZMETMass,"ZMETMass/F");
  outTree_.Branch("ZMassMET",&ZMassMET,"ZMassMET/F");
  outTree_.Branch("ZPtMET",&ZPtMET,"ZPtMET/F");
  outTree_.Branch("ZPMET",&ZPMET,"ZPMET/F");
  outTree_.Branch("dPhiMETZ",&dPhiMETZ,"dPhiMETZ/F");
  outTree_.Branch("dPhiMETL1",&dPhiMETL1,"dPhiMETL1/F");
  outTree_.Branch("dPhiMETL2",&dPhiMETL2,"dPhiMETL2/F");
  outTree_.Branch("METProjZ",METProjZ,"METProjZ[2]/F");
  outTree_.Branch("METProjL1",METProjL1,"METProjL1[2]/F");
  outTree_.Branch("METProjL2",METProjL2,"METProjL2[2]/F");
  outTree_.Branch("ZRecoil",&ZRecoil,"ZRecoil/F");
  outTree_.Branch("RecoilProjZ",RecoilProjZ,"RecoilProjZ[2]/F");
  outTree_.Branch("RecoilProjL1",RecoilProjL1,"RecoilProjL1[2]/F");
  outTree_.Branch("RecoilProjL2",RecoilProjL2,"RecoilProjL2[2]/F");
  outTree_.Branch("dPhiRecoilZ",&dPhiRecoilZ,"dPhiRecoilZ/F");
  outTree_.Branch("dPhiRecoilL1",&dPhiRecoilL1,"dPhiRecoilL1/F");
  outTree_.Branch("dPhiRecoilL2",&dPhiRecoilL2,"dPhiRecoilL2/F");

  outTree_.Branch("HT",HT,"HT[2]/F");

  //Jet related variables
  outTree_.Branch("dPhiJetMET",&dPhiJetMET,"dPhiJetMET/F");
  outTree_.Branch("dPhiJetZ",&dPhiJetZ,"dPhiJetZ/F");
  outTree_.Branch("dEtaJetZ",&dEtaJetZ,"dEtaJetZ/F");
  outTree_.Branch("ZPtJet",&ZPtJet,"ZPtJet/F");
  outTree_.Branch("ZJetMass",&ZJetMass,"ZJetMass/F");
  outTree_.Branch("dPhiJetL1",&dPhiJetL1,"dPhiJetL1/F");
  outTree_.Branch("dPhiJetL2",&dPhiJetL2,"dPhiJetL2 /F");
  outTree_.Branch("dEtaJetL1",&dEtaJetL1,"dEtaJetL1/F");
  outTree_.Branch("dEtaJetL2",&dEtaJetL2,"dEtaJetL2/F");
  
  //Photon related variables
  outTree_.Branch("dPhiPhotonMET",&dPhiPhotonMET,"dPhiPhotonMET/F");
  outTree_.Branch("dPhiPhotonZ",&dPhiPhotonZ,"dPhiPhotonZ/F");
  outTree_.Branch("dEtaPhotonZ",&dEtaPhotonZ,"dEtaPhotonZ/F");
  outTree_.Branch("ZPtPhoton",&ZPtPhoton,"ZPtPhoton/F");
  outTree_.Branch("ZPhotonMass",&ZPhotonMass,"ZPhotonMass/F");
  outTree_.Branch("dPhiPhotonL1",&dPhiPhotonL1,"dPhiPhotonL1/F");
  outTree_.Branch("dPhiPhotonL2",&dPhiPhotonL2,"dPhiPhotonL2/F");
  outTree_.Branch("dEtaPhotonL1",&dEtaPhotonL1,"dEtaPhotonL1/F");
  outTree_.Branch("dEtaPhotonL2",&dEtaPhotonL2,"dEtaPhotonL2/F");

  //MC truth
  outTree_.Branch("McNeutrino",McNeutrino,"MCNeutrino[2][5]/F");
  outTree_.Branch("McLepton",McLepton,"MCLepton[2][5]/F");
  outTree_.Branch("McZn",McZn,"MCZn[4]/F");
  outTree_.Branch("McZl",McZl,"MCZl[4]/F");

  //Vertices
   outTree_.Branch("NVertex",&Nvertex,"Nvertex/I");
  
}

void
WWCrossCheck::fillTree(Candidate* ZCand) {

  fillLeptonTree(ZCand);
  fillZTree(ZCand);
  fillMETTree(ZCand);
  fillJetTree(ZCand);
  fillPhotonTree(ZCand);
  fillEventTree(ZCand);
  ZFrameVariables(ZCand);
  fillMETPrime(ZCand);
  fillProjectedMET();
  
  GetElNNMETValue();

  //MCtruth
  //  if(SampleName == "ZZ_2l2n")
    FillGenEvent();

  Run =  _e.run();
  Event = _e.eventInFile();
  AbsEvent = _e.event();
  strcpy (sample , ((string)SampleName).c_str() );

  int n = (Config::dataPath).size() + SampleName.size() +1;
  TString tmpName = _e.fileName();
  strcpy (fileName, (((string)tmpName).substr(n, ((string)tmpName).size()-n )).c_str() );
  
  
  //Categ of event
  strcpy( EvtCateg, _e.categ().c_str() );


  outTree_.Fill();
}

bool
WWCrossCheck::SpikeCleaning() {


  bool spike[2]={false,false};

  for(int i=0;i<2;i++) {

    Candidate* dau_ = dau[i];

    //Load RecHits
  const CandMap& Sc_Bc_ = _e.candMap("SuperCluster-BasicCluster");
  const CandIdRecHitMap& Bc_rh_ = _e.ecalRecHitMap("BasicCluster-EcalRecHit");

  Candidate* superCluster = _e.getSecond("electron-SuperCluster", dau_ );
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
	    float swisscross = rh.SpikeVar.swCross;
	    
	    if(swisscross > 0.95 )
	      { spike[i]=true; }
	  }
	}
    }

  }

  bool Spike=false;
  if(spike[0] || spike[1])
    Spike=true;

  return Spike;

}




bool
WWCrossCheck::ConversionRejection() {

 //Rejection options
  bool firstPixel[2] ={true,true};
  bool expFirstPixel[2] ={true,true};
  bool EGamConv[2] = {true,true};
  
  bool Conv[2][3]={ {true,true,true},
		 {true,true,true} };


  for(int i=0;i<2;i++) {

    Candidate* dau_ = dau[i];
    if(abs(dau_->pdgCode()) != 11 ) continue;
  Candidate* gsfTrack_ = _e.getSecond( "electron-GsfTrack", dau_);
  CandInfo* info_ = gsfTrack_->info();
  //possibility, use first pixel hit
  if(firstPixel){
    bool FHP_ = info_->getBool("isFPHitTrk");
    if(!FHP_)
      Conv[i][0] = false;
  }

  //other Possibility , using expected hits in pixel layers
  if(expFirstPixel){
    int EIP_ = info_->getInt("expInnHits");
    int maxExpInHit =0;
    if(EIP_ > maxExpInHit)
      Conv[i][1] = false;
  }
  //other possibility, using EGamma rejection, not available with Summer09 datasets
  if(EGamConv) {
    CandInfo* elInfo = dau_->info();
    bool isConv = elInfo->getBool("isConv");
    
    if(isConv)
      Conv[i][2] = false;
  }
  

  }
  inHit1 = Conv[0][0];
  expInHit1 = Conv[0][1];
  ConvRej1 = Conv[0][2];
  inHit2 = Conv[1][0];
  expInHit2 = Conv[1][1];
  ConvRej2 = Conv[1][2];
 
  return (Conv[0][1] && Conv[1][1]);

}





void WWCrossCheck::fillLeptonTree(Candidate* ZCand)
{


  for(int ii=0;ii<2;ii++) 
    {
      
      CandInfo* info = dau[ii]->info();
      
      //ID and Iso ***********
      if((Zmode=="Zemu" && ii==0) || (Zmode=="Zmue" && ii==1)) {
	IDVar[ii][0] = info->getFloat("sigiEtaiEta");
	IDVar[ii][1] = info->getFloat("dEtaIn");
	IDVar[ii][2] = info->getFloat("dPhiIn");
	IDVar[ii][3] = info->getFloat("HOE");
      }
      else {

	bool MuID=true;
	if( !info->getBool("muidGlobalMuonPromptTight") ) MuID=false;
	if(  info->getInt("nTrkHits") < 10 ) MuID=false;
	//	if( fabs(dau[ii]->d0(Primary)) >= 0.2 )     MuID=false;
	if(info->getFloat("normChi2") >= 10 ) MuID=false;
	if(info->getInt("nMatch") <2 )    MuID=false; 
	if(info->getBool("muidTrackerMuonArbitrated")!= true)  MuID=false;
	if(info->getInt("nPixHits") <1 ) MuID=false;
	if(info->getInt("nMuonHits") <1 ) MuID=false;


	IDVar[ii][0] = (float)MuID;
	IDVar[ii][1] = (float)info->getInt("nMuonHits");
	IDVar[ii][2] = (float)info->getInt("nTrkHits");
	IDVar[ii][3] = dau[ii]->d0(Primary);
      }
	//Iso
      if((Zmode=="Zemu" && ii==0) || (Zmode=="Zmue" && ii==1)) {
	IsoVar[ii][0] = info->getFloat("dr03TrkIso");
	IsoVar[ii][1] = info->getFloat("dr03EcalIso");
	IsoVar[ii][2] = info->getFloat("dr03HcalD1Iso") + info->getFloat("dr03HcalD2Iso") ;
      }
      else {
	IsoVar[ii][0] = info->getFloat("isoR03sum"); //trck
	IsoVar[ii][1] = info->getFloat("isoR03strkVeto"); //ecal
	IsoVar[ii][2] = info->getFloat("isoR03semVeto"); //hcal
      }

      //Lepton *********
      Lepton[ii][0] = dau[ii]->pt();
      Lepton[ii][1] = dau[ii]->eta();
      Lepton[ii][2] = dau[ii]->phi();
      Lepton[ii][3] = dau[ii]->p();

      charge[ii] = (int)dau[ii]->charge();
      pdgId[ii] = dau[ii]->pdgCode();
    }


  //Differential variables
  dPhiLepton = KineUtils::dPhi(dau[0]->phi(),dau[1]->phi())*180/Constants::pi;
  dPhiLeptonZ = KineUtils::dPhi(dau[0]->phi(),ZCand->phi())*180/Constants::pi;
  dRLepton = KineUtils::dR(dau[0]->eta(),dau[1]->eta(),dau[0]->phi(),dau[1]->phi());
  dEtaLepton = dau[0]->eta() - dau[1]->eta();
  dPtLepton  = dau[0]->pt() - dau[1]->pt();
  dPLepton   = dau[0]->p() - dau[1]->p();
  
  //Special variables
  //conversion
 

    for(int ii=0;ii<2;ii++) 
      {

	if((Zmode=="Zemu" && ii==0) || (Zmode=="Zmue" && ii==1)) {
	  ConversionRejection();
	  
	  Candidate* superCluster = _e.getSecond("electron-SuperCluster", dau[ii]  );
	  if(superCluster != NULL) {
	    SCLepton[ii][0] = superCluster->Et() ;
	    SCLepton[ii][1] = superCluster->eta();
	    SCLepton[ii][2] = superCluster->phi();
	  }
	  else {
	    SCLepton[ii][0] = -100;
	    SCLepton[ii][1] = -100;
	    SCLepton[ii][2] = -100;
	  }
	  
	  
	}
	else
	  {
	    inHit1    = 1;
	    expInHit1 = 1;
	    ConvRej1  = 1;
	    inHit2    = 1;
	    expInHit2 = 1;
	    ConvRej2  = 1;
	    
	    
	    for(int ii=0;ii<2;ii++)
	      for(int jj=0;jj<3;jj++)
		SCLepton[ii][jj]=-100;
	  }
      }
 
  
  //Lepton multiplicity
  ElecMult =0;
  MuonMult =0;

  const CandList electrons = _e.electrons();
  const CandList muons = _e.muons();
  for(int unsigned iel=0;iel<electrons.size();iel++) {
    Candidate* el = electrons[iel];

    if(el->pt()> 15 && el->uid()!= dau[0]->uid() &&  el->uid()!= dau[1]->uid())
      ElecMult++;
  }
  for(int unsigned imu=0;imu<muons.size();imu++) {
    Candidate* mu = muons[imu];

    if(mu->pt()> 15 && mu->uid()!= dau[0]->uid() &&  mu->uid()!= dau[1]->uid())
      MuonMult++;
  }

}

void WWCrossCheck::fillZTree(Candidate* ZCand)
{
  //Z Variables
  Z[0] = ZCand->pt();
  Z[1] = ZCand->eta();
  Z[2] = ZCand->phi();
  Z[3] = ZCand->p();

  ZMass = ZCand->mass();

  ZMT = sqrt(2*dau[0]->pt()*dau[1]->pt()* (1-cos(KineUtils::dPhi(dau[0]->phi(),dau[1]->phi()))) );
  
  if(CatEta=="BB")
    ZCateg = 0;
  else if(CatEta=="EB")
    ZCateg = 1;
  else if(CatEta=="EE")
    ZCateg = 2;

  //FIXME
  ZMultiplicity = 0;


  //Pf Z
  Candidate* pfdau[2];
  pfdau[0] = _e.getSecond("electron-PfCand", dau[0] );
  pfdau[1] = _e.getSecond("electron-PfCand", dau[1] );

  if(pfdau[0] != NULL && pfdau[1]!=NULL) {
    Candidate* PfZCand = Candidate::create(pfdau[0], pfdau[1]);

    PfZ[0] = PfZCand->pt();
    PfZ[1] = PfZCand->eta();
    PfZ[2] = PfZCand->phi();
    PfZ[3] = PfZCand->p();
    PfZMass = PfZCand->mass();
  }
  else {
    for(int i=0;i<4;i++)
      PfZ[i]=-100;
    PfZMass = -100;
  }


}


void WWCrossCheck::fillMETTree(Candidate* ZCand)
{

  const Candidate* met = _e.met( EventManager::kPatMet , _vmet );
  TVector2 recoil = ComputeRecoil(met,ZCand);
  Double_t rphi = TVector2::Phi_0_2pi(recoil.Phi() );
  Candidate* ZMET = Candidate::create(ZCand,met);

  //MET Variables
  
  MET[0] = met->pt();
  MET[1] = met->phi();

  Recoil[0] = recoil.Mod();
  Recoil[1] = rphi;

  CandInfo* info = met->info();
  SumEt = info->getFloat("sumEt");

  //MET related variables
  _e.refreshRecoMET();
  const Candidate* trkmet = _e.met( EventManager::kPatMet, -1 );
  TrkMET[0] = trkmet->pt();
  TrkMET[1] = trkmet->phi();
  _e.refreshRecoMET();
  
  //Projections
  ComputeCBMs(ZCand);
  
  float MZPara = ZCBM_[0][0]*met->px() + ZCBM_[0][1]*met->py();
  float MZPerp = ZCBM_[1][0]*met->px() + ZCBM_[1][1]*met->py(); 
  float ML1Para = L1CBM_[0][0]*met->px() + L1CBM_[0][1]*met->py(); 
  float ML1Perp = L1CBM_[1][0]*met->px() + L1CBM_[1][1]*met->py(); 
  float ML2Para = L2CBM_[0][0]*met->px() + L2CBM_[0][1]*met->py();  
  float ML2Perp = L2CBM_[1][0]*met->px() + L2CBM_[1][1]*met->py(); 

  float RZPara = ZCBM_[0][0]*recoil.Px() + ZCBM_[0][1]*recoil.Py();
  float RZPerp = ZCBM_[1][0]*recoil.Px() + ZCBM_[1][1]*recoil.Py(); 
  float RL1Para = L1CBM_[0][0]*recoil.Px() + L1CBM_[0][1]*recoil.Py(); 
  float RL1Perp = L1CBM_[1][0]*recoil.Px() + L1CBM_[1][1]*recoil.Py(); 
  float RL2Para = L2CBM_[0][0]*recoil.Px() + L2CBM_[0][1]*recoil.Py();  
  float RL2Perp = L2CBM_[0][0]*recoil.Px() + L2CBM_[0][1]*recoil.Py(); 

  METProjZ[0] = MZPara;
  METProjZ[1] = MZPerp;
  METProjL1[0] = ML1Para;
  METProjL1[1] = ML1Perp;
  METProjL2[0] = ML2Para;
  METProjL2[1] = ML2Perp;

  RecoilProjZ[0] = RZPara;
  RecoilProjZ[1] = RZPerp;
  RecoilProjL1[0] = RL1Para;
  RecoilProjL1[1] = RL1Perp;
  RecoilProjL2[0] = RL2Para;
  RecoilProjL2[1] = RL2Perp;

  //Resolution based variables
  vector<float> MCRval = GetRValues(ZCand->pt(), RecoilProjZ , Nvertex);
  vector<float> DataRval = GetRValuesFromData(ZCand->pt(), RecoilProjZ , Nvertex);

  RValMC[0] = MCRval[0];
  RValMC[1] = MCRval[1];
  RValData[0] = DataRval[0];
  RValData[1] = DataRval[1];

  //Angles
  dPhiMETZ = KineUtils::dPhi(met->phi(), ZCand->phi() )*180/Constants::pi;
  dPhiMETL1 = KineUtils::dPhi(met->phi(), dau[0]->phi() )*180/Constants::pi;
  dPhiMETL2 = KineUtils::dPhi(met->phi(), dau[1]->phi() )*180/Constants::pi;

  dPhiRecoilZ = KineUtils::dPhi(rphi, ZCand->phi() )*180/Constants::pi;
  dPhiRecoilL1 = KineUtils::dPhi(rphi, dau[0]->phi() )*180/Constants::pi;
  dPhiRecoilL2 = KineUtils::dPhi(rphi, dau[1]->phi() )*180/Constants::pi;

  //Kinematics
  
  ZMETMass = ZMET->mass();
  ZMassMET = ZCand->mass() + met->pt();
  ZPtMET = ZCand->pt() + met->pt();
  ZPMET = ZCand->p() + met->pt();
  
  ZRecoil =  ZCand->pt() + recoil.Mod();



  //Special variables
  //Find the closest lepton
  int closest; float dphilep;
  if(dPhiMETL1 < dPhiMETL2)
    {
      closest = 0; dphilep = dPhiMETL1*Constants::pi/180.;
    }
  else 
    {
      closest = 1; dphilep = dPhiMETL2*Constants::pi/180.;
    }
  
  // cout<<" closest "<<closest<<"   "<<dphilep<<endl;
  //Now see if dPhilep < pi/2
  if(fabs(dphilep) < (Constants::pi/2.) )
    {
      SpecMET = met->pt()*sin(dphilep);
      
    }
  else
    SpecMET = met->pt();

  // cout<<" specMET "<<SpecMET<<"    "<<met->pt()<<endl;


}


void WWCrossCheck::fillJetTree(Candidate* ZCand)
{

  //Access to jet and MET
  const Candidate* met = _e.met( EventManager::kPatMet, _vmet );
  JetMultiplicity=0;

  TVector2 Ht(met->px(),met->py());
  TVector2 Jet(0,0);

  int ThrJetMult[10]={15,20,25,30,35,40,50,60,80,100};
  for(int i=0;i<10;i++) {
    MultiJetMult[i] =0;
  }

  Candidate* theJet;
  float pttmp=0;
  const CandList& Jets = _e.jetList( EventManager::kPatJet);
  for(int unsigned ij=0;ij<Jets.size();ij++) { 
    Candidate* jet = Jets[ij];
    {
      float dr1_ = KineUtils::dR(dau[0]->eta(), jet->eta(),
				 dau[0]->phi(), jet->phi() );
      float dr2_ = KineUtils::dR(dau[1]->eta(), jet->eta(),
				 dau[1]->phi(), jet->phi() );

      /*float drE_ = KineUtils::dR(LeadAddElec[1], jet->eta(),
				 LeadAddElec[2], jet->phi() );
      float drM_ = KineUtils::dR(LeadAddMuon[1], jet->eta(),
				 LeadAddMuon[2], jet->phi() );
      */
      if(jet->pt()> 5 && dr1_>0.2 && dr2_>0.2 /*&& drE_>0.2 && drM_>0.2*/ )
	{

	  Jet.SetMagPhi(jet->pt(), jet->phi() );
	  Ht += Jet;

	  if(jet->pt()> 10) {
	    if(pttmp<jet->pt() )
	      { theJet=jet; pttmp = jet->pt(); }
	  }
	  if(jet->pt()> 25) {
	    JetMultiplicity++;
	  }
	  for(int i=0;i<10;i++) {
	    if( ThrJetMult[i] < jet->pt() ){
	      MultiJetMult[i]++;
	    }
	  }
	  
	}
    }
  }

  //Leading Jet
 

  if(Jets.size()!=0 && pttmp!=0) {

    CandInfo* info = theJet->info();
    LJet[0] = theJet->pt();
    LJet[1] = theJet->eta();
    LJet[2] = theJet->phi();
    LJet[3] = 0;//FIXME//info->getFloat("chargedHadronEnergyFraction") + info->getFloat("neutralHadronEnergyFraction");

    Candidate* ZJet = Candidate::create(ZCand,theJet);

  dPhiJetMET = KineUtils::dPhi(theJet->phi(), met->phi() )*180/Constants::pi;
  dPhiJetZ = KineUtils::dPhi(theJet->phi(), ZCand->phi())*180/Constants::pi;
  dEtaJetZ = ZCand->eta() - theJet->eta();
  ZPtJet = ZCand->pt() + theJet->pt();
  ZJetMass = ZJet->mass();
  dPhiJetL1 = KineUtils::dPhi(theJet->phi(), dau[0]->phi())*180/Constants::pi;
  dPhiJetL2 = KineUtils::dPhi(theJet->phi(), dau[1]->phi())*180/Constants::pi;
  dEtaJetL1 = dau[0]->eta() - theJet->eta();
  dEtaJetL2 = dau[1]->eta() - theJet->eta();

  
  bTag[0] = info->getFloat("btagTkC");
  bTag[1] = info->getFloat("btagSoftM");

    }
  else {
    LJet[0] = -1000;
    LJet[1] = -1000;
    LJet[2] = -1000;
    LJet[3] = -1000;
    dPhiJetMET = -1000;
    dPhiJetZ = -1000;
    dEtaJetZ = -1000;
    ZPtJet = -1000;
    ZJetMass = -1000;
    dPhiJetL1 = -1000;
    dPhiJetL2 = -1000;
    dEtaJetL1 = -1000;
    dEtaJetL2 = -1000;

    bTag[0] = -1000;
    bTag[1] = -1000;

  }

  //HT
  HT[0] = Ht.Mod();
  HT[1] = Ht.Phi();


  //trackJet
  const CandList& TcJets = _e.jetList( EventManager::kPatJet);
  for(int unsigned ij=0;ij<TcJets.size();ij++) { 
    Candidate* tcjet = TcJets[ij];
    
    float dr1_ = KineUtils::dR(dau[0]->eta(), tcjet->eta(),
			       dau[0]->phi(), tcjet->phi() );
    float dr2_ = KineUtils::dR(dau[1]->eta(), tcjet->eta(),
			       dau[1]->phi(), tcjet->phi() );
    if(tcjet->pt()> 15 && dr1_>0.2 && dr2_>0.2  )
      {
	TcJetMultiplicity++;
      }
  }
  




}


void WWCrossCheck::fillPhotonTree(Candidate* ZCand)
{

  //Access to photon and MET
  const Candidate* met = _e.met( EventManager::kPatMet, _vmet );
  PhotonMultiplicity=0;
  Candidate* thePhoton;
  float pttmp=0;
  const CandList& Photons = _e.photons();
  for(int unsigned ij=0;ij<Photons.size();ij++) { 
    Candidate* photon = Photons[ij];
    {
      float dr1_ = KineUtils::dR(dau[0]->eta(), photon->eta(),
				dau[0]->phi(), photon->phi() );
      float dr2_ = KineUtils::dR(dau[1]->eta(), photon->eta(),
				 dau[1]->phi(), photon->phi() );
      
      if(photon->pt()> 15 && dr1_>0.2 && dr2_>0.2  )
	{
	  PhotonMultiplicity++;
	  if(pttmp<photon->pt() )
	    { thePhoton=photon; pttmp = photon->pt(); }
	}
    }
  }
  

  if(Photons.size()!=0 && pttmp!=0) {

    //Leading Photon
    LPhoton[0] = thePhoton->pt();
    LPhoton[1] = thePhoton->eta();
    LPhoton[2] = thePhoton->phi();
    
    CandInfo* info = thePhoton->info();
    float RelIso = (info->getFloat("ecalRecHitSumEtConeDR03") 
		    + info->getFloat("hcalDepth1TowerSumEtConeDR03")
		    + info->getFloat("hcalDepth2TowerSumEtConeDR03")
		    + info->getFloat("trkSumPtHollowConeDR03") )/ thePhoton->pt();

    LPhoton[3] = info->getFloat("sigmaIetaIeta");
    LPhoton[4] = RelIso;
    

    Candidate* ZPhoton = Candidate::create(ZCand,thePhoton);
    
  dPhiPhotonMET = KineUtils::dPhi(thePhoton->phi(), met->phi() )*180/Constants::pi;
  dPhiPhotonZ = KineUtils::dPhi(thePhoton->phi(), ZCand->phi())*180/Constants::pi;
  dEtaPhotonZ = ZCand->eta() - thePhoton->eta();
  ZPtPhoton = ZCand->pt() + thePhoton->pt();
  ZPhotonMass = ZPhoton->mass();
  dPhiPhotonL1 = KineUtils::dPhi(thePhoton->phi(), dau[0]->phi())*180/Constants::pi;
  dPhiPhotonL2 = KineUtils::dPhi(thePhoton->phi(), dau[1]->phi())*180/Constants::pi;
  dEtaPhotonL1 = dau[0]->eta() - thePhoton->eta();
  dEtaPhotonL2 = dau[1]->eta() - thePhoton->eta();

    }
  else {
    dPhiPhotonMET = 0;
    dPhiPhotonZ = 0;
    dEtaPhotonZ = 0;
    ZPtPhoton = 0;
    ZPhotonMass = 0;
    dPhiPhotonL1 = 0;
    dPhiPhotonL2 = 0;
    dEtaPhotonL1 = 0;
    dEtaPhotonL2 = 0;

    for(int kk=0;kk<5;kk++)
      LPhoton[kk] = -1000;

  }


  
}

void WWCrossCheck::fillEventTree(Candidate* ZCand)
{
  //  const VtxList& vertices = _e.vertices();
  //  Nvertex = (int)vertices.size();
}



void WWCrossCheck::ComputeCBMs(Candidate* ZCand) 
{
 //Compute ZCBM and LCBMs
  ZCBM_[0][0]=ZCand->cosPhi();
  ZCBM_[0][1]=ZCand->sinPhi();
  ZCBM_[1][0]=-ZCand->sinPhi();
  ZCBM_[1][1]=ZCand->cosPhi();

  L1CBM_[0][0]=dau[0]->cosPhi();
  L1CBM_[0][1]=dau[0]->sinPhi();
  L1CBM_[1][0]=-dau[0]->sinPhi();
  L1CBM_[1][1]=dau[0]->cosPhi();

  L2CBM_[0][0]=dau[1]->cosPhi();
  L2CBM_[0][1]=dau[1]->sinPhi();
  L2CBM_[1][0]=-dau[1]->sinPhi();
  L2CBM_[1][1]=dau[1]->cosPhi();
  //======
}


void WWCrossCheck::ZFrameVariables(Candidate* ZCand) {

  //  TLorentzVector p4Z = ZCand->p4();
  TLorentzVector p4l1 = dau[0]->p4();
  TLorentzVector p4l2 = dau[1]->p4();


  const Candidate* met = _e.met( EventManager::kPatMet, _vmet  );
  Candidate* ZMET = Candidate::create(ZCand,met);
  TLorentzVector p4Z = ZMET->p4();
  
  //Compute the Z boost
  TVector3 boost = p4Z.BoostVector();
  //remove the z component
  boost.SetZ(0);
  //and create daughters in the Z rest frame
  p4l1.Boost(-boost);
  p4l2.Boost(-boost);
  Candidate* daustar1 = Candidate::create(p4l1.Vect(), dau[0]->charge(), dau[0]->mass());
  Candidate* daustar2 = Candidate::create(p4l2.Vect(), dau[1]->charge(), dau[1]->mass());
  
  if(daustar1->charge()<0) {
    CosThetaM = daustar1->cosTheta();
    CosThetaP = daustar2->cosTheta();
  }
  else {
    CosThetaP = daustar1->cosTheta();
    CosThetaM = daustar2->cosTheta();
  }

  dPhiLeptonS = KineUtils::dPhi(daustar1->phi(),daustar2->phi())*180/Constants::pi;
  dRLeptonS = KineUtils::dR(daustar1->eta(),daustar2->eta(),daustar1->phi(),daustar2->phi());
  dEtaLeptonS = daustar1->eta() - daustar2->eta();
  dPtLeptonS  = daustar1->pt() - daustar2->pt();
  dPLeptonS   = daustar1->p() - daustar2->p();


  //Collins Soper frame

  //Decompose the boost
  /*  TVector3 longB;
  TVector3 transB;

  TVector3 boostZ = ZCand->p4().BoostVector();

  transB.SetXYZ( boostZ.Px(),boostZ.Py(),0 );
  longB.SetXYZ(0,0,boostZ.Pz() );

  TLorentzVector p4l1cs = dau[0]->p4(); 
  TLorentzVector p4l2cs = dau[1]->p4(); 
  
  cout<<"New event "<<endl;


  cout<<" Boost "<<boostZ.Pz()<<endl;
  cout<<" Boost2 "<<boostZ.Px()<<"  "<<boostZ.Py()<<"    "<<boostZ.Pt()<<endl;

  cout<<" L1 "<<p4l1cs.Px()<<"   "<<p4l1cs.Py()<<"   "<<p4l1cs.Pz()<<"   "<<p4l1cs.Pt()<<"   "<<endl;
  cout<<" L2 "<<p4l2cs.Px()<<"   "<<p4l2cs.Py()<<"   "<<p4l2cs.Pz()<<"   "<<p4l2cs.Pt()<<"   "<<endl;

  //First boost
  p4l1cs.Boost(-longB);
  p4l2cs.Boost(-longB);

  cout<<" First boost "<<endl;
  cout<<" L1 "<<p4l1cs.Px()<<"   "<<p4l1cs.Py()<<"   "<<p4l1cs.Pz()<<"   "<<p4l1cs.Pt()<<"   "<<endl;
  cout<<" L2 "<<p4l2cs.Px()<<"   "<<p4l2cs.Py()<<"   "<<p4l2cs.Pz()<<"   "<<p4l2cs.Pt()<<"   "<<endl;

  //Second boost
  p4l1cs.Boost(-transB);
  p4l2cs.Boost(-transB);

  cout<<" Second boost "<<endl;
  cout<<" L1 "<<p4l1cs.Px()<<"   "<<p4l1cs.Py()<<"   "<<p4l1cs.Pz()<<"   "<<p4l1cs.Pt()<<"   "<<endl;
  cout<<" L2 "<<p4l2cs.Px()<<"   "<<p4l2cs.Py()<<"   "<<p4l2cs.Pz()<<"   "<<p4l2cs.Pt()<<"   "<<endl<<endl;


  Candidate* daucs1 = Candidate::create(p4l1cs.Vect(), dau[0]->charge(), dau[0]->mass());
  Candidate* daucs2 = Candidate::create(p4l2cs.Vect(), dau[1]->charge(), dau[1]->mass());
  
  if(daustar1->charge()<0) {
      CosThetaCSM = daucs1->cosTheta();
      CosThetaCSP = daucs2->cosTheta();
    }
    else {
      CosThetaCSP = daucs1->cosTheta();
      CosThetaCSM = daucs2->cosTheta();
    }
  */

  float pp, pm, Ep, Em;

  if(dau[0]->charge() == 1) {
    pp = dau[0]->pz() ; Ep = dau[0]->E();
    pm = dau[1]->pz() ; Em = dau[1]->E();
  }
  else {
    pp = dau[1]->pz() ; Ep = dau[1]->E();
    pm = dau[0]->pz() ; Em = dau[0]->E(); 
  }


  CosThetaCSP = 2*(pm*Ep-pp*Em)/(ZCand->mass() * sqrt(pow(ZCand->mass(),2) + pow(ZCand->pt(),2) ) );
    
  CosThetaCSM = CosThetaCSP;

}


void WWCrossCheck::fillMETPrime(Candidate* ZCand) {

  //From d0 analysis, against mismeasurement


  //First, define the frame
  TVector2 TruthFrame(0,0);

  TruthFrame = dau[0]->p2() - dau[1]->p2();
  
  TVector2 tl = TruthFrame.Unit();
  TVector2 Zpt = ZCand->p2();
  
  TVector2 tt = tl.Rotate(Constants::pi/2);
  tt *= cos(tt.DeltaPhi(dau[1]->p2()))/fabs( cos(tt.DeltaPhi(dau[1]->p2())) ); 
  if( fabs(dPhiLepton) < Constants::pi/2.) {
    tt = Zpt.Unit();
    tl = tt.Rotate( -Constants::pi/2);
    tl *= -cos(tl.DeltaPhi(dau[1]->p2()))/fabs( cos(tl.DeltaPhi(dau[1]->p2())) );
  }

  //Now determine at and al for Z system
  float at,al;

  at =  Zpt*tt;
  al =  Zpt*tl;


  //Compute the maximum uncertainty on dilepton ET (2 sigma deviation)
  
  float sigma1=0; 
  float sigma2=0;
  if(Zmode=="Zee") {
    sigma1 = ComputeResolution(dau[0]->E() );
    sigma2 = ComputeResolution(dau[1]->E() );
  }
  else {
    sigma1=0;
    sigma2=0;
  }
  float dat,dal;

  TVector2 dTFtmp = (1-sigma1)*dau[0]->p2() - (1-sigma2)*dau[1]->p2();
  TVector2 dZpt = (1-sigma1)*dau[0]->p2() + (1-sigma2)*dau[1]->p2();
  TVector2 dtt = (dTFtmp.Unit()).Rotate(Constants::pi/2.);
  dtt *= cos(dtt.DeltaPhi(dau[1]->p2()))/fabs( cos(dtt.DeltaPhi(dau[1]->p2())) ); 

  dat = dZpt*dtt - at;//(dZpt.Proj(dtt)).Mod() * cos(dZpt.DeltaPhi(dtt))/fabs(cos(dZpt.DeltaPhi(dtt)))  - at;
  
  TVector2 dTFtmp2 = -sigma1*dau[0]->p2() + sigma2*dau[1]->p2();
  
  dal = dTFtmp2*tl;//(dTFtmp2.Proj(tl)).Mod() * cos(dTFtmp2.DeltaPhi(tl))/fabs(cos(dTFtmp2.DeltaPhi(tl)));

  //===============================================

  //Underlying event effects
  //to be sure we take the good one
  int v=-1;
  for(int unsigned i=0;i<_e.a().recoMET.size();i++) {
    if("metJESCorPFAK5" ==_e.a().recoMET[i])
      { v=i; break; }
  }
  const Candidate* met = _e.met( EventManager::kPatMet, _vmet );

  TVector2 recoil = ComputeRecoil(met,ZCand);

  float rt,rl;
  
  rt = recoil*tt;//(recoil.Proj(tt)).Mod() * cos(recoil.DeltaPhi(tt))/fabs(cos(recoil.DeltaPhi(tt)));
  rl = recoil*tl;//(recoil.Proj(tl)).Mod() * cos(recoil.DeltaPhi(tl))/fabs(cos(recoil.DeltaPhi(tl)));

  METPrime[0] = at ;
  METPrime[1] = al ;

  if(rt>0) rt =0;
  if(rl>0) rl =0;

  METPrimeRecoil[0] = rt;
  METPrimeRecoil[1] = rl;

  METPrimeCor[0] = dat;
  METPrimeCor[1] = dal;

  
  float Stjetm=0,Sljetm=0;
  float Stjetp=0,Sljetp=0;
  
  const CandList& Jets = _e.jetList( EventManager::kPatJet);
  for(int unsigned ij=0;ij<Jets.size();ij++) { 
    Candidate* jet = Jets[ij];
    
      float dr1_ = KineUtils::dR(dau[0]->eta(), jet->eta(),
				 dau[0]->phi(), jet->phi() );
      float dr2_ = KineUtils::dR(dau[1]->eta(), jet->eta(),
				 dau[1]->phi(), jet->phi() );
      
      float atjet = jet->p2() * tt;
      float aljet = jet->p2() * tl;

      // cout<<jet->pt()<<"   "<<dr1_<<"  "<<dr2_<<endl;

      if( dr1_>0.2 && dr2_>0.2 && jet->pt() > 10 )
	{
	  if(atjet < 0) {
	    Stjetm += atjet; }
	  else {
	    Stjetp += atjet;
	  }
	  
	  if(aljet<0) {
	    Sljetm += aljet; }
	  else {
	    Sljetp += aljet;
	  }

	}
  }

    float metal=0,metat=0;
    metat -= (met->p2()) * tt;
    metal -= (met->p2()) * tl;
    
   
    METPrimeMETM[0] = metat;
    METPrimeMETM[1] = metal;

    METPrimeJetM[0] =Stjetm;
    METPrimeJetM[1] =Sljetm;
    METPrimeJetP[0] =Stjetp;
    METPrimeJetP[1] =Sljetp;	
   
    
    //Compute correction with pfcand ===============================
    float Stcalom=0,Slcalom=0;
    float Stcalop=0,Slcalop=0;

    const CandList& pfcands_ = _e.pfCandidates();
    for(int unsigned i=0;i<pfcands_.size();i++) {
      Candidate* pfCand = pfcands_[i];
      float dR;
      bool noMatch=true;
      
      for(int unsigned ij=0;ij<Jets.size();ij++) { 
	Candidate* jet = Jets[ij];
	
	if(jet->pt()<10) continue;

	dR = KineUtils::dR(jet->eta(),pfCand->eta(),jet->phi(),pfCand->phi() );
	if(dR<0.51 && noMatch) {noMatch=false; continue;}
      }

      if(noMatch) {
	
	float dr1_ = KineUtils::dR(dau[0]->eta(), pfCand->eta(),
				   dau[0]->phi(), pfCand->phi() );
	float dr2_ = KineUtils::dR(dau[1]->eta(), pfCand->eta(),
				   dau[1]->phi(), pfCand->phi() );
	
	float atcalo = pfCand->p2() * tt;
	float alcalo = pfCand->p2() * tl;
	
	if( dr1_>0.1 && dr2_>0.1  )
	  {
	    if(atcalo < 0) {
	      Stcalom += atcalo; }
	    else {
	      Stcalop += atcalo;
	    }

	    if(alcalo<0) {
	      Slcalom += alcalo; }
	    else {
	      Slcalop += atcalo;
	    }
	  }
      }
    }
    
    _e.refreshPfCandidates(); //to be sure we take the good one
    int vpfc=1; //pfCand_lowPt
    _e.reloadPfCandidates( vpfc );
    const CandList& pfCandsLowPt_ = _e.candList( EventManager::kPfCand_lowPt );
    for(int unsigned i=0;i<pfCandsLowPt_.size();i++) {
      Candidate* pfCand = pfCandsLowPt_[i];
      float dR;
      bool noMatch=true;
      
      for(int unsigned ij=0;ij<Jets.size();ij++) { 
	Candidate* jet = Jets[ij];
	
	if(jet->pt()<10) continue;

	dR = KineUtils::dR(jet->eta(),pfCand->eta(),jet->phi(),pfCand->phi() );
	if(dR<0.51 && noMatch) {noMatch=false; continue;}
      }

      if(noMatch) {
	
	float dr1_ = KineUtils::dR(dau[0]->eta(), pfCand->eta(),
				   dau[0]->phi(), pfCand->phi() );
	float dr2_ = KineUtils::dR(dau[1]->eta(), pfCand->eta(),
				   dau[1]->phi(), pfCand->phi() );
	
	float atcalo = pfCand->p2() * tt;
	float alcalo = pfCand->p2() * tl;
	
	if( dr1_>0.1 && dr2_>0.1  )
	  {
	    if(atcalo < 0) {
	      Stcalom += atcalo; }
	    else {
	      Stcalop += atcalo;
	    }

	    if(alcalo<0) {
	      Slcalom += alcalo; }
	    else {
	      Slcalop += atcalo;
	    }
	  }
      }
    }

    METPrimeUnclusM[0] = Stcalom;
    METPrimeUnclusM[1] = Slcalom;
    METPrimeUnclusP[0] = Stcalop;
    METPrimeUnclusP[1] = Slcalop;
    
    //End correction computing with pfcand ===============================

    

    double Tt = at + rt + Stcalom + dat;
    double Tl = al + rl + Slcalom + dal;

    Tt = max(Tt,0.);
    Tl = max(Tl,0.);
  
}


void WWCrossCheck::fillProjectedMET() {

    const Candidate* met = _e.met( EventManager::kPatMet, _vmet );

  float dphi1 = fabs( KineUtils::dPhi(dau[0]->phi(), met->phi() ) );
  float dphi2 = fabs( KineUtils::dPhi(dau[1]->phi(), met->phi() ) );

  

  float dphimin = min(dphi1, dphi2);

  if(dphimin > (Constants::pi /2) )
    ProjMET = met->pt();
  else
    ProjMET = met->pt()*sin(dphimin);


  //Compute corrected ProjMET (Pileup)
  CorProjMET = (ProjMET -1.099*(Nvertex-1))/( 0.896 + 0.1056*Nvertex );
  

}



void WWCrossCheck::FillGenEvent() {


  MCToGen.clear();
  GenToMC.clear();
  MCok=true;
  
  MCTruthChain();
  if(!MCok) 
    MCTruthChain();
  
  float lepton[4][5];
  for(int i=0;i<4;i++)
    for(int l=0;l<5;l++)
      lepton[i][l]=-100;

  int kl=0;

  int nZn=-1, lZn=-1;

  int npdg[3]={12,14,16};
  int lpdg[3]={11,13,15};

  const CandList& mcCands = _e.mcTruthCandidates();

  for(int unsigned imc=0;imc<mcCands.size();imc++) {

    Candidate* mcCand = mcCands[imc];
    if(abs(mcCand->pdgCode()) >=  11  && abs(mcCand->pdgCode()) <= 18  ) {
      int pdgid = OriginOfMCCandidate(mcCand);
    
      if(abs(pdgid) == 23 )
	{
	  lepton[kl][0] = mcCand->pt();
	  lepton[kl][1] = mcCand->eta();
	  lepton[kl][2] = mcCand->phi();
	  lepton[kl][3] = mcCand->E();
	  lepton[kl][4] = mcCand->pdgCode();
	  kl++;

	  if(abs(pdgid)==23 && (mcCand->pdgCode()==npdg[0] ||
				mcCand->pdgCode()==npdg[1] ||
				mcCand->pdgCode()==npdg[2] )  )
	    nZn = imc;
	  if(abs(pdgid)==23 && (mcCand->pdgCode()==lpdg[0] ||
				mcCand->pdgCode()==lpdg[1] ||
				mcCand->pdgCode()==lpdg[2] )  )
	    lZn = imc;


	}
    }
  }


  int nn=0;
  int ln=0;

  for(int il=0;il<4;il++) {
    
    for(int ii=0;ii<3;ii++)
      if(fabs(lepton[il][4])==npdg[ii])
	{
	  for(int k=0;k<4;k++)
	    McNeutrino[nn][k] =lepton[il][k];
	  nn++;
	  break;
	}

     for(int ii=0;ii<3;ii++)
       if(fabs(lepton[il][4])==lpdg[ii])
	{
	  /* cout<<"  ";
	  for(int i=0;i<5;i++)
	    if( lepton[il][i]> 10000)
	      cout<< lepton[il][i]<<"  "<<i<<"   "<<il;
	      cout<<endl;*/

	  for(int k=0;k<4;k++)
	    McLepton[ln][k] =lepton[il][k];
	  ln++;
	  break;
	}
  }

  
  Candidate* Zneutrino = NULL;
  if(nZn != -1)
    if(mcCands[nZn] != NULL)
    Zneutrino = MotherOfMCCandidate(mcCands[nZn]);
  Candidate* Zlepton = NULL;
  if(lZn != -1)
    if(mcCands[lZn] != NULL)
   Zlepton = MotherOfMCCandidate(mcCands[lZn]);

  if( Zneutrino!=NULL) {

    McZn[0]=Zneutrino->pt();
    McZn[1]=Zneutrino->eta();
    McZn[2]=Zneutrino->phi();
    McZn[3]=Zneutrino->mass();
  }
  else {
    
    McZn[0]= -100;
    McZn[1]= -100;
    McZn[2]= -100;
    McZn[3]= -100;
  }

  if( Zlepton!=NULL) {
    McZl[0]=Zlepton->pt();
    McZl[1]=Zlepton->eta();
    McZl[2]=Zlepton->phi();
    McZl[3]=Zlepton->mass();
  }
  else {
    McZl[0]= -100;
    McZl[1]= -100;
    McZl[2]= -100;
    McZl[3]= -100;
    
  }



}


bool WWCrossCheck::MCTruthChain()
{

  PrintTree printTree;

  Candidate* fCand= _e.decayTree();
  if(fCand==NULL) return false;
  const CandList& MCList = _e.mcTruthCandidates();
  int pdgIds[11]={11,12,13,
		  14,15,16,
		  17,18,
		  22,23,24}; //Matching of all leptons and bosons
  for(int il=0;il<11;il++) {
   
    CandList CandToMatch;
      
    CandUtil::get( pdgIds[il] , fCand , CandToMatch );
    if(pdgIds[il] != 22 && pdgIds[il]!=23)
      CandUtil::get( -pdgIds[il] , fCand , CandToMatch );
     
    for(unsigned int ic=0;ic<CandToMatch.size();ic++) {
      Candidate* genCand = CandToMatch[ic];
      int pdgId = genCand->pdgCode();
      std::map< Candidate* , pair<float,float> > MCDrAssocMap;
      
      if( abs(pdgId) == pdgIds[il] && genCand->theMother()->pdgCode() != pdgId &&  genCand->p()>0.00001 ) { 
	
	if(!MCok)
	  cout<<" Particle "<<pdgId<<" from "<<genCand->theMother()->pdgCode()<<" ;  "<<genCand->pt()<<"  "<<genCand->p()<<endl;
	for(unsigned int imc=0;imc<MCList.size(); imc++) {
	  Candidate* MCCand = MCList[imc];
	 
	  if( MCCand->pdgCode() == pdgId ) {

	    float dR=KineUtils::dR(MCCand->eta(), genCand->eta(), MCCand->phi(), genCand->phi() );
	    float dPt = fabs( fabs(MCCand->pt()) - fabs(genCand->pt()) )/genCand->pt();
	    
	    MCDrAssocMap[ MCCand ].first= dR;
	    MCDrAssocMap[ MCCand ].second= dPt;
	 
	  }
	} //MC List
	
	if( MCDrAssocMap.size() != 0 ) {
	  Candidate* BestMC = FindBestdRCandidate(MCDrAssocMap);
	  if(BestMC==NULL)
	    /* cout<<" BestMC "<<BestMC<<"  "<<BestMC->pdgCode()<<endl;
	       else*/
	    cout<< " truc nul "<<_ievt<<"   "<<BestMC<<endl;
	  if(MCok) {
	    MCToGen[ BestMC ] = genCand;
	    GenToMC[ genCand ] = BestMC;
	  }
	}
	else
	  cout<<" Pas d'association avec  "<<pdgId<<" , "<<genCand->pt()<<"  "<<genCand->theMother()->pdgCode()<<endl;
	
      }// End cut on pdgID and pt for matching interesting particles
    }//Gen List
  }//Particles
  
  /*  MatchingIterator it;
  MatchingIterator itbegin = MCToGen.begin();
  MatchingIterator itend = MCToGen.end();
  cout<<" MCTogen size  "<< MCToGen.size()<<endl;
  for(it = itbegin; it!=itend; it++) {
   
    // int pdgidM = (*it).second->theMother()->pdgCode();
    int pdgid = (*it).first->pdgCode();
    
    if( abs(pdgid)==11 )
      fill("pdgId", "DeltaPdgId_electron",  (float)(*it).first->pdgCode() -  (float)(*it).second->pdgCode() );
    if( abs(pdgid)==23) {
      fill("mass","ZMCmass",(*it).first->mass() );
      Candidate * Zcand = (*it).second;
      

      int nd=Zcand->nDaughters();
      for(int i=0;i<nd;i++) {
	Candidate* dau= Zcand->daughter(i);
	MatchingIterator test = GenToMC.find( dau);
	if(test != GenToMC.end() ) {
	  // Candidate* mcdau = GenToMC.find( dau)->second;
	  nE++;
	}
	  else  {
	  if(MCok)
	    Nent++;
	  MCok=false;
	  cout<<" Lepton non trouve "<< Nent<<endl;
	}
      }
      
    }
  }*/			     
   return true;
}


Candidate* WWCrossCheck::FindBestdRCandidate(std::map<Candidate*, std::pair<float,float> > AssocMap ) {

  std::map<Candidate*, pair<float,float> >::const_iterator it;
  std::map<Candidate*, pair<float,float> >::const_iterator itbegin= AssocMap.begin();
  std::map<Candidate*, pair<float,float> >::const_iterator itend= AssocMap.end();
  
  float dRtmp=10;
  //float dPttmp;
  Candidate* BestCand;
  for(it = itbegin;it != itend; it++) {
    
    float dR = (*it).second.first;
    // float dPt =  (*it).second.second;
    if(dR < dRtmp) {
      dRtmp = dR;
      BestCand = (*it).first;
    }
  }
  
  return BestCand;
}



int WWCrossCheck::OriginOfMCCandidate(Candidate* cand) {

  MatchingIterator it = MCToGen.find( cand);
  if(it != MCToGen.end() )
    return (*it).second->theMother()->pdgCode();
  else
    return -100;
  
}




Candidate* WWCrossCheck::MotherOfMCCandidate(Candidate* cand) {

  MatchingIterator it = MCToGen.find( cand);
  if(it != MCToGen.end() )
    return (*it).second->theMother();
  else
    return NULL;
  
}




float WWCrossCheck::ComputeResolution(float ElE) {

  return sqrt(((pow(0.028/sqrt(ElE),2) + pow(0.12/ElE,2) + 0.000009 )));



}


vector<float> WWCrossCheck::GetRValues(float ZPt, float Recoil[2], int Nv) {

  vector<float> values(2,0);

  float Rcor = 1./GetResponseForPFT1(ZPt);
  
  float stotpara = sqrt( pow( sqrt(ZPt)*1.26114 + 0  ,2 )   + pow(4.9167 ,2 )*pow( Rcor, 2 ) + ( pow(3.51947,2)*Nv) *pow( Rcor, 2 ) );
  float stotperp = sqrt( pow( sqrt(ZPt)*0.665988 + 1.11022e-15  ,2 )   + pow(4.89571 ,2 )*pow( Rcor, 2 ) + ( pow(3.50427,2)*Nv) *pow( Rcor, 2 ) );
  
  values[0] = (Recoil[0]*Rcor-ZPt)/stotpara;
  values[1] = (Recoil[1]*Rcor)/stotperp;

  return values;
}


vector<float> WWCrossCheck::GetRValuesFromData(float ZPt, float Recoil[2], int Nv) {

  vector<float> values(2,0);

  float Rcor = 1./GetResponseForPFT1(ZPt);
  
  float stotpara = sqrt( pow( sqrt(ZPt)*0.892115 + 1.94915  ,2 )   + pow(4.99392 ,2 )*pow( Rcor, 2 ) + ( pow(3.40582,2)*Nv) *pow( Rcor, 2 ) );
  float stotperp = sqrt( pow( sqrt(ZPt)*0.612103 + 1.5691  ,2 )   + pow(4.77942 ,2 )*pow( Rcor, 2 ) + ( pow(3.66613,2)*Nv) *pow( Rcor, 2 ) );

  values[0] = (Recoil[0]*Rcor-ZPt)/stotpara;
  values[1] = (Recoil[1]*Rcor)/stotperp;

  return values;
}


float WWCrossCheck::GetResponseForPFT1(float ZPt) {
  if(ZPt>3)
    return 1.019 -36.8/(113.98 +  ZPt*ZPt);
  else 
    return 0.0146*ZPt + 0.676; 
}



bool WWCrossCheck::fillVertices() {

  const VtxList& vertices = _e.vertices();
 
  Nvertex = 0;

 
  bool findPV=false;
  //DZ vertices
  
  for(int unsigned iv=0;iv<vertices.size();iv++) {
    
    Vertex* Vert = vertices[iv];
    CandInfo* infoV = Vert->info();
   
    bool goodV=false;
    if(fabs((Vert->pos()).Z())<24 && (Vert->pos()).Perp()<2 && infoV->getFloat("ndof")> 4 && !(infoV->getBool("isFake")) )
      goodV=true;
  
    if(!findPV && goodV )
      {Primary = Vert; Nvertex++; findPV=true; continue;}
      
    if(goodV && findPV && vertices.size()>1) {
      Nvertex++;
    }
  }
 
  if(findPV)
    return true;
  else
    return false;

}




void WWCrossCheck::GetElNNMETValue() {

  vector<double> inputVecMET(9,0);

  inputVecMET[0] = (double)METPrime[0];
  inputVecMET[1] = (double)METPrime[1];
  //  inputVecMET[2] = (double)METPrimeCor[0];
  //  inputVecMET[3] = (double)METPrimeCor[1];
  inputVecMET[2] = (double)METPrimeMETM[0];
  inputVecMET[3] = (double)METPrimeMETM[1];
  inputVecMET[4] = (double)METPrimeJetM[0]/Nvertex;
  inputVecMET[5] = (double)METPrimeJetM[1]/Nvertex;
  // inputVecMET[8] = (double)METPrimeJetP[0]/Nvertex;
  //  inputVecMET[9] = (double)METPrimeJetP[1]/Nvertex;
  inputVecMET[6] = (double)METPrimeUnclusM[0]/Nvertex;
  inputVecMET[7] = (double)METPrimeUnclusM[1]/Nvertex;
  inputVecMET[8] = (double)METPrimeUnclusP[0]/Nvertex;
  //  inputVecMET[13] = (double)METPrimeUnclusP[1]/Nvertex;
  outMETNN = (float)(metNN_->GetMvaValue( inputVecMET ));


}


void WWCrossCheck::InitializeElNNMET() {


 //Définition des variables (comme pour le training)
  vector<string> inputVars;

  inputVars.push_back("METPrime[0]");
  inputVars.push_back("METPrime[1]");
  // inputVars.push_back("METPrimeCor[0]");
  //  inputVars.push_back("METPrimeCor[1]");
  inputVars.push_back("METPrimeMETM[0]");
  inputVars.push_back("METPrimeMETM[1]");
  inputVars.push_back("METPrimeJetM[0]/NVertex");
  inputVars.push_back("METPrimeJetM[1]/NVertex");
  //  inputVars.push_back("METPrimeJetP[0]/NVertex");
  //  inputVars.push_back("METPrimeJetP[1]/NVertex");
  inputVars.push_back("METPrimeUnclusM[0]/NVertex");
  inputVars.push_back("METPrimeUnclusM[1]/NVertex");
  inputVars.push_back("METPrimeUnclusP[0]/NVertex");
  // inputVars.push_back("METPrimeUnclusP[1]/NVertex");
  // on crée un neural net du type METNN
  metNN_ = new METNN( inputVars );
  
}

