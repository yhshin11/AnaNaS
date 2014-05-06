#include "Analysis/src/OfficialZAnalysis.hh"

#include <algorithm>
#include <cassert>
#include <sstream>

using namespace std;

int OfficialZAnalysis::imode = EventServer::kZmm;

ClassImp( OfficialZAnalysis )

OfficialZAnalysis::OfficialZAnalysis( Sample& sample, const std::string & collectionFileName ):
        SampleAnalysis( "Official"+EventServer::Zmode[imode] , sample, collectionFileName ),
	outTree_(0)
{
        Zmode  = EventServer::Zmode[imode];

        cout << "\t-------------------------------------------" << endl; 
        cout << "\t---             Z  analysis            ----" << endl; 
        cout << "\t-------------------------------------------" << endl; 

        _nDebug = 10;

        // init counters
        Nevt_ = 0;
        NAcceptance_ = 0;
        NevtCand_ = 0;
        NHLT_  = 0;
        NPresel_ = 0;
        NIsol_ = 0;
        NID_ = 0;
        NMassCut_ = 0;
        NZ_ = 0;

        outTree_ =new TTree("METCommVariables","Variables for MET Splots");

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

        TcUparaQpara = new TProfile("TcUparaQpara","TcUparaQpara",50,0,500);
        PFUparaQpara = new TProfile("PFUparaQpara","PFUparaQpara",50,0,500);
        RawCaloUparaQpara = new TProfile("RawCaloUparaQpara","RawCaloUparaQpara",50,0,500);
        CaloUparaQpara = new TProfile("CaloUparaQpara","CaloUparaQpara",50,0,500);

        _e.a().verbosity = 0;

}

OfficialZAnalysis::~OfficialZAnalysis()
{
}

bool
OfficialZAnalysis::analyzeEvent()
{
        // _outputFile->cd();

  CatEta="";
  dau.clear();
 
  if((_ievt+1)%500 == 0)
    cout<<" Event "<<_ievt+1<<endl;
 
  //Elements
  Nevt_++;
  
  if(_e.inAcceptance() )
    NAcceptance_++;
   
  //Need one Z candidate or more
  const CandList& zllList = MakeZCandidates();
  if( zllList.size()==0) return false;
  else  NevtCand_++;
  Candidate* ZCand = BestZCand(zllList);
  fill( "n", Zmode+"_ZCand", zllList.size() );
  fill( "mass",Zmode+"_Begin", ZCand->mass() );

  //Compute ZCBM
  ZCBM_[0][0]=ZCand->cosPhi();
  ZCBM_[0][1]=ZCand->sinPhi();
  ZCBM_[1][0]=-ZCand->sinPhi();
  ZCBM_[1][1]=+ZCand->cosPhi();
  //======

  
  //End preselection events ***********
  // identical to Official VBTF analysis *******

  //Vertices
  if(SampleName.substr(0,5)!="ttbar")
    if(!TreeFillVertices() )return false;


  //Analysis Way
  if(!HLTFiring() ) return false; //1st Step
  else NHLT_++;                   //Disabled for test

  if(!EventPreselection(ZCand)) return false; //2nd Step
  else  NPresel_++;

  if(!LeptonID() ) return false; //3rd step, actually combined with iso
  else NID_++;
  
  if(!LeptonIso() ) return false; //4th step, cf ID
  else NIsol_++;
  
  if(!(ZCand->mass() > 40) )return false;
  
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
  
  PrepareTree(ZCand);
  
  if(!MassSelection( ZCand) ) return false; //Just to give an approximation
  else NMassCut_++;                         // not implemented in Official analysis
  
  return true;
  
}

void
OfficialZAnalysis::bookHistograms()
{
        //
        // histogram booking
        //
        outTree_->SetDirectory( _outputFile->GetDirectory(NULL) );
        //outTree_->SetDirectory( _outputFile->GetDirectory("") );
        initRootTree();

}

void
OfficialZAnalysis::writeHistograms()
{
        outTree_->Write();

        /*  TcUparaQpara->Write();
            PFUparaQpara->Write();
            CaloUparaQpara->Write();
            RawCaloUparaQpara->Write();

            for(int i=0;i<10;i++) {
            dPhiRecZ[i]->Write();
            RecpVsPt[i]->Write();
            RecVsPt[i]->Write();
            UQTvsQT[i]->Write();
            UpQTvsQT[i]->Write();

            RecpCorVsPt[i]->Write();
            RecpeCorVsPt[i]->Write();
            RecpCorVsGenPt[i]->Write();
            RecpeCorVsGenPt[i]->Write();

            }
            */
        cout<<" Number of events analyzed "<<Nevt_<<endl;
        cout<<" Number of events containing one or more ZCand "<<NevtCand_<<endl;
        cout<<" Number of events in Acceptance "<<NAcceptance_<<"  "<<((double)NAcceptance_*100)/Nevt_<<" %"<<endl;
        cout<<" Number of events firing HLT "<<NHLT_<<"  "<<((double)NHLT_*100)/NAcceptance_<<" %;  "<<((double)NHLT_*100)/Nevt_<<endl;
        cout<<" Number of events preselected "<<NPresel_<<"  "<<((double)NPresel_*100)/NHLT_<<" %;  "<<((double)NPresel_*100)/Nevt_<<endl;
        cout<<" Number of events identified "<<NID_<<"  "<<((double)NID_*100)/NPresel_<<" %;  "<<((double)NID_*100)/Nevt_<<" %"<<endl;
        cout<<" Number of events Isolated "<<NIsol_<<"  "<<((double)NIsol_*100)/NID_<<" %;  "<<((double)NIsol_*100)/Nevt_<<" %"<<endl;
        cout<<" Number of events Passing";
        cout<<" trought the Mass window "<<NMassCut_<<"  "<<((double)NMassCut_*100)/NIsol_<<" %; "<<((double)NMassCut_*100)/Nevt_<<" %"<<endl;

}

bool
OfficialZAnalysis::HLTFiring() {

        /* cout<<" HLt fired? " <<_e.isFired("HLT_Ele15_SW_L1R")<<endl;
           cout<<" HLt fired? " <<_e.isFired("HLT_Photon15_Cleaned_L1R")<<endl;
           cout<<" HLt fired? " <<_e.isFired("HLT_Ele15_LW_L1R")<<endl;*/

        bool fired=false;

        int run = _e.run();

        if(Zmode=="Zee" && SampleName.substr(0,1)=="E") { //FIXME MC...

                // cout<<_e.isFired("HLT_Ele17_SW_CaloEleId_L1R")<<"    "<<_e.isFired("HLT_Ele17_SW_TightEleId_L1R")<<endl;

                if( ( (run <= 140401) && _e.isFired("HLT_Ele15_LW_L1R") ) ||
                    ( (run >= 140402 && run <= 143962 ) && _e.isFired("HLT_Ele15_SW_L1R") ) ||
                    ( (run >= 143963 && run <= 144114 ) && _e.isFired("HLT_Ele15_SW_CaloEleId_L1R") ) ||
                    ( (run >= 144115 && run <= 147116 ) && _e.isFired("HLT_Ele17_SW_CaloEleId_L1R") ) ||
                    ( (run >= 147117 && run <= 148058 ) && _e.isFired("HLT_Ele17_SW_TightEleId_L1R") ) ||
                    ( (run >= 148059 && run <= 149064 ) && _e.isFired("HLT_Ele17_SW_TighterEleIdIsol_L1R_v2") ) ||
                    ( (run >= 149065 && run <= 149442 ) && _e.isFired("HLT_Ele17_SW_TighterEleIdIsol_L1R_v3") ) )
                        fired =true;

                //    fired =true;//FIXME MM
                return fired;
        }
        else {

                return true;//_e.isFired("HLT_Mu9");
        }
}


bool
OfficialZAnalysis::EventPreselection(Candidate* ZCand) {

        bool accept=false;

        bool Presel[2]={false,false};
        string categ_[2];

        //CandInfo* ZInfo = ZCand->info();
        float pt_0      =  ZCand->daughter(0)->pt();//ZInfo->getFloat("pt_0"); 
        float pt_1      =  ZCand->daughter(1)->pt();//ZInfo->getFloat("pt_1"); 

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
        for( int idau=0; idau<2; idau++ )
        {
                Candidate* dau_ = dau[idau];

                CandInfo* info = dau_->info();

                float eta_; 
                float Et_;


                if(Zmode=="Zee") {
                        eta_ = info->getFloat("caloEta");
                        Et_ = info->getFloat("caloEnergy")/cosh(eta_);
                }
                else {
                        eta_ = dau_->eta();
                        Et_ = dau_->pt();
                }

                if( (Zmode=="Zee" && fabs(eta_)<1.479) || (Zmode=="Zmm" && fabs(eta_)<1.3 ) )
                        categ_[idau]="B";
                else if( (Zmode=="Zee" && fabs(eta_)>=1.479) || (Zmode=="Zmm" && fabs(eta_)>=1.3 ) )
                        categ_[idau]="E";

                //FIXME MM, used for DY validation
                if(Zmode=="Zee" && Et_> 20 && ( fabs(eta_)<1.4442
                                                || (fabs(eta_)>1.56 && fabs(eta_)<2.5) ) )
                        Presel[idau]=true;
                if(Zmode=="Zmm" && Et_> 15 && ( fabs(eta_)<2.4 ) )
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
OfficialZAnalysis::LeptonIso() {

        float TrkIso[2];
        float EcalIso[2];
        float HcalIso[2];

        float etaSep;
        if(Zmode=="Zee") {
                etaSep=1.479;

                //Electron Working Point
                TrkIso[0]= 0.15; TrkIso[1]= 0.08;
                EcalIso[0]= 2.0; EcalIso[1]= 0.06;
                HcalIso[0]= 0.12; HcalIso[1]= 0.05;
                //*****************

        }
        else
                etaSep=1.3;

        bool accept=false; 
        bool acc[2]={false,false};

        for( int idau=0; idau<2; idau++ )
        {
                Candidate* dau_ = dau[idau];
                CandInfo* info_ = dau_->info();

                if(Zmode=="Zee") {

                        //	IsoElec.ecalUseIsoDeposits = false;
                        //	IsoElec.computeIso( *dau_ );

                        float trkIso = info_->getFloat("dr03TrkIso")/dau_->pt();
                        float ecalIso = info_->getFloat("dr03EcalIso")/dau_->pt();
                        float hcalIso = (info_->getFloat("dr03HcalD1Iso") + info_->getFloat("dr03HcalD2Iso") )/dau_->pt() ;

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
                }//Z ee mode

        }//End daughter loop

        if(acc[0] && acc[1])
                accept=true;

        return accept;

}

bool
OfficialZAnalysis::LeptonID() {

        //Cut definition ***
        //0 barrel, 1 endcap
        float dEtaIn[2];
        float dPhiIn[2];
        float sigiEiE[2];
        float HOE[2];

        // Electron Working Point
        dEtaIn[0]= 0.007; dEtaIn[1]= 0.01; 
        dPhiIn[0]= 0.8; dPhiIn[1]= 0.7; 
        sigiEiE[0]= 0.01; sigiEiE[1]= 0.03;
        HOE[0]= 0.15; HOE[1]= 0.07;
        //****************

        /*
           float decor[2];
           float dpcor[2];
           if(Zmode=="Zee") {
           for( int idau=0; idau<2; idau++ ) {
           Candidate* dau_ = dau[idau];
           CandInfo* info_ = dau_->info();
           decor[idau] = vbtfSel.detaCorrections(info_->getFloat("caloEta"),dau_->phi());
           dpcor[idau] = vbtfSel.dphiCorrections(info_->getFloat("caloEta"),dau_->phi());		    
           }
           }
           else
           {
           for( int idau=0; idau<2; idau++ ) {
           decor[idau] =0;
           dpcor[idau] =0;
           }
           }
           */
        //    Candidate* ZCand = Candidate::create(dau[0],dau[1]);
        // PrepareTree(decor,dpcor,ZCand);

        bool accept=false;
        bool acc[2]={false,false};
        for( int idau=0; idau<2; idau++ ) {
                Candidate* dau_ = dau[idau];
                CandInfo* info_ = dau_->info();
                if(Zmode=="Zee") {

                        float etaSep = 1.479;

                        /* float detaCor = vbtfSel.detaCorrections(info_->getFloat("caloEta"),dau_->phi());
                           float dphiCor = vbtfSel.dphiCorrections(info_->getFloat("caloEta"),dau_->phi());		    

                           if(SampleName.substr(0,8)!="FullLumi" && SampleName.substr(0,2)!="EG" &&  SampleName.substr(0,2)!="SD")
                           {
                           detaCor = 0;
                           dphiCor =0;
                           }*/

                        float hoe = info_->getFloat("HOE");
                        float dphi = info_->getFloat("dPhiIn") ;
                        float deta = info_->getFloat("dEtaIn") ;
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
                        acc[idau] = info_->getBool("muidGlobalMuonPromptTight");
                        // PrepareTree(0,0,ZCand);
                }  
        }

        if(acc[0] && acc[1])
                accept=true;

        return accept;
}

bool 
OfficialZAnalysis::MassSelection(Candidate* ZCand) {

        bool accept=false;

        if(ZCand->mass()> 70 && ZCand->mass()< 120 )
                accept=true;

        return accept;

}

void
OfficialZAnalysis::FillingHistos(Candidate* ZCand, const string & type) {

        //Filling Z histos
        fill( "pt", Zmode+type+"_"+CatEta  , ZCand->pt() );
        fill( "mass", Zmode+type+"_"+CatEta , ZCand->mass() );
        fill( "eta", Zmode+type+"_"+CatEta, ZCand->eta() );
        fill( "phi", Zmode+type+"_"+CatEta, ZCand->phi( kDeg ) );

        fill( "pt", Zmode+type , ZCand->pt() );
        fill( "mass", Zmode+type , ZCand->mass() );
        fill( "eta", Zmode+type, ZCand->eta() );
        fill( "phi", Zmode+type, ZCand->phi( kDeg ) );

        //Filling daughters histos
        static string sdau[2]={"_dau_1","_dau_2"};
        for(int idau=0;idau<2;idau++) {
                Candidate* dau_ = dau[idau];

                fill(  "pt", Zmode+type+"_"+sdau[idau], dau_->pt()  );
                fill( "eta", Zmode+type+"_"+sdau[idau],  dau_->eta() );
                fill( "phi", Zmode+type+"_"+sdau[idau],  dau_->phi()*180/Constants::pi );
        }
        float dphi = KineUtils::dPhi(dau[0]->phi(),dau[1]->phi());
        fill( "dphi", Zmode+type, dphi*180/Constants::pi );

}

Candidate *
OfficialZAnalysis::BestZCand( const CandList& zllList ) {

        Candidate* BCand;
        vector<float> tmp(2,0);
        vector<vector<float> > Pts;
        for(int unsigned iz=0;iz<zllList.size();iz++) {

                Pts.push_back(tmp);

                Candidate* ZCand = zllList[iz];

                float pt_0      =  ZCand->daughter(0)->pt();//ZInfo->getFloat("pt_0"); 
                float pt_1      =  ZCand->daughter(1)->pt();//ZInfo->getFloat("pt_1"); 

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
        return BCand;
}


void
OfficialZAnalysis::FillUnderlyingPhoton(Candidate* ZCand ) {

        const CandList& Photons = _e.photons();
        Candidate* thePhoton_;
        bool found=false;
        float pttmp=0;
        for(int unsigned ip=0;ip<Photons.size();ip++) { 
                Candidate* photon = Photons[ip];
                {
                        float dr1_ = KineUtils::dR(dau[0]->eta(), photon->eta(),
                                                   dau[0]->phi(), photon->phi() ); 
                        float dr2_ = KineUtils::dR(dau[1]->eta(), photon->eta(),
                                                   dau[1]->phi(), photon->phi() );
                        if(pttmp<photon->pt() && dr1_>0.1 && dr2_>0.1 )
                        { thePhoton_=photon; pttmp = photon->pt(); found=true; }
                }
                fill( "pt", Zmode+CatEta+"_AllPhotons", photon->pt() );
        }

        fill( "n", Zmode+CatEta+"_NPhotons" , Photons.size() );
        if( Photons.size()!=0 && pttmp!=0) {
                fill( "phi", Zmode+"_LeadPhoton", thePhoton_->phi( kDeg ) );
                fill( "pt", Zmode+"_LeadPhoton" , thePhoton_->pt() );
                fill( "dphi", Zmode+"_LeadPhoton", KineUtils::dPhi(ZCand->phi(), thePhoton_->phi() )*180/3.1415 );

                //Including photon in the invariant mass
                /* TVector3 v3= dau[0]->p3() + dau[1]->p3() + thePhoton_->p3();
                   float Minv2= pow(dau[0]->E()+dau[1]->E()+thePhoton_->E(),2)
                   - pow(v3.Mag(),2);
                   fill( "mass", Zmode+"_InvPhoton",sqrt(Minv2));
                   fill( "massvsdpt", Zmode+"_MassVsPtPhoton", sqrt(Minv2), fabs( (thePhoton_->pt()-ZCand->pt())/thePhoton_->pt()) );*/
        }


}

void
OfficialZAnalysis::FillUnderlyingJets(Candidate* ZCand ) {

        const CandList& Jets = _e.jetList( EventManager::kPatJet);
        Candidate* theJet_;

        float pttmp=0;
        for(int unsigned ij=0;ij<Jets.size();ij++) { 
                Candidate* jet = Jets[ij];
                {
                        float dr1_ = KineUtils::dR(dau[0]->eta(), jet->eta(),
                                                   dau[0]->phi(), jet->phi() ); 
                        float dr2_ = KineUtils::dR(dau[1]->eta(), jet->eta(),
                                                   dau[1]->phi(), jet->phi() );
                        if(pttmp<jet->pt() && dr1_>0.05 && dr2_>0.05  )
                        { theJet_=jet; pttmp = jet->pt(); }
                }
                fill( "pt", Zmode+CatEta+"_AllJets", jet->pt() );
        }
        fill( "n", Zmode+CatEta+"_NJets" , Jets.size() );
        if( Jets.size()!=0 && pttmp!=0) {
                fill( "phi", Zmode+CatEta+"_LeadJet", theJet_->phi( kDeg ) );
                fill( "pt", Zmode+CatEta+"_LeadJet" , theJet_->pt() );
                fill( "dphi", Zmode+CatEta+"_LeadJet", KineUtils::dPhi(ZCand->phi(), theJet_->phi() )*180/3.1415 );
        }

}



void
OfficialZAnalysis::FillUnderlyingMET(Candidate* ZCand ) {

        const Candidate* met_  = _e.met( EventManager::kPatMet);
        fill( "pt", Zmode+CatEta+"_MET", met_->Et() );
        fill( "phi", Zmode+CatEta+"_MET", met_->phi() );

        float dPhi = KineUtils::dPhi(met_->phi(),  ZCand->phi() )*180/3.1415;
        fill( "dphi", Zmode+CatEta+"_MET", dPhi );

}


CandList OfficialZAnalysis::MakeZCandidates() {

        const CandList* leptons;

        if(Zmode=="Zee")
                leptons = &(_e.electrons());
        else
                leptons = &(_e.muons());

        //Making Z Candidates
        CandList ZCands;
        for(int unsigned il1=0;il1<leptons->size();il1++) {

                Candidate* l1 = (*leptons)[il1];

                for(int unsigned il2=0;il2<leptons->size();il2++) {

                        if(il1 != il2) {
                                Candidate* l2 = (*leptons)[il2];

                                Candidate* ZCand = Candidate::create(l1,l2);
                                ZCands.push_back(ZCand);
                        }
                }
        }
        return ZCands;

}



CandList OfficialZAnalysis::MakeZCandidatesWithCharge() {

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
OfficialZAnalysis::ComputeRecoil(const Candidate* met, Candidate* ZCand, double cor) {

        TVector2 Recoil(0,0);
        //  cout<<" MET "<<met->px()<<"   "<<met->py()<<"   "<<met->phi()<<endl;

        Recoil -= met->p2()*cor + ZCand->p2();
        return Recoil;

}


void
OfficialZAnalysis::FillRecoil(TVector2 Recoil, Candidate* ZCand, const string & type){
        /*if(type.substr(0,2)=="Tr") {
          cout<<"  event   "<<_e.eventInFile()<<endl;
          cout<<"Recoil "<<Recoil.Px()<<"   "<<Recoil.Py()<<"    "<<TVector2::Phi_0_2pi(Recoil.Phi() )<<endl;
          cout<<"Z "<<ZCand->px()<<"   "<<ZCand->py()<<"   "<<ZCand->phi()<<endl;
          }*/
        fill( "pt",type+"Recoil",Recoil.Mod() );
        //Recoil angle btw 0,2 pi



        Double_t rphi = TVector2::Phi_0_2pi(Recoil.Phi() ) ;
        fill( "dPhi",type+"Recoil", KineUtils::dPhi(rphi,ZCand->phi() )*180/3.141592 );
        fill( "dPhivsMag",type+"Recoil", Recoil.Mod() , KineUtils::dPhi(rphi,ZCand->phi() )*180/3.141592);

        fill( "dPhivsZMag",type+"Recoil", ZCand->pt() , KineUtils::dPhi(rphi,ZCand->phi() )*180/3.141592);

        //Recoil Projection
        float Rpara= ZCBM_[0][0]*Recoil.Px()+ZCBM_[0][1]*Recoil.Py();
        float Rperp= ZCBM_[1][0]*Recoil.Px()+ZCBM_[1][1]*Recoil.Py();
        // if(type.substr(0,2)=="Tr")
        //  cout<<" proj "<<Rpara<<"   "<<Rperp<<endl;

        //  float response = fabs(Rpara/ZCand->pt());

        fill("noAbsPt",type+"Recoil_ParaZCor", Rpara);
        fill("noAbsPt",type+"Recoil_PerpZCor", Rperp);

        fill("noAbsPt",type+"Recoil_ParaZ", Rpara);
        fill("noAbsPt",type+"Recoil_PerpZ", Rperp);

        //Recoil projection in XY frame
        fill("noAbsPt",type+"Recoil_X", Recoil.Px() );
        fill("noAbsPt",type+"Recoil_Y", Recoil.Py() );

        //Recoil Projection in leptons frame
        TVector2 pN0 = ZCand->daughter(0)->p2().Unit();
        TVector2 pN1 = ZCand->daughter(1)->p2().Unit();
        TVector2 Plep = pN0 + pN1;

        float Leptangle =  TVector2::Phi_0_2pi(Plep.Phi() );

        float RLpara = cos(Leptangle)*Recoil.Px() + sin(Leptangle)*Recoil.Py();
        float RLperp = -sin(Leptangle)*Recoil.Px() + cos(Leptangle)*Recoil.Py();

        float ZLpara = cos(Leptangle)*ZCand->px() + sin(Leptangle)*ZCand->py();
        float ZLperp = -sin(Leptangle)*ZCand->px() + cos(Leptangle)*ZCand->py();

        fill("noAbsPt",type+"Recoil_ParaLep", RLpara);
        fill("noAbsPt",type+"Recoil_PerpLep", RLperp);

        fill("ptVsPt",type+"Recoil_Z",  ZCand->pt(),  Recoil.Mod() );
        fill("ptVsPt",type+"RecoilBis_ZBis", ZLpara, RLpara );

        fill("ptVsPt",type+"Ut_vs_Qt", ZCand->pt(), Recoil.Mod() );
        fill("UQtVsQt",type+"UtQtratio_vs_Qt", ZCand->pt(),Recoil.Mod()/ZCand->pt() );

        fill( "dPhivsMag", type+"Z_Bis", ZCand->pt(), KineUtils::dPhi(Leptangle,ZCand->phi() )*180/Constants::pi );


        fill("ptVsPt",type+"Recoil_ParaLep_Z", ZLpara, RLpara );
        fill("ptVsPt",type+"Recoil_PerpLep_Z", ZLperp, RLperp );
        //*************

        //  Candidate* ZGenCand;


        NZ_++;

        ///Recoil Projection in lepton SCs frame
        if (type.substr(0,1)=="C") {
                Candidate* sC1_ = _e.getSecond("electron-SuperCluster",ZCand->daughter(0) );
                Candidate* sC2_ = _e.getSecond("electron-SuperCluster",ZCand->daughter(1) );
                if(sC1_!=NULL && sC2_!=NULL ) {
                        float SCangle =  TVector2::Phi_0_2pi(Plep.Phi() );

                        float RSCpara = cos(SCangle)*Recoil.Px() + sin(SCangle)*Recoil.Py();
                        float RSCperp = -sin(SCangle)*Recoil.Px() + cos(SCangle)*Recoil.Py();

                        fill("noAbsPt",type+"Recoil_ParaSCBis", RSCpara);
                        fill("noAbsPt",type+"Recoil_PerpSCBis", RSCperp);

                }
        }
        //Fill angles associated with recoil
        float RoPt = KineUtils::dPhi(rphi, ZCand->phi() )*180/(Constants::pi) ;
        float RoB = KineUtils::dPhi(rphi , Leptangle)*180/(Constants::pi) ;
        float PtoB = KineUtils::dPhi(ZCand->phi(), Leptangle)*180/(Constants::pi)  ;

        fill( "dPhi",type+"RoPt",RoPt);
        fill( "dPhi",type+"RoB",RoB );
        fill( "dPhi",type+"PtoB", PtoB);

}


void
OfficialZAnalysis::FillMETProjection(const Candidate* met, Candidate* ZCand, const string & type) {

        //MET Projection on Z axis
        float METpara= ZCBM_[0][0]*met->px()+ZCBM_[0][1]*met->py();
        float METperp= ZCBM_[1][0]*met->px()+ZCBM_[1][1]*met->py();

        fill("noAbsPt",type+"MET_ParaZ", METpara);
        fill("noAbsPt",type+"MET_PerpZ", METperp);

        //MET Projections in XY frame
        fill("noAbsPt",type+"MET_X", met->px() );
        fill("noAbsPt",type+"MET_Y", met->py() );

        //MET Projection in leptons frame
        float phi_0 = ZCand->daughter(0)->phi();
        float phi_1 = ZCand->daughter(1)->phi();
        float Leptangle = (phi_0 + phi_1)/2.;



        float METLpara = cos(Leptangle)*met->px() + sin(Leptangle)*met->py();
        float METLperp = -sin(Leptangle)*met->px() + cos(Leptangle)*met->py();

        if( fabs(Leptangle) > Constants::pi) {
                METLpara *= -1;
                METLperp *= -1;
        }

        fill("noAbsPt",type+"MET_ParaLep", METLpara);
        fill("noAbsPt",type+"MET_PerpLep", METLperp);

}



void
OfficialZAnalysis::FillMETVariables(const Candidate* met, Candidate* ZCand, const string & type) {


        Float_t Bins2[14] = {0,10,20,30,40,50,60,80,100,150,200,250,300,350};
        double CorrectionCalo[13]={0.354796,0.450207,0.520878,0.560771,0.594442,0.64499,0.666981,0.708246,0.720698,0.754153,0.80916,0.756484,0.8};
        double CorrectionCT1[13]={0.35537,0.585815,0.789873,0.902999,0.959416,1.01803,1.02837,1.05108,1.025,1.02194,1.03689,0.981633,1};
        double CorrectionCT2[13]={0.837067,0.978843,1.03213,1.0504,1.06541,1.09021,1.09502,1.1033,1.06841,1.0549,1.06635,0.970004,1};
        double CorrectionTc[13]={0.604812,0.705049,0.779874,0.816408,0.833187,0.86844,0.901517,0.91991,0.922474,0.918964,0.967263,0.901705,1};
        double CorrectionPF[13]={0.743238,0.819149,0.882104,0.911021,0.914964,0.933811,0.949826,0.950231,0.926152,0.918988,0.949833,0.931791,1};
        float correction =1;

        int bin=12;
        for(int i=0;i<13;i++) {
                if(ZCand->pt()<Bins2[i+1] && ZCand->pt()>=Bins2[i] )
                { bin = i; break;}
        }
        //  cout<<" ZCand "<<ZCand->pt()<<"   "<<bin<<endl;
        if(type.substr(0,3)=="Raw")
                correction = 1/CorrectionCalo[bin];
        else if(type.substr(0,2)=="PF")
                correction = 1/CorrectionPF[bin];
        else if(type.substr(0,3)=="Trk")
                correction = 1/CorrectionTc[bin];
        else if(type.substr(0,10)=="CaloTypeII")
                correction = 1/CorrectionCT2[bin];
        else
                correction = 1/CorrectionCT1[bin];
        // cout<<type <<"   "<<correction<<endl;
        fill( "dphi", type+"Dphi_"+CatEta, KineUtils::dPhi(met->phi(),ZCand->phi() )*180/3.141592 );

        correction = 1;
        // cout<<"  "<<correction<<endl;
        TVector2 GenRecoil = ComputeGenRecoil(met);
        genRecoil = GenRecoil;

        fill("pt", type , met->pt() );
        FillMETProjection(met, ZCand, type); 

        TVector2 Recoil = ComputeRecoil(met, ZCand,1);
        FillRecoil(Recoil*correction, ZCand, "Basic"+type);

        TVector2 CaloRecoil = ComputeRecoilWithCalo(met, ZCand,1);   


        //TVector2 GenRecoil = ComputeGenRecoil(met);

        if( GenRecoil.Px()!=0 && GenRecoil.Py()!=0 )
                FillRecoil(GenRecoil, ZCand, type+"Gen");

        TVector2 PFRecoil = ComputeRecoilWithPF(met, ZCand,1);

        if( type.substr(0,2)=="PF" && PFRecoil.Px()!=0 && PFRecoil.Py()!=0 ) {
                // FillRecoil(Recoil, ZCand, type);

                FillRecoil(PFRecoil*correction, ZCand, type+"PF");
        }
        else if( type.substr(0,2)=="Tr") {

                FillRecoil(Recoil*correction, ZCand, type);
        }
        else if( type.substr(0,2)!="PF" && type.substr(0,2)!="Tr"  && CaloRecoil.Px()!=0 && CaloRecoil.Py()!=0 ) {
                // FillRecoil(Recoil, ZCand, type);

                FillRecoil(CaloRecoil*correction, ZCand, type+"Calo");
                //  FillRecoil(CalRecHitRecoil, ZCand, type+"RecHit");
        }

}


void
OfficialZAnalysis::FillSeveralMET( Candidate* ZCand, const string & ALevel) {

        //string type;

        //With PF MET
        const Candidate* pfmet_ = _e.met( EventManager::kPfMet );
        FillMETVariables(pfmet_, ZCand, "PF"+ALevel);

        //With other MET : Calo Corrected T1 =======
        _e.refreshCaloMET(); //to be sure we take the good one
        const Candidate* calomet_  = _e.met( EventManager::kCaloMet);
        FillMETVariables(calomet_, ZCand, "Calo"+ALevel);
        //========================================
        //With Gen MET ========
        //const Candidate* genmet_ = _e.met( EventManager::kGenMet);
        //FillMETVariables(genmet_, "Gen"+ALevel);
        //========================================

        //With Track Corr Met
        const Candidate* trkcmet_  = _e.met( EventManager::kPatMet);
        FillMETVariables(trkcmet_, ZCand, "TrkC"+ALevel);
        //========================================

        //With Raw Calo MET
        _e.refreshCaloMET(); //to be sure we take the good one

        int v=-1;
        for(int unsigned i=0;i<_e.a().caloMET.size();i++) {
                if("met" ==_e.a().caloMET[i])
                { v=i; break; }
        }
        const Candidate* cmet_ = _e.met( EventManager::kCaloMet, v );
        FillMETVariables(cmet_, ZCand, "RawCalo"+ALevel);
        //========================================

        //With  Calo Corrected T2
        _e.refreshCaloMET(); //to be sure we take the good one

        v=-1;
        for(int unsigned i=0;i<_e.a().caloMET.size();i++) {
                if("metMuonJESCorAK5" ==_e.a().caloMET[i])
                { v=i; break; }
        }
        const Candidate* cmetT2_ = _e.met( EventManager::kCaloMet, v );
        FillMETVariables(cmetT2_, ZCand, "CaloTypeII"+ALevel);
        //========================================

} 


TVector2
OfficialZAnalysis::ComputeRecoilWithCaloRecHit(const Candidate* met, Candidate* ZCand) {

        TVector2 Recoil(0,0);;

        TVector2 sumRH_(0,0);

        //const CandAssoc& e_Sc_ =  _e.candAssoc("electron-SuperCluster");
        const CandMap& Sc_Bc_ = _e.candMap("SuperCluster-BasicCluster");
        const CandIdRecHitMap& Bc_rh_ = _e.ecalRecHitMap("BasicCluster-EcalRecHit");


        //looping on the two electrons
        for( int idau=0; idau<2; idau++ ) {
                Candidate* electron_ = ZCand->daughter(idau)->theBase();
                Candidate* superCluster = _e.getSecond("electron-SuperCluster", electron_ );
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
        }


        Recoil -=  met->p2() + sumRH_;
        return Recoil;

}

TVector2
OfficialZAnalysis::ComputeRecoilWithCalo(const Candidate* met, Candidate* ZCand, double cor) {

        TVector2 Recoil(0,0);;
        TVector2 sumSC_(0,0);
        bool test[2]={false,false};
        for( int idau=0; idau<2; idau++ ) {
                Candidate* lepton_ = ZCand->daughter(idau)->theBase();

                Candidate* superCluster_ = _e.getSecond("electron-SuperCluster",lepton_ );
                if(superCluster_!=NULL)
                { sumSC_ +=superCluster_->p2(); test[idau]=true; }
        }
        if(test[0] && test[1] )
                Recoil -=  met->p2()*cor + sumSC_;
        return Recoil;

}

TVector2
OfficialZAnalysis::ComputeRecoilWithPF(const Candidate* met, Candidate* ZCand, double cor) {

        TVector2 Recoil(0,0);
        TVector2 sumPF_(0,0);

        bool test[2]={false,false};

        for( int idau=0; idau<2; idau++ ) {
                Candidate* lepton_ = ZCand->daughter(idau)->theBase();
                Candidate* pfcand_ = _e.getSecond("electron-PfCand",lepton_ );
                if(pfcand_!=NULL)
                { sumPF_ +=pfcand_->p2(); test[idau]=true; }
        }
        if(test[0] && test[1] )
                Recoil -= met->p2()*cor + sumPF_;
        return Recoil;

}

TVector2
OfficialZAnalysis::ComputeRecoilWithTowers(const Candidate* met, Candidate* ZCand, double cor) {

        TVector2 Recoil(0,0);
        TVector2 sumCT_(0,0);

        const CandList& caloTowers_ = _e.caloTowers();

        for(int i=0;i<2;i++) {

                Candidate* lepton_ = ZCand->daughter(i)->theBase();

                float dR;
                for(int unsigned i=0;i<caloTowers_.size();i++) {
                        Candidate* caloTower = caloTowers_[i];

                        dR = KineUtils::dR(lepton_->eta(),caloTower->eta(),lepton_->phi(),caloTower->phi());
                        if(dR<0.3)
                                sumCT_ += caloTower->p2();
                }
        }
        Recoil -= met->p2()*cor + sumCT_;
        return Recoil;
}



TVector2
OfficialZAnalysis::ComputeEtVectorFromRecHit(ecalRecHit rh) {


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
        //  cout<<" rechit "<< rh.ix<<"<->"<<eta_<<"  "<< rh.iy<<"<->"<<phi_<<"  "
        //      << rh.iz<<"  "<<rh.e<<"  "<<et_<<"  <=>  vec "
        //      <<etVector.Px()<<"  "<<etVector.Py()<<"  "<<etVector.Mod()<<endl;
        return etVector;

}

TVector2
OfficialZAnalysis::ComputeGenRecoil(const Candidate* met) {

        TVector2 Recoil(0,0);

        TVector2 sumGen_(0,0);
        bool test[2]={false,false};

        const CandList& genList = _e.mcTruthCandidates();
        // cout<<genList.size()<<endl;
        for(int idau=0;idau<2;idau++) {
                Candidate* dau_ = dau[idau];
                //   cout<<dau_->pt()<<endl;
                for(int unsigned imc=0;imc<genList.size();imc++) {

                        Candidate* genCand = genList[imc];
                        float dR = KineUtils::dR(genCand->eta(), dau_->eta(), genCand->phi(), dau_->phi() );
                        if(abs(genCand->pdgCode())==abs(dau_->pdgCode()) && dR<0.1 )
                        {
                                //	cout<<genCand->pt()<<"  "<<genCand->pdgCode()<<endl;
                                sumGen_ += genCand->p2(); test[idau]=true; break; 
                        }
                }
        }
        // cout<< " test "<<test[0]<<"  "<<test[1]<<sumGen_.Px()<<endl;
        if(test[0] && test[1] )
                Recoil =  met->p2() + sumGen_;
        return sumGen_;

}



//*************************************=========================
//*************************************=========================



void 
OfficialZAnalysis::initRootTree() {

  
 //Event
  outTree_->Branch("Run",&Run,"Run/I");
  outTree_->Branch("Event",&Event,"Event/I");
  outTree_->Branch("AbsEvent",&AbsEvent,"AbsEvent/I");
  outTree_->Branch("sample",&sampleName,"sample[21]/C"); 
  outTree_->Branch("fileName",&fileName,"fileName[21]/C"); 
  outTree_->Branch("timestamp",&_time,"time/I");
 

  //Selection
  outTree_->Branch("InHit1",&inHit1,"inHit1/B");
  outTree_->Branch("ExpInHit1",&expInHit1,"expInHit1/B");
  outTree_->Branch("ConvRej1",&ConvRej1,"ConvRej1/B");

  outTree_->Branch("InHit2",&inHit2,"inHit2/B");
  outTree_->Branch("ExpInHit2",&expInHit2,"expInHit2/B");
  outTree_->Branch("ConvRej2",&ConvRej2,"ConvRej2/B");

  //Z
  outTree_->Branch("ZPt",&ZPt,"ZPt/F");
  outTree_->Branch("ZMass",&ZMass,"ZMass/F");
  outTree_->Branch("GenZPt",&GenZPt,"GenZPt/F");
  outTree_->Branch("GenZdPhi",&dPhiGenZ,"dPhiGenZ/F");
  
  outTree_->Branch("ZPtSC",&ZPtSC,"ZPtSC/F");
  outTree_->Branch("ZMassSC",&ZMassSC,"ZMassSC/F");
  outTree_->Branch("ZPtSCUnCor",&ZPtSCUnCor,"ZPtSCUnCor/F");
  outTree_->Branch("ZMassSCUnCor",&ZMassSCUnCor,"ZMassSCUnCor/F");
  outTree_->Branch("SCCor",SCCor,"SCCor[2][3]/F");
  outTree_->Branch("UnCorSC",UnCorSC,"UnCorSC[2][3]/F");

  //DeltaZ between vertices
  outTree_->Branch("DZvertices",&dZ,"dZ/F");

  //Jet
  outTree_->Branch("VetoJet",&VetoJet,"VetoJet/B");
  outTree_->Branch("JetMultiplicity",&JetMultiplicity,"JetMultiplicity/I");
  outTree_->Branch("LeadingJet",&LJetPt,"LJetPt/F");
  outTree_->Branch("dPhiJetZ",&dPhiJetZ,"dPhiJetZ/F");
  
  //Photon
  outTree_->Branch("PhotonMultiplicity",&PhotonMultiplicity,"PhotonMultiplicity/I");
  outTree_->Branch("LeadingPhoton",&LPhotonPt,"LPhotonPt/F");
  outTree_->Branch("dPhiPhotonZ",&dPhiPhotonZ,"dPhiPhotonZ/F");

  outTree_->Branch("NPhot1",&NPhot1,"NPhot1/I");
  outTree_->Branch("NPhot2",&NPhot2,"NPhot2/I");
  //  outTree_->Branch("",,"");

  //Angles
  outTree_->Branch("dPhiMETZ",dPhiMETZ,"dPhiMETZ[6]/F");
  outTree_->Branch("dPhiMETJet",dPhiMETJet,"dPhiMETJet[6]/F");

  outTree_->Branch("dPhiRecoilZ",dPhiRecoilZ,"dPhiRecoilZ[6]/F");
  outTree_->Branch("BasicdPhiRecoilZ",BasicdPhiRecoilZ,"BasicdPhiRecoilZ[6]/F");

  //Correction
  outTree_->Branch("Corrections",Corrections,"Corrections[6]/F");

  //Corrected MET variables
  outTree_->Branch("MetCor",MetCor,"MetCor[6][2]/F");
  outTree_->Branch("MetCorProj",MetCorProj,"MetCorProj[6][2]/F");
  outTree_->Branch("RecoilCor",RecoilCor,"RecoilCor[6][2]/F");
  outTree_->Branch("RecoilCorProj",RecoilCorProj,"RecoilCorProj[6][2]/F");
  outTree_->Branch("RecoilBisProj",BisectorRecoilProj,"RecoilBisProj[6][2]/F");

  //Lepton
  outTree_->Branch("IDVar",IDVar,"IDVar[2][6]/F");
  outTree_->Branch("IsoVar",IsoVar,"IsoVar[2][3]/F"); 
  outTree_->Branch("Lepton",Lepton,"Lepton[2][3]/F");
  outTree_->Branch("PFLepton",PFLepton,"PFLepton[2][3]/F");
  outTree_->Branch("SCLepton",SCLepton,"SCLepton[2][3]/F");
  outTree_->Branch("EtSC",&EtSC,"EtSC[2]/F");

  outTree_->Branch("Charge",&charge,"charge[2]/I");

  //SumET
  outTree_->Branch("SumET",SumET,"SumET[6]/F");
  
  //Pf MET variables
  outTree_->Branch("cEmFrac",&cEmFrac,"cEmFrac/F");
  outTree_->Branch("nEmFrac",&nEmFrac,"nEmFrac/F");
  outTree_->Branch("cHadFrac",&cHadFrac,"cHadFrac/F");
  outTree_->Branch("nHadFrac",&nHadFrac,"nhadFrac/F");
  outTree_->Branch("cMuFrac",&cMuFrac,"cMuFrac/F");
  
  //Decomposed MET
  outTree_->Branch("CaloSumETEB",&CaloSumETEB,"CaloSumETEB/F");
  outTree_->Branch("CaloSumETEE",&CaloSumETEE,"CaloSumETEE/F");
  outTree_->Branch("CaloSumETHB",&CaloSumETHB,"CaloSumETHB/F");
  outTree_->Branch("CaloSumETHE",&CaloSumETHE,"CaloSumETEE/F");
  outTree_->Branch("CaloSumETHF",&CaloSumETHF,"CaloSumETHF/F");

  outTree_->Branch("CaloMETEB",CaloMETEB,"CaloMETEB[2]/F");
  outTree_->Branch("CaloMETEE",CaloMETEE,"CaloMETEE[2]/F");
  outTree_->Branch("CaloMETHB",CaloMETHB,"CaloMETHB[2]/F");
  outTree_->Branch("CaloMETHE",CaloMETHE,"CaloMETEE[2]/F");
  outTree_->Branch("CaloMETHF",CaloMETHF,"CaloMETHF[2]/F");
  
  outTree_->Branch("PfSumETEB",&PfSumETEB,"PfSumETEB/F");
  outTree_->Branch("PfSumETEE",&PfSumETEE,"PfSumETEE/F");
  outTree_->Branch("PfSumETHB",&PfSumETHB,"PfSumETHB/F");
  outTree_->Branch("PfSumETHE",&PfSumETHE,"PfSumETEE/F");
  outTree_->Branch("PfSumETHF",&PfSumETHF,"PfSumETHF/F");
  
  outTree_->Branch("PfMETEB",PfMETEB,"PfMETEB[2]/F");
  outTree_->Branch("PfMETEE",PfMETEE,"PfMETEE[2]/F");
  outTree_->Branch("PfMETEEx",PfMETEEx,"PfMETEEx[2]/F");
  outTree_->Branch("PfMETHB",PfMETHB,"PfMETHB[2]/F");
  outTree_->Branch("PfMETHE",PfMETHE,"PfMETEE[2]/F");
  outTree_->Branch("PfMETHF",PfMETHF,"PfMETHF[2]/F");
  
  outTree_->Branch("PfEcalMET",PfEcalMET,"PfEcalMET[2]/F");
  outTree_->Branch("PfHcalMET",PfHcalMET,"PfHcalMET[2]/F");
  
  outTree_->Branch("CaloEcalMET",CaloEcalMET,"CaloEcalMET[2]/F");
  outTree_->Branch("CaloHcalMET",CaloHcalMET,"CaloHcalMET[2]/F");

  //Init MET
  std::stringstream buffer;
  buffer << "E"
         <<"[" << sizeof(Mets)/sizeof(Mets[0]) << "]"
         << "[" << sizeof(Mets[0])/sizeof(Mets[0][0]) << "]/F";
  outTree_->Branch("Mets", Mets, buffer.str().c_str());
  buffer.str("");
  buffer << "E"
         <<"[" << sizeof(MetProj)/sizeof(MetProj[0]) << "]"
         << "[" << sizeof(MetProj[0])/sizeof(MetProj[0][0]) << "]/F";
  outTree_->Branch("MetProj", MetProj, buffer.str().c_str());
  buffer.str("");
  //Init Recoil
  buffer << "E"
         <<"[" << sizeof(Recoils)/sizeof(Recoils[0]) << "]"
         << "[" << sizeof(Recoils[0])/sizeof(Recoils[0][0]) << "]/F";
  outTree_->Branch("Recoils", Recoils, buffer.str().c_str());
  buffer.str("");
  buffer << "E"
         <<"[" << sizeof(RecoilProj)/sizeof(RecoilProj[0]) << "]"
         << "[" << sizeof(RecoilProj[0])/sizeof(RecoilProj[0][0]) << "]/F";
  outTree_->Branch("RecoilProj", RecoilProj, buffer.str().c_str());
  buffer.str("");

  //Init Basic Recoil
  buffer << "E"
         <<"[" << sizeof(BasicRecoils)/sizeof(BasicRecoils[0]) << "]"
         << "[" << sizeof(BasicRecoils[0])/sizeof(BasicRecoils[0][0]) << "]/F";
  outTree_->Branch("BasicRecoils", BasicRecoils, buffer.str().c_str());
  buffer.str("");
  buffer << "E"
         <<"[" << sizeof(BasicRecoilProj)/sizeof(BasicRecoilProj[0]) << "]"
         << "[" << sizeof(BasicRecoilProj[0])/sizeof(BasicRecoilProj[0][0]) << "]/F";
  outTree_->Branch("BasicRecoilProj", BasicRecoilProj, buffer.str().c_str());
  buffer.str("");
  buffer << "E"
         <<"[" << sizeof(GenRecoil)/sizeof(GenRecoil[0]) << "]"
         << "[" << sizeof(GenRecoil[0])/sizeof(GenRecoil[0][0]) << "]/F";
  outTree_->Branch("GenRecoil", GenRecoil, buffer.str().c_str());
  buffer.str("");
  buffer << "E"
         <<"[" << sizeof(CTRecoil)/sizeof(CTRecoil[0]) << "]"
         << "[" << sizeof(CTRecoil[0])/sizeof(CTRecoil[0][0]) << "]/F";
  outTree_->Branch("CTRecoil", CTRecoil, buffer.str().c_str());
  buffer.str("");
  
  outTree_->Branch("NVertex",&Nvertex,"Nvertex/I");
 

}


void
OfficialZAnalysis::PrepareTree(/*float detacor[2], float dphicor[2],*/ Candidate* ZCand ) {

        float ID[2][6];
        float Iso[2][3];

        for(int i=0;i<2;i++) {
                Candidate* dau_ = dau[i];
                CandInfo* info = dau_->info();


                if(Zmode=="Zee") {
                        ID[i][0] = info->getFloat("sigiEtaiEta");
                        ID[i][1] = info->getFloat("dEtaIn");
                        ID[i][2] = info->getFloat("dPhiIn");
                        ID[i][3] = info->getFloat("HOE");
                        ID[i][4] = 0;//detacor[i];
                        ID[i][5] = 0;//dphicor[i];

                        Iso[i][0] = info->getFloat("dr03TrkIso") ;
                        Iso[i][1] = info->getFloat("dr03EcalIso");
                        Iso[i][2] = info->getFloat("dr03HcalD1Iso") + info->getFloat("dr03HcalD2Iso") ;
                }
                else {
                        ID[i][0] = info->getFloat("normChi2");
                        ID[i][1] = (float)info->getInt("nValMuHits");
                        ID[i][2] = (float)info->getInt("nValTrackerHits");
                        ID[i][3] = info->getFloat("dxy");
                        ID[i][4] = (float)info->getBool("muidGlobalMuonPromptTight");
                        ID[i][5] = (float)info->getBool("muidTrackerMuonArbitrated");

                        Iso[i][0] = info->getFloat("isoR03strk");
                        Iso[i][1] = info->getFloat("isoR03strkVeto");
                        Iso[i][2] = info->getFloat("isoR03sum");

                }


        }


        fillTree(Iso, ID, ZCand);

}

void
OfficialZAnalysis::fillTree(float Iso[2][3], float ID[2][6], Candidate* ZCand) {

        for(int j=0;j<2;j++) {
                for(int i=0;i<3;i++) {
                        IsoVar[j][i] = Iso[j][i];

                }
                for(int i=0;i<6;i++) {
                        IDVar[j][i] = ID[j][i];
                }
        }    

        TreeFillLepton();
        TreeFillUnderlyingEvent();

        if(SampleName.substr(0,3) == "Z_2")
                ComputeGenRecoil(ZCand);

        //Event
        Run = _e.run();
        Event = _e.eventInFile();

        AbsEvent = _e.event();
        // cout<<" sname "<<sampleName<<endl;
        strcpy (sampleName , (((string)SampleName)).c_str() );
	
	int n = (Config::dataPath).size() + SampleName.size() +1;
	TString tmpName = _e.fileName();
	strcpy (fileName, ((((string)tmpName).substr(n, ((string)tmpName).size()-n ))).c_str() );
        
	_time = _e.timestamp();


        //VetoJet
        VetoJet = VetoJetFunction();

        bool convrej =false;
        if(Zmode=="Zee")
                convrej = ConversionRejection();

        //  VetoE = SecondEcut();  //true c'ets que ça passe la coupure
        //  ConvRej = ConversionRejection(); //true c'ets que ça passe la coupure


        Float_t Bins2[14] = {0,10,20,30,40,50,60,80,100,150,200,250,300,350};
        double CorrectionCalo[13]={0.354796,0.450207,0.520878,0.560771,0.594442,0.64499,0.666981,0.708246,0.720698,0.754153,0.80916,0.756484,0.8};
        double CorrectionCT1[13]={0.35537,0.585815,0.789873,0.902999,0.959416,1.01803,1.02837,1.05108,1.025,1.02194,1.03689,0.981633,1};
        double CorrectionCT2[13]={0.837067,0.978843,1.03213,1.0504,1.06541,1.09021,1.09502,1.1033,1.06841,1.0549,1.06635,0.970004,1};
        double CorrectionTc[13]={0.604812,0.705049,0.779874,0.816408,0.833187,0.86844,0.901517,0.91991,0.922474,0.918964,0.967263,0.901705,1};
        double CorrectionPF[13]={0.743238,0.819149,0.882104,0.911021,0.914964,0.933811,0.949826,0.950231,0.926152,0.918988,0.949833,0.931791,1};
        //float correction =1;

        int bin=12;
        for(int i=0;i<13;i++) {
                if(ZCand->pt()<Bins2[i+1] && ZCand->pt()>=Bins2[i] )
                { bin = i; break;}
        }

        int METType=0;

        //With PF MET
        {
                METType=0;
                const Candidate* pfmet_ = _e.met( EventManager::kPfMet );
                TVector2 PFRecoil = ComputeRecoilWithPF(pfmet_,ZCand,1);
                TVector2 Recoil = ComputeRecoil(pfmet_,ZCand,1);
                TVector2 CTRecoil = ComputeRecoilWithTowers(pfmet_,ZCand,1);

                TreeFillMETProjection(pfmet_,METType,ZCand);
                TreeFillRecoil(PFRecoil,METType,ZCand);
                TreeFillUnderlyingJetsZ(pfmet_,METType,ZCand);
                TreeFillUnderlyingPhoton(pfmet_,METType,ZCand);
                TreeFillBasicRecoil(Recoil,METType,ZCand);
                TreeFillBisectorRecoil(Recoil,METType,ZCand);
                TreeFillCTRecoil(CTRecoil,METType,ZCand);

                Corrections[METType] = 1/CorrectionPF[bin];

                FillMETWithCorrections(pfmet_,METType,ZCand);

                CandInfo* info_=pfmet_->info();
                SumET[METType] = info_->getFloat("sumEt");

        }
        //With Track Corr Met
        {
                METType=1;
                const Candidate* trkcmet_  = _e.met( EventManager::kPatMet);
                TVector2 Recoil = ComputeRecoil(trkcmet_,ZCand,1);
                TVector2 CTRecoil = ComputeRecoilWithTowers(trkcmet_,ZCand,1);

                TreeFillMETProjection(trkcmet_,METType,ZCand);
                TreeFillRecoil(Recoil,METType,ZCand);
                TreeFillUnderlyingJetsZ(trkcmet_,METType,ZCand);
                TreeFillUnderlyingPhoton(trkcmet_,METType,ZCand);
                TreeFillBasicRecoil(Recoil,METType,ZCand); 
                TreeFillBisectorRecoil(Recoil,METType,ZCand);
                TreeFillCTRecoil(CTRecoil,METType,ZCand);

                Corrections[METType] = 1/CorrectionTc[bin];

                FillMETWithCorrections(trkcmet_,METType,ZCand);

                CandInfo* info_=trkcmet_->info();
                SumET[METType] = info_->getFloat("sumEt");

        }
        //With TypeI Calo =======
        {
                METType=3;
                _e.refreshCaloMET(); //to be sure we take the good one
                const Candidate* calomet_  = _e.met( EventManager::kCaloMet);
                TVector2 CaloRecoil = ComputeRecoilWithCalo(calomet_,ZCand,1);
                TVector2 Recoil = ComputeRecoil(calomet_,ZCand,1);
                TVector2 CTRecoil = ComputeRecoilWithTowers(calomet_,ZCand,1);

                TreeFillMETProjection(calomet_,METType,ZCand);
                TreeFillRecoil(CaloRecoil,METType,ZCand);
                TreeFillUnderlyingJetsZ(calomet_,METType,ZCand);
                TreeFillUnderlyingPhoton(calomet_,METType,ZCand);
                TreeFillBasicRecoil(Recoil,METType,ZCand);
                TreeFillBisectorRecoil(Recoil,METType,ZCand);
                TreeFillCTRecoil(CTRecoil,METType,ZCand);

                Corrections[METType] = 1/CorrectionCT1[bin];

                FillMETWithCorrections(calomet_,METType,ZCand);

                CandInfo* info_=calomet_->info();
                SumET[METType] = info_->getFloat("sumEt");

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
                TVector2 CaloRecoil = ComputeRecoilWithCalo(cmet_,ZCand,1);
                TVector2 Recoil = ComputeRecoil(cmet_,ZCand,1);
                TVector2 CTRecoil = ComputeRecoilWithTowers(cmet_,ZCand,1);

                TreeFillMETProjection(cmet_,METType,ZCand);
                TreeFillRecoil(CaloRecoil,METType,ZCand);
                TreeFillUnderlyingJetsZ(cmet_,METType,ZCand);
                TreeFillUnderlyingPhoton(cmet_,METType,ZCand);
                TreeFillBasicRecoil(Recoil,METType,ZCand);
                TreeFillBisectorRecoil(Recoil,METType,ZCand);
                TreeFillCTRecoil(CTRecoil,METType,ZCand);

                Corrections[METType] = 1/CorrectionCalo[bin];

                FillMETWithCorrections(cmet_,METType,ZCand);

                CandInfo* info_=cmet_->info();
                SumET[METType] = info_->getFloat("sumEt");

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
                TVector2 CaloRecoil = ComputeRecoilWithCalo(cmetT2_,ZCand,1);
                TVector2 Recoil = ComputeRecoil(cmetT2_,ZCand,1);
                TVector2 CTRecoil = ComputeRecoilWithTowers(cmetT2_,ZCand,1);

                TreeFillMETProjection(cmetT2_,METType,ZCand);
                TreeFillRecoil(CaloRecoil,METType,ZCand);
                TreeFillUnderlyingJetsZ(cmetT2_,METType,ZCand);
                TreeFillUnderlyingPhoton(cmetT2_,METType,ZCand);
                TreeFillBasicRecoil(Recoil,METType,ZCand);
                TreeFillBisectorRecoil(Recoil,METType,ZCand);
                TreeFillCTRecoil(CTRecoil,METType,ZCand);

                Corrections[METType] = 1/CorrectionCT2[bin];

                FillMETWithCorrections(cmetT2_,METType,ZCand);

                CandInfo* info_=cmetT2_->info();
                SumET[METType] = info_->getFloat("sumEt");
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
                TVector2 CaloRecoil = ComputeRecoilWithCalo(cmetPFT1_,ZCand,1);
                TVector2 Recoil = ComputeRecoil(cmetPFT1_,ZCand,1);
                TVector2 CTRecoil = ComputeRecoilWithTowers(cmetPFT1_,ZCand,1);

                TreeFillMETProjection(cmetPFT1_,METType,ZCand);
                TreeFillRecoil(CaloRecoil,METType,ZCand);
                TreeFillUnderlyingJetsZ(cmetPFT1_,METType,ZCand);
                TreeFillUnderlyingPhoton(cmetPFT1_,METType,ZCand);
                TreeFillBasicRecoil(Recoil,METType,ZCand);
                TreeFillBisectorRecoil(Recoil,METType,ZCand);
                TreeFillCTRecoil(CTRecoil,METType,ZCand);

                Corrections[METType] = 1/CorrectionCT2[bin];

                FillMETWithCorrections(cmetPFT1_,METType,ZCand);

                CandInfo* info_=cmetPFT1_->info();
                SumET[METType] = info_->getFloat("sumEt");

        }

        //Decompose the MET/SumET
        SeparateMETElements();
        SpecPfMETVariables();

        outTree_->Fill();
}


void
OfficialZAnalysis::TreeFillMETProjection(const Candidate* met, int metType, Candidate* ZCand) {


  //MET Projection on lepton axis
  float METpara= ZCBM_[0][0]*met->px()+ZCBM_[0][1]*met->py();
  float METperp= ZCBM_[1][0]*met->px()+ZCBM_[1][1]*met->py();

  MetProj[ metType ][0] = -METpara;
  MetProj[ metType ][1] = -METperp;

  Mets[ metType][0] = met->pt();
  Mets[ metType][1] = met->phi();
  
  //Z
  if(metType==0) {
    ZPt=ZCand->pt();
    ZMass=ZCand->mass();

    vector<Candidate*> Zcands =  RemoveLaserCor(ZCand);
    if(Zcands.size() != 0) {
      
      ZMassSC = Zcands[0]->mass();
      ZPtSC = Zcands[0]->pt();
      ZMassSCUnCor = Zcands[1]->mass();
      ZPtSCUnCor = Zcands[1]->pt();
      
      for(int id=0;id<2;id++) {
	
	SCCor[id][0] = Zcands[0]->daughter(id)->Et();
	SCCor[id][1] = Zcands[0]->daughter(id)->eta();
	SCCor[id][2] = Zcands[0]->daughter(id)->phi();
	UnCorSC[id][0] = Zcands[1]->daughter(id)->pt();
	UnCorSC[id][1] = Zcands[1]->daughter(id)->eta();
	UnCorSC[id][2] = Zcands[1]->daughter(id)->phi();
	
      }
      
    }
    else {
      ZMassSC = 0;
      ZPtSC = 0;
      ZMassSCUnCor =0;
      ZPtSCUnCor =0;
      
      for(int id=0;id<2;id++) {
	SCCor[id][0] = 0;
	SCCor[id][1] = 0;
	SCCor[id][2] = 0;
	UnCorSC[id][0] =0;
	UnCorSC[id][1] =0;
	UnCorSC[id][2] =0;
      }
      
    }

  }

}

void
OfficialZAnalysis::TreeFillRecoil(TVector2 Recoil, int metType, Candidate* ZCand){

        if(Recoil.Px()!=0 && Recoil.Py()!=0 ) {

                //Recoil angle btw 0,2 pi
                Double_t rphi = TVector2::Phi_0_2pi(Recoil.Phi() );

                //Recoil Projection
                float Rpara= ZCBM_[0][0]*Recoil.Px()+ZCBM_[0][1]*Recoil.Py();
                float Rperp= ZCBM_[1][0]*Recoil.Px()+ZCBM_[1][1]*Recoil.Py();


                Recoils[ metType][0] = Recoil.Mod();
                Recoils[ metType][1] = rphi;

                RecoilProj[ metType ][0] = Rpara;
                RecoilProj[ metType ][1] = Rperp;

                dPhiRecoilZ[ metType ] = KineUtils::dPhi(rphi,ZCand->phi() )*180/Constants::pi;

        }
        else {
                Recoils[ metType][0] = -100000;
                Recoils[ metType][1] = -100000;

                RecoilProj[ metType ][0] = -100000;
                RecoilProj[ metType ][1] = -100000;

                dPhiRecoilZ[ metType ] = -100000;
        }

}

void
OfficialZAnalysis::TreeFillBasicRecoil(TVector2 Recoil, int metType, Candidate* ZCand){

        Double_t rphi = TVector2::Phi_0_2pi(Recoil.Phi() );
        float Rpara= ZCBM_[0][0]*Recoil.Px()+ZCBM_[0][1]*Recoil.Py();
        float Rperp= ZCBM_[1][0]*Recoil.Px()+ZCBM_[1][1]*Recoil.Py();

        BasicRecoils[ metType][0] = Recoil.Mod();
        BasicRecoils[ metType][1] = rphi;

        BasicRecoilProj[ metType ][0] = Rpara;
        BasicRecoilProj[ metType ][1] = Rperp;

        BasicdPhiRecoilZ[ metType ] = KineUtils::dPhi(rphi,ZCand->phi() )*180/Constants::pi;

}

void
OfficialZAnalysis::TreeFillBisectorRecoil(TVector2 Recoil, int metType, Candidate* ZCand ) {

        //Recoil Projection in leptons frame
        TVector2 pN0 = ZCand->daughter(0)->p2().Unit();
        TVector2 pN1 = ZCand->daughter(1)->p2().Unit();
        TVector2 Plep = pN0 + pN1;

        float Leptangle =  TVector2::Phi_0_2pi(Plep.Phi() );

        float RLpara = cos(Leptangle)*Recoil.Px() + sin(Leptangle)*Recoil.Py();
        float RLperp = -sin(Leptangle)*Recoil.Px() + cos(Leptangle)*Recoil.Py();

        BisectorRecoilProj[metType][0] = RLpara;
        BisectorRecoilProj[metType][1] = RLperp;

}

void
OfficialZAnalysis::TreeFillCTRecoil(TVector2 Recoil, int metType, Candidate* ZCand){

        //Double_t rphi = TVector2::Phi_0_2pi(Recoil.Phi() );
        float Rpara= ZCBM_[0][0]*Recoil.Px()+ZCBM_[0][1]*Recoil.Py();
        float Rperp= ZCBM_[1][0]*Recoil.Px()+ZCBM_[1][1]*Recoil.Py();

        CTRecoil[ metType ][0] = Rpara;
        CTRecoil[ metType ][1] = Rperp;
}


void
OfficialZAnalysis::TreeFillUnderlyingJetsZ(const Candidate* met, int metType, Candidate* ZCand ) {

        const CandList& Jets = _e.jetList( EventManager::kPatJet);
        Candidate* theJet_;

        float pttmp=0;
        for(int unsigned ij=0;ij<Jets.size();ij++) { 
                Candidate* jet = Jets[ij];
                {
                        float dr1_ = KineUtils::dR(dau[0]->eta(), jet->eta(),
                                                   dau[0]->phi(), jet->phi() );
                        float dr2_ = KineUtils::dR(dau[1]->eta(), jet->eta(),
                                                   dau[1]->phi(), jet->phi() );

                        if(pttmp<jet->pt() && dr1_>0.2  && dr2_>0.2  )
                        { theJet_=jet; pttmp = jet->pt(); }
                }
        }
        if( Jets.size()!=0 && pttmp!=0) {
                LJetPt = theJet_->pt();
                dPhiJetZ =  KineUtils::dPhi(ZCand->phi(), theJet_->phi() )*180/Constants::pi;

                dPhiMETJet[metType] = KineUtils::dPhi(met->phi(), theJet_->phi() )*180/Constants::pi;
        }
        else {
                LJetPt = -1000;
                dPhiJetZ = -1000;
                dPhiMETJet[metType] = -1000;
        }

        //The Z
        if( Jets.size()!=0 && pttmp!=0 && metType==0) {
                dPhiJetZ = KineUtils::dPhi(ZCand->phi(), theJet_->phi() );
        }
        else
                dPhiJetZ = 0;

        dPhiMETZ[metType] = KineUtils::dPhi(met->phi(), ZCand->phi() )*180/Constants::pi;

}

void
OfficialZAnalysis::TreeFillUnderlyingPhoton(const Candidate* met, int metType, Candidate* ZCand ) {

        const CandList& Photons = _e.photons();
        Candidate* thePhoton_;

        NPhot1=0;
        NPhot2=0;

        float pttmp=0;
        for(int unsigned ij=0;ij<Photons.size();ij++) { 
                Candidate* photon = Photons[ij];
                {
                        float dr1_ = KineUtils::dR(dau[0]->eta(), photon->eta(),
                                                   dau[0]->phi(), photon->phi() );
                        float dr2_ = KineUtils::dR(dau[1]->eta(), photon->eta(),
                                                   dau[1]->phi(), photon->phi() );

                        if(pttmp<photon->pt() && dr1_>0.2 && dr2_>0.2   )
                        { thePhoton_=photon; pttmp = photon->pt(); }

                        if(dr1_>0.05 && dr1_<0.3 && photon->pt()>1)
                        {NPhot1++;}
                        if(dr2_>0.05 && dr2_<0.3 && photon->pt()>1)
                        {NPhot2++;}

                }
        }
        if( Photons.size()!=0 && pttmp!=0) {
                LPhotonPt = thePhoton_->pt();
        }
        else {
                LPhotonPt = -1000;
        }

        //The Z
        if( Photons.size()!=0 && pttmp!=0 && metType==0) {
                dPhiPhotonZ = KineUtils::dPhi(ZCand->phi(), thePhoton_->phi() );
        }
        else
                dPhiPhotonZ = 0;

}

bool OfficialZAnalysis::VetoJetFunction() {


        //Jet
        const CandList& Jets = _e.jetList( EventManager::kPatJet);
        //Candidate* theJet_;

        for(int unsigned ij=0;ij<Jets.size();ij++) { 
                Candidate* jet = Jets[ij];
                {
                        float dr1_ = KineUtils::dR(dau[0]->eta(), jet->eta(),
                                                   dau[0]->phi(), jet->phi() );
                        float dr2_ = KineUtils::dR(dau[1]->eta(), jet->eta(),
                                                   dau[1]->phi(), jet->phi() );

                        if(jet->pt()> 25 && dr1_>0.2 && dr2_>0.2  )
                                return true;
                }
        }
        return false;
}

void OfficialZAnalysis::TreeFillUnderlyingEvent(){

        //Photon

        PhotonMultiplicity=0;

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
                        }
                }
        }

        //Jet

        JetMultiplicity=0;

        const CandList& Jets = _e.jetList( EventManager::kPatJet);
        for(int unsigned ij=0;ij<Jets.size();ij++) { 
                Candidate* jet = Jets[ij];
                {
                        float dr1_ = KineUtils::dR(dau[0]->eta(), jet->eta(),
                                                   dau[0]->phi(), jet->phi() );
                        float dr2_ = KineUtils::dR(dau[1]->eta(), jet->eta(),
                                                   dau[1]->phi(), jet->phi() );
                        if(jet->pt()> 15 && dr1_>0.2 && dr2_>0.2  )
                        {
                                JetMultiplicity++;
                        }
                }
        }

}

void OfficialZAnalysis::TreeFillLepton() {

        for(int i=0;i<2;i++) {

                Candidate* dau_ = dau[i];

                Lepton[i][0] = dau_->pt();
                Lepton[i][1] = dau_->eta();
                Lepton[i][2] = dau_->phi();

                charge[i] = (int)dau_->charge();

                CandInfo* info = dau_->info();
                if(Zmode=="Zee")
                        EtSC[i] = info->getFloat("caloEnergy")/cosh(info->getFloat("caloEta"));
                else
                        EtSC[i] = dau_->pt();


                Candidate* pfCand = _e.getSecond("electron-PfCand",dau_);
                if(pfCand!=NULL)
                {
                        PFLepton[i][0] = pfCand->pt();
                        PFLepton[i][1] = pfCand->eta();
                        PFLepton[i][2] = pfCand->phi();
                }
                else
                {
                        PFLepton[i][0]=-1; PFLepton[i][1]=0; PFLepton[i][2]=0;
                }

                Candidate* SCCand(0);
                if(Zmode=="Zee")
                        SCCand= _e.getSecond("electron-SuperCluster",dau_);
                if(SCCand!=NULL)
                {
                        SCLepton[i][0] = SCCand->pt();
                        SCLepton[i][1] = SCCand->eta();
                        SCLepton[i][2] = SCCand->phi();
                }
                else {
                        SCLepton[i][0]=-1; SCLepton[i][1]=0; SCLepton[i][2]=0;
                }


        }

}


bool OfficialZAnalysis::TreeFillVertices() {

        const VtxList& vertices = _e.vertices();

        Nvertex = 0;

        Vertex* Primary; 
        bool findPV=false;
        //DZ vertices
        if(vertices.size() <2) 
                dZ = -1000;

        float dz=1000;

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



void OfficialZAnalysis::ComputeGenRecoil(Candidate* ZCand) {

        Candidate* Zmc;
        Candidate* lep1mc;
        Candidate* lep2mc;

        bool find[3]={false,false,false};

        //Find the Z
        const CandList& mc_ = _e.mcTruthCandidates();
        for(int unsigned i=0;i<mc_.size();i++) {
                Candidate* mc = mc_[i];

                if(abs(mc->pdgCode())==23)
                { Zmc = mc; find[0]=true; }

                float dR1= KineUtils::dR(mc->eta(), dau[0]->eta(),mc->phi(), dau[0]->phi() );
                float dR2= KineUtils::dR(mc->eta(), dau[1]->eta(),mc->phi(), dau[1]->phi() );

                if(abs(mc->pdgCode())==11 && dR1<0.1)
                {lep1mc = mc; find[1]=true;}
                if(abs(mc->pdgCode())==11 && dR2<0.1)
                {lep2mc = mc; find[2]=true;}

                if(find[0] && find[1] && find[2])
                        break;
        }

        /* float genZCBM_[2][2];

        //Lepton Base Matrix
        genZCBM_[0][0]=Zmc->cosPhi();
        genZCBM_[0][1]=Zmc->sinPhi();
        genZCBM_[1][0]=-Zmc->sinPhi();
        genZCBM_[1][1]=Zmc->cosPhi();

        float Rpara= genZCBM_[0][0]*Zmc->px()+genZCBM_[0][1]*Zmc->py();
        float Rperp= genZCBM_[1][0]*Zmc->px()+genZCBM_[1][1]*Zmc->py();

        GenRecoil[0][0] = -Rpara;
        GenRecoil[0][1] = -Rperp;


        float Rparar= ZCBM_[0][0]*Zmc->px()+ZCBM_[0][1]*Zmc->py();
        float Rperpr= ZCBM_[1][0]*Zmc->px()+ZCBM_[1][1]*Zmc->py();

        GenRecoil[1][0] = -Rparar;
        GenRecoil[1][1] = -Rperpr;
        */
        //W Gen Pt
        GenZPt = Zmc->pt();

        dPhiGenZ = KineUtils::dPhi(Zmc->phi() , ZCand->phi() );

}

bool
OfficialZAnalysis::SpikeCleaning() {


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
				  //                                        float swisscross = rh.SpikeVar.swCross;
				  float swisscross = 1.;

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
OfficialZAnalysis::ConversionRejection() {

        //Rejection options
        bool firstPixel[2] ={true,true};
        bool expFirstPixel[2] ={true,true};
        bool EGamConv[2] = {true,true};

        bool Conv[2][3]={ {true,true,true},
                {true,true,true} };


        for(int i=0;i<2;i++) {

                Candidate* dau_ = dau[i];

                Candidate* gsfTrack_ = _e.getSecond( "electron-GsfTrack", dau_);
                CandInfo* info_ = gsfTrack_->info();
                //possibility, use first pixel hit
                if(firstPixel[i]){
                        bool FHP_ = info_->getBool("isFPHitTrk");
                        if(!FHP_)
                                Conv[i][0] = false;
                }

                //other Possibility , using expected hits in pixel layers
                if(expFirstPixel[i]){
                        int EIP_ = info_->getInt("expInnHits");
                        int maxExpInHit =0;
                        if(EIP_ > maxExpInHit)
                                Conv[i][1] = false;
                }
                //other possibility, using EGamma rejection, not available with Summer09 datasets
                if(EGamConv[i]) {
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



void OfficialZAnalysis::FillMETWithCorrections(const Candidate* met,int metType, Candidate* ZCand) {

        float cor[2];
        for(int i=0;i<2;i++)
        {
                if(fabs(dau[i]->eta())>1.5)
                        cor[i]=0.025; //endcaps
                else
                        cor[i]=0.007; //barrel
        }

        TVector2 MetCorV2 = met->p2() - dau[0]->p2()*cor[0] -  dau[1]->p2()*cor[1];

        TVector2 ZptCor = dau[0]->p2()*(1+cor[0]) + dau[1]->p2()*(1+cor[1]);

        float ZBaseCor_[2][2];
        ZBaseCor_[0][0] = cos(ZptCor.Phi());
        ZBaseCor_[0][1] = sin(ZptCor.Phi());
        ZBaseCor_[1][0] = -sin(ZptCor.Phi());
        ZBaseCor_[1][1] = cos(ZptCor.Phi());

        float MetCpara= ZBaseCor_[0][0]*MetCorV2.Px()+ZBaseCor_[0][1]*MetCorV2.Py();
        float MetCperp= ZBaseCor_[1][0]*MetCorV2.Px()+ZBaseCor_[1][1]*MetCorV2.Py();


        //MET variables
        MetCor[metType][0] = MetCorV2.Mod();
        MetCor[metType][1] = TVector2::Phi_0_2pi(MetCorV2.Phi());
        MetCorProj[metType][0] = -MetCpara;
        MetCorProj[metType][1] = -MetCperp;

        TVector2 Recoil(0,0);
        Recoil -= MetCorV2 + dau[0]->p2()*(1+cor[0]) + dau[1]->p2()*(1+cor[1]);

        if(Recoil.Px()!=0 && Recoil.Py()!=0 ) {
                Double_t rphi = TVector2::Phi_0_2pi(Recoil.Phi() );
                float Rpara= ZBaseCor_[0][0]*Recoil.Px()+ZBaseCor_[0][1]*Recoil.Py();
                float Rperp= ZBaseCor_[1][0]*Recoil.Px()+ZBaseCor_[1][1]*Recoil.Py();

                RecoilCor[ metType][0] = Recoil.Mod();
                RecoilCor[ metType][1] = rphi;

                RecoilCorProj[ metType ][0] = Rpara;
                RecoilCorProj[ metType ][1] = Rperp;

        }

}



void OfficialZAnalysis::SeparateMETElements()
{

        //CaloTowers, raw calo MET
        {
                const CandList& caloTowers = _e.caloTowers();

                TVector2 caloMETEB(0,0);
                TVector2 caloMETEE(0,0);
                TVector2 caloMETHB(0,0);
                TVector2 caloMETHE(0,0);
                TVector2 caloMETHF(0,0);

                float caloSumEtEB = 0;
                float caloSumEtEE = 0;
                float caloSumEtHB = 0;
                float caloSumEtHE = 0;
                float caloSumEtHF = 0;

                float emEt=0;
                float hadEt=0;
                TVector2 ctP2(0,0);
                float outer=0;
                float dR[2];

                for(int unsigned i=0;i<caloTowers.size();i++) {
                        Candidate* caloTower = caloTowers[i];

                        if(caloTower->Et()<=0.3) continue;

                        for(int id=0;id<2;id++) {
                                dR[id] = KineUtils::dR(dau[id]->eta(),caloTower->eta(),dau[id]->phi(),caloTower->phi());
                        }

                        CandInfo* info = caloTower->info();

                        outer += info->getFloat("outerEt");
                        //ECAL
                        emEt = info->getFloat("emEt");
                        if(fabs(caloTower->eta())<1.5) //Barrel
                        {
                                ctP2.SetMagPhi(emEt, caloTower->phi() );
                                caloMETEB -= ctP2;
                                if(dR[0]>0.15 && dR[1]>0.15)
                                        caloSumEtEB += emEt;
                        }
                        else
                        {
                                ctP2.SetMagPhi(emEt, caloTower->phi() );
                                caloMETEE -=  ctP2;
                                if(dR[0]>0.15 && dR[1]>0.15)
                                        caloSumEtEE += emEt;
                        }

                        //HCAL
                        hadEt = info->getFloat("hadEt");
                        if(fabs(caloTower->eta())<1.25) //HB
                        {
                                ctP2.SetMagPhi(hadEt, caloTower->phi() );
                                caloMETHB -=  ctP2;
                                caloSumEtHB += hadEt;
                        }
                        else if(fabs(caloTower->eta())>2.7) //HF
                        {
                                ctP2.SetMagPhi(hadEt, caloTower->phi() );
                                caloMETHF -=  ctP2;
                                caloSumEtHF += hadEt;

                                if(fabs(caloTower->eta())>3) //containing also emEt ?
                                {
                                        ctP2.SetMagPhi(emEt, caloTower->phi() );
                                        caloMETHF -=  ctP2;
                                        caloSumEtHF += emEt;
                                }
                        }
                        else //HE
                        {
                                ctP2.SetMagPhi(hadEt, caloTower->phi() );
                                caloMETHE -=  ctP2;
                                caloSumEtHE += hadEt;
                        }
                }

                //Now Store in the Tree
                CaloSumETEB = caloSumEtEB;
                CaloSumETEE = caloSumEtEE;
                CaloSumETHB = caloSumEtHB;
                CaloSumETHE = caloSumEtHE;
                CaloSumETHF = caloSumEtHF;

                /*
                   CaloMETEB[0] = caloMETEB.Mod(); CaloMETEB[1] = caloMETEB.Phi();
                   CaloMETEE[0] = caloMETEE.Mod(); CaloMETEE[1] = caloMETEE.Phi();
                   CaloMETHB[0] = caloMETHB.Mod(); CaloMETHB[1] = caloMETHB.Phi();
                   CaloMETHE[0] = caloMETHE.Mod(); CaloMETHE[1] = caloMETHE.Phi();
                   CaloMETHF[0] = caloMETHF.Mod(); CaloMETHF[1] = caloMETHF.Phi();*/

                if(caloMETEB.Mod() != 0) {
                        CaloMETEB[0] = caloMETEB.Mod(); CaloMETEB[1] = caloMETEB.Phi(); }
                else {
                        CaloMETEB[0] = -1000; CaloMETEB[1] = -1000; }

                if(caloMETEE.Mod() != 0) {
                        CaloMETEE[0] = caloMETEE.Mod(); CaloMETEE[1] = caloMETEE.Phi(); }
                else {
                        CaloMETEE[0] = -1000; CaloMETEE[1] = -1000; }

                if(caloMETHB.Mod() != 0) {
                        CaloMETHB[0] = caloMETHB.Mod(); CaloMETHB[1] = caloMETHB.Phi(); }
                else {
                        CaloMETHB[0] = -1000; CaloMETHB[1] = -1000; }

                if(caloMETHE.Mod() != 0) {
                        CaloMETHE[0] = caloMETHE.Mod(); CaloMETHE[1] = caloMETHE.Phi(); }
                else {
                        CaloMETHE[0] = -1000; CaloMETHE[1] = -1000; }

                if(caloMETHF.Mod() != 0) {
                        CaloMETHF[0] = caloMETHF.Mod(); CaloMETHF[1] = caloMETHF.Phi(); }
                else {
                        CaloMETHF[0] = -1000; CaloMETHF[1] = -1000; }


                TVector2 caloEcalMET = caloMETEB + caloMETEE;
                TVector2 caloHcalMET = caloMETHB + caloMETHE + caloMETHF;

                if(caloEcalMET.Mod() != 0) {
                        CaloEcalMET[0] = caloEcalMET.Mod();  CaloEcalMET[1] = caloEcalMET.Phi(); }
                else {
                        CaloEcalMET[0] = -1000;  CaloEcalMET[1] = -1000; }

                if(caloHcalMET.Mod() != 0) {
                        CaloHcalMET[0] = caloHcalMET.Mod();  CaloHcalMET[1] = caloHcalMET.Phi(); }
                else {
                        CaloHcalMET[0] = -1000;  CaloHcalMET[1] = -1000; }

        }
        TVector2 METpf(0,0); 
        //Pf candidates, pf MET
        {
                TVector2 pfMETEB(0,0);
                TVector2 pfMETEE(0,0);
                TVector2 pfMETEEp(0,0);
                TVector2 pfMETEEm(0,0);
                TVector2 pfMETHB(0,0);
                TVector2 pfMETHE(0,0);
                TVector2 pfMETHF(0,0);

                float pfSumEtEB = 0;
                float pfSumEtEE = 0;
                float pfSumEtHB = 0;
                float pfSumEtHE = 0;
                float pfSumEtHF = 0;
                TVector2 pfP2(0,0);

                float emE=0,hadE=0,emEt=0,hadEt=0;

                //Get the PF electrons
                CandList pfdau; // cout<<" high pf cand "<<endl;
                for(int i=0;i<2;i++) {
                        Candidate* dau_ = dau[i];
                        Candidate* pfCand = _e.getSecond("electron-PfCand",dau_);
                        pfdau.push_back(pfCand);
                }

                TVector2 METpftmp(0,0);
                const CandList& pfCands = _e.pfCandidates();

                float dR[2]={1000,1000};

                for(int unsigned i=0;i<pfCands.size();i++) {
                        Candidate* pfCand = pfCands[i];
                        CandInfo* info = pfCand->info();

                        if(pfCand->Et()<=0) continue;

                        for(int id=0;id<2;id++) {
                                if(pfdau[id] != 0) {
                                        if( abs(pfCand->pdgCode()) == abs(pfdau[id]->pdgCode()) )
                                                dR[id] = KineUtils::dR(pfCand->eta(),pfdau[id]->eta(),pfCand->phi(),pfdau[id]->phi());
                                        else if(fabs(pfCand->pt()-pfdau[id]->pt())/pfCand->pt() <0.05)
                                                dR[id] = KineUtils::dR(pfCand->eta(),pfdau[id]->eta(),pfCand->phi(),pfdau[id]->phi());
                                        else
                                                dR[id] = 1000;
                                }
                        }

                        METpftmp.SetMagPhi(pfCand->pt(),pfCand->phi());
                        METpf -= METpftmp;
                        //Compute Et
                        emE = info->getFloat("ecalE");

                        if(info->getFloat("pS1E")>=0)
                                emE += info->getFloat("pS1E");
                        if(info->getFloat("pS2E")>=0)
                                emE += info->getFloat("pS2E");

                        hadE = info->getFloat("hcalE");
                        emEt = KineUtils::et(emE, pfCand->eta() );
                        hadEt = KineUtils::et(hadE, pfCand->eta() );
                        //=========
                        //      if(emEt<0)
                        //	cout<<" EM < 0  : "<<emE<<"  -> "<<info->getFloat("ecalE")<<"   "<<pfCand->eta()<<"   "<<pfCand->uid()<<"   "<<pfCand->pt()<<"   "<<pfCand->pdgCode()<<endl;
                        //ECAL
                        if(fabs(pfCand->eta())<1.5) //Barrel
                        {
                                pfP2.SetMagPhi(emEt, pfCand->phi() );
                                pfMETEB -= pfP2;
                                if(dR[0]>0.01 && dR[1]>0.01)
                                        pfSumEtEB += emEt;
                        }
                        else
                        {
                                pfP2.SetMagPhi(emEt, pfCand->phi() );
                                pfMETEE -=  pfP2;
                                if(dR[0]>0.01 && dR[1]>0.01)
                                        pfSumEtEE += emEt;

                                if(pfCand->eta()>0)
                                        pfMETEEp -= pfP2;
                                else
                                        pfMETEEm -= pfP2;

                        }

                        //HCAL
                        if(fabs(pfCand->eta())<1.25) //HB
                        {
                                pfP2.SetMagPhi(hadEt, pfCand->phi() );
                                pfMETHB -=  pfP2;
                                pfSumEtHB += hadEt;
                        }
                        else if(fabs(pfCand->eta())>2.7) //HF
                        {
                                pfP2.SetMagPhi(hadEt, pfCand->phi() );
                                pfMETHF -=  pfP2;
                                pfSumEtHF += hadEt;

                                if(fabs(pfCand->eta())>3) //containing also emEt ?
                                {
                                        pfP2.SetMagPhi(emEt, pfCand->phi() );
                                        pfMETHF -=  pfP2;
                                        pfSumEtHF += emEt;
                                }
                        }
                        else //HE
                        {
                                pfP2.SetMagPhi(hadEt, pfCand->phi() );
                                pfMETHE -=  pfP2;
                                pfSumEtHE += hadEt;
                        }
                }


                //Second passage with low pt pfCand
                // cout<<" low pf cand "<<endl;
                _e.refreshPfCandidates(); //to be sure we take the good one
                int v=1; //pfCand_lowPt
                _e.reloadPfCandidates( v );
                const CandList& pfCandsLowPt = _e.candList( EventManager::kPfCand_lowPt );

                for(int unsigned i=0;i<pfCandsLowPt.size();i++) {
                        Candidate* pfCand = pfCandsLowPt[i];
                        CandInfo* info = pfCand->info();

                        if(pfCand->Et()<=0) continue;
                        METpftmp.SetMagPhi(pfCand->pt(),pfCand->phi());
                        METpf -= METpftmp;
                        //Compute Et
                        emE = info->getFloat("ecalE");

                        if(info->getFloat("pS1E")>=0)
                                emE += info->getFloat("pS1E");
                        if(info->getFloat("pS2E")>=0)
                                emE += info->getFloat("pS2E");

                        hadE = info->getFloat("hcalE");
                        emEt = KineUtils::et(emE, pfCand->eta() );
                        hadEt = KineUtils::et(hadE, pfCand->eta() );
                        //=========
                        if(emEt<0 )
                                cout<<" EM < 0  : "<<emE<<"   "<<pfCand->eta()<<"   "<<pfCand->uid()<<"   "<<pfCand->pt()<<"   "<<pfCand->pdgCode()<<endl;
                        //ECAL
                        if(fabs(pfCand->eta())<1.5) //Barrel
                        {
                                pfP2.SetMagPhi(emEt, pfCand->phi() );
                                pfMETEB -= pfP2;
                                pfSumEtEB += emEt;
                        }
                        else
                        {
                                pfP2.SetMagPhi(emEt, pfCand->phi() );
                                pfMETEE -=  pfP2;
                                pfSumEtEE += emEt;
                        }

                        //HCAL
                        if(fabs(pfCand->eta())<1.25) //HB
                        {
                                pfP2.SetMagPhi(hadEt, pfCand->phi() );
                                pfMETHB -=  pfP2;
                                pfSumEtHB += hadEt;
                        }
                        else if(fabs(pfCand->eta())>2.7) //HF
                        {
                                pfP2.SetMagPhi(hadEt, pfCand->phi() );
                                pfMETHF -=  pfP2;
                                pfSumEtHF += hadEt;

                                if(fabs(pfCand->eta())>3) //containing also emEt ?
                                {
                                        pfP2.SetMagPhi(emEt, pfCand->phi() );
                                        pfMETHF -=  pfP2;
                                        pfSumEtHF += emEt;
                                }
                        }
                        else //HE
                        {
                                pfP2.SetMagPhi(hadEt, pfCand->phi() );
                                pfMETHE -=  pfP2;
                                pfSumEtHE += hadEt;
                        }
                }

                //Now Store in the Tree
                PfSumETEB = pfSumEtEB;
                PfSumETEE = pfSumEtEE;
                PfSumETHB = pfSumEtHB;
                PfSumETHE = pfSumEtHE;
                PfSumETHF = pfSumEtHF;

                if(pfMETEB.Mod() != 0) {
                        PfMETEB[0] = pfMETEB.Mod(); PfMETEB[1] = pfMETEB.Phi(); }
                else {
                        PfMETEB[0] = -1000; PfMETEB[1] = -1000; }

                if(pfMETEE.Mod() != 0) {
                        PfMETEE[0] = pfMETEE.Mod(); PfMETEE[1] = pfMETEE.Phi(); }
                else {
                        PfMETEE[0] = -1000; PfMETEE[1] = -1000; }

                if(pfMETEEp.Mod() != 0) {
                        PfMETEEx[0] = pfMETEEp.Phi(); }
                else {
                        PfMETEEx[0] = -1000; }
                if(pfMETEEm.Mod() != 0) {
                        PfMETEEx[1] = pfMETEEm.Phi(); }
                else {
                        PfMETEEx[1] = -1000; }

                if(pfMETHB.Mod() != 0) {
                        PfMETHB[0] = pfMETHB.Mod(); PfMETHB[1] = pfMETHB.Phi(); }
                else {
                        PfMETHB[0] = -1000; PfMETHB[1] = -1000; }

                if(pfMETHE.Mod() != 0) {
                        PfMETHE[0] = pfMETHE.Mod(); PfMETHE[1] = pfMETHE.Phi(); }
                else {
                        PfMETHE[0] = -1000; PfMETHE[1] = -1000; }

                if(pfMETHF.Mod() != 0) {
                        PfMETHF[0] = pfMETHF.Mod(); PfMETHF[1] = pfMETHF.Phi(); }
                else {
                        PfMETHF[0] = -1000; PfMETHF[1] = -1000; }


                TVector2 pfEcalMET = pfMETEB + pfMETEE;
                TVector2 pfHcalMET = pfMETHB + pfMETHE + pfMETHF;

                if(pfEcalMET.Mod() != 0) {
                        PfEcalMET[0] = pfEcalMET.Mod();  PfEcalMET[1] = pfEcalMET.Phi(); }
                else {
                        PfEcalMET[0] = -1000;  PfEcalMET[1] = -1000; }

                if(pfHcalMET.Mod() != 0) {
                        PfHcalMET[0] = pfHcalMET.Mod();  PfHcalMET[1] = pfHcalMET.Phi(); }
                else {
                        PfHcalMET[0] = -1000;  PfHcalMET[1] = -1000; }


        }


        //  cout<<" test "<<METpf.Mod()<<"   "<<Mets[0][0]<<"  "<<METpf.Phi()<<"   "<<Mets[0][1]<<endl;

        //cout<<" test calo "<<METcalo.Mod()<<"   "<<Mets[3][0]<<"  "<<METcalo.Phi()<<"   "<<Mets[3][1]<<endl;





}



void OfficialZAnalysis::SpecPfMETVariables() {


        const Candidate* pfmet_ = _e.met( EventManager::kPfMet );
        CandInfo* info = pfmet_->info();


        cEmFrac = info->getFloat("cEmEF");
        nEmFrac = info->getFloat("nEmEF");
        cHadFrac = info->getFloat("cHadEF");
        nHadFrac = info->getFloat("nHadEF");
        cMuFrac = info->getFloat("cMuEF");



}




vector<Candidate*>
OfficialZAnalysis::RemoveLaserCor(Candidate* ZCand) {
  
  //const CandAssoc& e_Sc_ =  _e.candAssoc("electron-SuperCluster");
  const CandMap& Sc_Bc_ = _e.candMap("SuperCluster-BasicCluster");
  const CandIdRecHitMap& Bc_rh_ = _e.ecalRecHitMap("BasicCluster-EcalRecHit");

  Candidate* scs[2];
  Candidate* scsuncor[2];
  vector<Candidate*> listZ;

  //looping on the two electrons
  for( int idau=0; idau<2; idau++ ) {

    float energy=0;
    float uncorenergy = 0;

    Candidate* electron_ = ZCand->daughter(idau)->theBase();
    Candidate* superCluster = _e.getSecond("electron-SuperCluster", electron_ );
  
    if(superCluster==NULL) return listZ;

    scs[idau] = superCluster;

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
	    energy +=  (*itrh).second.first.e;
	    uncorenergy +=  ((*itrh).second.first.e)/(((*itrh).second.first).LaserVar.laserCorr);
	  }
      }

    CandInfo* info= superCluster->info();
    float brem = info->getFloat("phiWidth") / info->getFloat("etaWidth") ;

    if( fabs(superCluster->eta())>1.479 ) //endcap
      {
	energy += info->getFloat("preshowerEnergy");
	uncorenergy += info->getFloat("preshowerEnergy");
	energy = (float)(electron_br1EE( (double)brem , (double)energy));
	energy =  KineUtils::e( (float)electron_br1_completeEE( KineUtils::et(energy, superCluster->eta()), superCluster->eta() ), superCluster->eta() );
	uncorenergy = (float)(electron_br1EE( (double)brem , (double)uncorenergy));
	uncorenergy =  KineUtils::e( (float)electron_br1_completeEE( KineUtils::et(uncorenergy, superCluster->eta()), superCluster->eta() ), superCluster->eta() );
      }
    else {
      energy = energy/f5x5( (double)((int)(fabs(superCluster->eta())*(5/0.087))));
	energy = (float)(electron_br1EB( (double)brem , (double)energy));
	energy =  KineUtils::e( (float)electron_br1_completeEB( KineUtils::et(energy, superCluster->eta()), superCluster->eta() ), superCluster->eta() );
	uncorenergy = (float)(electron_br1EB( (double)brem , (double)uncorenergy));
	uncorenergy =  KineUtils::e( (float)electron_br1_completeEB( KineUtils::et(uncorenergy, superCluster->eta()), superCluster->eta() ), superCluster->eta() );
    }


    TVector3 p3_(0, 0, 0);
    p3_.SetPtEtaPhi( KineUtils::et(uncorenergy, superCluster->eta()), superCluster->eta(), superCluster->phi() );
    scsuncor[idau] = Candidate::create( p3_, 22, const_cast<Vertex*>(&(_e.primaryVertex())) ); 

    cout<<" electron_   "<<superCluster->E()<<"    "<<energy<<"    "<<uncorenergy<<"    "<<superCluster->eta()<<endl;
  }
  
  //make the candidate with laser corrections
  Candidate* ZCandCor = Candidate::create(scs[0],scs[1]);

  //And now without
  Candidate* ZCandUnCor = Candidate::create(scsuncor[0],scsuncor[1]);

  listZ.push_back(ZCandCor);
  listZ.push_back(ZCandUnCor);
  
  return listZ;
}




//Default clustering corrections

double 
OfficialZAnalysis::electron_br1EB(double brLinear, double e) {
    //
    // first parabola (for brLinear < threshold)
    // p0*x^2 + p1*x + p2
    // second parabola (for brLinear >= threshold)
    // ax^2 + bx + c, make y and y' the same in threshold
    // y = p0*threshold^2 + p1*threshold + p2
    // yprime = p1 + 2*p0*threshold
    // a = p3
    // b = yprime - 2*a*threshold
    // c = y - a*threshold^2 - b*threshold
    // final result is multiplied by cos(p5/br + p6*br + p7)

    // make NO correction if brLinear is invalid!
    if ( brLinear == 0 ) return e;
    //

    // this is for 31X
    if ( brLinear < 1.1 ) brLinear = 1.1;
    if ( brLinear > 8 ) brLinear = 8.0;

    double p0 = -0.05185;
    double p1 = 0.1354;
    double p2 = 0.9165;
    double p3 = -0.0005626;
    double p4 = 1.385;    

    double threshold = p4;
    
    double y = p0*threshold*threshold + p1*threshold + p2;
    double yprime = 2*p0*threshold + p1;
    double a = p3;
    double b = yprime - 2*a*threshold;
    double c = y - a*threshold*threshold - b*threshold;

    double fCorr = 1;
    if ( brLinear < threshold ) 
      fCorr = p0*brLinear*brLinear + p1*brLinear + p2;
    else 
      fCorr = a*brLinear*brLinear + b*brLinear + c;
    
    return e/fCorr;

}

double
OfficialZAnalysis::electron_br1_completeEB(double et, double eta) {
    double fCorr = 0;

    // hybrid SC 
    // fBrem
    double c0 = 1.002;
    double c1 = -0.7424;
    double c2 = 0;

    // fEtEta
    double c4 = 0.5558;
    double c5 = 2.375; 
    double c6 = 0.1869;
 
    // final fitting
    double c7 = 1.081;  // curve point in eta distribution
    double c8 = 7.6;     // sharpness of the curve
    double c3 = -0.00181;


    double p0 = c0 + c1/(et + c2);
    double p1 = c4/(et + c5) + c6/(et*et);

    fCorr = p0 + p1*atan(c8*(c7 - fabs(eta))) + c3*fabs(eta);

    return et/fCorr;
}




double
OfficialZAnalysis::electron_br1EE(double brLinear, double e) {
    // These corrections are tuned on CMSSW 18x dynamicHybridSuperClusters!
    // YM: 02/05/2008
    //
    // first parabola (for brLinear < threshold)
    // p0*x^2 + p1*x + p2
    // second parabola (for brLinear >= threshold)
    // ax^2 + bx + c, make y and y' the same in threshold
    // y = p0*threshold^2 + p1*threshold + p2
    // yprime = p1 + 2*p0*threshold
    // a = p3
    // b = yprime - 2*a*threshold
    // c = y - a*threshold^2 - b*threshold
    // final result is multiplied by cos(p5/br + p6*br + p7)

    // make NO correction if brLinear is invalid!
    if ( brLinear == 0 ) return e;
    // make a flat correction if brLinear is too big (>9)

    // FM with preshower
    if ( brLinear > 6.5 ) brLinear = 6.5;
    if ( brLinear < 0.9 ) brLinear = 0.9;

    // ============= Fixed Matrix With Preshower SC
    double p0 = -0.1214;   //-0.04163;
    double p1 = 0.2362;    //0.08551;
    double p2 = 0.8847;     //0.95048;
    double p3 = -0.00193; //-0.002308;
    double p4 = 1.057;      //1.077;

    double threshold = p4;
    
    double y = p0*threshold*threshold + p1*threshold + p2;
    double yprime = 2*p0*threshold + p1;
    double a = p3;
    double b = yprime - 2*a*threshold;
    double c = y - a*threshold*threshold - b*threshold;

    double fCorr = 1;
    if ( brLinear < threshold ) 
      fCorr = p0*brLinear*brLinear + p1*brLinear + p2;
    else 
      fCorr = a*brLinear*brLinear + b*brLinear + c;
    
    return e/fCorr;

  }

double
OfficialZAnalysis::electron_br1_completeEE(double et, double eta) {
    double fCorr = 0.;

    double c0 = 2.213;
    double c1 = -17.29;

    double c2 = -0.599;
    double c3 = 8.874;

    double c4 = 0.09632;
    double c5 = -1.457;

    double c6 = -0.7584;
    double c7 = 10.29;

    double p0 = c0 + c1/sqrt(et);
    double p1 = c2 + c3/sqrt(et);
    double p2 = c4 + c5/sqrt(et);
    double p3 = c6 + c7/sqrt(et);

    fCorr = p0 + p1*fabs(eta) + p2*eta*eta + p3/fabs(eta);

    return et/fCorr;
}


double
OfficialZAnalysis::f5x5( double iEta ) {
    if ( iEta < 40.2198 ) return 1;
    return 1 - 3.03103e-6*(iEta - 40.2198)*(iEta - 40.2198);
  }
