#include "TagAndProbeAnalysis.hh"

#include <algorithm>
#include <cassert>
#include <sstream>

using namespace std;

#include "TH1F.h"
#include "TH2F.h"
#include "TH1I.h"
#include "TH2I.h"
#include "TProfile.h"
#include "TGraphErrors.h"

#include "MECore/src/ME.hh"
#include "MECore/src/MEEBDisplay.hh"
#include "MECore/src/MEEEDisplay.hh"
#include "Analysis/tools/PrintTree.hh"
#include "Analysis/selectors/TagAndProbeElectronSelector.hh"
#include "Analysis/utils/StatUtils.hh"

ClassImp( TagAndProbeAnalysis )

TagAndProbeAnalysis::TagAndProbeAnalysis( Sample& sample , const std::string & collectionFileName) : 
  SampleAnalysis( "TAPAnalysis", sample, collectionFileName ),
  outTreeID(0),outTreeIso(0)
 
{

  SampleName = sample.name();


  npt=400;
  neta=42;
  
  ptStep = 5;
  etaStep = 0.125;
  
  vector<vector<int> > tmp(npt,vector<int>(neta,0));
  Ntt_Iso=tmp;
  Ntp_Iso=tmp;
  Ntf_Iso=tmp;

  Ntt_ID=tmp;
  Ntp_ID=tmp;
  Ntf_ID=tmp;
  
  Ntt_Reco=tmp;
  Ntp_Reco=tmp;
  Ntf_Reco=tmp;


  SampleAnalysis::defineTemplate("EfficiencyPt", npt,0,npt*ptStep);
  SampleAnalysis::defineTemplate("EfficiencyEta", neta,-2.5 -etaStep,neta*etaStep-2.5 -etaStep);

  outTreeID  = new TTree("TAPVariablesID", "Variables for tag and probe ID");
  outTreeIso = new TTree("TAPVariablesIso", "Variables for tag and probe Iso");
  outTreeReco = new TTree("TAPVariablesReco", "Variables for tag and probe Reco");

}


TagAndProbeAnalysis::~TagAndProbeAnalysis()
{

}




bool
TagAndProbeAnalysis::analyzeEvent()
{

  if((_ievt+1)%500 == 0)
    cout<<"========== Event "<<_ievt+1<<endl;

  bool acc = _e.inAcceptance();
  if(!acc) return false;

  if( !HLTFiring())
    return false;

  GetNvertex();
  GetNJet();

  RecoEfficiency();
  IsoEfficiency();
  IDEfficiency();

  return true;

}


void
TagAndProbeAnalysis::bookHistograms()
{

  outTreeID->SetDirectory( _outputFile->GetDirectory(NULL) );
  outTreeIso->SetDirectory( _outputFile->GetDirectory(NULL) );
  initRootTree();



}

void
TagAndProbeAnalysis::writeHistograms()
{
  EndReconstructionEfficiency();
  EndIsolationEfficiency();
  EndIdentificationEfficiency();


  outTreeID->Write();
  outTreeIso->Write();
  outTreeReco->Write();

}

bool
TagAndProbeAnalysis::HLTFiring() {
  
  bool fired=false;
  
  int run = _e.run();
  
  if(SampleName.substr(0,1)=="E") { //FIXME MC...
    
    if( ( (run <= 140401) && _e.isFired("HLT_Ele15_LW_L1R") ) ||
	( (run >= 140402 && run <= 143962 ) && _e.isFired("HLT_Ele15_SW_L1R") ) ||
	( (run >= 143963 && run <= 144114 ) && _e.isFired("HLT_Ele15_SW_CaloEleId_L1R") ) ||
	( (run >= 144115 && run <= 147116 ) && _e.isFired("HLT_Ele17_SW_CaloEleId_L1R") ) ||
	( (run >= 147117 && run <= 148058 ) && _e.isFired("HLT_Ele17_SW_TightEleId_L1R") ) ||
	( (run >= 148059 && run <= 149064 ) && _e.isFired("HLT_Ele17_SW_TighterEleIdIsol_L1R_v2") ) ||
	( (run >= 149065 && run <= 149442 ) && _e.isFired("HLT_Ele17_SW_TighterEleIdIsol_L1R_v3") ) )
      fired =true;
    
    return fired;
  }
  else
    return true;

}


void
TagAndProbeAnalysis::GetNvertex() {

 const VtxList& vertices = _e.vertices();
 
 _NVertex = 0;
 
 for(int unsigned iv=0;iv<vertices.size();iv++) {
    
    Vertex* Vert = vertices[iv];
    CandInfo* infoV = Vert->info();
   
    bool goodV=false;
    if(fabs((Vert->pos()).Z())<24 && (Vert->pos()).Perp()<2 && infoV->getFloat("ndof")> 4 && !(infoV->getBool("isFake")) )
      goodV=true;
    
    if(goodV) {
      _NVertex++;
    }
  }
  
}


void
TagAndProbeAnalysis::GetNJet() {

  _NJet=0;

  const CandList& electrons = _e.electrons();
  const CandList& Jets = _e.jetList( EventManager::kPatJet);
  for(int unsigned ij=0;ij<Jets.size();ij++) { 
    Candidate* jet = Jets[ij];
    
    bool goodjet=true;

    if(jet->pt()> 25)
      {
	for(int unsigned ie=0;ie<electrons.size();ie++) {
	  Candidate* electron = electrons[ie];
	  
	  if(electron->pt()> 20) {
	    
	    float dr = KineUtils::dR(jet->eta(), electron->eta(),
				     jet->phi(), electron->phi() );
	    
	    if(dr<0.2) {
	      goodjet=false;
	    }
	    
	  }
	}
	
	if(goodjet)
	  _NJet++;
	
      }
  }
  
}


void
TagAndProbeAnalysis::IsoEfficiency() {

  //CandList Tag and Probe
  CandList _TagCandList;
  CandList _ProbeCandList;

  vector<int> _TagList;
  vector<int> _ProbeList;

  //Maps containing the different cuts 
  CutMap _TagCuts;
  CutMap _ProbeCuts;

  const CandList& electrons = _e.electrons();

  //Tag definition ======================================
  
  // Presel/acceptance
  vector<string> pteta;
  pteta.push_back( "20." );
  pteta.push_back( "fiducial" );


  // Identification
   vector<string> IDs;
   IDs.push_back("eid80" );
  

  vector<string> Isos;
  Isos.push_back("iso80");
 

  AddCondition( _TagCuts, "PtEta", pteta );
  AddCondition( _TagCuts, "ID", IDs );
  AddCondition( _TagCuts, "Iso", Isos );
  


  TagAndProbeElectronSelector TagSelector(EventManager::kElectron, _TagCuts );

  
  //Probe definition ======================================
  
  AddCondition( _ProbeCuts, "PtEta", pteta );
  AddCondition( _ProbeCuts, "ID", IDs );

  TagAndProbeElectronSelector ProbeSelector(EventManager::kElectron, _ProbeCuts );
  // cout<<endl<<" electron size  "<<electrons.size()<<endl;
   for(int unsigned i=0; i<electrons.size(); i++) {
     bool tag=false;  bool probe=false;
     // cout<<" el number "<<i<<endl;
     tag=TagSelector.accept( (*electrons[i]) );
     probe=ProbeSelector.accept( (*electrons[i]) );
     
     if(tag ) {
       
       //    cout<<"-----> electron tagged"<<endl;
       _TagCandList.push_back( (electrons[i]) );
       _TagList.push_back( (electrons[i])->uid() );
     }
     if(probe ) {
       //   cout<<"-----> electron probbed"<<endl;
       _ProbeCandList.push_back( (electrons[i]) );
       _ProbeList.push_back( (electrons[i])->uid() );
     }
     
   }

   EfficiencyMeasurement(electrons, electrons, _TagList, _ProbeList, "Iso", Isos, Ntt_Iso, Ntf_Iso, Ntp_Iso);

}




void
TagAndProbeAnalysis::RecoEfficiency() {

  //CandList Tag and Probe
  CandList _TagCandList;
  CandList _ProbeCandList;

  vector<int> _TagList;
  vector<int> _ProbeList;

  //Maps containing the different cuts 
  CutMap _TagCuts;
 
  const CandList& electrons = _e.electrons();

  //Tag definition ======================================
  
  // Presel/acceptance
  vector<string> pteta;
  pteta.push_back( "20" );
  pteta.push_back( "fiducial" );


  // Identification
   vector<string> IDs;
   IDs.push_back("eid80" );
  

  vector<string> Isos;
  Isos.push_back("iso80");
 

  AddCondition( _TagCuts, "PtEta", pteta );
  AddCondition( _TagCuts, "ID", IDs );
  AddCondition( _TagCuts, "Iso", Isos );
  
 
 

  TagAndProbeElectronSelector TagSelector(EventManager::kElectron, _TagCuts );

  
  //Probe definition ======================================

  CutMap _ProbeCutsEB;

  vector<string> hoeEB;
  hoeEB.push_back("hadronicOverEm");
  hoeEB.push_back("0.15");
  hoeEB.push_back("f");
  hoeEB.push_back("l");

  vector<string> sieieEB;
  sieieEB.push_back("sigmaIetaIeta");
  sieieEB.push_back("0.01");
  sieieEB.push_back("f");
  sieieEB.push_back("l");
  
  AddCondition( _ProbeCutsEB, "PtEta", pteta );
  AddCondition( _ProbeCutsEB, "Misc", hoeEB );
  AddCondition( _ProbeCutsEB, "Misc", sieieEB );


  TagAndProbeElectronSelector ProbeSelectorEB(EventManager::kElectron, _ProbeCutsEB );

  CutMap _ProbeCutsEE;

  vector<string> hoeEE;
  hoeEE.push_back("hadronicOverEm");
  hoeEE.push_back("0.07");
  hoeEE.push_back("f");
  hoeEE.push_back("l");

  vector<string> sieieEE;
  sieieEE.push_back("sigmaIetaIeta");
  sieieEE.push_back("0.03");
  sieieEE.push_back("f");
  sieieEE.push_back("l");
  
  AddCondition( _ProbeCutsEE, "PtEta", pteta );
  AddCondition( _ProbeCutsEE, "Misc", hoeEE );
  AddCondition( _ProbeCutsEE, "Misc", sieieEE );

  TagAndProbeElectronSelector ProbeSelectorEE(EventManager::kElectron, _ProbeCutsEE );
  
   for(int unsigned i=0; i<electrons.size(); i++) {
     bool tag=false; 
     
     tag=TagSelector.accept( (*electrons[i]) );
     
     if(tag ) {
       
       _TagCandList.push_back( (electrons[i]) );
       _TagList.push_back( (electrons[i])->uid() );
     }

   }
   
   const CandList& photons = _e.photons();
   for(int unsigned i=0; i<photons.size(); i++) {
     bool probe=false; 

     if( fabs(photons[i]->eta()) < 1.479)
       probe=ProbeSelectorEB.accept( (*photons[i]) );
     else
       probe=ProbeSelectorEE.accept( (*photons[i]) );

     if(probe ) {
       
       _ProbeCandList.push_back( (photons[i]) );
       _ProbeList.push_back( (photons[i])->uid() );
     }
     
   }

   vector<string> match;
   match.push_back("0.3");
   match.push_back("l");

   EfficiencyMeasurement(electrons, photons, _TagList, _ProbeList, "Reco", match, Ntt_Reco, Ntf_Reco, Ntp_Reco);

}


void
TagAndProbeAnalysis::IDEfficiency() {

  //CandList Tag and Probe
  CandList _TagCandList;
  CandList _ProbeCandList;

  vector<int> _TagList;
  vector<int> _ProbeList;

  //Maps containing the different cuts
  CutMap _TagCuts;
  CutMap _ProbeCuts;

  const CandList& electrons = _e.electrons();

  //Tag definition ======================================
  
  // Presel/acceptance
  vector<string> pteta;
  pteta.push_back( "20." );
  pteta.push_back( "fiducial" );
 

  // Identification
   vector<string> IDs;
   IDs.push_back("eid80" );
 

  vector<string> Isos;
  Isos.push_back("iso80");
 

  AddCondition( _TagCuts, "PtEta", pteta );
  AddCondition( _TagCuts, "ID", IDs );
  AddCondition( _TagCuts, "Iso", Isos );
  
  AddCondition( _ProbeCuts, "PtEta", pteta );
 

  TagAndProbeElectronSelector TagSelector(EventManager::kElectron, _TagCuts );

  
  //Probe definition ======================================
  
  vector<string> IDsprobe;
  IDsprobe.push_back("eid80");

  TagAndProbeElectronSelector ProbeSelector(EventManager::kElectron, _ProbeCuts );
  // cout<<endl<<" electron size  "<<electrons.size()<<endl;
   for(int unsigned i=0; i<electrons.size(); i++) {
     //     CandInfo* info = electrons[i]->info();

     bool tag=false;  bool probe=false;
     //cout<<" el number "<<i<<endl;
     tag=TagSelector.accept( (*electrons[i]) );
     probe=ProbeSelector.accept( (*electrons[i]) );
     
     if(tag ) {
       /*  
       cout<<"-----> electron tagged"
	   <<"   "<<info->getBool("eid80")
	   <<"    "<<info->getBool("iso80")<<endl;*/
       _TagCandList.push_back( (electrons[i]) );
       _TagList.push_back( (electrons[i])->uid() );
     }
     if(probe ) {

       /*   cout<<"-----> electron probed"
	   <<"   "<<info->getBool("eid80")
	   <<"    "<<info->getBool("iso80")<<endl;
       */
       //   cout<<"-----> electron probbed"<<endl;
       _ProbeCandList.push_back( (electrons[i]) );
       _ProbeList.push_back( (electrons[i])->uid() );
     }
     
   }

   EfficiencyMeasurement(electrons, electrons, _TagList, _ProbeList, "ID" , IDsprobe , Ntt_ID, Ntf_ID, Ntp_ID);

}



void TagAndProbeAnalysis::EndIsolationEfficiency() {
 
  vector<vector<float> > effeta(2,vector<float>(neta,0));
  vector<vector<float> > effpt(2,vector<float>(npt,0));
  
  float Num=0, Denom=0;

  int Ntt=0, Ntf=0, Ntp=0;

  for(int ipt=0;ipt<npt;ipt++) {
    for(int ieta=0;ieta<neta;ieta++) {
      effpt[0][ipt] += (float)(Ntt_Iso[ipt][ieta]+ Ntp_Iso[ipt][ieta]);
      effpt[1][ipt] += (float)(Ntt_Iso[ipt][ieta]+ Ntp_Iso[ipt][ieta]+ Ntf_Iso[ipt][ieta]);

      effeta[0][ieta] += (float)(Ntt_Iso[ipt][ieta]+ Ntp_Iso[ipt][ieta]);
      effeta[1][ieta] += (float)(Ntt_Iso[ipt][ieta]+ Ntp_Iso[ipt][ieta]+ Ntf_Iso[ipt][ieta]);
    
      Num += (float)(Ntt_Iso[ipt][ieta]+ Ntp_Iso[ipt][ieta]);
      Denom += (float)(Ntt_Iso[ipt][ieta]+ Ntp_Iso[ipt][ieta]+ Ntf_Iso[ipt][ieta]);
  
      Ntt += Ntt_Iso[ipt][ieta];
      Ntf += Ntf_Iso[ipt][ieta];
      Ntp += Ntp_Iso[ipt][ieta];

    }
  }
  cout<<" Fin calcul nombres isolation "<<endl;
  cout<<"  Ntagtag   "<<Ntt<<endl;
  cout<<"  Ntagpass   "<<Ntp<<endl;
  cout<<"  Ntagfail   "<<Ntf<<endl;
  cout<<" numerateur : "<<Num<<"   denominateur :  "<<Denom<<endl;
  cout<<" Eff : "<<Num/Denom<<" +- "<<StatUtils::EBinom(Num/Denom, (int)Denom)<<endl;
  _outputFile->cd();
  TGraphErrors EPt(npt);
  EPt.SetName("IsoPt"); 
  TGraphErrors EEta(neta);
  EEta.SetName("IsoEta");

  for(int ipt=0;ipt<npt;ipt++) {
    if(ipt<20)
      cout<<ipt*ptStep + ptStep/2.<<" Nevt "<<effpt[0][ipt]<<" : "<<effpt[1][ipt]<<endl;
    if(effpt[1][ipt]!=0) {
      if(ipt<20)
	cout<<ipt*ptStep + ptStep/2.<<"  "<<effpt[0][ipt]/effpt[1][ipt]<<" +-  "<<StatUtils::EBinom(effpt[0][ipt]/effpt[1][ipt],(int)effpt[1][ipt])<<endl;
      EPt.SetPoint(ipt,ipt*ptStep + ptStep/2., effpt[0][ipt]/effpt[1][ipt]);
      EPt.SetPointError(ipt, ptStep/2., StatUtils::EBinom(effpt[0][ipt]/effpt[1][ipt],(int)effpt[1][ipt]) );
      fill("EfficiencyPt","IsoPt", ipt*ptStep ,effpt[0][ipt]/effpt[1][ipt]);
    }
  }
  cout<<" Fin pt "<<endl;
  for(int ieta=0;ieta<neta;ieta++)
    if(effeta[1][ieta]!=0) {
      // cout<<ieta*etaStep -2.5 + etaStep/2.<<"  "<<effeta[0][ieta]/effeta[1][ieta]<<endl;
      EEta.SetPoint(ieta,ieta*etaStep -2.5 + etaStep/2., effeta[0][ieta]/effeta[1][ieta]);
      EEta.SetPointError(ieta, etaStep/2., StatUtils::EBinom(effeta[0][ieta]/effeta[1][ieta],(int)effeta[1][ieta]) );
      fill("EfficiencyEta","IsoEta", ieta*etaStep -2.5 ,effeta[0][ieta]/effeta[1][ieta]);
    }
  EPt.Write();
  EEta.Write();

  cout<<" Fin eta "<<endl;

}


void TagAndProbeAnalysis::EndReconstructionEfficiency() {
 
  vector<vector<float> > effeta(2,vector<float>(neta,0));
  vector<vector<float> > effpt(2,vector<float>(npt,0));
  
  float Num=0, Denom=0;

  int Ntt=0, Ntf=0, Ntp=0;

  for(int ipt=0;ipt<npt;ipt++) {
    for(int ieta=0;ieta<neta;ieta++) {
      effpt[0][ipt] += (float)(Ntt_Reco[ipt][ieta]+ Ntp_Reco[ipt][ieta]);
      effpt[1][ipt] += (float)(Ntt_Reco[ipt][ieta]+ Ntp_Reco[ipt][ieta]+ Ntf_Reco[ipt][ieta]);

      effeta[0][ieta] += (float)(Ntt_Reco[ipt][ieta]+ Ntp_Reco[ipt][ieta]);
      effeta[1][ieta] += (float)(Ntt_Reco[ipt][ieta]+ Ntp_Reco[ipt][ieta]+ Ntf_Reco[ipt][ieta]);
    
      Num += (float)(Ntt_Reco[ipt][ieta]+ Ntp_Reco[ipt][ieta]);
      Denom += (float)(Ntt_Reco[ipt][ieta]+ Ntp_Reco[ipt][ieta]+ Ntf_Reco[ipt][ieta]);
  
      Ntt += Ntt_Reco[ipt][ieta];
      Ntf += Ntf_Reco[ipt][ieta];
      Ntp += Ntp_Reco[ipt][ieta];

    }
  }
  cout<<" Fin calcul nombres reconstruction "<<endl;
  cout<<"  Ntagtag   "<<Ntt<<endl;
  cout<<"  Ntagpass   "<<Ntp<<endl;
  cout<<"  Ntagfail   "<<Ntf<<endl;
  cout<<" numerateur : "<<Num<<"   denominateur :  "<<Denom<<endl;
  cout<<" Eff : "<<Num/Denom<<" +- "<<StatUtils::EBinom(Num/Denom, (int)Denom)<<endl;
  _outputFile->cd();
  TGraphErrors EPt(npt);
  EPt.SetName("RecoPt"); 
  TGraphErrors EEta(neta);
  EEta.SetName("RecoEta");

  for(int ipt=0;ipt<npt;ipt++) {
    if(ipt<20)
      cout<<ipt*ptStep + ptStep/2.<<" Nevt "<<effpt[0][ipt]<<" : "<<effpt[1][ipt]<<endl;
    if(effpt[1][ipt]!=0) {
      if(ipt<20)
	cout<<ipt*ptStep + ptStep/2.<<"  "<<effpt[0][ipt]/effpt[1][ipt]<<" +-  "<<StatUtils::EBinom(effpt[0][ipt]/effpt[1][ipt],(int)effpt[1][ipt])<<endl;
      EPt.SetPoint(ipt,ipt*ptStep + ptStep/2., effpt[0][ipt]/effpt[1][ipt]);
      EPt.SetPointError(ipt, ptStep/2., StatUtils::EBinom(effpt[0][ipt]/effpt[1][ipt],(int)effpt[1][ipt]) );
      fill("EfficiencyPt","RecoPt", ipt*ptStep ,effpt[0][ipt]/effpt[1][ipt]);
    }
  }
  cout<<" Fin pt "<<endl;
  for(int ieta=0;ieta<neta;ieta++)
    if(effeta[1][ieta]!=0) {
      // cout<<ieta*etaStep -2.5 + etaStep/2.<<"  "<<effeta[0][ieta]/effeta[1][ieta]<<endl;
      EEta.SetPoint(ieta,ieta*etaStep -2.5 + etaStep/2., effeta[0][ieta]/effeta[1][ieta]);
      EEta.SetPointError(ieta, etaStep/2., StatUtils::EBinom(effeta[0][ieta]/effeta[1][ieta],(int)effeta[1][ieta]) );
      fill("EfficiencyEta","RecoEta", ieta*etaStep -2.5 ,effeta[0][ieta]/effeta[1][ieta]);
    }
  EPt.Write();
  EEta.Write();

  cout<<" Fin eta "<<endl;

}


void TagAndProbeAnalysis::EndIdentificationEfficiency() {
 
  vector<vector<float> > effeta(2,vector<float>(neta,0));
  vector<vector<float> > effpt(2,vector<float>(npt,0));
  
  float Num=0, Denom=0;

  int Ntt=0, Ntf=0, Ntp=0;

  for(int ipt=0;ipt<npt;ipt++) {
    for(int ieta=0;ieta<neta;ieta++) {
      effpt[0][ipt] += (float)(Ntt_ID[ipt][ieta]+ Ntp_ID[ipt][ieta]);
      effpt[1][ipt] += (float)(Ntt_ID[ipt][ieta]+ Ntp_ID[ipt][ieta]+ Ntf_ID[ipt][ieta]);

      effeta[0][ieta] += (float)(Ntt_ID[ipt][ieta]+ Ntp_ID[ipt][ieta]);
      effeta[1][ieta] += (float)(Ntt_ID[ipt][ieta]+ Ntp_ID[ipt][ieta]+ Ntf_ID[ipt][ieta]);
    
      Num += (float)(Ntt_ID[ipt][ieta]+ Ntp_ID[ipt][ieta]);
      Denom += (float)(Ntt_ID[ipt][ieta]+ Ntp_ID[ipt][ieta]+ Ntf_ID[ipt][ieta]);

      Ntt += Ntt_ID[ipt][ieta];
      Ntf += Ntf_ID[ipt][ieta];
      Ntp += Ntp_ID[ipt][ieta];
    }
  }
  cout<<" Fin calcul nombres identification "<<endl;
  cout<<"  Ntagtag   "<<Ntt<<endl;
  cout<<"  Ntagpass   "<<Ntp<<endl;
  cout<<"  Ntagfail   "<<Ntf<<endl;
  cout<<" numerateur : "<<Num<<"   denominateur :  "<<Denom<<endl;
  cout<<" Eff : "<<Num/Denom<<" +- "<<StatUtils::EBinom(Num/Denom, (int)Denom)<<endl;
  _outputFile->cd();
  TGraphErrors EPt(npt);
  EPt.SetName("IDPt"); 
  TGraphErrors EEta(neta);
  EEta.SetName("IDEta");

  for(int ipt=0;ipt<npt;ipt++) {
    if(ipt<20)
      cout<<ipt*ptStep + ptStep/2.<<" Nevt "<<effpt[0][ipt]<<" : "<<effpt[1][ipt]<<endl;
    if(effpt[1][ipt]!=0) {
      if(ipt<20)
	cout<<ipt*ptStep + ptStep/2.<<"  "<<effpt[0][ipt]/effpt[1][ipt]<<" +-  "<<StatUtils::EBinom(effpt[0][ipt]/effpt[1][ipt],(int)effpt[1][ipt])<<endl;
      EPt.SetPoint(ipt,ipt*ptStep + ptStep/2., effpt[0][ipt]/effpt[1][ipt]);
      EPt.SetPointError(ipt, ptStep/2., StatUtils::EBinom(effpt[0][ipt]/effpt[1][ipt],(int)effpt[1][ipt]) );
      fill("EfficiencyPt","IDPt", ipt*ptStep ,effpt[0][ipt]/effpt[1][ipt]);
    }
  }
  cout<<" Fin pt "<<endl;
  for(int ieta=0;ieta<neta;ieta++)
    if(effeta[1][ieta]!=0) {
      // cout<<ieta*etaStep -2.5 + etaStep/2.<<"  "<<effeta[0][ieta]/effeta[1][ieta]<<endl;
      EEta.SetPoint(ieta,ieta*etaStep -2.5 + etaStep/2., effeta[0][ieta]/effeta[1][ieta]);
      EEta.SetPointError(ieta, etaStep/2., StatUtils::EBinom(effeta[0][ieta]/effeta[1][ieta],(int)effeta[1][ieta]) );
      fill("EfficiencyEta","IDEta", ieta*etaStep -2.5 ,effeta[0][ieta]/effeta[1][ieta]);
    }
  EPt.Write();
  EEta.Write();

  cout<<" Fin eta "<<endl;

}







void
TagAndProbeAnalysis::AddCondition(CutMap& map, string name, vector<string> values) {

  std::pair<string,  vector<string> > tmp(name,values);
  map.push_back(tmp);
    
}




void
TagAndProbeAnalysis::EfficiencyMeasurement( const CandList& taglist, const CandList& probelist,
					    vector<int> _TagList,
					    vector<int> _ProbeList, string typeCut,
					    vector<string> valCut, NumTP& ntt, 
					    NumTP& ntf, NumTP& ntp) {

  //  cout<<typeCut<<endl;

  //Efficiency measurement
  CutMap dummy;
  TagAndProbeElectronSelector ProbeSelector(EventManager::kElectron, dummy );

  int cnttt=0;

  CandList Zcands;

  map<int, std::pair<vector<bool>, vector<float> > > Zs;

  for(int unsigned ie=0; ie<taglist.size();ie++) {
     
    Candidate* tagcand = taglist[ie];
    int uid = tagcand->uid();
    bool t=false;
     
    for(int unsigned ii=0;ii<_TagList.size();ii++)
      if(uid == _TagList[ii]) { t=true; break; }
     
    if(t) {
      for(int unsigned ie2=0; ie2<probelist.size();ie2++) {
	
	Candidate* probecand = probelist[ie2];
	 
	Candidate* Z= Candidate::create(tagcand,probecand);
	int uid2 = probecand->uid();

	int ptbin = (int)(probecand->pt()/ptStep);
	int etabin = (int)((probecand->eta()+2.5-etaStep)/etaStep );
	 
	//protection against very high/low pt/eta leptons overflow bins
	if(etabin<0)
	  etabin =0;
	if(etabin>=neta)
	  etabin=neta;
	if(ptbin>=npt)
	  ptbin=npt;

	float dr_ = KineUtils::dR(tagcand->eta(), probecand->eta(),
				  tagcand->phi(), probecand->phi() );

	//	cout<<dr_<<"    "<<t<<endl;

	if( uid != uid2 && dr_>0.1  && ( Z->mass()<120 &&  Z->mass()>60 ) ) {
	   
	  //  cout<<" pair "<<uid<<" <>  "<<uid2<<endl;

	  bool tt=false;
	  bool tp=false;
	  bool tf=false;

	  for(int unsigned ii=0;ii<_TagList.size();ii++)
	    if(uid2 == _TagList[ii]) { tt=true; break; }
	   
	  for(int unsigned ii=0;ii<_ProbeList.size();ii++)
	    if(uid2 == _ProbeList[ii]) { tp=true; break; }

	  if(typeCut!="Reco") {
	    if(tp) {
	      bool accCut = ProbeSelector.PassByType((const Candidate&) (*probecand), typeCut, valCut, tagcand );
	      if(!accCut)
		tf=true;
	    }
	  }
	  else {
	    tf=true;
	    for(int unsigned it=0; it<taglist.size();it++) {
	      bool accCut = ProbeSelector.PassByType((const Candidate&) (*probecand), typeCut, valCut, taglist[it] );
	      if(accCut && dr_>0.1 )
		{tf=false; break;}
	    }

	  }

	  //Filling numbers
	  
	  if(tt) //tag et probe sur les deux : les deux sonts comptés sur al double boucle
	    { ntt[ptbin][etabin]++; cnttt++; } 
	  else if(tp && tf) //coupure ratée
	    { ntf[ptbin][etabin]++; } 
	  else if(tp && !tf) //coupure passée
	    { ntp[ptbin][etabin]++; } 
	  
	  if(tp || tt) {
	    //   cout<<" t "<<t<<" tt  "<<tt<<"  tf   "<<(tp&&tf)<<"  tpass  "<<(tp&&!tf)<<"     "<<cnttt<<"   "<<ptbin<<"    "<<etabin<<"   "<<probecand->eta()<<endl;
	    
	    fillTree(typeCut, Z->mass(), tf, tp, tt, probecand->pt(), probecand->eta() );
	  }
	  
	  //Make all the Z cands
	  /*	  Zcands.push_back(Z);
	  
	  vector<bool> flags(3,0);
	  flags[0] = tt;
	  flags[1] = tp;
	  flags[2] = tf;
	  
	  vector<float> etapt(4,0);
	  etapt[0] =  probecand->pt();
	  etapt[1] =  probecand->eta();
	  etapt[2] = ptbin;
	  etapt[3] = etabin;
	  
	  std::pair<vector<bool>, vector<float> > tmp(flags, etapt);
	  Zs[ Z->uid() ] = tmp;
	  */

	}// condition masse
      }// 2e boucle
    } //condition tag
  } //1ere boucle
  

  /*
  //Filling numbers
 
  Candidate* bestZ;
  if(Zcands.size() > 1)
    bestZ = BestZCand(Zcands);
  else if(Zcands.size() == 1)
    bestZ = Zcands[0];
  else
    return;

  // vector<bool> flags(3,0);
  // vector<float> etapt(4,0);
 
  map<int, std::pair<vector<bool>, vector<float> > >::const_iterator itermap;

  itermap = Zs.find( bestZ->uid());
  if(itermap==Zs.end()) cout<<" gu"<<endl;
 

  bool tt = ((*itermap).second.first)[0];
  bool tp = ((*itermap).second.first)[1];
  bool tf = ((*itermap).second.first)[2];

  float pt = ((*itermap).second.second)[0];
  float eta = ((*itermap).second.second)[1];
  int ptbin = (int)(((*itermap).second.second)[2]);
  int etabin = (int)(((*itermap).second.second)[3]);
  
  if(tt) //tag et probe sur les deux : faut compter deux fois
    { ntt[ptbin][etabin]+=2; cnttt++;} 
  else if(tp && tf) //coupure ratée
    { ntf[ptbin][etabin]++; } 
  else if(tp && !tf) //coupure passée
      { ntp[ptbin][etabin]++; } 
  
  //cout<<typeCut<<"    "  <<" tt  "<<tt<<"  tf   "<<(tp&&tf)<<"  tpass  "<<(tp&&!tf)<<"     "<<cnttt<<"   "<<ptbin<<"    "<<etabin<<"   "<<pt<<"   "<<eta<<endl;
  
  fillTree(typeCut, bestZ->mass(), tf, tp, tt, pt, eta);
  */
	  
}



void 
TagAndProbeAnalysis::initRootTree() {


  outTreeID->Branch("ptprobeID",&PptID,"PptID/F");
  outTreeID->Branch("etaprobeID",&PetaID,"PetaID/F");
  outTreeID->Branch("ZMassID",&ZMassID,"ZMassID/F");
  outTreeID->Branch("IDPass",&IDPass,"IDPass/I");
  outTreeID->Branch("IDFail",&IDFail,"IDFail/I");
  outTreeID->Branch("IDTT",&IDTT,"IDTT/I");
  outTreeID->Branch("IDNVtx",&IDNvtx,"IDNvtx/I");
  outTreeID->Branch("IDNJet",&IDNJet,"IDNJet/I");

  outTreeIso->Branch("ptprobeIso",&PptIso,"PptIso/F");
  outTreeIso->Branch("etaprobeIso",&PetaIso,"PetaIso/F");
  outTreeIso->Branch("ZMassIso",&ZMassIso,"ZMassIso/F");
  outTreeIso->Branch("IsoPass",&IsoPass,"IsoPass/I");
  outTreeIso->Branch("IsoFail",&IsoFail,"IsoFail/I");
  outTreeIso->Branch("IsoTT",&IsoTT,"IsoTT/I");
  outTreeIso->Branch("IsoNVtx",&IsoNvtx,"IsoNvtx/I");
  outTreeIso->Branch("IsoNJet",&IsoNJet,"IsoNJet/I");
  
  outTreeReco->Branch("ptprobeReco",&PptReco,"PptReco/F");
  outTreeReco->Branch("etaprobeReco",&PetaReco,"PetaReco/F");
  outTreeReco->Branch("ZMassReco",&ZMassReco,"ZMassReco/F");
  outTreeReco->Branch("RecoPass",&RecoPass,"RecoPass/I");
  outTreeReco->Branch("RecoFail",&RecoFail,"RecoFail/I");
  outTreeReco->Branch("RecoTT",&RecoTT,"RecoTT/I");
  outTreeReco->Branch("RecoNVtx",&RecoNvtx,"RecoNvtx/I");
  outTreeReco->Branch("RecoNJet",&RecoNJet,"RecoNJet/I");

}



void
TagAndProbeAnalysis::fillTree(string type, float mass, bool tf, bool tp, bool tt, float pt, float eta)
{
  
  if(type=="ID") {
    ZMassID = mass;
    
    PptID = pt;
    PetaID = eta;

    IDPass = (int)tp;
    IDFail = (int)tf;
    IDTT   = (int)tt;

    IDNvtx = _NVertex;
    IDNJet = _NJet;

    outTreeID->Fill();
  }

  if(type=="Iso") {

    ZMassIso = mass;
    
    PptIso = pt;
    PetaIso = eta;

    IsoPass = (int)tp;
    IsoFail = (int)tf;
    IsoTT   = (int)tt;

    IsoNvtx = _NVertex;
    IsoNJet = _NJet;

    outTreeIso->Fill();

  }
  
  if(type=="Reco") {
    ZMassReco = mass;
    
    PptReco = pt;
    PetaReco = eta;

    RecoPass = (int)tp;
    RecoFail = (int)tf;
    RecoTT   = (int)tt;

    RecoNvtx = _NVertex;
    RecoNJet = _NJet;

    outTreeReco->Fill();

  }


}




Candidate *
TagAndProbeAnalysis::BestZCand( const CandList& zllList ) {

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

