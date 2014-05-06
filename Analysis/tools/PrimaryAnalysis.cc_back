#include "PrimaryAnalysis.hh"

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <map>
#include <math.h>

#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TH1F.h"
#include "TCut.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "THStack.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TTree.h"
#include "TChain.h"

#include <RooRealVar.h>
#include <RooFormulaVar.h>
#include <RooStringVar.h>
#include <RooDataSet.h>
#include <RooArgSet.h>
#include <RooArgList.h>
#include <RooPlot.h>
#include <RooAddPdf.h>
#include <RooFitResult.h>
#include <RooDataHist.h>
#include <RooRandom.h>
#include <RooGaussian.h>
#include <RooExponential.h>
#include <RooPolynomial.h>
#include <RooGExpModel.h>
#include <RooVoigtian.h>
#include <RooCBShape.h>
#include <RooArgusBG.h>
#include <RooBreitWigner.h>
#include <RooBifurGauss.h>
#include <RooFFTConvPdf.h>
#include <RooNumConvPdf.h>
#include <RooNovosibirsk.h>
#include <RooLandau.h>
#include <RooWorkspace.h>

using namespace std;

ClassImp(PrimaryAnalysis)

PrimaryAnalysis::PrimaryAnalysis( const char* configfile ) 
  :   _configfile(configfile) 
{  
  // initializations
  init();

  // configure
  configure();
}


PrimaryAnalysis::~PrimaryAnalysis()
{
}


void
PrimaryAnalysis::init()
{

  cout<<"================="<< endl;
  cout<<"Initializing ..."<< endl;

  ifstream fin;
  fin.open(_configfile);
  string datapath, newdatapath, outputpath, fitoutputpath, VBTFFile, weightsFile,ConfigFile,ConfigFullLumi;

  string nameRef;
  string nameReal;
  int k=0;

  while( !fin.eof() && k<8 ){    
    fin >> nameRef >> nameReal;
    if(nameRef=="datapath")datapath=nameReal;
    if(nameRef=="newdatapath") newdatapath=nameReal;
    if(nameRef=="outputpath") outputpath=nameReal;
    if(nameRef=="fitoutputpath") fitoutputpath=nameReal;
    if(nameRef=="VBTFFile") VBTFFile=nameReal;
    if(nameRef=="WeightsFile") weightsFile=nameReal;   
    if(nameRef=="ConfigFile") ConfigFile=nameReal;    
    if(nameRef=="ConfigFullLumi") ConfigFullLumi=nameReal;    
    k++;
  }
  cout<<" ... Done "<< endl;
  fin.close();

  _datapath=datapath;
  _newdatapath=newdatapath;
  _outputpath=outputpath;
  _VBTFFile=VBTFFile;
  _weightsFile=weightsFile;
  _ConfigFullLumi=ConfigFullLumi;

}


void
PrimaryAnalysis::configure(){

  cout<<"================="<< endl;
  cout<<"Configuring ..."<< endl;
  ReadWeights();
  ReadVBTF();
  ReadFullLumiFile();
  cout<<" ... Done "<< endl;
}

void 
PrimaryAnalysis::ReadWeights(){
    
  cout<<" "<< endl;
  cout<<"    reading weights in: "<<_weightsFile<<endl;

  
  ifstream fin;
  fin.open(_weightsFile.c_str());
  string nameSpecy;
  string name[5];
  double nEvt, xSect, filterEff, effXSect, lumi;
  int k=0;

  while( !fin.eof() ){
    {
      if(k>0){
	fin >> nameSpecy >> nEvt>> xSect>> filterEff>> effXSect>> lumi;
	//cout<<k<<" "<< nameSpecy<<" "<< nEvt<<" "<< effXSect<< endl;
	_Weights[nameSpecy]=(effXSect/nEvt);
      }else{
	fin>>nameSpecy>>name[0]>>name[1]>>name[2]>>name[3]>>name[4];
	k++;
      }
      
    }
  }
  cout<<"    ... done "<< endl;
  fin.close();
}


void 
PrimaryAnalysis::ReadVBTF(){
    
  cout<<" "<< endl;
  cout<<"    reading VBTF cuts in: "<<_VBTFFile<< endl;

  ifstream fin;
  fin.open(_VBTFFile.c_str());

  string namePart;
  string nameVar;
  int percent;
  double cutMax;
  int iPart, iPer, iVar, iType;

  while( !fin.eof()){
    
    fin >> namePart>> nameVar>> percent>> cutMax;    
    //    cout<< namePart<<" "<<nameVar<<" "<<percent<<" "<<cutMax<< endl;
  
    if(namePart=="EB") iPart=0;
    else iPart=1;
    
    if(nameVar=="trackIso" || nameVar=="ecalIso" || nameVar=="hcalIso" ) iType=0;
    else iType=1;
    
    if(nameVar=="trackIso" || nameVar=="sigma2ieta")iVar=0;
    else if (nameVar=="ecalIso" || nameVar=="dEta")iVar=1;
    else if (nameVar=="hcalIso" || nameVar=="dPhi")iVar=2;
    else iVar=3;
    
    if (percent ==95 ) iPer=3;
    else if (percent ==90 ) iPer=2;
    else if (percent ==80 ) iPer=1;
    else iPer=0;
    
    if(iType==0) _IsoCuts[iPart][iPer][iVar]=cutMax;
    else _IDCuts[iPart][iPer][iVar]=cutMax;
    
  }
  fin.close();

  cout<<"    ... done "<< endl;
}


void 
PrimaryAnalysis::ReadFullLumiFile(){
    
  cout<<" "<< endl;
  cout<<"    reading lumi files: "<<_ConfigFullLumi<<endl;

  
  ifstream fin;
  fin.open(_ConfigFullLumi.c_str());
  string nameFile;


  while( !fin.eof() )
    {
      fin >>nameFile;
      LumiFiles.push_back(nameFile);
    }
  cout<<" n Lumi Files=" <<LumiFiles.size()<< endl;
  
  cout<<"    ... done "<< endl;
  fin.close();
}

pair<float, float>
PrimaryAnalysis::VBTFSelection( int percent, const string & var, const string & ecalPart){
  
  float min=-9999.;
  float max=9999.;
  float cut=-9999;

  pair<float, float> output;
  
  int iEcalPart=0;
  if( ecalPart!= "EE" && ecalPart!= "EB"){ 
    cout <<" Wrong ecalPart: "<< ecalPart << endl;
    return output;
  }

  if( ecalPart== "EE") iEcalPart=1;
  
  int ipercent=4;

  //cout<< "percent: "<< percent<< endl;

  if (percent ==95 ) ipercent=3;
  else if (percent ==90 ) ipercent=2;
  else if (percent ==80 ) ipercent=1;
  else if (percent ==70 ) ipercent=0;
  
  
  //  cout<< "ipercent: "<< ipercent<<" "<< percent<< endl;
  
  int ivar=0; 
  string type="ISO";

  if(ipercent<4){
    if( var == "s2ieta" || var =="dPhi"|| var=="dEta" || var=="HoverE"){
      type="ID";
      
      if     (  var == "s2ieta"   ) ivar=0;
      else if(  var == "dEta"     ) ivar=1;
      else if(  var == "dPhi"     ) ivar=2;
      else if(  var == "HoverE"   ) ivar=3; 
      cut=_IDCuts[iEcalPart][ipercent][ivar];
      
    }else{
      
      type="ISO";
      if     (  var == "trackIso"    ) ivar=0;
      else if(  var == "ecalIso"     ) ivar=1;
      else if(  var == "hcalIso"     ) ivar=2;
      cut=_IsoCuts[iEcalPart][ipercent][ivar];
      
    }
    
    if     (  var == "s2ieta"   ){ min=0. ; max=cut;} 
    else if(  var == "dEta"     ){ min=-cut ; max=cut;}  
    else if(  var == "dPhi"     ){ min=-cut ; max=cut;}   
    else if(  var == "HoverE"   ){ min=0. ;  max=cut;}   
    else if(  var == "trackIso"    ){ min=-10. ; max=cut;}   
    else if(  var == "ecalIso"     ){ min=-10. ; max=cut;}   
    else if(  var == "hcalIso"     ){ min=-10. ; max=cut;}    
  }
  
  output.first=min;
  output.second=max;
  return output;

}


void PrimaryAnalysis::AddWeight(const string & filename){

  
  double weight=1.0;
  if(filename != "FullLumi"  
     && filename !="W_en" && filename!="W_tn" ) weight=_Weights[filename];

  TTree * tree = 0;


  if(filename!="FullLumi"){
    
    TFile file( (_datapath+"/"+filename+".root").c_str() ,"UPDATE");
    tree = (TTree*)file.Get("METCommVariables");
    
    if(tree->FindBranch("W")==NULL) {
      
      double W;
      TBranch* newBranch = tree->Branch("W",&W,"W/D");
      for(int k=0;k<tree->GetEntries();k++)
	{
	    W = weight;
	    newBranch->Fill();
	}
      tree->Write();
      
    }
  }else{
    for(unsigned i=0;i<LumiFiles.size();i++){

      TFile file( (_datapath+"/"+LumiFiles[i]+".root").c_str() ,"UPDATE");
      tree = (TTree*)file.Get("METCommVariables");
      
      if(tree->FindBranch("W")==NULL) {
	
	double W;
	TBranch* newBranch = tree->Branch("W",&W,"W/D");
	for(int k=0;k<tree->GetEntries();k++)
	  {
	    W = weight;
	    newBranch->Fill();
	  }
	tree->Write();
	
      }
    }
  }
}

void PrimaryAnalysis::AddWeights(const string & specy){

  if(specy=="GamJet"){
    static string names[5] = {"GamJet_PT1520","GamJet_PT2030","GamJet_PT3050","GamJet_PT5080","GamJet_PT80120"};
    
    for (int is=0;is<5;is++){
      AddWeight(names[is]);
    }
  }else if(specy=="QCD"){
    static string names[6] = {"QCD_EM2030","QCD_EM3080","QCD_EM80170","QCD_bc2030","QCD_bc3080","QCD_bc80170"};
    for (int is=0;is<6;is++){
      AddWeight(names[is]);
    }
  }else if(specy=="W_en" || specy== "W_tn" || specy== "Z_2e" || specy=="FullLumi"){
    AddWeight(specy);
  }else{
    cout<< "Unknown specy "<< endl;
  }
}


void PrimaryAnalysis::AddAllWeights(){

  cout<<"========================== "<< endl;
  cout<<"Adding weights to Trees ..."<< endl;
  AddWeights("GamJet");
  AddWeights("QCD");
  AddWeights("W_en");
  AddWeights("W_tn");
  AddWeights("Z_2e");
  AddWeights("FullLumi");
  cout<<"... Done"<< endl;
 
}


void PrimaryAnalysis::ModifyTree(const string & filename){

  int percent[3];
  percent[0]=95;
  percent[1]=80;
  percent[2]=70;
  
  TFile fIn( (_datapath +"/"+ filename+".root").c_str());
  TTree* inTree = (TTree*)fIn.Get("METCommVariables");
  
  //cout<< "Tree from: "<<(_datapath +"/"+ filename).c_str()<< " opened"<< endl;

  // Declaration of leaf types
  Char_t          VetoE;
  Char_t          VetoJet;
  Char_t          ExpInHit;
  Int_t           Charge;
  Int_t           NVertex;
  Int_t           Run;
  Float_t         dPhiJetLep;
  Float_t         dPhiMETLep[5];
  Float_t         dPhiMETJet[5];
  Float_t         IDVar[6];
  Float_t         IsoVar[3];
  Float_t         Lepton[3];
  Float_t         PFLepton[3];
  Float_t         SCLepton[3];
  Float_t         EtSC;
  Float_t         Mets[5][2];
  Float_t         MetProj[5][2];
  Float_t         Recoils[5][2];
  Float_t         RecoilProj[5][2];
  Float_t         MT[5];
  Double_t        W;

  // List of branches
  TBranch        *b_VetoE;   //!
  TBranch        *b_VetoJet;   //!
  TBranch        *b_ExpInHit;   //!
  TBranch        *b_Charge;   //!
  TBranch        *b_dPhiJetLep;   //!
  TBranch        *b_dPhiMETLep;   //!
  TBranch        *b_dPhiMETJet;   //!
  TBranch        *b_IDVar;   //!
  TBranch        *b_IsoVar;   //!
  TBranch        *b_Lepton;   //!
  TBranch        *b_PFLepton;   //!
  TBranch        *b_SCLepton;   //!
  TBranch        *b_EtSC;   //!
  TBranch        *b_Mets;   //!
  TBranch        *b_MetProj;   //!
  TBranch        *b_Recoils;   //!
  TBranch        *b_RecoilProj;   //!
  TBranch        *b_MT;   //!
  TBranch        *b_W;   //!
  TBranch        *b_NVertex;   //!
  TBranch        *b_Run;   //!

  inTree->SetBranchAddress("VetoE", &VetoE, &b_VetoE);
  inTree->SetBranchAddress("VetoJet", &VetoJet, &b_VetoJet);
  inTree->SetBranchAddress("ExpInHit", &ExpInHit, &b_ExpInHit);
  inTree->SetBranchAddress("Charge", &Charge, &b_Charge);
  inTree->SetBranchAddress("dPhiJetLep", &dPhiJetLep, &b_dPhiJetLep);
  inTree->SetBranchAddress("dPhiMETLep", dPhiMETLep, &b_dPhiMETLep);
  inTree->SetBranchAddress("dPhiMETJet", dPhiMETJet, &b_dPhiMETJet);
  inTree->SetBranchAddress("IDVar", IDVar, &b_IDVar);
  inTree->SetBranchAddress("IsoVar", IsoVar, &b_IsoVar);
  inTree->SetBranchAddress("Lepton", Lepton, &b_Lepton);
  inTree->SetBranchAddress("PFLepton", PFLepton, &b_PFLepton);
  inTree->SetBranchAddress("SCLepton", SCLepton, &b_SCLepton);
  inTree->SetBranchAddress("EtSC", &EtSC, &b_EtSC);
  inTree->SetBranchAddress("Mets", Mets, &b_Mets);
  inTree->SetBranchAddress("MetProj", MetProj, &b_MetProj);
  inTree->SetBranchAddress("Recoils", Recoils, &b_Recoils);
  inTree->SetBranchAddress("RecoilProj", RecoilProj, &b_RecoilProj);
  inTree->SetBranchAddress("MT", MT, &b_MT);
  inTree->SetBranchAddress("NVertex", &NVertex, &b_NVertex);
  inTree->SetBranchAddress("Run", &Run, &b_Run);
  inTree->SetBranchAddress("W", &W, &b_W);
  
  
  string foutname;
  foutname=_newdatapath +"/"+ filename+".root";

  TFile fOut( foutname.c_str(),"RECREATE" );
  
  //cout<< "Output file "<<(_newdatapath +"/"+ filename).c_str()<< " created"<< endl;


  TTree* outTree = new TTree("AnalysisTree","AnalysisTree");

  int ID[3], Iso[3], All[3];
  int ID70, ID80, ID95;
  int Iso70, Iso80, Iso95;
  int All70, All80, All95;
  int InitSel, EmaSel;
  float         phi, eta, dPhi, dEta, s2ieta, HoverE; 
  float         phiSC, etaSC;
  float         mT, MET, pT, recoil;
  float         dphiMETJet, dphiJetLep, dphiMETLep;
  float         iso, ecalIso, trackIso, hcalIso;
  float         METPara, METPerp;

  // List of branches

  TBranch* b_phi; TBranch* b_eta;TBranch* b_pT;TBranch* b_phiSC; TBranch* b_etaSC;
  TBranch* b_ID70,* b_ID80,* b_ID95;
  TBranch* b_Iso70,* b_Iso80,* b_Iso95;
  TBranch* b_All70,* b_All80 ,* b_All95;
  TBranch* b_InitSel,* b_EmaSel;
  TBranch* b_dPhi; TBranch* b_dEta; TBranch* b_s2ieta; TBranch* b_HoverE;
  TBranch* b_mT; TBranch* b_MET; TBranch* b_recoil;
  TBranch* b_dphiMETJet; TBranch* b_dphiJetLep; TBranch* b_dphiMETLep;
  TBranch* b_iso; TBranch* b_ecalIso; TBranch* b_trackIso; TBranch* b_hcalIso;
  TBranch* b_METPara; TBranch* b_METPerp;

  outTree->Branch("VetoE", &VetoE);
  outTree->Branch("VetoJet", &VetoJet);
  outTree->Branch("ExpInHit", &ExpInHit);
  outTree->Branch("Charge", &Charge);
  outTree->Branch("eta", &eta);
  outTree->Branch("phi", &phi);
  outTree->Branch("etaSC", &etaSC);
  outTree->Branch("phiSC", &phiSC);
  outTree->Branch("pT", &pT);
  outTree->Branch("dPhi", &dPhi);
  outTree->Branch("dEta", &dEta);
  outTree->Branch("s2ieta", &s2ieta);
  outTree->Branch("HoverE", &HoverE);
  outTree->Branch("mT", &mT);
  outTree->Branch("MET", &MET);
  outTree->Branch("METPara", &METPara);
  outTree->Branch("METPerp", &METPerp);
  outTree->Branch("recoil", &recoil); 
  outTree->Branch("EtSC",&EtSC);
  outTree->Branch("dphiMETJet", &dphiMETJet);
  outTree->Branch("dphiJetLep", &dphiJetLep);
  outTree->Branch("dphiMETLep", &dphiMETLep);
  outTree->Branch("iso", &iso);
  outTree->Branch("ecalIso", &ecalIso);
  outTree->Branch("trackIso", &trackIso);
  outTree->Branch("hcalIso", &hcalIso);
  outTree->Branch("NVertex", &NVertex);
  outTree->Branch("Run", &Run);
  outTree->Branch("W", &W);
  outTree->Branch("ID70", &ID70);
  outTree->Branch("ID80", &ID80);
  outTree->Branch("ID95", &ID95);
  outTree->Branch("Iso70", &Iso70);
  outTree->Branch("Iso80", &Iso80);
  outTree->Branch("Iso95", &Iso95);
  outTree->Branch("All70", &All70);
  outTree->Branch("All80", &All80);
  outTree->Branch("All95", &All95);
  outTree->Branch("InitSel", &InitSel);
  outTree->Branch("EmaSel", &EmaSel);

  outTree->SetBranchAddress("VetoE", &VetoE, &b_VetoE);
  outTree->SetBranchAddress("VetoJet", &VetoJet, &b_VetoJet);
  outTree->SetBranchAddress("ExpInHit", &ExpInHit, &b_ExpInHit );
  outTree->SetBranchAddress("pT", &pT, &b_pT );
  outTree->SetBranchAddress("phi", &phi, &b_phi );
  outTree->SetBranchAddress("etaSC", &etaSC, &b_etaSC );
  outTree->SetBranchAddress("phiSC", &phiSC, &b_phiSC );
  outTree->SetBranchAddress("eta", &eta, &b_eta );
  outTree->SetBranchAddress("dPhi", &dPhi, &b_dPhi );
  outTree->SetBranchAddress("dEta", &dEta, &b_dEta );
  outTree->SetBranchAddress("s2ieta", &s2ieta, &b_s2ieta );
  outTree->SetBranchAddress("HoverE", &HoverE, &b_HoverE );
  outTree->SetBranchAddress("mT", &mT, &b_mT );
  outTree->SetBranchAddress("MET", &MET, &b_MET );
  outTree->SetBranchAddress("METPara", &METPara, &b_METPara );
  outTree->SetBranchAddress("METPerp", &METPerp, &b_METPerp );
  outTree->SetBranchAddress("recoil", &recoil, &b_recoil );
  outTree->SetBranchAddress("EtSC", &EtSC, &b_EtSC );
  outTree->SetBranchAddress("dphiMETJet", &dphiMETJet, &b_dphiMETJet );
  outTree->SetBranchAddress("dphiJetLep", &dphiJetLep, &b_dphiJetLep );
  outTree->SetBranchAddress("dphiMETLep", &dphiMETLep, &b_dphiMETLep );
  outTree->SetBranchAddress("iso", &iso, &b_iso );
  outTree->SetBranchAddress("ecalIso", &ecalIso, &b_ecalIso );
  outTree->SetBranchAddress("trackIso", &trackIso, &b_trackIso );
  outTree->SetBranchAddress("hcalIso", &hcalIso, &b_hcalIso );
  outTree->SetBranchAddress("NVertex", &NVertex,&b_NVertex);
  outTree->SetBranchAddress("Run", &Run,&b_Run);
  outTree->SetBranchAddress("W", &W, &b_W );
  outTree->Branch("ID70", &ID70);
  outTree->Branch("ID80", &ID80);
  outTree->Branch("ID95", &ID95);
  outTree->Branch("Iso70", &Iso70);
  outTree->Branch("Iso80", &Iso80);
  outTree->Branch("Iso95", &Iso95);
  outTree->Branch("All70", &All70);
  outTree->Branch("All80", &All80);
  outTree->Branch("All95", &All95);
  outTree->Branch("InitSel", &InitSel);
  outTree->Branch("EmaSel", &EmaSel);
  outTree->SetBranchAddress("ID70", &ID70, &b_ID70 );
  outTree->SetBranchAddress("ID80", &ID80, &b_ID80 );
  outTree->SetBranchAddress("ID95", &ID95, &b_ID95 );
  outTree->SetBranchAddress("Iso70", &Iso70, &b_Iso70 );
  outTree->SetBranchAddress("Iso80", &Iso80, &b_Iso80 );
  outTree->SetBranchAddress("Iso95", &Iso95, &b_Iso95 );
  outTree->SetBranchAddress("All70", &All70, &b_All70 );
  outTree->SetBranchAddress("All80", &All80, &b_All80 );
  outTree->SetBranchAddress("All95", &All95, &b_All95 );
  outTree->SetBranchAddress("InitSel", &InitSel, &b_InitSel );
  outTree->SetBranchAddress("EmaSel", &EmaSel, &b_EmaSel );


  if (inTree == 0) return;
  
  float EtSCCut=20;
  string varCutName[7];
  varCutName[0]="dPhi";
  varCutName[1]="dEta";
  varCutName[2]="s2ieta";
  varCutName[3]="HoverE";
  varCutName[4]="ecalIso";
  varCutName[5]="trackIso";
  varCutName[6]="hcalIso";
  
  string ecalPart[2];
  ecalPart[0]="EB";
  ecalPart[1]="EE";

  pair<float,float> VBTFCuts[7][2][3];

  for (int kPer=0;kPer<3;kPer++){
    for (int i=0;i<7;i++){ // var
      for (int j=0;j<2;j++){ // ecalPart
	VBTFCuts[i][j][kPer]=VBTFSelection( percent[kPer], varCutName[i], ecalPart[j]);
	//cout<< varCutName[i]<<" "<<ecalPart[j]<<" "<<VBTFCuts[i][j][kPer].first<<" "<<VBTFCuts[i][j][kPer].second << endl;
	
      }
    }
  }

  int iFirstCutIso=4;
  int iLastCutIso=7;
  int iFirstCutID=0;
  int iLastCutID=4;
  
  Long64_t nentries = inTree->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    
    nb = inTree->GetEntry(jentry);   nbytes += nb;
    

    // ID Vars
    //=========
    

    s2ieta=IDVar[0];
    dEta=IDVar[1];
    dPhi=IDVar[2];

    //corrections for EE only:
    if(TMath::Abs(etaSC)>1.56 && TMath::Abs(etaSC)<2.5){
      dEta-=IDVar[4];
      dPhi-=IDVar[5];
    }
    HoverE=IDVar[3];

    // Iso Vars
    //===========

    trackIso=IsoVar[0];
    ecalIso=IsoVar[1];
    hcalIso=IsoVar[2];
    iso=Lepton[0]/(Lepton[0]+ecalIso+trackIso+hcalIso);

    // Other Kine Vars
    //================

    mT=MT[0];
    MET=sqrt(MetProj[0][0]*MetProj[0][0]+MetProj[0][1]*MetProj[0][1]);
    recoil=sqrt(RecoilProj[0][0]*RecoilProj[0][0]+RecoilProj[0][1]*RecoilProj[0][1]);;
    dphiMETJet=dPhiMETJet[0];
    dphiMETLep=dPhiMETLep[0];
    dphiJetLep=dPhiJetLep;
    pT=Lepton[0];
    eta=Lepton[1];
    phi=Lepton[2];
    etaSC=SCLepton[1];
    phiSC=SCLepton[2];
     
    METPara=MetProj[0][0];
    METPerp=MetProj[0][1];

    int jpart=0;
    string ecalPart="EB";
    if(TMath::Abs(etaSC)>1.56 && TMath::Abs(etaSC)<2.5){
      ecalPart="EE";
      jpart=1;
    }
    else if(TMath::Abs(etaSC)<1.44){
      ecalPart="EB";
      jpart=0;
    }else {
      cout<<" unknown ecal part"<<endl;
      return;
    }
    
    // EmaSel
    EmaSel=1;
    if(ExpInHit==0 || EtSC<EtSCCut || mT<20 ||  TMath::Abs(etaSC)>2.5 || HoverE>0.15 || ecalIso>4 || VetoJet==1 ) EmaSel=0; 
    
    // InitSel
    InitSel=1;
    if(VetoE==0 || ExpInHit==0 || EtSC<EtSCCut ) InitSel=0;
 

    float var[7];
    var[0]=dPhi;
    var[1]=dEta;
    var[2]=s2ieta;
    var[3]=HoverE;
    var[4]=ecalIso/pT;
    var[5]=trackIso/pT;
    var[6]=hcalIso/pT;

    for(int iPer=0;iPer<3;iPer++){
    
      All[iPer]=1;
      Iso[iPer]=1;
      ID[iPer]=1;
   
      for(int i=0;i<7;i++){
	if (var[i]<VBTFCuts[i][jpart][iPer].first || var[i]>VBTFCuts[i][jpart][iPer].second ){
	  All[iPer]=0;
	  if(i>=iFirstCutIso && i< iLastCutIso ) Iso[iPer]=0;
	  if(i>=iFirstCutID && i< iLastCutID ) ID[iPer]=0;
	}
      }
    }
    
    All70=All[2]; All80=All[1]; All95=All[0];
    ID70=ID[2]; ID80=ID[1]; ID95=ID[0];
    Iso70=Iso[2]; Iso80=Iso[1]; Iso95=Iso[0];
    
    outTree->Fill();
    
  }
  
  outTree->Write();
  fOut.Close();
}
  
double PrimaryAnalysis::errEff(double v1, double vtot){

  double N=vtot;
  double e=v1/N;
  return sqrt(e*(1-e)/N);

}
double PrimaryAnalysis::errRat(double v1, double v2, double s1, double s2){

  double rat=v1/v2;
  double sum= sqrt((s1/v1)*(s1/v1)+(s2/v2)*(s2/v2));

  return rat*sum;
}
double PrimaryAnalysis::errProd(double v1, double v2, double v3, double v4,double s1, double s2,double s3, double s4){

  double prod=v1*v2*v3*v4;
  double sum= sqrt((s1/v1)*(s1/v1)+(s2/v2)*(s2/v2)+(s3/v3)*(s3/v3)+(s4/v4)*(s4/v4));

  return prod*sum;
}


void PrimaryAnalysis::ModifyTrees(){

  cout<<"========================== "<< endl;
  cout<<"Generating new Trees..."<< endl;
  // GamJet
  string names[5] = {"GamJet_PT1520","GamJet_PT2030","GamJet_PT3050","GamJet_PT5080","GamJet_PT80120"};
  for (int is=0;is<5;is++){
    ModifyTree(names[is]);
  }
  // QCD
  string names2[6] = {"QCD_EM2030","QCD_EM3080","QCD_EM80170","QCD_bc2030","QCD_bc3080","QCD_bc80170"};
  for (int is=0;is<6;is++){
    ModifyTree(names2[is]);
  }

  ModifyTree("W_en");
  ModifyTree("W_tn");
  ModifyTree("Z_2e");

  ModifyTree("FullLumi");

  cout<<"... Done"<< endl;


}

TCut PrimaryAnalysis::GetCut( int percent, bool doIso, bool doID, bool doEmaAndNoInit, bool doFlipIso, bool doFlipID, bool isRoo) {

  stringstream percentString;
  percentString<<percent;
  
  string cutpercent;
  string cutFlip;

  if( doFlipID && doFlipIso){
    cutFlip="All"+percentString.str()+"==0";
    doIso=false;
    doID=false;
  }else if(doFlipID){
    cutFlip="ID"+percentString.str()+"==0";
    doID=false;
  }else if(doFlipIso){
    cutFlip="Iso"+percentString.str()+"==0";
    doIso=false;
  }else{
    cutFlip="";
  }
  
  if ( doIso && doID ){
    cutpercent="All"+percentString.str()+"==1";
  }else if (doIso){
    cutpercent="Iso"+percentString.str()+"==1";
  }else if (doID){
    cutpercent="ID"+percentString.str()+"==1";
  }else{
    cutpercent="";
  }
  
  //  cout<<"cutpercent: "<< cutpercent<< endl;

  string cutinit;
  if(doEmaAndNoInit){
    cutinit="EmaSel==1";
  }else{
    cutinit="InitSel==1";
  }

  string allcut;

  if(!isRoo){
    if(cutpercent!="" && cutFlip!="") allcut="W*("+cutpercent+" && "+cutinit+" && "+cutFlip+")";
    else if(cutpercent!="" && cutFlip=="") allcut="W*("+cutpercent+" && "+cutinit+")";
    else if(cutpercent=="" && cutFlip!="") allcut="W*("+cutinit+" && "+cutFlip+")";
    else allcut="W*("+cutinit+")";
  }else{

    if(cutpercent!="" && cutFlip!="") allcut=cutpercent+" && "+cutinit+" && "+cutFlip;
    else if(cutpercent!="" && cutFlip=="") allcut=cutpercent+" && "+cutinit;
    else if(cutpercent=="" && cutFlip!="") allcut=cutinit+" && "+cutFlip;
    else allcut=cutinit;
  }
  //cout<<allcut << endl;
  TCut cut(allcut.c_str());

  return cut;

}


string PrimaryAnalysis::CutName( int percent, bool doIso, bool doID, bool doEmaAndNoInit , bool doFlipIso, bool doFlipID) {

  string cutname;

  string cutinit;
  if(doEmaAndNoInit) cutinit="Ema";
  else cutinit="Init";

  stringstream percentString;
  percentString<<percent;
  string cutpercent;
  if ( doIso && doID ){
    cutpercent="All"+percentString.str();
  }else if (doIso){
    cutpercent="Iso"+percentString.str();
  }else if (doID){
    cutpercent="ID"+percentString.str();
  }else{
    cutpercent="";
  }  

  string cutFlip;
  if( doFlipID && doFlipIso){
    cutFlip="FlipAll"+percentString.str();
    doIso=false;
    doID=false;
  }else if(doFlipID){
    cutFlip="FlipID"+percentString.str();
    doID=false;
  }else if(doFlipIso){
    cutFlip="FlipIso"+percentString.str();
    doIso=false;
  }else{
    cutFlip="";
  }
  
  cutname=cutinit+cutpercent+cutFlip;
  return cutname; 

}


std::vector<double> PrimaryAnalysis::efficiency(const string & specy, TCut cut, TCut cutRef, double METParaMin, double METParaMax){
  

  TChain* chain=new TChain("AnalysisTree");
  
  if(specy=="GamJet"){
    string names[5] = {"GamJet_PT1520","GamJet_PT2030","GamJet_PT3050","GamJet_PT5080","GamJet_PT80120"};
    int nt=5;
 
    for (int i=0;i<nt;i++){
      chain->Add((_newdatapath +"/"+ names[i]+".root").c_str());  
      //cout<< "chain "<<chain->GetEntries()<<endl;
    }
    //cout<< " End chain "<<chain->GetEntries()<<endl;
  }else if(specy=="QCD"){
    string names[6] = {"QCD_EM2030","QCD_EM3080","QCD_EM80170","QCD_bc2030","QCD_bc3080","QCD_bc80170"};
    int nt=6;
    for (int i=0;i<nt;i++){
      chain->Add((_newdatapath+"/"+names[i]+".root").c_str());  
      //cout<< "chain "<<chain->GetEntries()<<endl;
    }
    //cout<< " End chain "<<chain->GetEntries()<<endl;
  }else if(specy=="W_en"){
    chain->Add((_newdatapath +"/"+"W_en"+".root").c_str());
  }else if(specy== "W_tn"){
    chain->Add((_newdatapath +"/"+"W_tn"+".root").c_str());
  }else if(specy== "Z_2e"){ 
    chain->Add((_newdatapath +"/"+"Z_2e"+".root").c_str());
  }else{
    cout<< "Unknown specy "<< endl;
  }
  TH1D *hTrue = new TH1D("hTrue","hTrue",500,METParaMin,METParaMax);
  TH1D *hRef = new TH1D("hRef","hRef",500,METParaMin,METParaMax);
  
  chain->Draw("METPara>>hTrue", cut, "goff");       
  chain->Draw("METPara>>hRef" , cutRef, "goff" );                         

  double nok =hTrue->Integral("w") ;
  double nref=hRef->Integral("w") ;

  double eff=nok/nref;
  double err=errEff(nok,nref);

  //cout<<" Eff= "<<eff<<" +/- "<< err<< endl;

  std::vector<double> res;
  res.push_back(eff);
  res.push_back(err);

  delete hTrue;
  delete hRef;


  return res;

}

void PrimaryAnalysis::AllEfficiencies(const string & specy){

  double METParaMin=-50.;
  double METParaMax=80.;

  int percent[3];
  bool doID[2];  
  bool doIso[2];
  bool doEmaAndNoInit[2];

  percent[0]=95;
  percent[1]=80;
  percent[2]=70; 
  doID[0]=false;
  doID[1]=true;
  doIso[0]=false;
  doIso[1]=true;
  doEmaAndNoInit[0]=false;
  doEmaAndNoInit[1]=true;

  for (int ipercent=0;ipercent<3;ipercent++){
    for (int jinit=0;jinit<2;jinit++){
      for (int jid=0;jid<2;jid++){ 
	for (int jiso=0;jiso<2;jiso++){
	  string cutName=CutName(percent[ipercent], doIso[jiso], doID[jid],doEmaAndNoInit[jinit],false,false);
	  //cout<<cutName<< endl;
	  TCut cut= GetCut(percent[ipercent], doIso[jiso], doID[jid],doEmaAndNoInit[jinit], false, false, false);
	  TCut cutRef= GetCut(percent[ipercent], false, false,doEmaAndNoInit[jinit], false, false, false);
	  std::vector<double> eff=efficiency(specy, cut, cutRef,METParaMin,METParaMax);

	  string EffName="Eff"+cutName;
	  RooRealVar Eff("Eff","Eff", eff[0]);
	  Eff.setRange( eff[0]-eff[1],eff[0]+eff[1]);
	  
	  RooArgList roolist(Eff);
	  RooArgSet* params= new RooArgSet(roolist);
	  params->writeToFile((_outputpath +"/Efficiency_"+specy+"_"+cutName+".txt").c_str());
	  
	}
      }
    }   
  }
  
}
void PrimaryAnalysis::AllEff(){

  cout<<"========================"<< endl;
  cout<<"Getting efficiencies ..."<< endl;

  AllEfficiencies("W_en");
  cout<<"  W_en done "<< endl;
  AllEfficiencies("W_tn");
  cout<<"  W_tn done "<< endl;
  AllEfficiencies("Z_2e");
  cout<<"  Z_2e done "<< endl;
  AllEfficiencies("QCD");
  cout<<"  QCD done "<< endl;
  AllEfficiencies("GamJet");
  cout<<"  GamJet done "<< endl;

  cout<< "... Done"<< endl;

}


std::vector<double> PrimaryAnalysis::WTauEEff(int percent, bool doID, bool doIso, bool doEmaAndNoInit, double METParaMin, double METParaMax, bool doFlipID, bool doFlipIso ){

  double eMat, eInitSel, ePercentSel;
  double sMat, sInitSel, sPercentSel;
  double PDGratio,sPDGratio;


  double eInitSelTau, ePercentSelTau;
  double sInitSelTau, sPercentSelTau;
  double eInitSelE, ePercentSelE;
  double sInitSelE, sPercentSelE;


  // PDG ratio:
  double PDGTau=11.25;
  double sPDGTau=0.20;
  double PDGE=10.75;
  double sPDGE=0.13;
  PDGratio=PDGTau/PDGE;
  sPDGratio=errRat(PDGTau,PDGE,sPDGTau,sPDGE);

  // Mat eff:
  double nMatTau=3579.;
  double nMatE=43195.;
  eMat=nMatTau/nMatE;
  sMat=errRat(nMatTau,nMatE,sqrt(nMatTau),sqrt(nMatE));


  // Get right numbers:
  


  TCut cutInit= GetCut(percent, false, false,doEmaAndNoInit, false, false, false);
  TCut noCut("");

  std::vector<double> effE=efficiency("W_en", cutInit, noCut, METParaMin, METParaMax);
  std::vector<double> effTau=efficiency("W_tn", cutInit, noCut, METParaMin, METParaMax);
  
  eInitSelTau=effTau[0];
  eInitSelE=effE[0];
  sInitSelTau=effTau[1];
  sInitSelE=effE[1];

  // Init Sel Eff


  eInitSel=eInitSelTau/eInitSelE;
  sInitSel=errRat(eInitSelTau,eInitSelE,sInitSelTau,sInitSelE);

  //cout<< " Init Sel Ratio: " <<eInitSel<<"+/-"<<sInitSel<< endl;

    
  // ID and Iso Sel Eff


  TCut cutSel= GetCut(percent, doIso, doID,doEmaAndNoInit,doFlipIso,doFlipID,false);
 
  std::vector<double> effESel=efficiency("W_en", cutSel, cutInit, METParaMin, METParaMax);
  std::vector<double> effTauSel=efficiency("W_tn", cutSel, cutInit, METParaMin, METParaMax);
  ePercentSelTau=effTauSel[0];
  ePercentSelE=effESel[0];
  sPercentSelTau=effTauSel[1];
  sPercentSelE=effESel[1];

  ePercentSel=ePercentSelTau/ePercentSelE;
  sPercentSel=errRat(ePercentSelTau,ePercentSelE,sPercentSelTau,sPercentSelE);

  //cout<< " Percent Sel Ratio: " <<ePercentSel<<"+/-"<<sPercentSel<< endl;

  double totalProd=PDGratio*eMat*eInitSel*ePercentSel;
  double totalErr=errProd(PDGratio,eMat,eInitSel,ePercentSel,sPDGratio,sMat,sInitSel,sPercentSel);

  std::vector<double> res;
  res.push_back(totalProd);
  res.push_back(totalErr);

  //cout<<" ============================================================" << endl;
  //cout<<" t/e ratio including efficiencies: " <<totalProd<<" +/- "<<totalErr<< endl;
  //cout<<" ============================================================" << endl;


  string cutName=CutName(percent, doIso, doID,doEmaAndNoInit,doFlipIso,doFlipID);

  RooRealVar TauERatio("TauERatio","TauERatio", totalProd);
  TauERatio.setRange( totalProd-totalErr,totalProd+totalErr);
  TauERatio.setConstant();

  RooArgList roolist(TauERatio);
  RooArgSet* params= new RooArgSet(roolist);
  params->writeToFile((_outputpath +"/"+ "TauERatio_"+cutName+".txt").c_str());
   
  return res;  

}

void PrimaryAnalysis::AllWTauEff(){

  cout<<"=============================="<< endl;
  cout<<"Getting W->tn/W->en ratios ..."<< endl;

  double METParaMin=-50.;
  double METParaMax=80.;

  int percent[3];
  bool doID[2];  
  bool doIso[2];
  bool doEmaAndNoInit[2];

  percent[0]=95;
  percent[1]=80;
  percent[2]=70; 
  doID[0]=false;
  doID[1]=true;
  doIso[0]=false;
  doIso[1]=true;
  doEmaAndNoInit[0]=false;
  doEmaAndNoInit[1]=true;

  for (int ipercent=0;ipercent<3;ipercent++){
    for (int jinit=0;jinit<2;jinit++){
      for (int jid=0;jid<2;jid++){ 
	for (int jiso=0;jiso<2;jiso++){
	  std::vector<double> rat=WTauEEff(percent[ipercent], doIso[jiso], doID[jid],doEmaAndNoInit[jinit],METParaMin,METParaMax);	  
	}
      }
    }   
  }

  std::vector<double> rat1=WTauEEff(95, true, false, false,METParaMin,METParaMax, false, true);	    
  
  cout<<"... Done "<<endl;
  
}

std::vector<double> PrimaryAnalysis::ZEEEEff(int percent, bool doID, bool doIso, bool doEmaAndNoInit, double METParaMin, double METParaMax, bool doFlipIso , bool doFlipID ){



  double eMat, eInitSel, ePercentSel;
  double sMat, sInitSel, sPercentSel;
  double PDGratio,sPDGratio;

  double eInitSelTau, ePercentSelTau;
  double sInitSelTau, sPercentSelTau;
  double eInitSelE, ePercentSelE;
  double sInitSelE, sPercentSelE;


  // PDG ratio:
  double PDGZEE=3.363;
  double sPDGZEE=0.40;
  double PDGE=10.75;
  double sPDGE=0.13;
  PDGratio=PDGZEE/PDGE;
  sPDGratio=errRat(PDGZEE,PDGE,sPDGZEE,sPDGE);

  // Mat eff:
  double nMatTau=6324.;
  double nMatE=43195.;
  eMat=nMatTau/nMatE;
  sMat=errRat(nMatTau,nMatE,sqrt(nMatTau),sqrt(nMatE));


  // Get right numbers:
  TCut cutInit= GetCut(percent, false, false,doEmaAndNoInit, false, false, false);
  TCut noCut("");

  std::vector<double> effE=efficiency("W_en", cutInit, noCut, METParaMin, METParaMax);
  std::vector<double> effTau=efficiency("Z_2e", cutInit, noCut, METParaMin, METParaMax);
  
  eInitSelTau=effTau[0];
  eInitSelE=effE[0];
  sInitSelTau=effTau[1];
  sInitSelE=effE[1];

  // Init Sel Eff


  eInitSel=eInitSelTau/eInitSelE;
  sInitSel=errRat(eInitSelTau,eInitSelE,sInitSelTau,sInitSelE);

  //  cout<< " Init Sel Ratio: " <<eInitSel<<"+/-"<<sInitSel<< endl;

    
  // ID and Iso Sel Eff
  
  TCut cutSel= GetCut(percent, doIso, doID,doEmaAndNoInit,doFlipIso,doFlipID, false);
 
  std::vector<double> effESel=efficiency("W_en", cutSel, cutInit, METParaMin, METParaMax);
  std::vector<double> effTauSel=efficiency("Z_2e", cutSel, cutInit, METParaMin, METParaMax);
  ePercentSelTau=effTauSel[0];
  ePercentSelE=effESel[0];
  sPercentSelTau=effTauSel[1];
  sPercentSelE=effESel[1];

  ePercentSel=ePercentSelTau/ePercentSelE;
  sPercentSel=errRat(ePercentSelTau,ePercentSelE,sPercentSelTau,sPercentSelE);

  //  cout<< " Percent Sel Ratio: " <<ePercentSel<<"+/-"<<sPercentSel<< endl;

  double totalProd=PDGratio*eMat*eInitSel*ePercentSel;
  double totalErr=errProd(PDGratio,eMat,eInitSel,ePercentSel,sPDGratio,sMat,sInitSel,sPercentSel);

  std::vector<double> res;
  res.push_back(totalProd);
  res.push_back(totalErr);

  //cout<<" ============================================================" << endl;
  //cout<<" Z->2e to W->en ratio including efficiencies: " <<totalProd<<" +/- "<<totalErr<< endl;
  //cout<<" ============================================================" << endl;


  string cutName=CutName(percent, doIso, doID,doEmaAndNoInit,doFlipIso,doFlipID);

  RooRealVar ZEERatio("ZEERatio","ZEERatio", totalProd);
  ZEERatio.setRange( totalProd-totalErr,totalProd+totalErr);
  ZEERatio.setConstant();

  RooArgList roolist(ZEERatio);
  RooArgSet* params= new RooArgSet(roolist);
  params->writeToFile((_outputpath +"/"+ "ZEERatio_"+cutName+".txt").c_str());
   
  return res;  

}
void PrimaryAnalysis::AllZEEEff(){

  cout<<"=============================="<< endl;
  cout<<"Getting Z->2e/W->en ratios ..."<< endl;

  double METParaMin=-50.;
  double METParaMax=80.;

  int percent[3];
  bool doID[2];  
  bool doIso[2];
  bool doEmaAndNoInit[2];

  percent[0]=95;
  percent[1]=80;
  percent[2]=70; 
  doID[0]=false;
  doID[1]=true;
  doIso[0]=false;
  doIso[1]=true;
  doEmaAndNoInit[0]=false;
  doEmaAndNoInit[1]=true;

  for (int ipercent=0;ipercent<3;ipercent++){
    for (int jinit=0;jinit<2;jinit++){
      for (int jid=0;jid<2;jid++){ 
	for (int jiso=0;jiso<2;jiso++){
	  std::vector<double> rat=ZEEEEff(percent[ipercent], doIso[jiso], doID[jid],doEmaAndNoInit[jinit],METParaMin,METParaMax);	  
	}
      }
    }   
  }
  
  std::vector<double> rat2=ZEEEEff(95, true, false, false, METParaMin,METParaMax, false, true);	 
  cout<<"... Done "<<endl;
}

void PrimaryAnalysis::fitMETPara(const string & specy, TCut cut, const string & cutName){



  TChain* chain=new TChain("AnalysisTree");
  
  if(specy=="GamJet" || specy=="MCBkg" ){
    string names[5] = {"GamJet_PT1520","GamJet_PT2030","GamJet_PT3050","GamJet_PT5080","GamJet_PT80120"};
    int nt=5;
    
    for (int i=0;i<nt;i++){
      chain->Add((_newdatapath +"/"+ names[i]+".root").c_str());  
      //cout<< "chain "<<chain->GetEntries()<<endl;
    }
    //cout<< " End chain "<<chain->GetEntries()<<endl;
  }else if(specy=="QCD"|| specy=="MCBkg" ){
    string names[6] = {"QCD_EM2030","QCD_EM3080","QCD_EM80170","QCD_bc2030","QCD_bc3080","QCD_bc80170"};
    int nt=6;
    for (int i=0;i<nt;i++){
      chain->Add((_newdatapath+"/"+names[i]+".root").c_str());  
      //cout<< "chain "<<chain->GetEntries()<<endl;
    }
    //cout<< " End chain "<<chain->GetEntries()<<endl;
  }else if(specy=="W_en"){
    chain->Add((_newdatapath +"/W_en.root").c_str());
  }else if(specy== "W_tn"){
    chain->Add((_newdatapath +"/W_tn.root").c_str());
  }else if(specy== "Z_2e" ){ 
    chain->Add((_newdatapath +"/Z_2e.root").c_str());
  }else if(specy=="FullLumi"){
    chain->Add((_newdatapath +"/FullLumi.root").c_str());
  } else{
    cout<< "Unknown specy "<< endl;
  }


  RooRealVar METPara("METPara","METPara",-50.,80.);
  RooRealVar etaSC("etaSC","etaSC",-2.5,2.5);
  RooRealVar W("W","W",0,2);
  

  RooArgList AllVars;
  AllVars.add(RooArgList(METPara,etaSC, W));


  TTree* chainSel = (TTree*) chain->CopyTree(cut);
  RooDataSet DataSet("data","data",chainSel,AllVars,"","W");
  
  //cout<<" Initial Tree entries:  "<<chain->GetEntries()<< endl;
  //cout<<" Selected Tree entries: "<<chainSel->GetEntries()<< endl;
  //cout<<" RooDataset entries:    "<<DataSet.numEntries()<< endl;
  
  // METParaBkg
  RooRealVar sigRBkg("sigmaRBkg","sigmaRBkg",18,0.0,100.);
  RooRealVar sigLBkg("sigmaLBkg","sigmaLBkg",22,0.0,100.);
  RooRealVar meanBifBkg("meanBifBkg","meanBif",0.,-100.0,100.0);
  RooRealVar sigBkg("sigmaBkg","sigmaBkg",2.,1.0,100.);
  RooRealVar meanBkg("meanBkg","meanBkg",4.,2.0,100.0);
  RooRealVar coefBkg("coefBkg","coefBkg",0.6,0.,1.0);

  RooBifurGauss BifGaussBkg("BifGaussBkg","BifGaussBkg",METPara,meanBifBkg,sigLBkg, sigRBkg);
  RooGaussian GaussBkg("GaussBkg","GaussBkg",METPara,meanBkg,sigBkg);  
  RooAddPdf METParaBkg("METBkg","METBkg",BifGaussBkg, GaussBkg,coefBkg);


  // METParaSig
  RooRealVar coefSig("coefSig","coefSig",0.752,0.,1.0);
  RooRealVar sigRSig("sigmaRSig","sigmaRSig",6.1,0.0,100.);
  RooRealVar sigLSig("sigmaLSig","sigmaLSig",9.7,0.0,100.);
  RooRealVar meanBifSig("meanBifSig","meanBifSig",36.,0.0,100.0);
  RooBifurGauss BifGaussSig("BifGausSig","BifGausSig",METPara,meanBifSig,sigLSig, sigRSig);
  RooRealVar sigSig("sigmaSig","sigmaSig",22.,0.0,100.);
  RooRealVar meanSig("meanSig","meanSig",35.5,0.0,100.);
  RooRealVar alphaSig("alphaSig","alphaSig",1., -100., 100.); 
  RooRealVar nSig("nSig","nSig",5); 


  RooCBShape METCBSignal("METCBSignal","METCBSignal",METPara,meanSig,sigSig,alphaSig,nSig);
  RooAddPdf METParaSignal("METSignal","METSignal",BifGaussSig,METCBSignal ,coefSig);
  
  if(specy=="W_tn"){
    
    coefSig.SetName("coefSigTau");
    sigRSig.SetName("sigmaRSigTau");
    sigLSig.SetName("sigmaLSigTau");
    meanBifSig.SetName("meanBifSigTau");
    sigSig.SetName("sigmaSigTau");
    meanSig.SetName("meanSigTau");
    alphaSig.SetName("alphaSigTau");
    nSig.SetName("nSigTau");
    BifGaussSig.SetName("BifGaussSigTau");
    METCBSignal.SetName("METCBSignalTau");
    METParaSignal.SetName("METSignalTau");
        
  }
  if(specy=="Z_2e"){

    sigRBkg.SetName("sigmaRZ2E");
    sigLBkg.SetName("sigmaLZ2E");
    meanBifBkg.SetName("meanBifZ2E");
    sigBkg.SetName("sigmaZ2E");
    meanBkg.SetName("meanZ2E");
    coefBkg.SetName("coefZ2E");
  }

  RooAbsPdf *pdfFit;
  if(specy=="W_en"|| specy=="W_tn"){
    pdfFit= (RooAbsPdf*) &METParaSignal;
  } else{
    pdfFit= (RooAbsPdf*) &METParaBkg;
  }

  TFile fOut( (_outputpath +"/"+"PDFMETPara_"+specy+"_"+cutName+".root").c_str(),"RECREATE" );  

  pdfFit->fitTo(DataSet,RooFit::PrintLevel(3),RooFit::SumW2Error(kTRUE), RooFit::Save(),RooFit::Verbose(kFALSE));
   
  coefSig.setConstant();
  sigRSig.setConstant();
  sigLSig.setConstant();
  meanBifSig.setConstant();
  sigSig.setConstant();
  meanSig.setConstant();
  alphaSig.setConstant();
  nSig.setConstant();

  if(specy=="Z_2e"){
    sigRBkg.setConstant();
    sigLBkg.setConstant();
    meanBifBkg.setConstant();
    sigBkg.setConstant();
    meanBkg.setConstant();
    coefBkg.setConstant();
    
  }
  double nsigma=2.0;

  sigRSig.setRange(sigRSig.getVal()-nsigma*sigRSig.getError(),sigRSig.getVal()+nsigma*sigRSig.getError());
  sigLSig.setRange(sigLSig.getVal()-nsigma*sigLSig.getError(),sigLSig.getVal()+nsigma*sigLSig.getError());
  meanBifSig.setRange(meanBifSig.getVal()-nsigma*meanBifSig.getError(),meanBifSig.getVal()+nsigma*meanBifSig.getError());
  sigSig.setRange(sigSig.getVal()-nsigma*sigSig.getError(),sigSig.getVal()+nsigma*sigSig.getError());
  meanSig.setRange(meanSig.getVal()-nsigma*meanSig.getError(),meanSig.getVal()+nsigma*meanSig.getError());
  alphaSig.setRange(alphaSig.getVal()-nsigma*alphaSig.getError(),alphaSig.getVal()+nsigma*alphaSig.getError());

  sigRBkg.setRange(sigRBkg.getVal()-nsigma*sigRBkg.getError(),sigRBkg.getVal()+nsigma*sigRBkg.getError());
  sigLBkg.setRange(sigLBkg.getVal()-nsigma*sigLBkg.getError(),sigLBkg.getVal()+nsigma*sigLBkg.getError());
  meanBifBkg.setRange(meanBifBkg.getVal()-0.5*nsigma*meanBifBkg.getError(),meanBifBkg.getVal()+0.5*nsigma*meanBifBkg.getError());
  sigBkg.setRange(sigBkg.getVal()-nsigma*sigBkg.getError(),sigBkg.getVal()+nsigma*sigBkg.getError());
  meanBkg.setRange(meanBkg.getVal()-0.5*nsigma*meanBkg.getError(),meanBkg.getVal()+0.5*nsigma*meanBkg.getError());
  
  string nameCan=( "can"+ specy);

  TCanvas *can = new TCanvas(nameCan.c_str(),nameCan.c_str());
  RooPlot* metframe = METPara.frame();
  DataSet.plotOn(metframe);
  pdfFit->plotOn(metframe); 
  metframe->SetTitle(specy.c_str());
  metframe->Draw();
  
  can->Write();
  fOut.Close();

  RooArgSet* params;
  params = pdfFit->getParameters( RooDataSet() );
  params->writeToFile((_outputpath +"/"+ "PDFMETPara_"+specy+"_"+cutName+".txt").c_str());
 
}


void PrimaryAnalysis::fitAllMETPara(const string & specy){
 
  int percent[3];
  bool doID[2];  
  bool doIso[2];
  bool doEmaAndNoInit[2];

  percent[0]=95;
  percent[1]=80;
  percent[2]=70; 
  doID[0]=false;
  doID[1]=true;
  doIso[0]=false;
  doIso[1]=true;
  doEmaAndNoInit[0]=false;
  doEmaAndNoInit[1]=true;

  for (int ipercent=0;ipercent<3;ipercent++){
    for (int jinit=0;jinit<2;jinit++){
      for (int jid=0;jid<2;jid++){ 
	for (int jiso=0;jiso<2;jiso++){
	  string cutName=CutName(percent[ipercent], doIso[jiso], doID[jid],doEmaAndNoInit[jinit],false,false);
	  TCut cut= GetCut(percent[ipercent], doIso[jiso], doID[jid],doEmaAndNoInit[jinit], false, false, false);
	  fitMETPara( specy, cut,  cutName);	  
	}
      }
    }   
  }
}


void PrimaryAnalysis::fitAllFlipCut(){
 
  cout<<"==============================="<< endl;
  cout<<"Fitting PDF with FlipCuts ..."<< endl;


  string specy="FullLumi";
  
  int percent[2];
  bool doEmaAndNoInit[2];
  
  bool doFlipIso[2];
  bool doFlipID[2];

  doFlipIso[0]=false;
  doFlipIso[1]=true;
  doFlipID[0]=false;
  doFlipID[1]=true;

  bool doID;  
  bool doIso;

  percent[0]=95;
  percent[1]=80;

  doEmaAndNoInit[0]=false;
  doEmaAndNoInit[1]=true;
  string cutName;TCut cut;
  for (int ipercent=0;ipercent<2;ipercent++){
    for (int jinit=0;jinit<2;jinit++){
      for (int jid=0;jid<2;jid++){ 
	for (int jiso=0;jiso<2;jiso++){

	  if(doFlipIso[jiso]==true && doFlipID[jid]==true){
	    doID=false; doIso=false;
	  }else if (doFlipIso[jiso]==true){
	    doID=true; doIso=false;
	  }else if (doFlipID[jid]==true){
	    doIso=true; doID=false;
	  }else{
	    continue;
	  }
	  
	  string cutName=CutName(percent[ipercent], doIso, doID, doEmaAndNoInit[jinit],doFlipIso[jiso],doFlipID[jid]);
	  TCut cut= GetCut(percent[ipercent], doIso, doID, doEmaAndNoInit[jinit],doFlipIso[jiso],doFlipID[jid],false);
	  fitMETPara( specy, cut,  cutName);	  
	}
      }
    }   
  }

  
  specy="W_en";
  cutName=CutName(95, true, false, false, false, true);
  cut= GetCut(95, true, false, false, false, true,false);
  fitMETPara( specy, cut,  cutName);	  
  specy="W_tn";
  fitMETPara( specy, cut,  cutName);	  
  specy="Z_2e";
  fitMETPara( specy, cut,  cutName);

  cout<<" ... Done "<< endl;

}

void PrimaryAnalysis::fitAll(){

  cout<<"==============================="<< endl;
  cout<<"Fitting PDF for all species ..."<< endl;

  string species[6];
  species[0]="W_en";
  species[1]="W_tn";
  species[2]="Z_2e";
  species[3]="QCD";
  species[4]="GamJet";
  for(int i=0;i<5;i++){
    fitAllMETPara(species[i]);
    cout<<"   "<< species[i]<<" done"<< endl;
  }
  cout<<" ... Done "<< endl;

}



