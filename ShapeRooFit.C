#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>

#include <TROOT.h>
#include <TTree.h>
#include <TChain.h>
#include <TFile.h>
#include <TMath.h>
#include <TH1F.h>
#include <TF1.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TVector2.h>
#include <TString.h>
#include <TLorentzVector.h>
#include <TMatrix.h>
#include <TLatex.h>

#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooFitResult.h>
#include <RooGaussian.h>
#include <RooNovosibirsk.h>
#include <RooExponential.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include <TStyle.h>
#include <RooVoigtian.h>
#include <RooFormulaVar.h>

#include <RooMsgService.h>


using namespace std;
using namespace RooFit;
vector<vector<float> > vf(3,vector<float>(2,0));

bool isNotIso=false;
float tauPtC=30;
float lepPtC=30;

float massC=0;

TString fName;

int NBin=50;

bool gen=false;

void SetVars(float tpT, float leppT, float mass, bool istauIso, TString file, int Nbin) {

  isNotIso = 1-istauIso;
  //tauPtC = tpT;
  lepPtC = leppT;

  fName = file;
  NBin = Nbin;
  massC= mass;

}

TTree* SkimTree() {

  //TFile* file=new TFile("/home/mmarionn/Documents/CMS/LeptoQuark/KKCode/ElTau/DYJets_tree.root","READ");

  TTree*  tree;
  if(fName=="MC") {
    TChain* chain = new TChain("mytree");
    chain->Add("/home/mmarionn/Documents/CMS/LeptoQuark/KKCode/ElTau/central/DYJets_tree.root");
    chain->Add("/home/mmarionn/Documents/CMS/LeptoQuark/KKCode/ElTau/central/WJets_tree.root");
    chain->Add("/home/mmarionn/Documents/CMS/LeptoQuark/KKCode/ElTau/central/TT_tree.root");
    chain->Add("/home/mmarionn/Documents/CMS/LeptoQuark/KKCode/ElTau/central/DIBOSON_tree.root");
    tree = (TTree*)chain;
  } 
  else if(fName=="DYComb") {
    TChain* chain = new TChain("mytree");
    chain->Add("/home/mmarionn/Documents/CMS/LeptoQuark/KKCode/ElTau/central/DYJets_tree.root");
    chain->Add("/home/mmarionn/Documents/CMS/LeptoQuark/MuTau/central/DYJets_tree.root");
    tree = (TTree*)chain;
  }
  else if(fName=="WComb") {
    TChain* chain = new TChain("mytree");
    chain->Add("/home/mmarionn/Documents/CMS/LeptoQuark/KKCode/ElTau/central/WJets_tree.root");
    chain->Add("/home/mmarionn/Documents/CMS/LeptoQuark/MuTau/central/WJets_tree.root");
    tree = (TTree*)chain;
  }
  else if (fName=="WZComb") {
    TChain* chain = new TChain("mytree");
    chain->Add("/home/mmarionn/Documents/CMS/LeptoQuark/KKCode/ElTau/central/WJets_tree.root");
    chain->Add("/home/mmarionn/Documents/CMS/LeptoQuark/MuTau/central/WJets_tree.root");
    chain->Add("/home/mmarionn/Documents/CMS/LeptoQuark/KKCode/ElTau/central/DYJets_tree.root");
    chain->Add("/home/mmarionn/Documents/CMS/LeptoQuark/MuTau/central/DYJets_tree.root");
    tree = (TTree*)chain;
  }
  else {
    TFile* file=new TFile( fName,"READ");
    tree =(TTree*)file->Get("mytree");
  }
  
  //Input Var
  float pt_lep;
  float pt_tau;
  float phi_lep;
  float phi_tau;
  int pass_trigger;
  float dz_lep;
  float d0_lep;
  int nbtagjets_tchel;
  int njets;
  int isoL_tau;
  int isoVL_tau;

  float iso_tau;

  float pt_jet1;
  float pt_jet1_1bjet;

  float pt_jet2;
  float pt_jet2_1bjet;

  float st_2jet;
  float st_1bjet;
  float st_2bjet;

  int charge_tau,charge_lep;

  tree->SetBranchAddress("pt_lep",&pt_lep);
  tree->SetBranchAddress("pt_tau",&pt_tau);
  tree->SetBranchAddress("phi_lep",&phi_lep);
  tree->SetBranchAddress("phi_tau",&phi_tau);
  tree->SetBranchAddress("iso_tau",&iso_tau);
  tree->SetBranchAddress("charge_lep",&charge_lep);
  tree->SetBranchAddress("charge_tau",&charge_tau);
  tree->SetBranchAddress("pass_trigger",&pass_trigger);
  tree->SetBranchAddress("dz_lep",&dz_lep);
  tree->SetBranchAddress("d0_lep",&d0_lep);
  tree->SetBranchAddress("nbtagjets_tchel",&nbtagjets_tchel);
  tree->SetBranchAddress("njets",&njets);
  tree->SetBranchAddress("isoL_tau",&isoL_tau);
  tree->SetBranchAddress("isoVL_tau",&isoVL_tau);
  tree->SetBranchAddress("pt_jet1",&pt_jet1);
  tree->SetBranchAddress("pt_jet1_1bjet",&pt_jet1_1bjet);
  tree->SetBranchAddress("pt_jet2",&pt_jet2);
  tree->SetBranchAddress("pt_jet2_1bjet",&pt_jet2_1bjet);

  tree->SetBranchAddress("st_2jet",&st_2jet);
  tree->SetBranchAddress("st_1bjet",&st_1bjet);
  tree->SetBranchAddress("st_2bjet",&st_2bjet);

  float m11,m12,m21,m22;
  float m1_reco_2bjet, m2_reco_2bjet;
  float m1_reco_2jet, m2_reco_2jet;
  tree->SetBranchAddress("m11",&m11);
  tree->SetBranchAddress("m12",&m12);
  tree->SetBranchAddress("m21",&m21);
  tree->SetBranchAddress("m22",&m22);
  tree->SetBranchAddress("m1_reco_2bjet",&m1_reco_2bjet);
  tree->SetBranchAddress("m2_reco_2bjet",&m2_reco_2bjet);
  tree->SetBranchAddress("m1_reco_2jet",&m1_reco_2jet);
  tree->SetBranchAddress("m2_reco_2jet",&m2_reco_2jet);

  float weight;
  tree->SetBranchAddress("Weight",&weight);

  float SFtrig;
  tree->SetBranchAddress("SFtrig",&SFtrig);

  float Wtag_2bjet;
  tree->SetBranchAddress("Wtag_2bjet",&Wtag_2bjet);

 //gen level
  float glep_pt,glep_eta,glep_phi;
  float gtau_pt,gtau_eta,gtau_phi;
  float gjet1_pt,gjet1_eta,gjet1_phi;
  float gjet2_pt,gjet2_eta,gjet2_phi;

  int gtau_pdgId, glep_pdgId, gjet1_pdgId, gjet2_pdgId;

  int comb;
  
  if(gen) {
    tree->SetBranchAddress("glep_pt",&glep_pt);
    tree->SetBranchAddress("glep_eta",&glep_eta);
    tree->SetBranchAddress("glep_phi",&glep_phi);

    tree->SetBranchAddress("gtau_pt",&gtau_pt);
    tree->SetBranchAddress("gtau_eta",&gtau_eta);
    tree->SetBranchAddress("gtau_phi",&gtau_phi);
 
    tree->SetBranchAddress("gjet1_pt",&gjet1_pt);
    tree->SetBranchAddress("gjet1_eta",&gjet1_eta);
    tree->SetBranchAddress("gjet1_phi",&gjet1_phi);
    tree->SetBranchAddress("gjet2_pt",&gjet2_pt);
    tree->SetBranchAddress("gjet2_eta",&gjet2_eta);
    tree->SetBranchAddress("gjet2_phi",&gjet2_phi);

    tree->SetBranchAddress("gtau_pdgId",&gtau_pdgId);
    tree->SetBranchAddress("glep_pdgId",&glep_pdgId);
    tree->SetBranchAddress("gjet1_pdgId",&gjet1_pdgId);
    tree->SetBranchAddress("gjet2_pdgId",&gjet2_pdgId);
   
    tree->SetBranchAddress("comb",&comb);
 }

  TLorentzVector tau(0,0,0,0),lep(0,0,0,0),jet1(0,0,0,0),jet2(0,0,0,0);
  TLorentzVector LQt(0,0,0,0);
  TLorentzVector LQl(0,0,0,0);
  int ST=0;

  //Weights
  TFile *wf=new TFile("../Weights.root","READ");
  TH1F* W1=(TH1F*)wf->Get("W1");
  TH1F* W2=(TH1F*)wf->Get("W2");
  vector<double> ratio1B,ratio2B;
  vector<double> ratio1Be,ratio2Be;
  for(int i=0;i<51;i++) {
    ratio1B.push_back( W1->GetBinContent(i+1) );
    ratio2B.push_back( W2->GetBinContent(i+1) );
    ratio1Be.push_back( W1->GetBinError(i+1) );
    ratio2Be.push_back( W2->GetBinError(i+1) );
  }

  //  OutputVar
  float o_ST;
  float o_tau_pt;
  float o_weight;
  float o_mLQ;
  float o_isot;

  TTree* oTree=new TTree("tree","tree");
  oTree->SetDirectory( 0);
  oTree->Branch("ST",&o_ST,"o_ST/F");
  oTree->Branch("pttau",&o_tau_pt,"tau_pt/F");
  oTree->Branch("weight",&o_weight,"o_weight/F");
  oTree->Branch("mLQ",&o_mLQ,"o_mLQ/F");
  //  oTree->Branch("isoT",&o_isot,"o_isot/F");

  //Now read the tree; 


  TVector2 t(0,0);
  TVector2 l(0,0);
  TVector2 tl(0,0);
  float w=1;

  int Nent = tree->GetEntries();
  int N=0;
  cout<< Nent<<endl;

  for(int ie=0;ie<Nent;ie++) {

    tree->GetEntry(ie);

    if( pt_lep <= lepPtC)  continue;
    if( pt_tau <= 30)  continue; //minimum cont
    if(pass_trigger != 1 ) continue;
    if(abs(dz_lep)>=0.2) continue;
    if(abs(d0_lep)>=0.045) continue;
    if(charge_tau*charge_lep!=-1) continue;
  
   
    if(njets<2) continue; //ok inf or equal to 1
    if(isoL_tau == isNotIso ) continue;
    //  if(!isoVL_tau ) continue;
   
    //Z/W+jet
    {
      if(nbtagjets_tchel!=0) continue;

      int bin = (int)(pt_jet1/10)+1;
      int bin2 = (int)(pt_jet2/10)+1;
      if(bin>51) bin =51;
      if(bin>51) bin =51;
      w = ratio1B[bin]*ratio2B[bin2]*weight;//*SFtrig*Wtag_2bjet
    }
    //TTbar
    // if(nbtagjets_tchel<0) continue;
    // if(fName=="MC")
    //   w = weight*SFtrig*Wtag_2bjet;
    // else
    //   w=1.;
    
    //if(m2_reco_2jet<massC) continue;
    
    t.SetMagPhi( pt_tau, phi_tau );
    l.SetMagPhi( pt_lep, phi_lep );
    
    tl = t+l;
    
    
   
    
    // float w1 = ratio1B[bin];
    // float w2 = ratio2B[bin2];

    // float w1e = ratio1Be[bin];
    // float w2e = ratio2Be[bin2];
    
    //ok, all values computed, now fill the output tree
    // cout<<st_2jet<<"   "<<st_2bjet<<"   "<<st_1bjet<<endl;
    o_ST = st_2jet;
    o_tau_pt = pt_tau;
    o_weight = w;
    o_mLQ = m2_reco_2jet;
    o_isot = iso_tau/pt_tau;

    //  cout<<o_isot<<endl;

    //Now generator level
    if(gen) {
      tau.SetPtEtaPhiM(gtau_pt,gtau_eta,gtau_phi,0.);
      lep.SetPtEtaPhiM(glep_pt,glep_eta,glep_phi,0.);
      jet1.SetPtEtaPhiM(gjet1_pt,gjet1_eta,gjet1_phi,0.);
      jet2.SetPtEtaPhiM(gjet2_pt,gjet2_eta,gjet2_phi,0.);
      if(comb==1) {
	LQt = tau + jet2;
	LQl = tau + jet1;
      }
      if(comb==2) {
	LQt = tau + jet1;
	LQl = tau + jet2;
      }  
      if(comb==-1) abort();
    
      ST = tau.Pt() + lep.Pt() + jet1.Pt() + jet2.Pt();


      //remake the pair
      if( fabs((tau+jet2).M()-(lep+jet1).M()) < fabs((lep+jet2).M()-(tau+jet1).M()) ) {
	LQt = tau + jet2;
	LQl = tau + jet1;
      }
      else {
	LQt = tau + jet1;
	LQl = tau + jet2;
      }

      o_ST = ST;
      o_tau_pt = pt_tau;
      o_weight = w;
      o_mLQ = LQt.M();
    }


    oTree->Fill();
    N++;
  }

  cout<<" Nevent "<<N<<endl;

  return oTree;

  

}




void MakeFit(TTree* tree,float taucut=30) {

  //Variables
  RooRealVar weight("weight","weight",-10,100);
  RooRealVar ST("ST","ST",0,1000);
  RooRealVar pttau("pttau","pttau",taucut,1000);

  //dataset
  RooDataSet data("data","data",tree,RooArgSet(weight,ST,pttau),"","weight");
  

  //Novosibirsk function==========
  RooRealVar peak("#mu","#mu",180,100,500);
  RooRealVar width("#sigma","#sigma",65,10,200);
  RooRealVar tail("#tau","#tau",0.5,0.05,1.5);
  
  RooNovosibirsk novo("novo","novo",ST,peak,width,tail);
  

  RooFitResult* result = novo.fitTo( data , RooFit::SumW2Error(kFALSE), RooFit::Save(kTRUE), RooFit::PrintLevel(-1) );

  // TCanvas * c = new TCanvas("cR","cR");

  RooPlot* frame = ST.frame() ;
  data.plotOn(frame);
  novo.plotOn(frame);
  novo.paramOn(frame);
  frame->Draw();

  cout<<" youpali "<<endl;

  //Get back the parameters
  // vector<vector<float> > vf(3,vector<float>(2,0));
  vf[0][0] = (float)(peak.getVal());  vf[0][1] = (float)peak.getError();
  vf[1][0] = (float)width.getVal();  vf[1][1] = (float)width.getError();
  vf[2][0] = (float)tail.getVal();  vf[2][1] = (float)tail.getError();
  cout<<" youpal1 "<<endl;
  // return vf;

}


void SingleFit(float taucut=30) {

  TTree* tree=SkimTree();
  MakeFit(tree, taucut);
 
}



void MultipleFits() {

  TTree* tree=SkimTree();
  
  int NP=4;

  vector<TGraphErrors*> vgraph;
  for(int i=0;i<3;i++) {
    TGraphErrors* tmpg=new TGraphErrors(NP);
    vgraph.push_back( tmpg );
  }

  TCanvas * cc=new TCanvas("cc","cc",1500,1000);
  cc->Divide(3,2);
   //TCanvas * cc=new TCanvas("cc","cc");

  // RooPlot* frame;
  int col[6]={kBlue,kGreen+1,kOrange-3,kRed+1,kViolet+7,kViolet+3};

  for(int iv=0;iv<NP;iv++) {
   
    float thr= 30+iv*5.;
    //float thr= 100+iv*10.;
    // float thr = 0.20-iv*0.04;

    //Variables
    TString ss="ST";
    ss+="(M>";
    ss+=thr;
    ss+=") (GeV)";

    RooRealVar weight("weight","weight",-10,100);
    RooRealVar ST("ST",ss,0,1000);
    RooRealVar pttau("pttau","pttau",thr,1000);
    RooRealVar mLQ("mLQ","mLQ",massC,1000);
    //    RooRealVar iso("isoT","isoT",-0.10,thr);

    //dataset
    RooDataSet data("data","data",tree,RooArgSet(weight,ST,pttau,mLQ),"","weight");
  
    //Novosibirsk function==========
    float mum=360;
    float sigmam=120;
    float tailm = 0.5;
    // if(massC> 100) {
    //   ST.setBins(NBin);
    //   mum = 400;
    //   sigmam = 100;
    //   tailm = 0.2;
    // }
    // if(iv!=0) {
    //   mum = vf[0][0];
    //   sigmam = vf[1][0];
    //   tailm = vf[2][0];
    // }
    ST.setBins(NBin);

    RooRealVar peak("#mu","#mu",mum,100,500);
    RooRealVar width("#sigma","#sigma",sigmam,10,200);
    RooRealVar tail("#tau","#tau",tailm,-0.20,1.5);
  
    RooNovosibirsk novo("novo","novo",ST,peak,width,tail);
  
    RooRealVar shift("shift","shift",230,150,400);
    RooFormulaVar ST_shift("STs","@0-@1",RooArgSet(ST,shift));
    RooRealVar slope("slope","slope",-0.007,-10,10);
    RooExponential expo("expo","expo",ST,slope);

    RooFitResult* result2 = novo.fitTo( data , RooFit::SumW2Error(kFALSE), RooFit::Save(kTRUE), RooFit::PrintLevel(-1) );
    // TCanvas * c = new TCanvas("cR","cR");

    RooFitResult* result = expo.fitTo( data , RooFit::SumW2Error(kFALSE), RooFit::Save(kTRUE), RooFit::Range(peak.getVal()+width.getVal()/3.,1000));//, RooFit::PrintLevel(-1) );

    cc->cd(iv+1);
    // if(iv==0) {
    RooPlot* frame = ST.frame() ;
      data.plotOn(frame);
      //}
    expo.plotOn(frame,LineColor(col[iv]));
    expo.paramOn(frame);
    novo.plotOn(frame,LineColor(col[iv]), LineStyle(iv+1));
    novo.paramOn(frame);
    // if(iv==NP-1)
      frame->Draw();

    vf[0][0] = (float)(peak.getVal());  vf[0][1] = (float)peak.getError();
    vf[1][0] = (float)width.getVal();  vf[1][1] = (float)width.getError();
    vf[2][0] = (float)tail.getVal();  vf[2][1] = (float)tail.getError();

    //  if(iv!=1)// && iv!=14)
    for(int j=0;j<3;j++) {
      cout<<vf[j][0]<<endl;
      vgraph[j]->SetPoint(iv,thr,vf[j][0]);
      vgraph[j]->SetPointError(iv,0,vf[j][1]);
    }
  }

  vector<TF1*> laws;
  gStyle->SetOptFit(0);
  TCanvas* c2= new TCanvas("c","c",1500,450);
  c2->Divide(3);

  string title[3]={"mean [GeV]","#sigma [GeV]","#tau  "};

  for(int i=0;i<3;i++) {
    c2->cd(i+1);
    vgraph[i]->Draw("AP");
    vgraph[i]->GetYaxis()->SetTitle(title[i].c_str() );
    vgraph[i]->GetXaxis()->SetTitle("p_{T}(#tau) thr [GeV]");
    TString p="p";
    p+=i;
    TF1* f = new TF1(p,"pol1",30,161);
   
    if(i==2)
      f = new TF1(p,"pol1",30,161);
    vgraph[i]->Fit(p,"R");
    laws.push_back(f);
  }


  //For Keti !!! below this step ==========================================================================



  //Check
  //Variables
  RooRealVar weight("weight","weight",-10,100);
  RooRealVar ST("ST","ST [GeV]",0,1000);
  RooRealVar pttau("pttau","pttau",50,1000); //50
  RooRealVar mLQ("mLQ","mLQ",massC,1000); //massC
  // RooRealVar iso("isoT","isoT",-0.10,0.15); //massC
  //dataset
  // TString Wstr="";
  // if(fName=="/home/mmarionn/Documents/CMS/LeptoQuark/KKCode/ElTau/central/DYJets_tree.root") Wstr="weight";

  RooDataSet data("data","data",tree,RooArgSet(weight,ST,pttau,mLQ),"","weight");//weight,,"","weight"

  ST.setBins(NBin);

  //Novosibirsk function==========
  RooRealVar peak("#mu","#mu",laws[0]->Eval(50),100,500);
  RooRealVar width("#sigma","#sigma",laws[1]->Eval(50),10,200);
  RooRealVar tail("#tau","#tau",laws[2]->Eval(50),-0.20,1.5);
  
  RooNovosibirsk novo("novo","novo",ST,peak,width,tail);

  RooFitResult* result = novo.fitTo( data , RooFit::SumW2Error(kFALSE), RooFit::Save(kTRUE));//, RooFit::PrintLevel(-1) 

  TMatrixDSym corrM = result->correlationMatrix();
  corrM.Print();
  //cout<<result->correlationMatrix()<<endl;

 //Novosibirsk function of reference==========
  RooRealVar peak_ref("#mu_ref","#mu",laws[0]->Eval(50)); //50
  RooRealVar width_ref("#sigma_ref","#sigma",laws[1]->Eval(50));
  RooRealVar tail_ref("#tau_ref","#tau",laws[2]->Eval(50));
  
  RooNovosibirsk novo_ref("novo_ref","novo_ref",ST,peak_ref,width_ref,tail_ref);

  //Expofit now
  
  RooRealVar shift("shift","shift",230,150,400);
  RooFormulaVar ST_shift("STs","@0-@1",RooArgSet(ST,shift));
  RooRealVar slope("slope","slope",-0.007,-10,10);
  RooExponential expo("expo","expo",ST,slope);
  RooFitResult* result2 = expo.fitTo( data , RooFit::SumW2Error(kFALSE), RooFit::Save(kTRUE), RooFit::Range(peak.getVal()+width.getVal()/3.,1000), RooFit::PrintLevel(-1) );


  RooRealVar n_Ex( "n_Ex","n_Ex",1.);
  RooRealVar n_No( "n_No","n_No",-1.);
  RooArgList listPdfVal( n_Ex, n_No );
  RooAddPdf diff("diff","diff",RooArgList(expo,novo),listPdfVal);



  TCanvas * c = new TCanvas("cR","cR");

  RooPlot* frameN = ST.frame() ;
  data.plotOn(frameN,Name("dati"));
  novo.plotOn(frameN,Range(0.0,1000),Name("novog"));
  //novo_ref.plotOn(frameN,Range(0.0,1000),LineColor(kRed),LineStyle(kDashed));
  //novo.paramOn(frameN );
  novo.plotOn(frameN,VisualizeError(*result,1,kFALSE),DrawOption("L"),LineWidth(2),LineColor(kRed),LineStyle(kDashed),Range(0.0,1000),Name("novoeg") );
  //expo.plotOn(frameN,/*LineStyle(kDashed),*/LineColor(kGreen+1),Range(500,1000),Name("expog"));
  //expo.paramOn(frameN);
  //  diff.plotOn(frameN,Range(500,1000),LineColor(kViolet+3));
  //expo.plotOn(frameN,VisualizeError(*result,1,kFALSE),DrawOption("L"),LineWidth(2),LineColor(kGreen),Range(500,1000) );
  // novo_ref.plotOn(frameN,LineStyle(kDashed),LineColor(kRed));
  frameN->Draw();

  //  frameN->Print("v");

  //==========================================================

  // cout<<frameN->chiSquare("novog","dati",3)<<endl;

  TGraph* novoPart = (TGraph*)(frameN->getCurve("novog"));
  TGraph* novoErrorPart = (TGraph*)(frameN->getCurve("novoeg"));
  TGraph* expoPart = (TGraph*)(frameN->getCurve("expog"));
  
  cout<<novoPart->GetN()<<"   "<<novoErrorPart->GetN()<<"   "<<expoPart->GetN()<<"   "<<endl;

  double x,y,xl,yl,xh,yh;
  int n=0;

  int NN=novoPart->GetN();

  TGraph* gUp=new TGraph(NN);
  TGraph* gLow=new TGraph(NN);

  vector<double> ratios;
  for(int i=0;i<NN;i++) {
    
    novoPart->GetPoint(i,x,y);
    //  cout<<i<<"  "<<x<<"   "<<y<<endl;
    novoErrorPart->GetPoint(i,xl,yl);
    novoErrorPart->GetPoint(i+NN,xh,yh);

    cout<<x<<" ---> "<<yl/y<<" ;  "<<yh/y<<"     =====       "<<y<<" - "<<yl<<" + "<<yh<<endl;

    gUp->SetPoint(NN-1-i,x,yh);
    gLow->SetPoint(i,x,yl);

    // if(x>=500) {
    //   expoPart->GetPoint(n,x2,y2);
    //   cout<<" --->  "<<i<<"  "<<x<<"   "<<y<<endl;
     
    //   //ratios
    //   n++;
    // }
  }

  cout<<gUp->Integral()/(1000./NBin)<<"   "<<gLow->Integral()/(1000./NBin)<<"   "<<novoPart->Integral()/(1000./NBin)<<"   "<<data.sumEntries()<<endl;

  // float corYieldUp=gUp->Integral()/novoPart->Integral();
  // float corYieldlow=gLow->Integral()/novoPart->Integral();

  float intUp = gUp->Integral()/(1000./NBin);
  float intCentral = novoPart->Integral()/(1000./NBin);
  float intLow = gLow->Integral()/(1000./NBin);

  TH1F* hUp=new TH1F("hUp","hUp",10,100,1100);
  TH1F* hCentral=new TH1F("hCentral","hCentral",10,100,1100);
  TH1F* hLow=new TH1F("hLow","hLow",10,100,1100);

  cout<<intUp<<"   "<<intLow<<endl;
  double xc,yc,xu,yu,xl2,yl2;

  for(int i=0;i<NN;i++) {
    
    gUp->GetPoint(i,xu,yu);
    novoPart->GetPoint(i,xc,yc);
    gLow->GetPoint(i,xl2,yl2);

    for(int j=0;j<11;j++) {
      if( fabs(xc-(j*100+50.))<0.1 ) {
  	cout<<xc<<"  "<<j<<endl;
  	hUp->SetBinContent(j,yu);
  	hCentral->SetBinContent(j,yc);
  	hLow->SetBinContent(j,yl2);
      }

    }
    
  }

  TFile* ofile=new TFile("ofile.root","RECREATE");
  hUp->Write();
  hCentral->Write();
  hLow->Write();
  ofile->Close();

 //For Keti !!! above this step ==========================================================================


  //====================================================

  //For comparison
 //  fName="MC";

 //  TTree* tree2=SkimTree();
 //  RooDataSet data2("data2","data2",tree2,RooArgSet(ST,pttau,mLQ,weight),"","weight" );
 // //Novosibirsk function==========
 //  RooRealVar peak2("#mu_{d}","#mu_{d}",230,100,500);
 //  RooRealVar width2("#sigma_{d}","#sigma_{d}",63,10,200);
 //  RooRealVar tail2("#tau_{d}","#tau_{d}",0.6,0.05,1.5);
  
 //  RooNovosibirsk novo2("novo2","novo_data",ST,peak2,width2,tail2);

 //  RooFitResult* result2d = novo2.fitTo( data2 , RooFit::SumW2Error(kFALSE), RooFit::Save(kTRUE)  );
 //  //cout<<result2d->correlationMatrix()<<endl;

  
 //  //Expofit now
 //  RooRealVar slope2("slope2","slope2",-0.007,-10,10);
 //  RooExponential expo2("expo2","expo2",ST,slope2);
 //  RooFitResult* result22 = expo2.fitTo( data2 , RooFit::SumW2Error(kFALSE), RooFit::Save(kTRUE), RooFit::Range(peak2.getVal()+width2.getVal()/3.,1000), RooFit::PrintLevel(-1) );
 //  //data2.plotOn(frameN,LineColor(kRed),PointColor(kRed));
 //  novo2.plotOn(frameN,LineStyle(kDashed),LineColor(kRed));
 //  //  novo2.paramOn(frameN);
 //  // expo2.plotOn(frameN,LineStyle(kDashed),LineColor(kRed));
 //  // expo2.paramOn(frameN);
 //  frameN->Draw();

 // // TCanvas * cdd = new TCanvas("cRdd","cR");

 //  RooPlot* dummyframe = ST.frame() ;
 //  data2.plotOn(dummyframe,Name("datiti"));
 //  novo2.plotOn(dummyframe,Range(0.0,1000),Name("novog2"));
 // //  dummyframe->Draw();

 //  c->cd();

  // double mu_mc = peak2.getVal();double emu_mc = peak2.getError();
  // double sig_mc = width2.getVal();double esig_mc = width2.getError();
  // double tau_mc = tail2.getVal();double etau_mc = tail2.getError();
  // double chi2_mc = dummyframe->chiSquare("novog2","datiti",3);

   double mu_d = peak.getVal();double emu_d = peak.getError();
   double sig_d = width.getVal();double esig_d = width.getError();
   double tau_d = tail.getVal();double etau_d = tail.getError();
   double chi2_d = frameN->chiSquare("novog","dati",3);

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.03);
  // latex.DrawLatex(0.624,0.600,Form("#mu_{MC} = %.2f #pm %.2f GeV", mu_mc, emu_mc) ); 
  // latex.DrawLatex(0.624,0.541,Form("#sigma_{MC} = %.2f #pm %.2f GeV", sig_mc, esig_mc ) );
  // latex.DrawLatex(0.624,0.482,Form("#tau_{MC} = %.2f #pm %.2f ", tau_mc, etau_mc ) );
  // latex.DrawLatex(0.624,0.423,Form("#chi^{2}/ndof (MC) = %.2f", chi2_mc) );


  latex.DrawLatex(0.614,0.893,Form("#mu = %.2f #pm %.2f GeV", mu_d, emu_d) );//  _{data}
   latex.DrawLatex(0.614,0.834,Form("#sigma = %.2f #pm %.2f GeV", sig_d, esig_d ) );//  _{data}
   latex.DrawLatex(0.614,0.765,Form("#tau = %.2f #pm %.2f ", tau_d, etau_d ) );//  _{data}
   // latex.DrawLatex(0.614,0.706,Form("#chi^{2}/ndof = %.2f", chi2_d) ); //(data)


  latex.Draw();

  
  TLegend* leg = new TLegend(0.607,0.508,0.919,0.656);
  leg->SetFillColor(0); leg->SetLineColor(1); leg->SetShadowColor(0);
  
  TLine* l1=new TLine(1,1,1,1);
  l1->SetLineColor(kBlue+1);l1->SetLineWidth(2);
  TLine* l2=new TLine(1,1,1,1);
  l2->SetLineColor(kRed);l2->SetLineWidth(2);l2->SetLineStyle(2);
  
  leg->AddEntry(novoPart,"W+jets","p");
  // leg->AddEntry(l1,"S_{T} shape (data)","l");
  // leg->AddEntry(l2,"S_{T} shape (MC)","l");
  leg->AddEntry(l1,"central fit","l");
  // leg->AddEntry(l2,"extrapolation ","l");
  leg->AddEntry(l2,"1 #sigma variation","l");
  leg->Draw();
}
