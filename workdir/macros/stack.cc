TCanvas*  ColorWeel()
{
  // Main color names
  // ================
  // kMagenta
  // kPink
  // kRed
  // kOrange
  // kYellow
  // kSpring
  // kGreen
  // kTeal
  // kCyan
  // kAzure
  // kBlue
  // kViolet
  // kGray
  // kWhite
  // kBlack
  TColorWheel *w = new TColorWheel();
  w->Draw();
  return w->GetCanvas();
}

void colors()
{
  ColorWeel()->Draw();
}

void fillHistFromSet( TH1* h_data, TFile* f_, TString sampleName, int lumi, TString varName )
{

  TString setName_ = "set_";
  setName_ += sampleName;
  setName_ += "_";
  setName_ += lumi;
  RooDataSet* set_ = (RooDataSet*)f_->Get(setName_);
  if( set_==0 ) continue;
  set_->Print();
  int jj(0);
  const RooArgSet* argset =  set_->get(jj);
  while( argset!=0 )
    {
      float val_ = ((RooAbsReal*)argset->find(varName.Data()) )->getVal();
      h_data->Fill( val_ );
      jj++;
      argset =  set_->get(jj);      
    }
}

void stack()
{
  //  ColorWeel();

  //  int fillColor[4] = { kYellow, kGreen, kCyan,  kBlue };
  //  int lineColor[4] = { kSpring, kTeal,  kAzure, kViolet };
//   int fillColor[10] = { kYellow-7, 
// 			kGreen-3, 
// 			kBlue-4,  k
// 			Red-3, 
// 			kYellow-7, 
// 			kGreen-3, 
// 			kBlue-4,  
// 			kRed-3, 
// 			kYellow+4, 
// 			kRed-5 };
//   int lineColor[10] = { kOrange-3, 
// 			kBlack,   
// 			kBlack,   
// 			kBlack,   
// 			kBlack,   
// 			kBlack,   
// 			kBlack,   
// 			kBlack,   
// 			kBlack ,   
// 			kBlack };
  int fillColor[50];  
  int lineColor[50];
  TH1*  h[50];

  // histogram file produced by src/runTest
  TString filename = "hist/stack_Zee_Zee_mass_100.root";
  TFile* f_ = TFile::Open( filename, "READ" );

  TString histName = "Zee_mass";
  //  TString histName = "mll_k11";

  

  int ipl=0;
  h[ipl++] = (TH1*)f_->Get("total_Zee_Wenu_WJet_GamJet_QCD")->Clone();
  bool reverse = true;
  if( reverse )
    {
      h[ipl++] = (TH1*)f_->Get("total_Zee_Wenu_WJet_GamJet")->Clone();
      h[ipl++] = (TH1*)f_->Get("total_Zee_Wenu_WJet")->Clone();
      h[ipl++] = (TH1*)f_->Get("total_Zee_Wenu")->Clone();
      h[ipl++] = (TH1*)f_->Get("total_Zee")->Clone();
    }
  else
    {
      h[ipl++] = (TH1*)f_->Get("total_Wenu_WJet_GamJet_QCD")->Clone();
      h[ipl++] = (TH1*)f_->Get("total_WJet_GamJet_QCD")->Clone();
      h[ipl++] = (TH1*)f_->Get("total_GamJet_QCD")->Clone();
      h[ipl++] = (TH1*)f_->Get("total_QCD")->Clone();
    }

  int n_    = ipl;
  cout << "Number of stacked histograms = " << ipl << endl;
  int scale_= 1;

  int mainColor_ = kRed;
  fillColor[0]    = mainColor_-7;
  lineColor[0]    = kBlack;
  fillColor[n_-1] = kYellow-7;
  lineColor[n_-1] = kBlack;
  for( int ii=1; ii<n_-1; ii++ )
    {
      fillColor[ii] = mainColor_-7+ii;
      lineColor[ii] = kBlack;
    }

  for ( int ii=0; ii<n_; ii++ )
    {
      //      TString name_ = TString("Stacked-")+histName+TString("_");
      //      name_ += (n_-ii-1);
      TH1* h_ = h[ii];
      TString name_ = h_ -> GetName();
      cout << name_ << endl;
      // name_ += ii;
      //      TH1* h_ = (TH1*)gROOT->FindObject(name_);
      //      TH1* h_ = (TH1*)f_->Get(name_)->Clone();
      h_->SetDirectory(0);
      h_->Scale( scale_ );
      h_->SetFillColor(fillColor[ii]);
      h_->SetLineColor(lineColor[ii]);
      //      h_->SetLineColor(kBlack);
      h_->SetLineWidth(1);
      if( ii==0 ) h_->SetLineWidth(2);
      if( ii==0 ) 
	{
	  h_->SetMinimum(0.2);
	  h_->Draw();
	}
      else
	{
	  h_->Draw("Same");
	}
      h[ii] = h_;
    }

   TAxis* xaxis = h[0]->GetXaxis();
   TH1* h_data = new TH1F( "h_data", "h_data",
 			  xaxis->GetNbins()/scale_,
 			  xaxis->GetXmin(), xaxis->GetXmax() );
   h_data->SetDirectory(0);

   TString sampleName[3] = {"QCD_EM2030","QCD_EM3080","QCD_EM80170" };
   TString varName="mll";
   int lumi = 100;

   for( int ii=0; ii<3; ii++ )
     {
       fillHistFromSet( h_data, f_, sampleName[ii], lumi, varName ); 
     }


   f_->Close();


   f_ = TFile::Open( "pseudoExps/exp_Zee_100_1.root", "READ" );
  if( f_!=0 )
    {
      int ii(0);
      const RooArgSet* argset =  dataset->get(ii);
      while( argset!=0 )
	{
	  float val_ = ((RooAbsReal*)argset->find("mll") )->getVal();
	  h_data->Fill( val_ );
	  ii++;
	  argset =  dataset->get(ii);      
	}
      
      f_->Close();
    }

  
  h_data->Draw("ESame");


  RooUtils::fixOverlay();
}

void checkToy( TString histName, TString sampleName, int lumi, TString varName )
{
  TString filename = "hist/stack_new_10000.root";

  // histogram file produced by src/runTest
  TFile* f_ = TFile::Open( filename, "READ" );

  histName += "_";
  histName += sampleName;
  histName += "_";
  histName += lumi;
  TH1* h_ = (TH1*)f_->Get(histName)->Clone();
  h_->SetDirectory(0);

  TAxis* xaxis = h_->GetXaxis();
  TH1* h_data = new TH1F( "h_data", "h_data",
			  xaxis->GetNbins(),
			  xaxis->GetXmin(), xaxis->GetXmax() );
  h_data->SetDirectory(0);
  fillHistFromSet( h_data, f_, sampleName, lumi, varName ); 

  f_->Close();

  h_->Draw();
  h_data->Draw("ESame");


  RooUtils::fixOverlay();
}
