void plotDataPlusFit( TString analysisName, TString sampleName, TString subSampleName,  TString histName, TString variableName, float xmin, float xmax, TString xtitle, TString ytitle, bool log=true )
{
  TString key = sampleName + "__" + analysisName + "__" + variableName + "__" + histName;
  TString key = 
    analysisName 
    + "__" 
    + sampleName 
    + "__" 
    + subSampleName 
    + "__" 
    + histName
    + "__" 
    + variableName 
    ;

  TString name_(key);
  if( log ) name_ += "_log";
  else      name_ += "_lin";

  TCanvas *c1 = new TCanvas(name_, name_, 5, 49, 600, 600);
  c1->SetLeftMargin(0.18);
  c1->SetRightMargin(0.04);

  TString filename= TString("hist/") + key + "__chi2.root";
  TFile* f = TFile::Open(filename);

  TString modelName = TString("pdf__")+key;

  TH1* model = (TH1*)((gROOT->FindObject(modelName))->Clone());
  TH1* data  = (TH1*)((gROOT->FindObject(histName))->Clone());

  data->SetAxisRange( 1.001*xmin, 0.999*xmax );

  model->GetXaxis()->SetTitle(xtitle);
  model->GetYaxis()->SetTitle(ytitle);
  model->SetLineColor( kRed );
  model->SetFillColor( kYellow );
  model->SetLineWidth( 2 );
  model->GetXaxis()->SetTitleOffset(1.2);
  model->GetXaxis()->SetLabelOffset(0.01);
  model->GetYaxis()->SetTitleOffset(2.2);
  model->GetYaxis()->SetLabelOffset(0.01);
    
  data->SetLineWidth( 2 );
  data->SetLineColor( kBlue );
  data->SetMarkerStyle( 20 );
  data->SetMarkerColor( kBlue );      
  
  // 
  c1->SetRightMargin(0.04);
  if( log )
    {
      c1->SetLogy();
      c1->SetGridy();
      //      model->SetMinimum(100);
      model->GetYaxis()->SetTitleOffset(1.7);
      data->SetMarkerSize( 0.5 );      
    }
  else
    {
      //      model->SetMaximum( 32000 );
      data->SetMarkerSize( 0.7 );      
    }

  c1->cd();
  data->Draw();
  model->Draw("LSame");
  data->Draw("ESame");

  //  f->Close();
  RooUtils::fixOverlay();
}

void plotZee( bool log=true )
{
  //  plotDataPlusFit( "Zee", "Zee", "mll", "Zee_mass","dilepton mass (GeV/c^{2 })", "events / (0.5 GeV/c^{2 })", log );

  // commented 20/11/09
  //  plotDataPlusFit( "Diboson", "Z_2e", "all", "stat__mll__0_all_nosel", "mll", 50, 120, "dilepton mass (GeV/c^{2 })", "events / (0.5 GeV/c^{2 })", log );

  plotDataPlusFit( "Zee", "Z_2e", "Z_2e", "mll__Zee", "mll", 50, 120, "dilepton mass (GeV/c^{2 })", "events / (0.5 GeV/c^{2 })", log );
}
