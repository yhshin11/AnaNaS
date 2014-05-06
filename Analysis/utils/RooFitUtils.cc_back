#include "Analysis/utils/RooFitUtils.hh"

#include <TH1F.h>
#include <TCanvas.h>

#include <RooGlobalFunc.h>
#include <RooRealVar.h>
#include <RooArgList.h>
#include <RooFormula.h>
#include <RooDataSet.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include <RooFitResult.h>
#include <RooRealIntegral.h>

#include <iostream>
using namespace std;

void
RooFitUtils::testRooFit()
{
  RooRealVar x("x","x",-1,1);
  RooRealVar y("y","y",-1,1);
  RooArgList list( x, y );
  RooFormula f("formula","x<0.90&&x>-0.90&&y<0.90&&y>-0.90",list);
  
  x.setVal(0);
  y.setVal(0);
  cout << f.eval() << endl;
}

RooPlot*
RooFitUtils::addPdfPlot( const RooDataSet* theDataSet, const RooRealVar* theVar, 
		 RooAddPdf* thePdf, const char* theName, int nBins )
{

  RooPlot* theFrame = theVar->frame( nBins ) ;
  theFrame->SetName( theName );
  RooPlot* dataplot = theDataSet -> plotOn(theFrame);
  dataplot -> Print();
  TAttMarker* dataMarker = dataplot -> getAttMarker();
  dataMarker -> SetMarkerStyle(1);
  TAttFill* dataFill = dataplot -> getAttFill();
  dataFill -> SetFillStyle(1001);
  dataFill -> SetFillColor(0);

  if( thePdf==0 ) return theFrame;
  
  // line 
  //  RooPlot* pdfparam = thePdf->paramOn(theFrame, 
  //				      theDataSet, "", 2, 
  //				      "NELU", 0.15, 0.50 );

  // plot sub-pdfs
  const RooArgList& thePdfList = thePdf->pdfList();
  RooArgSet  subset( thePdfList );
  TIterator* pdfIter = thePdfList.createIterator() ;
  RooAbsPdf* pdf ;

  Int_t iColor = 2;

  // remove the top pdf
  pdf = (RooAbsPdf*) pdfIter->Next();
  subset.remove( *pdf );
  while( ( pdf = (RooAbsPdf*) pdfIter->Next() ) ) 
    {
      //      RooPlot* subpdfplot =  thePdf->plotOn(theFrame,subset) ;
      // TO BE REIMPLEMENTED !!!!!
      RooPlot* subpdfplot(0);
      TAttLine* subpdfLine =   subpdfplot -> getAttLine();
      //      TAttFill* subpdfFill =   subpdfplot -> getAttFill();
      //      TAttText* subpdfText =   subpdfplot -> getAttText();
      
      subpdfLine->SetLineWidth(1);
      subpdfLine->SetLineColor(++iColor);
      subpdfLine->SetLineStyle(1);
     
      subset.remove( *pdf );
      
    }


  RooPlot*  pdfplot =   thePdf->plotOn(theFrame) ;
  TAttLine* pdfLine =   pdfplot -> getAttLine();
  //  TAttFill* pdfFill =   pdfplot -> getAttFill();
  //  TAttText* pdfText =   pdfplot -> getAttText();
      
  pdfLine->SetLineWidth(2);
  pdfLine->SetLineColor(2);
  pdfLine->SetLineStyle(1);
      
  return theFrame;
}

RooPlot* 
RooFitUtils::addPdfPlot( const RooDataSet* theDataSet, const RooRealVar* theVar, 
		 const RooArgList& thePdfList, 
		 const RooArgList& thePdfListCoef,
		 const char* theName, int nBins )
{
  RooAddPdf* thePdf = new RooAddPdf( "thePdf", "the Pdf", 
				     thePdfList, thePdfListCoef );
  RooPlot* thePlot =  addPdfPlot( theDataSet, theVar, thePdf, theName, nBins );
  delete thePdf;
  return thePlot;
}


RooPlot*
RooFitUtils::pdfPlot( const RooDataSet* theDataSet, const RooRealVar* theVar, 
	      RooAbsPdf* thePdf, const char* theName, int nBins )
{
  RooPlot* theFrame = theVar->frame( nBins ) ;
  theFrame->SetName( theName );
  RooPlot* dataplot = theDataSet -> plotOn(theFrame);
  dataplot -> Print();
  TAttMarker* dataMarker = dataplot -> getAttMarker();
  dataMarker -> SetMarkerStyle(1);
  TAttFill* dataFill = dataplot -> getAttFill();
  dataFill -> SetFillStyle(1001);
  dataFill -> SetFillColor(0);

  // line 
  //  RooPlot* pdfparam = thePdf->paramOn(theFrame, 
  //				     theDataSet, "", 2, 
  //				     "NELU", 0.15, 0.50 );

  RooPlot*  pdfplot =   thePdf->plotOn(theFrame) ;
  TAttLine* pdfLine =   pdfplot -> getAttLine();
  //  TAttFill* pdfFill =   pdfplot -> getAttFill();
  //  TAttText* pdfText =   pdfplot -> getAttText();
      
  pdfLine->SetLineWidth(2);
  pdfLine->SetLineColor(2);
  pdfLine->SetLineStyle(1);
      
  return theFrame;
}

TCanvas*
RooFitUtils::newAddPdfPlot( const RooAbsPdf& pdf, const RooDataSet& dataset, RooRealVar& var, int N_ )
{
  //
  // Gautier's version of RooPlotting of RooAddPdfs
  // (with correct normalization)
  // 
  TString varname = var.GetName();
  TH1* h_data = new TH1F( "h_data", "h_data",
			  var.getBins(),
			  var.getMin(),var.getMax() );
  int ii(0);
  const RooArgSet* argset =  dataset.get(ii);
  while( argset!=0 )
    {
      float val_ = ((RooAbsReal*)argset->find(varname.Data()) )->getVal();
      h_data->Fill( val_ );
      ii++;
      argset =  dataset.get(ii);      
    }
  Double_t n_ = h_data->GetEntries();

  TH1* h_pdf = new TH1F("h_pdf","h_pdf",
		       	N_*var.getBins(),
			var.getMin(),var.getMax() );


  TCanvas* c = new TCanvas();
  c->cd();

  h_data->Draw("0");

  TIterator* pdfIter(0);
  TIterator* coefIter(0);
  int size_pdf(0);
  int size_coef(0);
  const RooAddPdf* pdf_ = dynamic_cast< const RooAddPdf* >(&pdf);
  if( pdf_!=0 )
    {
      const RooArgList& thePdfList = pdf_->pdfList();
      const RooArgList& theCoefList = pdf_->coefList();
      pdfIter  = thePdfList.createIterator() ;
      coefIter = theCoefList.createIterator() ;
      size_pdf  = thePdfList.getSize();
      size_coef = theCoefList.getSize();

      thePdfList.Print();
      cout << "Number of sub pdfs " << size_pdf << endl; 
      cout << "Number of coefs    " << size_coef << endl;       
    }
  
  vector< TString > Name;
  vector< float > Integral;
  vector< float > Coef;
  vector< TH1*  > Hist;
  if( pdfIter!=0 )
    {
      float sum_f(0);
      RooAbsPdf* subpdf_;
      while( ( subpdf_ = (RooAbsPdf*) pdfIter->Next() ) ) 
	{
	  RooAbsReal* coef_ =  (RooAbsReal*) coefIter->Next();
	  float c_(0);
	  if( coef_!=0 )
	    {
	      c_ = coef_->getVal();
	      sum_f += c_;
	    }
	  else 
	    {
	      c_ = 1-sum_f;
	    }
	  float i_ = subpdf_->createIntegral( RooArgSet(var) )->getVal();
	  TString name_ = subpdf_->GetName();
	  TH1* h_  = (TH1*)h_pdf->Clone( TString("h_")+name_ );
	  
	  Integral.push_back(i_);
	  Name.push_back( name_ );
	  Hist.push_back( h_ );
	  Coef.push_back( c_ );
	}
    }

  Integral.push_back( pdf.createIntegral( RooArgSet(var) )->getVal() );
  Name.push_back( pdf.GetName() );
  Hist.push_back( (TH1*)h_pdf->Clone( TString("h_")+pdf.GetName() ) );
      
  TAxis* xaxis = h_pdf->GetXaxis();
  for( int ii=1; ii<=xaxis->GetNbins(); ii++ )
    {  
      float x_ = xaxis->GetBinCenter( ii );
      float b_ = xaxis->GetBinWidth( ii );
      float norm_ = b_*n_*N_;
      var.setVal( x_ );
      
      float sum_(0);
      if( pdfIter!=0 )
	{
	  pdfIter->Reset();
	  int jj(0);
	  RooAbsPdf* subpdf_;
	  while( ( subpdf_ = (RooAbsPdf*) pdfIter->Next() ) ) 
	    {
	      float i_ = Integral[jj];
	      float c_ = Coef[jj];
	      float pdfval_  = norm_*c_*subpdf_->getVal()/i_;
	      sum_ += pdfval_;
	      
	      TString name_ = Name[jj];
	      //	  cout << "Int[" << name_ << "]=" << i_ << endl;      
	      //	  cout << "Coef[" << name_ << "]=" << Integral[jj] << endl; 
	      
	      TH1* h_ = Hist[jj];
	      h_->SetBinContent( ii, pdfval_ );
	      
	      jj++;
	    }
	}
      else
	{
	  sum_ = norm_*pdf.getVal()/Integral[size_pdf];
	}
      Hist[size_pdf] -> SetBinContent( ii, sum_ );
    }
    
  Int_t iColor = 2;
  for( int jj=0; jj<size_pdf; jj++ )
    {
      Hist[jj] -> SetLineStyle(kDashed);
      Hist[jj] -> SetLineColor(iColor);
      Hist[jj] -> Draw("LSame");
      iColor++;
    }
  Hist[size_pdf] -> SetLineColor(kBlue);
  Hist[size_pdf] -> SetLineWidth(2);
  Hist[size_pdf] -> Draw("LSame");

  h_data->Draw("ESame");

  return c;
}

Double_t 
RooFitUtils::funcIntegral( const RooAbsReal& func, RooRealVar& var, 
		   Double_t min, Double_t max, int verbose )
{
  bool keepRange = ( min==0. && max==0. );

  RooArgSet* nset(0);
  
  double min0 = var.getMin();
  double max0 = var.getMax();
  double val0 = var.getVal();
  
  if( keepRange )
    {
      min = min0;
      max = max0;
    }
  else
    {
      assert( max>min );
      
      if( min>min0 || max<max0 ) var.setVal( 0.5*(min+max) );
      if( min>min0 ) var.setMin( min );
      if( max<max0 ) var.setMax( max );
      //      var.Print();
    }
  
  RooArgSet iset( var ) ;

  Double_t norm;
  norm = RooRealIntegral("i","i", func, iset, nset).getVal();

  if( verbose==1 )
    {
      cout << "integral(" << func.GetName() << "|" << min << "," << max << ") = " << norm << endl;
    }

  if( !keepRange )
    {
      var.setMin( min0 );
      var.setMax( max0 );
      var.setVal( val0 );
    }

  return norm;
}

Double_t 
RooFitUtils::funcFraction( const RooAbsReal& func, RooRealVar& var, 
	      Double_t min, Double_t max, int verbose )
{
  //  RooArgSet* nset(0);
  
  double min0 = var.getMin();
  double max0 = var.getMax();

  assert( max>=min );

  if( verbose==1 )   cout << endl;

  Double_t norm1, norm2;
  norm2 = funcIntegral( func, var, min, max, verbose );
  norm1 = funcIntegral( func, var, 0., 0., verbose );
  
  assert( norm1!=0 );

  double frac = norm2/norm1;

  if( verbose==1 )
    {
      cout << "----> frac(" << func.GetName() << "|" << min << "," << max << "|" << min0 << "," << max0 << ")=" << frac*100 << "\%" << endl;
    }

  return frac;
}

Double_t 
RooFitUtils::funcMaxVal( const RooAbsReal& func, RooRealVar& var, 
		 Double_t min, Double_t max, Int_t N, int verbose )
{
  bool keepRange = ( min==0. && max==0. );

  double min0 = var.getMin();
  double max0 = var.getMax();
  double val0 = var.getVal();
  
  if( keepRange )
    {
      min = min0;
      max = max0;
    }
  else
    {
      if( max<=min )
	{
	  cout << "warning max/min " << max << "/" << min << endl;
	  cout << "warning max0/min0 " << max0 << "/" << min0 << endl;
	  func.Print();
	  var.Print();
	  assert( 0 );
	}
      
      if( min<min0 ) min = min0;
      if( max>max0 ) max = max0;
    }

  double maxVal   = 0.;
  double step     = (max-min)/N;
  double x = min;
  double xAtMax = min;
  for( int i=0; i<=N; i++ )
    {
      x = min + i*step;
      var.setVal( x );
      double val = func.getVal();
      if( val>maxVal )
	{
	  maxVal = val;
	  xAtMax = x;
	}
    }

  if( verbose ==1 )
    {
      cout << "max(" << func.GetName() << "|" << min << "," << max << ") = " << maxVal << " for x=" << xAtMax << endl;
    }

  var.setVal( val0 );

  return maxVal;
}

RooFitResult* 
RooFitUtils::fit( RooAbsPdf& pdf, RooDataSet& d )
{
  //  RooDataSet* d_red = (RooDataSet*) d.reduce( _genVars );
  Int_t mark = time(0);
  
  //   RooNLLVar nll("nll","nll",*pdf(),*d_red);
  //   RooMinuit m(nll);
  //   m.optimizeConst(kTRUE);
  //   m.setProfile(kTRUE);
  //   m.setVerbose(kTRUE);
  //   m.migrad();
  //   m.hesse();
  //   RooFitResult* fit = m.save();
  
  RooFitResult* fit = pdf.fitTo( d, RooFit::FitOptions("mHrt") );
// Minimizer(type,algo)           -- Choose minimization package and algorithm to use. Default is MINUIT/MIGRAD through the RooMinuit
//                                    interface, but others can be specified (through RooMinimizer interface)

//                                           Type         Algorithm
//                                           ------       ---------
//                                           Minuit       migrad, simplex, minimize (=migrad+simplex), migradimproved (=migrad+improve)
//                                           Minuit2      migrad, simplex, minimize, scan
//                                           GSLMultiMin  conjugatefr, conjugatepr, bfgs, bfgs2, steepestdescent
//                                           GSLSimAn     -


//  InitialHesse(Bool_t flag)      -- Flag controls if HESSE before MIGRAD as well, off by default
//  Optimize(Bool_t flag)          -- Activate constant term optimization of test statistic during minimization (on by default)
//  Hesse(Bool_t flag)             -- Flag controls if HESSE is run after MIGRAD, on by default
//  Minos(Bool_t flag)             -- Flag controls if MINOS is run after HESSE, on by default
//  Minos(const RooArgSet& set)    -- Only run MINOS on given subset of arguments
//  Save(Bool_t flag)              -- Flac controls if RooFitResult object is produced and returned, off by default
//  Strategy(Int_t flag)           -- Set Minuit strategy (0 through 2, default is 1)
//  FitOptions(const char* optStr) -- Steer fit with classic options string (for backward compatibility). Use of this option
//                                    excludes use of any of the new style steering options.

//  RooFitResult* fit = pdf.fitTo( d, RooFit::Minimizer("Minuit","migrad") );

  cout << "Dataset fit in " << difftime(time(0),mark) << " seconds" << endl;
  
  string fitName( pdf.GetName() );
  fitName += "_fit";
  fit->SetName( fitName.c_str() );
  fit->Print("v");  

  //  delete d_red;
  return fit;
}

ClassImp( RooFitUtils )

