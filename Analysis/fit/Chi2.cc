#include "Analysis/fit/Chi2.hh"

#include <iomanip>
#include <iostream>
#include <fstream>

#include <TFile.h>

#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooMinuit.h>
#include <RooFitResult.h>

using namespace std;

ClassImp(Chi2);

using namespace std;

Chi2::Chi2( const char* name, 
	    const TH1* dataHist, RooRealVar& var, RooAbsPdf& pdf )
  :  RooAbsReal( name, name ), _dataHist(dataHist), _var(var), _pdf(pdf),
     _chi2(0.), _ndof(0), _called(0), 
     _theFitResult(0)
{
  cout << "Entering Chi2 constructor" << endl;

  _fitHist  = (TH1*) dataHist -> Clone("fitHist");
  _fitHist->Reset();
  _chi2Hist = (TH1*) dataHist -> Clone("chi2Hist");
  _chi2Hist->Reset();

  _params = _pdf.getParameters( RooArgSet(_var) );
  //  _var = (RooRealVar*) _params->find(varname);
  //  assert( _var!=0 );

  //  _params = _pdf.getParameters( RooArgSet() );
  //  _params = _pdf.getDependents( RooArgSet() );
  addServerList( *_params );

  _var.Print();
  _pdf.Print();
}

Chi2::~Chi2()
{
  //  delete _fitHist;
  //  delete _chi2Hist;
}

RooFitResult*
Chi2::fit()
{
  cout << "Start the MINUIT fit " << endl;

  // Start Minuit session on Chi2
  RooMinuit m2(*this) ;
  m2.migrad() ;
  m2.hesse() ;

  _theFitResult = m2.save() ;

  cout << "result of chi2 fit" << endl ;
  _theFitResult->Print("v") ;

  cout << "Ciao Bello (if Julie then Ciao Bella)" << endl;

  return _theFitResult;
}

Double_t 
Chi2::evaluate() const
{
  _called++;

  // now fill the signal data set 
  TAxis* xaxis = _dataHist->GetXaxis();
  int _nBin = (int) xaxis->GetNbins();

  for( int i=1; i<=_nBin; i++ )
    {
      double x_ = _dataHist->GetBinCenter(i);
      _var.setVal( x_ );
      double s_(0);
      s_ = _pdf.getVal(); 

      // reset bin	
      _chi2Hist -> SetBinContent( i, 0 );	 
      _fitHist  -> SetBinContent( i, s_ );

      //      cout << "i/x/pdf " << i << "/" << _var.getVal() << "/" << s_ << endl;
    }
  
  // now normalize the signal data set and compute the chi squared
  double n_data = getSumOfWeights(*_dataHist);
  double n_fit  = getSumOfWeights(*_fitHist);
  _fitHist->Scale( n_data/n_fit );

  _ndof=0;
  _chi2=0.;
  for( int i=0; i<_nBin; i++ )
    {
      double n_ = _dataHist -> GetBinContent(i);
      double s_ = _fitHist  -> GetBinContent(i);//  * n_data/n_fit;

      if( n_==0 ) continue;

      double chi2bin = pow( n_-s_, 2 )/n_;
      _chi2 += chi2bin;
      _ndof++;

      // fill control data sets
      _chi2Hist -> SetBinContent( i, chi2bin );
    }

  if( _called%1 == 1 )
    {
      cout << "." << flush;
    }
  if( _called%10 == 1 )
    {
      cout << endl;
      cout << _called << " ";
      cout << "chi2=" << _chi2 << " ";
      cout << "ndof=" << _ndof << " ";
      cout << "chi2/ndof=" << _chi2/_ndof << " ";
      cout << endl;
    }

  return _chi2;

}

void 
Chi2::writeHist( TFile* file ) const
{
  _dataHist->Clone()->Write();
  _fitHist->Clone()->Write();
  _chi2Hist->Clone()->Write();
  _theFitResult->Clone()->Write();

  float n_ = getSumOfWeights(*_fitHist);

  int scale_ = 10;
  TAxis* xaxis = _dataHist->GetXaxis();
  int nbin_ = int(xaxis->GetNbins())*scale_;
  double min_ = xaxis->GetXmin();
  double max_ = xaxis->GetXmax();
  double bin_ = (max_-min_)/nbin_;
  TString vname("pdf__");
  vname += GetName();
  
  TH1* h_model_ = new TH1F( vname, vname, nbin_, min_, max_ );

  for( int ii=0; ii<nbin_; ii++ )
    {
      double x_=min_+(ii+0.5)*bin_;
      _var.setVal( x_ );
      h_model_ -> Fill( x_, _pdf.getVal() ); 
    }
  float norm_ = scale_*n_/getSumOfWeights(*h_model_);
  h_model_->Scale( norm_ );

  h_model_ -> Write();
}

void
Chi2::write( ofstream& os ) const
{
  os << "// After " << _called << " calls ";
  os << " chi2/ndof=" << _chi2/_ndof << " " << endl;
  
  os << "// Observable" << endl;
  os << _var.GetName() << "\t=\t" ;
  _var.writeToStream(os,kFALSE) ; 
  os << endl;

  TIterator *iterator=0;
  RooAbsArg *next=0;
  
  for( int ii_=0; ii_<2; ii_++ )
    //  for( int ii_=0; ii_<1; ii_++ )  // write only the floating parameters
    {
      bool keepConstant = (ii_==1);
      if( keepConstant )
	{
	  os << "// Constant parameters" << endl;
	}
      else
	{
	  os << "// Floating parameters" << endl;
	}
      iterator = _params->createIterator();
      next=0;
      while((0 != (next= (RooAbsArg*)iterator->Next()))) 
	{
	  bool isConstant = next->isConstant();
	  if(  (isConstant && !keepConstant) 
	       || 
	       (!isConstant && keepConstant) ) continue;
	  os << next->GetName() << "\t=\t" ;
	  next->writeToStream(os,kFALSE) ;
	  os << endl ;
	}
      delete iterator;
    }  
}

double
Chi2::getSumOfWeights( const TH1& h ) const
{
  double n_=0.;
  float xmin_ = _var.getMin();
  float xmax_ = _var.getMax();
  for(int binx=1; binx<=h.GetNbinsX(); binx++) 
    {
      float low_ = h.GetBinLowEdge(binx);
      if( low_<xmin_ )  continue;
      if( low_>=xmax_ ) continue;
      n_ += h.GetBinContent(binx);
    }
  return n_;
}
