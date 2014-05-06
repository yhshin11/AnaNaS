/*********************************************************************
 *
 *   Original author:    Gautier Hamel de Monchenault
 *   Creation:           Jul 22, 2009
 *   Last modification:  Jul 22, 2009
 *
 *
 *********************************************************************/

#include "Analysis/fit/Model.hh"

#include "Analysis/utils/Constants.hh"
#include "Analysis/utils/Config.hh"

//ghm -- test
#include "Analysis/utils/HistoManager.hh"

#include <cmath>
#include <cassert>
#include <fstream>
using namespace std;

#include <TFile.h>
#include <TH1.h>

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
// models
#include <RooGaussian.h>
#include <RooExponential.h>
#include <RooPolynomial.h>
#include <RooGExpModel.h>
#include <RooVoigtian.h>
#include <RooCBShape.h>
#include <RooArgusBG.h>
#include <RooBreitWigner.h>
#include <RooFFTConvPdf.h>
#include <RooNumConvPdf.h>
#include <RooNovosibirsk.h>
using namespace RooFit;

#include "Analysis/core/Sample.hh"
#include "Analysis/fit/Chi2.hh"

ClassImp(Model)

Model::Model( const char* analysis, 
	      const char* sample, 
	      const char* subsample,
	      const char* hist, 
	      const char* variable 
	      ) 
  :   _analysis(analysis), 
      _sample(sample), 
      _subsample(subsample), 
      _hist(hist), 
      _variable(variable), 
      _type(0)
{  
  // initializations
  init();

  // create dependent variables
  createVariables();

  // building basic density functions
  buildPdfs();  

  // configure
  configure();
}

void
Model::init()
{
  RooRealVar::printScientific(1) ;
  _key        =  _analysis + "__" + _sample + "__" + _subsample + "__" + _hist + "__" + _variable; 
  _configFile = Config::confPath + _key + ".config";
}

void
Model::createVariables()
{
  // configuration variables
  _strVar["mode"]  = new RooStringVar( "mode", "Z decay mode", "Z -> mu+ mu-" );
  _strVar["type"]  = new RooStringVar( "type", "type of fit", "SignalAndBkg" );
  _var["nTot"]     = new RooRealVar("nTot","total number of events", 0 );
  _var["nInRange"] = new RooRealVar("nInRange","number of events in range", 0 );
  _var["eff"]      = new RooFormulaVar("eff","efficiency", "@1/@0",
				       RooArgList( var("nTot"), var("nInRange") ) );

  // constant variables
  _var["m_Z0"]     = new RooRealVar( "m_Z0",     "Z0 mass",  
				     Constants::Z0_m, "GeV/c^{2}" );
  _var["gamma_Z0"] = new RooRealVar( "gamma_Z0", "Z0 width",  
				     Constants::Z0_Gamma, "GeV/c^{2}" );
  // dependent variables
  _var["mll"] = new RooRealVar("mll","mll", 91, 50, 120, "GeV/c^{2}" );
  _var["met"]  = new RooRealVar("met","met", 0, 0, 1000, "GeV/c^{2}" );
}

void
Model::buildPdfs()
{
  if( _variable==TString("mll") )
    {
      //
      // parameters for background
      //
      _var["lambda"] = new RooRealVar( "lambda",    "exponent of background", 
				       0., -100., 0.);

      //
      // parameters for the double Voigtian
      //
      _var["biasNarrow"]  = new RooRealVar( "biasNarrow", "biasNarrow",   
					    0.0, "GeV/c^{2}" );  
      _var["widthNarrow"] = new RooRealVar( "widthNarrow", "widthNarrow", 
					    0.007, "GeV/c^{2}" );  
      _var["fracNarrow"]  = new RooRealVar( "fracNarrow", "fracNarrow",       
					    0.84 );  
      _var["biasWide"]    = new RooRealVar( "biasWide", "biasWide",       
					    0.0, "GeV/c^{2}" );  
      _var["widthWide"]   = new RooRealVar( "widthWide", "widthWide",     
					    0.015, "GeV/c^{2}" );  
      _var["meanNarrow"]  = new RooFormulaVar( "meanNarrow", "meanNarrow",  "@0+@1",
					       RooArgList(var("m_Z0"),var("biasNarrow")) );
      _var["meanWide"]    = new RooFormulaVar( "meanWide", "meanWide",  "@0+@1",
					       RooArgList(var("m_Z0"),var("biasWide")) );
  
      // yields
      _var["n_Z"]     = new RooRealVar( "n_Z","n_Z",     0., 1000000.);
      _var["n_bkg"]   = new RooRealVar( "n_bkg","n_bkg", 0., 100000.);
      _var["f_bkg"]   = new RooRealVar( "f_bkg","f_bkg", 0, 0., 1. );

      // Parameters for Crystal Ball Lineshape 
      _var["cb_bias"]   = new RooRealVar( "cb_bias", "bias", 
					  0.0, -1.0, 1.0,"GeV/c^{2}"); 
      _var["cb_width"]  = new RooRealVar( "cb_width","width", 
					  1.31,0.9,1.8,"GeV/c^{2}"); 
      _var["cb_alpha"]  = new RooRealVar( "cb_alpha","alpha", 1.1,0.6,2.0); 
      _var["cb_power"]  = new RooRealVar( "cb_power","power", 3.5, 0.5, 5.0); 
  
      // Crystal Ball Lineshape
      _pdf["cb_pdf"] = new RooCBShape( "cb_pdf", "A  Crystal Ball Lineshape", 
				       var("mll"), 
				       var("cb_bias"), var("cb_width"), 
				       var("cb_alpha"), var("cb_power") );
  
      // Breit-Wigner Distribution
      _pdf["bw_Z0"] = new RooBreitWigner( "bw_Z0","The true Z0 lineshape (BW)", 
					  var("mll") , 
					  var("m_Z0"), 
					  var("gamma_Z0") );

      _pdf["pdf_sig"] = new RooNumConvPdf("pdf_sig","Convolution",
					  realVar("mll"), 
					  pdf("bw_Z0"), 
					  pdf("cb_pdf") );
      ((RooNumConvPdf*)_pdf["pdf_sig"])
	->setConvolutionWindow( 
			       var("cb_bias"), var("cb_width"), 50 );
  
      _pdf["pdf_bkg"]    = new RooExponential( "pdf_bkg", 
					       "Exponential background", 
					       var("mll"), var("lambda") );
    }
  else if( _variable==TString("met") )
    {
      _var["peak"]   = new RooRealVar("peak","peak",10,0,50);
      _var["tail"]   = new RooRealVar("tail","tail",0.1,-1,5);
      _var["width"]  = new RooRealVar("width","width",10,0,50);
      
      _pdf["pdf_sig"]   = new RooNovosibirsk("pdf_sig", "pdf_sig", 
					     var("met"), 
					     var("peak"), 
					     var("width"), var("tail") );
    }
}

void
Model::configure()
{
  // check the configuration file
  FILE *test_;
  test_ = fopen( _configFile, "r" );
  if( !test_ )
    {
      cout << "Configuration file does not exist -- exit " << _configFile << endl;
      return;
    }
  fclose( test_ );
  cout << "\nconfiguration file --> " << _configFile << endl; 

  // configuration variables
  RooArgSet fitConfig( *_strVar["mode"], 
		       *_strVar["type"],
		       *_var["nTot"],
		       *_var["nInRange"] );

  fitConfig.readFromFile( _configFile.Data(), 0, "Config Data" );
  cout << "mode \t=\t" << _strVar["mode"]->getVal() << endl; 
  cout << "type \t=\t" << _strVar["type"]->getVal() << endl; 
  cout << "eff  \t=\t" << _var["eff"]->getVal()     << endl; 
  cout << "................................done"    << endl ;

  TString type_(_strVar["type"]->getVal());
  if
    ( type_==TString("SignalPlusBkg") ) _type=kSignalPlusBkg;
  else if
    ( type_==TString("SignalOnly") )    _type=kSignalOnly;
  else if
    ( type_==TString("BkgOnly") )       _type=kBkgOnly;

  if( _variable==TString("mll") )
    {
      if( _type==kSignalPlusBkg )
	{
	  _pdf["thePdf"] = new RooAddPdf( "thePdf", "thePdf", 
					  pdf("pdf_bkg"), 
					  pdf("pdf_sig"), 
					  var("f_bkg")    ); 
	}
      else if( _type==kSignalOnly )
	{
	  _pdf["thePdf"] = (RooAbsPdf*) pdf("pdf_sig").clone("thePdf");
	}
      else if( _type==kBkgOnly )
	{
	  _pdf["thePdf"] = (RooAbsPdf*) pdf("pdf_bkg").clone("thePdf");
	}
    }
  else if( _variable==TString("met") )
    {
      _pdf["thePdf"] = (RooAbsPdf*) pdf("pdf_sig").clone("thePdf");
    }

  // configure density function
  cout << "\nConfigure density function" << endl;
  thePdf().Print();
  RooArgSet* params;
  params = thePdf().getParameters( RooDataSet() );
  params->readFromFile( _configFile.Data(), "READ", "Fit Parameters" );
  cout << "Initial Parameters from '" << _configFile << "' : " ; 
  params->selectByAttrib("READ",kTRUE)->Print("v") ;
  cout << "Parameters NOT configured from '" << _configFile << "' : " ;
  params->selectByAttrib("READ",kFALSE)->Print("v") ;
  cout << "\n................................done\n" << endl ;
}

RooFitResult*
Model::fit()
{  
  cout << "\nfit" << endl;

  //  TH1* h_ = _sample.getHist( _analysis, _hist );
  TH1* h_ = Sample::getClone( _analysis, _sample, _subsample, _hist );
  h_->Print();
  float nTot = h_->GetEntries();
  realVar("nTot").setVal( nTot );

  //  float nInRange = h_->GetSumOfWeights();
  float nInRange=0.;
  float xmin_ = theVar().getMin();
  float xmax_ = theVar().getMax();
  cout << "nTot=" << nTot << endl;
  cout << "xmin/xmax " << xmin_ << "/" << xmax_ << endl;
  //  TAxis* xaxis_ = h_->GetXaxis();
  for(int binx=1; binx<=h_->GetNbinsX(); binx++) 
    {
      float low_ = h_->GetBinLowEdge(binx);
      if( low_<xmin_ )  continue;
      if( low_>=xmax_ ) continue;
      nInRange += h_->GetBinContent(binx);
    }
  cout << "nInRange=" << nInRange << endl;
 
  realVar("nInRange").setVal( nInRange ); 

  cout << "nTot/nInRange " << nTot << "/" << nInRange << endl;

  TString fileName_ = Config::histPath + _key + "__chi2.root";
  TFile* f_ = TFile::Open( fileName_, "RECREATE" );

  cout << "Create Chi2 object " << endl;
  //  Chi2 chi2( _key.Data(), h_, theVar().GetName(), thePdf() );
  Chi2 chi2( _key.Data(), h_, theVar(), thePdf() );
  cout << "Done. Now fit." << endl;
  chi2.fit();
  //  chi2.evaluate();
  cout << "Done." << endl;

  bool replace = false;
  TDatime datime;
  TString outfile_;
  outfile_ = Config::confPath;
  outfile_ += _key;
  outfile_ += ".config";
  if( !replace )
    {
      outfile_ += "_";
      outfile_ += "new";
      //      outfile_ += datime.Get();
    }
  ofstream os(outfile_.Data());

  os << "\n\n[Config Data]" << endl;
  os << "//===================" << endl;
  os << "\nmode\t=\t";     _strVar["mode"]  -> writeToStream(os,kFALSE) ;
  os << "\ntype\t=\t";     _strVar["type"]  -> writeToStream(os,kFALSE) ;
  os << "\nnTot\t=\t";     _var["nTot"]     -> writeToStream(os,kFALSE) ;
  os << "\nnInRange\t=\t"; _var["nInRange"] -> writeToStream(os,kFALSE) ;

  os << "\n\n[Fit Parameters]" << endl;
  os << "//===================" << endl;
  os << "// fit performed on " 
     << datime.AsString() << endl;
  chi2.write( os );
  chi2.writeHist(f_);

  f_->Close();
  delete f_;
  cout << "................................done" << endl ;
  return 0;
}

RooAbsReal&
Model::var( const char* varname )
{
  if( _var.count( varname )==0 )
    {
      cout << varname << endl;
      assert(0);
    }
  return *_var[varname];
}

RooRealVar&
Model::realVar( const char* varname )
{
  assert( _var.count( varname )!=0 );
  RooRealVar* realvar_ = dynamic_cast<RooRealVar*>(_var[varname]);
  assert( realvar_!=0 );
  return *realvar_;
}

RooAbsPdf&
Model::pdf( const char* pdfname )
{
  assert( _pdf.count( pdfname )!=0 );
  return *_pdf[pdfname];
}

RooAbsPdf&
Model::thePdf()
{
  return pdf("thePdf");
}

RooRealVar&
Model::theVar()
{
  return realVar(_variable.Data());
}

RooDataSet*
Model::generate( float weight )
{
  // fix me!
  //  int nEvt = ceil(_var["eff"]->getVal()*nTot);
  int nEvt = int( var("nInRange").getVal() * weight );
  RooDataSet* set_ = thePdf().generate( RooArgSet( theVar() ), nEvt );
  return set_;
}

Model::~Model()
{
  cout << "Destroying Model object " << endl;

  // cleanup
  for( map< TString, RooAbsReal* >::iterator it_=_var.begin();
       it_!=_var.end(); ++it_ )
    {
      delete it_->second;
    }
  for( map< TString, RooStringVar* >::iterator it_=_strVar.begin();
       it_!=_strVar.end(); ++it_ )
    {
      delete it_->second;
    }
  for( map< TString, RooAbsPdf* >::iterator it_=_pdf.begin();
       it_!=_pdf.end(); ++it_ )
    {
      delete it_->second;
    }
}

