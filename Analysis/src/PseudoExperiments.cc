#include "Analysis/src/PseudoExperiments.hh"

#include <cassert>

#include "Analysis/utils/Config.hh"
#include "Analysis/utils/RooFitUtils.hh"

using namespace std;

#include <RooRandom.h>
#include <RooDataSet.h>
#include <RooFit.h>
#include <RooCmdArg.h>
#include <RooFitResult.h>
#include <RooAbsReal.h>
#include <RooAbsPdf.h>
#include <RooPlot.h>
#include <RooGlobalFunc.h>
using namespace RooFit;

ClassImp( PseudoExperiments )

PseudoExperiments::PseudoExperiments( Sample& sample, float lumi ) : SampleAnalysis( "PseudoExp", sample), _model( "Diboson", "Z_2e", "all", "Zee_mass", "mll" ), _rooArgSet(), _dataset(0), _lumi(lumi), _n(0), _first(0), _last(0), _random(false)
{
  cout << "\t-------------------------------------------" << endl; 
  cout << "\t---     Pseudo experiments             ----" << endl; 
  cout << "\t-------------------------------------------" << endl; 

  _rooArgSet.add(_model.var("mll"));
  _rooArgSet.add(_model.var("met"));
  finalisePseudoExperiment();
}

PseudoExperiments::~PseudoExperiments()
{
}

bool
PseudoExperiments::analyzeEvent()
{
  if( _ievt>_last )
    {
      cout << "ievt/last " << _ievt << "/" << _last << endl; 
      finalisePseudoExperiment();
    }

  int imode  = EventServer::kZee;
  std::string Zmode  = EventServer::Zmode[imode];

  _e.loadZllCandidates( imode );

  // choose the best Z candidate in the event
  HCandidate* ZCand = _e.bestZCandidate( imode );

  // request at least one Z candidate in the event -- take the "best" one
  if( ZCand==0 ) return false;

  // get dilepton mass of candidate
  TLorentzVector Z_xyzt = ZCand->P4();
  float mass_     = Z_xyzt.Mag();
  if( mass_<_model.realVar("mll").getMin() 
      || 
      mass_>_model.realVar("mll").getMax() )
    {
      return false;
    }

  // angle of Z candidate to missing ET
  const TVector2& MET_xy = _e.caloMET();
  float met_ = MET_xy.Mod();

  _model.realVar("mll").setVal( mass_ );
  _model.realVar("met").setVal( met_ );
  _dataset->add( _rooArgSet );

  return true;
}

void 
PseudoExperiments::finalisePseudoExperiment()
{
  bool doFit = false;
  if( _dataset!=0 ) 
    {
      TString fname_ = Config::workPath + "pseudoExps/" + "exp_";
      fname_ += _sample.name();
      fname_ += "_";
      fname_ += int(_lumi);
      fname_ += "_";
      fname_ += _n;
      fname_ += ".root";
      TFile* f_ = TFile::Open( fname_, "RECREATE" );
      _dataset->Print();
      _dataset -> Write();
      if( doFit )
	{
	  RooRealVar& var_ = _model.theVar();
	  RooAbsPdf& pdf_  = _model.thePdf();
	  RooFitResult* result_ = RooFitUtils::fit( pdf_, *_dataset );
	  assert( result_!=0 ); result_->Print();
	  result_ -> Write();	  
	  // Plot the fit results
	  RooPlot* plot_ = var_.frame(ObjectName("frame"),Range(50,120),Bins(70));
	  plot_->SetName("rooPlot");
	  _dataset->plotOn(plot_);
	  pdf_.plotOn(plot_, 
		      Components("pdf_bkg"), LineStyle(kDashed), LineColor(kRed));
	  pdf_.plotOn(plot_);
	  plot_->Write();
	}
      f_->Close();
    }
  delete _dataset;

  _n++;
  TString setname_("dataset");
  //  setname_ += _n;
  _dataset = new RooDataSet( setname_, setname_, _rooArgSet );
  determineLastEvent();
}

void 
PseudoExperiments::determineLastEvent()
{
  float n_ =  _sample.n( _lumi );
  if( _random )
    {
      float g_ = RooRandom::gaussian();
      float sig_ = sqrt(n_);
      n_ += g_*sig_;
    }
  _first = _ievt;
  _last  = _first + int(n_) - 1;
  cout << "Experiment[" << _n << "]=(" << _first << "-" << _last << ")"
       << " --> " << _last-_first +1 << "/" << n_ << endl;
}

void
PseudoExperiments::bookHistograms()
{
}

void
PseudoExperiments::writeHistograms()
{
}


