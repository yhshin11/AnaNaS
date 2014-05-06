#include "Analysis/core/PseudoExp.hh"

#include <cassert>

#include "Analysis/utils/Config.hh"
#include "Analysis/utils/RooFitUtils.hh"

using namespace std;

#include <TFile.h>

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

ClassImp( PseudoExp )

PseudoExp::PseudoExp( Sample& sample, float lumi ) : _theSample(sample), _rooArgSet(), _dataset(0), _lumi(lumi), _n(0), _current(0), _first(0), _last(0), _random(false)
{
  cout << "\t-------------------------------------------" << endl; 
  cout << "\t---     Allows Pseudo Experiments      ----" << endl; 
  cout << "\t-------------------------------------------" << endl; 
}

PseudoExp::~PseudoExp()
{
}

 
void
PseudoExp::writeDatasets()
{
  cout << "Entering PseudoExp writeDatasets" << endl;
  for( size_t ii=0; ii<_datasetStore.size(); ii++ )
    {
      RooDataSet* set_ = _datasetStore[ii];
      TString fname_ = Config::workPath + "pseudoExps/" + "exp_";
      fname_ += _theSample.name();
      fname_ += "_";
      fname_ += int(_lumi);
      fname_ += "_";
      fname_ += ii;
      fname_ += ".root";
      TFile* f_ = TFile::Open( fname_, "RECREATE" );
      set_->Print();      
      set_ -> Write();
      f_->Close();
    }
}
  
void 
PseudoExp::createVar(  const char* name, float min, float max )
{
  assert( _var.count( name )==0 );
  RooRealVar* var_ = new RooRealVar(name,name, min, max );
  _var[name] = var_;
  _rooArgSet.add( *var_ );
}

void 
PseudoExp::setVarValue( const char* name, float val )
{
  assert( _var.count( name )==1 );
  _var[name]->setVal( val );
}

RooDataSet&
PseudoExp::dataset()
{
  if( _dataset==0 )
    {
      _n++;
      TString setname_("dataset");
      //  setname_ += _n;
      _dataset = new RooDataSet( setname_, setname_, _rooArgSet );
      determineLastEvent();
    }
  return *_dataset;
}

void
PseudoExp::checkEvent( int ievt )
{
  _current = ievt;
 
  if( _current>_last )
    {
      cout << "current/last " << _current << "/" << _last << endl; 
      finalisePseudoExperiment();
    }
}

void
PseudoExp::saveEvent()
{
  dataset().add( _rooArgSet );
}

void 
PseudoExp::finalisePseudoExperiment()
{
  if( _dataset!=0 ) 
    {

      //      TFile* f_ = TFile::Open( fname_, "RECREATE" );
      _dataset->Print();

      _datasetStore.push_back( _dataset );
      //      _dataset -> Write();
      //      f_->Close();
    }
  //  delete _dataset;
  _dataset = 0;
}

void 
PseudoExp::determineLastEvent()
{
  float n_ =  _theSample.n( _lumi );
  if( _random )
    {
      float g_ = RooRandom::gaussian();
      float sig_ = sqrt(n_);
      n_ += g_*sig_;
    }
  _first = _current;
  _last  = _first + int(n_) - 1;
  cout << "Experiment[" << _n << "]=(" << _first << "-" << _last << ")"
       << " --> " << _last-_first +1 << "/" << n_ << endl;
}



