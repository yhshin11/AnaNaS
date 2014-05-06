#ifndef TupleAnalysis_hh
#define TupleAnalysis_hh

#include <string>
#include <vector>
#include <map>

#include "Analysis/core/Sample.hh"
#include "Analysis/ghmzz/Weights.hh"
#include "Analysis/utils/HistoManager.hh"
#include "Analysis/utils/DatasetManager.hh"
#include "Analysis/utils/TupleManager.hh"

using namespace std;

typedef vector<double> vectorFloat; 

class TupleAnalysis
{
public:

  TupleAnalysis( int im=0, string ip="2011" ) 
    : imode(im), period(ip)
  {}

  virtual ~TupleAnalysis() {}
  
  void analyze( float lumi )
  {
    setOutputFile();
    defineSamples();
    init_();
    analyze_( lumi );
  }
			 
protected:

  virtual void init_();

  virtual void analyze_( float lumi ) =0;

  void defineTemplate( string templ, 
		       int nbin, float min, float max );

  void defineTemplate( string templ, 
		       int nbin, float* xbins );

  TH1* hist( string templ, string cut, string dir, string prefix );

  void fill( string templ, string cut, string dir, string prefix,
	     float val, float weight=1, bool overflow=false, bool norm=true );

  int imode;    // lepton mode: 0=electron, 1=muon
  string period;  // 2010, 2011A, 2011B, etc ...


  string analysis;
  string histfile;
  void setOutputFile();

  // input
  TupleManager   tm;

  // output
  HistoManager   hm;
  DatasetManager dm;
  TupleManager   tmsel;

  map< string, float >   xsect;
  map< string, float >   BR;
  map< string, float>    kf;
  map< string, Sample* > samples;

  void defineSamples();

  map< string, KFactor* > kFactors;
  void setVVKFactors();
  bool VVDynamicKFactors;

  map< string, MomentumWeight* > momWeights;
  void setMomentumWeights();

  // results of Novosibisrsk fits
  map< string, float> tail;
  map< string, float> alpha_0;
  map< string, float> alpha_PU;
  map< string, float> sigma_0;
  map< string, float> sigma_PU;
  
  ClassDef( TupleAnalysis, 0 )  
};

#endif
