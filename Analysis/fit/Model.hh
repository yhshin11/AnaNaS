/*********************************************************************
 * 
 *********************************************************************/

#ifndef Model_hh
#define Model_hh

#include <map>
#include <string>

#include <TString.h>
#include <TNtuple.h>

#include <RooArgSet.h>
#include <RooDataSet.h>
#include <RooStringVar.h>
#include <RooRealVar.h>
#include <RooFormulaVar.h>
#include <RooFitResult.h>
#include <RooAbsPdf.h>

class Model
{
public:

  enum { kSignalPlusBkg=0, kSignalOnly, kBkgOnly };

  Model( const char* analysis  = "Diboson",
	 const char* sample    = "Z_2e",
	 const char* subsample = "all",
	 const char* hist      = "stat__mll__0_ZZ_2e2n_presel",
	 const char* variable  = "mll" );

  virtual ~Model();
  
  RooFitResult* fit();

  RooDataSet* generate( float weight );

  RooAbsReal&     var( const char* );
  RooRealVar& realVar( const char* );
  RooAbsPdf&      pdf( const char* );

  RooRealVar& theVar();
  RooAbsPdf&  thePdf();

protected:

  void init();
  void createVariables();
  void buildPdfs();
  void configure();

  //  TString _sample;
  TString _analysis;
  TString _sample;
  TString _subsample;
  TString _hist;
  TString _variable;
  int _type;

  // configuration file
  TString _configFile;

  // key for configuration
  TString _key;

  // maps of (owned) variables and pdfs
  std::map< TString, RooAbsReal* >    _var;
  std::map< TString, RooStringVar* >  _strVar;
  std::map< TString, RooAbsPdf* >     _pdf;

  ClassDef(Model,0)
};


#endif
