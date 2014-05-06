#ifndef Chi2_hh
#define Chi2_hh

#include <fstream>

#include <RooRealVar.h>
#include <RooAbsPdf.h>
#include <TH1.h>

class RooFitResult;
class TFile;

class Chi2 : public RooAbsReal
{
public:

  Chi2(  const char* name, 
  	 const TH1* hist, RooRealVar& var, RooAbsPdf& pdf );

  virtual TObject* clone( const char* newname ) const 
  { 
    return 0; 
  }

  virtual ~Chi2();

  RooFitResult* fit();

  void writeHist( TFile* file ) const;
  void write( ofstream& ) const;
 
public:

  virtual Double_t evaluate() const;

protected:

  const TH1* _dataHist;
  RooRealVar&     _var;
  RooAbsPdf&      _pdf;
  
  mutable TH1* _fitHist;
  mutable TH1* _chi2Hist;

  mutable double _chi2;
  mutable int    _ndof;
  mutable int _called;

  RooFitResult* _theFitResult;

  RooArgSet* _params;

  double getSumOfWeights( const TH1& h ) const;

  ClassDef(Chi2,0)
};
#endif
