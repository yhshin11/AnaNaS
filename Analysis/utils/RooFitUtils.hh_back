#ifndef RooFitUtils_h
#define RooFitUtils_h

#include <TObject.h>
class TCanvas;
class RooPlot;
class RooAbsReal;
class RooRealVar;
class RooDataSet;
class RooAddPdf;
class RooAbsPdf;
class RooArgList;
class RooFitResult;

class RooFitUtils
{
public:

  virtual ~RooFitUtils(){}

  static void testRooFit();

  static 
  RooPlot* addPdfPlot( const RooDataSet* theDataSet, const RooRealVar* theVar, 
		       RooAddPdf* thePdf, const char* theName, int nBins );
  
  static 
  RooPlot* addPdfPlot( const RooDataSet* theDataSet, const RooRealVar* theVar, 
		       const RooArgList& thePdfList, 
		       const RooArgList& thePdfListCoef,
		       const char* theName, int nBins );
  
  static 
  RooPlot* pdfPlot( const RooDataSet* theDataSet, const RooRealVar* theVar, 
		    RooAbsPdf* thePdf, const char* theName, int nBins );
  
  static
  TCanvas*
  newAddPdfPlot( const RooAbsPdf& pdf, const RooDataSet& dataset, RooRealVar& var, int N_=10 );

  static
  RooFitResult* fit( RooAbsPdf& pdf, RooDataSet& d );
  
  static
  Double_t funcIntegral( const RooAbsReal& func, RooRealVar& var, 
			 Double_t min=0., Double_t max=0., int verbose=0 );
  
  static
  Double_t funcFraction( const RooAbsReal& func, RooRealVar& var, 
			 Double_t min, Double_t max, int verbose=0 );
  
  static
  Double_t funcMaxVal( const RooAbsReal& func, RooRealVar& var, 
		       Double_t min=0., Double_t max=0., Int_t N=10, 
		       int verbose=0 );

  ClassDef( RooFitUtils, 0 )
};

#endif
