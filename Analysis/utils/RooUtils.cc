#include "Analysis/utils/RooUtils.hh"

#include <cmath>
#include <iostream>
using namespace std;

#include <TPad.h>

#include "Analysis/utils/Constants.hh"

ClassImp( RooUtils )

void 
RooUtils::grid( TStyle* style, bool gridOn ) 
{
  style->SetPadGridX( gridOn );
  style->SetPadGridY( gridOn );
}

void 
RooUtils::fixOverlay() 
{
  gPad->RedrawAxis();
}

TStyle*
RooUtils::setTdr() 
{

  cout << "+++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "++++    setting up P-TDR style     ++++" << endl;
  cout << "+++++++++++++++++++++++++++++++++++++++" << endl;

  TStyle *tdrStyle = new TStyle( "tdrStyle", "Style for P-TDR" );

  // For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

  // For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

  // For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

  // For the histo:
  // tdrStyle->SetHistFillColor(1);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);

  tdrStyle->SetEndErrorSize(2);
  //  tdrStyle->SetErrorMarker(20); //Ghm -- does not exist???
  tdrStyle->SetErrorX(0.);
  
  tdrStyle->SetMarkerStyle(20);
  
  //For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

  //For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);

  // For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat("nemruo"); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);

  // Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.16);
  //  tdrStyle->SetPadRightMargin(0.02);
  tdrStyle->SetPadRightMargin(0.05);  // GHM

  // For the Global title:

  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);

  // For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.25);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

  // For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  //  tdrStyle->SetLabelFont(41, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "X");
  tdrStyle->SetLabelOffset(0.010, "Y");
  tdrStyle->SetLabelOffset(0.007, "Z");
  tdrStyle->SetLabelSize(0.04, "X");
  tdrStyle->SetLabelSize(0.04, "Y");
  tdrStyle->SetLabelSize(0.04, "Z");

  // For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  //  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetTickLength(0.02, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");

  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickX(1);  
  tdrStyle->SetPadTickY(1);

  // Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

  // Postscript options:
  tdrStyle->SetPaperSize(20.,20.);
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);

  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);

  tdrStyle->SetPalette(1);
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);

  tdrStyle->cd();

  return tdrStyle;
}

void
RooUtils::setTdrHistoStyle( TH1* h )
{
  if( h==0 ) return;
  
  h->SetLineColor(1);
  h->SetLineWidth(1);
  //  h->SetFillColor(38);
  TAxis* axis[3];
  axis[0] = h->GetXaxis();
  axis[1] = h->GetYaxis();
  axis[2] = h->GetZaxis();
  for( int ii=0; ii<3; ii++ )
    {
      TAxis* a = axis[ii];
      if( !a ) continue;
      a->SetLabelColor(1);
      a->SetLabelFont(42);
      a->SetLabelOffset(0.005);
      a->SetLabelSize(0.04);
      a->SetTitleFont(42);
      a->SetTitleOffset(1.2);
      a->SetTitleSize(0.045);
      a->SetNdivisions(510);
      a->SetTickLength(0.02);
    } 

  // special Y axis
  axis[1]->SetLabelOffset(0.010);
  axis[1]->SetTitleOffset(1.65);

  h->SetStats( kFALSE );  
}

void
RooUtils::welcome() 
{
  cout << "+++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "++++    AnaNaS  ... Welcome!       ++++" << endl;
  cout << "++++   Analyse non autorisee       ++++" << endl;
  cout << "++++    de Ntuples a Saclay        ++++" << endl;
  cout << "+++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "++++ Version:       00-01-00       ++++" << endl;
  cout << "++++     samedi 18 septembre 2009  ++++" << endl;
  cout << "+++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "+++ Team members:                   +++" << endl;
  cout << "+++  Federico Ferri                 +++" << endl;
  cout << "+++  Philippe Gras                  +++" << endl;
  cout << "+++  Gautier  Hamel de Monchenault  +++" << endl;
  cout << "+++  Julie    Malcles               +++" << endl;
  cout << "+++  Matthieu Marionneau            +++" << endl;
  cout << "+++++++++++++++++++++++++++++++++++++++" << endl;

}

bool 
RooUtils::toRhoEtaPhi( float x, float y, float z, float& rho, float& eta, float& phi )
{
  rho = sqrt( x*x + y*y );
  if( rho==0 ) return false;
  phi = atan2( y, x ); 
  float theta = Constants::pi/2.;
  if( z!=0 ) 
    {
      theta = atan( rho/z );
      if( theta<0 ) theta += 2*Constants::pi;
    }
  eta = -log( tan( theta/2.) );
  return true;
}

bool
RooUtils::toXYZ( float rho, float eta, float phi, float& x, float& y, float& z )
{
  x = rho*cos(phi);
  y = rho*sin(phi);
  float r, theta;
  if( !toRTheta( rho, eta, r, theta ) ) return false;
  float ttheta = tan(theta);
  if( ttheta==0 ) return false;
  z = rho/ttheta;
  return true;
}  

bool
RooUtils::toRTheta( float rho, float eta, float& r, float& theta )
{
  theta = 2*atan( exp( -eta ) );
  if( theta<0 ) theta+=2*Constants::pi;
  float stheta = sin(theta);
  if( stheta==0 ) return false;
  r = rho/stheta;
  return true;
}

int
RooUtils::color( float e, float E )
{
  int ncol = gStyle->GetNumberOfColors();
  float f = e/E;
  if( f>1 ) f=1;
  int jj = (int) (f*ncol);      
  if( jj<1 ) jj=1;
  if( jj==ncol ) jj--;
  int icol = gStyle->GetColorPalette( jj );
  return icol;
}  

