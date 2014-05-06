using namespace std;

#include "Display/src/MultiDisplay.hh"

#include <cassert>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <string>
#include <list>
#include <cmath>
#include <algorithm>

#include <TStyle.h>
#include <TArc.h>
#include <TArrow.h>
#include <TEllipse.h>
#include <THelix.h>
#include <TMarker.h>
#include <TText.h>
#include <TRootCanvas.h>

#include "MECore/src/ME.hh"

#include "Analysis/utils/Constants.hh"
#include "Analysis/utils/RooUtils.hh"
#include "Analysis/utils/KineUtils.hh"

TString MultiDisplay::canvName[ MultiDisplay::nCanv ] = 
  { "XY", "RZ", "XZ", "YZ", "VXY", "VRZ", "VXZ", "VYZ", "EtaPhi", "Pt", "Hist" };

int   MultiDisplay::showDet   = 1; 

// detector dimensions
float MultiDisplay::rPipe     =   3.0; 
float MultiDisplay::rPixel[3] =  { 4.4, 7.3, 10.2 };
float MultiDisplay::rTibIn    =  20.0;
float MultiDisplay::rPixOut   =  MultiDisplay::rTibIn;
float MultiDisplay::zPixOut   =  55.0; 
float MultiDisplay::rTibOut   =  55.0;
float MultiDisplay::rTobIn    = MultiDisplay::rTibOut; 
float MultiDisplay::rTobOut   = 116.0;
float MultiDisplay::zTobOut   = 118.0;
float MultiDisplay::zTidIn    =  75.0; 
float MultiDisplay::zTibOut   = MultiDisplay::zTidIn; 
float MultiDisplay::zTidOut   = MultiDisplay::zTobOut;
float MultiDisplay::zTecIn    = 118.0;
float MultiDisplay::zTecOut   = 282.0;
float MultiDisplay::rTecOut   = MultiDisplay::rTobOut;
float MultiDisplay::rCtOut    = MultiDisplay::rTobOut;
float MultiDisplay::zCtOut    = MultiDisplay::zTecOut;
float MultiDisplay::rEcalIn   = 129.0;
float MultiDisplay::rEcalOut  = 161.0;
float MultiDisplay::zEcalOut  = 400.0;
float MultiDisplay::rHcalIn   = 166.6;
float MultiDisplay::rHcalOut  = 254.4;
float MultiDisplay::rHcalExt  = 261.6;
float MultiDisplay::rCoilIn   = 286.7;
float MultiDisplay::rCoilOut  = 315.3;
float MultiDisplay::zCoilOut  = 800.0;  // fixme
float MultiDisplay::rYoke1In  =  458.5+1;
float MultiDisplay::rYoke1Out = 488.5-1;
float MultiDisplay::rYoke2In  =  534.5+1;
float MultiDisplay::rYoke2Out = 597.5-1;
float MultiDisplay::rYoke3In  =  647.0+1;
float MultiDisplay::rYoke3Out = 710.0-1;

float MultiDisplay::refRadOut = MultiDisplay::rHcalOut+10;
float MultiDisplay::refRadMid = MultiDisplay::rHcalIn;
float MultiDisplay::refRadIn  = MultiDisplay::rTobOut+1;
float MultiDisplay::caloRad   = MultiDisplay::rEcalIn-1;

ClassImp(MultiDisplay)

MultiDisplay::MultiDisplay() 
: _style(2), 
  fCanvas(0), fPad(0), 
  fTextTL(0), fTextTR(0), fTextBL(0), fTextBR(0)
{
  setSessionStyle();

  setXYZLimits();
  setEtaPhiLimits();
  setVertexLimits();
  setEnergyLimits();

  int ww = 800;
  int wh = 800;

  // style with white background
  _style[0].kCanvasFillStyle   = 1001;
  _style[0].kCanvasBorderSize  = 4;
  _style[0].kCanvasFillColor   = kOrange;
  _style[0].kPadFillStyle      = 1001;
  _style[0].kPadBorderSize     = 0;
  _style[0].kPadFrameFillColor = kWhite;
  _style[0].kPadFrameLineColor = kWhite;
  _style[0].kPadFillColor      = kWhite;
  _style[0].kBeamAxisColor     = kBlack;
  _style[0].kAxisColor         = kWhite;
  _style[0].kLineWidth         = 1;
  _style[0].kPipeFillColor     = kWhite;
  _style[0].kPipeLineColor     = kAzure+1;
  _style[0].kPixelFillColor    = kCyan-10;
  _style[0].kPixelLineColor    = kAzure+1;
  _style[0].kCtFillColor       = kCyan-9;
  _style[0].kCt2FillColor      = kCyan-4;
  _style[0].kCtLineColor       = kAzure+1;
  _style[0].kEcalFillColor     = kSpring+8;
  _style[0].kEcalLineColor     = kGreen;
  //  _style[0].kXtalLineColor     = kMagenta;
  //  _style[0].kXtalFillColor     = kMagenta-10;
  _style[0].kXtalLineColor     = kGreen;
  _style[0].kXtalFillColor     = kSpring+7;
  _style[0].kHcalFillColor     = kYellow-9;
  _style[0].kHcalLineColor     = kYellow+3;
  _style[0].kCoilFillColor     = kAzure-9;
  _style[0].kCoilLineColor     = kAzure+3;
  _style[0].kYokeFillColor     = kRed;
  _style[0].kYokeLineColor     = kRed;
  _style[0].kPtLineColor       = kGray+2;
  _style[0].kPtTextColor       = kBlack;
  _style[0].kPtFillColor       = kWhite;
  _style[0].kCTHadColor        = kBlue;
  _style[0].kCTEmColor         = kRed;

  // style with black background
  _style[1].kCanvasFillStyle   = 1001;
  _style[1].kCanvasBorderSize  = 4;
  _style[1].kCanvasFillColor   = kGray+2;
  _style[1].kPadFillStyle      = 1001;
  _style[1].kPadBorderSize     = 0;
  _style[1].kPadFrameFillColor = kBlack;
  _style[1].kPadFrameLineColor = kBlack;
  _style[1].kPadFillColor      = kBlack;
  _style[1].kBeamAxisColor     = kBlack;
  _style[1].kAxisColor         = kBlack;
  _style[1].kLineWidth         = 3;
  _style[1].kPipeFillColor     = kBlack;
  _style[1].kPipeLineColor     = kGray;
  _style[1].kPixelFillColor    = kGray+2;
  //  _style[1].kPixelFillColor    = kBlack;
  //  _style[1].kPixelLineColor    = kWhite;
  //  _style[1].kPipeFillColor     = kGray+3;
  //  _style[1].kPipeLineColor     = kGray+1;
  //  _style[1].kPixelFillColor    = kBlack;
  _style[1].kPixelLineColor    = kGray;
  _style[1].kCtFillColor       = kBlack;
  _style[1].kCt2FillColor      = kGray+3;
  //_style[1].kCt2FillColor      = kBlack;
  _style[1].kCtLineColor       = kGray;
  //_style[1].kCtLineColor       = kWhite;
  _style[1].kEcalFillColor     = kGreen+3;
  _style[1].kEcalLineColor     = kGreen;
  _style[1].kXtalLineColor     = kGreen+3;
  //  _style[1].kXtalFillColor     = kGreen;
  _style[1].kXtalFillColor     = kSpring+7;
  _style[1].kHcalFillColor     = kYellow-9;
  _style[1].kHcalLineColor     = kYellow+3;
  _style[1].kCoilFillColor     = kAzure-9;
  _style[1].kCoilLineColor     = kAzure+3;
  _style[1].kYokeFillColor     = kRed;
  _style[1].kYokeLineColor     = kRed;
  _style[1].kPtLineColor       = kGray+2;
  _style[1].kPtTextColor       = kWhite;
  _style[1].kPtFillColor       = kBlack;
  _style[1].kCTHadColor        = kYellow;
  _style[1].kCTEmColor         = kRed;

  size_t is=1;


  //
  // define the standard panels
  //

  ww = 800;
  wh = 800;
  setPanel( canvName[cXY], 0, 0, ww, wh, _style[is], true );

  ww = 1600;
  wh = 800;
  setPanel( canvName[cRZ], 0, 0, ww, wh, _style[is], true );

  ww = 1600;
  wh = 800;
  setPanel( canvName[cXZ], 0, 0, ww, wh, _style[is], false );

  ww = 1600;
  wh = 800;
  setPanel( canvName[cYZ], 0, 0, ww, wh, _style[is], false );

  ww = 800;
  wh = 800;
  setPanel( canvName[cVXY],  50, 100, ww, wh, _style[is], false );

  ww = 1600;
  wh = 800;
  setPanel( canvName[cVRZ], 50, 100, ww, wh, _style[is], false );

  ww = 1600;
  wh = 800;
  setPanel( canvName[cVXZ], 50, 100, ww, wh, _style[is], false );

  ww = 1600;
  wh = 800;
  setPanel( canvName[cVYZ], 50, 100, ww, wh, _style[is], false );

  ww = 800;
  wh = 800;
  setPanel( canvName[cPt], 800, 0, ww, wh, _style[0], true );

  ww = 1000;
  wh = 800;
  setPanel( canvName[cEtaPhi], 400, 0, ww, wh, _style[0], true );

  
  setCtLimits();

  //  for( int ii=nCanv-1; ii>=0; ii-- )
  //    {
  //      goToPanel( canvName[ii] );
  //    }
  //  displayAllPanels();

  //  setEvent( 180000, 9999999 );

  _globalEcalHist=0;
}

void
MultiDisplay::setPanel( TString str_, 
			unsigned x, unsigned y,
			unsigned w, unsigned h,
			Style style,
			bool keep )
{
  TString title = str_;
  if( _panels.count( str_ )!=0 ) return;
  //  cout << "creating panel " << str_ << endl;
  Panel& panel_ = _panels[str_];
  panel_.canv = 0;
  panel_.pad  = 0;
  panel_.hist = 0;
  panel_.title = title;
  panel_.x = x;
  panel_.y = y;
  panel_.w = w;
  panel_.h = h;
  panel_.style = style;
  panel_.keep = keep;
}

bool
MultiDisplay::refreshPanel( TString str_, bool force )
{
  
  // if the panel does not exist, exit
  if( _panels.count( str_ )==0 ) 
    {      
      return false;
    }

  //  cout << "entering refreshPanel for " << str_ << endl;

  // OK it exists
  Panel& panel_ = _panels[str_];
  
  if( panel_.canv !=0 )
    {
      assert( panel_.pad!=0 );
      if( force || !panel_.keep )
	{	  
	  //	  panel_.pad->Clear();
	  //	}
	  //!      else
	  //!	{
	  //	  cout << "Hello " << str_ << endl;
	  delete panel_.hist;
	  delete panel_.pad;
	  delete panel_.canv;
	  panel_.pad  = 0;
	  panel_.canv = 0;
	  panel_.hist = 0;
	}
    }
  else
    {
      assert( panel_.pad==0 && panel_.hist==0 );
    }
  //  fCanvas   = panel_.canv;
  //  fPad      = panel_.pad;
  //  fHist     = panel_.hist;
  //  if( fCanvas==0 ) return false;
  return true;
}

bool
MultiDisplay::goToPanel( TString str, bool create )
{
  TString str_=str;

  if( _panels.count( str_ )== 0 ) return false;

  Panel& panel_ = _panels[str_];
  
  TCanvas* canv_ = panel_.canv;
  if( canv_==0 )
    {
      if( !create ) return false;
      TPad* pad_  = panel_.pad;
      TH2*  hist_ = panel_.hist;
      assert( pad_==0 && hist_==0 );

      TString title = panel_.title;
      //      cout << "Create Panel " << title << endl;
      int x = panel_.x;
      int y = panel_.y;
      int w = panel_.w;
      int h = panel_.h;

      float _aspectRatio = float(h)/float(w);

      double fEmptyX     = 0.03*_aspectRatio; 
      double fEmptyY     = 0.03;
      double wpad        = 0.5*(1.0 - 4*fEmptyX);
      double fXLowPad1   = fEmptyX;
      double fXUpPad1    = fXLowPad1 + wpad;
      double fXLowPad2   = fXUpPad1 + 2*fEmptyX;
      double fXUpPad2    = fXLowPad2 + wpad;
      double fYLowPad    = fEmptyY;
      double fYUpPad     = 1.0 - fEmptyY;
  
      TString cname = title;
      fCanvas = new TCanvas( cname, cname, x, y, w, h );
      TRootCanvas* imp = (TRootCanvas*)fCanvas->GetCanvasImp();
      imp->DontCallClose();
      fCanvas->SetEditable();
      fCanvas->Range(0,0,1,1);

      Style style = panel_.style;

      fCanvas->SetFillStyle(  style.kCanvasFillStyle  );
      fCanvas->SetBorderSize( style.kCanvasBorderSize );
      fCanvas->SetFillColor(  style.kCanvasFillColor  );

      fCanvas->cd();

      TString pname = title;
      pname += "_pad";
      fPad = new TPad(pname,pname,fXLowPad1,fYLowPad,fXUpPad2,fYUpPad);

      fPad->SetGridx(0);
      fPad->SetGridy(0);
      fPad->SetLogx(0);
      fPad->SetLogy(0);
      fPad->SetCrosshair(0);
      fPad->SetFillStyle(  style.kPadFillStyle  );
      fPad->SetBorderSize( style.kPadBorderSize );
      fPad->SetNumber(1);

      if( str_.Contains("EtaPhi") )
	{
	  fPad->SetTopMargin( 0.08 );
	  fPad->SetRightMargin( 0.08 );
	  fPad->SetBottomMargin( 0.1 );
	  fPad->SetLeftMargin( 0.1 );
	  fPad->SetTickx();
	  fPad->SetTicky();
	}
      else if( str_.Contains("Hist") )
	{
	  fPad->SetTopMargin( 0.08 );
	  fPad->SetRightMargin( 0.1 );
	  fPad->SetBottomMargin( 0.2 );
	  fPad->SetLeftMargin( 0.2 );
	  fPad->SetTickx();
	  fPad->SetTicky();
	}
      else
	{
	  fPad->SetTopMargin( 0 );
	  fPad->SetRightMargin( 0 );
	  fPad->SetBottomMargin( 0 );
	  fPad->SetLeftMargin( 0 );
	}

      fPad->SetFrameFillColor( style.kPadFrameFillColor );
      fPad->SetFrameLineColor( style.kPadFrameLineColor );
      fPad->SetFillColor(      style.kPadFillColor      );        

      fPad->Draw();
      fPad->cd();

      // now the histogram

      if( str_.Contains("Hist") ) fHist=0;
      else
	{
	  title = "hist";
	  title += "_"; title += str_;
	  float xmin(-100), xmax(100), ymin(-100), ymax(100);
	  if( str_=="XY" )
	    {
	      if( _R<=0 ) _R = 100.;
	      xmin = _x0-_R;
	      xmax = _x0+_R;
	      ymin = _y0-_R;
	      ymax = _y0+_R;
	    }
	  else if( str_=="RZ" || str_=="XZ" || str_=="YZ" )
	    {
	      ymin = _r0-_R;
	      ymax = _r0+_R;
	      xmin = _z0-_L;
	      xmax = _z0+_L;
	    }
	  else if( str_=="VXY" )
	    {
	      xmin = _xPV-_RV;
	      xmax = _xPV+_RV;
	      ymin = _yPV-_RV;
	      ymax = _yPV+_RV;
	    }
	  else if( str_=="VRZ")
	    {
	      ymin = _rPV-_RV;
	      ymax = _rPV+_RV;
	      xmin = _zPV-_LV;
	      xmax = _zPV+_LV;
	    }
	  else if( str_=="VXZ" )
	    {
	      ymin = _xPV-_RV;
	      ymax = _xPV+_RV;
	      xmin = _zPV-_LV;
	      xmax = _zPV+_LV;
	    }
	  else if( str_=="VYZ" )
	    {
	      ymin = _yPV-_RV;
	      ymax = _yPV+_RV;
	      xmin = _zPV-_LV;
	      xmax = _zPV+_LV;
	    }
	  else if( str_=="Pt" )
	    {
	      xmin = -_ptmax;
	      xmax = +_ptmax;
	      ymin = -_ptmax;
	      ymax = +_ptmax;
	    }
	  else if( str_.Contains("EtaPhi") )
	    {
	      xmin = _eta0-_DE;
	      xmax = _eta0+_DE;
	      ymin = _phi0-_DPh;
	      ymax = _phi0+_DPh;
	    }
	  fHist = new TH2F(title,title,10,xmin,xmax,10,ymin,ymax);
	  fHist->SetDirectory(0);
	  fHist->SetStats( kFALSE );
	  fHist->SetTitleSize(0);
	  fHist->SetMinimum(0);
	  fHist->SetMaximum(_emax);
	  fHist->Fill(1000,1000,_emax); // fixme !!!
	  fHist->GetXaxis()->SetNdivisions( 2*100 + 6 );
	  fHist->GetYaxis()->SetNdivisions( 2*100 + 6 );
	  
	  if( str_=="EtaPhi" )
	    {
	      // 	  int kAxisColor = panel_.style.kAxisColor;
	      fHist->GetXaxis()->SetAxisColor(  15 );
	      // 	  fHist->GetXaxis()->SetLabelColor( kAxisColor);
	      fHist->GetXaxis()->SetLabelSize( 0 );
	      fHist->GetXaxis()->SetLabelOffset( 99 );
	      fHist->GetXaxis()->SetTickLength( 0.01 );
	      fHist->GetYaxis()->SetAxisColor(  15);
	      // 	  fHist->GetYaxis()->SetLabelColor( kAxisColor);
	      fHist->GetYaxis()->SetLabelSize( 0 );
	      fHist->GetYaxis()->SetLabelOffset( 99 );
	      fHist->GetYaxis()->SetTickLength( 0.01 );
	    }
	  else
	    {
	      int kAxisColor = panel_.style.kAxisColor;
	      fHist->GetXaxis()->SetAxisColor(  kAxisColor);
	      fHist->GetXaxis()->SetLabelColor( kAxisColor);
	      fHist->GetXaxis()->SetLabelSize( 0 );
	      fHist->GetYaxis()->SetAxisColor(  kAxisColor);
	      fHist->GetYaxis()->SetLabelColor( kAxisColor);
	      fHist->GetYaxis()->SetLabelSize( 0 );
	    }
	  fHist->Draw("COL");
	}

      if( str_=="EtaPhi" )
	{
	  float y0 = -0.32;
	  float x0 =  5.7;
	  int align = 22;
	  float posx[5] = { -4 , -2 , 0 , 2 , 4 };
	  TString labx[5] = {"-4","-2","0","2","4"};
	  float posy[5] = { 0 ,   0.5 ,   1 ,  1.5   ,  2   };
	  TString laby[5] = {"0","#pi/2","#pi","3#pi/2","2#pi"};
	  // draw labels along X
	  TLatex t;
	  t.SetTextAngle(0);
	  t.SetTextSize(0.03);
	  //	  t.SetTextFont(4);
	  t.SetTextColor(15);
	  for ( size_t i=0; i<5; i++ ) 
	    {	      
	      t.SetTextAlign(23);
	      t.DrawLatex( posx[i], y0, labx[i] );
	    }
	  for ( size_t i=0; i<5; i++ ) 
	    {
	      t.SetTextAlign(12);
	      t.DrawLatex(  x0, posy[i], laby[i] );
	      t.SetTextAlign(32);
	      t.DrawLatex( -x0, posy[i], laby[i] );
	    }
	}
      panel_.canv = fCanvas; 
      panel_.pad  = fPad;
      panel_.hist = fHist;      
    }

  assert( fCanvas!=0 && fPad!=0 );
  fCanvas = panel_.canv;
  fPad    = panel_.pad;
  fHist   = panel_.hist;
  curPanel = str_;

  fCanvas->cd();
  fPad->cd();

  return true;
}

void
MultiDisplay::refreshAllPanels()
{
  //  cout << "Refreshing all panels " << endl;

  for( map< TString, Panel >::iterator it=_panels.begin();
       it!=_panels.end(); it++ )
    //  for( size_t i=0; i<nCanv; i++ )
    {
      //      refreshPanel( canvName[i] );
      refreshPanel( it->first );
    }
}

void
MultiDisplay::displayAllPanels()
{
  //  Pt();
  //  RZ( "VRZ" );
  //  RZ( "VXZ" );
  //  RZ( "VYZ" );
  //  XY( "VXY" );
  //  RZ( "XZ"  );
  //  RZ( "YZ"  );
  RZ( "RZ"  );
  XY( "XY"  );
  setEvent( curRun, curEvent );
}


void
MultiDisplay::setSessionStyle()
{
  Color_t fStyleHistColor  = kWhite;
  Color_t fStyleTitleColor = kBlack;
  Color_t fStyleTFColor    = kGray+2;
  Color_t fStyleTTextColor = kBlack;

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(1111);
  gStyle->SetStatFont(22);
  gStyle->SetStatColor(18);
  gStyle->SetTitleFillColor(fStyleTFColor);
  gStyle->SetTitleTextColor(fStyleTTextColor);
  gStyle->SetTitleW( 0.76 );
  gStyle->SetHistFillColor(fStyleHistColor);
  gStyle->SetTitleFont(22);
  gStyle->SetTitleColor(173);
  gStyle->SetTitleFillColor(18);
  gStyle->SetTitleColor(fStyleTitleColor);
  gStyle->SetTitleTextColor(46);
  gStyle->SetLabelSize( 0.02, "XYZ" );

  // use plain black on white colors
  gStyle->SetFrameBorderMode(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetTitleBorderSize(0);

  // use bold lines and markers
  //  gStyle->SetMarkerStyle(8);
  gStyle->SetHistLineWidth(2);
  gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  // do not display any of the standard histogram decorations
  //  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(0);

  const int ncol = 60;
  gStyle->SetNumberContours(ncol);
  gStyle->SetPalette(1);

}

void
MultiDisplay::XY( TString str_ )
{
  assert( str_=="XY" || str_=="VXY" );

  refreshDet();

  if( !goToPanel(str_,true) ) return;
  Style& style = curStyle();

  TArc* arc;

  bool notInEvent(false);

  bool insideEcalOnly = true;

  if( !insideEcalOnly )
    {

      if( _R>rYoke3Out )
	{
	  drawPGon( 15, rYoke3In, rYoke3Out, 0, 1, 
		    style.kYokeLineColor, style.kYokeFillColor, notInEvent );
	}
      
      if( _R>rYoke2Out )
	{
	  drawPGon( 15, rYoke2In, rYoke2Out, 0, 1, 
		    style.kYokeLineColor, style.kYokeFillColor, notInEvent );
	}
      
      if( _R>rYoke1Out )
	{
	  drawPGon( 15, rYoke1In, rYoke1Out, 0, 1, style.kYokeLineColor, 
		    style.kYokeFillColor, notInEvent );
	}
      
      if( _R>rCoilOut )
	{
	  arc = new TArc( 0, 0, rCoilOut );
	  arc->SetLineColor(style.kCoilLineColor);  
	  arc->SetLineWidth(style.kLineWidth);  
	  arc->SetFillColor(style.kCoilFillColor);
	  arc->Draw();
	}
      
      if( _R>rCoilIn )
	{
	  arc = new TArc( 0, 0, rCoilIn );
	  arc->SetLineColor(style.kCoilLineColor);  
	  arc->SetLineWidth(style.kLineWidth);  
	  arc->SetFillColor(style.kBeamAxisColor);
	  arc->Draw();
	}
      
      if( _R>rHcalOut )
	{
	  drawPGon( 10, rHcalIn, rHcalOut, 0, style.kLineWidth, 
		    style.kHcalLineColor, style.kHcalFillColor, notInEvent );
	}
    }

  if( _R>rEcalOut )
    {
      drawPGon( 10, rEcalIn, rEcalOut, 0, style.kLineWidth, 
		style.kEcalLineColor, style.kEcalFillColor, notInEvent );

      MEEBDisplay::drawXy( style.kXtalLineColor, style.kXtalFillColor );
      
      //      TArc* circle = new TArc( 0, 0, refRad );
      //      circle->SetLineColor(kBlue);
      //      circle->Draw();
    }

  if( _R>rTobOut )
    {
      arc = new TArc( 0, 0, rTobOut );
      arc->SetLineColor(style.kCtLineColor);  
      arc->SetLineWidth(style.kLineWidth);  
      arc->SetFillColor(style.kCtFillColor);
      arc->Draw();

      arc = new TArc( 0, 0, rTibOut );
      arc->SetLineColor(style.kCt2FillColor);  
      //      arc->SetLineWidth(style.kLineWidth);  
      arc->SetLineWidth(1);  
      arc->SetFillColor(style.kCt2FillColor);
      arc->Draw();
    }

  //  if( _R>rTibIn )
    {
      arc = new TArc( 0, 0, rTibIn );
      arc->SetLineColor(style.kPixelLineColor);  
      //      arc->SetLineWidth(style.kLineWidth);  
      arc->SetLineWidth(1);  
      arc->SetFillColor(style.kPixelFillColor);
      arc->Draw();

      arc = new TArc( 0, 0, rPixel[2] );
      arc->SetLineColor(style.kPixelLineColor);  
      //      arc->SetLineWidth(style.kLineWidth);  
      arc->SetLineWidth(1);  
      arc->SetFillColor(style.kPixelFillColor);
      arc->Draw();

      arc = new TArc( 0, 0, rPixel[1] );
      arc->SetLineColor(style.kPixelLineColor);  
      //      arc->SetLineWidth(style.kLineWidth);  
      arc->SetLineWidth(1);  
      arc->SetFillColor(style.kPixelFillColor);
      arc->Draw();

      arc = new TArc( 0, 0, rPixel[0] );
      arc->SetLineColor(style.kPixelLineColor);  
      //      arc->SetLineWidth(style.kLineWidth);  
      arc->SetLineWidth(1);  
      arc->SetFillColor(style.kPixelFillColor);
      arc->Draw();

    }

  if( _R>rPipe )
    {
      arc = new TArc( 0, 0, rPipe );
      arc->SetLineColor(style.kPipeLineColor);  
      //      arc->SetLineWidth(style.kLineWidth);  
      arc->SetLineWidth(1);  
      arc->SetFillColor(style.kPipeFillColor);
      arc->Draw();
    }

  fCanvas->Modified();
  fCanvas->Update();
}

void
MultiDisplay::RZ( TString str_ )
{
  assert( str_=="RZ" || str_=="VRZ" ||
	  str_=="XZ" || str_=="VXZ" ||
	  str_=="YZ" || str_=="VYZ"    );

  refreshDet();

  if( !goToPanel(str_,true) ) return;
  Style& style = curStyle();

  //  if( rmax>rEcalOut )
    MEEBDisplay::drawRz( style.kXtalLineColor, style.kXtalFillColor );
    //  if( zmax>zEcalOut )  
    MEEEDisplay::drawRz( style.kXtalLineColor, style.kXtalFillColor );
 
    //  float rmin=_r0-_R; 
  float rmax = _r0+_R; 
  //  float zmin=_z0-_L; 
  float zmax=_z0+_L;

  // draw the central tracker
  if( rmax>rCtOut && zmax>zCtOut )
    {
      drawBox( -zCtOut, zCtOut, -rCtOut, rCtOut, style.kLineWidth, style.kCtLineColor, style.kCtFillColor, false );
    }
  // draw the Tec's
  if( zmax>zTecOut )
    {
      drawBox( zTecIn, zTecOut, -rTecOut, rTecOut, 1, style.kCtLineColor, style.kCtFillColor, false );
      drawBox( -zTecOut, -zTecIn, -rTecOut, rTecOut, 1, style.kCtLineColor, style.kCtFillColor, false );
    }
  // draw the pixels
  if( rmax>rPixOut && zmax>zPixOut )
    {
      drawBox( -zPixOut, zPixOut, -rPixOut, rPixOut, 1, style.kPixelLineColor, style.kPixelFillColor, false );
    }
  // draw the TIB
  if( rmax>rTibOut && zmax>zTibOut )
    {
      drawBox( -zTibOut, zTibOut, rTibIn, rTibOut, 1, style.kCtLineColor, style.kCt2FillColor, false );
      drawBox( -zTibOut, zTibOut, -rTibOut, -rTibIn, 1, style.kCtLineColor, style.kCt2FillColor, false );
    }
  // draw the TID 
  if( rmax>rTibOut && zmax>zTidOut )
    {
      drawBox( zTidIn, zTidOut, rTibIn, rTibOut, 1, style.kCtLineColor, style.kCt2FillColor, false  );
      drawBox( -zTidIn, -zTidOut, rTibIn, rTibOut, 1, style.kCtLineColor, style.kCt2FillColor, false  );
      drawBox( -zTidIn, -zTidOut, -rTibIn, -rTibOut, 1, style.kCtLineColor, style.kCt2FillColor, false  );
      drawBox(  zTidIn,  zTidOut, -rTibIn, -rTibOut, 1, style.kCtLineColor, style.kCt2FillColor, false  );
    }
  // draw the TOB
  if( rmax>rTobOut && zmax>zTobOut )
    {
      drawBox( -zTobOut, zTobOut, rTobIn, rTobOut, 1, style.kCtLineColor, style.kCtFillColor, false  );
      drawBox( -zTobOut, zTobOut, -rTobOut, -rTobIn, 1, style.kCtLineColor, style.kCtFillColor, false  );
    }
  //draw the beam pipe
  if( _L>50 )
    {
      drawBox( -_L, _L, -rPipe, rPipe, 1, style.kPipeLineColor, style.kPipeFillColor, false  );  
    }
  else
    {    
      drawBox( -_L, 0, -rPipe, rPipe, 2, style.kPipeLineColor, style.kPipeFillColor, false  );  
      drawBox(  0, _L, -rPipe, rPipe, 2, style.kPipeLineColor, style.kPipeFillColor, false  );  
      ///      drawLine( hRZ, 0, 0, -2*L, 0, 0, 2*L, 2*L, 4, 1, kWhite );
      TLine* line = new TLine( -2*_L, 0, 2*_L, 0 );
      //      registerTObject( line );
      line->SetLineColor( style.kBeamAxisColor );
      line->SetLineWidth( 1 );
      line->SetLineStyle( kDashed );
      line->Draw();    
    }
      //   //draw the coil pipe
//   float Lcoil = zCoilOut;
//   if( Lcoil>L ) Lcoil = L;
//   drawBox( -Lcoil, Lcoil, rCoilIn, rCoilOut, kLineWidth, kCoilLineColor, kCoilFillColor, false  );  
//   drawBox( -Lcoil, Lcoil, -rCoilIn, -rCoilOut, kLineWidth, kCoilLineColor, kCoilFillColor, false  );  
  
}

void
MultiDisplay::Pt()
{
  if( !goToPanel("Pt",true) ) return;
  Style& style = curStyle();

  TArc* arc;

  int iptcur = (int)( 2*_ptmax*10 );
  iptcur/=10;
  int ii(0);

  int idpt = 10;
  if( _ptmax>150 ) idpt = 50;
  while( iptcur>idpt )
    {
      ii++;
      iptcur -= idpt;
      int lineWidth = 2;
      if( ii%2 ) lineWidth = 1;
      arc = new TArc( 0, 0, iptcur );
      arc->SetLineColor(style.kPtLineColor);  
      arc->SetFillColor(style.kPtFillColor);  
      arc->SetLineWidth(lineWidth);  
      arc->Draw();
      //      registerTObject( arc );
    }

  int nrad = 2;
  for( int ii=0; ii<nrad; ii++ )
    {
      float xx_ = sqrt(2.)*_ptmax*cos( ii*Constants::pi/nrad );
      float yy_ = sqrt(2.)*_ptmax*sin( ii*Constants::pi/nrad );
      TLine* line = new TLine( -xx_, -yy_, xx_, yy_ );
      //      registerTObject( line );
      line->SetLineColor( style.kPtLineColor );
      line->SetLineWidth( 1 );
      int lineStyle = kDashed;
      if( ii==0 || ii==3 ) 
	{
	  lineStyle = kSolid;
	}
      line->SetLineStyle( lineStyle );
      //  line->SetLineStyle( kDashed );
      line->Draw();
    }

  float grad_ = idpt;
  int jjmax = (int)(1.3*_ptmax/grad_);
  if( jjmax<4 )
    {
      grad_ /= 2.;
      jjmax *=2;
    }
  int imod = 1;
  if( jjmax>15 ) imod = 3;
  else if( jjmax>10 ) imod = 2;
  else imod = 1;
  
  float size_ = 0.03;
  int align_ = 22;
  int color_ = style.kPtTextColor;
  float scale_ = 1;

  for( int jj=1; jj<2*jjmax; jj++ )
    {
      if( jj%imod!=0 ) continue;
      float pt_ = jj*grad_;
      if( jj!=1 && pt_>1.3*_ptmax ) break;
      ostringstream oss;
      oss << pt_;
      float xx_[4] = { pt_,   0, -pt_,    0 };
      float yy_[4] = {   0, pt_,    0, -pt_ };
      for( int kk=0; kk<4; kk++ )
	{
	  TText* text_ = new TText( xx_[kk], yy_[kk], oss.str().c_str() );
	  //	  registerTObject( text_ );
	  text_->SetTextSize( scale_*size_ );
	  text_->SetTextAlign( align_ );
	  text_->SetTextColor( color_ );
	  text_->Draw();
	}
    }
}

void
MultiDisplay::EtaPhi()
{
  if( !goToPanel("EtaPhi",true) ) return;
  Style& style = curStyle();

  refreshDet();
  
  drawEtaPhiDet();

}

void
MultiDisplay::drawEtaPhiDet()
{
    {
      {
	// fixme !!!
	for( int ieta=-41; ieta<=41; ieta++ )
	  {
	    for( int iphi=1; iphi<=72; iphi++ )
	      {
		drawCT( ieta, iphi, 1, 15, kWhite, false );
	      }
	  }
      }
    }

  RooUtils::fixOverlay();
}

void
MultiDisplay::drawEtaPhiEcal()
{
  {
    //    MEEEDisplay::bkgColor =  kWhite;
    MEEEDisplay::bkgColor =  -1;
    MEEEDisplay::lineColor = 15;
    MEEEDisplay::lineWidth = 1;
    
    //    MEEBDisplay::bkgColor =  kWhite;
    MEEBDisplay::bkgColor =  -1;
    MEEBDisplay::lineColor = 15;
    MEEBDisplay::lineWidth = 1;
  }
  {
    MEEEDisplay::drawEE( _etamin, _etamax, _phimin, _phimax );
    MEEBDisplay::drawEB( _etamin, _etamax, _phimin, _phimax ); 
  }
  RooUtils::fixOverlay();
}

void
MultiDisplay::drawCT( int ieta, int iphi, int lineWidth, int lineColor, int fillColor, bool isInEvent )
{
  if( ieta==0 || abs(ieta)>41 )
    {
      return; 
    }
  int ii= abs(ieta); 

  {
    float e0 = ct_eta_[ ii-1 ]; 
    float e1 = ct_eta_[ ii   ]; 
    
    // fixme
    if( e1>_DE ) return; 
    
    if( ieta<0 )
      {
	e0*=-1;  e1*=-1; 
      }
    float ph0 = (iphi-1)*2./72; 
    float dph = 2./72; 
    if( ii>20 ) 
      {
	if( (iphi-1)%2!=0 ) return; 
	dph *= 2; 
      }
    if( ii>39 ) 
      {
	if( (iphi-1)%4!=0 ) return;      
	dph *= 2; 
      }
    float ph1 = ph0 + dph; 
    
    drawBox( e0, e1, ph0, ph1 , lineWidth, lineColor, fillColor, isInEvent );      
  }
}

void
MultiDisplay::refresh()
{
  for( list<TObject*>::iterator it=_list.begin(); 
       it!=_list.end(); ++it )
    {
      delete (*it);
      (*it) = 0;
    }
  _list.clear();
}

void
MultiDisplay::refreshDet()
{
  //  MEEBDisplay::refresh();
  //  MEEEDisplay::refresh();
}

void
MultiDisplay::registerTObject( TObject* o )
{
  _list.push_back( o );
}

void 
MultiDisplay::drawBox( float xmin, float xmax, float ymin, float ymax, int lineWidth, int lineColor, int fillColor, bool isInEvent )
{
  double xx_[5];
  double yy_[5];
  xx_[0] = xmin; yy_[0] = ymin;
  xx_[1] = xmax; yy_[1] = ymin;
  xx_[2] = xmax; yy_[2] = ymax;
  xx_[3] = xmin; yy_[3] = ymax;
  xx_[4] = xx_[0]; yy_[4] = yy_[0];
  TPolyLine* pline = new TPolyLine( 5, xx_, yy_ );
  if( isInEvent ) registerTObject( pline );
  pline ->SetFillColor( fillColor );
  pline ->SetLineColor( lineColor );
  pline ->SetLineWidth( lineWidth );
  pline ->Draw("f");
  pline ->Draw();
} 

void
MultiDisplay::drawWedge( float dd, float rin, float rout, float phi0, int lineWidth, int lineColor, int fillColor, bool isInEvent )  
{
  double dd_  = dd;
  float scx0_ = rin;
  float scX0_ = scx0_ + (rout-rin);
  float scy0_ = scx0_*tan(dd_*Constants::pi/180);
  float scY0_ = scX0_*tan(dd_*Constants::pi/180.);
  float sphi_ = sin( phi0*Constants::pi/180. );
  float cphi_ = cos( phi0*Constants::pi/180. );
  double scxx[5];
  double scyy[5];
  scxx[0] =  scx0_*cphi_ + scy0_*sphi_;
  scyy[0] =  scx0_*sphi_ - scy0_*cphi_;
  scxx[1] =  scX0_*cphi_ + scY0_*sphi_;
  scyy[1] =  scX0_*sphi_ - scY0_*cphi_;
  scxx[2] =  scX0_*cphi_ - scY0_*sphi_;
  scyy[2] =  scX0_*sphi_ + scY0_*cphi_;
  scxx[3] =  scx0_*cphi_ - scy0_*sphi_;
  scyy[3] =  scx0_*sphi_ + scy0_*cphi_;
  scxx[4] = scxx[0];
  scyy[4] = scyy[0];
  TPolyLine* pline = new TPolyLine( 5, scxx, scyy );
  if( isInEvent ) registerTObject( pline );
  pline ->SetFillColor( fillColor );
  pline ->SetLineColor( lineColor );
  pline ->SetLineWidth( lineWidth );
  pline ->Draw("f");
  pline ->Draw();
}

void
MultiDisplay::drawPGon( float dd, float rin, float rout, float phi0, int lineWidth, int lineColor, int fillColor, bool isInEvent )  
{
  double dd_  = dd;
  float scx0_ = rin;
  float scX0_ = scx0_ + (rout-rin);
  float scy0_ = scx0_*tan(dd_*Constants::pi/180);
  float scY0_ = scX0_*tan(dd_*Constants::pi/180.);
  float phi_ = 0;
  while( phi_<360 )
    {
      float sphi_ = sin( (phi0+phi_)*Constants::pi/180. );
      float cphi_ = cos( (phi0+phi_)*Constants::pi/180. );
      double scxx[5];
      double scyy[5];
      scxx[0] =  scx0_*cphi_ + scy0_*sphi_;
      scyy[0] =  scx0_*sphi_ - scy0_*cphi_;
      scxx[1] =  scX0_*cphi_ + scY0_*sphi_;
      scyy[1] =  scX0_*sphi_ - scY0_*cphi_;
      scxx[2] =  scX0_*cphi_ - scY0_*sphi_;
      scyy[2] =  scX0_*sphi_ + scY0_*cphi_;
      scxx[3] =  scx0_*cphi_ - scy0_*sphi_;
      scyy[3] =  scx0_*sphi_ + scy0_*cphi_;
      scxx[4] = scxx[0];
      scyy[4] = scyy[0];
      TPolyLine* pline = new TPolyLine( 5, scxx, scyy );
      if( isInEvent )     registerTObject( pline );
      pline ->SetFillColor( fillColor );
      pline ->SetLineColor( lineColor );
      pline ->SetLineWidth( lineWidth );
      pline ->Draw("f");
      pline ->Draw();
      phi_ += 2*dd_;
    }  
}


void
MultiDisplay::setCtLimits()
{
  ct_eta_[0]=0.000; ct_etaout_[0]=0.000;
  ct_eta_[1]=0.087; ct_etaout_[1]=0.000;
  ct_eta_[2]=0.174; ct_etaout_[2]=0.087;
  ct_eta_[3]=0.261; ct_etaout_[3]=0.174;
  ct_eta_[4]=0.348; ct_etaout_[4]=0.261;
  ct_eta_[5]=0.435; ct_etaout_[5]=0.348;
  ct_eta_[6]=0.522; ct_etaout_[6]=0.435;
  ct_eta_[7]=0.609; ct_etaout_[7]=0.522;
  ct_eta_[8]=0.696; ct_etaout_[8]=0.609;
  ct_eta_[9]=0.783; ct_etaout_[9]=0.696;
  ct_eta_[10]=0.870; ct_etaout_[10]=0.783; 
  ct_eta_[11]=0.957; ct_etaout_[11]=0.870; 
  ct_eta_[12]=1.044; ct_etaout_[12]=0.957; 
  ct_eta_[13]=1.131; ct_etaout_[13]=1.044; 
  ct_eta_[14]=1.218; ct_etaout_[14]=1.131; 
  ct_eta_[15]=1.305; ct_etaout_[15]=1.218; 
  ct_eta_[16]=1.392; ct_etaout_[16]=1.305; 
  ct_eta_[17]=1.479; ct_etaout_[17]=1.392; 
  ct_eta_[18]=1.566; ct_etaout_[18]=1.479;
  ct_eta_[19]=1.653; ct_etaout_[19]=1.566;
  ct_eta_[20]=1.740; ct_etaout_[20]=1.653;
  ct_eta_[21]=1.830; ct_etaout_[21]=1.740;
  ct_eta_[22]=1.930; ct_etaout_[22]=1.830;
  ct_eta_[23]=2.043; ct_etaout_[23]=1.930;
  ct_eta_[24]=2.172; ct_etaout_[24]=2.043;
  ct_eta_[25]=2.332; ct_etaout_[25]=2.172;
  ct_eta_[26]=2.500; ct_etaout_[26]=2.332;
  ct_eta_[27]=2.650; ct_etaout_[27]=2.500;
  ct_eta_[28]=2.868; ct_etaout_[28]=2.650;
  ct_eta_[29]=3.000; ct_etaout_[29]=2.868;
  ct_eta_[30]=3.139; ct_etaout_[30]=2.964;
  ct_eta_[31]=3.314; ct_etaout_[31]=3.139;
  ct_eta_[32]=3.489; ct_etaout_[32]=3.314;
  ct_eta_[33]=3.664; ct_etaout_[33]=3.489;
  ct_eta_[34]=3.839; ct_etaout_[34]=3.664;
  ct_eta_[35]=4.013; ct_etaout_[35]=3.839;
  ct_eta_[36]=4.191; ct_etaout_[36]=4.013;
  ct_eta_[37]=4.363; ct_etaout_[37]=4.191;
  ct_eta_[38]=4.538; ct_etaout_[38]=4.363;
  ct_eta_[39]=4.716; ct_etaout_[39]=4.538;
  ct_eta_[40]=4.889; ct_etaout_[40]=4.716;
  ct_eta_[41]=5.191; ct_etaout_[41]=4.889;

  ct_rin_=143.136; ct_rout_=407.388; 
  ct_zin_[0]=0.000; ct_zout_[0]=0.000; 
  ct_zin_[1]=12.469; ct_zout_[1]=35.487; 
  ct_zin_[2]=25.032; ct_zout_[2]=71.244; 
  ct_zin_[3]=37.784; ct_zout_[3]=107.540; 
  ct_zin_[4]=50.823; ct_zout_[4]=144.650; 
  ct_zin_[5]=64.247; ct_zout_[5]=182.856; 
  ct_zin_[6]=78.157; ct_zout_[6]=222.446; 
  ct_zin_[7]=92.659; ct_zout_[7]=263.722; 
  ct_zin_[8]=107.863; ct_zout_[8]=306.995; 
  ct_zin_[9]=123.884; ct_zout_[9]=352.593; 
  ct_zin_[10]=140.843; ct_zout_[10]=400.861; 
  ct_zin_[11]=158.869; ct_zout_[11]=452.166; 
  ct_zin_[12]=178.098; ct_zout_[12]=506.895; 
  ct_zin_[13]=198.676; ct_zout_[13]=565.463; 
  ct_zin_[14]=220.759; ct_zout_[14]=628.314; 
  ct_zin_[15]=244.514; ct_zout_[15]=695.923; 
  ct_zin_[16]=270.120; ct_zout_[16]=768.804; 
  ct_zin_[17]=297.773; ct_zout_[17]=847.507; 

  ct_zin1_=320.000;  ct_zout1_=568.000;
  ct_rin1_[0]=153.821; ct_rout1_[0]=273.032; 
  ct_rin1_[1]=139.781; ct_rout1_[1]=248.112; 
  ct_rin1_[2]=127.208; ct_rout1_[2]=225.793; 
  ct_rin1_[3]=115.904; ct_rout1_[3]=205.729; 
  ct_rin1_[4]=105.376; ct_rout1_[4]=187.043; 
  ct_rin1_[5]=94.894;  ct_rout1_[5]=168.437; 
  ct_rin1_[6]=84.387;  ct_rout1_[6]=149.788; 
  ct_rin1_[7]=73.887;  ct_rout1_[7]=131.150; 
  ct_rin1_[8]=62.736;  ct_rout1_[8]=111.357; 
  ct_rin1_[9]=52.891;  ct_rout1_[9]=93.881; 
  ct_rin1_[10]=45.444; ct_rout1_[10]=80.662; 
  ct_rin1_[11]=36.478; ct_rout1_[11]=64.748; 
  ct_rin1_[12]=31.943; ct_rout1_[12]=56.699; 
 
  ct_zin2_=1100.000;  ct_zout2_=1265.000;
  ct_rin2_[0]=113.850; ct_rout2_[0]=130.927; 
  ct_rin2_[1]=95.497; ct_rout2_[1]=109.821; 
  ct_rin2_[2]=80.121; ct_rout2_[2]=92.139; 
  ct_rin2_[3]=67.232; ct_rout2_[3]=77.316; 
  ct_rin2_[4]=56.423; ct_rout2_[4]=64.886; 
  ct_rin2_[5]=47.355; ct_rout2_[5]=54.458; 
  ct_rin2_[6]=39.787; ct_rout2_[6]=45.755; 
  ct_rin2_[7]=33.296; ct_rout2_[7]=38.291; 
  ct_rin2_[8]=28.033; ct_rout2_[8]=32.238; 
  ct_rin2_[9]=23.531; ct_rout2_[9]=27.061; 
  ct_rin2_[10]=19.694; ct_rout2_[10]=22.648; 
  ct_rin2_[11]=16.565; ct_rout2_[11]=19.049; 
  ct_rin2_[12]=12.247; ct_rout2_[12]=14.084; 

  //   ct_etaLimits[0]=0.; 
  //   ct_etaLimits[1]=0.087; 
  //   ct_etaLimits[2]=0.174; 
  //   ct_etaLimits[3]=0.261; 
  //   ct_etaLimits[4]=0.348; 
  //   ct_etaLimits[5]=0.435; 
  //   ct_etaLimits[6]=0.522; 
  //   ct_etaLimits[7]=0.609; 
  //   ct_etaLimits[8]=0.696; 
  //   ct_etaLimits[9]=0.783; 
  //   ct_etaLimits[10]=0.870; 
  //   ct_etaLimits[11]=0.957; 
  //   ct_etaLimits[12]=1.044; 
  //   ct_etaLimits[13]=1.131; 
  //   ct_etaLimits[14]=1.218; 
  //   ct_etaLimits[15]=1.305; 
  //   ct_etaLimits[16]=1.392; 
  //   ct_etaLimits[17]=1.479; 
  //   ct_etaLimits[18]=1.566; 
  //   ct_etaLimits[19]=1.653; 
  //   ct_etaLimits[20]=1.740; 
  //   ct_etaLimits[21]=1.830; 
  //   ct_etaLimits[22]=1.930; 
  //   ct_etaLimits[23]=2.043; 
  //   ct_etaLimits[24]=2.172; 
  //   ct_etaLimits[25]=2.332; 
  //   ct_etaLimits[26]=2.500; 
  //   ct_etaLimits[27]=2.650; 
  //   ct_etaLimits[28]=2.868; 
  //   ct_etaLimits[29]=3.000; 
  //   ct_etaLimits[30]=3.139; 
  //   ct_etaLimits[31]=3.314; 
  //   ct_etaLimits[32]=3.489; 
  //   ct_etaLimits[33]=3.664; 
  //   ct_etaLimits[34]=3.839; 
  //   ct_etaLimits[35]=4.013; 
  //   ct_etaLimits[36]=4.191; 
  //   ct_etaLimits[37]=4.363; 
  //   ct_etaLimits[38]=4.538; 
  //   ct_etaLimits[39]=4.716; 
  //   ct_etaLimits[40]=4.889; 
  //   ct_etaLimits[41]=5.191; 
}

bool
MultiDisplay::isInView( TPolyLine* pl, float xmin, float xmax, float ymin, float ymax )
{
  return true;
}

bool
MultiDisplay::setEvent( int run, int event )
{
  curRun = run;
  curEvent = event;

  _runStr  = "Run";
  _runStr += curRun;
  _runStr += "_";
  _runStr += "Evt";
  _runStr += curEvent;
  _runStr += "_";
  
  _runStrLong  = "Run ";
  _runStrLong += curRun;
  _runStrLong += " ";
  _runStrLong += "Event ";
  _runStrLong += curEvent;
  
  cout << _runStrLong << endl;

  drawComments();

  refresh();

  return true;
}

// void
// MultiDisplay::setComments()
// {
//   TString comment1;
//   comment1 = "Run ";
//   comment1 += curRun;
//   comment1 += " Event ";
//   comment1 += curEvent;
//   setComments( "CMS Experiment at the LHC", 
// 	       "AnaNaS Event Display", 
// 	       comment1, "" );
// }


void
MultiDisplay::drawComments( )
{  
  drawComments( "CMS Experiment at the LHC", 
		"AnaNaS Event Display", 
		_runStrLong.Data(), "" );
}

void
MultiDisplay::drawComments( const char* textTL, const char* textBL, const char* textTR, const char* textBR )
{  
  TString fStringTL(textTL);
  TString fStringBL(textBL);
  TString fStringTR(textTR);
  TString fStringBR(textBR);

  Width_t    fWidthTex     = 2;

  double   fXTexTL         = 0.03;
  double   fYTexTL         = 0.975;
  double   fXTexBL         = 0.03;
  double   fYTexBL         = 0.030;
  double   fXTexBR         = 0.97;
  double   fYTexBR         = 0.030;
  double   fXTexTR         = 0.97;
  double   fYTexTR         = 0.975;

  //  Font_t     fFontTex      = 132;
  Font_t     fFontTex      = 62;
  int        fColorTex     = kWhite;
  Float_t    fSizeTex      = 0.03;

  delete fTextTL;
  delete fTextTR;
  delete fTextBL;
  delete fTextBR;
  fTextTL = new TLatex(fXTexTL,fYTexTL,fStringTL.Data());
  fTextTR = new TLatex(fXTexTR,fYTexTR,fStringTR.Data());
  fTextBL = new TLatex(fXTexBL,fYTexBL,fStringBL.Data());
  fTextBR = new TLatex(fXTexBR,fYTexBR,fStringBR.Data());
  // from DrawLabels()
  //Draws the 4 labels on the Canvas
  fTextTL->SetTextFont(fFontTex);
  fTextTL->SetTextSize(fSizeTex);
  fTextTL->SetTextColor(fColorTex);
  fTextTL->SetLineWidth(fWidthTex);
  fTextTL->SetTextAlign(11);
  fTextBL->SetTextFont(fFontTex);
  fTextBL->SetTextSize(fSizeTex);
  fTextBL->SetTextColor(fColorTex);
  fTextBL->SetLineWidth(fWidthTex);
  fTextBL->SetTextAlign(13);
  fTextBR->SetTextFont(fFontTex);
  fTextBR->SetTextSize(fSizeTex);
  fTextBR->SetTextColor(fColorTex);
  fTextBR->SetLineWidth(fWidthTex);
  fTextBR->SetTextAlign(33);
  fTextTR->SetTextFont(fFontTex);
  fTextTR->SetTextColor(fColorTex);
  fTextTR->SetTextSize(fSizeTex);
  fTextTR->SetLineWidth(fWidthTex);
  fTextTR->SetTextAlign(31);
  
  TString _curPanel = curPanel;
  map< TString, Panel >::iterator it;
  for( it=_panels.begin(); it!=_panels.end(); it++ )
    {
      //      cout << it->first << endl;
      if( !goToPanel( it->first ) ) continue;
      fCanvas->cd();
      fTextTL->Draw();
      fTextBL->Draw();
      //      fTextBR->Draw();
      fTextTR->Draw();
      TString title = _runStr;	
      title += curPanel;
      fCanvas->SetTitle( title );
      fCanvas->SetName( title );
      fCanvas->Modified();
      fCanvas->Update();
      fPad->cd();
    }  
  goToPanel( _curPanel );
}


void 
MultiDisplay::isolation( float eta_in, float phi_in, float DR_in, float Jhs,
			 float eta_out, float phi_out, float DR_out, 
			 int linestyle, int linewidth, int linecolor )
{
  typedef pair<float,float> point;
  float phi0_ = _phi0*Constants::pi;
  KineUtils::adjust( phi_in, phi0_ );
  KineUtils::adjust( phi_out, phi0_ );

  int N=120;

  if( Jhs>0 )
    {
      // compute the eight point coordinates
      point pt_[8];
      for( int ii=0; ii<4; ii++ )
	{
	  int sgn_ = 1;
	  if( ii==0 || ii==3 ) sgn_=-1;
	  float x0_ = eta_in;
	  float y0_ = phi_in;
	  float R_  = DR_in;
	  if( ii>1 )
	    {
	      x0_ = eta_out;
	      y0_ = phi_out;
	      R_  = DR_out;
	    }
	  float R2_  = R_*R_;
	  for( int jj=0; jj<2; jj++ )
	    {
	      int kk = ii + jj*4;
	      float xx_ = eta_in + (1-2*jj)*Jhs; 
	      float zz_ = R2_ - pow( xx_-x0_, 2 );
	      assert( zz_>0 );
	      float yy_ = y0_ + (1-2*jj)*sgn_*sqrt( zz_ );
	      pt_[kk] = make_pair( xx_, yy_ );
	    }
	}
  
      list<point> v[4];      
      for( int ii=0; ii<4; ii++ )
	{
	  float x0_ = eta_in;
	  float y0_ = phi_in;
	  float R_  = DR_in;
	  int sgnCos_ = (ii==0||ii==1)?-1:1;
	  int sgnSin_ = (ii==0||ii==3)?-1:1;
	  float xth_ = x0_ - sgnCos_*Jhs;
	  if( ii==1 || ii==3 )
	    {
	      x0_ = eta_out;
	      y0_ = phi_out;
	      R_  = DR_out;
	    }
	  for( int jj=0; jj<=N; jj++ )
	    {	  
	      float a_ = 2.*Constants::pi*(jj+0.5)/N;
	      float c_ = cos(a_);
	      float s_ = sin(a_);
	      float x_ = x0_ + sgnCos_*R_*c_;
	      float y_ = y0_ + sgnSin_*R_*s_;
	      if( sgnCos_*x_<sgnCos_*xth_ )
		{
		  v[ii].push_back( make_pair( x_, y_ ) );
		}
	    }
	  v[ii].push_front( pt_[2*ii] );
	  v[ii].push_back( pt_[2*ii+1] );	  
	}
      
      float xx[2][1000];
      float yy[2][1000];
      int ipt[2];
      for( int ii=0; ii<2; ii++ )
	{
	  //      cout << "Arc " << ii << endl;
	  ipt[ii]=0;
	  size_t npt_ = v[2*ii].size()+v[2*ii+1].size()+1;
	  assert( npt_<1000 );
	  for( int jj=0; jj<2; jj++ )
	    {
	      int kk = 2*ii+jj;
	      float curphi_=phi0_;
	      for( list<point>::iterator it=v[kk].begin(); it!=v[kk].end(); ++it ) 
		{
		  point pt = (*it);
		  float eta_ = pt.first;
		  float phi_ = pt.second;
		  KineUtils::adjust( phi_, curphi_ );
		  curphi_=phi_;
		  xx[ii][ipt[ii]] = eta_;
		  yy[ii][ipt[ii]] = phi_/Constants::pi;
		  ipt[ii]++;
		}
	    }      
	  xx[ii][ipt[ii]] = xx[ii][0];
	  yy[ii][ipt[ii]] = yy[ii][0];
	  ipt[ii]++;
	}  
      
      for( int ii=0; ii<2; ii++ )
	{
	  TPolyLine* pline = new TPolyLine( ipt[ii], xx[ii], yy[ii] );
	  registerTObject( pline );
	  //	  pline ->SetLineColor( kBlack );
	  //	  pline ->SetLineWidth( 2 );
	  pline ->SetLineColor( linecolor );
	  pline ->SetLineWidth( linewidth );
	  pline ->SetLineStyle( linestyle );
	  pline ->Draw();
	}
    }
  else
    {
      list<point> v[2];      
      for( int ii=0; ii<2; ii++ )
	{
	  float x0_ = eta_in;
	  float y0_ = phi_in;
	  float R_  = DR_in;
	  if( ii==1  )
	    {
	      x0_ = eta_out;
	      y0_ = phi_out;
	      R_  = DR_out;
	    }
	  for( int jj=0; jj<=N; jj++ )
	    {	  
	      float a_ = 2.*Constants::pi*(jj+0.5)/N;
	      float c_ = cos(a_);
	      float s_ = sin(a_);
	      float x_ = x0_ + R_*c_;
	      float y_ = y0_ + R_*s_;
	      v[ii].push_back( make_pair( x_, y_ ) );
	    }
	}

      float xx[2][1000];
      float yy[2][1000];
      int ipt[2];
      for( int ii=0; ii<2; ii++ )
	{
	  ipt[ii]=0;
	  size_t npt_ = v[ii].size()+1;
	  assert( npt_<1000 );
	  float curphi_=phi0_;
	  for( list<point>::iterator it=v[ii].begin(); it!=v[ii].end(); ++it ) 
	    {
	      point pt = (*it);
	      float eta_ = pt.first;
	      float phi_ = pt.second;
	      KineUtils::adjust( phi_, curphi_ );
	      curphi_=phi_;
	      xx[ii][ipt[ii]] = eta_;
	      yy[ii][ipt[ii]] = phi_/Constants::pi;
	      ipt[ii]++;
	    }
	  xx[ii][ipt[ii]] = xx[ii][0];
	  yy[ii][ipt[ii]] = yy[ii][0];
	  ipt[ii]++;
	}  
      
      for( int ii=0; ii<2; ii++ )
	{
	  TPolyLine* pline = new TPolyLine( ipt[ii], xx[ii], yy[ii] );
	  registerTObject( pline );
	  //	  pline ->SetLineColor( kBlack );
	  //	  pline ->SetLineWidth( 2 );
	  pline ->SetLineColor( linecolor );
	  pline ->SetLineWidth( linewidth );
	  pline ->SetLineStyle( linestyle );
	  pline ->Draw();
	}
    }
}

void
MultiDisplay::jurassic( float eta_in,  float phi_in,
			float DR_in, float Jhs, 
			float eta_out, float phi_out, 
			float DR_out, 
			int linestyle, int linewidth, int linecolor )
{
  assert( Jhs>0 && Jhs<DR_in );      

  // compute the eight point coordinates
  typedef pair<float,float> point;

  float phi0_ = _phi0*Constants::pi;
  KineUtils::adjust( phi_in, phi0_ );
  KineUtils::adjust( phi_out, phi0_ );

  point pt_[8];
  for( int ii=0; ii<4; ii++ )
    {
      int sgn_ = 1;
      if( ii==0 || ii==3 ) sgn_=-1;
      float x0_ = eta_in;
      float y0_ = phi_in;
      float R_  = DR_in;
      if( ii>1 )
	{
	  x0_ = eta_out;
	  y0_ = phi_out;
	  R_  = DR_out;
	}
      float R2_  = R_*R_;
      for( int jj=0; jj<2; jj++ )
	{
	  int kk = ii + jj*4;
	  float xx_ = eta_in + (1-2*jj)*Jhs; 
	  float zz_ = R2_ - pow( xx_-x0_, 2 );
	  assert( zz_>0 );
	  float yy_ = y0_ + (1-2*jj)*sgn_*sqrt( zz_ );
	  pt_[kk] = make_pair( xx_, yy_ );
	}
    }

  int N=120;
  
  list<point> v[4];

  for( int ii=0; ii<4; ii++ )
    {
      float x0_ = eta_in;
      float y0_ = phi_in;
      float R_  = DR_in;
      int sgnCos_ = (ii==0||ii==1)?-1:1;
      int sgnSin_ = (ii==0||ii==3)?-1:1;
      float xth_ = x0_ - sgnCos_*Jhs;
      if( ii==1 || ii==3 )
	{
	  x0_ = eta_out;
	  y0_ = phi_out;
	  R_  = DR_out;
	}
      for( int jj=0; jj<N; jj++ )
	{	  
	  float a_ = 2.*Constants::pi*(jj+0.5)/N;
	  float c_ = cos(a_);
	  float s_ = sin(a_);
	  float x_ = x0_ + sgnCos_*R_*c_;
	  float y_ = y0_ + sgnSin_*R_*s_;
	  if( sgnCos_*x_<sgnCos_*xth_ )
	    {
	      v[ii].push_back( make_pair( x_, y_ ) );
	    }
	}
      v[ii].push_front( pt_[2*ii] );
      v[ii].push_back( pt_[2*ii+1] );

    }

  float xx[2][1000];
  float yy[2][1000];
  int ipt[2];
  for( int ii=0; ii<2; ii++ )
    {
      //      cout << "Arc " << ii << endl;
      ipt[ii]=0;
      size_t npt_ = v[2*ii].size()+v[2*ii+1].size()+1;
      assert( npt_<1000 );
      for( int jj=0; jj<2; jj++ )
	{
	  int kk = 2*ii+jj;
	  for( list<point>::iterator it=v[kk].begin(); it!=v[kk].end(); ++it ) 
	    {
	      //	  cout << "eta=" << it->first << " phi=" << it->second << endl;
	      point pt = (*it);
	      float eta_ = pt.first;
	      float phi_ = pt.second;
	      KineUtils::adjust( phi_, phi0_ );
	      //	      float dphi_ = phi_-phi_in;
	      //	      if( dphi_<-Constants::pi ) phi_+=2*Constants::pi;
	      //	      else if( dphi_>+Constants::pi ) phi_-=2*Constants::pi;

	      xx[ii][ipt[ii]] = eta_;
	      yy[ii][ipt[ii]] = phi_/Constants::pi;
	      //	      cout << " eta=" <<  xx[ii][ipt[ii]]
	      //		   << " phi=" <<  yy[ii][ipt[ii]]
	      //		   << endl;

	      ipt[ii]++;
	    }
	}      
      xx[ii][ipt[ii]] = xx[ii][0];
      yy[ii][ipt[ii]] = yy[ii][0];
      ipt[ii]++;
    }  

  for( int ii=0; ii<2; ii++ )
    {
      TPolyLine* pline = new TPolyLine( ipt[ii], xx[ii], yy[ii] );
      registerTObject( pline );
      //	  pline ->SetLineColor( kBlack );
      //	  pline ->SetLineWidth( 2 );
      pline ->SetLineColor( linecolor );
      pline ->SetLineWidth( linewidth );
      pline ->SetLineStyle( linestyle );
      pline ->Draw();
    }
}

void 
MultiDisplay::drawLine( float x0, float y0, float z0, float x, float y, float z, float R, int lineStyle, int lineWidth, int lineColor )
{
  int sign_ = (R>0)?1:-1;
  R*=sign_;

  float cx0_(0), cy0_(0), cx1_(0), cy1_(0), cx2_(0), cy2_(0);
  float cr_(0), crmax_(0), vx_(0), vy_(0);

  // keep track of the current panel
  TString curPanel_ = curPanel;

  for( size_t iview=0; iview<nCanv; iview++ )
    {
      TString str_ = canvName[iview];
      if( str_=="Pt" )    continue;
      if( !goToPanel( str_ ) ) continue;
 

      if( str_.Contains("XY") )
	{
	  cx0_ = x0;
	  cy0_ = y0;
	  cx2_ = x;
	  cy2_ = y;
	  crmax_ = R;
	}
      else if( str_.Contains("RZ") )
	{
	  cx0_ = z0;
	  cy0_ = sqrt( x0*x0+y0*y0 );
	  cx2_ = z;
	  
	  //      cout << "sign= " << sign_ << endl;
	  
	  cy2_ = sqrt( x*x + y*y );
	  float th_ = atan( cy2_/fabs(z) );
	  if ( th_>0 )
	    crmax_ = R/sin(th_);
	  else 
	    crmax_ = R;
	  
	  cy2_ *= sign_;
	}
      else if( str_.Contains("XZ") )
	{
	  cx0_ = z0;
	  cy0_ = x0;
	  cx2_ = z;
	  cy2_ = x;
	  crmax_ = R;
	}
      else if( str_.Contains("YZ") )
	{
	  cx0_ = z0;
	  cy0_ = y0;
	  cx2_ = z;
	  cy2_ = y;
	  crmax_ = R;
	}
      else if( str_.Contains("PhiZ") )
	{
	  float phi0_ = _phi0 * Constants::pi;
	  float cphi0_ = cos( phi0_ );
	  float sphi0_ = sin( phi0_ );
	  int sgn_ = sphi0_>0 ? 1:-1;
	  cx0_ = z0;
	  cy0_ = sgn_*( cphi0_*x0 + sphi0_*y0 );
	  cx2_ = z;
	  cy2_ = sgn_*( cphi0_*x0 + sphi0_*y0 );
	  crmax_ = R;
	}
      else
	continue;
      vx_ = cx2_-cx0_;
      vy_ = cy2_-cy0_;
      cr_ = sqrt( pow( vx_,2 ) + pow( vy_,2 ) );
      
      //  cout << "cr/vx/vy " << cr_ << "/" << vx_ << "/" << vy_ << endl;
      
      if( cr_==0 ) continue;
      vx_ /= cr_;
      vy_ /= cr_;
      if( cr_>crmax_ ) cr_=crmax_;
      
      //  cout << "cr/vx/vy " << cr_ << "/" << vx_ << "/" << vy_ << endl;
      
      cx1_ = cx0_;
      cx1_ += cr_*vx_;
      cy1_ = cy0_;
      cy1_ += cr_*vy_;
      
      //  cout << "cx1=" << cx1_;
      //  cout << " cy1=" << cy1_ << endl;
      
      TLine* line = new TLine( cx0_, cy0_, cx1_, cy1_ );
      registerTObject( line );
      line->SetLineColor( lineColor );
      line->SetLineWidth( lineWidth );
      line->SetLineStyle( lineStyle );
      line->Draw();
    }
  curPanel = curPanel_;
}

void
MultiDisplay::drawHelix( float x0, float y0, float z0, 
			 float pt0, float eta0, float phi0, 
			 float rmax, float zmax, 
			 float B, 
			 int color, int linewidth, int linestyle )
{
  float rmax2 = rmax*rmax;

  // helix...
  float theta = 2*atan( exp( -eta0 ) );
  if( theta<0 ) theta += 2*Constants::pi;
  float lambda = Constants::pi/2 -theta;
  float tanl   = tan(lambda);
  float ctheta  = cos(theta);
  float stheta  = sin(theta);
  float ttheta  = tan( theta );
  if( ttheta==0 ) 
    {
      cout << "drawHelix : Warning, theta=" << theta << endl;
      return;
    }
  //  float phi0_ = _phi0*Constants::pi;
  float phi   = phi0;
  float cphi  = cos( phi );
  float sphi  = sin( phi );
  float pt    = fabs(pt0);
  assert( pt!=0 );

  float chg   = pt0/pt;
  float vx =  cphi;
  float vy =  sphi;
  float ux = -chg*sphi;
  float uy =  chg*cphi;

  if( B<=0 )
    {
      float x_ = x0 + cphi*stheta*10*rmax;
      float y_ = y0 + sphi*stheta*10*rmax;
      float z_ = z0 +      ctheta*10*rmax;

      int sign_ = ( sphi>0 ) ? 1:-1;

      drawLine( x0, y0, z0, x_, y_, z_, 
		sign_*rmax, linestyle, linewidth, color );
      return;
    }

  float R  =  pt/(0.003*B);
  float cx =  x0 - R*ux;
  float cy =  y0 - R*uy;
  float cz =  z0;
      
  // keep track of the current panel
  TString curPanel_ = curPanel;

  for( size_t iview=0; iview<nCanv; iview++ )
    {
      TString str_ = canvName[iview];
      if( str_=="Pt" ) continue;
      if( !goToPanel( str_ ) ) continue;

      const unsigned maxpt=5000;
      float xh[maxpt];
      float yh[maxpt];
      float d = 0.05/R;
      //float d = 1/R;
      float a = 0;
      int npt = 0;

      // OK, draw the Helix in that view
      for( unsigned ii=0; ii<maxpt; ii++ )
	{
	  //      if( ii+1==maxpt ) cout << "Warning--- Maximum of point along helix reached " << endl;
	  if( a>3*Constants::pi/2. ) break;  // against loopers
	  
	  float cosa = cos(a);
	  float sina = sin(a);
	  float uxa_ = cosa*ux + sina*vx;
	  float uya_ = cosa*uy + sina*vy;
	  
	  // the position
	  float x_ = cx + R*uxa_;
	  float y_ = cy + R*uya_;
	  float z_ = cz + a*R*tanl;
	  
	  float r2 = x_*x_ + y_*y_;
	  if( r2>=rmax2 ) break;
	  if( fabs(z_)>=zmax ) break;
	  

	  if( str_=="EtaPhi" )
	    {
	      return;
		{
		  // 	  // unit vector along the momentum
		  float vxa_ =  chg*uya_*stheta;
		  float vya_ = -chg*uxa_*stheta;
		  float vza_ =  ctheta;
		  
		  float Eta_(0);
		  float Phi_(0);
		  
		  KineUtils::straightToEnvelop( Eta_, Phi_, x_, y_, z_, vxa_, vya_, vza_, rmax, zmax );
		  
		  KineUtils::adjust( Phi_, phi0 ); 
		  //	  xh[ii] = sqrt(r2);
		  //	  yh[ii] = atan2( y_, x_ );
		  xh[ii] = Eta_;
		  yh[ii] = Phi_/Constants::pi;
		  //	  cout << "eta/phi " << xh[ii] << "/" << yh[ii] << endl;
		}
	    }
	  else if( str_.Contains("XY") )
	    {
	      xh[ii] = x_;
	      yh[ii] = y_;
	    }
	  else if( str_.Contains("RZ") )
	    {
	      int sign_ = 1;
	      if( sin( phi0 )<0 ) sign_ = -1;
	      xh[ii] = z_;
	      yh[ii] = sign_*sqrt(r2);
	    }
	  else if( str_.Contains("YZ") )
	    {
	      xh[ii] = z_;
	      yh[ii] = y_;
	    }
	  else if( str_.Contains("XZ") )
	    {
	      xh[ii] = z_;
	      yh[ii] = x_;
	    }
	  // FIXME GHM
	  // 	  else if( str_=="PhiZ" )
	  // 	    {
	  // 	      float phi0_ = _phi0 * Constants::pi;
	  // 	      float cphi0_ = cos( phi0_ );
	  // 	      float sphi0_ = sin( phi0_ );
	  // 	      int sgn_ = sphi0_>0 ? 1:-1;
	  // 	      xh[ii] = z_;
	  // 	      yh[ii] = sgn_ * ( cphi0_ * x_ + sphi0_ * y_ );
	  // 	    }
	  npt++;
	  a += d;
	}
  
      TPolyLine* phline = new TPolyLine( npt, xh, yh );
      registerTObject( phline );

      // FIXME GHM
//       if( str_ == "EtaPhi" )
// 	{
// 	  //      TMarker* marker; 
// 	  float eta1_ = xh[0];
// 	  float phi1_ = yh[0];
// 	  float eta2_ = xh[npt-1];
// 	  float phi2_ = yh[npt-1];

// 	  //       TLine* line = new TLine( eta0, phi1_, eta0, phi2_ );
// 	  //       registerTObject( line );
// 	  //       line->SetLineWidth( 2 );
// 	  //       line->SetLineColor( color );
// 	  //       line->SetLineStyle( linestyle );
// 	  //       line->Draw();
// 	  phline->SetLineColor( color );
// 	  phline->SetLineWidth( 2 );
// 	  phline->SetLineStyle( linestyle );
// 	  phline->Draw();

// 	  //       marker = new TMarker( eta1_, phi1_, 20 );
// 	  //       registerTObject( marker );
// 	  //       marker -> SetMarkerSize( 0.5 );
// 	  //       marker -> SetMarkerColor( kWhite );
// 	  //       marker -> Draw();
// 	  //       marker = new TMarker( eta1_, phi1_, 24 );
// 	  //       registerTObject( marker );
// 	  //       marker -> SetMarkerSize( 0.5 );
// 	  //       marker -> SetMarkerColor( color );
// 	  //       marker -> Draw();
// 	  //       marker = new TMarker( eta2_, phi2_, 21 );
// 	  //       registerTObject( marker );
// 	  //       marker -> SetMarkerSize( 0.5 );
// 	  //       marker -> SetMarkerColor( kWhite );
// 	  //       marker -> Draw();
// 	  //       marker = new TMarker( eta2_, phi2_, 25 );
// 	  //       registerTObject( marker );
// 	  //       marker -> SetMarkerSize( 0.5 );
// 	  //       marker -> SetMarkerColor( color );
// 	  //       marker -> Draw();
// 	  int mark_ = kFullTriangleDown;
// 	  if( chg<0 ) mark_ = kFullTriangleUp;

// 	  drawMarker( eta1_, phi1_, kOpenCircle, 1, color );
// 	  drawMarker( eta2_, phi2_, mark_, 1, color );

// 	  //      cout << "eta0/phi0/eta1/phi1/eta2/phi2 ";
// 	  //      cout << eta0  << "/" << phi0/Constants::pi << " / " ;
// 	  //      cout << eta1_  << "/" << phi1_ << " / " ;
// 	  //      cout << eta2_  << "/" << phi2_ << endl ;

// 	}
//       else
//	{
      phline->SetLineColor( color );
      phline->SetLineWidth( linewidth );
      phline->Draw();
      //	}      
    }  
  curPanel = curPanel_;
}

void 
MultiDisplay::drawVertex( float X_, float Y_, float Z_, int markerColor, bool isPrimary )
{
  // keep track of the current panel
  TString curPanel_ = curPanel;

  for( size_t iview=0; iview<nCanv; iview++ )
    {
      TString str_ = canvName[iview];
      if( str_=="Pt" || str_=="EtaPhi" ) continue;
      if( !goToPanel( str_ ) ) continue;

      if( !str_.Contains("V") && !isPrimary ) continue;

      float xx_(0);
      float yy_(0);
      if( str_.Contains("XY") )
	{
	  xx_ = X_;
	  yy_ = Y_;
	}
      else if( str_.Contains("RZ") )
	{
	  xx_ = Z_;
	  yy_ = sqrt( pow( X_ ,2 ) + pow( Y_, 2 ) );
	  if( Y_<0 ) yy_*=-1;
	}
      else if( str_.Contains("XZ") )
	{
	  xx_ = Z_;
	  yy_ = X_;
	}
      else if( str_.Contains("YZ") )
	{
	  xx_ = Z_;
	  yy_ = Y_;
	}
//       else if( str_.Contains("PhiZ") )
// 	{
// 	  float phi0_ = phi0*Constants::pi;
// 	  xx_ = Z_;
// 	  yy_ = cos(phi0_)*X_ + sin(phi0_)*Y_;
// 	}
      else
	continue;
      
      float vtxSize = 0.5;
      if( isPrimary ) vtxSize = 0.8;


      drawMarkerXY( xx_, yy_, vtxSize, kWhite, markerColor, kFullCircle, kOpenCircle, true );

    }
}

void
MultiDisplay::drawMarker( float x_, float y_, int style, float size, int icol2, int icol1, int shift )
{
//  Marker number         Marker shape          Marker name
//         1                    dot                  kDot
//         2                    +                    kPlus
//         3                    *                    kStar
//         4                    o                    kCircle
//         5                    x                    kMultiply
//         6                    small scalable dot   kFullDotSmall
//         7                    medium scalable dot  kFullDotMedium
//         8                    large scalable dot   kFullDotLarge
//         9 -->19              dot
//        20                    full circle          kFullCircle
//        21                    full square          kFullSquare
//        22                    full triangle up     kFullTriangleUp
//        23                    full triangle down   kFullTriangleDown
//        24                    open circle          kOpenCircle
//        25                    open square          kOpenSquare
//        26                    open triangle up     kOpenTriangleUp
//        27                    open diamond         kOpenDiamond
//        28                    open cross           kOpenCross
//        29                    open star            kOpenStar
//        30                    full star            kFullStar

  int style_1(0);
  int style_2(0);
  bool doublePass_(true);
  if( style==kOpenStar ) 
    {
      style_1 = kFullStar;
      style_2 = kOpenStar;      
    }
  else if( style==kOpenSquare )
    {
      style_1 = kFullSquare;
      style_2 = kOpenSquare;
    }
  else if( style==kOpenCircle )
    {
      style_1 = kFullCircle;
      style_2 = kOpenCircle;
    }
  else
    {
      style_2 = style;
      doublePass_ = false;
    }

  if( x_<_etamin || x_>_etamax ) return;
  y_ += shift;
  if( y_<_phimin )
    {
      if( shift==0 ) return drawMarker( x_, y_, style, size, icol2, icol1, 2 );
      return;
    }
  if( y_>_phimax )
    {
      if( shift==0 ) return drawMarker( x_, y_, style, size, icol2, icol1, -2 );
      return;
    }
  //  TMarker* marker = new TMarker( x_, y_, style );
  //  registerTObject( marker );
  //  marker->SetMarkerSize( size );
  //  marker->SetMarkerColor( icol );
  //  marker->Draw(); 

//   TMarker* marker;

//   if( doublePass_ )
//     {
//       marker = new TMarker( x_, y_, style_1 );
//       registerTObject( marker );
//       marker -> SetMarkerSize( size );
//       marker -> SetMarkerColor( icol1 );
//       marker -> Draw();
//     }

//   marker = new TMarker( x_, y_, style_2 );
//   registerTObject( marker );
//   marker -> SetMarkerSize( size );
//   marker -> SetMarkerColor( icol2 );

//   marker -> Draw();
  drawMarkerXY( x_, y_, size, icol1, icol2, style_1, style_2, doublePass_ );

}
 
void
MultiDisplay::drawMarkerXY( float x_, float y_, float size, int icol1, int icol2, int style_1,  int style_2, bool doublePass_ )
{
  TMarker* marker;

  if( doublePass_ )
    {
      marker = new TMarker( x_, y_, style_1 );
      registerTObject( marker );
      marker -> SetMarkerSize( size );
      marker -> SetMarkerColor( icol1 );
      marker -> Draw();
    }

  marker = new TMarker( x_, y_, style_2 );
  registerTObject( marker );
  marker -> SetMarkerSize( size );
  marker -> SetMarkerColor( icol2 );

  marker -> Draw();
}

void 
MultiDisplay::drawMarkerRhoEtaPhi( float rho, float eta, float phi, float markerSize, int markerColor )
{
  int markerStyle;
  TMarker* marker;

  float x_, y_, z_;
  bool ok = RooUtils::toXYZ( rho, eta, phi, x_, y_, z_ );
  if( !ok ) return;

  // keep track of the current panel
  TString curPanel_ = curPanel;

  for( size_t iview=0; iview<nCanv; iview++ )
    {
      TString str_ = canvName[iview];
      if( str_=="Pt" ) continue;
      if( !goToPanel( str_ ) ) continue;

      float X_(0), Y_(0);

      if( str_.Contains("XY" ) )
	{
	  X_ = x_;
	  Y_ = y_;
	}
      else if( str_.Contains("YZ" ) )
	{
	  Y_ = y_;
	  X_ = z_;
	}
      else if( str_.Contains("XZ" ) )
	{
	  Y_ = x_;
	  X_ = z_;
	}
      else if( str_.Contains("RZ" ) )
	{
	  Y_ = sqrt( x_*x_+y_*y_ );
	  Y_ *= (sin(phi)>0)?1:-1;
	  X_ = z_;
	}
      else if( str_.Contains("PhiZ" ) )
	{
	  float phi0_ = _phi0 * Constants::pi;
	  float cphi0_ = cos( phi0_ );
	  float sphi0_ = sin( phi0_ );
	  int sgn_ = sphi0_>0 ? 1:-1;
	  Y_ = sgn_*( cphi0_*x_ + sphi0_*y_ );
	  X_ = z_;
	}
      else
	continue;
	  
      markerStyle=21;
      marker = new TMarker( X_, Y_, markerStyle );
      registerTObject( marker );
      marker->SetMarkerSize( markerSize );
      marker->SetMarkerColor( kWhite );
      marker->Draw();   

      markerStyle=25;
      marker = new TMarker( X_, Y_, markerStyle );
      registerTObject( marker );
      marker->SetMarkerSize( markerSize );
      marker->SetMarkerColor( markerColor );
      marker->Draw();   
    }
  curPanel = curPanel_;
}

void 
MultiDisplay::drawMET( float et, float phi, float scale, float ptmax, int icol )
{
  if( !goToPanel( "XY" ) ) return;
    {
      double rin(0);
      if( _R==0 )
	{
	  rin = 0.;
	}
      else
	{
	  //	  if( _R>refRadOut )      rin = refRadOut;
	  //	  else if( _R>refRadMid ) rin = refRadMid;
	  //	  else if( _R>refRadIn )  rin = refRadIn;
	  if( _R>refRadOut )      rin = refRadIn;
	  else if( _R>refRadMid ) rin = refRadIn;
	  else if( _R>refRadIn )  rin = refRadIn;
	  else return;
	}

      float scx0_ = rin;
      float scy0_ = 0.;
      float scX0_ = scx0_ + et*scale;
      float scY0_ = 0.;
      float sphi_ = sin( phi );
      float cphi_ = cos( phi );
      double scxx[2];
      double scyy[2];
      scxx[0] =  scx0_*cphi_ + scy0_*sphi_;
      scyy[0] =  scx0_*sphi_ - scy0_*cphi_;
      scxx[1] =  scX0_*cphi_ + scY0_*sphi_;
      scyy[1] =  scX0_*sphi_ - scY0_*cphi_;
      TArrow* metarrow = new TArrow( scxx[0], scyy[0], scxx[1], scyy[1] );
      registerTObject( metarrow );
      metarrow->SetLineColor( icol );
      metarrow->SetLineWidth(4);
      metarrow ->Draw();

      // new !!!
      //      drawLine( hXY, 0, 0, 0, scxx[0], scyy[0], 0, 
      //		rin, kDashed, 2, icol );
      // new !!!
      
      float dphi_ = 5.;
      float phi_  = phi*Constants::radToDeg;
      TArc* metarc = new TArc( 0, 0, rin, phi_-dphi_/2., phi_+dphi_/2. );
      registerTObject( metarc );
      metarc->SetLineColor(icol);
      metarc->SetLineWidth(2);
      metarc->SetFillStyle(0);
      metarc->SetFillColor(0);
      metarc->Draw("only");

    }
//   else if( proj==hEtaPhi )
//     {
//       double x_ = 0;
//       double y_ = phi/Constants::pi;
//       if( y_<0 ) y_+=2;
//       drawMarker( x_, y_, kOpenStar, int( et*scale/50. ), icol );
//     }  
}

void 
MultiDisplay::drawScales( float xpos, float ypos )
{
  // keep track of the current panel
  TString curPanel_ = curPanel;

  for( size_t iview=0; iview<nCanv; iview++ )
    {
      TString str_ = canvName[iview];
      if( str_=="Pt" ) continue;
      if( !str_.Contains("V") ) continue;
      if( !goToPanel( str_ ) ) continue;
      Style& style = curStyle();

      // il faut que la longueur de la scale soit comprise entre
      // 10 et 20% de la largeur 2*_RV
      //      float test1_ =   _RV/5.;
      float test2_ = 2*_RV/5.;
      int ll = 5*int( test2_/5. );

      //      int ll = 10;
      float fact = 1.;
      TString unit = " cm";
      if( ll==0 )
	{
	  // en dessous de 5 cm, de cm en cm
	  ll = int( test2_ );
	  //	  cout << ll << " " << test2_ << endl;
	  if( ll==0 )
	    {
	      // on passe aux millimetres
	      ll = 5*int( 10.*test2_/5. );	      
	      if( ll==0 ) ll=1;
	      fact = 0.1;
	      unit = " mm";
	    }
	}
      else if( ll>100 )
	{
	  // on passe aux metres
	  ll = 5*int( 0.1*test2_/5. );
	  fact = 10.;
	  unit = " m";
	}
      TString txt_;
      txt_ += ll;
      txt_ += unit;
      float length = fact*ll;          

      TLatex ltx_;
      ltx_.SetTextColor( kWhite );
      ltx_.SetTextFont( 42 );   // note 10*font+precision=2
      ltx_.SetTextSize( 0.02 ); // fraction of the current pad size
      ltx_.SetTextAlign( 21 ); 
      
      TLine line_;
      line_.SetLineColor( kWhite );
      line_.SetLineWidth( 2 );

      float xmax;
      float xmin;
      float ymax;
      float ymin;
      if( str_.Contains("XY") )
	{
	  xmax = _xPV+_RV;
	  xmin = _xPV-_RV;
	  ymax = _yPV+_RV;
	  ymin = _yPV-_RV;
	}
      else if( str_.Contains("XZ") )
	{
	  ymin = _xPV-_RV;
	  ymax = _xPV+_RV;
	  xmin = _zPV-_LV;
	  xmax = _zPV+_LV;
	}
      else if( str_.Contains("YZ") )
	{
	  ymin = _yPV-_RV;
	  ymax = _yPV+_RV;
	  xmin = _zPV-_LV;
	  xmax = _zPV+_LV;
	}
      else if( str_.Contains("RZ") )
	{
	  ymin = _rPV-_RV;
	  ymax = _rPV+_RV;
	  xmin = _zPV-_LV;
	  xmax = _zPV+_LV;
	}
      else 
	continue;

      float xm_ = xmin + xpos*(xmax-xmin) - length/2.;
      float xl_ = xm_ - length/2.;
      float xr_ = xm_ + length/2.;
      float ym_ = ymin + ypos*(ymax-ymin);
      float yu_ = ym_ + 0.0075*(ymax-ymin);
      float yd_ = ym_ - 0.0075*(ymax-ymin);
      float yt_ = ym_ + 0.01*(ymax-ymin);

      line_.DrawLine(xl_,ym_,xr_,ym_);
      line_.DrawLine(xl_,yu_,xl_,yd_);
      line_.DrawLine(xr_,yu_,xr_,yd_);
      ltx_.DrawLatex(xm_,yt_,txt_);
      
    }

  curPanel = curPanel_;
}


void 
MultiDisplay::drawCaloTowerHistos( const std::vector< TH1* >& ct_h )
{
  // keep track of the current panel
  TString curPanel_ = curPanel;

  for( size_t iview=0; iview<nCanv; iview++ )
    {
      TString str_ = canvName[iview];
      if( str_=="Pt" ) continue;
      if( str_.Contains("V") ) continue;
      if( !goToPanel( str_ ) ) continue;

      Style& style = curStyle();

      if( str_.Contains("XY") )
	{
	  TH1* ct_pt_h = ct_h[0]; 
	  TH1* ct_em_h = ct_h[1]; 
	  drawPhiHisto( ct_pt_h, _CTScale, 1, style.kCTHadColor, style.kCTHadColor, false, true );
	  drawPhiHisto( ct_em_h, _CTScale, 1, style.kCTEmColor, style.kCTEmColor, false, true );
	}
      else if( str_.Contains("Z") )
	{
	  for( size_t kk=0; kk<4; kk++ )
	    {
	      TH1* ct_pt_h = ct_h[2+2*kk];
	      TH1* ct_em_h = ct_h[2+2*kk+1];
	      drawEtaHisto( ct_pt_h, kk, _CTScale, 1, style.kCTHadColor, style.kCTHadColor, false, true );
	      drawEtaHisto( ct_em_h, kk, _CTScale, 1, style.kCTEmColor, style.kCTEmColor, false, true );
	    }
	}
      else 
	continue;
    }
  curPanel = curPanel_;
}

void
MultiDisplay::drawPhiHisto( TH1* h, float scale, int lineWidth, int lineColor, int fillColor,bool midBin,  bool isInEvent )
{
  assert( h!=0 );
  int nbins = h->GetXaxis()->GetNbins();
  assert( nbins>0 );
  float dphi = 360./nbins;
  float phi0 = 0;
  if( !midBin ) phi0 = dphi/2.;

  float R_(0);
//   if( _R>refRadOut )      R_ = refRadOut;
//   else if( _R>refRadMid ) R_ = refRadMid;
//   else if( _R>refRadIn )  R_ = refRadIn;
//   else return;
  R_ = refRadMid;
  
  for( int ii=1; ii<=nbins; ii++ )
    {
      float val_ = h->GetBinContent( ii );
      if( val_<=0 ) continue;
      drawWedge( dphi/2, R_, R_+val_*scale, phi0+(ii-1)*dphi, lineWidth, lineColor, fillColor, isInEvent ); 
    }
  TArc* arc = new TArc( 0, 0, R_ );
  arc->SetLineColor(kBlack);  
  arc->SetLineWidth(1);  
  arc->SetFillStyle(0);
  arc->SetFillColor(0);
  arc->Draw();
  registerTObject( arc );
}

void
MultiDisplay::drawEtaHisto( TH1* h, int jj, float scale, int lineWidth, int lineColor, int fillColor, bool midBin,  bool isInEvent )
{
  assert( jj>=0 && jj<4 );
  assert( h!=0 );
  int nbins = h->GetXaxis()->GetNbins();
  assert( nbins>0 );

  //  cout << "jj= " << jj << endl;
  int sign1_ =  (jj<2)?   1:-1;
  int sign2_ = (jj%2==0)? 1:-1;

  for( int ii=1; ii<=nbins; ii++ )
    {
      float val_ = h->GetBinContent( ii );
      if( val_<=0 ) continue;
      //      cout << "val[" << ii << "]=" << val_ << endl; 

      //      float scale_ = scale; 
      float scale_ = (val_/200.)*scale; 
      //float scale_ = 0.1;

      float r_[5];
      float z_[5];

      float r_in_ = rEcalIn+40; //!!!
      float z_in_ = ct_zin1_+30; //!!!
      if( ii<18 )
	{
	  float a1 = r_in_/ct_rin_;
	  r_[0] = a1*sign2_*ct_rin_;
	  r_[1] = a1*sign2_*ct_rin_*(1+scale_);
	  r_[2] = r_[1];
	  r_[3] = r_[0];
	  r_[4] = r_[0];
	  z_[0] = a1*sign1_*ct_zin_[ii-1];
	  z_[1] = a1*sign1_*ct_zin_[ii-1]*( 1 + scale_ );
	  z_[2] = a1*sign1_*ct_zin_[ii]*( 1 + scale_ );
	  z_[3] = a1*sign1_*ct_zin_[ii];
	  z_[4] = z_[0];
	}
      else if( ii<30 )
	//      else
	{
	  float a1 = z_in_/ct_zin1_;
	  //	  float a2 = ct_zout_[29]*rHcalOut/ct_rout_/ct_zout1_;
	  z_[0] = a1*sign1_*ct_zin1_;
	  z_[1] = z_[0];
	  z_[2] = a1*sign1_*ct_zin1_*(1+scale_);
	  z_[3] = z_[2];
	  z_[4] = z_[0];
	  r_[0] = a1*sign2_*ct_rin1_[ii-17];
	  r_[1] = a1*sign2_*ct_rin1_[ii-17-1];
	  r_[2] = a1*sign2_*ct_rin1_[ii-17-1]*(1+scale_);
	  r_[3] = a1*sign2_*ct_rin1_[ii-17]*(1+scale_);
	  r_[4] = r_[0];
	}
      else
	{
	  float a1 = z_in_/ct_zin2_;
	  z_[0] = a1*sign1_*ct_zin2_;
	  z_[1] = z_[0];
	  z_[2] = a1*sign1_*ct_zin2_*(1+scale_);
	  z_[3] = z_[2];
	  z_[4] = z_[0];
	  r_[0] = a1*sign2_*ct_rin2_[ii-29];
	  r_[1] = a1*sign2_*ct_rin2_[ii-29-1];
	  r_[2] = a1*sign2_*ct_rin2_[ii-29-1]*(1+scale_);
	  r_[3] = a1*sign2_*ct_rin2_[ii-29]*(1+scale_);
	  r_[4] = r_[0];
	}
      TPolyLine* pline = new TPolyLine( 5, z_, r_ );
      if( isInEvent ) registerTObject( pline );
      pline ->SetFillColor( fillColor );
      pline ->SetLineColor( lineColor );
      pline ->SetLineWidth( lineWidth );
      pline ->Draw("f");
      pline ->Draw();

    }
}

void
MultiDisplay::drawRadialArrow( float l, float phi, 
			       float dphi, float r1, float r2,
			       int lineColor, 
			       int lineWidth, 
			       int fillColor, 
			       int fillStyle,
			       float x0, float y0
			       )
{
  // warning, phi and dphi are in radians
  //  dphi *= Constants::pi/180.;
  double xx[7];
  double yy[7];
  xx[0] = 0; yy[0] = 0;
  xx[1] = r1*l*cos(dphi*(1-r2)/2.);        
  yy[1] = r1*l*sin(dphi*(1-r2)/2.);
  xx[2] = r1*l*cos(dphi/2.); 
  yy[2] = r1*l*sin(dphi/2.);
  xx[3] = l; yy[3] = 0;
  xx[4] = xx[2]; yy[4] = -yy[2];
  xx[5] = xx[1]; yy[5] = -yy[1];
  xx[6] = xx[0]; yy[6] = yy[0];
  for( int ii=0; ii<7; ii++ )
    {
      float xx_ = xx[ii]; float yy_ = yy[ii];
      xx[ii] = x0 + cos(phi) * xx_ - sin(phi) * yy_;
      yy[ii] = y0 + sin(phi) * xx_ + cos(phi) * yy_;
    }
  //  TLine* line = new TLine( xx[0], yy[0], xx[3], yy[3] );
  //  line->SetLineColor( lineColor );
  //  line->SetLineStyle( kSolid );
  //  line->Draw();
  //  RegisterTObject( line );
  TPolyLine* pline = new TPolyLine( 7, xx, yy );
  pline->SetLineColor( lineColor );
  pline->SetLineStyle( kSolid );
  pline->SetLineWidth( lineWidth );
  pline->SetFillStyle( fillStyle );
  pline->SetFillColor( fillColor );
  pline->Draw("f");
  pline->Draw();
  registerTObject( pline );
}

void
MultiDisplay::print( TString fileType )
{
  for( size_t i=0; i<nCanv; i++ )
    {
      Panel& panel_ = _panels[ canvName[i] ];
  
      TCanvas* canv_ =  panel_.canv;
      if( canv_ ==0 ) continue;
      
      canv_->Print( fileType );
      
    }
}


void 
MultiDisplay::refreshGlobalEcalHist( float Emax )
{
  delete _globalEcalHist;
  _globalEcalHist = MEGeom::getGlobalHist("Global_Ecal_Hist");
  RooUtils::setTdrHistoStyle( _globalEcalHist );
  _globalEcalHist->SetMinimum( 0    );
  _globalEcalHist->SetMaximum( Emax );
}

void 
MultiDisplay::fillGlobalEcalHist( int ix, int iy, int iz, float val )
{
  assert( _globalEcalHist!=0 );
  MEGeom::setBinGlobalHist( _globalEcalHist, ix, iy, iz, val );
}

void 
MultiDisplay::drawGlobalEcalHist( string opt )
{
  _globalEcalHist->Draw( opt.c_str() );
  MEGeom::drawGlobalBoundaries( kGray );
}

void 
MultiDisplay::drawElectronEcalHist( float pt, float eta, float phi,
				    vector<int>& v_ix, 
				    vector<int>& v_iy,
				    vector<int>& v_iz,
				    vector<float>& v_e,
				    int hist_size, string opt  )
{
  if( hist_size<2 ) hist_size=2;
  bool isLego = strstr( opt.c_str(),"LEGO" )!=0;

  refreshGlobalEcalHist( pt );

  TH2* hist_ = (TH2*) _globalEcalHist->Clone();

  size_t nh_ = v_e.size();
  assert( v_ix.size()==nh_ );
  assert( v_iy.size()==nh_ );
  assert( v_iz.size()==nh_ );
  float dRmin=0.5;
  float Emax=0.;
  size_t iimin = 0;
  list<int>   v_binx;
  list<int>   v_biny;
  for( size_t ii=0; ii<nh_; ii++ )
    {
      int ix = v_ix[ii];
      int iy = v_iy[ii];
      int iz = v_iz[ii];
      float e  = v_e[ii];
      
      //      fillGlobalEcalHist( ix, iy, iz, e );
      MEGeom::setBinGlobalHist( hist_, ix, iy, iz, e );

      pair<float, float> etaPhiPoint;
      if( iz==0 ) 
	{
	  etaPhiPoint = MEEBDisplay::center( ix, iy );
	}
      else
	{
	  etaPhiPoint = MEEEDisplay::center( ix, iy, iz );
	}
      float etaXtal = etaPhiPoint.first;
      float phiXtal = etaPhiPoint.second * Constants::pi ; //!!!!!
      float dR = KineUtils::dR( eta, etaXtal, phi, phiXtal );

//       if( e>10 )
// 	{
// 	  cout << "ix/iy/iz/e " << ix << "/" << iy << "/" << iz << "/" << e;
// 	  cout << " ---> eta/phi/etaXtal/phiXtal/dR " << etaXtal << "/" << phiXtal << "/" << dR;
// 	  cout << endl;
// 	}

      if( dR>dRmin ) continue;
      int ibinx(0);
      int ibiny(0);
      MEGeom::getBinGlobalHist( hist_, ix, iy, iz, ibinx, ibiny );

//       if( dR<dRmin )
// 	{
// 	  dRmin = dR;
// 	  iimin = ii;
// 	  v_binx.push_front( ibinx );
// 	  v_biny.push_front( ibiny );
// 	}
//       else
// 	{
// 	  v_binx.push_back( ibinx );
// 	  v_biny.push_back( ibiny );
// 	}
      if( e>Emax)
	{
	  Emax = e;
	  iimin = ii;
	  v_binx.push_front( ibinx );
	  v_biny.push_front( ibiny );
	}
      else
	{
	  v_binx.push_back( ibinx );
	  v_biny.push_back( ibiny );
	}
    }
  int ibinx_0 = v_binx.front();
  int ibiny_0 = v_biny.front();
  int iz = v_iz[iimin];

  MEGeom::getBinGlobalHist( hist_, v_ix[iimin], v_iy[iimin], v_iz[iimin], ibinx_0, ibiny_0 );

  if( Emax>0 )
    {
      //      cout << "ix/iy/iz/e " << v_ix[iimin] << "/" << v_iy[iimin] << "/" << v_iz[iimin] << "/" << v_e[iimin] << " dR=" << dRmin << endl;
    }
  else
    {
      cout << "No closeby xtal signal found " << endl;
      return;
    }      


  //  float Emax = pt;
  //  if( v_e[iimin]>Emax ) Emax = v_e[iimin];
  hist_->SetMaximum( Emax );
  
  //  int hist_size = 20;

  TAxis* ax = hist_->GetXaxis();
  TAxis* ay = hist_->GetYaxis();
  TAxis* az = hist_->GetZaxis();

  int ibinx_min = ibinx_0 - hist_size;
  int ibinx_max = ibinx_min + 2*hist_size+1;
  if( ibinx_min<1 ) 
    {
      ibinx_min = 1;
      ibinx_max = ibinx_min + 2*hist_size+1;
    }
  if( ibinx_max>=ax->GetNbins() )
    {
      ibinx_max = ax->GetNbins();
      ibinx_min = ibinx_max - 2*hist_size-1;
    }
  float ax_min = ax->GetBinLowEdge( ibinx_min );
  float ax_max = ax->GetBinUpEdge( ibinx_max );
  ax->SetNdivisions( 3 );
  ax->SetRangeUser( ax_min, 
		    ax_max );

  int ibiny_min = ibiny_0 - hist_size;
  int ibiny_max = ibiny_min + 2*hist_size+1;
  if( ibiny_min<1 ) 
    {
      ibiny_min = 1;
      ibiny_max = ibiny_min + 2*hist_size+1;
    }
  if( ibiny_max>=ay->GetNbins() )
    {
      ibiny_max = ay->GetNbins();
      ibiny_min = ibiny_max - 2*hist_size-1;
    }
  float ay_min = ay->GetBinLowEdge( ibiny_min );
  float ay_max = ay->GetBinUpEdge( ibiny_max );
  ay->SetNdivisions( 3 );
  ay->SetRangeUser( ay_min, 
		    ay_max );
  if( iz==0 )
    {
      ax->SetTitle( "ieta" );
      ay->SetTitle( "iphi" );
    }
  else 
    {
      ay->SetTitle( "iy" );
      ax->SetTitle( "ix" );
    }
  ax->CenterTitle();
  ay->CenterTitle();

//   cout << ibinx_0 << "/" << ibiny_0 << endl;
//   cout << "x-->" << ibinx_min << "/" << ibinx_max << endl;
//   cout 
//     << "y-->" << ibiny_min << "/" << ibiny_max 
//     << "-->" << ay_min << "/" << ay_max 
//     << endl;

  if( isLego )
    {
      az->SetTitle( "energy (GeV)" );
      az->CenterTitle();
      az->SetNdivisions( 5 );
      az->SetLabelOffset( 0.02 );
      az->SetTitleOffset( 1.30 );

      ax->SetLabelSize(0.02);
      ay->SetLabelSize(0.02);
      ax->SetLabelColor(kGray+1);
      ay->SetLabelColor(kGray+1);
      ax->SetNdivisions( 5 );
      ay->SetNdivisions( 5 );

      ax->SetTitleOffset( 1.2 );
      ay->SetTitleOffset( 1.2 );
    }
  hist_->Draw( opt.c_str() );
  
  if( !isLego ) MEGeom::drawGlobalBoundaries( kGray );

}

void
MultiDisplay::setEtaPhiLimits( float eta0, float phi0, float DE, float DPh )
{
  _DE      = DE;
  _DPh     = DPh;
  _eta0    = eta0;
  _phi0    = phi0;
  _etamin  = _eta0 - _DE;
  _etamax  = _eta0 + _DE;
  _phimin  = _phi0 - _DPh;
  _phimax  = _phi0 + _DPh;
}

void
MultiDisplay::setXYZLimits( float x0, float y0, float z0, float R )
{
  _x0 = x0;
  _y0 = y0;
  _z0 = z0;
  _r0 = sqrt( _x0*_x0+_y0*_y0 );
  _R = R;
  _L = 2*_R;
}

void
MultiDisplay::setVertexLimits( float xPV, float yPV, float zPV, float RV )
{
  _xPV = xPV;
  _yPV = yPV;
  _zPV = zPV;
  _rPV = sqrt( xPV*xPV + yPV*yPV );
  _RV = RV;
  _LV = 2*_RV;
}

void
MultiDisplay::setEnergyLimits( float ptmax, float Emax, float CTScale )
{
  _CTScale = CTScale;
  _ptmax   = ptmax;
  _emax    = Emax;
}

void
MultiDisplay::drawEcalRecHit( int ix, int iy, int iz, float e_, float Emax, float emin, float etmin )
{ 
  int det=-1;
  if( iz==0 )
    {
      if( ix<0 )
	det = ME::iEBM;
      if( iy>0 )
	det = ME::iEBP;
    }
  else if( iz==-1 ) det = ME::iEEM;
  else if( iz==1  ) det = ME::iEEP;
  assert( det!=-1 );

  float eta_(0), phi_(0), theta_(0);
  float et_=e_;
  if( det==ME::iEBM || det==ME::iEBP )
    {
      int ieta = ix;
      int iphi = iy;
      MEEBGeom::EtaPhiPoint point_ = MEEBDisplay::center( ieta, iphi );
      eta_ = point_.first;
      phi_ = point_.second;
      theta_ = 2*atan( exp( -eta_ ) );
      et_ = e_*sin(theta_);
    }
  else
    {
      MEEEGeom::EtaPhiPoint point_ = MEEEDisplay::center( ix, iy, iz );
      eta_ = point_.first;
      phi_ = point_.second;
      theta_ = 2*atan( exp( -eta_ ) );
      et_ = e_*sin(theta_);
    }

  int icol = 19;
  if( e_>0 ) 
    {
      if( e_>emin && et_>etmin )
	{	      
	  icol = RooUtils::color( e_, Emax );
	}
      else
	{
	  icol = kGray;
	}
    }
  
  if( det==ME::iEBM || det==ME::iEBP )
    {
      int ieta = ix;
      int iphi = iy;
      //      cout << "\t\tieta/iphi " << ieta << "/" << iphi
      //      	       << fixed  << " E=" << e_ << endl;
      drawEBXtal( ieta, iphi, icol, 0. , _etamin, _etamax, _phimin, _phimax );
    }
  else
    {
      //      cout << "\t\tix/iy/iz " << ix << "/" << iy << "/" << iz
      //	   << fixed  << " E=" << e_ << endl;
      drawEEXtal( ix, iy, iz, icol, 0. , _etamin, _etamax, _phimin, _phimax );
    }	  	  
}

void 
MultiDisplay::drawEBXtal( MEEBGeom::EBGlobalCoord ieta, 
		       MEEBGeom::EBGlobalCoord iphi, 
		       int color, float shift, float etamin, float etamax, float phimin, float phimax )
{
  TPolyLine* pline = MEEBDisplay::getXtalPolyLine( ieta, iphi, shift, _etamin, _etamax, _phimin, _phimax );
  if( pline!=0 ) 
    {
      registerTObject( pline );
      
      pline->SetFillColor( color  );
      pline->SetLineColor( color );
      pline->SetLineWidth( 0 );
      pline->Draw("f"); pline->Draw();
      
    }
  if( shift==0 )
    {
      //      int iSM = MEEBGeom::sm( ieta, iphi );
      //      if( iSM==1 || iSM==19 )  
      drawEBXtal( ieta, iphi, color, 2.,  _etamin, _etamax, _phimin, _phimax );
      drawEBXtal( ieta, iphi, color, -2., _etamin, _etamax, _phimin, _phimax );
    }
}

void 
MultiDisplay::drawEEXtal( int ix, int iy, int iz, int color, float shift, float etamin, float etamax, float phimin, float phimax )
{
  TPolyLine* pline = MEEEDisplay::getXtalPolyLine( ix, iy, iz, shift, _etamin, _etamax, _phimin, _phimax );
  if( pline!=0 ) 
    {
      registerTObject( pline );
      
      pline->SetFillColor( color  );
      pline->SetLineColor( color );
      pline->SetLineWidth( 0 );
      pline->Draw("f"); pline->Draw();
      
    }
  if( shift==0 )
    {
      //      int iX = (ix-1)/5+1;
      //      int iY = (iy-1)/5+1;
      //      if( iY==10 && iX>=11 ) 
      drawEEXtal( ix, iy, iz, color, -2., _etamin, _etamax, _phimin, _phimax ); 
      //      if( iY==11 && iX>=11 ) 
      drawEEXtal( ix, iy, iz, color,  2., _etamin, _etamax, _phimin, _phimax );
    }
}


void
MultiDisplay::drawCaloTower_( int ieta, int iphi, float energy, float emEt, float hadEt, float Emax )
{
  float et_ = emEt+hadEt;
  if( et_<=0 ) return;
  if( energy<=0 ) return;

  if( ieta==0 || abs(ieta)>41 )
    {
      return; 
    }
  int ii= abs(ieta); 

  float e0 = ct_eta_[ ii-1 ]; 
  float e1 = ct_eta_[ ii   ];       
  if( ieta<0 )
    {
      e0*=-1;  e1*=-1; 
    }

  float ph0 = (iphi-1)*2./72; 
  float dph = 2./72; 
  if( ii>20 ) 
    {
      if( (iphi-1)%2!=0 ) return; 
      dph *= 2; 
    }
  if( ii>39 ) 
    {
      if( (iphi-1)%4!=0 ) return;      
      dph *= 2; 
    }
  float ph1 = ph0 + dph; 

  //  int lineWidth = 2;
  //  int lineColor = kRed;
  //  int fillColor = 0;
  //  drawBox( e0, e1, ph0, ph1 , lineWidth, lineColor, fillColor, true );
  float eta_ = 0.5*(e0+e1);
  float phi_ = 0.5*(ph0+ph1);

  drawCaloTower( eta_, phi_, energy, emEt, hadEt, Emax );
}

void
MultiDisplay::drawCaloTower( float eta, float phi, float energy, float emEt, float hadEt, float Emax )
{
  float eta_ = eta;
  float phi_ = phi;

  int icol1 = RooUtils::color( energy* hadEt/(emEt+hadEt), Emax );
  int icol2 = RooUtils::color( energy* emEt/(emEt+hadEt), Emax );

  int isize1 = 3;
  int isize2 = 1;
  
  if( hadEt>0 )
    drawMarker( eta_, phi_, kOpenSquare, isize1, kBlack, icol1 );      
  if( emEt>0 )
    drawMarker( eta_, phi_, kOpenCircle, isize2, kBlack, icol2 );      
}
