#ifndef MultiDisplay_hh
#define MultiDisplay_hh

//
// Author  : Gautier Hamel de Monchenault
//

//
// AnaNaS Event Display 
//

// std headers
#include <string>
#include <vector>
#include <map>
#include <list>

// root headers
#include <TROOT.h>
#include <TString.h>

#include <TCanvas.h>
#include <TPad.h>
#include <TVirtualPadEditor.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <TColor.h>

// MusEcal headers
#include "MECore/src/MEGeom.hh"
#include "MECore/src/MEEBDisplay.hh"
#include "MECore/src/MEEEDisplay.hh"

// Canvases

class MultiDisplay
{
public:

  enum Canv { cXY=0, cRZ, cXZ, cYZ, cVXY, cVRZ, cVXZ, cVYZ, cEtaPhi, cPt, cHist, nCanv };
  static TString canvName[nCanv];  

  struct Style;
  struct Panel;

  // contructor/destructor
  MultiDisplay();
  virtual ~MultiDisplay() {}

  // refresh drawing
  void refresh();
  void refreshDet();

  // set run & event number
  bool setEvent( int run, int event );

  // general panel management
  bool refreshPanel(   TString, bool force=false );
  bool goToPanel(      TString, bool create=false );

  void refreshAllPanels();
  void displayAllPanels();

  static int   showDet;
  static float refRadOut;
  static float refRadMid;
  static float refRadIn;
  static float caloRad;

  // detector drawing and dimension setting
  void XY( TString view="XY" );
  void RZ( TString view="RZ" );
  void Pt();
  void EtaPhi();
  //  void EtaPhi();

  // central tracker
  static float rPipe   ;
  static float rPixOut ;
  static float zPixOut ;
  static float rPixel[3]  ;
  static float rTibIn ;
  static float rTibOut ;  
  static float rTobIn ;
  static float rTobOut ;
  static float zTobOut ;
  static float zTidIn ;
  static float zTidOut ;
  static float zTibOut ;
  static float zTecIn ;
  static float zTecOut ;
  static float rTecOut ;
  static float rCtOut ;
  static float zCtOut ;
  
  // ecal
  static float rEcalIn ;
  static float rEcalOut;
  static float zEcalOut;

  // hcal
  static float rHcalIn ;
  static float rHcalOut;
  static float rHcalExt;

  // coil
  static float rCoilIn ;
  static float rCoilOut;
  static float zCoilOut;

  // yoke
  static float rYoke1In ;
  static float rYoke1Out;
  static float rYoke2In ;
  static float rYoke2Out;
  static float rYoke3In ;
  static float rYoke3Out;

  // primitives
  void drawBox( float xmin, float xmax, float ymin, float ymax, 
		int lineWidth, int lineColor, int fillColor,
		bool isInEvent=true  ); 
  void drawWedge( float dd, float rin, float rout, float phi0, 
		  int lineWidth, int lineColor, int fillColor,
		  bool isInEvent=true  ); 
  void drawPGon( float dd, float rin, float rout, float phi0, 
		 int lineWidth, int lineColor, int fillColor,
		 bool isInEvent=true  ); 
  bool isInView( TPolyLine* pl, float xmin, float xmax, float ymin, float ymax );
  void drawEtaPhiDet();
  void drawEtaPhiEcal();

  //
  // Drawing functions
  //

  // isolation in eta-phi
  void isolation( float eta_in, float phi_in, float DRin, float Jhs,
		  float eta_out, float phi_out, float DR_out,
		  int linestyle=0, int linewidth=1, int linecolor=1 );

  void jurassic( float eta_in, float phi_in, float DRin, float Jhs,
		 float eta_out, float phi_out, float DR_out,
		 int linestyle=0, int linewidth=1, int linecolor=1 );

  void drawHelix( float x0, float y0, float z0, 
		  float pt0, float eta0, float phi0, 
		  float rmax, float zmax, float B, 
		  int color, int linewidth, int linestyle );

  void drawVertex( float x0, float y0, float z0, int markerColor, bool isPrimary );

  void drawEBXtal( int ieta, int iphi, int color, float shift=0, float etamin=-1000, float etamax=1000, float phimin=-1000, float phimax=1000 );

  void drawEEXtal( int ix, int iy, int iz, int color, float shift=0, float etamin=-1000, float etamax=1000, float phimin=-1000, float phimax=1000 );

  void drawEcalRecHit( int ix, int iy, int iz, float e, float Emax=30, float emin=-1000, float etmin=-1000 );

  void drawCT( int ieta, int iphi, int lineWidth, int lineColor, int fillColor, bool isInEvent=true );  

  void drawCaloTower( float eta, float phi, float energy, float emEt, float hadEt, float Emax );  

  void drawCaloTower_( int ieta, int iphi, float energy, float emEt, float hadEt, float Emax );  
  
  void drawCaloTowerHistos( const std::vector< TH1* >& ct_h );

  void drawLine( float x0, float y0, float z0, float x, float y, float z, float R, int lineStyle, int lineWidth, int lineColor );

  void drawPhiHisto( TH1* h, float scale, int lineWidth, int lineColor, int fillColor, bool midBin=false, bool isInEvent=true );

  void drawEtaHisto( TH1* h, int jj, float scale, int lineWidth, int lineColor, int fillColor, bool midBin=false, bool isInEvent=true );

  void drawRadialArrow( float l, float phi,
			float dphi=10, float ratio1=0.80, float ratio2=0.60,
			int lineColor=kRed, 
			int lineWidth=4, 			
			int fillColor=kWhite, 
			int fillStyle=3001,
			float x0=0, float y0=0
			);

//   void drawArrow( float x1, float y1, float x2, float y2, 
// 		  int lineColor=kRed, 
// 		  int lineWidth=4, 			
// 		  int lineStyle=kSolid,
// 		  float frac=0.3
// 		  );

  void drawMarkerXY( float x_, float y_, float size, int icol1, int icol2, int style_1, int style_2, bool doublePass=true );

  void drawMarkerRhoEtaPhi( float rho, float eta, float phi, float markerSize, int markerColor );

  void drawMarker( float x_, float y_, int style, float size, int icol1, int icol2=kWhite, int shift=0 );

  void drawMET( float et, float phi, float scale, float ptmax, int col=kBlue );

  void drawScales( float xpos=0.95, float ypos=0.05 );

//   void drawJet( float et, float eta, float phi, float ptmax, int icol );

//   void drawCaloJet( float et, float eta, float phi, float area, float emf, float scale, float ptmax );

//   // ecal hists
  void refreshGlobalEcalHist( float Emax );
  void fillGlobalEcalHist( int ix, int iy, int iz, float val );
  void drawGlobalEcalHist( std::string opt );
  void drawElectronEcalHist( float pt, float eta, float phi,
 			     std::vector<int>& v_ix, 
 			     std::vector<int>& v_iy,
 			     std::vector<int>& v_iz,
 			     std::vector<float>& v_e,
 			     int hist_size, std::string opt  );
  


  // calo towers
  float ct_eta_[42];
  float ct_etaout_[42];
  float ct_rin_;
  float ct_rout_;
  float ct_zin_[18];
  float ct_zout_[18];
  float ct_zin1_;
  float ct_zout1_;
  float ct_rin1_[13];
  float ct_rout1_[13];
  float ct_zin2_;
  float ct_zout2_;
  float ct_rin2_[13];
  float ct_rout2_[13];
  void setCtLimits();



  void print( TString fileType=".pdf" );

public:

  struct Style
  {
    int kCanvasFillStyle;
    int kCanvasBorderSize;
    int kCanvasFillColor;
    int kPadFillStyle;
    int kPadBorderSize;
    int kPadFillColor;
    int kPadFrameFillColor;
    int kPadFrameLineColor;
    int kBeamAxisColor;
    int kAxisColor;
    int kLineWidth;
    int kPipeFillColor;
    int kPipeLineColor;
    int kPixelFillColor;
    int kPixelLineColor;
    int kCtFillColor;
    int kCt2FillColor;
    int kCtLineColor;
    int kEcalFillColor;
    int kEcalLineColor;
    int kXtalLineColor;
    int kXtalFillColor;
    int kHcalFillColor;
    int kHcalLineColor;
    int kCoilFillColor;
    int kCoilLineColor;
    int kYokeFillColor;
    int kYokeLineColor;
    int kPtLineColor;
    int kPtTextColor;
    int kPtFillColor;
    int kCTHadColor;
    int kCTEmColor;
  };

  // styles
  std::vector< Style > _style;
  
  struct Panel
  {
    TCanvas*    canv;          // the canvas
    TPad*       pad;           // the pad
    TH2*        hist;          // the histogram
    TString     title;         // the title
    unsigned    x;             // x top left
    unsigned    y;             // y top left
    unsigned    w;             // width
    unsigned    h;             // height
    Style       style;         // white or black
    bool        keep;          // if true, do not close at each event 
  };

  //  std::vector< Panel > _panels;
  std::map< TString, Panel > _panels;

  void setPanel( TString title, 
		 unsigned x, unsigned y,
		 unsigned w, unsigned h,
		 Style style,
		 bool keep
		 );

  void setSessionStyle();

  TCanvas* fCanvas;
  TPad* fPad;
  TH2*  fHist;
  TString curPanel;
  Style& curStyle() { return _panels[curPanel].style; }

  void setXYZLimits( float x0=0, 
		     float y0=0, 
		     float z0=0,
		     float R=260. );
  float _x0;
  float _y0;
  float _r0;
  float _z0;
  float _R;
  float _L;

  void setVertexLimits( float xPV=0,
			float yPV=0,
			float zPV=0,
		        float RV=30 );
  float _xPV;
  float _yPV;
  float _rPV;
  float _zPV;
  float _RV;
  float _LV;

  void setEtaPhiLimits( float eta0=0, 
			float phi0=1, 
			float DE=5.5, 
			float DPh=1.3 );
  float _DE;
  float _DPh;
  float _eta0;
  float _etamin;
  float _etamax;
  float _phi0;
  float _phimin;
  float _phimax;

  void setEnergyLimits( float ptmax=110., float Emax=20., float CTScale=1. );
  float _CTScale;
  float _ptmax;
  float _emax;

  int curRun;
  int curEvent;
  TString _runStr;
  TString _runStrLong;

  TLatex *fTextTL, *fTextTR, *fTextBL, *fTextBR;
  void drawComments();
  void drawComments( const char* textTL, 
		     const char* textBL, 
		     const char* textTR, 
		     const char* textBR );

  TH2* _globalEcalHist;

  // garbage collection
  std::list< TObject* > _list;
  void registerTObject( TObject* );

ClassDef(MultiDisplay,0) // MultiDisplay -- Ecal event display

};

#endif

