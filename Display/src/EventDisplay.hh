#ifndef EventDisplay_hh
#define EventDisplay_hh

//
// Authors  : Gautier Hamel de Monchenault and Julie Malcles, Saclay 
//

//
// Simple CMS Event Display from Reco Collections
//   based on ROOT input files created by CMSSW module GhmNtupleMaker
//   (source code available in CVS in:
//         UserCode/GautierHdeM/src/GhmAnalysis/GhmNtupleMaker)
//

// std headers
#include <string>
#include <vector>
#include <map>
#include <list>

// root headers
#include <TROOT.h>
#include <TChain.h>
#include <TH2I.h>
#include <TH2F.h>
#include <TFile.h>
#include <TString.h>

// MusEcal headers
#include "MECore/src/MECanvasHolder.hh"
#include "MECore/src/MEGeom.hh"
#include "MECore/src/MEEBDisplay.hh"
#include "MECore/src/MEEEDisplay.hh"

// views/projections
enum View { hXY=1, hYZ, hXZ, hPhiZ, hRZ, hPt, hEtaPhi, hEcal };  

class EventDisplay : public MECanvasHolder
{

public:

  // contructor/destructor
  EventDisplay( const char* name="AnaNaS", int ww=800, int wh=800 );
  virtual ~EventDisplay();

  // event navigation
  bool event( int run, int event );

  static int   showDet;
  static float Ptcut;
  static float McPtcut;
  static float McEtacut;
  static float Wghcut;
  static float JetEtcut;
  static float delta;
  static float refRadOut;
  static float refRadMid;
  static float refRadIn;
  static float caloRad;
  static float length;

  void primer();

  void print( const char* name_=0 );

  // isolation in eta-phi
  void isolation( float eta_in, float phi_in, float DRin, float Jhs,
		  float eta_out, float phi_out, float DR_out,
		  int linestyle=0, int linewidth=1, int linecolor=1 );
  void jurassic( float eta_in, float phi_in, float DRin, float Jhs,
		 float eta_out, float phi_out, float DR_out,
		 int linestyle=0, int linewidth=1, int linecolor=1 );

  // refresh drawing
  void refresh();
  void refreshDet();

  // predefined views
  void xyView( float R=350,  float x0=0, float y0=0, float Emax=30 ); 
  void rzView( float L=600, float aspectRatio=1, float r0=0, float z0=0, float Emax=30 );
  void xzView( float L=600, float aspectRatio=1, float x0=0, float z0=0, float Emax=30 );
  void yzView( float L=600, float aspectRatio=1, float y0=0, float z0=0, float Emax=30 );
  void etaPhiView( float Deta=3.2, float Dphi=1.2, float eta0=0, float phi0=1, float Emax=30 );
  void ptView( float ptmax=100., float Emax=30 );

  // detector drawing and dimension setting
  void xyPlane( float R=350, float x0=0, float y0=0, float Emax=30 );
  void rzPlane( float L=350, float aspectRatio=1, float r0=0, float z0=0, float Emax=30 );
  void etaPhiPlane( float Deta=3.2, float Dphi=1.2, float eta0=0, float phi0=1, float Emax =30 );
  void drawEtaPhiDet();

  // event in the transverse momentum plane
  void ptPlane( float ptmax = 100., float Emax =30 );

  // B Field
  static float BField;
 
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

  // colors & widths
  static int kBeamAxisColor;
  static int kAxisColor;
  static int kLineWidth;
  static int kPipeFillColor;
  static int kPipeLineColor;
  static int kPixelFillColor;
  static int kPixelLineColor;
  static int kCtFillColor;
  static int kCt2FillColor;
  static int kCtLineColor;
  static int kEcalFillColor;
  static int kEcalLineColor;
  static int kHcalFillColor;
  static int kHcalLineColor;
  static int kCoilFillColor;
  static int kCoilLineColor;
  static int kYokeFillColor;
  static int kYokeLineColor;

  static int kPtFillColor;
  static int kPtTextColor;
  static int kPtLineColor;
  
  // helper
  void pzView( int proj=hYZ, 
	       float L=600, float aspectRatio=1, float r0=0, float z0=0, float Emax=30 );

  // primitives
  void drawHelix( int proj, float x0, float y0, float z0, float pt0, float eta0, float phi0, float rmax, float zmax, float B, 
		  int color, int linewidth, int linestyle );
  void drawBox( float xmin, float xmax, float ymin, float ymax, 
		int lineWidth, int lineColor, int fillColor,
		bool isInEvent=true  ); 
  void drawWedge( float dd, float rin, float rout, float phi0, 
		  int lineWidth, int lineColor, int fillColor,
		  bool isInEvent=true  ); 
  void drawPGon( float dd, float rin, float rout, float phi0, 
		 int lineWidth, int lineColor, int fillColor,
		 bool isInEvent=true  ); 
  void drawEBXtal( int ieta, int iphi, int color, float shift=0, float etamin=-1000, float etamax=1000, float phimin=-1000, float phimax=1000 );
  void drawEEXtal( int ix, int iy, int iz, int color, float shift=0, float etamin=-1000, float etamax=1000, float phimin=-1000, float phimax=1000 );
  void drawEcalRecHit( int ix, int iy, int iz, float e, float Emax=30, float emin=-1000, float etmin=-1000 );
  void drawCT( int proj, int ieta, int iphi, int lineWidth, int lineColor, int fillColor, bool isInEvent=true );  
  void drawCaloTower( float eta, float phi, float energy, float emEt, float hadEt, float Emax );  
  void drawCaloTower_( int ieta, int iphi, float energy, float emEt, float hadEt, float Emax );  

  void drawLine( int proj, float x0, float y0, float z0, float x, float y, float z, float R, int lineStyle, int lineWidth, int lineColor );
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

  void drawArrow( float x1, float y1, float x2, float y2, 
		  int lineColor=kRed, 
		  int lineWidth=4, 			
		  int lineStyle=kSolid,
		  float frac=0.3
		  );

  void drawMarkerRhoEtaPhi( int proj, float rho, float eta, float phi, float markerSize, int markerColor );
  void drawMarker( float x_, float y_, int style, int size, int icol1, int icol2=kWhite, int shift=0 );
  void drawMarkerXY( float x_, float y_, float size, int icol1, int icol2, int style_1, int style_2, bool doublePass=true );

  void drawMET( int proj, float et, float phi, float scale, float ptmax, int col=kBlue );
  void drawJet( int proj, float et, float eta, float phi, float ptmax, int icol );
  void drawCaloJet( int proj, float et, float eta, float phi, float area, float emf, float scale, float ptmax );

  // ecal hists
  TH2* _globalEcalHist;
  void refreshGlobalEcalHist( float Emax );
  void fillGlobalEcalHist( int ix, int iy, int iz, float val );
  void drawGlobalEcalHist( std::string opt );
  void drawElectronEcalHist( float pt, float eta, float phi,
			     std::vector<int>& v_ix, 
			     std::vector<int>& v_iy,
			     std::vector<int>& v_iz,
			     std::vector<float>& v_e,
			     int hist_size, std::string opt  );
  

  void setAtlas( bool atlas=true );

  // color scale

  // view dimensions
  float _Z, _R, _aspectRatio;
  float _r0, _x0, _y0, _z0, _eta0, _phi0;
  float _etamin, _etamax, _phimin, _phimax; 
  float _E, _P;
  void setEtaPhiLimits( float eta0, float phi0, float DE, float DPh );

  //  class TPolyLine;
  bool isInView( TPolyLine* pl, float xmin, float xmax, float ymin, float ymax );

  // run and event number
  void setCanvasName( const char* name=0 );
  TString _runAndEvt;

  // calo towers  FIXME !!!
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

  void registerTObject( TObject* );

  int curView;
  int curRun;
  int curEvent;

protected:

  // current histogram
  TH2* h_; 

  // garbage collection
  std::list< TObject* > _list;
  
ClassDef(EventDisplay,0) // EventDisplay -- Ecal event display

};

#endif

