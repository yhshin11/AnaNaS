#ifndef DisplayManager_h
#define DisplayManager_h

#include <list>

#include "Analysis/core/EventManager.hh"
#include "Display/src/MultiDisplay.hh"

class DisplayManager : public EventManager
{
public:

  // options & parameters
  static bool  showEcalRH;
  static bool  showEcalSC;
  static bool  showCaloMET;
  static bool  showCaloJets;
  static bool  showCaloTowers;
  static bool  showPfCands;
  static bool  showPfJets;
  static bool  showPfLeptons;
  static bool  showPfMet;
  static bool  showVertices;
  static bool  showTracks;
  static bool  showTrackHits;
  static bool  showLeptonHits;
  static bool  showElectrons;
  static bool  showMuons;
  static bool  showPhotons;
  static bool  showMcTruth;
  static bool  showTrkMET;
  static bool  showOnlyPVTracks;
  static bool  showOnlyPVTracks_vertex;

  static bool  showElectronIso;
  static bool  showMuonIso;
  static bool  showPhotonIso;
  static float electronDRin; 
  static float electronDRout;
  static float muonDRin;     
  static float muonDRout;    
  static float photonDRin;   
  static float photonDRout;  

  static float PtMax;
  static float TrackPtMin;
  static float TrackRMin;
  static float TrackRMax;
  static float TrackZMax;
  static float LeptonPtMin;
  static float CaloJetEtMin;
  static float PfJetEtMin;
//  static float McPtcut;
//  static float McEtacut;
//  static float Wghcut;
//  static float JetEtcut;
  static float Emax;
  static float length;
  static float delta;
  static float BField;
  static float r0;
  static float x0;
  static float y0;
  static float z0;
  static float R;
  static float L;
  static float AR;
  static float DE;
  static float DPh;
  static float eta0;
  static float phi0;
  static float CTScale;

  static int lowPtTrackColor;
  static int trackColor;
  //  static int hitColor;
  static int CTHadColor;
  static int CTEmColor;

  static int  col_electron      ;
  static int  col_muon          ;
  static int  col_pvTrack_clu   ;
  static int  col_pvTrack_unclu ;
  static int  col_puTrack_clu   ;
  static int  col_puTrack_unclu ;

  static int firstRun;
  static int lastRun;
  static int firstEvent;
  static int lastEvent;

  DisplayManager( const char* filename, const char * collectionFileName );
  DisplayManager( std::vector< std::string > filenames, const char * collectionFileName ) ;
  virtual ~DisplayManager( );

  virtual bool isSelected();

  virtual bool nextEvent();
  virtual bool goToEvent( int evt );

//   void setView( int v ) { _view=v; }
//   void XYView()     { setView(hXY); view(); }
//   void YZView()     { setView(hYZ); view(); } 
//   void XZView()     { setView(hXZ); view(); }
//   void RZView()     { setView(hRZ); view(); }
//   void PhiZView()   { setView(hPhiZ); view(); }
//   void PtView()     { setView(hPt); view(); }
//   void EtaPhiView() { setView(hEtaPhi); view(); }

  bool view();
  void vertexViews();
  void projectionViews();
  void ptView();
  void etaPhiView();

  void zoomXYZ( float zoom=1, float x0=0, float y0=0, float z0=0 );
  void zoomVertex( float zoom=1 );
  void energyScales( float pt_scale=1, float E_scale=1, float CT_scale=1 );

  //  void setPhi0( float );
  //  void setPhi0( const Candidate* );
  
  //private:

  void displayElectrons( float zoom=1., float scaleE=1. );
  bool displayElectron( size_t ii, float zoom=1., float scaleE=1. );
  bool displayLepton( int ilepton, size_t ii, float zoom=1., float scaleE=1. );

  //  bool drawLepton( int ilepton, size_t ii, float scale=1., float scaleE=1. );

  //  void drawNextElectron(  float scale=1., float scaleE=1. );
  //  bool redrawElectron( float scale=1., float scaleE=1. );
  //  bool drawElectron( size_t ii, float scale=1., float scaleE=1. );

  //  void drawNextMuon(  float scale=1., float scaleE=1. );
  //  bool redrawMuon( float scale=1., float scaleE=1. );
  //  bool drawMuon( size_t ii, float scale=1., float scaleE=1. );

  //  void drawNextLepton(  float scale=1., float scaleE=1. );

  //  void drawNextMuon(      float scale=1., float scaleE=1. );

  void drawTrack( const Candidate&, int trajCol=1, int trajWidth=1,
		  int trajStyle=1 );
  void drawCaloTower( const Candidate&  );
  void drawVertex( const Vertex& vtx, int markerColor, bool isPrimary );
  void drawPrimaryVertex();
  void drawVertices();
  //  void drawTracks( int color=kBlack, int width=1, int style=1 );
  void drawTracks();
  void drawElectrons();
  void drawMuons();
  //  void drawTrackHits( Candidate* trk, float size, int color );
  void drawPfCandidates();
  //  void drawPfJets();
  void drawCaloTowers();
  //  void drawCaloJets();
  //  void drawCaloMet();
  //  void drawTrkMet();
  void drawEcalRecHits( float Emax=30., bool applyThresh=true );

  void drawJet( size_t ii, int color=kYellow );
  //  void drawEcalSuperClusters();

  //  void drawGlobalEcalHist( std::string opt="COLZ" );

  //  void drawElectronEcalHist( int hist_size=50, std::string opt="COLZ" );

  //  void drawElectronHits( int ii );

  void draw2DVectors( float scale, float dphi, 
 		      int lineWidth, int fillStyle, 
 		      int color1, int color2=kWhite, 
 		      bool fillSameColor=false );

  //  void setAtlas( float ptmin=0.5 );

  //  void drawTexts();

  MultiDisplay* _ed;
  // public forward functions
  bool goToPanel( TString str_ )
  {
    return _ed->goToPanel( str_ );
  }

  TPad* goToPad( TString str_ )
  {
    TPad* fPad_(0);
    if( _ed->goToPanel(str_) ) fPad_ = _ed->fPad;
    return fPad_;
  }

  //  int _view;

  //  int _currentElectron;
  //  int _currentMuon;

  struct ThresAndColor
  {
    float pt;
    int color;
    int width;
    int style;
  };
  struct SortThresAndColor
  {
    bool operator()( const ThresAndColor& x1, const ThresAndColor& x2 ) const
    {
      float th1 = x1.pt;
      float th2 = x2.pt;
      return fabs( th1 ) < fabs( th2 );
    }    
  };
  std::vector<ThresAndColor> _thres;
  void setThreshold( float val, int color, int width=1, int style=1 );
  void setPtMin( float ptmin );
  void setupDisplay();
    
  // class definition
  ClassDef( DisplayManager, 0 )    
     
};

#endif
