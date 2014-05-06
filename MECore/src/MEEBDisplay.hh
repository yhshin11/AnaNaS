#ifndef MEEBDisplay_hh
#define MEEBDisplay_hh

//
// Authors  : Gautier Hamel de Monchenault and Julie Malcles, Saclay 
//

#include <list>
#include <map>

#include "MEEBGeom.hh"

class MEEBDisplay
{
  // static functions
public:

  static MEEBGeom::EtaPhiPoint getNode( int iSM, 
					MEEBGeom::EBTTLocalCoord iX, 
					MEEBGeom::EBTTLocalCoord iY, 
					int jx, int jy );
  static void drawEB(  float etamin=-1000, float etamax=1000,
		       float phimin=-1000, float phimax=1000
		       );
  static void drawSM(  int iSM, float shift=0 );
  static void drawTT(  int iSM, 
		       MEEBGeom::EBTTLocalCoord iX, 
		       MEEBGeom::EBTTLocalCoord iY, float shift=0,
		       float etamin=-1000, float etamax=1000,
		       float phimin=-1000, float phimax=1000
		       );
  static void drawXtal( MEEBGeom::EBGlobalCoord ieta, 
			MEEBGeom::EBGlobalCoord iphi, 
			int color=kBlue, float shift=0,
			float etamin=-1000, float etamax=1000,
			float phimin=-1000, float phimax=1000
			);
  static MEEBGeom::EtaPhiPoint center( MEEBGeom::EBGlobalCoord ieta, 
				       MEEBGeom::EBGlobalCoord iphi  ); 

  static void drawEBGlobal();
  static void drawEBLocal();

  static TPolyLine* getXtalPolyLine( MEEBGeom::EBGlobalCoord ieta, 
				     MEEBGeom::EBGlobalCoord iphi, 
				     float shift=0,
				     float etamin=-1000, float etamax=1000,
				     float phimin=-1000, float phimax=1000
				     );    
  static TPolyLine* getTTPolyLine( int iSM, 
				   MEEBGeom::EBTTLocalCoord iX, 
				   MEEBGeom::EBTTLocalCoord iY, 
				   float shift=0,
				   float etamin=-1000, float etamax=1000,
				   float phimin=-1000, float phimax=1000
); 
  static TPolyLine* getSMPolyLine( int iSM, float shift=0 );
  static void drawRz( int lc=kGreen+3, int fc=kGreen-9 );
  // GHM feb12
  static void drawXy( int lc=kGreen+3, int fc=kGreen-9 );
  // GHM feb12
  static int bkgColor;
  static int lineColor;
  static int lineWidth;

  static void refresh();

  virtual ~MEEBDisplay() {}

private:

  static void setPhiLimits( int iSM,  
			    MEEBGeom::EBLocalCoord iy, 
			    MEEBGeom::EBGlobalCoord iphi, 
			    float phioverpi_0, float phioverpi_1 );
  static void setEtaLimits( int iSM,  
			    MEEBGeom::EBLocalCoord ix, 
			    MEEBGeom::EBGlobalCoord ieta, 
			    float eta_0, float eta_1 );
  static void setSM_2_and_20();

  static std::map< int, std::pair<float,float>  > _phiLimits;
  static std::map< int, std::pair<float,float>  > _etaLimits;

  static void setRzXtals();
  static std::map< int, TPolyLine* > _rzXtals;
  static std::map< int, TPolyLine* > _xyXtals;

  static std::list<TObject*> _list;
  static void registerTObject( TObject* );
  
  ClassDef(MEEBDisplay,0) // MEEBDisplay -- Monitoring utility for survey of Ecal
};

#endif

