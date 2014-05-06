#ifndef RooUtils_h
#define RooUtils_h

#include <TStyle.h>
#include <TH1.h>

class RooUtils
{
public:

  virtual ~RooUtils() {}
  
  // grid: Turns the grid lines on (true) or off (false)
  static void grid( TStyle*, bool gridOn );

  // fixOverlay: Redraws the axis
  static void fixOverlay();

  // set style
  static TStyle* setTdr();

  static void setTdrHistoStyle( TH1* h );

  // welcome banner
  static void welcome();

  // conversion functions
  static bool toXYZ(       float rho, float eta, float phi, float& x, float& y, float& z );
  static bool toRTheta(    float rho, float eta, float& r, float& theta );
  static bool toRhoEtaPhi( float x, float y, float z, float& rho, float& eta, float& phi );
  static int color( float e, float Emax );
  
  ClassDef( RooUtils, 0 )
};

#endif
