/*********************************************************************
 * 
 *********************************************************************/

#ifndef PrimaryAnalysis_hh
#define PrimaryAnalysis_hh

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <map>
#include <TCut.h>
#include <TString.h>


using namespace std;

class PrimaryAnalysis
{
public:


  PrimaryAnalysis( const char* configfile="/home/usr201/mnt/jmalcles/Bosons/AnaNaS/forGautier/configPA.txt" );

  virtual ~PrimaryAnalysis();



public:

  pair<float, float> VBTFSelection( int percent, const string & var, const string & ecalPart);
  void AddWeight(const string & filename);
  void AddWeights(const string & specy);
  void AddAllWeights();
  void ModifyTree(const string & filename); 

  double errEff(double v1, double vtot);
  double errRat(double v1, double v2, double s1, double s2);
  double errProd(double v1, double v2, double v3, double v4,double s1, double s2,double s3, double s4);


  void ModifyTrees();

  TCut GetCut( int percent, bool doIso, bool doID, bool doEmaAndNoInit, bool doFlipIso=false, bool doFlipID=false, bool isRoo=false );
  string CutName( int percent, bool doIso, bool doID, bool doEmaAndNoInit , bool doFlipIso=false, bool doFlipID=false);
  std::vector<double> efficiency(const string & specy, TCut cut, TCut cutRef, double METParaMin, double METParaMax);

  
  void AllEfficiencies(const string & specy);
  void AllEff();

  std::vector<double> WTauEEff(int percent=95, bool doID=true, bool doIso=true, bool doEmaAndNoInit=true, double METParaMin=-50., double METParaMax=80., bool doFlipID=false, bool doFlipIso=false );
  void AllWTauEff();
  std::vector<double> ZEEEEff(int percent=95, bool doID=true, bool doIso=true, bool doEmaAndNoInit=true, double METParaMin=-50., double METParaMax=80., bool doFlipIso=false , bool doFlipID=false );


  void AllZEEEff();
  void fitMETPara(const string & specy, TCut cut, const string & cutName);

  void fitAllMETPara(const string & specy);


  void fitAllFlipCut();
  void fitAll();


private:

  void init();
  void configure();
  void ReadWeights();
  void ReadVBTF();
  void ReadFullLumiFile();

  TString _configfile;
  string _datapath;
  string _newdatapath;
  string _outputpath;
  string _VBTFFile;
  string _weightsFile;
  string _ConfigFullLumi;

  float _IsoCuts[2][4][3];
  float _IDCuts[2][4][4];
  
  map <string, double> _VBTF;
  map <string, double> _Weights;

  vector<string> LumiFiles;

  ClassDef(PrimaryAnalysis,0)

};


#endif
