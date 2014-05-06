#include <iostream>
#include <errno.h>
#include <dirent.h>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cassert>
#include <string>
#include <map>

#include <dirent.h>
#include <sys/types.h>
using namespace std;

#include <TString.h>
#include <TFile.h>
#include <TH1.h>

#include "Analysis/tools/PrimaryAnalysis.hh"
#include "Analysis/core/EventServer.hh"

int main(int argc, char **argv)
{
  
  PrimaryAnalysis TreeAnalysis("/home/usr201/mnt/jmalcles/Bosons/AnaNaS/forGautier/configPA.txt");

  //TreeAnalysis.AddAllWeights();
  //TreeAnalysis.ModifyTrees();
  //TreeAnalysis.AllEff();
  //TreeAnalysis.AllWTauEff();
  //TreeAnalysis.AllZEEEff();
  //TreeAnalysis.fitAll();
  //TreeAnalysis.fitAllFlipCut();

  TreeAnalysis.AddWeight("FullLumi");

  //TreeAnalysis.fitAllMETPara("W_en");

  return 0;
}
