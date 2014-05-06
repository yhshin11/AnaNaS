//
// Dalitz Analysis
// 
           
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <unistd.h>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cassert>
#include <cmath>

#include "Analysis/fit/Model.hh"
#include "Analysis/core/Sample.hh"

using namespace std;

int main( int argc, char *argv[] )
{
  cout << "Fitting program" 
       << endl;
  cout << "Author : Gautier Hamel de Monchenault " 
       << endl;
  cout << "-- July 2009" 
       << endl;
  //
  // Read the options from the command line
  //
  TString analysis("Diboson");
  TString sample("Z_2e");
  TString subsample("all");
  TString hist("stat__mll__0_all_nosel");
  TString variable("mll");

  int c;
  while ( (c = getopt(argc, argv, "S:s:a:v:h:")) != EOF ) {
    switch (c) {
    case 'a': analysis  = TString(optarg); break;
    case 's': sample    = TString(optarg); break;
    case 'S': subsample = TString(optarg); break;
    case 'h': hist      = TString(optarg); break;
    case 'v': variable  = TString(optarg); break;
    }
  }

  Model model( analysis, sample, subsample, hist, variable );
  cout << "fit=" << model.fit() << endl;

  // be polite, say good-bye
  cout << "\nEnd of program -- take care" 
       << endl;

  return 0;

}

