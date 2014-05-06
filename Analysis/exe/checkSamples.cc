#include <iostream>
#include <string>
#include <vector>
#include <map>

#include <TString.h>

using namespace std;

#include "Analysis/core/Sample.hh"
#include "Analysis/core/SampleAnalysis.hh"

#include <iostream>
#include <string>
#include <vector>
using namespace std;
#include <TFile.h>
#include <TTree.h>

#include "Analysis/utils/Config.hh"

int main(int argc, char **argv)
{  
  Sample::verbose = true;
  Sample::check   = true;
  TString sampleFile = Config::confPath + "samples_test.txt";
  Sample::set( sampleFile );
  
  return(0);
}
