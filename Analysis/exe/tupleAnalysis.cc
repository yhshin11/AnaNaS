#include <unistd.h> 
#include <iostream>
#include <sstream>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <typeinfo>
#include <string>
#include <cmath>
#include <algorithm>
#include <functional>
#include <vector>
#include <map>

#include <TFile.h>
#include <TBranch.h>
#include <TKey.h>
#include <TTree.h>
#include <TString.h>
#include <TBits.h>
#include <Math/VectorUtil.h>
#include <TMath.h>
#include <TVector3.h>
#include <TLorentzVector.h>

#include "Analysis/core/Sample.hh"
#include "Analysis/utils/Constants.hh"
#include "Analysis/utils/Config.hh"
#include "Analysis/utils/RooUtils.hh"

#include "Analysis/ghmzz/GhmTupleAnalysis.hh"

int main(int argc, char **argv)
{
  //
  // Histograms
  //
  RooUtils::setTdr();

  //
  // Samples
  //
  Sample::verbose = true;
  Sample::check   = false;

  string sampleSetName("data2011_DoubleElectron");  
  int n_max = -1;

  int c;
  while ( (c = getopt( argc, argv, "s:n:h" ) ) != EOF ) 
    {
      switch (c) 
	{
	case 's': sampleSetName = string(optarg); break;
	case 'n': n_max = atoi(optarg); break;
	case 'h':
	  {
	    cout << "tupleAnalysis -s [sampleSet]" << endl;
	    return 0;
	  }
	}
    }

  GhmTupleAnalysis tupleAnalysis( sampleSetName.c_str(), n_max );
  tupleAnalysis.go();

  return 0;
}

