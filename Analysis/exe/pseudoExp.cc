#include <iostream>
#include <string>
#include <vector>
#include <map>

#include <TTimeStamp.h>
#include <TString.h>

using namespace std;

#include "Analysis/utils/Config.hh"
#include "Analysis/utils/Constants.hh"
#include "Analysis/core/Sample.hh"
#include "Analysis/src/PseudoExperiments.hh"

int main(int argc, char **argv)
{  
  TString sampleFile = Config::confPath + "samples.txt";

  int n_(-1);
  float lumi_(100);
  TString sample_("Zee");

  int c;
  while ( (c = getopt( argc, argv, "n:l:s:" ) ) != EOF ) 
    {
      switch (c) 
	{
	case 's': sample_ = TString(optarg); break;
	case 'n': n_      = atoi(optarg); break;
	case 'l': lumi_   = atof(optarg); break;
	case 'h':
	  {
	    cout << "--- usage ---" << endl;
	    cout << "workdir\% pseudoExp -l<Lumi=\"" << lumi_ << " nb^-1\">" << endl;
	    return 0;
	  }
	}
    }
  
  Sample* s_ = Sample::get( sample_, sampleFile );
  if( s_==0 ) 
    {
      cout << "Sample not found " << endl;
      return -1;
    }

  // start
  TTimeStamp date;
  date.Set();
  cout << "Start at date " << date.AsString() << endl;

  // instanciate pseudo-experiment analysis object
  PseudoExperiments pe_( *s_, lumi_ );

  // analyze n_ events
  pe_.analyze( n_ );

  // end
  date.Set();
  cout << "End at date   " << date.AsString() << endl;

  return(0);
}
