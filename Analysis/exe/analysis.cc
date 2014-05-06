#include <iostream>
#include <string>
#include <vector>
#include <map>

#include <TTimeStamp.h>

#include <string>
using namespace std;

#include "Analysis/utils/Config.hh"
#include "Analysis/core/Sample.hh"
#include "Analysis/src/AnalysisFactory.hh"

int main(int argc, char **argv)
{  
  string sample("ZZ");
  string analysis("Diboson");
  string sampleFile("samples.txt");
  string configFile("collections.txt");
  string eventFile("");
  string outpath("");
  string version("");

  string filename("");

  int n_(-1);  // number of events
  int k_(-1);   // number of events to skip

  int c;
  while ( (c = getopt( argc, argv, "a:n:s:f:k:c:F:o:v:h" ) ) != EOF ) 
    {
      switch (c) 
	{
	case 'a': analysis = string(optarg); break;
	case 'n': n_     = atoi(optarg); break; 
	case 'k': k_     = atoi(optarg); break;
	case 's': sample = string(optarg); break;
	case 'f': filename = string(optarg); break;
	case 'c': configFile = string(optarg); break;
	case 'F': eventFile = string(optarg); break;
	case 'o': outpath = string(optarg); break;
	case 'v': version = string(optarg); break;
	case 'h':
	  {
	    cout << "--- usage ---" << endl;
	    cout << "workdir\% analysis -a<Analysis=\"Diboson\"> -s<Sample=\"ZZ\"> -c <ConfigFile=\"collections.txt\"> -n<NumberOfEvents=-1>" << endl;
	    cout << "--- to check samples: --- " << endl;
	    cout << "workdir\% checkSamples" << endl;
	    return 0;
	  }
	}
    }
  
  // start
  TTimeStamp date;
  date.Set();
  cout << "Start at date " << date.AsString() << endl;

  string samplePath  = Config::confPath + sampleFile;
  string configPath  = Config::confPath + configFile;

  // additionnal outpath, useful for gpfs redirection
  if(outpath!="")
    {Config::rootPath = outpath;}

  //additionnal version, useful when samples are not back-compatible
   if(version!="")
    {Config::version = version;}


  // sample
  Sample* s_(0);
  if( !eventFile.empty() ) 
    {
      s_ = new Sample( eventFile.c_str() );
      s_->setName( eventFile.c_str() );
    }
  else if( filename!="" )
    {
      s_ = new Sample( filename.c_str() );
    }
  else
    {
      s_ = Sample::get( sample.c_str(), samplePath.c_str() );
    }
  if( s_==0 )
    {
      cout << "Sample \"" << sample << "\" does not exist!" << endl;
      return 1;
    }

  EventManager* e_ = new EventManager( s_->filenames, configPath.c_str() );

  // instanciate analysis object
  SampleAnalysis* a_ = AnalysisFactory::get( analysis, *s_, *e_ );
  if( a_==0 )
    {
      cout << "Analysis \"" << analysis << "\" does not exist!" << endl;
      return 2;
    }
  
  // analyze n_ events
  a_->analyze( n_, k_ );

  // end
  date.Set();
  cout << "End at date   " << date.AsString() << endl;

  return(0);
}
