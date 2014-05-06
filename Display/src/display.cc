#include <iostream>
#include <cstdlib>

#include "TRint.h"
//#include "TGClient.h"

using namespace std;
#include "Display/src/EventDisplay.hh"

#include "Display/src/DisplayManager.hh"

#include "Analysis/utils/Config.hh"

int main(int argc, char **argv)
{
  TString sample("Z_2e");
  int i(0);  // file number
  int n(1);  // event number in file
  bool minBias2009 = false;

  int r(0);  // first run number
  int R(0);  // last run number
  int e(-1); // first event number
  int E(-1); // last event number
  
  TString filename;
  TString configFile("collections.txt");

  bool styleAtlas(false);

  int c;
  while ( (c = getopt( argc, argv, "s:f:i:n:r:R:e:c:MA" ) ) != EOF ) 
    {
      switch (c) 
  	{
  	case 's': sample  = optarg;    break;
  	case 'f': filename  = optarg;    break;
  	case 'c': configFile  = optarg;    break;
  	case 'i': i       = atoi(optarg);    break;
  	case 'n': n       = atoi(optarg);    break;
  	case 'r': r       = atoi(optarg);    break;
  	case 'R': R       = atoi(optarg);    break;
  	case 'e': e       = atoi(optarg);    break;
  	case 'M': minBias2009 = true;    break;
  	case 'A': styleAtlas = true;    break;
  	}
    }

  if( r>0 )
    {
      DisplayManager::firstRun = r;
      DisplayManager::lastRun  = r;
      if( R>r ) DisplayManager::lastRun  = R;
      if( e>0 )
  	{
  	  DisplayManager::firstEvent = e;
  	  DisplayManager::lastEvent  = e;
  	  if( E>e ) DisplayManager::lastEvent  = E;
  	}	
    }

  if( filename=="" )
    {
      filename = Config::workPath + TString("data/") + sample;
      if( i!=0 )
  	{
  	  filename += TString("/Ntuple_");      
  	  filename += i;
  	  filename += TString(".root");
  	}
    }
  cout << "Display for file=" << filename << endl;

  //  TRint *theApp = new TRint("App", &argc, argv, NULL, 0);
  TRint *theApp = new TRint("App", 0, 0, NULL, 0);

  TString command("DisplayManager m(\""); 
  command += filename;
  command += "\", \"";
  command += Config::confPath;
  //  command +=  "collections.txt\"";
  command +=  configFile;
  command +=  "\"";
  command += "); ";
  cout << command << endl;
  theApp->ProcessLine(command);  

  //  if( n!=0 )
  //    {
  //      command = "m.goToEvent("; command += n; command += ");";
  //    }
  //  else
    {
      command = "m.nextEvent();"; 
    }
  cout << command << endl;
  theApp->ProcessLine(command);
  
  //  theApp->ProcessLine(command);

  if( styleAtlas )
    {
      command = "m.setAtlas(); ";
      theApp->ProcessLine(command);
    }

  theApp->Run();

  return 0;
}
