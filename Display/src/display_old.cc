#include <iostream>
#include <cstdlib>

#include "TRint.h"
#include "TGClient.h"

using namespace std;
#include "Display/src/EventDisplay.hh"

int main(int argc, char **argv)
{
  TString sample("Zee");
  int i(1);
  
  int c;
  while ( (c = getopt( argc, argv, "s:i:" ) ) != EOF ) 
    {
      switch (c) 
	{
	case 's': sample  = optarg;    break;
	case 'i': i       = atoi(optarg);    break;
	}
    }

  TString filename = TString("data/") + sample + TString("/Ntuple_");
  filename += i;
  filename += TString(".root");
  cout << "Display for file=" << filename << endl;

  TRint *theApp = new TRint("App", &argc, argv, NULL, 0);
  TString command("EventDisplay d(\""); command += filename; command += "\")";
  theApp->ProcessLine(command);
  theApp->ProcessLine("d.nextEvent()");
  theApp->ProcessLine("d.primer()");
  theApp->Run();

  return 0;
}
