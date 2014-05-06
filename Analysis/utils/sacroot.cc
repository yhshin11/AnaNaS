#include <iostream>
using namespace std;

#include <TString.h>
#include <TRint.h>

int main(int argc, char **argv)
{
  TString command;
  TRint *theApp = new TRint("App", &argc, argv, NULL, 0);
  command="RooUtils::welcome();";
  theApp->ProcessLine(command);
  command="RooUtils::setTdr();";
  theApp->ProcessLine(command);
  command="gROOT->SetStyle(\"tdrStyle\");";
  theApp->ProcessLine(command);
  theApp->Run();

  return 0;
}
