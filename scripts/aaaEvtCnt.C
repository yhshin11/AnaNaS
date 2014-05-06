// ROOT macro.
// Gets number of events contains in a Ntuple_*.root AnaNaS file.
// The number of events is obtained from the number of entries of tree
// 'eventSummary'
// author: Ph Gras CEA/IRFU Saclay 2010-04-28

#include <TFile.h>
#include <TTree.h>
#include <iostream>

using namespace std;

void aaaEvtCnt(const char* filename){
  TFile* f = new TFile(filename);
  if(f->IsZombie()){
    cerr << "Failed to open " << filename << endl;
    return;
  }

  TTree* t;
  const char* sumTreeName = "eventSummary";
  f->GetObject(sumTreeName, t);
  if(t == 0){
    cerr << "ROOT tree sumTreeName" << "was not found." <<endl;
    return;
  }

  int nevts = t->GetEntries();
  cout << nevts << endl;
}
