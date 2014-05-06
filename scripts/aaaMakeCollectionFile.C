/** ROOT macro to produce a collections.txt template file for ananas
 * from an AnaNaS ntuple file.
 * Note the letter case of the section might be incorrected.
 * Author: Ph. Gras CEA/Saclay 2010-04-28
 */

#include <iostream>
#include <fstream>
#include <ctype.h>
#include <map>
#include <vector>
#include "TFile.h"
#include "TKey.h"
#include "TString.h"
#include "TROOT.h"
#include "TPRegexp.h"
#include "TClass.h"

using namespace std;

void aaaMakeCollectionFile(const char* filename = "Ntuple_1.root",
                           const char* outFileName = "_collections.txt"){
  TFile* file =  new TFile(filename);
  if(file->IsZombie()){
    cerr << "Failed to find or open file " << filename << endl;
    return;
  }

  TIter next(file->GetListOfKeys());
  TKey* key;
  std::map<TString, std::vector<TString> > m;
  while ((key = (TKey*)next())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TTree")) continue;
    const char* tname = key->GetName();

    TStringToken tok(tname, "_");
    TString pref = "";
    TString cand = "";
    TString algo = "";
    if(tok.NextToken()) pref = tok;
    if(tok.NextToken()) cand = tok;
    if(tok.NextToken()) algo = tok;

    if(pref.CompareTo("ntu")!=0) continue;
    if(cand.Length()==0){
      cerr << tok.Data() << ": no candidate name?" << endl;
      continue;
    }
    
    if(algo.Length()==0){
      cerr << tok.Data() << ": no algo name?" << endl;
      continue;
    }

    //capitalize candidate name:
    cand[0] = toupper(cand[0]);

    m[cand].push_back(algo);
  }

  ofstream f(outFileName);

  for(std::map<TString, std::vector<TString> >::const_iterator itCand = m.begin();
      itCand != m.end();
      ++itCand){
    if(itCand != m.begin()) f << "\n";
    f << "[" << itCand->first.Data() << "]\n";
    const vector<TString>& algos = itCand->second;
    for(vector<TString>::const_iterator itAlgo = algos.begin();
        itAlgo != algos.end();
        ++itAlgo){
      f << "\t" << itAlgo->Data() <<"\n";
    }
  }
}
