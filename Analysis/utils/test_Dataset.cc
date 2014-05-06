#include <cassert>
#include <map>
#include <string>
#include <typeinfo>
#include <cmath>

#include <sstream>
#include <iostream>

#include <TFile.h>

#include "Config.hh"
#include "DatasetManager.hh"


using namespace std;

class Test {
public:
  static int value();
};

int Test::value()
{
  return 2;
}

int main()
{
  {
    cout << "Test reading" << endl;  
    TString file_ = Config::histPath;
    file_ += "ZZ_2l2n_cs.root";
    TFile* inputFile  = new TFile( file_, "READ" );
    inputFile->Print();

    DatasetManager dsm("dataset");
   
    cout << "set file " << endl;
    dsm.setFile( inputFile, DatasetManager::kRead );

    string dir_ = "data/q30_jet_met_bal_phj_phz_btg_dib";

    cout << dir_ << endl;
    RooDataSet* ds = dsm.dataset( dir_ );
    ds->Print();
    int numEntries = ds->numEntries();
    for( int jj=0; jj<numEntries; jj++ )
      {
	cout << "jj=" << jj 
	     << " Mll=" << ds->get(jj)->getRealValue("mll")
	     << endl;
      }
    return 0;
  }
  if(0)
    {
      {
	cout << "Test writing" << endl;  
	TFile* outputFile  = new TFile( "test_ds.root", "RECREATE" );
	DatasetManager dsm("METDataset");
   
	cout << "set file " << endl;
	dsm.setFile( outputFile, DatasetManager::kWrite );
	cout << "OK " << endl;
	float MET(0.1);
	float Mll(90.);
	vector<string> dir(2);
	dir[0] = "data";
	dir[1] = "Z_2e";
	for( size_t idir=0; idir<2; idir++ )
	  {
	    cout << "dir[" << idir << "]=" << dir[idir] << endl;
	    for( size_t ii=0; ii<(idir+1)*10; ii++ )
	      {
		MET += (idir+1)*ii*0.05;
		Mll += (idir+1)*ii*0.001;

		cout << "ii=" << ii << " MET=" << MET << " Mll=" << Mll << endl;

		dsm.add( "MET", MET, 0., 100. );
		dsm.add( "Mll", Mll, 0., 7000. );
	   
		dsm.flush( dir[idir] );
	      }
	  }
	dsm.save();
	outputFile->Close();
   
      }
      {
	cout << "Test reading" << endl;  
	TFile* inputFile  = new TFile( "test_ds.root", "READ" );
	DatasetManager dsm("METDataset");
   
	cout << "set file " << endl;
	dsm.setFile( inputFile, DatasetManager::kRead );
	vector<string> dir(2);
	dir[0] = "data";
	dir[1] = "Z_2e";

	for( size_t idir=0; idir<2; idir++ )
	  {
	    cout << "dir[" << idir << "]=" << dir[idir] << endl;
	    RooDataSet* ds = dsm.dataset( dir[idir] );
	    ds->Print();
	    int numEntries = ds->numEntries();
	    for( int jj=0; jj<numEntries; jj++ )
	      {
		cout << "jj=" << jj 
		     << " MET=" << ds->get(jj)->getRealValue("MET")
		     << " Mll=" << ds->get(jj)->getRealValue("Mll")
		     << endl;
	      }
	  }   
      }
    }
  return 0;
}
