#include <cassert>
#include <map>
#include <string>
#include <typeinfo>

#include <iostream>

#include <TFile.h>
#include <TBranch.h>
#include <TKey.h>
#include <TTree.h>

#include "TreeManager.hh"

#include <iostream>
#include <TBits.h>

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
  // simple variable
  if ( true )
    {
      std::string ofilename("test.root");
      std::cout << "Will write ntuple in output file " << ofilename << std::endl;
      TFile* outputFile  = new TFile( ofilename.c_str(), "RECREATE" );

      std::cout << "--------------- variable I/O\n";
      //      TreeManager tm( "test.root", TreeManager::kWrite );
      TreeManager tm( outputFile, TreeManager::kWrite );
      tm.addTree( "test_tree" );

      int wval = 3;
      tm.branchio<int>( "test_tree", "pluto", (Test::value()) );
      float wfval = 3;
      tm.branchio<float>( "test_tree", "papero", wfval );
      tm.sync( "test_tree" );

      wval = 5;
      tm.branchio<int>( "test_tree", "pluto", wval );
      wfval = 500;
      tm.branchio<float>( "test_tree", "papero", wfval );
      tm.sync( "test_tree" );
      tm.scan();
      tm.save();

      TreeManager tm_read;
      tm_read.open( "test.root", TreeManager::kRead );

      int rval = -1;
      tm_read.branchio<int>( "test_tree", "pluto", &rval );
      float rfval = -1;
      tm_read.branchio<float>( "test_tree", "papero", &rfval );
      {
	int nentries = tm_read.t( "test_tree"  )->GetEntries();
	for ( int ientry = 0; ientry < nentries; ++ientry ) {
	  tm_read.entry( ientry, "test_tree" );
	  std::cout << "pluto: " << rval << "\n";
	  std::cout << "papero: " << rfval << "\n";
	}
      }

      outputFile->Close();

      std::cout << "--------------- variable I/O finished\n";
    }
  // std::vector
  if ( false )
    {
      std::cout << "--------------- std::vector<> I/O\n";
      TreeManager tm( "test.root", TreeManager::kWrite );
      tm.addTree( "test_tree" );
      std::vector<int> * v = new std::vector<int>;
      v->push_back(0);
      v->push_back(1);
      v->push_back(2);
      v->push_back(3);
      v->push_back(50);
      v->push_back(27);
      v->push_back(81);
      tm.branchio< std::vector<int> * >( "test_tree", "vect", &v );
      tm.sync( "test_tree" );
      v->clear();
      v->push_back(50);
      v->push_back(27);
      v->push_back(81);
      tm.branchio< std::vector<int> * >( "test_tree", "vect", &v );
      tm.sync( "test_tree" );
      v->clear();
      for ( size_t i = 0; i < 100; ++i ) {
	v->push_back(i*i);
      }
      tm.branchio< std::vector<int> * >( "test_tree", "vect", &v );
      tm.sync( "test_tree" );
      tm.save(true);
      TreeManager tm_read;
      tm_read.open( "test.root", TreeManager::kRead );
      std::vector<int> *w = 0;
      std::cout << "**> address: " << &w << "\n";
      tm_read.branchio< std::vector<int> * >( "test_tree", "vect", &w );
      {
	int nentries = tm_read.t( "test_tree"  )->GetEntries();
	std::cout << nentries << " entries found here.\n";
	for ( int ientry = 0; ientry < nentries; ++ientry ) {
	  tm_read.entry( ientry, "test_tree" );
	  std::cout << w << " size: " << w->size() << "\n";;
	  for ( size_t i=0; i<w->size(); ++i) {
	    std::cout << "w " << i << " " << (*w)[i] << "\n";
	  }
	}
      }

      std::cout << "--------------- std::vector<> I/O finished\n";
    }
  // std::string
  if ( false )
    {
      std::cout << "--------------- std::string<> I/O\n";
      TreeManager tm( "test.root", TreeManager::kWrite );
      tm.addTree( "test_tree" );
      std::string s = "ciao!!!!!!!!!!!!!!!!";
      tm.branchio< std::string >( "test_tree", "s", s );
      tm.sync( "test_tree" );
      tm.save();
      TreeManager tm_read;
      tm_read.open( "test.root", TreeManager::kRead );
      std::string * rs = 0;
      std::cout << "**> address: " << &rs << "\n";
      tm_read.branchio< std::string * >( "test_tree", "s", &rs );
      {
	int nentries = tm_read.t( "test_tree"  )->GetEntries();
	std::cout << nentries << " entries found here.\n";
	for ( int ientry = 0; ientry < nentries; ++ientry ) {
	  tm_read.entry( ientry, "test_tree" );
	  std::cout << "the string was: " << (*rs) << "\n";
	}
      }
      std::cout << "--------------- std::string<> I/O finished\n";
    }
  // object
  if ( false )

    {
      std::cout << "--------------- object I/O\n";
      TreeManager tm( "test.root", TreeManager::kWrite );
      tm.addTree( "test_bits" );
      TBits *wbits = new TBits(256);
      wbits->SetBitNumber(2, kTRUE);
      wbits->SetBitNumber(3, kTRUE);
      wbits->SetBitNumber(4, kTRUE);
      wbits->SetBitNumber(15, kTRUE);
      wbits->SetBitNumber(16, kTRUE);
      wbits->SetBitNumber(17, kTRUE);
      tm.branchio<TBits*>( "test_bits", "bits", &wbits );
      std::cout << "bits: " << (*wbits) << "\n";
      tm.sync( "test_bits" );
      tm.save();

      TreeManager tm_read;
      tm_read.open( "test.root", TreeManager::kRead );
      //TBits *rbits = new TBits(256); // one or the other are fine
      TBits *rbits = 0;
      tm_read.branchio<TBits*>( "test_bits", "bits", &rbits );
      {
	int nentries = tm_read.t( "test_bits" )->GetEntries();
	for ( int ientry = 0; ientry < nentries; ++ientry ) {
	  tm_read.entry( ientry, "test_bits" );
	  std::cout << "bits: " << (*rbits) << "\n";
	}
      }
      std::cout << "--------------- object I/O finished\n";
    }

  TreeManager tm_read;
  tm_read.open( "data/Z_2e/Ntuple_101_1_Zdq.root", TreeManager::kRead );
  
  int rval = -1;
  tm_read.branchio<int>( "ntu_tracks", "event", &rval );

  float rfval = -1;
  tm_read.branchio<float>( "ntu_tracks", "pt", &rfval );

  //TBits *rbits = new TBits(256); // one or the other are fine
  TBits *rbits = 0;
  tm_read.branchio<TBits*>( "ntu_triggerResults", "hltBits", &rbits );

  std::string * rs = 0;
  std::cout << "**> address: " << &rs << "\n";
  tm_read.branchio< std::string * >( "runSummary", "electron", &rs );

  std::vector<float> *w = 0;
  std::cout << "**> address: " << &w << "\n";
  tm_read.branchio< std::vector<float> * >( "ntu_electrons", "E_ECALdep", &w );

  if( false )
  {
    int nentries = tm_read.t( "ntu_tracks"  )->GetEntries();
    for ( int ientry = 0; ientry < nentries; ++ientry ) 
      {
	tm_read.entry( ientry, "ntu_tracks" );
	std::cout << "event: " << rval << "\t";
	std::cout << "pt: " << rfval << "\n";
      }
  }

  if( false )
  {
    int nentries = tm_read.t( "ntu_electrons"  )->GetEntries();
    std::cout << nentries << " entries found here.\n";
    for ( int ientry = 0; ientry < nentries; ++ientry ) 
      {
	tm_read.entry( ientry, "ntu_electrons" );
	std::cout << w << " size: " << w->size() << "\n";;
	for ( size_t i=0; i<w->size(); ++i) 
	  {
	    std::cout << "w " << i << " " << (*w)[i] << "\n";
	  }
      }
  }

  if( false )
  {
    int nentries = tm_read.t( "ntu_triggerResults" )->GetEntries();
    for ( int ientry = 0; ientry < nentries; ++ientry ) 
      {
	tm_read.entry( ientry, "ntu_triggerResults" );
	std::cout << ientry << "/";
	std::cout << "bits: " << (*rbits) << "\n";
      }
  }

  if( false )
  {
    int nentries = tm_read.t( "runSummary"  )->GetEntries();
    std::cout << nentries << " entries found here.\n";
    for ( int ientry = 0; ientry < nentries; ++ientry ) 
      {
	tm_read.entry( ientry, "runSummary" );
	std::cout << "the string was: " << (*rs) << "\n";
      }
  }


  std::cout << "--------------- variable I/O finished\n";
  return 0;
}
