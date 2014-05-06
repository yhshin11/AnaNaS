#include <cassert>
#include <map>
#include <string>
#include <typeinfo>
#include <cmath>

#include <sstream>
#include <iostream>

#include <TFile.h>
#include <TBranch.h>
#include <TKey.h>
#include <TTree.h>
#include <TString.h>

#include "TupleManager.hh"

#include <iostream>
#include <TBits.h>

using namespace std;

// class Test {
// public:
//   static int value();
// };

// int Test::value()
// {
//   return 2;
// }

class Test
{
public:
  Test(): i(0), x(0), b(false) {}
  void test()
  {
    file  = new TFile( "test.root", "RECREATE" );
    tm_out.setFile( file, TupleManager::kWrite );    
    
    tm_out.setTree( "theTree", "theDir" );    
    tm_out.declare<int>("i",&i);
    tm_out.declare<float>("x",&x);
    tm_out.declare<bool>("b",&b);
    
    tm_out.setTree( "theOtherTree", "theDir" );    
    tm_out.declare<int>("i",&i);
    tm_out.declare<float>("x",&x);
    tm_out.declare<bool>("b",&b);
    
    for( int ii=0; ii<100; ii++ )
      {
	i = ii;
	x = 3.14*ii;
	b = (ii%2==0);
	
	tm_out.setTree( "theTree", "theDir" );    
	tm_out.fill();      
	
	i = 100-ii;
	x = (1-2*(ii%2))*3.14*ii;
	b = x>0;
	
	tm_out.setTree( "theOtherTree", "theDir" );    
	tm_out.fill();      
      }
  }
  ~Test() 
  {
    tm_out.save();
    file->Close();    
  }
private:
  TFile* file;
  TupleManager tm_out;
  int i;
  float x;
  bool b;
};

int main()
{
  Test t;
  t.test();
  return 0;

  int i(0);
  float x(0);
  bool b(true);

  TupleManager tm_out;
  TFile* file  = new TFile( "test.root", "RECREATE" );
  tm_out.setFile( file, TupleManager::kWrite );    

  tm_out.setTree( "theTree", "theDir" );    
  tm_out.declare<int>("i",&i);
  tm_out.declare<float>("x",&x);
  tm_out.declare<bool>("b",&b);

  tm_out.setTree( "theOtherTree", "theDir" );    
  tm_out.declare<int>("i",&i);
  tm_out.declare<float>("x",&x);
  tm_out.declare<bool>("b",&b);

  for( int ii=0; ii<100; ii++ )
    {
      i = ii;
      x = 3.14*ii;
      b = (ii%2==0);

      tm_out.setTree( "theTree", "theDir" );    
      tm_out.fill();      

      i = 100-ii;
      x = (1-2*(ii%2))*3.14*ii;
      b = x>0;

      tm_out.setTree( "theOtherTree", "theDir" );    
      tm_out.fill();      
    }
  tm_out.save();
  file->Close();

  return 0;

  typedef double Real;
  typedef vector<Real> VectorReal;
  {
    cout << "Test writing" << endl;
    VectorReal* v = new VectorReal;
    TFile* outputFile  = new TFile( "test.root", "RECREATE" );
    TupleManager tm;
    tm.setFile( outputFile, TupleManager::kWrite );
    for( size_t jj=0; jj<3;    jj++ )
      {
	tm.setTree( "tagada", "dir" );
	v->clear();
	for( size_t ii=0; ii<jj+1; ii++ )
	  {
	    v->push_back((jj+1)+0.1*(ii+1) );      
	  }
	tm.add< VectorReal* >( "v", &v );
	tm.flush();
      }
    tm.save();
    outputFile->Close();
  }
  {
    cout << "Test reading" << endl;
    TFile* inputFile  = new TFile( "test.root", "READ" );
    TupleManager tm;
    tm.setFile( inputFile, TupleManager::kRead );
    TTree* tree = tm.getTree( "tagada", "dir" );
    VectorReal* v = new VectorReal;
    tm.add< VectorReal* >( "v", &v );
    int n = tree->GetEntriesFast();
    for( int ii=0; ii<n; ii++ )
      {
	cout << "Entry " << ii << endl;
	tree->GetEntry( ii );
	for ( size_t i=0; i<v->size(); ++i) 
	  {
	    cout << "\tv[" << i << "]=" << (*v)[i] << endl;
	  }
      }
    inputFile->Close();
  }
  if(0)
  {
    cout << "Test reading 2" << endl;
    TFile* inputFile  = new TFile( "root/Diboson/WMultijet.txt.root", "READ" );
    TupleManager tm;
    tm.setFile( inputFile, TupleManager::kRead );
    TTree* tree = tm.getTree( "Wtuple", "data/W_en" );
    VectorReal* jEt     = new VectorReal;
    VectorReal* jDRL    = new VectorReal;
    VectorReal* jDPhiNu = new VectorReal;
    int run;  
    int event;
    float LPt;
    float LEta;
    float MET;
    float MT;
    float jRecoilPt;
    float jRecoilDphi; 
    float jRecoilDpt;
    tm.add< int >(    "run",   &run       );
    tm.add< int >(    "event", &event     );
    tm.add< float >(  "LPt",   &LPt       );
    tm.add< float >(  "LEta",  &LEta      );
    tm.add< float >(  "MET",   &MET        );
    tm.add< float >(  "MT",    &MT        );
    tm.add< float >(  "jRecPt",  &jRecoilPt   );
    tm.add< float >(  "jRecDphi",&jRecoilDphi );
    tm.add< float >(  "jRecDpt", &jRecoilDpt  );
    tm.add< VectorReal* >( "jEt", &jEt    );
    tm.add< VectorReal* >( "jDRL", &jDRL    );
    tm.add< VectorReal* >( "jDPhiNu", &jDPhiNu );
    int n = tree->GetEntriesFast();
    for( int ii=0; ii<n; ii++ )
      {
	tree->GetEntry( ii );
	cout << "Entry " << ii <<  " run/event " << run << "/" << event << endl;
	cout << "\tLPt=" << LPt << endl;
	cout << "\tLEta=" << LEta << endl;
	cout << "\tMET=" << MET << endl;
	cout << "\tMT=" << MT << endl;
	for ( size_t i=0; i<jEt->size(); ++i) 
	  {
	    cout << "\tJet " << i << endl;
	    cout << "\t\tjEt[" << i << "]=" << (*jEt)[i] << endl;
	    cout << "\t\tjDRL[" << i << "]=" << (*jDRL)[i] << endl;
	    cout << "\t\tjDPhiNu[" << i << "]=" << (*jDPhiNu)[i] << endl;
	  }
      }
    inputFile->Close();
  }

  if ( false )
    {
      TupleManager tm;
      
      std::string ofilename("test.root");
      std::cout << "Will write ntuple in output file " << ofilename << std::endl;
      TFile* outputFile  = new TFile( ofilename.c_str(), "RECREATE" );

      std::cout << "--------------- Tuple creation\n";
      //      TreeManager tm( "test.root", TreeManager::kWrite );
      tm.setFile( outputFile, TupleManager::kWrite );

      int ix;
      float x;
      std::string str;
      
      ix  = 1;
      x   = 3.1416;
      str = "/home/gautier/cms/AnaNaS/workdir/data/EG_skim/Ntuple_100_1_4wV_skim.root";

      std::string iName("BirthYear");
      std::string fName("FavoriteNumber");
      std::string sName("FirstName");

      string tname1 = "tagada";
      string dir1   = "top/dir1/dir2";

      string tname2 = "tagada";
      string dir2   = "top/dir1/dir3";

      for( size_t ii=0; ii<10; ii++ )
	{	  
	  tm.add< int   >( tname1, "Birth", &ix , dir1 );
	  tm.add< float >( tname1, "Name",  &x  , dir1 );

	  int dum = ii*ii+1;
	  tm.add< int >( tname1, "sqare", &dum, dir1 );

	  // was in sync
	  tm.flush( tname1, dir1 );

	  ix *= ii+1;
	  x += 1;
	}

      //      std::vector<int> V(1000);
      VectorReal* v = new VectorReal;
      //std::vector<int>* v = &V;
      for( size_t ii=0; ii<8; ii++ )
	{
	  tm.setTree( tname2, dir2 );


	  TBits *wbits = new TBits(8);
	  wbits->SetBitNumber(ii, kTRUE);
	  tm.add< TBits* >( "bits", &wbits );      

	  v->clear();
	  for( size_t jj=0; jj<ii; jj++ )
	    {
	      v->push_back(jj*111.111);
	    }
	  tm.add< VectorReal* >( "vect_real", &v );      
      
	  ostringstream o;
	  o << str
	    << "_" << ii;
	  string str_(o.str());
	  string* strptr = &str_;
	  tm.add< std::string* >( "strptr", &strptr );   
	  

	  //	  TString tstr(str.c_str());
	  // 	  tstr += "_";
	  // 	  tstr += ii;
	  //	  TString* tstrptr
	  // 	  tm.add< TString* >( "tstr", &tstrptr );
	    

	  tm.flush();

	}
      delete v;

      tm.save();
      
      tm.scan( tname1, dir1 );
      
      outputFile->Close();
      
      std::cout << "--------------- Tuple creation finished\n";
      
    }
  
  if( false )
    {
      TupleManager tm;

      std::string ifilename("test.root");
      std::cout << "Will read ntuple from input file " << ifilename << std::endl;
      TFile* inputFile  = new TFile( ifilename.c_str(), "READ" );

      std::cout << "--------------- Tuple reading\n";
      //      TreeManager tm( "test.root", TreeManager::kWrite );
      tm.setFile( inputFile, TupleManager::kRead );

      string tname1 = "tagada";
      string dir1   = "top/dir1/dir2";

      tm.setTree( tname1, dir1 );
      tm.scan( tname1,dir1 );

      int ix;
      float x;
      int dum;
      tm.add< int   >( "Birth", &ix );
      tm.add< float >( "Name", &x );
      tm.add< int   >( "square", &dum );

      TTree* tree;
      int n;

      tree = tm.tree();
      n = tree->GetEntriesFast();
      for( int ii=0; ii<n; ii++ )
	{
	  tree->GetEntry( ii );
	  cout << "*** " << ii << "" << " ***" << endl;
	  cout << "\tBirth="  << ix << endl;
	  cout << "\tName="   << x << endl;
	  cout << "\tsquare=" << dum << endl;
	}
      

      string tname2 = "tagada";
      string dir2   = "top/dir1/dir3";

      tm.setTree( tname2, dir2 );
      tree = tm.tree();
      tree->Print();


//       //      double  gold(0);
//       //      string  str(256,'0');
      TBits*  bits(0);
      string* strptr(0);
      VectorReal* vect_real = new VectorReal();
      TString tstr;

      tm.add< TBits*   >( "bits", &bits );
//       //      tm.add< string   >( "str", &str );
      tm.add< string*  >( "strptr", &strptr );
      tm.add< VectorReal* >( "vect_real", &vect_real );
//       //      tm.add< TString >( "tstr", &tstr );

      n = tree->GetEntriesFast();
      cout << "n= " << n << endl;
      for( int ii=0; ii<n; ii++ )
 	{
	  tree->GetEntry( ii );
// 	  //	  cout << "str=" << str << endl;
 	  cout << "bits=" << (*bits) << endl;
 	  cout << "strptr=" << (*strptr) << endl;
// 	  //	  cout << "tstr=" << tstr << endl;
 	  for ( size_t i=0; i<vect_real->size(); ++i) {
	    std::cout << "w " << i << " " << (*vect_real)[i] << "\n";
	  }
	}
    }

  return 0;
}
