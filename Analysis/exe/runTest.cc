#include <iostream>
#include <errno.h>
#include <dirent.h>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cassert>
#include <string>
#include <map>

#include <dirent.h>
#include <sys/types.h>
using namespace std;

#include <TString.h>
#include <TFile.h>
#include <TH1.h>

#include "Analysis/utils/ObjectStore.hh"
#include "Analysis/utils/HistoManager.hh"

#include "Analysis/utils/Constants.hh"
#include "Analysis/utils/Config.hh"
#include "Analysis/core/Sample.hh"
#include "Analysis/core/Candidate.hh"
#include "Analysis/core/EventManager.hh"
#include "Analysis/tools/CandUtil.hh"

#include "Analysis/core/EventServer.hh"

int main(int argc, char **argv)
{

  float BR_ZZ_4e   = Constants::Z0_br_ee       * Constants::Z0_br_ee;
  float BR_ZZ_4m   = Constants::Z0_br_mumu       * Constants::Z0_br_mumu;
  float BR_ZZ_2e2m = 2 * Constants::Z0_br_ee       * Constants::Z0_br_mumu;
  float BR_ZZ_2e2n = 2 * Constants::Z0_br_ee       * Constants::Z0_br_nunu;
  float BR_ZZ_2m2n = 2 * Constants::Z0_br_mumu     * Constants::Z0_br_nunu;
  float BR_ZZ_2l2n = 2 * Constants::Z0_br_ll       * Constants::Z0_br_nunu;
  float BR_ZZ_4L   = BR_ZZ_4e + BR_ZZ_4m + BR_ZZ_2e2m; 
  float BR_ZZ_4l = 
    BR_ZZ_4L
    +     Constants::Z0_br_tautau   * Constants::Z0_br_tautau
    + 2 * Constants::Z0_br_ee       * Constants::Z0_br_tautau 
    + 2 * Constants::Z0_br_mumu     * Constants::Z0_br_tautau;

  float BR_ZZ = 7.406;
  
  cout << "ZZ " << BR_ZZ << endl;
  cout << "ZZ_2e2m " << BR_ZZ * BR_ZZ_2e2m << endl;
  cout << "ZZ_2e2n " << BR_ZZ * BR_ZZ_2e2n << endl;
  cout << "ZZ_2m2n " << BR_ZZ * BR_ZZ_2m2n << endl;
  cout << "ZZ_2l2n " << BR_ZZ * BR_ZZ_2l2n << endl;
  cout << "ZZ_4e   " << BR_ZZ * BR_ZZ_4e << endl;
  cout << "ZZ_4m   " << BR_ZZ * BR_ZZ_4m << endl;

  
  if(1) return 0;

  HistoManager _histoManager;
  TFile* _f  = new TFile( "test.root", "RECREATE" );
  _histoManager.setFile( _f );
  _histoManager.addTemplate<TH1F>( "essai",
				   new TH1F( "essai", "essai", 
					     10,0,10) );

//   string dir("top/dir1/dir2");
//   size_t size =  dir.size();
//   if( dir[size-1]!='/' )
//     {
//       // no slash at the end of string -- add one
//       dir += "/";
//       size++;
//     }

//   vector<string> dirs;
//   const char* dir_ = dir.c_str();
//   while( char *slash_ = strchr(dir_,'/') )
//     {
//       size_t size_ = size_t(slash_-dir_);
//       char* workdir_ = new char[size_+1];
//       strncpy(workdir_, dir_, size_);
//       workdir_[size_] = 0; // end string
//       dirs.push_back(workdir_);
//       delete[] workdir_;
//       dir_ = slash_+1;
//     }

//   TDirectory* tdir_(_f);
//   for( size_t idir=0; idir<dirs.size(); idir++ )
//     {
//       string dir_ =  dirs[idir];
//       cout << "directory : " << dir_ << endl;
//       TDirectory* tmpdir = tdir_->GetDirectory( dir_.c_str() );
//       if( tmpdir==0 ) tmpdir = tdir_->mkdir( dir_.c_str() );
//       tdir_ = tmpdir;
//       tdir_->cd();
//     }
  
  TH1* h_ = (TH1*)_histoManager.h<TH1F>( "essai", "cut", "top/dir1/dir2", "hello" );
  
  //  tdir_ = _f->GetDirectory( dir.c_str(), false );
  //  assert( tdir_!=0 );
  
  //  TH1* h_ = new TH1F("h","h",10,0,10);
  //  h_->DirectoryAutoAdd( tdir_ );
  
  for( int ii=0; ii<10; ii++ )
    h_->Fill(ii+0.5, ii );

  //    h_->Write();

  _histoManager.save();
  _f->Close();
  if(1) return 0;

  //  unsigned char isFile   =0x8;
  //  unsigned char isFolder =0x4;

  
  //  cout << "Hello, World" << endl;

  //   string run = argv[1];
  //   ostringstream out;
  //   out << run;
  //   string rep = "run" + out.str() +"/";

//   DIR *Dir;
//   struct dirent *DirEntry;
//   Dir = opendir(Root.c_str());

//   while(DirEntry=readdir(Dir))
//       {
	
// 	if ( DirEntry->d_type == isFile)
// 	{
// 	  cout <<"Found a File : " << DirEntry->d_name << endl;
// 	}
// 	else if(  DirEntry->d_type == isFolder )
// 	{
// 	  cout <<"Found a Directory : " << DirEntry->d_name << endl;
// 	}
// 	cout <<"----> " << DirEntry->d_name << " ===> " << DirEntry->d_type <<  endl;
//       }

//   if( 1 ) return 0;
 
  string Root ="/home/gpfs/manip/mnt/cms/mmarionn/data/MinBias/";


  vector<string> dirnames;
  vector<string> filenames;
 
  DIR *topdir;
  DIR *rundir;
  struct dirent *dirname;
  struct dirent *fichier;
  string file;
  //  cout<< Root+rep <<endl;
  topdir=opendir( Root.c_str() );
  if (topdir == NULL)
    {
      perror("erreur chemin\n" ); //Répertoire non valide
      return(-1);
    }
  
  while( (dirname=readdir(topdir) ) != NULL )
    {
      if (   strcmp(dirname->d_name, "."  ) == 0 
	     || strcmp(dirname->d_name, ".." ) == 0 )
	continue;
      if(  strstr( dirname->d_name, "run" ) == 0 ) continue;
      cout << "----> " << dirname->d_name << endl;
      string runpath = Root+dirname->d_name;
      rundir = opendir( runpath.c_str() );      
      while( (fichier=readdir(rundir)) != NULL )
	{
	  string path = runpath;
	  path += "/";
	  file= fichier->d_name;
	  if( strstr( file.c_str(), ".root" )==0 ) continue;
	  path += file.c_str(); //Chemin global du fichier	 
	  cout << path << endl;	  
	  filenames.push_back( path );// cout<<path<<endl; //Dans le vecteur
	}
    }


  if( 1 ) return 0;

  //  TString sample("QCD_EM2030");
  TString sample("Z_2e");
  TString sampleFile = Config::confPath + "samples.txt";
  TString collectionFile = Config::confPath + "collections.txt";

  // sample
  Sample* s_ = Sample::get( sample, sampleFile );
  if( s_==0 )
    {
      cout << "Sample \"" << sample << "\" does not exist!" << endl;
      return 1;
    }
  {
    EventServer::verbosity = 0;
    EventServer _e( s_->filenames, collectionFile );
    
    
    cout << "OK -- now loop " << endl;

    int ievt_(0);
    while( _e.nextEvent() )
      {
	ievt_++;
	if( ievt_>10 ) break;

	TBits hltBits_;
	assert( _e.get( "hlt", "hltBits", hltBits_ ) );
	//	cout << "\n" << ievt_ << " hlt:" <<  hltBits_ << "\n";

	// loop over electrons
	size_t nElectrons = _e.n("el");
	cout << "number of electrons " << nElectrons << endl;
	for( size_t ii=0; ii<nElectrons; ii++ )
	  {
	    assert( _e.load( "el", ii ) );
	    float pt_  = _e.get_f("el","pt");
	    float eta_ = _e.get_f("el","eta");
	    float phi_ = _e.get_f("el","phi");
	    size_t nEcalDep_ = (size_t) _e.get_i("el","nECALdep");
	    vector<float>& E_ECALdep_ = _e.get_vf("el","E_ECALdep");

	    cout << "el[" << ii << "]=";
	    cout << pt_ << "/";
	    cout << eta_ << "/";
	    cout << phi_ << " -- ";
	    cout << nEcalDep_ << "/" << E_ECALdep_.size() << endl;
// 	    for( size_t jj=0; jj<nEcalDep_; jj++ )
// 	      {
// 		cout << jj << "--> ";
// 		cout << _e.get_vf("el","dRECALdep",jj) << "/";
// 		cout << _e.get_vf("el","E_ECALdep",jj) << "/";
// 		cout << _e.get_vf("el","etaECALdep",jj) << "/";
// 		cout << _e.get_vf("el","phiECALdep",jj) << "/";
// 		cout << endl;
//	      }
	    TBits eidBits_;
	    assert( _e.get( "el", "electronIDs", eidBits_ ) );
	    cout << " eid:" <<  eidBits_ << "\n";
	  }
      }
  }
  
  cout << "\nEnd of test" << endl;

  cout << "Good bye, World" << endl;
  
  return 0;
}
