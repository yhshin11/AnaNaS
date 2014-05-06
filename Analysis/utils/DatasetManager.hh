#ifndef DatasetManager_hh
#define DatasetManager_hh

//
// Author: Gautier Hamel de Monchenault 
//

#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TBits.h>

#include <map>
#include <string>
#include <vector>
#include <iostream>

#include <cassert>
using namespace std;

#include <RooDataSet.h>
#include <RooAbsReal.h>
#include <RooRealVar.h>
//using namespace RooFit;

class DatasetManager 
{
public:

  typedef map< string, RooRealVar* > MapVar;  
  typedef map< string, RooDataSet* > MapSet;

  enum Mode { kRead, kWrite, kSaved, kUnknown };

  // constructor
  DatasetManager( string name="dataset" ) 
    : _name(name), _f(0), _mode(kUnknown), _dataset(0), _rooArgSet() 
  {
  }

  ~DatasetManager() {}

  DatasetManager(  TFile* f, Mode mode = kRead )
  {
    setFile( f, mode );
  }

  void setFile( TFile* f, Mode mode = kRead )
  {
    assert( _f==0 );
    _f = f;
    _mode = mode;
  }

  TDirectory* getDirectory( string dir )
  {
    TDirectory* tdir_(_f);
    if( _f!=0 && dir!="" )
      {
	string pwd = gDirectory->GetPath();
	vector<string> dirs;
	const char* dir_ = dir.c_str();
	while( const char *slash_ = strchr(dir_,'/') )
	  {
	    size_t size_ = size_t(slash_-dir_);
	    char* workdir_ = new char[size_+1];
	    strncpy(workdir_, dir_, size_);
	    workdir_[size_] = 0; // end string
	    dirs.push_back(workdir_);
	    delete[] workdir_;
	    dir_ = slash_+1;
	  }
	
	for( size_t idir=0; idir<dirs.size(); idir++ )
	  {
	    string dir_ =  dirs[idir];
	    TDirectory* tmpdir = tdir_->GetDirectory( dir_.c_str(), false );
	    if( tmpdir==0 )
	      {
		tmpdir = tdir_->mkdir( dir_.c_str() );
		cout << "create directory with path " << tmpdir->GetPath() << endl;
	      }
	    tdir_ = tmpdir;
	    tdir_->cd();
	  }
	gDirectory->cd( pwd.c_str() );
      }
    return tdir_;
  }

  void createVar( const char* name, float min, float max )
  {
    assert( m_var.count( name )==0 );
    RooRealVar* var_ = new RooRealVar(name,name, min, max );
    m_var[name] = var_;
    _rooArgSet.add( *var_ );
  }

  void add( const char* name, float val, 
	    float min=-100000., float max=1000000. )
  {
    if( m_var.count( name )==0 )
      {
	createVar( name, min, max );
      }
    assert( m_var.count( name )==1 );
    m_var[name]->setVal( val );
  }

  void
  flush( string dir="" )
  {
    dataset( dir )->add( _rooArgSet );
  }
  
  RooDataSet*
  dataset( string dir )
  {
    size_t size =  dir.size();
    if( dir!="" && dir[size-1]!='/' )
      {
	dir += "/";
	size++;
      }
    RooDataSet* set;
    if( m_set.find(dir)==m_set.end() )
      {
	if( _mode==kWrite )
	  {
	    // the dataset does not exist, create it
	    string pwd = gDirectory->GetPath();
	    TDirectory* tdir_ = getDirectory( dir );
	    tdir_->cd();
	    cout << "create dataset " << _name << " in directory " 
		 << tdir_->GetName() << endl;    	    
	    set = new RooDataSet( _name.c_str(), _name.c_str(), _rooArgSet );
	    cout << "OK" << endl;
	    gDirectory->cd( pwd.c_str() );
	  }
	else
	  {
	    string key = dir;
	    key += _name;
	    _f->GetObject( key.c_str(), set );
	  }
	m_set[dir] = set; 
      }
    return m_set[dir];
  }

  void save()
  {
    string pwd = gDirectory->GetPath();
    MapSet::const_iterator it;
    for ( it = m_set.begin(); it != m_set.end(); ++it ) 
      {
        RooDataSet* t_ = it->second;
	TDirectory* tdir_ = getDirectory( it->first.c_str() );
	//        TDirectory* tdir_ = t_->GetDirectory();
        tdir_->cd();
        t_->Write();
      }
    gDirectory->cd( pwd.c_str() );
  }


private:

  // the name
  string _name;

  // the file
  TFile* _f;

  // the variables
  MapVar   m_var;

  // the datasets
  MapSet   m_set;

  // the mode
  Mode _mode;
  
  // the current dataset
  RooDataSet* _dataset;

  // the set of variables
  RooArgSet _rooArgSet;

};


#endif
