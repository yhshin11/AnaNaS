#ifndef TupleManager_hh
#define TupleManager_hh

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

class TupleManager 
{
public:

  typedef map< string, TTree* > Map;  

  enum Mode { kRead, kWrite, kSaved, kUnknown };

  // constructor
  TupleManager() : _f(0), _mode(kUnknown), _tree(0) {}

  TupleManager(  TFile* f, Mode mode = kRead )
  {
    setFile( f, mode );
  }

  void setFile( TFile* f, Mode mode = kRead )
  {
    assert( _f==0 );
    _f = f;
    _mode = mode;
  }

  ~TupleManager() {}

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

  string treeKey( string tname, string dir="" )
  {
    if( dir!="" && dir[dir.size()]!='/' )
      {
	dir += "/";
      }
    return dir+tname;
  } 
    
  string varKey( string tname, string vname, 
		      string dir="" )
  {
    string str = treeKey( tname, dir );
    str += "_";
    str += vname;
    return str;
  } 
    
  TTree* getTree( string tname, 
		  string dir="" )
  {
    size_t size =  dir.size();
    if( dir!="" && dir[size-1]!='/' )
      {
	dir += "/";
	size++;
      }
    string key = treeKey(tname,dir); 
    TTree* tree; 
    if( m_trees.find( key ) == m_trees.end() )
      {
	if( _mode==kWrite )
	  {
	    // the tree does not exist, create it
	    tree = new TTree( tname.c_str(), tname.c_str() );
	    TDirectory* tdir_ = getDirectory( dir );
	    tree->DirectoryAutoAdd( tdir_ );
	  }
	else
	  {
	    _f->GetObject( key.c_str(), tree );
	  }

	m_trees[key] = tree; 
      }
    else 
      {
	tree = m_trees[key];
      }
    _tree = tree;
    return tree;
  }

  void flush( string tname, string dir )
  {
    getTree( tname, dir )->SetEntries();
  }

  void flush()
  {
    assert( _tree!=0 );
    _tree->SetEntries();
    _tree=0;
  }
  
  void scan( string tname, string dir, Long64_t i=10 )
  {
    getTree( tname, dir )->Scan("","","",i);
  }

  void print( string tname, string dir )
  {
    getTree( tname, dir )->Print();
  }
  
  void save()
  {
    std::cout << "TupleManager -- Saving Tuples" << std::endl; 
    string pwd = gDirectory->GetPath();
    Map::const_iterator it;
    for ( it = m_trees.begin(); it != m_trees.end(); ++it ) 
      {
        TTree* t_ = it->second;
        TDirectory* tdir_ = t_->GetDirectory();
        tdir_->cd();
	//t_->Print();
        t_->Write();
      }
    gDirectory->cd( pwd.c_str() );
  }

  template < class T >
  TBranch* add( string tname, string vname, T* vaddress, string dir )
  {
    return add<T>( getTree( tname, dir ), vname, vaddress );
  }

  void 
  setTree(  string tname, string dir )
  {
    getTree( tname, dir );
  }

  TTree* 
  tree()
  {
    return _tree;
  }

  template < class T >
  TBranch* add( string vname, T* vaddress )
  {
    assert( _tree!=0 );
    return add<T>( _tree, vname, vaddress );
  }

  template < class T >
  TBranch* add( TTree* tree, string vname, T* vaddress )
  {
    TBranch* b = tree->FindBranch( vname.c_str() );
    if( b!=0 )
      {
	b->SetAddress( vaddress );
	if( _mode == kWrite )
	{
	  b->Fill();
	  b->ResetAddress();
	}
      }
    else
      {
	if ( _mode == kWrite ) 
	  {
	    b = tree->Branch( vname.c_str(), vaddress );
	    b->SetFile( _f );
	    b->Fill();
	    b->ResetAddress();
	  } 
      }
    return b;
  }

  template < class T >
  bool declare( string vname, T* vaddress )
  {
    assert( _tree!=0 );
    TBranch* b = _tree->FindBranch( vname.c_str() );
    if( b!=0 ) return false;
    if ( _mode != kWrite ) return false; 
    b = _tree->Branch( vname.c_str(), vaddress );
    b->SetFile( _f );
    return true;
  }

  bool
  fill()
  {
    assert( _tree!=0 );
    if( _mode != kWrite ) return false;
    _tree->Fill();
    _tree->SetEntries();
    return true;
  }

  void
  reset()
  {
    _f = 0;
    _mode = kUnknown;    
    m_trees.clear();
    _tree=0;
  }

private:

  TFile* _f;
  Map m_trees;
  Mode _mode;
  
  // current tree
  TTree* _tree;
};


#endif
