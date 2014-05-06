#include <cassert>
#include <map>
#include <string>
#include <typeinfo>

#include <iostream>
#include <cstdlib>

#include "TFile.h"
#include "TBranch.h"
#include "TKey.h"
#include "TTree.h"

class TreeManager {
public:

  typedef std::map< std::string, TTree * > Map;

  enum Mode { kRead, kWrite, kSaved, kUnknown };

  TreeManager( TFile* f=0, TreeManager::Mode mode = TreeManager::kRead ) : _f(f), _mode(mode), _own(false)
  {
  }

  TreeManager( const char * fileName, TreeManager::Mode mode = TreeManager::kRead ) : _f(0)
  {
    open( fileName, mode );
  }

  ~TreeManager()
  {
    save();
    close();
  }

  // GHM
  void close()
  {
    if( _f==0 ) return;
    if( !_own ) return;
    if( _mode != kRead ) return;
    _f->Close();
    _f=0;
  }

  bool open( const char * fileName, TreeManager::Mode mode )
  {
    if ( _f == 0 ) {
      std::string option;
      if ( mode == kRead ) {
	_mode = kRead;
	option = "READ";
      } else if ( mode == kWrite ) {
	_mode = kWrite;
	option = "RECREATE";
      } else {
	std::cerr << mode << ": unkown mode to open file.\n";
	return false;
      }
      _f = TFile::Open( fileName, option.c_str() );
      _own = true;
      if ( _f->IsOpen() ) {
	if ( mode == kRead ) {
	  readTrees( _f );
	  return true;
	}
      } 
    } else {
      std::cerr << "A file has already been open. Please check your code.\n";
      _own = false;
      return false;
    }
    return true; //GHM
  }

  void save( bool printTrees = false )
  {
    if ( _mode != kWrite || _mode == kSaved ) return;
    std::string pwd = gDirectory->GetPath();
    _f->cd();
    std::cout<<pwd<<"   "<<_f->GetName()<<std::endl;
    if ( printTrees ) {
      for ( Map::const_iterator it = _trees.begin(); it != _trees.end(); ++it ) {
	it->second->Print();
      }
    }
    _f->Write();
    gDirectory->cd( pwd.c_str() );
    _mode = kSaved;
  }

  // get all the trees in a directory, recursively
  void readTrees( TDirectory *d )
  {
    TIter nextkey = d->GetListOfKeys();
    TKey * k = 0;
    while ( (k = (TKey*)nextkey()) ) {
      if ( std::string( k->GetClassName() ) == "TDirectoryFile" ) {
	readTrees( (TDirectory*)d->Get( k->GetName() ) );
      } else if ( std::string( k->GetClassName() ) == "TTree" ) {
	_trees[ k->GetName() ] = (TTree*)d->Get( k->GetName() );
      }
    }
  }

  void addTree( const char * name )
  {
    assert( _mode == kWrite );
    std::string pwd = gDirectory->GetDirectory( "" )->GetPath();
    _f->cd();
    std::string tname = name;
    _trees[ tname ] = new TTree( name, name );
    gDirectory->cd( pwd.c_str() );
  }

  // if mode == kWrite fill the branch `branchName'
  //    with the object at the address `address' 
  // if mode == kRead  set the address for the 
  //    branch `branchName' to `address'
  template < class T >
  void branchio( const char * treeName, const char * branchName, T * address )
  {
    Map::const_iterator it = _trees.find( treeName );
    if( it == _trees.end() )
      {
	//	std::cout << treeName << "/" << branchName << std::endl; 
	//	assert(0);
	return;
      }
    // following line modified by GHM, 14/05/11
    if( _mode==kRead && it->second->GetNbranches()==0 ) return;
    //    if( it->second->GetNbranches()==0 ) return;

    TBranch * b = it->second->FindBranch( branchName );
    if ( b != 0 ) {
      b->SetAddress( address );
      if ( _mode == kWrite ) {
	b->Fill();
	b->ResetAddress();
      }
    } else {
      // assert ( _mode == TreeManager::kWrite );
      if ( _mode == kWrite ) {
	TBranch * b = it->second->Branch( branchName, address );
	b->SetFile( _f );
	b->Fill();
	b->ResetAddress();
      } else {
	std::cerr << "Cannot find branch " << branchName 
		  << " in tree " << treeName
                  << ".\n";
          //		  << ". Abort.\n";
	// exit(1);
	return;
      }
    }
  }

  // fill the branch `branchName' with `value' for T
  template < class T >
  void branchio( const char * treeName, const char * branchName, T value )
  {
    if ( _mode != kWrite ) {
      std::cerr << "This method is usable only for writing."
	" Please pass the pointer to your object. Abort.\n";
      //GHM                    exit(2);
      return;
    }
    branchio<T>( treeName, branchName, &value );
  }

  // synchronize the tree `treeName' and check that 
  // all the branches have been filled
  void sync( const char * treeName )
  {
    Map::const_iterator it = _trees.find( treeName );
    assert( it != _trees.end() );
    TTree * t = it->second;
    t->SetEntries();
  }

  // get the entry `ientry' for the tree `treeName'
  // TreeManager::branchio has to be previously called for all the
  // branches to be retrieved
  void entry( size_t ientry, const char * treeName = 0 )
  {
    if ( treeName != 0 ) {
      // read a specific tree
      Map::const_iterator it = _trees.find( treeName );
      //      assert( it != _trees.end() );
      if( it== _trees.end() ){
	std::cerr << "Tree " << treeName << " not found!\n";
        abort();
      }
      it->second->GetEntry( ientry );
    } else {
      // read all the trees
      for ( Map::const_iterator it = _trees.begin(); it != _trees.end(); ++it ) {
	((TTree*)it->second)->GetEntry( ientry );
      }
    }
  }

  // get the number of entries for the tree `treeName'
  size_t entries( const char * treeName )
  {
    Map::const_iterator it = _trees.find( treeName );
    assert( it != _trees.end() );
    return it->second->GetEntriesFast();
  }

  // scan the tree `treeName' 
  // or all the trees if `treeName' is null
  void scan( const char * treeName = 0 )
  {
    if ( treeName != 0 ) {
      Map::const_iterator it = _trees.find( treeName );
      assert( it != _trees.end() );
      it->second->Scan();
    } else {
      for ( Map::const_iterator it = _trees.begin(); it != _trees.end(); ++it ) {
	it->second->Scan();
      }
    }
  }

  // return a pointer to the tree `treeName'
  TTree * t( const char * treeName )
  {
    Map::const_iterator it = _trees.find( treeName );
    
    if ( it != _trees.end() ) {
      return it->second;
    } else {
      return 0;
    }     
  }

private:

  Map _trees;
  TFile *_f;
  Mode _mode;
  bool _own;
};
