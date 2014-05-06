#ifndef HistoManager_hh
#define HistoManager_hh

//
// Author: Federico Ferri (Saclay) 2009
//

// 26/08/09 -- GHM: Directory management
// 29/09/09 -- GHM: pass strings instead of const char*
//  9/05/11 -- GHM: improve directory management (subdirectories)

// HistoManager histos;
// add a template histogram called ZZmass
// histos.addTemplate<TH1F>( "ZZmass", new TH1F( "ZZmass", "ZZmass", 100, 50., 150.) );
// [...]
// histos.h<TH1F>( "ZZmass", "mass_after_eleID" )->Fill( 4lMass );>
// [..]
// histos.h<TH1F>( "ZZmass", "mass_after_isolation" )->Fill( 4lMass );
// [...]
// std::string outPlotFile_ = "toto.root";
// histos.save( outPlotFile_.c_str() );

#include <TROOT.h>
#include <TH1.h>
#include <TFile.h>

#include <map>
#include <string>
#include <vector>
#include <iostream>

#include <cassert>

class HistoManager 
{
public:

  typedef std::map< std::string, TH1* > Map;  //ghm: TObject* -> TH1*

  HistoManager() : _f(0) {}

  void setFile( TFile* f ) { _f=f; }

  ~HistoManager()  {}

  template < class T > void addTemplate( std::string type, T* templ )
  {
    assert( m_templates.find( type ) == m_templates.end() );
    m_templates[ type ] = templ;
  }

  template < class T > T* h( std::string type, std::string name, 
			     std::string dir="", std::string prefix="" )
  {
    // check slash at end of directory name
    size_t size =  dir.size();
    if( dir[size-1]!='/' )
      {
	dir += "/";
	size++;
      }
    
    if( m_templates.find( type ) == m_templates.end() ) 
      { std::cout << type << std::endl; assert(0); }
    std::string hname;
    if( prefix!="" ) hname += ( prefix + "__" );
    hname += ( type + "__" );
    hname += name;
    std::string id(hname);
    if( dir!="" ) id += ( std::string("__") + dir );
    if ( m_histos.find( id ) != m_histos.end() ) 
      {
	return (T*)m_histos[id];
      } 
    else 
      {
	m_histos[id] = (T*)m_templates[type]->Clone(hname.c_str());

	TDirectory* tdir_(_f);
	if( _f!=0 && dir!="" )
	  {
	    std::vector<std::string> dirs;
	    const char* dir_ = dir.c_str();
	    while( char *slash_ = strchr(dir_,'/') )
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
		std::string dir_ =  dirs[idir];
		TDirectory* tmpdir = tdir_->GetDirectory( dir_.c_str(), false );
		if( tmpdir==0 )
		  {
		    tmpdir = tdir_->mkdir( dir_.c_str() );
		    std::cout << "create directory with path " << tmpdir->GetPath() << std::endl;
		  }
		tdir_ = tmpdir;
		tdir_->cd();
	      }
	  }
	m_histos[id] -> DirectoryAutoAdd( tdir_ ); 
	
	return (T*)(m_histos[id]);
      }
  }

  void save()
  {
    Map::const_iterator it;
    for ( it = m_histos.begin(); it != m_histos.end(); ++it ) 
      {
        TH1* h_ = it->second;
        TDirectory* tdir_ = h_->GetDirectory();
        tdir_->cd();
        h_->Write();
      }
  }

private:

  TFile* _f;
  Map m_templates;
  Map m_histos;

};


#endif
