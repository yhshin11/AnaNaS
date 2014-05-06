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
	
	    std::string pwd = gDirectory->GetPath();
	    std::vector<std::string> dirs;
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
	
	    gDirectory->cd( pwd.c_str() );
	
	  }
	m_histos[id] -> DirectoryAutoAdd( tdir_ ); 
	return (T*)(m_histos[id]);
      }
  }

  void save()
  {
    std::string pwd = gDirectory->GetPath();
    Map::const_iterator it;
    for ( it = m_histos.begin(); it != m_histos.end(); ++it ) 
      {
        TH1* h_ = it->second;
        TDirectory* tdir_ = h_->GetDirectory();
        tdir_->cd();
        h_->Write();
      }
    gDirectory->cd( pwd.c_str() );
  }

private:

  TFile* _f;
  Map m_templates;
  Map m_histos;

};


class HBool
{
  // an histogram and a reference to a boolean
  TH1*   _h;
  float* _v;  
  int*   _i;
  bool*  _b;
  bool   _u; // underflow
  bool   _o; // overflow

  // histogram 
  TAxis* _ax;
  float _min;
  float _max;
  int   _N;


public:

  void fill( float weight=1. )
  {
    assert( _h!=0 && _b!=0 );
    assert( _v!=0 || _i!=0 );
    
    if( !(*_b) ) return;

    float v_(0);
    if( _v ) v_ = (*_v);
    else     v_ = (*_i);

    if( _ax==0 )
      {
	_ax  =  _h->GetXaxis();
	_min = _ax->GetXmin();
	_max = _ax->GetXmax();
	_N   = _ax->GetNbins();
      }

    if( _u && v_<_min ) 
      _h->SetBinContent(  1, _h->GetBinContent(1)  + weight );
    else if( _o && v_>_max )
      _h->SetBinContent( _N, _h->GetBinContent(_N) + weight );
    else
      // keep it simple for the moment...
      _h->Fill( v_, weight );

  }

  HBool( TH1* h, float* v, bool* b, bool u=false, bool o=false ) 
    : _h(h), _v(v), _i(0), _b(b), _u(u), _o(o),
      _ax(0), _min(0), _max(0), _N(0) {}
  HBool( TH1* h,   int* i, bool* b, bool u=false, bool o=false ) 
    : _h(h), _v(0), _i(i), _b(b), _u(u), _o(o),
      _ax(0), _min(0), _max(0), _N(0) {}
  ~HBool() {}
};


#endif
