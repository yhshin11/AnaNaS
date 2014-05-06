#include <iostream>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <cassert>
#include <map>

using namespace std;

#include <TFile.h>
#include <TH1.h>
#include <TTree.h>

#include "Analysis/utils/Constants.hh"
#include "Analysis/utils/Config.hh"
#include "Analysis/utils/FileUtils.hh"
#include "Analysis/core/Sample.hh"
#include "Analysis/fit/Model.hh"

using namespace std;

ClassImp( Sample )

bool Sample::check   = false;
bool Sample::verbose = false;
bool Sample::useRootFile = false;
map< TString, Sample* > Sample::_samples = map< TString, Sample* >();

Sample::Sample( const char* name, const char* coll, const char* comments,
		float sig_tot, 
		float filter_eff, float br, float mean_kf, float n_tot )
  : _name(name), _coll(coll), _comments( comments), _n(0), _n_proc(0), _sampleList(0), _weights(0), _subsample("all"), _isData(false)
{
  assert( !this->isComposite() );

  cout << "Create Sample : " << _name << endl;

  _sig_tot    = sig_tot;
  _filter_eff = filter_eff;
  _br         = br;
  _mean_kf    = mean_kf;
  _n_tot      = n_tot;
  if( _sig_tot>0 )
    {
      assert( _filter_eff>0 && _filter_eff<=1 );
      assert( _br>0 && _br<=1 );
      assert( _mean_kf>0 );
      assert( _n_tot>0 );
    }

  setDataFiles_();

  // set the number of event in ntuples
  if( check ) setNFromFiles_();
  else 
    _n_proc = _n_tot;
  if( this->n()==0 )
    {
      if( useRootFile )
	{
	  cout << "n is not set " << endl;
	  cout << "-- try to read the root file " << endl;      
	  // FIXME !!!!  Temporary ugly solution :(
	  float n_ = checkRootFile( "Diboson", this->name(), "all", "n_electrons" );
	  cout << "WARNING " << this->name() << " set N=" << n_ << endl;
	  _n = n_;
	}
      else
	_n = _n_tot;
    }
}

Sample::Sample( const char* name, SampleList& sampleList ) 
  : _name(name), _coll("composite_sample"), _comments("no comment"),
    _sampleList(sampleList), _isData(false)
{
  // sort( _sampleList.begin(), _sampleList.end(), SortBySigma() );

  assert( this->isComposite() );
  _sig_tot    = 0;
  _filter_eff = 1;
  _br         = 1;
  _mean_kf    = 1;
  _n_tot      = 0;
  _n          = 0;
  _n_proc     = 0;

  cout << "Create composite sample " << endl;
  for( size_t ii=0; ii<_sampleList.size(); ii++ )
    {
      Sample* s_ = _sampleList[ii];
      assert( s_!=0 );
      float w_ = s_->w( integratedLumi() ); 
      _weights.push_back( w_ ); 
      cout << "Sample " << s_->name() << " w=" << w_ << endl; 
    }
}

Sample::Sample( const char* filename ) 
  : _name("test"), _coll("from_file"), _comments("no comment"), _isData(false)
{
  assert( !this->isComposite() );
  _sig_tot    = 0;
  _filter_eff = 1;
  _br         = 1;
  _mean_kf    = 1;
  _n_tot      = 0;
  _n          = 0;
  //  cout<<" The file is "<<filename<<endl;
  if( filename!=0 )
    filenames.push_back( filename );

}

Sample::Sample( const Sample& o )
  : _name(o._name), _coll(o._coll), _comments( o._comments),_n_tot(o._n_tot), _n(o._n), _sig_tot(o._sig_tot), _filter_eff(o._filter_eff), _br(o._br), _mean_kf(o._mean_kf), _sampleList(o._sampleList), _weights(o._weights), _subsample(o._subsample), _isData(o._isData)
{
}

bool
Sample::isComposite() const
{
  return _sampleList.size()>0;
}

void 
Sample::setNFromFiles_()
{
  if( _n>0 )
    {
      cout << "[" << this->name() << "] already set, n=" << _n << endl;
      return;
    }
  _n = checkDataFiles( this->name() );
  _n_proc = checkNProcFromDataFiles( this->name() );
}


void
Sample::print( ostream& os ) const
{
  char line_[512];
  float ratio_=1.;
  if( n_proc()>0 ) ratio_ = n()/n_proc();
  //  sprintf( line_,
  //	   "[%-10s] %7.0f/%7.0f/%7.0f sig(pb)=%-9.2f sig_tot(pb)=%-9.2f sum(L)(/pb)=%-9.2f ",
  //	   _name.Data(), n(), n_proc(), n_tot(), sigma(), sigma_tot(), integratedLumi() );
  sprintf( line_,
	   "%-10s %7.0f %6.2f %5.3f data data ",
	   _name.Data(), n(), sigma(), ratio_ );
  os << line_ << endl;
  //  if( verbose )
  //    {
  //     sprintf(line_,"    collection:%-s", _coll.Data() );
  //     os << line_ << endl;
  //     sprintf(line_,"    comments  :%-s", _comments.Data() );
  //     os << line_ << endl;
  //   }
}

float
Sample::trig_eff() const
{
  if( isComposite() ) return 1.;
  if( n_tot()==0 )    return 0.;
  if( n()==0     )    return 1.;
  return n()/n_tot();
}

float 
Sample::br() const
{
  return _br;
}

float 
Sample::filter_eff() const
{
  return _filter_eff;
}

float 
Sample::sigma_tot() const 
{ 
  return _mean_kf * _sig_tot; 
}

float 
Sample::sigma() const
{
  if( !isComposite() )
    {
      return sigma_tot() * br() * filter_eff(); 
    }
  else
    {
      return n_tot() / integratedLumi();
    }
}

float 
Sample::integratedLumi() const 
{ 
  if( !isComposite() )
    {
      return n_proc() / sigma(); 
    }
  else
    {
      return _sampleList[0]->integratedLumi();
    }
}

float 
Sample::n_tot( float lumi ) const
{
  return lumi * sigma();
}

float 
Sample::n( float lumi ) const
{
  if( _sig_tot==0 ) return _n_tot; 
  return n_tot( lumi ) * trig_eff();
}

float 
Sample::w( float lumi ) const
{
  if( _sig_tot==0 ) return 1;
  if( _isData ) return 1;
  return lumi / integratedLumi();
} 

float 
Sample::n_tot()       const 
{
  if( _n_tot==0 && isComposite() )
    {
      float n_tot_(0.);
      for( size_t ii=0; ii<_sampleList.size(); ii++ )
	{
	  n_tot_ += _weights[ii] * _sampleList[ii]->n_tot();
	}
      return n_tot_;
    }
  return _n_tot;
}

float 
Sample::n()       const 
{ 
  return (float)_n; 
}

float 
Sample::n_proc()       const 
{ 
  return (float)_n_proc; 
}

void 
Sample::getListOfStackedHists( vector< TH1* >&listOfStackedHists,
			       map< Sample*, TH1* >& mapOfHists,
			       float lumi,
			       const char* analysisName, 
			       const char* subSampleName, 
			       const char* histName
			       ) 
{
  //  if( !isComposite() ) return;

  //
  // First treat the sub-samples
  //
  for( size_t ii=0; ii<_sampleList.size(); ii++ )
    {
      Sample* s_ = _sampleList[ii];
      s_->getListOfStackedHists( listOfStackedHists, mapOfHists,
 				 lumi, analysisName, subSampleName, histName );
    }
  
  //
  // then get weighted histograms in the map
  //
  this->getMapOfHists( mapOfHists, 
		       lumi, analysisName, subSampleName, histName );

  //
  // get the combinatorics of histograms
  //
  vector< vector<size_t> > vec = combinatorics_( _sampleList.size() );

  //
  //
  //
  //  cout << "Sample " << name() << endl;

  TH1* h_stack(0);
  for( size_t ii=0; ii<vec.size(); ii++ )
    {
      TString stackName_ = TString(this->name());
      h_stack=0;

      vector<size_t>& vec_ = vec[ii];
      for( size_t jj=0; jj<vec_.size(); jj++ )
 	{
 	  size_t kk = vec_[jj];
	  Sample* s_ = _sampleList[ kk ];
	  assert( s_!=0 );
	  cout << "\t" << s_->name();
	  if( mapOfHists.count(s_)==0 )
	    {
	      cout << "-";
	      continue;
	    }
	  TH1* h_ =  mapOfHists[s_];
	  cout << "+";
	  
	  stackName_+="_";
	  stackName_+=s_->name();
	  if( h_stack==0 )
	    {
	      h_stack = (TH1*)h_->Clone();
	      h_stack -> SetDirectory(0);
	    }
	  else
	    {
	      TAxis* xaxis = h_stack->GetXaxis();
	      for( int iBin=1; iBin<=xaxis->GetNbins(); iBin++ )
		{
		  h_stack -> AddBinContent( iBin, h_->GetBinContent( iBin ) );
		}        
	    }    
 	}
      h_stack->SetName( stackName_ );
      listOfStackedHists.push_back( h_stack );
      //      cout << endl;
    }
  //
  // by constuction, the last histograms contains the full weighted sum
  // -> add it to the map
  //
  if( h_stack!=0 )
    {
      TString histName_ = TString(histName)+"_"+this->name()+"_";
      histName_ += int(lumi);
      TH1* h_ = (TH1*)h_stack->Clone( histName_ );
      if( mapOfHists.count( this )!=0 )
	{
	  assert( histName_==TString(mapOfHists[this]->GetName() ) );
	  return;
	}
      //      cout << this->name() << " add hist to the map with name " << histName_ << endl;
      //      h_stack->Print();
      //      h_->Print();
      mapOfHists[ this ] = h_;
    }
}

void 
Sample::getMapOfHists( std::map< Sample*, TH1* >& mapOfHists,
		       float lumi,
		       const char* analysisName, 
		       const char* subSampleName, 
		       const char* histName
		       )
{
  if( isComposite() )
    {
      for( size_t ii=0; ii<_sampleList.size(); ii++ )
	{
	  Sample* s_ = _sampleList[ii];
	  s_->getMapOfHists( mapOfHists,
			     lumi,
			     analysisName, 
			     subSampleName, 
			     histName );
	} 
    }
  else
    {
      TString sampleName_ = this->name();
      TString histName_ = TString(histName)+"_"+sampleName_+"_";
      histName_ += int(lumi);
      if( mapOfHists.count( this )!=0 )
	{
	  assert( histName_==TString(mapOfHists[this]->GetName() ) );
	  return;
	}
      cout << "analysis=" << analysisName << endl;
      cout << "sample=" << sampleName_ << endl;
      cout << "subsample=" << subSampleName << endl;
      cout << "hist=" << histName << endl;
      TH1* h_ = getClone( analysisName, sampleName_, subSampleName, histName );
      assert( h_!=0 );
      h_->SetName( histName_ );
      float w_ = this->w( lumi );
      cout << "w[" << sampleName_ << ",L=" << lumi << "]=" << w_ << endl;
      h_->Scale( w_ );
      //      h_->Print();
      mapOfHists[ this ] = h_;
    }
}

vector< vector<size_t> > 
Sample::combinatorics_( size_t n )
{
  vector< vector<size_t> > vec;
  for( size_t ii=0; ii<n; ii++ )
    {
      vector< vector<size_t> > vec_( vec );
      vec.push_back( vector<size_t>(1,n-ii-1) );
      for( size_t jj=0; jj<vec_.size(); jj++ )
	{
	  vector<size_t> vec__( vec_[jj] );
	  vec__.push_back( n-ii-1 );
	  vec.push_back( vec__ );
	}
    }
  return vec;
}

//
// static functions
//
int 
Sample::checkDataFiles( const char* sampleName )
{
  TString files_ = Config::dataPath + sampleName + "/Ntuple_*.root";
  vector< string > filenames_ = FileUtils::getFileList( files_.Data() );
  int n_(0); 
  for( size_t ii=0; ii<filenames_.size(); ii++ )
    {
      string filename_ = filenames_[ii];
      if( verbose )      cout << "\t" << filename_ << endl;
      try 
	{	  
	  int m_ = FileUtils::checkDataFile( filename_.c_str() ); 
	  //	  if( m_<=0 ) abort();  // m_=0 est possible (skims)
	  n_ += m_;
	}
      catch( string& badfile )
	{
	  cout << "Bad ntuple file " << badfile << endl;
	  abort();
	}
    }
  return n_;
}

int 
Sample::checkNProcFromDataFiles( const char* sampleName )
{
  TString files_ = Config::dataPath + sampleName + "/Ntuple_*.root";
  vector< string > filenames_ = FileUtils::getFileList( files_.Data() );
  int n_proc_(0); 
  for( size_t ii=0; ii<filenames_.size(); ii++ )
    {
      string filename_ = filenames_[ii];
      if( verbose )      cout << "\t" << filename_ << endl;
      try 
	{
	  int m_ = FileUtils::checkNProc( filename_.c_str() ); 
	  //	  if( m_<=0 ) abort();
	  n_proc_ += m_;
	}
      catch( string& badfile )
	{
	  cout << "Bad ntuple file " << badfile << endl;
	  abort();
	}
    }
  return n_proc_;
}

void
Sample::setDataFiles_()
{
  filenames = getDataFiles( this->name() );
}

vector< string > 
Sample::getDataFiles( const char* sampleName )
{
  TString files_ = Config::dataPath + sampleName + "/Ntuple_*.root";
  return FileUtils::getFileList( files_.Data() );
}


float
Sample::checkRootFile(   const char* analysisName,
			 const char* sampleName, 
			 const char* subSampleName,
			 const char* histName )
{
  // FIXME !!!!  Temporary ugly solution :(
  TH1* h_ = getClone( analysisName, sampleName, subSampleName, histName );
  float n_(0);
  if( h_!=0 )
    {
      n_ = h_->GetEntries();
      delete h_;
    }
  return n_;
}

bool
Sample::exists( const char* sampleName )
{
  TString dataPath = Config::dataPath + sampleName;
  if( fopen( dataPath.Data(), "r" )==0 )
    {
      if( verbose )
	cout << "Directory " << dataPath << " does not exists!" << endl;
      return false;
    }
  return true;
}

void
Sample::set( const char* sampleFile )
{  
  if( verbose ) cout << "Read the Sample file " << sampleFile << endl;
  ifstream fin;
  fin.open(sampleFile);
  
  int nline(0);
  char line_[512];
  char c_;
  while( fin.good() )
    {
      c_=fin.peek();
      //      cout << "c=" << c_ << endl;
      if( c_=='/' )
	{      
	  fin.getline( line_, 512 );
	  continue;
	} 
      nline++;
      
      TString sampleName;
      TString collectionName;
      int   nEvt_(0);
      float sig_tot_(0);
      float filter_eff_(1);
      float mean_kf_(1);
      float br_(1);
      
      fin >> sampleName;
      fin >> nEvt_;
      fin >> sig_tot_;
      fin >> filter_eff_;
      fin >> collectionName;
      fin.getline( line_, 512 );
      TString comments(line_);
      
      if( sampleName.Sizeof()==0 ) continue;
      if( nEvt_==0 ) continue;
      
      //GHM!!!      if( !Sample::exists( sampleName.Data() ) ) continue;
      
      if( sampleName==TString("Z_2e") )
	{
	  br_ = Constants::Z0_br_ee;
	  sig_tot_ /= br_;
	}
      else if( sampleName==TString("Z_2m") )
	{
	  br_ = Constants::Z0_br_mumu;
	  sig_tot_ /= br_;	  
	}
      else if( sampleName==TString("Z_2t") )
	{
	  br_ = Constants::Z0_br_tautau;
	  sig_tot_ /= br_;	  
	}
      else if( sampleName==TString("W_en") )
	{
	  br_ = Constants::W_br_enu;
	  sig_tot_ /= br_;	  
	}
      else if( sampleName==TString("W_mn") )
	{
	  br_ = Constants::W_br_munu;
	  sig_tot_ /= br_;	  
	}
      else if( sampleName==TString("ZZ_2l2n") )  
	{
	br_ = 2*Constants::Z0_br_nunu*Constants::Z0_br_ll;
	sig_tot_ /= br_;
	} 
      else if( sampleName==TString("ZZ_4l") ) 
	{
	  br_ =      
	    Constants::Z0_br_ee       * Constants::Z0_br_ee 
	    +     Constants::Z0_br_tautau   * Constants::Z0_br_tautau
	    +     Constants::Z0_br_mumu     * Constants::Z0_br_mumu
	    + 2 * Constants::Z0_br_ee       * Constants::Z0_br_mumu 
	    + 2 * Constants::Z0_br_ee       * Constants::Z0_br_tautau 
	    + 2 * Constants::Z0_br_mumu     * Constants::Z0_br_tautau;
	  sig_tot_ /= br_;
	}
      else if( sampleName==TString("WW_2l2n") )
	{
	  br_ = 
	    Constants::W_br_enu   * Constants::W_br_enu 
	    +     Constants::W_br_munu  * Constants::W_br_munu
	    +     Constants::W_br_taunu * Constants::W_br_taunu
	    + 2 * Constants::W_br_enu   * Constants::W_br_munu 
	    + 2 * Constants::W_br_enu   * Constants::W_br_taunu 
	    + 2 * Constants::W_br_munu  * Constants::W_br_taunu; 
	  sig_tot_ /= br_;
	}
      else if( sampleName==TString("WZ_3ln") )
	{
	  br_ = Constants::Z0_br_ll * Constants::W_br_lnu;
	  sig_tot_ /= br_;
	}
      
      cout << nline <<"\t";
      cout << sampleName << "\t";
      cout << nEvt_  << "\t";
      cout << sig_tot_ << "\t";
      cout << filter_eff_   << "\t";
      //      cout << generatorName << "\t";
      //      cout << collectionName << "\t";
      cout << endl;
      
      Sample* s_ = 
	new Sample( sampleName, collectionName, comments,
		    sig_tot_, filter_eff_, br_, mean_kf_, nEvt_ );
      _samples[sampleName] = s_; 
      //      s_->print( cout );

    }
  fin.close();
  
  for( map< TString, Sample* >::iterator it_=_samples.begin(); 
       it_!=_samples.end(); ++it_ )
    {
      Sample* s_ = it_->second;
      s_ -> print( cout );
    }
}

map< TString, Sample* >& 
Sample::samples( const char* sampleFile )
{
  if( _samples.size()==0 )
    {
      assert( sampleFile!=0 );

      set( sampleFile );
    }
  return _samples;
}

Sample*
Sample::get( const char* sampleName, const char* sampleFile )
{
  if( verbose ) cout << "Read the Sample file " << sampleFile << endl;
  ifstream fin;
  fin.open(sampleFile);
  
  int nline(0);
  char line_[512];
  char c_;
  while( fin.good() )
    {
      c_=fin.peek();
      if( c_=='/' )
	{      
	  fin.getline( line_, 512 );
	  continue;
	} 
      nline++;
      
      TString sampleName_;
      TString collectionName;
      int   nEvt_(0);
      float sig_tot_(0);
      float filter_eff_(1);
      float mean_kf_(1);
      float br_(1);
      
      fin >> sampleName_;
      fin >> nEvt_;
      fin >> sig_tot_;
      fin >> filter_eff_;
      fin >> collectionName;
      fin.getline( line_, 512 );
      TString comments(line_);
      
      if( sampleName_.Sizeof()==0 ) continue;
      if( sampleName_==TString(sampleName) )
	cout << "Sample " << sampleName << " exists " << endl;
      else 
	continue;
      cout << nline <<"\t";
      cout << sampleName << "\t";
      cout << nEvt_  << "\t";
      cout << sig_tot_ << "\t";
      cout << filter_eff_   << "\t";
      //      cout << generatorName << "\t";
      //      cout << collectionName << "\t";
      cout << endl;
      
      return 
	new Sample( sampleName, collectionName, comments,
		    sig_tot_, filter_eff_, br_, mean_kf_, nEvt_ );

    }
  return 0;
}

TH1* 
Sample::getHist( float L, 
		 const char* analysisName, 
		 const char* subSampleName, 
		 const char* histName )
{
  vector< TH1* > listOfStackedHists_;
  map< Sample*, TH1* > mapOfHists_;
  getListOfStackedHists( listOfStackedHists_, mapOfHists_, L, analysisName, subSampleName, histName );
  assert( mapOfHists_.count(this)==1 );
  TH1* h_ = (TH1*) mapOfHists_[this] -> Clone();
  h_->SetDirectory(0);
  return h_;
}

TH1* 
Sample::getHist( const char* analysisName, 
		 const char* subSampleName, 
		 const char* histName )
{
  return getHist( integratedLumi(), analysisName, subSampleName, histName );
}

TH1*
Sample::getClone( const char* analysisName, 
		  const char* sampleName, 
		  const char* subsample, 
		  const char* histName )
{
  TString rootFileName_ = Config::rootPath +
    analysisName + "/" + sampleName + ".root";
  
  TFile* rootFile_ = TFile::Open( rootFileName_, "READ" );
  if( rootFile_==0 ) return 0;
  TString hname_ = TString(subsample) + "/" + histName;
  const TH1* href_ = (TH1*) rootFile_->Get( hname_ );
  if( href_==0 )
    {
      cout << rootFile_ << " " << hname_ << endl;
      return 0;
    }
  TH1* h_ = (TH1*) href_-> Clone();
  h_->SetDirectory(0);
  rootFile_->Close();
  return h_;
}
