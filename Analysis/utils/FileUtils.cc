#include "Analysis/utils/Config.hh"
#include "Analysis/utils/FileUtils.hh"

#include <iostream>
#include <cassert>
#include <cstdlib>
#include <glob.h>
#include <algorithm>
#include <fstream>
using namespace std;

#include <TFile.h>
#include <TTree.h>

vector< string > 
FileUtils::getFileList( const char* str )
{
  glob_t files;
  glob( str, GLOB_TILDE, NULL, &files);
  vector< string > v;
  char **f;
  int fcnt;
  for( f = files.gl_pathv, fcnt = files.gl_pathc; fcnt; f++, fcnt-- ) 
    {
      v.push_back( *f );
    }
  return v;
}

void 
FileUtils::checkDataFile( const char* filename, int& analyzed_, int& selected_, int& nEventInJob_ )
{
  string _filename_(filename);
  if(_filename_.find("eos")!=(size_t)-1) _filename_ = "root://eoscms/"+_filename_;
  TFile* _file_ = TFile::Open( _filename_.c_str() );
  if( _file_==0 ) 
    {
      cout << "file " << _filename_ << " not found " << endl;
      abort();
    }
  
  // run summary
  static string str_ = "eventSummary";
  TTree* _ttree_ = (TTree*) _file_->Get( str_.c_str() );
  int n_ = (int)(_ttree_->GetEntriesFast());

  TTree* _runSum_ = (TTree*) _file_->Get( "runSummary" );
  //  int selected_;
  _runSum_->SetBranchAddress("selected",&selected_);
  //  int analyzed_;
  _runSum_->SetBranchAddress("analyzed",&analyzed_); 
  //  int nEventInJob_;
  _runSum_->SetBranchAddress("nEventInJob",&nEventInJob_);
  
  _runSum_->GetEntry(0);
  if( n_ != selected_ )
    cout << " Warning, number of selected events is not equal to filled event in file "<< n_ << " <-> " << selected_ << endl;
  if( analyzed_ != selected_ )
    cout << " Warning, number of selected events is not equal to the number of analyzed events " << analyzed_ << " <-> " << selected_ << endl;
  if( analyzed_ > nEventInJob_ )
    cout << " Warning, number of analyzed events is larger to the number of events in job " << analyzed_ << " <-> " << nEventInJob_ << endl;
  
  _file_->Close();

  //  return n_;
}  

int
FileUtils::checkDataFile( const char* filename )
{
  string _filename_(filename);
  if(_filename_.find("eos")!=(size_t)-1) _filename_ = "root://eoscms/"+_filename_;
  TFile* _file_ = TFile::Open( _filename_.c_str() );
  if( _file_==0 ) 
    {
      cout << "file " << _filename_ << " not found " << endl;
      abort();
    }
  
  // run summary
  static string str_ = "eventSummary";
  TTree* _ttree_ = (TTree*) _file_->Get( str_.c_str() );
  int n_ = (int)(_ttree_->GetEntriesFast());

  TTree* _runSum_ = (TTree*) _file_->Get( "runSummary" );
  int selected_;
  _runSum_->SetBranchAddress("selected",&selected_);
  
  _runSum_->GetEntry(0);
  if( n_ != selected_ )
    cout << " Warning, number of selected events is not equal to filled event in file "<< n_ << " <-> " << selected_ << endl;
  
  _file_->Close();

  return n_;
}  

int 
FileUtils::checkNProc( const char* filename )
{
  string _filename_(filename);
  if(_filename_.find("eos")!=(size_t)-1) _filename_ = "root://eoscms/"+_filename_;
  TFile* _file_ = TFile::Open( _filename_.c_str() );
  if( _file_==0 ) 
    {
      cout << "file " << _filename_ << " not found " << endl;
      return 0;
    }
  
  // run summary
  static string str_ = "eventSummary";
  
  TTree* _runSum_ = (TTree*) _file_->Get( "runSummary" );
  
  int nInJob_;
  _runSum_->SetBranchAddress("nEventInJob",&nInJob_);
  // _runSum_->SetBranchAddress("selected",&nInJob_); //!!!GHM 10/05/11 !!!
  _runSum_->GetEntry(0);  
 
  _file_->Close();

  return nInJob_;
}

void
FileUtils::getListOfHeaders( const char* listOfEvents, HeaderList& headers )
{

  // open the file 
  ifstream infile;
  infile.open( listOfEvents, ifstream::in );
  if( infile.fail() )
    {
      cout << "Event file does not exist -- ignore " << endl;
      return;
    }
  
  // parse items
  TString dataset;
  int currun(0), ievt(0);
  EventHeader header;

  char item[50];
  while(1)
    {
      infile >> item;
      if( infile.eof() ) break;

      TString str(item);
      if( isalpha(item[0]) ) 
	{
	  // new data set
	  dataset = str;
	  if( currun!=0 )
	    {
	      currun = 0;
	      ievt = 0;
	    }
	  header.sample = dataset;
	  header.run    = currun;
	  header.event  = 0;
	}
      else
	{
	  header.run = str.Atoi();
	  if( header.run!=currun )
	    {
	      currun = header.run;
	      ievt = 0;
	    }
	  infile >> item;
	  TString str_(item);
	  assert( str_.IsDigit() );
	  header.event = str_.Atoi();
	  headers.push_back(header);
	  ievt++;
	}
    }

  headers.sort();
  headers.unique();
  cout << "Selection: " << headers.size() << " events." << endl;

  HeaderList::iterator it;
  for (it=headers.begin(); it!=headers.end(); ++it)
    it->print();
  
  infile.close();
}

void
FileUtils::getListOfEvents( const char* listOfEvents, EventList& events )
{
  HeaderList headers;
  FileUtils::getListOfHeaders( listOfEvents, headers );
  if( headers.empty() ) return;
  
  TFile* file;
  TTree* tree;

  bool ok(true);
  while( ok )
  {
    EventHeader header = headers.front();  

    string sample = header.sample;
    TString files_ = Config::dataPath;
    files_ += sample;
    files_ += "/Ntuple_*.root";
    vector< string > filenames = FileUtils::getFileList( files_ );
    int nevt(0);
    for( size_t ii=0; ii<filenames.size(); ii++ )
      {
	string filename_ = filenames[ii];
	if(filename_.find("eos")!=(size_t)-1) filename_ = "root://eoscms/"+filename_;
	file = TFile::Open( filename_.c_str(), "READ" );
	file->GetObject( "eventSummary", tree );
	tree->FindBranch( "run"   )->SetAddress( &header.run );
	tree->FindBranch( "event" )->SetAddress( &header.event );
	int n = tree->GetEntriesFast();
	for( int jj=0; jj<n; jj++ )
	  {
	    tree->GetEntry( jj );
	    bool found = 
	      binary_search( headers.begin(), headers.end(), header );
	    if( found ) 
	      {
		headers.remove(header);
		nevt++;
		events[filename_].push_back(jj);
	      }
	  }
	if( file ) file->Close();
	if( headers.size()==0 ) break;
	if( headers.front().sample!=sample ) break;
      }
    if( headers.size()!=0 )
      {
	while( headers.front().sample==sample )
	  {
	    cout << "Event not found -- ";
	    headers.front().print();
	    headers.pop_front();
	    if( headers.size()==0 )
	      {
		ok = false;
		break;
	      }
	  }
      }    
    else 
      {
	ok = false;
      }
  }
}

void
FileUtils::squeezerScriptFromListOfEvents( const char* listOfEvents )
{
  ofstream o("dum.sh" );
  o << "#!/bin/bash" << endl;
  o << "outbase=\"$ANANAS/workdir/data\"" << endl;
  o << "outcoll=$outbase/$1_`date +%m%d%y`" << endl;
  o << "mkdir $outcoll" << endl;
  o << "squeezer_(){" << endl;
  o << "TMPFILE=`mktemp -q /tmp/eventList.XXXXXX`" << endl;
  o << "cat /dev/null > \"$TMPFILE\"" << endl;
  o << "echo \"$2 $3\" >> $TMPFILE" << endl;
  o << "echo \"temporary event file: $TMPFILE\"" << endl;
  o << "cat $TMPFILE" << endl;
  o << "$ANANAS/Squeezer/squeezer --event-list-from $TMPFILE $1 $outcoll/Ntuple_$2_$3.root" << endl;
  o << "}" << endl;
  
  HeaderList headers;
  FileUtils::getListOfHeaders( listOfEvents, headers );
  if( headers.empty() ) return;
  
  TFile* file;
  TTree* tree;

  char line_[512];
  vector< string > events_found;

  bool ok(true);
  while( ok )
  {
    EventHeader header = headers.front();  

    string sample = header.sample;
    TString files_ = Config::dataPath;
    files_ += sample;
    files_ += "/Ntuple_*.root";
    vector< string > filenames = FileUtils::getFileList( files_ );
    int nevt(0);
    for( size_t ii=0; ii<filenames.size(); ii++ )
      {
	string filename_ = filenames[ii];
	if(filename_.find("eos")!=(size_t)-1) filename_ = "root://eoscms/"+filename_;
	file = TFile::Open( filename_.c_str(), "READ" );
	file->GetObject( "eventSummary", tree );
	tree->FindBranch( "run"   )->SetAddress( &header.run );
	tree->FindBranch( "event" )->SetAddress( &header.event );
	int n = tree->GetEntriesFast();
	for( int jj=0; jj<n; jj++ )
	  {
	    tree->GetEntry( jj );
	    bool found = 
	      binary_search( headers.begin(), headers.end(), header );
	    if( found ) 
	      {
		headers.remove(header);
		nevt++;
		sprintf( line_, "%s\t%10d\t%10d\n", filename_.c_str(), header.run, header.event );
		events_found.push_back(line_);
	      }
	  }
	if( file ) file->Close();
	if( headers.size()==0 ) break;
	if( headers.front().sample!=sample ) break;
      }
    if( headers.size()!=0 )
      {
	while( headers.front().sample==sample )
	  {
	    cout << "Event not found -- ";
	    headers.front().print();
	    headers.pop_front();
	    if( headers.size()==0 )
	      {
		ok = false;
		break;
	      }
	  }
      }    
    else 
      {
	ok = false;
      }
  }
  cout << "\n\n";
  cout << "events found " << endl;
  for( size_t ii=0; ii<events_found.size(); ii++ )
    {
      o << "squeezer_ " << events_found[ii];
    }
}
