#include <cassert>
#include <string>
#include <glob.h>
#include <algorithm>
#include <fstream>
#include <iostream>

#include <vector>
#include <map>

#include <TFile.h>
#include <TBranch.h>
#include <TTree.h>

#include "Analysis/utils/Config.hh"
#include "Analysis/utils/FileUtils.hh"
#include "Analysis/core/EventServer.hh"

using namespace std;

int main()
{
  if ( true )
    {
      cout << "Hello, World" << endl;

      string listOfEvents = "ZZ_2e2n_2.sel";

      EventList  events;

      //      FileUtils::getListOfEvents( listOfEvents.c_str(), events );
      FileUtils::squeezerScriptFromListOfEvents( listOfEvents.c_str() );
      
      if( 1 ) return 0;
	
      EventList::iterator it;
      for ( it=events.begin() ; it != events.end(); it++ )
	{
	  string filename_ = (*it).first;
	  cout << filename_ << endl;

	  vector<int> vect_ = (*it).second;
	  for( size_t ii=0; ii<vect_.size(); ii++ )	    
	    {
	      cout << "\t" << vect_[ii] << endl;
	    }
	}

      if( 1) return 0;

      cout << "Create the Event Server " << endl;
      
      string configPath = Config::confPath + string("collections.txt"); 
      EventServer es( events, configPath.c_str()  );
      
      while( es.nextEvent() )
	{ 
	  es.oneLine(cout);
	}
      

//       EventList::iterator it;
//       int ifile_ = 0;
//       for ( it=events.begin() ; it != events.end(); it++ )
// 	{
// 	  string filename_ = (*it).first;
// 	  cout << filename_ << endl;
	  
// 	  //	  es.openFile( filename_.c_str() );
	  
// 	  vector<int> vect_ = (*it).second;
// 	  for( size_t ii=0; ii<vect_.size(); ii++ )	    
// 	    {
// 	      cout << "\t" << vect_[ii] << endl;
// 	      assert( es.eventInFile( ifile_, vect_[ii]+1 ) );
// 	      es.oneLine(cout);
// 	    }
// 	  ifile_++;
// 	}    
    }
  return 0;
}
