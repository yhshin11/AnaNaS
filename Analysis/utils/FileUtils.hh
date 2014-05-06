#ifndef FileUtils_hh
#define FileUtils_hh

#include <vector>
#include <list>
#include <string>
#include <map>
#include <iostream>

#include <TObject.h>

// an event is defined by its collection, run and event number
class EventHeader
{
public:
  
  std::string sample;
  int run;
  int event;

  EventHeader() : sample(""), run(0), event(0) {}
  bool operator<( const EventHeader& h2 ) const
  {
    return this->operator()( *this, h2 );
  }
  bool operator==( const EventHeader& h2 ) const
  {
    return h2.sample==sample && h2.run==run && h2.event==event;
  }
  bool operator()( const EventHeader& h1, const EventHeader& h2 ) const
  {
    if( h1.sample == h2.sample ) 
      if( h1.run == h2.run )
	return h1.event < h2.event;
      else
	return h1.run < h2.run; 
    else
      return h1.sample < h2.sample;
  }
  void print() const
  {
    std::cout << sample << " " 
	      << run    << " " 
	      << event  << std::endl;
  }
};

//typedef std::vector<EventHeader>          HeaderList;
typedef std::list<EventHeader>          HeaderList;
typedef std::map< std::string, std::vector<int> > EventList;

class FileUtils 
{
public:


  virtual ~FileUtils(){}
  
  static std::vector< std::string > getFileList( const char* str );

  static int checkDataFile( const char* filename );

  static void checkDataFile( const char* filename, int& analyzed, int& selected, int& nEventInJob );

  static int checkNProc( const char* filename );

  static void getListOfHeaders( const char* listOfEvents, HeaderList& headers );

  static void getListOfEvents( const char* listOfEvents, EventList& events );

  static void squeezerScriptFromListOfEvents( const char* listOfEvents );

  ClassDef( FileUtils, 0 )

};
#endif
