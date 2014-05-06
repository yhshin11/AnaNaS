#ifndef CandInfo_hh
#define CandInfo_hh

#include <map>

#include <TNamed.h>
#include <string>

#include "Analysis/utils/Counter.hh"
#include "Analysis/utils/ObjectStore.hh"

class CandInfo : public Counter<CandInfo>, public ObjectStore
{
public :

  // get the cand info for a given candidate
  static CandInfo* create()
  { return new CandInfo; }
  
  // set functions
  void setBool(  const std::string &,  bool  ); 
  void setInt(   const std::string &,   int  ); 
  void setFloat( const std::string &, float  ); 

  // accessors that abort if the values do not exist
  bool  getBool(  const std::string &  ) const;
  int   getInt(   const std::string &  ) const;
  float getFloat( const std::string &  ) const;

  // accessors that tests if the values exists
  bool getBool(   const std::string &, bool&  ) const;
  bool getInt(    const std::string &, int&   ) const;
  bool getFloat(  const std::string &, float& ) const;

  // print
  virtual void print( ostream& o ) const;
  
private:

  CandInfo();

  virtual ~CandInfo();

  std::map< std::string, bool  > _bool;
  std::map< std::string, int   > _int;
  std::map< std::string, float > _float;

  std::multimap< std::string, CandInfo* > _map;

ClassDef(CandInfo,0) 
};

#endif

