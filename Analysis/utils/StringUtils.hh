#ifndef StringUtils_hh
#define StringUtils_hh

#include <TObject.h>

#include <string>

class StringUtils 
{
public:
  
  virtual ~StringUtils(){}
  
  // returns a string with appended "_i", where i is a positive integer
  static std::string indexedString( std::string, size_t i );
  
  // returns the appendix "_i"
  static std::string stringIndex( size_t i );
  
  ClassDef( StringUtils, 0 )
};
#endif
