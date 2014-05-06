#include "Analysis/utils/StringUtils.hh"

#include <sstream>
using namespace std;

ClassImp( StringUtils )

string
StringUtils::stringIndex( size_t ii )
{
  stringstream ss;
  ss << "_" << ii;
  return ss.str();
}

string
StringUtils::indexedString( string str, size_t ii )
{
  string str_(str);
  str_ += stringIndex( ii );
  return str_;
}
