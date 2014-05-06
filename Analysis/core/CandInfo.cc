#include <cassert>
#include <iostream>
#include <algorithm>
using namespace std;

#include "Analysis/core/CandInfo.hh"

// instanciate the counter
#include "Analysis/utils/Counter.cc"

ClassImp(CandInfo)

CandInfo::CandInfo() 
{
}

CandInfo::~CandInfo()
{
}

void 
CandInfo::setBool( const string & key, bool bval )
{
  _bool[key] = bval;
}

void 
CandInfo::setInt( const string & key, int ival )
{
  _int[key] = ival;
}

void 
CandInfo::setFloat( const string & key, float val )
{
  _float[key] = val;  
}

bool
CandInfo::getBool( const string & key ) const
{
  bool val;
  //  assert( getBool( key, val ) );
  if( !getBool( key, val ) )
    {
      cout << "key=" << key << endl;
      abort();
    }
  return val;
}

int
CandInfo::getInt( const string & key ) const
{
  int val;
  if( !getInt( key, val ) )
    {
      cout << "key=" << key << endl;
      assert(0);
    }
  return val;
}

float
CandInfo::getFloat( const string & key ) const
{
  float val;
  if( !getFloat( key, val ) )
    {
      cout << key << endl;
      assert(0);
    }
  return val;
}

bool
CandInfo::getBool( const string & key, bool& val ) const
{
  if( _bool.count( key ) != 1 ) return false;
  val = _bool.find(key)->second;
  return true;
}

bool
CandInfo::getInt( const string & key, int& val ) const
{
  if( _int.count( key ) != 1 ) return false;
  val = _int.find(key)->second;
  return true;
}

bool
CandInfo::getFloat( const string & key, float& val ) const
{
  if( _float.count( key ) !=1 ) return false;
  val = _float.find(key)->second;
  return true;
}

void
CandInfo::print( ostream& o ) const
{
  {
    map< string, float >::const_iterator it;
    for( it=_float.begin(); it!=_float.end(); it++ )
      {
	o << it->first << "\t=\t" << it->second << endl;
	//	o << "Float[" << it->first << "]\t=\t" << it->second << endl;
      }
  }
  {
    map< string, int >::const_iterator it;
    for( it=_int.begin(); it!=_int.end(); it++ )
      {
	o << it->first << "\t=\t" << it->second << endl;
	//	o << "Int[" << it->first << "]\t=\t" << it->second << endl;
      }
  }
  {
    map< string, bool >::const_iterator it;
    for( it=_bool.begin(); it!=_bool.end(); it++ )
      {
	o << it->first << "\t=\t" << it->second << endl;
	//	o << "Bool[" << it->first << "]\t=\t" << it->second << endl;
      }
  }
}
