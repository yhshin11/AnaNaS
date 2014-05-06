#ifndef CutUtils_hh
#define CutUtils_hh

#include <vector>
#include <sstream>
#include <iostream>
using namespace std;

#include <TBits.h>

class AbsCut
{
protected:

  string _name;
  string _type;

public:

  AbsCut( const string& name, const string& type )
    : _name(name),_type(type) {}

  virtual ~AbsCut() {}

  virtual bool operator==( const AbsCut& c2 ) const
  {
    return c2._name==_name && c2.cutType()==cutType();
  }

  virtual bool operator()()    const =0;

  virtual string cutType()     const =0;

  virtual void print()         const =0;

  virtual const string& name() const =0;

  virtual const string& type() const =0;

  //  virtual int getBits( TBits& bits, int ibit=0 ) const =0; 
  virtual int getBits( TBits& bits, int ibit=0 ) const
  {
    bits.SetBitNumber( ibit, this->operator()() );
    return ++ibit;
  } 

  virtual void getString( string& str ) const 
  { 
    str += name();
  }

  virtual bool ignore( AbsCut& cut ) const
  {
    if( this==&cut ) return true;
    return this->operator()();
  }

  virtual bool invert( AbsCut& cut ) const
  {
    bool out = this->operator()();
    if( this==&cut ) 
      {
	if( out ) out=false;
	else      out=true;
      }
    return out;
  }

  bool sequential( AbsCut& cut ) const
  {
    bool ok=true;
    return ignoreAfter( cut, ok );
  }

  virtual bool ignoreAfter( AbsCut& cut, bool& ok ) const
  {
    if( !ok ) return true;
    bool out = this->operator()();
    if( this==&cut )
      {
	if(ok) ok=false; else ok=true;
      }
    return out;
  }

};

template < class T >
class Cut : public AbsCut
{
public:

  Cut( const string& name, const string& type )
    : AbsCut(name,type) {}

  virtual ~Cut() {}

//   int getBits( TBits& bits, int ibit=0 ) const
//   {
//     bits.SetBitNumber( ibit, (*this)() );
//     return ++ibit;
//   } 

};

template < class T >
class OneSidedCut : public Cut<T>
{
public:

  OneSidedCut( const string& name, T& t, T val, const string& type )
    : Cut<T>(name,type), _t(t), _val(val) {}

  ~OneSidedCut() {}

  const string& name() const { return Cut<T>::_name; }
  const string& type() const { return Cut<T>::_type; }
  string cutType() const { return "OneSidedCut"; }

  void print() const
  {
    std::cout << cutType() << "[" << name() << "] " 
	      << _t << " " << type() << " " << _val << " -- "; 
    std::cout << ( this->operator()() ?  "true" : "false" )
	      << std::endl;
  }

  bool operator()() const
  {
    if( type()=="<"   ) return _t <  _val;
    if( type()=="<="  ) return _t <= _val;
    if( type()==">"   ) return _t >  _val;
    if( type()==">="  ) return _t >= _val;
    if( type()=="=="  ) return _t == _val;
    return false;
  }

private:

  T& _t;
  T  _val;
};

template < class T >
class TwoSidedCut : public virtual Cut<T>
{
public:

  TwoSidedCut( const string& name, T& t, T val1, T val2, const string& type )
    : Cut<T>(name,type), _t(t), _val1(val1), _val2(val2) {}

  ~TwoSidedCut() {}

  const string& name() const { return Cut<T>::_name; }
  const string& type() const { return Cut<T>::_type; }
  string cutType() const { return "TwoSidedCut"; }

  void print() const
  {
    char t1 = (type())[0];
    char t2 = (type())[1];
    std::cout << cutType() << "[" << name() << "] " 
	      << _t << " in " << t1 << _val1 << "," << _val2 << t2 << " --" ; 
    std::cout << ( this->operator()() ?  "true" : "false" )
	      << std::endl;
  }

  bool operator()() const
  {    
    if( type()=="[]"  ) return _t >= _val1 && _t <= _val2;
    if( type()=="[["  ) return _t >= _val1 && _t <  _val2;
    if( type()=="]]"  ) return _t >  _val1 && _t <= _val2;
    if( type()=="]["  ) return _t >  _val1 && _t <  _val2;
    if( type()=="![]" ) return _t <  _val1 || _t >  _val2;
    if( type()=="![[" ) return _t <  _val1 || _t <= _val2;
    if( type()=="!]]" ) return _t <= _val1 || _t >  _val2;
    if( type()=="!][" ) return _t <= _val1 || _t >= _val2;
    return false;
  }

private:

  T& _t;
  T  _val1;
  T  _val2;
};

class CompositeCut : public AbsCut
{
public:

  CompositeCut( const string& name, const string& type )
    : AbsCut(name,type) {}

  ~CompositeCut() {}

  const string& name() const { return AbsCut::_name; }
  const string& type() const { return AbsCut::_type; }
  string cutType() const { return "CompositeCut"; }

  bool add( AbsCut& cut )
  {
    _cuts.push_back( &cut );
    return true;
  }
  
  void print() const
  {
    std::cout << cutType() << "[" << name() << "] " << std::endl;
    for( size_t ii=0; ii<_cuts.size(); ii++ )
      {
	std::cout << "\t";
	if( ii>0 ) std::cout << type() << " ";
	_cuts[ii]->print();
      }
    std::cout << "---> " << ( this->operator()() ?  "true" : "false" )
	 << std::endl;
  }

  bool operator()() const
  {
    return isSelected(); 
  }

  bool isSelected( TBits* mask=0 ) const
  {
    bool out(false);
    if( type()=="&&" || type()=="and" )
      {
	out = true;
	for( size_t ii=0; ii<_cuts.size(); ii++ )
	  {
	    if( mask && !(*mask)[ii] ) 
	      {
		continue;
	      }
	    if( !(*_cuts[ii])() ) { out=false; break; }
	  }	
      }
    else 
      if( type()=="||" || type()=="or" )
	{
	  out = false;
	  for( size_t ii=0; ii<_cuts.size(); ii++ )
	    {
	      if( mask && !(*mask)[ii] ) 
		{
		  continue;
		}
	      if( (*_cuts[ii])() ) { out=true; break; }
	    }	
	}
    return out;
  }

  int getBits( TBits& bits, int ibit=0 ) const
  {
    for( size_t ii=0; ii<_cuts.size(); ii++ )
      {
	ibit = _cuts[ii]->getBits( bits, ibit );
      }
    bits.Compact();
    return ibit;
  } 

  void getString( string& str ) const 
  {
    str += "( ";
    for( size_t ii=0; ii<_cuts.size(); ii++ )
      {
	if( ii>0 ) 
	  {
	    str += " ";
	    str += type();
	    str += " ";
	  }
	_cuts[ii]->getString( str );
      }    
    str += " )";
  }

  bool invert( AbsCut& cut_ ) const
  {
    if( this==&cut_ ) return AbsCut::invert( cut_ );
    bool out(false);
    if( type()=="&&" || type()=="and" )
      {
	out = true;
	for( size_t ii=0; ii<_cuts.size(); ii++ )
	  {
	    if( !(_cuts[ii]->invert( cut_)) ) { out=false; break; }
	  }	
      }
    else 
      if( type()=="||" || type()=="or" )
	{
	  out = false;
	  for( size_t ii=0; ii<_cuts.size(); ii++ )
	    {
	      if( _cuts[ii]->invert( cut_ ) ) 
		//	      if( (_cuts[ii]->ignore( cut_ )) ) 
		{ out=true; break; }
	    }	
	}
    return out;
  } 

  bool ignore( AbsCut& cut_ ) const
  {
    if( this==&cut_ ) return true;
    bool out(false);
    if( type()=="&&" || type()=="and" )
      {
	out = true;
	for( size_t ii=0; ii<_cuts.size(); ii++ )
	  {
	    if( !(_cuts[ii]->ignore( cut_)) ) { out=false; break; }
	  }	
      }
    else 
      if( type()=="||" || type()=="or" )
	{
	  out = false;
	  for( size_t ii=0; ii<_cuts.size(); ii++ )
	    {
	      if( _cuts[ii]!=&cut_ && (_cuts[ii]->ignore( cut_ )) ) 
		//	      if( (_cuts[ii]->ignore( cut_ )) ) 
		{ out=true; break; }
	    }	
	}
    return out;
  } 

  bool ignoreAfter( AbsCut& cut_, bool& ok ) const
  {
    if( !ok ) return true;
    bool out(false);
    if( type()=="&&" || type()=="and" )
      {
	out = true;
	for( size_t ii=0; ii<_cuts.size(); ii++ )
	  {
	    if( !(_cuts[ii]->ignoreAfter( cut_, ok )) ) { out=false; break; }
	  }	
      }
    else 
      if( type()=="||" || type()=="or" )
	{
	  out = false;
	  for( size_t ii=0; ii<_cuts.size(); ii++ )
	    {
	      if( _cuts[ii]!=&cut_ && (_cuts[ii]->ignoreAfter( cut_, ok )) ) 
		//	      if( (_cuts[ii]->ignore( cut_ )) ) 
		{ out=true; break; }
	    }	
	}
    if( this==&cut_ )
      {
	if(ok) ok=false; else ok=true;
      }
    return out;
  } 

  bool NMinusOne( AbsCut& cut_ ) const
  {
    TBits mask;
    for( size_t ii=0; ii<_cuts.size(); ii++ )
      {	
	AbsCut* ptr1 = _cuts[ii];
	AbsCut* ptr2 = &cut_;
	if( ptr1 != ptr2 ) mask.SetBitNumber(ii);
      }        
    mask.Compact();
    //    std::cout << "mask" << mask << std::endl;
    return isSelected( &mask );
  } 

  bool NMinusOne( size_t icut ) const
  {
    TBits mask;
    for( size_t ii=0; ii<_cuts.size(); ii++ )
      {	
	if( ii != icut ) mask.SetBitNumber(ii);
      }        
    mask.Compact();
    return isSelected( &mask );
  } 

  size_t size() const
  {
    return _cuts.size();
  }

  string name( size_t icut ) const
  {
    if( icut>=_cuts.size() ) return name();
    ostringstream o;
    o << icut << "_" << _cuts[icut]->name();
    return o.str();
  }

  AbsCut* cut( size_t icut ) const
  {
    if( icut>=_cuts.size() ) return 0;
    return _cuts[icut];    
  }

  bool cutUpTo( size_t icut ) 
  {
    if( icut>= _cuts.size() ) return true;
    TBits mask;
    bool ok(true);
    for( size_t ii=0; ii<_cuts.size(); ii++ )
      {	
	if( ii > icut ) ok=false;
	if(ok) mask.SetBitNumber(ii);
      }        
    mask.Compact();
    return isSelected( &mask );
  }

private:

  std::vector< AbsCut* > _cuts;
};

#endif
