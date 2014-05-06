#include "Analysis/core/Vertex.hh"

#include <cassert>
#include <algorithm>
using namespace std;

// instanciate static counter
#include "Analysis/utils/Counter.cc"

#include "Analysis/utils/Constants.hh"

ClassImp(Vertex)

Vertex::Vertex() : _pos(0,0,0), _isPrimary(false), _info(0)
{
}

Vertex::Vertex( const TVector3& pos ) : _pos(pos), _isPrimary(false), _info(0)
{
}

Vertex::Vertex( const Vertex& o )
  : _pos(   o._pos ), _isPrimary(o._isPrimary), _info( o._info )
{
}

Vertex::~Vertex()
{
}

void
Vertex::setXYZ( float x, float y, float z )
{
  _pos.SetXYZ( x, y, z );
}

void
Vertex::print( ostream& o ) const
{
  o 
    << _pos.X() <<"/"
    << _pos.Y() <<"/"
    << _pos.Z() <<"/"
    << endl;
}
      
bool  
Vertex::operator==( const Vertex& o ) const
{
  return _pos == o._pos;
}

void 
Vertex::setIncomingCand( Candidate*  cand )
{
  _incoming.push_back( cand );
}

void 
Vertex::setOutgoingCand( Candidate*  cand )
{
  _outgoing.push_back( cand );
}

Candidate*
Vertex::incomingCand()
{
  if( _incoming.size()==0 ) return 0;
  return _incoming[0];
}

size_t 
Vertex::nOutgoingCand() const
{
  return _outgoing.size();
}

Candidate*
Vertex::outgoingCand( size_t ii )
{
  return _outgoing[ii];
}

Candidate*
Vertex::theMother() const
{
  if( _incoming.size()>0 )  return _incoming[0];
  return 0;
}

void
Vertex::setInfo( CandInfo* info )
{
  assert( _info==0 );
  _info = info;
}

int
Vertex::index() const
{
  int id_(-1);
  CandInfo* info_ = info();
  if( info_!=0 ) info_->getInt( "index", id_ );
  return id_;
}
