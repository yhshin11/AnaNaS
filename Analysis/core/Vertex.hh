#ifndef Vertex_hh
#define Vertex_hh

#include "Analysis/utils/ObjectStore.hh"
#include "Analysis/utils/Counter.hh"
#include "Analysis/core/CandInfo.hh"

#include <iostream>

#include <TVector3.h>

class Candidate;

class Vertex : public ObjectStore, public Counter<Vertex>
{

public:

  // default construction
  static Vertex* create() 
  { return new Vertex; }

  // construction from a point in space
  static Vertex* create( const TVector3& pos ) 
  { return new Vertex(pos); }

  // cloning
  Vertex* clone() 
  { return new Vertex(*this); }
  
  // modifiers
  void setXYZ( float x, float y, float z );

  // operators
  bool  operator==( const Vertex& ) const;
  
  const TVector3&  pos() const { return _pos; }
    
  // Prints
  void print( ostream& o ) const;

  void setIncomingCand( Candidate*  cand );
  void setOutgoingCand( Candidate*  cand );

  Candidate* incomingCand();
  size_t nOutgoingCand() const;
  Candidate* outgoingCand( size_t ii );
  
  Candidate* theMother() const;

  bool isPrimary() const { return _isPrimary; }
  void setAsPrimary()   { _isPrimary = true; }
  void unsetAsPrimary() { _isPrimary = false; }

  // extra information
  CandInfo* info() const { return _info; }
  CandInfo* info()       { return _info; }

  // assign info block
  void setInfo( CandInfo* info );

  // vertex index
  int index() const;

private:

  Vertex();

  // constructor from momentum, charge and mass
  Vertex( const TVector3& pos );

  // copy constructor
  Vertex( const Vertex& );

  // destructor
  ~Vertex();

private:

  TVector3 _pos;

  bool _isPrimary;

  std::vector< Candidate* >  _incoming;
  std::vector< Candidate* >  _outgoing;

  // info block
  CandInfo* _info;

  ClassDef(Vertex,0) 
};
#endif




