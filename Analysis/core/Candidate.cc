#include "Analysis/core/Candidate.hh"

#include <cassert>
#include <iostream>
#include <algorithm>
using namespace std;

// instanciate static counter
#include "Analysis/utils/Counter.cc"

#include "Analysis/utils/Constants.hh"
#include "Analysis/utils/KineUtils.hh"

ClassImp(Candidate)

size_t Candidate::_curUid=0;
map<size_t,const Candidate*> Candidate::_baseCand;

void
Candidate::reset()
{
  _curUid=0;
  _baseCand.clear();
}

Candidate::Candidate()
{
  init();
}

Candidate::Candidate( const TVector3& mom,
		      float charge,
		      float mass,  
		      Vertex* vtx )
{
  init();
  setP3( mom );
  _q = charge;
  _m = mass;  
  _vtx = vtx;
}

Candidate::Candidate( const TVector3& mom,
		      const TParticlePDG& pdtEntry,
		      Vertex* vtx )
{
  init();
  setP3( mom ); 
  setPdtEntry( pdtEntry );
  _vtx = vtx;
}

Candidate::Candidate( const TVector2& tmom,
		      Vertex* vtx )
{
  init();
  _type = kTransverse;
  _pt   = tmom.Mod();
  _phi  = tmom.Phi();
  _eta  = 0.;
  _vtx  = vtx;
}

Candidate::Candidate( Vertex* vtx )  
{
  init();
  _vtx = vtx;
  assert( _vtx!=0 );
  
  for( size_t ii=0; ii<_vtx->nOutgoingCand(); ii++ )
    {
      addDaughter( _vtx->outgoingCand( ii ) );
    }
}

//Candidate::Candidate( const vector< const Candidate* >& listOfDau )
Candidate::Candidate( const CandList& listOfDau )
{
  init();
  for( size_t idau=0; idau<listOfDau.size(); idau++ )
    {
      const Candidate* dau = listOfDau[idau];
      if( dau->isTransverse() ) _type = kTransverse;
      addDaughter( dau->clone() );
    }
  // now loop on the daughters
  TVector2 p2_(0,0);
  float pz_(0);
  float E_(0);
  _q=0;
  Vertex* vtx_(0);
  bool sameVtx(false);
  for( size_t idau=0; idau<nDaughters(); idau++ )
    {
      const Candidate* dau = daughter(idau);
      if( vtx_==0 ) 
	{
	  sameVtx=true;
	  vtx_=dau->vertex();
	}
      else if( dau->vertex()!=vtx_ )
	{
	  sameVtx=false;
	}
      _q  += dau->charge(); 
      p2_ += dau->p2();
      if( isTransverse() )
	{
	  E_  += dau->Et();
	}
      else
	{
	  E_  += dau->E();
	  pz_ += dau->pz();
	}
    }
  _pt  = p2_.Mod(); 
  _phi = p2_.Phi();
  TVector3 p3_( p2_.X(), p2_.Y(), pz_ );
  _eta = 0;
  if( !isTransverse() && _pt>0 ) 
    { 
      _eta = p3_.Eta();
    }
  float m2_= pow(E_,2)-p3_.Mag2();
  if( m2_<0 ) m2_=0;
  _m = sqrt(m2_);
  if( sameVtx ) setVertex( vtx_ );
}

Candidate::Candidate( const Candidate& o )
  : 
  _uid(    o._uid  ),
  _type(   o._type ),
  _name(   o._name ),
  _q(      o._q    ),
  _m(      o._m    ),
  _pt(     o._pt   ),
  _eta(    o._eta  ),
  _phi(    o._phi  ),
  _pdtEntry(    o._pdtEntry  ),
  _vtx(    o._vtx  ),
  _info(   o._info )
{
  _status = kUnlocked;

  // make sure the original is in the bank of base candidates
  if( _baseCand.count(_uid)==0 ) _baseCand[_uid]=&o;

  for( size_t idau=0; idau<o.nDaughters(); idau++ )
    {
      addDaughter( o.daughter(idau)->clone() );
    } 
}

Candidate::~Candidate()
{
}

Candidate* 
Candidate::create()
{ 
  return new Candidate; 
}

Candidate* 
Candidate::create( const TVector3& mom,
		   float charge,
		   float mass,  
		   Vertex* vtx        )
{ 
  return new Candidate( mom, charge, mass, vtx ); 
} 

Candidate* 
Candidate::create( const TVector3& mom,
		   const TParticlePDG& pdtEntry,
		   Vertex* vtx )
{ 
  return new Candidate( mom, pdtEntry, vtx ); 
}

Candidate* 
Candidate::create( const TVector3& mom,
		   int pdgId,
		   Vertex* vtx )
{ 
  //FIXME MM
  if(pdgId==4414 || pdgId==4424 || pdgId==4434 || pdgId==5242)
    {pdgId=445;}
  // cout<<pdgId<<endl;
  const TParticlePDG* pdt_ 
    = ParticleData::Table()->GetParticle( pdgId );
  // assert( pdt_!=0 );
  //FIXME YHS
  if (pdt_ == 0)
	  cout << "pdgId of " << pdgId << " not recognized" << endl;
  return create( mom, *pdt_, vtx ); 
}

Candidate* 
Candidate::create( const TVector2& tmom,
		   Vertex* vtx )
{ 
  return new Candidate( tmom, vtx ); 
}
	     
Candidate* 
Candidate::create( Vertex* vtx )
{ 
  return new Candidate( vtx ); 
}
	     
Candidate*
//Candidate::create( const vector< const Candidate* >& listOfDau )
Candidate::create( const CandList& listOfDau )
{
  return new Candidate( listOfDau );
}

Candidate*
Candidate::create( const Candidate* dau1, const Candidate* dau2 )
{
  //  vector< const Candidate* > listOfDau;
  CandList listOfDau;
  listOfDau.push_back( (Candidate*)dau1 ); 
  listOfDau.push_back( (Candidate*)dau2 );
  return new Candidate( listOfDau );
}

Candidate* 
Candidate::clone() const
{ 
  return new Candidate(*this); 
}

Candidate* 
Candidate::transverseClone() const
{ 
  Candidate* out = new Candidate(*this);
  if( out->isTransverse() )  return out;
  out->_type = kTransverse;
  out->_eta  = 0;
  return out;
}

void
Candidate::setP3( const TVector3& mom )
{
  assert( !isLocked() );
  _pt  = mom.Pt();
  _eta = 1.e09;
  if( _pt!=0 ) _eta = mom.Eta();
  if( _eta>1.e09 )
    {
      mom.Print();
      //assert(0);
    }
  _phi = mom.Phi();
}


void
Candidate::setP2( const TVector2& mom )
{
  assert( !isLocked() );
  _pt  = mom.Mod();
  _phi = mom.Phi();
}


void
Candidate::setPdtEntry( const TParticlePDG& pdtEntry )
{
  assert( !isLocked() );
  _pdtEntry = &pdtEntry;
  _name = _pdtEntry->GetName();
  _q    = _pdtEntry->Charge();
  _m    = _pdtEntry->Mass();
}

void
Candidate::setPdgCode( int pdgId )
{
  const TParticlePDG* pdt_ 
    = ParticleData::Table()->GetParticle( pdgId );
  assert( pdt_!=0 );
  setPdtEntry( *pdt_ );
}    

void
Candidate::setPxPyPz( float px, float py, float pz )
{
  assert( !isLocked() );
  TVector3 mom( px, py, pz );
  setP3( mom );
}

void
Candidate::setPtEtaPhi( float pt, float eta, float phi )
{
  assert( !isLocked() );
  _pt  = pt;
  _eta = eta;
  _phi = phi;
}

void
Candidate::setMass( float m )
{
  assert( !isLocked() );
  _m = m;
}

void
Candidate::setVertex( Vertex* v )
{
  assert( !isLocked() );
  _vtx = v;
}

void
Candidate::init()
{
  _name = "";
  _status = kUnlocked;
  _uid  = ++_curUid;
  _type = kFull;
  _q    = 0.;
  _m    = 0.;
  _pt   = 0.;
  _eta  = 0.;
  _phi  = 0.;
  _pdtEntry  = 0;
  _vtx  = 0;

  _mother = 0;
  _info   = 0;
}

void
Candidate::lock()
{
  _status=kLocked;   
  assert( _baseCand.count(_uid)==0 );
  _baseCand[_uid]=this;
}

float          
Candidate::charge()               const
{
  return _q;
}

float          
Candidate::mass()                 const
{
  return _m;
}

TVector3       
Candidate::pos()                  const
{
  if( _vtx!=0 ) return _vtx->pos();
  return TVector3(0,0,0);
}

TVector3       
Candidate::p3()                   const
{
  TVector3 p3_;
  p3_.SetPtEtaPhi( _pt, _eta, _phi );
  return p3_;
}

TVector2       
Candidate::p2()                   const
{
  return TVector2( _pt*cos(_phi), _pt*sin(_phi) );
}

TLorentzVector 
Candidate::p4()                   const
{
  TLorentzVector p4_;
  p4_.SetPtEtaPhiM( _pt, _eta, _phi, _m );
  return p4_;
}

float          
Candidate::px()                    const
{
  return p3().X();
}

float          
Candidate::py()                    const
{
  return p3().Y();
}

float
Candidate::pz()                    const
{
  return p3().Z();
}

float          
Candidate::E()                    const
{
  return sqrt( pow(_m,2) + pow(p(),2) );
}

float          
Candidate::Et()                    const
{
  return sqrt( pow(_m,2) + pow(_pt,2) );
}

float          
Candidate::p()                   const
{
  return p3().Mag();
}

float          
Candidate::pt()                   const
{
  return _pt;
}

float          
Candidate::eta()                  const
{
  return _eta;
}

float          
Candidate::phi( int unit )   const
{
  float phi_ = _phi;
  if( phi_<0     ) phi_ += 2*Constants::pi;
  if( unit==kDeg ) phi_ *= Constants::radToDeg;
  return phi_;
}

float          
Candidate::sinPhi()               const
{
  return sin(_phi);
}

float          
Candidate::cosPhi()               const
{
  return cos(_phi);
}

float          
Candidate::tanPhi()               const
{
  return tan(_phi);
}

float          
Candidate::theta( int unit ) const
{
  float th_ = p3().Theta();
  if( unit==kDeg ) th_ *= Constants::radToDeg;
  return th_;
}

float          
Candidate::sinTheta()               const
{
  return sin(theta());
}

float          
Candidate::cosTheta()               const
{
  return cos(theta());
}

float          
Candidate::tanTheta()               const
{
  return tan(theta());
}

float          
Candidate::tanDip()               const
{
  return tan(Constants::pi/2.-theta());
}

Vertex*          
Candidate::vertex()               const
{
  return _vtx;
}

const TParticlePDG* 
Candidate::pdtEntry() const
{
  return _pdtEntry;
}

int
Candidate::pdgCode() const
{
  if( _pdtEntry==0 ) return 0;
  return _pdtEntry->PdgCode();
}

const string &
Candidate::name() const
{
  return _name;
}

void
Candidate::print( ostream& o ) const
{
  o << "\n";
  oneLine(o);
  o << "Candidate -- ";
  o << name();
  //  o << "\tuid=" << uid() << "\tptr=" << this << endl;
  o << "-- (Px=" << px() << ", Py=" << py() << ", Pz=" << pz() << "; E=" << E() <<  ") q=" << charge() << " M=" << mass() << " Gev/c2 ";    
  o << endl;
  // print daughter links
  if( _vtx ) { o << " vtx= "; _vtx->print(o); }
  o << "pt="  << pt() << " GeV, ";
  o << "\teta=" << eta() << ", ";
  o << "\tphi=" << phi( kDeg ) << " deg";
  o << endl;
  if( isComposite() )   
    {
      o << "M="   << mass() << " GeV, ";
      o << "E="   << E()    << " GeV ";
      o << "d: "<< nDaughters();
      for( size_t ii=0; ii<nDaughters(); ii++ ) 
	{
	o << " - " << daughter(ii)->uid();
	o << endl;
	daughter(ii)->oneLine( o );
      }
      o << endl;
    }						
  //  o << "mother=" << theMother() << endl;
  if( _info!=0 ) _info->print( o );
}

void
Candidate::oneLine( ostream& o, bool ltx ) const
{
  const CandInfo* info_ = info();

  char line_[512];
  //  o << _pdtEntry->PdgCode() << " ";
  //  o << name();
  //  o << " UID=" << uid(); 
  //  o << " pt="  << pt() << " GeV,";
  //  o << " eta=" << eta() << ",";
  //  o << " phi=" << phi( kDeg ) << " deg, ";
  if( ltx )
    {
      string str_;
      if( name()=="e-" )  str_="e^{-}"; 
      if( name()=="e+" )  str_="e^{+}"; 
      if( name()=="mu-" ) str_="#mu^{-}"; 
      if( name()=="mu+" ) str_="#mu^{+}"; 
      sprintf( line_, "%s : P_{T}=%7.2f GeV / #eta=%5.2f / #phi=%6.1f deg ",str_.c_str(),pt(),eta(), phi(kDeg) );
      o << line_;
      return;
    }
  else
    {
      sprintf( line_, "%s : P_T=%7.2f GeV / eta=%5.2f / phi=%6.1f deg ",name().c_str(),pt(),eta(), phi(kDeg) );
      o << line_;
    }
  if( info_!=0 ) 
    {
      int elID_(0);

      bool ok_;
      if( info_->getBool( "eidLoose", ok_ ) )
	o << (( ok_ ) ? "1":"0");
      else o << "x";
      if( info_->getBool( "eidRobustHighEnergy", ok_ ) )
	o << (( ok_ ) ? "1":"0");
      else o << "x";
      if( info_->getBool( "eidRobustLoose", ok_ ) )
	o << (( ok_ ) ? "1":"0");
      else o << "x";
      if( info_->getBool( "eidRobustTight", ok_ ) )
	o << (( ok_ ) ? "1":"0");
      else o << "x";
      if( info_->getBool( "eidTight", ok_ ) )
	o << (( ok_ ) ? "1":"0");
      else o << "x";
      o << " ";

      if( info_->getInt( "elID", elID_ ) )
	o << " ID=" << elID_ << ",";
      float iso_;
      if( info_->getFloat( "iso", iso_ ) )
	o << " iso=" << iso_ << ",";
      bool b_(false);
      if(      info_->getBool("VeryTight", b_ ) && b_ ) o << " VT";
      else if( info_->getBool("Tight", b_ ) && b_ )     o << " T";
      else if( info_->getBool("Loose", b_ ) && b_ )     o << " L";
      else if( info_->getBool("VeryLoose", b_ ) && b_ ) o << " VL";
    }
  if( isComposite() )   
    {
      sprintf( line_, " M=%7.2f GeV ", mass() );
      //      o << "\tM="   << mass() << " GeV" << endl;
      o << line_;
    }

  int vtxid = vertexIndex();
  o << " -- Vtx=" << vtxid;
  if( !isAtPrimaryVertex() )
    {
      if( vtxid!=-1 ) o << " (not PV)";
    }
  else
    {
      o << " (PV)";
    }
  o << endl;
  if( isComposite() )   
    {
      //      o << "\tM="   << mass() << " GeV" << endl;
      for( size_t idau=0; idau<nDaughters(); idau++ )
	{
	  const Candidate* dau_ = daughter(idau);
	  o << "\t--> ";
	  dau_->oneLine( o );
	}
    }
}
      
bool  
Candidate::operator==( const Candidate& o ) const
{
  return _uid == o._uid;
}

size_t
Candidate::nDaughters() const
{
  return _daughter.size();
}

bool 
Candidate::isComposite() const
{
  return nDaughters()>0;
}

const Candidate* 
Candidate::theBase() const
{
  //  if( isComposite() ) return this;
  //  if( isLocked()    ) return this;
   const Candidate* cand = _baseCand[_uid];
   if( cand==0 ) cand = this;
   return cand;
} 
Candidate* 
Candidate::theBase()
{
  //  if( isComposite() ) return this;
  //  if( isLocked()    ) return this;
  Candidate* cand =  const_cast<Candidate*>(_baseCand[_uid]);
  if( cand==0 ) cand = this;
  return cand;
} 

const Candidate* 
Candidate::daughter( size_t idau ) const
{
  assert( idau<nDaughters() );
  return _daughter[idau];
}

Candidate* 
Candidate::daughter( size_t idau )
{
  assert( idau<nDaughters() );
  return _daughter[idau];
}

void 
Candidate::setMother( Candidate* mo )
{
  assert( !isLocked() );
  assert( mo!=0 );
  _mother = mo;
  //  mo->addDaughter(this);
}

void 
Candidate::addDaughter( Candidate* dau )
{
  assert( !isLocked() );  
  assert( dau!=0 );
  _daughter.push_back( dau );
  dau->setMother(this);
}

void
Candidate::setInfo( CandInfo* info )
{
  assert( _info==0 );
  _info = info;
}

void
Candidate::basePrint( ostream& o )
{
  for( map<size_t,const Candidate*>::const_iterator it=_baseCand.begin();
       it!=_baseCand.end(); ++it )
    {
      it->second->print(o);
    }
}

const Candidate* 
Candidate::base( size_t uid )
{
  map<size_t,const Candidate*>::const_iterator it = _baseCand.find(uid);
  if( it==_baseCand.end() ) return 0;
  return it->second; 
}



float
Candidate::dz( const Vertex* pv ) const
{

  float dz_=0;
  
  float vx_= _vtx->pos().X();
  float vy_= _vtx->pos().Y();
  float vz_= _vtx->pos().Z();

  float pvx_= pv->pos().X();
  float pvy_= pv->pos().Y();
  float pvz_= pv->pos().Z();

  
  dz_ = (vz_-pvz_) - ((vx_-pvx_)*px()+(vy_-pvy_)*py())/pt() * pz()/pt(); 


  return dz_;

}


float
Candidate::d0( const Vertex* pv ) const
{

  float d0_=0;

  float vx_= _vtx->pos().X();
  float vy_= _vtx->pos().Y();

  float pvx_= pv->pos().X();
  float pvy_= pv->pos().Y();

  d0_ = ( -(vx_-pvx_)*py() + (vy_-pvy_)*px() ) / pt();

  return d0_;

}


float
Candidate::dR( const Candidate* cand) const
{

  return KineUtils::dR( eta(),  cand->eta(), phi() , cand->phi() );

}

float
Candidate::d3r( const Vertex* pv ) const
{
  
  float vx_= _vtx->pos().X();
  float vy_= _vtx->pos().Y();
  float vz_= _vtx->pos().Z();

  float pvx_= pv->pos().X();
  float pvy_= pv->pos().Y();
  float pvz_= pv->pos().Z();

  return (vz_-pvz_)*pt()/p() - ((vx_-pvx_)*px()+(vy_-pvy_)*py())/pt() *pz()/p(); 

}

bool
Candidate::isAtPrimaryVertex() const
{
  Vertex* vtx_ = vertex();
  if( vtx_==0 ) return false;
  return vtx_->isPrimary();
}

int
Candidate::vertexIndex() const
{
  int vtxid =-1;
  Vertex* vtx_=vertex();
  if( vtx_!=0 )
    {
      vtxid = vtx_->index();
    }
  return vtxid;
}
