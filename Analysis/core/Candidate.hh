#ifndef Candidate_hh
#define Candidate_hh

#include "Analysis/utils/ObjectStore.hh"
#include "Analysis/utils/Counter.hh"

#include <vector>

#include <TLorentzVector.h>
#include <TVector3.h>
#include <TVector2.h>

#include "Analysis/utils/ParticleData.hh"
#include "Analysis/core/Vertex.hh"
#include "Analysis/core/CandInfo.hh"

#include "Analysis/core/CandFwd.hh"

class Candidate : public ObjectStore, public Counter<Candidate>
{
public:

  // states
  enum Type   { kFull=0,     kTransverse };
  enum Status { kUnlocked=0, kLocked     };

  // default construction
  static Candidate* create();

  // construction from momentum, charge and mass
  static Candidate* create( const TVector3& mom,
			    float charge=0,
			    float mass=0,  
			    Vertex* vtx=0         );

  // construction from momentum and pggId
  static Candidate* create( const TVector3& mom,
			    const TParticlePDG& pdtEntry,
			    Vertex* vtx=0 );

  // construction from momentum and pggId
  static Candidate* create( const TVector3& mom,
			    int pdgId,
			    Vertex* vtx=0 );

  // construction for transverse candidates
  static Candidate* create( const TVector2& tmom,
			    Vertex* vtx=0 );
	     
  // construction for composite candidates 
  //  static Candidate* create( const std::vector< const Candidate* >& listOfDau );
  static Candidate* create( const CandList& listOfDau );
  static Candidate* create( const Candidate* dau1, const Candidate* dau2 );
	     
  // cloning
  Candidate* clone() const;

  // get a transverse clone 
  Candidate* transverseClone() const;

  // return the base candidate
  const Candidate* theBase() const; 
  Candidate* theBase(); 

  // accessors
  const TParticlePDG*  pdtEntry()        const;
  int                  pdgCode()         const;
  const std::string &  name()            const;
  float                charge()          const;
  float                mass()            const;
  TVector3             pos()             const;
  TVector2             p2()              const;
  TVector3             p3()              const;
  TLorentzVector       p4()              const;
  float                px()              const;
  float                py()              const;
  float                pz()              const;
  float                E()               const;
  float                Et()              const;
  float                p()               const;
  float                pt()              const;
  float                eta()             const;
  float                phi(int u=kRad)   const;
  float                sinPhi()          const;
  float                cosPhi()          const;
  float                tanPhi()          const;
  float                theta(int u=kRad) const;
  float                sinTheta()        const;
  float                cosTheta()        const;
  float                tanTheta()        const;
  float                tanDip()          const;
  Vertex*              vertex()          const;
  int                  vertexIndex()     const;
  bool                 isAtPrimaryVertex() const;

  //Impacts parameter w.r.t one vertex
  float                d0( const Vertex* pv )    const;
  float                dz( const Vertex* pv )    const;
  float                d3r( const Vertex* pv )   const;

  //Useful dR function
  float                 dR( const Candidate* cand ) const; 

  // state
  size_t  uid()       const { return _uid; }
  Status  status()    const { return _status; }
  Type    type()      const { return _type;   }
  bool isTransverse() const { return type()==kTransverse; }
  bool isLocked(    ) const { return status()==kLocked;  }

  // extra information
  CandInfo* info() const { return _info; }
  CandInfo* info()       { return _info; }

  // navigate throught decay tree
  bool             isComposite() const;
  size_t           nDaughters()  const;
  const CandList&  theDaughters() const { return _daughter; }
  const Candidate* daughter( size_t idau ) const;
  Candidate*       daughter( size_t idau );
  const Candidate* theMother()   const { return _mother; }
  Candidate*       theMother() { return _mother; }

  // modifiers
  void setName( std::string name ) { _name=name; }
  void setP3(       const TVector3& mom            );
  void setP2(       const TVector2& mom            );
  void setPxPyPz(   float px, float py, float pz   );
  void setPtEtaPhi( float pt, float eta, float phi );
  void setPdtEntry( const TParticlePDG& pdtEntry   );
  void setPdgCode(  int pdgCode   );
  void setMass(     float m  );
  void setVertex(   Vertex* v );

  // assign info block
  void setInfo( CandInfo* info );

  // set mother link
  void setMother( Candidate* );

  // set daughter link
  void addDaughter( Candidate* );

  // lock candidate -- cannot be unlocked 
  void lock();
      
  // equality operator (based on uids)
  bool  operator==( const Candidate& )  const;
  
  // resets uids
  static void reset();

  // prints
  void print(   ostream& o ) const;
  void oneLine( ostream& o, bool ltx=false ) const;
  
  // global prints
  static void basePrint( ostream& o );
  static const Candidate* base( size_t uid );

  // sorting
  enum Sorting { kSortDefault=0, kSortByPt, kSortByEnergy, kSortByEta }; 
  
  struct SortByPt
  {
    bool operator()( const Candidate* c1, const Candidate* c2 ) const
    {
      float pt1 = (c1!=0) ? c1->pt() : 0;
      float pt2 = (c2!=0) ? c2->pt() : 0;
      return pt1>pt2;
    }
  };
  
  struct SortByEnergy
  {
    bool operator()( const Candidate* c1, const Candidate* c2 ) const
    {
      float e1 = (c1!=0) ? c1->E() : 0;
      float e2 = (c2!=0) ? c2->E() : 0;
      return e1>e2;
    }
  };

  struct SortByEta
  {
    bool operator()( const Candidate* c1, const Candidate* c2 ) const
    {
      float eta1 = (c1!=0) ? c1->eta() : 10;
      float eta2 = (c2!=0) ? c2->eta() : 10;
      return eta1<eta2;  // warning, from minus to plus
    }
  };

  
private:

  // constructors and destructor are private 
  // -- use "create" functions

  // default constructor
  Candidate();

  // constructor from momentum, charge and mass
  Candidate( const TVector3& mom,
	     float charge,
	     float mass,  
	     Vertex* vtx         );

  // constructor from momentum and pggId
  Candidate( const TVector3& mom,
	     const TParticlePDG& pdtEntry,
	     Vertex* vtx );

  // constructor for transverse candidates
  Candidate( const TVector2& tmom,
	     Vertex* vtx );
	     
  // constructor from vertex
  Candidate( Vertex* vtx );
	     
  // construction from a list of candidates
  //  Candidate( const std::vector< const Candidate* >& listOfDau );
  Candidate( const CandList& listOfDau );

  // copy constructor
  Candidate( const Candidate& );
  
  // construction for a vertex
  static Candidate* create( Vertex* vtx );
	     
  // destructor
  ~Candidate();
  
private:

  // must be called by all ctor, except copy ctor
  void init();
  
  // identifier; for clones, identifier of the original
  size_t _uid; 

  // state
  Status   _status;
  Type     _type;

  // data 
  std::string   _name;
  float            _q;
  float            _m;
  float           _pt;
  float          _eta;
  float          _phi;
  const TParticlePDG* _pdtEntry;
  Vertex*       _vtx;

  // info block
  CandInfo* _info;

  // family
  Candidate* _mother;
  CandList _daughter;

  // candidate bookkeeping
  static size_t _curUid;
  static std::map< size_t, const Candidate* > _baseCand; 
      
  ClassDef(Candidate,0) 
};

#endif




