#ifndef GhmTupleAnalysis_hh
#define GhmTupleAnalysis_hh

#include "Analysis/selectors/SelectorBase.hh"
#include "Analysis/ghmzz/BaseTupleAnalysis.hh"


// WZ candidates
class VVCandidate
{
public:
  VVCandidate() {}
  virtual ~VVCandidate() {}
  virtual void print() const {}
};

class WZCandidate : public VVCandidate
{
  WZCandidate( Candidate* Zboson, Candidate* lepton, Candidate* etmiss );
  void print() const;
};

// W candidate
class WCandidate
{
public:

  WCandidate( Candidate* lepton, Candidate* etmiss );
  ~WCandidate() {}

  float mT() const { return _mT; }
  Candidate* l()   const { return _l;  }
  Candidate* n()  const { return _n; }
  Candidate* ln() const { return _ln; }
  Candidate* W() const { if( n() ) return _Wsol[0]; return 0; }
  const CandList& Wsol() const { return _Wsol; }

private:

  // the leptons
  Candidate* _l;
  Candidate* _n;

  // the transverse W candidate
  Candidate* _ln;
    
  // the sorted list of W candidates
  CandList _Wsol;

  TVector2 _pT_l_vec;
  TVector2 _pT_n_vec;
  float _m, _mT;
  float _m2, _mT2;
  float _p_l, _pT_l, _pL_l; 
  float _pT_n;
  float _A, _B; 
};

class GhmTupleAnalysis : public BaseTupleAnalysis
{
public:

  GhmTupleAnalysis( const char* sampleSet=0, int n=-1 );
  virtual ~GhmTupleAnalysis() {}

  // selection & sorting
  class ZSelector : public SelectorBase
  {
  public:
    bool accept( const Candidate& cand ) const;
  };

  class ThirdLeptonSelector : public SelectorBase
  {
  public:
    bool accept( const Candidate& cand ) const;
  };

  struct ZCompare
  {
    bool operator()( const Candidate* Z1, const Candidate* Z2 ) const;
  };

protected:
			 
  virtual void init();
  virtual bool analyzeEvent();
  virtual void finalize();

  static std::string indexedStr( std::string, size_t ii ); 

  ClassDef( GhmTupleAnalysis, 0 )  
};

#endif
