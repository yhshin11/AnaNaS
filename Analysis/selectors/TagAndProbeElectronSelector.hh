#ifndef TAPESelector_h
#define TAPESelector_h

#include "Analysis/selectors/SelectorBase.hh"
#include "Analysis/selectors/VBTFElectronSelector.hh"
#include "Analysis/core/EventManager.hh"

#include <vector>
#include <string>
#include <map>

//Maps containing cuts to determine Tags and Probes
typedef std::vector< std::pair< std::string, std::vector<std::string> > > CutMap;
//typedef CutMap::const_iterator CutMapIterator;


class TagAndProbeElectronSelector : public SelectorBase     
{

public:
 //  TagAndProbeElectronSelector();
  TagAndProbeElectronSelector( unsigned settings , CutMap cMap);
 //  TagAndProbeElectronSelector(const CandList& cand, CutMap cMap);
  virtual ~TagAndProbeElectronSelector() {}

  
  // return true if the candidate is selected 
  virtual bool accept( const Candidate& cand );

  // sorting
  //virtual void sortList( CandList& list, int sorting ) const;

  

public:

  //Map
  CutMap _cMap;
  unsigned _settings;

  //All possibles cuts
  bool PtEtaCut(const Candidate& cand, float ptCut, std::vector<float> etaCut);
  bool HLTCut(const Candidate& cand, std::vector<std::string> Lines) ; 
  bool IsoCut(const Candidate& cand, std::vector<std::string> Isos ) ; 
  bool IDCut(const Candidate& cand, std::vector<std::string> IDs); 
  bool MatchingDR(const Candidate& cand1, const Candidate& cand2, float cut, std::string cuttype );
  bool GeneralCut(const Candidate& cand, std::string typevar, float val,
		  std::string typvar, std::string cuttype );
  bool QualityCutByGreaterVal( float val, float cut)
  { return  val > cut;  };
  bool QualityCutByLowerVal( float val, float cut)
  { return  val < cut;  };
  bool QualityCutByEgalVal( float val, float cut)
  { return  val == cut;  };
  bool QualityCutByNotEgalVal( float val, float cut)
  { return  val != cut;  };
  
  //Match the cut
  bool PassByType(const Candidate& cand, std::string type,
		  std::vector<std::string> vals,  Candidate* addc=NULL);
  //std::vector<std::string>

  void LoadVBTFSel(const Candidate& cand);

private:

  bool _debug;

  std::vector<float> _FiduEtaCuts;
  VBTFElectronSelector VBTFSel;

  ClassDef( TagAndProbeElectronSelector, 0 ) 
};

#endif
