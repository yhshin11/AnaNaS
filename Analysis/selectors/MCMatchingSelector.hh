#ifndef MCMatchingSelector_h
#define MCMatchingSelector_h

#include <vector>
#include <map>

#include "Analysis/core/Candidate.hh"

#include "Analysis/utils/KineUtils.hh"


class MCMatchingSelector
{

public:

  //Contructors and destructors

  MCMatchingSelector();
  MCMatchingSelector(int mcStatus, int pdgId, 
		     bool chargeMatching,
		     float dRcut, float dPtcut);
  MCMatchingSelector(std::vector<int> mcStatus, std::vector<int> pdgId, 
		     bool chargeMatching,
		     float dRcut, float dPtcut);
  virtual ~MCMatchingSelector() {};

  //Predefined official values
  void loadMCMatchPhoton();
  void loadMCMatchElectron();
  void loadMCMatchMuon();

  Candidate* MCMatch(Candidate* cand, const CandList& );
  bool MCMatchTruth(Candidate* cand, const CandList& );
  std::map<Candidate*, Candidate*> MCMatch(const CandList&, const CandList&);

private:

  std::vector<int> _mcStatus;
  std::vector<int> _pdgIds;
  
  bool _chargeMatch;

  float _dRcut;
  float _dPtcut;

  //Private functions
  bool Preselection(float charge, Candidate* mcCand);
  Candidate* FindMcCandidate(Candidate* cand, const CandList&);
  

  ClassDef( MCMatchingSelector, 0)
};



#endif
