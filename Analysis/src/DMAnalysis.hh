#ifndef DMAnalysis_h
#define DMAnalysis_h

#include <vector>
#include <map>
#include <string>

class TH1;
class TH2;

#include "Analysis/core/SampleAnalysis.hh"

#include "Analysis/selectors/MuonSelector.hh"
#include "Analysis/selectors/ElectronSelector.hh"

// Note: kLeptonic excludes tau+nu decay mode.
namespace WDecayChannel {enum decayChannel { kInternal=-100, kUnclassified=-1, kLeptonic=0, kTauUnclassified=1, kTauLeptonic=2, kTauHadronic=3, kHadronic=10 };}
namespace TauDecayChannel {enum decayChannel { kInternal=-100, kUnclassified=-1, kLeptonic=0, kHadronic=10 };}
 
class DMAnalysis : public SampleAnalysis
{
public:

  //  DMAnalysis( Sample& sample, const std::string & collectionFileName );
  DMAnalysis( Sample& sample, EventManager& manager );
  virtual ~DMAnalysis();

public:

  virtual void bookHistograms();
  virtual bool analyzeEvent();
  virtual void writeHistograms();

private:
  // insert private members here
  //XXX TH1* hist_mll;
  //XXX double mll_tree;
  float getRapidity(const Candidate* zCandidate);
	bool checkVertices();
	bool isLepton(const Candidate* jet);
	void fillGenEvent();
	Candidate* findAncestorOfMcCandidate(Candidate* mcCand, const CandList& mcCands, CandList& motherCands);
	Candidate* findMotherOfMcCandidate(Candidate* mcCand, const CandList& mcCands);
	bool isInList(const Candidate* cand, const CandList& candlist);
	void fillWChan();
	WDecayChannel::decayChannel classifyWDaughters(Candidate*);
	TauDecayChannel::decayChannel classifyTauDaughters(Candidate*);
	Candidate* createGenCand(Candidate*);
	void findDaughters(Candidate*, CandList&);

  MuonSelector MuonSel;
  ElectronSelector ElectronSel;

  bool HLTFiredMultiChan();

  //For matching
  typedef std::map< Candidate* , Candidate* > Matching;
  typedef Matching::const_iterator MatchingIterator;

  //MCTruth
  Matching MCToGen;
  Matching GenToMC;

  int _zChan; // 0 for Z->mumu; 1 for Z->ee; 2 for Z->mue

	// Variables to be filled.
	int Nvertex; // Not filled right now.

  ClassDef( DMAnalysis, 0 )  
};

#endif

     
