#include "Analysis/selectors/MCMatchingSelector.hh"

#include <stdlib.h>

using namespace std;


ClassImp( MCMatchingSelector )

//Constructors

MCMatchingSelector::MCMatchingSelector() {

  //By default, no mcStatus requirement, no charge matching, no pdgId matching
  _chargeMatch = false;

  //By default, dR=0.5, dPt=0.5
  _dRcut  = 0.5;
  _dPtcut = 0.5;

}


MCMatchingSelector::MCMatchingSelector(int mcStatus, int pdgId, 
				       bool chargeMatching,
				       float dRcut, float dPtcut) {

  _mcStatus.push_back(mcStatus);
  _pdgIds.push_back( abs(pdgId) );
  
  _chargeMatch = chargeMatching;

  _dRcut  = dRcut;
  _dPtcut = dPtcut;

}


MCMatchingSelector::MCMatchingSelector(std::vector<int> mcStatus,
				       std::vector<int> pdgId, 
				       bool chargeMatching,
				       float dRcut, float dPtcut) {
  
  _mcStatus = mcStatus;
  _pdgIds = pdgId;
  
  _chargeMatch = chargeMatching;
  
  _dRcut  = dRcut;
  _dPtcut = dPtcut;
  
}


//Predefined official values
void
MCMatchingSelector::loadMCMatchPhoton() {

  _mcStatus.push_back(1);
  _pdgIds.push_back(22);

  _chargeMatch = true;

  _dRcut  = 0.2;
//   _dPtcut = 1.0;
  _dPtcut = 0.1;

}

void
MCMatchingSelector::loadMCMatchElectron() {

  _mcStatus.push_back(1);
  _pdgIds.push_back(11);

  _chargeMatch = true;

  _dRcut  = 0.5;
  _dPtcut = 0.5;

}

void
MCMatchingSelector::loadMCMatchMuon() {

  _mcStatus.push_back(1);
  _pdgIds.push_back(13);

  _chargeMatch = true;

  _dRcut  = 0.5;
  _dPtcut = 0.5;
  
}

//Matching functions

Candidate*
MCMatchingSelector::MCMatch(Candidate* cand, const CandList& mcCandList) {
  
  return FindMcCandidate(cand, mcCandList);
}

bool
MCMatchingSelector::MCMatchTruth(Candidate* cand, const CandList& mcCandList) {
  
  if( FindMcCandidate(cand, mcCandList) ) return true;
  else                                    return false;
  
}

map<Candidate*, Candidate*>
MCMatchingSelector::MCMatch(const CandList& candlist, const CandList& mcCandList) {

  map<Candidate*, Candidate*> mcMatchMap;

  for(unsigned int ic=0;ic<candlist.size();ic++) {
    
    Candidate* cand = candlist[ic];
    Candidate* mcCand = FindMcCandidate(cand, mcCandList);
    
    mcMatchMap[ cand ] = mcCand;
  }

  return mcMatchMap;

}

//Private functions

bool
MCMatchingSelector::Preselection(float charge, Candidate* mcCand) {

  bool checkStatus=false;
  bool pdgId=false;
  bool checkCharge =false;

  CandInfo* info = mcCand->info();
  int stat = info->getInt("status") ;
  checkStatus = (stat==1) ;
  
  if(_pdgIds.size() == 0 ) pdgId = true; 
  else {
    for (unsigned int ii=0;ii<_pdgIds.size();ii++) {
      if( abs(mcCand->pdgCode()) == _pdgIds[ii] ) {
	  pdgId = true; break;
      }
    }
  }
  
  if(_chargeMatch && (charge==mcCand->charge() ) )
    checkCharge = true;
  else if(!_chargeMatch)
    checkCharge = true;

  
  if(checkStatus && pdgId && checkCharge)
    return true;
  else
    return false;
  
}


Candidate*
MCMatchingSelector::FindMcCandidate(Candidate* cand, const CandList& mcCandList)  {

  float dR  = 10;
  float dPt = 10;
  float dRtmp = 10;
  Candidate* BestMcCand = 0;

  for(unsigned int imc=0;imc<mcCandList.size();imc++) {
    
    Candidate* mcCand = mcCandList[imc];
    if(!Preselection( cand->charge(), mcCand) ) continue;

    dR  = KineUtils::dR(cand->eta(), mcCand->eta(), cand->phi(), mcCand->phi() );
    dPt = fabs( cand->pt() - mcCand->pt() ) / cand->pt() ;

    if(dR>_dRcut)   continue;
    if(dPt>_dPtcut) continue;
    
    if(dR<dRtmp) {
      dRtmp = dR;
      BestMcCand = mcCand;
    }
  }
  return BestMcCand;
}




