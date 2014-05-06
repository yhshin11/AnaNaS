#include "Analysis/tools/JetSmearing.hh"

using namespace std;



JetSmearing::JetSmearing(EventManager* em) {
  
  LoadJetSmearingDB();
  LoadJetResolutionDB();

  smearBy_ = ;
  shiftBy_ = ;
  _e = em;
}


JetSmearing::~JetSmearing() {
}


void JetSmearing::LoadJetSmearingDB() {

  string db="";


}


void JetSmearing::LoadJetResolutionDB() {

  string db="";


}


float JetSmearing::jetResolutionExtractor(const Candidate* jet) {
  
  float res;


  return res;
}


Candidate* JetSmearing::smearing(jetCol jets) {

  // jet : corrected or raw jet entered by the user
  // rawjet : matched raw jet...

  const Candidate* rawJet=jets.rawJet;
  const Candidate* genJet=jets.genJet;
  const Candidate* jet=jets.corJet;
  const Candidate* smJet=jets.smearJet;
  
  double smearFactor = 1.;   
  double newSF = 1.;   
  double x = TMath::Abs(jet.eta());
  double y = jet.pt();
  if ( x > lut_->GetXaxis()->GetXmin() && x < lut_->GetXaxis()->GetXmax() &&
       y > lut_->GetYaxis()->GetXmin() && y < lut_->GetYaxis()->GetXmax() ) {
    int binIndex = lut_->FindBin(x, y);
	
    if ( smearBy_ > 0. ) 
      newSF += smearBy_*(lut_->GetBinContent(binIndex) - 1.);
      smearFactor += smearBy_*(lut_->GetBinContent(binIndex) - 1.);
    double smearFactorErr = lut_->GetBinError(binIndex);
   
    if ( shiftBy_ != 0. ) {
      smearFactor += (shiftBy_*smearFactorErr);
      newSF+= (shiftBy_*smearFactorErr);
    }
  }


  // New Smearing =============================================================== 

  float smearedjetEn = smJet.E()/jet.E();

  float Max=max(rawJet.E(),jet.E() );

  if(genJet==NULL) {

    if( smearFactor> 1.) {
    
      float gausSmear = (smearedjetEn/rawJet.E() - 1)*Max;
    
      float fPol = smearFactor*smearFactor -2*smearFactor + 1;
      float fPolp = newSF*newSF -2*newSF + 1;
      float newGausS = exp(log(gausSmear)*(fPol/fPolp));
    
      newSmearEn = rawJet.E()*(1+newGausS)/Max;
    }
  
  }
  else {
 
    if(smJet.E()==jet.E()) continue; //no smearing because too far from genJet

    double dEn = jet->energy() - genJet->energy();
    
    float newSmearEn = rawJet.energy()*(1. + (newSF - 1.)*dEn/TMath::Max(rawJet.E(), jet.E()));


  }
  //============================================================================
 

  // CV: keep minimum jet energy, in order not to loose direction information
  const double minJetEn = 1.e-2;
  if ( newSmearEn < minJetEn ) newSmearEn = minJetEn;

  // CV: skip smearing in case either "raw" or "corrected" jet energy is very low
  //     or jet passes selection configurable via python
  //    (allows for protection against "pathological cases",
  //     cf. PhysicsTools/PatUtils/python/tools/metUncertaintyTools.py)
  TLorentzVector smearedJetP4 = jet->p4();
  if ( !( abs(jet->energy() - rawJet->energy() ) > (5.*min(jet->energy(), rawJet->energy() )) 
	   //(skipJetSelection_ && (*skipJetSelection_)(jet)) || 
	 rawJet.pt()  < skipRawJetPtThreshold_          ||
	 jet.pt() < skipCorrJetPtThreshold_         ) ) {
    smearedJetP4 *= (newSmearEn/jet.energy());
  }

      
  Candidate* smearJet = Candidate::create( smearedJetP4->P3(), jet->charge(), 
					   smearedJetP4->Mass(), jet->vertex() );
      
  return smearJet;

}


const Candidate* 
JetSmearing::genJetMatch(const Candidate* jet) {

  const CandList& mcJets = _e->jetList( EventManager::kGenJet);
  
  bool matched=false;
  for(size_t ij=0;ij<mcJets.size();ij++) {
    
    const Candidate* mcjet = mcJets[ij];    
    if( jet->dR( mcjet ) < min(0.5, 0.1 + 0.3*exp(-0.05*(mcJet->pt() - 10.))) )
      return mcJet;
  }

  return NULL;

}


const Candidate* 
JetSmearing::rawJetMatch(const Candidate* jet) {

  const CandList& rawJets = _e->jetList( EventManager::kPfJet, "ak5PFJets" ); //default rawPfJets
  
  bool matched=false;
  for(size_t ij=0;ij<rawJets.size();ij++) {
    
    const Candidate* rawjet = rawJets[ij];    
    if( jet->dR( rawjet ) < 0.15 )
      return rawjet;
  }

  return jet;

}

const Candidate* 
JetSmearing::rawJetMatch(const Candidate* jet) {

  const CandList& rawJets = _e->jetList( EventManager::kPfJet, "smearedAK5PFJets" ); //default rawPfJets
  
  bool matched=false;
  for(size_t ij=0;ij<rawJets.size();ij++) {
    
    const Candidate* rawjet = rawJets[ij];    
    if( jet->dR( rawjet ) < 0.15 )
      return rawjet;
  }

  return jet;

}
