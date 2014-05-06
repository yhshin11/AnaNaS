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
  
  double smearFactor = 1.;      
  double x = TMath::Abs(jet.eta());
  double y = jet.pt();
  if ( x > lut_->GetXaxis()->GetXmin() && x < lut_->GetXaxis()->GetXmax() &&
       y > lut_->GetYaxis()->GetXmin() && y < lut_->GetYaxis()->GetXmax() ) {
    int binIndex = lut_->FindBin(x, y);
	
    if ( smearBy_ > 0. ) smearFactor += smearBy_*(lut_->GetBinContent(binIndex) - 1.);
    double smearFactorErr = lut_->GetBinError(binIndex);
   
    if ( shiftBy_ != 0. ) {
      smearFactor += (shiftBy_*smearFactorErr);
    }
  }

  double smearedJetEn = jet.energy();
  double sigmaEn = jetResolutionExtractor(rawJet)*TMath::Sqrt(smearFactor*smearFactor - 1.);
  bool isGenMatched = false;
  if ( genJet ) {
    double dEn = jet->energy() - genJet->energy();
    if ( dEn < (sigmaMaxGenJetMatch_*sigmaEn) ) {
      //--- case 1: reconstructed jet matched to generator level jet, 
      //            smear difference between reconstructed and "true" jet energy

      smearedJetEn = rawJet.energy()*(1. + (smearFactor - 1.)*dEn/TMath::Max(rawJet.E(), jet.E()));
      isGenMatched = true;
    }
  }
  if ( !isGenMatched ) {
    //--- case 2: reconstructed jet **not** matched to generator level jet, 
    //            smear jet energy using MC resolution functions implemented in PFMEt significance algorithm (CMS AN-10/400)

    if ( smearFactor > 1. ) {
      // CV: MC resolution already accounted for in reconstructed jet,
      //     add additional Gaussian smearing of width = sqrt(smearFactor^2 - 1) 
      //     to account for Data/MC **difference** in jet resolutions.
      //     Take maximum(rawJetEn, corrJetEn) to avoid pathological cases
      //    (e.g. corrJetEn << rawJetEn, due to L1Fastjet corrections)

      smearedJetEn = rawJet.energy()*(1. + rnd_.Gaus(0., sigmaEn)/TMath::Max(rawJet.E(), jet.E()));
    }
  }

  // CV: keep minimum jet energy, in order not to loose direction information
  const double minJetEn = 1.e-2;
  if ( smearedJetEn < minJetEn ) smearedJetEn = minJetEn;

  // CV: skip smearing in case either "raw" or "corrected" jet energy is very low
  //     or jet passes selection configurable via python
  //    (allows for protection against "pathological cases",
  //     cf. PhysicsTools/PatUtils/python/tools/metUncertaintyTools.py)
  TLorentzVector smearedJetP4 = jet->p4();
  if ( !( abs(jet->energy() - rawJet->energy() ) > (5.*min(jet->energy(), rawJet->energy() )) 
	   //(skipJetSelection_ && (*skipJetSelection_)(jet)) || 
	 rawJet.pt()  < skipRawJetPtThreshold_          ||
	 jet.pt() < skipCorrJetPtThreshold_         ) ) {
    smearedJetP4 *= (smearedJetEn/jet.energy());
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
