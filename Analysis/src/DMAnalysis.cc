#include "Analysis/src/DMAnalysis.hh"
 
#include <cassert>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <math.h>
using namespace std;

#include <TH1I.h>
#include <TH1F.h>
#include <TH2I.h>
#include <TH2F.h>

#include "Analysis/utils/Config.hh"
#include "Analysis/utils/Constants.hh"
#include "Analysis/utils/KineUtils.hh"
#include "Analysis/utils/CutUtils.hh"
#include "Analysis/core/Candidate.hh"
#include "Analysis/core/CandInfo.hh"
#include "Analysis/core/EventManager.hh"
#include "Analysis/tools/CandUtil.hh"
#include "Analysis/selectors/ElectronSelector.hh"
#include "Analysis/selectors/MuonSelector.hh"
#include "Analysis/selectors/LeptonSelector.hh"
#include "Analysis/selectors/ConeSelector.hh"

ClassImp( DMAnalysis )

typedef vector<float> vectorFloat; 


DMAnalysis::DMAnalysis( Sample& sample, EventManager& manager ) : 
SampleAnalysis( "DM", sample, manager ),
MuonSel(MuonSelector::kPt20,MuonSelector::kTight),
ElectronSel(ElectronSelector::kPt20,ElectronSelector::kMedium)
{
  cout << "\t-------------------------------------------" << endl; 
  cout << "\t---  DM preselection analysis     ----" << endl; 
  cout << "\t-------------------------------------------" << endl; 

  _hlt = false; // apply HLT selection
  _nDebug = 0;
  //XXX hist_mll = new TH1F("m_ll", "m_ll", 100, 0, 200);

}

DMAnalysis::~DMAnalysis()
{
}

void
DMAnalysis::bookHistograms()
{
//  defineTemplate( "Eproj", 400, -200, 200 );
//  defineTemplate( "Eproj2D", 200, -200, 200, 200, -200, 200 );
//  defineTemplate( "METVsMT", 80, 0, 80, 80, 0, 160 );
//  defineTemplate( "sigMET",  200, 0, 20    );
//  defineTemplate( "balance", 200, 0, 20    );
//  defineTemplate( "jmult", 20, -0.5, 19.5 );
//  defineTemplate( "LP", 100, -2., 2. );
  defineTemplate("m_ll", 100, 0, 200);
  tm.setTree("NTuple", "");
	cout << "DMAnalysis::bookHistograms()" << endl;
  
}

void
DMAnalysis::writeHistograms()
{
//  for( map< string, TH1* >::const_iterator it=h_.begin(); 
//       it!=h_.end(); it++ )
//    {
//      it->second->Write();
//    }
//  for( map< string, TH2* >::const_iterator it=h2_.begin(); 
//       it!=h2_.end(); it++ )
//    {
//      it->second->Write();
//    }  
//
  //XXX hist_mll->Write();

	// Write out histograms.
  tm.flush();
}

bool
DMAnalysis::analyzeEvent()
{
  int eventIndex = _ievt;
  if(eventIndex % 1000 == 0) {
	  cout << "Analyzing event " << eventIndex << endl;
  }

  // build lists of leptons and prepare the MC matching map
  buildLeptonLists();


  // Retrieve leptons.
  const CandList& Electrons = _e.electrons();
  const CandList& Muons = _e.muons();
  int nElectrons, nMuons;
  nElectrons = Electrons.size();
  nMuons = Muons.size();
  Candidate *l1, *l2;
  // Lepton variables to store.
  bool boolSelector_l1 = false;
  int pdgCode_l1 = 0;
  float pt_l1, eta_l1, phi_l1, charge_l1;
  pt_l1 = eta_l1 = phi_l1 = charge_l1 = 0;
  bool boolSelector_l2 = false;
  int pdgCode_l2 = 0;
  float pt_l2, eta_l2, phi_l2, charge_l2;
  pt_l2 = eta_l2 = phi_l2 = charge_l2 = 0;

  Candidate* ZCand;
  // Z candidate variables to store.
  bool boolValidZ = false;
  float m_ll, pt_ll, phi_ll, y_ll;

  if ( nElectrons >= 2 && !boolValidZ ) {
	  l1 = Electrons[0];
	  l2 = Electrons[1];
	  // XXX fill("m_ll", "ee", m_ll, "");
	  boolSelector_l1 = ElectronSel.accept(*l1);
	  boolSelector_l2 = ElectronSel.accept(*l2);
	  if (boolSelector_l1 && boolSelector_l2) {
		  boolValidZ = true;
		  _zChan = 1;
	  }
  }
  if ( nMuons >= 2 && !boolValidZ ) {
	  l1 = Muons[0];
	  l2 = Muons[1];
	  // XXX fill("m_ll", "mumu", m_ll, "");
	  boolSelector_l1 = MuonSel.accept(*l1);
	  boolSelector_l2 = MuonSel.accept(*l2);
	  if (boolSelector_l1 && boolSelector_l2) {
		  boolValidZ = true;
		  _zChan = 0;
	  }
  }
  if (!boolValidZ) {
	  //cout << "no Z candidate in this event" << endl;
  }
  if (boolValidZ) {
		pdgCode_l1 = l1->pdgCode();
		pt_l1 = l1->pt();
		eta_l1 = l1->eta();
		phi_l1 = l1->phi();
		charge_l1 = l1->charge();
		pdgCode_l2 = l2->pdgCode();
		pt_l2 = l2->pt();
		eta_l2 = l2->eta();
		phi_l2 = l2->phi();
		charge_l2 = l2->charge();
		// Create Z candidate with 2 leptons.
		ZCand = Candidate::create(l1, l2);
		m_ll = ZCand->mass();
		pt_ll = ZCand->pt();
		phi_ll = ZCand->phi();
		y_ll = getRapidity(ZCand);

  }

  if(!HLTFiredMultiChan() ) return false;

  // Check Z-candidate validity
  if (!boolValidZ) {return false;}

  // Check Z-candidate pt (threshold at 30 GeV)
  if (!pt_ll > 30) {return false;}

  // Count the number of "good" vertices. Check that there are at least one.
  if(!checkVertices()) return false;

	int run = _e.run();
  int lumisec = _e.lumisec();
  int event = _e.event();
  int eventInFile = _e.eventInFile();
  tm.add<int>("run", &run);
  tm.add<int>("lumisec", &lumisec);
  tm.add<int>("event", &event);
  tm.add<int>("eventInFile", &eventInFile);

	// cout << "Filling lepton info." << endl;
  tm.add<bool>("boolSelector_l1", &boolSelector_l1);
  tm.add<int>("pdgCode_l1", &pdgCode_l1);
  tm.add<float>("pt_l1", &pt_l1);
  tm.add<float>("eta_l1", &eta_l1);
  tm.add<float>("phi_l1", &phi_l1);
  tm.add<float>("charge_l1", &charge_l1);
  tm.add<bool>("boolSelector_l2", &boolSelector_l2);
  tm.add<int>("pdgCode_l2", &pdgCode_l2);
  tm.add<float>("pt_l2", &pt_l2);
  tm.add<float>("eta_l2", &eta_l2);
  tm.add<float>("phi_l2", &phi_l2);
  tm.add<float>("charge_l2", &charge_l2);

  tm.add<bool>("boolValidZ", &boolValidZ);
  tm.add<float>("m_ll", &m_ll);
  tm.add<float>("pt_ll", &pt_ll);
  tm.add<float>("phi_ll", &phi_ll);
  tm.add<float>("y_ll", &y_ll);
  //XXX hist_mll->Fill(m_ll);

  // Save additional leptons information. (PDG ID and Loose boolean)
  vector<int> *pdgCode_addLeptons = new vector<int>;
  vector<float> *pt_addLeptons = new vector<float>;
  vector<float> *eta_addLeptons = new vector<float>;
  vector<float> *phi_addLeptons = new vector<float>;
  vector<bool> *Loose_addLeptons = new vector<bool>;
  for(int unsigned il=0;il<Electrons.size();il++) {
	  Candidate* lepton = Electrons[il];
	  if (lepton != l1 && lepton != l2 && lepton->pt() > 10) {
			ElectronSel.accept(*lepton);
			pdgCode_addLeptons->push_back(lepton->pdgCode());
			pt_addLeptons->push_back(lepton->pt());
			eta_addLeptons->push_back(lepton->eta());
			phi_addLeptons->push_back(lepton->phi());
			bool boolLoose = lepton->info()->getBool("Loose");
			Loose_addLeptons->push_back(boolLoose);
	  }
  }
  for(int unsigned il=0;il<Muons.size();il++) {
	  Candidate* lepton = Muons[il];
	  if (lepton != l1 && lepton != l2 && lepton->pt() > 10) {
			MuonSel.accept(*lepton);
			pdgCode_addLeptons->push_back(lepton->pdgCode());
			pt_addLeptons->push_back(lepton->pt());
			eta_addLeptons->push_back(lepton->eta());
			phi_addLeptons->push_back(lepton->phi());
			bool boolLoose = lepton->info()->getBool("Loose");
			Loose_addLeptons->push_back(boolLoose);
	  }
  }
	// // Insert dummy values if no additional leptons are present.
	// if (pt_addLeptons->size()==0)
	// {
	// 		pdgCode_addLeptons->push_back(0);
	// 		pt_addLeptons->push_back(-100.0);
	// 		eta_addLeptons->push_back(-100.0);
	// 		phi_addLeptons->push_back(-100.0);
	// 		Loose_addLeptons->push_back(false);
	// }
	// cout << "pt_addLeptons.size(): " << pt_addLeptons->size() << endl;
	// cout << "Filling additional lepton info." << endl;
  tm.add<vector<int>*>("pdgCode_addLeptons", &pdgCode_addLeptons);
  tm.add<vector<float>*>("pt_addLeptons", &pt_addLeptons);
  tm.add<vector<float>*>("eta_addLeptons", &eta_addLeptons);
  tm.add<vector<float>*>("phi_addLeptons", &phi_addLeptons);
  tm.add<vector<bool>*>("Loose_addLeptons", &Loose_addLeptons);

  // Retrieve jets.
  const CandList& Jets = _e.jetList( EventManager::kPatJet);
  // Jet variables to store.
  vector<float> *pt_jets = new vector<float>;
  vector<float> *eta_jets = new vector<float>;
  vector<float> *phi_jets = new vector<float>;
  vector<bool> *mvaIdMedium_jets = new vector<bool>;
  vector<float> *CSVBT_jets	= new vector<float>;
  vector<float> *CSVMVA_jets = new vector<float>;
  for(int unsigned ij=0;ij<Jets.size();ij++) { 
	  Candidate* jet = Jets[ij];
	  if( jet->pt() > 30 && !isLepton(jet)) {
		  pt_jets->push_back(jet->pt());
		  eta_jets->push_back(jet->eta());
		  phi_jets->push_back(jet->phi());
		  CandInfo* info = jet->info();
		  bool mvaIdMedium = info->getBool("mvaIdMedium");
		  mvaIdMedium_jets->push_back(mvaIdMedium);
		  float CSVBT = info->getFloat("btagCSVBT");
		  float CSVMVA = info->getFloat("btagCSVMVA");
		  CSVBT_jets->push_back(CSVBT);
		  CSVMVA_jets->push_back(CSVMVA);
	  }
  }
	// cout << "Filling jet info." << endl;
  tm.add<vector<float>*>("pt_jets", &pt_jets);
  tm.add<vector<float>*>("eta_jets", &eta_jets);
  tm.add<vector<float>*>("phi_jets", &phi_jets);
  tm.add<vector<bool>*>("mvaIdMedium_jets", &mvaIdMedium_jets);
  tm.add<vector<float>*>("CSVBT_jets", &CSVBT_jets);
  tm.add<vector<float>*>("CSVMVA_jets", &CSVMVA_jets);

  //for(int unsigned ij=0;ij<Jets.size();ij++) { 
  //  Candidate* jet = Jets[ij];
  //  if((fabs(jet->eta())<1.4442 ||
  //  (fabs(jet->eta())<2.5 && fabs(jet->eta())>1.56 ) ) )
  //    {
  //  if(pttmp<jet->pt())
  //    { theJet_=jet; pttmp = jet->pt(); }
  //    }
  //}

  // Retrieve met candidates.
  //PATType1CorrectedPFMet
  //MVAPFMet
  //MVANoPUMET
  Candidate* pfmet;
  Candidate* pfmetMva;
  if( !(_e.a().isMC) )  { //data 
	  pfmet = const_cast<Candidate*>(_e.met( "pat", "patType1CorrectedPFMetNoSmear" ));
	  pfmetMva = const_cast<Candidate*>(_e.met( "pat", "patPFMetMVANoSmear" ));
  }
  else {
	  pfmet = const_cast<Candidate*>(_e.met( "pat", "patType1CorrectedPFMetSmeared" ));
	  pfmetMva = const_cast<Candidate*>(_e.met( "pat", "patPFMetMVASmeared" ));
  }
  // MET variables to store.
  float pt_met, phi_met, sig_met;
  pt_met = phi_met = sig_met = 0;
  {
	  CandInfo* info = pfmet->info();
	  pt_met = pfmet->pt();
	  phi_met = pfmet->phi();
	  sig_met = info->getFloat("mEtSig");
  }
  tm.add<float>("pt_met", &pt_met);
  tm.add<float>("phi_met", &phi_met);
  tm.add<float>("sig_met", &sig_met);
  // MVA MET variables to store.
  float pt_metMva, phi_metMva, sig_metMva;
  pt_metMva = phi_metMva = sig_metMva = 0;
  {
	  CandInfo* info = pfmetMva->info();
	  pt_metMva = pfmetMva->pt();
	  phi_metMva = pfmetMva->phi();
	  sig_metMva = info->getFloat("mEtSig");
  }
  tm.add<float>("pt_metMva", &pt_metMva);
  tm.add<float>("phi_metMva", &phi_metMva);
  tm.add<float>("sig_metMva", &sig_metMva);

  // Get true number of interactions.
  float trueI = _e.getTrueNInt();
  tm.add<float>( "trueNint", &trueI ); 

	fillGenEvent();
	fillWChan();

  // stat histograms
  //fillStatHistograms();

  // debug printouts
  //debugPrintouts( cout );

  delete pdgCode_addLeptons;
  delete pt_addLeptons;
  delete eta_addLeptons;
  delete phi_addLeptons;
  delete Loose_addLeptons;

  delete pt_jets;
  delete eta_jets;
  delete phi_jets;
  delete mvaIdMedium_jets;
  delete CSVBT_jets;
  delete CSVMVA_jets;

  //return nSelected>0;
  return true;
}

float 
DMAnalysis::getRapidity(const Candidate* zCandidate)
{
	if (zCandidate) {
		float pz = zCandidate->pz();
		float E = zCandidate->E();
		float rapidity = 0.5*log( (E + pz)/(E - pz) );
		return rapidity;
	}
	else return 0;
}

bool
DMAnalysis::HLTFiredMultiChan()
{

  if( (_zChan==0 && _e.isFired("HLT_Mu17_Mu8",1)) || 
      (_zChan==1 && _e.isFired("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL",1)) )
      // (_zChan==2 && _e.isFired("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL",1)) )
    return true;
  else  
    return false;
}

bool
DMAnalysis::checkVertices()
{

  const VtxList& vertices = _e.vertices();
 
  Nvertex = 0;

	Vertex* Primary;
   
  bool findPV=false;
  //DZ vertices
  
  for(int unsigned iv=0;iv<vertices.size();iv++) {
    
    Vertex* Vert = vertices[iv];
    CandInfo* infoV = Vert->info();
   
    bool goodV=false;
    if(fabs((Vert->pos()).Z())<24 && (Vert->pos()).Perp()<2 && infoV->getFloat("ndof")> 4 && !(infoV->getBool("isFake")) )
      goodV=true;
  
    if(!findPV && goodV )
      {Primary = Vert; Nvertex++; findPV=true; continue;}
      
    if(goodV && findPV && vertices.size()>1) {
      Nvertex++;
    }
  }
 
  if(findPV)
    return true;
  else
    return false;

}

bool
DMAnalysis::isLepton(const Candidate* jet)
{
  const CandList& Electrons = _e.electrons();
  const CandList& Muons = _e.muons();

  // if( jet->dR( _ZCand->daughter(0) )<0.2 && jet->pt()/_ZCand->daughter(0)->pt() >0.8 ) return true;
  // if( jet->dR( _ZCand->daughter(1) )<0.2 && jet->pt()/_ZCand->daughter(1)->pt() >0.8 ) return true;

  for(int unsigned il=0;il<Electrons.size();il++) {
	  Candidate* lepton = Electrons[il];
	  if (lepton->pt() > 10) {
			ElectronSel.accept(*lepton);
			bool boolLoose = lepton->info()->getBool("Loose");
			if (boolLoose) { 
				if (jet->dR(lepton) < 0.2) return true;
			}
	  }
  }

  return false;
}

void 
DMAnalysis::fillGenEvent()
{

  const CandList& mcCands = _e.mcTruthCandidates();
	CandList motherCands;

	// Find valid mother particles (Z or W).
  for(int unsigned imc=0;imc<mcCands.size();imc++) {
    Candidate* mcCand = mcCands[imc];
		if (mcCand->pdgCode() == 23 || abs(mcCand->pdgCode()) == 24) {
			motherCands.push_back(mcCand);
		}
	}

	// Variables to store
  vector<int> *pdgCode_genLeptons = new vector<int>;
  vector<float> *pt_genLeptons = new vector<float>;
  vector<float> *eta_genLeptons = new vector<float>;
  vector<float> *phi_genLeptons = new vector<float>;
  vector<int> *motherPdgCode_genLeptons = new vector<int>;

	// Loop through list of valid mother particles and record the index and pdg code.
	// imc is equal to "index" in mcTruth tree from NTupleMaker.
  for(int unsigned imc=0; imc<mcCands.size(); imc++) {
    Candidate* mcCand = mcCands[imc];// cout<<" debug "<<imc<<"  "<<mcCands.size()<<endl;

		bool isLeptonAndStable = mcCand->info()->getInt("status") == 1
													&& abs(mcCand->pdgCode()) >= 11 && abs(mcCand->pdgCode()) <= 16 ;
		if (isLeptonAndStable && mcCand->pt()>10) {
			int pdgCode = mcCand->pdgCode();
			float pt = mcCand->pt();
			float eta = mcCand->eta();
			float phi = mcCand->phi();
			int motherPdgCode = 0; 

			Candidate* ancestor = findAncestorOfMcCandidate(mcCand, mcCands, motherCands);
			if(ancestor!=NULL) motherPdgCode = ancestor->pdgCode();

			pdgCode_genLeptons->push_back(pdgCode);
			pt_genLeptons->push_back(pt);
			eta_genLeptons->push_back(eta);
			phi_genLeptons->push_back(phi);
			motherPdgCode_genLeptons->push_back(motherPdgCode);

		}
	}

  tm.add<vector<int>*>("pdgCode_genLeptons", &pdgCode_genLeptons);
  tm.add<vector<float>*>("pt_genLeptons", &pt_genLeptons);
  tm.add<vector<float>*>("eta_genLeptons", &eta_genLeptons);
  tm.add<vector<float>*>("phi_genLeptons", &phi_genLeptons);
  tm.add<vector<int>*>("motherPdgCode_genLeptons", &motherPdgCode_genLeptons);

  delete pdgCode_genLeptons;
  delete pt_genLeptons;
  delete eta_genLeptons;
  delete phi_genLeptons;
  delete motherPdgCode_genLeptons;

  const CandList& mcJets = _e.jetList( EventManager::kGenJet);
	// Variables to store
  vector<float> *pt_genJets = new vector<float>;
  vector<float> *eta_genJets = new vector<float>;
  vector<float> *phi_genJets = new vector<float>;
	for (int i=0; i<mcJets.size(); i++) {
		Candidate* genJet = mcJets[i]; 
		float pt = genJet->pt();
		float eta = genJet->eta();
		float phi = genJet->phi();
		pt_genJets->push_back(pt);
		eta_genJets->push_back(eta);
		phi_genJets->push_back(phi);
	}
  tm.add<vector<float>*>("pt_genJets", &pt_genJets);
  tm.add<vector<float>*>("eta_genJets", &eta_genJets);
  tm.add<vector<float>*>("phi_genJets", &phi_genJets);

	delete pt_genJets;
	delete eta_genJets;
	delete phi_genJets;
}

Candidate*
DMAnalysis::findAncestorOfMcCandidate(Candidate* mcCand, const CandList& mcCands, CandList& motherCands)
{
	int pdgCode = mcCand->pdgCode();
	Candidate* previousAncestor = findMotherOfMcCandidate(mcCand, mcCands);
	// while( previousAncestor->pdgCode() == pdgCode ) {
	while( !isInList(previousAncestor, motherCands) && previousAncestor!=NULL ) {
		previousAncestor = findMotherOfMcCandidate(previousAncestor, mcCands);
	}
	return previousAncestor;
}

Candidate*
DMAnalysis::findMotherOfMcCandidate(Candidate* mcCand, const CandList& mcCands)
{
	assert(mcCand!=NULL);
	int motherId = mcCand->info()->getInt("motherIndex");
	Candidate* mother;
	if (motherId != -1) {
		mother = mcCands[motherId];
	}
	else {
		mother = NULL;
	}
	return mother;
}

bool
DMAnalysis::isInList(const Candidate* cand, const CandList& candlist)
{
	bool isInList_ = false;
	for (int i=0; i<candlist.size(); i++) {
		if (cand == candlist[i]) {
			isInList_ = true;
			break;
		}
	}
	return isInList_;
}

void 
DMAnalysis::fillWChan() 
{
  const CandList& mcCands = _e.mcTruthCandidates();
	CandList WCands;
	vector<int> *decayChannels = new vector<int>;

	// Find valid mother particles (W).
  for(int unsigned imc=0; imc<mcCands.size(); imc++) { 
    Candidate* mcCand = mcCands[imc];
		if (abs(mcCand->pdgCode()) == 24) {
			WCands.push_back(mcCand);
			WDecayChannel::decayChannel decayChan = classifyWDaughters(mcCand);
			decayChannels->push_back(decayChan);
		}
	}

  tm.add<vector<int>*>("WDecayChannels", &decayChannels);
	delete decayChannels;
}

WDecayChannel::decayChannel
DMAnalysis::classifyWDaughters(Candidate* mcCand) 
{
	assert(mcCand!=NULL);
	int nLep = 0, nLeptonicTau = 0, nHadronicTau = 0, nUnclassifiedTau = 0, nHad = 0;
	CandList daughters;
	findDaughters(mcCand, daughters);
	int nDau = daughters.size();
	for (int i=0; i<nDau; i++) {
		Candidate* daughter = daughters[i];
		assert(daughter!=NULL);
		int pdgCode = daughter->pdgCode();
		if (abs(pdgCode)==24) {
			// If daughter is another W, this means the W is intermediate.
			// Check its daughter instead.
			return classifyWDaughters(daughter);
		}
		else if (abs(pdgCode)>=11 && abs(pdgCode)<=14) {
			nLep++;
		}
		else if (abs(pdgCode)==15) {
			if (classifyTauDaughters(daughter)==TauDecayChannel::kHadronic) {
				nHadronicTau++;
			}
			else if (classifyTauDaughters(daughter)==TauDecayChannel::kLeptonic) {
				nLeptonicTau++;
			}
			else if (classifyTauDaughters(daughter)==TauDecayChannel::kUnclassified) {
				nUnclassifiedTau++;
			}
		}
		else if (abs(pdgCode)>=1 && abs(pdgCode)<=6) {
			nHad++;
		}
	}

	WDecayChannel::decayChannel decayChannel;
	// Hadronic decay channel.
	if (nHad >= 2) {
		decayChannel = WDecayChannel::kHadronic;
	}
	// Leptonic decay channel.
	else if (nLep >=2) {
		decayChannel = WDecayChannel::kLeptonic;
	}
	// Leptonic tau decay channel.
	else if (nLeptonicTau >= 1) {
		decayChannel = WDecayChannel::kTauLeptonic;
	}
	// Hadronic tau decay channel.
	else if (nHadronicTau >= 1) {
		decayChannel = WDecayChannel::kTauHadronic;
	}
	// Unclassified tau decay channel.
	else if (nUnclassifiedTau >= 1) {
		decayChannel = WDecayChannel::kTauUnclassified;
	}
	// Unclassified decay channel.
	else {
		int run = _e.run();
		int lumisec = _e.lumisec();
		int event = _e.event();
		cout << "Unclassified W decay at run, lumi, eventnumber: " << run << ", " << lumisec << ", " << event << endl;
	}
	decayChannel = WDecayChannel::kUnclassified;
	return decayChannel;
}

TauDecayChannel::decayChannel 
DMAnalysis::classifyTauDaughters(Candidate* mcCand) 
{
	assert(mcCand!=NULL);
	int nLep = 0, nHad = 0;
	CandList daughters;
	findDaughters(mcCand, daughters);
	int nDau = daughters.size();
	for (int i=0; i<nDau; i++) {
		Candidate* daughter = daughters[i];
		assert(daughter!=NULL);
		int pdgCode = daughter->pdgCode();
		if (abs(pdgCode)==15 | abs(pdgCode)==24) {
			// If daughter is an intermediate tau or W, check its daughter instead.
			return classifyTauDaughters(daughter);
		}
		else if (abs(pdgCode)==11 | abs(pdgCode)==13) {
			nLep++;
		}
		else if (abs(pdgCode)>=100) {
			nHad++;
		}
	}

	TauDecayChannel::decayChannel decayChannel;
	// Hadronic decay channel.
	if (nHad >= 1) {
		decayChannel = TauDecayChannel::kHadronic;
	}
	// Leptonic decay channel.
	else if (nLep >=1) {
		decayChannel = TauDecayChannel::kLeptonic;
	}
	// Unclassified decay channel.
	else {
		decayChannel = TauDecayChannel::kUnclassified;
		int run = _e.run();
		int lumisec = _e.lumisec();
		int event = _e.event();
		cout << "Unclassified tau decay at run, lumi, eventnumber: " << run << ", " << lumisec << ", " << event << endl;
	}
	return decayChannel;
}

Candidate*
DMAnalysis::createGenCand(Candidate* cand)
{
	assert(cand!=NULL);
	CandList daughters;
  const CandList& mcCands = _e.mcTruthCandidates();
  for(int unsigned imc=0;imc<mcCands.size();imc++) {
    Candidate* mcDaughter = mcCands[imc];
		assert(mcDaughter!=NULL);
		int motherId = mcDaughter->info()->getInt("motherIndex");
		if (cand==mcCands[motherId]) {
			Candidate* genDaughter = createGenCand(mcDaughter);
			assert(genDaughter!=NULL);
			daughters.push_back(genDaughter);
			MCToGen[mcDaughter] = genDaughter;
			GenToMC[genDaughter] = mcDaughter;
		}
	}
	Candidate* genCand = Candidate::create(daughters);
	// MCToGen[cand] = genCand;
	// GenToMC[genCand] = cand;
	return genCand;
}

// Returns dummy Candidate with daughters.
void
DMAnalysis::findDaughters(Candidate* cand, CandList& daughters)
{
  const CandList& mcCands = _e.mcTruthCandidates();
  for(int unsigned imc=0;imc<mcCands.size();imc++) {
    Candidate* mcDaughter = mcCands[imc];
		assert(mcDaughter!=NULL);
		int motherId = mcDaughter->info()->getInt("motherIndex");
		if (cand==mcCands[motherId]) {
			assert(mcDaughter!=NULL);
			daughters.push_back(mcDaughter);
		}
	}
}
