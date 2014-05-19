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
  bool boolTightIso_l1 = false;
  int pdgCode_l1 = 0;
  float pt_l1, eta_l1, phi_l1, charge_l1;
  pt_l1 = eta_l1 = phi_l1 = charge_l1 = 0;
  bool boolTightIso_l2 = false;
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
	  boolTightIso_l1 = ElectronSel.accept(*l1);
	  boolTightIso_l2 = ElectronSel.accept(*l2);
	  if (boolTightIso_l1 && boolTightIso_l2) {
		  boolValidZ = true;
		  _zChan = 1;
	  }
  }
  if ( nMuons >= 2 && !boolValidZ ) {
	  l1 = Muons[0];
	  l2 = Muons[1];
	  // XXX fill("m_ll", "mumu", m_ll, "");
	  boolTightIso_l1 = MuonSel.accept(*l1);
	  boolTightIso_l2 = MuonSel.accept(*l2);
	  if (boolTightIso_l1 && boolTightIso_l2) {
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

  tm.add<bool>("boolTightIso_l1", &boolTightIso_l1);
  tm.add<int>("pdgCode_l1", &pdgCode_l1);
  tm.add<float>("pt_l1", &pt_l1);
  tm.add<float>("eta_l1", &eta_l1);
  tm.add<float>("phi_l1", &phi_l1);
  tm.add<float>("charge_l1", &charge_l1);
  tm.add<bool>("boolTightIso_l2", &boolTightIso_l2);
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
	  if( jet->pt() > 30 ) {
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

  // stat histograms
  //fillStatHistograms();

  // debug printouts
  //debugPrintouts( cout );


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
DMAnalysis::HLTFiredMultiChan() {

  if( (_zChan==0 && _e.isFired("HLT_Mu17_Mu8",1)) || 
      (_zChan==1 && _e.isFired("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL",1)) )
      // (_zChan==2 && _e.isFired("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL",1)) )
    return true;
  else  
    return false;
}
