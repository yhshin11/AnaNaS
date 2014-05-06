
#include "Analysis/src/TnP.hh"

#include <algorithm>
#include <cassert>
#include <sstream>

using namespace std;

int TnP::imode = EventServer::kWm;

ClassImp( TnP )

TnP::TnP( Sample& sample, const std::string & collectionFileName ) :
  SampleAnalysis( "TnP" , sample, collectionFileName ),
  outTree_("tnpVariables","Variables for Diphoton Analysis")
{
  debug       = false;
  _nDebug     = 0;
  SampleName  = sample.name();

  // set counters
  Nevt_  = 0;
  Nnosc_ = 0;

  // set constants
  _minPt = 20;
  _minM  = 50;
  _maxM  = 120;
}

TnP::~TnP()
{
}

void 
TnP::bookHistograms()
{
  initRootTree();
}

bool
TnP::analyzeEvent()
{
  if( debug )            { cout << "Event " << _ievt << endl; if(_ievt<1) return false; }
  else if( _ievt%500 == 0) cout << "Event " << _ievt << endl;
  Nevt_++;
  
  // isData
  if( SampleName.substr(0,6)=="Photon" || SampleName.substr(0,2)=="EG" ) isData = true;
  else                                                                   isData = false;

  // find tag_ and probe_
  if( !preselection() ) return false;

  // fill tree
  fillTree();
  return true;
}

bool
TnP::preselection()
{
  tags.clear();
  probes.clear();
  eprobes.clear();

  // looking for tags
  const CandList& electrons_ = _e.electrons();
  for( int e=0; e<(int)electrons_.size(); e++) {
    Candidate* ele = electrons_[e];
    if( ele->pt()<_minPt ) continue;
    CandInfo* elei = ele->info();
    float caloEta = elei->getFloat("caloEta");
    if( caloEta>2.5 ) continue;
    if( caloEta>1.4442 && caloEta<1.560 ) continue;
    // hardcore tag cuts
    if( !isTag( ele ) ) continue;
    if( isData && isSpike( ele ) ) continue;
    tags.push_back(ele);
  }
  if( tags.size()==0 ) return false;
  else {
    /*tag_ = tags[0];
      for( int t=0; t<(int)tags.size(); t++ ) if( tags[t]->pt()>tag_->pt() ) tag_=tags[t];*/
    int di = (int)round(rg.Rndm()*tags.size()-.5);
    tag_ = tags[di];
  }
  if( debug ) cout << tags.size() << " tags found." << endl;

  // looking for probes
  for( int e=0; e<(int)electrons_.size(); e++) {
    Candidate* ele = electrons_[e];
    if( ele->pt()<_minPt ) continue;
    CandInfo* elei = ele->info();
    float caloEta = elei->getFloat("caloEta");
    if( caloEta>2.5 ) continue;
    if( caloEta>1.4442 && caloEta<1.560 ) continue;
    // spike cleaning
    if( isData && isSpike( ele ) ) continue;
    // inv mass cut
    Candidate* zcand = Candidate::create( tag_,ele );
    float zcandm = zcand->mass();
    if( zcandm<_minM || zcandm>_maxM ) continue;
    // some loose cuts
    if     ( caloEta<1.45 && elei->getFloat("HOE")>0.15 ) continue;
    else if( caloEta>1.45 && elei->getFloat("HOE")>0.07 ) continue;
    if     ( caloEta<1.45 && elei->getFloat("sigiEtaiEta")>0.015 ) continue;
    else if( caloEta>1.45 && elei->getFloat("sigiEtaiEta")>0.035 ) continue;
    // match to a photon and keep it
    ConeSelector coneSel(*ele, 0.05);
    const CandList& photons_ = _e.photons();
    for( int p=0; p<(int)photons_.size(); p++) if( coneSel.accept(*photons_[p]) ) { 
      probes.push_back(photons_[p]); 
      eprobes.push_back(ele);
    }
  }
  if( debug ) cout << probes.size() << " probes found." << endl;
  
  // decide which probe to keep

  /*if( probes.size()>0 ) {
    if     ( probes.size()>2 ) return false;
    else if( probes.size()==2 ) {
      bool iT0 = isTag(eprobes[0]);
      bool iT1 = isTag(eprobes[1]);
      if     (  iT0 &&  iT1 ) { probe_ = probes[rg.Rndm()>.5];}
      else if( !iT0 &&  iT1 ) { probe_ = probes[1]; }
      else if(  iT0 && !iT1 ) { probe_ = probes[0]; }
      else if( !iT0 && !iT1 ) { probe_ = probes[rg.Rndm()>.5];}
    }
    else probe_ = probes[0];
    return true;
  }*/

  if( probes.size() ) {
    int divineintervention = (int)round(rg.Rndm()*probes.size()-.5);
    probe_ = probes[divineintervention];
    return true;
  }

  else return false;
}

bool
TnP::isTag( Candidate* cand )
{
  bool pass = true;
  CandInfo*  candi  = cand->info();
  if( candi->getFloat("caloEta")<1.4442 ) {
    //if( candi->getBool("ecalDriven")!=true ) return false;
    if( candi->getFloat("dr03TrkIso")/cand->pt()>0.09 ) return false;
    if( candi->getFloat("dr03EcalIso")/cand->pt()>0.07 ) return false;
    if( (candi->getFloat("dr04HcalD1Iso")+candi->getFloat("dr04HcalD2Iso"))/cand->pt()>.1 ) return false;
    if( candi->getFloat("sigiEtaiEta")>0.01 ) return false;
    if( fabs( candi->getFloat("dPhiIn") )>0.06 ) return false;
    if( fabs( candi->getFloat("dEtaIn") )>0.004 ) return false;
    if( candi->getFloat("HOE")>0.04 ) return false;
  }
  else if( candi->getFloat("caloEta")<1.560 ) return false;
  else if( candi->getFloat("caloEta")<2.5 ) {
    //if( candi->getBool("ecalDriven")!=true ) return false;
    if( candi->getFloat("dr03TrkIso")/cand->pt()>0.04 ) return false;
    if( candi->getFloat("dr03EcalIso")/cand->pt()>0.05 ) return false;
    if( (candi->getFloat("dr04HcalD1Iso")+candi->getFloat("dr04HcalD2Iso"))/cand->pt()>0.025 ) return false;
    if( candi->getFloat("sigiEtaiEta")>0.03 ) return false;
    if( fabs( candi->getFloat("dPhiIn") )>0.7 ) return false;
    if( fabs( candi->getFloat("dEtaIn") )>99. ) return false;
    if( candi->getFloat("HOE")>0.025 ) return false;
  }
  if( debug ) cout << "isTag done." << endl;
  return pass;
}

void
TnP::writeHistograms() // write all the output
{ 

  TDirectory* savedDir = gDirectory;
  _outputFile->cd();
  outTree_.Write(); // do not touch
  savedDir->cd();
}

void 
TnP::initRootTree() 
{
  // general stuff
  outTree_.Branch("run",               &run,               "run/I");
  outTree_.Branch("event",             &event,             "event/I");
  outTree_.Branch("nProbes",           &nProbes,           "nProbes/I");
  outTree_.Branch("nGoodVtx",          &nGoodVtx,          "nGoodVtx/I");
  // probe variables
  outTree_.Branch("p_eta",             &p_eta,             "p_eta/F");
  outTree_.Branch("p_phi",             &p_phi,             "p_phi/F");
  outTree_.Branch("p_pt",              &p_pt,              "p_pt/F");
  outTree_.Branch("p_sigIeta",         &p_sigIeta,         "p_sigIeta/F");
  outTree_.Branch("p_hcalIso",         &p_hcalIso,         "p_hcalIso/F");
  outTree_.Branch("p_ecalIso",         &p_ecalIso,         "p_ecalIso/F");
  outTree_.Branch("p_trackIso",        &p_trackIso,        "p_trackIso/F");
  outTree_.Branch("p_r9",              &p_r9,              "p_r9/F");
  outTree_.Branch("p_hoe",             &p_hoe,             "p_hoe/F");
  // common variables
  outTree_.Branch("mass",              &mass,              "mass/F");
}

void
TnP::fillTree() {

  // general stuff
  run     = _e.run();
  event   = _e.event();
  nProbes = probes.size();
  
  // vertex
  const VtxList& vertices = _e.vertices();
  nGoodVtx = 0;
  for( int v=0; v<(int)vertices.size(); v++) {
    Vertex* ver = vertices[v];
    CandInfo* infop = ver->info();
    bool vergo = infop->getBool("isGoodVtx");
    if( vergo ) nGoodVtx++;
  }
  Vertex* Primary = vertices[0];
  vX = Primary->pos().X();
  vY = Primary->pos().Y();
  vZ = Primary->pos().Z();
  
  // trigger bits
  HLT_Photon10_Cleaned_L1R             = _e.isFired("HLT_Photon10_Cleaned_L1R");
  HLT_Photon15_Cleaned_L1R             = _e.isFired("HLT_Photon15_Cleaned_L1R");
  HLT_Photon17_SC17HE_L1R_v1           = _e.isFired("HLT_Photon17_SC17HE_L1R_v1");
  HLT_Photon17_Isol_SC17HE_L1R_v1      = _e.isFired("HLT_Photon17_Isol_SC17HE_L1R_v1");
  HLT_Photon20_Cleaned_L1R             = _e.isFired("HLT_Photon20_Cleaned_L1R");
  HLT_Photon20_Isol_Cleaned_L1R_v1     = _e.isFired("HLT_Photon20_Isol_Cleaned_L1R_v1");
  HLT_Photon22_SC22HE_L1R              = _e.isFired("HLT_Photon22_SC22HE_L1R");
  HLT_Photon25_Cleaned_L1R             = _e.isFired("HLT_Photon25_Cleaned_L1R");
  HLT_Photon30_Cleaned_L1R             = _e.isFired("HLT_Photon30_Cleaned_L1R");
  HLT_Photon35_Isol_Cleaned_L1R_v1     = _e.isFired("HLT_Photon35_Isol_Cleaned_L1R_v1");
  HLT_Photon40_Isol_Cleaned_L1R_v1     = _e.isFired("HLT_Photon40_Isol_Cleaned_L1R_v1");
  HLT_Photon40_CaloId_Cleaned_L1R_v1   = _e.isFired("HLT_Photon40_CaloId_Cleaned_L1R_v1");
  HLT_Photon50_Cleaned_L1R             = _e.isFired("HLT_Photon50_Cleaned_L1R");
  HLT_Photon70_Cleaned_L1R             = _e.isFired("HLT_Photon70_Cleaned_L1R");
  HLT_DoublePhoton5_L1R                = _e.isFired("HLT_DoublePhoton5_L1R");
  HLT_DoublePhoton10_L1R               = _e.isFired("HLT_DoublePhoton10_L1R");
  HLT_DoublePhoton15_L1R               = _e.isFired("HLT_DoublePhoton15_L1R");
  HLT_DoublePhoton17_L1R               = _e.isFired("HLT_DoublePhoton17_L1R");
  HLT_DoublePhoton17_SingleIsol_L1R_v1 = _e.isFired("HLT_DoublePhoton17_SingleIsol_L1R_v1");
  HLT_DoublePhoton20_L1R               = _e.isFired("HLT_DoublePhoton20_L1R");
  HLT_DoublePhoton22_L1R_v1            = _e.isFired("HLT_DoublePhoton22_L1R_v1");
  
  // tag variables
  
  
  // probe variables
  CandInfo*  pinfo = probe_->info();
  Candidate* psc  = _e.getSecond("photon-SuperCluster",probe_);
  CandInfo*  psci = psc->info();
  p_eta      = fabs(probe_->eta());
  p_phi      = probe_->phi();
  p_pt       = probe_->pt();
  p_sigIeta  = pinfo->getFloat("sigmaIetaIeta");
  p_hcalIso  = pinfo->getFloat("hcalTowerSumEtConeDR04");
  p_trackIso = pinfo->getFloat("trkSumPtHollowConeDR04");
  p_ecalIso  = pinfo->getFloat("ecalRecHitSumEtConeDR04");
  p_r9       = psci->getFloat("R9");
  p_hoe      = pinfo->getFloat("hadronicOverEm");

  // common variables
  Candidate* zc = Candidate::create( tag_,probe_ );  
  mass = zc->mass();
  
  // finish it
  TDirectory* savedDir = gDirectory;
  _outputFile->cd();
  outTree_.Fill();
  savedDir->cd();

}

bool
TnP::isSpike( Candidate* cand ) {

  bool spike=false;
   
  //Load RecHits
  const CandMap&         Sc_Bc_ = _e.candMap("SuperCluster-BasicCluster");
  const CandIdRecHitMap& Bc_rh_ = _e.ecalRecHitMap("BasicCluster-EcalRecHit");

  Candidate* superCluster = _e.getSecond("electron-SuperCluster", cand );
  CandMapIterator itbc_lower_ = Sc_Bc_.lower_bound( superCluster );
  CandMapIterator itbc_upper_ = Sc_Bc_.upper_bound( superCluster );
  CandMapIterator itbc;
  for (itbc=itbc_lower_; itbc!=itbc_upper_; itbc++) { 
    CandIdRecHitMapIterator itrh_lower_ = Bc_rh_.lower_bound( (*itbc).second->uid() );
    CandIdRecHitMapIterator itrh_upper_ = Bc_rh_.upper_bound( (*itbc).second->uid() );
    CandIdRecHitMapIterator itrh;
    for (itrh=itrh_lower_; itrh!=itrh_upper_; itrh++) { 
      ecalRecHit rh = (*itrh).second.first;
      if ( rh.e > 10 && rh.iz==0 ) {	    
        float swisscross = rh.SpikeVar.swCross;
        if( swisscross > 0.95 ) { spike=true; }
	int severitylevel = rh.SpikeVar.sevLvl;
	if( severitylevel > 0 ) { spike=true; }
      }
    }
  }
  if( debug ) cout << "Spikes cleaned." << endl;
  return spike;
}
