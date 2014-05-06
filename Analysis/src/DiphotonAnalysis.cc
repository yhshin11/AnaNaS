
#include "Analysis/src/DiphotonAnalysis.hh"

#include <algorithm>
#include <cassert>
#include <sstream>

using namespace std;

int DiphotonAnalysis::imode = EventServer::kWm;

ClassImp( DiphotonAnalysis )

DiphotonAnalysis::DiphotonAnalysis( Sample& sample, const std::string & collectionFileName ) :
  SampleAnalysis( "Diphotons" , sample, collectionFileName ),
  outTree_("PhotonVariables","Variables for Diphoton Analysis")
{
  // switches
  debug       = false;
  doDiphotons = false;
  addSigCut   = 0;
  noPixSeed   = 0;
  doShuffle   = true  && doDiphotons;
  doFSel      = false && doDiphotons;
  _nDebug     = 0;
  SampleName  = sample.name();

  // init counters
  Nevt_        = 0;
  Nacceptance_ = 0;
  NHLT_        = 0;
  Npresel_     = 0;

  // histograms
  SampleAnalysis::defineTemplate( "All",    500,0,2);
  SampleAnalysis::defineTemplate( "Pt",     500,0,2000);
  SampleAnalysis::defineTemplate( "dR",     1000,0,1);

  // basic selection
  pho1_minEt = 21.;
  pho2_minEt = 20.;
  pho_maxHoe = 0.05;

  // eta limits
  barrel_maxEta = 1.4442;
  endcap_minEta = 1.5600;
  endcap_maxEta = 2.5000;

  // further selections
  pho_maxHcal    = doFSel ? 3.     : 999. ;
  pho_maxTrack   = doFSel ? 2.     : 999. ;
  pho_maxSietaEB = doFSel ? 0.011  : 999. ;
  pho_maxSietaEE = doFSel ? 0.03   : 999. ; 
  pho_maxSieta=0.;

}

DiphotonAnalysis::~DiphotonAnalysis()
{
}

void 
DiphotonAnalysis::bookHistograms()
{
  //outTree_.SetDirectory(0);
  initRootTree();
}

bool
DiphotonAnalysis::analyzeEvent()
{
  if( (_ievt)%1000 == 0) cout<<" Event "<<_ievt<<endl;
  if( debug==true ) {
    cout<<" Event "<<_ievt<<endl;
    //if( _ievt < 4740 ) return false;
  }
    
  Nevt_++;
  fill( "All", "Nevt", 1.  );

  if(_e.inAcceptance()) Nacceptance_++;

  float pttmp=0;
  const CandList& Jets = _e.jetList( EventManager::kPfJet);
  for(int unsigned ij=0;ij<Jets.size();ij++) {
    Candidate* jet = Jets[ij];
    if((fabs(jet->eta())<1.4442 ||
	(fabs(jet->eta())<20.5 && fabs(jet->eta())>1.56 ) ) )
      {
	if(pttmp<jet->pt()) { hardestJet_=jet; pttmp = jet->pt(); }
      }
  }

  // HLT selection
  //if( SampleName.substr(0,1)=="W" || SampleName.substr(0,1)=="Z" || SampleName.substr(0,6)=="GamGam" || SampleName.substr(0,2)=="BO" || SampleName.substr(0,7)=="QCDPT20" ) NHLT_++;
  //else { if( !HLTFiring() ) return false; else NHLT_++; }
  fill( "All", "NHLT", 1.  );

  if( !EventPreselection()) return false; else Npresel_++;
  fill( "Pt", "Npresel", thePhoton1_->pt()  );
  if( debug ) cout << "Here." << endl;

  // if Data:
  if( SampleName.substr(0,2)=="EG" || SampleName.substr(0,6)=="Photon" ) {
    int run,event;
    for(int i=0;i<(int) Events.size();i++) {
      run = Events[i].first;
      event = Events[i].second;
      if(run == _e.run() && event==_e.event()) return false;
    }
    std::pair<int,int> tmp(_e.run(),_e.event());
    Events.push_back(tmp);
    hasGenMatch1 = hasGenMatch2 = false;
  }
  // if Monte Carlo
  else { 
    MCtruthMatching();
  }
  JetMatching();
  
  // fill the actual tree
  if( debug ) cout << "Start filling the tree." << endl;
  fillTree();
  return true;
}

void
DiphotonAnalysis::writeHistograms() // write all the output
{ 

  TDirectory* savedDir = gDirectory;
  _outputFile->cd();
  outTree_.Write(); // do not touch
  savedDir->cd();
}

bool
DiphotonAnalysis::HLTFiring() 
{
  bool HLTfiring=false;
  if( _e.isFired("HLT_Photon10_Cleaned_L1R") ||
      _e.isFired("HLT_Photon15_Cleaned_L1R") || 
      _e.isFired("HLT_Photon17_SC17HE_L1R_v1") ||
      _e.isFired("HLT_Photon17_Isol_SC17HE_L1R_v1") ||
      _e.isFired("HLT_Photon20_Cleaned_L1R") ||
      _e.isFired("HLT_Photon20_Isol_Cleaned_L1R_v1") ||
      _e.isFired("HLT_Photon22_SC22HE_L1R") ||
      _e.isFired("HLT_Photon25_Cleaned_L1R") ||
      _e.isFired("HLT_Photon30_Cleaned_L1R") ||
      _e.isFired("HLT_Photon35_Isol_Cleaned_L1R_v1") ||
      _e.isFired("HLT_Photon40_Isol_Cleaned_L1R_v1") ||
      _e.isFired("HLT_Photon40_CaloId_Cleaned_L1R_v1") ||
      _e.isFired("HLT_Photon50_Cleaned_L1R") ||
      _e.isFired("HLT_Photon70_Cleaned_L1R") ||
      _e.isFired("HLT_DoublePhoton5_L1R") ||
      _e.isFired("HLT_DoublePhoton10_L1R") ||
      _e.isFired("HLT_DoublePhoton15_L1R") ||
      _e.isFired("HLT_DoublePhoton17_L1R") ||
      _e.isFired("HLT_DoublePhoton17_SingleIsol_L1R_v1") ||
      _e.isFired("HLT_DoublePhoton20_L1R") ||
      _e.isFired("HLT_DoublePhoton22_L1R_v1")
    ) HLTfiring = true;
  else if( debug ) cout << "HLT fail." << endl;
  return HLTfiring;

}

bool
DiphotonAnalysis::EventPreselection() 
{
  // return and store
  bool presel=false;
  hasPhoton1=false;
  hasPhoton2=false;

  // changing variables
  float et=0.;
  float hovere=0.;
  float eta=0.;
  float sigeta=0.;
  float hiso=0.;
  float tiso=0.;
  bool  haspix=false;
  bool  isEB=false;

  // photon selection
  vector<int> firstPhos, secndPhos; // indices
  vector<int> firstPacc, secndPacc; // whether they pass further cuts
  firstPhos.clear();
  firstPacc.clear();
  secndPhos.clear();
  secndPacc.clear();

  // run on photon collection
  const CandList& photons_ = _e.photons();

  for (int unsigned ig=0; ig<photons_.size(); ig++) {
    Candidate* photon_ = photons_[ig];
    CandInfo* info = photon_->info();
    
    et       = photon_->pt();
    eta      = info->getFloat("caloEta");
    hovere   = info->getFloat("hadronicOverEm");
    sigeta   = info->getFloat("sigmaIetaIeta");
    haspix   = info->getBool("hasPixelSeed");
    hiso     = info->getFloat("hcalTowerSumEtConeDR04");
    tiso     = info->getFloat("trkSumPtHollowConeDR04");
    isEB     = info->getBool("isEB");
    pho_maxSieta = isEB? pho_maxSietaEB : pho_maxSietaEE ;

    if( (fabs(eta)<=barrel_maxEta||(fabs(eta)>=endcap_minEta && fabs(eta)<endcap_maxEta)) && hovere<pho_maxHoe ) {

      if( addSigCut==1 && sigeta>pho_maxSieta )  continue;
      if( addSigCut==2 && sigeta<pho_maxSieta )  continue;
      if( noPixSeed==1 && haspix==true)          continue; 
      if( noPixSeed==2 && haspix==false)         continue; 

      if( SpikeCleaning( photon_ ) ) continue;
      bool fcut = sigeta<pho_maxSieta && hiso<pho_maxHcal && tiso<pho_maxTrack;

      if( et>pho1_minEt ) { 
	firstPhos.push_back(ig); 
	firstPacc.push_back(fcut); 
      }
      if( et>pho2_minEt ) { 
	secndPhos.push_back(ig); 
	secndPacc.push_back(fcut); 
      }
    }
  }

  // photon1
  TVector3 nullPt(0,.1,10);
  if( firstPhos.size() == 0 ) return false;
  else {
    tmpPhoton1_ = Candidate::create(nullPt);
    for( int p=0; p<(int)firstPhos.size(); p++ ) {
      Candidate* photon_ = photons_[ firstPhos[p] ];
      if( photon_->pt()>tmpPhoton1_->pt() && firstPacc[p]==true) { 
	tmpPhoton1_ = photon_ ; 
	hasPhoton1 = true; 
      }
    }
    if( !hasPhoton1 ) {
      for( int p=0; p<(int)firstPhos.size(); p++ ) {
	Candidate* photon_ = photons_[ firstPhos[p] ];
	if( photon_->pt()>tmpPhoton1_->pt() ) { 
	  tmpPhoton1_ = photon_ ; 
	  hasPhoton1 = true; 
	}
      }
    }
  }
  if( !hasPhoton1 ) return false;

  // photon2
  if( secndPhos.size() > 0 ) {
    tmpPhoton2_ = Candidate::create(nullPt);
    for( int p=0; p<(int)secndPhos.size(); p++ ) {
      Candidate* photon_ = photons_[ secndPhos[p] ];
      if( photon_->pt()>tmpPhoton2_->pt() && photon_->pt()<tmpPhoton1_->pt() && secndPacc[p]==true) { 
	tmpPhoton2_ = photon_ ; 
	hasPhoton2 = true; 
      }
    }
    if( !hasPhoton2 ) {
      for( int p=0; p<(int)secndPhos.size(); p++ ) {
        Candidate* photon_ = photons_[ secndPhos[p] ];
        if( photon_->pt()>tmpPhoton2_->pt() && photon_->pt()<tmpPhoton1_->pt() ) { 
	  tmpPhoton2_ = photon_ ; 
	  hasPhoton2 = true; 
	}
      }
    }
  }

  // shuffle the photons
  float rndout = random.Integer(2)<0.5;
  if( doShuffle && hasPhoton2 && rndout<0.5 ) {
    thePhoton1_ = tmpPhoton2_;
    thePhoton2_ = tmpPhoton1_;
    isShuffled  = true; 
  }
  else if( hasPhoton2 ) {
    thePhoton1_ = tmpPhoton1_;
    thePhoton2_ = tmpPhoton2_;
    isShuffled   = false;
    if( debug ) cout << "hasPhoton2" << endl;
  }
  else {
    thePhoton1_ = tmpPhoton1_;
  }
  // diphoton selection
  if( hasPhoton1 ) presel=true;
  if( doDiphotons && hasPhoton2==false ) presel=false;
  if( debug ) cout << "Preselection done." << endl;
  return presel;
  
}

void 
DiphotonAnalysis::initRootTree() 
{
  //Event
  outTree_.Branch("Run",               &Run,               "Run/I");
  outTree_.Branch("Event",             &Event,             "Event/I");
  outTree_.Branch("LumiSec",           &LumiSec,           "LumiSec/I");
  outTree_.Branch("nPhotons",          &nPhotons,          "nPhotons/I");  
  outTree_.Branch("nVertex",           &nVertex,           "nVertex/I"); 
  outTree_.Branch("nGoodVertex",       &nGoodVertex,       "nGoodVertex/I"); 
  outTree_.Branch("nJets",             &nJets,             "nJets/I");
  outTree_.Branch("nTracks",           &nTracks,           "nTracks/I");
  outTree_.Branch("nPixSeeds",         &nPixSeeds,         "nPixSeeds/I");
  outTree_.Branch("nClusters",         &nClusters,         "nClusters/I");
  outTree_.Branch("nSC",               &nSC,               "nSC/I");
  outTree_.Branch("nRecHits",          &nRecHits,          "nRecHits/I");
  outTree_.Branch("met",               &met,               "met/F");
  outTree_.Branch("met_phi",           &met_phi,           "met_phi/F");
  outTree_.Branch("mtrans",            &mtrans,            "mtrans/F");
  outTree_.Branch("ptHat",             &ptHat,             "ptHat/F");
  outTree_.Branch("isShuffled",        &isShuffled,        "isShuffled/B");
  
  //Vertex
  outTree_.Branch("vX",             &vX,         "vX/F");
  outTree_.Branch("vY",             &vY,         "vY/F");
  outTree_.Branch("vZ",             &vZ,         "vZ/F");
  outTree_.Branch("vNtracks",       &vNtracks,   "vNtracks/I");
  outTree_.Branch("vNout",          &vNout,      "vNout/I");
  outTree_.Branch("vIsGood",        &vIsGood,    "vIsGood/B");

  //Trigger bits
  outTree_.Branch("HLT_Photon10_Cleaned_L1R",     &HLT_Photon10_Cleaned_L1R,   "HLT_Photon10_Cleaned_L1R/B");
  outTree_.Branch("HLT_Photon15_Cleaned_L1R",     &HLT_Photon15_Cleaned_L1R,   "HLT_Photon15_Cleaned_L1R/B");
  outTree_.Branch("HLT_Photon17_Isol_SC17HE_L1R_v1",&HLT_Photon17_Isol_SC17HE_L1R_v1,"HLT_Photon17_Isol_SC17HE_L1R_v1/B");
  outTree_.Branch("HLT_Photon17_SC17HE_L1R_v1",   &HLT_Photon17_SC17HE_L1R_v1, "HLT_Photon17_SC17HE_L1R_v1/B");
  outTree_.Branch("HLT_Photon20_Cleaned_L1R",     &HLT_Photon20_Cleaned_L1R,   "HLT_Photon20_Cleaned_L1R/B");
  outTree_.Branch("HLT_Photon20_Isol_Cleaned_L1R_v1",&HLT_Photon20_Isol_Cleaned_L1R_v1,"HLT_Photon20_Isol_Cleaned_L1R_v1/B");
  outTree_.Branch("HLT_Photon22_SC22HE_L1R",      &HLT_Photon22_SC22HE_L1R,    "HLT_Photon22_SC22HE_L1R/B");
  outTree_.Branch("HLT_Photon25_Cleaned_L1R",     &HLT_Photon25_Cleaned_L1R,   "HLT_Photon25_Cleaned_L1R/B");
  outTree_.Branch("HLT_Photon30_Cleaned_L1R",     &HLT_Photon30_Cleaned_L1R,   "HLT_Photon30_Cleaned_L1R/B");
  outTree_.Branch("HLT_Photon35_Isol_Cleaned_L1R_v1",&HLT_Photon35_Isol_Cleaned_L1R_v1,"HLT_Photon35_Isol_Cleaned_L1R_v1/B");
  outTree_.Branch("HLT_Photon40_Isol_Cleaned_L1R_v1",&HLT_Photon40_Isol_Cleaned_L1R_v1,"HLT_Photon40_Isol_Cleaned_L1R_v1/B");
  outTree_.Branch("HLT_Photon40_CaloId_Cleaned_L1R_v1",&HLT_Photon40_CaloId_Cleaned_L1R_v1,"HLT_Photon40_CaloId_Cleaned_L1R_v1/B");
  outTree_.Branch("HLT_Photon50_Cleaned_L1R",     &HLT_Photon50_Cleaned_L1R,   "HLT_Photon50_Cleaned_L1R/B");
  outTree_.Branch("HLT_Photon70_Cleaned_L1R",     &HLT_Photon70_Cleaned_L1R,   "HLT_Photon70_Cleaned_L1R/B");
  outTree_.Branch("HLT_DoublePhoton5_L1R",        &HLT_DoublePhoton5_L1R,      "HLT_DoublePhoton5_L1R/B");
  outTree_.Branch("HLT_DoublePhoton10_L1R",       &HLT_DoublePhoton10_L1R,     "HLT_DoublePhoton10_L1R/B");
  outTree_.Branch("HLT_DoublePhoton15_L1R",       &HLT_DoublePhoton15_L1R,     "HLT_DoublePhoton15_L1R/B");
  outTree_.Branch("HLT_DoublePhoton17_L1R",       &HLT_DoublePhoton17_L1R,     "HLT_DoublePhoton17_L1R/B");
  outTree_.Branch("HLT_DoublePhoton17_SingleIsol_L1R_v1",&HLT_DoublePhoton17_SingleIsol_L1R_v1,"HLT_DoublePhoton17_SingleIsol_L1R_v1/B");
  outTree_.Branch("HLT_DoublePhoton20_L1R",       &HLT_DoublePhoton20_L1R,     "HLT_DoublePhoton20_L1R/B");
  outTree_.Branch("HLT_DoublePhoton22_L1R_v1",    &HLT_DoublePhoton22_L1R_v1,  "HLT_DoublePhoton22_L1R_v1/B");

  //Photon1
  outTree_.Branch("p1",             &p1,         "p1/F");
  outTree_.Branch("px1",            &px1,        "px1/F");
  outTree_.Branch("py1",            &py1,        "py1/F");
  outTree_.Branch("pz1",            &pz1,        "pz1/F");
  outTree_.Branch("et1",            &et1,        "et1/F");
  outTree_.Branch("energy1",        &energy1,    "energy1/F");
  outTree_.Branch("eta1",           &eta1,       "eta1/F");
  outTree_.Branch("phi1",           &phi1,       "phi1/F");
  outTree_.Branch("isEBGap1",       &isEBGap1,   "isEBGap1/B"); 
  outTree_.Branch("isEEGap1",       &isEEGap1,   "isEEGap1/B"); 
  outTree_.Branch("isEBEEGap1",     &isEBEEGap1, "isEBEEGap1/B"); 
  outTree_.Branch("isEB1",          &isEB1,      "isEB1/B"); 
  outTree_.Branch("isEE1",          &isEE1,      "isEE1/B"); 
  
  outTree_.Branch("sigmaEtaEta1",   &sigmaEtaEta1,   "sigmaEtaEta1/F");
  outTree_.Branch("sigmaIetaIeta1", &sigmaIetaIeta1, "sigmaIetaIeta1/F");
  outTree_.Branch("e1x51",          &e1x51,          "e1x51/F");
  outTree_.Branch("e2x51",          &e2x51,          "e2x51/F");
  outTree_.Branch("ecalRecHitSumEtConeDR031",        &ecalRecHitSumEtConeDR031,       "ecalRecHitSumEtConeDR031/F");
  outTree_.Branch("ecalRecHitSumEtConeDR041",        &ecalRecHitSumEtConeDR041,       "ecalRecHitSumEtConeDR041/F");
  outTree_.Branch("hadronicDepth1OverEm1",           &hadronicDepth1OverEm1,          "hadronicDepth1OverEm1/F");
  outTree_.Branch("hadronicDepth2OverEm1",           &hadronicDepth2OverEm1,          "hadronicDepth2OverEm1/F");
  outTree_.Branch("hadronicOverEm1",                 &hadronicOverEm1,                "hadronicOverEm1/F");
  outTree_.Branch("e2hcalDepth1TowerSumEtConeDR031", &hcalDepth1TowerSumEtConeDR031,  "hcalDepth1TowerSumEtConeDR031/F");
  outTree_.Branch("hcalDepth1TowerSumEtConeDR041",   &hcalDepth1TowerSumEtConeDR041,  "hcalDepth1TowerSumEtConeDR041/F");
  outTree_.Branch("hcalDepth2TowerSumEtConeDR031",   &hcalDepth2TowerSumEtConeDR031,  "hcalDepth2TowerSumEtConeDR031/F");
  outTree_.Branch("hcalDepth2TowerSumEtConeDR041",   &hcalDepth2TowerSumEtConeDR041,  "hcalDepth2TowerSumEtConeDR041/F");
  outTree_.Branch("hcalTowerSumEtConeDR031",         &hcalTowerSumEtConeDR031,        "hcalTowerSumEtConeDR031/F");
  outTree_.Branch("hcalTowerSumEtConeDR041",         &hcalTowerSumEtConeDR041,        "hcalTowerSumEtConeDR041/F");
  outTree_.Branch("nTrkHollowConeDR031",             &nTrkHollowConeDR031,            "nTrkHollowConeDR031/F");
  outTree_.Branch("nTrkHollowConeDR041",             &nTrkHollowConeDR041,            "nTrkHollowConeDR041/F");
  outTree_.Branch("nTrkSolidConeDR031",              &nTrkSolidConeDR031,             "nTrkSolidConeDR031/F");
  outTree_.Branch("nTrkSolidConeDR041",              &nTrkSolidConeDR041,             "nTrkSolidConeDR041/F");
  outTree_.Branch("trkSumPtHollowConeDR031",         &trkSumPtHollowConeDR031,        "trkSumPtHollowConeDR031/F");
  outTree_.Branch("trkSumPtHollowConeDR041",         &trkSumPtHollowConeDR041,        "trkSumPtHollowConeDR041/F");
  outTree_.Branch("trkSumPtSolidConeDR031",          &trkSumPtSolidConeDR031,         "trkSumPtSolidConeDR031/F");
  outTree_.Branch("trkSumPtSolidConeDR041",          &trkSumPtSolidConeDR041,         "trkSumPtSolidConeDR041/F");
  outTree_.Branch("patSumIso1",                      &patSumIso1,                     "patSumIso1/F");
  outTree_.Branch("ecalIso1",                        &ecalIso1,                       "ecalIso1/F");
  outTree_.Branch("ecalIso1_280",                    &ecalIso1_280,                   "ecalIso1_280/F");
  outTree_.Branch("ecalIso1_320",                    &ecalIso1_320,                   "ecalIso1_320/F");
  outTree_.Branch("ecalIso1_350",                    &ecalIso1_350,                   "ecalIso1_350/F");
  outTree_.Branch("ecalIso1_400",                    &ecalIso1_400,                   "ecalIso1_400/F");
  outTree_.Branch("ecalIso1_n",                      &ecalIso1_n,                     "ecalIso1_n/I");
  outTree_.Branch("ecalIso1_ntk",                    &ecalIso1_ntk,                   "ecalIso1_ntk/I");
  outTree_.Branch("ecalIso1_nTkGinV",                &ecalIso1_nTkGinV,               "ecalIso1_nTkGinV/I");
  outTree_.Branch("ecalIso1_nTkG",                   &ecalIso1_nTkG,                  "ecalIso1_nTkG/I");
  outTree_.Branch("ecalIso1_nTkG_28",                &ecalIso1_nTkG_28,               "ecalIso1_nTkG_28/I");
  outTree_.Branch("ecalIso1_nTkG_32",                &ecalIso1_nTkG_32,               "ecalIso1_nTkG_32/I");
  outTree_.Branch("ecalIso1_nTkG_40",                &ecalIso1_nTkG_40,               "ecalIso1_nTkG_40/I");
  outTree_.Branch("ecalIso1_nTkGnoZ",                &ecalIso1_nTkGnoZ,               "ecalIso1_nTkGnoZ/I");
  outTree_.Branch("ecalIso1_Vieta",                  &ecalIso1_Vieta);
  outTree_.Branch("ecalIso1_Veta",                   &ecalIso1_Veta);
  outTree_.Branch("ecalIso1_Viphi",                  &ecalIso1_Viphi);
  outTree_.Branch("ecalIso1_Vphi",                   &ecalIso1_Vphi);
  outTree_.Branch("ecalIso1_Vet",                    &ecalIso1_Vet);
  outTree_.Branch("ecalIso1_Ve",                     &ecalIso1_Ve);
  outTree_.Branch("ecalIso1_Vieta_nv",               &ecalIso1_Vieta_nv);
  outTree_.Branch("ecalIso1_Veta_nv",                &ecalIso1_Veta_nv);
  outTree_.Branch("ecalIso1_Viphi_nv",               &ecalIso1_Viphi_nv);
  outTree_.Branch("ecalIso1_Vphi_nv",                &ecalIso1_Vphi_nv);
  outTree_.Branch("ecalIso1_Vet_nv",                 &ecalIso1_Vet_nv);
  outTree_.Branch("ecalIso1_Ve_nv",                  &ecalIso1_Ve_nv);
  outTree_.Branch("ecalIso1_Vtkpt",                  &ecalIso1_Vtkpt);
  outTree_.Branch("ecalIso1_VtkCnt",                 &ecalIso1_VtkCnt);
  outTree_.Branch("ecalIso1_Vtkd0",                  &ecalIso1_Vtkd0);
  outTree_.Branch("ecalIso1_Vtkz0",                  &ecalIso1_Vtkz0);
  outTree_.Branch("ecalIso1_VtkEtaH",                &ecalIso1_VtkEtaH);
  outTree_.Branch("ecalIso1_VtkPhiH",                &ecalIso1_VtkPhiH);
  outTree_.Branch("ecalIso1_VtkDEta",                &ecalIso1_VtkEtaH);
  outTree_.Branch("ecalIso1_VtkDPhi",                &ecalIso1_VtkEtaH);
  outTree_.Branch("ecalIso1_VtkDR",                  &ecalIso1_VtkEtaH);
  outTree_.Branch("ecalIso1_VtkNh",                  &ecalIso1_VtkNh);
  outTree_.Branch("ecalIso1_VtkIh",                  &ecalIso1_VtkIh);
  outTree_.Branch("ecalIso1_VtkChi2",                &ecalIso1_VtkChi2);
  outTree_.Branch("ecalIso1_VtkNdof",                &ecalIso1_VtkNdof);
  outTree_.Branch("ecalIso1_VtkQ",                   &ecalIso1_VtkQ);
  outTree_.Branch("ecalIso1_VtkFP",                  &ecalIso1_VtkFP);
  outTree_.Branch("ecalIso1_VtkSV",                  &ecalIso1_VtkSV);
  outTree_.Branch("ecalIso1_VtkInV",                 &ecalIso1_VtkInV);
  outTree_.Branch("ecalIso1_VtkVid",                 &ecalIso1_VtkVid);
  outTree_.Branch("ecalIso1_VtkDepE",                &ecalIso1_VtkDepE);
  outTree_.Branch("ecalIso1_VtkDepEt",               &ecalIso1_VtkDepEt);
  outTree_.Branch("hcalIso1",                        &hcalIso1,                       "hcalIso1/F");
  outTree_.Branch("hcalIso1_n",                      &hcalIso1_n,                     "hcalIso1_n/I");
  outTree_.Branch("hcalIsoH1",                       &hcalIsoH1,                      "hcalIsoH1/F");
  outTree_.Branch("hcalIsoH1_n",                     &hcalIsoH1_n,                    "hcalIsoH1_n/I");
  outTree_.Branch("trackIso1",                       &trackIso1,                      "trackIso1/F");
  outTree_.Branch("trackIsoRd1",                     &trackIsoRd1,                    "trackIsoRd1/F");
  outTree_.Branch("trackIso1_n",                     &trackIso1_n,                    "trackIso1_n/I");
  outTree_.Branch("sumIso1",                         &sumIso1,                        "sumIso1/F");

  outTree_.Branch("scPhi1",     &scPhi1,     "scPhi1/F");
  outTree_.Branch("scEta1",     &scEta1,     "scEta1/F");
  outTree_.Branch("scRawE1",    &scRawE1,    "scRawE1/F");
  outTree_.Branch("nXtalsSC1",  &nXtalsSC1,  "nXtalsSC1/I");
  outTree_.Branch("scPhiWidth1",&scPhiWidth1,"scPhiWidth1/F");
  outTree_.Branch("scEtaWidth1",&scEtaWidth1,"scEtaWidth1/F");
  outTree_.Branch("scIxSeed1",  &scIxSeed1,  "scIxSeed1/I"); 
  outTree_.Branch("scIySeed1",  &scIySeed1,  "scIySeed1/I"); 
  outTree_.Branch("scIzSeed1",  &scIzSeed1,  "scIzSeed1/I"); 
  outTree_.Branch("r91",        &r91,        "r91/F");
  outTree_.Branch("fBrem1",     &fBrem1,     "fBrem1/F");
  outTree_.Branch("r41",        &r41,        "r41/F");

  outTree_.Branch("isConv1",    &isConv1,    "isConv1/B"); 
  outTree_.Branch("hasPixSeed1",&hasPixSeed1,"hasPixSeed1/B"); 
  outTree_.Branch("inCrack1",   &inCrack1,   "inCrack1/B");
  
  outTree_.Branch("hasGenMatch1",&hasGenMatch1,"hasGenMatch1/B"); 
  outTree_.Branch("genDR1",      &genDR1,      "genDR1/F"); 
  outTree_.Branch("genEnergy1",  &genEnergy1,  "genEnergy1/F"); 
  outTree_.Branch("genEta1",     &genEta1,     "genEta1/F"); 
  outTree_.Branch("genPhi1",     &genPhi1,     "genPhi1/F"); 
  outTree_.Branch("genEt1",      &genEt1,      "genEt1/F"); 
  outTree_.Branch("genIso1",     &genIso1,     "genIso1/F"); 
  outTree_.Branch("genIsoAll1",  &genIsoAll1,  "genIsoAll1/F"); 
  outTree_.Branch("genMother1",  &genMother1,  "genMother1/I");
  outTree_.Branch("isFrag1",     &isFrag1,     "isFrag1/B"); 
  
  outTree_.Branch("hasJetMatch1",&hasJetMatch1,"hasJetMatch1/B"); 
  outTree_.Branch("jetDR1",      &jetDR1,      "jetDR1/F"); 
  outTree_.Branch("jetEnergy1",  &jetEnergy1,  "jetEnergy1/F"); 
  outTree_.Branch("jetEta1",     &jetEta1,     "jetEta1/F"); 
  outTree_.Branch("jetPhi1",     &jetPhi1,     "jetPhi1/F"); 
  outTree_.Branch("jetEt1",      &jetEt1,      "jetEt1/F"); 
  
//   outTree_.Branch("preshE1",    &preshE1,    "preshE1/F"); // handle
//   outTree_.Branch("nESclu1",    &nESclu1,    "nESclu1/F"); // handle
  
  //Photon2
  outTree_.Branch("hasPhoton2",     &hasPhoton2, "hasPhoton2/B");
  
  outTree_.Branch("p2",             &p2,         "p2/F");
  outTree_.Branch("px2",            &px2,        "px2/F");
  outTree_.Branch("py2",            &py2,        "py2/F");
  outTree_.Branch("pz2",            &pz2,        "pz2/F");
  outTree_.Branch("et2",            &et2,        "et2/F");
  outTree_.Branch("energy2",        &energy2,    "energy2/F");
  outTree_.Branch("eta2",           &eta2,       "eta2/F");
  outTree_.Branch("phi2",           &phi2,       "phi2/F");
  outTree_.Branch("isEBGap2",       &isEBGap2,   "isEBGap2/B"); 
  outTree_.Branch("isEEGap2",       &isEEGap2,   "isEEGap2/B"); 
  outTree_.Branch("isEBEEGap2",     &isEBEEGap2, "isEBEEGap2/B"); 
  outTree_.Branch("isEB2",          &isEB2,      "isEB2/B"); 
  outTree_.Branch("isEE2",          &isEE2,      "isEE2/B"); 
  
  outTree_.Branch("sigmaEtaEta2",   &sigmaEtaEta2,   "sigmaEtaEta2/F");
  outTree_.Branch("sigmaIetaIeta2", &sigmaIetaIeta2, "sigmaIetaIeta2/F");
  outTree_.Branch("e1x52",          &e1x52,          "e1x52/F");
  outTree_.Branch("e2x52",          &e2x52,          "e2x52/F");
  outTree_.Branch("ecalRecHitSumEtConeDR032",        &ecalRecHitSumEtConeDR032,      "ecalRecHitSumEtConeDR032/F");
  outTree_.Branch("ecalRecHitSumEtConeDR042",        &ecalRecHitSumEtConeDR042,      "ecalRecHitSumEtConeDR042/F");
  outTree_.Branch("hadronicDepth1OverEm2",           &hadronicDepth1OverEm2,         "hadronicDepth1OverEm2/F");
  outTree_.Branch("hadronicDepth2OverEm2",           &hadronicDepth2OverEm2,         "hadronicDepth2OverEm2/F");
  outTree_.Branch("hadronicOverEm2",                 &hadronicOverEm2,               "hadronicOverEm2/F");
  outTree_.Branch("e2hcalDepth1TowerSumEtConeDR032", &hcalDepth1TowerSumEtConeDR032, "hcalDepth1TowerSumEtConeDR032/F");
  outTree_.Branch("hcalDepth1TowerSumEtConeDR042",   &hcalDepth1TowerSumEtConeDR042, "hcalDepth1TowerSumEtConeDR042/F");
  outTree_.Branch("hcalDepth2TowerSumEtConeDR032",   &hcalDepth2TowerSumEtConeDR032, "hcalDepth2TowerSumEtConeDR032/F");
  outTree_.Branch("hcalDepth2TowerSumEtConeDR042",   &hcalDepth2TowerSumEtConeDR042, "hcalDepth2TowerSumEtConeDR042/F");
  outTree_.Branch("hcalTowerSumEtConeDR032",         &hcalTowerSumEtConeDR032,       "hcalTowerSumEtConeDR032/F");
  outTree_.Branch("hcalTowerSumEtConeDR042",         &hcalTowerSumEtConeDR042,       "hcalTowerSumEtConeDR042/F");
  outTree_.Branch("nTrkHollowConeDR032",             &nTrkHollowConeDR032,           "nTrkHollowConeDR032/F");
  outTree_.Branch("nTrkHollowConeDR042",             &nTrkHollowConeDR042,           "nTrkHollowConeDR042/F");
  outTree_.Branch("nTrkSolidConeDR032",              &nTrkSolidConeDR032,            "nTrkSolidConeDR032/F");
  outTree_.Branch("nTrkSolidConeDR042",              &nTrkSolidConeDR042,            "nTrkSolidConeDR042/F");
  outTree_.Branch("trkSumPtHollowConeDR032",         &trkSumPtHollowConeDR032,       "trkSumPtHollowConeDR032/F");
  outTree_.Branch("trkSumPtHollowConeDR042",         &trkSumPtHollowConeDR042,       "trkSumPtHollowConeDR042/F");
  outTree_.Branch("trkSumPtSolidConeDR032",          &trkSumPtSolidConeDR032,        "trkSumPtSolidConeDR032/F");
  outTree_.Branch("trkSumPtSolidConeDR042",          &trkSumPtSolidConeDR042,        "trkSumPtSolidConeDR042/F");
  outTree_.Branch("patSumIso2",                      &patSumIso2,                    "patSumIso2/F");
  outTree_.Branch("ecalIso2",                        &ecalIso2,                      "ecalIso2/F");
  outTree_.Branch("ecalIso2_280",                    &ecalIso2_280,                   "ecalIso2_280/F");
  outTree_.Branch("ecalIso2_320",                    &ecalIso2_320,                   "ecalIso2_320/F");
  outTree_.Branch("ecalIso2_350",                    &ecalIso2_350,                   "ecalIso2_350/F");
  outTree_.Branch("ecalIso2_400",                    &ecalIso2_400,                   "ecalIso2_400/F");
  outTree_.Branch("ecalIso2_n",                      &ecalIso2_n,                    "ecalIso2_n/I");
  outTree_.Branch("ecalIso2_ntk",                    &ecalIso2_ntk,                  "ecalIso2_ntk/I");
  outTree_.Branch("ecalIso2_nTkGinV",                &ecalIso2_nTkGinV,              "ecalIso2_nTkGinV/I");
  outTree_.Branch("ecalIso2_nTkG",                   &ecalIso2_nTkG,                 "ecalIso2_nTkG/I");
  outTree_.Branch("ecalIso2_nTkG_28",                &ecalIso2_nTkG_28,               "ecalIso2_nTkG_28/I");
  outTree_.Branch("ecalIso2_nTkG_32",                &ecalIso2_nTkG_32,               "ecalIso2_nTkG_32/I");
  outTree_.Branch("ecalIso2_nTkG_40",                &ecalIso2_nTkG_40,               "ecalIso2_nTkG_40/I");
  outTree_.Branch("ecalIso2_nTkGnoZ",                &ecalIso2_nTkGnoZ,              "ecalIso2_nTkGnoZ/I");
  outTree_.Branch("ecalIso2_Vieta",                  &ecalIso2_Vieta);
  outTree_.Branch("ecalIso2_Veta",                   &ecalIso2_Veta);
  outTree_.Branch("ecalIso2_Viphi",                  &ecalIso2_Viphi);
  outTree_.Branch("ecalIso2_Vphi",                   &ecalIso2_Vphi);
  outTree_.Branch("ecalIso2_Vet",                    &ecalIso2_Vet);
  outTree_.Branch("ecalIso2_Ve",                     &ecalIso2_Ve);
  outTree_.Branch("ecalIso2_Vieta_nv",               &ecalIso2_Vieta_nv);
  outTree_.Branch("ecalIso2_Veta_nv",                &ecalIso2_Veta_nv);
  outTree_.Branch("ecalIso2_Viphi_nv",               &ecalIso2_Viphi_nv);
  outTree_.Branch("ecalIso2_Vphi_nv",                &ecalIso2_Vphi_nv);
  outTree_.Branch("ecalIso2_Vet_nv",                 &ecalIso2_Vet_nv);
  outTree_.Branch("ecalIso2_Ve_nv",                  &ecalIso2_Ve_nv);
  outTree_.Branch("ecalIso2_Vtkpt",                  &ecalIso2_Vtkpt);
  outTree_.Branch("ecalIso2_VtkCnt",                 &ecalIso2_VtkCnt);
  outTree_.Branch("ecalIso2_Vtkd0",                  &ecalIso2_Vtkd0);
  outTree_.Branch("ecalIso2_Vtkz0",                  &ecalIso2_Vtkz0);
  outTree_.Branch("ecalIso2_VtkEtaH",                &ecalIso2_VtkEtaH);
  outTree_.Branch("ecalIso2_VtkPhiH",                &ecalIso2_VtkPhiH);
  outTree_.Branch("ecalIso2_VtkDEta",                &ecalIso2_VtkEtaH);
  outTree_.Branch("ecalIso2_VtkDPhi",                &ecalIso2_VtkEtaH);
  outTree_.Branch("ecalIso2_VtkDR",                  &ecalIso2_VtkEtaH);
  outTree_.Branch("ecalIso2_VtkNh",                  &ecalIso2_VtkNh);
  outTree_.Branch("ecalIso2_VtkIh",                  &ecalIso2_VtkIh);
  outTree_.Branch("ecalIso2_VtkChi2",                &ecalIso2_VtkChi2);
  outTree_.Branch("ecalIso2_VtkNdof",                &ecalIso2_VtkNdof);
  outTree_.Branch("ecalIso2_VtkQ",                   &ecalIso2_VtkQ);
  outTree_.Branch("ecalIso2_VtkFP",                  &ecalIso2_VtkFP);
  outTree_.Branch("ecalIso2_VtkSV",                  &ecalIso2_VtkSV);
  outTree_.Branch("ecalIso2_VtkInV",                 &ecalIso2_VtkInV);
  outTree_.Branch("ecalIso2_VtkVid",                 &ecalIso2_VtkVid);
  outTree_.Branch("ecalIso2_VtkDepE",                &ecalIso2_VtkDepE);
  outTree_.Branch("ecalIso2_VtkDepEt",               &ecalIso2_VtkDepEt);
  outTree_.Branch("eIso2_RECO",                      &eIso2_RECO,                    "eIso2_RECO/F");
  outTree_.Branch("hcalIso2",                        &hcalIso2,                      "hcalIso2/F");
  outTree_.Branch("hcalIsoH2",                       &hcalIsoH2,                     "hcalIsoH2/F");
  outTree_.Branch("trackIso2",                       &trackIso2,                     "trackIso2/F");
  outTree_.Branch("trackIsoRd2",                     &trackIsoRd2,                   "trackIsoRd2/F");
  outTree_.Branch("sumIso2",                         &sumIso2,                       "sumIso2/F");

  outTree_.Branch("scPhi2",     &scPhi2,     "scPhi2/F");
  outTree_.Branch("scEta2",     &scEta2,     "scEta2/F");
  outTree_.Branch("scRawE2",    &scRawE2,    "scRawE2/F");
  outTree_.Branch("nXtalsSC2",  &nXtalsSC2,  "nXtalsSC2/I");
  outTree_.Branch("scPhiWidth2",&scPhiWidth2,"scPhiWidth2/F");
  outTree_.Branch("scEtaWidth2",&scEtaWidth2,"scEtaWidth2/F");
  outTree_.Branch("scIxSeed2",  &scIxSeed2,  "scIxSeed2/I"); 
  outTree_.Branch("scIySeed2",  &scIySeed2,  "scIySeed2/I"); 
  outTree_.Branch("scIzSeed2",  &scIzSeed2,  "scIzSeed2/I"); 
  outTree_.Branch("r92",        &r92,        "r92/F");
  outTree_.Branch("fBrem2",     &fBrem2,     "fBrem2/F");
  outTree_.Branch("r42",        &r42,        "r42/F");

  outTree_.Branch("isConv2",    &isConv2,    "isConv2/B"); 
  outTree_.Branch("hasPixSeed2",&hasPixSeed2,"hasPixSeed2/B"); 
  outTree_.Branch("inCrack2",   &inCrack2,   "inCrack2/B");
  
  outTree_.Branch("hasGenMatch2",&hasGenMatch2,"hasGenMatch2/B"); 
  outTree_.Branch("genDR2",      &genDR2,      "genDR2/F"); 
  outTree_.Branch("genEnergy2",  &genEnergy2,  "genEnergy2/F"); 
  outTree_.Branch("genEta2",     &genEta2,     "genEta2/F"); 
  outTree_.Branch("genPhi2",     &genPhi2,     "genPhi2/F"); 
  outTree_.Branch("genEt2",      &genEt2,      "genEt2/F"); 
  outTree_.Branch("genIso2",     &genIso2,     "genIso2/F"); 
  outTree_.Branch("genIsoAll2",  &genIsoAll2,  "genIsoAll2/F"); 
  outTree_.Branch("genMother2",  &genMother2,  "genMother2/I"); 
  outTree_.Branch("isFrag2",     &isFrag2,     "isFrag2/B"); 
  
  outTree_.Branch("hasJetMatch2",&hasJetMatch2,"hasJetMatch2/B"); 
  outTree_.Branch("jetDR2",      &jetDR2,      "jetDR2/F"); 
  outTree_.Branch("jetEnergy2",  &jetEnergy2,  "jetEnergy2/F"); 
  outTree_.Branch("jetEta2",     &jetEta2,     "jetEta2/F"); 
  outTree_.Branch("jetPhi2",     &jetPhi2,     "jetPhi2/F"); 
  outTree_.Branch("jetEt2",      &jetEt2,      "jetEt2/F"); 
  
  //Diphoton
  outTree_.Branch("diphoMass",   &diphoMass,   "diphoMass/F"); 
  outTree_.Branch("diphoEta",    &diphoEta,    "diphoEta/F"); 
  outTree_.Branch("diphoY",      &diphoY,      "diphoY/F"); 
  outTree_.Branch("deltaY",      &deltaY,      "deltaY/F"); 
  outTree_.Branch("diphoQt",     &diphoQt,     "diphoQt/F"); 
  outTree_.Branch("diphoPhi",    &diphoPhi,    "diphoPhi/F"); 
  outTree_.Branch("deltaEta",    &deltaEta,    "deltaEta/F"); 
  outTree_.Branch("deltaPhi",    &deltaPhi,    "deltaPhi/F"); 
  outTree_.Branch("deltaR",      &deltaR,      "deltaR/F");  
  outTree_.Branch("cosTheta",    &cosTheta,    "cosTheta/F");
  outTree_.Branch("cosThetaH",   &cosThetaH,   "cosThetaH/F");
  outTree_.Branch("cosThetaS",   &cosThetaS,   "cosThetaS/F");
  outTree_.Branch("cosThetaP",   &cosThetaP,   "cosThetaP/F");

  //Hardest Jet
  outTree_.Branch("hardestJetPt",  &hardestJetPt,        "hardestJetPt/F");
  outTree_.Branch("secondJetPt",   &secondJetPt,         "secondJetPt/F");
  
  //Random Cone
  outTree_.Branch("cc1_eta",                 &cc1_eta,      "cc1_eta/F");
  outTree_.Branch("cc1_phi",                 &cc1_phi,      "cc1_phi/F");
  outTree_.Branch("cc1_hcalIso",             &cc1_hcalIso,  "cc1_hcalIso/F");
  outTree_.Branch("cc1_Nhcal",               &cc1_Nhcal,    "cc1_Nhcal/I");
  outTree_.Branch("cc1_trackIso",            &cc1_trackIso, "cc1_trackIso/F");
  outTree_.Branch("cc1_Ntrack",              &cc1_Ntrack,   "cc1_Ntrack/I");
  outTree_.Branch("cc1_sumIso",              &cc1_sumIso,   "cc1_sumIso/F");
  outTree_.Branch("cc1_hasJet",              &cc1_hasJet,   "cc1_hasJet/B");

  outTree_.Branch("cc1_ecalIso",             &cc1_ecalIso,    "cc1_ecalIso/F");
  outTree_.Branch("cc1_ecalIso_n",           &cc1_ecalIso_n,  "cc1_ecalIso_n/I");
  outTree_.Branch("cc1_ecalIso_ntk",         &cc1_ecalIso_ntk,"cc1_ecalIso_n/I");
  outTree_.Branch("cc1_ecalIso_nTkG",        &cc1_ecalIso_nTkG,"cc1_ecalIso_nTkG/I");
  outTree_.Branch("cc1_ecalIso_Vieta",       &cc1_ecalIso_Vieta);
  outTree_.Branch("cc1_ecalIso_Veta",        &cc1_ecalIso_Veta);
  outTree_.Branch("cc1_ecalIso_Viphi",       &cc1_ecalIso_Viphi);
  outTree_.Branch("cc1_ecalIso_Vphi",        &cc1_ecalIso_Vphi);
  outTree_.Branch("cc1_ecalIso_Vet",         &cc1_ecalIso_Vet);
  outTree_.Branch("cc1_ecalIso_Ve",          &cc1_ecalIso_Ve);
  outTree_.Branch("cc1_ecalIso_Vtkpt",       &cc1_ecalIso_Vtkpt);
  outTree_.Branch("cc1_ecalIso_Vtkd0",       &cc1_ecalIso_Vtkd0);
  outTree_.Branch("cc1_ecalIso_Vtkz0",       &cc1_ecalIso_Vtkz0);
  outTree_.Branch("cc1_ecalIso_VtkNh",       &cc1_ecalIso_VtkNh);
  outTree_.Branch("cc1_ecalIso_VtkIh",       &cc1_ecalIso_VtkIh);
  outTree_.Branch("cc1_ecalIso_VtkChi2",     &cc1_ecalIso_VtkChi2);
  outTree_.Branch("cc1_ecalIso_VtkNdof",     &cc1_ecalIso_VtkNdof);
  outTree_.Branch("cc1_ecalIso_VtkQ",        &cc1_ecalIso_VtkQ);
  outTree_.Branch("cc1_ecalIso_VtkFP",       &cc1_ecalIso_VtkFP);
  outTree_.Branch("cc1_ecalIso_VtkSV",       &cc1_ecalIso_VtkSV);
  outTree_.Branch("cc1_ecalIso_VtkEtaH",     &cc1_ecalIso_VtkEtaH);
  outTree_.Branch("cc1_ecalIso_VtkPhiH",     &cc1_ecalIso_VtkPhiH);
  outTree_.Branch("cc1_ecalIso_VtkInV",      &cc1_ecalIso_VtkInV);
  outTree_.Branch("cc1_ecalIso_VtkVid",      &cc1_ecalIso_VtkVid);
  outTree_.Branch("cc1_ecalIso_VtkDepE",     &cc1_ecalIso_VtkDepE);
  outTree_.Branch("cc1_ecalIso_VtkDepEt",    &cc1_ecalIso_VtkDepEt);

  outTree_.Branch("cc2_eta",     &cc2_eta,     "cc2_eta/F");
  outTree_.Branch("cc2_phi",     &cc2_phi,     "cc2_phi/F");
  outTree_.Branch("cc2_ecalIso", &cc2_ecalIso, "cc2_ecalIso/F");
  outTree_.Branch("cc2_Necal",   &cc2_Necal,   "cc2_Necal/I");
  outTree_.Branch("cc2_hcalIso", &cc2_hcalIso, "cc2_hcalIso/F");
  outTree_.Branch("cc2_Nhcal",   &cc2_Nhcal,   "cc2_Nhcal/I");
  outTree_.Branch("cc2_trackIso",&cc2_trackIso,"cc2_trackIso/F");
  outTree_.Branch("cc2_Ntrack",  &cc2_Ntrack,  "cc2_Ntrack/I");
  outTree_.Branch("cc2_sumIso",  &cc2_sumIso,  "cc2_sumIso/F");
  outTree_.Branch("cc2_hasJet",  &cc2_hasJet,  "cc2_hasJet/B");
  
  outTree_.Branch("cc3_eta",     &cc3_eta,     "cc3_eta/F");
  outTree_.Branch("cc3_phi",     &cc3_phi,     "cc3_phi/F");
  outTree_.Branch("cc3_ecalIso", &cc3_ecalIso, "cc3_ecalIso/F");
  outTree_.Branch("cc3_Necal",   &cc3_Necal,   "cc3_Necal/I");
  outTree_.Branch("cc3_hcalIso", &cc3_hcalIso, "cc3_hcalIso/F");
  outTree_.Branch("cc3_Nhcal",   &cc3_Nhcal,   "cc3_Nhcal/I");
  outTree_.Branch("cc3_trackIso",&cc3_trackIso,"cc3_trackIso/F");
  outTree_.Branch("cc3_Ntrack",  &cc3_Ntrack,  "cc3_Ntrack/I");
  outTree_.Branch("cc3_sumIso",  &cc3_sumIso,  "cc3_sumIso/F");
  outTree_.Branch("cc3_hasJet",  &cc3_hasJet,  "cc3_hasJet/B");
  outTree_.Branch("cc3_Vieta",   &cc3_Vieta);
  outTree_.Branch("cc3_Viphi",   &cc3_Viphi);
  outTree_.Branch("cc3_Vet",     &cc3_Vet);
  outTree_.Branch("cc3_Ve",      &cc3_Ve);
  
  outTree_.Branch("cc4_eta",     &cc4_eta,     "cc4_eta/F");
  outTree_.Branch("cc4_phi",     &cc4_phi,     "cc4_phi/F");
  outTree_.Branch("cc4_ecalIso", &cc4_ecalIso, "cc4_ecalIso/F");
  outTree_.Branch("cc4_Necal",   &cc4_Necal,   "cc4_Necal/I");
  outTree_.Branch("cc4_hcalIso", &cc4_hcalIso, "cc4_hcalIso/F");
  outTree_.Branch("cc4_Nhcal",   &cc4_Nhcal,   "cc4_Nhcal/I");
  outTree_.Branch("cc4_trackIso",&cc4_trackIso,"cc4_trackIso/F");
  outTree_.Branch("cc4_Ntrack",  &cc4_Ntrack,  "cc4_Ntrack/I");
  outTree_.Branch("cc4_sumIso",  &cc4_sumIso,  "cc4_sumIso/F");
  outTree_.Branch("cc4_hasJet",  &cc4_hasJet,  "cc4_hasJet/B");
  outTree_.Branch("cc4_Vieta",   &cc4_Vieta);
  outTree_.Branch("cc4_Viphi",   &cc4_Viphi);
  outTree_.Branch("cc4_Vet",     &cc4_Vet);
  outTree_.Branch("cc4_Ve",      &cc4_Ve);
  
  outTree_.Branch("ecalIso1_Vieta", &ecalIso1_Vieta);
  outTree_.Branch("ecalIso1_Viphi", &ecalIso1_Viphi);
  outTree_.Branch("ecalIso1_Vet",   &ecalIso1_Vet);
  outTree_.Branch("ecalIso1_Ve",    &ecalIso1_Ve);
  
}


void
DiphotonAnalysis::fillTree() {

  TreeFillPhoton();
  TreeFillJet();
  TreeFillVertices();
  
  Run     = _e.run();
  Event   = _e.event();
  LumiSec = _e.lumisec();
  
  // missing et and transverse mass
  met     = _e.met(EventManager::kPfMet)->pt();
  met_phi = _e.met(EventManager::kPfMet)->phi();
  TLorentzVector  vE, vN, vW;
  vE.SetPtEtaPhiE(et1,0.,phi1,   et1);
  vN.SetPtEtaPhiE(met,0.,met_phi,met);
  vW = vE + vN;
  mtrans = vW.M();
  ptHat  = _e.ptHat();
  //weight = _e.weight();
  
  nTracks   = _e.tracks().size();
  nClusters = 0;//_e.clusters().size();
  nPixSeeds = 0;//_e.pixelSeeds().size();
  nSC       = 0;//_e.superClusters().size();
  nRecHits  = (int) _e.a().n("rh");//_e.recHits().size();

  TDirectory* savedDir = gDirectory;
  _outputFile->cd();
  outTree_.Fill();
  savedDir->cd();

}

void 
DiphotonAnalysis::TreeFillJet()
{
  nJets=0;
  hardestJetPt = secondJetPt = -1000 ;
  
  const CandList& Jets = _e.jetList( EventManager::kPatJet);
  for(int unsigned ij=0;ij<Jets.size();ij++) {
    Candidate* jet = Jets[ij];
    
    if(jet->pt() > 15) {
      if (ij==0 ) hardestJetPt = jet->pt() ; // leading jet pt
      if (ij==1 )  secondJetPt = jet->pt() ; // trailing jet pt
      nJets++;
    }
  }
}

void
DiphotonAnalysis::TreeFillPhoton() {

  const CandList& photons_ = _e.photons();
  nPhotons = photons_.size();
  float pi = TMath::Pi();

  // Trigger bits
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
  
  //Photon1
  CandInfo* infopho1 = thePhoton1_->info();  
  p1                            = thePhoton1_->p();
  px1                           = thePhoton1_->px();
  py1                           = thePhoton1_->py();
  pz1                           = thePhoton1_->pz();
  et1                           = thePhoton1_->pt();
  energy1                       = thePhoton1_->E();
  eta1                          = thePhoton1_->eta();
  phi1                          = thePhoton1_->phi();
  isEBGap1                      = infopho1->getBool("isEBGap");
  isEEGap1                      = infopho1->getBool("isEEGap");
  isEBEEGap1                    = infopho1->getBool("isEBEEGap");
  isEB1                         = infopho1->getBool("isEB");
  isEE1                         = infopho1->getBool("isEE");
  sigmaEtaEta1                  = infopho1->getFloat("sigmaEtaEta");
  sigmaIetaIeta1                = infopho1->getFloat("sigmaIetaIeta");
  e1x51                         = infopho1->getFloat("e1x5");
  e2x51                         = infopho1->getFloat("e2x5");
  ecalRecHitSumEtConeDR031      = infopho1->getFloat("ecalRecHitSumEtConeDR03");
  ecalRecHitSumEtConeDR041      = infopho1->getFloat("ecalRecHitSumEtConeDR04");
  hadronicDepth1OverEm1         = infopho1->getFloat("hadronicDepth1OverEm");
  hadronicDepth2OverEm1         = infopho1->getFloat("hadronicDepth2OverEm");
  hadronicOverEm1               = infopho1->getFloat("hadronicOverEm");
  hcalDepth1TowerSumEtConeDR031 = infopho1->getFloat("hcalDepth1TowerSumEtConeDR03");
  hcalDepth1TowerSumEtConeDR041 = infopho1->getFloat("hcalDepth1TowerSumEtConeDR04");
  hcalDepth2TowerSumEtConeDR031 = infopho1->getFloat("hcalDepth2TowerSumEtConeDR03");
  hcalDepth2TowerSumEtConeDR041 = infopho1->getFloat("hcalDepth2TowerSumEtConeDR04");
  hcalTowerSumEtConeDR031       = infopho1->getFloat("hcalTowerSumEtConeDR03");
  hcalTowerSumEtConeDR041       = infopho1->getFloat("hcalTowerSumEtConeDR04");
  nTrkHollowConeDR031           = infopho1->getFloat("nTrkHollowConeDR03");
  nTrkHollowConeDR041           = infopho1->getFloat("nTrkHollowConeDR04");
  nTrkSolidConeDR031            = infopho1->getFloat("nTrkSolidConeDR03");
  nTrkSolidConeDR041            = infopho1->getFloat("nTrkSolidConeDR04");
  trkSumPtHollowConeDR031       = infopho1->getFloat("trkSumPtHollowConeDR03");
  trkSumPtHollowConeDR041       = infopho1->getFloat("trkSumPtHollowConeDR04");
  trkSumPtSolidConeDR031        = infopho1->getFloat("trkSumPtSolidConeDR03");
  trkSumPtSolidConeDR041        = infopho1->getFloat("trkSumPtSolidConeDR04");
  patSumIso1                    = ecalRecHitSumEtConeDR041 + hcalTowerSumEtConeDR041 + trkSumPtHollowConeDR041;
  if( debug ) cout << "Basics1 done." << endl;

  Candidate* pho1sc_ = _e.getSecond("photon-SuperCluster",thePhoton1_);
  CandInfo*  infosc  = pho1sc_->info();
  
  // compute the isolation arount thePhoton1_
  keepEta.push_back(thePhoton1_->eta());
  keepPhi.push_back(thePhoton1_->phi());
  
  scPhi1      = pho1sc_->phi();
  scEta1      = pho1sc_->eta();
  scRawE1     = infosc->getFloat("rawEnergy");
  nXtalsSC1   = infosc->getInt("nHits");
  scPhiWidth1 = infosc->getFloat("phiWidth");
  scEtaWidth1 = infosc->getFloat("etaWidth");
  scIxSeed1   = infosc->getInt("ixSeed");
  scIySeed1   = infosc->getInt("iySeed");
  scIzSeed1   = infosc->getInt("izSeed");
  r91         = infosc->getFloat("R9");
  r41         = infosc->getFloat("R4");
  if( debug ) cout << "SC1 done." << endl;

  // get the equivalent electron
  
  
  float r_ecal1 = 129;
  float z_ecal1 = 325;
    
  const CandList& electrons_ = _e.electrons();
  float maxiDRsel = 0.1, elTkPt, impEta, impPhi, deltRA, deltRC;
  ConeSelector coneSel1(*thePhoton1_, maxiDRsel);
  vector<Candidate*> electrons1;
  Candidate*         phoElec1;
  Candidate*         elTk;
  Candidate*         elGsf;
  const CandAssoc&     El_Tk_  = _e.candAssoc("electron-Track");
  const CandAssoc&     El_Gsf_ = _e.candAssoc("electron-GsfTrack");
  for( int e=0; e<(int)electrons_.size(); e++) {
    if( coneSel1.accept(*electrons_[e]) ) electrons1.push_back(electrons_[e]);
  }
  if     ( electrons1.size()==0 ) {
    fBrem1 = -1000.;
    elTkPt = -1000.;
  }
  else if( electrons1.size()==1 ) {
    phoElec1 = electrons1[0];
    CandInfo* phoElecInfo1 = phoElec1->info();
    fBrem1 = phoElecInfo1->getFloat("fBrem");
    elTk  = El_Tk_.getSecond( phoElec1 ) ;
    elGsf = El_Gsf_.getSecond( phoElec1 ) ; 
    if( elTk==NULL ) elTkPt = impEta = impPhi = -1000;
    else {
      CandInfo* elInfo = elTk->info();
      elTkPt = elTk->pt();
      impEta = elInfo->getFloat("impEta");
      impPhi = elInfo->getFloat("impPhi");
      float hiteta, hitphi;
      TVector3 v3_ = elTk->pos();
      float x0_ = v3_.X();
      float y0_ = v3_.Y();
      float z0_ = v3_.Z();
      int qtk_ = (int) elTk->charge();
      KineUtils::helixToEnvelop( hiteta, hitphi, qtk_*elTkPt, elTk->eta(), elTk->phi(), x0_, y0_, z0_, r_ecal1, z_ecal1, 3.8);
      deltRA =  KineUtils::dR(hiteta, phoElec1->eta(), hitphi, phoElec1->phi() );
      deltRC =  KineUtils::dR(impEta, phoElec1->eta(), impPhi, phoElec1->phi() );
    }
  }
  else {
    phoElec1 = electrons1[0];
    for( int e=0; e<(int)electrons1.size(); e++ ) {
      if( fabs(electrons1[e]->pt()-et1)<fabs(phoElec1->pt()-et1) ) phoElec1 = electrons1[e];
    }
    CandInfo* phoElecInfo1 = phoElec1->info();
    fBrem1 = phoElecInfo1->getFloat("fBrem");
    elTk = El_Tk_.getSecond( phoElec1 ) ;
    if( elTk==NULL ) elTkPt = impEta = impPhi = -1000; 
    else {
      CandInfo* elInfo = elTk->info();
      elTkPt = elTk->pt();
      impEta = elInfo->getFloat("impEta");
      impPhi = elInfo->getFloat("impPhi");
      float hiteta, hitphi;
      TVector3 v3_ = elTk->pos();
      float x0_ = v3_.X();
      float y0_ = v3_.Y();
      float z0_ = v3_.Z();
      int qtk_ = (int) elTk->charge();
      KineUtils::helixToEnvelop( hiteta, hitphi, qtk_*elTkPt, elTk->eta(), elTk->phi(), x0_, y0_, z0_, r_ecal1, z_ecal1, 3.8);
      deltRA =  KineUtils::dR(hiteta, scEta1, hitphi, scPhi1 );
      deltRC =  KineUtils::dR(impEta, scEta1, impPhi, scPhi1 );
    }
  }
  
  fill( "dR", "dR_A", deltRA  );
  fill( "dR", "dR_C", deltRC  );
  if( debug ) cout << "Electron match business done." << endl;
  
  isConv1     = infopho1->getBool("isConverted");
  hasPixSeed1 = infopho1->getBool("hasPixelSeed");
  inCrack1    = inCrack(scIxSeed1, scIySeed1, scIzSeed1);
  
  nXtlSC1 = -1000;
  preshE1 = -1000;
  nESclu1 = -1000;

  if( debug ) cout << "Basics done." << endl;
  
  const CandList& genPlist = _e.mcTruthCandidates();
  CandInfo* infogm1 = hasGenMatch1 ? theMCpho1_->info() : 0 ;
  genDR1     = hasGenMatch1 ? KineUtils::dR(theMCpho1_->eta(), thePhoton1_->eta(), theMCpho1_->phi(), thePhoton1_->phi() ) : -1000;
  genEnergy1 = hasGenMatch1 ? theMCpho1_->E()           : -1000;
  genEta1    = hasGenMatch1 ? theMCpho1_->eta()         : -1000;
  genPhi1    = hasGenMatch1 ? theMCpho1_->phi()         : -1000;
  genEt1     = hasGenMatch1 ? theMCpho1_->pt()          : -1000;
  int moIdx  = hasGenMatch1 ? infogm1->getInt("motherIndex") : 0; 
  genMother1 = (moIdx==-1)? -1000 : (hasGenMatch1 ? genPlist[moIdx]->pdgCode() : -1000);
  genIso1    = hasGenMatch1 ? getGenIso(theMCpho1_,false)  : -1000;
  genIso1    = hasGenMatch1 ? getGenIso(theMCpho1_,true)  : -1000;
  isFrag1    = hasGenMatch1 && (abs(genMother1)<23) && (genIso1<5) ;
  
  const CandList& listTracks = _e.tracks();
  for( size_t ii=0; ii<listTracks.size(); ii++ ) {
  
    float r_ecal1 = 129;
    float z_ecal1 = 325;
    
    const Candidate* trk = listTracks[ii];
    const CandInfo* tkin = trk->info();
    
    impEta = tkin->getFloat("impEta");
    impPhi = tkin->getFloat("impPhi");
    
    float pt_  = trk->pt();
    float hiteta, hitphi;
    
    TVector3 v3_ = trk->pos();
    float x0_ = v3_.X();
    float y0_ = v3_.Y();
    float z0_ = v3_.Z();
    
    int qtk_ = (int) trk->charge();
    KineUtils::helixToEnvelop( hiteta, hitphi, qtk_*pt_, trk->eta(), trk->phi(), x0_, y0_, z0_, r_ecal1, z_ecal1, 3.8);
    
//     cout << "extEta(A/C) " << hiteta << "/" << impEta << " extPhi(A/C) " << hitphi << "/" << impPhi << endl;
    float delR =  KineUtils::dR(hiteta, impEta, hitphi, impPhi );
    
    fill( "dR", "tkImp", (_ievt%2==0)?delR:-delR  );
  }
  
  
  
  if( debug ) cout << "Gen1 done." << endl;

  jetDR1     = hasJetMatch1 ? KineUtils::dR(theJet1_->eta(), thePhoton1_->eta(), theJet1_->phi(), thePhoton1_->phi() ) : -1000;
  jetEnergy1 = hasJetMatch1 ? theJet1_->E()   : -1000;
  jetEta1    = hasJetMatch1 ? theJet1_->eta() : -1000;
  jetPhi1    = hasJetMatch1 ? theJet1_->phi() : -1000;
  jetEt1     = hasJetMatch1 ? theJet1_->pt()  : -1000;

  if( debug ) cout << "Jet1 done." << endl;

  //Photon2
  CandInfo* infopho2;
  if (hasPhoton2) infopho2 = thePhoton2_->info();  
  p2                            = hasPhoton2 ? thePhoton2_->p()   : -1000;
  px2                           = hasPhoton2 ? thePhoton2_->px()  : -1000;
  py2                           = hasPhoton2 ? thePhoton2_->py()  : -1000;
  pz2                           = hasPhoton2 ? thePhoton2_->pz()  : -1000;
  et2                           = hasPhoton2 ? thePhoton2_->pt()  : -1000;
  energy2                       = hasPhoton2 ? thePhoton2_->E()   : -1000;
  eta2                          = hasPhoton2 ? thePhoton2_->eta() : -1000;
  phi2                          = hasPhoton2 ? thePhoton2_->phi() : -1000;
  /*isEBGap2                      = hasPhoton2 ? infopho2->getBool("isEBGap")   : false;
  isEEGap2                      = hasPhoton2 ? infopho2->getBool("isEEGap")   : false;
  isEBEEGap2                    = hasPhoton2 ? infopho2->getBool("isEBEEGap") : false;
  isEB2                         = hasPhoton2 ? infopho2->getBool("isEB")      : false;
  isEE2                         = hasPhoton2 ? infopho2->getBool("isEE")      : false;*/
  sigmaEtaEta2                  = hasPhoton2 ? infopho2->getFloat("sigmaEtaEta")   : -1000;
  sigmaIetaIeta2                = hasPhoton2 ? infopho2->getFloat("sigmaIetaIeta") : -1000;
  e1x52                         = hasPhoton2 ? infopho2->getFloat("e1x5") : -1000;
  e2x52                         = hasPhoton2 ? infopho2->getFloat("e2x5") : -1000;
  ecalRecHitSumEtConeDR032      = hasPhoton2 ? infopho2->getFloat("ecalRecHitSumEtConeDR03") : -1000;
  ecalRecHitSumEtConeDR042      = hasPhoton2 ? infopho2->getFloat("ecalRecHitSumEtConeDR04") : -1000;
  hadronicDepth1OverEm2         = hasPhoton2 ? infopho2->getFloat("hadronicDepth1OverEm") : -1000;
  hadronicDepth2OverEm2         = hasPhoton2 ? infopho2->getFloat("hadronicDepth2OverEm") : -1000;
  hadronicOverEm2               = hasPhoton2 ? infopho2->getFloat("hadronicOverEm") : -1000;
  hcalDepth1TowerSumEtConeDR032 = hasPhoton2 ? infopho2->getFloat("hcalDepth1TowerSumEtConeDR03") : -1000;
  hcalDepth1TowerSumEtConeDR042 = hasPhoton2 ? infopho2->getFloat("hcalDepth1TowerSumEtConeDR04") : -1000;
  hcalDepth2TowerSumEtConeDR032 = hasPhoton2 ? infopho2->getFloat("hcalDepth2TowerSumEtConeDR03") : -1000;
  hcalDepth2TowerSumEtConeDR042 = hasPhoton2 ? infopho2->getFloat("hcalDepth2TowerSumEtConeDR04") : -1000;
  hcalTowerSumEtConeDR032       = hasPhoton2 ? infopho2->getFloat("hcalTowerSumEtConeDR03") : -1000;
  hcalTowerSumEtConeDR042       = hasPhoton2 ? infopho2->getFloat("hcalTowerSumEtConeDR04") : -1000;
  nTrkHollowConeDR032           = hasPhoton2 ? infopho2->getFloat("nTrkHollowConeDR03") : -1000;
  nTrkHollowConeDR042           = hasPhoton2 ? infopho2->getFloat("nTrkHollowConeDR04") : -1000;
  nTrkSolidConeDR032            = hasPhoton2 ? infopho2->getFloat("nTrkSolidConeDR03") : -1000;
  nTrkSolidConeDR042            = hasPhoton2 ? infopho2->getFloat("nTrkSolidConeDR04") : -1000;
  trkSumPtHollowConeDR032       = hasPhoton2 ? infopho2->getFloat("trkSumPtHollowConeDR03") : -1000;
  trkSumPtHollowConeDR042       = hasPhoton2 ? infopho2->getFloat("trkSumPtHollowConeDR04") : -1000;
  trkSumPtSolidConeDR032        = hasPhoton2 ? infopho2->getFloat("trkSumPtSolidConeDR03") : -1000;
  trkSumPtSolidConeDR042        = hasPhoton2 ? infopho2->getFloat("trkSumPtSolidConeDR04") : -1000;
  patSumIso2                    = hasPhoton2 ? ecalRecHitSumEtConeDR042 + hcalTowerSumEtConeDR042 + trkSumPtHollowConeDR042 : -1000;
  if( debug ) cout << "Basic2 done." << endl;

  Candidate* pho2sc_;
  CandInfo*  infosc2;
  if (hasPhoton2) {
    pho2sc_  = _e.getSecond("photon-SuperCluster",thePhoton2_);
    infosc2  = pho2sc_->info();
  }
  
  // compute the isolation arount thePhoton2_
  if( hasPhoton2 ) {
    keepEta.push_back(thePhoton2_->eta());
    keepPhi.push_back(thePhoton2_->phi());
  }
  
  scPhi2      = hasPhoton2 ? pho2sc_->phi()                 : -1000;
  scEta2      = hasPhoton2 ? pho2sc_->eta()                 : -1000;
  scRawE2     = hasPhoton2 ? infosc2->getFloat("rawEnergy") : -1000;
  nXtalsSC2   = hasPhoton2 ? infosc2->getInt("nHits")       : -1000;
  scPhiWidth2 = hasPhoton2 ? infosc2->getFloat("phiWidth")  : -1000;
  scEtaWidth2 = hasPhoton2 ? infosc2->getFloat("etaWidth")  : -1000;
  scIxSeed2   = hasPhoton2 ? infosc2->getInt("ixSeed")      : -1000;
  scIySeed2   = hasPhoton2 ? infosc2->getInt("iySeed")      : -1000;
  scIzSeed2   = hasPhoton2 ? infosc2->getInt("izSeed")      : -1000;
  r92         = hasPhoton2 ? infosc2->getFloat("R9")        : -1000;
  r42         = hasPhoton2 ? infosc2->getFloat("R4")        : -1000;

  if( debug ) cout << "Basics2 done." << endl;
  
  // get the equivalent electron
  if( hasPhoton2 ) {
    ConeSelector coneSel2(*thePhoton2_, maxiDRsel);
    vector<Candidate*> electrons2;
    Candidate*         phoElec2;
    for( int e=0; e<(int)electrons_.size(); e++) {
      if( coneSel2.accept(*electrons_[e]) ) electrons2.push_back(electrons_[e]);
    }
    if     ( electrons2.size()==0 ) fBrem2 = -1000.;
    else if( electrons2.size()==1 ) {
      phoElec2 = electrons2[0];
      CandInfo* phoElecInfo2 = phoElec2->info();
      fBrem2 = phoElecInfo2->getFloat("fBrem");
    }
    else {
      phoElec2 = electrons2[0];
      for( int e=0; e<(int)electrons2.size(); e++ ) {
	if( fabs(electrons2[e]->pt()-et2)<fabs(phoElec2->pt()-et2) ) phoElec2 = electrons2[e];
      }
      CandInfo* phoElecInfo2 = phoElec2->info();
      fBrem2 = phoElecInfo2->getFloat("fBrem");
    }
    if( debug && fBrem2>-100 ) cout << "phoEt " << et2 << " eleEt " << phoElec2->Et() << " fBrem " << fBrem2 << endl;
  }
  
  if( debug ) cout << "Electron2 done." << endl;
  
  isConv2     = hasPhoton2 ? infopho2->getBool("isConverted")         : false;
  hasPixSeed2 = hasPhoton2 ? infopho2->getBool("hasPixelSeed")        : false;
  inCrack2    = hasPhoton2 ? inCrack(scIxSeed2, scIySeed2, scIzSeed2) : false;
  
  CandInfo* infogm2 = hasGenMatch2 ? theMCpho2_->info() : 0 ;
  genDR2     = hasGenMatch2 ? KineUtils::dR(theMCpho2_->eta(), thePhoton2_->eta(), theMCpho2_->phi(), thePhoton2_->phi() ) : -1000;
  genEnergy2 = hasGenMatch2 ? theMCpho2_->E()           : -1000;
  genEta2    = hasGenMatch2 ? theMCpho2_->eta()         : -1000;
  genPhi2    = hasGenMatch2 ? theMCpho2_->phi()         : -1000;
  genEt2     = hasGenMatch2 ? theMCpho2_->pt()          : -1000;
  moIdx      = hasGenMatch2 ? infogm2->getInt("motherIndex") : 0;
  genMother2 = (moIdx==-1)? -1000 : (hasGenMatch2 ? genPlist[moIdx]->pdgCode() : -1000);
  genIso2    = hasGenMatch2 ? getGenIso(theMCpho2_,false)  : -1000;
  genIsoAll2 = hasGenMatch2 ? getGenIso(theMCpho2_,true)   : -1000;
  isFrag2    = hasGenMatch2 && (abs(genMother2)<23) && (genIso2<5) ;  

  if( debug ) cout << "Gen2 done." << endl;
  
  jetDR2     = hasJetMatch2 ? KineUtils::dR(theJet2_->eta(), thePhoton2_->eta(), theJet2_->phi(), thePhoton2_->phi() ) : -1000;
  jetEnergy2 = hasJetMatch2 ? theJet2_->E() : -1000;
  jetEta2    = hasJetMatch2 ? theJet2_->eta() : -1000;
  jetPhi2    = hasJetMatch2 ? theJet2_->phi() : -1000;
  jetEt2     = hasJetMatch2 ? theJet2_->pt() : -1000;
  
  nXtlSC2 = -1000;
  preshE2 = -1000;
  nESclu2 = -1000;

  if( debug ) cout << "Jet2 done." << endl;

  // Diphoton variables
  
  Candidate* diphoton;
  if (hasPhoton2) diphoton = Candidate::create(thePhoton1_,thePhoton2_);
  diphoMass = hasPhoton2 ? diphoton->mass() : -1000;
  diphoQt   = hasPhoton2 ? diphoton->pt() : -1000;
  diphoEta  = hasPhoton2 ? diphoton->eta()  : -1000;
  
  diphoPhi  = hasPhoton2 ? diphoton->phi()  : -1000;
  deltaEta  = hasPhoton2 ? fabs(eta1-eta2)  : -1000;
  float protoDphi = fabs(phi1-phi2) ; 
  if( protoDphi>TMath::Pi() ) protoDphi = 2*TMath::Pi()-protoDphi;
  deltaPhi  = hasPhoton2 ? protoDphi  : -1000;
  
  deltaY    = hasPhoton2 ? 0.25*log((energy1+pz1)/(energy2+pz2)*(energy2-pz2)/(energy1-pz1)) : -1000; 

  diphoY    = hasPhoton2 ? 0.5*log( (diphoton->E()+diphoton->pz())/(diphoton->E()-diphoton->pz()) ) : -1000 ; 

  protoDphi = min(fabs(deltaPhi),2*pi-fabs(deltaPhi));
  deltaR    = hasPhoton2 ? sqrt( deltaEta*deltaEta + protoDphi*protoDphi ) : -1000;
  cosTheta  = hasPhoton2 ? ROOT::Math::VectorUtil::CosTheta(thePhoton1_->p3(), diphoton->p3()) : -1000;
  
  TLorentzVector p1, p2;
  TLorentzVector pL, pT;
  if( hasPhoton2 ) {
    if( !isShuffled ) {
      p1 = thePhoton1_->p4();
      p2 = thePhoton2_->p4();
      pL = thePhoton1_->p4();
      pT = thePhoton2_->p4();
    }
    else {
      p1 = thePhoton2_->p4();
      p2 = thePhoton1_->p4();
      pL = thePhoton2_->p4();
      pT = thePhoton1_->p4();
    }
    ROOT::Math::LorentzVector< ROOT::Math::PxPyPzE4D<double> > lv(diphoton->px(),diphoton->py(),diphoton->pz(),diphoton->E());
    ROOT::Math::DisplacementVector3D< ROOT::Math::Cartesian3D<double> > bvt = lv.BoostToCM();
    TVector3 bv(bvt.x(), bvt.y(), bvt.z());
    //cout << "Boost " << bv.X() << " " << bv.Y() << " " << bv.Z() << endl;
    pL.Boost(bv);
    pT.Boost(bv);
  }

  cosThetaH = hasPhoton2 ? cosThetaStar(p1,p2) : -1000;
  
  cosThetaS = hasPhoton2 ? ROOT::Math::VectorUtil::CosTheta(pL.Vect(), diphoton->p3()) : -1000;  

  cosThetaP = hasPhoton2 ? TMath::Abs( tanh( 0.25*log( (energy1+pz1)/(energy2+pz2)*(energy2-pz2)/(energy1-pz1) ) ) ) : -1000;

  if( debug ) cout << "Diphoton business done." << endl;

  // Random Cone variables
  
  cc1_eta      = thePhoton1_->eta() ;
  float ccTmp  = (2*random.Rndm()-1)*pi/2 ; 
  int ccTmpNeg = (ccTmp<0);
  cc1_phi      = thePhoton1_->phi() + TMath::Power(-1,ccTmpNeg)*pi/4 + ccTmp; // random in [45,135] away from photon
  
  cc2_eta      = cc1_eta;
  cc2_phi      = (cc1_phi+pi<2*pi)? cc1_phi+pi : cc1_phi-pi;

  if( hasPhoton2 ) {
    cc3_eta      = thePhoton2_->eta();
    ccTmp        = (2*random.Rndm()-1)*pi/2 ; 
    ccTmpNeg     = (ccTmp<0);
    cc3_phi      = thePhoton2_->phi() + TMath::Power(-1,ccTmpNeg)*pi/4 + ccTmp; // random in [45,135] away from photon
    
    cc4_eta      = cc3_eta;
    cc4_phi      = (cc3_phi+pi<2*pi)? cc3_phi+pi : cc3_phi-pi;
  }
  
  // create the 2/4 cc candidates
  Candidate* cc1 = Candidate::create();  cc1->setPtEtaPhi(10.,cc1_eta,cc1_phi);
  Vertex* ccVer = _e.vertices()[0];
  cc1->setVertex( ccVer ); cc1->lock();
  
  Candidate* cc2 = Candidate::create();  cc2->setPtEtaPhi(10.,cc2_eta,cc2_phi);
  cc2->setVertex( ccVer ); cc2->lock();
  
  Candidate* cc3 = Candidate::create();
  Candidate* cc4 = Candidate::create();
  
  if( hasPhoton2 ) {
    cc3->setPtEtaPhi(10.,cc3_eta,cc3_phi);
    cc3->setVertex( ccVer ); cc3->lock();
    cc4->setPtEtaPhi(10.,cc4_eta,cc4_phi);
    cc4->setVertex( ccVer ); cc4->lock();
  }
  
  // compute the isolation around cc1to4
  keepEta.push_back(cc1->eta());
  keepPhi.push_back(cc1->phi());
  keepEta.push_back(cc2->eta());
  keepPhi.push_back(cc2->phi());
  if( hasPhoton2 ) {
    keepEta.push_back(cc3->eta());
    keepPhi.push_back(cc3->phi());
    keepEta.push_back(cc4->eta());
    keepPhi.push_back(cc4->phi());
  }

  if( debug ) cout << "Randcom cones created." << endl;

  // Compute all the isolations, quicker with the keepEta, keepPhi method
  makeRHvector();
  
  RandomConeSelector isoSel1(rhVector);
  isoSel1.computeIso(*thePhoton1_,*pho1sc_,false); // boolean to print

  ecalIsoCMS1   = infopho1->getFloat("ecalIso");
  char sochn[200];
  sprintf(sochn,"eta %7.3f phi %7.3f  ecalIso   hand %7.3f   cms %7.3f  diff %7.3f",scEta1,scPhi1,ecalIso1,ecalIsoCMS1,ecalIso1-ecalIsoCMS1);
  //cout << sochn << endl;
  ecalIso1_n    = isoSel1.getNecal();
  ecalIso1_ntk  = isoSel1.getNetk();

  ecalIso1_Vieta = isoSel1.getVieta();
  ecalIso1_Viphi = isoSel1.getViphi();
  ecalIso1_Veta  = isoSel1.getVeta();
  ecalIso1_Vphi  = isoSel1.getVphi();
  ecalIso1_Vet   = isoSel1.getVet();
  ecalIso1_Ve    = isoSel1.getVe();

  ecalIso1_Vtkpt  = isoSel1.getVtkpt();
  ecalIso1_VtkCnt = isoSel1.getVtkCnt();
  ecalIso1_Vtkd0  = isoSel1.getVtkd0();
  ecalIso1_Vtkz0  = isoSel1.getVtkz0();
  ecalIso1_VtkNh  = isoSel1.getVtkNh();
  ecalIso1_VtkIh  = isoSel1.getVtkIh(); // expInnHits
  ecalIso1_VtkChi2  = isoSel1.getVtkChi2();
  ecalIso1_VtkNdof  = isoSel1.getVtkNdof();
  ecalIso1_VtkQ     = isoSel1.getVtkQ();
  ecalIso1_VtkEtaH  = isoSel1.getVtkEtaH();
  ecalIso1_VtkPhiH  = isoSel1.getVtkPhiH();
  ecalIso1_VtkDEta  = isoSel1.getVtkDEta();
  ecalIso1_VtkDPhi  = isoSel1.getVtkDPhi();
  ecalIso1_VtkDR    = isoSel1.getVtkDR();
  ecalIso1_VtkFP    = isoSel1.getVtkFP();
  ecalIso1_VtkSV    = isoSel1.getVtkSV();
  ecalIso1_VtkInV   = isoSel1.getVtkInVeto();
  ecalIso1_VtkVid   = isoSel1.getVtkVid();
  ecalIso1_VtkDepE  = isoSel1.getVtkDepE();
  ecalIso1_VtkDepEt = isoSel1.getVtkDepEt();
  if( debug ) cout << "Simple iso done." << endl;
  
  // how many good tracks in ecalIso zone?
  float niceTkPt = 3.;
  ecalIso1_nTkGinV = 0;
  ecalIso1_nTkGnoZ = 0;
  ecalIso1_nTkG    = 0;
  for( int t=0; t<(int)ecalIso1_Vtkpt.size(); t++) {
    if( ecalIso1_VtkSV[t]==false || ecalIso1_VtkFP[t]==false || ecalIso1_Vtkpt[t]<niceTkPt ) continue;
    ecalIso1_nTkG++;
    if( ecalIso1_VtkInV[t]==true ) ecalIso1_nTkGinV++;
    if( fabs(elTkPt-ecalIso1_Vtkpt[t])/elTkPt>0.01 && elTkPt/et1>0.6) { ecalIso1_nTkGnoZ++; }
  }
  // other cuts
  ecalIso1_nTkG_28 = 0;
  ecalIso1_nTkG_32 = 0;
  ecalIso1_nTkG_40 = 0;
  for( int t=0; t<(int)ecalIso1_Vtkpt.size(); t++) {
    if( ecalIso1_VtkSV[t]==false || ecalIso1_VtkFP[t]==false ) continue;
    if( ecalIso1_Vtkpt[t]>2.8 ) ecalIso1_nTkG_28++;
    if( ecalIso1_Vtkpt[t]>3.2 ) ecalIso1_nTkG_32++;
    if( ecalIso1_Vtkpt[t]>4.0 ) ecalIso1_nTkG_40++;
  }

  ecalIso1     = 0;
  ecalIso1_280 = 0;
  ecalIso1_320 = 0;
  ecalIso1_350 = 0;
  ecalIso1_400 = 0;
  for( int x=0; x<(int)ecalIso1_Vet.size(); x++ ) {
    if( ecalIso1_Vet[x]>0.300 ) ecalIso1     += ecalIso1_Vet[x];
    if( ecalIso1_Vet[x]>0.280 ) ecalIso1_280 += ecalIso1_Vet[x];
    if( ecalIso1_Vet[x]>0.320 ) ecalIso1_320 += ecalIso1_Vet[x];
    if( ecalIso1_Vet[x]>0.350 ) ecalIso1_350 += ecalIso1_Vet[x];
    if( ecalIso1_Vet[x]>0.400 ) ecalIso1_400 += ecalIso1_Vet[x];
  }

  // more isolation
  hcalIso1    = isoSel1.getHcalIso();
  hcalIsoH1   = isoSel1.getHcalIsoH();
  trackIso1   = isoSel1.getTrackIso();
  float tmpTrackIso = 0.;
  for( int t=0; t<(int)ecalIso1_Vtkpt.size(); t++ ) tmpTrackIso += (ecalIso1_VtkCnt[t]) ? ecalIso1_Vtkpt[t] : 0.;
  trackIsoRd1 = trkSumPtHollowConeDR041 - tmpTrackIso; // reduced trackIso
  sumIso1     = isoSel1.getSumIso();

  hcalIso1_n   = isoSel1.getNhcal();
  hcalIsoH1_n  = isoSel1.getNhcalH();
  trackIso1_n  = isoSel1.getNtrack();
  if( debug ) cout << "Track and hcalIso done." << endl;
  
  isoSel1.setIC(0.);
  isoSel1.setStrip(0.);
  isoSel1.computeEcalIso(*thePhoton1_,*pho1sc_,false);
  
  ecalIso1_Vieta_nv = isoSel1.getVieta();
  ecalIso1_Viphi_nv = isoSel1.getViphi();
  ecalIso1_Veta_nv  = isoSel1.getVeta();
  ecalIso1_Vphi_nv  = isoSel1.getVphi();
  ecalIso1_Vet_nv   = isoSel1.getVet();
  ecalIso1_Ve_nv    = isoSel1.getVe();
  if( debug ) cout << "noVeto done." << endl;
 
  // when I will finally change the cut values
  eIso1_RECO = -1000.;
  eIso2_RECO = -1000.;
  
  RandomConeSelector isoSel2(rhVector); 
  if( hasPhoton2 ) isoSel2.computeIso(*thePhoton2_,*pho2sc_);
  //ecalIso2      = hasPhoton2 ? isoSel2.getEcalIso()  : -1000;
  ecalIso2_n    = hasPhoton2 ? isoSel2.getNecal()    : -1000;
  ecalIso2_ntk  = hasPhoton2 ? isoSel2.getNetk()     : -1000;

  vector<int> tmpi; vector<float> tmpf; vector<bool> tmpb;
  ecalIso2_Vieta = hasPhoton2 ? isoSel2.getVieta()   : tmpi;
  ecalIso2_Viphi = hasPhoton2 ? isoSel2.getViphi()   : tmpi;
  ecalIso2_Veta  = hasPhoton2 ? isoSel2.getVeta()    : tmpf;
  ecalIso2_Vphi  = hasPhoton2 ? isoSel2.getVphi()    : tmpf;
  ecalIso2_Vet   = hasPhoton2 ? isoSel2.getVet()     : tmpf;
  ecalIso2_Ve    = hasPhoton2 ? isoSel2.getVe()      : tmpf;

  ecalIso2_Vtkpt  = hasPhoton2 ? isoSel2.getVtkpt()        : tmpf;
  ecalIso2_VtkCnt = hasPhoton2 ? isoSel2.getVtkCnt()       : tmpb;
  ecalIso2_Vtkd0  = hasPhoton2 ? isoSel2.getVtkd0()        : tmpf;
  ecalIso2_Vtkz0  = hasPhoton2 ? isoSel2.getVtkz0()        : tmpf;
  ecalIso2_VtkNh  = hasPhoton2 ? isoSel2.getVtkNh()        : tmpi;
  ecalIso2_VtkIh  = hasPhoton2 ? isoSel2.getVtkIh()        : tmpi; // expInnHits
  ecalIso2_VtkChi2  = hasPhoton2 ? isoSel2.getVtkChi2()   : tmpf;
  ecalIso2_VtkNdof  = hasPhoton2 ? isoSel2.getVtkNdof()   : tmpi;
  ecalIso2_VtkQ     = hasPhoton2 ? isoSel2.getVtkQ()      : tmpi;
  ecalIso2_VtkEtaH  = hasPhoton2 ? isoSel2.getVtkEtaH()   : tmpf;
  ecalIso2_VtkPhiH  = hasPhoton2 ? isoSel2.getVtkPhiH()   : tmpf;
  ecalIso2_VtkDEta  = hasPhoton2 ? isoSel2.getVtkDEta()   : tmpf;
  ecalIso2_VtkDPhi  = hasPhoton2 ? isoSel2.getVtkDPhi()   : tmpf;
  ecalIso2_VtkDR    = hasPhoton2 ? isoSel2.getVtkDR()     : tmpf;
  ecalIso2_VtkFP    = hasPhoton2 ? isoSel2.getVtkFP()     : tmpb;
  ecalIso2_VtkSV    = hasPhoton2 ? isoSel2.getVtkSV()     : tmpb;
  ecalIso2_VtkInV   = hasPhoton2 ? isoSel2.getVtkInVeto() : tmpb;
  ecalIso2_VtkVid   = hasPhoton2 ? isoSel2.getVtkVid()    : tmpi;
  ecalIso2_VtkDepE  = hasPhoton2 ? isoSel2.getVtkDepE()   : tmpf;
  ecalIso2_VtkDepEt = hasPhoton2 ? isoSel2.getVtkDepEt()  : tmpf;
  if( debug ) cout << "Simple iso done." << endl;
  
  // how many good tracks in ecalIso zone?
  ecalIso2_nTkGinV = 0;
  ecalIso2_nTkG    = 0;
  if( hasPhoton2) {
    for( int t=0; t<(int)ecalIso2_Vtkpt.size(); t++) {
      if( ecalIso2_VtkSV[t]==false || ecalIso2_VtkFP[t]==false || ecalIso2_Vtkpt[t]<niceTkPt ) continue;
      ecalIso2_nTkG++;
      if( ecalIso2_VtkInV[t]==true ) ecalIso1_nTkGinV++;
    }
  }
  // other cuts
  ecalIso2_nTkG_28 = hasPhoton2 ? 0 : -1000;
  ecalIso2_nTkG_32 = hasPhoton2 ? 0 : -1000;
  ecalIso2_nTkG_40 = hasPhoton2 ? 0 : -1000;
  if( hasPhoton2 ) for( int t=0; t<(int)ecalIso2_Vtkpt.size(); t++) {
    if( ecalIso2_VtkSV[t]==false || ecalIso2_VtkFP[t]==false ) continue;
    if( ecalIso2_Vtkpt[t]>2.8 ) ecalIso2_nTkG_28++;
    if( ecalIso2_Vtkpt[t]>3.2 ) ecalIso2_nTkG_32++;
    if( ecalIso2_Vtkpt[t]>4.0 ) ecalIso2_nTkG_40++;
  }

  ecalIso2     = hasPhoton2 ? 0 : -1000;
  ecalIso2_280 = hasPhoton2 ? 0 : -1000;
  ecalIso2_320 = hasPhoton2 ? 0 : -1000;
  ecalIso2_350 = hasPhoton2 ? 0 : -1000;
  ecalIso2_400 = hasPhoton2 ? 0 : -1000;
  if( hasPhoton2 ) for( int x=0; x<(int)ecalIso2_Vet.size(); x++ ) {
    if( ecalIso2_Vet[x]>0.300 ) ecalIso2     += ecalIso2_Vet[x];
    if( ecalIso2_Vet[x]>0.280 ) ecalIso2_280 += ecalIso2_Vet[x];
    if( ecalIso2_Vet[x]>0.320 ) ecalIso2_320 += ecalIso2_Vet[x];
    if( ecalIso2_Vet[x]>0.350 ) ecalIso2_350 += ecalIso2_Vet[x];
    if( ecalIso2_Vet[x]>0.400 ) ecalIso2_400 += ecalIso2_Vet[x];
  }

  // other cuts
  hcalIso2    = hasPhoton2 ? isoSel2.getHcalIso()   : -1000;
  hcalIsoH2   = hasPhoton2 ? isoSel2.getHcalIsoH()  : -1000;  
  trackIso2   = hasPhoton2 ? isoSel2.getTrackIso()  : -1000;
  tmpTrackIso = 0.;
  for( int t=0; t<(int)ecalIso2_Vtkpt.size(); t++ ) tmpTrackIso += ecalIso2_VtkCnt[t]*ecalIso2_Vtkpt[t];
  trackIsoRd2 = trkSumPtHollowConeDR042 - tmpTrackIso; // reduced trackIso
  sumIso2     = hasPhoton2 ? isoSel2.getSumIso()    : -1000;
  //vector<int>   cc1_ecalIso1_Vieta_nv, cc1_ecalIso1_Viphi_nv;
  //vector<float> cc1_ecalIso1_Veta_nv,  cc1_ecalIso1_Vphi_nv;
  //vector<float> cc1_ecalIso1_Vet_nv,   cc1_ecalIso1_Ve_nv;

  if( debug ) cout << "RandomConeSelector business1 done." << endl;

  RandomConeSelector ccSel(rhVector);
  ccSel.computeIso(*cc1,*cc1);
  cc1_hcalIso  = ccSel.getHcalIso();
  cc1_trackIso = ccSel.getTrackIso();
  cc1_sumIso   = ccSel.getSumIso();
  cc1_hasJet   = jetInCone(cc1);
  cc1_Nhcal     = ccSel.getNhcal();
  cc1_Ntrack    = ccSel.getNtrack();

  cc1_ecalIso       = ccSel.getEcalIso();
  cc1_ecalIso_n     = ccSel.getNecal();
  cc1_ecalIso_ntk   = ccSel.getNetk();
  cc1_ecalIso_Vieta = ccSel.getVieta();
  cc1_ecalIso_Viphi = ccSel.getViphi();
  cc1_ecalIso_Veta  = ccSel.getVeta();
  cc1_ecalIso_Vphi  = ccSel.getVphi();
  cc1_ecalIso_Vet   = ccSel.getVet();
  cc1_ecalIso_Ve    = ccSel.getVe();
   
  cc1_ecalIso_Vtkpt     = ccSel.getVtkpt();
  cc1_ecalIso_Vtkd0     = ccSel.getVtkd0();
  cc1_ecalIso_Vtkz0     = ccSel.getVtkz0();
  cc1_ecalIso_VtkNh     = ccSel.getVtkNh();
  cc1_ecalIso_VtkIh     = ccSel.getVtkIh();
  cc1_ecalIso_VtkChi2   = ccSel.getVtkChi2();
  cc1_ecalIso_VtkNdof   = ccSel.getVtkNdof();
  cc1_ecalIso_VtkQ      = ccSel.getVtkQ();
  cc1_ecalIso_VtkFP     = ccSel.getVtkFP();
  cc1_ecalIso_VtkSV     = ccSel.getVtkSV();
  cc1_ecalIso_VtkEtaH   = ccSel.getVtkEtaH();
  cc1_ecalIso_VtkPhiH   = ccSel.getVtkPhiH();
  cc1_ecalIso_VtkInV    = ccSel.getVtkInVeto();
  cc1_ecalIso_VtkVid    = ccSel.getVtkVid();
  cc1_ecalIso_VtkDepE   = ccSel.getVtkDepE();
  cc1_ecalIso_VtkDepEt  = ccSel.getVtkDepEt();

  cc1_ecalIso_nTkG   = 0;
  for( int t=0; t<(int)cc1_ecalIso_Vtkpt.size(); t++) {
    if( cc1_ecalIso_VtkSV[t]==false || cc1_ecalIso_VtkFP[t]==false || cc1_ecalIso_Vtkpt[t]<niceTkPt ) continue;
    cc1_ecalIso_nTkG++;
  }
  
  ccSel.computeIso(*cc2,*cc2);
  cc2_ecalIso  = ccSel.getEcalIso();
  cc2_hcalIso  = ccSel.getHcalIso();
  cc2_trackIso = ccSel.getTrackIso();
  cc2_sumIso   = ccSel.getSumIso();
  cc2_Necal    = ccSel.getNecal();
  cc2_Nhcal    = ccSel.getNhcal();
  cc2_Ntrack   = ccSel.getNtrack();
  cc2_hasJet   = jetInCone(cc2);
  
  if( hasPhoton2 ) ccSel.computeIso(*cc3,*cc3);
  cc3_ecalIso  = (!hasPhoton2)? -1000 : ccSel.getEcalIso();
  cc3_hcalIso  = (!hasPhoton2)? -1000 : ccSel.getHcalIso();
  cc3_trackIso = (!hasPhoton2)? -1000 : ccSel.getTrackIso();
  cc3_sumIso   = (!hasPhoton2)? -1000 : ccSel.getSumIso();
  cc3_Necal    = (!hasPhoton2)? -1000 : ccSel.getNecal();
  cc3_Nhcal    = (!hasPhoton2)? -1000 : ccSel.getNhcal();
  cc3_Ntrack   = (!hasPhoton2)? -1000 : ccSel.getNtrack();
  cc3_hasJet   = (!hasPhoton2)? false : jetInCone(cc3);
  
  if( hasPhoton2 ) ccSel.computeIso(*cc4,*cc4);
  cc4_ecalIso  = (!hasPhoton2)? -1000 : ccSel.getEcalIso();
  cc4_hcalIso  = (!hasPhoton2)? -1000 : ccSel.getHcalIso();
  cc4_trackIso = (!hasPhoton2)? -1000 : ccSel.getTrackIso();
  cc4_sumIso   = (!hasPhoton2)? -1000 : ccSel.getSumIso();
  cc4_Necal    = (!hasPhoton2)? -1000 : ccSel.getNecal();
  cc4_Nhcal    = (!hasPhoton2)? -1000 : ccSel.getNhcal();
  cc4_Ntrack   = (!hasPhoton2)? -1000 : ccSel.getNtrack();
  cc4_hasJet   = (!hasPhoton2)? false : jetInCone(cc4);
  if( debug ) cout << "Rest RC done." << endl;
  
}

void 
DiphotonAnalysis::TreeFillVertices() {

  const VtxList& vertices = _e.vertices();
 
  nVertex = (int)vertices.size();
  nGoodVertex = 0;
  for( int v=0; v<(int)vertices.size(); v++) {
    Vertex* ver = vertices[v];
    CandInfo* infop = ver->info();
    bool vergo = infop->getBool("isGoodVtx");
    if( vergo ) nGoodVertex++;
  }
  Vertex* Primary = vertices[0];
  vX = Primary->pos().X();
  vY = Primary->pos().Y();
  vZ = Primary->pos().Z();
  
  CandInfo* infov = Primary->info();
  vNtracks = infov->getInt("nTracks");
  vNout    = Primary->nOutgoingCand();
  vIsGood  = infov->getBool("isGoodVtx");
}

bool
DiphotonAnalysis::SpikeCleaning(Candidate* cand) {

  bool spike=false;
  
  if( SampleName.substr(0,2)!="EG" && SampleName.substr(0,6)!="Photon" ) return false;
 
  //Load RecHits
  const CandMap&         Sc_Bc_ = _e.candMap("SuperCluster-BasicCluster");
  const CandIdRecHitMap& Bc_rh_ = _e.ecalRecHitMap("BasicCluster-EcalRecHit");

  //Leading photon
  Candidate* superCluster     = _e.getSecond("photon-SuperCluster", cand );
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
	//GHM        float swisscross = rh.SpikeVar.swCross;
	float swisscross = 1;
        if( swisscross > 0.95 ) { spike=true; }
	int severitylevel = 0;
	if( severitylevel > 0 ) { spike=true; }
      }
    }
  }
  return spike;
}

void
DiphotonAnalysis::MCtruthMatching() {
  const CandList& genParticles_ = _e.mcTruthCandidates();
   
  MCMatchingSelector matcher;
  matcher.loadMCMatchPhoton();
  hasGenMatch1 = hasGenMatch2 = false;
  
  // mc truth matching for the first photon
  hasGenMatch1 = matcher.MCMatchTruth(thePhoton1_, genParticles_);
  if (hasGenMatch1) theMCpho1_ = matcher.MCMatch(thePhoton1_, genParticles_);

  // mc truth matching for the second photon
  if (hasPhoton2) {
    hasGenMatch2 = matcher.MCMatchTruth(thePhoton2_, genParticles_);
    if (hasGenMatch2) theMCpho2_ = matcher.MCMatch(thePhoton2_, genParticles_);
  }
}

void
DiphotonAnalysis::JetMatching() {
  
  float minJetPt_ = 20.;
  float maxJetDR_ =  2.;

  const CandList& Jets = _e.jetList( EventManager::kPfJet);
  hasJetMatch1 = hasJetMatch2 = false;
  
  float dRtmp = 20.;
  float dRmin = 10.;
  
  // jet matching for the first photon
  for(int unsigned ij=0;ij<Jets.size();ij++) {
    Candidate* jet = Jets[ij];
    dRtmp = KineUtils::dR(jet->eta(), thePhoton1_->eta(), jet->phi(), thePhoton1_->phi() );
    if (jet->pt()>minJetPt_ && dRtmp>maxJetDR_ && dRtmp<dRmin) {
      dRmin = dRtmp;
      hasJetMatch1 = true;
      theJet1_ = jet;
    }
  }
   
  dRtmp = 20.;
  dRmin = 10.;
   
  // mc truth matching for the second photon
  if (hasPhoton2) {
    for (int unsigned ij=0;ij<Jets.size();ij++) {
      Candidate* jet = Jets[ij];
      dRtmp = KineUtils::dR(jet->eta(), thePhoton2_->eta(), jet->phi(), thePhoton2_->phi() );
      if (jet->pt()>minJetPt_ && dRtmp>maxJetDR_ && dRtmp<dRmin) {
        dRmin = dRtmp;
        hasJetMatch2 = true;
        theJet2_ = jet;
      }
    }
  }
}

float
DiphotonAnalysis::getGenIso(Candidate* candMatch, bool doAll)
{
  const float etMin  = 0.0;
  const float dRout  = 0.4;
  float genIsoSum = 0.0;
  const CandList& mcCandList = _e.mcTruthCandidates();
  CandInfo*  infomatch = candMatch->info();
  
  for(unsigned int imc=0;imc<mcCandList.size();imc++) {
    
    Candidate* mcCand = mcCandList[imc];
    CandInfo*  infomc = mcCand->info();
    
    // do not count the matched particle itself in the isolation
    if (mcCand->uid() == candMatch->uid()) continue; 
    // collisionID should be the same as matched photon
    if ( infomc->getInt("collId") != infomatch->getInt("collId") ) continue;
    // do not count unstable particles
    if ( infomc->getInt("status")!=1 ) continue;
    // do not keep muons and neutrinos 
    float pdgIdrun = fabs(mcCand->pdgCode());
    if (!doAll && pdgIdrun>11 && pdgIdrun<20) continue;
    
    float etrun = mcCand->pt();
    if (etrun<etMin) continue;
    
    float dRrun = KineUtils::dR(mcCand->eta(), candMatch->eta(), mcCand->phi(), candMatch->phi() );
    if (dRrun>dRout) continue;

    //cout << "pt " <<  mcCand->pt() << " uid " << mcCand->uid() << " pdg " << mcCand->pdgCode() << " status " << infomc->getInt("status") << " ptCand " << candMatch->pt() << " uidCand " << candMatch->uid() << " pidCand " << candMatch->pdgCode() << " statusCand " << infomatch->getInt("status") << endl;
        
    genIsoSum += etrun;
  }
  
  return genIsoSum;
}


bool
DiphotonAnalysis::inCrack(int ieta, int iphi, int iz)
{
  int c_eta[8] = {1,25,    26,45,   46,65,   66,85}; 
  int c_phi[36]= {1,20,    21,40,   41,60,   61,80,  
                  81,100,  101,120, 121,140, 141,160,
                  161,180, 181,200, 201,220, 221,240,
                  241,260, 261,280, 281,300, 301,320, 
                  321,340, 341,360};
  
  bool crackresult = false;
  if (iz != 0) return false;
  
  for (int eta=0; eta<8; eta++)  if (abs(ieta)==c_eta[eta]) crackresult = true;
  for (int phi=0; phi<36; phi++) if (iphi     ==c_phi[phi]) crackresult = true;
  
  return crackresult;
}

void
DiphotonAnalysis::makeRHvector()
    // keep the indices of the interesting crystals
    // cut on |xtalE| and on DR with all interesting points
{
  rhVector.clear();
  bool keep = false;
  float eRH, etaRH, phiRH, etaK, phiK, dRtmp;
  int ietaRH, iphiRH, izRH;
  if (keepEta.size() != keepPhi.size()) {cout << "[makeRHvector] Size problem." << endl; return;}
  size_t nRecHits = (size_t) _e.a().n("rh");
  for(size_t irh=0;irh<nRecHits;irh++) {
    keep = false;
    assert( _e.a().load("rh", irh) );
    eRH = _e.a().get_f("rh","energy");
    ietaRH = _e.a().get_i("rh","ix");
    iphiRH = _e.a().get_i("rh","iy");
    izRH   = _e.a().get_i("rh","iz");
    // cut on energy
    if ( fabs(eRH)<0.080 ) continue;
    // cut on deltaR
    if(izRH==0) { //Barrel
      MEEBGeom::EtaPhiPoint point_ = MEEBDisplay::center( ietaRH, iphiRH );
      etaRH = point_.first;
      phiRH = point_.second * Constants::pi;
    }
    else { //Endcaps
      MEEEGeom::EtaPhiPoint point_ = MEEEDisplay::center( ietaRH, iphiRH, izRH );
      etaRH = point_.first;
      phiRH = point_.second * Constants::pi;
    }
    for (int k=0; k<(int)keepEta.size(); k++) {
      etaK = keepEta[k];
      phiK = keepPhi[k];
      dRtmp = KineUtils::dR( etaRH, etaK, phiRH, phiK );
      if (dRtmp<0.6) keep = true;
    }
    if (keep) rhVector.push_back(irh);
  }
  keepEta.clear();
  keepPhi.clear();
  return;
}


bool
DiphotonAnalysis::jetInCone(Candidate* can)
{
  bool jic = false;
  float jetptmin = 10.;
  float drcut = .5;
  float dRtmp = 20.;
  
  const CandList& Jets = _e.jetList( EventManager::kPfJet);
  for(int unsigned ij=0;ij<Jets.size();ij++) {
    Candidate* jet = Jets[ij];
    dRtmp = KineUtils::dR(jet->eta(), can->eta(), jet->phi(), can->phi() );
    if( jet->pt()>jetptmin && dRtmp<drcut ) { jic=true; break; }
  }
  return jic;
}

float 
DiphotonAnalysis::cosThetaStar(TLorentzVector P1, TLorentzVector P2)
{
  float costhetastar=0;
  
  TLorentzVector Pgg = P1+P2;
  
  double e_gg = Pgg.Energy();
  double m_gg = Pgg.Mag();
  double p_gg = sqrt(e_gg*e_gg-m_gg*m_gg);
  
  double px_gg = Pgg.Px();
  double py_gg = Pgg.Py();
  double pz_gg = Pgg.Pz();
  
  double beta_b = p_gg/e_gg;
  double gamma_b = e_gg/m_gg;
  
  double xnx = px_gg/p_gg;
  double xny = py_gg/p_gg;
  double xnz = pz_gg/p_gg;
  
  double xwx = P1.Py()*xnz - P1.Pz()*xny;
  double xwy = P1.Pz()*xnx - P1.Px()*xnz;
  double xwz = P1.Px()*xny - P1.Py()*xnx;
  
  double sin_theta_3 = sqrt(xwx*xwx+xwy*xwy+xwz*xwz)/P1.Energy();
  double cos_theta_3 = (P1.Px()*xnx+P1.Py()*xny+P1.Pz()*xnz)/P1.Energy();
  
  double tg_thetas_lhc = sin_theta_3/(gamma_b*(cos_theta_3-beta_b));
  
  costhetastar = 1/sqrt(1+tg_thetas_lhc*tg_thetas_lhc);
  
  return costhetastar;
}
 


