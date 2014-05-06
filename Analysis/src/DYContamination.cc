#include "Analysis/src/DYContamination.hh"

#include <cassert>
#include <algorithm>
#include <sstream>

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
#include "Analysis/selectors/RandomConeSelector.hh"

ClassImp( DYContamination );

DYContamination::DYContamination( Sample& sample, const std::string & collectionFileName ) :
    SampleAnalysis( "Diboson", sample, collectionFileName )
{
    cout << "\t-------------------------------------------" << endl;
    cout << "\t---     DY contamination analysis      ----" << endl;
    cout << "\t-------------------------------------------" << endl;

    _hlt = false; // apply HLT selection
    _nDebug = 10;
    _nevts = 0;
    _nSelectedZ = 0;
    _tm.open("ntuple.root", TreeManager::kWrite);
    _tm.addTree("ntu");
}

DYContamination::~DYContamination()
{
    _tm.save();
}

void
DYContamination::bookHistograms()
{
    defineTemplate( "Eproj", 400, -200, 200 );
    defineTemplate( "Eproj2D", 200, -200, 200, 200, -200, 200 );
    defineTemplate( "METVsMT", 80, 0, 80, 80, 0, 160 );
    defineTemplate( "sigMET",  200, 0, 20    );
    defineTemplate( "balance", 200, 0, 20    );
    defineTemplate( "jmult", 20, -0.5, 19.5 );
    defineTemplate( "LP", 100, -2., 2. );
}

void
DYContamination::writeHistograms()
{
    printf("%d event(s) analysed.\n", _nevts);
    printf("%d Z candidate(s) found.\n", _nSelectedZ);
}


bool
DYContamination::analyzeEvent()
{
    ++_nevts;
    if (_nevts % 1000 == 0) printf("%d\n", _nevts);

    int nSel = 0;
    // build lists of leptons and prepare the MC matching map
    buildLeptonLists();

    // build the combinatorics of Z candidates
    CandList zCand;

    initRHV();
    makeZCandList(zCand);

    if (zCand.size() == 0) return 0;

    ++_nSelectedZ;

    //_tm.sync("ntu");
    return nSel>0;
}


void
DYContamination::initRHV()
{
    _rhV.clear();
//    bool keep = false;
//    CandList photons = _e.photons();
//    float eRH, etaRH, phiRH, etaK, phiK, dRtmp;
//    int ietaRH, iphiRH, izRH;
//    size_t nRecHits = (size_t)_e.a().n("rh");
//    printf("--> rh %d\n", nRecHits);
//    for(size_t irh = 0; irh < nRecHits; ++irh)
//    {
//        keep = false;
//        assert( _e.a().load("rh", irh) );
//        eRH = _e.a().get_f("rh", "energy");
//        ietaRH = _e.a().get_i("rh", "ix");
//        iphiRH = _e.a().get_i("rh", "iy");
//        izRH   = _e.a().get_i("rh", "iz");
//        // cut on energy
//        //if ( fabs(eRH)<0.080 ) continue;
//        // cut on deltaR
//        if(izRH==0) { //Barrel
//            MEEBGeom::EtaPhiPoint point_ = MEEBDisplay::center( ietaRH, iphiRH );
//            etaRH = point_.first;
//            phiRH = point_.second * Constants::pi;
//        }
//        else { //Endcaps
//            MEEEGeom::EtaPhiPoint point_ = MEEEDisplay::center( ietaRH, iphiRH, izRH );
//            etaRH = point_.first;
//            phiRH = point_.second * Constants::pi;
//        }
//        for (size_t k = 0; k < photons.size(); ++k) {
//            etaK = photons[k]->eta();
//            phiK = photons[k]->phi();
//            printf("--> %f %f\n", etaK, phiK);
//            dRtmp = KineUtils::dR( etaRH, etaK, phiRH, phiK );
//            if (dRtmp < 0.6) keep = true;
//        }
//        if (keep) _rhV.push_back(irh);
//    }
}


void
DYContamination::makeZCandList(CandList & cands)
{
    CandList eles = _e.photons();
    for (size_t i = 0; i < eles.size(); ++i) {
        Candidate * ei = eles[i];
        //if (!selectElectron(ei)) continue;
        for (size_t j = i + 1; j < eles.size(); ++j) {
            Candidate * ej = eles[j];
            _tm.branchio<float>("ntu", "run", _e.run());
            _tm.branchio<float>("ntu", "event", _e.event());
            //if (!selectElectron(ej)) continue;
            //printf("new cand: %d %d --> %f %f\n", i, j, ei->pt(), ej->pt());
            //if (ei->charge() * ej->charge() >= 0) continue;
            Candidate * Z = Candidate::create(ei, ej);
            //printf("    mass: %f\n", Z->mass());
            cands.push_back(Z);
            _tm.branchio<float>("ntu", "mZ", Z->mass());
            _tm.branchio<float>("ntu", "ptZ", Z->pt());
            saveElectron(ei, 1, i);
            saveElectron(ej, 2, j);
            _tm.sync("ntu");
        }
    }
}


    void
DYContamination::saveElectron(Candidate * c, int ionetwo, int iele)
{
    char name[64];
    sprintf(name, "index%d", ionetwo);
    _tm.branchio<int>("ntu", name, ionetwo);

    CandInfo * ci = c->info();

    // acceptance
    //if (c->pt() < 20) return false;
    sprintf(name, "pt%d", ionetwo);
    _tm.branchio<float>("ntu", name, c->pt());
    sprintf(name, "eta%d", ionetwo);
    _tm.branchio<float>("ntu", name, c->eta());
    sprintf(name, "phi%d", ionetwo);
    _tm.branchio<float>("ntu", name, c->phi());

    float eta    = ci->getFloat("caloEta");
    //if (fabs(eta) > 2.5 || ((eta) > 1.442 && fabs(eta) < 1.556)) return false;
    sprintf(name, "ecalEta%d", ionetwo);
    _tm.branchio<float>("ntu", name, eta);

    float phi    = ci->getFloat("caloPhi");
    sprintf(name, "ecalPhi%d", ionetwo);
    _tm.branchio<float>("ntu", name, phi);

    // pixel seed
    bool haspix  = ci->getBool("hasPixelSeed");
    //if (haspix) return false;
    sprintf(name, "hasPixelSeed%d", ionetwo);
    _tm.branchio<int>("ntu", name, haspix);


    bool isEB    = ci->getBool("isEB");
    sprintf(name, "isEB%d", ionetwo);
    _tm.branchio<int>("ntu", name, isEB);

    float sigeta = ci->getFloat("sigmaIetaIeta");
    //if (sigeta >= (isEB ? 0.011 : 0.028)) return false;
    sprintf(name, "sieie%d", ionetwo);
    _tm.branchio<float>("ntu", name, sigeta);

    float hovere = ci->getFloat("hadronicOverEm");
    sprintf(name, "hoe%d", ionetwo);
    _tm.branchio<float>("ntu", name, hovere);

    float eiso   = ci->getFloat("ecalTowerSumEtConeDR04");
    sprintf(name, "eiso%d", ionetwo);
    _tm.branchio<float>("ntu", name, eiso);

    float hiso   = ci->getFloat("hcalTowerSumEtConeDR04");
    //if (hiso >= (isEB ? 2. : 4.)) return false;
    sprintf(name, "hiso%d", ionetwo);
    _tm.branchio<float>("ntu", name, hiso);

    float tiso   = ci->getFloat("trkSumPtHollowConeDR04");
    //if (tiso >= (isEB ? 2. : 4.)) return false;
    sprintf(name, "tiso%d", ionetwo);
    _tm.branchio<float>("ntu", name, tiso);

    // niceTracks
    ///////////ConeSelector sel(*c, 0.1);
    ///////////CandList tkL;
    /////////////const CandList & tk_list = _e.tracks();
    ///////////sel.getList(_e.tracks(), tkL);
    ///////////printf("--> tkL = %d\n", tkL.size());
    ///////////float niceTkPt = 3.;
    ///////////for (size_t i = 0; i < tkL.size(); ++i) {
    ///////////    if (tkL[i]->pt() < niceTkPt) continue;
    ///////////    const CandInfo * ci = tkL[i]->info();
    ///////////    float tk_eta = ci->getFloat("impEta");
    ///////////    float tk_phi = ci->getFloat("impPhi");
    ///////////    float dR  = KineUtils::dR(tk_eta, eta, tk_phi, phi);
    ///////////    printf("--> dR = %f %f -- %f %f --> %f\n", tk_eta, tk_phi, eta, phi, dR);
    ///////////    if (dR > 0.4) continue;
    ///////////    printf("--> passo qui: %f %f\n", c->pt(), tkL[i]->pt());
    ///////////}
    //////RandomConeSelector isoSel(_rhV);
    //////Candidate * sc = _e.getSecond("photon-SuperCluster", c);
    //////isoSel.computeIso(*c, *sc, false);
    //////std::vector<float> ecalIso_Vtkpt = isoSel.getVtkpt();
    //////std::vector<bool> ecalIso_VtkSV = isoSel.getVtkSV();
    //////std::vector<bool> ecalIso_VtkFP = isoSel.getVtkFP();
    //////float niceTkPt = 0.;
    //////int nImpTk = 0;
    //////printf("--> %d\n", ecalIso_Vtkpt.size());
    //////for(size_t t = 0; t < ecalIso_Vtkpt.size(); ++t) {
    //////    if(ecalIso_VtkSV[t] == false || ecalIso_VtkFP[t] == false || ecalIso_Vtkpt[t] < niceTkPt) continue;
    //////    ++nImpTk;
    //////}
    //////sprintf(name, "impTk%d", ionetwo);
    //////_tm.branchio<int>("ntu", name, nImpTk);

    //return true;
}

    bool
DYContamination::passHlt()
{
    //
    // This is data
    //
    bool fired = false;
    int run = _e.run();

    // electron triggers, 2010
    if( ( (run <= 140401)                   && _e.isFired("HLT_Ele15_LW_L1R") ) ||
        ( (run >= 140402 && run <= 143962 ) && _e.isFired("HLT_Ele15_SW_L1R") ) ||
        ( (run >= 143963 && run <= 144114 ) && _e.isFired("HLT_Ele15_SW_CaloEleId_L1R") ) ||
        ( (run >= 144115 && run <= 147116 ) && _e.isFired("HLT_Ele17_SW_CaloEleId_L1R") ) ||
        ( (run >= 147117 && run <= 148058 ) && _e.isFired("HLT_Ele17_SW_TightEleId_L1R") ) ||
        ( (run >= 148059 && run <= 149064 ) && _e.isFired("HLT_Ele17_SW_TighterEleIdIsol_L1R_v2") ) ||
        ( (run >= 149065 && run <= 149442 ) && _e.isFired("HLT_Ele17_SW_TighterEleIdIsol_L1R_v3") ) )
        fired =true;

    // muon triggers, 2010

    if( ( (run <= 147116) &&  _e.isFired("HLT_Mu9") ) ||
        ( (run >= 147117) && (  _e.isFired("HLT_Mu15_v1") || _e.isFired("HLT_Mu15_v2") ) ) )
        fired =true;

    if( EventServer::isMC )
    { if( _e.isFired("HLT_Ele15_LW_L1R") || _e.isFired("HLT_Ele15_SW_L1R") || _e.isFired("HLT_Ele15_SW_CaloEleId_L1R") || _e.isFired("HLT_Ele17_SW_CaloEleId_L1R") || _e.isFired("HLT_Ele17_SW_TightEleId_L1R") ||
          _e.isFired("HLT_Ele17_SW_TighterEleIdIsol_L1R_v2") ||
          _e.isFired("HLT_Ele17_SW_TighterEleIdIsol_L1R_v3") ||
          _e.isFired("HLT_Mu9") ||
          _e.isFired("HLT_Mu15_v1") ||
          _e.isFired("HLT_Mu15_v2") )
        fired=true;
    }

    return fired;
}
