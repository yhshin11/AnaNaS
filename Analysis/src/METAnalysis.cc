#include "Analysis/src/METAnalysis.hh"

#include <sstream>

using namespace std;



ClassImp( METAnalysis )

  typedef vector<float> vectorFloat;


METAnalysis::METAnalysis( Sample& sample, EventManager& manager ):
			  // Sample& sample, const std::string& collectionFileName ):
SampleAnalysis( "MET", sample, manager ),
MuonSel(MuonSelector::kPt20,MuonSelector::kTight),
ElectronSel(MuonSelector::kPt20,MuonSelector::kMedium)
{
  cout << "\t-------------------------------------------" << endl; 
  cout << "\t---            MET  analysis           ----" << endl; 
  cout << "\t-------------------------------------------" << endl; 
  
  mainVtx= NULL;
  _ZCand= NULL;

  Nevt_ = 0;
  SampleName = sample.name();
  _nDebug = 0;


  LoadJetSmearingDB();
  smearBy_=1.0; //FIXME 1
  shiftBy_=0.;

  _useMu=false;

  _verbose=true;

}

METAnalysis::~METAnalysis() {
}


void
METAnalysis::bookHistograms() {
}

void
METAnalysis::writeHistograms() {
}

bool
METAnalysis::analyzeEvent() {
  Nevt_++;
  if(Nevt_%1000==0) cout<<"Nevt = "<<Nevt_<<endl;

  

  //Yeah!!!
  Candidate* ZCand = BuildBestZCandMultiChan();  //BuildBestZCand();
  
  if(ZCand==NULL) return false;
  mainVtx = ZCand->daughter(0)->vertex();
   _ZCand = ZCand;
   //Now MET
  


   //FIXME : atlas looks for events with no jets
   //  const CandList& jets = _e.jetList( EventManager::kPatJet);
   //  for(size_t ij=0;ij<jets.size();ij++) {
   
   //    const Candidate* j=jets[ij];
   //    if(j->pt()>20) return false;
   //  }

   //if(!HLTFired() ) return false;
   if(!HLTFiredMultiChan() ) return false;
  

  // count the number of "good" vertices
  const VtxList& vertices = _e.vertices(); 
  bool findPV=false;
  int nVertex=0;
  for( size_t iv=0; iv<vertices.size(); iv++ ) 
    {
      Vertex* vert = vertices[iv];
      CandInfo* infoV = vert->info();
      bool goodV=false;
      if( fabs((vert->pos()).Z())<24 && 
	  (vert->pos()).Perp()<2     && 
	  infoV->getFloat("ndof")> 4 && 
	  !(infoV->getBool("isFake"))    )
	goodV=true;
      
      if(!findPV && goodV )
	{ nVertex++; findPV=true; continue; }
      
      if(goodV && findPV && vertices.size()>1) {
	nVertex++;
      }
    }
      
  //cout<<findPV<<endl;
  if(!findPV) return false;
  

  TLorentzVector cL1 = pfConeMatch( ZCand->daughter(0) );
  TLorentzVector cL2 = pfConeMatch( ZCand->daughter(1) );
  //if(cL1==NULL || cL2==NULL) return false;

  //FIXME
  //const Candidate* MET = _e.met( "pat", "patPFMetNoPileUp" );
  //if(MET->pt()<150) return false;

  //store useful variables
  tm.setTree( "metTree", "" );
 
  int run = _e.run();
  int event = _e.event();
  int ls = _e.lumisec();
  tm.add< int >(    "run",   &run       );
  tm.add< int >(    "event", &event     );
  tm.add< int >(    "ls", &ls     );
  
  float mass=ZCand->mass();
  float pt = ZCand->pt();
  float phi =ZCand->phi();
  tm.add<float>( "massZ", &mass );
  tm.add<float>( "ptZ", &pt );
  tm.add<float>( "phiZ", &phi );
  
  float pt_m1 = ZCand->daughter(0)->pt();
  float eta_m1 = ZCand->daughter(0)->eta();
  float phi_m1 = ZCand->daughter(0)->phi();

  float pt_m2 = ZCand->daughter(1)->pt();
  float eta_m2 = ZCand->daughter(1)->eta();
  float phi_m2 = ZCand->daughter(1)->phi();
  
  tm.add<float>( "pt_m1", &pt_m1 );
  tm.add<float>( "eta_m1", &eta_m1 );
  tm.add<float>( "phi_m1", &phi_m1 );

  tm.add<float>( "pt_m2", &pt_m2 );
  tm.add<float>( "eta_m2", &eta_m2 );
  tm.add<float>( "phi_m2", &phi_m2 );
 
  int ch_m1=ZCand->daughter(0)->charge();
  int ch_m2=ZCand->daughter(1)->charge();
  tm.add< int >( "ch_m1", &ch_m1 );
  tm.add< int >( "ch_m2", &ch_m2 );

  int pdgCode_m1=ZCand->daughter(0)->pdgCode();
  int pdgCode_m2=ZCand->daughter(1)->pdgCode();
  tm.add< int >( "pdgId_m1", &pdgCode_m1 );
  tm.add< int >( "pdgId_m2", &pdgCode_m2 );

  
  float pt_cm1 =cL1.Pt();
  float eta_cm1 =cL1.Eta();
  float phi_cm1 =cL1.Phi();
  
  float pt_cm2 =cL2.Pt();
  float eta_cm2 =cL2.Eta();
  float phi_cm2 =cL2.Phi();
  
  
  //cout<<pt_cm1<<"   "<<pt_cm2<<endl;

  tm.add<float>( "pt_cm1", &pt_cm1 );
  tm.add<float>( "eta_cm1", &eta_cm1 );
  tm.add<float>( "phi_cm1", &phi_cm1 );

  tm.add<float>( "pt_cm2", &pt_cm2 );
  tm.add<float>( "eta_cm2", &eta_cm2 );
  tm.add<float>( "phi_cm2", &phi_cm2 );

  int njet=0,nbtag=0;
  float ptj1,ptj2,etaj1,etaj2,phij1,phij2;
  float tmp=0,tmp2=-1;
  
  
  vectorFloat* jetPts=new vectorFloat;
  vectorFloat* jetEtas=new vectorFloat;
  vectorFloat* jetPhis=new vectorFloat;

  TVector2 v(0,0);
  TVector2 vl(0,0);
  TVector2 vvl(0,0);

  const CandList& jets = _e.jetList( EventManager::kPatJet);
  for(unsigned int ijet=0;ijet<jets.size();ijet++) {
    const Candidate* jet = jets[ijet];

    if(isLepton(jet)) continue;

    if(jet->pt()>5) {
      jetPts->push_back(jet->pt());
      jetEtas->push_back(jet->eta());
      jetPhis->push_back(jet->phi());

      if(jet->pt()<30) {
	vl += jet->p2();
      }

    }
    else {
      vvl += jet->p2();
    }

    if(jet->pt()>30) {
      njet++;

      CandInfo* info = jet->info();
      float tkC = info->getFloat("btagCSVBT");
      v += jet->p2();

      if(tkC>0.244) { nbtag++;
          
	if(jet->pt()>tmp) {
	  tmp = jet->pt();
	  ptj1 = jet->pt();
	  etaj1 = jet->eta();
	  phij1 = jet->phi();
	}
	if(jet->pt()>tmp2 && jet->pt()<tmp) {
	  tmp2 = jet->pt();
	  ptj2 = jet->pt();
	  etaj2 = jet->eta();
	  phij2 = jet->phi();
	}
      }
    }

 } 
  
  tm.add<int>( "njet", &njet );
  tm.add<int>( "nbtag", &nbtag );

  tm.add<float>( "pt_j1", &ptj1 );
  tm.add<float>( "eta_j1", &etaj1 );
  tm.add<float>( "phi_j1", &phij1 );
  tm.add<float>( "pt_j2", &ptj2 );
  tm.add<float>( "eta_j2", &etaj2 );
  tm.add<float>( "phi_j2", &phij2 );

  float vj_pt=v.Mod();
  float vj_phi=v.Phi();
  float vlj_pt=vl.Mod();
  float vlj_phi=vl.Phi();
  float vvlj_pt=vvl.Mod();
  float vvlj_phi=vvl.Phi();

  
  tm.add<float>( "vj_pt",  &vj_pt); 
  tm.add<float>( "vj_phi", &vj_phi); 
  tm.add<float>( "vlj_pt",  &vlj_pt); 
  tm.add<float>( "vlj_phi", &vlj_phi); 
  tm.add<float>( "vvlj_pt",  &vvlj_pt); 
  tm.add<float>( "vvlj_phi", &vvlj_phi); 

  tm.add<vectorFloat*>("jetPts", & jetPts);
  tm.add<vectorFloat*>("jetEtas", & jetEtas);
  tm.add<vectorFloat*>("jetPhis", & jetPhis);

  
  //now pfcandidates
  Candidate* pf1 = pfMatch( ZCand->daughter(0) );
  Candidate* pf2 = pfMatch( ZCand->daughter(1) );

  float pfpt_m1 = pf1?pf1->pt():-1;
  float pfeta_m1 = pf1?pf1->eta():-1;
  float pfphi_m1 = pf1?pf1->phi():-1;
  
  float pfpt_m2 = pf2?pf2->pt():-1;
  float pfeta_m2 = pf2?pf2->eta():-1;
  float pfphi_m2 = pf2?pf2->phi():-1;


  // and total cone
  TLorentzVector Zpart = pfConeMatch( ZCand->daughter(0) );
  Zpart += pfConeMatch( ZCand->daughter(1) );

  float ptZpart=Zpart.Pt();
  float phiZpart=Zpart.Phi();
  
  
  tm.add<float>( "pfpt_m1", &pfpt_m1 );
  tm.add<float>( "pfeta_m1", &pfeta_m1 );
  tm.add<float>( "pfphi_m1", &pfphi_m1 );

  tm.add<float>( "pfpt_m2", &pfpt_m2 );
  tm.add<float>( "pfeta_m2", &pfeta_m2 );
  tm.add<float>( "pfphi_m2", &pfphi_m2 );

  tm.add<float>( "ptZpart", &ptZpart );
  tm.add<float>( "phiZpart", &phiZpart );
  
  if( _e.a().isMC ) {
    const Candidate* genMET = _e.met( EventManager::kGenMet );
    float genMETPt = genMET->pt() ;
    float genMETPhi = genMET->phi() ;
    
    tm.add<float>("genMETPt",&genMETPt );
    tm.add<float>("genMETPhi",&genMETPhi );
 
  }

//   float factVtx = 1./_e.nGoodVertices();
//   float factAtlas = 1./AtlasFactor();
   int nVtx = _e.nGoodVertices();

//   tm.add<float>( "fVtx", &factVtx );
//   tm.add<float>( "fAtlas", &factAtlas );
  tm.add<int>( "nVtx", &nVtx );
  tm.add<int>( "nGQVtx", &nVertex );

  float trueI = _e.getTrueNInt();
  tm.add<float>( "trueNint", &trueI ); 

  string algo[4]={"mva","atlasAssoc","atlas","basic"};
  string factLvl[3]={"","vtx","atlas"};
  string redLvl[3]={"","RedUnc","RedAll"};
  string idLvl[3]={"mvaIdLoose","mvaIdMedium","mvaIdTight"};
  float pts[6]={10.,15.,20.,25.,30.,35.};

 //  for(size_t ia=0;ia<4;ia++) {
//     for(size_t ir=0;ir<3;ir++) {
//       for(size_t ic=0;ic<3;ic++) {
	
// 	if( (ir==0 && ic!=0 && ia!=2) || (ir!=0 && ic==0 && ia!=2) ) continue;

// 	if(ia==2) {
// 	  	  if(ic==0 && ir==0)
// 	   for(int it=0;it<6;it++)
// 	    FillMET(algo[ia],"","atlas",pts[it],"");
// 	}
// 	else
// 	  if(ia!=0)
// 	    FillMET(algo[ia],redLvl[ir],factLvl[ic],0,"");
// 	  else {
// 	    for(int id=0;id<3;id++)
// 	      FillMET(algo[ia],redLvl[ir],factLvl[ic],0,idLvl[id]);
//       }
//       }
//     }
//   }
  
//   for(int it=0;it<6;it++) {
//     FillMET("atlasred","mva","atlas",pts[it]);
//     FillMET("atlasred","atlas","atlas",pts[it]);
//     FillMET("atlasred","mva","vtx",pts[it]);
//     FillMET("atlasred","atlas","vtx",pts[it]);
//   }


  //cleanPUMET
  //FillMET("clean","","vtx");

  //cout<<" new event "<<endl;
  
  //FillMET("Calo","","");
  FillMET("Calo","caloType1CorrectedMet","");
  FillMET("Calo", "caloType1CorrectedMetJetEnDown", "");
  FillMET("Calo", "caloType1CorrectedMetJetEnUp", "");
  FillMET("Calo", "caloType1CorrectedMetUnclusteredEnDown", "");
  FillMET("Calo", "caloType1CorrectedMetUnclusteredEnUp", "");
  //FillMET("pat", "caloType1p2CorrectedMet", "");

  FillMET("pf","","");
  FillMET("pf","pfMet","");
  
  //cleaned mets

  // FillMET("pf","smearedPFType1CorrectedMetRECO","");
  //if(SampleName.substr(0,5)!="TTbar") {
 //  FillMET("pf","pfType1PhiCorrectedMet","");
//   FillMET("pf","pfType1p0CorrectedMet","");
//   FillMET("pf","pfType1p0PhiCorrectedMet","");
//   FillMET("pf","pfType1p0p2CorrectedMet","");
//   FillMET("pf","pfType1p0p2PhiCorrectedMet","");
//   FillMET("pf","smearedPFType1PhiCorrectedMetRECO","");
//   FillMET("pf","smearedPFType1p0CorrectedMetRECO","");
//   FillMET("pf","smearedPFType1p0PhiCorrectedMetRECO","");
  
//   FillMET("pat","patPFMetNoPileUp","");
//   FillMET("pat","patPFMetMVA","");
  
//   FillMET("pat","patPFMetMVAJetEnDown","");
//   FillMET("pat","patPFMetMVAJetEnUp","");
//   FillMET("pat","patPFMetMVAJetResDown","");
//   FillMET("pat","patPFMetMVAJetResUp","");
//   FillMET("pat","patPFMetMVAUnclusteredEnDown","");
//   FillMET("pat","patPFMetMVAUnclusteredEnUp","");
//   FillMET("pat","patPFMetMVAMuonEnDown","");
//   FillMET("pat","patPFMetMVAMuonEnUp","");

//   FillMET("pat","patPFMetNoPileUpJetEnDown","");
//   FillMET("pat","patPFMetNoPileUpJetEnUp","");
//   FillMET("pat","patPFMetNoPileUpJetResDown","");
//   FillMET("pat","patPFMetNoPileUpJetResUp","");
//   FillMET("pat","patPFMetNoPileUpMuonEnDown","");
//   FillMET("pat","patPFMetNoPileUpMuonEnUp","");
//   FillMET("pat","patPFMetNoPileUpUnclusteredEnDown","");
//   FillMET("pat","patPFMetNoPileUpUnclusteredEnUp","");

//     // }
 

//   FillMET("pat", "patType1CorrectedPFMet", "");
//   FillMET("pat", "patType1p2CorrectedPFMet", "");
//   FillMET("pat", "patType1CorrectedPFMetElectronEnDown", "");
//   FillMET("pat", "patType1CorrectedPFMetElectronEnUp", "");
//   FillMET("pat", "patType1CorrectedPFMetJetEnDown", "");
//   FillMET("pat", "patType1CorrectedPFMetJetEnUp", "");
//   FillMET("pat", "patType1CorrectedPFMetJetResDown", "");
//   FillMET("pat", "patType1CorrectedPFMetJetResUp", "");
//   FillMET("pat", "patType1CorrectedPFMetMuonEnUp", "");
//   FillMET("pat", "patType1CorrectedPFMetMuonEnDown", "");
//   FillMET("pat", "patType1CorrectedPFMetTauEnDown", "");
//   FillMET("pat", "patType1CorrectedPFMetTauEnUp", "");
//   FillMET("pat", "patType1CorrectedPFMetUnclusteredEnDown", "");
//   FillMET("pat", "patType1CorrectedPFMetUnclusteredEnUp", "");


 //  FillMET("pat", "patPFMetNoSmear", "");
//   FillMET("pat", "patType1CorrectedPFMetNoSmear", "");
//   FillMET("pat", "patType1CorrectedPFMetElectronEnDownNoSmear", "");
//   FillMET("pat", "patType1CorrectedPFMetElectronEnUpNoSmear", "");
//   FillMET("pat", "patType1CorrectedPFMetJetEnDownNoSmear", "");
//   FillMET("pat", "patType1CorrectedPFMetJetEnUpNoSmear", "");
//   FillMET("pat", "patType1CorrectedPFMetJetResDownNoSmear", "");
//   FillMET("pat", "patType1CorrectedPFMetJetResUpNoSmear", "");
//   FillMET("pat", "patType1CorrectedPFMetMuonEnDownNoSmear", "");
//   FillMET("pat", "patType1CorrectedPFMetMuonEnUpNoSmear", "");
//   FillMET("pat", "patType1CorrectedPFMetTauEnDownNoSmear", "");
//   FillMET("pat", "patType1CorrectedPFMetTauEnUpNoSmear", "");
//   FillMET("pat", "patType1CorrectedPFMetUnclusteredEnDownNoSmear", "");
//   FillMET("pat", "patType1CorrectedPFMetUnclusteredEnUpNoSmear", "");
//   FillMET("pat", "patPFMetSmeared", "");
//   FillMET("pat", "patType1CorrectedPFMetSmeared", "");
//   FillMET("pat", "patType1CorrectedPFMetElectronEnDownSmeared", "");
//   FillMET("pat", "patType1CorrectedPFMetElectronEnUpSmeared", "");
//   FillMET("pat", "patType1CorrectedPFMetJetEnDownSmeared", "");
//   FillMET("pat", "patType1CorrectedPFMetJetEnUpSmeared", "");
//   FillMET("pat", "patType1CorrectedPFMetJetResDownSmeared", "");
//   FillMET("pat", "patType1CorrectedPFMetJetResUpSmeared", "");
//   FillMET("pat", "patType1CorrectedPFMetMuonEnDownSmeared", "");
//   FillMET("pat", "patType1CorrectedPFMetMuonEnUpSmeared", "");
//   FillMET("pat", "patType1CorrectedPFMetTauEnDownSmeared", "");
//   FillMET("pat", "patType1CorrectedPFMetTauEnUpSmeared", "");
//   FillMET("pat", "patType1CorrectedPFMetUnclusteredEnDownSmeared", "");
//   FillMET("pat", "patType1CorrectedPFMetUnclusteredEnUpSmeared", "");
//   FillMET("pat", "patPFMetNoPileUpSmeared", "");
//   FillMET("pat", "patPFMetMVASmeared", "");
//   FillMET("pat", "patPFMetMVAJetEnDownSmeared", "");
//   FillMET("pat", "patPFMetMVAJetEnUpSmeared", "");
//   FillMET("pat", "patPFMetMVAJetResDownSmeared", "");
//   FillMET("pat", "patPFMetMVAJetResUpSmeared", "");
//   FillMET("pat", "patPFMetMVAMuonEnDownSmeared", "");
//   FillMET("pat", "patPFMetMVAMuonEnUpSmeared", "");
//   FillMET("pat", "patPFMetMVAUnclusteredEnDownSmeared", "");
//   FillMET("pat", "patPFMetMVAUnclusteredEnUpSmeared", "");
//   FillMET("pat", "patPFMetNoPileUpJetEnDownSmeared", "");
//   FillMET("pat", "patPFMetNoPileUpJetEnUpSmeared", "");
//   FillMET("pat", "patPFMetNoPileUpJetResDownSmeared", "");
//   FillMET("pat", "patPFMetNoPileUpJetResUpSmeared", "");
//   FillMET("pat", "patPFMetNoPileUpMuonEnDownSmeared", "");
//   FillMET("pat", "patPFMetNoPileUpMuonEnUpSmeared", "");
//   FillMET("pat", "patPFMetNoPileUpUnclusteredEnDownSmeared", "");
//   FillMET("pat", "patPFMetNoPileUpUnclusteredEnUpSmeared", "");
//   FillMET("pat", "patPFMetNoPileUpNoSmear", "");
//   FillMET("pat", "patPFMetMVANoSmear", "");
//   FillMET("pat", "patPFMetMVAJetEnDownNoSmear", "");
//   FillMET("pat", "patPFMetMVAJetEnUpNoSmear", "");
//   FillMET("pat", "patPFMetMVAJetResDownNoSmear", "");
//   FillMET("pat", "patPFMetMVAJetResUpNoSmear", "");
//   FillMET("pat", "patPFMetMVAMuonEnDownNoSmear", "");
//   FillMET("pat", "patPFMetMVAMuonEnUpNoSmear", "");
//   FillMET("pat", "patPFMetMVAUnclusteredEnDownNoSmear", "");
//   FillMET("pat", "patPFMetMVAUnclusteredEnUpNoSmear", "");
//   FillMET("pat", "patPFMetNoPileUpJetEnDownNoSmear", "");
//   FillMET("pat", "patPFMetNoPileUpJetEnUpNoSmear", "");
//   FillMET("pat", "patPFMetNoPileUpJetResDownNoSmear", "");
//   FillMET("pat", "patPFMetNoPileUpJetResUpNoSmear", "");
//   FillMET("pat", "patPFMetNoPileUpMuonEnDownNoSmear", "");
//   FillMET("pat", "patPFMetNoPileUpMuonEnUpNoSmear", "");
//   FillMET("pat", "patPFMetNoPileUpUnclusteredEnDownNoSmear", "");
//   FillMET("pat", "patPFMetNoPileUpUnclusteredEnUpNoSmear", "");

//   FillMET("pat", "patPFMetMVAEMuSmeared", "");
//   FillMET("pat", "patPFMetMVAElectronEnDownEMuSmeared", "");
//   FillMET("pat", "patPFMetMVAElectronEnUpEMuSmeared", "");
//   FillMET("pat", "patPFMetMVAJetEnDownEMuSmeared", "");
//   FillMET("pat", "patPFMetMVAJetEnUpEMuSmeared", "");
//   FillMET("pat", "patPFMetMVAJetResDownEMuSmeared", "");
//   FillMET("pat", "patPFMetMVAJetResUpEMuSmeared", "");
//   FillMET("pat", "patPFMetMVAMuonEnDownEMuSmeared", "");
//   FillMET("pat", "patPFMetMVAMuonEnUpEMuSmeared", "");
//   FillMET("pat", "patPFMetMVATauEnDownEMuSmeared", "");
//   FillMET("pat", "patPFMetMVATauEnUpEMuSmeared", "");
//   FillMET("pat", "patPFMetMVAUnclusteredEnDownEMuSmeared", "");
//   FillMET("pat", "patPFMetMVAUnclusteredEnUpEMuSmeared", "");
//   FillMET("pat", "patPFMetMVAEleSmeared", "");
//   FillMET("pat", "patPFMetMVAElectronEnDownEleSmeared", "");
//   FillMET("pat", "patPFMetMVAElectronEnUpEleSmeared", "");
//   FillMET("pat", "patPFMetMVAJetEnDownEleSmeared", "");
//   FillMET("pat", "patPFMetMVAJetEnUpEleSmeared", "");
//   FillMET("pat", "patPFMetMVAJetResDownEleSmeared", "");
//   FillMET("pat", "patPFMetMVAJetResUpEleSmeared", "");
//   FillMET("pat", "patPFMetMVAMuonEnDownEleSmeared", "");
//   FillMET("pat", "patPFMetMVAMuonEnUpEleSmeared", "");
//   FillMET("pat", "patPFMetMVATauEnDownEleSmeared", "");
//   FillMET("pat", "patPFMetMVATauEnUpEleSmeared", "");
//   FillMET("pat", "patPFMetMVAUnclusteredEnDownEleSmeared", "");
//   FillMET("pat", "patPFMetMVAUnclusteredEnUpEleSmeared", "");
//   FillMET("pat", "patPFMetNoPileUpEMuSmeared", "");
//   FillMET("pat", "patPFMetNoPileUpElectronEnDownEMuSmeared", "");
//   FillMET("pat", "patPFMetNoPileUpElectronEnUpEMuSmeared", "");
//   FillMET("pat", "patPFMetNoPileUpJetEnDownEMuSmeared", "");
//   FillMET("pat", "patPFMetNoPileUpJetEnUpEMuSmeared", "");
//   FillMET("pat", "patPFMetNoPileUpJetResDownEMuSmeared", "");
//   FillMET("pat", "patPFMetNoPileUpJetResUpEMuSmeared", "");
//   FillMET("pat", "patPFMetNoPileUpMuonEnDownEMuSmeared", "");
//   FillMET("pat", "patPFMetNoPileUpMuonEnUpEMuSmeared", "");
//   FillMET("pat", "patPFMetNoPileUpTauEnDownEMuSmeared", "");
//   FillMET("pat", "patPFMetNoPileUpTauEnUpEMuSmeared", "");
//   FillMET("pat", "patPFMetNoPileUpUnclusteredEnDownEMuSmeared", "");
//   FillMET("pat", "patPFMetNoPileUpUnclusteredEnUpEMuSmeared", "");
//   FillMET("pat", "patPFMetNoPileUpEleSmeared", "");
//   FillMET("pat", "patPFMetNoPileUpElectronEnDownEleSmeared", "");
//   FillMET("pat", "patPFMetNoPileUpElectronEnUpEleSmeared", "");
//   FillMET("pat", "patPFMetNoPileUpJetEnDownEleSmeared", "");
//   FillMET("pat", "patPFMetNoPileUpJetEnUpEleSmeared", "");
//   FillMET("pat", "patPFMetNoPileUpJetResDownEleSmeared", "");
//   FillMET("pat", "patPFMetNoPileUpJetResUpEleSmeared", "");
//   FillMET("pat", "patPFMetNoPileUpMuonEnDownEleSmeared", "");
//   FillMET("pat", "patPFMetNoPileUpMuonEnUpEleSmeared", "");
//   FillMET("pat", "patPFMetNoPileUpTauEnDownEleSmeared", "");
//   FillMET("pat", "patPFMetNoPileUpTauEnUpEleSmeared", "");
//   FillMET("pat", "patPFMetNoPileUpUnclusteredEnDownEleSmeared", "");
//   FillMET("pat", "patPFMetNoPileUpUnclusteredEnUpEleSmeared", "");
  

  FillMET("pat", "patMETsPF", "");
 
  

  if( !(_e.a().isMC) )  { //data
    
    FillMET("pat", "patPFMetNoSmear", "");
    FillMET("pat", "patType1CorrectedPFMetNoSmear", "");
    FillMET("pat", "patType1CorrectedPFMetElectronEnDownNoSmear", "");
    FillMET("pat", "patType1CorrectedPFMetElectronEnUpNoSmear", "");
    FillMET("pat", "patType1CorrectedPFMetJetEnDownNoSmear", "");
    FillMET("pat", "patType1CorrectedPFMetJetEnUpNoSmear", "");
    FillMET("pat", "patType1CorrectedPFMetJetResDownNoSmear", "");
    FillMET("pat", "patType1CorrectedPFMetJetResUpNoSmear", "");
    FillMET("pat", "patType1CorrectedPFMetMuonEnDownNoSmear", "");
    FillMET("pat", "patType1CorrectedPFMetMuonEnUpNoSmear", "");
    FillMET("pat", "patType1CorrectedPFMetTauEnDownNoSmear", "");
    FillMET("pat", "patType1CorrectedPFMetTauEnUpNoSmear", "");
    FillMET("pat", "patType1CorrectedPFMetUnclusteredEnDownNoSmear", "");
    FillMET("pat", "patType1CorrectedPFMetUnclusteredEnUpNoSmear", "");
    FillMET("pat", "patType1p2CorrectedPFMetNoSmear", "");
    FillMET("pat", "patType1p2CorrectedPFMetElectronEnDownNoSmear", "");
    FillMET("pat", "patType1p2CorrectedPFMetElectronEnUpNoSmear", "");
    FillMET("pat", "patType1p2CorrectedPFMetJetEnDownNoSmear", "");
    FillMET("pat", "patType1p2CorrectedPFMetJetEnUpNoSmear", "");
    FillMET("pat", "patType1p2CorrectedPFMetJetResDownNoSmear", "");
    FillMET("pat", "patType1p2CorrectedPFMetJetResUpNoSmear", "");
    FillMET("pat", "patType1p2CorrectedPFMetMuonEnDownNoSmear", "");
    FillMET("pat", "patType1p2CorrectedPFMetMuonEnUpNoSmear", "");
    FillMET("pat", "patType1p2CorrectedPFMetTauEnDownNoSmear", "");
    FillMET("pat", "patType1p2CorrectedPFMetTauEnUpNoSmear", "");
    FillMET("pat", "patType1p2CorrectedPFMetUnclusteredEnDownNoSmear", "");
    FillMET("pat", "patType1p2CorrectedPFMetUnclusteredEnUpNoSmear", "");
    FillMET("pat", "patPFMetMVANoSmear", "");
    FillMET("pat", "patPFMetMVAElectronEnDownNoSmear", "");
    FillMET("pat", "patPFMetMVAElectronEnUpNoSmear", "");
    FillMET("pat", "patPFMetMVAJetEnDownNoSmear", "");
    FillMET("pat", "patPFMetMVAJetEnUpNoSmear", "");
    FillMET("pat", "patPFMetMVAJetResDownNoSmear", "");
    FillMET("pat", "patPFMetMVAJetResUpNoSmear", "");
    FillMET("pat", "patPFMetMVAMuonEnDownNoSmear", "");
    FillMET("pat", "patPFMetMVAMuonEnUpNoSmear", "");
    FillMET("pat", "patPFMetMVATauEnDownNoSmear", "");
    FillMET("pat", "patPFMetMVATauEnUpNoSmear", "");
    FillMET("pat", "patPFMetMVAUnclusteredEnDownNoSmear", "");
    FillMET("pat", "patPFMetMVAUnclusteredEnUpNoSmear", "");
    
    FillMET("pat", "patPFMetMVAUnityNoSmear", "");
    FillMET("pat", "patPFMetMVAElectronEnDownUnityNoSmear", "");
    FillMET("pat", "patPFMetMVAElectronEnUpUnityNoSmear", "");
    FillMET("pat", "patPFMetMVAJetEnDownUnityNoSmear", "");
    FillMET("pat", "patPFMetMVAJetEnUpUnityNoSmear", "");
    FillMET("pat", "patPFMetMVAJetResDownUnityNoSmear", "");
    FillMET("pat", "patPFMetMVAJetResUpUnityNoSmear", "");
    FillMET("pat", "patPFMetMVAMuonEnDownUnityNoSmear", "");
    FillMET("pat", "patPFMetMVAMuonEnUpUnityNoSmear", "");
    FillMET("pat", "patPFMetMVATauEnDownUnityNoSmear", "");
    FillMET("pat", "patPFMetMVATauEnUpUnityNoSmear", "");
    FillMET("pat", "patPFMetMVAUnclusteredEnDownUnityNoSmear", "");
    FillMET("pat", "patPFMetMVAUnclusteredEnUpUnityNoSmear", "");
    
    FillMET("pat", "patPFMetNoPileUpNoSmear", "");
    FillMET("pat", "patPFMetNoPileUpElectronEnDownNoSmear", "");
    FillMET("pat", "patPFMetNoPileUpElectronEnUpNoSmear", "");
    FillMET("pat", "patPFMetNoPileUpJetEnDownNoSmear", "");
    FillMET("pat", "patPFMetNoPileUpJetEnUpNoSmear", "");
    FillMET("pat", "patPFMetNoPileUpJetResDownNoSmear", "");
    FillMET("pat", "patPFMetNoPileUpJetResUpNoSmear", "");
    FillMET("pat", "patPFMetNoPileUpMuonEnDownNoSmear", "");
    FillMET("pat", "patPFMetNoPileUpMuonEnUpNoSmear", "");
    FillMET("pat", "patPFMetNoPileUpTauEnDownNoSmear", "");
    FillMET("pat", "patPFMetNoPileUpTauEnUpNoSmear", "");
    FillMET("pat", "patPFMetNoPileUpUnclusteredEnDownNoSmear", "");
    FillMET("pat", "patPFMetNoPileUpUnclusteredEnUpNoSmear", "");
    FillMET("pat", "patPFMetMVAEMuNoSmear", "");
    FillMET("pat", "patPFMetMVAElectronEnDownEMuNoSmear", "");
    FillMET("pat", "patPFMetMVAElectronEnUpEMuNoSmear", "");
    FillMET("pat", "patPFMetMVAJetEnDownEMuNoSmear", "");
    FillMET("pat", "patPFMetMVAJetEnUpEMuNoSmear", "");
    FillMET("pat", "patPFMetMVAJetResDownEMuNoSmear", "");
    FillMET("pat", "patPFMetMVAJetResUpEMuNoSmear", "");
    FillMET("pat", "patPFMetMVAMuonEnDownEMuNoSmear", "");
    FillMET("pat", "patPFMetMVAMuonEnUpEMuNoSmear", "");
    FillMET("pat", "patPFMetMVATauEnDownEMuNoSmear", "");
    FillMET("pat", "patPFMetMVATauEnUpEMuNoSmear", "");
    FillMET("pat", "patPFMetMVAUnclusteredEnDownEMuNoSmear", "");
    FillMET("pat", "patPFMetMVAUnclusteredEnUpEMuNoSmear", "");
    FillMET("pat", "patPFMetMVAEleNoSmear", "");
    FillMET("pat", "patPFMetMVAElectronEnDownEleNoSmear", "");
    FillMET("pat", "patPFMetMVAElectronEnUpEleNoSmear", "");
    FillMET("pat", "patPFMetMVAJetEnDownEleNoSmear", "");
    FillMET("pat", "patPFMetMVAJetEnUpEleNoSmear", "");
    FillMET("pat", "patPFMetMVAJetResDownEleNoSmear", "");
    FillMET("pat", "patPFMetMVAJetResUpEleNoSmear", "");
    FillMET("pat", "patPFMetMVAMuonEnDownEleNoSmear", "");
    FillMET("pat", "patPFMetMVAMuonEnUpEleNoSmear", "");
    FillMET("pat", "patPFMetMVATauEnDownEleNoSmear", "");
    FillMET("pat", "patPFMetMVATauEnUpEleNoSmear", "");
    FillMET("pat", "patPFMetMVAUnclusteredEnDownEleNoSmear", "");
    FillMET("pat", "patPFMetMVAUnclusteredEnUpEleNoSmear", "");
    FillMET("pat", "patPFMetNoPileUpEMuNoSmear", "");
    FillMET("pat", "patPFMetNoPileUpElectronEnDownEMuNoSmear", "");
    FillMET("pat", "patPFMetNoPileUpElectronEnUpEMuNoSmear", "");
    FillMET("pat", "patPFMetNoPileUpJetEnDownEMuNoSmear", "");
    FillMET("pat", "patPFMetNoPileUpJetEnUpEMuNoSmear", "");
    FillMET("pat", "patPFMetNoPileUpJetResDownEMuNoSmear", "");
    FillMET("pat", "patPFMetNoPileUpJetResUpEMuNoSmear", "");
    FillMET("pat", "patPFMetNoPileUpMuonEnDownEMuNoSmear", "");
    FillMET("pat", "patPFMetNoPileUpMuonEnUpEMuNoSmear", "");
    FillMET("pat", "patPFMetNoPileUpTauEnDownEMuNoSmear", "");
    FillMET("pat", "patPFMetNoPileUpTauEnUpEMuNoSmear", "");
    FillMET("pat", "patPFMetNoPileUpUnclusteredEnDownEMuNoSmear", "");
    FillMET("pat", "patPFMetNoPileUpUnclusteredEnUpEMuNoSmear", "");
    FillMET("pat", "patPFMetNoPileUpEleNoSmear", "");
    FillMET("pat", "patPFMetNoPileUpElectronEnDownEleNoSmear", "");
    FillMET("pat", "patPFMetNoPileUpElectronEnUpEleNoSmear", "");
    FillMET("pat", "patPFMetNoPileUpJetEnDownEleNoSmear", "");
    FillMET("pat", "patPFMetNoPileUpJetEnUpEleNoSmear", "");
    FillMET("pat", "patPFMetNoPileUpJetResDownEleNoSmear", "");
    FillMET("pat", "patPFMetNoPileUpJetResUpEleNoSmear", "");
    FillMET("pat", "patPFMetNoPileUpMuonEnDownEleNoSmear", "");
    FillMET("pat", "patPFMetNoPileUpMuonEnUpEleNoSmear", "");
    FillMET("pat", "patPFMetNoPileUpTauEnDownEleNoSmear", "");
    FillMET("pat", "patPFMetNoPileUpTauEnUpEleNoSmear", "");
    FillMET("pat", "patPFMetNoPileUpUnclusteredEnDownEleNoSmear", "");
    FillMET("pat", "patPFMetNoPileUpUnclusteredEnUpEleNoSmear", "");

    //FIXME
    FillMET("pat", "patType1CorrectedPFMetEleNoSmear", "");
    FillMET("pat", "patType1CorrectedPFMetEMuNoSmear", "");
    
  }
  else { //MC

    //Brian, temporary
  //   FillMET("pat", "patType1CorrectedPFMetNoSmear", "");
//     FillMET("pat", "patType1CorrectedPFMetElectronEnDownNoSmear", "");
//     FillMET("pat", "patType1CorrectedPFMetElectronEnUpNoSmear", "");
//     FillMET("pat", "patType1CorrectedPFMetJetEnDownNoSmear", "");
//     FillMET("pat", "patType1CorrectedPFMetJetEnUpNoSmear", "");
//     FillMET("pat", "patType1CorrectedPFMetJetResDownNoSmear", "");
//     FillMET("pat", "patType1CorrectedPFMetJetResUpNoSmear", "");
//     FillMET("pat", "patType1CorrectedPFMetMuonEnDownNoSmear", "");
//     FillMET("pat", "patType1CorrectedPFMetMuonEnUpNoSmear", "");
//     FillMET("pat", "patType1CorrectedPFMetTauEnDownNoSmear", "");
//     FillMET("pat", "patType1CorrectedPFMetTauEnUpNoSmear", "");
//     FillMET("pat", "patType1CorrectedPFMetUnclusteredEnDownNoSmear", "");
//     FillMET("pat", "patType1CorrectedPFMetUnclusteredEnUpNoSmear", "");

    

    FillMET("pat", "patPFMetSmeared", "");
    FillMET("pat", "patType1CorrectedPFMetSmeared", "");
    FillMET("pat", "patType1CorrectedPFMetElectronEnDownSmeared", "");
    FillMET("pat", "patType1CorrectedPFMetElectronEnUpSmeared", "");
    FillMET("pat", "patType1CorrectedPFMetJetEnDownSmeared", "");
    FillMET("pat", "patType1CorrectedPFMetJetEnUpSmeared", "");
    FillMET("pat", "patType1CorrectedPFMetJetResDownSmeared", "");
    FillMET("pat", "patType1CorrectedPFMetJetResUpSmeared", "");
    FillMET("pat", "patType1CorrectedPFMetMuonEnDownSmeared", "");
    FillMET("pat", "patType1CorrectedPFMetMuonEnUpSmeared", "");
    FillMET("pat", "patType1CorrectedPFMetTauEnDownSmeared", "");
    FillMET("pat", "patType1CorrectedPFMetTauEnUpSmeared", "");
    FillMET("pat", "patType1CorrectedPFMetUnclusteredEnDownSmeared", "");
    FillMET("pat", "patType1CorrectedPFMetUnclusteredEnUpSmeared", "");
    FillMET("pat", "patType1p2CorrectedPFMetSmeared", "");
    FillMET("pat", "patType1p2CorrectedPFMetElectronEnDownSmeared", "");
    FillMET("pat", "patType1p2CorrectedPFMetElectronEnUpSmeared", "");
    FillMET("pat", "patType1p2CorrectedPFMetJetEnDownSmeared", "");
    FillMET("pat", "patType1p2CorrectedPFMetJetEnUpSmeared", "");
    FillMET("pat", "patType1p2CorrectedPFMetJetResDownSmeared", "");
    FillMET("pat", "patType1p2CorrectedPFMetJetResUpSmeared", "");
    FillMET("pat", "patType1p2CorrectedPFMetMuonEnDownSmeared", "");
    FillMET("pat", "patType1p2CorrectedPFMetMuonEnUpSmeared", "");
    FillMET("pat", "patType1p2CorrectedPFMetTauEnDownSmeared", "");
    FillMET("pat", "patType1p2CorrectedPFMetTauEnUpSmeared", "");
    FillMET("pat", "patType1p2CorrectedPFMetUnclusteredEnDownSmeared", "");
    FillMET("pat", "patType1p2CorrectedPFMetUnclusteredEnUpSmeared", "");
    FillMET("pat", "patPFMetMVASmeared", "");
    FillMET("pat", "patPFMetMVAElectronEnDownSmeared", "");
    FillMET("pat", "patPFMetMVAElectronEnUpSmeared", "");
    FillMET("pat", "patPFMetMVAJetEnDownSmeared", "");
    FillMET("pat", "patPFMetMVAJetEnUpSmeared", "");
    FillMET("pat", "patPFMetMVAJetResDownSmeared", "");
    FillMET("pat", "patPFMetMVAJetResUpSmeared", "");
    FillMET("pat", "patPFMetMVAMuonEnDownSmeared", "");
    FillMET("pat", "patPFMetMVAMuonEnUpSmeared", "");
    FillMET("pat", "patPFMetMVATauEnDownSmeared", "");
    FillMET("pat", "patPFMetMVATauEnUpSmeared", "");
    FillMET("pat", "patPFMetMVAUnclusteredEnDownSmeared", "");
    FillMET("pat", "patPFMetMVAUnclusteredEnUpSmeared", "");
    
    FillMET("pat", "patPFMetMVAUnitySmeared", "");
    FillMET("pat", "patPFMetMVAElectronEnDownUnitySmeared", "");
    FillMET("pat", "patPFMetMVAElectronEnUpUnitySmeared", "");
    FillMET("pat", "patPFMetMVAJetEnDownUnitySmeared", "");
    FillMET("pat", "patPFMetMVAJetEnUpUnitySmeared", "");
    FillMET("pat", "patPFMetMVAJetResDownUnitySmeared", "");
    FillMET("pat", "patPFMetMVAJetResUpUnitySmeared", "");
    FillMET("pat", "patPFMetMVAMuonEnDownUnitySmeared", "");
    FillMET("pat", "patPFMetMVAMuonEnUpUnitySmeared", "");
    FillMET("pat", "patPFMetMVATauEnDownUnitySmeared", "");
    FillMET("pat", "patPFMetMVATauEnUpUnitySmeared", "");
    FillMET("pat", "patPFMetMVAUnclusteredEnDownUnitySmeared", "");
    FillMET("pat", "patPFMetMVAUnclusteredEnUpUnitySmeared", "");
    
    FillMET("pat", "patPFMetNoPileUpSmeared", "");
    FillMET("pat", "patPFMetNoPileUpElectronEnDownSmeared", "");
    FillMET("pat", "patPFMetNoPileUpElectronEnUpSmeared", "");
    FillMET("pat", "patPFMetNoPileUpJetEnDownSmeared", "");
    FillMET("pat", "patPFMetNoPileUpJetEnUpSmeared", "");
    FillMET("pat", "patPFMetNoPileUpJetResDownSmeared", "");
    FillMET("pat", "patPFMetNoPileUpJetResUpSmeared", "");
    FillMET("pat", "patPFMetNoPileUpMuonEnDownSmeared", "");
    FillMET("pat", "patPFMetNoPileUpMuonEnUpSmeared", "");
    FillMET("pat", "patPFMetNoPileUpTauEnDownSmeared", "");
    FillMET("pat", "patPFMetNoPileUpTauEnUpSmeared", "");
    FillMET("pat", "patPFMetNoPileUpUnclusteredEnDownSmeared", "");
    FillMET("pat", "patPFMetNoPileUpUnclusteredEnUpSmeared", "");
    FillMET("pat", "patPFMetMVAEMuSmeared", "");
    FillMET("pat", "patPFMetMVAElectronEnDownEMuSmeared", "");
    FillMET("pat", "patPFMetMVAElectronEnUpEMuSmeared", "");
    FillMET("pat", "patPFMetMVAJetEnDownEMuSmeared", "");
    FillMET("pat", "patPFMetMVAJetEnUpEMuSmeared", "");
    FillMET("pat", "patPFMetMVAJetResDownEMuSmeared", "");
    FillMET("pat", "patPFMetMVAJetResUpEMuSmeared", "");
    FillMET("pat", "patPFMetMVAMuonEnDownEMuSmeared", "");
    FillMET("pat", "patPFMetMVAMuonEnUpEMuSmeared", "");
    FillMET("pat", "patPFMetMVATauEnDownEMuSmeared", "");
    FillMET("pat", "patPFMetMVATauEnUpEMuSmeared", "");
    FillMET("pat", "patPFMetMVAUnclusteredEnDownEMuSmeared", "");
    FillMET("pat", "patPFMetMVAUnclusteredEnUpEMuSmeared", "");
    FillMET("pat", "patPFMetMVAEleSmeared", "");
    FillMET("pat", "patPFMetMVAElectronEnDownEleSmeared", "");
    FillMET("pat", "patPFMetMVAElectronEnUpEleSmeared", "");
    FillMET("pat", "patPFMetMVAJetEnDownEleSmeared", "");
    FillMET("pat", "patPFMetMVAJetEnUpEleSmeared", "");
    FillMET("pat", "patPFMetMVAJetResDownEleSmeared", "");
    FillMET("pat", "patPFMetMVAJetResUpEleSmeared", "");
    FillMET("pat", "patPFMetMVAMuonEnDownEleSmeared", "");
    FillMET("pat", "patPFMetMVAMuonEnUpEleSmeared", "");
    FillMET("pat", "patPFMetMVATauEnDownEleSmeared", "");
    FillMET("pat", "patPFMetMVATauEnUpEleSmeared", "");
    FillMET("pat", "patPFMetMVAUnclusteredEnDownEleSmeared", "");
    FillMET("pat", "patPFMetMVAUnclusteredEnUpEleSmeared", "");
    FillMET("pat", "patPFMetNoPileUpEMuSmeared", "");
    FillMET("pat", "patPFMetNoPileUpElectronEnDownEMuSmeared", "");
    FillMET("pat", "patPFMetNoPileUpElectronEnUpEMuSmeared", "");
    FillMET("pat", "patPFMetNoPileUpJetEnDownEMuSmeared", "");
    FillMET("pat", "patPFMetNoPileUpJetEnUpEMuSmeared", "");
    FillMET("pat", "patPFMetNoPileUpJetResDownEMuSmeared", "");
    FillMET("pat", "patPFMetNoPileUpJetResUpEMuSmeared", "");
    FillMET("pat", "patPFMetNoPileUpMuonEnDownEMuSmeared", "");
    FillMET("pat", "patPFMetNoPileUpMuonEnUpEMuSmeared", "");
    FillMET("pat", "patPFMetNoPileUpTauEnDownEMuSmeared", "");
    FillMET("pat", "patPFMetNoPileUpTauEnUpEMuSmeared", "");
    FillMET("pat", "patPFMetNoPileUpUnclusteredEnDownEMuSmeared", "");
    FillMET("pat", "patPFMetNoPileUpUnclusteredEnUpEMuSmeared", "");
    FillMET("pat", "patPFMetNoPileUpEleSmeared", "");
    FillMET("pat", "patPFMetNoPileUpElectronEnDownEleSmeared", "");
    FillMET("pat", "patPFMetNoPileUpElectronEnUpEleSmeared", "");
    FillMET("pat", "patPFMetNoPileUpJetEnDownEleSmeared", "");
    FillMET("pat", "patPFMetNoPileUpJetEnUpEleSmeared", "");
    FillMET("pat", "patPFMetNoPileUpJetResDownEleSmeared", "");
    FillMET("pat", "patPFMetNoPileUpJetResUpEleSmeared", "");
    FillMET("pat", "patPFMetNoPileUpMuonEnDownEleSmeared", "");
    FillMET("pat", "patPFMetNoPileUpMuonEnUpEleSmeared", "");
    FillMET("pat", "patPFMetNoPileUpTauEnDownEleSmeared", "");
    FillMET("pat", "patPFMetNoPileUpTauEnUpEleSmeared", "");
    FillMET("pat", "patPFMetNoPileUpUnclusteredEnDownEleSmeared", "");
    FillMET("pat", "patPFMetNoPileUpUnclusteredEnUpEleSmeared", "");
    
    //FIXME
    FillMET("pat", "patType1CorrectedPFMetEleSmeared", "");
    FillMET("pat", "patType1CorrectedPFMetEMuSmeared", "");
  }


      //MET smearing
  //  string n="ak5PFJetsL1FastL2L3";
//    const CandList& sjets = _e.jetList( EventManager::kPfJet, n );
//    for(size_t ij=0;ij<sjets.size();ij++) {
//      if(!isLepton(sjets[ij])) {
//        smearing( sjets[ij] );
//      }

//    }



  //  FillMET("decMET","core","",0,"mvaIdTight");
  //  FillMET("decMET","pu","",0,"mvaIdTight");
  
  //  missEDensity();


  //test "switch"
  // const Candidate* pfmet0 = _e.met( "pf" );
  // const Candidate* pfmet1 = _e.met( "pat" );
  // const Candidate* pfmet2 = _e.met( "pat", "patType1CorrectedPFMet" );
  // const Candidate* pfmet3 = _e.met( "pat", "patType1p2CorrectedPFMet" );

  // cout<<pfmet1->pt()<<"   "<<pfmet2->pt()<<"   "<<pfmet3->pt()<<"   "<<pfmet0->pt()<<endl;
  
  tm.flush();

  return true;
}


Candidate* METAnalysis::BuildBestZCand() {

  const CandList muons = _useMu?_e.muons():_e.electrons();

  CandList ZCands;

  float pt1=0, pt2=0;
  float pt_1=0, pt_2=0;
  int nZ=-1;

  for(int unsigned il1=0;il1<muons.size();il1++) {

    Candidate* l1 = muons[il1];
 
    if(_useMu?!MuonSel.accept( *l1 ):!ElectronSel.accept( *l1) ) continue;
 
    CandInfo* info1 = l1->info();
    if( info1->getFloat("fBrem") > 0.1 ) continue;

    for(int unsigned il2=il1+1;il2<muons.size();il2++) {
   
      Candidate* l2 = muons[il2];
    
      if(_useMu?!MuonSel.accept( *l2 ):!ElectronSel.accept( *l2)) continue;      
    
      CandInfo* info2 = l2->info();
      if( info2->getFloat("fBrem") > 0.1 ) continue;

      if(l1->vertex() != l2->vertex() ) continue; 

      Candidate* ZCand = Candidate::create(l1,l2);
      ZCands.push_back(ZCand);
      
      pt_1 = l1->pt()>l2->pt()?l1->pt():l2->pt();
      pt_2 = l1->pt()>l2->pt()?l2->pt():l1->pt();

      if(pt_1 > pt1 ) {
	pt1 = l1->pt();
	nZ = ZCands.size()-1;
      }
      if(pt_2 > pt2 && pt_1 == pt1) {
	pt2 = l2->pt();
	nZ = ZCands.size()-1;
      }
    
    }
  }
  
  if(nZ==-1)
    return NULL;

   return ZCands[nZ];
}



Candidate* METAnalysis::BuildGSFZCand() {

  const CandList gsfTrks = _e.gsfTracks();

  CandList ZCands;

  float pt1=0, pt2=0;
  float pt_1=0, pt_2=0;
  int nZ=-1;

  for(int unsigned it1=0;it1<gsfTrks.size();it1++) {

    Candidate* t1 = gsfTrks[it1];
    //cout<<t1->pt()<<"  "<<t1->eta()<<"  "<<t1->phi()<<"   "<<t1->vertex()<<endl;
    if(fabs(t1->eta() )>2.5 || t1->pt()<20) continue;
 
    for(int unsigned it2=it1+1;it2<gsfTrks.size();it2++) {
   
      Candidate* t2 = gsfTrks[it2];
    
      if(fabs(t2->eta() )>2.5 || t2->pt()<20) continue;

      if(t1->vertex() != t2->vertex() ) continue; 

      Candidate* ZCand = Candidate::create(t1,t2);
      ZCands.push_back(ZCand);
      
      pt_1 = t1->pt()>t2->pt()?t1->pt():t2->pt();
      pt_2 = t1->pt()>t2->pt()?t2->pt():t1->pt();

      if(pt_1 > pt1 ) {
	pt1 = t1->pt();
	nZ = ZCands.size()-1;
      }
      if(pt_2 > pt2 && pt_1 == pt1) {
	pt2 = t2->pt();
	nZ = ZCands.size()-1;
      }
    
    }
  }
 
  if(nZ==-1)
    return NULL;
  _zChan=1;
   return ZCands[nZ];
}

Candidate* METAnalysis::BuildBestZCandMultiChan() {

  _zChan=-1;

  const CandList muons = _e.muons();//:_e.electrons();

  CandList ZCands;

  float pt1=0, pt2=0;
  float pt_1=0, pt_2=0;
  int nZ=-1;
  int nE=0;
  int nM=0;


  for(int unsigned il1=0;il1<muons.size();il1++) {

    Candidate* l1 = muons[il1];
 
    if(!MuonSel.accept( *l1 ) ) continue;
    
    for(int unsigned il2=il1+1;il2<muons.size();il2++) {
   
      Candidate* l2 = muons[il2];
    
      if(!MuonSel.accept( *l2 )) continue;      
    
      if(l1->vertex() != l2->vertex() ) continue; 

      Candidate* ZCand = Candidate::create(l1,l2);
      ZCands.push_back(ZCand);
      
      pt_1 = l1->pt()>l2->pt()?l1->pt():l2->pt();
      pt_2 = l1->pt()>l2->pt()?l2->pt():l1->pt();

      if(pt_1 > pt1 ) {
	pt1 = l1->pt();
	nZ = ZCands.size()-1;
      }
      if(pt_2 > pt2 && pt_1 == pt1) {
	pt2 = l2->pt();
	nZ = ZCands.size()-1;
      }
    
    }
  }
  
  if(nZ!=-1) {
    _zChan=0;
    return ZCands[nZ];
  }

  //we did not have found any Z->mm, let's try a Z->ee

  const CandList elecs = _e.electrons();
  for(int unsigned il1=0;il1<elecs.size();il1++) {

    Candidate* l1 = elecs[il1];
 
    if(!ElectronSel.accept( *l1) ) continue;
 
    // CandInfo* info1 = l1->info();
    // if( info1->getFloat("fBrem") > 0.1 ) continue;

    nE++;

    for(int unsigned il2=il1+1;il2<elecs.size();il2++) {
   
      Candidate* l2 = elecs[il2];
    
      if(!ElectronSel.accept( *l2)) continue;      
    
      // CandInfo* info2 = l2->info();
      // if( info2->getFloat("fBrem") > 0.1 ) continue;

      nE++;

      if(l1->vertex() != l2->vertex() ) continue; 

      Candidate* ZCand = Candidate::create(l1,l2);
      ZCands.push_back(ZCand);
      
      pt_1 = l1->pt()>l2->pt()?l1->pt():l2->pt();
      pt_2 = l1->pt()>l2->pt()?l2->pt():l1->pt();

      if(pt_1 > pt1 ) {
	pt1 = l1->pt();
	nZ = ZCands.size()-1;
      }
      if(pt_2 > pt2 && pt_1 == pt1) {
	pt2 = l2->pt();
	nZ = ZCands.size()-1;
      }
    
    }
  }
  
  if(nZ!=-1) {
    _zChan=1;
    return ZCands[nZ];
  } 

  //and now... no Z, let's try the emu channel
   for(int unsigned il1=0;il1<muons.size();il1++) {

    Candidate* l1 = muons[il1];
 
    if(!MuonSel.accept( *l1 ) ) continue;
    for(int unsigned il2=0;il2<elecs.size();il2++) {

      Candidate* l2 = elecs[il2];
 
      if(!ElectronSel.accept( *l2) ) continue;
 
 
      if(l1->vertex() != l2->vertex() ) continue; 

      Candidate* ZCand = Candidate::create(l1,l2);
      ZCands.push_back(ZCand);
      
      pt_1 = l1->pt()>l2->pt()?l1->pt():l2->pt();
      pt_2 = l1->pt()>l2->pt()?l2->pt():l1->pt();

      if(pt_1 > pt1 ) {
	pt1 = l1->pt();
	nZ = ZCands.size()-1;
      }
      if(pt_2 > pt2 && pt_1 == pt1) {
	pt2 = l2->pt();
	nZ = ZCands.size()-1;
      }
    
    }
   }

   if(nZ!=-1) {
     _zChan=2;
     return ZCands[nZ];
   } 

   return NULL;

}


Candidate* METAnalysis::MVAJetMET(string type,float fact, string lvl) {
  
 TVector2 met(0,0);
 const CandList& jets = _e.jetList( EventManager::kPatJet);
 CandList uncPfc_ = _e.uncPfCandidates();

 for(unsigned int ijet=0;ijet<jets.size();ijet++) {

   const Candidate* jet = jets[ijet];
   const CandInfo* info = jet->info();

   if(info->getBool(lvl)|| isLepton(jet)) met -= jet->p2();
   else {
     if(type=="RedAll") met -= jet->p2()/fact;
   }
 }
 
 for(unsigned int ipfc=0;ipfc<uncPfc_.size();ipfc++) {
   
   const Candidate* pfc = uncPfc_[ipfc];

   if(mainVtx == pfc->vertex()  && pfc->charge()!=0 )  
     met -= pfc->p2();
   else {
     if(type=="RedAll" || type=="RedUnc") 
       met -= pfc->p2()/fact;
   }
 }

 Candidate* cmet = Candidate::create(met);
 return cmet;
}


Candidate* METAnalysis::AtlasMETAssoc(string type,  float fact) {

 TVector2 met(0,0);
 const CandList& jets = _e.jetList( EventManager::kPatJet);
 CandList uncPfc_ = _e.uncPfCandidates();

 for(unsigned int ijet=0;ijet<jets.size();ijet++) {

   const Candidate* jet = jets[ijet];
   //const CandInfo* info = jet->info();

   Vertex* vtxJet = JetVtxAssoc(jet, "atlas");
   //cout<<vtxJet<<"  "<<mainVtx<<"  "<< (vtxJet == mainVtx)<<endl;
   if(vtxJet == mainVtx|| isLepton(jet)) met -= jet->p2();
   else {
     if(type=="RedAll") met -= jet->p2()/fact;
   }
 }

 for(unsigned int ipfc=0;ipfc<uncPfc_.size();ipfc++) {
   
   const Candidate* pfc = uncPfc_[ipfc];
   //cout<<pfc->pt()<<"  "<<pfc->eta()<<endl;
   if(mainVtx == pfc->vertex()  && pfc->charge()!=0 )  
     met -= pfc->p2();
   else {
     if(type=="RedAll" || type=="RedUnc") 
       met -= pfc->p2()/fact;
   }
 }
 
 Candidate* cmet = Candidate::create(met);
 return cmet;
 
}


Candidate* METAnalysis::AtlasMET(float jpt, float fact) {

 TVector2 met(0,0);
 const CandList& jets = _e.jetList( EventManager::kPatJet);
 CandList uncPfc_ = _e.uncPfCandidates();

 // cout<<"======================== "<<jpt<<endl;

 for(unsigned int ijet=0;ijet<jets.size();ijet++) {

   const Candidate* jet = jets[ijet];
   // const CandInfo* info = jet->info();

   //cout<<jet->pt()<<endl;

   if(jet->pt()>jpt|| isLepton(jet)) met -= jet->p2();
   else met -= jet->p2()/fact;
  
 }

 for(unsigned int ipfc=0;ipfc<uncPfc_.size();ipfc++) {
   
   const Candidate* pfc = uncPfc_[ipfc];
   met -= pfc->p2()/fact;
  
 }
 
 //cout<<" mod "<<met.Mod()<<endl;

 Candidate* cmet = Candidate::create(met);
 //cout<<cmet->pt()<<endl;
 
 return cmet;
 
}


Candidate* METAnalysis::AtlasMETPt(string type, float jpt, float fact) {

 TVector2 met(0,0);
 const CandList& jets = _e.jetList( EventManager::kPatJet);
 CandList uncPfc_ = _e.uncPfCandidates();

 // cout<<"======================== "<<jpt<<endl;

 for(unsigned int ijet=0;ijet<jets.size();ijet++) {

   const Candidate* jet = jets[ijet];
   const CandInfo* info = jet->info();

   //cout<<jet->pt()<<endl;
    bool asJ=false;
    Vertex* vtxJet = JetVtxAssoc(jet, "atlas");
    if(type=="atlas" && vtxJet==mainVtx ) asJ=true;
    if(type=="mva" && info->getBool("mvaIdMedium")) asJ=true;
    if( isLepton(jet) ) asJ=true;

   if(jet->pt()>jpt && asJ) met -= jet->p2();
   else met -= jet->p2()/fact;
  
 }

 for(unsigned int ipfc=0;ipfc<uncPfc_.size();ipfc++) {
   
   const Candidate* pfc = uncPfc_[ipfc];

   if(pfc->charge()!=0 && mainVtx == pfc->vertex() )
     met -= pfc->p2();
   else
     met -= pfc->p2()/fact;
  
 }
 
 //cout<<" mod "<<met.Mod()<<endl;

 Candidate* cmet = Candidate::create(met);
 //cout<<cmet->pt()<<endl;
 
 return cmet;
 
}




Candidate* METAnalysis::AssocMET(string type,float fact) {

  TVector2 met(0,0);
  const CandList& jets = _e.jetList( EventManager::kPatJet);
  CandList uncPfc_ = _e.uncPfCandidates();

 for(unsigned int ijet=0;ijet<jets.size();ijet++) {

   const Candidate* jet = jets[ijet];
   // const CandInfo* info = jet->info();

   Vertex* vtxJet = JetVtxAssoc(jet, "basic");

   if(vtxJet == mainVtx|| isLepton(jet)) met -= jet->p2();
   else {
     if(type=="RedAll") met -= jet->p2()/fact;
   }
 }

 for(unsigned int ipfc=0;ipfc<uncPfc_.size();ipfc++) {
   
   const Candidate* pfc = uncPfc_[ipfc];
   
   if(mainVtx == pfc->vertex()  && pfc->charge()!=0 )  
     met -= pfc->p2();
   else {
     if(type=="RedAll" || type=="RedUnc") 
       met -= pfc->p2()/fact;
   }
 }
 
 Candidate* cmet = Candidate::create(met);
 return cmet;
 

}



// Jet-vertex association

Vertex* METAnalysis::JetVtxAssoc(const Candidate* jet, string opt) {

  Vertex* mVertex(0);

  map< Vertex* , float> sumPtAssoc;
  map< Vertex* , float>::const_iterator itmap;
  float sumPt;

  const CandMap& jetConst_ = _e.candMap("Jet-PfCand");	 

  CandMapIterator jc_lower_ = jetConst_.lower_bound( const_cast<Candidate*>(jet) );
  CandMapIterator jc_upper_ = jetConst_.upper_bound( const_cast<Candidate*>(jet) );
  CandMapIterator itjc;

  for( itjc=jc_lower_; itjc!=jc_upper_; itjc++ )
    { 
       Candidate* cand = itjc->second;
       if(cand->charge()==0) continue;

       sumPt += cand->pt();
       itmap = sumPtAssoc.find( cand->vertex() );
	
       if( itmap != sumPtAssoc.end() ) {
	  sumPtAssoc[ cand->vertex() ] += cand->pt();
       }
       else {
	 sumPtAssoc[ cand->vertex() ] = cand->pt();
       }
    }

  //now different options
  if(opt=="basic") {
    float sPttmp=0;
    for(itmap = sumPtAssoc.begin(); itmap != sumPtAssoc.end(); itmap++) {
      if(sPttmp < (*itmap).second ) {
	sPttmp = (*itmap).second;
	mVertex =  (*itmap).first;
      }
    }
  }
  if(opt=="atlas") {
    
    for(itmap = sumPtAssoc.begin(); itmap != sumPtAssoc.end(); itmap++) {
      if((*itmap).second/sumPt>0.75)
	{
	  mVertex = (*itmap).first;
	  break;
	}
    }
  }

//   cout<<" vertex for jet "<<jet->pt()<<"  "<<jet->eta()<<" with option "<<opt<<"  ---> ";
//   if(mVertex!=NULL) {
//     mVertex->print( cout );
//   }
//   else
//     cout<<" NUUUUUUUUUUUUUUlllll "<<endl;

  return mVertex;
}



void METAnalysis::FillMET(string METType, string redLvl, string factType, float pTVal, string lvlId) {

  float fact=1;
  bool bmet=false;

  Candidate* met;

  //cout<<METType<<"   "<<lvlId<<endl;

  if(factType=="vtx") {
    fact = _e.nGoodVertices();
  }
  if(factType=="atlas") {
    fact = AtlasFactor();
  }
  
  if(METType=="mva") {
    met = MVAJetMET( redLvl, fact, lvlId);
  }
  else if(METType=="atlasAssoc") {
    met = AtlasMETAssoc( redLvl, fact );
  }
  else if(METType=="atlas") {
    met = AtlasMET(pTVal, fact );
  }
  else if(METType=="atlasred") {
    met = AtlasMETPt(redLvl, pTVal, fact );
  }
  else if(METType=="basic") {
    met = AssocMET( redLvl, fact );
  }
  else if(METType=="clean") {
    met = const_cast<Candidate*>(_e.cleanPUmet(mainVtx, _ZCand, fact) );
  }
  else if(METType=="decMET") {
    float sumEt =0;
    if(redLvl=="core")
      met = coreMET(lvlId, sumEt);
    if(redLvl=="pu")
      met = puMET(lvlId, sumEt);

    tm.add<float>( redLvl+"_sumEt", &sumEt );

  }
  else {
    bmet=true;
    if(redLvl=="")
      met = const_cast<Candidate*>(_e.met( METType ) ); //ugly fix
    else
      met = const_cast<Candidate*>(_e.met( METType, redLvl) );

    const CandInfo* info = met->info();
    float sumEt = info->getFloat("sumEt");
    tm.add<float>( METType+"_"+redLvl+"__sumEt", &sumEt );
  }
  
  if(met==NULL) return;

  //Filling
  string name = METType+"_"+redLvl+"_"+factType;
  float pt = met->pt();
  float phi = met->phi();

  // cout<<redLvl<<" ->  "<<pt<<"   "<<phi<<endl;

  if(METType=="atlas") {
    ostringstream os; os << (int)pTVal;
    name = METType+"_"+os.str();
  }
  if(METType=="atlasred") {
    ostringstream os; os << (int)pTVal;
    name += "_"+os.str();
  }
  if(METType=="mva") name+="_"+lvlId;

  tm.add<float>( name+"_pt", &pt );
  tm.add<float>( name+"_phi", &phi );


  //Now smearing
//   if(!bmet || SampleName.find("Double")!=(size_t)-1 ) return;
//   Candidate* smMet = metJetSmearing( met );

  
//   float smPt = smMet->pt();
//   float smPhi = smMet->phi();

//   const CandInfo* info = smMet->info();
//   float smSumEt = info->getFloat("sumEt");
//   tm.add<float>( name+"Smear_pt", &smPt );
//   tm.add<float>( name+"Smear_phi", &smPhi );
//   tm.add<float>( name+"Smear_sumEt", &smSumEt );

}


float METAnalysis::AtlasFactor() {

  CandList tracks=_e.tracks();

  float sumPt=0;
  float sumPtMV=0;

  for(unsigned int it=0;it<tracks.size();it++) {

    const Candidate* trk= tracks[it];
    sumPt += trk->pt();
    
    if(mainVtx == trk->vertex())
      sumPtMV = trk->pt();
  }

  return sumPt/sumPtMV;

}


bool METAnalysis::isLepton(const Candidate* jet) {

  if( jet->dR( _ZCand->daughter(0) )<0.1 && jet->pt()/_ZCand->daughter(0)->pt() >0.8 ) return true;
  if( jet->dR( _ZCand->daughter(1) )<0.1 && jet->pt()/_ZCand->daughter(1)->pt() >0.8 ) return true;

  return false;
}



Candidate* METAnalysis::coreMET(string lvl, float& sumEt) {

  vectorFloat* chPts=new vectorFloat;
  vectorFloat* chEtas=new vectorFloat;
  vectorFloat* chPhis=new vectorFloat;
  vectorFloat* nePts=new vectorFloat;
  vectorFloat* neEtas=new vectorFloat;
  vectorFloat* nePhis=new vectorFloat;
  vectorFloat* phPts=new vectorFloat;
  vectorFloat* phEtas=new vectorFloat;
  vectorFloat* phPhis=new vectorFloat;

  vectorFloat* jetPts=new vectorFloat;
  vectorFloat* jetEtas=new vectorFloat;
  vectorFloat* jetPhis=new vectorFloat;

  TVector2 met(0,0);
  const CandList& jets = _e.jetList( EventManager::kPatJet);
  CandList uncPfc_ = _e.uncPfCandidates();
  
  for(unsigned int ijet=0;ijet<jets.size();ijet++) {
    
    const Candidate* jet = jets[ijet];
    const CandInfo* info = jet->info();
    
    bool isMatched=false;

    if(info->getBool(lvl)|| isLepton(jet)) {
      met -= jet->p2();
      sumEt += jet->pt();

      if(isLepton(jet)) isMatched=isMcMatched( jet );
      else isMatched = isJetMcMatched( jet );
   
      if(!isMatched) { //this jet is no matched with a genjet
	
	jetPts->push_back(jet->pt());
	jetEtas->push_back(jet->eta());
	jetPhis->push_back(jet->phi());
      }//matched
    }//jet id
  }//loop
  
  for(unsigned int ipfc=0;ipfc<uncPfc_.size();ipfc++) {
    
    const Candidate* pfc = uncPfc_[ipfc];
    
    if(mainVtx == pfc->vertex()  && pfc->charge()!=0 ) {
      met -= pfc->p2();
      sumEt += pfc->pt();
    
      bool isMatched = isMcMatched( pfc );
  
      if(!isMatched) { //particle does not exists

	if(pfc->pdgCode()==22) {
	  phPts->push_back( pfc->pt() );
	  phEtas->push_back( pfc->eta() );
	  phPhis->push_back( pfc->phi() );      

	}

	else if(pfc->charge()!=0 && fabs(pfc->eta())<2.45) {
	  chPts->push_back( pfc->pt() );
	  chEtas->push_back( pfc->eta() );
	  chPhis->push_back( pfc->phi() );      
	}
	  
	else {
	  nePts->push_back( pfc->pt() );
	  neEtas->push_back( pfc->eta() );
	  nePhis->push_back( pfc->phi() );      

	}
	
      }//matched
    }//id as core
    
  }//loop
    

  tm.add<vectorFloat*>("add_chPts", & chPts);
  tm.add<vectorFloat*>("add_chEtas", & chEtas);
  tm.add<vectorFloat*>("add_chPhis", & chPhis);
  tm.add<vectorFloat*>("add_nePts", & nePts);
  tm.add<vectorFloat*>("add_neEtas", & neEtas);
  tm.add<vectorFloat*>("add_nePhis", & nePhis);
  tm.add<vectorFloat*>("add_phPts", & phPts);
  tm.add<vectorFloat*>("add_phEtas", & phEtas);
  tm.add<vectorFloat*>("add_phPhis", & phPhis);
  
  tm.add<vectorFloat*>("add_jetPts", & jetPts);
  tm.add<vectorFloat*>("add_jetEtas", & jetEtas);
  tm.add<vectorFloat*>("add_jetPhis", & jetPhis);

  Candidate* cmet = Candidate::create(met);
  return cmet;
}



Candidate* METAnalysis::puMET(string lvl, float& sumEt) {

 TVector2 met(0,0);
 const CandList& jets = _e.jetList( EventManager::kPatJet);
 CandList uncPfc_ = _e.uncPfCandidates();

 for(unsigned int ijet=0;ijet<jets.size();ijet++) {

   const Candidate* jet = jets[ijet];
   const CandInfo* info = jet->info();

   if(!(info->getBool(lvl)|| isLepton(jet))) {
     met -= jet->p2();
     sumEt += jet->pt();
   }
}
 
 for(unsigned int ipfc=0;ipfc<uncPfc_.size();ipfc++) {
   
   const Candidate* pfc = uncPfc_[ipfc];

   if(!(mainVtx == pfc->vertex()  && pfc->charge()!=0) ) {
     met -= pfc->p2();
     sumEt += pfc->pt();
   }
 }

 Candidate* cmet = Candidate::create(met);
 return cmet;


}


bool METAnalysis::isJetMcMatched(const Candidate* jet) {

  const CandList& mcJets = _e.jetList( EventManager::kGenJet);

  //cout<<" New jet "<<endl;

  bool matched=false;
  for(size_t ij=0;ij<mcJets.size();ij++) {

    const Candidate* mcjet = mcJets[ij];    
    // cout<<jet->E()<<"  "<<jet->pt()<<"   "<<jet->eta()<<"   "<<jet->phi()<<endl;
    // cout<<"\t--->"<<mcjet->pt()<<"   "<<mcjet->eta()<<"   "<<mcjet->phi()<<"  "<<jet->dR( mcjet )<<endl;

    if( jet->dR( mcjet ) < 0.1 && mcjet->pt()/jet->pt() > 0.3 )
      {matched=true;// cout<<jet->pt()<<"  "<<mcjet->pt()<<endl;
	break;}
  }

  return matched;
  
}



bool METAnalysis::isMcMatched(const Candidate* pfc) {

  const CandList& mcPfcs = _e.mcTruthCandidates();

  //cout<<" New Part "<<endl;

  bool matched=false;
  for(size_t ij=0;ij<mcPfcs.size();ij++) {

    const Candidate* mcpfc = mcPfcs[ij];    

    //cout<<mcpfc->pt()<<"  "<<mcpfc->eta()<<"  "<<mcpfc->phi()<<"   "<<mcpfc->pdgCode()<<endl;
    //cout<<"\t ----> "<<pfc->pt()<<" <> "<<mcpfc->pt()<<"  "<< pfc->dR( mcpfc )<<endl;
    if( pfc->dR( mcpfc ) < 0.1 && mcpfc->pt()/pfc->pt() > 0.3 )
      {matched=true; break;
	//cout<<"\t ----> "<<pfc->pt()<<"  "<<pfc->pdgCode()<<"  <>  "<<mcpfc->pdgCode()<<endl;
      }
  }
  
  return matched;
}


void METAnalysis::missEDensity() {

  const CandList& mcPfcs = _e.mcTruthCandidates();
  CandList uncPfc_ = _e.uncPfCandidates();  
  const CandList& jets = _e.jetList( EventManager::kPatJet);
  //const CandList& genJets = _e.jetList( EventManager::kGenJet);

  vectorFloat* chPts=new vectorFloat;
  vectorFloat* chEtas=new vectorFloat;
  vectorFloat* chPhis=new vectorFloat;
  vectorFloat* nePts=new vectorFloat;
  vectorFloat* neEtas=new vectorFloat;
  vectorFloat* nePhis=new vectorFloat;
  vectorFloat* phPts=new vectorFloat;
  vectorFloat* phEtas=new vectorFloat;
  vectorFloat* phPhis=new vectorFloat;

  vectorFloat* chDrJ=new vectorFloat;
  vectorFloat* neDrJ=new vectorFloat;
  vectorFloat* phDrJ=new vectorFloat;

  for(size_t ip=0;ip<mcPfcs.size();ip++) {
    
    const Candidate* mcpfc = mcPfcs[ip];    
    
    int pdgId=mcpfc->pdgCode();

    if(abs(pdgId)<10 || abs(pdgId)==23 || abs(pdgId)==21 || abs(pdgId)==24 ) continue;

    bool notInJet=true;
    bool corePart=false;

    float tmp=1000;

    for(unsigned int ijet=0;ijet<jets.size();ijet++) {
      
      const Candidate* jet = jets[ijet];
      const CandInfo* info = jet->info();
      
      if(info->getBool("mvaIdTight")|| isLepton(jet)) {
	
	if(!isLepton(jet) )
	  if(mcpfc->dR(jet)<0.45) { 
	    //particle in a core jet
	    notInJet=false; break;
	  }
	  else {
	    if(tmp>mcpfc->dR(jet)) tmp=mcpfc->dR(jet);
	  }
      }
    }//jets


    if(!notInJet) continue; //we did not missed the particle
    
    for(unsigned int ipfc=0;ipfc<uncPfc_.size();ipfc++) {
      
      const Candidate* pfc = uncPfc_[ipfc];
      
      if(mainVtx == pfc->vertex()  && pfc->charge()!=0 )
	{
	  if(mcpfc->dR( pfc )< 0.05 && mcpfc->pt()/pfc->pt()<0.3 ) {
	    //particle is a core particle
	    corePart=true; break;
	  }
	  
	}
    }//uncpfcs

    if(corePart) continue;


    //now we have only missing particles

    if(abs(pdgId)==22 || abs(pdgId)==111) {

      phPts->push_back( mcpfc->pt() );
      phEtas->push_back( mcpfc->eta() );
      phPhis->push_back( mcpfc->phi() );      
      phDrJ->push_back( tmp );
    }
    
    else if(abs(mcpfc->charge())!=0 && fabs(mcpfc->eta())<2.5) { //charged
      
      chPts->push_back( mcpfc->pt() );
      chEtas->push_back( mcpfc->eta() );
      chPhis->push_back( mcpfc->phi() );
      chDrJ->push_back( tmp );  
    }
    
    else {

      nePts->push_back( mcpfc->pt() );
      neEtas->push_back( mcpfc->eta() );
      nePhis->push_back( mcpfc->phi() );
      neDrJ->push_back( tmp );
    }

  }//mc pfcs



  tm.add<vectorFloat*>("miss_chPts", & chPts);
  tm.add<vectorFloat*>("miss_chEtas", & chEtas);
  tm.add<vectorFloat*>("miss_chPhis", & chPhis);
  tm.add<vectorFloat*>("miss_nePts", & nePts);
  tm.add<vectorFloat*>("miss_neEtas", & neEtas);
  tm.add<vectorFloat*>("miss_nePhis", & nePhis);
  tm.add<vectorFloat*>("miss_phPts", & phPts);
  tm.add<vectorFloat*>("miss_phEtas", & phEtas);
  tm.add<vectorFloat*>("miss_phPhis", & phPhis);
  
  tm.add<vectorFloat*>("miss_chDrJ", & chDrJ);
  tm.add<vectorFloat*>("miss_neDrJ", & neDrJ);
  tm.add<vectorFloat*>("miss_phDrJ", & phDrJ);

  
}







////////////
//////////// Jet Smearing =====================================
////////////



void 
METAnalysis::LoadJetSmearingDB() {

  string db="/afs/cern.ch/user/m/mmarionn/workspace/private/AnaNaS/workdir/database/pfJetResolutionMCtoDataCorrLUT.root";
  TFile* file=new TFile(db.c_str(), "READ");

  lut_ = (TH2D*)file->Get("pfJetResolutionMCtoDataCorrLUT");
}


Candidate* 
METAnalysis::oldToNewsmearing(const Candidate* smJet) {

  // jet : corrected or raw jet entered by the user
  // rawjet : matched raw jet...

  bool debug_=true;

  const Candidate* jet=corJetMatch(smJet);
  const Candidate* rawJet=rawJetMatch(smJet);
  const Candidate* genJet=genJetMatch(smJet);
  //const Candidate* smJet=smearedJetMatch(jet);
  
  if(debug_) {
    cout<<endl<<" ========================================================================="<<endl;
    cout<< " jet "<<jet->pt()<<"  "<<jet->eta()<<"  "<<jet->phi()<<"  "<<jet->E()<<endl;
    cout<< " raw "<<rawJet->pt()<<"  "<<rawJet->eta()<<"  "<<rawJet->phi()<<"  "<<rawJet->E()<<endl;
    cout<< " smj "<<smJet->pt()<<"  "<<smJet->eta()<<"  "<<smJet->phi()<<"  "<<smJet->E()<<endl;
    if(genJet!=NULL)
      cout<< " gen "<<genJet->pt()<<"  "<<genJet->eta()<<"  "<<genJet->phi()<<"  "<<genJet->E()<<endl;
  }

  double smearFactor = 1.;   
  double newSF = 1.;   
  double x = TMath::Abs(jet->eta());
  double y = jet->pt();

  if(debug_)
    cout<<" lut_ exists? "<<lut_<<endl;

  if ( x > lut_->GetXaxis()->GetXmin() && x < lut_->GetXaxis()->GetXmax() &&
       y > lut_->GetYaxis()->GetXmin() && y < lut_->GetYaxis()->GetXmax() ) {
    int binIndex = lut_->FindBin(x, y);
	
    if(debug_) {
      cout<<" ---> "<<x<<"  "<<y<<"  "<<binIndex<<endl;
    }

      if ( smearBy_ > 0. ) {
	newSF += smearBy_*(lut_->GetBinContent(binIndex) - 1.);
	smearFactor += 1.0*(lut_->GetBinContent(binIndex) - 1.);
      
	if(debug_) {
	  cout<<" smear factor : "<<smearBy_<<"  "<<newSF<<"  "<<smearFactor<<endl;
	}
	
      }
    double smearFactorErr = lut_->GetBinError(binIndex);
   
    if ( shiftBy_ != 0. ) {
      smearFactor += (shiftBy_*smearFactorErr);
      newSF+= (shiftBy_*smearFactorErr);

      if(debug_) {
	cout<<" shift factor : "<<shiftBy_<<"  "<<newSF<<"  "<<smearFactor<<endl;
      }
    }
  }


  // New Smearing =============================================================== 

  float smearedjetEn = smJet->E()/jet->E();
  float newSmearEn = jet->E();
  float Max=max(rawJet->E(),jet->E() );


  if(debug_) {
    cout<<" old smeared energy : "<<smearedjetEn<<" max: "<<Max<<endl;
  }


  if(genJet==NULL) {

    if( smearFactor> 1.) {
    
      float gausSmear = (smearedjetEn/rawJet->E() - 1)*Max;
    
      float fPol = smearFactor*smearFactor -2*smearFactor + 1;
      float fPolp = newSF*newSF -2*newSF + 1;
      float newGausS = exp(log(fabs(gausSmear))*(fPol/fPolp));
    
      newSmearEn = rawJet->E()*(1+ (gausSmear>0)?newGausS:(newGausS*-1 ) )/Max;

      if(debug_) {
	cout<<" remove old smearing : gS="<<gausSmear<<" fPol="<<fPol<<" fPolp="<<fPolp<<" ngS="<<newGausS<<endl;
	cout<<" ---> new smearing "<<newSmearEn<<endl;
      }
    }
  
  }
  else {
 
    if(smJet->E()!=jet->E()) { //no smearing because too far from genJet

      double dEn = jet->E() - genJet->E();
      
      newSmearEn = rawJet->E()*(1. + (newSF - 1.)*dEn/TMath::Max(rawJet->E(), jet->E()));
      
      if(debug_) {
	cout<<" New Smearing :  dEn="<<dEn<<" newSmearing="<<newSmearEn<<endl;
      }
    }

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
  if ( !( abs(jet->E() - rawJet->E() ) > (5.*min(jet->E(), rawJet->E() )) 
	   //(skipJetSelection_ && (*skipJetSelection_)(jet)) || 
	|| rawJet->pt()  < skipRawJetPtThreshold_          ||
	 jet->pt() < skipCorrJetPtThreshold_         ) ) {
    smearedJetP4 *= (newSmearEn/jet->E());

    if(debug_) {
      cout<<" Default jet     : "<<jet->pt()<<"  "<<jet->eta()<<"  "<<jet->phi()<<"  "<<jet->E()<<endl;
      cout<<" New Smeared Jet : "<<smearedJetP4.Pt()<<"  "<<smearedJetP4.Eta()<<"  "<<smearedJetP4.Phi()<<"  "<<smearedJetP4.E()<<endl;
      cout<<" Old Smeared Jet : "<<smJet->pt()<<"  "<<smJet->eta()<<"  "<<smJet->phi()<<"  "<<smJet->E()<<endl;
      
    } 
  }

      
  Candidate* smearJet = Candidate::create( smearedJetP4.Vect(), jet->charge(), 
					   smearedJetP4.M(), jet->vertex() );
      
  return smearJet;

}



Candidate* 
METAnalysis::smearing(const Candidate* jet) {

  // jet : corrected or raw jet entered by the user
  // rawjet : matched raw jet...

  CandInfo* jetInfo=jet->info();

  bool debug_=false;

  if(debug_) {
    cout<<"        =================== New Jet ==================== "<<endl;
    cout<<"initial jet : "<<jet->pt()<<"   "<<jet->eta()<<"   "<<jet->phi()<<endl;
    }

  //const Candidate* jet=corJetMatch(smJet);
  const Candidate* rawJet=rawJetMatch(jet);
  const Candidate* genJet=genJetMatch(jet);
  //const Candidate* smJet=smearedJetMatch(jet);

  double smearFactor = 1.;      
  double x = TMath::Abs(jet->eta());
  double y = jet->pt();// cout<<" ---------------- > "<<lut_<<endl;
  if ( x > lut_->GetXaxis()->GetXmin() && x < lut_->GetXaxis()->GetXmax() &&
       y > lut_->GetYaxis()->GetXmin() && y < lut_->GetYaxis()->GetXmax() ) {
    int binIndex = lut_->FindBin(x, y);
	
    if ( smearBy_ > 0. ) smearFactor += smearBy_*(lut_->GetBinContent(binIndex) - 1.);
    double smearFactorErr = lut_->GetBinError(binIndex);
    if ( debug_ ) std::cout << "smearFactor = " << smearFactor << " +/- " << smearFactorErr << std::endl;
	
    if ( shiftBy_ != 0. ) {
      smearFactor += (shiftBy_*smearFactorErr);
      if ( debug_ ) std::cout << "smearFactor(shifted) = " << smearFactor << std::endl;
    }
  }

  double smearedJetEn = jet->E();
  double sigmaEn = jetInfo->getFloat("jetResolution")*TMath::Sqrt(smearFactor*smearFactor - 1.);
  bool isGenMatched = false;
  if ( genJet ) {
    if ( debug_ ) {
      std::cout << "genJet: Pt = " << genJet->pt() << ", eta = " << genJet->eta() << ", phi = " << genJet->phi() << std::endl;
    }
    double dEn = jet->E() - genJet->E();
    if ( dEn < (3.0*sigmaEn) ) {
      //--- case 1: reconstructed jet matched to generator level jet, 
      //            smear difference between reconstructed and "true" jet energy

      if ( debug_ ) {
	std::cout << " successfully matched to genJet" << std::endl;	
	std::cout << "corrJetEn = " << jet->E() << ", genJetEn = " << genJet->E() << " --> dEn = " << dEn << std::endl;
      }

      smearedJetEn = jet->E()*(1. + (smearFactor - 1.)*dEn/TMath::Max(rawJet->E(), jet->E()));
      isGenMatched = true;
    }
  }
  if ( !isGenMatched ) {
    //--- case 2: reconstructed jet **not** matched to generator level jet, 
    //            smear jet energy using MC resolution functions implemented in PFMEt significance algorithm (CMS AN-10/400)

    if ( debug_ ) {
      std::cout << " not matched to genJet" << std::endl;
      std::cout << "corrJetEn = " << jet->E() << ", sigmaEn = " << sigmaEn << std::endl;
    }

    if ( smearFactor > 1. ) {
      // CV: MC resolution already accounted for in reconstructed jet,
      //     add additional Gaussian smearing of width = sqrt(smearFactor^2 - 1) 
      //     to account for Data/MC **difference** in jet resolutions.
      //     Take maximum(rawJetEn, corrJetEn) to avoid pathological cases
      //    (e.g. corrJetEn << rawJetEn, due to L1Fastjet corrections)

      smearedJetEn = jet->E()*(1. + rnd_.Gaus(0., sigmaEn)/TMath::Max(rawJet->E(), jet->E()));
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
  if ( !( abs(jet->E() - rawJet->E() ) > (5.*min(jet->E(), rawJet->E() ) ) || 
	  //(skipJetSelection_ && (*skipJetSelection_)(jet)) ||
	  rawJet->pt()  < skipRawJetPtThreshold_          ||
	  jet->pt() < skipCorrJetPtThreshold_         ) ) {
    if ( debug_ ) {
      std::cout << " smearing jetP4 by factor = " << (smearedJetEn/jet->E()) << " --> smearedJetEn = " << smearedJetEn << std::endl;
    }
    smearedJetP4 *= (smearedJetEn/jet->E());
  }
      
  if ( debug_ ) {
    std::cout << "smearedJet: Pt = " << smearedJetP4.Pt() << ", eta = " << smearedJetP4.Eta() << ", phi = " << smearedJetP4.Phi() << std::endl;
    std::cout << " dPt = " << (smearedJetP4.Pt() - jet->pt()) 
	      << " (Px = " << (smearedJetP4.Px() - jet->px()) << ", Py = " << (smearedJetP4.Py() - jet->py()) << ")" << std::endl;
  }
    
      
  Candidate* smearJet = Candidate::create( smearedJetP4.Vect(), jet->charge(), 
					   smearedJetP4.M(), jet->vertex() );
      
  return smearJet;
 
}



const Candidate* 
METAnalysis::genJetMatch(const Candidate* jet) {

  const CandList& mcJets = _e.jetList( EventManager::kGenJet);
  
  //  bool matched=false;
  for(size_t ij=0;ij<mcJets.size();ij++) {
    
    const Candidate* mcjet = mcJets[ij];    
  //   if(mcjet->pt()> 20)
//       cout<<jet->dR( mcjet )<<"   "<<min(0.5, 0.1 + 0.3*exp(-0.05*(mcjet->pt() - 10.)))<<endl;
    if( jet->dR( mcjet ) < min(0.5, 0.1 + 0.3*exp(-0.05*(mcjet->pt() - 10.))) )
      return mcjet;
  }

  return NULL;

}


const Candidate* 
METAnalysis::rawJetMatch(const Candidate* jet) {

  string s= "ak5PFJets";
  const CandList& rawJets = _e.jetList( EventManager::kPfJet, s ); //default rawPfJets
  
  //  bool matched=false;
  for(size_t ij=0;ij<rawJets.size();ij++) {
    
    const Candidate* rawjet = rawJets[ij];    
    //cout<<" raw ---> "<<jet->dR( rawjet )<<"  "<<jet->pt()<<"   "<<rawjet->pt()<<endl; 
    if( jet->dR( rawjet ) < 0.02 )
      return rawjet;
  }

  return jet;

}

const Candidate* 
METAnalysis::smearedJetMatch(const Candidate* jet) {

  string s= "smearedAK5PFJets";
  const CandList& rawJets = _e.jetList( EventManager::kPfJet, s ); //default rawPfJets
  
  //  bool matched=false;
  for(size_t ij=0;ij<rawJets.size();ij++) {
    
    const Candidate* rawjet = rawJets[ij];    
    if( jet->dR( rawjet ) < 0.02 )
      return rawjet;
  }

  return jet;
}


const Candidate* 
METAnalysis::corJetMatch(const Candidate* jet) {

  string s= "ak5PFJetsL1FastL2L3";
  const CandList& corJets = _e.jetList( EventManager::kPfJet, s ); //default rawPfJets
  
  //  bool matched=false;
  for(size_t ij=0;ij<corJets.size();ij++) {
    
    const Candidate* corjet = corJets[ij];    
    if( jet->dR( corjet ) < 0.02 )
      return corjet;
  }

  return jet;
}






Candidate*
METAnalysis::metJetSmearing(Candidate* met) {

  TVector2 smMet= met->p2();

  CandInfo* info = met->info();
  float sumEt = info->getFloat("sumEt");

  string n="ak5PFJetsL1FastL2L3";
  const CandList& sjets = _e.jetList( EventManager::kPfJet, n );
  for(size_t ij=0;ij<sjets.size();ij++) {
    if(!isLepton(sjets[ij])) {
      Candidate* smearJet =  smearing( sjets[ij] );

      smMet += smearJet->p2() - sjets[ij]->p2();
      sumEt += smearJet->pt() - sjets[ij]->pt();
    }
    
  }

  Candidate* cmet = Candidate::create(smMet);

  CandInfo* cInfo = CandInfo::create();
  cmet->setInfo( cInfo  );
  
  cInfo->setFloat("sumEt", sumEt );
  
  return cmet;
  
}



bool 
METAnalysis::HLTFired() {

  if( _useMu? _e.isFired("HLT_Mu17_Mu8",1):
      _e.isFired("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL",1) )
    return true;
  else  
    return false;
}

bool 
METAnalysis::HLTFiredMultiChan() {

  if( (_zChan==0 && _e.isFired("HLT_Mu17_Mu8",1)) || 
      (_zChan==1 && _e.isFired("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL",1) ) ||
      (_zChan==2 && _e.isFired("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL",1)) )
    return true;
  else  
    return false;
}



// bool
// METAnalysis::ElectronSel(const Candidate* c) {


//   const CandInfo* info = c->info();

//   float pt= info->getFloat("caloEnergy");
//   float eta= info->getFloat("caloEta");

//   if( pt < 20){// cout<<" you1 "<<pt<<endl;
//     return false;}

//   if(fabs( eta )>2.45 || 
//      (fabs( eta )>1.4444 &&
//       fabs( eta )<1.56) ) {//cout<<" you2 "<<endl; 
//     return false;}

//   if(!info->getBool("mediumID")) return false;

//   return true;


//   bool isEB=true;
//   if( fabs(eta)>1.5) isEB=false;


//   if( fabs(info->getFloat("dEtaIn"))>(isEB?0.004:0.005)) {//cout<<" you3 "<<info->getFloat("dEtaIn")<<endl; 
//     return false;}
//   if( fabs(info->getFloat("dPhiIn"))>(isEB?0.03:0.02)) {//cout<<" you4 "<<info->getFloat("dPhiIn")<<endl;
//     return false;}
//   if( info->getFloat("HOE")>(isEB?0.12:0.10 )) {//cout<<" you5 "<<endl;
//     return false;}
//   if( info->getFloat("sigiEtaiEta")>(isEB?0.01:0.03) ) {//cout<<" you6 "<<endl;
//     return false;}
  
//   if( c->d0( c->vertex() )>0.02 ) {//cout<<" you7 "<<endl;
//     return false;}

//   if( info->getBool("isConv") ) {//cout<<" you8 "<<endl;
//     return false;}
  
//   float iso = (info->getFloat("dr03TrkIso") + 
// 	       info->getFloat("dr03EcalIso") +
// 	       info->getFloat("dr03HcalD1Iso") +
// 	       info->getFloat("dr03HcalD2Iso") -
// 	       0.10*_e.getRhoFastJet())/pt;

//   if( iso > (isEB?0.15:0.12) ) {//cout<<" you9 ->"<<iso<<endl;
//     return false;}

//   return true;

//}

Candidate* METAnalysis::pfMatch(Candidate* el) {


  const CandList pfcands = _e.pfCandidates();


  for(size_t ip=0; ip<pfcands.size();ip++) {

    Candidate* pfc=pfcands[ip];

    if(pfc->pt() < 5) continue;

    if( el->dR( pfc ) <0.1 && fabs(pfc->pdgCode())==13)
      return pfc;
  }

  return 0;
}


TLorentzVector METAnalysis::pfConeMatch(Candidate* el) {


  const CandList pfcands = _e.pfCandidates();

  TLorentzVector Zpart(0.00001,0,0,0);

  for(size_t ip=0; ip<pfcands.size();ip++) {

    Candidate* pfc=pfcands[ip];

    if( el->dR( pfc ) <0.2 )
      //	&& (fabs(pfc->pdgCode())==13
      //    || (fabs(pfc->pdgCode())==22 && pfc->pt()>0.5) ))
      Zpart += pfc->p4();
  }

 //  _e.refreshPfCandidates();
//   int v=1; //pfCand_lowPt
//   _e.reloadPfCandidates( v ); //protection
//   const CandList pfcandslow = _e.pfCandidates();
  
//   for(size_t ip=0; ip<pfcandslow.size();ip++) {
    
//     Candidate* pfc=pfcandslow[ip];
    
//     if( el->dR( pfc ) <0.1)
//       Zpart += pfc->p2();
//   }
  

  // if(el->pt() > Zpart.Pt() ) {
  //   TLorentzVector Z(0,0,0,0);
  //   cout<<" main candidate : "<<el->pt()<<"  ==> "<<Zpart.Pt()<<endl;
  //   cout<<" detail : "<<endl;
  //   for(size_t ip=0; ip<pfcands.size();ip++) {

  //     Candidate* pfc=pfcands[ip];

  //     if( el->dR( pfc ) <0.1 ) {
  // 	Z += pfc->p4();
  // 	cout<<" -> "<<pfc->pdgCode()<<" :: "<<pfc->pt()<<"   "<<Z.Pt()<<endl;
  //     }
  //   }
  // }

  return Zpart;
  
}


Candidate* 
METAnalysis::lepTrkMap(Candidate* lep) {

  Candidate* trk(0);

  if(abs(lep->pdgCode())==11) {
    trk = _e.getSecond( "electron-GsfTrack", lep );
  }
  else {
    trk = _e.getSecond( "muon-Track", lep );
  }
  
  return trk;
}
