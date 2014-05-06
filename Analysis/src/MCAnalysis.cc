#include "Analysis/src/MCAnalysis.hh"

#include <algorithm>
#include <cassert>
#include <sstream>

using namespace std;

#include <TH1I.h>
#include <TH1F.h>
#include <TH2I.h>
#include <TH2F.h>

#include "Analysis/utils/Constants.hh"
#include "Analysis/utils/Config.hh"
#include "Analysis/core/CandInfo.hh"
#include "Analysis/core/EventManager.hh"
#include "Analysis/selectors/IsolationSelector.hh"
#include "Analysis/selectors/OfficialIsolationSelector.hh"
#include "Analysis/selectors/ConeSelector.hh"
#include "Analysis/utils/StatUtils.hh"
#include "Analysis/utils/CutUtils.hh"
#include "Analysis/tools/PrintTree.hh"
#include "Analysis/ghmzz/Weights.hh"


ClassImp( MCAnalysis )
  
  MCAnalysis::MCAnalysis( Sample& sample, const std::string & collectionFileName ):
  SampleAnalysis( "MCAnalysis" , sample, collectionFileName )
{
  cout << "\t-------------------------------------------" << endl; 
  cout << "\t---            MC  analysis            ----" << endl; 
  cout << "\t-------------------------------------------" << endl; 

  _verbose = true;

  cout << "\nK-factors for ZZ" << endl;
  string kfDB_ZZ_jv30 = Config::confPath + "kFactorDatabaseZZ_jv30.txt";
  _kFactors["ZZ_jv30"] = new KFactor(kfDB_ZZ_jv30);
  string kfDB_ZZ = Config::confPath + "kFactorDatabaseZZ.txt";
  _kFactors["ZZ"] = new KFactor(kfDB_ZZ);

  cout << "\nK-factors for WZ" << endl;
  string kfDB_WZ_jv30 = Config::confPath + "kFactorDatabaseWZ_jv30.txt";
  _kFactors["WZ_jv30"] = new KFactor(kfDB_WZ_jv30);
  string kfDB_WZ = Config::confPath + "kFactorDatabaseWZ.txt";
  _kFactors["WZ"] = new KFactor(kfDB_WZ);
  
  cout << "\nK-factors for WW" << endl;
  string kfDB_WW_jv30 = Config::confPath + "kFactorDatabaseWW_jv30.txt";
  _kFactors["WW_jv30"] = new KFactor(kfDB_WW_jv30);
  string kfDB_WW = Config::confPath + "kFactorDatabaseWW.txt";
  _kFactors["WW"] = new KFactor(kfDB_WW);

  _N=0;

}

MCAnalysis::~MCAnalysis()
{
}

bool
MCAnalysis::analyzeEvent()
{
  if( _N<100 ) _verbose = true;
  else         _verbose = false;

  float  mVV(0);

  float  mV[2] = {0,0};
  float ptV[2] = {0,0};
  float pzV[2] = {0,0};
  float  EV[2] = {0,0};
  float  yV[2] = {0,0};
  
  // acceptance
  bool leptonEta   = true;
  bool leptonPt    = true;

  vector<float> ptL(4,0);
  vector<float> etaL(4,0);
  vector<float> phiL(4,0);
  vector<float> EL(4,0);
  vector<int> pdgL(4,0);

  Candidate* genCand_ = _e.decayTree();
  assert( genCand_!=0 );
  string categ_ = _e.categ();

  bool isZZ =  categ_.find("ZZ")!=string::npos;
  bool isWZ =  categ_.find("WZ")!=string::npos;

  KFactor* vvkf;
  KFactor* vvkf_jv30;
  //  if( _e.categ() == "ZZ_2e2n" || _e.categ() == "ZZ_2m2n" )     
  if(      isZZ )  
    {
      vvkf      = _kFactors["ZZ"];
      vvkf_jv30 = _kFactors["ZZ_jv30"];
    }
  else if( isWZ )  
    {
      vvkf      = _kFactors["WZ"];
      vvkf_jv30 = _kFactors["WZ_jv30"];
    }
	
  // get the two Z bosons
  CandList Z_;
  CandUtil::get(  23, genCand_, Z_ );

  // get the W bosons
  CandList W_;
  CandUtil::get(  24, genCand_, W_ );
  CandUtil::get( -24, genCand_, W_ );

  // sanity
  assert( ( isZZ && Z_.size()==2 ) || 
	  ( isWZ && Z_.size()==1 && W_.size()==1 ) );
    
  Candidate* V_[2] = {0,0};
  if( isZZ )
    {
      // search for neutrinos (12:nu_e; 14:nu_mu; 16:nu_tau)
      int iV2 = 2;
      for( int iZ=0; iZ<2; iZ++ )
	{
	  if( Z_[iZ]->nDaughters()!=2 ) continue;
	  int pdgCode_ =  Z_[iZ]->daughter(0)->pdgCode();
	  if(    pdgCode_ == 12 || pdgCode_ == -12 ||
		 pdgCode_ == 14 || pdgCode_ == -14 ||
		 pdgCode_ == 16 || pdgCode_ == -16   )
	    {
	      iV2 = iZ;
	      break;
	    }
	}
      assert( iV2<2 );
      int iV1 = 1-iV2;
      V_[0] = Z_[iV1];
      V_[1] = Z_[iV2];
    }
  else if( isWZ )
    {
      V_[0] = Z_[0];
      V_[1] = W_[0];
    }

  Candidate* VV_ = Candidate::create( V_[0], V_[1] );
  
  mVV  = VV_->mass();
  
  for( size_t ii=0; ii<2 ; ii++ )
    {
      mV[ii]  = V_[ii]->mass();
      ptV[ii] = V_[ii]->pt();
      pzV[ii] = V_[ii]->pz();
      EV[ii]  = V_[ii]->E();
      yV[ii]  = KineUtils::y( EV[ii], pzV[ii] );
    }

  for(int il=0;il<4;il++) 
    {
      int iV = ( il/2 == 0 ) ? 0:1;
      ptL[  il ]  = ( V_[ iV ]->daughter(il%2)->pt() );
      etaL[ il ]  = ( V_[ iV ]->daughter(il%2)->eta() );
      phiL[ il ]  = ( V_[ iV ]->daughter(il%2)->phi() );
      pdgL[ il ]  = ( V_[ iV ]->daughter(il%2)->pdgCode() );
      EL[   il ]  = ( V_[ iV ]->daughter(il%2)->E() );
    }      
  
  bool isEEorMM = true;
  for( int il=0; il<2; il++ )
    {
      if( fabs(ptL[il])<20 ) 
	{
	  leptonPt = false;
	  break;
	}
      if( abs(pdgL[il]) == 11 )
	{
	  if( fabs(etaL[il])>2.5 ) 
	    {
	      leptonEta = false;
	      break;
	    }	      
	}	 
      else if( abs(pdgL[il]) == 13 )
	{
	  if( fabs(etaL[il])>2.4 ) 
	    {
	      leptonEta = false;
	      break;
	    }	      
	}
      else
	{
	  leptonEta = false;
	  isEEorMM = false;
	}
    }
  
  float kfactor       = (*vvkf)( ptV[0] );	
  float kfactor_jv30  = (*vvkf_jv30)( ptV[0] );	
  
  fill( "mll", "mZ", mV[0], 1. );
  fill( "mll", "mZ_KF", mV[0], kfactor );
    
  bool dileptonM   = mV[0]>=60 && mV[0]<=120;
  bool inAcceptance  = dileptonM && leptonEta  && leptonPt;
  if( isZZ && !dileptonM ) return false;
  if( isWZ ) 
    {
      if( !isEEorMM ) return false;
      if( !dileptonM ) return false;
    }
  fill( "mll", "mZ_60-120", mV[0], 1. );
  fill( "mll", "mZ_60-120_KF", mV[0], kfactor );

  if( inAcceptance )
    {
      fill( "mll", "mZ_acc", mV[0], 1. );
      fill( "mll", "mZ_acc_KF", mV[0], kfactor );
    }
      
  // if( inAcceptance )
  //	{
  //	  cout << "*****************" << endl;
  //	  cout << "MC Analysis -- Event " << _ievt << " -- categ " << _e.categ() << "  k=" << kfactor << endl;
  //	}
  _counters[categ_+"__N"] += 1;
  _counters[categ_+"__N_KF"] += kfactor;
  if( leptonEta  )     _counters[categ_+"__N_KF_lepEta"] += kfactor;
  if( leptonPt   )     _counters[categ_+"__N_KF_lepPt" ] += kfactor;
  if( inAcceptance )   _counters[categ_+"__N_KF_Acc"   ] += kfactor;

  buildLeptonLists();
  
  size_t nEleVL  = _e.userCandList("ElectronVeryLoose").size();
  size_t nEleL   = _e.userCandList("ElectronLoose").size();
  size_t nEleT   = _e.userCandList("ElectronTight").size();
  size_t nEleVT  = _e.userCandList("ElectronVeryTight").size();

  size_t nMuoVL  = _e.userCandList("MuonVeryLoose").size();
  size_t nMuoL   = _e.userCandList("MuonLoose").size();
  size_t nMuoT   = _e.userCandList("MuonTight").size();
  size_t nMuoVT  = _e.userCandList("MuonVeryTight").size();

  size_t nLepL(0);
  
// 	      const Candidate* mcEl1  = _e.mcCand( *el1->theBase() );	
// 	      const Candidate* mcMot1(0);
// 	      if( mcEl1 ) mcMot1 = mcEl1->theMother();
// 	      const Candidate* mcEl2  = _e.mcCand( *el2->theBase() );
// 	      const Candidate* mcMot2(0);
// 	      if( mcEl2 ) mcMot2 = mcEl2->theMother();

// 	      bool goodMatch = false;
// 	      if( mcMot1!=0 && mcMot1==mcMot2 )
// 		{
// 		  goodMatch = true;
// 		} 

  CandList& elListL = _e.userCandList("ElectronLoose");
  CandList& elListT = _e.userCandList("ElectronTight");

  // two electrons from the tight list
  CandList& dielList = _e.userCandList( "Dielectron" );
  for( size_t iel1=0; iel1<elListT.size(); iel1++ )       
    {
      Candidate* el1 = elListT[iel1];
      for( size_t iel2=iel1+1; iel2<elListT.size(); iel2++ )
	{
	  Candidate* el2 = elListT[iel2];
	  Candidate* diel = Candidate::create( el1, el2 );
	  if( el1->charge()*el2->charge()>=0 ) continue; 
	  float mll = diel->mass();
	  if( mll<60 || mll>120 ) continue;
	  
	  dielList.push_back( diel );	  
	}
    }

  CandList& muListL = _e.userCandList("MuonLoose");
  CandList& muListT = _e.userCandList("MuonTight");

  // one muon from the tight list, one muon from the loose list
  CandList& dimuList = _e.userCandList( "Dimuon" );
  for( size_t imu1=0; imu1<muListL.size(); imu1++ )       
    {
      Candidate* mu1 = muListL[imu1];
      for( size_t imu2=imu1+1; imu2<muListL.size(); imu2++ )
	{
	  Candidate* mu2 = muListL[imu2];
	  if( CandUtil::overlap( *mu1, *mu2 ) ) continue;
	  Candidate* dimu = Candidate::create( mu1, mu2 );
	  bool isTight1 =  CandUtil::overlaps( *mu1, muListT );
	  bool isTight2 =  CandUtil::overlaps( *mu2, muListT );
	  if( !( isTight1 || isTight2 ) ) continue;  // at least one tight
	  if( mu1->charge()*mu2->charge()>=0 ) continue; 
	  float mll = dimu->mass();
	  if( mll<60 || mll>120 ) continue;
	  
	  dimuList.push_back( dimu );	  
	}
    }

  sort( dielList.begin(), dielList.end(), ZllCompare() );
  sort( dimuList.begin(), dimuList.end(), ZllCompare() );
  CandList& ZList = _e.userCandList("ZList");
  CandUtil::mergeLists( dielList, dimuList, ZList );
  sort( ZList.begin(), ZList.end(), ZllCompare() );
  if( ZList.size()==0 ) return false; 

  Candidate* ZCand = ZList[0];
  //  int nZCand = 1;
  bool isEE = CandUtil::overlaps( *ZCand, dielList );
  bool isMM = CandUtil::overlaps( *ZCand, dimuList );

  if( !( isEE || isMM ) ) return false;

  _N++;

  if( isEE ) nLepL = nEleL;
  if( isMM ) nLepL = nMuoL;

  if( _verbose )
    {
      PrintTree prtTree;
      cout << "*****************" << endl;
      cout << "MC Analysis -- Event " << _ievt << " -- categ " << _e.categ() << "  k=" << kfactor << endl;
      cout << prtTree.Print( *genCand_ );
      
      cout << "Number of electron VL/L/T/VT = " 
	   << nEleVL << "/"
	   << nEleL << "/"
	   << nEleT << "/"
	   << nEleVT << endl;
      cout << "Number of muons VL/L/T/VT = " 
	   << nMuoVL << "/"
	   << nMuoL << "/"
	   << nMuoT << "/"
	   << nMuoVT << endl;
      
      CandUtil::printList( cout, dielList );
      CandUtil::printList( cout, dimuList );
      ZCand->print(cout);
      if( isEE ) cout << "---> electron channel " << endl;
      if( isMM ) cout << "---> muon channel " << endl;
    }


  float mll = ZCand->mass();
  const Candidate* metCand =  _e.met( EventManager::kPatMet );
  assert( metCand!=0 );
  Candidate* ZnnCand = metCand->clone();
  //  Candidate* ZZCand = Candidate::create( ZCand, ZnnCand );
  float MET    = metCand->pt();
  //  float METPhi = metCand->phi();


  CandList leptonList[2] = { elListL, muListL };
  int nLoose = 0;

  vector< Candidate* > LCand;
  list< Candidate* > LCandAll;
  for( size_t ilep=0; ilep<2; ilep++ )
    {
      for( size_t ii=0; ii<leptonList[ilep].size(); ii++ )
	{
	  Candidate* LCand_ = leptonList[ilep][ii];
	  CandInfo* info = LCand_->info();
	  if( CandUtil::overlap( *LCand_, *ZCand ) ) 
	    {
	      LCand.push_back( LCand_ );
	    }
	  else
	    {
	      if( info->getBool("Loose") ) nLoose++;
	      LCandAll.push_back( LCand_ );
	    }
	}
    }
  assert( LCand.size()==2 );
  LCandAll.push_front( LCand[1] );
  LCandAll.push_front( LCand[0] );

  vector< float > LPt;
  vector<int> LPdg;
  vector< float > LEta;
  vector< float > LPhi;
  vector<int> LID;
  vector< float > LTrkIso;
  vector< float > LEcalIso;
  vector< float > LHcalIso;
  vector< float > LCombIso;
  vector<int> LVtx;
  float rhoFJ = _e.getRhoFastJet();

  for( list< Candidate* >::iterator it=LCandAll.begin(); 
       it!=LCandAll.end(); ++it )
    {
      Candidate* LCand_ = *it;  
      Vertex* vertex = LCand_->vertex();
      int vtxIndex_ = -1;
      if( vertex!=0 )
	{
	  CandInfo* vtxInfo_ = vertex->info();
	  if( vtxInfo_ )
	    {
	      vtxIndex_ = vtxInfo_->getInt("index");
	    }
	}
      int pdgCode_ = LCand_->pdgCode();
      LPt  .push_back( LCand_->charge() * LCand_->pt()  );
      LPdg .push_back( pdgCode_ );
      LEta .push_back( LCand_->eta() );
      LPhi .push_back( LCand_->phi() );
      LVtx .push_back( vtxIndex_ );
      int id_(0);
      CandInfo* info = LCand_->info();
      if( info->getBool("VeryLoose") ) id_ += 1;
      if( info->getBool("Loose") )     id_ += 2;
      if( info->getBool("Tight") )     id_ += 4;
      if( info->getBool("VeryTight") ) id_ += 8;
      LID  .push_back( id_ );

      // isolation variables
      float trkIso_(0);
      float ecalIso_(0);
      float hcalIso_(0);
      float combIso_(0);
     
      if( abs(pdgCode_)==13 )
	{
	  // muons
	  trkIso_ = info->getFloat("isoR03strk");
	  ecalIso_ = info->getFloat("isoR03sem");
	  hcalIso_ = info->getFloat("isoR03shad");
	  combIso_ = (trkIso_ + ecalIso_ + hcalIso_ - rhoFJ*0.11 )/ LCand_->pt();
	}
      else if( abs(pdgCode_)==11 )
	{
	  // electrons
	  trkIso_ = info->getFloat("dr03TrkIso");
	  ecalIso_ = info->getFloat("dr04EcalIso");
	  hcalIso_ = 
	    info->getFloat("dr04HcalD1Iso") + 
	    info->getFloat("dr04HcalD2Iso");	 

	  combIso_ = info->getFloat("combIso");
	}
      else
	assert(0);

      LTrkIso.push_back( trkIso_ );
      LEcalIso.push_back( ecalIso_ );
      LHcalIso.push_back( hcalIso_ );
      LCombIso.push_back( combIso_ );
      
    }

  float dPhiMin = 1000;
  for( size_t ii=0; ii<2; ii++ )
    {
      Candidate* LCand_ = LCand[ii];
      float dPhi_ = CandUtil::dPhi( ZnnCand, LCand_ ); 
      if( fabs( dPhi_ )<dPhiMin ) dPhiMin = fabs(dPhi_);
    }
  assert( dPhiMin<=Constants::pi );  

  // count the number of "good" vertices
  const VtxList& vertices = _e.vertices(); 
  bool findPV=false;
  int nVertex = 0;
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

  //The cleaned MET //MM
  Vertex* vtx = ZCand->daughter(0)->vertex();
  const Candidate* cpuMETVtxC  = _e.cleanPUmet(vtx , ZCand, nVertex );
  float cpuMETVtx = cpuMETVtxC->pt();

  float cpuMET = cpuMETVtx;
  float cpuMETMin = MET;
  if( cpuMET<MET ) cpuMETMin = cpuMET;
  float theMET = cpuMETMin;	      

  float alpha_0_;  
  float alpha_PU_; 
  float sigma_0_;  
  float sigma_PU_; 
  float alpha_ref_;
  float sigma_ref_;
  if( isEE )
    {
      alpha_ref_  = 4.997;
      sigma_ref_  = 3.962;
      alpha_0_    = 3.807;
      alpha_PU_   = 0.305;
      sigma_0_    = 2.924;
      sigma_PU_   = 1.267;
    }
  else
    {
      alpha_ref_  = 4.900;
      sigma_ref_  = 3.773;
      alpha_0_    = 3.776;
      alpha_PU_   = 0.315;
      sigma_0_    = 2.908;
      sigma_PU_   = 1.262;
    }

  float sigma_T  = 
    sqrt( pow(sigma_0_,2) + (nVertex-1)*pow(sigma_PU_,2) );
  float peak_T   = alpha_0_ + (nVertex-1)*alpha_PU_;
  float ghmSig = ( theMET - peak_T ) / sigma_T;
  float ghmMET = sigma_ref_ * ghmSig + alpha_ref_;
  



  string modestr("");
  if( isEE ) modestr = "e";
  else       modestr = "m";

  string selstr("");
  
  selstr = "presel";
  fill( "mll", modestr+"_"+selstr, mll, 1. );
  fill( "mll", modestr+"_"+selstr+"_KF", mll, kfactor );
  fill( "mll", modestr+"_"+selstr+"_KFJV", mll, kfactor_jv30 );

  // jet veto
  vector< float > jEt;
  vector< float > jEta;
  vector< float > jPhi;
  vector< float > jDR1;
  vector< float > jDR2;
  vector< float > jDPhiMET;
  vector< float > jBtagTkC;
  vector< float > jBtagSoftM;
  vector<int> jVtx; 


  const CandList& jets = _e.jetList( EventManager::kPatJet );
  CandList jList10; 
  CandList jList30; 

  int jMult10(0);
  int jMult30(0);
  float jBtagTkC_max   = -100.;
  float jBtagSoftM_max = -100.;

  for( size_t ij=0; ij<jets.size(); ij++ ) 
    { 
      Candidate* jet = jets[ij];
      Vertex* vertex = jet->vertex();
      int vtxIndex_ = -1;

      bool ismatch=false;
      for( list< Candidate* >::iterator it=LCandAll.begin();
	   it!=LCandAll.end(); ++it )
	{
	  Candidate* LCand_ = *it;
	  float dpt_ = (jet->pt()-LCand_->pt())/LCand_->pt();
	  if( jet->dR( LCand_ ) <0.2 && dpt_>0 && dpt_<0.5 )
	    { ismatch=true; break;}
	}
      if(ismatch) continue;

      //now photons
      //      for(size_t ip=0;ip<SelPh.size(); ip++) {
      //	Candidate* pCand_ = SelPh[ip];
      //	float dpt_ = (jet->pt()-pCand_->pt())/pCand_->pt();
      //	if( jet->dR( pCand_ ) <0.2 && dpt_>0 && dpt_<0.5 )
      //	    { ismatch=true; break; }
      //      }
      // if(ismatch) continue; //!!!
      
      if( vertex!=0 )
	{
	  CandInfo* vtxInfo_ = vertex->info();
	  if( vtxInfo_ )
	    {
	      vtxIndex_ = vtxInfo_->getInt("index");
	    }
	}
      float jEta_ = jet->eta();
      float jEt_  = jet->Et();
      float jPhi_ = jet->phi();
      float jDPhiMET_ = CandUtil::dPhi( jet, ZnnCand )*Constants::radToDeg ;
      float jDR1_ = CandUtil::dPhi( jet, LCand[0] )*Constants::radToDeg ;
      float jDR2_ = CandUtil::dPhi( jet, LCand[1] )*Constants::radToDeg ;
      CandInfo* info_ = jet->info();
      float jBtagTkC_   = info_->getFloat("btagTkC");
      float jBtagSoftM_ = info_->getFloat("btagSoftM");

      
      if( fabs(jEta_)<=5 )
	if( jEt_>=10 )
	{	  
	  // b-tagging
	  if( jBtagTkC_>jBtagTkC_max ) 
	    jBtagTkC_max = jBtagTkC_;
	  if( jBtagSoftM_>jBtagSoftM_max ) 
	    jBtagSoftM_max = jBtagSoftM_;

	  jPhi    .push_back(jPhi_ );
	  jEta    .push_back(jEta_ );
	  jEt     .push_back(jEt_ );
	  jDR1    .push_back(jDR1_ ); 
	  jDR2    .push_back(jDR2_ ); 
	  jDPhiMET.push_back(jDPhiMET_ ); 
	  jBtagTkC.push_back(jBtagTkC_);
	  jBtagSoftM.push_back(jBtagSoftM_);
	  jVtx    .push_back(vtxIndex_);
	  	  	  
	  //	  if( jEt_<30 )
	    {
	      jList10.push_back( jet );
	      jMult10++;
	    }
	    if( jEt_>=30 )
	    {
	      jMult30++;
	      jList30.push_back( jet );
	    }
	}      
    }
  sort( jList10.begin(), jList10.end(), Candidate::SortByPt() );
  sort( jList30.begin(), jList30.end(), Candidate::SortByPt() );
    
  float dPhiJMET = -190.;
  if( jList10.size()>0 ) dPhiJMET = jDPhiMET[0];
  
  TVector2 ZCand2D   = ZCand->p2();
  TVector2 ZnnCand2D = ZnnCand->p2();
  TVector2 minusZnnCand2D = ZnnCand2D;
  minusZnnCand2D *= -1;
  float dPhiZMET = KineUtils::dPhi( ZCand2D.Phi(), minusZnnCand2D.Phi() )*Constants::radToDeg;

  float qTll = ZCand->pt();
  float bal = MET / qTll;

  // presel
  //  OneSidedCut<int> c_presel   ( "presel",    nLepL, 2, ">=" );

  // Z quality
  bool ZQual = true;
  OneSidedCut<bool>  c_ZQual  ( "ZLQual",     ZQual, true,   "==" );

  // Z qT
  OneSidedCut<float> c_qT30   ( "qT30", qTll, 30., ">=" );

  // jet veto
  OneSidedCut<int>   c_j30Veto ( "j30Veto",    jMult30, 0,      "==" );

  // MET cut
  OneSidedCut<float> c_MET     ( "MET",  ghmMET,  50.,   ">"  );

  // balance
  TwoSidedCut<float> c_bal  ( "bal",    bal, 0.4, 1.8, "[]" ); 

  // jet/MET quality
  TwoSidedCut<float> c_dPhiJMET( "dPhiJMET", dPhiJMET, -20.,20, "![]");

  // Z/MET quality
  TwoSidedCut<float> c_dPhiZMET ( "dPhiZMET", dPhiZMET, -60, 60, "[]" ); //!!! GHM

  // B-tagging
  OneSidedCut<float> c_BtagTkCVeto   ( "BtagTkCVeto",   jBtagTkC_max, 2.5, "<=" );
  OneSidedCut<float> c_BtagSoftMVeto ( "BtagSoftMVeto", jBtagSoftM_max, 0.30, "<=" );
  CompositeCut c_BtagVeto( "BTagVeto", "&&" );
  c_BtagVeto.add( c_BtagTkCVeto    );
  c_BtagVeto.add( c_BtagSoftMVeto  );	  

  // diboson reduction
  OneSidedCut<int>   c_lVeto    ( "lVeto",       nLoose, 0,      "==" );

  TwoSidedCut<float> c_massCut  ( "massCut", mll, 80., 100, "[]");

    // final selection
  CompositeCut c_finalSelection( "final", "&&" );
  c_finalSelection.add( c_ZQual       );
  c_finalSelection.add( c_qT30        );
  c_finalSelection.add( c_j30Veto     );
  c_finalSelection.add( c_lVeto       );
  c_finalSelection.add( c_MET         ); 
  c_finalSelection.add( c_bal         );
  c_finalSelection.add( c_dPhiJMET    );
  c_finalSelection.add( c_dPhiZMET    );
  c_finalSelection.add( c_BtagVeto    );
  c_finalSelection.add( c_massCut     );

  for( size_t ic=0; ic<c_finalSelection.size(); ic++ )
    {
      string name_ = c_finalSelection.name(ic);
      bool ok = c_finalSelection.cutUpTo( ic );
      if( ok )
	{
	  selstr = name_;
	  fill( "mll", modestr+"_"+selstr, mll, 1. );
	  fill( "mll", modestr+"_"+selstr+"_KF", mll, kfactor );
	  fill( "mll", modestr+"_"+selstr+"_KFJV", mll, kfactor_jv30 );
	}
    }


  return true;
  
}

void
MCAnalysis::bookHistograms()
{
  //
  // histogram booking
  //

}

void
MCAnalysis::writeHistograms()
{

  vector<string>  categ;
  categ.push_back( "ZZ_2e2n" );
  categ.push_back( "ZZ_2m2n" );
  categ.push_back( "ZZ_2e2t6n" );
  categ.push_back( "ZZ_2m2t6n" );
  categ.push_back( "WZ_3en" );
  categ.push_back( "WZ_2emn" );
  categ.push_back( "WZ_e2mn" );
  categ.push_back( "WZ_3mn" );
  categ.push_back( "WZ_3e3t7n" );
  categ.push_back( "WZ_2em3t7n" );
  categ.push_back( "WZ_e2m3t7n" );
  categ.push_back( "WZ_3m3t7n" );
  for( size_t ic=0; ic<categ.size(); ic++ )
    {
      string categ_ = categ[ic];
      float n_  = _counters[categ_ + "__N"];
      cout << categ_ << "-Number of events in sample " << n_ << endl;
      if( n_==0 ) continue;
      float nTot_ = _counters[categ_ + "__N_KF"];
      vector<string> selstr;
      selstr.push_back(categ_ + "__N_KF_lepEta");
      selstr.push_back(categ_ + "__N_KF_lepPt");
      selstr.push_back(categ_ + "__N_KF_Acc");
      
      cout << categ_ << "-Number of events (NLO) " << nTot_ << endl;
      for( size_t ii=0; ii<selstr.size(); ii++ )
	{
	  const string & str_ = selstr[ii];
	  float nSelected_ = _counters[str_];
	  float nAnalyzed_ = nTot_;
	  if( nAnalyzed_<=0 ) continue;
	  double p_   = (1.*nSelected_)/nAnalyzed_;
	  double var_ = nAnalyzed_*p_*(1-p_);
	  double sig_ = sqrt(var_)/nAnalyzed_;
	  printf( "Sel/Tot[%-30s]=%8.2f/%8.2f --> (%6.2f+/-%-5.2f)%1s\n", str_.c_str(), nSelected_, nAnalyzed_, p_*100, sig_*100, "%" );
	}
    }
  
}
