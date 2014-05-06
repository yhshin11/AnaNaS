#include "Analysis/src/WlnAnalysis.hh"

#include <cassert>
#include <algorithm>
using namespace std;

#include <TH1I.h>
#include <TH1F.h>
#include <TH2I.h>
#include <TH2F.h>

#include "Analysis/utils/Config.hh"
#include "Analysis/utils/Constants.hh"
#include "Analysis/core/Candidate.hh"
#include "Analysis/core/CandInfo.hh"
#include "Analysis/core/EventManager.hh"
#include "Analysis/tools/CandUtil.hh"
#include "Analysis/selectors/ElectronSelector.hh"
#include "Analysis/selectors/MuonSelector.hh"
#include "Analysis/selectors/LeptonSelector.hh"
#include "Analysis/selectors/IsolationSelector.hh"
#include "Analysis/selectors/ConeSelector.hh"

ClassImp( WlnAnalysis )

string WlnAnalysis::modeName[WlnAnalysis::kNModes] = 
{
  "all",
  "W_en",      "W_mn"
};

string WlnAnalysis::levelName[WlnAnalysis::kNLevels] = 
{
  "nosel", "presel", "sel"   
};

WlnAnalysis::WlnAnalysis( Sample& sample, const std::string & collectionFileName ) : 
  SampleAnalysis( "Wln", sample, collectionFileName ), 
  PseudoExp( sample )
{
  cout << "\t-------------------------------------------" << endl; 
  cout << "\t---  W -> lnu preselection analysis    ----" << endl; 
  cout << "\t-------------------------------------------" << endl; 

  _hlt = false; // apply HLT selection

  //  _hltLines.push_back("HLT_DoubleEle10_SW_L1R");
  //  _hltLines.push_back("HLT_DoubleMu3");
  //  _hltLines.push_back("HLT_Mu9");
  //  _hltLines.push_back("HLT_Mu15");
  _hltLines.push_back("HLT_Ele15_SW_L1R");   // 20 or 25

  // lines at startup 2010@7TeV
  //  _hltLines.push_back("HLT_L2Mu9");
  //  _hltLines.push_back("HLT_L1SingleEG8");

  _nDebug = 1;

  //
  // pseudo experiments
  //
  setLumi( 39000.0 );

  // define variables
  createVar( "pt",  -100., 100. );
  createVar( "eta", -3, 3. );
  createVar( "phi", -360, 360. );//?
  createVar( "pf_mt",  0., 1000. );
  createVar( "pf_met", 0., 1000. );
  createVar( "pf_met_lx", -1000., 1000. );
  createVar( "pf_met_ly", -1000., 1000. );
}

WlnAnalysis::~WlnAnalysis()
{
}

void
WlnAnalysis::bookHistograms()
{
  // projection histograms
  defineTemplate( "Eproj", 400, -200, 200 );
  defineTemplate( "Eproj2D", 200, -200, 200, 200, -200, 200 );
  defineTemplate( "METVsMT", 80, 0, 80, 80, 0, 160 );
}

void
WlnAnalysis::writeHistograms()
{
  for( map< string, TH1* >::const_iterator it=h_.begin(); 
       it!=h_.end(); it++ )
    {
      it->second->Write();
    }
  for( map< string, TH2* >::const_iterator it=h2_.begin(); 
       it!=h2_.end(); it++ )
    {
      it->second->Write();
    }  

  // temporary
  string categ;
  int nall(0), nin(0);
  for( map< pair< string, string >, int >::const_iterator it=_sel.begin();
       it!=_sel.end(); it++ )
    {
      const string & categ_ = it->first.first;
      const string & mode_  = it->first.second;
      if( categ_!=categ )
	{
	  categ = categ_; 
	  map< string, int >::const_iterator it_;
	  it_ = _all.find( categ );
	  nall = it_->second;
	  it_ = _in.find( categ_ );
	  nin =  it_->second;
	  cout << nin << "/" << nall << endl;
	  float fin = (1.*nin)/nall;
	  printf( "\nn[%-6s]=%6d, %6.2f%1s in acceptance\n", categ.c_str(), nall, fin*100., "%" );
	}

      int n_ = nall;
      if( categ==mode_ ) n_ = nin;
      int nsel = it->second;
      float p_ = (1.*nsel)/n_;
      float var_ = n_*p_*(1-p_);
      float sig_ = sqrt(var_)/n_;

      printf( "n[%-6s, %-6s]=%6d --> (%6.2f+/-%-5.2f)%1s\n", categ.c_str(), mode_.c_str(), nsel, p_*100., sig_*100., "%" );
    }
  cout << endl;

  writeDatasets();
}

bool
WlnAnalysis::analyzeEvent()
{
  // event initialisations
  for( size_t lev=kAll; lev<kNLevels; lev++ )
    {
      _nMode[lev].clear();
    }
  _nMode[kNosel][kAll]=1;

  // build lists of leptons and prepare the MC matching map
  //  buildLeptonLists();

  // prepare the lists of veryloose electrons
  ElectronSelector elSelVL( LeptonSelector::kVeryLoose ); 
  CandList& veryLooseElectrons = _e.userCandList("ElectronVeryLoose");
  elSelVL.getList( _e.electrons(), veryLooseElectrons );
  sort( veryLooseElectrons.begin(), veryLooseElectrons.end(), LeptonCompare() );

  // loop on pre-selected hypotheses
  int nSelected(0);

  if( W_en__analysis() ) nSelected++;

  // debug printouts
  debugPrintouts( cout );

  return nSelected>0;
}


bool
WlnAnalysis::W_en__analysis()
{
  // pseudo experiments
  checkEvent( _ievt );

  int imode(0);
  string prefix(modeName[kW_en]); 

  TString leptonListName( "ElectronVeryLoose" );
  CandList& leptonList_ = _e.userCandList( leptonListName.Data() );
  
  if( leptonList_.size() == 0 ) return false;

  // take the "best" lepton of the list
  Candidate& cand = *leptonList_[0];

  // isolation
  unsigned settings;
  if( imode==0 ) settings=EventManager::kElectron;  // !!!GHM
  if( imode==1 ) settings=EventManager::kMuon;
  //  IsolationSelector isoSel_(settings);
  //  isoSel_.computeIso( cand );

  CandInfo* info = cand.info();
  float iso(0.), S_trk(0.), S_ecal(0.), S_hcal(0.);
  info->getFloat("iso", iso );
  info->getFloat("S_trk",  S_trk  );
  info->getFloat("S_ecal", S_ecal );
  info->getFloat("S_hcal", S_hcal );
//   int id(-1);
//   for( int lepSel=0; lepSel<LeptonSelector::kNSel; lepSel++ )
//     {
//       string str_ = LeptonSelector::selection[lepSel];
//       if( info->getBool(str_.c_str()) )  id = lepSel;
//     }

  float pt  = cand.pt();
  //  float eta = cand.eta();

  //  float caloE   = info->getFloat("caloEnergy");
  float caloEta = info->getFloat("caloEta");

  // should be the SC eta
  bool isEB = fabs(caloEta)<=1.44;
  bool isEE = fabs(caloEta)>=1.56 && fabs(caloEta)<=2.5;
  bool isEcalFiducial = isEB || isEE;

  

  float dEtaIn      = info->getFloat("dEtaIn");
  float dPhiIn      = info->getFloat("dPhiIn");
  float sigiEtaiEta = info->getFloat("sigiEtaiEta");
  float HOE         = info->getFloat("HOE");

  // get the unit vectors along and perpendicular to LCand
  TVector2 ux_ = cand.p2().Unit();
  TVector2 uy_( -ux_.Y(), ux_.X() ); 

  // charge 
  int chg = ( cand.charge()>0 ) ? 1:-1;

  // pdgCode
  int LPdgCode_  = cand.pdgCode();
  int nuPdgCode_ = -1*(LPdgCode_-chg);

  //
  // from Georgios Daskalakis
  //
  // http://www.hep.ph.ic.ac.uk/~nrompoti/ForOctX/7TeVEleSelectionWenuZeePlotCollection/WENU/
  // http://www.hep.ph.ic.ac.uk/~nrompoti/ForOctX/7TeVEleSelectionWenuZeePlotCollection/ZEE/
  // http://www.hep.ph.ic.ac.uk/~nrompoti/ForOctX/7TeVEleSelectionWenuZeePlotCollection/WENU/wenuSummary.txt
  // http://www.hep.ph.ic.ac.uk/~nrompoti/ForOctX/7TeVEleSelectionWenuZeePlotCollection/ZEE/zeeSummary.txt

  bool is70(true);
  bool is80(true);
  bool is90(true);
  bool is95(true);
  bool isAnti(true);
  {
    if( isEB )
      {
	// identification
	if(       sigiEtaiEta  >0.01  ) is70=false;
	if( fabs( dPhiIn      )>0.02  ) is70=false;
	if( fabs( dEtaIn      )>0.006 ) is70=false;
	if(       HOE          >0.02  ) is70=false;
	
	// isolation
	if( S_trk >2.5 ) is70=false;
	if( S_ecal>3.0 ) is70=false;
	if( S_hcal>5.0 ) is70=false;
      }
    else if( isEE )
      {
	// identification
	if(       sigiEtaiEta  >0.03   ) is70=false;
	if( fabs( dPhiIn      )>0.02   ) is70=false;
	if( fabs( dEtaIn      )>0.003  ) is70=false;
	if(       HOE          >0.0025 ) is70=false;
	
	// isolation
	if( S_trk >0.8  ) is70=false;
	if( S_ecal>2.5  ) is70=false;
	if( S_hcal>0.25 ) is70=false;
      }
  }
  if( !is70 ) {
    if( isEB )
      {
	// identification
	if(       sigiEtaiEta  >0.01  ) is80=false;
	if( fabs( dPhiIn      )>0.02  ) is80=false;
	if( fabs( dEtaIn      )>0.006 ) is80=false;
	if(       HOE          >0.05  ) is80=false;
	
	// isolation
	if( S_trk >3.0 ) is80=false;
	if( S_ecal>4.0 ) is80=false;
	if( S_hcal>5.0 ) is80=false;
      }
    else if( isEE )
      {
	// identification
	if(       sigiEtaiEta  >0.03   ) is80=false;
	if( fabs( dPhiIn      )>0.02   ) is80=false;
	if( fabs( dEtaIn      )>0.006  ) is80=false;
	if(       HOE          >0.025  ) is80=false;
	
	// isolation
	if( S_trk >1.5  ) is80=false;
	if( S_ecal>2.5  ) is80=false;
	if( S_hcal>0.7  ) is80=false;
      }
  }
  if( !is80 ) {
    if( isEB )
      {
	// identification
	if(       sigiEtaiEta  >0.01  ) is90=false;
	if( fabs( dPhiIn      )>0.04  ) is90=false;
	if( fabs( dEtaIn      )>0.006 ) is90=false;
	if(       HOE          >0.05  ) is90=false;
	
	// isolation
	if( S_trk >6.0 ) is90=false;
	if( S_ecal>5.0 ) is90=false;
	if( S_hcal>5.0 ) is90=false;
      }
    else if( isEE )
      {
	// identification
	if(       sigiEtaiEta  >0.03   ) is90=false;
	if( fabs( dPhiIn      )>0.025  ) is90=false;
	if( fabs( dEtaIn      )>0.008  ) is90=false;
	if(       HOE          >0.025  ) is90=false;
	
	// isolation
	if( S_trk >6.0  ) is90=false;
	if( S_ecal>2.5  ) is90=false;
	if( S_hcal>1.5  ) is90=false;
      }
  }
  if( !is95 ) {
    if( isEB )
      {
	// identification
	if(       sigiEtaiEta  >0.01  ) is95=false;
	if( fabs( dPhiIn      )>0.80  ) is95=false;
	if( fabs( dEtaIn      )>0.006 ) is95=false;
	if(       HOE          >0.05  ) is95=false;
	
	// isolation
	if( S_trk >7.0 ) is95=false;
	if( S_ecal>5.0 ) is95=false;
	if( S_hcal>5.0 ) is95=false;
      }
    else if( isEE )
      {
	// identification
	if(       sigiEtaiEta  >0.03   ) is95=false;
	if( fabs( dPhiIn      )>0.7    ) is95=false;
	if( fabs( dEtaIn      )>0.008  ) is95=false;
	if(       HOE          >0.04   ) is95=false;
	
	// isolation
	if( S_trk >8.0  ) is95=false;
	if( S_ecal>3.0  ) is95=false;
	if( S_hcal>2.0  ) is95=false;
      }
  }

  {
    if( isEB )
      {
	// identification
	if(       sigiEtaiEta  >0.01   ) isAnti=false;
	//	if( fabs( dPhiIn      )>0.02   ) isAnti=false;
	//	if( fabs( dEtaIn       )>0.006 ) isAnti=false;

	if( fabs( dPhiIn      )<0.04   ) isAnti=false;
	if( fabs( dEtaIn      )<0.006  ) isAnti=false;

	if(       HOE          >0.02   ) isAnti=false;
	
	// isolation
	if( S_trk >2.5 ) isAnti=false;
	if( S_ecal>3.0 ) isAnti=false;
	if( S_hcal>5.  ) isAnti=false;
      }
    else if( isEE )
      {
	// identification
	if(       sigiEtaiEta  >0.03   ) isAnti=false;

	//	if( fabs( dPhiIn      )>0.02   ) isAnti=false;
	//	if( fabs( dEtaIn      )>0.003  ) isAnti=false;
	if( fabs( dPhiIn      )<0.025  ) isAnti=false;
	if( fabs( dEtaIn      )<0.008  ) isAnti=false;

	if(       HOE          >0.0025 ) isAnti=false;
	
	// isolation
	if( S_trk >0.8  ) isAnti=false;
	if( S_ecal>2.5  ) isAnti=false;
	if( S_hcal>0.25 ) isAnti=false;
      }
  }

  bool is25 = pt>=25;
  bool is27 = pt>=27;
  bool is30 = pt>=30;
  bool is32 = pt>=32;

  //  int ipt=2;
  for( int ipt=0; ipt<5; ipt++ )
    {
      string suffix1("sel");
      if( ipt==1 )
	{
	  suffix1+="_pt25";
	  if( !is25 ) continue;
	} 
      else if( ipt==2 )
	{
	  //	  continue;
	  suffix1+="_pt27";
	  if( !is27 ) continue;
	} 
      else if( ipt==3 )
	{
	  suffix1+="_pt30";
	  if( !is30 ) continue;
	} 
      else if( ipt==4 )
	{
	  suffix1+="_pt32";
	  if( !is32 ) continue;
	} 
      else
	{
	  suffix1+="_all";
	}
	 
      int ireg=1;
      //      for( int ireg=0; ireg<4; ireg++ )
	{
	  string suffix2(suffix1);
	  if( ireg==1 )
	    {
	      suffix2+="_Ecal";
	      if( !isEcalFiducial ) continue;
	    }
	  else if( ireg==2 )
	    {
	      suffix2+="_EB";
	      if( !isEB ) continue;
	    }
	  else if( ireg==3 )
	    {
	      suffix2+="_EE";
	      if( !isEE ) continue;
	    }
	  else
	    {
	      continue;
	      //	      suffix2+="_all";
	    }


	  int icut=0;
	  //	  for( int icut=0; icut<5; icut++ )
	    {
	      
	      string suffix(suffix2);
	      if( icut==1 )
		{
		  suffix += "_cut90";
		  if( !is90 ) continue;
		}
	      else if( icut==2 ) 
		{
		  suffix += "_cut80";
		  if( !is80 ) continue;
		}
	      else if( icut==3 ) 
		{
		  suffix += "_cut70";
		  if( !is70 ) continue;
		}
	      else if( icut==4 )
		{
		  suffix += "_anticut";
		  if( !isAnti ) continue;
		}
	      else
		{
		  suffix+="_presel";
		}

	      if( imode==0 )
		{
		  //		  fillElectronHistograms( cand, suffix, prefix );
		}
	      else if( imode==1 )
		{
		  //		  fillMuonHistograms( cand, suffix, prefix );
		}
      
	      TVector2 neuP2 = _e.neutrinosFromBosonsP4().Vect().XYvector();
	      string sgn_ = "neg";
	      if( chg==1 ) sgn_ = "pos"; 

	      map< string, const Candidate* > metCand;
	      map< string, Candidate* >       LNuCand;
	      map< string, Candidate* >         WCand;
	      metCand["calo"] =  _e.met("Calo");
	      metCand["paf"]  =  _e.met("Pat");  // this is the tcMET
	      metCand["pf"]   =  _e.met("ParticleFlow");
	      //	      for( map<string,const Candidate*>::const_iterator it_ = metCand.begin();
	      //		   it_!=metCand.end(); ++it_ )
		{
		  string str_ = "pf";
		  //		  if( str_!="pf" ) continue; //!!!
		  //		  const Candidate* metCand_ = it_->second;
		  const Candidate* metCand_ = metCand[str_];
		  assert( metCand_!=0 );
		  Candidate* nuCand_ = metCand_->clone();
		  nuCand_->setPdgCode( nuPdgCode_ ); 
		  TVector2 nuP2_ = nuCand_->p2();
		  Candidate* LNuCand_ = Candidate::create( &cand, nuCand_ );
		  LNuCand_->setName("LNu");
		  LNuCand[str_] = LNuCand_;
		  float MET_   = nuCand_->pt();
		  //		  float MTref_ = mT( &cand, nuCand_ );
		  float MT_    = LNuCand_->mass();

		  bool isMT40 = MT_>=40;
		  bool isMET20 = MET_>=20;
		  bool isMET27 = MET_>=27;

		  float dphi_  = nuCand_->p2().DeltaPhi( neuP2 ) * Constants::radToDeg;
		  float dpt_   = nuCand_->pt()-neuP2.Mod();
		  fill( "met"    , suffix+"__"+str_,           MET_,    prefix );
		  //		  fill( "mt"     , suffix+"__"+str_+"__ref",   MTref_,  prefix );
		  fill( "mt"     , suffix+"__"+str_,           MT_,     prefix );
		  //		  fill( "dphi"   , suffix+"__"+str_+"__"+sgn_, dphi_,   prefix );
		  //		  fill( "dpt"    , suffix+"__"+str_+"__"+sgn_, dpt_,    prefix );
		  fill( "METVsMT", suffix+"__"+str_, MET_, MT_, prefix );
		  if( isMT40 )
		    {
		      fill( "met"    , suffix+"__"+str_+"__MT40",           MET_,    prefix );		      
		    }
		  if( isMET20 )
		    {
		      fill( "mt"    , suffix+"__"+str_+"__MET20",      MT_,    prefix );		      
		      if( isMET27 )
			{
			  fill( "mt"    , suffix+"__"+str_+"__MET27",  MT_,    prefix );		      
			}
		    }
		  
		  // projections
		  float proj_x_ = nuP2_ * ux_;
		  float proj_y_ = nuP2_ * uy_;
		  fill( "Eproj"  , suffix+"__"+str_+"__ux",      proj_x_ , prefix );
		  fill( "Eproj"  , suffix+"__"+str_+"__uy",      proj_y_ , prefix );
		  //		  fill( "Eproj2D" , suffix+"__"+str_, proj_x_, proj_y_ , prefix );
		  //		  fill( "Eproj"  , suffix+"__"+str_+"__ux"+"__"+sgn_,      proj_x_ , prefix );
		  //		  fill( "Eproj"  , suffix+"__"+str_+"__uy"+"__"+sgn_,      proj_y_ , prefix );
		  //		  fill( "Eproj2D" , suffix+"__"+str_+"__"+sgn_, proj_x_, proj_y_ , prefix );

		  Candidate* WCand_ = LNuCand_->clone();
		  WCand_->setPdgCode( chg*24 );

		  //		  if( str_=="pf" && ipt==1 && ireg==1 && icut==2 )
		  if( str_=="pf" )
		    {
		      // cout << suffix << endl;
		      //		      select_ = true;
		      //		      TString toto_(str_);
		      setVarValue( "pt",      cand.charge()*cand.pt()     );
		      setVarValue( "eta",     cand.eta()     );
		      setVarValue( "phi",     cand.phi()     );
		      setVarValue( "pf_met",    MET_    );
		      setVarValue( "pf_met_lx", proj_x_ );
		      setVarValue( "pf_met_ly", proj_y_ );
		      saveEvent();
		    }

		  if( _verbose )
		    {
		      cout << "W Candidate with " << str_ << " met " << endl;
		      WCand_->oneLine( cout );
		      WCand_->theBase()->oneLine( cout );
		    }
		}
	    }
	}
    }
  //  hypothesis is selected
  //  if( select_ ) select( hyp, kmode, kSel );

  return true;
}


void
WlnAnalysis::debugPrintouts( ostream& o )
{
  // counters
  string _categ = _e.categ();
  bool _inAcceptance = _e.inAcceptance();
  _all[_categ]++;
  if( _inAcceptance ) _in[_categ]++;
  ModeStatIterator it = _nMode[kSel].begin();
  for( ; it!=_nMode[kSel].end(); it++ )
    {
      string mode_ = modeName[it->first];
      if( _categ==mode_ && !_inAcceptance ) mode_+="_out";
      pair< string, string >  p_( _categ, mode_ );
      _sel[p_]++;
    }
  
}

float
WlnAnalysis::mT( const Candidate* cand, const Candidate* met )
{
  //   returns  the transverse mass for a candidate and a missing momentum
  TVector2 metP2  =  met->p2();
  TVector2 candP2 = cand->p2();
  
  float mT2 =  2*( candP2.Mod() * metP2.Mod() -  candP2 * metP2 );

  if( mT2<0 ) mT2=0;
  
  return sqrt( mT2 );
}


