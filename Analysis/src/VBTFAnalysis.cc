 #include "Analysis/src/VBTFAnalysis.hh"

 #include <cassert>
 #include <algorithm>
 using namespace std;

 #include <TH1I.h>
 #include <TH1F.h>
 #include <TH2I.h>
 #include <TH2F.h>

 #include "Analysis/utils/Config.hh"
 #include "Analysis/utils/Constants.hh"
 #include "Analysis/utils/KineUtils.hh"
 #include "Analysis/core/Candidate.hh"
 #include "Analysis/core/CandInfo.hh"
 #include "Analysis/core/EventManager.hh"
 #include "Analysis/tools/CandUtil.hh"
 #include "Analysis/selectors/ElectronSelector.hh"
 #include "Analysis/selectors/VBTFElectronSelector.hh"
 #include "Analysis/selectors/MuonSelector.hh"
 #include "Analysis/selectors/LeptonSelector.hh"
 #include "Analysis/selectors/IsolationSelector.hh"
 #include "Analysis/selectors/ConeSelector.hh"

 ClassImp( VBTFAnalysis )

 string VBTFAnalysis::lmodeName[VBTFAnalysis::kNLModes] = 
 {
   "electron", "muon"
 };

 string VBTFAnalysis::analysisName[VBTFAnalysis::kNLModes] = 
 {
   "VBTF_e", "VBTF_m"
 };

 string VBTFAnalysis::modeName[VBTFAnalysis::kNModes] = 
 {
   "all",
   "W_en", "Z_2e",  
   "W_mn", "Z_2m"
 };

 string VBTFAnalysis::levelName[VBTFAnalysis::kNLevels] = 
 {
   "nosel", "HLT", "presel", "sel"   
 };

 VBTFAnalysis::VBTFAnalysis( Sample& sample, unsigned imode, const std::string & collectionFileName ) : 
   SampleAnalysis( analysisName[imode], sample, collectionFileName ), PseudoExp( sample ), _pseudo(false), _imode(imode), _theTuple(0)
 {
   cout << "\t-------------------------------------------" << endl; 
   cout << "\t---  VBTF analysis                     ----" << endl; 
   cout << "\t-------------------------------------------" << endl; 

   _hlt = true; // apply HLT selection

   if( _imode==kElectron )
     {
       // electron line
       _hltLines.push_back("HLT_L1SingleEG8");
       _hltLines.push_back("HLT_DoubleEle10_SW_L1R");
       _hltLines.push_back("HLT_Ele15_SW_L1R");   // 20 or 25
     }
   else if( _imode==kMuon )
     {
       // muon lines
       _hltLines.push_back("HLT_DoubleMu3");
       _hltLines.push_back("HLT_Mu9");
       _hltLines.push_back("HLT_Mu15");
       _hltLines.push_back("HLT_L2Mu9");
     }
   else
     abort();

   _nDebug = 10;
 }

 VBTFAnalysis::~VBTFAnalysis()
 {
 }

 void
 VBTFAnalysis::bookHistograms()
 {
   // projection histograms
   defineTemplate( "Eproj", 400, -200, 200 );
   defineTemplate( "Eproj2D", 200, -200, 200, 200, -200, 200 );

   TString name = "SelectedCandidates";
   _indexVar.push_back("run");
   _indexVar.push_back("event");
   _indexVar.push_back("e");
   _indexVar.push_back("et");
   _indexVar.push_back("fid");
   _indexVar.push_back("sel");
   _indexVar.push_back("eid");
   _indexVar.push_back("iso");
   _indexVar.push_back("pt");
   _indexVar.push_back("eta");
   _indexVar.push_back("phi");
   _indexVar.push_back("eop");
   _indexVar.push_back("hoe");
   _indexVar.push_back("sigieie");
   _indexVar.push_back("dphi");
   _indexVar.push_back("deta");
   _indexVar.push_back("strk");
   _indexVar.push_back("secal");
   _indexVar.push_back("shcal");
   _indexVar.push_back("calo_met");
   _indexVar.push_back("calo_mt");
   _indexVar.push_back("calo_dph");
   _indexVar.push_back("calo_dpt");
   _indexVar.push_back("calo_lproj_x");
   _indexVar.push_back("calo_lproj_y");
   _indexVar.push_back("calo_wproj_x");
   _indexVar.push_back("calo_wproj_y");
   _indexVar.push_back("pf_met");
   _indexVar.push_back("pf_mt");
   _indexVar.push_back("pf_dph");
   _indexVar.push_back("pf_dpt");
   _indexVar.push_back("pf_lproj_x");
   _indexVar.push_back("pf_lproj_y");
   _indexVar.push_back("pf_wproj_x");
   _indexVar.push_back("pf_wproj_y");

   //  map< string, float >::const_iterator it_ = _tupleVar.begin();
   //  TString varNames_(it_->first); it_++;
   //  for( ; it_!=_tupleVar.end(); it_++ )
   TString varNames_ = _indexVar[0];
   for( size_t ii=1; ii<_indexVar.size(); ii++ )
     {
       varNames_ += ":"; varNames_ += _indexVar[ii];
     }
   cout << varNames_ << endl;
   _theTuple = 
     new TNtuple( name, name, varNames_ );  

 }

 void
 VBTFAnalysis::writeHistograms()
 {
   _theTuple->Write();

   // temporary
   string & categ;
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
 }

bool
VBTFAnalysis::analyzeEvent()
{
  // event initialisations
  for( size_t lev=kAll; lev<kNLevels; lev++ )
    {
      _nMode[lev].clear();
    }
  _nMode[kNosel][kAll]=1;
  

  // build VBTF lepton lists
  buildLeptonLists();

  buildVBTFLeptonLists();

  // build the combinatorics of Z candidates 
  // (Z-hypothesis = sorted lists of non-overlapping Z candidates and leptons)
  nZHypos = buildZllLists();

  unsigned nSelected(0);

  if( nZHypos>0 )
    {
      // there is at least one Z candidate
    }
  else
    {
      if( W_ln__analysis() ) nSelected++;
    }

  // stat histograms
  fillStatHistograms();

  // debug printouts
  debugPrintouts( cout );

  return nSelected>0;
}

void
VBTFAnalysis::select( Mode mode, Level lev )
{
  _nMode[lev][mode]++;
}

size_t
VBTFAnalysis::buildZllLists()
{
  // get the lepton list
  CandList& leptonList = _e.userCandList( leptonListName() );
  
//   for( size_t ii=0; ii< leptonList.size(); ii++ )
//     {
//       leptonList[ii]->oneLine(cout);
//     }

  const CandList& initialList = _e.compositeCandList( EventServer::Zmode[_imode] );

  // the list of selected Z candidates
  CandList& list = _e.userCandList( ZListName() );
  for( size_t iZ=0; iZ<initialList.size(); iZ++ )
    {
      Candidate* cand = initialList[iZ];

      // kinematics of the candidate
      float mass_ = cand->mass();
      
      // mass cut
      bool massOK_ =  mass_>=50.;

      bool idOK_(true);
      for( size_t idau=0; idau<cand->nDaughters(); idau++ )
	{
	  const Candidate* dau = cand->daughter(idau);
	  if( !CandUtil::overlaps( *dau, leptonList ) ) idOK_=false;
	}

      if( massOK_ && idOK_ )
	{
	  list.push_back( cand );
	}
    }

  return list.size();
}

void
VBTFAnalysis::print( ostream& o ) const
{
  o << "*****************" << endl;
  o << "VBTF-like Analysis -- Event " << _ievt << " -- categ " << _e.categ() << endl;
  o << "-+-+-+-+-+" << endl;
  
  const CandList& leptonList 
    = EventManager::e()->userCandList( leptonListName() );
  CandUtil::printList( o, leptonList );

  cout << "\nelectrons"  << endl;
  CandUtil::printList( cout, _e.electrons() ); 

  cout << "\nsuper clusters"  << endl;
  CandUtil::printList( cout, _e.candList( EventManager::kEcalSC ) ); 



  for( size_t ii=0; ii<_e.electrons().size(); ii++ )
    {
      Candidate* el_ = _e.electrons()[ii];
      el_->oneLine(cout);

      for( size_t isc_=0; isc_<2; isc_++ )
	{
	  string sc_mapname_ = "electron-SuperCluster";
	  string bc_mapname_ = "SuperCluster-BasicCluster";
	  string rh_mapname_ = "BasicCluster-EcalRecHit";
	  string prefix_     = "E-Gamma";
	  if( isc_==1 )
	    {
	      sc_mapname_ = "electron-PfSuperCluster";
	      bc_mapname_ = "PfSuperCluster-PfBasicCluster";
	      rh_mapname_ = "PfBasicCluster-EcalRecHit";
	      prefix_     = "Particle-Flow";
	    }

	  const CandAssoc& sc_assoc_ = _e.candAssoc( sc_mapname_ );
	  Candidate* sc_ = sc_assoc_.getSecond( el_ );
	  if( sc_==0 )
	    {
	      cout << "Warning -- No " << prefix_ 
		   << " Super-Cluster associated: " << endl;
	      continue;
	    }
	  cout << prefix_ << " Super-Cluster" << endl;
	  cout << "\t";
	  sc_->oneLine( cout );
	  cout << "\t... composed of the following basic clusters:" << endl;
	  //GHM const CandIdEcalMap& bc_map_ = ecalMap(bc_mapname_);
	  const CandMap& bc_map_ = _e.candMap(bc_mapname_);
	  CandMapIterator bc_it_lower_ = bc_map_.lower_bound( sc_ );
	  CandMapIterator bc_it_upper_ = bc_map_.upper_bound( sc_ );
	  CandMapIterator bc_it;
	  for( bc_it=bc_it_lower_; bc_it!=bc_it_upper_; bc_it++ )
	    {	    
	      Candidate* bc_ = bc_it->second;
	      cout << "\t\t";
	      bc_->oneLine( cout );
	      cout << "\t\t... composed of the following RecHits:" << endl;
	      const CandIdRecHitMap& rh_map_ = _e.ecalRecHitMap(rh_mapname_);
	      CandIdRecHitMapIterator rh_it_lower_ = rh_map_.lower_bound( bc_->uid() );
	      CandIdRecHitMapIterator rh_it_upper_ = rh_map_.upper_bound( bc_->uid() );
	      CandIdRecHitMapIterator rh_it;
	      size_t nrh_(0);
	      for( rh_it=rh_it_lower_; rh_it!=rh_it_upper_; rh_it++ )
		{	    
		  ecalRecHit rh_   = rh_it->second.first;
		  float      frac_ = rh_it->second.second;
		  float rh_e_ = rh_.e;
		  nrh_++;
		  if( rh_e_>0.5 )
		    {
		      cout << "\t\t\t";
		      cout << "\tix/iy/iz=" 
			   << rh_.ix << "/"
			   << rh_.iy << "/"
			   << rh_.iz;
		      cout << "\te=" << rh_.e;
		      cout << "\tfrac=" << frac_;
		      cout << endl;
		    }
		} // rechits
	      if( nrh_>0 ) cout << "\t\t..." << nrh_ << " RecHits " << endl;
	    } // basic clusters		  		  
	} // e/g -> PF
      
      // tracks
      for( size_t itrk_=0; itrk_<2; itrk_++ ) {	  
	string trk_mapname_ = "electron-Track";
	if( itrk_==1 )
	  {
	    trk_mapname_ = "electron-GsfTrack";
	    cout << "GSF ";
	  }
	cout << "track(s)" << endl;
	const CandAssoc& trk_map_ = _e.candAssoc( trk_mapname_ );
	Candidate* trk_ = trk_map_.getSecond( el_ );
	cout << "\t";
	trk_->oneLine( cout );
      } // track <-> GSF track      
    }// electrons
  
  cout << "Super Clusters" << endl;
  const CandList& ecalSCList = _e.candList(EventManager::kEcalSC);
  const CandAssoc& sc_assoc_ = _e.candAssoc( "electron-SuperCluster" );
  for( size_t ii=0; ii<ecalSCList.size(); ii++ )
    {      
      Candidate* egsc = ecalSCList[ii];
      egsc->oneLine( cout );
      Candidate* el_ = sc_assoc_.getFirst( egsc );
      if( el_!=0 )
	{
	  cout << "---->";
	  el_->oneLine( cout );
	}
    }

  cout << "Tracks" << endl;
  const CandList& trackList = _e.tracks();
  const CandAssoc& trk_assoc_ = _e.candAssoc( "electron-Track" );
  const CandAssoc& gsf_assoc_ = _e.candAssoc( "electron-GsfTrack" );
  for( size_t ii=0; ii<trackList.size(); ii++ )
    {      
      Candidate* trk_ = trackList[ii];
      //      trk_->oneLine( cout );
      Candidate* el_ = trk_assoc_.getFirst( trk_ );
      if( el_!=0 )
	{
	  cout << "---->";
	  el_->oneLine( cout );
	  Candidate* gsf_ = gsf_assoc_.getSecond( el_ );
	  assert( gsf_!=0 );
	  cout << "---->";
	  gsf_->oneLine( cout );
	}
    }

  if( _e.electrons().size()!=0 )
    {
      Candidate* el = _e.electrons()[0];
      el->oneLine(cout);
      const CandList& pfList = _e.pfCandidates();
      cout << "PFCandidates ... " << pfList.size() << endl;
      const CandAssoc& pf_assoc_ = _e.candAssoc( "electron-PfCand" );
      Candidate* pf = pf_assoc_.getSecond( el );
      float drmin(0);
      float ptmin(5.);
      if( pf==0 )
	{
	  cout << "PF electron not found, find the closest electron" << endl;
	  pf = CandUtil::closest( el, pfList, drmin, ptmin, true ); 
	}
      if( pf==0 )
	{
	  cout << "PF electron not found, find the closest hadron" << endl;
	  pf = CandUtil::closest( el, pfList, drmin, ptmin, false, true ); 
	}
      if( pf!=0 )
	{
	  cout << "at dR=" << drmin << "...";
	  pf->oneLine( cout );
	}
      else
	cout << "No hard PF candidate with same charge in vicinity" << endl;
    }

}

bool
VBTFAnalysis::Z_2l__analysis()
{
  return false;
}

bool
VBTFAnalysis::W_ln__analysis()
{
  //GHM!!!
  if( _imode!=0 ) return false;

  // pseudo experiments
  if( _pseudo ) checkEvent( _ievt );

  string & prefix = modeName[kW_en]; 


  if( _verbose )
    {
      for( unsigned sel_=0; sel_<VBTFElectronSelector::kNWP; sel_++ )
	{
	  string & level_ =  VBTFElectronSelector::level[sel_];
	  CandList& list_ 
	    = _e.userCandList( leptonListName( level_.c_str() ) );
	  cout << level_ << endl;
	  CandUtil::printList( cout, list_ );
	}    
    }  
  
  CandList& leptonList_ = _e.userCandList( leptonListName("Presel") );
  if( leptonList_.size() == 0 ) return false;

  int nCand=0;
  int nPresel=0;
  //  for( unsigned ii=0; ii< leptonList_.size(); ii++ )
  for( unsigned ii=0; ii<1; ii++ )
    {
      Candidate& cand = *leptonList_[ii];
      //      cand.oneLine(cout);
      
      CandInfo* info = cand.info();
      if( !info->getBool("isEcalFiducial") ) continue;
      
      bool isEB = info->getBool("isEB");
      bool isEE = info->getBool("isEE");
      assert( isEB || isEE );
      assert( !(isEB && isEE) );

      bool isEt20 = info->getBool("isEt20");
      if( !isEt20 ) continue;

      nPresel++;

      bool isEt25 = info->getBool("isEt25");
      bool isEt30 = info->getBool("isEt30");
     
      float iso(0.), S_trk(0.), S_ecal(0.), S_hcal(0.);
      info->getFloat("iso", iso );
      info->getFloat("S_trk",  S_trk  );
      info->getFloat("S_ecal", S_ecal );
      info->getFloat("S_hcal", S_hcal );

      bool eid70   = info->getBool( "eid70" );
      bool eid80   = info->getBool( "eid80" );
      bool eid90   = info->getBool( "eid90" );
      bool eid95   = info->getBool( "eid95" );
      bool eidAnti = info->getBool( "eidAnti" );
      
      bool iso70   = info->getBool( "iso70" );
      bool iso80   = info->getBool( "iso80" );
      bool iso90   = info->getBool( "iso90" );
      bool iso95   = info->getBool( "iso95" );
      bool isoAnti = info->getBool( "isoAnti" );
      float sumTrk  = info->getFloat("sumTrk"  );
      float sumEcal = info->getFloat("sumEcal" );
      float sumHcal = info->getFloat("sumHcal" );
      
      bool sel70   = info->getBool( "sel70" );
      bool sel80   = info->getBool( "sel80" );
      bool sel90   = info->getBool( "sel90" );
      bool sel95   = info->getBool( "sel95" );
      bool selAnti = info->getBool( "selAnti" );
      
      float pt  = cand.pt();
      float eta = cand.eta();
      float phi = cand.phi();
      // charge 
      int chg = ( cand.charge()>0 ) ? 1:-1;


      float caloEta = info->getFloat("caloEta");

      float caloEnergy  = info->getFloat("caloEnergy");
      float dEtaIn      = info->getFloat("dEtaIn");
      float dPhiIn      = info->getFloat("dPhiIn");
      float sigiEtaiEta = info->getFloat("sigiEtaiEta");
      float HOE         = info->getFloat("HOE");
      float EOP         = info->getFloat("EOP");

      _tupleVar["run"]=_e.run();
      _tupleVar["event"]=_e.event();
      _tupleVar["e"]=caloEnergy;
      _tupleVar["et"]=KineUtils::et( caloEnergy, caloEta );
      _tupleVar["fid"]= ( isEB ? 1:-1 );
      {
	int isel_=0;
	if( sel70 ) isel_ = 1;
	else if( sel80 ) isel_ = 2;
	else if( sel90 ) isel_ = 3;
	else if( sel95 ) isel_ = 4;
	if( selAnti  ) isel_ += 10;
	_tupleVar["sel"]=isel_;
      }
      {
	int ieid_=0;
	if( eid70 ) ieid_ = 1;
	else if( eid80 ) ieid_ = 2;
	else if( eid90 ) ieid_ = 3;
	else if( eid95 ) ieid_ = 4;
	if( eidAnti  ) ieid_ += 10;
	_tupleVar["eid"]=ieid_;
      }
      {
	int iiso_=0;
	if( iso70 ) iiso_ = 1;
	else if( iso80 ) iiso_ = 2;
	else if( iso90 ) iiso_ = 3;
	else if( iso95 ) iiso_ = 4;
	if( isoAnti  ) iiso_ += 10;
	_tupleVar["iso"]=iiso_;
      }
      _tupleVar["pt"]=chg*pt;
      _tupleVar["eta"]=eta;
      _tupleVar["phi"]=phi;
      _tupleVar["eop"]=EOP;
      _tupleVar["hoe"]=HOE;
      _tupleVar["sigieie"]=sigiEtaiEta;
      _tupleVar["dphi"]=dPhiIn;
      _tupleVar["deta"]=dEtaIn;
      _tupleVar["strk"]=sumTrk;
      _tupleVar["secal"]=sumEcal;
      _tupleVar["shcal"]=sumHcal;
      
      // get the unit vectors along and perpendicular to LCand 
      TVector2 ux_ = cand.p2().Unit();
      TVector2 uy_( -ux_.Y(), ux_.X() ); 
      
      // pdgCode
      int LPdgCode_  = cand.pdgCode();
      int nuPdgCode_ = -1*(LPdgCode_-chg);
      
      TVector2 neuP2 = _e.neutrinosFromBosonsP4().Vect().XYvector();

      Candidate* WCand_(0);
      for( int imet=0; imet<3; imet++ )
	{
	  string suffix1;      
	  const Candidate* metCand_(0);
	  if( imet==0 )
	    {
	      suffix1 = "_calo";
	      metCand_ = _e.met("Calo");
	    }
	  else if( imet==1 )
	    {
	      suffix1 = "_paf";
	      metCand_ = _e.met("Pat");
	    }
	  else if( imet==2 )
	    {
	      suffix1 = "_pf";
	      metCand_ = _e.met("ParticleFlow");
	    }
	  else
	    continue;

	  assert( metCand_!=0 );
	  Candidate* nuCand_ = metCand_->clone();
	  nuCand_->setPdgCode( nuPdgCode_ ); 
	  TVector2 nuP2_ = nuCand_->p2();
	  Candidate* LNuCand_ = Candidate::create( &cand, nuCand_ );
	  LNuCand_->setName("LNu");
	  
	  // get the unit vectors along and perpendicular to LCand 
	  TVector2 lnux_ = LNuCand_->p2().Unit();
	  TVector2 lnuy_( -lnux_.Y(), lnux_.X() ); 
	  
	  float MET_   = nuCand_->pt();
	  float MT_    = LNuCand_->mass();
	  float dphi_  = nuCand_->p2().DeltaPhi( neuP2 ) * Constants::radToDeg;
	  float dpt_   = nuCand_->pt()-neuP2.Mod();
	  
	  if( _verbose )
	    {
	      cand.print(cout);
	      LNuCand_->oneLine(cout);
	      cout << string("Met")+suffix1 << "=" << MET_ << "  MT=" << MT_ << endl;
	      cout << "dphi/dpt " << dphi_ << "/" << dpt_ << endl;
	    }
	  
	  if( sel70 && isEt30 && imet==2 )
	    {
	      WCand_ = LNuCand_->clone();
	      WCand_->setPdgCode( chg*24 );
	    }

	  // projections
	  float lproj_x_ = -( nuP2_ * ux_ );
	  float lproj_y_ = -( nuP2_ * uy_ );

	  float wproj_x_ = -( nuP2_ * lnux_ );
	  float wproj_y_ = -( nuP2_ * lnuy_ );

	  if( imet==0 )
	    {	      
	      _tupleVar["calo_met"]=MET_;
	      _tupleVar["calo_mt"]=MT_;
	      _tupleVar["calo_dph"]=dphi_;
	      _tupleVar["calo_dpt"]=dpt_;
	      _tupleVar["calo_lproj_x"]=lproj_x_;
	      _tupleVar["calo_lproj_y"]=lproj_y_;
	      _tupleVar["calo_wproj_x"]=wproj_x_;
	      _tupleVar["calo_wproj_y"]=wproj_y_;
	    }
	  else if( imet==2 )
	    {
	      _tupleVar["pf_met"]=MET_;
	      _tupleVar["pf_mt"]=MT_;
	      _tupleVar["pf_dph"]=dphi_;
	      _tupleVar["pf_dpt"]=dpt_;
	      _tupleVar["pf_lproj_x"]=lproj_x_;
	      _tupleVar["pf_lproj_y"]=lproj_y_;
	      _tupleVar["pf_wproj_x"]=wproj_x_;
	      _tupleVar["pf_wproj_y"]=wproj_y_;
	    }
	  
	  for( int icut=-1; icut<11; icut++ )
	    {		  
	      string suffix2(suffix1);
	      if( icut==0 )
		{
		  suffix2 += "_sel95";
		  if( !sel95 ) continue;
		}
	      else if( icut==1 )
		{
		  suffix2 += "_sel90";
		  if( !sel90 ) continue;
		}
	      else if( icut==2 ) 
		{
		  suffix2 += "_sel80";
		  if( !sel80 ) continue;
		}
	      else if( icut==3 ) 
		{
		  suffix2 += "_sel70";
		  if( !sel70 ) continue;
		}
	      else if( icut==4 )
		{
		  suffix2 += "_eid80";
		  if( !eid80 ) continue;
		}
	      else if( icut==5 )
		{
		  suffix2 += "_eid80_iso70";
		  if( !(eid80 && iso70) ) continue;
		}
	      else if( icut==6 )
		{
		  suffix2 += "_eid80_iso90";
		  if( !(eid80 && iso90) ) continue;
		}
	      else if( icut==7 )
		{
		  suffix2 += "_iso80";
		  if( !iso80 ) continue;
		}
	      else if( icut==8 )
		{
		  suffix2 += "_eid70_iso80";
		  if( !(eid70 && iso80) ) continue;
		}
	      else if( icut==9 )
		{
		  suffix2 += "_eid90_iso80";
		  if( !(eid90 && iso80) ) continue;
		}
	      else if( icut==10 )
		{
		  suffix2 += "_anticut";
		  if( !selAnti ) continue;
		}
	      else
		{
		  suffix2 += "_presel";
		}

	      for( int ipt=0; ipt<4; ipt++ )
		{
		  string suffix3(suffix2);
		  if( ipt==1 )
		    {
		      suffix3+="_pt20";
		      if( !isEt20 ) continue;
		    } 
		  else if( ipt==2 )
		    {
		      suffix3+="_pt25";
		      if( !isEt25 ) continue;
		    } 
		  else if( ipt==3 )
		    {
		      suffix3+="_pt30";
		      if( !isEt30 ) continue;
		    } 
		  else
		    {
		      continue;
		    }
	  
		  for( int ichg=-1; ichg<2; ichg++ )
		    {
		      string suffix5(suffix3);
		      if( ichg==-1 )
			{
			  suffix5 += "_neg";
			  if( chg!=-1 ) continue; 
			}
		      else if( ichg==1 )
			{
			  suffix5 += "_pos";
			  if( chg!=1 )  continue; 
			}

		      for( int ireg=0; ireg<3; ireg++ )
			{
			  string suffix4(suffix5);
			  if( ireg==0 )
			    {
			      suffix4+="_EB";
			      if( !isEB ) continue;
			    }
			  else if( ireg==1 )
			    {
			      suffix4+="_EE";
			      if( !isEE ) continue;
			    }
	      	      
			  string suffix(suffix4);
			  fill( "pt"     , suffix,           pt,    prefix );
			  fill( "eta"    , suffix,           eta,    prefix );
			  fill( "met"    , suffix,           MET_,    prefix );
			  fill( "mt"     , suffix,           MT_,     prefix );
			  fill( "dphi"   , suffix,         dphi_,   prefix );
			  fill( "dpt"    , suffix,          dpt_,    prefix );
			  
			  fill( "Eproj"  , string("_lx")+suffix, lproj_x_ , prefix );
			  fill( "Eproj"  , string("_ly")+suffix, lproj_y_ , prefix );
			  fill( "Eproj2D" , string("l")+suffix, lproj_x_, lproj_y_ , prefix );
			  fill( "Eproj"  , string("_wx")+suffix, wproj_x_ , prefix );
			  fill( "Eproj"  , string("_wy")+suffix, wproj_y_ , prefix );
			  fill( "Eproj2D" , string("w")+suffix, wproj_x_, wproj_y_ , prefix );
			}
		    }
		}
	    }
	  
	  if( WCand_!=0 )
	    {
	      if( _pseudo )
		{
		  {
		    // cout << suffix << endl;
		    //		      select_ = true;
		    //		      TString toto_(str_);
		    setVarValue( "pf_mt",     MT_     );
		    setVarValue( "pf_met",    MET_    );
		    setVarValue( "pf_met_lx", lproj_x_ );
		    setVarValue( "pf_met_ly", lproj_y_ );
		    setVarValue( "pf_met_wx", wproj_x_ );
		    setVarValue( "pf_met_wy", wproj_y_ );
		    saveEvent();
		  }
		}
	      
	      if( _verbose )
		{
		  cout << "W Candidate " << endl;
		  WCand_->oneLine( cout );
		  WCand_->theBase()->oneLine( cout );	      
		}
	      nCand++;
	    }
	}

      {
	float val_[100];
	//	size_t ival_(0);
	//	map< string, float >::const_iterator it_ = _tupleVar.begin();
	//	for( ; it_!=_tupleVar.end(); it_++ )
	for( size_t ii=0; ii<_indexVar.size(); ii++ )
	  {
	    val_[ii] = _tupleVar[_indexVar[ii]];
	  }
	_theTuple->Fill( val_ );
      }      
    }
  
  if( nCand==0 ) return false;
  
  return true;
}

void
VBTFAnalysis::fillStatHistograms()
{
  // stat histograms contain exactly the number of events selected at
  // each step of the selection, irrespective of the number of hypotheses 
  // --> they can be used for statistics
  string prefix("stat");

  // Just for stat histograms, take the "best" Z candidates
  // (=first Z candidate of first Z-hypothesis, if any)
  // This should be unbiased 
  // Mass=zero means no candidate.
  //  const Candidate* bestZ[2] = {0,0};
  float mll;

  const CandList& list = _e.userCandList( ZListName() );
  if( list.size()>0 )
    {
      const Candidate* cand_ = list[0];
      assert( cand_!=0 );
      //	  bestZ[ii] = cand_;
      mll = cand_->mass();
    }

  // missing et  !!!FIXME!!!
  //  const TVector2& pafMet  = _e.met( EventManager::kPatMet  );
  //  const TVector2& caloMet = _e.met( EventManager::kCaloMet );
  //  const TVector2& pfMet   = _e.met( EventManager::kPfMet   );
  //  const TVector2& genMet  = _e.met( EventManager::kGenMet  );

  for( size_t lev=kAll; lev<kNLevels; lev++ )
    {
      ModeStatIterator it = _nMode[lev].begin();
      for( ; it!=_nMode[lev].end(); it++ )
	{
	  string modeName_ = modeName[it->first];
	  string suffix = modeName_;
	  // if the selected mode corresponds to the true mode,
	  //  this is signal. 
	  bool isSignal = ( _e.categ() == modeName_ );
	  // check if the true event is within acceptance
	  if( isSignal && !_e.inAcceptance() )
	    {
	      suffix += "_";
	      suffix += "out";		
	    }
	  suffix += "_";
	  suffix += levelName[lev];
	  // !!!FIXME!!!
	  //	    fill( "met", "calo_"+suffix, caloMet.Mod(), prefix );  
	  //	    fill( "met", "pf_"+suffix,     pfMet.Mod(), prefix );  
	  //	    fill( "met", "paf_"+suffix,   pafMet.Mod(), prefix );  
	  //	    fill( "met", "gen_"+suffix,   genMet.Mod(), prefix );  
	  for( int ii=0; ii<2; ii++ )
	    {
	      TString str_;
	      str_ += ii;
	      str_ += "_";
	      str_ += suffix;
	      fill( "mll", str_.Data(), mll, prefix );
	    }
	}
    }
}

void
VBTFAnalysis::debugPrintouts( ostream& o )
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


  static int nDbg_=0;
  for( unsigned ii=0; ii<kNModes; ii++ ) if( _categ==modeName[ii] ) 
					   {
					     nDbg_++;
					     if( nDbg_>10 ) continue;

					     //      o << "*****************" << endl;
					     //      o << "Diboson Analysis -- Event " << _ievt << " -- categ " << _categ << endl;
					     print( o );
					   }

  static int nSel_=0;
  nSel_++;
  if( nSel_<=1000 )
    {
      o << "Selected event -- ";
      o << " ievt/file/evtInFile ";
      o << _ievt << "/"; 
      o << _e.file() << "/";
      o << _e.eventInFile();
      o << "\t[" << _categ << "]";
      ModeStatIterator it = _nMode[kSel].begin();
      for( ; it!=_nMode[kSel].end(); it++ )
	{
	  o << "/" << modeName[it->first] << "=" << it->second; 
	}
      o << endl;
    }
}

void
VBTFAnalysis::buildVBTFLeptonLists()
{
  if( _imode!=kElectron ) return;
  {
    unsigned etcut_ = VBTFElectronSelector::kEt20;
    unsigned sel_;

    // prepare an ordered list of preselected electrons 
    // starting from the electron list
    sel_= VBTFElectronSelector::kPresel;		       
    const string & level_ =  VBTFElectronSelector::level[sel_];
    VBTFElectronSelector elPresel(etcut_,sel_); 
    CandList& preselElectrons 
      = _e.userCandList( leptonListName( level_.c_str() ) );
    elPresel.getList( _e.electrons(), preselElectrons );

    if( preselElectrons.size()==0 ) return;

    // sort the list
    sort( preselElectrons.begin(), preselElectrons.end(), CandCompare() );    

    // get other lists
    for( sel_++; sel_<VBTFElectronSelector::kNWP; sel_++ )
      {
	level_ =  VBTFElectronSelector::level[sel_];
	CandList& list_ 
	  = _e.userCandList( leptonListName( level_.c_str() ) );
	VBTFElectronSelector selector_(etcut_,sel_); 
	selector_.getList( preselElectrons, list_ );
      }    
  }
}

const char*
VBTFAnalysis::leptonListName( const char* level ) const
{
  string lname = lmodeName[_imode];
  lname += "_";
  lname += string(level);
  return lname.c_str();
}

const char*
VBTFAnalysis::ZListName() const
{
  string lname = "ZList";
  lname += "_presel";
  return lname.c_str();
}
