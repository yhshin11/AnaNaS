#include "Analysis/src/DibosonAnalysis.hh"
 
#include <cassert>
#include <algorithm>
#include <sstream>
#include <fstream>
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

ClassImp( DibosonAnalysis )

typedef vector<double> vectorFloat; 

string DibosonAnalysis::modeName[DibosonAnalysis::kNModes] = 
{
  "all",
  "ZZ_2e2n",   "ZZ_2m2n",   
  "ZZ_4e",     "ZZ_2e2m",     "ZZ_4m",
  "WZ_3en",    "WZ_2emn",   "WZ_e2mn",  "WZ_3mn",
  "WW_2e2n",   "WW_em2n",   "WZ_2m2n",
  "Z_2e",      "Z_2m",
  "W_en",      "W_mn"
};

string DibosonAnalysis::levelName[DibosonAnalysis::kNLevels] = 
{
  "nosel", "presel", "sel"   
};

DibosonAnalysis::DibosonAnalysis( Sample& sample, EventManager& manager ) : 
  SampleAnalysis( "Diboson", sample, manager )
{
  cout << "\t-------------------------------------------" << endl; 
  cout << "\t---  Diboson preselection analysis     ----" << endl; 
  cout << "\t-------------------------------------------" << endl; 

  _hlt = false; // apply HLT selection
  _nDebug = 0;
}

DibosonAnalysis::~DibosonAnalysis()
{
}

void
DibosonAnalysis::bookHistograms()
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
DibosonAnalysis::writeHistograms()
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
}

string 
DibosonAnalysis::ZmodeName(          int imode )
{
  assert( imode>=00 && imode<kNZmode );
  switch ( imode )
    {
    case kZeeMode:
      return "Zee";
    case kZmmMode:
      return "Zmm";
    case kZemMode:
      return "Zem";
    }    
  return "";
}

string 
DibosonAnalysis::rejectedZListName(  int imode )
{
  return ZmodeName(imode)+"Rejected";
}

string 
DibosonAnalysis::ZListName(          int imode, size_t ihyp )
{
  string ZmodeName_( ZmodeName(imode) );
  ZmodeName_ += "_";
  ZmodeName_ += ihyp;
  return ZmodeName_;
}

string 
DibosonAnalysis::leptonName(         int ilepton )
{
  assert( ilepton==0 || ilepton==1 );
  switch ( ilepton )
    {
    case 0: 
      return "Electron";
    case 1:
      return "Muon";
    }
  return 0;
}

string 
DibosonAnalysis::leptonListName(     int ilepton )
{
  return leptonName(ilepton)+"VeryLoose";
}

string 
DibosonAnalysis::leptonListName(     int ilepton, size_t ihyp )
{
  string leptonListName_( leptonListName( ilepton ) );
  leptonListName_ += "_";
  leptonListName_ += ihyp;
  return leptonListName_;
}

bool
DibosonAnalysis::analyzeEvent()
{ 
  
  // event initialisations
  for( size_t lev=kAll; lev<kNLevels; lev++ )
    {
      _hypoList[lev].clear();
      _nMode[lev].clear();
      _nHypo[lev].clear();
    }
  _nMode[kNosel][kAll]=1;

  // build lists of leptons and prepare the MC matching map
  buildLeptonLists();

  // build the combinatorics of Z candidates 
  // (Z-hypothesis = sorted lists of non-overlapping Z candidates and leptons)
  for( int imode=0; imode<3; imode++ )
    {
      nZHypos[imode] = buildZllLists( imode );
    }
 
  // preselection
  for( size_t ihyp=0; ihyp<nZHypos[0]; ihyp++ )
    {
      for( size_t jhyp=0; jhyp<nZHypos[1]; jhyp++ )
	{
	  const CandList& ZeeList = _e.userCandList( ZListName(0,ihyp)      );
	  const CandList& eList   = _e.userCandList( leptonListName(0,ihyp) );
	  const CandList& ZmmList = _e.userCandList( ZListName(1,jhyp)      );
	  const CandList& mList   = _e.userCandList( leptonListName(1,jhyp) );
	  size_t nZee = ZeeList.size();
	  size_t nZmm = ZmmList.size();
	  size_t ne   = eList.size();
	  size_t nm   = mList.size();
	  Hypo h_ = make_pair( ihyp, jhyp );
	  if( nZee==1 && nZmm==0 )                    select( h_, kZZ_2e2n, kPresel );
	  if( nZee==0 && nZmm==1 )                    select( h_, kZZ_2m2n, kPresel );
	  if( nZee>=2 && nZmm==0 )                    select( h_, kZZ_4e  , kPresel );
	  if( nZee>=1 && nZmm>=1 )                    select( h_, kZZ_2e2m, kPresel );
	  if( nZee==0 && nZmm>=2 )                    select( h_, kZZ_4m  , kPresel );
	  if( nZee==0 && nZmm==1 && ne>=1 )           select( h_, kWZ_e2mn, kPresel );
	  if( nZee==0 && nZmm==1 && nm>=1 )           select( h_, kWZ_3mn , kPresel );
	  if( nZee==1 && nZmm==0 && ne>=1 )           select( h_, kWZ_3en , kPresel );
	  if( nZee==1 && nZmm==0 && nm>=1 )           select( h_, kWZ_2emn, kPresel );
	  if( nZee==0 && nZmm==0 && nm>=1 && ne>=1 )  select( h_, kWW_em2n, kPresel );
	  if( nZee==0 && nZmm==0 && ne==0 && nm>=2 )  select( h_, kWW_2m2n, kPresel );
	  if( nZee==0 && nZmm==0 && ne>=2 && nm==0 )  select( h_, kWW_2m2n, kPresel );
          if( nZee>=1 )                               select( h_, kZ_2e,    kPresel );
          if( nZmm>=1 )                               select( h_, kZ_2m,    kPresel );
          if( nZee==0 && nZmm==0 && ne>=1 && nm>=0 )  select( h_, kW_en,    kPresel );
          if( nZee==0 && nZmm==0 && ne>=0 && nm>=1 )  select( h_, kW_mn,    kPresel );
 	}
    }

  // loop on pre-selected hypotheses
  int nSelected(0);
  int nSelZZ(0);
  int hypo;
  for( HypoListIterator ithyp=_hypoList[kPresel].begin(); 
       ithyp!=_hypoList[kPresel].end(); ithyp++ )
    {
      int kmode   = ithyp->second;

      if( kmode==kZZ_2e2n || kmode==kZZ_2m2n )
	{
	  hypo =kmode;
	  if( ZZ_2l2n__analysis( ithyp ) ) {nSelected++; nSelZZ++;}
	}
      else if( kmode==kZZ_4e || kmode==kZZ_2e2m || kmode==kZZ_4m )
	{	  
	  if( ZZ_4l__analysis( ithyp ) ) nSelected++;
	}
      else if( kmode==kZ_2e || kmode==kZ_2m )
        {
          if( Z_2l__analysis( ithyp ) ) nSelected++;
        }
//       else if( kmode==kW_en || kmode==kW_mn )
//         {
//           if( W_ln__analysis( ithyp ) ) nSelected++;
//         }
//       else if( kmode==kWW_em2n )
// 	{
// 	  if( WW_em2n__analysis( ithyp ) ) nSelected++;
// 	}
    }

  // stat histograms
  fillStatHistograms();

  // debug printouts
  debugPrintouts( cout );


  return nSelected>0;
}

void
DibosonAnalysis::select( Hypo hyp, Mode mode, Level lev )
{
  _hypoList[lev].insert( make_pair( hyp, mode ) );
  _nMode[lev][mode]++;
  _nHypo[lev][hyp]++;
}

size_t
DibosonAnalysis::buildZllLists( int imode )
{
  if( imode==kZemMode ) return buildZemList();
  
  // get the lepton list
  CandList& leptonList = _e.userCandList( leptonListName(imode) );

  const CandList& initialList = _e.compositeCandList( ZmodeName(imode) );

  // create the list of rejected candidates
  CandList& rejectedList = 
    EventManager::e()->userCandList( rejectedZListName(imode) );
  assert( rejectedList.size()==0 );

  // the list of selected Z candidates
  CandList list;
  for( size_t iZ=0; iZ<initialList.size(); iZ++ )
    {
      Candidate* cand = initialList[iZ];

      // select candidates in the Z0 window 
      //  with at least two Loose lepton candidates
      //      if( selectZCand( *cand, Z0_window[0], Z0_window[1], 0, 0, 2, 2 ) )
      if( selectZCand( *cand, Z0_fullDY[0], Z0_fullDY[1], 0, 0, 0, 2 ) ) // !!!! to be checked !!!
	{
	  list.push_back( cand );
	}
      else
	{
	  rejectedList.push_back( cand );
	}
    }

  // is the list is empty, there is only one trivial hypothesis
  if( list.size()==0 ) 
    {
      CandList& ZList_ = _e.userCandList( ZListName(imode,0) );
      ZList_ = list;
      CandList& leptonList_ = _e.userCandList( leptonListName(imode,0) );
      leptonList_ = leptonList;      
      return 1;
    }

  // sort and dispatch the list
  sort( list.begin(), list.end(), ZllCompare() );
  vector< CandList > Zhypos;
  CandUtil::dispatchList( list, Zhypos );
  size_t nZhypos_ = Zhypos.size();

  // create the hypotheses
  for( size_t ihyp=0; ihyp<nZhypos_; ihyp++ )
    {
      CandList& ZList_ = _e.userCandList( ZListName(imode,ihyp) );
      ZList_ = Zhypos[ihyp];
      CandList& leptonList_ = _e.userCandList( leptonListName(imode,ihyp) );
      CandUtil::pruneList( leptonList, leptonList_, ZList_ );
    }
  
  // return the number of hypotheses
  return nZhypos_;
}

size_t
DibosonAnalysis::buildZemList()
{
  // get the lepton list
  CandList& electronList = _e.userCandList( leptonListName(0) );
  CandList&     muonList = _e.userCandList( leptonListName(1) );

  CandList list;
  for( size_t ie=0; ie<electronList.size(); ie++ )
    {
      Candidate* electron = electronList[ie];
      for( size_t im=0; im<muonList.size(); im++ )
	{
	  Candidate* muon = muonList[im];

	  // create an electron-muon candidate
	  Candidate* cand = Candidate::create( electron, muon );
	  
	  if( selectZCand( *cand, Z0_fullDY[0], Z0_fullDY[1], 0, 0, 2, 2 ) )
	    {
	      list.push_back( cand );
	    }
	}
    }

  // sort and dispatch the list
  sort( list.begin(), list.end(), ZllCompare() );
  CandList& ZList_ = _e.userCandList( "Zem" );
  ZList_ = list;
  return list.size();
}

void
DibosonAnalysis::print( ostream& o, HypoStatIterator it, Level lev ) const
{
  Hypo hyp = it->first;
  size_t ihyp = hyp.first;
  size_t jhyp = hyp.second;
  o << "======================" << endl;
  o << "Event hypothesis [" << ihyp << "," << jhyp << "] " << endl;
  HypoListIterator it_lower = _hypoList[lev].lower_bound( hyp );
  HypoListIterator it_upper = _hypoList[lev].upper_bound( hyp );
  HypoListIterator it_;
  if( it_==it_upper ) 
    {
      o << "Not selected at level " << lev << endl;
      return; 
    }
  for( it_=it_lower; it_!=it_upper; it_++ )
    {
      assert( hyp==it_->first  );
      Mode kmode = it_->second;
      o << "[" << ihyp << "," << jhyp << "] --> " << modeName[kmode] << endl;
    }
  const CandList& ZeeList = EventManager::e()->userCandList( ZListName(0,ihyp) );
  const CandList& ZmmList = EventManager::e()->userCandList( ZListName(1,jhyp) );
  const CandList& ZemList = EventManager::e()->userCandList( ZListName(2,jhyp) );
  const CandList& eList   = EventManager::e()->userCandList( leptonListName(0,ihyp) );
  const CandList& mList   = EventManager::e()->userCandList( leptonListName(1,jhyp) );
  CandUtil::printList( o, ZeeList, ZmodeName(0).c_str() );
  CandUtil::printList( o, ZmmList, ZmodeName(1).c_str() );
  CandUtil::printList( o, ZemList, ZmodeName(2).c_str() );
  CandUtil::printList( o, eList, leptonName(0).c_str() );
  CandUtil::printList( o, mList, leptonName(1).c_str() );
}

void
DibosonAnalysis::print( ostream& o ) const
{
  o << "*****************" << endl;
  o << "Diboson Analysis -- Event " << _ievt << " -- categ " << _e.categ() << endl;
  for( int imode=0; imode<2; imode++ )
    {
      o << "-+-+-+-+-+" << endl;
      if( nZHypos[imode]==0 )
	{
	  o << "No selected " << ZmodeName(imode) << " candidate" << endl;
	  const CandList& leptonList = EventManager::e()->userCandList( leptonListName(imode) );
	  CandUtil::printList( o, leptonList, leptonName(imode).c_str() );
	}
      else
	{
	  o << "Number of " << ZmodeName(imode) << " hypotheses = " 
	    << nZHypos[imode] <<  endl;
	  for( size_t ihyp=0; ihyp<nZHypos[imode]; ihyp++ )
	    {
	      o << "++++ " ;
	      o << ZmodeName(imode) << " hypothesis " << ihyp << endl;
	      
	      const CandList& ZList 
		= EventManager::e()->userCandList( ZListName(imode,ihyp) );
	      CandUtil::printList( o, ZList, ZmodeName(imode).c_str() );
	      
	      const CandList& leptonList 
		= EventManager::e()->userCandList( leptonListName(imode,ihyp) );
	      CandUtil::printList( o, leptonList, leptonName(imode).c_str() );
	    }
	}
      o << "-+-+-+-+-+" << endl;
      const CandList& rejectedZList 
	= EventManager::e()->userCandList( rejectedZListName(imode) );
      CandUtil::printList( o, rejectedZList, (string("rejected ")+ZmodeName(imode)).c_str() );
    }
  o << "-+-+-+-+-+" << endl;
  //  o << "Diboson Analysis -- End of event " << _ievt << endl;
}

bool
DibosonAnalysis::ZZ_4l__analysis(  HypoListIterator ithyp )
{
  int run   = _e.run();
  int event = _e.event();

  // declare the main analysis variables
  int   nZCand(0);        // number of Z candidates
  //  int nVertex(0);         // number of good vertices

  float mll_1(0);           // di-lepton mass of first pair
  float qTll_1(0);          // di-lepton transverse momentum
  float yll_1(0);           // di-lepton rapidity
  float hll_1(0);           // di-lepton helicity

  float mll_0(0);           // di-lepton mass of first pair
  float qTll_0(0);          // di-lepton transverse momentum
  float yll_0(0);           // di-lepton rapidity
  float hll_0(0);           // di-lepton helicity

  float mTZZ(0);

  int nLoose(0);
  //
  // Cuts
  //

  Hypo hyp    = ithyp->first;
  size_t ihyp = hyp.first;
  size_t jhyp = hyp.second;
  Mode kmode   = ithyp->second;

  string categ = _e.categ();
  string* categ_ptr = & categ;

  string prefix(modeName[kmode]); 

  CandList ZList[2] = 
    { EventManager::e()->userCandList( ZListName(0,ihyp) ),
      EventManager::e()->userCandList( ZListName(1,jhyp) )  };
  CandList leptonList[2] = 
    { EventManager::e()->userCandList( leptonListName(0,ihyp) ),
      EventManager::e()->userCandList( leptonListName(1,jhyp) ) };
  
  // require no other "loose" lepton in the rest of the event
  for( size_t ilep=0; ilep<2; ilep++ )
    {
      for( size_t ii=0; ii<leptonList[ilep].size(); ii++ )
	{
	  CandInfo* info = leptonList[ilep][ii]->info();
	  if( info->getBool("Loose") ) nLoose++;       	 
	}
    }
  //  fill( "n10", "loose", nLoose, prefix );
  //  if( nLoose>0 ) return false;
  
  // require two and only two Z candidate
  size_t nZ = ZList[0].size() + ZList[1].size();
  //  fill( "n10", "Z", nZ, prefix );
  if( nZ!=2 ) return false;

  Candidate* ZCand[2];
  size_t iZ(0);
  for( size_t imode=0; imode<2; imode++ )
    {
      for( size_t ii=0; ii<ZList[imode].size(); ii++ )
	{
	  ZCand[iZ++] = ZList[imode][ii];
	}
    }
  map< Candidate*, Candidate* > other;
  other[ ZCand[0] ] = ZCand[1];
  other[ ZCand[1] ] = ZCand[0];

  // require that one of the Z boson contains at least one Tight lepton 
  // and two Loose leptons and lives in the veto window
  int nTightZ(0);
  for( size_t imode=0; imode<2; imode++ )
    {
      for( size_t ii=0; ii<ZList[imode].size(); ii++ )
	{
	  Candidate* cand1 = ZList[imode][ii];
	  Candidate* cand2 = other[cand1];
	  
	  if( selectZCand( *cand2, Z0_veto[0], Z0_veto[1], 0, 1, 2, 2 ) ) 
	    {
	      TString str_;
	      str_ += imode;
	      nTightZ++;
	      // plot the other Z boson mass
	      //	      fill( "mll", str_.Data(), cand1->mass(), prefix );
	    }
	}
    }

  //  fill( "n10", "tightZ", nTightZ, prefix );
  if( nTightZ==0 ) return false;

  Candidate* ZZCand = Candidate::create( ZCand[0], ZCand[1] );
  mTZZ = ZZCand->mass();

  mll_0 = ZCand[0]->mass();
  qTll_0 = ZCand[0]->pt();
  CandInfo* ZInfo_0 = ZCand[0]->info();
  yll_0 = ZInfo_0->getFloat("rapidity");
  hll_0 = ZInfo_0->getFloat("helicity");

  mll_1 = ZCand[1]->mass();
  qTll_1 = ZCand[1]->pt();
  CandInfo* ZInfo_1 = ZCand[1]->info();
  yll_1 = ZInfo_1->getFloat("rapidity");
  hll_1 = ZInfo_1->getFloat("helicity");


  {
    // tuple
    //    string dir_ = _e.categ();
    string dir_ = prefix;
    tm.setTree( "ZZ4ltuple", dir_ );
    tm.add< int >(    "run",   &run       );
    tm.add< int >(    "event", &event     );
    tm.add< string*>( "categ", &categ_ptr );
    tm.add< int >(  "nZCand",   &nZCand        );
    tm.add< float >(  "mll_0",   &mll_0        );
    tm.add< float >(  "qTll_0",   &qTll_0        );
    tm.add< float >(  "yll_0",   &yll_0        );
    tm.add< float >(  "hll_0",   &hll_0        );
    tm.add< float >(  "mll_1",   &mll_1        );
    tm.add< float >(  "qTll_1",   &qTll_1        );
    tm.add< float >(  "yll_1",   &yll_1        );
    tm.add< float >(  "hll_1",   &hll_1        );
    tm.add< float >(  "mTZZ",   &mTZZ      );
    tm.add< int >(  "nLoose",&nLoose    );
    tm.flush();
  }    


  //  hypothesis is selected
  select( hyp, kmode, kSel );

  return true;
}

bool
DibosonAnalysis::ZZ_2l2n__analysis(  HypoListIterator ithyp )
{
  int run   = _e.run();
  int event = _e.event();

  // declare the main analysis variables
  int   nZCand(0);        // number of Z candidates
  int nVertex(0);         // number of good vertices
  float mll(0);           // di-lepton mass
  float qTll(0);          // di-lepton transverse momentum
  float phill(0);         // di-lepton phi
  float yll(0);           // di-lepton rapidity
  float hll(0);           // di-lepton helicity
  float mTZZ(0);          // Z+MET transverse mass
  float MET(0);           // PF-MET type-I
  float METPhi(0);        // PF-MET type-I - phi
  float sumEt(0);         // sum of transverse energies
  float mEtSig(0);        // sigma of MET
  float projMET(0);       // projected MET
  float corProjMET(0);    // PU-vertex corrected MET
  float sigMET(0);        // corrected MET significance
  float dPhiMin(0);       // minimum angle of MET to a lepton
  float balZMET;          // "balance" Z/MET
  bool  isTight(false);   // is the Z candidate tight?
  int   nLoose(0);        // number of additional loose leptons
  int   jMult40(0);       // number of jets above 40 GeV-Et
  int   jMult35(0);       // number of jets above 35 GeV-Et
  int   jMult30(0);       // number of jets above 30 GeV-Et
  int   jMult25(0);       // number of jets above 25 GeV-Et
  int   jMult20(0);       // number of jets above 25 GeV-Et
  int   jMult15(0);       // number of jets above 15 GeV-Et
  int   jMult10(0);       // number of jets above 10 GeV-Et
  float dPhiJMET(0);      // angle between leading jet and MET

  //MM
  float cpuMETVtx(0);  
  float cpuMETSum(0);  
  float cpuMETSqrt(0);  

  // create the jet vectors to be stored in the ntuple -- MUST be destroyed!
  // (warning: no return before the end)
  // Note to myself: implement garbage collection for these...
  vectorFloat* jEt        = new vectorFloat;
  vectorFloat* jEta       = new vectorFloat;
  vectorFloat* jPhi       = new vectorFloat;
  vectorFloat* jDR1       = new vectorFloat;
  vectorFloat* jDR2       = new vectorFloat;
  vectorFloat* jDPhiMET   = new vectorFloat;
  vectorFloat* jBtagTkC   = new vectorFloat;
  vectorFloat* jBtagSoftM = new vectorFloat;
  vector<int>* jVtx       = new vector<int>;
 
  //Photon vectors
  vectorFloat* pEt        = new vectorFloat;
  vectorFloat* pEta       = new vectorFloat;
  vectorFloat* pPhi       = new vectorFloat;
  vectorFloat* pDR1       = new vectorFloat;
  vectorFloat* pDR2       = new vectorFloat;
  vectorFloat* pDPhiMET   = new vectorFloat;

  // leptons -- the first two are from the Z
  vectorFloat* LPt   = new vectorFloat;
  vector<int>* LPdg = new vector<int>;
  vectorFloat* LEta  = new vectorFloat;
  vectorFloat* LPhi  = new vectorFloat;
  vectorFloat* LMT   = new vectorFloat;
  vector<int>* LID   = new vector<int>;
  vectorFloat* LTrkIso   = new vectorFloat;
  vectorFloat* LEcalIso  = new vectorFloat;
  vectorFloat* LHcalIso  = new vectorFloat;
  vectorFloat* LCombIso   = new vectorFloat;
  vector<int>* LVtx   = new vector<int>;

  float rhoFJ = _e.getRhoFastJet();


  //
  // Cuts
  //

  // Z quality
  OneSidedCut<int>   c_nZCand    ( "nZCand",      nZCand,    1,   "==" );
  TwoSidedCut<float> c_ZWindow   ( "ZWindow",     mll, 60., 120., "[]" );
  OneSidedCut<bool>  c_ZLQual    ( "ZLQual",     isTight, true,   "==" );
  CompositeCut        c_ZQual     ( "ZQual", "&&" );
  c_ZQual.add( c_nZCand );
  c_ZQual.add( c_ZWindow );
  c_ZQual.add( c_ZLQual );

  // MET cuts
  OneSidedCut<float> c_balZMET   ( "balZMET",    balZMET, 2.,     "<=" );
  OneSidedCut<float> c_MET0p50   ( "MET0p50",    sigMET,  0.50,   ">"  );
  OneSidedCut<float> c_MET1p00   ( "MET1p00",    sigMET,  1.00,   ">"  );
  OneSidedCut<float> c_MET1p50   ( "MET1p50",    sigMET,  1.50,   ">"  );
  OneSidedCut<float> c_MET2p00   ( "MET2p00",    sigMET,  2.00,   ">"  );
  OneSidedCut<float> c_MET2p50   ( "MET2p50",    sigMET,  2.50,   ">"  );
  OneSidedCut<float> c_MET3p00   ( "MET3p00",    sigMET,  3.00,   ">"  );
  OneSidedCut<float> c_MET3p50   ( "MET3p50",    sigMET,  3.50,   ">"  );
  OneSidedCut<float> c_MET4p00   ( "MET4p00",    sigMET,  4.00,   ">"  );
  OneSidedCut<float> c_MET4p20   ( "MET4p20",    sigMET,  4.20,   ">"  );
  CompositeCut       c_MET       ( "MET", "&&" );
  c_MET.add( c_balZMET );
  c_MET.add( c_MET0p50 );
  c_MET.add( c_MET1p00 );
  c_MET.add( c_MET1p50 );
  c_MET.add( c_MET2p00 );
  c_MET.add( c_MET2p50 );
  c_MET.add( c_MET3p00 );
  c_MET.add( c_MET3p50 );
  c_MET.add( c_MET4p00 );
  c_MET.add( c_MET4p20 );

  // jet veto
  OneSidedCut<int>   c_j30Veto   ( "j30Veto",    jMult30, 0,      "==" );
  OneSidedCut<int>   c_j25Veto   ( "j25Veto",    jMult25, 0,      "==" );
  OneSidedCut<int>   c_j20Veto   ( "j20Veto",    jMult20, 0,      "==" );
  OneSidedCut<int>   c_j15Veto   ( "j15Veto",    jMult15, 0,      "==" );
  TwoSidedCut<float> c_dPhiJMET  ( "dPhiJMET", dPhiJMET, -20.,20, "![]");
  CompositeCut c_j20Veto_OR_dPhiMET( "j20Veto_OR_dPhiMET", "||" );
  c_j20Veto_OR_dPhiMET.add( c_j20Veto  );
  c_j20Veto_OR_dPhiMET.add( c_dPhiJMET );  
  CompositeCut c_jVeto( "jVeto", "&&" );
  c_jVeto.add( c_j25Veto );
  c_jVeto.add( c_j20Veto_OR_dPhiMET );
  
  // diboson separation
  OneSidedCut<int>   c_LVeto     ( "LVeto",       nLoose, 0,      "==" );
  CompositeCut c_dibosonVeto( "dibosonSeparation", "&&" );
  c_dibosonVeto.add( c_LVeto );

  // final selection
  CompositeCut c_finalSelection( "finalSelection", "&&" );
  c_finalSelection.add( c_ZQual       );
  c_finalSelection.add( c_jVeto       );
  c_finalSelection.add( c_MET         );
  c_finalSelection.add( c_dibosonVeto );

  TBits   c_bits;
  TBits*  c_bits_ptr = &c_bits;
  string  c_str;
  string* c_str_ptr  = &c_str;
  c_finalSelection.getString(c_str);
  string categ = _e.categ();
  string* categ_ptr = & categ;

  //  cout << "categ=" << categ << " ptr=" << categ_ptr << endl;

  // constants
  float sigma_0  = 7.34;
  float sigma_PU = 3.9;
  float sigma_T(0);
  float alpha = 0.856;

  // the hypothesis
  Hypo hyp    = ithyp->first;
  size_t ihyp = hyp.first;
  size_t jhyp = hyp.second;
  Mode kmode   = ithyp->second;

  // the lepton mode for this hypothesis
  int imode(0);
  if( kmode==kZZ_2e2n ) imode=0;
  else if( kmode==kZZ_2m2n ) imode=1;
  else assert(0);
  string prefix(modeName[kmode]); 

  // The Z boson:take the best (=first) in the list
  CandList ZList[2] = 
    { EventManager::e()->userCandList( ZListName(0,ihyp) ),
      EventManager::e()->userCandList( ZListName(1,jhyp) )  };    
  Candidate* ZCand;
  assert( ZList[imode].size()>0 );
  nZCand = ZList[imode].size();
  ZCand = ZList[imode][0];
  assert( ZCand!=0 );

  // get the vertex of the Z and make it the primary vertex
  Vertex* ZVertex = ZCand->vertex();
  assert( ZVertex!=0 );
  _e.setPrimaryVertex( ZVertex );

  // the di-lepton invariant mass
  mll  = ZCand->mass();
  qTll = ZCand->pt();
  phill = ZCand->phi();
  CandInfo* ZInfo = ZCand->info();
  yll = ZInfo->getFloat("rapidity");
  hll = ZInfo->getFloat("helicity");

  // a Z boson candidate is tight if it is made out of at least one Tight 
  // and one Loose lepton
  isTight =  selectZCand( *ZCand, Z0_window[0], Z0_window[1], 0, 1, 2, 2 ); 

  // the missing transverse momentum 
  const Candidate* metCand =  _e.met( EventManager::kPfMet );
  assert( metCand!=0 );
  MET = metCand->pt();
  METPhi = metCand->phi();
  CandInfo* metInfo = metCand->info();
  sumEt = metInfo->getFloat( "sumEt" );
  mEtSig = metInfo->getFloat( "mEtSig" );

  Candidate* ZnnCand = metCand->clone();

  Candidate* ZZCand = Candidate::create( ZCand, ZnnCand );
  mTZZ = ZZCand->mass();
  
  // the leptons in the event (two of them come from the Z boson)  
  CandList leptonList[2] = 
    { EventManager::e()->userCandList( leptonListName(0) ),
      EventManager::e()->userCandList( leptonListName(1) ) };

  // - count the number of "loose" leptons in the rest of the event
  // - get the isolation of the two leptons from the Z candidate
  // - determine the minimal angle with the missing pT in the transverse plane
  nLoose = 0;
  //  vector< float > iso;
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
	      //	      float iso_(0);
	      //	      info->getFloat( "iso",iso_ );
	      //	      iso.push_back( iso_ );
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
      LPt  -> push_back( LCand_->charge() * LCand_->pt()  );
      LPdg -> push_back( pdgCode_ );
      LEta -> push_back( LCand_->eta() );
      LPhi -> push_back( LCand_->phi() );
      LMT  -> push_back( mT(LCand_,ZnnCand) );
      LVtx -> push_back( vtxIndex_ );
      int id_(0);
      CandInfo* info = LCand_->info();
      if( info->getBool("VeryLoose") ) id_ += 1;
      if( info->getBool("Loose") )     id_ += 2;
      if( info->getBool("Tight") )     id_ += 4;
      if( info->getBool("VeryTight") ) id_ += 8;
      LID  -> push_back( id_ );

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

      LTrkIso->push_back( trkIso_ );
      LEcalIso->push_back( ecalIso_ );
      LHcalIso->push_back( hcalIso_ );
      LCombIso->push_back( combIso_ );
      
    }

  dPhiMin = 1000;
  for( size_t ii=0; ii<2; ii++ )
    {
      Candidate* LCand_ = LCand[ii];
      float dPhi_ = CandUtil::dPhi( ZnnCand, LCand_ ); 
      if( fabs( dPhi_ )<dPhiMin ) dPhiMin = fabs(dPhi_);
    }
  assert( dPhiMin<=Constants::pi );  

  // the projected MET 
  projMET = MET;
  if( dPhiMin < Constants::pi/2. )
    {
      projMET *= sin( dPhiMin );
    }

  // count the number of "good" vertices
  const VtxList& vertices = _e.vertices(); 
  bool findPV=false;
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
      
  // the corrected MET
  sigma_T    = sqrt( pow(sigma_0,2) + ( nVertex-1 )*pow(sigma_PU,2) );
  corProjMET = ( projMET - alpha * ( nVertex-1 ) ) * sigma_0 / sigma_T;
  sigMET     = corProjMET / sigma_0;

  balZMET = 0;
  //  if( fabs(qTll)>0 ) balZMET    = corProjMET / qTll; 
  if( fabs(qTll)>0 ) balZMET    = MET / qTll;  // !!! change GHM 29/6


  //MM : Photon list, needed for Zgamma control
  const CandList& photons = _e.photons();
  CandList SelPh;
  for( size_t ip=0; ip<photons.size(); ip++ ) 
    { 
      Candidate* ph = photons[ip];
          
      bool ismatch=false;
      for( list< Candidate* >::iterator it=LCandAll.begin();
	   it!=LCandAll.end(); ++it )
	{
	  Candidate* LCand_ = *it;
	  float dpt_ = (ph->pt()-LCand_->pt())/LCand_->pt();
	  if( ph->dR( LCand_ ) <0.2 && dpt_>-0.2 && dpt_<0.5 )
	    { ismatch=true; break;}
	}
      if(ismatch) { continue; }
      
      float pEta_ = ph->eta();
      float pEt_  = ph->Et();
      float pPhi_ = ph->phi();
      float pDPhiMET_ = CandUtil::dPhi( ph, ZnnCand )*Constants::radToDeg ;
      float pDR1_ = CandUtil::dPhi( ph, LCand[0] )*Constants::radToDeg ;
      float pDR2_ = CandUtil::dPhi( ph, LCand[1] )*Constants::radToDeg ;
    
      if( !isPhotonIDIso(ph) ) continue;

      SelPh.push_back( ph );

      if( fabs(pEta_)<=2.5 )
	if( pEt_>=20 )
	{	  
	  pPhi    ->push_back(pPhi_ );
	  pEta    ->push_back(pEta_ );
	  pEt     ->push_back(pEt_ );
	  pDR1    ->push_back(pDR1_ ); 
	  pDR2    ->push_back(pDR2_ ); 
	  pDPhiMET->push_back(pDPhiMET_ ); 
	  
	}
    }	  
	  
  // NEW GHM (April 2012) : use list of clean jets
  jMult40=0;
  jMult35=0;
  jMult30=0;
  jMult25=0;
  jMult20=0;
  jMult15=0;
  jMult10=0;
  CandList jList10; 
  CandList jList15; 
  CandList jList20; 
  CandList jList25; 
  CandList jList30; 
  CandList jList35; 
  CandList jList40; 
  CandList& cleanJets = _e.userCandList("cleanJets");
  for( size_t ij=0; ij<cleanJets.size(); ij++ ) 
    { 
      Candidate* jet = cleanJets[ij];

      float jEta_ = jet->eta();
      float jEt_  = jet->Et();
      float jPhi_ = jet->phi();
      float jDPhiMET_ = CandUtil::dPhi( jet, ZnnCand )*Constants::radToDeg ;
      float jDR1_ = CandUtil::dPhi( jet, LCand[0] )*Constants::radToDeg ;
      float jDR2_ = CandUtil::dPhi( jet, LCand[1] )*Constants::radToDeg ;
      CandInfo* info_ = jet->info();
      float jBtagTkC_   = info_->getFloat("btagTkC");
      float jBtagSoftM_ = info_->getFloat("btagSoftM");
      int vtxIndex_ = jet->vertexIndex();


      if( fabs(jEta_)<=5 )
	if( jEt_>=10 )
	{	  
	  jPhi    ->push_back(jPhi_ );
	  jEta    ->push_back(jEta_ );
	  jEt     ->push_back(jEt_ );
	  jDR1    ->push_back(jDR1_ ); 
	  jDR2    ->push_back(jDR2_ ); 
	  jDPhiMET->push_back(jDPhiMET_ ); 
	  jBtagTkC->push_back(jBtagTkC_);
	  jBtagSoftM->push_back(jBtagSoftM_);
	  jVtx    ->push_back(vtxIndex_);
	  	  
	  jList10.push_back( jet );
	  jMult10++;
	  if( jEt_>=15 ) 
	    {
	      jMult15++;
	      jList15.push_back( jet );
	    }
	  if( jEt_>=20 ) 
	    {
	      jMult20++;
	      jList20.push_back( jet );
	    }
	  if( jEt_>=25 ) 
	    {
	      jMult25++;
	      jList25.push_back( jet );
	    }
	  if( jEt_>=30 ) 
	    {
	      jMult30++;
	      jList30.push_back( jet );
	    }
	  if( jEt_>=35 ) 
	    {
	      jMult35++;
	      jList35.push_back( jet );
	    }
	  if( jEt_>=40 ) 
	    {
	      jMult40++;
	      jList40.push_back( jet );
	    }
	}      
    }

//   // jets: Particle flow with energy scale correction and lepton cleaning
//   //MM : added second leptonphoton cleaning (PF2PAT :s )
//   jMult40=0;
//   jMult35=0;
//   jMult30=0;
//   jMult25=0;
//   jMult20=0;
//   jMult15=0;
//   jMult10=0;
//   const CandList& jets = _e.jetList( EventManager::kPatJet );
//   CandList jList10; 
//   CandList jList15; 
//   CandList jList20; 
//   CandList jList25; 
//   CandList jList30; 
//   CandList jList35; 
//   CandList jList40; 

//   for( size_t ij=0; ij<jets.size(); ij++ ) 
//     { 
//       Candidate* jet = jets[ij];
//       Vertex* vertex = jet->vertex();
//       int vtxIndex_ = -1;
      
//       bool ismatch=false;
//       for( list< Candidate* >::iterator it=LCandAll.begin();
// 	   it!=LCandAll.end(); ++it )
// 	{
// 	  Candidate* LCand_ = *it;
// 	  float dpt_ = (jet->pt()-LCand_->pt())/LCand_->pt();
// 	  if( jet->dR( LCand_ ) <0.2 && dpt_>0 && dpt_<0.5 )
// 	    { ismatch=true; break;}
// 	}
//       if(ismatch) { continue; }
      
//       //now photons
//       for(size_t ip=0;ip<SelPh.size(); ip++) {
// 	Candidate* pCand_ = SelPh[ip];
// 	float dpt_ = (jet->pt()-pCand_->pt())/pCand_->pt();
// 	if( jet->dR( pCand_ ) <0.2 && dpt_>0 && dpt_<0.5 )
// 	    { ismatch=true; break;}

//       }

//       if( vertex!=0 )
// 	{
// 	  CandInfo* vtxInfo_ = vertex->info();
// 	  if( vtxInfo_ )
// 	    {
// 	      vtxIndex_ = vtxInfo_->getInt("index");
// 	    }
// 	}
//       float jEta_ = jet->eta();
//       float jEt_  = jet->Et();
//       float jPhi_ = jet->phi();
//       float jDPhiMET_ = CandUtil::dPhi( jet, ZnnCand )*Constants::radToDeg ;
//       float jDR1_ = CandUtil::dPhi( jet, LCand[0] )*Constants::radToDeg ;
//       float jDR2_ = CandUtil::dPhi( jet, LCand[1] )*Constants::radToDeg ;
//       CandInfo* info_ = jet->info();
//       float jBtagTkC_   = info_->getFloat("btagTkC");
//       float jBtagSoftM_ = info_->getFloat("btagSoftM");
      
//       if( fabs(jEta_)<=5 )
// 	if( jEt_>=10 )
// 	{	  
// 	  jPhi    ->push_back(jPhi_ );
// 	  jEta    ->push_back(jEta_ );
// 	  jEt     ->push_back(jEt_ );
// 	  jDR1    ->push_back(jDR1_ ); 
// 	  jDR2    ->push_back(jDR2_ ); 
// 	  jDPhiMET->push_back(jDPhiMET_ ); 
// 	  jBtagTkC->push_back(jBtagTkC_);
// 	  jBtagSoftM->push_back(jBtagSoftM_);
// 	  jVtx    ->push_back(vtxIndex_);
	  	  
// 	  jList10.push_back( jet );
// 	  jMult10++;
// 	  if( jEt_>=15 ) 
// 	    {
// 	      jMult15++;
// 	      jList15.push_back( jet );
// 	    }
// 	  if( jEt_>=20 ) 
// 	    {
// 	      jMult20++;
// 	      jList20.push_back( jet );
// 	    }
// 	  if( jEt_>=25 ) 
// 	    {
// 	      jMult25++;
// 	      jList25.push_back( jet );
// 	    }
// 	  if( jEt_>=30 ) 
// 	    {
// 	      jMult30++;
// 	      jList30.push_back( jet );
// 	    }
// 	  if( jEt_>=35 ) 
// 	    {
// 	      jMult35++;
// 	      jList35.push_back( jet );
// 	    }
// 	  if( jEt_>=40 ) 
// 	    {
// 	      jMult40++;
// 	      jList40.push_back( jet );
// 	    }
// 	}      
//     }
  sort( jList10.begin(), jList10.end(), Candidate::SortByPt() );
  sort( jList15.begin(), jList15.end(), Candidate::SortByPt() );
  sort( jList20.begin(), jList20.end(), Candidate::SortByPt() );
  sort( jList25.begin(), jList25.end(), Candidate::SortByPt() );
  sort( jList30.begin(), jList30.end(), Candidate::SortByPt() );
  sort( jList35.begin(), jList35.end(), Candidate::SortByPt() );
  sort( jList40.begin(), jList40.end(), Candidate::SortByPt() );
  
  
  //The cleaned MET //MM
  Vertex* vtx = ZCand->daughter(0)->vertex();
  const Candidate* cpuMETVtxC  = _e.cleanPUmet(vtx , ZCand, nVertex );
  const Candidate* cpuMETSumC  = _e.cleanPUmet(vtx , ZCand, sumEt-ZCand->pt());
  const Candidate* cpuMETSqrtC = _e.cleanPUmet(vtx , ZCand, sqrt(sumEt-ZCand->pt()) );
  cpuMETVtx = cpuMETVtxC->pt();
  cpuMETSum = cpuMETSumC->pt();
  cpuMETSqrt = cpuMETSqrtC->pt();

  
  // the "jet-recoil" candidate
  float JZB15(-1000.);
  float JZB20(-1000.);
  float JZB25(-1000.);
  float qT = ZCand->pt();

  if( jList15.size()>0 )
    {
      Candidate* jRecoil15   = Candidate::create( jList15 );
      JZB15 = jRecoil15->pt() - qT;
    }
  if( jList20.size()>0 )
    {
      Candidate* jRecoil20 = Candidate::create( jList20 );
      JZB20 = jRecoil20->pt() - qT;
    }
  if( jList25.size()>0 )
    {
      Candidate* jRecoil25 = Candidate::create( jList25 );
      JZB25 = jRecoil25->pt() - qT;
    }

  dPhiJMET = -190.;
  if( jList10.size()>0 ) dPhiJMET = (*jDPhiMET)[0];

  // selection cuts
  for( size_t ic=0; ic<=c_finalSelection.size(); ic++ )
    {
      string name_ = c_finalSelection.name(ic);
      for( int ii=0; ii<2; ii++ )
	{ 
	  string suffix_ = name_;
	  bool ok(true);
	  if( ii==0 )
	    {
	      if( ic==c_finalSelection.size() ) continue;	  
	      suffix_ += "/Nm1";
	      ok = c_finalSelection.NMinusOne( ic );
	    }
	  else
	    {
	      if( ic!=c_finalSelection.size() )
		suffix_ += "/Seq";
	      ok = c_finalSelection.cutUpTo( ic );
	    }
	  if( ok )
	    {
	      string dir_ = prefix;
	      dir_ += "/";
	      dir_ += suffix_;
	      fill( "mll",     "",        mll,        dir_ );
	      fill( "met",     "MET",     MET,        dir_ );
	      fill( "met",     "projMET", projMET,    dir_ );
	      fill( "met",     "corProj", corProjMET, dir_ );
	      fill( "balance", "ZMET",    balZMET,     dir_ );
	      fill( "sigMET",  "",        sigMET,     dir_ );
	      fill( "jmult",   "25GeV", jMult25,    dir_ );
	      fill( "jmult",   "20GeV", jMult20,    dir_ );
	      fill( "jmult",   "15GeV", jMult15,    dir_ );
	      fill( "dphi",    "jMET", dPhiJMET, dir_ );
	      fill( "dphi",    "LMET", dPhiMin,  dir_ );
	    }
	}
    }

  bool isSelected = c_finalSelection(); 


  if( isSelected )
    {
      //  hypothesis is selected
      select( hyp, kmode, kSel );
    }

  c_finalSelection.getBits(c_bits);
  if( _verbose )
    {
      c_finalSelection.print();
      cout << c_str << ": " << c_bits << endl;
    }

  {
    // tuple
    //    string dir_ = _e.categ();
    string dir_ = prefix;
    tm.setTree( "ZZtuple", dir_ );
    tm.add< int >(    "run",   &run       );
    tm.add< int >(    "event", &event     );
    tm.add< string*>( "categ", &categ_ptr );
    tm.add< string*>( "c_str", &c_str_ptr  );      
    tm.add< TBits* >( "c_bits", &c_bits_ptr  );      
    tm.add< int >(  "nZCand",   &nZCand        );
    tm.add< bool  >(  "isTight",&isTight    );
    tm.add< float >(  "mll",   &mll        );
    tm.add< float >(  "qTll",   &qTll        );
    tm.add< float >(  "phill",   &phill        );
    tm.add< float >(  "yll",   &yll        );
    tm.add< float >(  "hll",   &hll        );
    tm.add< float >(  "mTZZ",   &mTZZ      );
    tm.add< int >(  "nLoose",&nLoose    );
    tm.add< int >(  "nVertex",&nVertex    );
    tm.add< vectorFloat* >(  "LPt",   &LPt   );
    tm.add< vectorFloat* >(  "LEta",  &LEta  );
    tm.add< vectorFloat* >(  "LPhi",  &LPhi  );
    tm.add< vector<int>* >(  "LID",   &LID   );
    tm.add< vector<int>* >(  "LPdg",   &LPdg   );
    tm.add< vectorFloat* >(  "LMT",   &LMT   );
    tm.add< vectorFloat* >( "LTrkIso", &LTrkIso );
    tm.add< vectorFloat* >( "LEcalIso", &LEcalIso );
    tm.add< vectorFloat* >( "LHcalIso", &LHcalIso );
    tm.add< vectorFloat* >( "LCombIso", &LCombIso );
    tm.add< vector<int>* >(  "LVtx",   &LVtx   );
    tm.add< float >( "rhoFJ" , &rhoFJ );
    tm.add< float >(  "MET",   &MET        );
    tm.add< float >(  "METPhi",&METPhi    );
    tm.add< float >(  "sumEt", &sumEt );
    tm.add< float >(  "mEtSig", &mEtSig );
    tm.add< float >(  "projMET",&projMET    );
    tm.add< float >(  "corProjMET",&corProjMET    );
    tm.add< float >(  "sigMET",&sigMET    );
    tm.add< float >(  "dPhiMin",&dPhiMin    );
    tm.add< float >(  "balZMET",&balZMET    );

    tm.add< float >(  "cpuMETVtx",&cpuMETVtx    );
    tm.add< float >(  "cpuMETSum",&cpuMETSum    );
    tm.add< float >(  "cpuMETSqrt",&cpuMETSqrt  );

    tm.add< vectorFloat* > ( "jEt", &jEt );
    tm.add< vectorFloat* > ( "jEta", &jEta );
    tm.add< vectorFloat* > ( "jPhi", &jPhi );
    tm.add< vectorFloat* > ( "jDR1", &jDR1 );
    tm.add< vectorFloat* > ( "jDR2", &jDR2 );
    tm.add< vectorFloat* > ( "jDPhiMET", &jDPhiMET );
    tm.add< vectorFloat* > ( "jBtagTkC", &jBtagTkC );
    tm.add< vectorFloat* > ( "jBtagSoftM", &jBtagSoftM );
    tm.add< vector<int>* > ( "jVtx", &jVtx );
    tm.add< float > ( "JZB15", &JZB15 );
    tm.add< float > ( "JZB20", &JZB20 );
    tm.add< float > ( "JZB25", &JZB25 );
    tm.add< int >(  "jMult40",&jMult40  );
    tm.add< int >(  "jMult35",&jMult35  );
    tm.add< int >(  "jMult30",&jMult30  );
    tm.add< int >(  "jMult25",&jMult25  );
    tm.add< int >(  "jMult20",&jMult20  );
    tm.add< int >(  "jMult15",&jMult15  );
    tm.add< int >(  "jMult10",&jMult10  );
    tm.add< float >(  "dPhiJMET",&dPhiJMET );

    //MM : photons
    tm.add< vectorFloat* > ( "pEt", &pEt );
    tm.add< vectorFloat* > ( "pEta", &pEta );
    tm.add< vectorFloat* > ( "pPhi", &pPhi );
    tm.add< vectorFloat* > ( "pDR1", &pDR1 );
    tm.add< vectorFloat* > ( "pDR2", &pDR2 );
    tm.add< vectorFloat* > ( "pDPhiMET", &pDPhiMET );
 
    tm.flush();

    for( int ilep=0; ilep<2; ilep++ )
      {
	tm.setTree( "ZZLeptons", dir_ );
	tm.add< int >(    "run",   &run       );
	tm.add< int >(    "event", &event     );
	Candidate* LCand_ = LCand[ilep];
	CandInfo* info_ = LCand_->info();
	float pt_ = LCand_->charge()*LCand_->pt();
	float eta_ = LCand_->eta();
	float phi_ = LCand_->phi() * Constants::radToDeg;
	tm.add<float>(    "pt", &pt_ ); 
	tm.add<float>(    "eta", &eta_ );
	tm.add<float>(    "phi", &phi_ );

	if( imode==0 )
	  {
	    float fBrem  = info_->getFloat("fBrem"   );
	    float EOP    = info_->getFloat("EOP"     );
	    float dPhiIn = info_->getFloat("dPhiIn"  );
	    float dEtaIn = info_->getFloat("dEtaIn"  );

	    float dr03TrkIso = info_->getFloat("dr03TrkIso");
	    float dr04EcalIso = info_->getFloat("dr04EcalIso");
	    float dr04HcalD1Iso = info_->getFloat("dr04HcalD1Iso"); 
	    float dr04HcalD2Iso = info_->getFloat("dr04HcalD2Iso");	  

	    bool isConv = info_->getBool("isConv");
	    Candidate* gsfTrk_ = _e.getSecond( "electron-GsfTrack", LCand_ );
	    CandInfo* gsfInfo_ = gsfTrk_->info();
	    assert( gsfInfo_ );
	    bool isFPHitTrk = gsfInfo_->getBool("isFPHitTrk");
	    int  expInnHits = gsfInfo_->getInt("expInnHits");

	    tm.add< float >( "fBrem",  &fBrem  );
	    tm.add< float >( "EOP",    &EOP    );
	    tm.add< float >( "dPhiIn", &dPhiIn );
	    tm.add< float >( "dEtaIn", &dEtaIn );
	    tm.add< float >( "dr03TrkIso", &dr03TrkIso );
	    tm.add< float >( "dr04EcalIso", &dr04EcalIso );
	    tm.add< float >( "dr04HcalD1Iso", &dr04HcalD1Iso );
	    tm.add< float >( "dr04HcalD2Iso", &dr04HcalD2Iso );
	    tm.add< bool  >( "isConv", &isConv );
	    tm.add< bool  >( "isFPHitTrk", &isFPHitTrk );
	    tm.add< int   >( "expInnHits", &expInnHits );
	  }
	else
	  {
	    bool isGlobal = info_->getBool("muidGlobalMuonPromptTight");
	    bool muIdTracker = info_->getBool("muidTrackerMuonArbitrated");
	    float normChi2 = info_->getFloat("normChi2");
	    int nVHits = info_->getInt("nTrkHits"); 
	    int nMatch = info_->getInt("nMatch");
	    int nPixHit = info_->getInt("nPixHits");
	    int nMuonHit = info_->getInt("nMuonHits");
	    float isoR03strk = info_->getFloat("isoR03strk");
	    float isoR03sem  = info_->getFloat("isoR03sem");
	    float isoR03shad = info_->getFloat("isoR03shad");
	    float dPV = info_->getFloat("dPV"); //MM not dPV but d0
	    
	    tm.add< bool >( "isGlobal",  &isGlobal  );
	    tm.add< bool >( "muIdTracker",  &muIdTracker  );
	    tm.add< float >( "normChi2",    &normChi2    );
	    tm.add< int >( "nVHits", &nVHits );
	    tm.add< int >( "nMatch", &nMatch );
	    tm.add< int >( "nPixHit", &nPixHit );
	    tm.add< int >( "nMuonHit", &nMuonHit );
	    tm.add< float >( "isoR03strk",    &isoR03strk   );
	    tm.add< float >( "isoR03sem",    &isoR03sem    );
	    tm.add< float >( "isoR03shad",    &isoR03shad    );
	    tm.add< float >( "dPV",    &dPV    );
	  }
	tm.flush();
      }
  }    
 
  // garbage
  delete jEt;
  delete jEta;
  delete jPhi;
  delete jDR1;
  delete jDR2;
  delete jDPhiMET;
  delete jBtagTkC;
  delete jBtagSoftM;
  delete jVtx;
  delete pEt;
  delete pEta;
  delete pPhi;
  delete pDR1;
  delete pDR2;
  delete pDPhiMET;
  delete LPt;
  delete LEta;
  delete LPhi;
  delete LID;
  delete LPdg;
  delete LMT;
  delete LTrkIso;
  delete LEcalIso;
  delete LHcalIso;
  delete LCombIso;
  delete LVtx;

  float mVV(0);
  
  float  mV1(0);
  float ptV1(0);
  float pzV1(0);
  float  EV1(0);
  float  yV1(0);
  
  float  mV2(0);
  float ptV2(0);
  float pzV2(0);
  float  EV2(0);
  float  yV2(0);

  vectorFloat* ptL= new vectorFloat(4,0);
  vectorFloat* etaL= new vectorFloat(4,0);
  vectorFloat* phiL= new vectorFloat(4,0);
  vectorFloat* EL= new vectorFloat(4,0);
  vectorFloat* pdgL= new vectorFloat(4,0);

  bool storeMcTruth(false);
  if( _e.categ().substr(0,2) == "ZZ" ) //_2e2n" || _e.categ() == "ZZ_2m2n" )     
    {      
      storeMcTruth = true;
      Candidate* genCand_ = _e.decayTree();
      assert( genCand_!=0 );
	
      // get the two Z bosons
      CandList Z_;
      CandUtil::get(  23, genCand_, Z_ );
      if( Z_.size()==2 )
	{
      
	  Candidate* ZZ_ = Candidate::create( Z_[0], Z_[1] );
      
	  // search for neutrinos (12:nu_e; 14:nu_mu; 16:nu_tau)
	  //or quarks!!! MM (pdgId<6) 
	  //quarks are less interesting than neutrinos -> always 2nd Z
	  //after, taus
	  //then muons
	  int iV2 = 2; 
	  bool fZ=true;
	  int lvl= 0;
	  for( int iZ=0; iZ<2; iZ++ )
	    {
	      if( Z_[iZ]->nDaughters()!=2 ) continue;
	      int pdgCode_ =  Z_[iZ]->daughter(0)->pdgCode();


	      if( fabs(pdgCode_) == 11 && ( lvl<1 || fZ) ){
		iV2 = iZ; lvl=1;
	      }
	      else if( fabs(pdgCode_) == 13 && ( lvl<2 || fZ) ) {
		iV2 = iZ; lvl=2;
	      }
	      else if( fabs(pdgCode_) == 15 && ( lvl<3 || fZ) ) {
		iV2 = iZ; lvl=3;
	      }
	      else if( ( pdgCode_ == 12 || pdgCode_ == -12 ||
			 pdgCode_ == 14 || pdgCode_ == -14 ||
			 pdgCode_ == 16 || pdgCode_ == -16 ) && ( lvl<4 || fZ) ) {
		iV2 = iZ; lvl=4;
	      }
	      else if( fabs(pdgCode_) < 6 && ( lvl<5 || fZ) ) {
		iV2 = iZ;  lvl=5;
	      }
 
	      fZ =false;
	    }

	  assert( iV2<2 );
	  int iV1 = 1-iV2;
      
	  mVV  = ZZ_->mass();

	  mV1 = Z_[iV1]->mass();
	  ptV1 = Z_[iV1]->pt();
	  pzV1 = Z_[iV1]->pz();
	  EV1 = Z_[iV1]->E();
	  yV1 = KineUtils::y( EV1, pzV1 );

	  mV2 = Z_[iV2]->mass();
	  ptV2 = Z_[iV2]->pt();
	  pzV2 = Z_[iV2]->pz();
	  EV2 = Z_[iV2]->E();
	  yV2 = KineUtils::y( EV2, pzV2 );

	  for(int il=0;il<4;il++) {
	    int iZ = ( il/2 == iV1 ) ? 0:1;
	    if( Z_[iZ]->nDaughters()!=2 ) continue;
	    (*ptL)[ il] = ( Z_[ iZ ]->daughter(il%2)->pt() );
	    (*etaL)[ il] = ( Z_[ iZ ]->daughter(il%2)->eta() );
	    (*phiL)[ il] = ( Z_[ iZ ]->daughter(il%2)->phi() );
	    (*pdgL)[ il] = ( Z_[ iZ ]->daughter(il%2)->pdgCode() );
	    (*EL)[ il] = ( Z_[ iZ ]->daughter(il%2)->E() );
	  }
	}    
    }

  if( _e.categ().substr(0,2) == "WW" )//_2e2n" || _e.categ() == "WW_2m2n" )     
    {      
      storeMcTruth = true;
      Candidate* genCand_ = _e.decayTree();
      assert( genCand_!=0 );
	
      // get the two Z bosons
      CandList W_;
      CandUtil::get(  24, genCand_, W_ );
      CandUtil::get( -24, genCand_, W_ );
      if( W_.size()==2 )
	{
      
	  Candidate* WW_ = Candidate::create( W_[0], W_[1] );
      
	  mVV  = WW_->mass();

	  mV1 = W_[0]->mass();
	  ptV1 = W_[0]->pt();
	  pzV1 = W_[0]->pz();
	  EV1 = W_[0]->E();
	  yV1 = KineUtils::y( EV1, pzV1 );

	  mV2 = W_[1]->mass();
	  ptV2 = W_[1]->pt();
	  pzV2 = W_[1]->pz();
	  EV2 = W_[1]->E();
	  yV2 = KineUtils::y( EV2, pzV2 );

	  for(int il=0;il<4;il++) {

	    if( W_[ il/2 ]->nDaughters()!=2 ) continue;
	    (*ptL)[ il] = ( W_[ il/2 ]->daughter(il%2)->pt() );
	    (*etaL)[ il] = ( W_[ il/2 ]->daughter(il%2)->eta() );
	    (*phiL)[ il] = ( W_[ il/2 ]->daughter(il%2)->phi() );
	    (*pdgL)[ il] = ( W_[ il/2 ]->daughter(il%2)->pdgCode() );
	    (*EL)[ il] = ( W_[ il/2 ]->daughter(il%2)->E() );
	  }
	}
    }

  if( _e.categ().substr(0,2) == "WZ" )//_3en" || _e.categ() == "WZ_3mn" )     
    {      
      storeMcTruth = true;
      Candidate* genCand_ = _e.decayTree();
      assert( genCand_!=0 );
	
      // get the two Z bosons
      CandList V_;
      CandUtil::get(  24, genCand_, V_ );
      CandUtil::get( -24, genCand_, V_ );
      CandUtil::get(  23, genCand_, V_ );
      if( V_.size()==2 )
	{
      
	  Candidate* WZ_ = Candidate::create( V_[0], V_[1] );
      
	  mVV  = WZ_->mass();

	  mV1 = V_[0]->mass();
	  ptV1 = V_[0]->pt();
	  pzV1 = V_[0]->pz();
	  EV1 = V_[0]->E();
	  yV1 = KineUtils::y( EV1, pzV1 );

	  mV2 = V_[1]->mass();
	  ptV2 = V_[1]->pt();
	  pzV2 = V_[1]->pz();
	  EV2 = V_[1]->E();
	  yV2 = KineUtils::y( EV2, pzV2 );

	  for(int il=0;il<4;il++) {
	    if( V_[ il/2 ]->nDaughters()!=2 ) continue;
	    (*ptL)[ il] = ( V_[ il/2 ]->daughter(il%2)->pt() );
	    (*etaL)[ il] = ( V_[ il/2 ]->daughter(il%2)->eta() );
	    (*phiL)[ il] = ( V_[ il/2 ]->daughter(il%2)->phi() );
	    (*pdgL)[ il] = ( V_[ il/2 ]->daughter(il%2)->pdgCode() );
	    (*EL)[ il] = ( V_[ il/2 ]->daughter(il%2)->E() );
	  }
	}
    }
  
  if( storeMcTruth )
    {
      // tuple
      //      string dir_ = _e.categ();
      string dir_ = prefix;
      tm.setTree( "VV_mcTruth", dir_ );
      tm.add< int   >(    "run",   &run   );
      tm.add< int   >(    "event", &event );
      tm.add< string*>( "categ", &categ_ptr );
      tm.add< float >(    "mVV",   &mVV   );
      tm.add< float >(    "mV1",   &mV1   );
      tm.add< float >(    "ptV1",  &ptV1  );
      tm.add< float >(    "pzV1",  &pzV1  );
      tm.add< float >(    "EV1",   &EV1   );
      tm.add< float >(    "yV1",   &yV1   );
      tm.add< float >(    "mV2",   &mV2   );
      tm.add< float >(    "ptV2",  &ptV2  );
      tm.add< float >(    "pzV2",  &pzV2  );
      tm.add< float >(    "EV2",   &EV2   );
      tm.add< float >(    "yV2",   &yV2   );
      tm.add< vectorFloat* > ( "ptL", &ptL );
      tm.add< vectorFloat* > ( "etaL", &etaL );
      tm.add< vectorFloat* > ( "phiL", &phiL );
      tm.add< vectorFloat* > ( "pdgL", &pdgL );
      tm.add< vectorFloat* > ( "EL", &EL );
      
      tm.flush();
    
    }

  //garbage
  delete ptL;
  delete etaL;
  delete phiL;
  delete pdgL;
  delete EL;


  return isSelected;
}

bool
DibosonAnalysis::WW_em2n__analysis(  HypoListIterator ithyp )
{
  int run   = _e.run();
  int event = _e.event();

  // declare the main analysis variables
  int   nZCand(0);        // number of Z candidates
  int nVertex(0);         // number of good vertices
  float mll(0);           // di-lepton mass
  float qTll(0);          // di-lepton transverse momentum
  float phill(0);         // di-lepton phi
  float yll(0);           // di-lepton rapidity
  float hll(0);           // di-lepton helicity
  float mTZZ(0);          // Z+MET transverse mass
  float MET(0);           // PF-MET type-I
  float METPhi(0);        // PF-MET type-I - phi
  float sumEt(0);         // sum of transverse energies
  float mEtSig(0);        // sigma of MET
  float projMET(0);       // projected MET
  float corProjMET(0);    // PU-vertex corrected MET
  float sigMET(0);        // corrected MET significance
  float dPhiMin(0);       // minimum angle of MET to a lepton
  float balZMET;          // "balance" Z/MET
  bool  isTight(false);   // is the Z candidate tight?
  int   nLoose(0);        // number of additional loose leptons
  int   jMult30(0);       // number of jets above 30 GeV-Et
  int   jMult25(0);       // number of jets above 25 GeV-Et
  int   jMult20(0);       // number of jets above 25 GeV-Et
  int   jMult15(0);       // number of jets above 15 GeV-Et
  int   jMult10(0);       // number of jets above 10 GeV-Et
  float dPhiJMET(0);      // angle between leading jet and MET

  // create the jet vectors to be stored in the ntuple -- MUST be destroyed!
  // (warning: no return before the end)
  // Note to myself: implement garbage collection for these...
  vectorFloat* jEt        = new vectorFloat;
  vectorFloat* jEta       = new vectorFloat;
  vectorFloat* jPhi       = new vectorFloat;
  vectorFloat* jDR1       = new vectorFloat;
  vectorFloat* jDR2       = new vectorFloat;
  vectorFloat* jDPhiMET   = new vectorFloat;
  vectorFloat* jBtagTkC   = new vectorFloat;
  vectorFloat* jBtagSoftM = new vectorFloat;
  vector<int>* jVtx       = new vector<int>;

  // leptons -- the first two are from the Z
  vectorFloat* LPt   = new vectorFloat;
  vector<int>* LPdg = new vector<int>;
  vectorFloat* LEta  = new vectorFloat;
  vectorFloat* LPhi  = new vectorFloat;
  vectorFloat* LMT   = new vectorFloat;
  vector<int>* LID   = new vector<int>;
  vectorFloat* LTrkIso   = new vectorFloat;
  vectorFloat* LEcalIso  = new vectorFloat;
  vectorFloat* LHcalIso  = new vectorFloat;
  vector<int>* LVtx   = new vector<int>;

 
  //
  // Cuts
  //

  // Z quality
  OneSidedCut<int>   c_nZCand    ( "nZCand",      nZCand,    1,   "==" );
  TwoSidedCut<float> c_ZWindow   ( "ZWindow",     mll, 60., 120., "[]" );
  OneSidedCut<bool>  c_ZLQual    ( "ZLQual",     isTight, true,   "==" );
  CompositeCut        c_ZQual     ( "ZQual", "&&" );
  c_ZQual.add( c_nZCand );
  c_ZQual.add( c_ZWindow );
  c_ZQual.add( c_ZLQual );

  // MET cuts
  OneSidedCut<float> c_balZMET   ( "balZMET",    balZMET, 2.,     "<=" );
  OneSidedCut<float> c_MET0p50   ( "MET0p50",    sigMET,  0.50,   ">"  );
  OneSidedCut<float> c_MET1p00   ( "MET1p00",    sigMET,  1.00,   ">"  );
  OneSidedCut<float> c_MET1p50   ( "MET1p50",    sigMET,  1.50,   ">"  );
  OneSidedCut<float> c_MET2p00   ( "MET2p00",    sigMET,  2.00,   ">"  );
  OneSidedCut<float> c_MET2p50   ( "MET2p50",    sigMET,  2.50,   ">"  );
  OneSidedCut<float> c_MET3p00   ( "MET3p00",    sigMET,  3.00,   ">"  );
  OneSidedCut<float> c_MET3p50   ( "MET3p50",    sigMET,  3.50,   ">"  );
  OneSidedCut<float> c_MET4p00   ( "MET4p00",    sigMET,  4.00,   ">"  );
  OneSidedCut<float> c_MET4p20   ( "MET4p20",    sigMET,  4.20,   ">"  );
  CompositeCut       c_MET       ( "MET", "&&" );
  c_MET.add( c_balZMET );
  c_MET.add( c_MET0p50 );
  c_MET.add( c_MET1p00 );
  c_MET.add( c_MET1p50 );
  c_MET.add( c_MET2p00 );
  c_MET.add( c_MET2p50 );
  c_MET.add( c_MET3p00 );
  c_MET.add( c_MET3p50 );
  c_MET.add( c_MET4p00 );
  c_MET.add( c_MET4p20 );

  // jet veto
  OneSidedCut<int>   c_j30Veto   ( "j30Veto",    jMult30, 0,      "==" );
  OneSidedCut<int>   c_j25Veto   ( "j25Veto",    jMult25, 0,      "==" );
  OneSidedCut<int>   c_j20Veto   ( "j20Veto",    jMult20, 0,      "==" );
  OneSidedCut<int>   c_j15Veto   ( "j15Veto",    jMult15, 0,      "==" );
  TwoSidedCut<float> c_dPhiJMET  ( "dPhiJMET", dPhiJMET, -20.,20, "![]");
  CompositeCut c_j20Veto_OR_dPhiMET( "j20Veto_OR_dPhiMET", "||" );
  c_j20Veto_OR_dPhiMET.add( c_j20Veto  );
  c_j20Veto_OR_dPhiMET.add( c_dPhiJMET );  
  CompositeCut c_jVeto( "jVeto", "&&" );
  c_jVeto.add( c_j25Veto );
  c_jVeto.add( c_j20Veto_OR_dPhiMET );
  
  // diboson separation
  OneSidedCut<int>   c_LVeto     ( "LVeto",       nLoose, 0,      "==" );
  CompositeCut c_dibosonVeto( "dibosonSeparation", "&&" );
  c_dibosonVeto.add( c_LVeto );

  // final selection
  CompositeCut c_finalSelection( "finalSelection", "&&" );
  c_finalSelection.add( c_ZQual       );
  c_finalSelection.add( c_jVeto       );
  c_finalSelection.add( c_MET         );
  c_finalSelection.add( c_dibosonVeto );

  TBits   c_bits;
  TBits*  c_bits_ptr = &c_bits;
  string  c_str;
  string* c_str_ptr  = &c_str;
  c_finalSelection.getString(c_str);
  string categ = _e.categ();
  string* categ_ptr = & categ;

  // constants
  float sigma_0  = 7.34;
  float sigma_PU = 3.9;
  float sigma_T(0);
  float alpha = 0.856;

  // the hypothesis
  Hypo hyp    = ithyp->first;
  //  size_t ihyp = hyp.first;
  //  size_t jhyp = hyp.second;
  Mode kmode   = ithyp->second;

  // the lepton mode for this hypothesis
  //  int imode(0);
  assert( kmode==kWW_em2n );
  string prefix(modeName[kmode]); 

  nZCand = buildZemList();
  if( nZCand==0 ) return false;
  
  CandList ZemList = EventManager::e()->userCandList( "Zem" );
  Candidate* ZCand = ZemList[0];
  assert( ZCand!=0 );

  // the di-lepton invariant mass
  mll  = ZCand->mass();
  qTll = ZCand->pt();
  phill = ZCand->phi();
  CandInfo* ZInfo = ZCand->info();
  //  yll = ZInfo->getFloat("rapidity");
  //  hll = ZInfo->getFloat("helicity");
  yll = ZInfo->getFloat(0); //!!!
  hll = ZInfo->getFloat(0); //!!!

  // a Z boson candidate is tight if it is made out of at least one Tight 
  // and one Loose lepton
  isTight =  selectZCand( *ZCand, Z0_window[0], Z0_window[1], 0, 1, 2, 2 ); 

  // the missing transverse momentum 
  const Candidate* metCand =  _e.met( EventManager::kPfMet );
  assert( metCand!=0 );
  MET = metCand->pt();
  METPhi = metCand->phi();
  CandInfo* metInfo = metCand->info();
  sumEt = metInfo->getFloat( "sumEt" );
  mEtSig = metInfo->getFloat( "mEtSig" );

  Candidate* ZnnCand = metCand->clone();

  Candidate* ZZCand = Candidate::create( ZCand, ZnnCand );
  mTZZ = ZZCand->mass();
  
  // the leptons in the event (two of them come from the Z boson)  
  CandList leptonList[2] = 
    { EventManager::e()->userCandList( leptonListName(0) ),
      EventManager::e()->userCandList( leptonListName(1) ) };

  // - count the number of "loose" leptons in the rest of the event
  // - get the isolation of the two leptons from the Z candidate
  // - determine the minimal angle with the missing pT in the transverse plane
  nLoose = 0;
  //  vector< float > iso;
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
	      //	      float iso_(0);
	      //	      info->getFloat( "iso",iso_ );
	      //	      iso.push_back( iso_ );
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

  for( list< Candidate* >::iterator it=LCandAll.begin();
       it!=LCandAll.end(); ++it )
    {
      Candidate* LCand_ = *it;
      Vertex* vertex = LCand_->vertex();
      int vtxIndex_ = -1;
      //      if( vertex!=0 )
      //	vtxIndex_ = vertex->info()->getInt("index");
      if( vertex!=0 )
	{
	  CandInfo* vtxInfo_ = vertex->info();
	  if( vtxInfo_ )
	    {
	      vtxIndex_ = vtxInfo_->getInt("index");
	    }
	}
      int pdgCode_ = LCand_->pdgCode();
      LPt  -> push_back( LCand_->charge() * LCand_->pt()  );
      LPdg -> push_back( pdgCode_ );
      LEta -> push_back( LCand_->eta() );
      LPhi -> push_back( LCand_->phi() );
      LMT  -> push_back( mT(LCand_,ZnnCand) );
      LVtx -> push_back( vtxIndex_ );
      int id_(0);
      CandInfo* info = LCand_->info();
      if( info->getBool("VeryLoose") ) id_ += 1;
      if( info->getBool("Loose") )     id_ += 2;
      if( info->getBool("Tight") )     id_ += 4;
      if( info->getBool("VeryTight") ) id_ += 8;
      LID  -> push_back( id_ );

      // isolation variables
      float trkIso_(0);
      float ecalIso_(0);
      float hcalIso_(0);
      if( abs(pdgCode_)==13 )
	{
	  // muons
	  trkIso_ = info->getFloat("isoR03strk");
	  ecalIso_ = info->getFloat("isoR03sem");
	  hcalIso_ = info->getFloat("isoR03shad");
	}
      else if( abs(pdgCode_)==11 )
	{
	  // electrons
	  trkIso_ = info->getFloat("dr03TrkIso");
	  ecalIso_ = info->getFloat("dr04EcalIso");
	  hcalIso_ = 
	    info->getFloat("dr04HcalD1Iso") + 
	    info->getFloat("dr04HcalD2Iso");	  
	}
      else
	assert(0);
      LTrkIso->push_back( trkIso_ );
      LEcalIso->push_back( ecalIso_ );
      LHcalIso->push_back( hcalIso_ );
    }

  dPhiMin = 1000;
  for( size_t ii=0; ii<2; ii++ )
    {
      Candidate* LCand_ = LCand[ii];
      float dPhi_ = CandUtil::dPhi( ZnnCand, LCand_ ); 
      if( fabs( dPhi_ )<dPhiMin ) dPhiMin = fabs(dPhi_);
    }
  assert( dPhiMin<=Constants::pi );  

  // the projected MET 
  projMET = MET;
  if( dPhiMin < Constants::pi/2. )
    {
      projMET *= sin( dPhiMin );
    }

  // count the number of "good" vertices
  const VtxList& vertices = _e.vertices(); 
  bool findPV=false;
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
      
  // the corrected MET
  sigma_T    = sqrt( pow(sigma_0,2) + ( nVertex-1 )*pow(sigma_PU,2) );
  corProjMET = ( projMET - alpha * ( nVertex-1 ) ) * sigma_0 / sigma_T;
  sigMET     = corProjMET / sigma_0;
  if( fabs(qTll)>0 ) balZMET    = corProjMET / qTll;

  // jets: Particle flow with energy scale correction and lepton cleaning
  jMult30=0;
  jMult25=0;
  jMult20=0;
  jMult15=0;
  jMult10=0;
  const CandList& jets = _e.jetList( EventManager::kPatJet );
  CandList jList10; 
  CandList jList15; 
  CandList jList20; 
  CandList jList25; 
  CandList jList30; 
  for( size_t ij=0; ij<jets.size(); ij++ ) 
    { 
      Candidate* jet = jets[ij];
      Vertex* vertex = jet->vertex();
      int vtxIndex_ = -1;
      //      if( vertex!=0 )
      //	vtxIndex_ = vertex->info()->getInt("index");
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
      float jDR2_ = CandUtil::dPhi( jet, LCand[0] )*Constants::radToDeg ;
      CandInfo* info_ = jet->info();
      float jBtagTkC_   = info_->getFloat("btagTkC");
      float jBtagSoftM_ = info_->getFloat("btagSoftM");
      
      if( fabs(jEta_)<=5 )
	if( jEt_>=10 )
	{	  
	  jPhi    ->push_back(jPhi_ );
	  jEta    ->push_back(jEta_ );
	  jEt     ->push_back(jEt_ );
	  jDR1    ->push_back(jDR1_ ); 
	  jDR2    ->push_back(jDR2_ ); 
	  jDPhiMET->push_back(jDPhiMET_ ); 
	  jBtagTkC->push_back(jBtagTkC_);
	  jBtagSoftM->push_back(jBtagSoftM_);
	  jVtx    ->push_back(vtxIndex_);
	  	  
	  jList10.push_back( jet );
	  jMult10++;
	  if( jEt_>=15 ) 
	    {
	      jMult15++;
	      jList15.push_back( jet );
	    }
	  if( jEt_>=20 ) 
	    {
	      jMult20++;
	      jList20.push_back( jet );
	    }
	  if( jEt_>=25 ) 
	    {
	      jMult25++;
	      jList25.push_back( jet );
	    }
	  if( jEt_>=30 ) 
	    {
	      jMult30++;
	      jList30.push_back( jet );
	    }
	}      
    }
  sort( jList10.begin(), jList10.end(), Candidate::SortByPt() );
  sort( jList15.begin(), jList15.end(), Candidate::SortByPt() );
  sort( jList20.begin(), jList20.end(), Candidate::SortByPt() );
  sort( jList25.begin(), jList25.end(), Candidate::SortByPt() );
  sort( jList30.begin(), jList30.end(), Candidate::SortByPt() );

  // the "jet-recoil" candidate
  float JZB15(-1000.);
  float JZB20(-1000.);
  float JZB25(-1000.);
  float qT = ZCand->pt();

  if( jList15.size()>0 )
    {
      Candidate* jRecoil15   = Candidate::create( jList15 );
      JZB15 = jRecoil15->pt() - qT;
    }
  if( jList20.size()>0 )
    {
      Candidate* jRecoil20 = Candidate::create( jList20 );
      JZB20 = jRecoil20->pt() - qT;
    }
  if( jList25.size()>0 )
    {
      Candidate* jRecoil25 = Candidate::create( jList25 );
      JZB25 = jRecoil25->pt() - qT;
    }

  dPhiJMET = -190.;
  if( jList10.size()>0 ) dPhiJMET = (*jDPhiMET)[0];

  // selection cuts
  for( size_t ic=0; ic<=c_finalSelection.size(); ic++ )
    {
      string name_ = c_finalSelection.name(ic);
      for( int ii=0; ii<2; ii++ )
	{ 
	  string suffix_ = name_;
	  bool ok(true);
	  if( ii==0 )
	    {
	      if( ic==c_finalSelection.size() ) continue;	  
	      suffix_ += "/Nm1";
	      ok = c_finalSelection.NMinusOne( ic );
	    }
	  else
	    {
	      if( ic!=c_finalSelection.size() )
		suffix_ += "/Seq";
	      ok = c_finalSelection.cutUpTo( ic );
	    }
	  if( ok )
	    {
	      string dir_ = prefix;
	      dir_ += "/";
	      dir_ += suffix_;
	      fill( "mll",     "",        mll,        dir_ );
	      fill( "met",     "MET",     MET,        dir_ );
	      fill( "met",     "projMET", projMET,    dir_ );
	      fill( "met",     "corProj", corProjMET, dir_ );
	      fill( "balance", "ZMET",    balZMET,     dir_ );
	      fill( "sigMET",  "",        sigMET,     dir_ );
	      fill( "jmult",   "25GeV", jMult25,    dir_ );
	      fill( "jmult",   "20GeV", jMult20,    dir_ );
	      fill( "jmult",   "15GeV", jMult15,    dir_ );
	      fill( "dphi",    "jMET", dPhiJMET, dir_ );
	      fill( "dphi",    "LMET", dPhiMin,  dir_ );
	    }
	}
    }

  bool isSelected = c_finalSelection(); 


  if( isSelected )
    {
      //  hypothesis is selected
      select( hyp, kmode, kSel );
    }

  c_finalSelection.getBits(c_bits);
  if( _verbose )
    {
      c_finalSelection.print();
      cout << c_str << ": " << c_bits << endl;
    }

  {
    // tuple
    //    string dir_ = _e.categ();
    string dir_ = prefix;
    tm.setTree( "WWtuple", dir_ );
    tm.add< int >(    "run",   &run       );
    tm.add< int >(    "event", &event     );
    tm.add< string*>( "categ", &categ_ptr );
    tm.add< string*>( "c_str", &c_str_ptr  );      
    tm.add< TBits* >( "c_bits", &c_bits_ptr  );      
    tm.add< int >(  "nZCand",   &nZCand        );
    tm.add< bool  >(  "isTight",&isTight    );
    tm.add< float >(  "mll",   &mll        );
    tm.add< float >(  "qTll",   &qTll        );
    tm.add< float >(  "phill",   &phill        );
    tm.add< float >(  "yll",   &yll        );
    tm.add< float >(  "hll",   &hll        );
    tm.add< float >(  "mTZZ",   &mTZZ      );
    tm.add< int >(  "nLoose",&nLoose    );
    tm.add< int >(  "nVertex",&nVertex    );
    tm.add< vectorFloat* >(  "LPt",   &LPt   );
    tm.add< vectorFloat* >(  "LEta",  &LEta  );
    tm.add< vectorFloat* >(  "LPhi",  &LPhi  );
    tm.add< vector<int>* >(  "LID",   &LID   );
    tm.add< vector<int>* >(  "LPdg",   &LPdg   );
    tm.add< vectorFloat* >(  "LMT",   &LMT   );
    tm.add< vectorFloat* >( "LTrkIso", &LTrkIso );
    tm.add< vectorFloat* >( "LEcalIso", &LEcalIso );
    tm.add< vectorFloat* >( "LHcalIso", &LHcalIso );
    tm.add< vector<int>* >(  "LVtx",   &LVtx   );
    tm.add< float >(  "MET",   &MET        );
    tm.add< float >(  "METPhi",&METPhi    );
    tm.add< float >(  "sumEt", &sumEt );
    tm.add< float >(  "mEtSig", &mEtSig );
    tm.add< float >(  "projMET",&projMET    );
    tm.add< float >(  "corProjMET",&corProjMET    );
    tm.add< float >(  "sigMET",&sigMET    );
    tm.add< float >(  "dPhiMin",&dPhiMin    );
    tm.add< float >(  "balZMET",&balZMET    );
    tm.add< vectorFloat* > ( "jEt", &jEt );
    tm.add< vectorFloat* > ( "jEta", &jEta );
    tm.add< vectorFloat* > ( "jPhi", &jPhi );
    tm.add< vectorFloat* > ( "jDR1", &jDR1 );
    tm.add< vectorFloat* > ( "jDR2", &jDR2 );
    tm.add< vectorFloat* > ( "jDPhiMET", &jDPhiMET );
    tm.add< vectorFloat* > ( "jBtagTkC", &jBtagTkC );
    tm.add< vectorFloat* > ( "jBtagSoftM", &jBtagSoftM );
    tm.add< vector<int>* > ( "jVtx", &jVtx );
    tm.add< float > ( "JZB15", &JZB15 );
    tm.add< float > ( "JZB20", &JZB20 );
    tm.add< float > ( "JZB25", &JZB25 );
    tm.add< int >(  "jMult30",&jMult30  );
    tm.add< int >(  "jMult25",&jMult25  );
    tm.add< int >(  "jMult20",&jMult20  );
    tm.add< int >(  "jMult15",&jMult15  );
    tm.add< int >(  "jMult10",&jMult10  );
    tm.add< float >(  "dPhiJMET",&dPhiJMET );
    tm.flush();

//     for( int ilep=0; ilep<2; ilep++ )
//       {
// 	tm.setTree( "WWLeptons", dir_ );
// 	tm.add< int >(    "run",   &run       );
// 	tm.add< int >(    "event", &event     );
// 	Candidate* LCand_ = LCand[ilep];
// 	CandInfo* info_ = LCand_->info();
// 	float pt_ = LCand_->charge()*LCand_->pt();
// 	float eta_ = LCand_->eta();
// 	float phi_ = LCand_->phi() * Constants::radToDeg;
// 	tm.add<float>(    "pt", &pt_ ); 
// 	tm.add<float>(    "eta", &eta_ );
// 	tm.add<float>(    "phi", &phi_ );

// 	if( imode==0 )
// 	  {
// 	    float fBrem  = info_->getFloat("fBrem"   );
// 	    float EOP    = info_->getFloat("EOP"     );
// 	    float dPhiIn = info_->getFloat("dPhiIn"  );
// 	    float dEtaIn = info_->getFloat("dEtaIn"  );

// 	    float dr03TrkIso = info_->getFloat("dr03TrkIso");
// 	    float dr04EcalIso = info_->getFloat("dr04EcalIso");
// 	    float dr04HcalD1Iso = info_->getFloat("dr04HcalD1Iso"); 
// 	    float dr04HcalD2Iso = info_->getFloat("dr04HcalD2Iso");	  

// 	    bool isConv = info_->getBool("isConv");
// 	    Candidate* gsfTrk_ = _e.getSecond( "electron-GsfTrack", LCand_ );
// 	    CandInfo* gsfInfo_ = gsfTrk_->info();
// 	    assert( gsfInfo_ );
// 	    bool isFPHitTrk = gsfInfo_->getBool("isFPHitTrk");
// 	    int  expInnHits = gsfInfo_->getInt("expInnHits");

// 	    tm.add< float >( "fBrem",  &fBrem  );
// 	    tm.add< float >( "EOP",    &EOP    );
// 	    tm.add< float >( "dPhiIn", &dPhiIn );
// 	    tm.add< float >( "dEtaIn", &dEtaIn );
// 	    tm.add< float >( "dr03TrkIso", &dr03TrkIso );
// 	    tm.add< float >( "dr04EcalIso", &dr04EcalIso );
// 	    tm.add< float >( "dr04HcalD1Iso", &dr04HcalD1Iso );
// 	    tm.add< float >( "dr04HcalD2Iso", &dr04HcalD2Iso );
// 	    tm.add< bool  >( "isConv", &isConv );
// 	    tm.add< bool  >( "isFPHitTrk", &isFPHitTrk );
// 	    tm.add< int   >( "expInnHits", &expInnHits );
// 	  }
// 	else
// 	  {
// 	    bool isGlobal = info_->getBool("muidGlobalMuonPromptTight");
// 	    bool muIdTracker = info_->getBool("muidTrackerMuonArbitrated");
// 	    float normChi2 = info_->getFloat("normChi2");
// 	    int nVHits = info_->getInt("nTrkHits"); 
// 	    int nMatch = info_->getInt("nMatch");
// 	    int nPixHit = info_->getInt("nPixHits");
// 	    int nMuonHit = info_->getInt("nMuonHits");
// 	    float isoR03strk = info_->getFloat("isoR03strk");
// 	    float isoR03sem  = info_->getFloat("isoR03sem");
// 	    float isoR03shad = info_->getFloat("isoR03shad");
// 	    float dPV = info_->getFloat("dPV");
	    
// 	    tm.add< bool >( "isGlobal",  &isGlobal  );
// 	    tm.add< bool >( "muIdTracker",  &muIdTracker  );
// 	    tm.add< float >( "normChi2",    &normChi2    );
// 	    tm.add< int >( "nVHits", &nVHits );
// 	    tm.add< int >( "nMatch", &nMatch );
// 	    tm.add< int >( "nPixHit", &nPixHit );
// 	    tm.add< int >( "nMuonHit", &nMuonHit );
// 	    tm.add< float >( "isoR03strk",    &isoR03strk   );
// 	    tm.add< float >( "isoR03sem",    &isoR03sem    );
// 	    tm.add< float >( "isoR03shad",    &isoR03shad    );
// 	    tm.add< float >( "dPV",    &dPV    );
// 	  }
// 	tm.flush();
//       }
  }    
  
  // garbage
  delete jEt;
  delete jEta;
  delete jPhi;
  delete jDR1;
  delete jDR2;
  delete jDPhiMET;
  delete jBtagTkC;
  delete jBtagSoftM;
  delete jVtx;
  delete LPt;
  delete LEta;
  delete LPhi;
  delete LID;
  delete LPdg;
  delete LMT;
  delete LTrkIso;
  delete LEcalIso;
  delete LHcalIso;
  delete LVtx;

  float mVV(0);
  
  float  mV1(0);
  float ptV1(0);
  float pzV1(0);
  float  EV1(0);
  float  yV1(0);
  
  float  mV2(0);
  float ptV2(0);
  float pzV2(0);
  float  EV2(0);
  float  yV2(0);

  bool storeMcTruth(false);
  if( _e.categ() == "ZZ_2e2n" || _e.categ() == "ZZ_2m2n" )     
    {      
      storeMcTruth = true;
      Candidate* genCand_ = _e.decayTree();
      assert( genCand_!=0 );
	
      // get the two Z bosons
      CandList Z_;
      CandUtil::get(  23, genCand_, Z_ );
      assert( Z_.size()==2 );
      
      Candidate* ZZ_ = Candidate::create( Z_[0], Z_[1] );
      
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
      
      mVV  = ZZ_->mass();

      mV1 = Z_[iV1]->mass();
      ptV1 = Z_[iV1]->pt();
      pzV1 = Z_[iV1]->pz();
      EV1 = Z_[iV1]->E();
      yV1 = KineUtils::y( EV1, pzV1 );

      mV2 = Z_[iV2]->mass();
      ptV2 = Z_[iV2]->pt();
      pzV2 = Z_[iV2]->pz();
      EV2 = Z_[iV2]->E();
      yV2 = KineUtils::y( EV2, pzV2 );
    }

  if( _e.categ() == "WW_2e2n" || _e.categ() == "WW_2m2n" || _e.categ()=="WW_em2n" )     
    {      
      storeMcTruth = true;
      Candidate* genCand_ = _e.decayTree();
      assert( genCand_!=0 );
	
      // get the two Z bosons
      CandList W_;
      CandUtil::get(  24, genCand_, W_ );
      CandUtil::get( -24, genCand_, W_ );
      assert( W_.size()==2 );
      
      Candidate* WW_ = Candidate::create( W_[0], W_[1] );
      
      mVV  = WW_->mass();

      mV1 = W_[0]->mass();
      ptV1 = W_[0]->pt();
      pzV1 = W_[0]->pz();
      EV1 = W_[0]->E();
      yV1 = KineUtils::y( EV1, pzV1 );

      mV2 = W_[1]->mass();
      ptV2 = W_[1]->pt();
      pzV2 = W_[1]->pz();
      EV2 = W_[1]->E();
      yV2 = KineUtils::y( EV2, pzV2 );
    }

  if( _e.categ() == "WZ_3en" || _e.categ() == "WZ_3mn" )     
    {      
      storeMcTruth = true;
      Candidate* genCand_ = _e.decayTree();
      assert( genCand_!=0 );
	
      // get the two Z bosons
      CandList V_;
      CandUtil::get(  24, genCand_, V_ );
      CandUtil::get( -24, genCand_, V_ );
      CandUtil::get(  23, genCand_, V_ );
      assert( V_.size()==2 );
      
      Candidate* WZ_ = Candidate::create( V_[0], V_[1] );
      
      mVV  = WZ_->mass();

      mV1 = V_[0]->mass();
      ptV1 = V_[0]->pt();
      pzV1 = V_[0]->pz();
      EV1 = V_[0]->E();
      yV1 = KineUtils::y( EV1, pzV1 );

      mV2 = V_[1]->mass();
      ptV2 = V_[1]->pt();
      pzV2 = V_[1]->pz();
      EV2 = V_[1]->E();
      yV2 = KineUtils::y( EV2, pzV2 );
    }
  
  if( storeMcTruth )
    {
      // tuple
      //      string dir_ = _e.categ();
      string dir_ = prefix;
      tm.setTree( "VV_mcTruth", dir_ );
      tm.add< int   >(    "run",   &run   );
      tm.add< int   >(    "event", &event );
      tm.add< string*>( "categ", &categ_ptr );
      tm.add< float >(    "mVV",   &mVV   );
      tm.add< float >(    "mV1",   &mV1   );
      tm.add< float >(    "ptV1",  &ptV1  );
      tm.add< float >(    "pzV1",  &pzV1  );
      tm.add< float >(    "EV1",   &EV1   );
      tm.add< float >(    "yV1",   &yV1   );
      tm.add< float >(    "mV2",   &mV2   );
      tm.add< float >(    "ptV2",  &ptV2  );
      tm.add< float >(    "pzV2",  &pzV2  );
      tm.add< float >(    "EV2",   &EV2   );
      tm.add< float >(    "yV2",   &yV2   );
      
      tm.flush();
    }

  return isSelected;
}

bool
DibosonAnalysis::Z_2l__analysis(  HypoListIterator ithyp )
{
  int run   = _e.run();
  int event = _e.event();

  // declare the main analysis variables
  int   nZCand(0);        // number of Z candidates
  //  int nVertex(0);         // number of good vertices

  float mll(0);           // di-lepton mass of first pair
  float qTll(0);          // di-lepton transverse momentum
  float yll(0);           // di-lepton rapidity
  float hll(0);           // di-lepton helicity

  int nLoose(0);

  Hypo hyp    = ithyp->first;
  size_t ihyp = hyp.first;
  size_t jhyp = hyp.second;
  Mode kmode   = ithyp->second;

  int imode(0);
  if( kmode==kZ_2e ) imode=0;
  else if( kmode==kZ_2m ) imode=1;
  else assert(0);

  string prefix(modeName[kmode]); 

  CandList ZList[2] = 
    { EventManager::e()->userCandList( ZListName(0,ihyp) ),
      EventManager::e()->userCandList( ZListName(1,jhyp) )  };
    
  // take the Z candidate as the best in the list
  Candidate* ZCand;
  assert( ZList[imode].size()>0 );
  ZCand = ZList[imode][0];
  assert( ZCand!=0 );

  // A Z boson candidate is tight if it is made out of at least one Tight 
  // and one Loose leptons
  bool isTight =  selectZCand( *ZCand, Z0_fullDY[0], Z0_fullDY[1], 0, 1, 2, 2 );
   
  CandList leptonList[2] = 
    { EventManager::e()->userCandList( leptonListName(0) ),
      EventManager::e()->userCandList( leptonListName(1) ) };
  
  // count the number of leptons in the rest of the event
  //   and get the isolation of the two leptons from the Z candidate
  map<string, int> nLepSel[2];
  vector<float> iso;
  vector< Candidate* > LCand;
  for( size_t ilep=0; ilep<2; ilep++ )
    {
      for( size_t ii=0; ii<leptonList[ilep].size(); ii++ )
	{
	  Candidate* LCand_ = leptonList[ilep][ii];
	  
	  CandInfo* info = LCand_->info();
	  if( CandUtil::overlap( *LCand_, *ZCand ) ) 
	    {
	      LCand.push_back( LCand_ );
	      float iso_(0.);
	      info->getFloat( "iso", iso_ );
	      iso.push_back( iso_ );
	    }
	  else
	    {
	      for( int lepSel=0; lepSel<LeptonSelector::kNSel; lepSel++ )
		{
		  string str_ = LeptonSelector::selection[lepSel];
		  if( info->getBool(str_.c_str()) )  nLepSel[ilep][str_]++;
		}
	    }
	}
    }
  assert( iso.size()==2 );
  sort( iso.begin(), iso.end() );

  mll = ZCand->mass();

  if( !isTight ) return false;

  string categ = _e.categ();
  string* categ_ptr = & categ;

  // require no other "loose" lepton in the rest of the event
  for( size_t ilep=0; ilep<2; ilep++ )
    {
      for( size_t ii=0; ii<leptonList[ilep].size(); ii++ )
	{
	  CandInfo* info = leptonList[ilep][ii]->info();
	  if( info->getBool("Loose") ) nLoose++;       	 
	}
    }
  //  fill( "n10", "loose", nLoose, prefix );
  //  if( nLoose>0 ) return false;
  
  mll = ZCand->mass();
  qTll = ZCand->pt();
  CandInfo* ZInfo = ZCand->info();
  yll = ZInfo->getFloat("rapidity");
  hll = ZInfo->getFloat("helicity");

  {
    // tuple
    //    string dir_ = _e.categ();
    string dir_ = prefix;
    tm.setTree( "Ztuple", dir_ );
    tm.add< int >(    "run",   &run       );
    tm.add< int >(    "event", &event     );
    tm.add< string*>( "categ", &categ_ptr );
    tm.add< int >(  "nZCand",   &nZCand        );
    tm.add< bool >(  "isTight",   &isTight        );
    tm.add< float >(  "mll",   &mll        );
    tm.add< float >(  "qTll",   &qTll        );
    tm.add< float >(  "yll",   &yll        );
    tm.add< float >(  "hll",   &hll        );
    tm.add< int >(  "nLoose",&nLoose    );
    tm.flush();

    for( int ilep=0; ilep<2; ilep++ )
      {
	tm.setTree( "ZLeptons", dir_ );
	tm.add< int >(    "run",   &run       );
	tm.add< int >(    "event", &event     );
	Candidate* LCand_ = LCand[ilep];
	CandInfo* info_ = LCand_->info();
	float pt_ = LCand_->charge()*LCand_->pt();
	float eta_ = LCand_->eta();
	float phi_ = LCand_->phi() * Constants::radToDeg;
	tm.add<float>(    "pt", &pt_ ); 
	tm.add<float>(    "eta", &eta_ );
	tm.add<float>(    "phi", &phi_ );

	if( imode==0 )
	  {
	    float fBrem  = info_->getFloat("fBrem"   );
	    float EOP    = info_->getFloat("EOP"     );
	    float dPhiIn = info_->getFloat("dPhiIn"  );
	    float dEtaIn = info_->getFloat("dEtaIn"  );

	    float dr03TrkIso = info_->getFloat("dr03TrkIso");
	    float dr04EcalIso = info_->getFloat("dr04EcalIso");
	    float dr04HcalD1Iso = info_->getFloat("dr04HcalD1Iso"); 
	    float dr04HcalD2Iso = info_->getFloat("dr04HcalD2Iso");	  

	    tm.add< float >( "fBrem",  &fBrem  );
	    tm.add< float >( "EOP",    &EOP    );
	    tm.add< float >( "dPhiIn", &dPhiIn );
	    tm.add< float >( "dEtaIn", &dEtaIn );
	    tm.add< float >( "dr03TrkIso", &dr03TrkIso );
	    tm.add< float >( "dr04EcalIso", &dr04EcalIso );
	    tm.add< float >( "dr04HcalD1Iso", &dr04HcalD1Iso );
	    tm.add< float >( "dr04HcalD2Iso", &dr04HcalD2Iso );
	  }
	else
	  {
	    bool isGlobal = info_->getBool("muidGlobalMuonPromptTight");
	    bool muIdTracker = info_->getBool("muidTrackerMuonArbitrated");
	    float normChi2 = info_->getFloat("normChi2");
	    int nVHits = info_->getInt("nTrkHits"); 
	    int nMatch = info_->getInt("nMatch");
	    int nPixHit = info_->getInt("nPixHits");
	    int nMuonHit = info_->getInt("nMuonHits");
	    float isoR03strk = info_->getFloat("isoR03strk");
	    float isoR03sem  = info_->getFloat("isoR03sem");
	    float isoR03shad = info_->getFloat("isoR03shad");
	    float dPV = info_->getFloat("d0"); //MM not dPV but d0
	    
	    tm.add< bool >( "isGlobal",  &isGlobal  );
	    tm.add< bool >( "muIdTracker",  &muIdTracker  );
	    tm.add< float >( "normChi2",    &normChi2    );
	    tm.add< int >( "nVHits", &nVHits );
	    tm.add< int >( "nMatch", &nMatch );
	    tm.add< int >( "nPixHit", &nPixHit );
	    tm.add< int >( "nMuonHit", &nMuonHit );
	    tm.add< float >( "isoR03strk",    &isoR03strk   );
	    tm.add< float >( "isoR03sem",    &isoR03sem    );
	    tm.add< float >( "isoR03shad",    &isoR03shad    );
	    tm.add< float >( "dPV",    &dPV    );
	  }
	tm.flush();
      }
  }    

  //  hypothesis is selected
  select( hyp, kmode, kSel );

  return true;
}

bool
DibosonAnalysis::W_ln__analysis(  HypoListIterator ithyp )
{
  return false;

  Hypo hyp    = ithyp->first;
  size_t ihyp = hyp.first;
  size_t jhyp = hyp.second;
  Mode kmode   = ithyp->second;

  int imode(0);
  if( kmode==kW_en ) imode=0;
  else if( kmode==kW_mn ) imode=1;
  else assert(0);

  string subdir(modeName[kmode]); 

  CandList ZList[2] = 
    { EventManager::e()->userCandList( ZListName(0,ihyp) ),
      EventManager::e()->userCandList( ZListName(1,jhyp) )  };
    
  // exit if there is a Z candidate in the event
  if( ZList[0].size()!=0 || ZList[1].size()!=0 ) return false;

  CandList leptonList[2] = 
    { EventManager::e()->userCandList( leptonListName(0) ),
      EventManager::e()->userCandList( leptonListName(1) ) };
  
  if( leptonList[imode].size() == 0 ) return false;

  // take the "best" lepton of the list
  Candidate* LCand = leptonList[imode][0];
  assert( LCand!=0 );

  // charge 
  int chg = ( LCand->charge()>0 ) ? 1:-1;

  // pdgCode
  int LPdgCode_  = LCand->pdgCode();
  int nuPdgCode_ = -1*(LPdgCode_-chg);

  CandInfo* info = LCand->info();
  float iso(0.);
  info->getFloat("iso", iso );
  int id(-1);
  for( int lepSel=0; lepSel<LeptonSelector::kNSel; lepSel++ )
    {
      const string & str_ = LeptonSelector::selection[lepSel];
      if( info->getBool(str_.c_str()) )  id = lepSel;
    }

  if( id<LeptonSelector::kTight ) return false;

  string sgn_ = "m";
  if( chg==1 ) sgn_ = "p"; 

  // get the particle flow MET
  const Candidate* metCand =  _e.met( EventManager::kPfMet );
  assert( metCand!=0 );

  // create the neutrino candidate
  Candidate* nuCand = metCand->clone();
  nuCand->setPdgCode( nuPdgCode_ ); 
  TVector2 nuP2_ = nuCand->p2();

  // create the lepton+nu candidate
  Candidate* LNuCand = Candidate::create( LCand, nuCand );
  LNuCand->setName("LNu");

  // create the W candidate
  Candidate* WCand = LNuCand->clone();
  WCand->setPdgCode( chg*24 );

  // the missing transverse energy
  float MET   = nuCand->pt();
  float METPhi = nuCand->phi();

  // the transverse mass
  float MT    = LNuCand->mass();

  // the W pT
  float WPt   = LNuCand->pt();
  float WPhi  = LNuCand->phi();

  // the LP variable
  float LP = -1000;
  if( fabs( WPt )>0. )
    {
      LP  = LNuCand->p2() * LCand->p2() ;
      LP /= pow( WPt, 2 );
    }
  
  // projections
  // get the unit vectors along and perpendicular to LCand
  TVector2 ux_ = LCand->p2().Unit();
  TVector2 uy_( -ux_.Y(), ux_.X() ); 
  float proj_x = nuP2_ * ux_;
  float proj_y = nuP2_ * uy_;

  float dphi_(0);
  float dpt_(0);
  if( EventServer::isMC ) 
    {
      TVector2 neuP2 = _e.neutrinosFromBosonsP4().Vect().XYvector();
      dphi_  = nuCand->p2().DeltaPhi( neuP2 ) * Constants::radToDeg;
      dpt_   = nuCand->pt()-neuP2.Mod();
    }

  // jets: Particle flow with energy scale correction and lepton cleaning
  int jmult(0);
  float jEtThres = 15;
  float jEtaMax  = 2.4;
  const CandList& jets = _e.jetList( EventManager::kPatJet);

  // the "jet-recoil" candidate
  Candidate* jRecoil = Candidate::create( jets );
  float jRecoilPt(-1000);
  float jRecoilDphi(-1000);
  float jRecoilDpt(-1000);
  if( jets.size()>0 ) 
    {
      jRecoilPt = jRecoil->pt();
      jRecoilDphi  = CandUtil::dPhi( jRecoil, LNuCand ) * Constants::radToDeg;
      jRecoilDpt   = jRecoilPt-WPt;     
    }

  CandList jList; 
  for( size_t ij=0; ij<jets.size(); ij++ ) 
    { 
      Candidate* jet = jets[ij];
      float jeta_ = jet->eta();
      float jEt_  = jet->Et();
      if( jEt_>=jEtThres && fabs(jeta_)<=jEtaMax )
	{
	  jList.push_back( jet );
	  jmult++;
	}      
    }
  sort( jList.begin(), jList.end(), Candidate::SortByPt() );
  int jMult = (int)jList.size();
  
  // create the vectors to be stored in the ntuple -- MUST be destroyed!
  vectorFloat* jEt     = new vectorFloat;
  vectorFloat* jEta    = new vectorFloat;
  vectorFloat* jDRL    = new vectorFloat;
  vectorFloat* jDPhiL  = new vectorFloat;
  vectorFloat* jDPhiNu = new vectorFloat;
  vectorFloat* jDPhiW  = new vectorFloat;
  for( size_t ij=0; ij<jList.size(); ij++ ) 
    { 
      Candidate* jet = jets[ij];
      float jEta_    = jet->eta() ;
      float jEt_     = jet->Et()  ;
      float jDRL_    = CandUtil::dR(   jet, LCand ) ;
      float jDPhiL_  = CandUtil::dPhi( jet, LCand  )*Constants::radToDeg ;
      float jDPhiNu_ = CandUtil::dPhi( jet, nuCand )*Constants::radToDeg ;
      float jDPhiW_  = CandUtil::dPhi( jet, WCand  )*Constants::radToDeg ;

      if( _verbose )
	{
	  cout << "Jet " << ij << endl;
	  cout << "\tEt=" << jEt_;
	  cout << " Eta=" << jEta_;
	  cout << " DRL=" << jDRL_;
	  cout << " DPhiL=" << jDPhiL_;
	  cout << " DPhiNu=" << jDPhiNu_;
	  cout << " DPhiW=" << jDPhiW_;
	  cout << endl;
	}

      jEta    ->push_back(jEta_ );
      jEt     ->push_back(jEt_ );
      jDRL    ->push_back(jDRL_ );
      jDPhiL  ->push_back(jDPhiL_ );
      jDPhiNu ->push_back(jDPhiNu_ );
      jDPhiW  ->push_back(jDPhiW_ );
    }

  // histograms

  bool isInEtaBin[7];  
  bool isInPtRange[4];  
  bool tightness[2];
  float LEta(0.);
  float LPt(0.);
  float LPhi(0.);

  tightness[0] = id>=LeptonSelector::kMedium;
  tightness[1] = id>=LeptonSelector::kTight;
  if( imode==0 )
    {
      LEta = info->getFloat("caloEta"); 
      float E_  = info->getFloat("caloEnergy");
      LPt = E_/cosh(LEta);
      LPhi = info->getFloat("caloPhi"); 

      // electron
      isInEtaBin[0] =  true;
      isInEtaBin[1] =  (LEta>=0  && LEta<=0.4);
      isInEtaBin[2] =  (LEta>0.4 && LEta<=0.8); 
      isInEtaBin[3] =  (LEta>0.8 && LEta<=1.2); 
      isInEtaBin[4] =  (LEta>1.2 && LEta<=1.4); 
      isInEtaBin[5] =  (LEta>1.6 && LEta<=2.0); 
      isInEtaBin[6] =  (LEta>2.0 && LEta<=2.4); 
      
      //      isInPtRange[0] = true;
      //      isInPtRange[1] = info->getBool("isEt20"); 
      //      isInPtRange[2] = info->getBool("isEt25"); 
      //      isInPtRange[3] = info->getBool("isEt30"); 

    }
  else
    {
      LEta = LCand->eta();
      LPt  = LCand->pt();
      LPhi = LCand->phi();
      
      // electron
      isInEtaBin[0] =  true;
      isInEtaBin[1] =  (LEta>=0  && LEta<=0.4);
      isInEtaBin[2] =  (LEta>0.4 && LEta<=0.8); 
      isInEtaBin[3] =  (LEta>0.8 && LEta<=1.2); 
      isInEtaBin[4] =  (LEta>1.2 && LEta<=1.5); 
      isInEtaBin[5] =  (LEta>1.5 && LEta<=1.8); 
      isInEtaBin[6] =  (LEta>1.8 && LEta<=2.1); 

//       isInPtRange[0] = true; 
//       isInPtRange[1] = info->getBool("isPt20"); 
//       isInPtRange[2] = info->getBool("isPt25"); 
//       isInPtRange[3] = info->getBool("isPt30"); 
    }
  isInPtRange[0] = true;
  isInPtRange[1] = LPt>=27;
  isInPtRange[2] = LPt>=30;
  isInPtRange[3] = LPt>=32;

  bool isMT40  = MT>=40.;
  bool isMET20 = MET>=20.;
  
  TBits* idBits = new TBits(4);
  if( id>=LeptonSelector::kTight ) idBits->SetBitNumber(3, kTRUE);
  if( id>=LeptonSelector::kMedium     ) idBits->SetBitNumber(2, kTRUE);
  if( id>=LeptonSelector::kLoose     ) idBits->SetBitNumber(1, kTRUE);
  if( id>=LeptonSelector::kVeryLoose ) idBits->SetBitNumber(0, kTRUE);

  for( int kk=0; kk<7; kk++ )
    // eta
    {
      if( !isInEtaBin[kk] ) continue;
      for( int jj=0; jj<4; jj++ )
	// pt 
	{
	  if( !isInPtRange[jj] ) continue;
	  for( int ii=0; ii<2; ii++ )
	    // loop on charge
	    {
	      for( int ll=0; ll<2; ll++ )
		// tightness
		{
		  if( !tightness[ll] ) continue;

		  ostringstream o;
		  o << subdir << "/" << "eta" << kk << "/" << "pt" << jj << "/";
		  string dir_(o.str());

		  string dumstr_("W");
		  if( ii==1 ) dumstr_+=sgn_;
		  if( ll==1 ) dumstr_+="__VT";
		  else        dumstr_+="__T";
		  
		  string suffix = (dumstr_);

		  // isolation
		  fill( "iso"    , suffix, iso, dir_ );
		  
		  if( EventServer::isMC ) 
		    {
		      fill( "dphi"   , suffix, dphi_,   dir_ );
		      fill( "dpt"    , suffix, dpt_,    dir_ );
		    }
		  
		  fill( "met"    , suffix, MET,    dir_ );
		  fill( "mt"     , suffix, MT,     dir_ );
		  fill( "METVsMT", suffix, MET, MT, dir_ );
		  if( isMT40 )
		    {
		      fill( "met"    , suffix+"__"+"MT40", MET, dir_ );
		    }
		  if( isMET20 )
		    {
		      fill( "mt", suffix+"__"+"MET20", MT, dir_ );
		    }
		  
		  fill( "Eproj", suffix+"__ux",      proj_x , dir_ );
		  fill( "Eproj", suffix+"__uy",      proj_y , dir_ );
		  
		  fill( "pt",    suffix,      WPt , dir_ );

		  // jets
		  
		  fill( "pt"  , suffix+"__recoil", jRecoilPt, dir_ );
		  fill( "dphi"  , suffix+"__recoil", jRecoilDphi, dir_ );
		  fill( "dpt"  , suffix+"__recoil", jRecoilDpt, dir_ );
		  fill( "jmult", suffix,  jMult, dir_ );
		  if( MT>40 )
		    {
		      if( WPt>30. )
			fill( "LP", suffix+"__Wpt30",  LP, dir_ );
		      if( WPt>50. )
			fill( "LP", suffix+"__Wpt50",  LP, dir_ );
		    }		  
		}
	    }
	}
    }

  
  if( _verbose )
    {
      cout << "W Candidate " << endl;
      WCand->oneLine( cout );
      WCand->theBase()->oneLine( cout );
    }

  {
    // tuple
    string dir_ = _e.categ();
    dir_ += "/"; 
    dir_ += subdir;
    int run_   = _e.run();
    int event_ = _e.event();
    string fname_(_e.fileName().Data());
    string* fnameptr_ = &fname_;
    tm.setTree( "Wtuple", dir_ );
    tm.add< int >(    "run",   &run_       );
    tm.add< int >(    "event", &event_     );
    tm.add< string* >("fname", &fnameptr_  );  
    tm.add< TBits* >( "idBits", &idBits    );      
    tm.add< float >(  "LPt",   &LPt       );
    tm.add< float >(  "LEta",  &LEta      );
    tm.add< float >(  "LPhi",  &LPhi     );
    tm.add< float >(  "MET",   &MET        );
    tm.add< float >(  "METPhi",&METPhi    );
    tm.add< float >(  "METx",  &proj_x    );
    tm.add< float >(  "METy",  &proj_y    );
    tm.add< float >(  "MT",    &MT        );
    tm.add< float >(  "WPt",   &WPt       );
    tm.add< float >(  "WPhi",  &WPhi      );
    tm.add< int   >(  "jmult", &jMult     );
    tm.add< float >(  "LP",    &LP        );
    tm.add< float >(  "jRecPt",  &jRecoilPt   );
    tm.add< float >(  "jRecDphi",&jRecoilDphi );
    tm.add< float >(  "jRecDpt", &jRecoilDpt  );
    tm.add< vectorFloat* > ( "jEt", &jEt );
    tm.add< vectorFloat* > ( "jEta", &jEta );
    tm.add< vectorFloat* > ( "jDRL", &jDRL );
    tm.add< vectorFloat* > ( "jDPhiL", &jDPhiL );
    tm.add< vectorFloat* > ( "jDPhiNu", &jDPhiNu );
    tm.add< vectorFloat* > ( "jDPhiW", &jDPhiW );
    tm.flush();
  }    

  delete idBits;
  delete jEt;
  delete jEta;
  delete jDRL;
  delete jDPhiL;
  delete jDPhiNu;
  delete jDPhiW;
   
  //  hypothesis is selected
  select( hyp, kmode, kSel );

  return true;
}

void
DibosonAnalysis::fillStatHistograms()
{
  // stat histograms contain exactly the number of events selected at
  // each step of the selection, irrespective of the number of hypotheses 
  // --> they can be used for statistics
  static string prefix("stat");

  // Just for stat histograms, take the "best" Z candidates
  // (=first Z candidate of first Z-hypothesis, if any)
  // This should be unbiased 
  // Mass=zero means no candidate.
  //  const Candidate* bestZ[2] = {0,0};
  float mll[2] = {0.,0.};
  for( int ii=0; ii<2; ii++ )
    {
      const CandList& list = _e.userCandList( ZListName(ii,0) );
      if( list.size()>0 )
	{
	  const Candidate* cand_ = list[0];
	  assert( cand_!=0 );
	  //	  bestZ[ii] = cand_;
	  mll[ii] = cand_->mass();
	}
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
		fill( "mll", str_.Data(), mll[ii], prefix );
	      }
	  }
    }
}

void
DibosonAnalysis::debugPrintouts( ostream& o )
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
      for( HypoStatIterator it = _nHypo[kPresel].begin();
	   it!=_nHypo[kPresel].end(); it++ )
	{
	  print( o, it, kPresel );
	}
    }

  // that's all if no hypothesis is selected
  if( _nHypo[kSel].size()==0 ) return;

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
  if( nSel_<=100 )
    {
      for( HypoStatIterator it = _nHypo[kSel].begin();
	   it!=_nHypo[kSel].end(); it++ )
	{
	  print( o, it, kSel );
	}
    }
}

float
DibosonAnalysis::mT( const Candidate* cand, const Candidate* met )
{
  //   returns  the transverse mass for a candidate and a missing momentum
  TVector2 metP2  =  met->p2();
  TVector2 candP2 = cand->p2();
  
  float mT2 =  2*( candP2.Mod() * metP2.Mod() -  candP2 * metP2 );

  if( mT2<0 ) mT2=0;
  
  return sqrt( mT2 );
}

Candidate*
DibosonAnalysis::correctedMet( HypoListIterator ithyp )
{
  // clone the uncorrected caloMet
  Candidate* corrMet = _e.met("Calo")->clone(); 
  TVector2 vec = corrMet->p2();

  Hypo hyp    = ithyp->first;
  size_t ihyp = hyp.first;
  size_t jhyp = hyp.second;
  
  CandList ZList[2] = 
    { EventManager::e()->userCandList( ZListName(0,ihyp) ),
      EventManager::e()->userCandList( ZListName(1,jhyp) )  };
  
  CandList lList;
  
  // first, get the leptons from the Z candidates
  for( int imode=0; imode<2; imode++ )
    {
      for( size_t iZ=0; iZ<ZList[imode].size(); iZ++ )
 	{
 	  Candidate* ZCand = ZList[imode][iZ];
	  
 	  // add the Z transverse momentum contribution
 	  vec -= ZCand->p2();
	  
 	  for( size_t idau=0; idau<ZCand->nDaughters(); idau++ )
 	    {
 	      lList.push_back( ZCand->daughter(idau) );
 	    }
 	}            
    }

  // then get other leptons in the event
  CandList leptonList[2] = 
    { EventManager::e()->userCandList( leptonListName(0,ihyp) ),
      EventManager::e()->userCandList( leptonListName(1,jhyp) ) };

  for( int imode=0; imode<2; imode++ )
    {
      for( size_t il=0; il<leptonList[imode].size(); il++ )
	{
	  Candidate* lCand = leptonList[imode][il];
	  lList.push_back( leptonList[imode][il] );

	  // add the lepton transverse momentum contribution
	  vec -= lCand->p2();
	}
    }

  // list of Uids of CaloJets associated with leptons
  const CandList& list_cj = _e.jetList( "Calo" );
  vector< int > veto_cjUid;  

  for( size_t il=0; il<lList.size(); il++ )
    {
      Candidate* lCand = lList[il];
      int lUid = lCand->uid();

      int pdgId =  lCand->pdgCode();
      if( pdgId==11 || pdgId==-11 )
	{
	  // remove the transverse energy associated with Ecal clusters

	  // get the map of cluster iso object for this candidate
	  const CandIdIsoMap&  cluMap      = _e.isoMap("electron-EcalCluster");
	  CandIdIsoMapIterator cluIt_lower = cluMap.lower_bound( lUid );
	  CandIdIsoMapIterator cluIt_upper = cluMap.upper_bound( lUid );
	  CandIdIsoMapIterator cluIt;

	  bool dbg = false;
	  if( dbg )
	    {
	      lCand->oneLine( cout );
	    }
	  for( cluIt=cluIt_lower; cluIt!=cluIt_upper; cluIt++ )
	    {
	      float e_   = (*cluIt).second.val;
	      float eta_ = (*cluIt).second.eta;
	      float phi_ = (*cluIt).second.phi;
	      float theta_ = 2*atan( exp( -eta_ ) );
	      float et_    = e_*sin(theta_);
	      TVector2 p2_(0,0);
	      p2_.SetMagPhi( fabs( et_), phi_ );
	      vec += p2_;
	      if( dbg )
		{		  
		  cout << "--> e/et/eta/phi ";
		  cout << e_   << "/";
		  cout << et_  << "/";
		  cout << eta_ << "/";
		  cout << phi_ << endl;
		}
	    }
	}
      else if( pdgId==13 || pdgId==-13 )
	{
	  ConeSelector sel_( *lCand, 0.1  );
	  CandList l_cjList;
	  sel_.getList( list_cj, l_cjList );
	  if( l_cjList.size()>1 )
	    {
	      cout << "WARNING -- More than one CaloJet associated "
		   << endl;
	    }
	  size_t icj=0;
	  for( icj=0; icj<l_cjList.size(); icj++ )
	    {
	      const Candidate* cjCand = l_cjList[icj];
	      veto_cjUid.push_back( cjCand->uid() );
	    }	  
	}
      else
	assert(0);
    }  
  
  //
  // then remove the contribution of CaloTowers associated 
  //  with jets in the veto list
  //
  for( size_t icj=0; icj<list_cj.size(); icj++ )
    {
      const Candidate* cjCand = list_cj[icj];
      int cjUid = cjCand->uid();

      vector<int>::const_iterator vetoIt = 
	find( veto_cjUid.begin(), veto_cjUid.end(), cjUid );
      if( vetoIt != veto_cjUid.end() ) continue;

      // get the map of iso object for this candidate
      const CandIdIsoMap&  cjMap      = _e.isoMap("caloJet-caloTower");
      CandIdIsoMapIterator cjIt_lower = cjMap.lower_bound( cjUid );
      CandIdIsoMapIterator cjIt_upper = cjMap.upper_bound( cjUid );
      CandIdIsoMapIterator cjIt;	
      for( cjIt=cjIt_lower; cjIt!=cjIt_upper; cjIt++ )
	{
	  float et_  = (*cjIt).second.val;
	  float eta_ = (*cjIt).second.eta;
	  float phi_ = (*cjIt).second.phi;
	  TVector3 p3_(0,0,0);
	  p3_.SetPtEtaPhi( fabs( et_), eta_, phi_ );
	  vec += p3_.XYvector();
	}
    }
  
  corrMet -> setPxPyPz( vec.X(), vec.Y(), 0 );
  return corrMet;
}

bool
DibosonAnalysis::passHlt()
{
  //
  // This is data
  //
  bool fired = false;

  // to be fixed
  if( 1 ) return fired;
  
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
    {
      if( _e.isFired("HLT_Ele15_LW_L1R") ||
	  _e.isFired("HLT_Ele15_SW_L1R") ||
	  _e.isFired("HLT_Ele15_SW_CaloEleId_L1R") ||
	  _e.isFired("HLT_Ele17_SW_CaloEleId_L1R") ||
	  _e.isFired("HLT_Ele17_SW_TightEleId_L1R") ||
	  _e.isFired("HLT_Ele17_SW_TighterEleIdIsol_L1R_v2") ||
	  _e.isFired("HLT_Ele17_SW_TighterEleIdIsol_L1R_v3") ||
	  _e.isFired("HLT_Mu9") ||
	  _e.isFired("HLT_Mu15_v1") ||
	  _e.isFired("HLT_Mu15_v2") )      
	fired=true;
    }
  
  return fired;
}



bool
DibosonAnalysis::isPhotonIDIso(Candidate* ph) {

  CandInfo* info_ = ph->info();

  //SC
  Candidate* SC = _e.getSecond("photon-SuperCluster", ph );
  CandInfo* infoSC = SC->info();

  float sigphi = infoSC->getFloat("phiWidth");
  float sigieie = info_->getFloat("sigmaIetaIeta");
  float HoE = info_->getFloat("hadronicOverEm");
  float trkIso = info_->getFloat("trkSumPtHollowConeDR04");
  float ecalIso = info_->getFloat("ecalRecHitSumEtConeDR04");
  float hcalIso = info_->getFloat("hcalTowerSumEtConeDR04");
  bool pixel = info_->getBool("hasPixelSeed");

  if(pixel) return false;

  if( fabs(ph->eta() ) < 1.44 ) {
      if(sigieie > 0.01 ) return false;
      if(sigphi < 0.001 ) return false;
      if(HoE > 0.05) return false;
      if(trkIso/ph->pt() > 0.05) return false;
      if(ecalIso/ph->pt() > 0.05) return false;
      if(hcalIso/ph->pt() > 0.05) return false;

    }
    else if  (fabs(ph->eta() ) > 1.56 && fabs(ph->eta()) < 2.5 ) {
      if(sigieie > 0.01 ) return false;
      if(HoE > 0.05) return false;
      if(trkIso/ph->pt() > 0.03) return false;
      if(ecalIso/ph->pt() > 0.03) return false;
      if(hcalIso/ph->pt() > 0.03) return false;
    }
    else { return false;}
    
  return true;
}


