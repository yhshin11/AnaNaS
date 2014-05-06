#include <unistd.h> 
#include <iostream>
#include <sstream>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <typeinfo>
#include <string>
#include <cmath>
#include <algorithm>
#include <functional>
#include <vector>
#include <map>

#include <TFile.h>
#include <TBranch.h>
#include <TKey.h>
#include <TTree.h>
#include <TString.h>
#include <TBits.h>
#include <Math/VectorUtil.h>
#include <TMath.h>
#include <TVector3.h>
#include <TLorentzVector.h>

#include "Analysis/ghmzz/GhmTupleAnalysis.hh"

#include "Analysis/core/Sample.hh"
#include "Analysis/core/Candidate.hh"
#include "Analysis/core/CandInfo.hh"
#include "Analysis/core/SampleAnalysis.hh"
#include "Analysis/tools/CandUtil.hh"
#include "Analysis/tools/PrintTree.hh"
#include "Analysis/utils/Constants.hh"
#include "Analysis/utils/Config.hh"
#include "Analysis/utils/KineUtils.hh"
#include "Analysis/utils/RooUtils.hh"
#include "Analysis/utils/HistoManager.hh"
#include "Analysis/utils/TupleManager.hh"
#include "Analysis/utils/DatasetManager.hh"
#include "Analysis/utils/CutUtils.hh"
#include "Analysis/utils/StringUtils.hh"
#include "Analysis/selectors/ConeSelector.hh"

ClassImp( GhmTupleAnalysis )

GhmTupleAnalysis::GhmTupleAnalysis( const char* sampleSetName, int n ) : 
  BaseTupleAnalysis( sampleSetName, n )
{
}

void
GhmTupleAnalysis::init()
{
  defineTemplate( "mll", 14800, 20., 1500  );
  defineTemplate( "pt",   4000, 0., 1000.  );
  defineTemplate( "y",    6000, -3., 3.    );
  defineTemplate( "h",    3600, -90., 90.  );

//   tm_out.setTree("secondTuple","ghmDir");
//   tm_out.declare<int>( "run", &irun );
//   tm_out.declare<int>( "event", &ievent );
//   tm_out.declare<float>( "mll", &_float["mll"] );

}

bool
GhmTupleAnalysis::ThirdLeptonSelector::accept( const Candidate& cand ) const
{
  float pt_ = cand.pt();
  if( pt_<20 )   return false;
  bool isTight = cand.info()->getBool("isVeryTight");
  if( !isTight ) return false;
  return true;
}

bool
GhmTupleAnalysis::ZSelector::accept( const Candidate& cand ) const
{
  float mll_  = cand.mass();
  bool isInWindow = mll_>=60 && mll_<=120;
  bool isTight = cand.info()->getBool("isTight");
  return ( isInWindow && isTight );
}

bool 
GhmTupleAnalysis::ZCompare::operator()( const Candidate* Z1, const Candidate* Z2 ) const
{
  const float mZ = Constants::Z0_m;
  float m1 = Z1->mass();
  float m2 = Z2->mass();
  
  float dm1 = fabs( m1-mZ );
  float dm2 = fabs( m2-mZ );

  bool out = dm1<dm2;
  return out;
}

bool 
GhmTupleAnalysis::analyzeEvent()
{
  _verbose = false;
  //  static int N_= 0;
  //  if( ++N_<10 ) _verbose = true;

  bool isSelected=false;

  // loop over Z bosons at vertices
  for( map<int,CandList>::iterator itVtx=ZAtVertex.begin(); 
       itVtx != ZAtVertex.end(); ++itVtx )
    {
      // get the vertex
      int ivtx_ = itVtx->first;
      Vertex* vtx_ = vertices[ ivtx_ ];
      assert( vtx_!=0 );

      // get the Z list at vertex
      CandList ZListAll_ = itVtx->second;

      // number of Z candidates (including overlaps)
      int nZAll_ = (int) ZListAll_.size();
      assert( nZAll_>0 );

      // get the number of leptons at vertex (including overlaps)
      CandList leptonAtVertexAll_ = leptonAtVertex[ivtx_];
      int nLepAll_ = (int) leptonAtVertexAll_.size(); 
      assert( nLepAll_>=2 );

      // get the number of jets at vertex
      CandList jetAtVertex_ = jetAtVertex[ivtx_];
      int nJet_ = (int) jetAtVertex_.size();

      // dispatch the list into lists of non-overlapping candidates
      vector< CandList > hypLists_;
      CandUtil::dispatchList( ZListAll_, hypLists_ );

      // number of Z configuration hypotheses at that vertex
      int nHyp_ = (int) hypLists_.size();
      
      string hdir = dirName;
      hdir += "/";
      hdir += _sampleName;

      // loop over the hypotheses
      for( size_t iHyp=0; iHyp<hypLists_.size(); iHyp++ )
	{	  
	  int iHyp_ = (int) iHyp;
	
	  // list of non-overlapping Z candidates
	  CandList ZList_ = hypLists_[iHyp];
	  int nZ_ = (int) ZList_.size();
	  assert( nZ_>0 );

	  // list of non-overlapping leptons
	  CandList leptonAtVertex_;
	  CandUtil::pruneList( leptonAtVertexAll_, leptonAtVertex_, ZList_ );
	  int nLep_ = (int) leptonAtVertex_.size();

	  // list of selected lepton
	  CandList selectedLeptonAtVertex_;
	  ThirdLeptonSelector lepSel_;
	  lepSel_.getSortedList( leptonAtVertex_, selectedLeptonAtVertex_ );
	  int nSelLep_ = (int) selectedLeptonAtVertex_.size();

	  // just for the check, list of leptons from Z bosons
	  CandList leptonFromZ_;
	  CandUtil::pruneList( leptonAtVertexAll_, leptonFromZ_, leptonAtVertex_ );
	  int nLepFromZ_ = (int) leptonFromZ_.size();
	  assert( nLepFromZ_ == 2*nZ_ );

	  if( nZ_==1 )
	    {
	      Candidate* Z_ = ZList_[0];
	
	      CandInfo*  Zinfo_ = Z_->info();
	      bool isTight = Zinfo_->getBool("isTight");
	      
	      float mll  = Z_->mass();
	      int lmode  = Zinfo_->getInt(   "mode" );
	      float qTll = Zinfo_->getFloat( "qTll" );
	      float yll  = Zinfo_->getFloat( "yll"  );
	      float hll  = Zinfo_->getFloat( "hll"  );

	      tm_out.setTree("ZTuple", categ );
	      tm_out.declare<int>( "run", &irun );
	      tm_out.declare<int>( "event", &ievent );
	      tm_out.declare<float>( "mll", &mll );
	      tm_out.declare<float>( "qTll", &qTll );
	      tm_out.declare<float>( "yll", &yll );
	      tm_out.declare<float>( "hll", &hll );
	      tm_out.declare<bool>( "isTight", &isTight );
	      tm_out.declare<int>( "nLep", &nLep_ );
	      tm_out.declare<int>( "nSelLep", &nSelLep_ );
	      tm_out.declare<int>( "nJet", &nJet_ );
	      tm_out.declare<int>( "nHyp", &nHyp_ );
	      tm_out.declare<int>( "hyp", &iHyp_ );
	      tm_out.fill();

	      string hname =  ZllName[lmode];
	      if( !isTight ) hname+="_fail";
	  
	      // OK
	      hist( "mll", hname, hdir, "" )->Fill( mll,  1. );
	      hist( "pt",  hname, hdir, "" )->Fill( qTll, 1. );
	      hist( "y",   hname, hdir, "" )->Fill( yll,  1. );
	      hist( "h",   hname, hdir, "" )->Fill( hll,  1. );

	      ZSelector ZSel_;

	      if( ZSel_.accept(*Z_) && nSelLep_>0 )
		{
		  print();

		  Candidate* thirdLepton_ = selectedLeptonAtVertex_[0];
		  assert( thirdLepton_!=0 );

		  WCandidate Wcand_( thirdLepton_, MetCand );
		  Candidate* W_ = Wcand_.W();
		  assert( W_!=0 );

		  float DRmax_( 10.0 );
		  ConeSelector coneSel_( *thirdLepton_, DRmax_ );

		  CandList DRsortedLeptonFromZ_;
		  coneSel_.getSortedList( leptonFromZ_, DRsortedLeptonFromZ_, Candidate::kSortDefault );
		  cout << "Distance to leptons from Z" << endl;
		  for( size_t ii=0; ii<DRsortedLeptonFromZ_.size(); ii++ )
		    {
		      Candidate* c_=DRsortedLeptonFromZ_[ii];
		      cout << "DR=" << coneSel_.dr( *c_ ) << " -- ";
		      c_->oneLine(cout);
		    }
		  
		  CandList DRsortedLeptonAtVertex_;
		  coneSel_.getSortedList( leptonAtVertex_, DRsortedLeptonAtVertex_, Candidate::kSortDefault );
		  if( DRsortedLeptonAtVertex_.size()>0 )
		    {
		      cout << "Distance to other leptons" << endl;
		      for( size_t ii=0; ii<DRsortedLeptonAtVertex_.size(); ii++ )
			{
			  Candidate* c_=DRsortedLeptonAtVertex_[ii];
			  cout << "DR=" << coneSel_.dr( *c_ ) << " -- ";
			  c_->oneLine(cout);
			}
		    } 

		}

// 	      if( nLep_==0 )
// 		{
// 		}
// 	      else
// 		{
// 		}
	    }
	  else 
	    {
	      continue;
	      // sort the list of Z bosons
	      sort( ZList_.begin(), ZList_.end(), ZCompare() );
	      if( nZ_==2 )
		{
		  cout << "\nRun/Event " << irun << "/" << ievent << endl;
		  cout << "!!! HYPOTHESIS WITH 2 Z CANDIDATES !!!" << endl;		  
		}
	      else
		{
		  cout << "\nRun/Event " << irun << "/" << ievent << endl;
		  cout << "!!! HYPOTHESIS WITH MORE THAN 2 Z CANDIDATES !!!" << endl;
		}

	      //	      print();

	      // take only the first 2 Z bosons
	      tm_out.setTree("ZZTuple", categ );
	      tm_out.declare<int>( "run", &irun );
	      tm_out.declare<int>( "event", &ievent );
	      for( size_t iZZ=0; iZZ<2; iZZ++ )
		{
		  Candidate* Z_ = ZList_[iZZ];

		  CandInfo*  Zinfo_ = Z_->info();

		  float mll = Z_->mass();
		  bool isTight = Zinfo_->getBool("isTight");
		  int lmode  = Zinfo_->getInt(   "mode" );
		  float qTll = Zinfo_->getFloat( "qTll" );
		  float yll  = Zinfo_->getFloat( "yll"  );
		  float hll  = Zinfo_->getFloat( "hll"  );
	      
		  {
		    string istr_;

		    istr_ = StringUtils::indexedString( "mll", iZZ );		      
		    _float[istr_] = mll;
		    tm_out.declare<float>( istr_, &_float[istr_] );

		    istr_ = StringUtils::indexedString( "qTll", iZZ );		      
		    _float[istr_] = qTll;
		    tm_out.declare<float>( istr_, &_float[istr_] );

		    istr_ = StringUtils::indexedString( "yll", iZZ );		      
		    _float[istr_] = yll;
		    tm_out.declare<float>( istr_, &_float[istr_] );

		    istr_ = StringUtils::indexedString( "hll", iZZ );		      
		    _float[istr_] = hll;
		    tm_out.declare<float>( istr_, &_float[istr_] );

		    istr_ = StringUtils::indexedString( "mode", iZZ );		      
		    _int[istr_] = lmode;
		    tm_out.declare<int>( istr_, &_int[istr_] );

		    istr_ = StringUtils::indexedString( "isTight", iZZ );		      
		    _bool[istr_] = isTight;
		    tm_out.declare<bool>( istr_, &_bool[istr_] );
		    
		    
		    
		    // 		      tm_out.declare<float>( "qTll", &_float["qTll"] );
		    // 		      tm_out.declare<float>( "yll", &yll );
		    // 		      tm_out.declare<float>( "hll", &hll );		      
		    // 		      tm_out.declare<int>( "mode", &lmode );
		    // 		      tm_out.declare<bool>( "isTight", &isTight );
		  }
		}
	      tm_out.declare<int>( "nZ", &nZ_ );
	      tm_out.declare<int>( "nLep", &nLep_ );
	      tm_out.declare<int>( "nJet", &nJet_ );
	      tm_out.declare<int>( "nHyp", &nHyp_ );
	      tm_out.declare<int>( "hyp", &iHyp_ );
	      tm_out.fill();

	    }
	}

//       size_t iZ=0;
//       //      for( size_t iZ=0; iZ<nZ_; iZ++ )
//       {
// 	//	  if( iZ!=0 ) continue;
	
// 	Candidate* Zcand_ = ZList_[iZ];
	
// 	CandInfo*  Zinfo_ = Zcand_->info();
// 	bool isTight = Zinfo_->getBool("isTight");
	  
// 	// simple event selection...
// 	if( isTight ) isSelected = true;
	  
// 	float mll  = Zcand_->mass();
// 	int lmode  = Zinfo_->getInt(   "mode" );
// 	float qTll = Zinfo_->getFloat( "qTll" );
// 	float yll  = Zinfo_->getFloat( "yll"  );
// 	float hll  = Zinfo_->getFloat( "hll"  );
	
// 	tm_out.setTree("ZllTuple", categ );
// 	tm_out.declare<int>( "run", &irun );
// 	tm_out.declare<int>( "event", &ievent );
// 	tm_out.declare<float>( "mll", &mll );
// 	tm_out.declare<float>( "qTll", &qTll );
// 	tm_out.declare<float>( "yll", &yll );
// 	tm_out.declare<float>( "hll", &hll );
// 	tm_out.declare<bool>( "isTight", &isTight );
// 	tm_out.declare<int>( "nLep", &nLep_ );
// 	tm_out.declare<int>( "nJet", &nJet_ );
// 	tm_out.fill();


// 	string hname =  ZllName[lmode];
// 	if( !isTight ) hname+="_fail";
	  
// 	// OK
// 	hist( "mll", hname, hdir, "" )->Fill( mll,  1. );
// 	hist( "pt",  hname, hdir, "" )->Fill( qTll, 1. );
// 	hist( "y",   hname, hdir, "" )->Fill( yll,  1. );
// 	hist( "h",   hname, hdir, "" )->Fill( hll,  1. );
//       }
      
//     }


//   // fill simple histograms
//   for( int lmode=0; lmode<2; lmode++ )
//     {
//       string hdir = dirName;
//       hdir += "/";
//       hdir += _sampleName;
//       for( size_t iZ=0; iZ<Zll[lmode].size(); iZ++ )
// 	{
// 	  if( iZ!=0 ) continue;
// 	  Candidate* Zcand_ = Zll[lmode][iZ];
// 	  int vtx_ = Zcand_->vertexIndex();
// 	  _int["nLep"] = (int) leptonAtVertex[vtx_].size();
// 	  _int["nJet"] = (int) jetAtVertex[vtx_].size();

// 	  CandInfo*  Zinfo_ = Zcand_->info();
// 	  _bool["isTight"] = Zinfo_->getBool("isTight");
	  
// 	  // simple event selection...
// 	  if( _bool["isTight"] ) isSelected = true;
	  
// 	  _float["mll"]  = Zcand_->mass();
// 	  _float["qTll"] = Zinfo_->getFloat( "qTll" );
// 	  _float["yll"]  = Zinfo_->getFloat( "yll"  );
// 	  _float["hll"]  = Zinfo_->getFloat( "hll"  );

// 	  tm_out.setTree("ZllTuple", categ );
// 	  tm_out.declare<int>( "run", &irun );
// 	  tm_out.declare<int>( "event", &ievent );
// 	  tm_out.declare<float>( "mll", &_float["mll"] );
// 	  tm_out.declare<float>( "qTll", &_float["qTll"] );
// 	  tm_out.declare<float>( "yll", &_float["yll"] );
// 	  tm_out.declare<float>( "hll", &_float["hll"] );
// 	  tm_out.declare<bool>( "isTight", &_bool["isTight"] );
// 	  tm_out.declare<int>( "nLep", &_int["nLep"] );
// 	  tm_out.declare<int>( "nJet", &_int["nJet"] );
// 	  tm_out.fill();

// 	  string hname =  ZllName[lmode];
// 	  if( !_bool["isTight"] ) hname+="_fail";
	  
// 	  // OK
// 	  hist( "mll", hname, hdir, "" )->Fill( _float["mll"],  1. );
// 	  hist( "pt",  hname, hdir, "" )->Fill( _float["qTll"], 1. );
// 	  hist( "y",   hname, hdir, "" )->Fill( _float["yll"],  1. );
// 	  hist( "h",   hname, hdir, "" )->Fill( _float["hll"],  1. );
// 	}
      
    }

  return isSelected;
}

void
GhmTupleAnalysis::finalize()
{
}


// 	  if( hbool.size()==0 )
// 	    {
// 	      // bool histograms
// 	      hbool.push_back(
// 			      HBool(
// 				    hist(
// 					 "mll", "mll", dirName, "" 
// 					 ),            // the histogram
// 				    &_float["Zee_mll"],    // pointer to the variable
// 				    &_bool["Zee_isTight"], // pointer to the boolean
// 				    true,              // underflow bin?
// 				    true               // overflow bin?
// 				    )
// 			      );		
// 	    }
// 	  float w_ = 1;
// 	  for( size_t ih=0; ih<hbool.size(); ih++ )
// 	    {
// 	      hbool[ih].fill( w_ );
// 	    }

WCandidate::WCandidate( Candidate* lepton, Candidate* etmiss )
  : _l( lepton ), _n( etmiss ), _ln(0)
{
  assert( _l!=0 && _n!=0 );
  assert( etmiss->isTransverse() );
  
  _ln = Candidate::create( _l, _n );

  // the vertex
  Vertex* vtx_ = _l->vertex();

  // the lepton, neutrino and PDG code
  int lid_ = _l->pdgCode();
  int s1_ = lid_>0 ?  1 : -1;
  int s2_ = -1*s1_;
  int nid_ = s2_*(s1_*lid_+1);
  int Wid_ = s2_*24;

  // the lepton momentum
  _p_l   = _l->p();

  // the lepton transvere momentum
  _pT_l_vec = _l->p2();
  _pT_l  = _l->pt();

  // the lepton longitudinal momentum
  _pL_l  = _l->pz();

  // the neutrino  transverse momentum
  _pT_n_vec = _n->p2();
  _pT_n = _n->pt();

  // the true W mass
  _m = Constants::W_m;
  _m2 = pow(_m,2);

  // the transverse mass
  _mT2 = 2*( _pT_l*_pT_n - _pT_l_vec*_pT_n_vec );
  assert( _mT2>=0 );
  _mT = sqrt( _mT2 );

  // find the two W solutions
  _A = _m2 + 2 * _pT_l_vec * _pT_n_vec;
  float den_ = 2*pow(_pT_l,2);
  if( _m2>_mT2 )
    {
      _B = ( _m2 - _mT2 )*( _m2 + 2*( _pT_l*_pT_n + _pT_l_vec*_pT_n_vec ) );
      assert( _B>=0 );
      float sqrtB_ = sqrt(_B);
      // the two solutions for the longitudinal momentum
      float pL_n_p_ = ( _pL_l*_A + _p_l*sqrtB_ )/den_;
      float pL_n_m_ = ( _pL_l*_A - _p_l*sqrtB_ )/den_;
      // choose the solution with the smallest abs value
      
      vector< float > pL_n_(2);
      pL_n_[0] = pL_n_p_;
      pL_n_[1] = pL_n_m_;
      if( fabs(pL_n_p_)>fabs(pL_n_m_) )
	{
	  pL_n_[0] = pL_n_m_;
	  pL_n_[1] = pL_n_p_;
	}
      
      // loop over solutions
      for( size_t ii=0; ii<2; ii++ )
	{
	  // create the two neutrino candidates
	  TVector3 p_n_( _pT_n_vec.X(), _pT_n_vec.Y(), pL_n_[ii] );
	  Candidate* n_ = Candidate::create( p_n_, nid_, vtx_ );
      
	  // the W solution
	  Candidate* W_ =  Candidate::create( _l, n_ );
	  CandInfo* infoW_ = CandInfo::create();
	  W_->setInfo( infoW_ );

	  float mW_ = W_->mass();
	  infoW_->setFloat("mW",mW_ );
	  infoW_->setFloat("mT",_mT );
	  infoW_->setFloat("A",_A );
	  infoW_->setFloat("B",_B );

	  W_->setPdgCode( Wid_ );

	  _Wsol.push_back( W_ );
	}
    }
  else
    {
      _B = 0;
      float pL_n_ = _pL_l*_A/den_;
      TVector3 p_n_( _pT_n_vec.X(), _pT_n_vec.Y(), pL_n_ );
      Candidate* n_ = Candidate::create( p_n_, nid_, vtx_ );

      // the W solution
      Candidate* W_ =  Candidate::create( _l, n_ );
      CandInfo* infoW_ = CandInfo::create();
      W_->setInfo( infoW_ );
      
      float mW_ = W_->mass();
      infoW_->setFloat("mW",mW_ );
      infoW_->setFloat("mT",_mT );
      infoW_->setFloat("A",_A );
      infoW_->setFloat("B",_B );
      
      W_->setPdgCode( Wid_ );

      _Wsol.push_back( W_ );
    }
}

//		  cout << "Hypothesis: WZ candidate " << endl;
//		  cout << "Z candidate -- ";
//		  Z_->oneLine(cout);
//		  cout << "3d lepton -- ";
//		  thirdLepton_->oneLine(cout);
//		  W_->oneLine(cout);
