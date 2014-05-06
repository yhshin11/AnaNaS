#include "Analysis/core/SampleAnalysis.hh"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;

#include <TDatime.h>

#include "Analysis/utils/Constants.hh"
#include "Analysis/utils/Config.hh"
#include "Analysis/tools/CandUtil.hh"
#include "Analysis/core/EventManager.hh"
#include "Analysis/tools/PrintTree.hh"
#include "Analysis/selectors/ElectronSelector.hh"
#include "Analysis/selectors/VBTFElectronSelector.hh"
#include "Analysis/selectors/MuonSelector.hh"

ClassImp( SampleAnalysis )

// DY analysis binning: 15-20-30-40-50-60-76-86-96-106-120-150-200-600
//        large window            |++++++++++++++++++++++++++++++|   
//              window                  |----------------| 
//                veto                     |xxxxxxxxx|
//          tight veto                        |oo| 
//
const float SampleAnalysis::Z0_fullDY[2]      = { 20., 10000. };
const float SampleAnalysis::Z0_largeWindow[2] = { 40., 200. };
const float SampleAnalysis::Z0_window[2]      = { 60., 120. };
const float SampleAnalysis::Z0_veto[2]        = { 76., 106. };
const float SampleAnalysis::Z0_tightVeto[2]   = { 86., 96. };

SampleAnalysis::SampleAnalysis( const string & name, Sample& sample, const std::string & collectionFileName  ) 
  : _name(name), _sample(sample), _e( *(new EventManager( sample.filenames, collectionFileName.c_str() )) ),
    _verbose(false), _nDebug(0), _hlt(false), 
    _ievt(0), _outputFile(0), 
    _isSelected(false)
{

  assert( 0 );
}

SampleAnalysis::SampleAnalysis( const string & name, Sample& sample, EventManager& eventManager ) 
  : _name(name), _sample(sample), 
    //    _e( _sample.filenames, collectionFileName.c_str() ),
    _e( eventManager ),
    _verbose(false), _nDebug(0), _hlt(false), 
    _ievt(0), _outputFile(0), 
    _isSelected(false)
{
  cout << endl;
  cout << "\t-------------------------------------------" << endl; 
  cout << "\t--- Welcome to the Analysis of Sample  ";
  cout << _sample.name()  << endl;
  cout << "\t-------------------------------------------" << endl; 
}

void
SampleAnalysis::analyze( int nevt, int nskip, bool random )
{
  // book histograms
  bookHistograms_();

  // event loop
  int iskip=0;

  while( _e.nextEvent() )
    {
      if( iskip++<nskip ) continue;

      // event initialisations
      initEvent_();

      // stop event loop when number of events reached
      if( nevt>0 && _ievt>nevt ) break;
      //*****
      //****      if( _ievt%10000!=0 ) continue;

      // event analysis
      _isSelected = analyzeEvent_();

      // update counters
      updateCounters_();
    }

  // write histograms
  writeHistograms_();

  // final printings
  printStats_();
}

bool
SampleAnalysis::nextEvent()
{
  _e.nextEvent();
  initEvent_();
  return sameEvent();
}

bool
SampleAnalysis::sameEvent()
{
  bool out_ =  analyzeEvent_();
  print( cout );
  return out_;
}

void
SampleAnalysis::initEvent_()
{
  // increment event number
  _ievt++;
  
  // verbosity
  _verbose = false; if( _ievt<=_nDebug ) _verbose = true;

}

bool
SampleAnalysis::analyzeEvent_()
{
  if( _verbose ) cout << "\nStart of event " << _ievt << endl;

  if( _verbose )
    {
      PrintTree prtTree;
      Candidate* genCand_ = _e.decayTree();
      if( genCand_!=0 )
	{
	  cout << prtTree.Print( *genCand_ );
	  cout << "Event Category (long):  " << _e.categ_long() << endl;
	  cout << "Event Category (short): " << _e.categ() << endl;

	  CandList ewBosons;
	  CandUtil::get(  24, genCand_, ewBosons );
	  CandUtil::get( -24, genCand_, ewBosons );
	  CandUtil::get(  23, genCand_, ewBosons );

	  CandList lepFromBoson;
	  CandList qrkFromBoson;
	  
	  cout << "Electroweak bosons (MC truth)" << endl;
	  for( size_t ibos=0; ibos<ewBosons.size(); ibos++ )
	    {
	      Candidate* bos_ = ewBosons[ibos];
	      bos_->oneLine( cout );
	      //	      prtTree.Print( *bos_ );

	      for( int ii=11; ii<16; ii+=2 )
		{
		  CandUtil::get(  ii, bos_, lepFromBoson );
		  CandUtil::get( -ii, bos_, lepFromBoson );
		}

	      // search for neutrinos (12:nu_e; 14:nu_mu; 16:nu_tau)
	      for( int ii=12; ii<17; ii+=2 )
		{
		  CandUtil::get(  ii, bos_, lepFromBoson );
		  CandUtil::get( -ii, bos_, lepFromBoson );
		}

	      // search for quarks (1:d, 2:u, 3:s, 4:c, 5:b)
	      for( int ii=1; ii<6; ii++ )
		{
		  CandUtil::get(  ii, bos_, qrkFromBoson );
		  CandUtil::get( -ii, bos_, qrkFromBoson );
		}
	    }

	  cout << "Leptons from electroweak bosons (MC truth)" << endl;
	  for( size_t ilep=0; ilep<lepFromBoson.size(); ilep++ )
	    {
	      Candidate* lep_ = lepFromBoson[ilep];
	      lep_->oneLine( cout );
	    }	      
	}
    }
      
  // output flag
  bool selected_(false);
  
  string prefix("evt");
  string suffix("");
  //   for( int ii=0; ii<EventManager::kNType; ii++ )
  //     {
  //       fill( "n", EventManager::objectName[ii],    
  // 	    _e.candList(ii).size()   , prefix );
  //     }
  
  //   for( int ii=0; ii<EventManager::kNJet; ii++ )
  //     {
  //       fill( "n", EventManager::jetName[ii],    
  // 	    _e.jetList(ii).size()   , prefix );
  //     }

  fillMultHists( prefix, suffix );

  // trigger decision
  _trigHlt = passHlt();
  if( _trigHlt )
    {
      suffix = "_hlt";
      fillMultHists( prefix, suffix );
    }

  // analyze event if hlt trigged or ignored
  if( _trigHlt || !_hlt )
    {
      // subclass: analyze event
      selected_ = analyzeEvent();
    }
      
  // counters after selection
  if( selected_ )
    {
      suffix = "_sel";
      fillMultHists( prefix, suffix );
      //      string prefix("evt");
      //       for( int ii=0; ii<EventManager::kNType; ii++ )
      // 	{
      // 	  fill( "n", EventManager::objectName[ii]+"_sel",    
      // 		_e.candList(ii).size()   , prefix );
      // 	}
      //       for( int ii=0; ii<EventManager::kNJet; ii++ )
      // 	{
      // 	  fill( "n", EventManager::jetName[ii]+"_sel",    
      // 		_e.jetList(ii).size()   , prefix );
      // 	}
    }  
  if( _verbose ) cout << "End of event " << _ievt << endl;

  return selected_;
}

bool
SampleAnalysis::passHlt()
{
  // simple version -- should be redone in subclasses
  for( size_t ihlt=0; ihlt<_hltLines.size(); ihlt++ )
    {
      if( _e.isFired( _hltLines[ihlt] ) )  return true;
    }
  return false;
}

void
SampleAnalysis::updateCounters_()
{
  string _categ = _e.categ(); 
//   static int ipass=0;
//   if( ++ipass==0 )
//     {
//       _nAnalyzed["all"]=0;
//       _nAnalyzed[_categ]=0;      
//       _nHlt["all"]=0;
//       _nHlt[_categ]=0;      
//       _nSelected["all"]=0;
//       _nSelected[_categ]=0;            
//     }
  _nAnalyzed["all"]++;
  _nAnalyzed[_categ]++;      
  if( _trigHlt ) 
    {
      _nHlt["all"]++;
      _nHlt[_categ]++;
    }
  if( _isSelected ) 
    {
      _nSelected["all"]++;
      _nSelected[_categ]++;
    }
}

void
SampleAnalysis::printStats_()
{
  assert( _nAnalyzed.count("all")!=0 );
  int nTot_ = _nAnalyzed["all"];
  assert( nTot_!=0 );
  
  for( map< string, int >::iterator it_=_nAnalyzed.begin();
       it_!=_nAnalyzed.end(); it_++ )
    {
      const string & str_ = it_->first;
      int nAnalyzed_ = it_->second;
      float frac_ = (1.*nAnalyzed_)/nTot_;
      if( nAnalyzed_<=0 ) continue;
      int nHlt_      = _nHlt[str_];
      int nSelected_ = _nSelected[str_];
      double p_   = (1.*nSelected_)/nAnalyzed_;
      double var_ = nAnalyzed_*p_*(1-p_);
      double sig_ = sqrt(var_)/nAnalyzed_;
      printf( "%6.2f%1s - Sel/Hlt/Tot[%-12s]=%6d/%6d/%6d --> (%6.2f+/-%-5.2f)%1s\n", frac_*100, "%", str_.c_str(), nSelected_, nHlt_, nAnalyzed_, p_*100, sig_*100, "%" );
    }
}

void
SampleAnalysis::bookHistograms_()
{
  cout << "Histogram booking" << endl;

  //
  // output file for histograms
  //
  FILE *test_;

  TString dirname_ = Config::rootPath + _name;
  test_=fopen( dirname_.Data(), "r" );
  if( test_==0 )
    {
      TString command_ = TString("mkdir ") + dirname_; 
      assert( system( command_.Data() )==0 );
    }
  else
    fclose( test_ );

  TString ofilename_ = dirname_ + "/" + _sample.name() + TString(".root");
  test_=fopen( ofilename_.Data(), "r" );
  if( test_!=0 )
    {
      fclose( test_ );
      TDatime datime_;
      cout << "File " << ofilename_ << " already exists, save it." << endl;;
      TString command_ = TString("mv ") + ofilename_ + " " 
	+ ofilename_ + "_"; 
      command_ += datime_.Get();
      assert( system( command_ )==0 );
    }
  cout << "Will write histograms in output file " << ofilename_;
  _outputFile  = new TFile( ofilename_.Data(), "RECREATE" );

  hm.setFile( _outputFile );
  tm.setFile( _outputFile, TupleManager::kWrite );

  //
  // template for object counter
  //
  int n_=1000;
  defineTemplate(       "n",  n_+1,  -0.5, n_+0.5 );
  n_=1;
  defineTemplate(       "bool",  n_+1,  -0.5, n_+0.5 );
  n_=10;
  defineTemplate(       "n10",  n_+1,  -0.5, n_+0.5 );
  n_=20;
  defineTemplate(       "n20",  n_+1,  -0.5, n_+0.5 );

  //counter for number of processed events
  n_=1;
  defineTemplate(  "nEvtProc",  n_, -0.5, 0.5 );

  //
  // define general templates
  //
  defineTemplate(     "mll", 5920, 20., 1500. );
  defineTemplate(      "pt", 4000,  0., 1000. );
  defineTemplate(     "dpt", 400, -100., 100. );
  defineTemplate(     "met", 4000,  0., 1000. );
  defineTemplate(      "mt", 4000,  0., 1000. );
  defineTemplate( "met_mll",  360, 20., 200., 200,  0., 400. );
  defineTemplate(     "eta", 3000, -3., 3.    );
  defineTemplate(       "y", 3000, -3., 3.    );
  defineTemplate(     "phi",  720,    0., 360. );
  defineTemplate(    "dphi",  720, -180., 180. );
  defineTemplate(      "dr",  200,   -1., 1.   );
  defineTemplate(   "theta",  360, 0., 180. );
  defineTemplate( "eta_phi",   90, 0., 360., 50, -5., 5. );
  defineTemplate(     "hel",  180, 0, 90.   );

  // templates for electrons ID
  defineTemplate(    "elID", 12, -0.5, 5.5    );
  defineTemplate(   "fBrem", 200, -1., 1.    );
  defineTemplate(     "EOP", 200, 0., 5.    );
  defineTemplate(     "HOE", 200, 0., 0.1    );
  defineTemplate(  "dPhiIn", 200, -0.2, 0.2 );
  defineTemplate(  "dEtaIn", 200, -0.05, 0.05 );
  defineTemplate(  "sigiEtaiEta", 200, 0, 0.1 );
  defineTemplate(    "nTrk", 11, -0.5, 10.5 );
  defineTemplate("sumOfPt2", 200, 0., 200. );
  defineTemplate( "sumOfPt", 200, 0., 10.  );
  defineTemplate(     "iso", 200, 0., 1. );
  defineTemplate(  "dPhiIn_dEtaIn", 50, -0.05, 0.05, 50, -0.05, 0.05 );
 
  // templates for muon ID
  defineTemplate(    "d0", 600, -0.3, 0.3    );
  defineTemplate(    "normChi2", 200, 0, 20 );                                 
 
  //
  // subclass: book histograms
  // 
  bookHistograms();
}

void
SampleAnalysis::writeHistograms_()
{
  _outputFile->cd();

  // subclass: write histograms
  writeHistograms();

  //produce and save the default normalization counter
  fill( "nEvtProc", "run", 0, _e.nEventProc() , "all" );
  cout<<" ====================>>>>>>>>>>>> "<<_e.nEventProc()<<endl;

  // save histograms managed by the histo manager
  hm.save();

  // save ntuples managed by the tuple manager
  tm.save();

  // close file
  if( _outputFile!=0 ) _outputFile->Close();
  delete _outputFile;
  cout << " ....done." << endl;
}

void
SampleAnalysis::defineTemplate( const string & templ, 
				int nbin, float min, float max )
{
  cout << "1D-Template[" << templ << "]-->(nbin=" << nbin << "/min=" << min << "/max=" << max << ")" << endl; 
  hm.addTemplate<TH1F>( templ,
				   new TH1F( templ.c_str(), templ.c_str(), 
					     nbin, min, max ) );
}

void
SampleAnalysis::defineTemplate( const string & templ, 
				int nbin1, float min1, float max1,
				int nbin2, float min2, float max2 )
{
  cout << "2D-Template[" << templ << "]-->(nbin_x=" << nbin1 << "/min_x=" << min1 << "/max_x=" << max1 << ")" << "(nbin_y=" << nbin2 << "/min_y=" << min2 << "/max_y=" << max2 << ")" << endl; 
  hm.addTemplate<TH2F>( templ,
				   new TH2F( templ.c_str(), templ.c_str(), 
					     nbin1, min1, max1,
					     nbin2, min2, max2 ) );
}

TH1* 
SampleAnalysis::hist( const string & templ, const string & cut, const string & dir, const string & prefix )
{
  TH1* h_ = (TH1*)hm.h<TH1F>( templ, cut, dir, prefix );
  assert( h_!=0 );
  return h_;
}

void
SampleAnalysis::fill( const string & templ, const string & suffix, float val, const string & subdir )
{
//   string dir_ = _e.categ();
//   dir_ += "/";
//   dir_ += subdir;

  string dir_ = subdir; 
  dir_ += "/";
  dir_ += _e.categ();
  //  hist( templ, cut, _e.categ(), prefix )->Fill( val );  
  hist( templ, suffix, dir_, "" )->Fill( val );  

  if( EventServer::isMC )
    {
      //  cout<<dir_<<" 1  "<<endl;
      dir_ = "all";
      dir_ += "/";
      dir_ += subdir;
      hist( templ, suffix, dir_, "" )->Fill( val );  
    }

}

void
SampleAnalysis::fill( const string & templ, const string & suffix, float val1, float val2, const string & subdir )
{
//   string dir_ = _e.categ();
//   dir_ += "/";
//   dir_ += subdir;

  string dir_ = subdir; 
  dir_ += "/";
  dir_ += _e.categ();

  hist( templ, suffix, dir_, "" )->Fill( val1, val2 );  

  if( EventServer::isMC )
    {
      dir_ = "all";
      dir_ += "/";
      dir_ += subdir;
      hist( templ, suffix, dir_, "" )->Fill( val1, val2 );  
    }
}

SampleAnalysis::~SampleAnalysis()      
{
}

void SampleAnalysis::bookHistograms()  {}
bool SampleAnalysis::analyzeEvent()    { return false; }
void SampleAnalysis::writeHistograms() {}

bool
SampleAnalysis::CandCompare::operator()( Candidate* c1, Candidate* c2 ) const
{
  float pt1 = (c1!=0) ? c1->pt() : 0;
  float pt2 = (c2!=0) ? c2->pt() : 0;
  return pt1>pt2;
}

bool
SampleAnalysis::LeptonCompare::operator()( Candidate* c1, Candidate* c2 ) const
{
  CandInfo* info1 = c1->info();
  CandInfo* info2 = c2->info();
    
  bool isTight1 = info1->getBool("Tight");
  bool isTight2 = info2->getBool("Tight");
  if(  isTight1 && !isTight2 ) return true;
  if( !isTight1 &&  isTight2 ) return false;
    
  bool isMedium1 = info1->getBool("Medium");
  bool isMedium2 = info2->getBool("Medium");
  if(  isMedium1 && !isMedium2 ) return true;
  if( !isMedium1 &&  isMedium2 ) return false;
    
  bool isLoose1 = info1->getBool("Loose");
  bool isLoose2 = info2->getBool("Loose");
  if(  isLoose1 && !isLoose2 ) return true;
  if( !isLoose1 &&  isLoose2 ) return false;
    
  bool isVeryLoose1 = info1->getBool("VeryLoose");
  bool isVeryLoose2 = info2->getBool("VeryLoose");
  if(  isVeryLoose1 && !isVeryLoose2 ) return true;
  if( !isVeryLoose1 &&  isVeryLoose2 ) return false;
    
  float iso1(0.), iso2(0.);
  info1->getFloat("iso",iso1);
  info2->getFloat("iso",iso2);
  if( iso1!=iso2 ) return iso1>iso2;
    
  float pt1 = (c1!=0) ? c1->pt() : 0;
  float pt2 = (c2!=0) ? c2->pt() : 0;
  return pt1>pt2;
}

bool 
SampleAnalysis::ZllCompare::operator()( Candidate* Z1, Candidate* Z2 ) const
{
  LeptonCompare comp; // lepton comparator
    
  Candidate* Zcand_[2] = { Z1, Z2 };
  Candidate* dau_[2][2] = {{0,0},{0,0}};
  float m_[2]  = {0,0};
  float dm_[2] = {0,0};
  bool isInZWindow_[2]    = { false, false };
  bool isInVetoRegion_[2] = { false, false };
    
  for( int ii=0; ii<2; ii++ )
    {
      Candidate* dum1 = Zcand_[ii]->daughter(0);
      Candidate* dum2 = Zcand_[ii]->daughter(1);
      if( dum1->charge()>0 ) 
	{
	  dau_[ii][0] = dum1;
	  dau_[ii][1] = dum2;
	}
      else
	{
	  dau_[ii][0] = dum2;
	  dau_[ii][1] = dum1;
	}
      m_[ii]  = Zcand_[ii]->mass();

      // is the candidate in the Z mass region ?
      if( m_[ii]>Z0_window[0] && m_[ii]< Z0_window[1] ) 
	isInZWindow_[ii] = true;

      // is the candidate in the Z veto region ?
      if( m_[ii]>Z0_veto[0] && m_[ii]< Z0_veto[1] ) 
	isInVetoRegion_[ii] = true;

      // distance to the Z pole
      dm_[ii] = fabs( m_[ii]-Constants::Z0_m );
    }
    
  // by default take the candidate with the mass closest to the Z mass
  if( !isInVetoRegion_[0] &&  isInVetoRegion_[1] ) return false;
  if(  isInVetoRegion_[0] && !isInVetoRegion_[1] ) return true;
  if( !isInZWindow_[0]    &&  isInZWindow_[1]    ) return false;
  if(  isInZWindow_[0]    && !isInZWindow_[1]    ) return true;
    
  // check whether the two candidates have a common daughter
  int common_ = -1;
  for( int jj=0; jj<2; jj++ )
    {
      if( dau_[0][jj]->uid() == dau_[1][jj]->uid() )
	common_ = jj;
    }
    
  // if the two Z candidates have one lepton in common, 
  //  select the other lepton with best identification/isolation
  if( common_>=0 )
    {
      //	  cout << "Two candidates with a common dauthter" << endl;
      int other_ = 1-common_;
      // out_ = fabs(dau_[0][other_]->Pt())>fabs(dau_[1][other_]->Pt()); 
      return comp( dau_[0][other_], dau_[1][other_] );
    }
    
  // else select the Z candidate with the best lepton
  bool out_ = false;
  for( int jj=0; jj<2; jj++ )
    {
      if( comp( dau_[0][jj], dau_[1][0] )
	  && comp( dau_[0][jj], dau_[1][1] ) )
	{
	  out_ = true;
	  break;
	}
    }      
  return out_;
}

void
SampleAnalysis::buildLeptonLists()
{
  
  // prepare the standard lists of electrons
  buildElectronLists();

  // prepare the standard lists of muons
  buildMuonLists();

  // do the MC matching
  _e.mcMatching();
  {
    CandList matched   = _e.userCandList("ElectronVeryLoose_Matched");
    CandList unmatched = _e.userCandList("ElectronVeryLoose_Unmatched");
    const CandList& list = _e.userCandList("ElectronVeryLoose");
    for( size_t ii=0; ii< list.size(); ii++ )
      {
	Candidate* cand_ = list[ii];
	if( _e.mcCand( *list[ii] ) )   matched.push_back( cand_ );
	else                      unmatched.push_back( cand_ );
      }
    fillElectronHistograms( matched, "matched" );  
    fillElectronHistograms( unmatched, "unmatched" );  
  }
  {
    CandList matched   = _e.userCandList("MuonVeryLoose_Matched");
    CandList unmatched = _e.userCandList("MuonVeryLoose_Unmatched");
    const CandList& list = _e.userCandList("MuonVeryLoose");
    for( size_t ii=0; ii< list.size(); ii++ )
      {
	Candidate* cand_ = list[ii];
	if( _e.mcCand( *list[ii] ) )   matched.push_back( cand_ );
	else                      unmatched.push_back( cand_ );
      }
    fillMuonHistograms( matched, "matched" );  
    fillMuonHistograms( unmatched, "unmatched" );  
  }
}

void
SampleAnalysis::buildElectronLists()
{
  // prepare the lists of loose and tight electron and muon candidates
  ElectronSelector elSelVL( VBTFElectronSelector::kEt5, LeptonSelector::kVeryLoose );
  CandList& veryLooseElectrons = _e.userCandList("ElectronVeryLoose");
  elSelVL.getList( _e.electrons(), veryLooseElectrons );
  sort( veryLooseElectrons.begin(), veryLooseElectrons.end(), LeptonCompare() );

  ElectronSelector elSelL(  VBTFElectronSelector::kEt5, LeptonSelector::kLoose );
  CandList& looseElectrons = _e.userCandList("ElectronLoose");
  elSelL.getList( veryLooseElectrons, looseElectrons );

  ElectronSelector elSelM(  VBTFElectronSelector::kEt5, LeptonSelector::kMedium );
  CandList& mediumElectrons = _e.userCandList("ElectronMedium");
  elSelM.getList( looseElectrons, mediumElectrons );

  ElectronSelector elSelT( VBTFElectronSelector::kEt5, LeptonSelector::kTight );
  CandList& tightElectrons = _e.userCandList("ElectronTight");
  elSelT.getList( tightElectrons, tightElectrons );

  fillElectronHistograms( _e.electrons(),    "ref"  );
  fillElectronHistograms( veryLooseElectrons, "VL"  );
  fillElectronHistograms( looseElectrons,      "L"  );
  fillElectronHistograms( mediumElectrons,      "M"  );
  fillElectronHistograms( tightElectrons, "T"  );
}

void
SampleAnalysis::buildMuonLists()
{
  // prepare the lists of loose and tight muon and muon candidates
  MuonSelector selVL( MuonSelector::kPt5, LeptonSelector::kVeryLoose ); 
  CandList& veryLooseMuons = _e.userCandList("MuonVeryLoose");
  selVL.getList( _e.muons(), veryLooseMuons );
  sort( veryLooseMuons.begin(), veryLooseMuons.end(), LeptonCompare() );

  MuonSelector selL( MuonSelector::kPt5, LeptonSelector::kLoose );
  CandList& looseMuons = _e.userCandList("MuonLoose");
  selL.getList( veryLooseMuons, looseMuons );

  MuonSelector selM( MuonSelector::kPt5, LeptonSelector::kMedium );
  CandList& mediumMuons = _e.userCandList("MuonMedium");
  selM.getList( looseMuons, mediumMuons );

  MuonSelector selT( MuonSelector::kPt5, LeptonSelector::kTight );
  CandList& tightMuons = _e.userCandList("MuonTight");
  selT.getList( tightMuons, tightMuons );

  fillMuonHistograms( _e.muons(),    "ref"  );
  fillMuonHistograms( veryLooseMuons, "VL"  );
  fillMuonHistograms( looseMuons,      "L"  );
  fillMuonHistograms( mediumMuons,      "M"  );
  fillMuonHistograms( tightMuons, "T"  );
}

void
SampleAnalysis::fillElectronHistograms( const CandList& list, 
					const string & suffix    )
{
  for( size_t iel=0; iel<list.size(); iel++ )
    {
      Candidate* cand = list[iel];

      fillElectronHistograms( *cand, suffix, "el" );
    }
}

void
SampleAnalysis::fillElectronHistograms( const Candidate& cand, const string & suffix, const string & prefix )
{

  CandInfo* info_ = cand.info();

  float eta_ = /*cand.eta();*/ info_->getFloat("caloEta");
  float E_ = info_->getFloat("caloEnergy");
  float pt_  = /*cand.pt();*/  E_/cosh(eta_);
  float phi_ = cand.phi( kDeg );

  //     string prefix("el");

  fill( "pt",       suffix, pt_ , prefix );
  fill( "eta",      suffix, eta_ , prefix );
  fill( "phi",      suffix, phi_ , prefix );
  fill( "eta_phi",  suffix, phi_, eta_ , prefix );

  //      fill( "elID",     suffix, info_->getInt("elID"   ) , prefix );
  fill( "nTrk",     suffix, info_->getInt("nTrk"   ) , prefix );
  fill( "fBrem",    suffix, info_->getFloat("fBrem"   ) , prefix );
  fill( "EOP" ,     suffix, info_->getFloat("EOP"     ) , prefix );
  fill( "HOE" ,     suffix, info_->getFloat("HOE"     ) , prefix );
  fill( "dPhiIn" ,  suffix, info_->getFloat("dPhiIn"  ) , prefix );
  fill( "dEtaIn" ,  suffix, info_->getFloat("dEtaIn"  ) , prefix );
  fill( "sumOfPt2", suffix, info_->getFloat("sumOfPt2") , prefix );
  
  fill( "sumOfPt",  suffix, info_->getFloat("sumOfPt" ) , prefix );      
  float iso_(0.);
  if( info_->getFloat("iso",iso_) )
    fill( "iso",      suffix, iso_ , prefix );      
  fill( "dPhiIn_dEtaIn", 
	suffix, info_->getFloat("dEtaIn"), info_->getFloat("dPhiIn") , prefix );
  fill( "sigiEtaiEta", suffix, info_->getFloat("sigiEtaiEta"), prefix );
}

void
SampleAnalysis::fillMuonHistograms( const CandList& list, 
				    const string & suffix    )
{
  for( size_t imu=0; imu<list.size(); imu++ )
    {
      Candidate* cand = list[imu];

      fillMuonHistograms( *cand, suffix, "mu" );

    }
}

void
SampleAnalysis::fillMuonHistograms( const Candidate& cand, const string & suffix, const string & prefix )
{
  float pt_  = cand.pt();
  float eta_ = cand.eta();
  float phi_ = cand.phi( kDeg );

  float d0_ = cand.d0( &(_e.primaryVertex()) );

  CandInfo* info_ = cand.info();
  bool isGlobal = info_->getBool("muidGlobalMuonPromptTight");
  //  bool isGlobal = true;
  bool muIdTracker = info_->getBool("muidTrackerMuonArbitrated");
  //  bool muIdTracker = true;
  float normChi2 = info_->getFloat("normChi2");
  int nvHits = info_->getInt("nTrkHits");
  int nMatch = info_->getInt("nMatch");
  int nPixHit = info_->getInt("nPixHits");
  int nMuonHit = info_->getInt("nMuonHits");

  //  string prefix("mu");

  fill( "pt",       suffix, pt_, prefix );
  fill( "eta",      suffix, eta_, prefix );
  fill( "phi",      suffix, phi_, prefix );
  fill( "eta_phi",  suffix, phi_, eta_, prefix );
  float iso_(0.);
  if( info_->getFloat("iso",iso_) )
    fill( "iso",      suffix, iso_, prefix );      
  fill( "n20", string("nvHits_")+suffix, nvHits, prefix );
  fill( "n20", string("nMatch_")+suffix, nMatch, prefix );
  fill( "n20", string("nPixHits_")+suffix, nPixHit, prefix );
  fill( "n20", string("nMuonHit_")+suffix, nMuonHit, prefix );
  fill( "normChi2", suffix, normChi2, prefix );
  fill( "d0", suffix, d0_, prefix );
  fill( "bool",string("isGlobal_")+suffix, isGlobal, prefix );
  fill( "bool",string("muIdTracker_")+suffix, muIdTracker, prefix );
}

bool
SampleAnalysis::selectZCand( const Candidate& cand, float mmin, float mmax, 
			     unsigned nVT, unsigned nT, unsigned nL, unsigned nVL )
{      
  // kinematics of the candidate
  float mass_ = cand.mass();
      
  // mass cut
  bool massOK =  mass_>=mmin && mass_<=mmax;
  if( !massOK ) return false;

  assert( nVL>=nL && nL>=nT && nT>=nVT );
  if( nVL==0 ) return true;

  // require both leptons from same vertex
  int vtxid_ = cand.vertexIndex();
  if( vtxid_<0 ) return false;
      
  unsigned nVeryLoose(0), nLoose(0), nMedium(0), nTight(0);
  for( size_t idau=0; idau<cand.nDaughters(); idau++ )
    {
      const Candidate* dau = cand.daughter(idau);

      //*****
      //      dau->oneLine(cout);
      //*****

      CandInfo* info = dau->info();
      if( info->getBool("Tight") )
	{
	  nTight++;
	  nMedium++;
	  nLoose++;
	  nVeryLoose++;
	}
      else if( info->getBool("Medium") ) 
	{    
	  nMedium++;
	  nLoose++;
	  nVeryLoose++;
	}
      else if( info->getBool("Loose") )    
	{
	  nLoose++; 
	  nVeryLoose++;
	}
      else if( info->getBool("VeryLoose") )       
	{
	  nVeryLoose++;
	}      
    }
  if( nTight<nVT ) return false;
  if( nMedium<nT )      return false;
  if( nLoose<nL )      return false;
  if( nVeryLoose<nVL ) return false;

  return true;
}

void
SampleAnalysis::fillMultHists( const string & prefix, const string & suffix )
{
  EventServer& a_ = _e.a();

  vector<string> ntu_;
  ntu_.push_back("trk");
  ntu_.push_back("vtx");
  ntu_.push_back("ct");
  ntu_.push_back("el");
  ntu_.push_back("mu");
  ntu_.push_back("ph");
  ntu_.push_back("j");
  ntu_.push_back("sc");
  ntu_.push_back("rh");
  //  ntu_.push_back("pfj");
  for( size_t ii=0; ii<ntu_.size(); ii++ )
    {
      fill( "n", a_.getName( ntu_[ii] ) + suffix,   
	    a_.n( ntu_[ii] )   , prefix );
    }
}
  
