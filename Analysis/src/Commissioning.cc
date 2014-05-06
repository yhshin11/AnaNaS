#include "Analysis/src/Commissioning.hh"

#include <cassert>
#include <sstream>
#include <algorithm>
using namespace std;

#include <TH1I.h>
#include <TH1F.h>
#include <TH2I.h>
#include <TH2F.h>
#include <TMath.h>

#include "Analysis/utils/Config.hh"
#include "Analysis/utils/Constants.hh"
#include "Analysis/core/Candidate.hh"
#include "Analysis/core/CandInfo.hh"
#include "Analysis/core/EventManager.hh"
#include "Analysis/tools/CandUtil.hh"
#include "Analysis/selectors/ElectronSelector.hh"
#include "Analysis/selectors/MuonSelector.hh"
#include "Analysis/selectors/LeptonSelector.hh"
#include "Analysis/selectors/ConeSelector.hh"


//kinematical cuts
static const float phot_et_cut = 0.300;
static const float diphot_et_cut = 0.900;

static double deltar(const Candidate& a, const Candidate& b){
  double dphi = fabs(a.phi(kRad) - b.phi(kRad));
  dphi = min(dphi, 2*TMath::Pi()-dphi);
  double deta = fabs(a.eta() - b.eta());
  return sqrt(dphi*dphi+deta*deta);
}

ClassImp( Commissioning )

Commissioning::Commissioning( Sample& sample, const std::string & collectionFileName ) : 
  SampleAnalysis( "Commissioning", sample, collectionFileName )
{
  cout << "\t-------------------------------------------------" << endl; 
  cout << "\t---  Commissioning preselection analysis     ----" << endl; 
  cout << "\t-------------------------------------------------" << endl; 

  _hlt = false; // do not apply HLT selection

  //_hltLines.push_back("HLT_DoubleEle10_SW_L1R");
  //_hltLines.push_back("HLT_DoubleMu3");
  //_hltLines.push_back("HLT_Mu15");
  //_hltLines.push_back("HLT_Ele15_SW_L1R");   // 20 or 25

  _nDebug = 10;
}

Commissioning::~Commissioning()
{
}

void
Commissioning::bookHistograms()
{
  defineTemplate(    "Cmgg", 100,   0., 0.5    );
  defineTemplate(   "Cnvtx", 100,   0,  100    );
  defineTemplate(     "Cvx", 400,  -2.,   2.   );
  defineTemplate(     "Cvy", 400,  -2.,   2.   );
  defineTemplate(     "Cvz", 800, -20.,  20.   );
  defineTemplate(     "Cpt", 500,   0.,  50.   );
  defineTemplate(    "Cmet", 500,   0., 100.   );
  defineTemplate(    "Ceta", 240,  -3.,   3.   );
  defineTemplate(    "Cphi", 126,   0.,  6.3   );

  defineTemplate("CEt", 100, 0, 300);
  defineTemplate("CDr", 100, -4, 4);

  defineTemplate("CR4", 100, 0, 1);
}

void
Commissioning::writeHistograms()
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

  // print out runs
  for( map<int,int>::const_iterator it=_runs.begin(); it!=_runs.end(); it++ )
    {
      cout << "n[" << it->first << "]=" << it->second << endl;
    }

}

bool
Commissioning::analyzeEvent()
{
  int run_ = _e.run();
  ostringstream oss;
  oss << run_;
  string str_ = oss.str();

  string suffix("Com");


  static int i = 0;

  ++i;

  if(i%1000==1) cout << "--> " << i << "\n";
  
  
  // int evt_ = _e.event();

    // get the candidate list
    CandList listPhot = _e.candList( EventManager::kEcalSC );

    if ( listPhot.size() >= 2 )
      {
	for ( size_t i = 0; i < listPhot.size()-1; i++ ) {
	  Candidate* dau1 = listPhot[i];
	  CandInfo* info1 = dau1->info();

	  float r4_1 = info1->getFloat("R4");
	  
	  fill("CEt", "PhotEt", dau1->Et(), suffix);
	  fill("CEt", "PhotE",  dau1->E(), suffix);	  
	  fill("eta", "PhotEta", dau1->eta());
	  fill("Cphi", "PhotPhi", dau1->phi());
	  fill("CR4", "R4_1", r4_1);
	  
	  if(dau1->Et() < phot_et_cut) continue;

	  if(r4_1<0.85) continue;
	  
	  for ( size_t j = i+1; j < listPhot.size(); j++ ) {	    
            Candidate* dau2 = listPhot[j];
	    CandInfo* info2 = dau2->info();
	    
	    float r4_2 = info2->getFloat("R4");


	    if(dau2->Et() < phot_et_cut) continue;	    

            Candidate* c = Candidate::create( dau1, dau2 );

	    fill("CEt", "DiPhotEt", diphot_et_cut, suffix);
	    fill("eta", "DiPhotEta", c->eta());
	    fill("Cphi", "DiPhotPhi", c->phi());

	    fill("CDr", "drPhots", deltar(*dau1, *dau2)); 
	    fill("CR4", "R4_2", r4_2);
	    
	    if(c->Et() < diphot_et_cut) continue;

	    if(r4_2<0.85) continue;
	    
            //c->print( std::cout );
            fill("Cmgg", str_, c->mass(), "EcalSC" );
            fill("Cmgg", "all", c->mass(), "EcalSC" );
	  }
	}
      }

  const Vertex& vtx_ = _e.primaryVertex();
  const CandInfo* info = vtx_.info();
  int vtx_nTracks(0);

  if( !info->getInt( "nTracks", vtx_nTracks ) ) return false;
  if( vtx_nTracks==0 ) return false;
  //  cout << "n[tracks] at vertex=" << vtx_nTracks << endl;

  const TVector3& pos_ = vtx_.pos();

  fill( "Cnvtx", str_,  vtx_nTracks , suffix );
  fill(   "Cvx", str_,  pos_.X()    , suffix );
  fill(   "Cvy", str_,  pos_.Y()    , suffix );
  fill(   "Cvz", str_,  pos_.Z()    , suffix );
  fill( "Cnvtx", "all",  vtx_nTracks , suffix );
  fill(   "Cvx", "all",  pos_.X()    , suffix );
  fill(   "Cvy", "all",  pos_.Y()    , suffix );
  fill(   "Cvz", "all",  pos_.Z()    , suffix );

  size_t nTracks    = _e.tracks().size();
  size_t nElectrons = _e.electrons().size();
  size_t nMuons     = _e.muons().size();
  size_t nPhotons   = _e.photons().size();
  if(1) 
    {
      float TrackPtMin = 0.7;
      if( nTracks==0 ) return false;  
      //      cout << "n[tracks]=" << nTracks << endl;
      if( nTracks>1000 ) return false;  
      //      cout << "Tracks above pt=" << TrackPtMin << "GeV" << endl;
      for( size_t ii=0; ii<nTracks; ii++ )
	{
	  const Candidate* c_ = _e.tracks()[ii];
	  if( c_->pt()<TrackPtMin ) continue;
	  fill( "pt",  str_, c_->pt() , suffix );
	  fill( "eta", str_, c_->eta() , suffix );
	  fill( "phi", str_, c_->phi( kDeg ) , suffix );
	  fill( "eta_phi", str_, c_->phi( kDeg ), c_->eta() , suffix );
	  fill( "pt",  "all", c_->pt() , suffix );
	  fill( "eta", "all", c_->eta() , suffix );
	  fill( "phi", "all", c_->phi( kDeg ) , suffix );
	  fill( "eta_phi", "all", c_->phi( kDeg ), c_->eta() , suffix );
	  //	  c_->oneLine( cout );
	}
    }
  
  if(0)
    {
      if( nElectrons==0 ) return false;
    }
  
  if(0)
    {
      if( nMuons==0 ) return false;
    }
  
  if(0)
    {
      if( nPhotons==0 ) return false;
    }
  
  if(0)
    {
      if( nElectrons==0 && nMuons==0 && nPhotons==0 ) return false;
    }
  
  if( _runs.count(run_)==0 ) _runs.insert( make_pair( run_, 0 ) ); 
  int n_ = _runs[run_];
  _runs[run_] = ++n_;
    
    // combine them to see a mass
    return true;
}
