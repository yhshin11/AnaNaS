#include "Analysis/selectors/TagAndProbeElectronSelector.hh"

#include "Analysis/selectors/PtSelector.hh"
#include "Analysis/selectors/IsolationSelector.hh"
#include "Analysis/core/EventManager.hh"

#include "Analysis/utils/KineUtils.hh"

#include <cassert>

using namespace std;


ClassImp( TagAndProbeElectronSelector )
  
  
TagAndProbeElectronSelector::TagAndProbeElectronSelector(unsigned settings , CutMap cMap ) :
  VBTFSel( VBTFElectronSelector::kEt20 , VBTFElectronSelector::kPresel)
{
  _cMap = cMap;
  _settings = settings;


  _FiduEtaCuts.push_back(2.5);
  _FiduEtaCuts.push_back(1.442);
  _FiduEtaCuts.push_back(1.560);


  _debug =false;

}

// TagAndProbeElectronSelector::TagAndProbeElectronSelector(const CandList& cand, CutMap cMap ) :
//   VBTFSel( VBTFElectronSelector::kEt20 , VBTFElectronSelector::kPresel)
// {
  
// }

// TagAndProbeElectronSelector::TagAndProbeElectronSelector() :
//   VBTFSel( VBTFElectronSelector::kEt20 , VBTFElectronSelector::kPresel)
// {
  
// }


void
TagAndProbeElectronSelector::LoadVBTFSel(const Candidate& cand) {

  VBTFSel.computeFlags(cand);
  
}

bool 
TagAndProbeElectronSelector::accept( const Candidate& cand )
{
  
  if(abs(cand.pdgCode())==11)
    { VBTFSel.computeFlags(cand); }

  if(_debug) {
    cout<<endl<<" Begin accept of electron "<<endl;
  }

  string stmp;
  
  for(int unsigned it=0; it<_cMap.size(); it++) {
  
    bool Cpass_ = PassByType(cand, _cMap[it].first , _cMap[it].second  );
    if(!Cpass_)
      return false;
  }
    return true;
}

bool
TagAndProbeElectronSelector::PtEtaCut(const Candidate& cand, float ptCut, vector<float> etaCuts)
{

 if(_debug) {
   cout<<" Begin PtEtaCut "<<ptCut<<endl;
  }
 

 PtSelector PtSelec_(0, 1000 );
 PtSelec_.setPtmin( ptCut );
  
  if(etaCuts.size() == 1)
    PtSelec_.setEtamax( etaCuts[0] );
  else if(etaCuts.size() == 3) {
    PtSelec_.setEtamax( etaCuts[0] );
    PtSelec_.vetoRangeInEta( etaCuts[1] , etaCuts[2] );
  }
  else if(etaCuts.size() != 0)
    cerr<< " Warning on Eta cuts during Tag And Probe Analysis. \n Check etaCuts vector. Default cuts applied "<<endl;

   if(_debug) {
     cout<<" PtEta accepted ? "<<PtSelec_.accept( cand,  true)<<endl;
   }

   return PtSelec_.accept( cand,  true);
}

bool
TagAndProbeElectronSelector::HLTCut(const Candidate& cand, vector<string> Lines)
{
  EventManager& _e = *(EventManager::e());

  if(_debug) {
    cout<<" Begin HLT "<<Lines[0]<<endl;
  }

  // ****** /!\ FIXME MM /!\ ******
  bool fired =false;
  for(int unsigned i=0; i< Lines.size() ;i++) {
    
    if( _e.isFired(Lines[i]) )
      {  fired =true; }
  }

    if(_debug) {
      cout<<" HLT accepted ? "<<fired<<endl;
   }
  
  return fired;
}
 
bool
TagAndProbeElectronSelector::IsoCut(const Candidate& cand, vector<string> Isos) 
{

  if(_debug) {
    cout<<" Begin Iso "<<Isos[0]<<endl;
  }

  bool accept=true;
  for(int unsigned i=0; i<Isos.size(); i++) {
    bool iso_ = VBTFSel.isFlag(cand, Isos[i]);
    
if(!iso_)
      accept = false;
  }

  if(_debug) {
    cout<<" Iso accepted ? "<<accept<<endl;
   }

  return accept;

}

bool
TagAndProbeElectronSelector::IDCut(const Candidate& cand, vector<string> IDs)
{

  if(_debug) {
    cout<<" Begin ID "<<IDs[0]<<endl;
  }

  bool accept=true;
  for(int unsigned i=0; i<IDs.size(); i++) {
    bool eid_ = VBTFSel.isFlag(cand, IDs[i] );
   
    if(!eid_)
      accept = false;
  }

  if(_debug) {
    cout<<" ID accepted ? "<<accept<<endl;
   }

  return accept;
}

bool
TagAndProbeElectronSelector::GeneralCut(const Candidate& cand, std::string typevar,
				float val, std::string typvar, std::string cuttype )
{

  CandInfo* info =  cand.info();

  float var_;
  if(typvar=="f")
    var_ = (float)info->getFloat( typevar.c_str() );
  else if(typvar=="i")
    var_ = (float)info->getInt( typevar.c_str() );
  else {
    cerr<<" Bad type of Variable "<<endl;
    return false;
  }
  // cout<<" Variable générale "<<var_<<"   "<<val<<endl;
  if(cuttype=="g") 
    return QualityCutByGreaterVal(var_,val);
  else if(cuttype=="l") 
    return QualityCutByLowerVal(var_,val);
  else if(cuttype=="e") 
    return QualityCutByEgalVal(var_,val);
  else if(cuttype=="ne") 
    return QualityCutByNotEgalVal(var_,val);
  else
    return false;
  
  return true;
}

bool
TagAndProbeElectronSelector::MatchingDR(const Candidate& cand1, const Candidate& cand2, float cut, std::string cuttype )
{

  float dr = KineUtils::dR(cand1.eta(), cand2.eta(), cand1.phi(), cand2.phi() );

  //  cout<<dr<<"   " <<cuttype<<"    "<<cut<<endl;

   if(cuttype=="g") 
    return QualityCutByGreaterVal(dr,cut);
  else if(cuttype=="l") 
    return QualityCutByLowerVal(dr,cut);
  else if(cuttype=="e") 
    return QualityCutByEgalVal(dr,cut);
  else if(cuttype=="ne") 
    return QualityCutByNotEgalVal(dr,cut);
  else
    return false;
  

}


bool
TagAndProbeElectronSelector::PassByType(const Candidate& cand, string type,
					vector<string> vals, Candidate* addc ) 
{

    if(_debug) {
    cout<<" Begin Type "<<type<<endl;
  }
 
    bool pass_=true;

  //PtEta Cuts
  if(type=="PtEta") {
    
    float PtMin_ =0;
    vector<float> EtaCuts_(3,0);
    if(vals.size()==1) 
      PtMin_ = atof(vals[0].c_str());
    else if(vals.size()==2) {
      PtMin_ = atof(vals[0].c_str());
      if(vals[1]=="fiducial") {
	for(int ie=0;ie<3;ie++) {
	  EtaCuts_[ie] = _FiduEtaCuts[ie];
	}
      }
      else {
	EtaCuts_.push_back( atof(vals[1].c_str()) );
      }
    }
    else if(vals.size()==3) {
	for(int ie=0;ie<3;ie++) {
	  EtaCuts_[ie] = ( atof(vals[ie].c_str()) );
	}
    }
    else if(vals.size()==4) {
      PtMin_ = atof(vals[0].c_str());
      for(int ie=0;ie<3;ie++) {
	EtaCuts_[ie] = ( atof(vals[ie+1].c_str()) );
      }
    }
    else {
      cerr<<" Error in TagAndProbe PtEta Values , bad vector of values"<<endl;
      return false;
    }
    
    if( !PtEtaCut(cand, PtMin_, EtaCuts_ ) )
      {pass_ = false;}
  }

  //HLT
  else if(type=="HLT") {
    if( !HLTCut(cand, vals ) )
      { pass_ = false; }
  }
  
  //Leptons ID
  else if(type=="ID") {
   
    if( !IDCut(cand, vals ) )
      { pass_ = false; }
  }
  
  //Isolation
  else if(type=="Iso") {
    
    if( !IsoCut(cand, vals ) )
      { pass_ = false; }
  }

  //Matching Reco
  else if(type=="Reco") {
    float valcut_ = atof(vals[0].c_str());
    string cuttype_ = vals[1];
    if(  !MatchingDR(cand, (*addc), valcut_, cuttype_  ) )
      { pass_ = false; }
  }
  
  //Other
  else if(type=="Misc") {
    if(vals.size() ==4 ) {
      string typevar_ = vals[0];
      float valcut_ = atof(vals[1].c_str());
      string typv_ = vals[2];
      string cuttype_ = vals[3];
      if( !GeneralCut( cand, typevar_, valcut_, typv_, cuttype_) )
	{ pass_ = false; }
      
    }
    else {
      cout<<" En construction  --> "<<type<<"   "<<vals[0]<<endl;
    return false;
    }
  }
  else {
    cout<<" En construction  --> "<<type<<"   "<<vals[0]<<endl;
    return false;
  }
  if(_debug) {
    cout<<" Pass type ? "<<type<<"   "<<pass_<<endl;
  }

  
  return pass_;

}


