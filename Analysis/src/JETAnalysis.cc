#include "Analysis/src/JETAnalysis.hh"

#include <sstream>

using namespace std;



ClassImp( JETAnalysis )

  typedef vector<float> vectorFloat;


JETAnalysis::JETAnalysis( Sample& sample, EventManager& manager ):
			  // Sample& sample, const std::string& collectionFileName ):
SampleAnalysis( "JET", sample, manager ),
MuonSel(MuonSelector::kPt20,MuonSelector::kTight)
{

  cout << "\t-------------------------------------------" << endl; 
  cout << "\t---            JET  analysis            ----" << endl; 
  cout << "\t-------------------------------------------" << endl; 
  
  mainVtx= NULL;
  _ZCand= NULL;

  Nevt_ = 0;
  SampleName = sample.name();
  _nDebug = 0;

  SampleAnalysis::defineTemplate("etapt",200,0,500,50,-5,5);

}

JETAnalysis::~JETAnalysis() {
}


void
JETAnalysis::bookHistograms() {
}

void
JETAnalysis::writeHistograms() {
}

bool
JETAnalysis::analyzeEvent() {
  Nevt_++;
  if(Nevt_%100==0) cout<<"Nevt = "<<Nevt_<<endl;

  //Yeah!!!
  Candidate* ZCand = BuildBestZCand();

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


  //store useful variables

  FillJets();
  


  return true;
}


Candidate* JETAnalysis::BuildBestZCand() {

  const CandList muons = _e.muons();

  CandList ZCands;

  float pt1=0, pt2=0;
  float pt_1=0, pt_2=0;
  int nZ=-1;

  for(int unsigned il1=0;il1<muons.size();il1++) {

    Candidate* l1 = muons[il1];
 
    if(!MuonSel.accept( *l1 )) continue;
 
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
  
  if(nZ==-1)
    return NULL;

   return ZCands[nZ];
}


void JETAnalysis::FillJets() {


  const CandList& jets = _e.jetList( EventManager::kPatJet);

  for(unsigned int ijet=0;ijet<jets.size();ijet++) {
    
    const Candidate* jet = jets[ijet];
    const CandInfo* info = jet->info();

    if(isLepton(jet)) continue;

    bool matched=isJetMcMatched(jet);

    if(matched)
      fill( "etapt", "allHSJets", jet->pt(), jet->eta() );
    else
      fill( "etapt", "allPUJets", jet->pt(), jet->eta() );
    

    Vertex* vtxAJet = JetVtxAssoc(jet, "atlas");
    Vertex* vtxJet = JetVtxAssoc(jet, "basic");

    if(vtxJet == mainVtx) {
      if(matched)
	fill( "etapt", "JVCHSJets", jet->pt(), jet->eta() );
      else
	fill( "etapt", "JVCPUJets", jet->pt(), jet->eta() );
    }

    if(vtxAJet == mainVtx) {
      if(matched)
	fill( "etapt", "JVFHSJets", jet->pt(), jet->eta() );
      else
	fill( "etapt", "JVFPUJets", jet->pt(), jet->eta() );
    }

    if(info->getBool("mvaIdMedium") ) {
      if(matched)
	fill( "etapt", "mvaHSJets", jet->pt(), jet->eta() );
      else
	fill( "etapt", "mvaPUJets", jet->pt(), jet->eta() );
    }

  }

}




// Jet-vertex association

Vertex* JETAnalysis::JetVtxAssoc(const Candidate* jet, string opt) {

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

  return mVertex;
}



bool JETAnalysis::isLepton(const Candidate* jet) {

  if( jet->dR( _ZCand->daughter(0) )<0.1 && jet->pt()/_ZCand->daughter(0)->pt() >0.8 ) return true;
  if( jet->dR( _ZCand->daughter(1) )<0.1 && jet->pt()/_ZCand->daughter(1)->pt() >0.8 ) return true;

  return false;
}


bool JETAnalysis::isJetMcMatched(const Candidate* jet) {

  const CandList& mcJets = _e.jetList( EventManager::kGenJet);

  //cout<<" New jet "<<endl;

  bool matched=false;
  for(size_t ij=0;ij<mcJets.size();ij++) {

    const Candidate* mcjet = mcJets[ij];    

    if( jet->dR( mcjet ) < 0.1 && mcjet->pt()/jet->pt() > 0.3 )
      {matched=true;
	break;}
  }

  return matched;
  
}

