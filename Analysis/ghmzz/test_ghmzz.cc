#include <cassert>
#include <string>
#include <typeinfo>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <cassert>
#include <vector>
#include <map>
#include <algorithm>

using namespace std;

#include <TVector2.h>
#include <TVector3.h>
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

#include "Analysis/utils/Config.hh"
#include "Analysis/utils/DatasetManager.hh"
#include "Analysis/utils/TupleManager.hh"
#include "Analysis/utils/Constants.hh"
#include "Analysis/utils/KineUtils.hh"


using namespace std;

typedef vector<double> vectorFloat; 
typedef pair<string,float> Line; 
struct SortedLine
{
  bool operator()( const Line& p1, const Line& p2 ) const 
  {
    double dr1 = p1.second;
    double dr2 = p2.second;	      
    return dr1 < dr2;
  } 
};

bool tex = false;
bool firstLine = true;

TString infname = Config::histPath + "ZZ_2l2n.root";
TString oufname = Config::histPath + "ZZ_2l2n_ds.root";
string tupname = "ZZtuple";
string dirname = "data/q30_jet_met_bal_phj_phz_btg_dib";

void 
eventPrintouts()
{

  TFile* inputFile  = new TFile( infname, "READ" );

  TFile* outputFile  = new TFile( oufname, "RECREATE" );
  DatasetManager dsm("dataset");
  dsm.setFile( outputFile, DatasetManager::kWrite );

  TupleManager tm;
  tm.setFile( inputFile );
  TTree* tree = tm.getTree( tupname, dirname );

  TBits   c_bits;
  TBits*  c_bits_ptr = &c_bits;
  string  c_str;
  string* c_str_ptr  = &c_str;
  string  categ;
  string* categ_ptr = & categ;

  // declare the main analysis variables
  int   nZCand(0);        // number of Z candidates
  int   nVertex(0);       // number of vertices
  float mll(0);           // di-lepton mass
  float qTll(0);          // di-lepton transverse momentum
  float phill(0);         // di-lepton phi
  float yll(0);           // di-lepton rapidity
  float hll(0);           // di-lepton helicity
  float MET(0);           // PF-MET type-I
  float METPhi(0);        // PF-MET type-I - phi
  float sumEt(0);         // sum of transverse energies
  float mEtSig(0);        // sigma of MET
  float mTZZ(0);
  float projMET(0);       // projected MET
  float corProjMET(0);    // PU-vertex corrected MET
  float sigMET(0);        // corrected MET significance
  float ghmMET(0);
  float trueMET(0);
  float unclMET(0);
  float ghmSig(0);
  float dPhiMin(0);       // minimum angle of MET to a lepton
  float balZMET;          // "balance" Z/MET
  bool  isTight(false);   // is the Z candidate tight?
  int   nLoose(0);        // number of additional loose leptons
  int   jMult25(0);       // number of jets above 25 GeV-Et
  int   jMult20(0);       // number of jets above 25 GeV-Et
  int   jMult15(0);       // number of jets above 15 GeV-Et
  int   jMult10(0);       // number of jets above 10 GeV-Et
  float JZB15(0);
  float JZB20(0);
  float JZB25(0);
  float dPhiJMET(0);      // angle between leading jet and MET
  
  float jEt_1(0);
  float jEt_2(0);
  float jEt_3(0);
  float jEt_4(0);
  float dPhiSumJMET(0);
  float dPhiZMET(0);
  float thirdL_MT(0);
  float jBtagTkC_max(0);
  float jBtagSoftM_max(0);
  
  // The JZB variable (see CMS-PAS-SUS-10-010)
  //    JZB = |−∑pTjets| − |pTZ|
  
  // create the jet vectors to be stored in the ntuple -- MUST be destroyed!
  // (warning: no return before the end)
  // Note to myself: implement garbage collection for these...
  vectorFloat* jEt      = new vectorFloat;
  vectorFloat* jEta     = new vectorFloat;
  vectorFloat* jPhi     = new vectorFloat;
  vectorFloat* jDR1     = new vectorFloat;
  vectorFloat* jDR2     = new vectorFloat;
  vectorFloat* jDPhiMET = new vectorFloat;
  vectorFloat* jBtagTkC   = new vectorFloat;
  vectorFloat* jBtagSoftM = new vectorFloat;
  
  vectorFloat* LPt  = new vectorFloat;
  vector<int>* LPdg = new vector<int>;
  vectorFloat* LEta = new vectorFloat;
  vectorFloat* LPhi = new vectorFloat;
  vectorFloat* LMT   = new vectorFloat;
  vector<int>* LID  = new vector<int>;
  vectorFloat* LTrkIso  = new vectorFloat;
  vectorFloat* LEcalIso = new vectorFloat;
  vectorFloat* LHcalIso = new vectorFloat;

  int run;
  int event;
  tm.add< int >(    "run",   &run       );
  tm.add< int >(    "event", &event     );
  tm.add< string*>( "categ", &categ_ptr );
  tm.add< string*>( "c_str", &c_str_ptr  );      
  tm.add< TBits* >( "c_bits", &c_bits_ptr  );      
  tm.add< int >(  "nZCand",   &nZCand        );
  tm.add< int >(  "nVertex",   &nVertex      );
  tm.add< bool  >(  "isTight",&isTight    );
  tm.add< float >(  "mll",   &mll        );
  tm.add< float >(  "qTll",   &qTll        );
  tm.add< float >(  "phill",  &phill        );
  tm.add< float >(  "yll",   &yll        );
  tm.add< float >(  "hll",   &hll        );
  tm.add< int >(  "nLoose",&nLoose    );
  tm.add< vectorFloat* >(  "LPt",   &LPt   );
  tm.add< vector<int>* >(  "LPdg",   &LPdg   );
  tm.add< vectorFloat* >(  "LEta",  &LEta  );
  tm.add< vectorFloat* >(  "LPhi",  &LPhi  );
  tm.add< vector<int>* >(  "LID",   &LID   );
  tm.add< vectorFloat* >(  "LMT",   &LMT   );
  tm.add< vectorFloat* >( "LTrkIso", &LTrkIso );
  tm.add< vectorFloat* >( "LEcalIso", &LEcalIso );
  tm.add< vectorFloat* >( "LHcalIso", &LHcalIso );
  tm.add< float >(  "MET",   &MET        );
  tm.add< float >(  "METPhi",&METPhi    );
  tm.add< float >(  "sumEt", &sumEt );
  tm.add< float >(  "mEtSig", &mEtSig );
  tm.add< float >(  "mTZZ", &mTZZ );
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
  tm.add< int >(  "jMult25",&jMult25  );
  tm.add< int >(  "jMult20",&jMult20  );
  tm.add< int >(  "jMult15",&jMult15  );
  tm.add< int >(  "jMult10",&jMult10  );
  tm.add< float >(  "JZB25",&JZB25  );
  tm.add< float >(  "JZB20",&JZB20  );
  tm.add< float >(  "JZB15",&JZB15  );
  tm.add< float >(  "dPhiJMET",&dPhiJMET );	  

 
  vector<Line> lines; 
  int n = tree->GetEntriesFast();
  int n_win=0;

  for( int ii=0; ii<n; ii++ )
    {
      tree->GetEntry( ii );

      if( mll>=80 && mll<=100 ) n_win++;

      if( firstLine ) continue;
      
      TVector2 MET2D( MET*cos( METPhi), MET*sin( METPhi ) );
      TVector2 ll2D( qTll*cos(phill), qTll*sin(phill) );	     
      TVector2 minusMET2D = MET2D;
      minusMET2D *= -1;
      float dPhiZMET    = KineUtils::dPhi( ll2D.Phi(), minusMET2D.Phi() )
	* Constants::radToDeg;	      

      ostringstream o;
      o << fixed << setw(7) << run;
      if( tex ) o << fixed << setw(2) << "&";
      o << fixed << setw(12) << event;
      if( tex ) o << fixed << setw(2) << "&";
      o << fixed << setw(3) << nVertex; 
      if( tex ) o << fixed << setw(2) << "&";
      if( abs((*LPdg)[0])==11 ) 
	{
	  if( tex )
	    o << setw(2) << "$\\Pe$";
	  else
	    o << setw(2) << "e";
	}
      else
	{
	  if( tex )
	    o << setw(2) << "$\\Pgm$";     
	  else
	    o << setw(2) << "m";     
	}
      if( tex ) o << fixed << setw(2) << "&";
      o.precision(1);
      o << fixed << setw(7) << mll;
      if( tex ) o << fixed << setw(2) << "&";
      o << fixed << setw(7) << qTll;
      if( tex ) o << fixed << setw(2) << "&";
      o << fixed << setw(7) << MET;       
      if( tex ) o << fixed << setw(2) << "&";
      o << fixed << setw(5) << balZMET;       
      if( tex ) o << fixed << setw(2) << "&";
      o.precision(0);
      if( dPhiJMET<-180 ) 
	o << fixed << setw(7) << "--";       
      else
	o << fixed << setw(7) << dPhiJMET;       
      if( tex ) o << fixed << setw(2) << "&";
      o << fixed << setw(7) << dPhiZMET;       
      if( tex ) o << fixed << setw(2) << "&";
      o.precision(1);
      //      o << fixed << setw(7) << sumEt;       
      //      if( tex ) o << fixed << setw(2) << "&";
      o << fixed << setw(7) << mTZZ;       
      if( tex ) o << fixed << setw(2) << "&";
      o << fixed << setw(7) << (*LPt)[0];       
      if( tex ) o << fixed << setw(2) << "&";
      o << fixed << setw(7) << (*LPt)[1];       
      if( tex ) o << fixed << setw(3) << "\\\\";
      o << endl;

      Line line;
      line.first = o.str();
      line.second = mll;
      lines.push_back( line );

      //      lines.push_back(pair(o.str(), mll));
    }
  sort( lines.begin(), lines.end(), SortedLine() );
  //  lines.sort();

  if(firstLine )
    {
      firstLine = false;
      if( tex )
	{
	  cout << "\\begin{table}" << endl;
	  cout << "\\begin{center}" << endl;
	  cout << "\\small" << endl;
	  cout << "\\caption{  List of the " 
	       << n 
	       << " events in the final sample (of which " 
	       << n_win 
	       << " are in the reduced signal window $[\\,80$-$100\\,]~\\GeV$). All masses, energies, momenta are expressed in \\GeV; all angles, in degrees. The \\pt\\ of the leptons is multiply by the lepton electric charge.}" << endl;
	  cout << "\\begin{tabular} {|c|c|c|c|c|c|c|c|c|c|c|c|c|}" << endl;
	  cout << "\\hline" << endl;
	}
      else 
	{
	  cout << "Number of events in sample: " << n << endl;
	  
	}
      ostringstream o;      
      if( tex )
	{
	  o << fixed << setw(7) << "run";
	  o << fixed << setw(3) << "&";
	  o << fixed << setw(12)<< "event";
	  o << fixed << setw(3) << "&";
	  o << fixed << setw(3) << "nV"; 
	  o << fixed << setw(3) << "&";
	  o << fixed << setw(3) << "$\\ell$";
	  o << fixed << setw(3) << "&";
	  o << fixed << setw(7) << "$m(\\ell\\ell)$";
	  o << fixed << setw(3) << "&";
	  o << fixed << setw(7) << "$\\qt(\\ell\\ell)$";
	  o << fixed << setw(3) << "&";
	  o << fixed << setw(7) << "$\\met$";       
	  o << fixed << setw(3) << "&";
	  o << fixed << setw(5) << "$\\BAL$";       
	  o << fixed << setw(3) << "&";
	  o << fixed << setw(7) << "$\\Delta\\Phi({\\mathrm J})$";       
	  o << fixed << setw(3) << "&";
	  o << fixed << setw(7) << "$\\Delta\\Phi(\\PZ)$";       
	  o << fixed << setw(3) << "&";
	  //	  o << fixed << setw(7) << "$\\Sigma\\et$";       
	  //	  o << fixed << setw(3) << "&";
	  o << fixed << setw(7) << "$m_{\\mathrm T}(\\PZ\\PZ)$";       
	  o << fixed << setw(3) << "&";
	  o << fixed << setw(7) << "$\\pt(\\ell_{1})$";
	  o << fixed << setw(3) << "&";
	  o << fixed << setw(7) << "$\\pt(\\ell_{2})$";
	  o << fixed << setw(3) << "\\\\";
	  o << endl;
	  o << "\\hline" << endl;
	  o << "\\hline" << endl;
	}
      else
	{

	  o << fixed << setw(7) << "run";
	  o << fixed << setw(12)<< "event";
	  o << fixed << setw(3) << "nV"; 
	  o << fixed << setw(2) << "l";
	  o << fixed << setw(7) << "mll";
	  o << fixed << setw(7) << "qTll";
	  o << fixed << setw(7) << "MET";       
	  o << fixed << setw(5) << "Bal";       
	  o << fixed << setw(7) << "dPhJ";       
	  o << fixed << setw(7) << "dPhZ";       
	  //	  o << fixed << setw(7) << "SigET";       
	  o << fixed << setw(7) << "mTZZ";       
	  o << fixed << setw(7) << "pTl1";
	  o << fixed << setw(7) << "pTl2";
	  o << endl;
	}
      cout << o.str();
      eventPrintouts();
      return;
    }

  bool b_ = false;
  for( vector<Line>::iterator it = lines.begin(); it!=lines.end(); it++ )
    {
      Line line = (*it);
      float mll_ = line.second;

      dsm.add( "mll", mll_, 60, 120 );
      dsm.flush( dirname );

      if( tex ) 
	{
	  if( mll_<=100 )
	    {
	      if( !b_ && mll_>=80 ) 
		{ 
		  b_ = true;
		  cout << "\\hline" << endl;
		}
	    }
	  else if( b_ ) 
	    { 
	      b_ = false;
	      cout << "\\hline" << endl;
	    }
	}
      cout << line.first;
    }

  dsm.save();
  outputFile->Close();

  if( tex )
    {
      cout << "\\hline" << endl;
      cout << "\\end{tabular}" << endl;
      cout << "\\end{center}" << endl;
      cout << "\\end{table}" << endl;
    }

}



int main()
{
  bool write=true;
  if( write )
    {
      eventPrintouts();
    }
  else
    {
      cout << "Test reading" << endl;  
      TFile* inputFile  = new TFile( oufname, "READ" );
      inputFile->Print();
      
      DatasetManager dsm("dataset");
      
      cout << "set file " << endl;
      dsm.setFile( inputFile, DatasetManager::kRead );
      
      RooDataSet* ds = dsm.dataset( dirname );
      ds->Print();
      int numEntries = ds->numEntries();
      for( int jj=0; jj<numEntries; jj++ )
	{
	  cout << "jj=" << jj 
	     << " Mll=" << ds->get(jj)->getRealValue("mll")
	       << endl;
	}
    }
  return 0;
}
