#include <iostream>
#include <errno.h>
#include <dirent.h>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include <string>
#include <map>
using namespace std;

#include <TString.h>
#include <TFile.h>
#include <RooDataSet.h>

#include "Analysis/utils/Constants.hh"
#include "Analysis/utils/Config.hh"
#include "Analysis/fit/Model.hh"
#include "Analysis/core/Sample.hh"

int main(int argc, char **argv)
{
  cout << "Hello, World" << endl;

  TString sampleFile = Config::confPath + "samples.txt";

  map< TString, Sample* >& samples = Sample::samples( sampleFile );

  //  Sample::SampleList QCDSampleList;
  //  QCDSampleList.push_back( samples["QCD_EM2030"] );
  //  QCDSampleList.push_back( samples["QCD_EM80170"] );
  //  QCDSampleList.push_back( samples["QCD_EM3080"] );
  //  samples["QCD"] = new Sample("QCD", QCDSampleList );
  //  samples["QCD"]->print(cout);

  //  Sample::SampleList GamJetSampleList;
  //  GamJetSampleList.push_back( samples["GamJet1520"] );
  //  GamJetSampleList.push_back( samples["GamJet2025"] );
  //  GamJetSampleList.push_back( samples["GamJet2530"] );
  //  GamJetSampleList.push_back( samples["GamJet3035"] );
  //  GamJetSampleList.push_back( samples["GamJet35"] );
  //  samples["GamJet"] = new Sample("GamJet", GamJetSampleList );
  //  samples["GamJet"]->print(cout);

  // FIXME !!!! temporary !!!
  //  Sample::SampleList ZZSampleList;
  //  ZZSampleList.push_back( samples["ZZ_2l2n"] );
  //  ZZSampleList.push_back( samples["ZZ_4l"]   );
  //  samples["ZZ"] = new Sample("ZZ", ZZSampleList );
  //  samples["ZZ"]->print(cout);

  Sample::SampleList sampleList;
  //  sampleList.push_back( samples["QCD"]     );
  //  sampleList.push_back( samples["GamJet"]  );
  sampleList.push_back( samples["WW_2l2n"] );
  sampleList.push_back( samples["WZ_3ln"]  );
  sampleList.push_back( samples["ZZ_2l2n"] );
  sampleList.push_back( samples["ZZ_4l"]   );
  //  sampleList.push_back( samples["WJet"]    );
  sampleList.push_back( samples["W_en"]    );
  sampleList.push_back( samples["ttbar"]   );
  sampleList.push_back( samples["Z_2t"] );
  //  sampleList.push_back( samples["Zmm"] );
  sampleList.push_back( samples["Z_2e"] );
  samples["total"] = new Sample("total", sampleList );
  samples["total"]->print(cout);  
    
  float lumi=100;  // in pb-1 
  TString analysisName  = "Diboson";
  TString variableName  = "mll";
  TString histName      = "stat__mll__0_all_nosel";
  // TString histName = "mll_k11";

  map< Sample*, TH1* > mapOfHists;
  vector< TH1* >    listOfStackedHists;
  samples["total"]->getListOfStackedHists( listOfStackedHists, mapOfHists, 
					   lumi, 
					   analysisName, "all", histName );

  TString outname = Config::histPath + "stack_" 
    + analysisName + "_" + "_" + histName + "_";
  outname += int(lumi);
  outname += ".root";

  TFile* f_out = TFile::Open(outname,"RECREATE");

  cout << "\nDump map of hists " << endl;
  for( map< Sample*, TH1* >::iterator it_=mapOfHists.begin(); 
       it_!=mapOfHists.end(); ++it_ )
    {
      Sample* s_ = it_->first;
      TH1* h_ = it_->second;
      cout << ">>>>>>> " << s_->name() << endl;
      assert( h_!=0 ); 
      h_->Print();
      h_->Write();
    }
  cout << "\nDump list of stacked hists" << endl;
  for( size_t ii=0; ii<listOfStackedHists.size(); ii++ )
    {
      TH1* h_ = listOfStackedHists[ii];
      assert( h_!=0 );
      h_->Print();
      h_->Write();
    }
  
//    for( map< TString, Sample* >::iterator it_=samples.begin();
//         it_!=samples.end(); ++it_ )
//      {
//        Sample* s_ = it_->second;
      
//        float n_ = s_->n( lumi );
//        cout << "n[" << s_->name() << ",L=" << lumi << "]=" << n_ << endl;
//        float w_ = s_->w( lumi );
      
//        // FIXME !!!
//        //       Model* model_ = s_->getModel( analysisName, variableName, histName ); 
//        Model* model_ = new Model( analysisName, s_->name(), s_->subsample(), histName, variableName ); 
//        if( model_==0 ) continue;

//        RooDataSet* set_ = model_->generate( w_ );
      
//        TString setName_ = "set_";
//        setName_ += s_->name();
//        setName_ += "_";
//        setName_ += int(lumi);
//        set_->SetName( setName_ );
      
//        set_->Write();      
//      }
  
  f_out->Close();


  cout << "Good bye, World" << endl;
  
  return 0;
}
