// Implemantation of the Squeezer main class
// Author: Ph. Gras CEA/IRFU Saclay
// Jan. 3, 11

#include <iostream>
#include <algorithm>
#include <iomanip>
#include <sys/time.h>
			 
#include "TObjArray.h"
#include "TBranch.h"
#include "TFile.h"
#include "TKey.h"
#include "TIterator.h"

#include "CopyAnaTuple.h"

using std::cout;
using std::cerr;
using std::endl;


void CopyAnaTuple::listEvents(std::ostream& o, const char* inputDataFile){
    
  fin_ = std::auto_ptr<TFile>(new TFile(inputDataFile));
    
  if(fin_->IsZombie()){ cout << "Failed to open input ntuple file "
                             << inputDataFile << endl;
    return;
  }
  setEventSummaryTree();
  for(Long64_t i = 0; i < eventSummary_->GetEntries(); ++i){
    eventSummary_->GetEntry(i);
    o << std::setw(8) << runNum_ << std::setw(8) << " " << eventNum_ << endl;
  }
}

void CopyAnaTuple::listBranches(std::ostream& o, const char* inputDataFile){
    
  init(inputDataFile);

  TIter itTree(&trees_);
  while(itTree.Next()){
    TTree* tree = ((TreeRcd*) (*itTree))->tree;
    //TO CHECK: is list of leaves fine or should we use GetListOfBranches
    //recursively?
    TObjArray* ls = tree->GetListOfLeaves(); 
    TIter it(ls);
    while(it.Next()) o << tree->GetName() << "\t" << (*it)->GetName() << endl;
  }
}

void CopyAnaTuple::fillRunSummary(){
  TTree* tree = dynamic_cast<TTree*>(fin_->FindObject("runSummary"));

  if(tree==0){
    cout << "Error: run summary tree was not found!\n";
    return;
  }

  Int_t selected;
  tree->SetBranchAddress("selected", &selected);
    
  fout_->cd();
  TTree* treeOut = tree->CloneTree(0);
  tree->GetEntry(0);
  if(selected != eventSummary_->GetEntries()){
    cout << "Error: selected variable contains a different run "
      "summary than event summary tree entry count!\n";
  }
  selected = nCopied_;
  treeOut->Fill();
  treeOut->AutoSave();
  fin_->cd();
}
  
void CopyAnaTuple::run(const char* inputDataFile, const char* outputDataFile){
  timeval start;
    
  gettimeofday(&start, 0);

  init(inputDataFile, outputDataFile);
        
  fin_->cd();

  ////----    
  if(verbose_){
    cout << "Number of trees to copy: " << trees_.GetEntries() << endl;
    cout << "List of trees to copy: " << endl;
    TIter itT(&trees_);
    while(itT.Next()){
      cout << "\t" << ((TreeRcd*)(*itT))->tree->GetName() << endl;
    }
  }

  Long64_t nevts = eventSummary_->GetEntries();
  if(maxEvents_ >=0 && nevts > maxEvents_) nevts = maxEvents_;

  timeval t;
  for(Long64_t i = 1; i <= nevts; ++i){
    nextEvent();
    if(filterEvent()){
      copyEvent();
      ++nCopied_;
    }
    const static int step = 100;
    if(i%step==0 || i == nevts){
      cout << "\rNumber of events. Read: " << std::setw(8) << i << " Copied "
           << std::setw(8) << nCopied_
           << " Remaining: " << std::setw(8) << (nevts-i)
           << " Total: " << std::setw(8) << nevts; // << std::flush;
      timeval t0;
      if(i==step) gettimeofday(&t0, 0);
      else if(i>=2 * step) {
        gettimeofday(&t, 0);
        double remaining = double(nevts - step) / (i - step)
          * ((t.tv_sec - t0.tv_sec) + 1.e-6 * (t.tv_usec - t0.tv_usec));
        time_t eat = int(t0.tv_sec +  1.e-6 * t0.tv_usec + remaining + 0.5);
        char buf[256];
        strftime(buf, sizeof(buf),  "%a, %d %b %Y %T", localtime(&eat));
        cout << " EAT: " << std::setw(16) << buf;
      }
      cout << std::flush;
    }
  } //next event
  cout << "\n";
    
  fillRunSummary();
    
  TIter itTree(&trees_);
  while(itTree.Next()){ ((TreeRcd*) (*itTree))->outTree->AutoSave(); }
  fout_->Close();

  gettimeofday(&t, 0);

  double dt = (t.tv_sec - start.tv_sec) + 1.e-6 * (t.tv_usec - start.tv_usec);

  printf("Elapsed time: %2f sec. %2f event/sec.\n",  dt, nevts/dt);
  fflush(stdout);
}

void CopyAnaTuple::readEventList(const char* fileName){
  allEvent_ = false;
  FILE* f = fopen(fileName, "r");
  if(f==0) { cout << "Failed to open event list file " << fileName << endl; abort(); }
  eventList_.clear();
  int nerr = 0;
  while(!feof(f)){
    //      while(0==fscanf(f, " # ")) {/*NOP*/}
    int run;
    int event;
    int n = fscanf(f, " %d %d ", &run, &event);
    if(n!=2){
      ++nerr;
    } else{
      eventList_.push_back((((Long64_t)run) <<RUN_OFFSET) + event);
    }
  }
  sort(eventList_.begin(), eventList_.end());
  if(nerr!=0) cout << nerr << " lines of " << fileName << " were skipped\n";
  if(verbose_>1){
    cout << "Event to select: \n";
    for(size_t i = 0; i < eventList_.size(); ++i){
      cout << "Run " << (eventList_[i] >>32) << " event "
           << (eventList_[i] & 0xFFFFFFFF) << "\n";
    }
  }
}

void CopyAnaTuple::readExludedBranchList(const char* fileName){
  FILE* f = fopen(fileName, "r");
  if(f==0) {
    cout << "Failed to open exluded tree list file "
         << fileName << endl; abort();
  }
  excludedTreeList_.clear();
  excludedBranchList_.clear();
  int nerr = 0;
  while(!feof(f)){
    //      while(0==fscanf(f, " # ")) {/*NOP*/}
    char treeName[256];
    char branchName[256];
    int n = fscanf(f, " %255s %255s ", treeName, branchName);
    if(n!=2){
      ++nerr;
    } else{
      if(strcmp(branchName, "*")==0){//full tree exclusion
        excludedTreeList_.push_back(treeName);
      } else{
        excludedBranchList_.push_back(std::pair<std::string, std::string>(treeName,
                                                                          branchName));
      }
    }
  }
    
  sort(excludedTreeList_.begin(), excludedTreeList_.end());
  sort(excludedBranchList_.begin(), excludedBranchList_.end());
  if(nerr!=0) cout << nerr << " lines of " << fileName << " were skipped\n";
  if(verbose_>1){
    cout << "Tree to exclude from copy: \n";
    for(size_t i = 0; i < eventList_.size(); ++i){
      cout << "Run " << (eventList_[i] >>32) << " event "
           << (eventList_[i] & 0xFFFFFFFF) << "\n";
    }
  }
}

void CopyAnaTuple::setEventSummaryTree(){
  fin_->GetObject("eventSummary", eventSummary_);
  if(eventSummary_==0) { cout << "Event summary tree not found!" << endl; }
  eventSummary_->SetBranchAddress("run", &runNum_);
  eventSummary_->SetBranchAddress("event", &eventNum_);
}

/** Initialize input and output files and trees
 */
void CopyAnaTuple::init(const char* inputDataFile, const char* outputDataFile){
    
  fin_ = std::auto_ptr<TFile>(new TFile(inputDataFile));
    
  if(fin_->IsZombie()){
    cout << "Failed to open input ntuple file "
         << inputDataFile << endl;
    abort();
  }
    
  setEventSummaryTree();
    
  TObjArray* summaryBranches = eventSummary_->GetListOfBranches();
    
  TList* fcontent = fin_->GetListOfKeys();
  if(fcontent==0) {
    cout << "File " << inputDataFile << " is empty! " << endl;
    abort();
  }

  if(outputDataFile){
    fout_ = std::auto_ptr<TFile>(new TFile(outputDataFile, "RECREATE"));
      
    if(fout_->IsZombie()){
      cout << "Failed to open output ntuple file " << outputDataFile << endl;
      abort();
    }
  }

  fin_->cd();
    
  TIter it(fcontent);

  while(it.Next()){
    TKey* key = (TKey*) (*it);
      
    if(strcmp(key->GetClassName(), "TTree")!=0) continue;
    TreeRcd* treeRcd = new TreeRcd;
    treeRcd->tree = (TTree*)  key->ReadObj();

    if(treeRcd->tree==0) { cout << "Failed to read object " << key->GetName() << endl; continue; }

    //eventSummary tree must be clone in last once its branch statuses are set
    if(strcmp(treeRcd->tree->GetName(), "eventSummary")==0) continue;

    //runSummary will be updated and copied at end of the processing
    if(strcmp(treeRcd->tree->GetName(), "runSummary")==0) continue;
      
    bool toCopy = filterTree(treeRcd->tree->GetName());

    if(toCopy){
      //check branch to be copied
      TObjArray* bs = treeRcd->tree->GetListOfBranches();
      TIter itBranches(bs);
      while(itBranches.Next()){
        TBranch* b = (TBranch*)(*itBranches);
        //disable copy of filered-out branches:
        if(!filterBranch(treeRcd->tree->GetName(), b->GetName())){
          if(verbose_) cout << "Disable branch " << b->GetName() << " of tree "
                            << treeRcd->tree->GetName() << "\n";
          treeRcd->tree->SetBranchStatus(b->GetName(), 0);
        }
      }
    }
      
    if(strncmp(treeRcd->tree->GetName(), "ntu_", 4)==0){
      char buf[256];
      strcpy(buf, "n_");
      strncpy(buf + 2, treeRcd->tree->GetName() + 4, sizeof(buf) - 2);
      buf[sizeof(buf)-1] = 0;

      const TBranch* b = dynamic_cast<TBranch*>(summaryBranches->FindObject(buf));

      if(!toCopy && b){
        //disable branch processing and copy:
        if(verbose_) cout << "Disable branch " << buf << " of tree "
                          << eventSummary_->GetName() << "\n";
        eventSummary_->SetBranchStatus(buf, 0);
      } else{
        if(b==0){
          cout << "Tree " << treeRcd->tree->GetName()
               << " not found in event summary " << endl;
          toCopy = false;
        } else{
          eventSummary_->SetBranchAddress(buf, &(treeRcd->nelts));
        }
      }
    }
      
    if(toCopy){
      trees_.Add(treeRcd);
      if(fout_.get()){
        fout_->cd();    
        treeRcd->outTree = treeRcd->tree->CloneTree(0);
        fin_->cd();
      }
    }
  } //next tree

    //Clone eventSummary tree
  if(fout_.get()){
    fout_->cd();
    TreeRcd* treeRcd = new TreeRcd;
    treeRcd->tree = eventSummary_;
    treeRcd->outTree = treeRcd->tree->CloneTree(0);
    trees_.Add(treeRcd);
    fin_->cd();
  }
}
  
void CopyAnaTuple::copyEvent(){
  if(verbose_ > 1) cout << "Copying event " << eventNum_ << " of run " << runNum_ << endl;
  TIter itTree(&trees_);
  while(itTree.Next()){
    TreeRcd* tree = (TreeRcd*) *itTree;
    if(verbose_ > 1) cout << "\ttree " << tree->tree->GetName() << endl;
    for(Int_t ielt = 0; ielt < tree->nelts; ++ielt){
      if(verbose_ > 1) cout << "\t\tcopy object " << (ielt+1) << endl;
      fin_->cd();
      tree->tree->GetEntry(tree->begin + ielt);
      fout_->cd();
      tree->outTree->Fill();
    }
  }
}
  
bool CopyAnaTuple::nextEvent(){
  ++ievent;
  //update tree event startpointers:
  if(ievent>0){
    TIter it(&trees_);
    while(it.Next()) ((TreeRcd*)(*it))->begin += ((TreeRcd*)(*it))->nelts;
  }
  return eventSummary_->GetEntry(ievent)?true:false;
}
  
void CopyAnaTuple::processTree(TTree* tree){

  fout_->cd();    
  TTree* outTree = tree->CloneTree(0);
  fin_->cd();
  const Long64_t nentries = tree->GetEntriesFast();
  for(Long64_t i = 0; i < nentries && i < 2; ++i){
    tree->GetEntry(i);
    outTree->Fill();
  }
  outTree->AutoSave();
}
