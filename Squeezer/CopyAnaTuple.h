// Definition of the Squeezer main class
// Author: Ph. Gras CEA/IRFU Saclay
// Jan. 3, 11

#ifndef COPYANATUPLE_H
#define COPYANATUPLE_H

#include <algorithm>
#include <iomanip>
#include <memory>
#include <vector>
			 
#include "TObject.h"
#include "TTree.h"

class TFile;

/**
 * Utility class to copy and skim AnaNaS n-tuple files.
 * User can either use the standalone application squeezer
 * which is based on this class or write its own code
 * derving a customized class from this one.
 */
class CopyAnaTuple{
protected:

  struct TreeRcd: public TObject{
    
    TreeRcd(): nelts(1), begin(0) {}
    TreeRcd(const TreeRcd& r): TObject(r) {}
    virtual ~TreeRcd() {} 
    
    TTree* tree;
    TTree* outTree;
    Int_t  nelts;
    //index in the tree of the first element of current events (CopyAnaTuple::ievent)
    Long64_t begin;
  };

  std::auto_ptr<TFile> fin_;
  std::auto_ptr<TFile> fout_;
  
  TObjArray trees_;

  TTree* eventSummary_;

  Int_t maxEvents_;
  
  Int_t ievent;

  Int_t runNum_;

  Int_t eventNum_;
  
  std::vector<Long64_t> eventList_;

  std::vector<std::string> excludedTreeList_;

  std::vector<std::pair<std::string, std::string> > excludedBranchList_;
  
  static const int RUN_OFFSET = 32;

  bool allEvent_;

  //numbr of copied events:
  int nCopied_;
  
public:
  
  int verbose_;  

  CopyAnaTuple(): maxEvents_(-1), ievent(-1), allEvent_(true), nCopied_(0), verbose_(0){
    trees_.SetOwner(kTRUE);
  }

  virtual ~CopyAnaTuple(){
  }

  void setMaxEvents(int val) { maxEvents_ = val; }

  void listEvents(std::ostream& o, const char* inputDataFile);

  void listBranches(std::ostream& o, const char* inputDataFile);
  
  void fillRunSummary();
  
  void run(const char* inputDataFile, const char* outputDataFile);

  virtual void readEventList(const char* fileName);

  virtual void readExludedBranchList(const char* fileName);
  
protected:
  /** Check if tree must be copied. Return false if tree should not be copied.
   * Can be overwritten in a derived class to customize the selection.
   */
  virtual bool filterTree(const char* treeName){
    return std::find(excludedTreeList_.begin(), excludedTreeList_.end(),
                     std::string(treeName)) == excludedTreeList_.end();
  }

  /** Check if a branch must be copied. Return false the branch must excluded
   * from the copy,
   * Can be overwritten in a derived class to customize the selection.
   */
  virtual bool filterBranch(const char* treeName, const char* branchName){
    return std::find(excludedBranchList_.begin(), excludedBranchList_.end(),
                     std::pair<std::string, std::string>(treeName, branchName))
                     == excludedBranchList_.end();
  }
  
  /** Check if curren event must be used. Typically use eventNum_ and runNum_
   * fields which contain the event and the run numbers. Return true if event
   * must be copied
   */
  virtual bool filterEvent(){
    if(allEvent_) return true; //copy every event
    bool selected = (std::find(eventList_.begin(), eventList_.end(),
                               (((Long64_t) runNum_) <<RUN_OFFSET) + eventNum_)
                     != eventList_.end());
    return selected;
  }

private:
  /** initialized event summary tree. This method is already called by
   * init() method.
   */
  void setEventSummaryTree();

  /** Initialize input and output files and trees
   */
  void init(const char* inputDataFile, const char* outputDataFile = 0);
  
  void copyEvent();
  
  bool nextEvent();
  
private:
  void processTree(TTree* tree);
};

#endif //COPYANATUPLE_H not defined
