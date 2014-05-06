#ifndef DYContamination_h
#define DYContamination_h

#include <vector>
#include <map>
#include <string>

#include "Analysis/core/SampleAnalysis.hh"
#include "Analysis/utils/TreeManager.hh"

class DYContamination : public SampleAnalysis
{
public:

  DYContamination( Sample& sample, const std::string & collectionFileName );
  virtual ~DYContamination();
 
private:

  virtual void bookHistograms();
  virtual bool analyzeEvent();
  virtual void writeHistograms();

  void makeZCandList(CandList &);
  void saveElectron(Candidate *, int, int);

  // HLT decision
  virtual bool passHlt();

  // build the Z->ll solutions
  //  imode=0 (Z->ee) or 1 (Z->mumu)
  size_t buildZeeLists();
  void initRHV();

  size_t _nevts;
  size_t _nSelectedZ;

  std::vector<int> _rhV;

  TreeManager _tm;

  ClassDef( DYContamination, 0 );
};

#endif
