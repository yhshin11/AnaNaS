#ifndef Sample_h
#define Sample_h

#include <cassert>
#include <iostream>
#include <list>
#include <vector>
#include <map>
#include <string>

#include <TString.h>
#include <TH1.h>

class Model;

//
// a data sample
//

class Sample      
{
public:

  typedef std::vector<Sample*> SampleList;

  struct SortByLumi
  {
    bool operator()( const Sample* s1, const Sample* s2 ) const
    {
      float l1 = (s1!=0) ? s1->integratedLumi() : 0;
      float l2 = (s2!=0) ? s2->integratedLumi() : 0;
      return l1>l2;
    }
  };

  struct SortBySigma
  {
    bool operator()( const Sample* s1, const Sample* s2 ) const
    {
      float l1 = (s1!=0) ? s1->sigma() : 0;
      float l2 = (s2!=0) ? s2->sigma() : 0;
      return l1>l2;
    }
  };

  // constructor for base sample
  Sample( const char* name, const char* coll, const char* comments=0,
	  float sig_tot=0., float filter_eff=0., 
	  float br=0., float mean_kf=0, float n_tot=0 );

  // constructor for composite sample
  Sample( const char* name, SampleList& sampleList );

  // constructor from a file name
  Sample( const char* filename=0 );

  // copy constructor
  Sample( const Sample& );

  virtual ~Sample() {}

  // name of the sample
  const char* name() const { return _name; }
  void setName( const char* name_ ) { _name=name_; }

  // this is data
  bool isData() const { return _isData; }
  void setAsData() { _isData=true; }

  // name of the current subsample
  const char* subsample() const { return _subsample; }

  // total number of events in the original collection
  float n_tot()            const;

  // filter efficiency at generation
  float filter_eff()       const;

  // number of events in ntuples
  float n()                const;

  // number of processed events in ntuples
  float n_proc()                const;

  // the ratio of n ever n_tot
  float trig_eff()         const;

  // if applicable, branching ratio factor
  float br()               const;

  // total cross section 
  float sigma_tot()        const;

  // total cross section * filter efficiency * BR factor 
  float sigma()            const;

  // equivalent integrated luminosity (= n_tot/sigma)
  float integratedLumi()   const;
  
  // weight for a given integrated luminosity (=L/integratedLumi)
  float w(     float L ) const;

  // number of events for a given integrated luminosity
  float n_tot( float L ) const;
  float n(     float L ) const;

  // one-line printing
  void print( ostream& ) const;

  // list of rootfiles to feed the analyses
  std::vector< std::string > filenames;  

  void setSubsample( const char* subsample ) { _subsample = subsample; }

  //
  // for a given histogram (histName) in a given analysis (analysisNAme),
  // for a given luminosity (L)
  //

  // fill a map of normalized histograms indexed by sample
  // (contains the sample and all its descendants)
  void getMapOfHists( std::map< Sample*, TH1* >& mapOfHists,
		      float L,
		      const char* analysisName, 
		      const char* subSampleName,
		      const char* histName
		      );

  // fill the map described above
  //    plus the list of all possible stacked histograms
  void getListOfStackedHists( std::vector< TH1* >& listOfStackedHists,
			      std::map< Sample*, TH1* >& mapOfHists,
			      float L,
			      const char* analysisName, 
			      const char* subSampleName,
			      const char* histName
			      );
  TH1* getHist( float L, 
		const char* analysisName, 
		const char* subSampleName,
		const char* histName
		);

  TH1* getHist( const char* analysisName,
		const char* subSampleName,
		const char* histName 
		);

  // fill a map of normalized histograms indexed by sample
  // (contains the sample and all its descendants)
  void getMapOfHists( std::map< Sample*, TH1* >& mapOfHists,
		      float L,
		      const char* analysisName, 
		      const char* histName
		      ) 
  { 
    return getMapOfHists( mapOfHists, L, analysisName, subsample(), histName ); 
  }

  // fill the map described above
  //    plus the list of all possible stacked histograms
  void getListOfStackedHists( std::vector< TH1* >& listOfStackedHists,
			      std::map< Sample*, TH1* >& mapOfHists,
			      float L,
			      const char* analysisName, 
			      const char* histName
			      )
  { 
    return getListOfStackedHists( listOfStackedHists, mapOfHists, L, analysisName, subsample(), histName ); 
  }

  TH1* getHist( float L, 
		const char* analysisName, 
		const char* histName
		)
  {
    return getHist( L, analysisName, subsample(), histName ); 
  }

  TH1* getHist( const char* analysisName,
		const char* histName 
		)
  {
    return getHist( analysisName, subsample(), histName );
  }

  //
  // static helpers
  //
  static bool exists( const char* sampleName );
  static int checkDataFiles( const char* sampleName );  
  static std::vector< std::string > getDataFiles( const char* sampleName );  
  static int checkNProcFromDataFiles( const char* sampleName );  

  static float checkRootFile( const char* analysisName,
			      const char* sampleName, 
			      const char* subSampleName,
			      const char* histName     );

  static std::map< TString, Sample* >& samples( const char* sampleFile=0 );
  static void    set( const char* sampleFile );
  static Sample* get( const char* sampleName, const char* sampleFile=0 );

  static TH1* getClone( const char* analysisName,
			const char* sampleName, 		        
			const char* subSampleName,
			const char* histName
			);

  // for checks
  static bool verbose;
  static bool check;
  static bool useRootFile;
  
protected:

  TString _name;         // simple string 
  TString _coll;         // collection name
  TString _comments;     // collection name
  float   _n_tot;        // total number of events in collection
  float   _n;            // total number of events in ntuples
  float   _n_proc;       // total number of processed events in ntuples
  float   _sig_tot;      // total cross-section (pb)
  float   _filter_eff;   // filter efficiency
  float   _br;           // product of all branching fractions
  float   _mean_kf;      // mean K factor

  bool isComposite() const;
  SampleList _sampleList;
  std::vector< float > _weights;

  TString _subsample;    // name of the subsample

  void setDataFiles_();

  // determine the number of event from the root files
  void setNFromFiles_();

  static std::map< TString, Sample* > _samples;

  // helper
  // returns a vector of indices vectors
  static std::vector< std::vector<size_t> > combinatorics_( size_t n );

  // data sample ?
  bool _isData;

  ClassDef( Sample, 0 )  
};


#endif

     
