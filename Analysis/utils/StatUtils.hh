#ifndef StatUtils_hh
#define StatUtils_hh

#include <TObject.h>

class StatUtils 
{
public:
  virtual ~StatUtils(){}
  
  //return a binomial error
  static float EBinom(float eff, int N);

  ClassDef( StatUtils, 0 )

};
#endif
