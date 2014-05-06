#include "Analysis/utils/StatUtils.hh"

#include <cmath>
#include <iostream>
#include <cassert>
using namespace std;

#include "Analysis/utils/Constants.hh"



float 
StatUtils::EBinom(float eff, int N) {

  float error;
  if(eff> 1)
    error = sqrt(eff*(100-eff)/N) ;
  else
    error = sqrt(eff*(1-eff)/N) ;
  
  return error;

}
