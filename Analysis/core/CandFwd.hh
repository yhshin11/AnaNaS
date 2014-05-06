#ifndef CandFwd_h
#define CandFwd_h

#include <cassert>
#include <vector>
#include <map>

//
// useful global names
//
class Candidate;

enum AngleUnit   { kRad=0,      kDeg };

typedef std::vector< Vertex* > VtxList;
typedef VtxList::const_iterator VtxListIterator;

typedef std::vector< Candidate* > CandList;
typedef CandList::const_iterator CandListIterator;

typedef std::multimap< Candidate*, Candidate* > CandMap;
typedef CandMap::const_iterator CandMapIterator;

typedef std::vector< size_t > CandIdList;
typedef CandIdList::const_iterator CandIdListIterator;

typedef std::multimap< size_t, size_t > CandIdMap;
typedef CandIdMap::const_iterator CandIdMapIterator;

struct isoDep
{
  float dr;
  float val;
  float eta;
  float phi;
};

typedef std::multimap< Candidate*, isoDep > CandIsoMap;
typedef CandIsoMap::const_iterator CandIsoMapIterator;

typedef std::multimap< size_t, isoDep > CandIdIsoMap;
typedef CandIdIsoMap::const_iterator CandIdIsoMapIterator;

//typedef std::multimap< size_t, Candidate* > CandIdEcalMap;
//typedef CandIdEcalMap::const_iterator CandIdEcalMapIterator;

struct antiSpike
{
  int rFlag;
  float chi2;
  float ootChi2;
  float ootE;
  //  int sevLvl;
  //  float E1oE9;
  //  float swCross;
  //  float dbStatus;  
};

struct laser
{
  float laserCorr;
  float apdpnref;
  float apdpn_p2;
  float apdpn_t2;
  float alpha;
};

struct ecalRecHit
{
  float e;
  int ix;
  int iy;
  int iz;
  int index;
  int time;

  antiSpike SpikeVar;
  laser LaserVar;

};

struct ecalPSRecHit
{
  float e;
  int ix;
  int iy;
  int iz;
  int index;
  int time;
  int plane;
  int strip;

};

typedef std::vector< ecalRecHit > ecalRecHitList;
typedef ecalRecHitList::const_iterator ecalRecHitListIterator;

typedef std::multimap< size_t, std::pair<ecalRecHit, float> > CandIdRecHitMap;
typedef CandIdRecHitMap::const_iterator CandIdRecHitMapIterator;

typedef std::vector< ecalPSRecHit > ecalPSRecHitList;
typedef ecalPSRecHitList::const_iterator ecalPSRecHitListIterator;

typedef std::multimap< size_t, std::pair<ecalPSRecHit, float> > CandIdPSRecHitMap;
typedef CandIdPSRecHitMap::const_iterator CandIdPSRecHitMapIterator;

typedef std::pair<size_t,float> indexDRPair;

struct SortIndexDRPairs
{
  bool operator()( const indexDRPair& p1, const indexDRPair& p2 ) const 
  {
    double dr1 = p1.second;
    double dr2 = p2.second;
    if( dr1>0 && dr2>0 )
      return dr1 < dr2;
    else if( dr1>0 ) return false;
    else if( dr2>0 ) return true;
    return fabs(dr1) < fabs(dr2);
  } 
};
  
typedef std::multimap< size_t, indexDRPair > CandIdDrMap; 
typedef CandIdDrMap::const_iterator CandIdDrMapIterator;


#endif
