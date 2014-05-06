#ifndef RandomConeSelector_h
#define RandomConeSelector_h

#include <vector>
#include <stdlib.h>
#include "Analysis/core/EventManager.hh"
#include "Analysis/selectors/SelectorBase.hh"
#include "Analysis/utils/Constants.hh"

#include "MECore/src/ME.hh"
#include "MECore/src/MEEBDisplay.hh"
#include "MECore/src/MEEEDisplay.hh"

using namespace std;

class RandomConeSelector : public SelectorBase
{
public:

  RandomConeSelector(vector<int> interest);

  virtual ~RandomConeSelector() {}

  // compute isolation variables
  void computeIso( const Candidate&, const Candidate& , bool ptk=false);
  
  void computeEcalIso  ( const Candidate&, const Candidate& , bool btp=false);
  void computeHcalIso  ( const Candidate& );
  void computeHcalIsoH ( const Candidate& );
  void computeTrackIso ( const Candidate& );

  float getEcalIso()  { return S_ecal; }
  float getHcalIso()  { return S_hcal; }
  float getHcalIsoH() { return S_hcalH; } // hcal with hole
  float getTrackIso() { return S_trk; }
  float getSumIso()   { return S_sum; }

  bool  countsInTkIso( const Candidate*,const Candidate& );

  
  // getting functions
  int   getNecal()    { return n_ecal; }
  int   getNhcal()    { return n_hcal; }
  int   getNhcalH()   { return n_hcalH; }
  int   getNtrack()   { return n_trk; }
  int   getNetk()     { return n_tkecal; }
  
  vector<int>   getVieta() { return ieta; }
  vector<float> getVeta()  { return veta; }
  vector<int>   getViphi() { return iphi; }
  vector<float> getVphi()  { return vphi; }
  vector<float> getVe()    { return ve; }
  vector<float> getVet()   { return vet; }
  
  vector<float> getVtkpt()  { return vtkpt; }
  vector<bool>  getVtkCnt() { return vtkCnt; }
  vector<float> getVtkd0()  { return vtkd0; }
  vector<float> getVtkz0()  { return vtkz0; }
  vector<int>   getVtkNh()  { return vtkNhits; }
  vector<int>   getVtkIh()  { return vtkNxIhits; }
  vector<float> getVtkChi2()  { return vtkChi2; }
  vector<int>   getVtkNdof()  { return vtkNdof; }
  vector<int>   getVtkQ()     { return vtkQ; }
  vector<bool>  getVtkFP()    { return vtkFP; }
  vector<bool>  getVtkSV()    { return vtkSV; }
  vector<float> getVtkEtaH()  { return vtkEtaH; }
  vector<float> getVtkPhiH()  { return vtkPhiH; }
  vector<bool>  getVtkInVeto(){ return vtkInVeto; }
  vector<int>   getVtkVid()   { return vtkVid; }
  vector<float> getVtkDepE()  { return vtkDepE; }
  vector<float> getVtkDepEt() { return vtkDepEt; }
  vector<float> getVtkDEta()  { return vtkDEta; }
  vector<float> getVtkDPhi()  { return vtkDPhi; }
  vector<float> getVtkDR()    { return vtkDR; }
  
  // setting functions
  void  setEBzs  ( float ebZSthresh ) { e0_ecal[0] = ebZSthresh; }
  void  setAbsCut( bool ab )          { absCut = ab; }
  void  setSCCut ( bool scc )         { doSCcut = scc; }
  void  setIC    ( float innerC )     { dR_ecal_in[0] = innerC; }
  void  setStrip ( float strip  )     { strip_ecal[0] = strip; }
  void  loadXtalDatabase();
  int   getTrkGenPdgId(const Candidate*, bool);
  float getTkDeposit( float, float, int );
  void  genList( float, float );
  
public:

  vector<int> RHinterest;
  
  bool  debug;
  bool  absCut;
  bool  useGeometryDatabase;
  bool  doSCcut;
  float r_EB, z_EB, e_EE, r_EE, z_EE, r_HB, z_HE;

  vector<int>   ieta;
  vector<int>   iphi;
  vector<float> veta;
  vector<float> vphi;
  vector<int>   iz;
  vector<float> ve;
  vector<float> vet;
  vector<float> vtkpt;
  vector<bool>  vtkCnt;
  vector<float> vtkd0;
  vector<float> vtkz0;
  vector<int>   vtkQ;
  vector<int>   vtkNhits;
  vector<int>   vtkNxIhits;
  vector<float> vtkChi2;
  vector<int>   vtkNdof;
  vector<bool>  vtkFP;
  vector<bool>  vtkSV;
  vector<float> vtkEtaH;
  vector<float> vtkPhiH;
  vector<float> vtkDEta;
  vector<float> vtkDPhi;
  vector<float> vtkDR;
  vector<bool>  vtkInVeto;
  vector<int>   vtkVid;
  vector<float> vtkDepE;
  vector<float> vtkDepEt;
  vector<int>   usedTrkMatch;
  
  float pt0_trk, dR_trk_in,  dR_trk_out;
  float d0_trk, dz0_trk;
  float strip_trk;
  
  bool  ecalUseIsoDeposits;
  bool  isBarrel;
  int   iecal;
  float et0_ecal[2];
  float dR_ecal_in[2];
  float dR_ecal_out[2];
  float e0_ecal[2];
  float strip_ecal[2];

  float dR_hcal_in, dR_hcal_out;
  float et0_hcal;
  float e0_hcal;

  float S_trk;
  float S_ecal;
  float S_hcal;
  float S_hcalH;
  float S_sum;
  
  int   n_trk;
  int   n_hcal;
  int   n_hcalH;
  int   n_ecal;
  int   n_tkecal;

  bool  isInVeto;
  float niceTkCut;
  
private:

  RandomConeSelector* me() const { return (RandomConeSelector*)this; }

  bool _dbLoaded;
  std::map<int, std::pair<float, float> > MapXtal;
  std::map<int, std::pair<float, float> >::const_iterator itMapXtal;

  ClassDef( RandomConeSelector, 0 )  
};

#endif

     
