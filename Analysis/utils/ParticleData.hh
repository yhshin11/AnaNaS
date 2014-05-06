#ifndef ParticleData_hh
#define ParticleData_hh
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// ParticleData								//
//                                                                      //
// Read the BABAR style PDT.						//
//                                                                      //
// Author: Marcel Kunze, RUB, 1999-2000					//
// Copyright (C) 1999-2001, Ruhr-University Bochum.			//
//                                                                      //
// Adapted for AnaNaS by Gautier HAMEL de MONCHENAULT (09-2009)         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <TDatabasePDG.h>

class ParticleData : public TDatabasePDG
{
public:

  static ParticleData* Table();


  ParticleData() {}
  ParticleData(const char *file) { ReadPDGTable(file); }
  virtual ~ParticleData() {}
  virtual void ReadPDGTable(const char *filename); //Read the BABAR style pdt.table
  virtual TParticlePDG*   AddParticle(const char*  Name, 
				      const char*  Title,
				      Double_t     Mass, 
				      Bool_t       Stable,
				      Double_t     DecayWidth, 
				      Double_t     Charge, 
				      const char*  ParticleClass,
				      Int_t        PdgCode,
				      Int_t        Anti=-1,
				      Int_t        TrackingCode=0);
  
private:

  static ParticleData* _instance;

  ClassDef(ParticleData, 0)   // BABAR Particle Data Table
};

#endif
