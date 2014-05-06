#include "Analysis/utils/ParticleData.hh"

#include <sstream>
#include <iostream>

using namespace std;

#include <THashList.h>
#include <TString.h>

#include "Analysis/utils/Config.hh"

ClassImp(ParticleData)

ParticleData* ParticleData::_instance = 0;

ParticleData*
ParticleData::Table()
{
  if( _instance==0 )
    {
      TString _pdtTable = Config::workPath + "PDT/pdt.table";
      cout << "Instantiating TDatabasePDG from " << _pdtTable << endl;
      _instance = new ParticleData(_pdtTable);
    }
  return _instance;
}

// Read the BABAR table

void ParticleData::ReadPDGTable(const char *filename)
{
   // read list of particles from a file
   // if the particle list does not exist, it is created, otherwise
   // particles are added to the existing list
   // See $ROOTSYS/tutorials/pdg.dat to see the file format

  if (fParticleList == 0) {
     fParticleList = new THashList();
  }

  const Float_t HBARC = 197.327*1.e-3*1.e-13; // GeV*cm

  FILE *ifl = fopen(filename,"r");
  if (ifl == 0) {
    Error("ReadPDGTable","Could not open PDG particle file %s",filename);
    return;
  }

  char line[512];
  while ( fgets(line,512,ifl) ) {
    if (strlen(line) >= 511) {
       Error("ReadPDGTable","input line is too long");
       return;
    }
    istringstream linestr(line);
    TString opcode;
    char subcode;
    linestr >> opcode >> subcode;

    if( opcode == "*" )
      continue;

    if( opcode == "end" )
      break;

    else if( opcode == "add" ) {
      switch (subcode) {
	case 'p':
	  {
	    TString classname;
	    linestr >> classname;
	    // if (classname == "Collision" || classname == "Parton")
	    if (classname == "Collision" )
	      continue;
	
	    TString name;
	    int type;
	    float mass, width, cut, charge, spin, lifetime;
	
	    linestr >> name >> type;
	    linestr >> mass >> width >> cut >> charge;
	    linestr >> spin >> lifetime;
	
	    charge /= 3.0;
	    if (classname != "Meson")
	      spin /= 2.0;
	
	    // lifetime is c*tau (mm)
	    if (lifetime > 0.0 && width < 1e-10)
	      width = HBARC / (lifetime/10.0);

	    Bool_t stable = (lifetime <= 0);

            TParticlePDG *p = new TParticlePDG(name, name, mass, stable, width,
                                      charge, classname.Data(), type, -1, 0);
	    fParticleList->Add(p);
	    break;
	  }
	
	case 'c':
	  {
	    int     ptype, nchild;
	    float   bf;
	    TString decayer;
	
	    linestr >> ptype >> bf >> decayer >> nchild;
	    TParticlePDG *parent = GetParticle(ptype);
	    if (parent == 0) continue;
	
	    TList *kids = new TList();
	
	    int i;
	    for(i=0; i<nchild; i++ )
	    {
	      int ctype;
	      linestr >> ctype;
	      TParticlePDG* secondary = GetParticle(ctype);
	      if( secondary ==0 ) break;
	      kids->Add(secondary);
	    }
	
	    //parent->AddDecay(bf, kids ); // Not yet implemented
	    break;
	  }
	
	case 'd':
	  break;
	
	default:
	  Error("ReadPDGTable","unknown subcode %d for operation add",subcode);
	  break;
	}
     }
  }

  fclose(ifl);
}

//______________________________________________________________________________
TParticlePDG* ParticleData::AddParticle(const char *name, const char *title,
					Double_t mass, Bool_t stable,
					Double_t width, Double_t charge,
					const char* ParticleClass,
					Int_t PDGcode,
					Int_t Anti,
					Int_t TrackingCode)
{
  //
  //  Particle definition normal constructor. If the particle is set to be
  //  stable, the decay width parameter does have no meaning and can be set to
  //  any value. The parameters granularity, LowerCutOff and HighCutOff are
  //  used for the construction of the mean free path look up tables. The
  //  granularity will be the number of logwise energy points for which the
  //  mean free path will be calculated.
  //

  //TParticlePDG* old = GetParticle(PDGcode);

  //if (old) {
  //  printf(" *** TDatabasePDG::AddParticle: particle with PDGcode=%d already defined\n",PDGcode);
  //  return 0;
  //}

  TParticlePDG* p = new TParticlePDG(name, title, mass, stable, width,
				     charge, ParticleClass, PDGcode, Anti,
				     TrackingCode);
  fParticleList->Add(p);
/*
  TParticleClassPDG* pclass = GetParticleClass(ParticleClass);

  if (!pclass) {
    pclass = new TParticleClassPDG(ParticleClass);
    fListOfClasses->Add(pclass);
  }

  pclass->AddParticle(p);
*/
  return p;
}
