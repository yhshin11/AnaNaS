//==========================================================================
// File and Version Information:
// 	$Id: PrintTree.cc,v 1.1 2009/11/24 12:11:31 ghm Exp $
//
//--------------------------------------------------------------------------
// Description:
//	See PrintTree.h
//
//--------------------------------------------------------------------------
// Sample User Code:
//
// Author List:
//	Abi Soffer              (Original author)
//
//--------------------------------------------------------------------------
// Copyright Information:
//	Copyright (C) 1998	Colorado State University
//
// ROOT Version by Marcel Kunze, RUB
// HIGGS Version by Sergey Ganzhur, CEA-Saclay, Dapnia/SPP
//
// SacAnalysis version by Gautier HAMEL de MONCHENAULT, IRFU/SPP
//
//==========================================================================

//---------------
// C++ Headers --
//---------------
#include <sstream>
#include <iomanip>

using namespace std;

#include <TList.h>
#include <TObjString.h>

//-----------------------
// This Class's Header --
//-----------------------

#include "Analysis/core/PrintTree.hh"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "Analysis/utils/ParticleData.hh"
#include "Analysis/core/Candidate.hh"

//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------

ClassImp(PrintTree)

//-----------------------------------------------
//-- Static Data & Function Member Definitions --
//-----------------------------------------------
const char* 
PrintTree::PrintName(const Candidate & cand) {
    TString result;
    result += cand.name();
    return result.Data();
}

//------------------------------------------------------------------
const char* 
PrintTree::PrintP4(const Candidate & cand) {
    ostringstream stream;
    TLorentzVector q = cand.p4();
    // Reproduce the LorentzVector printout since the effect of setw(..) is lost
    // after each floating point output:
    stream << "(" 
	<<  q.X() << "," 
	<<  q.Y() << ","
	<<  q.Z() << ";" 
	<<  q.T() << ")" << ends;
    
    return stream.str().c_str();
}

//------------------------------------------------------------------
const char* 
PrintTree::PrintVertex(const Candidate & cand) {
    TString result;
    TVector3 q = cand.pos();
    ostringstream stream;
    stream << "(" 
	   << q.X() << "," 
	   << q.Y() << "," 
	   << q.Z() << ")" << ends;
    
    result += stream.str().c_str();
    return result.Data();
}

//----------------------------------------
//-- Public Function Member Definitions --
//----------------------------------------

//----------------
// Constructors --
//----------------
//-----------------------------------------------------------
PrintTree::PrintTree() {
    Initialize();
    fTerm = 0;
}

//--------------
// Destructor --
//--------------
PrintTree::~PrintTree() {
    _lines.Delete();  
}


//-------------
// Accessors --
//-------------
//-----------------------------------------------------------
TStringLong
PrintTree::Print(const Candidate & cand) {
    MakeLines(cand);
    return PrintAgain();
}

//--------------------------------------------------------------------  
TStringLong 
PrintTree::PrintAgain() const {
    TStringLong result;
    for (int l = 0; l < _lines.GetSize(); ++l){
	TObjString *string = (TObjString *) _lines.At(l);
	result = result + string->GetString()  + TString("\n");
    }
    result = result + TString("\0");
    return result;
}

//-------------
// Modifiers --
//-------------
void 
PrintTree::SetSeparator(const char * div) {
    _separator = div;
}

//-----------------------------------------------------------
void 
PrintTree::SetLastMark(const char * div) {
    _lastMark = div;
}

//-----------------------------------------------------------
void 
PrintTree::SetAligner(const char * div) {
    _aligner = div;
}

//-----------------------------------------------------------
void 
PrintTree::SetSpaceFiller(char spa) {
    _spaceFiller = spa;
}

//-----------------------------------------------------------
void 
PrintTree::SetMissingInfo(const char * info) {
    _missingInfo = info;
}

//-------------------------------------------------------------
void 
PrintTree::AddTerminal(Int_t particle) {
    if (fTerm<100) _terminals[fTerm++] = particle;
}

void 
PrintTree::AddTerminal(TString particles) {
    // Don't print the daughters of terminal particles:
    istringstream terminalsStream((char*)particles.Data());
    char terminal[128];
    while (terminalsStream >> terminal){
      const TParticlePDG * entry = ParticleData::Table()->GetParticle(terminal);
	if (entry != 0){
	    AddTerminal(entry->PdgCode());
	}
    }
}

//-------------------------------------------------------------
void 
PrintTree::RemoveTerminal(Int_t particle) {
}


//-----------
// Globals --
//-----------

//-------------------------------------------
//-- Protected Function Member Definitions --
//-------------------------------------------

void 
PrintTree::MakeLines(const Candidate & cand) {
    _lines.Delete();  
    _daughtersLineNum = 0;
    MakeLines(cand, 0, kFALSE, kTRUE, 0);
}

//-------------------------------------------------------------
void
PrintTree::MakeLines(const Candidate & cand, int lineNumber,
			Bool_t firstDaughter, Bool_t lastDaughter,
			int numInitialSpaces) {
    
    // See if the line number exists in _lines. If not, new it:
    TObjString *line = 0;
    if (_lines.GetSize() > lineNumber){
	line = (TObjString*)_lines.At(lineNumber);
    }
    else{
	line = new TObjString();
    }
    
    TString stringLine = line->GetString();
    
    // If cand is the first daughter, align it with is parent.
    // Otherwise, push it beyond the longest line below it:
    if (kTRUE == firstDaughter){
	// Regardless of _spaceFiller, always use ' ' for first daughter:
	char localSpaceFiller = ' ';
	
	stringLine = stringLine + SpaceFillerString(numInitialSpaces - stringLine.Length(), 
	    localSpaceFiller);
    }
    else {
	// Use spaceFillers only of cand needs space for daughters in next line:
	if (cand.nDaughters() > 0){
	    int longestLineLength = 0;
	    for (int l = lineNumber; l < _lines.GetSize(); ++l){
		if (longestLineLength < (((TObjString*)_lines.At(l))->GetString()).Length()){
		    longestLineLength = (((TObjString*)_lines.At(l))->GetString()).Length();
		}
	    }
	    stringLine = stringLine + SpaceFillerString(longestLineLength - stringLine.Length(), 
		_spaceFiller);
	}
    }
    
    // If first daughter, indicate decay from parent using aligner:
    if (kTRUE == firstDaughter){
	stringLine = stringLine + TString(_aligner);
    }
    
    // numInitialSpaces for next recursive calls:
    const int nextNumInitialSpaces = stringLine.Length();
    
    // print the information of interest on the line:
    TString funcOutput;
    if (funcOutput==""){
	///		funcOutput += _missingInfo;
      funcOutput += cand.name();
    }
    
    // Add separators:
    if (kTRUE == lastDaughter) {
	stringLine = stringLine + funcOutput + TString(_lastMark);
    }
    else {
	stringLine = stringLine + funcOutput + TString(_separator);
    }
    
    // If line is not in _lines yet, append it to the list:
    if (_lines.GetSize() <= lineNumber){
	_lines.Add(new TObjString(stringLine));
    }
    else
	line->SetString((char*)stringLine.Data());
    
    // If cand is a terminal particle, ignore its daughters:
    const TParticlePDG * candEntry = cand.pdtEntry();
    if (0 != candEntry){
	const int candType = candEntry->PdgCode();
	for (UInt_t t = 0; t < fTerm; ++t){
	    if (_terminals[t] == candType){
		return;
	    }
	}
    }
    
    // Otherwise, recursively call this function for daughters:
    Bool_t nextFirstDaughter = kTRUE;
    int nextLineNumber = lineNumber + 1;

    size_t sisterNum = 0;
    
    if(cand.nDaughters()>0) {
      for ( size_t idau=0; idau<cand.nDaughters(); idau++ )
	{
	  const Candidate* daughter = cand.daughter(idau);
	  Bool_t nextLastDaughter = kFALSE;
	  if (++sisterNum == cand.nDaughters()){
	    nextLastDaughter = kTRUE;
	  }
	  
	  MakeLines(*daughter, nextLineNumber, nextFirstDaughter, 
		    nextLastDaughter, nextNumInitialSpaces);
	  nextFirstDaughter = kFALSE;
	}
    }
}


//-----------------------------------------------------------
TString 
PrintTree::SpaceFillerString(int length, char spaceFiller) {
    TString result;
    if (length <= 0){
	return result;
    }
    else {
	char * spaces = new char[length + 1];  
	for (int s = 0; s < length; ++s){
	    spaces[s] = spaceFiller;
	}
	spaces[length] = '\0';
	result += spaces;
	delete[] spaces;
	return result;
    }
}

//-----------------------------------------
//-- Private Function Member Definitions --
//-----------------------------------------
// Do work common to different constructors:
void 
PrintTree::Initialize() {
    _aligner = "|_ ";
    _separator = " ";
    _lastMark = " ";
    _spaceFiller = '-';
    _missingInfo = "NoInfo";
    _daughtersLineNum = 0;
}

//-----------------------------------
//-- Internal Function Definitions --
//-----------------------------------
