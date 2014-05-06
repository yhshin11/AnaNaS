//==========================================================================
// File and Version Information:
// 	$Id: PrintTree.hh,v 1.1 2009/11/24 12:11:31 ghm Exp $
//
//--------------------------------------------------------------------------
// Description:
//	Class PrintTree can print the decay tree of a MC Candidate
//      to a TString (which can then be printed to an ostream or ErrStream).
//
//--------------------------------------------------------------------------
// Collaborating classes:
//
//--------------------------------------------------------------------------
// Sample User Code:
//  Simplest: print aCandidate with defaults: will print the names of all
//  particles in the decay tree with the default format:
//        PrintTree printTree;
//        cout          << printTree.print(aCandidate);
//        anotherStream << printTree.printAgain();    // reprint. saves CPU 
//
// More options: You can change what is printed for each node (Candidate)
// in the decay tree, by specifying the function to do this. The default 
// PrintTree uses the function PrintTree::printName, which prints the
// pdtEntry()->name() of the Candidate. Other functions which are provided
// are PrintTree::printP4, which prints the 4-momentum, and 
// PrintTree::printVertex, which prints the vertex point. There are two 
// ways to specify the output function:
//       PrintTree(anOutputFunc);            // specify in constructor
//       prinTree.setOutputFunc(anOutputFunc);  // or modifier
//
// The typedef PrintTree::OutputFunc defines pointers of such functions:
//       typedef TString (*OutputFunc)(const Candidate & cand, 
//				         const ComIOManip & manip);
//
// The second argument can be used to set the printout format. See 
// CommonUtils/ComIOManip.hh and ComFloatIOManip.hh, for example. The 
// ComIOManip to be used in a PrintTree is specified using
//       void setFormat(const ComIOManip & format);
//
// If you write your own OutputFunc, follow the guidelines given below.
//
// It is often useful to not print out daughters of some particles which would
// just crowd the output. In that case you let the PrintTree know that such
// particles types are "terminal":
//       prinTree.addTerminal(PdtLund::pi0);
//       prinTree.addTerminal(PdtLund::eta);
//
// To undo the last call and show eta decays, call
//       prinTree.removeTerminal(PdtLund::eta);
//
// Additional decay formatting:
// Change the character that fills spaces between sisters (white spaces make
// it hard to keep track of long straight lines, hiding sisterhood relations):
//       printTree.setSpaceFiller('.');  // default is '-'
//
// Change the string that aligns mothers and daughters to demonstrate their 
// relationship:
//       printTree.setAligner("-> ");    // default is "|_ "
//
// Change the string that follows the OutputFunc output for each particle:
//       printTree.setSeparator(",");    // default is " "
//
// Change the string that ends a daughter list:
//       printTree.setLastMark(".");     // default is " "
//
// Here is what the default settings look like in a complete printout with 
// default settings:
//
// B- 
// |_ D0 -------rho- 
//    |_ K- pi+ |_ pi- pi0 
//
// The "|_ " is the aligner, the "---" before the rho- are spaceFillers, the 
// " " following each particle is the separator, and the empty string
// at the end of each daughter list is the lastMark.
// 
// If a particle in the decay tree is missing the information needed, the 
// string missingInfo() appears instead of that information. The default is 
// "NoInfo", and you can change it using
//       printTree.setMissingInfo("?")
//
//--------------------------------------------------------------------------
// OutputFunc Writing Guidelines:
//   If you write your own OutputFunc:
//   1) If an OutputFunc gets a Candidate that lacks the information needed
//      to print it out, the OutputFunc must return an empty string. The 
//      PrintTree will take care of appending the _missingInfo string to 
//      that empty string.
//   2) When using a ComIOManip, note that the width of floating point output,
//      set by the manipulator setw(int) lasts for only one floating point
//      number printout, and is reset afterwards. See how vectors and points
//      are printed in PrintTree:printP4(..) and printVertex(..).
//
//--------------------------------------------------------------------------
// Compiler Default Member Functions:
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
//==========================================================================

#ifndef HPRINTTREE_H
#define HPRINTTREE_H


//---------------
// C++ Headers --
//---------------
#include "TObject.h"
#include "TStringLong.h"
#include "TParticlePDG.h"

//----------------------
// Base Class Headers --
//----------------------

//-------------------------------
// Collaborating Class Headers --
//-------------------------------

//--------------------------------------------
// Collaborating Class Forward Declarations --
//--------------------------------------------
class Candidate;


//		---------------------
// 		-- Class Interface --
//		---------------------

class PrintTree : public TObject {
    //----------------
    // Static members:
    //----------------
public:
    
    // Some pre-provided OutputFunc functions.
    // The function printName is set by the defult constructor:
    static const char* PrintName(const Candidate &);
    
    static const char* PrintP4(const Candidate &);
    
    static const char* PrintVertex(const Candidate &);
    
    
    //------------------
    // Instance members:
    //------------------
public:
    // Constructors:
    // Takes a function that prints a Candidate to a TString:
    PrintTree();
    
    // Destructor
    virtual ~PrintTree();
    
    PrintTree(const PrintTree &x) {}
    
    // Print operations:
    // print(cand) traverses the decay tree of cand:
    TStringLong Print(const Candidate & cand);
    // printAgain() just returns the string created by the last print(cand) call:
    TStringLong PrintAgain() const;
    
    // Accessors (const)
    inline const char* Aligner() const {return _aligner;}
    inline const char* Separator() const {return _separator;}
    inline const char* LastMark() const {return _lastMark;}
    inline char SpaceFiller() const {return _spaceFiller;}
    
    // Accessors to formatting information:
    inline const char* MissingInfo() const {return _missingInfo;}
    
    // Modifiers
    void SetAligner(const char * aligner);
    void SetSeparator(const char * separator);
    void SetLastMark(const char * lastMark);
    void SetSpaceFiller(char spaceFiller);
    
    void AddTerminal(Int_t particle);
    void AddTerminal(TString particles);
    void RemoveTerminal(Int_t particle);
    
    // Modifiers of formatting information:
    void SetMissingInfo(const char * missingInfo);  
    
protected:
    // helpers:
    void MakeLines(const Candidate & cand);
    
    void MakeLines(const Candidate & cand, int lineNumber,
	Bool_t firstDaughter, Bool_t lastDaughter, 
	int numInitialSpaces);
    
    TString SpaceFillerString(int length, char spaceFiller);
    
    
private:
    // Functions:
    void Initialize(); // carries out work for constructors
    
    // Data members:
    // Things the user sets (or uses defaults):
    TString _aligner; //!Do not stream
    TString _separator; //!Do not stream
    TString _lastMark; //!Do not stream
    char _spaceFiller; //!Do not stream
    TString _missingInfo; //!Do not stream
    Int_t _terminals[100]; //!Do not stream
    UInt_t fTerm; //!Do not stream
    
    // Hidden variables used to keep track of printing:
    int _daughtersLineNum; //!Do not stream
    TList _lines; //!Do not stream
    
public:
    ClassDef(PrintTree,0) // Print a tree of candidates
};

#endif
