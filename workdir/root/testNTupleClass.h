//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu May  1 12:10:21 2014 by ROOT version 5.34/01
// from TTree testNTuple/testNTuple
// found on file: Toy/test.root
//////////////////////////////////////////////////////////

#ifndef testNTupleClass_h
#define testNTupleClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class testNTupleClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Bool_t          boolTightIso_l1;
   Int_t           pdgCode_l1;
   Float_t         pt_l1;
   Float_t         eta_l1;
   Float_t         phi_l1;
   Float_t         charge_l1;
   Bool_t          boolTightIso_l2;
   Int_t           pdgCode_l2;
   Float_t         pt_l2;
   Float_t         eta_l2;
   Float_t         phi_l2;
   Float_t         charge_l2;
   Bool_t          boolValidZ;
   Float_t         m_ll;
   Float_t         pt_ll;
   Float_t         phi_ll;
   Float_t         y_ll;
   vector<float>   *pt_jets;
   vector<float>   *eta_jets;
   vector<float>   *phi_jets;
   vector<bool>    *mvaIdMedium_jets;
   vector<float>   *CSVBT_jets;
   vector<float>   *CSVMVA_jets;
   Float_t         pt_met;
   Float_t         phi_met;
   Float_t         sig_met;

   // List of branches
   TBranch        *b_boolTightIso_l1;   //!
   TBranch        *b_pdgCode_l1;   //!
   TBranch        *b_pt_l1;   //!
   TBranch        *b_eta_l1;   //!
   TBranch        *b_phi_l1;   //!
   TBranch        *b_charge_l1;   //!
   TBranch        *b_boolTightIso_l2;   //!
   TBranch        *b_pdgCode_l2;   //!
   TBranch        *b_pt_l2;   //!
   TBranch        *b_eta_l2;   //!
   TBranch        *b_phi_l2;   //!
   TBranch        *b_charge_l2;   //!
   TBranch        *b_boolValidZ;   //!
   TBranch        *b_m_ll;   //!
   TBranch        *b_pt_ll;   //!
   TBranch        *b_phi_ll;   //!
   TBranch        *b_y_ll;   //!
   TBranch        *b_pt_jets;   //!
   TBranch        *b_eta_jets;   //!
   TBranch        *b_phi_jets;   //!
   TBranch        *b_mvaIdMedium_jets;   //!
   TBranch        *b_CSVBT_jets;   //!
   TBranch        *b_CSVMVA_jets;   //!
   TBranch        *b_pt_met;   //!
   TBranch        *b_phi_met;   //!
   TBranch        *b_sig_met;   //!

   testNTupleClass(TTree *tree=0);
   virtual ~testNTupleClass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef testNTupleClass_cxx
testNTupleClass::testNTupleClass(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Toy/test.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("Toy/test.root");
      }
      f->GetObject("testNTuple",tree);

   }
   Init(tree);
}

testNTupleClass::~testNTupleClass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t testNTupleClass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t testNTupleClass::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void testNTupleClass::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   pt_jets = 0;
   eta_jets = 0;
   phi_jets = 0;
   mvaIdMedium_jets = 0;
   CSVBT_jets = 0;
   CSVMVA_jets = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("boolTightIso_l1", &boolTightIso_l1, &b_boolTightIso_l1);
   fChain->SetBranchAddress("pdgCode_l1", &pdgCode_l1, &b_pdgCode_l1);
   fChain->SetBranchAddress("pt_l1", &pt_l1, &b_pt_l1);
   fChain->SetBranchAddress("eta_l1", &eta_l1, &b_eta_l1);
   fChain->SetBranchAddress("phi_l1", &phi_l1, &b_phi_l1);
   fChain->SetBranchAddress("charge_l1", &charge_l1, &b_charge_l1);
   fChain->SetBranchAddress("boolTightIso_l2", &boolTightIso_l2, &b_boolTightIso_l2);
   fChain->SetBranchAddress("pdgCode_l2", &pdgCode_l2, &b_pdgCode_l2);
   fChain->SetBranchAddress("pt_l2", &pt_l2, &b_pt_l2);
   fChain->SetBranchAddress("eta_l2", &eta_l2, &b_eta_l2);
   fChain->SetBranchAddress("phi_l2", &phi_l2, &b_phi_l2);
   fChain->SetBranchAddress("charge_l2", &charge_l2, &b_charge_l2);
   fChain->SetBranchAddress("boolValidZ", &boolValidZ, &b_boolValidZ);
   fChain->SetBranchAddress("m_ll", &m_ll, &b_m_ll);
   fChain->SetBranchAddress("pt_ll", &pt_ll, &b_pt_ll);
   fChain->SetBranchAddress("phi_ll", &phi_ll, &b_phi_ll);
   fChain->SetBranchAddress("y_ll", &y_ll, &b_y_ll);
   fChain->SetBranchAddress("pt_jets", &pt_jets, &b_pt_jets);
   fChain->SetBranchAddress("eta_jets", &eta_jets, &b_eta_jets);
   fChain->SetBranchAddress("phi_jets", &phi_jets, &b_phi_jets);
   fChain->SetBranchAddress("mvaIdMedium_jets", &mvaIdMedium_jets, &b_mvaIdMedium_jets);
   fChain->SetBranchAddress("CSVBT_jets", &CSVBT_jets, &b_CSVBT_jets);
   fChain->SetBranchAddress("CSVMVA_jets", &CSVMVA_jets, &b_CSVMVA_jets);
   fChain->SetBranchAddress("pt_met", &pt_met, &b_pt_met);
   fChain->SetBranchAddress("phi_met", &phi_met, &b_phi_met);
   fChain->SetBranchAddress("sig_met", &sig_met, &b_sig_met);
   Notify();
}

Bool_t testNTupleClass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void testNTupleClass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t testNTupleClass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef testNTupleClass_cxx
