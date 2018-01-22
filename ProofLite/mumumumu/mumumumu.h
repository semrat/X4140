//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Jun 18 14:15:53 2017 by ROOT version 5.34/26
// from TTree X_data/X(4140) Data
// found on file: MuOniaParked_Run2012C_mumumumuPAT_merged0.root
//////////////////////////////////////////////////////////

#ifndef mumumumu_h
#define mumumumu_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

#include <TSystem.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TBranch.h>
//#include <TCint.h>
#include <TRandom.h>
#include <TMath.h>
#include <TDirectory.h>
#include "TEnv.h"
#include <TString.h>
#include <TSelector.h>
#include <TProof.h>
#include <TProofOutputFile.h>
#include <TLorentzVector.h>
#include "TPoint.h"
#include <TH1.h>
#include <TH2.h>
#include <TH2F.h>
#include <TF1.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <tuple>
#include <map>

// Header file for the classes stored in the TTree if any.
//#include "/build/bellenot/source/root_v5.34.26/root/cint/cint/lib/prec_stl/vector"
//#include "/build/bellenot/source/root_v5.34.26/root/cint/cint/lib/prec_stl/vector"
//#include "/build/bellenot/source/root_v5.34.26/root/cint/cint/lib/prec_stl/vector"
//#include "/build/bellenot/source/root_v5.34.26/root/cint/cint/lib/prec_stl/vector"
//#include "/build/bellenot/source/root_v5.34.26/root/cint/cint/lib/prec_stl/vector"
//#include "/build/bellenot/source/root_v5.34.26/root/cint/cint/lib/prec_stl/vector"

// Fixed size dimensions of array or collections stored in the TTree if any.

class mumumumu : public TSelector {
public :
  ////////////////////////////////////////////////////////////
  //INPUT TREE
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   // Declaration of leaf types

   UInt_t          evtNum;
   UInt_t          runNum;
   UInt_t          lumiNum;
   UInt_t          trigFlag;

   Float_t         xM;
   Float_t         phiM;
   Float_t         jPsiM;
   Float_t         vProb;
   Float_t         cosAlpha;
   Float_t         lxy;
   Float_t         lxyErr;
   Float_t         vNChi2;
   Float_t         dz;

   UInt_t          phiTrigger;
   UInt_t          jpsiTrigger;

   TLorentzVector xP4;

   // List of branches

   TBranch*         b_evtNum;
   TBranch*         b_runNum;
   TBranch*         b_lumiNum;
   TBranch*         b_trigFlag;

   TBranch*         b_xM;
   TBranch*         b_phiM;
   TBranch*         b_jPsiM;
   TBranch*         b_vProb;
   TBranch*         b_cosAlpha;
   TBranch*         b_lxy;
   TBranch*         b_lxyErr;
   TBranch*         b_vNChi2;
   TBranch*         b_dz;
   TBranch*         b_xP4;
   TBranch*         b_phiTrigger;
   TBranch*         b_jpsiTrigger;


   ////////////////////////////////////////////////////////////
   mumumumu(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~mumumumu() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   TProofOutputFile *OutFile;
   TFile            *fOut;

   // TTree *outTree;
   TNtuple *outTuple;
   //
   Float_t          run_out;
   Float_t          evt_out;
   Float_t          lum_out;
   Float_t         X_mass;
   Float_t         kk_mass;
   Float_t         mumu_mass;
   Float_t         X_LFly;
   Float_t         X_pt;
   Float_t         X_eta;
   Float_t         X_vtx;
   Float_t         X_cosAlpha;
   Float_t         X_hlt;

  //
  // TBranch*      X_mass_b;
  // TBranch*      kk_mass_b;
  // TBranch*      mumu_mass_b;
  // TBranch*      X_LFly_b;
  // TBranch*      X_pt_b;
  // TBranch*      X_eta_b;
  // TBranch*      X_vtx_b;
  // TBranch*      X_cosAlpha_b;
  // TBranch*      X_hlt_b;
  // TBranch*      X_run_b;
  // TBranch*      X_evt_b;
  // TBranch*      X_lum_b;

   double JPsi_mass;
   double Phi_mass,Phi_sigma,Phi_mean;
   double Bs0_Low_Mass;
   double Bs0_High_Mass;
   double Y_High_Mass;

  //  int muonQual[4] = {1,3,4,12};


   ClassDef(mumumumu,0);
};

#endif

#ifdef mumumumu_cxx
void mumumumu::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   evtNum        = 0;
   runNum        = 0;
   lumiNum       = 0;
   trigFlag      = 0;

   xM          = 0;
   phiM        = 0;
   jPsiM         = 0;
   vProb         = 0;
   cosAlpha        = 0;
   lxy         = 0;
   lxyErr        = 0;
   vNChi2        = 0;
   dz        = 0;

   phiTrigger = 0;
   jpsiTrigger = 0;

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event", &evtNum, &b_evtNum);
   // fChain->SetBranchAddress("run", &runNum, &b_runNum);
   // fChain->SetBranchAddress("lumiblock", &lumiNum, &b_lumiNum);
   // fChain->SetBranchAddress("trigger", &trigFlag, &b_trigFlag);
   // fChain->SetBranchAddress("xM", &xM, &b_xM);
   // fChain->SetBranchAddress("phi_M", &phiM, &b_phiM);
   // fChain->SetBranchAddress("jpsi_M", &jPsiM, &b_jPsiM);
   // fChain->SetBranchAddress("vProb", &vProb, &b_vProb);
   // fChain->SetBranchAddress("cosAlpha", &cosAlpha, &b_cosAlpha);
   // fChain->SetBranchAddress("l_xy", &lxy, &b_lxy);
   // fChain->SetBranchAddress("lErr_xy", &lxyErr, &b_lxyErr);
   // fChain->SetBranchAddress("vNChi2", &vNChi2, &b_vNChi2);
   // fChain->SetBranchAddress("dz", &dz, &b_dz);
   // fChain->SetBranchAddress("x_p4", &xP4);
   // fChain->SetBranchAddress("phi_trigger", &phiTrigger, &b_phiTrigger);
   // fChain->SetBranchAddress("jpsi_trigger", &jpsiTrigger, &b_jpsiTrigger);



}

Bool_t mumumumu::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif
