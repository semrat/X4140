//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jan 22 13:42:17 2018 by ROOT version 6.08/07
// from TTree xTree/Tree of xs
// found on file: /lustre/cms/store/user/adiflori/MuOnia/phiJpsiTriggersBCDEF.root
//////////////////////////////////////////////////////////

#ifndef mumumumu_h
#define mumumumu_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector
#include "TLorentzVector.h"

#include "Math/GenVector/PositionVector3D.h"

#include "DataFormats/VertexReco/interface/Vertex.h"



class mumumumu : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<UInt_t> lumiblock = {fReader, "lumiblock"};
   TTreeReaderValue<UInt_t> run = {fReader, "run"};
   TTreeReaderValue<ULong64_t> event = {fReader, "event"};
   TTreeReaderValue<TLorentzVector> x_p4 = {fReader, "x_p4"};
   TTreeReaderValue<Double_t> xM = {fReader, "xM"};
   TTreeReaderValue<UInt_t> trigger = {fReader, "trigger"};
   TTreeReaderValue<UInt_t> jpsi_i = {fReader, "jpsi_i"};
   TTreeReaderValue<UInt_t> phi_i = {fReader, "phi_i"};
   TTreeReaderValue<TLorentzVector> jpsi_p4 = {fReader, "jpsi_p4"};
   TTreeReaderValue<Double_t> jpsi_M = {fReader, "jpsi_M"};
   TTreeReaderValue<TLorentzVector> muonM_jpsi_p4 = {fReader, "muonM_jpsi_p4"};
   TTreeReaderValue<TLorentzVector> muonP_jpsi_p4 = {fReader, "muonP_jpsi_p4"};
   TTreeReaderValue<Int_t> jpsi_muonM_type = {fReader, "jpsi_muonM_type"};
   TTreeReaderValue<Int_t> jpsi_muonP_type = {fReader, "jpsi_muonP_type"};
   TTreeReaderValue<TLorentzVector> phi_p4 = {fReader, "phi_p4"};
   TTreeReaderValue<Double_t> phi_M = {fReader, "phi_M"};
   TTreeReaderValue<TLorentzVector> muonM_phi_p4 = {fReader, "muonM_phi_p4"};
   TTreeReaderValue<TLorentzVector> muonP_phi_p4 = {fReader, "muonP_phi_p4"};
   TTreeReaderValue<Int_t> phi_muonM_type = {fReader, "phi_muonM_type"};
   TTreeReaderValue<Int_t> phi_muonP_type = {fReader, "phi_muonP_type"};
   TTreeReaderValue<UInt_t> numPrimaryVertices = {fReader, "numPrimaryVertices"};
   TTreeReaderValue<Double_t> dz = {fReader, "dz"};
   TTreeReaderValue<Double_t> dz_jpsi = {fReader, "dzjpsi"};
   TTreeReaderValue<Double_t> dz_phi = {fReader, "dzphi"};
   TTreeReaderValue<Double32_t> fCoordinates_fX = {fReader, "fCoordinates.fX"};
   TTreeReaderValue<Double32_t> fCoordinates_fY = {fReader, "fCoordinates.fY"};
   TTreeReaderValue<Double32_t> fCoordinates_fZ = {fReader, "fCoordinates.fZ"};
   TTreeReaderValue<Float_t> chi2_ = {fReader, "chi2_"};
   TTreeReaderValue<Float_t> ndof_ = {fReader, "ndof_"};
   TTreeReaderValue<Double_t> position__fCoordinates_fX = {fReader, "position_.fCoordinates.fX"};
   TTreeReaderValue<Double_t> position__fCoordinates_fY = {fReader, "position_.fCoordinates.fY"};
   TTreeReaderValue<Double_t> position__fCoordinates_fZ = {fReader, "position_.fCoordinates.fZ"};
   TTreeReaderArray<Float_t> covariance_ = {fReader, "covariance_[10]"};
   TTreeReaderArray<UChar_t> refittedTracks__hitPattern__hitCount = {fReader, "refittedTracks_.hitPattern_.hitCount"};
   TTreeReaderArray<UChar_t> refittedTracks__hitPattern__beginTrackHits = {fReader, "refittedTracks_.hitPattern_.beginTrackHits"};
   TTreeReaderArray<UChar_t> refittedTracks__hitPattern__endTrackHits = {fReader, "refittedTracks_.hitPattern_.endTrackHits"};
   TTreeReaderArray<UChar_t> refittedTracks__hitPattern__beginInner = {fReader, "refittedTracks_.hitPattern_.beginInner"};
   TTreeReaderArray<UChar_t> refittedTracks__hitPattern__endInner = {fReader, "refittedTracks_.hitPattern_.endInner"};
   TTreeReaderArray<UChar_t> refittedTracks__hitPattern__beginOuter = {fReader, "refittedTracks_.hitPattern_.beginOuter"};
   TTreeReaderArray<UChar_t> refittedTracks__hitPattern__endOuter = {fReader, "refittedTracks_.hitPattern_.endOuter"};
   TTreeReaderArray<Double32_t> refittedTracks__vertex__fCoordinates_fX = {fReader, "refittedTracks_.vertex_.fCoordinates.fX"};
   TTreeReaderArray<Double32_t> refittedTracks__vertex__fCoordinates_fY = {fReader, "refittedTracks_.vertex_.fCoordinates.fY"};
   TTreeReaderArray<Double32_t> refittedTracks__vertex__fCoordinates_fZ = {fReader, "refittedTracks_.vertex_.fCoordinates.fZ"};
   TTreeReaderArray<Double32_t> refittedTracks__momentum__fCoordinates_fX = {fReader, "refittedTracks_.momentum_.fCoordinates.fX"};
   TTreeReaderArray<Double32_t> refittedTracks__momentum__fCoordinates_fY = {fReader, "refittedTracks_.momentum_.fCoordinates.fY"};
   TTreeReaderArray<Double32_t> refittedTracks__momentum__fCoordinates_fZ = {fReader, "refittedTracks_.momentum_.fCoordinates.fZ"};
   TTreeReaderArray<edm::RefCoreWithIndex> refittedTracks__extra__product_ = {fReader, "refittedTracks_.extra_.product_"};
   TTreeReaderArray<unsigned char> weights_ = {fReader, "weights_"};
   TTreeReaderValue<Bool_t> validity_ = {fReader, "validity_"};
   TTreeReaderValue<Double_t> time_ = {fReader, "time_"};
   TTreeReaderValue<UInt_t> jpsi_trigger = {fReader, "jpsi_trigger"};
   TTreeReaderValue<UInt_t> phi_trigger = {fReader, "phi_trigger"};
   TTreeReaderValue<UInt_t> jpsi_deltaR = {fReader, "jpsi_deltaR"};
   TTreeReaderValue<UInt_t> phi_deltaR = {fReader, "phi_deltaR"};
   TTreeReaderValue<UInt_t> countTksOfPV = {fReader, "countTksOfPV"};
   TTreeReaderValue<Double_t> vertexWeight = {fReader, "vertexWeight"};
   TTreeReaderValue<Double_t> sumPTPV = {fReader, "sumPTPV"};
   TTreeReaderValue<Double_t> vProb = {fReader, "vProb"};
   TTreeReaderValue<Double_t> vNChi2 = {fReader, "vNChi2"};
   TTreeReaderValue<Double_t> ctauBS = {fReader, "ctauBS"};
   TTreeReaderValue<Double_t> ctauErrBS = {fReader, "ctauErrBS"};
   TTreeReaderValue<Double_t> ctauPV = {fReader, "ctauPV"};
   TTreeReaderValue<Double_t> ctauErrPV = {fReader, "ctauErrPV"};
   TTreeReaderValue<Double_t> ctauPVMuLess = {fReader, "ctauPVMuLess"};
   TTreeReaderValue<Double_t> ctauErrPVMuLess = {fReader, "ctauErrPVMuLess"};
   TTreeReaderValue<Double_t> cosAlpha = {fReader, "cosAlpha"};
   TTreeReaderValue<Double_t> cosAlphaMuLess = {fReader, "cosAlphaMuLess"};
   TTreeReaderValue<Double_t> cosAlphaBS = {fReader, "cosAlphaBS"};
   TTreeReaderValue<Double_t> cosAlpha3D = {fReader, "cosAlpha3D"};
   TTreeReaderValue<Double_t> cosAlphaBS3D = {fReader, "cosAlphaBS3D"};
   TTreeReaderValue<Double_t> l_xy = {fReader, "l_xy"};
   TTreeReaderValue<Double_t> l_xyBS = {fReader, "l_xyBS"};
   TTreeReaderValue<Double_t> l_xyz = {fReader, "l_xyz"};
   TTreeReaderValue<Double_t> l_xyzBS = {fReader, "l_xyzBS"};
   TTreeReaderValue<Double_t> lErr_xy = {fReader, "lErr_xy"};
   TTreeReaderValue<Double_t> lErr_xyBS = {fReader, "lErr_xyBS"};
   TTreeReaderValue<Double_t> lErr_xyz = {fReader, "lErr_xyz"};
   TTreeReaderValue<Double_t> lErr_xyzBS = {fReader, "lErr_xyzBS"};
   TTreeReaderValue<Int_t> x_rank = {fReader, "x_rank"};


   mumumumu(TTree * /*tree*/ =0) { }
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

   ClassDef(mumumumu,0);

};

#endif

#ifdef mumumumu_cxx
void mumumumu::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
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


#endif // #ifdef mumumumu_cxx
~
