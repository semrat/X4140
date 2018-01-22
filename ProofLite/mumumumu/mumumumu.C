#define mumumumu_cxx
// The class definition in Y4140.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("Y4140.C")
// Root > T->Process("Y4140.C","some options")
// Root > T->Process("Y4140.C+")
//

#include "mumumumu.h"
#include <TH2.h>
#include <TStyle.h>

/// SEMRA added
#include <TCanvas.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TSystem.h>
#include <TTree.h>
#include <TBranch.h>
//#include <TCint.h>
#include <TRandom.h>
#include <TMath.h>
#include <TDirectory.h>
#include "TEnv.h"
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TString.h>
#include <TProof.h>
#include <TProofOutputFile.h>
#include "TLorentzVector.h"
#include "TPoint.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <tuple>
#include <map>
#include "TSelectorCint.h"

#include <iostream>
#include <fstream>

//for flags
#include <bitset>

void mumumumu::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

}


void mumumumu::SlaveBegin(TTree * /*tree*/)
{
  //std::cout << "SLAVE BEGIN" << std::endl;

  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

  std::string outputString = "mumumumu_tree.root";
  OutFile = new TProofOutputFile( outputString.data() );
  fOut = OutFile->OpenFile("RECREATE");
  if (!(fOut=OutFile->OpenFile("RECREATE")))
  {
    Warning("SlaveBegin","Problems opening file: %s%s", OutFile->GetDir(), OutFile->GetFileName() );
  }

  ////////////////// Histograms //////////////////
  JPsi_mass = 3.096916; /// pdg mass
  Phi_mass = 1.019455; /// pdg mass
  Phi_mean = 1.019723;
  Phi_sigma = 2.35607e-03;//2.28400e-03;

  outTuple = new TNtuple("outuple","outuple","run:evt:lum:xHlt:xM:phiM:jPsiM:xL:xProb:xCos:xP4:phiTrigger:jpsiTrigger");


}

bool mumumumu::Process(Long64_t entry)
{
  //std::cout << "PROCESS" << std::endl;

  GetEntry(entry); /// ROOT event loop

  // The entry argument specifies which entry in the currently loaded tree is to be processed.
  // It can be passed to either Y4140::GetEntry() or TBranch::GetEntry()
  // to read either all or the required parts of the data. When processing

  Float_t          run_out;
  Float_t          evt_out;
  Float_t          lum_out;
  Float_t         X_mass;
  Float_t         phi_mass;
  Float_t         jpsi_mass;
  Float_t         X_LFly;
  Float_t         X_pt;
  Float_t           X_eta;
  Float_t            X_vtx;
  Float_t            X_cosAlpha;
  Float_t            X_chi2;
  Float_t            X_dZ;
  UInt_t             X_hlt;
  UInt_t             p_hlt;
  UInt_t             j_hlt;

  TLorentzVector X_p4;

  // #Phi
  // 'HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi',
  // 'HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi',
  // 'HLT_Mu20_TkMu0_Phi',
  // 'HLT_Dimuon14_Phi_Barrel_Seagulls',
  // 'HLT_Mu25_TkMu0_Phi',
  // 'HLT_Dimuon24_Phi_noCorrL1',
  // #JPsi
  // 'HLT_DoubleMu4_JpsiTrkTrk_Displaced',
  // 'HLT_DoubleMu4_JpsiTrk_Displaced',
  // 'HLT_DoubleMu4_Jpsi_Displaced',
  // 'HLT_DoubleMu4_3_Jpsi_Displaced',
  // 'HLT_Dimuon20_Jpsi_Barrel_Seagulls',
  // 'HLT_Dimuon25_Jpsi',
  // 'HLT_Dimuon0_Jpsi'

  run_out = runNum;
  evt_out = evtNum;
  lum_out = lumiNum;

  X_mass = xM;

  // std::bitset<13> bitsXHLT(0x0);
  X_hlt = trigFlag;

  jpsi_mass = jPsiM;
  phi_mass = phiM;
  X_cosAlpha = cosAlpha;
  X_LFly = lxy / lxyErr;
  X_chi2 = vNChi2;
  X_dZ = dz;
  X_p4 = xP4;

  X_pt = X_p4.Pt();
  X_eta = X_p4.Eta();

  p_hlt = phiTrigger;
  j_hlt = jpsiTrigger;

  outTuple->Fill(run_out,evt_out,lum_out,X_mass,X_hlt,jpsi_mass,phi_mass,X_cosAlpha,X_LFly,X_chi2,X_dZ,X_pt,X_eta,p_hlt,j_hlt);

  return kTRUE;
}


void mumumumu::SlaveTerminate()
{

  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
  TDirectory *savedir = gDirectory;
  if (fOut)
  {
    fOut->cd();
    gStyle->SetOptStat(111111) ;


    outTuple->Write();
    OutFile->Print();
    fOutput->Add(OutFile);
    gDirectory = savedir;
    fOut->Close();

  }


}


void mumumumu::Terminate()
{


}
