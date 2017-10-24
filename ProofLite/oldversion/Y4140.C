#define Y4140_cxx
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

#include "Y4140.h"
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
#include "TSelectorCint.h"


void Y4140::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

}


void Y4140::SlaveBegin(TTree * /*tree*/)
{
  //std::cout << "SLAVE BEGIN" << std::endl;

  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

  /// SEMRA added
  OutFile = new TProofOutputFile( "Y4140.root" );
  fOut = OutFile->OpenFile("RECREATE");
  if (!(fOut=OutFile->OpenFile("RECREATE")))
  {
    Warning("SlaveBegin","Problems opening file: %s%s", OutFile->GetDir(), OutFile->GetFileName() );
  }

  ////////////////// Histograms //////////////////
  JPsi_mass = 3.0969; /// pdg mass
  Phi_mass = 1.020; /// pdg mass

  /// refit muons & kaons
  Hist_Muon1_pt         = new TH1F ("Muon1_pT", "Muon1_pT; pT(#mu^{-}) [GeV];Entries",100, 3, 20);
  Hist_Muon1_Eta        = new TH1F ("Muon1_Eta","Muon1_Eta; Eta(#mu^{-});Entries",100, -5, 5);
  Hist_Muon1_Phi        = new TH1F ("Muon1_Phi","Muon1_Phi; Phi(#mu^{-});Entries ",100, -5, 5);
  Hist_Muon2_pt         = new TH1F ("Muon2_pT", "Muon2_pT; pT(#mu^{+}) [GeV];Entries",100, 3, 20);
  Hist_Muon2_Eta        = new TH1F ("Muon2_Eta","Muon2_Eta; Eta(#mu^{+});Entries",100, -5, 5);
  Hist_Muon2_Phi        = new TH1F ("Muon2_Phi","Muon2_Phi; Phi(#mu^{+});Entries",100, -5, 5);
  Hist_Muon1PtvsMuon2Pt = new TH2F ("Muon1PtvsMuon2Pt","Muon1PtvsMuon2Pt; pT(#mu^{-}) [GeV];pT(#mu^{+}) [GeV]",100, 3, 20,100, 3, 20); /// scatter plot
  Hist_Kaon1_pt         = new TH1F ("Kaon1_pT","Kaon1_pT; pT(K^{-}) [GeV];Entries",100, 0, 10);
  Hist_Kaon1_Eta        = new TH1F ("Kaon1_Eta","Kaon1_Eta; Eta(K^{-});Entries",100, -5, 5);
  Hist_Kaon1_Phi        = new TH1F ("Kaon1_Phi","Kaon1_Phi; Phi(K^{-});Entries",100, -5, 5);
  Hist_Kaon2_pt         = new TH1F ("Kaon2_pT","Kaon2_pT; pT(K^{+}) [GeV];Entries",100, 0, 10);
  Hist_Kaon2_Eta        = new TH1F ("Kaon2_Eta","Kaon2_Eta; Eta(K^{+});Entries",100, -5, 5);
  Hist_Kaon2_Phi        = new TH1F ("Kaon2_Phi","Kaon2_Phi; Phi(K^{+});Entries",100, -5, 5);
  Hist_Kaon1PtvsKaon2Pt = new TH2F ("Kaon1PtvsKaon2Pt","Kaon1PtvsKaon2Pt; pT(K^{-}) [GeV];pT(K^{+}) [GeV];",100, 0, 10,100, 0, 10); /// scatter plot

  /// muon pairs from original muons
  Hist_mumu_mass = new TH1F ("mumu_mass","mumu_mass;M(#mu^{-}#mu^{+}) [GeV];Entries",100, 2.8, 3.4);
  Hist_mumu_pt   = new TH1F ("mumu_pT", "mumu_pT; pT(#mu^{-}#mu^{+}) [GeV];Entries",100, 0, 50);
  Hist_mumu_Eta  = new TH1F ("mumu_Eta","mumu_Eta; Eta(#mu^{-}#mu^{+});Entries",100, -5, 5);
  Hist_mumu_Phi  = new TH1F ("mumu_Phi","mumu_Phi; Phi(#mu^{-}#mu^{+});Entries ",100, -5, 5);
  Hist_kk_mass   = new TH1F ("kk_mass","kk_mass;M(K^{-}K^{+}) [GeV];Entries",100, 0.97, 1.07);
  Hist_kk_pt     = new TH1F ("kk_pT", "kk_pT; pT(K^{-}K^{+}) [GeV];Entries",100, 0, 50);
  Hist_kk_Eta    = new TH1F ("kk_Eta","kk_Eta; Eta(K^{-}K^{+});Entries",100, -5, 5);
  Hist_kk_Phi    = new TH1F ("kk_Phi","kk_Phi; Phi(K^{-}K^{+});Entries ",100, -5, 5);

  /// refit muons from JPsi & refit kaons from Phi
  Hist_JPsi_mass         = new TH1F ("JPsi_mass","JPsi_mass;M(#mu^{-}#mu^{+}) [GeV];Entries",100,2.8, 3.4);
  Hist_JPsi_pt           = new TH1F ("JPsi_pT", "JPsi_pT; pT(#mu^{-}#mu^{+}) [GeV];Entries",100, 6, 50);
  Hist_JPsi_Eta          = new TH1F ("JPsi_Eta","JPsi_Eta; Eta(#mu^{-}#mu^{+});Entries",100, -5, 5);
  Hist_JPsi_Phi          = new TH1F ("JPsi_Phi","JPsi_Phi; Phi(#mu^{-}#mu^{+});Entries ",100, -5, 5);
  Hist_JPsi_PtvsEta      = new TH2F ("JPsi_PtvsEta","JPsi_PtvsEta;Eta(#mu^{-}#mu^{+});pT(#mu^{-}#mu^{+}) [GeV]",100, -2.6, 2.6, 100, 5, 35);
  Hist_JPsi_mass_HLT     = new TH1F ("JPsi_mass_HLT","JPsi_mass_HLT;M(#mu^{-}#mu^{+}) [GeV];Entries",100,2.8, 3.4);
  Hist_JPsi_mass_SoftMu  = new TH1F ("JPsi_mass_SoftMu","JPsi_mass_SoftMu;M(#mu^{-}#mu^{+}) [GeV];Entries",100,2.8, 3.4);
  Hist_JPsi_PtvsEta_JPsi = new TH2F ("JPsi_PtvsEta_JPsi","JPsi_PtvsEta_JPsi;Eta(#mu^{-}#mu^{+});pT(#mu^{-}#mu^{+}) [GeV]",100, -2.6, 2.6, 100, 5, 35);
  Hist_Phi_mass          = new TH1F ("Phi_mass","Phi_mass;M(K^{-}K^{+}) [GeV];Entries",100, 0.97, 1.07);
  Hist_Phi_pt            = new TH1F ("Phi_pT", "Phi_pT; pT(K^{-}K^{+}) [GeV];Entries",100, 0, 50);
  Hist_Phi_Eta           = new TH1F ("Phi_Eta","Phi_Eta; Eta(K^{-}K^{+});Entries",100, -5, 5);
  Hist_Phi_Phi           = new TH1F ("Phi_Phi","Phi_Phi; Phi(K^{-}K^{+});Entries ",100, -5, 5);
  Hist_Phi_PtvsEta       = new TH2F ("Phi_PtvsEta","Phi_PtvsEta;Eta(K^{-}K^{+});pT(K^{-}K^{+}) [GeV]",100, -2.6, 2.6, 100, 5, 35);
  Hist_Phi_PtvsEta_JPsi  = new TH2F ("Phi_PtvsEta_JPsi","Phi_PtvsEta_JPsi;Eta(K^{-}K^{+});pT(K^{-}K^{+}) [GeV]",100, -2.6, 2.6, 100, 5, 35);

  /// Bs0 & Some Cuts
  Hist_Bs0_Y_mass_fromNTuple = new TH1F ("Bs0_Y_mass_fromNTuple","Bs0_Y_mass_fromNTuple;M(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV];Entries",100,4., 6.);
  Hist_Bs0_Y_mass            = new TH1F ("Bs0_Y_mass","Bs0_Y_mass;M(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV];Entries",100,4., 6.);
  Hist_Bs0_Y_pt              = new TH1F ("Bs0_Y_pT","Bs0_Y_pT; pT(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV];Entries",100, 6, 50);
  Hist_Bs0_Y_Eta             = new TH1F ("Bs0_Y_Eta","Bs0_Y_Eta; Eta(#mu^{-}#mu^{+}K^{-}K^{+});Entries",100, -5, 5);
  Hist_Bs0_Y_Phi             = new TH1F ("Bs0_Y_Phi","Bs0_Y_Phi; Phi (#mu^{-}#mu^{+}K^{-}K^{+});Entries",100, -5, 5);

  /// Y4140 & Some Cuts (HLT, SoftMu, JPsi, Bs0, Phi, Prompt, Middle, Non-Prompt)
  Hist_Y4140_mass           = new TH1F ("Y4140_mass","Y4140_mass;M(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV];Entries",100, 4.0, 4.9);
  Hist_Y4140_pt             = new TH1F ("Y4140_pT","Y4140_pT; pT(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV];Entries",100, 8, 50);
  Hist_Y4140_Eta            = new TH1F ("Y4140_Eta","Y4140_Eta; Eta(#mu^{-}#mu^{+}K^{-}K^{+});Entries",100, -5, 5);
  Hist_Y4140_Phi            = new TH1F ("Y4140_Phi","Y4140_Phi; Phi (#mu^{-}#mu^{+}K^{-}K^{+});Entries",100, -5, 5);
  Hist_Y4140_PtvsEta        = new TH2F ("Y4140_PtvsEta","Y4140_PtvsEta;Eta(#mu^{-}#mu^{+}K^{-}K^{+});pT(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV]",100, -2.6, 2.6, 100, 5, 35);
  Hist_Y4140_mass_HLT       = new TH1F ("Y4140_mass_HLT","Y4140_mass_HLT;M(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV];Entries",100, 4.0, 4.9);
  Hist_Y4140_mass_SoftMu    = new TH1F ("Y4140_mass_SoftMu","Y4140_mass_SoftMu;M(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV];Entries",100,4.0, 4.9);
  Hist_Y4140_mass_JPsi      = new TH1F ("Y4140_mass_JPsi","Y4140_mass_JPsi;M(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV];Entries",100,4.0, 4.9);
  Hist_Y4140_PtvsEta_JPsi   = new TH2F ("Y4140_PtvsEta_JPsi","Y4140_PtvsEta_JPsi;Eta(#mu^{-}#mu^{+}K^{-}K^{+});pT(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV]",100, -2.6, 2.6, 100, 5, 35);
  Hist_Y4140_mass_Prompt    = new TH1F ("Y4140_mass_Prompt","Y4140_mass_Prompt;M(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV];Entries",100,4.0, 4.9);
  Hist_Y4140_mass_Mixed     = new TH1F ("Y4140_mass_Mixed","Y4140_mass_Mixed;M(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV];Entries",100,4.0, 4.9);
  Hist_Y4140_mass_NonPrompt = new TH1F ("Y4140_mass_NonPrompt","Y4140_mass_NonPrompt;M(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV];Entries",100,4.0, 4.9);
  Hist_Y4140_mass_Track     = new TH1F ("Y4140_mass_Track","Y4140_mass_Track;M(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV];Entries",100,4.0, 4.9);
  Hist_Y4140_mass_TrackP    = new TH1F ("Y4140_mass_TrackP","Y4140_mass_TrackP;M(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV];Entries",100,4.0, 4.9);
  Hist_Y4140_mass_TrackM    = new TH1F ("Y4140_mass_TrackM","Y4140_mass_TrackM;M(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV];Entries",100,4.0, 4.9);
  Hist_Y4140_mass_TrackNP   = new TH1F ("Y4140_mass_TrackNP","Y4140_mass_TrackNP;M(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV];Entries",100,4.0, 4.9);
  Hist_Y4140_mass_Phi       = new TH1F ("Y4140_mass_Phi","Y4140_mass_Phi;M(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV];Entries",100,4.0, 4.9);
  Hist_Y4140_mass_PhiP      = new TH1F ("Y4140_mass_PhiP","Y4140_mass_PhiP;M(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV];Entries",100,4.0, 4.9);
  Hist_Y4140_mass_PhiM      = new TH1F ("Y4140_mass_PhiM","Y4140_mass_PhiM;M(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV];Entries",100,4.0, 4.9);
  Hist_Y4140_mass_PhiNP     = new TH1F ("Y4140_mass_PhiNP","Y4140_mass_PhiNP;M(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV];Entries",100,4.0, 4.9);
  Hist_Y4140_mass_Bs0       = new TH1F ("Y4140_mass_Bs0","Y4140_mass_Bs0;M(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV];Entries",100,4.0, 4.9);

  Hist_Y4140_mass_Bs0P      = new TH1F ("Y4140_mass_Bs0P","Y4140_mass_Bs0P;M(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV];Entries",100,4.0, 4.9);
  Hist_SC_JPsi_mass_Bs0P    = new TH1F ("SC_JPsi_mass_Bs0P","SC_JPsi_mass_Bs0P;M(#mu^{-}#mu^{+}) [GeV];Entries",100,2.8, 3.4);
  Hist_SC_Phi_mass_Bs0P     = new TH1F ("SC_Phi_mass_Bs0P","SC_Phi_mass_Bs0P;M(K^{-}K^{+}) [GeV];Entries",100, 0.97, 1.07);

  Hist_Y4140_mass_Bs0M      = new TH1F ("Y4140_mass_Bs0M","Y4140_mass_Bs0M;M(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV];Entries",100,4.0, 4.9);
  Hist_SC_JPsi_mass_Bs0M    = new TH1F ("SC_JPsi_mass_Bs0M","SC_JPsi_mass_Bs0M;M(#mu^{-}#mu^{+}) [GeV];Entries",100,2.8, 3.4);
  Hist_SC_Phi_mass_Bs0M     = new TH1F ("SC_Phi_mass_Bs0M","SC_Phi_mass_Bs0M;M(K^{-}K^{+}) [GeV];Entries",100, 0.97, 1.07);

  Hist_Y4140_mass_Bs0NP     = new TH1F ("Y4140_mass_Bs0NP","Y4140_mass_Bs0NP;M(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV];Entries",100,4.0, 4.9);
  Hist_SC_JPsi_mass_Bs0NP   = new TH1F ("SC_JPsi_mass_Bs0NP","SC_JPsi_mass_Bs0NP;M(#mu^{-}#mu^{+}) [GeV];Entries",100,2.8, 3.4);
  Hist_SC_Phi_mass_Bs0NP    = new TH1F ("SC_Phi_mass_Bs0NP","SC_Phi_mass_Bs0NP;M(K^{-}K^{+}) [GeV];Entries",100, 0.97, 1.07);

  /// Bs0 & Some Cuts (HLT, SoftMu, JPsi, Bs0, Phi, Prompt, Middle, Non-Prompt)
  Hist_Bs0_mass           = new TH1F ("Bs0_mass","Bs0_mass;M(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV];Entries",100,5.15, 5.55);
  Hist_Bs0_pt             = new TH1F ("Bs0_pT","Bs0_pT; pT(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV];Entries",100, 8, 50);
  Hist_Bs0_Eta            = new TH1F ("Bs0_Eta","Bs0_Eta; Eta(#mu^{-}#mu^{+}K^{-}K^{+});Entries",100, -5, 5);
  Hist_Bs0_Phi            = new TH1F ("Bs0_Phi","Bs0_Phi; Phi (#mu^{-}#mu^{+}K^{-}K^{+});Entries",100, -5, 5);
  Hist_Bs0_PtvsEta        = new TH2F ("Bs0_PtvsEta","Bs0_PtvsEta;Eta(#mu^{-}#mu^{+}K^{-}K^{+});pT(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV]",100, -2.6, 2.6, 100, 5, 35);
  Hist_Bs0_mass_HLT       = new TH1F ("Bs0_mass_HLT","Bs0_mass_HLT;M(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV];Entries",100,5.15, 5.55);
  Hist_Bs0_mass_SoftMu    = new TH1F ("Bs0_mass_SoftMu","Bs0_mass_SoftMu;M(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV];Entries",100,5.15, 5.55);
  Hist_Bs0_mass_JPsi      = new TH1F ("Bs0_mass_JPsi","Bs0_mass_JPsi;M(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV];Entries",100,5.15, 5.55);
  Hist_Bs0_PtvsEta_JPsi   = new TH2F ("Bs0_PtvsEta_JPsi","Bs0_PtvsEta_JPsi;Eta(#mu^{-}#mu^{+}K^{-}K^{+});pT(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV]",100, -2.6, 2.6, 100, 5, 35);
  Hist_Bs0_mass_Prompt    = new TH1F ("Bs0_mass_Prompt","Bs0_mass_Prompt;M(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV];Entries",100,5.15, 5.55);
  Hist_Bs0_mass_Mixed     = new TH1F ("Bs0_mass_Mixed","Bs0_mass_Mixed;M(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV];Entries",100,5.15, 5.55);
  Hist_Bs0_mass_NonPrompt = new TH1F ("Bs0_mass_NonPrompt","Bs0_mass_NonPrompt;M(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV];Entries",100,5.15, 5.55);
  Hist_Bs0_mass_Track     = new TH1F ("Bs0_mass_Track","Bs0_mass_Track;M(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV];Entries",100,5.15, 5.55);
  Hist_Bs0_mass_TrackP    = new TH1F ("Bs0_mass_TrackP","Bs0_mass_TrackP;M(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV];Entries",100,5.15, 5.55);
  Hist_Bs0_mass_TrackM    = new TH1F ("Bs0_mass_TrackM","Bs0_mass_TrackM;M(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV];Entries",100,5.15, 5.55);
  Hist_Bs0_mass_TrackNP   = new TH1F ("Bs0_mass_TrackNP","Bs0_mass_TrackNP;M(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV];Entries",100,5.15, 5.55);
  Hist_Bs0_mass_Phi       = new TH1F ("Bs0_mass_Phi","Bs0_mass_Phi;M(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV];Entries",100,5.15, 5.55);
  Hist_Bs0_mass_PhiP      = new TH1F ("Bs0_mass_PhiP","Bs0_mass_PhiP;M(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV];Entries",100,5.15, 5.55);
  Hist_Bs0_mass_PhiM      = new TH1F ("Bs0_mass_PhiM","Bs0_mass_PhiM;M(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV];Entries",100,5.15, 5.55);
  Hist_Bs0_mass_PhiNP     = new TH1F ("Bs0_mass_PhiNP","Bs0_mass_PhiNP;M(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV];Entries",100,5.15, 5.55);
  Hist_Bs0_mass_Bs0       = new TH1F ("Bs0_mass_Bs0","Bs0_mass_Bs0;M(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV];Entries",101,5.15, 5.55);

  Hist_Bs0_mass_Bs0P      = new TH1F ("Bs0_mass_Bs0P","Bs0_mass_Bs0P;M(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV];Entries",80,5.15, 5.55);
  Hist_CC_JPsi_mass_Bs0P  = new TH1F ("CC_JPsi_mass_Bs0P","CC_JPsi_mass_Bs0P;M(#mu^{-}#mu^{+}) [GeV];Entries",100,2.8, 3.4);
  Hist_CC_Phi_mass_Bs0P   = new TH1F ("CC_Phi_mass_Bs0P","CC_Phi_mass_Bs0P;M(K^{-}K^{+}) [GeV];Entries",100, 0.97, 1.07);

  Hist_Bs0_mass_Bs0M      = new TH1F ("Bs0_mass_Bs0M","Bs0_mass_Bs0M;M(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV];Entries",80,5.15, 5.55);
  Hist_CC_JPsi_mass_Bs0M  = new TH1F ("CC_JPsi_mass_Bs0M","CC_JPsi_mass_Bs0M;M(#mu^{-}#mu^{+}) [GeV];Entries",100,2.8, 3.4);
  Hist_CC_Phi_mass_Bs0M   = new TH1F ("CC_Phi_mass_Bs0M","CC_Phi_mass_Bs0M;M(K^{-}K^{+}) [GeV];Entries",100, 0.97, 1.07);

  Hist_Bs0_mass_Bs0NP     = new TH1F ("Bs0_mass_Bs0NP","Bs0_mass_Bs0NP;M(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV];Entries",80,5.15, 5.55);
  Hist_CC_JPsi_mass_Bs0NP = new TH1F ("CC_JPsi_mass_Bs0NP","CC_JPsi_mass_Bs0NP;M(#mu^{-}#mu^{+}) [GeV];Entries",100,2.8, 3.4);
  Hist_CC_Phi_mass_Bs0NP  = new TH1F ("CC_Phi_mass_Bs0NP","CC_Phi_mass_Bs0NP;M(K^{-}K^{+}) [GeV];Entries",100, 0.97, 1.07);


  Hist_CTau_CTauE_PV = new TH1F ("CTau_CTauE_PV","CTau_CTauE_PV;Entries", 100, -2, 16);
  Hist_Lxy_LxyE_PV   = new TH1F ("Lxy_LxyE_PV","Lxy_LxyE_PV;Entries", 100, -2, 16);
  Hist_CTauVsLxy_PV  = new TH2F ("CTauVsLxy_PV","CTauVsLxy_PV;LxyPV/LxyPVE;CTauPV/CTauPVE",100, -2, 16, 100, -2, 16);

  /*
  Hist_CTau_CTauE_BS = new TH1F ("CTau_CTauE_BS","CTau_CTauE_BS;Entries", 100, -2, 16);
  Hist_Lxy_LxyE_BS   = new TH1F ("Lxy_LxyE_BS","Lxy_LxyE_BS;Entries", 100, -2, 16);
  Hist_CTauVsLxy_BS  = new TH2F ("CTauVsLxy_BS","CTauVsLxy_BS;LxyBS/LxyBSE;CTauBS/CTauBSE",100, -2, 16, 100, -2, 16);

  Hist_CTau_CTauE_XLessPV = new TH1F ("CTau_CTauE_XLessPV","CTau_CTauE_XLessPV;Entries", 100, -2, 16);
  Hist_Lxy_LxyE_XLessPV   = new TH1F ("Lxy_LxyE_XLessPV","Lxy_LxyE_XLessPV;Entries", 100, -2, 16);
  Hist_CTauVsLxy_XLessPV  = new TH2F ("CTauVsLxy_XLessPV","CTauVsLxy_XLessPV;LxyXLessPV/LxyXLessPVE;CTauXLessPV/CTauXLessPVE",100, -2, 16, 100, -2, 16);

  Hist_CTau_CTauE_PVCosAlpha = new TH1F ("CTau_CTauE_PVCosAlpha","CTau_CTauE_PVCosAlpha;Entries", 100, -2, 16);
  Hist_Lxy_LxyE_PVCosAlpha   = new TH1F ("Lxy_LxyE_PVCosAlpha","Lxy_LxyE_PVCosAlpha;Entries", 100, -2, 16);
  Hist_CTauVsLxy_PVCosAlpha  = new TH2F ("CTauVsLxy_PVCosAlpha","CTauVsLxy_PVCosAlpha;LxyPVCosAlpha/LxyPVCosAlphaE;CTauPVCosAlpha/CTauPVCosAlphaE",100, -2, 16, 100, -2, 16);

  Hist_CTau_CTauE_PVCosAlpha3D = new TH1F ("CTau_CTauE_PVCosAlpha3D","CTau_CTauE_PVCosAlpha3D;Entries", 100, -2, 16);
  Hist_Lxy_LxyE_PVCosAlpha3D   = new TH1F ("Lxy_LxyE_PVCosAlpha3D","Lxy_LxyE_PVCosAlpha3D;Entries", 100, -2, 16);
  Hist_CTauVsLxy_PVCosAlpha3D  = new TH2F ("CTauVsLxy_PVCosAlpha3D","CTauVsLxy_PVCosAlpha3D;LxyPVCosAlpha3D/LxyPVCosAlpha3DE;CTauPVCosAlpha3D/CTauPVCosAlpha3DE",100, -2, 16, 100, -2, 16);

  Hist_CTau_CTauE_XLessPVCosAlpha = new TH1F ("CTau_CTauE_XLessPVCosAlpha","CTau_CTauE_XLessPVCosAlpha;Entries", 100, -2, 16);
  Hist_Lxy_LxyE_XLessPVCosAlpha   = new TH1F ("Lxy_LxyE_XLessPVCosAlpha","Lxy_LxyE_XLessPVCosAlpha;Entries", 100, -2, 16);
  Hist_CTauVsLxy_XLessPVCosAlpha  = new TH2F ("CTauVsLxy_XLessPVCosAlpha","CTauVsLxy_XLessPVCosAlpha;LxyXLessPVCosAlpha/LxyXLessPVCosAlphaE;CTauXLessPVCosAlpha/CTauXLessPVCosAlphaE",100, -2, 16, 100, -2, 16);

  Hist_CTau_CTauE_XLessPVCosAlpha3D = new TH1F ("CTau_CTauE_XLessPVCosAlpha3D","CTau_CTauE_XLessPVCosAlpha3D;Entries", 100, -2, 16);
  Hist_Lxy_LxyE_XLessPVCosAlpha3D   = new TH1F ("Lxy_LxyE_XLessPVCosAlpha3D","Lxy_LxyE_XLessPVCosAlpha3D;Entries", 100, -2, 16);
  Hist_CTauVsLxy_XLessPVCosAlpha3D  = new TH2F ("CTauVsLxy_XLessPVCosAlpha3D","CTauVsLxy_XLessPVCosAlpha;LxyXLessPVCosAlpha3D/LxyXLessPVCosAlpha3DE;CTauXLessPVCosAlpha3D/CTauXLessPVCosAlpha3DE",100, -2, 16, 100, -2, 16);
  */


  Xmultiplicity = new TH1F ("Xmultiplicity","Xmultiplicity",11, -0.5, 10.5); /// for counting how many X
  nIndex = new TH1F ("nIndex","nIndex",11, -0.5, 10.5);
  psiIndex = new TH1F ("psiIndex","psiIndex",11, -0.5, 10.5);
  psiIndexNx = new TH1F ("psiIndexNx","psiIndexNx",11, -0.5, 10.5);
}

Bool_t Y4140::Process(Long64_t entry)
{
  //std::cout << "PROCESS" << std::endl;

  GetEntry(entry); /// ROOT event loop

  // The entry argument specifies which entry in the currently loaded tree is to be processed.
  // It can be passed to either Y4140::GetEntry() or TBranch::GetEntry()
  // to read either all or the required parts of the data. When processing

  Double_t muon_mass = 0.1056;
  Double_t kaonCh_mass = 0.4936;
  Double_t Bs0_Low_Mass = 5.1;
  Double_t Bs0_High_Mass = 5.6;
  Double_t Y_High_Mass = 4.35;

  ////////////////// Bs0 & X(4140) Loop //////////////////
  Bool_t HLT_4_v9 = false, HLT_4_v10 = false, HLT_4_v11 = false, HLT_4_v12 = false ;
  Bool_t HLT_8_v3 = false, HLT_8_v4 = false, HLT_8_v5 = false, HLT_8_v6 = false, HLT_8_v7 = false ;
  Bool_t HLT_4_vAll = false;
  Bool_t HLT_8_vAll = false;

  for (Int_t i = 0; i != abs((int)(TrigRes->size())); ++i)
  {
    if ( TrigNames->at(i).find("HLT_DoubleMu4_Jpsi_Displaced_v9") != string::npos  &&  TrigRes->at(i) == 1) HLT_4_v9 = true;
    if ( TrigNames->at(i).find("HLT_DoubleMu4_Jpsi_Displaced_v10") != string::npos  &&  TrigRes->at(i) == 1) HLT_4_v10 = true;
    if ( TrigNames->at(i).find("HLT_DoubleMu4_Jpsi_Displaced_v11") != string::npos  &&  TrigRes->at(i) == 1) HLT_4_v11 = true;
    if ( TrigNames->at(i).find("HLT_DoubleMu4_Jpsi_Displaced_v12") != string::npos  &&  TrigRes->at(i) == 1) HLT_4_v12 = true;

    if ( TrigNames->at(i).find("HLT_Dimuon8_Jpsi_v3") != string::npos  &&  TrigRes->at(i) == 1) HLT_8_v3 = true;
    if ( TrigNames->at(i).find("HLT_Dimuon8_Jpsi_v4") != string::npos  &&  TrigRes->at(i) == 1) HLT_8_v4 = true;
    if ( TrigNames->at(i).find("HLT_Dimuon8_Jpsi_v5") != string::npos  &&  TrigRes->at(i) == 1) HLT_8_v5 = true;
    if ( TrigNames->at(i).find("HLT_Dimuon8_Jpsi_v6") != string::npos  &&  TrigRes->at(i) == 1) HLT_8_v6 = true;
    if ( TrigNames->at(i).find("HLT_Dimuon8_Jpsi_v7") != string::npos  &&  TrigRes->at(i) == 1) HLT_8_v7 = true;
  }

  if (HLT_4_v9 || HLT_4_v10 || HLT_4_v11 || HLT_4_v12) HLT_4_vAll = true;
  if (HLT_8_v3 || HLT_8_v4 || HLT_8_v5 || HLT_8_v6 || HLT_8_v7) HLT_8_vAll = true;

  int muonQual[] = {1,3,4,12};



  ////////////////// for loop //////////////////
  Xmultiplicity->Fill(nX); // added by me for counting how many X
  nIndex->Fill(sizeof(*XMuMuIdx)/sizeof(Int_t));
  std::cout<<nX<<std::endl;
  std::cout<<sizeof(*XMuMuIdx)/sizeof(Int_t)<<std::endl;
  for (unsigned int myXIdx=0; myXIdx<sizeof(*XMuMuIdx)/sizeof(Int_t); myXIdx++) // X
  {
    psiIndex->Fill((*XMuMuIdx)[myXIdx]);

  }
  //for (unsigned int myXIdx=0; myXIdx<nX; myXIdx++) // X
  // if(sizeof(*XMuMuIdx)/sizeof((*XMuMuIdx)[0])>=nX)
  //   return kTRUE;


  //for (unsigned int myXIdx=0; myXIdx<nX; myXIdx++)
  int nnn = XMuMuIdx->size();
  for (unsigned int myXIdx=0; myXIdx<nnn; myXIdx++)
  {

    Int_t myJPsiIdx = (*XMuMuIdx)[myXIdx];
    psiIndexNx->Fill(myJPsiIdx);
    int mymupIdx = (*mu1Idx)[myJPsiIdx] ; // define for original muon1
    int mymumIdx = (*mu2Idx)[myJPsiIdx] ; // define for original muon2
    int kaon1Idx = (*ka1Idx)[myXIdx] ; // define for original kaon1
    int kaon2Idx = (*ka2Idx)[myXIdx] ; // define for original kaon2
    //int kaon1Idx = (*XKaon1Idx)[myXIdx] ;
    //int kaon2Idx = (*XKaon2Idx)[myXIdx] ;

    Float_t mu1_E = 0., mu2_E = 0., K1_E = 0., K2_E = 0.;


    //////////// original muons & kaons ////////////
    TLorentzVector mu1, mu2 ;
    mu1_E = sqrt( pow((*muPx)[mymupIdx], 2) + pow((*muPy)[mymupIdx], 2) + pow((*muPz)[mymupIdx], 2) + pow(muon_mass, 2) ) ;
    mu2_E = sqrt( pow((*muPx)[mymumIdx], 2) + pow((*muPy)[mymumIdx], 2) + pow((*muPz)[mymumIdx], 2) + pow(muon_mass, 2) ) ;
    mu1.SetPxPyPzE( (*muPx)[mymupIdx], (*muPy)[mymupIdx], (*muPz)[mymupIdx], mu1_E) ;
    mu2.SetPxPyPzE( (*muPx)[mymumIdx], (*muPy)[mymumIdx], (*muPz)[mymumIdx], mu2_E) ;

    TLorentzVector kaon1;
    K1_E = sqrt( pow((*trackPx)[kaon1Idx],2)+pow((*trackPy)[kaon1Idx],2)+pow((*trackPz)[kaon1Idx],2)+pow(kaonCh_mass,2)) ;
    kaon1.SetPxPyPzE((*trackPx)[kaon1Idx],(*trackPy)[kaon1Idx],(*trackPz)[kaon1Idx], K1_E) ;
    TLorentzVector kaon2;
    K2_E = sqrt( pow((*trackPx)[kaon2Idx],2)+pow((*trackPy)[kaon2Idx],2)+pow((*trackPz)[kaon2Idx],2)+pow(kaonCh_mass,2)) ;
    kaon2.SetPxPyPzE((*trackPx)[kaon2Idx],(*trackPy)[kaon2Idx],(*trackPz)[kaon2Idx], K2_E) ;

    /// crate original muon & kaon pairs
    TLorentzVector mu1mu2;
    mu1mu2 = mu1 + mu2;
    TLorentzVector ka1ka2;
    ka1ka2 = kaon1 + kaon2;
    TLorentzVector mu1mu2ka1ka2;
    mu1mu2ka1ka2 = mu1mu2 + ka1ka2;


    //////////// refit muons from Jpsi & refit kaons from Phi ////////////
    double Kaon1Energy, Kaon2Energy;
    TLorentzVector mymuon1;
    mymuon1.SetPxPyPzE((*Muon1Px_MuMuKK)[myXIdx],(*Muon1Py_MuMuKK)[myXIdx],(*Muon1Pz_MuMuKK)[myXIdx],(*Muon1E_MuMuKK)[myXIdx]);
    TLorentzVector mymuon2;
    mymuon2.SetPxPyPzE((*Muon2Px_MuMuKK)[myXIdx],(*Muon2Py_MuMuKK)[myXIdx],(*Muon2Pz_MuMuKK)[myXIdx],(*Muon2E_MuMuKK)[myXIdx]);

    TLorentzVector mykaon1;
    Kaon1Energy=sqrt(pow((*Kaon1Px_MuMuKK)[myXIdx],2)+pow((*Kaon1Py_MuMuKK)[myXIdx],2)+pow((*Kaon1Pz_MuMuKK)[myXIdx],2)+pow(kaonCh_mass,2));
    mykaon1.SetPxPyPzE((*Kaon1Px_MuMuKK)[myXIdx],(*Kaon1Py_MuMuKK)[myXIdx],(*Kaon1Pz_MuMuKK)[myXIdx],Kaon1Energy);
    TLorentzVector mykaon2;
    Kaon2Energy=sqrt(pow((*Kaon2Px_MuMuKK)[myXIdx],2)+pow((*Kaon2Py_MuMuKK)[myXIdx],2)+pow((*Kaon2Pz_MuMuKK)[myXIdx],2)+pow(kaonCh_mass,2));
    mykaon2.SetPxPyPzE((*Kaon2Px_MuMuKK)[myXIdx],(*Kaon2Py_MuMuKK)[myXIdx],(*Kaon2Pz_MuMuKK)[myXIdx],Kaon2Energy);

    /// crate refit muon & kaon pairs
    TLorentzVector JPsi;
    JPsi = mymuon1 + mymuon2;
    TLorentzVector Phi;
    Phi = mykaon1 + mykaon2;
    TLorentzVector Bs0_Y;
    Bs0_Y = JPsi + Phi;



    ////////////////// Fill Histograms //////////////////
    /// refit muons & kaons
    Hist_Muon1_pt->Fill(mymuon1.Pt());
    Hist_Muon1_Eta->Fill(mymuon1.Eta());
    Hist_Muon1_Phi->Fill(mymuon1.Phi());
    Hist_Muon2_pt->Fill(mymuon2.Pt());
    Hist_Muon2_Eta->Fill(mymuon2.Eta());
    Hist_Muon2_Phi->Fill(mymuon2.Phi());
    Hist_Muon1PtvsMuon2Pt->Fill(mymuon1.Pt(),mymuon2.Pt());
    Hist_Kaon1_pt->Fill(mykaon1.Pt());
    Hist_Kaon1_Eta->Fill(mykaon1.Eta());
    Hist_Kaon1_Phi->Fill(mykaon1.Phi());
    Hist_Kaon2_pt->Fill(mykaon2.Pt());
    Hist_Kaon2_Eta->Fill(mykaon2.Eta());
    Hist_Kaon2_Phi->Fill(mykaon2.Phi());
    Hist_Kaon1PtvsKaon2Pt->Fill(mykaon1.Pt(),mykaon2.Pt());

    /// muon pairs from original muons
    Hist_mumu_mass->Fill(mu1mu2.M());
    Hist_mumu_pt->Fill(mu1mu2.Pt());
    Hist_mumu_Eta->Fill(mu1mu2.Eta());
    Hist_mumu_Phi->Fill(mu1mu2.Phi());
    Hist_kk_mass->Fill(ka1ka2.M());
    Hist_kk_pt->Fill(ka1ka2.Pt());
    Hist_kk_Eta->Fill(ka1ka2.Eta());
    Hist_kk_Phi->Fill(ka1ka2.Phi());

    /// refit muons from Jpsi & refit kaons from Phi & Bs0_Y
    Hist_JPsi_mass->Fill(JPsi.M());
    Hist_JPsi_pt->Fill(JPsi.Pt());
    Hist_JPsi_Eta->Fill(JPsi.Eta());
    Hist_JPsi_Phi->Fill(JPsi.Phi());
    Hist_JPsi_PtvsEta->Fill(JPsi.Eta(),JPsi.Pt());
    Hist_Phi_mass->Fill(Phi.M());
    Hist_Phi_pt->Fill(Phi.Pt());
    Hist_Phi_Eta->Fill(Phi.Eta());
    Hist_Phi_Phi->Fill(Phi.Phi());
    Hist_Phi_PtvsEta->Fill(Phi.Eta(),Phi.Pt());
    Hist_Bs0_Y_mass_fromNTuple->Fill((*XMass)[myXIdx]);
    Hist_Bs0_Y_mass->Fill(Bs0_Y.M());
    Hist_Bs0_Y_pt->Fill(Bs0_Y.Pt());
    Hist_Bs0_Y_Eta->Fill(Bs0_Y.Eta());
    Hist_Bs0_Y_Phi->Fill(Bs0_Y.Phi());

    Hist_CTau_CTauE_PV->Fill((*XCTauPV)[myXIdx] / (*XCTauPVE)[myXIdx]);
    Hist_Lxy_LxyE_PV->Fill((*XLxyPV)[myXIdx] / (*XLxyPVE)[myXIdx]);
    Hist_CTauVsLxy_PV->Fill(((*XLxyPV)[myXIdx] / (*XLxyPVE)[myXIdx]),((*XCTauPV)[myXIdx] / (*XCTauPVE)[myXIdx]));                        ;
    /*
    Hist_CTau_CTauE_BS->Fill((*XCTauBS)[myXIdx] / (*XCTauBSE)[myXIdx]);
    Hist_Lxy_LxyE_BS->Fill((*XLxyBS)[myXIdx] / (*XLxyBSE)[myXIdx]);
    Hist_CTauVsLxy_BS->Fill(((*XLxyBS)[myXIdx] / (*XLxyBSE)[myXIdx]),((*XCTauBS)[myXIdx] / (*XCTauBSE)[myXIdx]));

    Hist_CTau_CTauE_XLessPV->Fill((*XCTauXLessPV)[myXIdx] / (*XCTauXLessPVE)[myXIdx]);
    Hist_Lxy_LxyE_XLessPV->Fill((*XLxyXLessPV)[myXIdx] / (*XLxyXLessPVE)[myXIdx]);
    Hist_CTauVsLxy_XLessPV->Fill(((*XLxyXLessPV)[myXIdx] / (*XLxyXLessPVE)[myXIdx]),((*XCTauXLessPV)[myXIdx] / (*XCTauXLessPVE)[myXIdx]));

    Hist_CTau_CTauE_PVCosAlpha->Fill((*XCTauPVCosAlpha)[myXIdx] / (*XCTauPVCosAlphaE)[myXIdx]);
    Hist_Lxy_LxyE_PVCosAlpha->Fill((*XLxyPVCosAlpha)[myXIdx] / (*XLxyPVCosAlphaE)[myXIdx]);
    Hist_CTauVsLxy_PVCosAlpha->Fill(((*XLxyPVCosAlpha)[myXIdx] / (*XLxyPVCosAlphaE)[myXIdx]),((*XCTauPVCosAlpha)[myXIdx] / (*XCTauPVCosAlphaE)[myXIdx]));

    Hist_CTau_CTauE_PVCosAlpha3D->Fill((*XCTauPVCosAlpha3D)[myXIdx] / (*XCTauPVCosAlpha3DE)[myXIdx]);
    Hist_Lxy_LxyE_PVCosAlpha3D->Fill((*XLxyPVCosAlpha3D)[myXIdx] / (*XLxyPVCosAlpha3DE)[myXIdx]);
    Hist_CTauVsLxy_PVCosAlpha3D->Fill(((*XLxyPVCosAlpha3D)[myXIdx] / (*XLxyPVCosAlpha3DE)[myXIdx]),((*XCTauPVCosAlpha3D)[myXIdx] / (*XCTauPVCosAlpha3DE)[myXIdx]));

    Hist_CTau_CTauE_XLessPVCosAlpha->Fill((*XCTauXLessPVCosAlpha)[myXIdx] / (*XCTauXLessPVCosAlphaE)[myXIdx]);
    Hist_Lxy_LxyE_XLessPVCosAlpha->Fill((*XLxyXLessPVCosAlpha)[myXIdx] / (*XLxyXLessPVCosAlphaE)[myXIdx]);
    Hist_CTauVsLxy_XLessPVCosAlpha->Fill(((*XLxyXLessPVCosAlpha)[myXIdx] / (*XLxyXLessPVCosAlphaE)[myXIdx]),((*XCTauXLessPVCosAlpha)[myXIdx] / (*XCTauXLessPVCosAlphaE)[myXIdx]));

    Hist_CTau_CTauE_XLessPVCosAlpha3D->Fill((*XCTauXLessPVCosAlpha3D)[myXIdx] / (*XCTauXLessPVCosAlpha3DE)[myXIdx]);
    Hist_Lxy_LxyE_XLessPVCosAlpha3D->Fill((*XLxyXLessPVCosAlpha3D)[myXIdx] / (*XLxyXLessPVCosAlpha3DE)[myXIdx]);
    Hist_CTauVsLxy_XLessPVCosAlpha3D->Fill(((*XLxyXLessPVCosAlpha3D)[myXIdx] / (*XLxyXLessPVCosAlpha3DE)[myXIdx]),((*XCTauXLessPVCosAlpha3D)[myXIdx] / (*XCTauXLessPVCosAlpha3DE)[myXIdx]));

    */


    if ((Bs0_Y.M() > 4.05) && (Bs0_Y.M() < 4.8)) {
      Hist_Y4140_mass->Fill((*XMass)[myXIdx]);
      Hist_Y4140_pt->Fill(Bs0_Y.Pt());
      Hist_Y4140_Eta->Fill(Bs0_Y.Eta());
      Hist_Y4140_Phi->Fill(Bs0_Y.Phi());
      Hist_Y4140_PtvsEta->Fill(Bs0_Y.Eta(),Bs0_Y.Pt());
    }

    else if ((Bs0_Y.M() > 5.15) && (Bs0_Y.M() < 5.55)) {
      Hist_Bs0_mass->Fill((*XMass)[myXIdx]);
      Hist_Bs0_pt->Fill(Bs0_Y.Pt());
      Hist_Bs0_Eta->Fill(Bs0_Y.Eta());
      Hist_Bs0_Phi->Fill(Bs0_Y.Phi());
      Hist_Bs0_PtvsEta->Fill(Bs0_Y.Eta(),Bs0_Y.Pt());
    }



    ////////////////// Some cuts for Y4140 & Bs0 //////////////////

    //////////// HLT8 Cut ////////////
    /// We use HLT_4(displaced) for Bs0: non-prompt, HLT_8 for Y(4140): prompt
    if (HLT_8_vAll) {
      Hist_JPsi_mass_HLT->Fill(mu1mu2.M());
      if ((Bs0_Y.M() > 4.05) && (Bs0_Y.M() < 4.8)) {
        Hist_Y4140_mass_HLT->Fill((*XMass)[myXIdx]);
      }
      else if ((Bs0_Y.M() > 5.15) && (Bs0_Y.M() < 5.55)) {
        Hist_Bs0_mass_HLT->Fill((*XMass)[myXIdx]);
      }


      //////////// SoftMuID Cut ////////////
      if ( 1
        && ( ((*muQual)[mymupIdx]) & (1 << muonQual[3]) ) && ( ((*muQual)[mymumIdx]) & (1 << muonQual[3]) )
        && ( ( (*muChi2)[mymupIdx] / (*muNDF)[mymupIdx] ) < 3 ) &&( ( (*muChi2)[mymumIdx] / (*muNDF)[mymumIdx] ) < 3 )
        && fabs((*muDzVtx)[mymupIdx]) < 20.0 && fabs((*muDzVtx)[mymumIdx]) < 20.0
        && fabs((*muDxyVtx)[mymupIdx]) < 3.0 && fabs((*muDxyVtx)[mymumIdx]) < 3.0
        && (*muPhits)[mymupIdx] > 0 && (*muPhits)[mymumIdx] > 0
        && (*muShits)[mymupIdx] > 5 && (*muShits)[mymumIdx] > 5
      ) {
        Hist_JPsi_mass_SoftMu->Fill(mu1mu2.M());

        if ((Bs0_Y.M() > 4.05) && (Bs0_Y.M() < 4.8)) {
          Hist_Y4140_mass_SoftMu->Fill((*XMass)[myXIdx]);
        }
        else if ((Bs0_Y.M() > 5.15) && (Bs0_Y.M() < 5.55)) {
          Hist_Bs0_mass_SoftMu->Fill((*XMass)[myXIdx]);
        }


        Bool_t mu1_PtSel = false;
        Bool_t mu2_PtSel = false;
        if ((((fabs(mu1.Eta()) < 1.2) && (mu1.Pt() > 4.))) || ((mu1.Eta() >= 1.2 || mu1.Eta() <= -1.2 ) && (mu1.Pt() > 3.3))) mu1_PtSel = true;
        if ((((fabs(mu2.Eta()) < 1.2) && (mu2.Pt() > 4.))) || ((mu2.Eta() >= 1.2 || mu2.Eta() <= -1.2 ) && (mu2.Pt() > 3.3))) mu2_PtSel = true;

        //////////// JPsi Cut ////////////
        if (1
          && JPsi.Pt() > 8.0
          && fabs(JPsi.M() - JPsi_mass) < 0.12
          && ((*MuMuVtx_CL)[myJPsiIdx]) > 0.01
          && fabs(mu1.Eta()) < 2.4 && fabs(mu2.Eta()) < 2.4
          && mu1.Pt() > 3.3 && mu2.Pt() > 3.3
          && mu1_PtSel && mu2_PtSel
        ) {

          Hist_JPsi_PtvsEta_JPsi->Fill(JPsi.Eta(),JPsi.Pt());
          Hist_Phi_PtvsEta_JPsi->Fill(Phi.Eta(),Phi.Pt());

          if ((Bs0_Y.M() > 4.05) && (Bs0_Y.M() < 4.8)) {
            Hist_Y4140_mass_JPsi->Fill((*XMass)[myXIdx]);
            Hist_Y4140_PtvsEta_JPsi->Fill(Bs0_Y.Eta(),Bs0_Y.Pt());
          }
          else if ((Bs0_Y.M() > 5.15) && (Bs0_Y.M() < 5.55)) {
            Hist_Bs0_mass_JPsi->Fill((*XMass)[myXIdx]);
            Hist_Bs0_PtvsEta_JPsi->Fill(Bs0_Y.Eta(),Bs0_Y.Pt());
          }

          if (((*XLxyPV)[myXIdx] / (*XLxyPVE)[myXIdx]) < 2) {
            if ((Bs0_Y.M() > 4.05) && (Bs0_Y.M() < 4.8)) {
              Hist_Y4140_mass_Prompt->Fill((*XMass)[myXIdx]);
            }
            else if ((Bs0_Y.M() > 5.15) && (Bs0_Y.M() < 5.55)) {
              Hist_Bs0_mass_Prompt->Fill((*XMass)[myXIdx]);
            }
          } /// Prompt cuts

          else if (((*XLxyPV)[myXIdx] / (*XLxyPVE)[myXIdx]) > 2 && ((*XLxyPV)[myXIdx] / (*XLxyPVE)[myXIdx])  < 3.5) {
            if ((Bs0_Y.M() > 4.05) && (Bs0_Y.M() < 4.8)) {
              Hist_Y4140_mass_Mixed->Fill((*XMass)[myXIdx]);
            }
            else if ((Bs0_Y.M() > 5.15) && (Bs0_Y.M() < 5.55)) {
              Hist_Bs0_mass_Mixed->Fill((*XMass)[myXIdx]);
            }
          } /// Mixed cuts

          else if (((*XLxyPV)[myXIdx] / (*XLxyPVE)[myXIdx]) > 3.5) {
            if ((Bs0_Y.M() > 4.05) && (Bs0_Y.M() < 4.8)) {
              Hist_Y4140_mass_NonPrompt->Fill((*XMass)[myXIdx]);
            }
            else if ((Bs0_Y.M() > 5.15) && (Bs0_Y.M() < 5.55)) {
              Hist_Bs0_mass_NonPrompt->Fill((*XMass)[myXIdx]);
            }
          } /// NonPrompt cuts


          //////////// Tracks Cut ////////////
          if (1
            && ((*trackChi2)[kaon1Idx] / (*trackNDF)[kaon1Idx]) < 5.0
            && (*trackPhits)[kaon1Idx] > 0
            && (*trackShits)[kaon1Idx] >= 7
            && ((*trackChi2)[kaon2Idx] / (*trackNDF)[kaon2Idx]) < 5.0
            && (*trackPhits)[kaon2Idx] > 0
            && (*trackShits)[kaon2Idx] >= 7
          ) {

            if ((Bs0_Y.M() > 4.05) && (Bs0_Y.M() < 4.8)) {
              Hist_Y4140_mass_Track->Fill((*XMass)[myXIdx]);
            }
            else if ((Bs0_Y.M() > 5.15) && (Bs0_Y.M() < 5.55)) {
              Hist_Bs0_mass_Track->Fill((*XMass)[myXIdx]);
            }

            if (((*XLxyPV)[myXIdx] / (*XLxyPVE)[myXIdx]) < 2) {
              if ((Bs0_Y.M() > 4.05) && (Bs0_Y.M() < 4.8)) {
                Hist_Y4140_mass_TrackP->Fill((*XMass)[myXIdx]);
              }
              else if ((Bs0_Y.M() > 5.15) && (Bs0_Y.M() < 5.55)) {
                Hist_Bs0_mass_TrackP->Fill((*XMass)[myXIdx]);
              }
            } /// Prompt cuts

            else if (((*XLxyPV)[myXIdx] / (*XLxyPVE)[myXIdx]) > 2 && ((*XLxyPV)[myXIdx] / (*XLxyPVE)[myXIdx])  < 3.5) {
              if ((Bs0_Y.M() > 4.05) && (Bs0_Y.M() < 4.8)) {
                Hist_Y4140_mass_TrackM->Fill((*XMass)[myXIdx]);
              }
              else if ((Bs0_Y.M() > 5.15) && (Bs0_Y.M() < 5.55)) {
                Hist_Bs0_mass_TrackM->Fill((*XMass)[myXIdx]);
              }
            } /// Mixed cuts

            else if (((*XLxyPV)[myXIdx] / (*XLxyPVE)[myXIdx]) > 3.5) {
              if ((Bs0_Y.M() > 4.05) && (Bs0_Y.M() < 4.8)) {
                Hist_Y4140_mass_TrackNP->Fill((*XMass)[myXIdx]);
              }
              else if ((Bs0_Y.M() > 5.15) && (Bs0_Y.M() < 5.55)) {
                Hist_Bs0_mass_TrackNP->Fill((*XMass)[myXIdx]);
              }
            } /// NonPrompt cuts


            //////////// Phi Cut ////////////
            if ((ka1ka2.M()) > 1.013 && (ka1ka2.M()) < 1.026) {  /// 1.008-1.035 first cut // 1.0125-1.0265 CMS B+ // 1.012-1.029 D0 // 1.013-1.026 last cut
              if ((Bs0_Y.M() > 4.05) && (Bs0_Y.M() < 4.8)) {
                Hist_Y4140_mass_Phi->Fill((*XMass)[myXIdx]);
              }
              else if ((Bs0_Y.M() > 5.15) && (Bs0_Y.M() < 5.55)) {
                Hist_Bs0_mass_Phi->Fill((*XMass)[myXIdx]);
              }

              if (((*XLxyPV)[myXIdx] / (*XLxyPVE)[myXIdx]) < 2) {
                if ((Bs0_Y.M() > 4.05) && (Bs0_Y.M() < 4.8)) {
                  Hist_Y4140_mass_PhiP->Fill((*XMass)[myXIdx]);
                }
                else if ((Bs0_Y.M() > 5.15) && (Bs0_Y.M() < 5.55)) {
                  Hist_Bs0_mass_PhiP->Fill((*XMass)[myXIdx]);
                }
              } /// Prompt cuts

              else if (((*XLxyPV)[myXIdx] / (*XLxyPVE)[myXIdx]) > 2 && ((*XLxyPV)[myXIdx] / (*XLxyPVE)[myXIdx])  < 3.5) {
                if ((Bs0_Y.M() > 4.05) && (Bs0_Y.M() < 4.8)) {
                  Hist_Y4140_mass_PhiM->Fill((*XMass)[myXIdx]);
                }
                else if ((Bs0_Y.M() > 5.15) && (Bs0_Y.M() < 5.55)) {
                  Hist_Bs0_mass_PhiM->Fill((*XMass)[myXIdx]);
                }
              } /// Mixed cuts

              else if (((*XLxyPV)[myXIdx] / (*XLxyPVE)[myXIdx]) > 3.5) {

                if ((Bs0_Y.M() > 4.05) && (Bs0_Y.M() < 4.8)) {
                  Hist_Y4140_mass_PhiNP->Fill((*XMass)[myXIdx]);
                }
                else if ((Bs0_Y.M() > 5.15) && (Bs0_Y.M() < 5.55)) {
                  Hist_Bs0_mass_PhiNP->Fill((*XMass)[myXIdx]);
                }
              } /// NonPrompt cuts


              //////////// Bs0 Cut ////////////
              if (1
                && fabs((*XCosAlphaPV)[myXIdx]) > 0.9
                && ((*XVtx_CL)[myXIdx]) > 0.01
              ) {

                if ((Bs0_Y.M() > 4.05) && (Bs0_Y.M() < 4.8)) {
                  Hist_Y4140_mass_Bs0->Fill((*XMass)[myXIdx]);
                }
                else if ((Bs0_Y.M() > 5.15) && (Bs0_Y.M() < 5.55)) {
                  Hist_Bs0_mass_Bs0->Fill((*XMass)[myXIdx]);
                }

                if (((*XLxyPV)[myXIdx] / (*XLxyPVE)[myXIdx]) < 2) {
                  if ((Bs0_Y.M() > 4.05) && (Bs0_Y.M() < 4.8)) {
                    //if ((ka1ka2.M()) > 1.013 && (ka1ka2.M()) < 1.026) {
                    Hist_Y4140_mass_Bs0P->Fill((*XMass)[myXIdx]);
                    //}
                    Hist_SC_JPsi_mass_Bs0P->Fill(mu1mu2.M());
                    Hist_SC_Phi_mass_Bs0P->Fill(Phi.M());
                  }
                  else if ((Bs0_Y.M() > 5.15) && (Bs0_Y.M() < 5.55)) {
                    //if ((ka1ka2.M()) > 1.013 && (ka1ka2.M()) < 1.026) {
                    Hist_Bs0_mass_Bs0P->Fill((*XMass)[myXIdx]);
                    //}
                    Hist_CC_JPsi_mass_Bs0P->Fill(mu1mu2.M());
                    Hist_CC_Phi_mass_Bs0P->Fill(Phi.M());
                  }
                } /// Prompt cuts

                else if (((*XLxyPV)[myXIdx] / (*XLxyPVE)[myXIdx]) > 2 && ((*XLxyPV)[myXIdx] / (*XLxyPVE)[myXIdx])  < 3) {
                  if ((Bs0_Y.M() > 4.05) && (Bs0_Y.M() < 4.8)) {
                    //if ((ka1ka2.M()) > 1.013 && (ka1ka2.M()) < 1.026) {
                    Hist_Y4140_mass_Bs0M->Fill((*XMass)[myXIdx]);
                    //}
                    Hist_SC_JPsi_mass_Bs0M->Fill(mu1mu2.M());
                    Hist_SC_Phi_mass_Bs0M->Fill(Phi.M());
                  }
                  else if ((Bs0_Y.M() > 5.15) && (Bs0_Y.M() < 5.55)) {
                    //if ((ka1ka2.M()) > 1.013 && (ka1ka2.M()) < 1.026) {
                    Hist_Bs0_mass_Bs0M->Fill((*XMass)[myXIdx]);
                    //}
                    Hist_CC_JPsi_mass_Bs0M->Fill(mu1mu2.M());
                    Hist_CC_Phi_mass_Bs0M->Fill(Phi.M());
                  }
                } /// Mixed cuts

                else if (((*XLxyPV)[myXIdx] / (*XLxyPVE)[myXIdx]) >= 3.0) {
                  if ((Bs0_Y.M() > 4.05) && (Bs0_Y.M() < 4.8)) {
                    //if ((ka1ka2.M()) > 1.013 && (ka1ka2.M()) < 1.026) {
                    Hist_Y4140_mass_Bs0NP->Fill((*XMass)[myXIdx]);
                    //}
                    Hist_SC_JPsi_mass_Bs0NP->Fill(mu1mu2.M());
                    Hist_SC_Phi_mass_Bs0NP->Fill(Phi.M());
                  }
                  else if ((Bs0_Y.M() > 5.15) && (Bs0_Y.M() < 5.55)) {
                    //if ((ka1ka2.M()) > 1.013 && (ka1ka2.M()) < 1.026) {
                    Hist_Bs0_mass_Bs0NP->Fill((*XMass)[myXIdx]);
                    //}
                    Hist_CC_JPsi_mass_Bs0NP->Fill(mu1mu2.M());
                    Hist_CC_Phi_mass_Bs0NP->Fill(Phi.M());
                  }
                } /// NonPrompt cuts


              } /// Bs0 cut
            } /// Phi cut
          } /// Track cut
        } /// JPsi cut
      } /// softMu cut
    } /// HLT cut
  } /// for loop


  return kTRUE;
}


void Y4140::SlaveTerminate()
{
  //std::cout << "SLAVE TERMINATE" << std::endl;


  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

  /// SEMRA added
  TDirectory *savedir = gDirectory;
  if (fOut)
  {
    fOut->cd();
    gStyle->SetOptStat(111111) ;


    ////////////////// Write Histograms //////////////////
    /// refit muons & kaons
    Hist_Muon1_pt->Write();
    Hist_Muon1_Eta->Write();
    Hist_Muon1_Phi->Write();
    Hist_Muon2_pt->Write();
    Hist_Muon2_Eta->Write();
    Hist_Muon2_Phi->Write();
    Hist_Muon1PtvsMuon2Pt->Write();
    Hist_Kaon1_pt->Write();
    Hist_Kaon1_Eta->Write();
    Hist_Kaon1_Phi->Write();
    Hist_Kaon2_pt->Write();
    Hist_Kaon2_Eta->Write();
    Hist_Kaon2_Phi->Write();
    Hist_Kaon1PtvsKaon2Pt->Write();

    /// muon pairs from original muons
    Hist_mumu_mass->Write();
    Hist_mumu_pt->Write();
    Hist_mumu_Eta->Write();
    Hist_mumu_Phi->Write();
    Hist_kk_mass->Write();
    Hist_kk_pt->Write();
    Hist_kk_Eta->Write();
    Hist_kk_Phi->Write();

    /// refit muons from Jpsi & refit kaons from Phi & Bs0_Y
    Hist_JPsi_mass->Write();
    Hist_JPsi_pt->Write();
    Hist_JPsi_Eta->Write();
    Hist_JPsi_Phi->Write();
    Hist_JPsi_PtvsEta->Write();
    Hist_JPsi_mass_HLT->Write();
    Hist_JPsi_mass_SoftMu->Write();
    Hist_JPsi_PtvsEta_JPsi->Write();
    Hist_Phi_mass->Write();
    Hist_Phi_pt->Write();
    Hist_Phi_Eta->Write();
    Hist_Phi_Phi->Write();
    Hist_Phi_PtvsEta->Write();
    Hist_Phi_PtvsEta_JPsi->Write();

    /// Bs0_Y
    Hist_Bs0_Y_mass_fromNTuple->Write();
    Hist_Bs0_Y_mass->Write();
    Hist_Bs0_Y_pt->Write();
    Hist_Bs0_Y_Eta->Write();
    Hist_Bs0_Y_Phi->Write();

    /// Y4140 & Some Cuts
    Hist_Y4140_mass->Write();
    Hist_Y4140_pt->Write();
    Hist_Y4140_Eta->Write();
    Hist_Y4140_Phi->Write();
    Hist_Y4140_PtvsEta->Write();

    Hist_Y4140_mass_HLT->Write();
    Hist_Y4140_mass_SoftMu->Write();
    Hist_Y4140_mass_JPsi->Write();
    Hist_Y4140_PtvsEta_JPsi->Write();
    Hist_Y4140_mass_Prompt->Write();
    Hist_Y4140_mass_Mixed->Write();
    Hist_Y4140_mass_NonPrompt->Write();
    Hist_Y4140_mass_Track->Write();
    Hist_Y4140_mass_TrackP->Write();
    Hist_Y4140_mass_TrackM->Write();
    Hist_Y4140_mass_TrackNP->Write();
    Hist_Y4140_mass_Phi->Write();
    Hist_Y4140_mass_PhiP->Write();
    Hist_Y4140_mass_PhiM->Write();
    Hist_Y4140_mass_PhiNP->Write();
    Hist_Y4140_mass_Bs0->Write();
    Hist_Y4140_mass_Bs0P->Write();
    Hist_SC_JPsi_mass_Bs0P->Write();
    Hist_SC_Phi_mass_Bs0P->Write();
    Hist_Y4140_mass_Bs0M->Write();
    Hist_SC_JPsi_mass_Bs0M->Write();
    Hist_SC_Phi_mass_Bs0M->Write();
    Hist_Y4140_mass_Bs0NP->Write();
    Hist_SC_JPsi_mass_Bs0NP->Write();
    Hist_SC_Phi_mass_Bs0NP->Write();

    /// Bs0 & Some Cuts
    Hist_Bs0_mass->Write();
    Hist_Bs0_pt->Write();
    Hist_Bs0_Eta->Write();
    Hist_Bs0_Phi->Write();
    Hist_Bs0_PtvsEta->Write();

    Hist_Bs0_mass_HLT->Write();
    Hist_Bs0_mass_SoftMu->Write();
    Hist_Bs0_mass_JPsi->Write();
    Hist_Bs0_PtvsEta_JPsi->Write();
    Hist_Bs0_mass_Prompt->Write();
    Hist_Bs0_mass_Mixed->Write();
    Hist_Bs0_mass_NonPrompt->Write();
    Hist_Bs0_mass_Track->Write();
    Hist_Bs0_mass_TrackP->Write();
    Hist_Bs0_mass_TrackM->Write();
    Hist_Bs0_mass_TrackNP->Write();
    Hist_Bs0_mass_Phi->Write();
    Hist_Bs0_mass_PhiP->Write();
    Hist_Bs0_mass_PhiM->Write();
    Hist_Bs0_mass_PhiNP->Write();
    Hist_Bs0_mass_Bs0->Write();
    Hist_Bs0_mass_Bs0P->Write();
    Hist_CC_JPsi_mass_Bs0P->Write();
    Hist_CC_Phi_mass_Bs0P->Write();
    Hist_Bs0_mass_Bs0M->Write();
    Hist_CC_JPsi_mass_Bs0M->Write();
    Hist_CC_Phi_mass_Bs0M->Write();
    Hist_Bs0_mass_Bs0NP->Write();
    Hist_CC_JPsi_mass_Bs0NP->Write();
    Hist_CC_Phi_mass_Bs0NP->Write();

    Hist_CTau_CTauE_PV->Write();
    Hist_Lxy_LxyE_PV->Write();
    Hist_CTauVsLxy_PV->Write();

    /*
    Hist_CTau_CTauE_BS->Write();
    Hist_Lxy_LxyE_BS->Write();
    Hist_CTauVsLxy_BS->Write();

    Hist_CTau_CTauE_XLessPV->Write();
    Hist_Lxy_LxyE_XLessPV->Write();
    Hist_CTauVsLxy_XLessPV->Write();

    Hist_CTau_CTauE_PVCosAlpha->Write();
    Hist_Lxy_LxyE_PVCosAlpha->Write();
    Hist_CTauVsLxy_PVCosAlpha->Write();

    Hist_CTau_CTauE_PVCosAlpha3D->Write();
    Hist_Lxy_LxyE_PVCosAlpha3D->Write();
    Hist_CTauVsLxy_PVCosAlpha3D->Write();

    Hist_CTau_CTauE_XLessPVCosAlpha->Write();
    Hist_Lxy_LxyE_XLessPVCosAlpha->Write();
    Hist_CTauVsLxy_XLessPVCosAlpha->Write();

    Hist_CTau_CTauE_XLessPVCosAlpha3D->Write();
    Hist_Lxy_LxyE_XLessPVCosAlpha3D->Write();
    Hist_CTauVsLxy_XLessPVCosAlpha3D->Write();
    */

    // Xmultiplicity->Write();
    // nIndex->Write();
    // psiIndex->Write();
    // psiIndexNx->Write();

    OutFile->Print();
    fOutput->Add(OutFile);
    gDirectory = savedir;
    fOut->Close();

  }
}


void Y4140::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

}
