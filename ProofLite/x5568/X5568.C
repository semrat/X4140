#define X5568_cxx
// The class definition in X5568.h has been generated automatically
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
// Root > T->Process("X5568.C")
// Root > T->Process("X5568.C","some options")
// Root > T->Process("X5568.C+")
//

#include "X5568.h"
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


void X5568::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

}


void X5568::SlaveBegin(TTree * /*tree*/)
{

  TString option = GetOption();

  /// SEMRA added
  OutFile = new TProofOutputFile( "X5568.root" );
  fOut = OutFile->OpenFile("RECREATE");
  if (!(fOut=OutFile->OpenFile("RECREATE")))
  {
    Warning("SlaveBegin","Problems opening file: %s%s", OutFile->GetDir(), OutFile->GetFileName() );
  }

  ////////////////// Histograms //////////////////
  JPsi_mass = 3.0969; /// pdg mass
  Phi_mass = 1.020; /// pdg mass

  Hist_Bs0_mass_Bs0NP  = new TH1F ("Bs0_mass_Bs0NP","Bs0_mass_Bs0NP;M(#mu^{-}#mu^{+}K^{-}K^{+}) [GeV];Entries",80,5.15, 5.55);

  Xmultiplicity = new TH1F ("Xmultiplicity","Xmultiplicity",11, -0.5, 10.5); /// for counting how many X
  nIndex = new TH1F ("nIndex","nIndex",11, -0.5, 10.5);
  psiIndex = new TH1F ("psiIndex","psiIndex",11, -0.5, 10.5);
  psiIndexNx = new TH1F ("psiIndexNx","psiIndexNx",11, -0.5, 10.5);

}

Bool_t X5568::Process(Long64_t entry)
{
  //std::cout << "PROCESS" << std::endl;

  GetEntry(entry); /// ROOT event loop

  std::vector<bool> cuts;
  Double_t muon_mass = 0.1056;
  Double_t kaonCh_mass = 0.4936;
  Double_t Bs0_Low_Mass = 5.1;
  Double_t Bs0_High_Mass = 5.6;
  Double_t Y_High_Mass = 4.35;

  Bool_t HLT_4_v9 = false, HLT_4_v10 = false, HLT_4_v11 = false, HLT_4_v12 = false ;
  Bool_t HLT_8_v3 = false, HLT_8_v4 = false, HLT_8_v5 = false, HLT_8_v6 = false, HLT_8_v7 = false ;
  Bool_t HLT_4_vAll = false, HLT_4_vAny = false;
  Bool_t HLT_8_vAll = false, HLT_8_vAny = false;
  Bool_t HLT_All = false, HLT_Any = false;

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

  if (HLT_4_v9 || HLT_4_v10 || HLT_4_v11 || HLT_4_v12) HLT_4_vAny = true;
  if (HLT_8_v3 || HLT_8_v4 || HLT_8_v5 || HLT_8_v6 || HLT_8_v7) HLT_8_vAny = true;
  if (HLT_4_v9 && HLT_4_v10 && HLT_4_v11 && HLT_4_v12) HLT_4_vAll = true;
  if (HLT_8_v3 && HLT_8_v4 && HLT_8_v5 && HLT_8_v6 && HLT_8_v7) HLT_8_vAll = true;
  if (HLT_4_vAny || HLT_8_vAny) HLT_Any = true;

  int muonQual[] = {1,3,4,12};




  return kTRUE;
}


void X5568::SlaveTerminate()
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

    /// X5568 & Some Cuts
    Hist_X5568_mass->Write();
    Hist_X5568_pt->Write();
    Hist_X5568_Eta->Write();
    Hist_X5568_Phi->Write();
    Hist_X5568_PtvsEta->Write();

    Hist_X5568_mass_HLT->Write();
    Hist_X5568_mass_SoftMu->Write();
    Hist_X5568_mass_JPsi->Write();
    Hist_X5568_PtvsEta_JPsi->Write();
    Hist_X5568_mass_Prompt->Write();
    Hist_X5568_mass_Mixed->Write();
    Hist_X5568_mass_NonPrompt->Write();
    Hist_X5568_mass_Track->Write();
    Hist_X5568_mass_TrackP->Write();
    Hist_X5568_mass_TrackM->Write();
    Hist_X5568_mass_TrackNP->Write();
    Hist_X5568_mass_Phi->Write();
    Hist_X5568_mass_PhiP->Write();
    Hist_X5568_mass_PhiM->Write();
    Hist_X5568_mass_PhiNP->Write();
    Hist_X5568_mass_Bs0->Write();
    Hist_X5568_mass_Bs0P->Write();
    Hist_SC_JPsi_mass_Bs0P->Write();
    Hist_SC_Phi_mass_Bs0P->Write();
    Hist_X5568_mass_Bs0M->Write();
    Hist_SC_JPsi_mass_Bs0M->Write();
    Hist_SC_Phi_mass_Bs0M->Write();
    Hist_X5568_mass_Bs0NP->Write();
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


void X5568::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

}
