#define JPsiCount_cxx
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

#include "JPsiCount.h"
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

void JPsiCount::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

}


void JPsiCount::SlaveBegin(TTree * /*tree*/)
{
  //std::cout << "SLAVE BEGIN" << std::endl;

  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();
  OutFile = new TProofOutputFile( "JPsiCount.root" );
  fOut = OutFile->OpenFile("RECREATE");
  if (!(fOut=OutFile->OpenFile("RECREATE")))
  {
    Warning("SlaveBegin","Problems opening file: %s%s", OutFile->GetDir(), OutFile->GetFileName() );
  }

  ////////////////// Histograms //////////////////
  JPsi_mass = 3.0969; /// pdg mass
  Phi_mass = 1.020; /// pdg mass

  /// refit muons & kaons
  JPsi_Mass      = new TH1F ("JPsi_Mass_NoHLT", "JPsi_Mass; M(#mu^{-}#mu^{+}) [GeV];Entries",200, 2.5, 3.5);
  JPsi_Mass_HLT8 = new TH1F ("JPsi_Mass_HLT8","JPsi_Mass_HLT8; M(#mu^{-}#mu^{+}) [GeV];Entries",200, 2.5, 3.5);
  JPsi_Mass_HLT4 = new TH1F ("JPsi_Mass_HLT4","JPsi_Mass_HLT4; M(#mu^{-}#mu^{+}) [GeV];Entries",200, 2.5, 3.5);

  JPsi_Mass_ref      = new TH1F ("JPsi_Mass_NJPsi_Mass_refoHLT", "JPsi_Mass_ref; M(#mu^{-}#mu^{+}) [GeV];Entries",200, 2.5, 3.5);
  JPsi_Mass_ref_HLT8 = new TH1F ("JPsi_Mass_ref_HLT8","JPsi_Mass_ref_HLT8; M(#mu^{-}#mu^{+}) [GeV];Entries",200, 2.5, 3.5);
  JPsi_Mass_ref_HLT4 = new TH1F ("JPsi_Mass_ref_HLT4","JPsi_Mass_ref_HLT4; M(#mu^{-}#mu^{+}) [GeV];Entries",200, 2.5, 3.5);

  B0_Mass       = new TH1F ("B0_Mass","B0_Mass",200, 5.1, 5.6);
  B0_Mass_HLT8  = new TH1F ("B0_Mass_HLT8","B0_Mass_HLT8",200, 5.1, 5.6);
  B0_Mass_HLT4  = new TH1F ("B0_Mass_HLT4","B0_Mass_HLT4",200, 5.1, 5.6);

  B0_Mass_ref       = new TH1F ("B0_Mass_ref","B0_Mass_ref",200, 5.1, 5.6);
  B0_Mass_HLT8_ref  = new TH1F ("B0_Mass_HLT8_ref","B0_Mass_HLT8_ref",200, 5.1, 5.6);
  B0_Mass_HLT4_ref  = new TH1F ("B0_Mass_HLT4_ref","B0_Mass_HLT4_ref",200, 5.1, 5.6);

  JPsi_vs_run       = new TH1F ("JPsi_vs_run", "JPsi_vs_run; Run[#];J/Psi[#]",20000, 190000, 210000);
  JPsi_vs_run_HLT4  = new TH1F ("JPsi_vs_run_HLT4", "JPsi_vs_run_HLT4;  Run[#];J/Psi[#]",20000, 190000, 210000);
  JPsi_vs_run_HLT8  = new TH1F ("JPsi_vs_run_HLT8", "JPsi_vs_run_HLT8;  Run[#];J/Psi[#]",20000, 190000, 210000);

  JPsi_vs_run_ref       = new TH1F ("JPsi_vs_run_ref", "JPsi_vs_run_ref; Run[#];J/Psi[#]",20000, 190000, 210000);
  JPsi_vs_run_ref_HLT4  = new TH1F ("JPsi_vs_run_ref_HLT4", "JPsi_vs_run_ref_HLT4;  Run[#];J/Psi[#]",20000, 190000, 210000);
  JPsi_vs_run_ref_HLT8  = new TH1F ("JPsi_vs_run_ref_HLT8", "JPsi_vs_run_ref_HLT8;  Run[#];J/Psi[#]",20000, 190000, 210000);

  B0s_vs_run       = new TH1F ("B0s_vs_run", "B0s_vs_run; Run[#];J/Psi[#]",20000, 190000, 210000);
  B0s_vs_run_HLT4  = new TH1F ("B0s_vs_run_HLT4", "B0s_vs_run_HLT4;  Run[#];J/Psi[#]",20000, 190000, 210000);
  B0s_vs_run_HLT8  = new TH1F ("B0s_vs_run_HLT8", "B0s_vs_run_HLT8;  Run[#];J/Psi[#]",20000, 190000, 210000);

  B0s_vs_run_ref       = new TH1F ("B0s_vs_run_ref", "B0s_vs_run_ref; Run[#];J/Psi[#]",20000, 190000, 210000);
  B0s_vs_run_ref_HLT4  = new TH1F ("B0s_vs_run_ref_HLT4", "B0s_vs_run_ref_HLT4;  Run[#];J/Psi[#]",20000, 190000, 210000);
  B0s_vs_run_ref_HLT8  = new TH1F ("B0s_vs_run_ref_HLT8", "B0s_vs_run_ref_HLT8;  Run[#];J/Psi[#]",20000, 190000, 210000);

  JPsi_B0_vs_run       = new TH1F ("JPsi_B0_vs_run", "JPsi_B0_vs_run; Run[#];J/Psi[#]",20000, 190000, 210000);
  JPsi_B0_vs_run_HLT4  = new TH1F ("JPsi_B0_vs_run_HLT4", "JPsi_B0_vs_run_HLT4;  Run[#];J/Psi[#]",20000, 190000, 210000);
  JPsi_B0_vs_run_HLT8  = new TH1F ("JPsi_B0_vs_run_HLT8", "JPsi_B0_vs_run_HLT8;  Run[#];J/Psi[#]",20000, 190000, 210000);

  JPsi_B0_vs_run_ref       = new TH1F ("JPsi_B0_vs_run_ref", "JPsi_B0_vs_run_ref; Run[#];J/Psi[#]",20000, 190000, 210000);
  JPsi_B0_vs_run_ref_HLT4  = new TH1F ("JPsi_B0_vs_run_ref_HLT4", "JPsi_B0_vs_run_ref_HLT4;  Run[#];J/Psi[#]",20000, 190000, 210000);
  JPsi_B0_vs_run_ref_HLT8  = new TH1F ("JPsi_B0_vs_run_ref_HLT8", "JPsi_B0_vs_run_ref_HLT8;  Run[#];J/Psi[#]",20000, 190000, 210000);

  runLumi = new TH2F("runLumi","runLumi",20000, 190000, 210000,3000,0,3000);

}

bool JPsiCount::Process(Long64_t entry)
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
  bool HLT_4_v9 = false, HLT_4_v10 = false, HLT_4_v11 = false, HLT_4_v12 = false ;
  bool HLT_8_v3 = false, HLT_8_v4 = false, HLT_8_v5 = false, HLT_8_v6 = false, HLT_8_v7 = false ;
  bool HLT_4_vAny = false;
  bool HLT_8_vAny = false;
  bool HLT_Any = false;

  std::vector<bool> hltsFlags;

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
  if (HLT_8_vAny || HLT_4_vAny) HLT_Any = true;

  hltsFlags.push_back(true);
  hltsFlags.push_back(HLT_Any);
  hltsFlags.push_back(HLT_4_vAny);
  hltsFlags.push_back(HLT_8_vAny);
  hltsFlags.push_back(HLT_4_vAny && HLT_8_vAny);

  std::map<int,int> doneJPsi;
  std::map<int,int>::iterator doneJPsiIt;

  int muonQual[] = {1,3,4,12};

  runLumi->Fill(runNum,lumiNum);

  bool trueFalse[2] = {true,false};

  TH1F* JPsi_Masses[2]      = {JPsi_Mass,JPsi_Mass_ref};
  TH1F* JPsi_Masses_HLT8[2] = {JPsi_Mass_HLT8,JPsi_Mass_ref_HLT8};
  TH1F* JPsi_Masses_HLT4[2] = {JPsi_Mass_HLT4,JPsi_Mass_ref_HLT4};

  TH1F* B0_Masses[2]      = {B0_Mass,B0_Mass_ref};
  TH1F* B0_Masses_HLT8[2] = {B0_Mass_HLT8,B0_Mass_HLT8_ref};
  TH1F* B0_Masses_HLT4[2] = {B0_Mass_HLT4,B0_Mass_HLT4_ref};

  TH1F* JPsi_vs_runs[2]      = {JPsi_vs_run,JPsi_vs_run_ref};
  TH1F* JPsi_vs_runs_HLT4[2] = {JPsi_vs_run_HLT4,JPsi_vs_run_ref_HLT4};
  TH1F* JPsi_vs_runs_HLT8[2] = {JPsi_vs_run_HLT8,JPsi_vs_run_ref_HLT8};

  TH1F* JPsi_B0_vs_runs[2]      = {JPsi_B0_vs_run,JPsi_B0_vs_run_ref};
  TH1F* JPsi_B0_vs_runs_HLT4[2] = {JPsi_B0_vs_run_HLT4,JPsi_B0_vs_run_ref_HLT4};
  TH1F* JPsi_B0_vs_runs_HLT8[2] = {JPsi_B0_vs_run_HLT8,JPsi_B0_vs_run_ref_HLT8};

  TH1F* B0s_vs_runs[2]       = {B0s_vs_run,B0s_vs_run_ref};
  TH1F* B0s_vs_runs_HLT4[2]  = {B0s_vs_run_HLT4,B0s_vs_run_ref_HLT4};
  TH1F* B0s_vs_runs_HLT8[2]  = {B0s_vs_run_HLT8,B0s_vs_run_ref_HLT8};

  for(iRef = 0; iRef < 2; ++iRef)
  {
    for(int iX=0; iX<nX; ++iX)
    {
      bool refitMuons = trueFalse[iRef], refitKaons = true;

      std::map<std::string,bool> allCuts;
      std::map<std::string,bool> hltCuts;
      std::map<std::string,bool> lxyCuts;

      std::vector<bool> cutsFlags, winsFlags, regsFlags;

      bool muonQualityCut = false, muonChiCut = false, muonPhitsCut = false,           muonShitsCut = false, muonDZPVCut= false, muonDXYPVCut = false, muonSoftCuts = false;

      bool jPsiPtCut = false,jPsiMassCut = false, jPsiVtxCut = false, jPsiMuEtaPtCut = false,           jPsiMusPtCut = false, jPsiMu1SelCut  = false, jPsiMu2SelCut  = false, jPsiCuts = false;

      bool kaonOneChiCut = false, kaonOnePhitsCut = false, kaonOneShitsCut = false, kaonTwoChiCut = false, kaonTwoPhitsCut = false, kaonTwoShitsCut = false;
      bool phiCut = false;

      bool kaonOneCuts = false, kaonTwoCuts = false, kaonsCuts = false;
      bool cosAlphaCut = false, vtxCLCut = false;

      bool CWMass = false, SWMass = false;
      bool promptRegion = false, mixedRegion = false, nonPromptRegion = false;

      bool x5568_muOnePt = false, x5568_muOneEta = false, x5568_muTwoPt = false, x5568_muTwoEta = false;
      bool x5568_mumuPt = false, x5568_mumuVtx = false, x5568_Lxy = false, x5568_cosAlpha = false;
      bool x5568_mumuMass = false, x5568_muDxy = false ,x5568_muDz = false;
      bool x5568_kkMass = false, x5568_kkPt = false;
      bool x5568_B0smass = false, x5568_B0sVtx = false, x5568_B0scosAlpha = false;

      bool x5568_B0sCut = false, x5568_mumuCut = false, x5568_muCut = false, x5568_kkCut = false;

      int iJPsi = (*XMuMuIdx)[iX];

      doneJPsiIt = doneJPsi.find(iJPsi);

      // if(doneJPsiIt!=doneJPsi.end())
      //   continue;
      // else
      //   doneJPsi[iJPsi] = 1.0;

      int iMu1 = (*mu1Idx)[iJPsi] ; // define for original muon1
      int iMu2 = (*mu2Idx)[iJPsi] ; // define for original muon2
      int iK1 = (*ka1Idx)[iX] ; // define for original kaon1
      int iK2 = (*ka2Idx)[iX] ;

      double mu1_E = 0., mu2_E = 0., K1_E = 0., K2_E = 0.;

      TLorentzVector mu1, mu2;

      mu1_E = sqrt( pow((*muPx)[iMu1], 2) + pow((*muPy)[iMu1], 2) + pow((*muPz)[iMu1], 2) + pow(muon_mass, 2) ) ;
      mu2_E = sqrt( pow((*muPx)[iMu2], 2) + pow((*muPy)[iMu2], 2) + pow((*muPz)[iMu2], 2) + pow(muon_mass, 2) ) ;

      if(refitMuons)
      {
        mu1.SetPxPyPzE((*Muon1Px_MuMuKK)[iX],(*Muon1Py_MuMuKK)[iX],(*Muon1Pz_MuMuKK)[iX],(*Muon1E_MuMuKK)[iX]);
        mu2.SetPxPyPzE((*Muon2Px_MuMuKK)[iX],(*Muon2Py_MuMuKK)[iX],(*Muon2Pz_MuMuKK)[iX],(*Muon2E_MuMuKK)[iX]);
      }
      else
      {
        mu1.SetPxPyPzE( (*muPx)[iMu1], (*muPy)[iMu1], (*muPz)[iMu1], mu1_E) ;
        mu2.SetPxPyPzE( (*muPx)[iMu2], (*muPy)[iMu2], (*muPz)[iMu2], mu2_E) ;
      }

      TLorentzVector JPsi;
      JPsi = mu1 + mu2;

      TLorentzVector kaon1, kaon2;

      if(refitKaons)
      {
        K1_E=sqrt(pow((*Kaon1Px_MuMuKK)[iX],2)+pow((*Kaon1Py_MuMuKK)[iX],2)+pow((*Kaon1Pz_MuMuKK)[iX],2)+pow(kaonCh_mass,2));
        kaon1.SetPxPyPzE((*Kaon1Px_MuMuKK)[iX],(*Kaon1Py_MuMuKK)[iX],(*Kaon1Pz_MuMuKK)[iX],K1_E);
        K2_E=sqrt(pow((*Kaon2Px_MuMuKK)[iX],2)+pow((*Kaon2Py_MuMuKK)[iX],2)+pow((*Kaon2Pz_MuMuKK)[iX],2)+pow(kaonCh_mass,2));
        kaon2.SetPxPyPzE((*Kaon2Px_MuMuKK)[iX],(*Kaon2Py_MuMuKK)[iX],(*Kaon2Pz_MuMuKK)[iX],K2_E);
      }
      else
      {
        K1_E = sqrt( pow((*trackPx)[iK1],2)+pow((*trackPy)[iK1],2)+pow((*trackPz)[iK1],2)+pow(kaonCh_mass,2)) ;
        kaon1.SetPxPyPzE((*trackPx)[iK1],(*trackPy)[iK1],(*trackPz)[iK1], K1_E) ;
        K2_E = sqrt( pow((*trackPx)[iK2],2)+pow((*trackPy)[iK2],2)+pow((*trackPz)[iK2],2)+pow(kaonCh_mass,2)) ;
        kaon2.SetPxPyPzE((*trackPx)[iK2],(*trackPy)[iK2],(*trackPz)[iK2], K2_E);
      }

      TLorentzVector Phi;
      Phi = kaon1 + kaon2;

      Muon1_Mass->Fill(mu1.M());
      Muon2_Mass->Fill(mu2.M());

      TLorentzVector XCand;
      XCand = JPsi + Phi;

      // SWMass = (((XCand.M() > 4.05) && (XCand.M() < 4.8)));
      // CWMass = ((XCand.M() > 5.2) && (XCand.M() < 5.55));

      x5568_muOnePt   = (mu1.Pt() > 4.0);
      x5568_muOneEta  = (fabs(mu1.Eta()) < 2.2);
      x5568_muTwoPt   = (mu2.Pt() > 4.0);
      x5568_muTwoEta  = (fabs(mu2.Eta()) < 2.2);
      x5568_muDz = (fabs((*muDzVtx)[iMu1]) < 20.0 && fabs((*muDzVtx)[iMu2]) < 20.0);
      x5568_muDxy = (fabs((*muDxyVtx)[iMu1]) < 0.3 && fabs((*muDxyVtx)[iMu2]) < 0.3);

      x5568_muCut = (x5568_muOnePt && x5568_muTwoPt && x5568_muTwoEta && x5568_muOneEta && x5568_muDxy && x5568_muDz);
      //x5568_Lxy       = (
      //5568_cosAlpha = false;
      x5568_mumuPt    = (JPsi.Pt() > 7.0);
      x5568_mumuVtx   = ((*MuMuVtx_CL)[iJPsi]) > 0.1;
      x5568_mumuMass = (JPsi.M()<3.15 && JPsi.M()>3.04);

      x5568_mumuCut = (x5568_mumuPt && x5568_mumuVtx && x5568_mumuMass);

      x5568_kkMass = (fabs(Phi.M()-Phi_mass)<0.01);
      x5568_kkPt = ((kaon1.Pt()>0.7) && (kaon2.Pt()>0.7));

      x5568_kkCut = (x5568_kkMass && x5568_kkPt);

      x5568_B0smass = ((XCand.M() > 5.2) && (XCand.M() < 5.55));
      x5568_B0sVtx = (((*XVtx_CL)[iX]) > 0.01);
      x5568_B0scosAlpha = (fabs((*XCosAlphaPV)[iX]) > 0.99);

      x5568_B0sCut = (x5568_B0smass && x5568_B0sVtx && x5568_B0scosAlpha);

      if(x5568_mumuMass && x5568_B0smass)
      {

        if(doneJPsiIt==doneJPsi.end())
        {
          JPsi_vs_runs[iRef]->Fill(runNum);
          JPsi_Masses[iRef]->Fill(JPsi.M());

          if(HLT_8_vAny && ((*XLxyPV)[iX]/(*XLxyPVE)[iX]>3.) && JPsi.Pt()>8.0)
            {
              JPsi_vs_runs_HLT8[iRef]->Fill(runNum);
              JPsi_Masses_HLT8[iRef]->Fill(JPsi.M());
            }

          if(HLT_4_vAny && JPsi.Pt()>8.0 && ((*XLxyPV)[iX]/(*XLxyPVE)[iX]>3.))
          {
              JPsi_vs_runs_HLT4[iRef]->Fill(runNum);
              JPsi_Masses_HLT4[iRef]->Fill(JPsi.M());
          }

        }

        if(x5568_muCut && x5568_kkCut && x5568_mumuCut && x5568_B0sCut)
        {
          if(doneJPsiIt==doneJPsi.end())
          {
            JPsi_B0_vs_runs[iRef]->Fill(runNum);

            if(HLT_8_vAny && ((*XLxyPV)[iX]/(*XLxyPVE)[iX]>3.) && JPsi.Pt()>8.0)
              {
                JPsi_B0_vs_runs_HLT8[iRef]->Fill(runNum);
              }

            if(HLT_4_vAny && JPsi.Pt()>8.0 && ((*XLxyPV)[iX]/(*XLxyPVE)[iX]>3.))
            {
                JPsi_B0_vs_runs_HLT4[iRef]->Fill(runNum);
            }

          }

          B0s_vs_runs[iRef]->Fill(runNum);

          if(HLT_8_vAny && ((*XLxyPV)[iX]/(*XLxyPVE)[iX]>3.) && JPsi.Pt()>8.0)
            {
              B0s_vs_runs_HLT8[iRef]->Fill(runNum);
              B0_Masses_HLT8[iRef]->Fill(XCand.M());
            }

          if(HLT_4_vAny && JPsi.Pt()>8.0 && ((*XLxyPV)[iX]/(*XLxyPVE)[iX]>3.))
          {
              B0s_vs_runs_HLT4[iRef]->Fill(runNum);
              B0_Masses_HLT4[iRef]->Fill(XCand.M());
          }

          B0_Masses[iRef]->Fill(XCand.M());

        }


      }

    }
  }
  return kTRUE;
}


void JPsiCount::SlaveTerminate()
{

  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
  TDirectory *savedir = gDirectory;
  if (fOut)
  // if(false)
  {
    fOut->cd();
    gStyle->SetOptStat(111111) ;

    runLumi->Write();
    ////////////////// Write Histograms //////////////////
    /// refit muons & kaons
    JPsi_vs_run->Write();
    JPsi_vs_run_HLT8->Write();
    JPsi_vs_run_HLT4->Write();

    JPsi_vs_run_ref->Write();
    JPsi_vs_run_ref_HLT8->Write();
    JPsi_vs_run_ref_HLT4->Write();

    JPsi_B0_vs_run->Write();
    JPsi_B0_vs_run_HLT8->Write();
    JPsi_B0_vs_run_HLT4->Write();

    B0s_vs_run->Write();
    B0s_vs_run_HLT8->Write();
    B0s_vs_run_HLT4->Write();

    JPsi_Mass->Write();
    JPsi_Mass_HLT8->Write();
    JPsi_Mass_HLT4->Write();

    B0_Mass->Write();
    B0_Mass_HLT8->Write();
    B0_Mass_HLT4->Write();

    OutFile->Print();
    fOutput->Add(OutFile);
    gDirectory = savedir;
    fOut->Close();

  }


}


void JPsiCount::Terminate()
{


}
