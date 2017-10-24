#define X4140_MuMuKK_cxx
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

#include "X4140_MuMuKK.h"
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

void X4140_MuMuKK::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

}


void X4140_MuMuKK::SlaveBegin(TTree * /*tree*/)
{
  //std::cout << "SLAVE BEGIN" << std::endl;

  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

  std::string outputString = "X4140_MuMuKK_KRe_MuRef_EtaX1.1_NP3.5_Alpha99_CW5.2-5.55.root";
  OutFile = new TProofOutputFile( outputString.data() );
  fOut = OutFile->OpenFile("RECREATE");
  if (!(fOut=OutFile->OpenFile("RECREATE")))
  {
    Warning("SlaveBegin","Problems opening file: %s%s", OutFile->GetDir(), OutFile->GetFileName() );
  }

  ////////////////// Histograms //////////////////
  JPsi_mass = 3.096916; /// pdg mass
  Phi_mass = 1.019455; /// pdg mass

  /// refit muons & kaons
  // Muon1_Mass     = new TH1F ("Muon1_mass", "Muon1_mass; m(#mu^{+}) [GeV];Entries",1000, 0.0, 2.0);
  // Muon2_Mass     = new TH1F ("Muon2_mass", "Muon2_mass; m(#mu^{-}) [GeV];Entries",1000, 0.0, 2.0);

  // MuMu_Mass      = new TH1F ("MuMu_Mass_NoHLT", "MuMu_Mass; M(#mu^{-}#mu^{+}) [GeV];Entries",200, 2.5, 3.5);
  // MuMu_Mass_HLT8 = new TH1F ("MuMu_Mass_HLT8", "MuMu_Mass_HLT8; M(#mu^{-}#mu^{+}) [GeV];Entries",200, 2.5, 3.5);
  // MuMu_Mass_HLT4 = new TH1F ("MuMu_Mass_HLT4", "MuMu_Mass_HLT4; M(#mu^{-}#mu^{+}) [GeV];Entries",200, 2.5, 3.5);

  // XCand_mass = new TH1F ("XCand_mass","XCand_mass",400,4.0,6.0);
  // Y4140_mass = new TH1F ("Y4140_mass","Y4140_mass",400,4.0,6.0);
  // Bs0_mass = new TH1F ("Bs0_mass","Bs0_mass",400,4.0,6.0);
  // Lxy_LxyE_PV = new TH1F ("Lxy_LxyE_PV","Lxy_LxyE_PV",500,0.0,5.0);

  std::string xcandHisto = "Xcand_histo_";
  std::string jpsiHisto  = "JPsi_histo_";
  std::string phiHisto   = "Phi_hist_";


  for (size_t l = 0; l < 5; l++)
  {
    for (size_t j = 0; j < 3; ++j)
    {
      for (size_t k = 0; k < 4; ++k)
      {
        for (size_t i = 0; i < 6; ++i)
        {
          std::tuple < std::string, std::string, std::string, std::string > histokey(hlts[l],windows[j],regions[k],clabels[i]);
          std::string suffix = hlts[l] + "_" + windows[j] + "_" + regions[k] + "_" + clabels[i];

          XCandMassHistos[histokey] = new TH1F ((xcandHisto + suffix).data(),(xcandHisto + suffix).data(),2000, 4.0, 6.0);
          JPsiMassHistos[histokey]  = new TH1F ((jpsiHisto + suffix).data(),(jpsiHisto + suffix).data(),1000, 2.5, 3.5);
          PhiMassHistos[histokey]   = new TH1F ((phiHisto + suffix).data(),(phiHisto + suffix).data(),1000, 0.5, 1.5);

        }
      }
    }
  }

  for (size_t l = 0; l < 5; l++)
  {
    for (size_t j = 0; j < 3; ++j)
    {
      for (size_t k = 0; k < 4; ++k)
      {
        for (size_t i = 0; i < 6; ++i)
        {
          std::tuple < std::string, std::string, std::string, std::string > histokey(hlts[l],windows[j],regions[k],clabels[i]);
          std::string suffix = hlts[l] + "_" + windows[j] + "_" + regions[k] + "_" + clabels[i] + "_ref";

          XCandMassHistos_Ref[histokey] = new TH1F ((xcandHisto + suffix).data(),(xcandHisto + suffix).data(),2000, 4.0, 6.0);
          JPsiMassHistos_Ref[histokey]  = new TH1F ((jpsiHisto + suffix).data(),(jpsiHisto + suffix).data(),1000, 2.5, 3.5);
          PhiMassHistos_Ref[histokey]   = new TH1F ((phiHisto + suffix).data(),(phiHisto + suffix).data(),1000, 0.5, 1.5);

        }
      }
    }
  }

  X5568_Cand_Mass     = new TH1F ("X5568_Cand_Mass","X5568_Cand_Mass",500, 5.1, 5.6);
  X5568_Cand_Mass_Ref = new TH1F ("X5568_Cand_Mass_Ref","X5568_Cand_Mass_Ref",500, 5.1, 5.6);
  X5568_Cand_Mass_NoM = new TH1F ("X5568_Cand_Mass_Ref_NoMult","X5568_Cand_Mass_Ref_NoMult",500, 5.1, 5.6);
  // runLumi        = new TH2F("runLumi","runLumi",20000, 190000, 210000,3000,0,3000);
  // runLumiB0s     = new TH2F("runLumiB0s","runLumiB0s",20000, 190000, 210000,3000,0,3000);
  // runLumiB0sJPsi = new TH2F("runLumiB0sJPsi","runLumiB0sJPsi",20000, 190000, 210000,3000,0,3000);
}

bool X4140_MuMuKK::Process(Long64_t entry)
{
  //std::cout << "PROCESS" << std::endl;

  GetEntry(entry); /// ROOT event loop

  // The entry argument specifies which entry in the currently loaded tree is to be processed.
  // It can be passed to either Y4140::GetEntry() or TBranch::GetEntry()
  // to read either all or the required parts of the data. When processing

  double muon_mass = 0.1056583715;
  double kaonCh_mass = 0.493677;
  double Bs0_Low_Mass = 5.1;
  double Bs0_High_Mass = 5.6;
  double Y_High_Mass = 4.35;


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

  int X4140_MuMuKKer = 0,jPsiMuMuCounter = 0, jPsis = 0;
  int muonQual[] = {1,3,4,12};

  // for (std::map < std::pair <int,int>, int>::iterator lumiMapIt=lumiMap.begin(); lumiMapIt!=lumiMap.end(); ++lumiMapIt)
  //  runlumi << lumiMapIt->first.first << " : " << lumiMapIt->first.second << std::endl;


  bool trueFalse[2] = {false,true};

  TH1F* x5586_hists[2] = {X5568_Cand_Mass,X5568_Cand_Mass_Ref};
  std::map< std::tuple < std::string, std::string, std::string, std::string >,TH1F*> xCand_hists[2] = {XCandMassHistos,XCandMassHistos_Ref};
  std::map< std::tuple < std::string, std::string, std::string, std::string >,TH1F*> jPsis_hists[2] = {JPsiMassHistos,JPsiMassHistos_Ref};
  std::map< std::tuple < std::string, std::string, std::string, std::string >,TH1F*> phi_hists[2]   = {PhiMassHistos,PhiMassHistos_Ref};

  // runLumi->Fill(runNum,lumiNum);
  std::vector<double> x5568CandsB0sMasses;

  for(int iRef=1; iRef<2; ++iRef)
    {

      bool refitMuons = trueFalse[iRef], refitKaons = true;
      x5568CandsB0sMasses.clear();
      for(int iX=0; iX<nX; ++iX)
      {

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

        bool extraCuts = false;

        int iJPsi = (*XMuMuIdx)[iX];

        doneJPsiIt = doneJPsi.find(iJPsi);

        //if(doneJPsiIt!=doneJPsi.end())
          //continue;
        //else
          //doneJPsi[iJPsi] = 1.0;

        ++jPsis;

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

        // Muon1_Mass->Fill(mu1.M());
        // Muon2_Mass->Fill(mu2.M());

        TLorentzVector XCand;
        XCand = JPsi + Phi;

        SWMass = (((XCand.M() > 4.05) && (XCand.M() < 4.8)));
        CWMass = ((XCand.M() > 5.2) && (XCand.M() < 5.55));

        mixedRegion     = (((*XLxyPV)[iX] / (*XLxyPVE)[iX]) >= 2. && ((*XLxyPV)[iX] / (*XLxyPVE)[iX])  <= 3.5);
        nonPromptRegion = (((*XLxyPV)[iX] / (*XLxyPVE)[iX]) > 3.5);
        promptRegion    = (((*XLxyPV)[iX] / (*XLxyPVE)[iX]) < 2.);

        winsFlags.push_back(true);
        winsFlags.push_back(CWMass);
        winsFlags.push_back(SWMass);

        regsFlags.push_back(true);
        regsFlags.push_back(promptRegion);
        regsFlags.push_back(mixedRegion);
        regsFlags.push_back(nonPromptRegion);

        // XCand_mass->Fill(XCand.M());

        // if ((XCand.M() > 4.05) && (XCand.M() < 4.8))
        //   Y4140_mass->Fill(XCand.M());
        //
        // else if ((XCand.M() > 5.15) && (XCand.M() < 5.55))
        //   Bs0_mass->Fill(XCand.M());

        cutsFlags.push_back(true);

        //x5568 cuts

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

        if(x5568_muCut && x5568_kkCut && x5568_mumuCut && x5568_B0sCut && HLT_4_vAny)
          {
            x5586_hists[iRef]->Fill(XCand.M());
            x5568CandsB0sMasses.push_back(XCand.M());
          }

        ////////////////////////////////////
        //Muon Cuts

        muonQualityCut = ( ((*muQual)[iMu1]) & (1 << muonQual[3]) ) && ( ((*muQual)[iMu2]) & (1 << muonQual[3]) );
        muonChiCut     = (( ( (*muChi2)[iMu1] / (*muNDF)[iMu1] ) < 3 ) && ( ( (*muChi2)[iMu2] / (*muNDF)[iMu2] ) < 3 ));
        muonPhitsCut   = ((*muPhits)[iMu1] > 0 && (*muPhits)[iMu2] > 0);
        muonShitsCut   = ((*muShits)[iMu1] > 5 && (*muShits)[iMu2] > 5);
        muonDZPVCut    = (fabs((*muDzVtx)[iMu1]) < 20.0 && fabs((*muDzVtx)[iMu2]) < 20.0);
        muonDXYPVCut   = (fabs((*muDxyVtx)[iMu1]) < 3.0 && fabs((*muDxyVtx)[iMu2]) < 3.0);

        muonSoftCuts = (muonQualityCut && muonChiCut && muonPhitsCut && muonShitsCut && muonDXYPVCut && muonDXYPVCut);

        cutsFlags.push_back((muonSoftCuts));

        ////////////////////////////////////
        //JPsi cuts

        jPsiPtCut      = (JPsi.Pt() > 8.0);
        jPsiMassCut    = (fabs(JPsi.M() - JPsi_mass) < 0.1); //
        jPsiVtxCut     = (((*MuMuVtx_CL)[iJPsi]) > 0.01);
        jPsiMuEtaPtCut = (fabs(mu1.Eta()) < 2.4 && fabs(mu2.Eta()) < 2.4);
        jPsiMusPtCut   = (mu1.Pt() > 3.3 && mu2.Pt() > 3.3);
        jPsiMu1SelCut  = ((((fabs(mu1.Eta()) < 1.2) && (mu1.Pt() > 4.))) || ((mu1.Eta() >= 1.2 || mu1.Eta() <= -1.2 ) && (mu1.Pt() > 3.3)));
        jPsiMu2SelCut  = ((((fabs(mu2.Eta()) < 1.2) && (mu2.Pt() > 4.))) || ((mu2.Eta() >= 1.2 || mu2.Eta() <= -1.2 ) && (mu2.Pt() > 3.3)));

        jPsiCuts = jPsiPtCut && jPsiMassCut && jPsiVtxCut  && jPsiMuEtaPtCut && jPsiMuEtaPtCut && jPsiMusPtCut && jPsiMu1SelCut && jPsiMu2SelCut;

        cutsFlags.push_back((jPsiCuts));

        ////////////////////////////////////
        //Kaon Track Cuts
        kaonOneChiCut    = (((*trackChi2)[iK1] / (*trackNDF)[iK1]) < 5.0);
        kaonOnePhitsCut  = ((*trackPhits)[iK1] > 0);
        kaonOneShitsCut  = ((*trackShits)[iK1] >= 7);
        kaonTwoChiCut    = (((*trackChi2)[iK2] / (*trackNDF)[iK2]) < 5.0);
        kaonTwoPhitsCut  = ((*trackPhits)[iK2] > 0);
        kaonTwoShitsCut  = ((*trackShits)[iK2] >= 7);

        //Phi cut
        phiCut = (fabs(Phi.M()-Phi_mass)<0.007);
        //((Phi.M()) > 1.013 && (Phi.M()) < 1.026);

        kaonOneCuts = kaonOneChiCut && kaonOnePhitsCut && kaonOneShitsCut;
        kaonTwoCuts = kaonTwoChiCut && kaonTwoPhitsCut && kaonTwoShitsCut;
        kaonsCuts = kaonOneCuts && kaonTwoCuts;

        cutsFlags.push_back((kaonsCuts));
        cutsFlags.push_back((phiCut));

        ////////////////////////////////////
        //Vtx and cos(alpha) cuts
        cosAlphaCut = (fabs((*XCosAlphaPV)[iX]) > 0.99);
        vtxCLCut =  (((*XVtx_CL)[iX]) > 0.01);

        ////Extra cuts
        bool etaXCut = fabs(XCand.Eta()) < 1.1;

        extraCuts = etaXCut;

        cutsFlags.push_back((cosAlphaCut && vtxCLCut && extraCuts));

        std::string window,region,hlt,cut;

        for (size_t l = 0; l < hltsFlags.size(); l++)
        {
          if(hltsFlags[l])
          {
            hlt = hlts[l];
            for (size_t j = 0; j < winsFlags.size(); ++j)
            {
              if(winsFlags[j])
              {
                window = windows[j];
                for (size_t k = 0; k < regsFlags.size(); ++k)
                {
                  if(regsFlags[k])
                  {
                    region = regions[k];
                    bool cut = true;
                    for (size_t i = 0; i < cutsFlags.size(); ++i)
                    {
                      cut = cutsFlags[i] && cut;
                      if(cut)
                      {
                        std::tuple < std::string, std::string, std::string, std::string > histokey(hlts[l],windows[j],regions[k],clabels[i]);
                        xCand_hists[iRef][histokey]->Fill(XCand.M());
                        jPsis_hists[iRef][histokey]->Fill(JPsi.M());
                        phi_hists[iRef][histokey]->Fill(Phi.M());
                      }
                    }

                  }
                }
              }
            }
          }
        }


        cutsFlags.clear();
        winsFlags.clear();
        regsFlags.clear();
      }


  }

  if(x5568CandsB0sMasses.size()==1)
    X5568_Cand_Mass_NoM->Fill(x5568CandsB0sMasses[0]);

  return kTRUE;
}


void X4140_MuMuKK::SlaveTerminate()
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

    for (size_t l = 0; l < 5; l++)
    {
      for (size_t j = 0; j < 3; ++j)
      {
        for (size_t k = 0; k < 4; ++k)
        {
          for (size_t i = 0; i < 6; ++i)
          {
            std::tuple < std::string, std::string, std::string, std::string > histokey(hlts[l],windows[j],regions[k],clabels[i]);
            std::string suffix = hlts[l] + "_" + windows[j] + "_" + regions[k] + "_" + clabels[i];

            XCandMassHistos[histokey]->Write();
            JPsiMassHistos[histokey]->Write();
            PhiMassHistos[histokey]->Write();
            XCandMassHistos_Ref[histokey]->Write();
            JPsiMassHistos_Ref[histokey]->Write();
            PhiMassHistos_Ref[histokey]->Write();

          }
        }
      }
    }

    X5568_Cand_Mass->Write();
    X5568_Cand_Mass_Ref->Write();
    X5568_Cand_Mass_NoM->Write();

    OutFile->Print();
    fOutput->Add(OutFile);
    gDirectory = savedir;
    fOut->Close();

  }


}


void X4140_MuMuKK::Terminate()
{


}
