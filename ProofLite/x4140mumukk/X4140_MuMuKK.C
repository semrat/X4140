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

  std::string outputString = "X4140_MuMuKK_KRe_MuRef_Sidebands5-7_NP3.0_B0Cuts_CW5.15-5.55.root";
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


  for (size_t l = 0; l < 4; l++)
  {
    for (size_t j = 0; j < 3; ++j)
    {
      for (size_t k = 0; k < 4; ++k)
      {
        for (size_t i = 0; i < 5; ++i)
        {
          std::tuple < std::string, std::string, std::string, std::string > histokey(hlts[l],windows[j],regions[k],clabels[i]);
          std::string suffix = hlts[l] + "_" + windows[j] + "_" + regions[k] + "_" + clabels[i];

          XCandMassHistos[histokey] = new TH1F ((xcandHisto + suffix).data(),(xcandHisto + suffix).data(),2000, 4.0, 6.0);
          JPsiMassHistos[histokey]  = new TH1F ((jpsiHisto + suffix).data(),(jpsiHisto + suffix).data(),1000, 2.5, 3.5);
          PhiMassHistos[histokey]   = new TH1F ((phiHisto + suffix).data(),(phiHisto + suffix).data(),1000, 0.5, 1.5);
          //XCandMassHistosDeltaM[histokey] = new TH1F ((xcandHistoDM + suffix).data(),(xcandHistoDM + suffix).data(),2000, 4.0, 6.0);
        }
      }
    }
  }


  B0_Cand_Mass           = new TH1F ("B0_Cand_Mass","B0_Cand_Mass",500, 5.1, 5.6);
  B0_Cand_Mass_L_Side = new TH1F ("B0_Cand_Mass_L_Side","B0_Cand_Mass_L_Side",500, 5.1, 5.6);
  B0_Cand_Mass_R_Side = new TH1F ("B0_Cand_Mass_R_Side","B0_Cand_Mass_R_Side",500, 5.1, 5.6);

  B0_Cand_Mass_NoM       = new TH1F ("B0_Cand_Mass_NoM","B0_Cand_Mass_NoM",500, 5.1, 5.6);
  B0_Cand_Mass_L_Side_NoM = new TH1F ("B0_Cand_Mass_L_Side_NoM","B0_Cand_Mass_L_Side_NoM",500, 5.1, 5.6);
  B0_Cand_Mass_R_Side_NoM = new TH1F ("B0_Cand_Mass_R_Side_NoM","B0_Cand_Mass_R_Side_NoM",500, 5.1, 5.6);

  B0_PhiMassHisto        = new TH1F ("B0_PhiMassHisto","B0_PhiMassHisto",1000, 0.5, 1.5);
  B0_PhiMassHisto_NoM    = new TH1F ("B0_PhiMassHisto_NoM","B0_PhiMassHisto_NoM",1000, 0.5, 1.5);

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
  //hltsFlags.push_back(HLT_4_vAny && HLT_8_vAny);

  std::map<int,int> doneJPsi;
  std::map<int,int>::iterator doneJPsiIt;

  int X4140_MuMuKKer = 0,jPsiMuMuCounter = 0, jPsis = 0;
  int muonQual[] = {1,3,4,12};

  // for (std::map < std::pair <int,int>, int>::iterator lumiMapIt=lumiMap.begin(); lumiMapIt!=lumiMap.end(); ++lumiMapIt)
  //  runlumi << lumiMapIt->first.first << " : " << lumiMapIt->first.second << std::endl;


  // runLumi->Fill(runNum,lumiNum);
  std::vector<double> xCandsB0sMasses,B0_Cand_Masses_L_Sides,B0_Cand_Masses_R_Sides,B0_Cand_Masses,B0_phi_masses;
  B0_phi_masses.clear();
  xCandsB0sMasses.clear();
  B0_Cand_Masses_R_Sides.clear();
  B0_Cand_Masses_L_Sides.clear();
  B0_Cand_Masses.clear();

  for(int iX=0; iX<nX; ++iX)
  {

    std::map<std::string,bool> allCuts;
    std::map<std::string,bool> hltCuts;
    std::map<std::string,bool> lxyCuts;

    std::vector<bool> cutsFlags, winsFlags, regsFlags;

    bool muonQualityCut = false, muonChiCut = false, muonPhitsCut = false, muonShitsCut = false;
    bool muonDZPVCut= false, muonDXYPVCut = false, muonSoftCuts = false, muonsCuts = false;

    bool jPsiPtCut = false,jPsiMassCut = false, jPsiVtxCut = false, jPsiMuEtaPtCut = false, jPsiMusPtCut = false, jPsiCuts = false;

    bool kaonOneChiCut = false, kaonOnePhitsCut = false, kaonOneShitsCut = false, kaonTwoChiCut = false;
    bool kaonTwoPhitsCut = false, kaonTwoShitsCut = false, kaonsPt = false;

    bool kaonOneCuts = false, kaonTwoCuts = false, kaonsCuts = false;
    bool cosAlphaCut = false, vtxCLCut = false;

    bool CWMass = false, SWMass = false;
    bool promptRegion = false, mixedRegion = false, nonPromptRegion = false;

    bool new_muOnePt = false, new_muOneEta = false, new_muTwoPt = false, new_muTwoEta = false;
    bool new_mumuPt = false, new_mumuVtx = false, new_Lxy = false, new_cosAlpha = false;
    bool new_mumuMass = false, new_muDxy = false ,new_muDz = false;
    bool new_kkMass = false, new_kkPt = false;
    bool new_B0smass = false, new_B0sVtx = false, new_B0scosAlpha = false;

    bool new_B0sCut = false, new_mumuCut = false, new_muCut = false, new_kkCut = false;

    bool extraCuts = false;

    bool b0_side = false, b0_signal = false, b0_side_l = false, b0_side_r = false;

    int iJPsi = (*XMuMuIdx)[iX];

    //doneJPsiIt = doneJPsi.find(iJPsi);

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

    TLorentzVector mu1, mu2, refMu1, refMu2;

    mu1_E = sqrt( pow((*muPx)[iMu1], 2) + pow((*muPy)[iMu1], 2) + pow((*muPz)[iMu1], 2) + pow(muon_mass, 2) ) ;
    mu2_E = sqrt( pow((*muPx)[iMu2], 2) + pow((*muPy)[iMu2], 2) + pow((*muPz)[iMu2], 2) + pow(muon_mass, 2) ) ;
    refMu1.SetPxPyPzE((*Muon1Px_MuMuKK)[iX],(*Muon1Py_MuMuKK)[iX],(*Muon1Pz_MuMuKK)[iX],(*Muon1E_MuMuKK)[iX]);
    refMu2.SetPxPyPzE((*Muon2Px_MuMuKK)[iX],(*Muon2Py_MuMuKK)[iX],(*Muon2Pz_MuMuKK)[iX],(*Muon2E_MuMuKK)[iX]);
    mu1.SetPxPyPzE( (*muPx)[iMu1], (*muPy)[iMu1], (*muPz)[iMu1], mu1_E) ;
    mu2.SetPxPyPzE( (*muPx)[iMu2], (*muPy)[iMu2], (*muPz)[iMu2], mu2_E) ;

    TLorentzVector JPsi;
    JPsi = refMu1 + refMu2;

    TLorentzVector JPsiOriginal;
    JPsiOriginal = mu1 + mu2;

    TLorentzVector kaon1,kaon2;

    K1_E=sqrt(pow((*Kaon1Px_MuMuKK)[iX],2)+pow((*Kaon1Py_MuMuKK)[iX],2)+pow((*Kaon1Pz_MuMuKK)[iX],2)+pow(kaonCh_mass,2));
    kaon1.SetPxPyPzE((*Kaon1Px_MuMuKK)[iX],(*Kaon1Py_MuMuKK)[iX],(*Kaon1Pz_MuMuKK)[iX],K1_E);
    K2_E=sqrt(pow((*Kaon2Px_MuMuKK)[iX],2)+pow((*Kaon2Py_MuMuKK)[iX],2)+pow((*Kaon2Pz_MuMuKK)[iX],2)+pow(kaonCh_mass,2));
    kaon2.SetPxPyPzE((*Kaon2Px_MuMuKK)[iX],(*Kaon2Py_MuMuKK)[iX],(*Kaon2Pz_MuMuKK)[iX],K2_E);

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

    // muonQualityCut = ( ((*muQual)[iMu1]) & (1 << muonQual[3]) ) && ( ((*muQual)[iMu2]) & (1 << muonQual[3]) );
    // muonChiCut     = (( ( (*muChi2)[iMu1] / (*muNDF)[iMu1] ) < 3 ) && ( ( (*muChi2)[iMu2] / (*muNDF)[iMu2] ) < 3 ));
    // muonPhitsCut   = ((*muPhits)[iMu1] > 0 && (*muPhits)[iMu2] > 0);
    // muonShitsCut   = ((*muShits)[iMu1] > 5 && (*muShits)[iMu2] > 5);
    muonDZPVCut    = (fabs((*muDzVtx)[iMu1]) < 20.0 && fabs((*muDzVtx)[iMu2]) < 20.0);
    muonDXYPVCut   = (fabs((*muDxyVtx)[iMu1]) < 0.3 && fabs((*muDxyVtx)[iMu2]) < 0.3);

    //muonsCuts = muonQualityCut && muonChiCut && muonShitsCut && muonPhitsCut && muonDZPVCut && muonDXYPVCut;
    muonsCuts = muonDZPVCut && muonDXYPVCut;

    cutsFlags.push_back((muonsCuts));

    jPsiPtCut      = (JPsi.Pt() > 7.0);
    jPsiMassCut    =  (JPsi.M()<3.15 && JPsi.M()>3.04);
    jPsiVtxCut     = ((*MuMuVtx_CL)[iJPsi]) > 0.1;
    jPsiMuEtaPtCut = (fabs(mu1.Eta()) < 2.2) && (fabs(mu2.Eta()) < 2.2);
    jPsiMusPtCut   = ((mu1.Pt() > 4.0) && (mu2.Pt() > 4.0));

    jPsiCuts = jPsiPtCut && jPsiMassCut && jPsiVtxCut  && jPsiMuEtaPtCut && jPsiMuEtaPtCut && jPsiMusPtCut;

    cutsFlags.push_back((jPsiCuts));

    // kaonOneChiCut    = (((*trackChi2)[iK1] / (*trackNDF)[iK1]) < 5.0);
    // kaonOnePhitsCut  = ((*trackPhits)[iK1] > 0);
    // kaonOneShitsCut  = ((*trackShits)[iK1] >= 7);
    // kaonTwoChiCut    = (((*trackChi2)[iK2] / (*trackNDF)[iK2]) < 5.0);
    // kaonTwoPhitsCut  = ((*trackPhits)[iK2] > 0);
    // kaonTwoShitsCut  = ((*trackShits)[iK2] >= 7);
    kaonsPt = ((kaon1.Pt()>0.7) && (kaon2.Pt()>0.7));

    kaonOneCuts = kaonOneChiCut && kaonOnePhitsCut && kaonOneShitsCut;
    kaonTwoCuts = kaonTwoChiCut && kaonTwoPhitsCut && kaonTwoShitsCut;
    kaonsCuts = kaonOneCuts && kaonTwoCuts && kaonsPt;
    kaonsCuts = kaonsPt;

    cutsFlags.push_back((kaonsCuts));

    cosAlphaCut = (fabs((*XCosAlphaPV)[iX]) > 0.99);
    vtxCLCut =  (((*XVtx_CL)[iX]) > 0.01);

    extraCuts = vtxCLCut && cosAlphaCut;

    cutsFlags.push_back((extraCuts));

    mixedRegion     = (((*XLxyPV)[iX] / (*XLxyPVE)[iX]) >= 2. && ((*XLxyPV)[iX] / (*XLxyPVE)[iX])  <= 3.0);
    nonPromptRegion = (((*XLxyPV)[iX] / (*XLxyPVE)[iX]) > 3.0);
    promptRegion    = (((*XLxyPV)[iX] / (*XLxyPVE)[iX]) < 2.);

    //new_kkMass = (fabs(Phi.M()-Phi_mean)<3.0*Phi_sigma);
    b0_side = ((fabs((Phi.M()-Phi_mean)) <= 7.0*Phi_sigma) && ((fabs((Phi.M()-Phi_mean)) >= 5.0*Phi_sigma)));
    b0_side_l = ((fabs((Phi.M()-Phi_mean)) <= 7.0*Phi_sigma) && ((fabs((Phi.M()-Phi_mean)) >= 5.0*Phi_sigma)));
    b0_side_r = ((fabs((Phi.M()-Phi_mean)) <= 7.0*Phi_sigma) && ((fabs((Phi.M()-Phi_mean)) >= 5.0*Phi_sigma)));
    b0_signal = ((fabs((Phi.M()-Phi_mean)) <= 3.0*Phi_sigma));


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
                    XCandMassHistos[histokey]->Fill(XCand.M());
                    JPsiMassHistos[histokey]->Fill(JPsiOriginal.M());
                    PhiMassHistos[histokey]->Fill(Phi.M());
                  }
                }

              }
            }
          }
        }
      }
    }

    if(muonsCuts && kaonsCuts && jPsiCuts && extraCuts && HLT_4_vAny && nonPromptRegion && CWMass)
    {

      B0_PhiMassHisto->Fill(Phi.M());

      if(b0_signal)
        {
          B0_Cand_Mass->Fill(XCand.M());
          B0_Cand_Masses.push_back(XCand.M());
        }
      if(b0_side_l)
      {
        B0_Cand_Mass_L_Side->Fill(XCand.M());
        B0_Cand_Masses_L_Sides.push_back(XCand.M());
      }
      if(b0_side)
      {
        B0_Cand_Mass_R_Side->Fill(XCand.M());
        B0_Cand_Masses_R_Sides.push_back(XCand.M());
      }

      B0_phi_masses.push_back(Phi.M());

    }


    cutsFlags.clear();
    winsFlags.clear();
    regsFlags.clear();

  }




  if(B0_phi_masses.size()==1)
    B0_PhiMassHisto_NoM->Fill(B0_phi_masses[0]);
  if(B0_Cand_Masses.size()==1)
    B0_Cand_Mass_NoM->Fill(B0_Cand_Masses[0]);
  if(B0_Cand_Masses_R_Sides.size()==1)
    B0_Cand_Mass_L_Side_NoM->Fill(B0_Cand_Masses_R_Sides[0]);
  if(B0_Cand_Masses_L_Sides.size()==1)
    B0_Cand_Mass_R_Side_NoM->Fill(B0_Cand_Masses_L_Sides[0]);

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

          }
        }
      }
    }

    B0_PhiMassHisto->Write();
    B0_PhiMassHisto_NoM->Write();

    B0_Cand_Mass_L_Side->Write();
    B0_Cand_Mass_L_Side_NoM->Write();

    B0_Cand_Mass_R_Side->Write();
    B0_Cand_Mass_R_Side_NoM->Write();

    B0_Cand_Mass->Write();
    B0_Cand_Mass_NoM->Write();

    OutFile->Print();
    fOutput->Add(OutFile);
    gDirectory = savedir;
    fOut->Close();

  }


}


void X4140_MuMuKK::Terminate()
{


}
