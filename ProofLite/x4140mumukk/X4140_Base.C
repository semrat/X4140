#define X4140_Base_cxx
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

#include "X4140_Base.h"
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

void X4140_Base::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

}


void X4140_Base::SlaveBegin(TTree * /*tree*/)
{
  //std::cout << "SLAVE BEGIN" << std::endl;

  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

  std::string outputString = "X4140_Base_Ref_NoQual_Side4-6_NP3.0_CW5.15-5.55_BinWise_XnP.root";
  OutFile = new TProofOutputFile( outputString.data() );
  fOut = OutFile->OpenFile("RECREATE");
  if (!(fOut=OutFile->OpenFile("RECREATE")))
  {
    Warning("SlaveBegin","Problems opening file: %s%s", OutFile->GetDir(), OutFile->GetFileName() );
  }

  clabels.push_back("all"); clabels.push_back("muons"); clabels.push_back("jpsi"); clabels.push_back("kaons"); clabels.push_back("extra");
  regions.push_back("all"); regions.push_back("prompt"); regions.push_back("mixed"); regions.push_back("nonprompt");
  windows.push_back("all"); windows.push_back("cw"); windows.push_back("sw");
  hlts.push_back("all");hlts.push_back("any");hlts.push_back("hlt4");hlts.push_back("hlt8");

  ////////////////// Histograms //////////////////
  JPsi_mass = 3.096916; /// pdg mass
  Phi_mass = 1.019455; /// pdg mass

  CW_Phi_mean = 1.0195879;
  CW_Phi_sigma = 4.25298e-03;//2.35607e-03;//2.28400e-03;

  SW_Phi_mean = 1.0194159;
  SW_Phi_sigma = 5.10795e-03; //2.35607e-03;//2.28400e-03;
  SW_NP_Phi_sigma = 5.26254e-03;
  SW_NP_Phi_mean = 1.01955223318;

  std::string xcandHisto        = "Xcand_histo_";
  std::string xcandHistosignal  = "Xcand_histo_signal_";
  std::string xcandHistoLSide   = "Xcand_histo_l_side_";
  std::string xcandHistoRSide   = "Xcand_histo_r_side_";
  std::string jpsiHisto         = "JPsi_histo_";
  std::string jpsiHistoOriginal = "JPsi_histo_original_";
  std::string phiHisto          = "Phi_hist_";

  std::string bufferstring;

  for (size_t l = 0; l < hlts.size(); l++)
  {
    for (size_t j = 0; j < windows.size(); ++j)
    {
      for (size_t k = 0; k < regions.size(); ++k)
      {
        for (size_t i = 0; i < clabels.size(); ++i)
        {
          std::tuple < std::string, std::string, std::string, std::string > histokey(hlts[l],windows[j],regions[k],clabels[i]);
          std::string suffix = hlts[l] + "_" + windows[j] + "_" + regions[k] + "_" + clabels[i];

          XCandMass[histokey]  = new TH1F ((xcandHisto + suffix).data(),(xcandHisto + suffix).data(),4000, 4.0, 6.0);
          JPsiMass[histokey]   = new TH1F ((jpsiHisto + suffix).data(),(jpsiHisto + suffix).data(),2000, 2.5, 3.5);
          PhiMass[histokey]    = new TH1F ((phiHisto + suffix).data(),(phiHisto + suffix).data(),2000, 0.5, 1.5);
          XCandMass_Signal[histokey] = new TH1F ((xcandHistosignal + suffix).data(),(xcandHistosignal + suffix).data(),4000, 4.0, 6.0);
          JPsiMass_Original[histokey]   = new TH1F ((jpsiHistoOriginal + suffix).data(),(jpsiHistoOriginal + suffix).data(),2000, 2.5, 3.5);
          XCandMass_R_Side[histokey]   = new TH1F ((xcandHistoRSide + suffix).data(),(xcandHistoRSide + suffix).data(),4000, 4.0, 6.0);
          XCandMass_L_Side[histokey]   = new TH1F ((xcandHistoLSide + suffix).data(),(xcandHistoLSide + suffix).data(),4000, 4.0, 6.0);
          //XCandMassDeltaM[histokey] = new TH1F ((xcandHistoDM + suffix).data(),(xcandHistoDM + suffix).data(),2000, 4.0, 6.0);
        }
      }
    }
  }


  CW_Mass = new TH1F ("CW_Cand_Mass","CW_Cand_Mass",1000, 5.1, 5.6);
  CW_Mass_NoM= new TH1F ("CW_Cand_Mass_NoM","CW_Cand_Mass_NoM",1000, 5.1, 5.6);
  CW_Mass_Signal = new TH1F ("CW_Mass_Signal","CW_Mass_Signal",1000, 5.1, 5.6);
  CW_Mass_Signal_NoM = new TH1F ("CW_Mass_Signal_NoM","CW_Mass_Signal_NoM",1000, 5.1, 5.6);
  CW_JPsi = new TH1F ("CW_JPsi","CW_JPsi",2000, 2.5, 3.5);
  CW_JPsi_NoM = new TH1F ("CW_JPsi_NoM","CW_JPsi_NoM",2000, 2.5, 3.5);
  CW_PhiMass = new TH1F ("CW_PhiMass","CW_JPsi_NoM",2000, 0.5, 1.5);
  CW_PhiMass_NoM = new TH1F ("CW_PhiMass_NoM","CW_JPsi_NoM",2000, 0.5, 1.5);
  CW_Mass_R_Side = new TH1F ("CW_Mass_R_Side","CW_Mass_R_Side",1000, 5.1, 5.6);
  CW_Mass_R_Side_NoM = new TH1F ("CW_Mass_R_Side_NoM","CW_Mass_R_Side_NoM",1000, 5.1, 5.6);
  CW_Mass_L_Side = new TH1F ("CW_Mass_L_Side","CW_Mass_L_Side",1000, 5.1, 5.6);
  CW_Mass_L_Side_NoM = new TH1F ("CW_Mass_L_Side_NoM","CW_Mass_L_Side_NoM",1000, 5.1, 5.6);

  SW_Mass = new TH1F ("SW_Cand_Mass","SW_Cand_Mass",2000, 4.0, 5.0);
  SW_Mass_NoM = new TH1F ("SW_Cand_Mass_NoM","SW_Cand_Mass_NoM",2000, 4.0, 5.0);
  SW_Mass_Signal = new TH1F ("SW_Mass_Signal","SW_Mass_Signal",2000, 4.0, 5.0);
  SW_Mass_Signal_NoM= new TH1F ("SW_Mass_Signal_NoM","SW_Mass_Signal_NoM",2000, 4.0, 5.0);
  SW_JPsi = new TH1F ("SW_JPsi","SW_JPsi",2000, 2.5, 3.5);
  SW_JPsi_NoM = new TH1F ("SW_JPsi_NoM","SW_JPsi_NoM",2000, 2.5, 3.5);
  SW_PhiMass = new TH1F ("SW_PhiMass","SW_JPsi_NoM",2000, 0.5, 1.5);
  SW_PhiMass_NoM = new TH1F ("SW_PhiMass_NoM","SW_JPsi_NoM",2000, 0.5, 1.5);
  SW_Mass_R_Side = new TH1F ("SW_Mass_R_Side","SW_Mass_R_Side",2000, 4.0, 5.0);
  SW_Mass_R_Side_NoM = new TH1F ("SW_Mass_R_Side_NoM","SW_Mass_R_Side_NoM",2000, 4.0, 5.0);
  SW_Mass_L_Side = new TH1F ("SW_Mass_L_Side","SW_Mass_L_Side",2000, 4.0, 5.0);
  SW_Mass_L_Side_NoM = new TH1F ("SW_Mass_L_Side_NoM","SW_Mass_L_Side_NoM",2000, 4.0, 5.0);

  SW_NP_Mass = new TH1F ("SW_NP_Cand_Mass","SW_NP_Cand_Mass",2000, 4.0, 5.0);
  SW_NP_Mass_NoM = new TH1F ("SW_NP_Cand_Mass_NoM","SW_NP_Cand_Mass_NoM",2000, 4.0, 5.0);
  SW_NP_Mass_Signal = new TH1F ("SW_NP_Mass_Signal","SW_NP_Mass_Signal",2000, 4.0, 5.0);
  SW_NP_Mass_Signal_NoM= new TH1F ("SW_NP_Mass_Signal_NoM","SW_NP_Mass_Signal_NoM",2000, 4.0, 5.0);
  SW_NP_JPsi = new TH1F ("SW_NP_JPsi","SW_NP_JPsi",2000, 2.5, 3.5);
  SW_NP_JPsi_NoM = new TH1F ("SW_NP_JPsi_NoM","SW_NP_JPsi_NoM",2000, 2.5, 3.5);
  SW_NP_PhiMass = new TH1F ("SW_NP_PhiMass","SW_NP_JPsi_NoM",2000, 0.5, 1.5);
  SW_NP_PhiMass_NoM = new TH1F ("SW_NP_PhiMass_NoM","SW_NP_JPsi_NoM",2000, 0.5, 1.5);
  SW_NP_Mass_R_Side = new TH1F ("SW_NP_Mass_R_Side","SW_NP_Mass_R_Side",2000, 4.0, 5.0);
  SW_NP_Mass_R_Side_NoM = new TH1F ("SW_NP_Mass_R_Side_NoM","SW_NP_Mass_R_Side_NoM",2000, 4.0, 5.0);
  SW_NP_Mass_L_Side = new TH1F ("SW_NP_Mass_L_Side","SW_NP_Mass_L_Side",2000, 4.0, 5.0);
  SW_NP_Mass_L_Side_NoM = new TH1F ("SW_NP_Mass_L_Side_NoM","SW_NP_Mass_L_Side_NoM",2000, 4.0, 5.0);

  double cwsteps = (5.0 - 4.0)/bwwidth;
  double swsteps = (5.6 - 5.1)/bwwidth;

  for (int i = 0; i < bwwidth; i++) {
    bufferstring = "SW_PhiMass_BinWise_" + std::to_string(i);
    SW_PhiMass_BinWise.push_back(new TH1F (bufferstring.data(),bufferstring.data(),2000, 0.5, 1.5));
    bufferstring = "CW_PhiMass_BinWise_" + std::to_string(i);
    CW_PhiMass_BinWise.push_back(new TH1F (bufferstring.data(),bufferstring.data(),2000, 0.5, 1.5));

    bufferstring = "SW_PhiMass_BinWise_NoM_" + std::to_string(i);
    SW_PhiMass_BinWise_NoM.push_back(new TH1F (bufferstring.data(),bufferstring.data(),2000, 0.5, 1.5));
    bufferstring = "CW_PhiMass_BinWise_NoM_" + std::to_string(i);
    CW_PhiMass_BinWise_NoM.push_back(new TH1F (bufferstring.data(),bufferstring.data(),2000, 0.5, 1.5));

    cw_binwise.push_back(CW_Mass->GetBinLowEdge(1) + double(double(i) * cwsteps));
    sw_binwise.push_back(SW_Mass->GetBinLowEdge(1) + double(double(i) * swsteps));

  }

  cw_binwise.push_back(CW_Mass->GetBinLowEdge(1) + double(double(i) * bwwidth));
  sw_binwise.push_back(SW_Mass->GetBinLowEdge(1) + double(double(i) * bwwidth));

}

bool X4140_Base::Process(Long64_t entry)
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

  int X4140_Baseer = 0,jPsiMuMuCounter = 0, jPsis = 0;
  int muonQual[] = {1,3,4,12};

  // for (std::map < std::pair <int,int>, int>::iterator lumiMapIt=lumiMap.begin(); lumiMapIt!=lumiMap.end(); ++lumiMapIt)
  //  runlumi << lumiMapIt->first.first << " : " << lumiMapIt->first.second << std::endl;


  // runLumi->Fill(runNum,lumiNum);
  std::vector<double> xCandsB0sMasses,B0_Cand_Masses_L_Sides,B0_Cand_Masses_R_Sides,B0_Cand_Masses,B0_phi_masses;

  std::vector<double> CW_Mass_NoMVec;
  std::vector<double> CW_Mass_Signal_NoMVec;
  std::vector<double> CW_JPsi_NoMVec;
  std::vector<double> CW_PhiMass_NoMVec;
  std::vector<double> CW_Mass_R_Side_NoMVec;
  std::vector<double> CW_Mass_L_Side_NoMVec;

  std::vector<double> SW_Mass_NoMVec;
  std::vector<double> SW_Mass_Signal_NoMVec;
  std::vector<double> SW_JPsi_NoMVec;
  std::vector<double> SW_PhiMass_NoMVec;
  std::vector<double> SW_Mass_R_Side_NoMVec;
  std::vector<double> SW_Mass_L_Side_NoMVec;


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

    bool extraCuts = false;

    bool b0_side = false, b0_signal = false, b0_side_l = false, b0_side_r = false;
    bool y_side = false, y_signal = false, y_side_l = false, y_side_r = false;
    bool y_NP_side = false, y_NP_signal = false, y_NP_side_l = false, y_NP_side_r = false;

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

    TLorentzVector mu1, mu2, oMu1, oMu2;

    mu1.SetPxPyPzE((*Muon1Px_MuMuKK)[iX],(*Muon1Py_MuMuKK)[iX],(*Muon1Pz_MuMuKK)[iX],(*Muon1E_MuMuKK)[iX]);
    mu2.SetPxPyPzE((*Muon2Px_MuMuKK)[iX],(*Muon2Py_MuMuKK)[iX],(*Muon2Pz_MuMuKK)[iX],(*Muon2E_MuMuKK)[iX]);

    mu1_E = sqrt( pow((*muPx)[iMu1], 2) + pow((*muPy)[iMu1], 2) + pow((*muPz)[iMu1], 2) + pow(muon_mass, 2) ) ;
    mu2_E = sqrt( pow((*muPx)[iMu2], 2) + pow((*muPy)[iMu2], 2) + pow((*muPz)[iMu2], 2) + pow(muon_mass, 2) ) ;
    oMu1.SetPxPyPzE( (*muPx)[iMu1], (*muPy)[iMu1], (*muPz)[iMu1], mu1_E) ;
    oMu2.SetPxPyPzE( (*muPx)[iMu2], (*muPy)[iMu2], (*muPz)[iMu2], mu2_E) ;

    TLorentzVector JPsi;
    JPsi = mu1 + mu2;

    TLorentzVector JPsiOriginal;
    JPsiOriginal = oMu1 + oMu2;

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

    SWMass = (((XCand.M() > 4.0) && (XCand.M() < 5.0)));
    CWMass = ((XCand.M() > 5.15) && (XCand.M() < 5.55));

    mixedRegion     = (((*XLxyPV)[iX] / (*XLxyPVE)[iX]) >= 2. && ((*XLxyPV)[iX] / (*XLxyPVE)[iX])  <= 3.0);
    nonPromptRegion = (((*XLxyPV)[iX] / (*XLxyPVE)[iX]) > 3.0);
    promptRegion    = (((*XLxyPV)[iX] / (*XLxyPVE)[iX]) < 2.);

    // muonQualityCut = ( ((*muQual)[iMu1]) & (1 << muonQual[3]) ) && ( ((*muQual)[iMu2]) & (1 << muonQual[3]) );
    // muonChiCut     = (( ( (*muChi2)[iMu1] / (*muNDF)[iMu1] ) < 3 ) && ( ( (*muChi2)[iMu2] / (*muNDF)[iMu2] ) < 3 ));
    // muonPhitsCut   = ((*muPhits)[iMu1] > 0 && (*muPhits)[iMu2] > 0);
    // muonShitsCut   = ((*muShits)[iMu1] > 5 && (*muShits)[iMu2] > 5);
    muonDZPVCut    = (fabs((*muDzVtx)[iMu1]) < 20.0 && fabs((*muDzVtx)[iMu2]) < 20.0);
    muonDXYPVCut   = (fabs((*muDxyVtx)[iMu1]) < 0.3 && fabs((*muDxyVtx)[iMu2]) < 0.3);

    muonsCuts = muonQualityCut && muonChiCut && muonShitsCut && muonPhitsCut && muonDZPVCut && muonDXYPVCut;
    muonsCuts = muonDZPVCut && muonDXYPVCut;

    jPsiPtCut      = (JPsi.Pt() > 7.0);
    jPsiMassCut    =  (JPsiOriginal.M()<3.4 && JPsiOriginal.M()>2.8);
    jPsiVtxCut     = ((*MuMuVtx_CL)[iJPsi]) > 0.1;
    jPsiMuEtaPtCut = (fabs(mu1.Eta()) < 2.2) && (fabs(mu2.Eta()) < 2.2);
    jPsiMusPtCut   = ((mu1.Pt() > 4.0) && (mu2.Pt() > 4.0));

    jPsiCuts = jPsiPtCut && jPsiMassCut && jPsiVtxCut  && jPsiMuEtaPtCut && jPsiMuEtaPtCut && jPsiMusPtCut;

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

    cosAlphaCut = (fabs((*XCosAlphaPV)[iX]) > 0.99);
    vtxCLCut =  (((*XVtx_CL)[iX]) > 0.01);

    extraCuts = vtxCLCut && cosAlphaCut;

    mixedRegion     = (((*XLxyPV)[iX] / (*XLxyPVE)[iX]) >= 2. && ((*XLxyPV)[iX] / (*XLxyPVE)[iX])  <= 3.0);
    nonPromptRegion = (((*XLxyPV)[iX] / (*XLxyPVE)[iX]) > 3.0);
    promptRegion    = (((*XLxyPV)[iX] / (*XLxyPVE)[iX]) < 2.);

    //new_kkMass = (fabs(Phi.M()-Phi_mean)<3.0*Phi_sigma);
    b0_side = ((fabs((Phi.M()-CW_Phi_mean)) <= 6.0*CW_Phi_sigma) && (fabs((Phi.M()-CW_Phi_mean)) >= 4.0*CW_Phi_sigma));
    b0_side_r = (((Phi.M()-CW_Phi_mean) <= 6.0*CW_Phi_sigma) && ((Phi.M()-CW_Phi_mean) >= 4.0*CW_Phi_sigma));
    b0_side_l = (((Phi.M()-CW_Phi_mean) >= -6.0*CW_Phi_sigma) && ((Phi.M()-CW_Phi_mean) <= -4.0*CW_Phi_sigma));
    b0_signal = ((fabs((Phi.M()-CW_Phi_mean)) <= 3.0*CW_Phi_sigma));

    y_side = ((fabs((Phi.M()-SW_Phi_mean)) <= 6.0*SW_Phi_sigma) && (fabs((Phi.M()-SW_Phi_mean)) >= 4.0*SW_Phi_sigma));
    y_side_r = (((Phi.M()-SW_Phi_mean) <= 6.0*SW_Phi_sigma) && ((Phi.M()-SW_Phi_mean) >= 4.0*SW_Phi_sigma));
    y_side_l = (((Phi.M()-SW_Phi_mean) >= -6.0*SW_Phi_sigma) && ((Phi.M()-SW_Phi_mean) <= -4.0*SW_Phi_sigma));
    y_signal = ((fabs((Phi.M()-SW_Phi_mean)) <= 3.0*SW_Phi_sigma));

    y_NP_side_r = (((Phi.M()-SW_NP_Phi_mean) <= 6.0*SW_NP_Phi_sigma) && ((Phi.M()-SW_NP_Phi_mean) >= 4.0*SW_NP_Phi_sigma));
    y_NP_side_l = (((Phi.M()-SW_NP_Phi_mean) >= -6.0*SW_NP_Phi_sigma) && ((Phi.M()-SW_NP_Phi_mean) <= -4.0*SW_NP_Phi_sigma));
    y_NP_signal = ((fabs((Phi.M()-SW_NP_Phi_mean)) <= 3.0*SW_NP_Phi_sigma));

    winsFlags.push_back(true);
    winsFlags.push_back(CWMass);
    winsFlags.push_back(SWMass);

    regsFlags.push_back(true);
    regsFlags.push_back(promptRegion);
    regsFlags.push_back(mixedRegion);
    regsFlags.push_back(nonPromptRegion);

    cutsFlags.push_back(true);
    cutsFlags.push_back((muonsCuts));
    cutsFlags.push_back((jPsiCuts));
    cutsFlags.push_back((kaonsCuts));
    cutsFlags.push_back((extraCuts));

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
                    XCandMass[histokey]->Fill(XCand.M());
                    JPsiMass[histokey]->Fill(JPsi.M());
                    JPsiMass_Original[histokey]->Fill(JPsiOriginal.M());
                    PhiMass[histokey]->Fill(Phi.M());
                    if(b0_signal)
                      XCandMass_Signal[histokey]->Fill(XCand.M());
                    if(b0_side_r)
                      XCandMass_R_Side[histokey]->Fill(XCand.M());
                    if(b0_side_l)
                      XCandMass_L_Side[histokey]->Fill(XCand.M());

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

      CW_PhiMass->Fill(Phi.M());
      CW_PhiMass_NoMVec.push_back(Phi.M());

      CW_Mass->Fill(XCand.M());
      CW_Mass_NoMVec.push_back(XCand.M());

      CW_JPsi->Fill(JPsiOriginal.M());
      CW_JPsi_NoMVec.push_back(JPsiOriginal.M());

      if(b0_signal)
        {
          CW_Mass_Signal->Fill(XCand.M());
          CW_Mass_Signal_NoMVec.push_back(XCand.M());
        }
      if(b0_side_l)
      {
        CW_Mass_L_Side->Fill(XCand.M());
        CW_Mass_L_Side_NoMVec.push_back(XCand.M());
      }
      if(b0_side_r)
      {
        CW_Mass_R_Side->Fill(XCand.M());
        CW_Mass_R_Side_NoMVec.push_back(XCand.M());
      }

      for (int i = 0; i < bwwidth; i++)
      {
        if(XCand.M() >= cw_binwise[i] && XCand.M() < cw_binwise[i+1])
          CW_PhiMass_BinWise[i]->Fill(Phi.M());
      }

    }

    if(muonsCuts && kaonsCuts && jPsiCuts && extraCuts && HLT_8_vAny && promptRegion && SWMass)
    {

      SW_PhiMass->Fill(Phi.M());
      SW_PhiMass_NoMVec.push_back(Phi.M());

      SW_Mass->Fill(XCand.M());
      SW_Mass_NoMVec.push_back(XCand.M());

      SW_JPsi->Fill(JPsiOriginal.M());
      SW_JPsi_NoMVec.push_back(JPsiOriginal.M());

      if(y_signal)
        {
          SW_Mass_Signal->Fill(XCand.M());
          SW_Mass_Signal_NoMVec.push_back(XCand.M());
        }
      if(y_side_l)
      {
        SW_Mass_L_Side->Fill(XCand.M());
        SW_Mass_L_Side_NoMVec.push_back(XCand.M());
      }
      if(y_side_r)
      {
        SW_Mass_R_Side->Fill(XCand.M());
        SW_Mass_R_Side_NoMVec.push_back(XCand.M());
      }

      for (int i = 0; i < bwwidth; i++)
      {
        if(XCand.M() >= sw_binwise[i] && XCand.M() < sw_binwise[i+1])
          SW_PhiMass_BinWise[i]->Fill(Phi.M());
      }

    }

    if(muonsCuts && kaonsCuts && jPsiCuts && extraCuts && HLT_4_vAny && nonPromptRegion && SWMass)
    {

      SW_NP_PhiMass->Fill(Phi.M());
      // SW_NP_PhiMass_NoMVec.push_back(Phi.M());

      SW_NP_Mass->Fill(XCand.M());
      // SW_NP_Mass_NoMVec.push_back(XCand.M());

      SW_NP_JPsi->Fill(JPsiOriginal.M());
      // SW_NP_JPsi_NoMVec.push_back(JPsiOriginal.M());

      if(y_NP_signal)
        {
          SW_NP_Mass_Signal->Fill(XCand.M());
          //SW_NP_Mass_Signal_NoMVec.push_back(XCand.M());
        }
      if(y_NP_side_l)
      {
        SW_NP_Mass_L_Side->Fill(XCand.M());
        //SW_NP_Mass_L_Side_NoMVec.push_back(XCand.M());
      }
      if(y_NP_side_r)
      {
        SW_NP_Mass_R_Side->Fill(XCand.M());
        //SW_NP_Mass_R_Side_NoMVec.push_back(XCand.M());
      }

    }

    cutsFlags.clear();
    winsFlags.clear();
    regsFlags.clear();

  }

  if(SW_Mass_NoMVec.size()==1)
  {

    SW_PhiMass_NoM->Fill(SW_PhiMass_NoMVec[0]);

    SW_Mass_NoM->Fill(SW_Mass_NoMVec[0]);

    SW_JPsi_NoM->Fill(SW_JPsi_NoMVec[0]);

    if(SW_Mass_Signal_NoMVec.size()==1)
      SW_Mass_Signal_NoM->Fill(SW_Mass_Signal_NoMVec[0]);
    if(SW_Mass_L_Side_NoMVec.size()==1)
      SW_Mass_L_Side_NoM->Fill(SW_Mass_L_Side_NoMVec[0]);
    if(SW_Mass_R_Side_NoMVec.size()==1)
      SW_Mass_R_Side_NoM->Fill(SW_Mass_R_Side_NoMVec[0]);

    for (int i = 0; i < bwwidth; i++)
    {
      if(SW_Mass_NoMVec[0] >= sw_binwise[i] && SW_Mass_NoMVec[0] < sw_binwise[i+1])
        SW_PhiMass_BinWise_NoM[i]->Fill(SW_PhiMass_NoM[i]);
    }
  }

  if(CW_Mass_NoMVec.size()==1)
  {

    CW_PhiMass_NoM->Fill(CW_PhiMass_NoMVec[0]);

    CW_Mass_NoM->Fill(CW_Mass_NoMVec[0]);

    CW_JPsi_NoM->Fill(CW_JPsi_NoMVec[0]);

    if(CW_Mass_Signal_NoMVec.size()==1)
      CW_Mass_Signal_NoM->Fill(CW_Mass_Signal_NoMVec[0]);
    if(CW_Mass_L_Side_NoMVec.size()==1)
      CW_Mass_L_Side_NoM->Fill(CW_Mass_L_Side_NoMVec[0]);
    if(CW_Mass_R_Side_NoMVec.size()==1)
      CW_Mass_R_Side_NoM->Fill(CW_Mass_R_Side_NoMVec[0]);

    for (int i = 0; i < bwwidth; i++)
    {
      if(CW_Mass_NoMVec[0] >= cw_binwise[i] && CW_Mass_NoMVec[0] < cw_binwise[i+1])
          CW_PhiMass_BinWise_NoM[i]->Fill(CW_PhiMass_NoMVec[0]);
    }

  }


  // if(B0_phi_masses.size()==1)
  //   B0_PhiMassHisto_NoM->Fill(B0_phi_masses[0]);
  // if(B0_Cand_Masses.size()==1)
  //   B0_Cand_Mass_NoM->Fill(B0_Cand_Masses[0]);
  // if(B0_Cand_Masses_R_Sides.size()==1)
  //   B0_Cand_Mass_R_Side_NoM->Fill(B0_Cand_Masses_R_Sides[0]);
  // if(B0_Cand_Masses_L_Sides.size()==1)
  //   B0_Cand_Mass_L_Side_NoM->Fill(B0_Cand_Masses_L_Sides[0]);

  return kTRUE;
}


void X4140_Base::SlaveTerminate()
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

    for (size_t l = 0; l < hlts.size(); l++)
    {
      for (size_t j = 0; j < windows.size(); ++j)
      {
        for (size_t k = 0; k < regions.size(); ++k)
        {
          for (size_t i = 0; i < clabels.size(); ++i)
          {
            std::tuple < std::string, std::string, std::string, std::string > histokey(hlts[l],windows[j],regions[k],clabels[i]);
            std::string suffix = hlts[l] + "_" + windows[j] + "_" + regions[k] + "_" + clabels[i];

            XCandMass[histokey]->Write();
            JPsiMass[histokey]->Write();
            PhiMass[histokey]->Write();
            JPsiMass_Original[histokey]->Write();
            XCandMass_Signal[histokey]->Write();
            XCandMass_R_Side[histokey]->Write();
            XCandMass_L_Side[histokey]->Write();

          }
        }
      }
    }

    CW_Mass->Write(); CW_Mass_NoM->Write();
    CW_Mass_Signal->Write(); CW_Mass_Signal_NoM->Write();
    CW_JPsi->Write(); CW_JPsi_NoM->Write();
    CW_PhiMass->Write(); CW_PhiMass_NoM->Write();
    CW_Mass_R_Side->Write(); CW_Mass_R_Side_NoM->Write();
    CW_Mass_L_Side->Write(); CW_Mass_L_Side_NoM->Write();

    SW_Mass->Write(); SW_Mass_NoM->Write();
    SW_JPsi->Write(); SW_JPsi_NoM->Write();
    SW_Mass_Signal->Write(); SW_Mass_Signal_NoM->Write();
    SW_PhiMass->Write(); SW_PhiMass_NoM->Write();
    SW_Mass_R_Side->Write(); SW_Mass_R_Side_NoM->Write();
    SW_Mass_L_Side->Write(); SW_Mass_L_Side_NoM->Write();

    SW_NP_Mass->Write(); SW_NP_Mass_NoM->Write();
    SW_NP_JPsi->Write(); SW_NP_JPsi_NoM->Write();
    SW_NP_Mass_Signal->Write(); SW_NP_Mass_Signal_NoM->Write();
    SW_NP_PhiMass->Write(); SW_NP_PhiMass_NoM->Write();
    SW_NP_Mass_R_Side->Write(); SW_NP_Mass_R_Side_NoM->Write();
    SW_NP_Mass_L_Side->Write(); SW_NP_Mass_L_Side_NoM->Write();

    for (size_t i = 0; i < bwwidth; i++) {
      SW_PhiMass_BinWise[i]->Write();
      CW_PhiMass_BinWise[i]->Write();

      SW_PhiMass_BinWise_NoM[i]->Write();
      CW_PhiMass_BinWise_NoM[i]->Write();
    }
    OutFile->Print();
    fOutput->Add(OutFile);
    gDirectory = savedir;
    fOut->Close();

  }


}


void X4140_Base::Terminate()
{


}
