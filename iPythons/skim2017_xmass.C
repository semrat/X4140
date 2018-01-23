#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TDirectory.h>
#include <iostream>
#include <bitset>
#include <TCanvas.h>
#include <algorithm>
#include <TLegend.h>
#include <TStyle.h>
#include <string>

int noHlts = 13;

int colors[17] = {kRed,kBlue,kGreen,kMagenta,kBlack,kRed + 1,kBlue + 1,kGreen + 1,kMagenta + 1,kBlack + 1,
              kRed - 1,kBlue - 1,kGreen - 1,kMagenta - 1,kBlack - 1, kOrange, kOrange + 1};

std::string hltsName[13] = {"HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi","HLT_Mu20_TkMu0_Phi",
                            "HLT_Dimuon14_Phi_Barrel_Seagulls","HLT_Mu25_TkMu0_Phi","HLT_Dimuon24_Phi_noCorrL1",
                            "HLT_DoubleMu4_JpsiTrkTrk_Displaced",
                            "HLT_DoubleMu4_JpsiTrk_Displaced",
                            "HLT_DoubleMu4_Jpsi_Displaced",
                            "HLT_DoubleMu4_3_Jpsi_Displaced",
                            "HLT_Dimuon20_Jpsi_Barrel_Seagulls",
                            "HLT_Dimuon25_Jpsi",
                            "HLT_Dimuon0_Jpsi"};

int skimXTree(std::string path, std::string filename, std::string treename = "xTree", std::string dirname = "rootuple")
{

   TFile *oldfile = TFile::Open((path+filename).data());
   TDirectory *directory = (TDirectory*)oldfile->Get(dirname.data());
   TTree *oldtree = (TTree*)directory->Get(treename.data());
   Long64_t nentries = oldtree->GetEntries();
   ULong64_t event   = 0;
   oldtree->SetBranchAddress("event",&event);
   //Create a new file + a clone of old tree in new file
   TFile *newfile = new TFile((treename + "_skim_" + filename).data(),"RECREATE");
   TTree *newtree = oldtree->CloneTree();

   newtree->Print();
   newtree->Write();

   return 0;

}



int selectXTree()
{

   TFile *oldfile = TFile::Open("/Users/adrianodiflorio/Documents/Git/X4140/iPythons/skimmedNP.root");

   TTree *oldtree = (TTree*)oldfile->Get("xTree");
   Long64_t nentries = oldtree->GetEntries();
   Double_t xyl   = 0.0;
   Double_t xylErr   = 0.0;
   Double_t cosA  = 0.0;

   Double_t phiM  = 0.0;
   Double_t jPsiM  = 0.0;

   Double_t ctau  = 0.0;
   Double_t ctauErr  = 0.0;

   Double_t vProb  = 0.0;

   Int_t phiMType = 0, phiPType = 0;
   UInt_t phi_trigger = 0, jpsi_trigger = 0;

   oldtree->SetBranchAddress("vProb",&vProb);
   oldtree->SetBranchAddress("l_xy",&xyl);
   oldtree->SetBranchAddress("lErr_xy",&xylErr);
   oldtree->SetBranchAddress("cosAlpha",&cosA);
   oldtree->SetBranchAddress("phi_M",&phiM);
   oldtree->SetBranchAddress("jpsi_M",&jPsiM);

   oldtree->SetBranchAddress("phi_muonM_type",&phiMType);
   oldtree->SetBranchAddress("phi_muonP_type",&phiPType);

   oldtree->SetBranchAddress("ctauPV",&ctau);
   oldtree->SetBranchAddress("ctauErrPV",&ctauErr);

   oldtree->SetBranchAddress("phi_trigger",&phi_trigger);
   oldtree->SetBranchAddress("jpsi_trigger",&jpsi_trigger);
   //Create a new file + a clone of old tree in new file
   TFile *newfile = new TFile("skimmedNPCos.root","RECREATE");
   TTree *newtree = oldtree->CloneTree(0);

   for (Long64_t i=0;i<nentries; i++) {
      oldtree->GetEntry(i);
      std::bitset<16> pM(phiMType);
      std::bitset<16> pP(phiPType);
      // std::cout << phiMType << "-" << binary << std::endl;
      // std::cout<<vProb <<std::endl;
      //if (vProb > 0.0) newtree->Fill();
      if (jPsiM > 2.8 && phiM < 1.1 && cosA > 0.995 && phiM > 0.95 && pP.test(1) && pM.test(1) && vProb > 0.1  ) newtree->Fill();
   }

   newtree->Print();
   newtree->Write();

   return 0;

}

int drawXTree(std::string path = "/Users/adrianodiflorio/Documents/Git/X4140/iPythons/xTree.root")
{

   TFile *oldfile = TFile::Open(path.data());
   TTree *oldtree = (TTree*)oldfile->Get("xTree");

   Long64_t nentries = oldtree->GetEntries();
   Double_t xyl   = 0.0;
   Double_t xylErr   = 0.0;
   Double_t cosA  = 0.0;

   Double_t phiM  = 0.0;
   Double_t jPsiM  = 0.0;

   Double_t ctau  = 0.0;
   Double_t ctauErr  = 0.0;

   Double_t vProb  = 0.0;

   Int_t phiMType = 0, phiPType = 0;
   UInt_t phi_trigger = 0, jpsi_trigger = 0, trigger = 0;

   oldtree->SetBranchAddress("vProb",&vProb);
   oldtree->SetBranchAddress("trigger",&trigger);
   oldtree->SetBranchAddress("l_xy",&xyl);
   oldtree->SetBranchAddress("lErr_xy",&xylErr);
   oldtree->SetBranchAddress("cosAlpha",&cosA);
   oldtree->SetBranchAddress("phi_M",&phiM);
   oldtree->SetBranchAddress("jpsi_M",&jPsiM);

   oldtree->SetBranchAddress("phi_muonM_type",&phiMType);
   oldtree->SetBranchAddress("phi_muonP_type",&phiPType);

   oldtree->SetBranchAddress("ctauPV",&ctau);
   oldtree->SetBranchAddress("ctauErrPV",&ctauErr);

   oldtree->SetBranchAddress("phi_trigger",&phi_trigger);
   oldtree->SetBranchAddress("jpsi_trigger",&jpsi_trigger);
   //Create a new file + a clone of old tree in new file
   TCanvas c("c","c",1200,1200);

   TFile *newfile = new TFile("drawSkim.root","RECREATE");


   // for(int j = 0; j < 1; j++)
   // {

     TTree *newtree = oldtree->CloneTree(0);
     // TH1F* phi_triggrHist = new TH1F("phi_triggrHist","phi_triggrHist",600,0.6,1.2);
     TH1F* phi_allphiHist = new TH1F("phi_allphiHist","phi_allphiHist",1000,0.5,1.5);
     std::vector<TH1F*> phiHists;
     for (size_t i = 0; i < noHlts; i++)
     {
        //std::string name = "phi_allphiHist_" + std::to_string(i);

        phiHists.push_back(new TH1F(hltsName[i].data(),hltsName[i].data(),1000,0.5,1.5));

     }
     float max = -1.0;
     for (Long64_t i=0;i<nentries; i++) {
        oldtree->GetEntry(i);
        std::bitset<16> tB(trigger);
        std::bitset<16> pM(phiMType);
        std::bitset<16> pP(phiPType);
        max = std::max(max,float(trigger));
        for (int j = 0; j < noHlts; j++)
          if (tB.test(j) && pM.test(1) && pP.test(1)) phiHists[j]->Fill(phiM);

        phi_allphiHist->Fill(phiM);
     }

     std::cout << max << std::endl;
     //newtree->Draw("phi_M","","same");
   // }
   phi_allphiHist->SetMinimum(0.01);
   //oldtree->Draw("phi_M");
   // TH1F* phi_triggrHist = (TH1F*)gDirectory->Get("phi_triggrHist");

   phi_allphiHist->SetLineColor(kBlue);
   phi_allphiHist->Draw();

   TLegend leg(0.1,0.7,0.48,0.9);

   for (size_t i = 0; i < noHlts; i++)
    {
      phiHists[i]->SetLineColor(colors[i]);
      phiHists[i]->Draw("same");
      leg.AddEntry(phiHists[i],(phiHists[i]->GetName()),"l");
    }


   // phi_triggrHist->Draw("same");
   leg.Draw();
   c.SetLogy(1);
   c.SaveAs("triggerCheckTrakerMu.eps");

   return 0;

}

int selectXTreeHLT()
{

   TFile *oldfile = TFile::Open("/Users/adrianodiflorio/Documents/Git/X4140/iPythons/skimmed_vProb.root");

   TTree *oldtree = (TTree*)oldfile->Get("xTree");
   Long64_t nentries = oldtree->GetEntries();
   Double_t vProb   = 0.0;
   UInt_t trigger = 0;
   oldtree->SetBranchAddress("vProb",&vProb);
   oldtree->SetBranchAddress("trigger",&trigger);
   //Create a new file + a clone of old tree in new file
   TFile *newfile = new TFile("skimmed.root","RECREATE");
   TTree *newtree = oldtree->CloneTree(0);

   for (Long64_t i=0;i<nentries; i++) {
      oldtree->GetEntry(i);
      // std::bitset<16> binary(trigger);
      // std::cout << trigger << "-" << binary << std::endl;
      if(trigger>0)
      newtree->Fill();
   }

   newtree->Print();
   newtree->Write();

   return 0;

}

/*
# coding: utf-8

# In[1]:


import ROOT
from ROOT import TFile,TH1,TH1F,TCanvas,TNtuple,TTreeReader,TTreeReaderValue
from ROOT import RooFit
from ROOT.RooFit import Layout
from ROOT import RooStats
from ROOT import RooAbsData
RooAbsData.setDefaultStorageType ( RooAbsData.Tree )
from array import array
import sys


# In[3]:


from ROOT import RooRealVar,RooAbsPdf,RooChebychev,RooExponential,RooGaussian,RooAbsPdf,RooPlot,RooAddPdf,RooDataHist,RooArgSet,RooArgList
from ROOT import kGreen,kRed,kBlack,kBlue,kDashed,kDotted,kMagenta,RooVoigtian
from ROOT.RooFit import Components,LineColor,LineStyle,Name,Normalization,Range,SelectVars
from ROOT import RooDataSet,RooFormulaVar,RooLinkedList


# In[7]:


no_hlts = 13


# In[8]:


#rootfile = "/Users/adrianodiflorio/Desktop/mmkk2017/09Jan2017.root"
#rootfile = "/Users/adrianodiflorio/Desktop/mmkk2017/allphi_DataF.root"
rootfile = "/Users/adrianodiflorio/Desktop/mmkk2017/phiJpsiTriggersBCDEF.root"
inputfile = TFile(rootfile,"READ")
inputfile.ls()
xTupleDir = (inputfile.Get("rootuple"))
xTupleDir.ls()
#pTuple = (xTupleDir.Get("pTree"))
xTuple = (xTupleDir.Get("xTree"))
#jTuple = (xTupleDir.Get("jTree"))

event = 2

xTuple.SetBranchAddress("event",event);
newfile = TFile("small.root","recreate");
newtree = xTuple.CloneTree();
newtree.CopyEntries(xTuple);

newtree.Write()

sys.exit()

# In[16]:

file = TFile("newFile.root","RECREATE")
canvas = TCanvas("canvas","canvas",1200,1000)
mass = RooRealVar("xM","M(#mu#mu#mu#mu)[GeV]",5.15,5.55)
trigger = RooRealVar("trigger","trigger",0.0,10000)
vProb = RooRealVar("vProb","vProb",-1.0,1.0)
alldata = RooDataSet("alldata","alldata",xTuple,RooArgSet(mass), RooFormulaVar("vProb","vProb","vProb>0.01",RooArgList(vProb)))#,cutFormula)
frame = mass.frame(Range(5.15,5.55))
alldata.plotOn(frame,RooLinkedList())
alldata.Write()
frame.Draw()


# In[ ]:


canvas.SaveAs("testCanvas.eps")
*/
