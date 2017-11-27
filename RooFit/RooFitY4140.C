#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TROOT.h"
#include "TMath.h"

#include "RooAddPdf.h"
#include "RooAddition.h"
#include "RooExponential.h"
#include "RooChebychev.h"
#include "RooMinuit.h"
#include "RooHistPdf.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooArgList.h"
#include "TH1F.h"

#include <sys/time.h> // for timeval
#include <sys/times.h> // for tms
#include "TSystem.h" // to get number of CPUs
#include "TStyle.h" // to use gStyle
#include <TFile.h>
#include <TNtupleD.h>
#include "TPaveText.h"
#include "TPaletteAxis.h"

#include <string>

using namespace RooFit ;

timeval start, stop;
clock_t startCPU, stopCPU;
tms startProc, stopProc;

Float_t TH2_offset = 1.6;


void Y4140RooFit(TFile *inputfile,std::string histoname = "X5568_Cand_Mass_Ref")
{

  RooRealVar mass("mass","M(#mu#muKK)[GeV]",5.15,5.55);
  mass.setBins(80);
  RooRealVar mean("mean","mean of gaussian",5.38,5.31,5.41);
  RooRealVar sigma1("sigma1","width of gaussian1",0.01,0.001,0.05);
  RooRealVar sigma2("sigma2","width of gaussian2",0.005,0.001,0.01);

  // Build gaussian p.d.f in terms of x,mean and sigma
  RooGaussian gauss1("gauss1","gaussian PDF 1",mass,mean,sigma1);
  RooGaussian gauss2("gauss2","gaussian PDF 2",mass,mean,sigma2);

  RooRealVar a0("a0","a0",0.001,-1.,1.);
  RooRealVar a1("a1","a1",0.001,-0.5,0.5);
  RooRealVar a2("a2","a2",0.0001,-2.,2.);
  RooRealVar a3("a3","a3",0.00001,-2.,2.);
  RooRealVar a4("a4","a4",0.00001,-2.,2.);
  RooArgSet set(a0,a1,a2);
  RooChebychev cheb("cheb","Background",mass,set);//,a3,a4)) ;
  RooRealVar alpha("alpha","alpha",-0.5,-2.0,0.1);
  RooExponential exp("exp","exp",mass,alpha);

  RooRealVar gaussFrac("sig1frac","fraction of component 1 in signal",0.8,0.,1.);
  RooAddPdf  model("model","g1+g2",RooArgList(gauss1,gauss2),gaussFrac) ;

  RooRealVar nSig("nSig","nSig",5000,0.,10E6);
  RooRealVar nBkg("nBkg","nBkg",5000,0.,10E6);
  RooAddPdf  tot("tot","g1+g2+cheb",RooArgList(model,cheb),RooArgList(nSig,nBkg)) ;

  //TFile *inputFile = TFile::Open("../X4140_MuMuKK_KRe_MuMixed_NP3.5_Alpha99_CW5.2-5.55.root");
  //TH1F* hist = (TH1F*)inputFile->Get("Xcand_histo_hlt8_cw_nonprompt_cosalpha");
  TH1F* hist = (TH1F*)inputfile->Get(histoname.data());
  hist->Rebin(5);
  RooDataHist dh("dh","dh",RooArgList(mass),hist);
  // Construct plot frame in 'x'

  tot.fitTo(dh) ;

  RooPlot* massFrame = mass.frame(Title("B_{0}^{s} signal")) ;
  dh.plotOn(massFrame);
  tot.plotOn(massFrame);
  tot.plotOn(massFrame,Components(gauss1),LineColor(kGreen),LineStyle(kDashed),Name("gauss1"));
  tot.plotOn(massFrame,Components(gauss2),LineColor(kMagenta),LineStyle(kDashed),Name("gauss2"));
  //tot.plotOn(massFrame,Components(cheb),LineColor(kRed),LineStyle(kDotted),Name("exp"));
  tot.plotOn(massFrame,Components(cheb),LineColor(kRed),LineStyle(kDotted),Name("exp"));

  // gauss1.plotOn(massFrame);
  // gauss2.plotOn(massFrame);
  tot.paramOn(massFrame);

  TCanvas* c = new TCanvas("rf101_basics","rf101_basics",800,400) ;
  massFrame->Draw();
  //gauss1.plotOn(massFrame,LineColor(kRed)) ;

}
