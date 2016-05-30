using namespace RooFit;

void b_peak() {
  gStyle->SetOptStat(0);
  //gStyle->SetTitleOffset(1.3, "XY");
  gStyle->SetTitleOffset(1.2, "Y"); gStyle->SetTitleOffset(1.4, "X");
  gStyle->SetLabelSize(0.04) ;
  TGaxis::SetMaxDigits(3) ;

  TFile *f = TFile::Open("../26Jan_histos.root");
  TH1F* hmyPsiPKPiMassSelAltZoom = (TH1F*) f->Get("hmyPsiPKPiMassSelAltZoom") ;

  TCanvas *Bpeak_C = new TCanvas("Bpeak_C","Bpeak_C", 800, 600) ;  

  Double_t xMin = hmyPsiPKPiMassSelAltZoom->GetXaxis()->GetXmin();
  Double_t xMax = hmyPsiPKPiMassSelAltZoom->GetXaxis()->GetXmax();
  Int_t nBins = hmyPsiPKPiMassSelAltZoom->GetNbinsX();
  Double_t bin_width = (xMax-xMin)/nBins ;
  
  RooFit::RooRealVar xVar("xVar", hmyPsiPKPiMassSelAltZoom->GetXaxis()->GetTitle(), xMin, xMax) ;  xVar.setBins(nBins) ;  
  RooDataHist *MuMuPiKHist = new RooDataHist("MuMuPi_hist", hmyPsiPKPiMassSelAltZoom->GetTitle(), RooArgSet(xVar), Import(*hmyPsiPKPiMassSelAltZoom, kFALSE)) ;
  
  Double_t fitMax = xMax, fitMin = xMin ;

  Double_t meanval[3], sigma1val[3] ;
  meanval[0] = 5.28 ;  meanval[1] = 5.25 ;   meanval[2] = 5.31 ;
  sigma1val[0] = 0.01 ; sigma1val[1] = 0.0001 ; sigma1val[2] = 0.05 ;
  
  RooRealVar mean("m", "mean", meanval[0], meanval[1], meanval[2]) ;
  RooRealVar sigma1("#sigma_{1}", "sigma1", sigma1val[0], sigma1val[1], sigma1val[2]);
  //RooRealVar sigma2("#sigma_{2}", "sigma2", 0.02, 0.01, 0.2);
  
  //RooGaussian sig1PDF("sig1PDF", "Signal component 1", xVar, mean1, sigma1);
  //RooGaussian sig2PDF("sig2PDF", "Signal component 2", xVar, mean1, sigma2);
  //RooRealVar ratio("sig1frac", "fraction of component 1 in signal", 0.5, 0, 1);
  //RooAddPdf sigPDF("sigPDF", "Signal component", RooArgList(sig1PDF, sig2PDF), ratio);
  RooRealVar nSig1("N", "Number of 1^{st} signal", 1e+5, 1e+2, 1e+7);
  //RooRealVar nSig2("N_{2}", "Number of 2^{nd} signal", 1e+5, 1e+3, 1e+7);

  //RooAddPdf sigPDF("sigPDF", "Signal component", RooArgList(sig1PDF, sig2PDF), RooArgList(nSig1, nSig2)) ;
  RooGaussian sigPDF("sigPDF", "Signal component", xVar, mean, sigma1); mean.SetNameTitle("m","mean"); sigma1.SetNameTitle("#sigma","sigma");
  
  RooRealVar nSig("N_{Sig}", "Number of signal candidates", 1e+6, 1e+2, 1e+7);
  /*// not with Kai BF
  RooGaussian bumpPDF("bumpPDF", "Bump", xVar, mean3, sigma3);
  RooRealVar nBump("N_{3}", "Number of bump candidates", 1e+5, 1e+3, 1e+7);
  */
  RooRealVar c1("c_{1}", "c1", +0.02, -10, 10);
  //RooRealVar c1("c_{1}", "c1", -0.02, -10, 10);
  //RooChebychev bkgPDF("bkgPDF", "bkgPDF", xVar, RooArgSet(c1));  // with 5 sigma
  RooRealVar c2("c_{2}", "c2", -0.001, -10, 10);
  RooChebychev bkgPDF("bkgPDF", "bkgPDF", xVar, RooArgSet(c1,c2));
  //RooRealVar c3("c_{3}", "c3", -0.001, -10, 10);
  //RooChebychev bkgPDF("bkgPDF", "bkgPDF", xVar, RooArgSet(c1,c2,c3));
  //RooRealVar c4("c_{4}", "c4", -0.001, -10, 10);
  //RooChebychev bkgPDF("bkgPDF", "bkgPDF", xVar, RooArgSet(c1,c2,c3,c4));
  /* // not with Kai BF
  RooRealVar gamma("#Gamma", "gamma", -0.01, -5, -0.0005);
  RooExponential bkgPDF("bkgPDF", "Background", xVar, gamma) ;
  */
  RooRealVar nBkg("nBkg", "Bkg component", 1e+6, 1e+3, 1e+8);     

  RooAddPdf *totalPDF = new RooAddPdf("totalPDF", "totalPDF", RooArgList(sigPDF, bkgPDF), RooArgList(nSig, nBkg));  // with 5 sigma
  //RooAddPdf *totalPDF = new RooAddPdf("totalPDF", "totalPDF", RooArgList(sig1PDF, sig2PDF, bkgPDF), RooArgList(nSig1, nSig2, nBkg));     
  //RooAddPdf *totalPDF = new RooAddPdf("totalPDF", "totalPDF", RooArgList(sigPDF, bumpPDF, bkgPDF), RooArgList(nSig, nBump, nBkg));  // not with Kai BF 
  
  cout <<"\n\n---------------  Fitting histo " <<MuMuPiKHist->GetTitle() <<"  ---------------\n" <<endl;
  
  totalPDF->fitTo(*MuMuPiKHist, Range(fitMin, fitMax), Extended(kTRUE));
  
  RooPlot* xframe = xVar.frame();
  //xframe->SetAxisRange(fitMin, xMax) ;
  xframe->SetTitle( hmyPsiPKPiMassSelAltZoom->GetTitle() ); xframe->SetYTitle( TString::Format("Candidates / %.1f MeV/c^{2}", bin_width*1000) ) ;
  
  MuMuPiKHist->plotOn(xframe);
  totalPDF->plotOn(xframe);
  //totalPDF->plotOn(xframe, Components(RooArgSet(sig1PDF)), LineColor(kRed));
  //totalPDF->plotOn(xframe, Components(RooArgSet(sig2PDF)), LineColor(kMagenta));
  totalPDF->plotOn(xframe, Components(RooArgSet(sigPDF)), LineColor(kRed));  // not with Kai BF 
  //totalPDF->plotOn(xframe, Components(RooArgSet(sig2PDF)), LineColor(kMagenta));
  //totalPDF->plotOn(xframe, Components(RooArgSet(bumpPDF)), LineColor(kGreen));  // not with Kai BF 
  totalPDF->plotOn(xframe, Components(RooArgSet(bkgPDF)), LineColor(kBlue), LineStyle(2));  // with Kai BF

  totalPDF->paramOn(xframe, Parameters(RooArgSet(nSig,mean,sigma1)), Layout(0.41,0.74,0.85)) ;
  /*
  if (fitMin < 4.8) {
    //totalPDF->paramOn(xframe, Parameters(RooArgSet(nSig1,mean1,sigma1)), Layout(0.41,0.74,0.85)) ;
    //totalPDF->paramOn(xframe, Parameters(RooArgSet(nSig1,nSig2,mean1,mean2,sigma1,sigma2)), Layout(0.41,0.74,0.92)) ;
    totalPDF->paramOn(xframe, Parameters(RooArgSet(nSig,mean1,sigma1,sigma2)), Layout(0.41,0.74,0.85)) ;
  } else if (fixedMassWind != "_5sigma") {
    //totalPDF->paramOn(xframe, Parameters(RooArgSet(nSig1,mean1,sigma1)), Layout(0.17,0.52,0.8)) ;
    //totalPDF->paramOn(xframe, Parameters(RooArgSet(nSig,mean1,sigma1,sigma2)), Layout(0.17,0.52,0.8)) ;
    //} else {
    totalPDF->paramOn(xframe, Parameters(RooArgSet(nSig,mean1,sigma1,sigma2)), Layout(0.6,0.95,0.9)) ;
  }
  */

  xframe->Draw();
  
  //Bpeak_C->SaveAs("/cmshome/cristella/work/Z_analysis/b_peak/plots/"+trigger+"_"+data+fixedMassWind+".png");
  //Bpeak_C->SaveAs("/cmshome/cristella/work/Z_analysis/b_peak/plots/twoGauss/"+trigger+"_"+data+fixedMassWind+".png");
  //Bpeak_C->SaveAs("/cmshome/cristella/work/Z_analysis/exclusive/CMSSW_5_3_7_patch5/src/UserCode/PsiPrimePiKPAT/test/selector/b_peak/plots/twoGauss/"+trigger+"_"+data+fixedMassWind+".png");
  //
  //Bpeak_C->SaveAs("/cmshome/cristella/work/Z_analysis/exclusive/CMSSW_5_3_7_patch5/src/UserCode/PsiPrimePiKPAT/test/selector/b_peak/plots/twoGauss/fromKstar.png");
  //Bpeak_C->SaveAs("/cmshome/cristella/work/Z_analysis/exclusive/CMSSW_5_3_7_patch5/src/UserCode/PsiPrimePiKPAT/test/selector/b_peak/plots/twoGauss/fromKstar_ex.png");
  //Bpeak_C->SaveAs("/cmshome/cristella/work/Z_analysis/exclusive/CMSSW_5_3_7_patch5/src/UserCode/PsiPrimePiKPAT/test/selector/b_peak/plots/twoGauss/Bplus.png");
  Bpeak_C->SaveAs("/plots/B0_mass.png");
}
