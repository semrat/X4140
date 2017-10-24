#!/usr/bin/env python
import sys, os, os.path, re
import commands,string,getopt

from math import *
#from ROOT import gROOT, TCanvas, TF1
from ROOT import TString
from ROOT import TFile, TCanvas
import ROOT

#ROOT.gROOT.SetBatch(True)


usage="\n Usage: python JPsivsrun.py <options> \n Options: \n --trigger= \t\t trigger path like: Jpsi6p5 Jpsi6p5_Barrel \n --period \t\t reco period like: CertMay6  May10 \n  --selection \t\t selection like: Barrel Barrel_nocuts Barrel_vtx Barrel_muon Barrel_trigger"
valid = ['trigger=','period=','selection=']
try:
    opts, args = getopt.getopt(sys.argv[1:], "", valid)
except getopt.GetoptError, ex:
    print usage
    sys.exit(1)

trigger = None
period = None
selection = None

selection = "0"
period = "0"
trigger = "0"

# for opt, arg in opts:
#     if opt == "--trigger":
#        trigger = arg
#     if opt == "--period":
#        period = arg
#     if opt == "--selection":
#        selection = arg
#
# if trigger == None:
#     print "--trigger option not provided"
#     print usage
#     sys.exit()
# if period == None:
#     print "--period option not provided"
#     print usage
#     sys.exit()
# if selection == None:
#     print "--selection option not provided"
#     print usage
#     sys.exit()


filename="JPsiCount.root"
fileIn = ROOT.TFile(filename)


#jsonfile1 = "/gpfs_data/local/cms/fanfani/Collision11/rootfiles/X3872/CMSSW422/Certified4Nov/all.json"
jsonfile1 = "/afs/cern.ch/work/a/adiflori/CMSSW_5_3_22/src/X4140/MuMuKKPAT/test/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON_MuonPhys.txt"


# dirIn  = fileIn.Get("JPsi")
#
# dirIn.ls()

#histo = dirIn.Get("JPsi_vs_run_%s"%selection) # read histogram filled with number of Jpsi vs run number

# read histogram filled with number of Jpsi vs run number
binmap = {}
strangeruns = {}

# ###############################
def tdrStyle():
  ROOT.gStyle.SetCanvasBorderMode(0);
  ROOT.gStyle.SetCanvasColor(0);
  ROOT.gStyle.SetCanvasDefH(600); #Height of canvas
  ROOT.gStyle.SetCanvasDefW(600); #Width of canvas
  ROOT.gStyle.SetCanvasDefX(0);   #POsition on screen
  ROOT.gStyle.SetCanvasDefY(0);
  # For the Pad:
  ROOT.gStyle.SetPadBorderMode(0);
  # ROOT.gStyle.SetPadBorderSize(Width_t size = 1);
  ROOT.gStyle.SetPadColor(0);
  ROOT.gStyle.SetPadGridX(False);
  ROOT.gStyle.SetPadGridY(False);
  ROOT.gStyle.SetGridColor(0);
  ROOT.gStyle.SetGridStyle(3);
  ROOT.gStyle.SetGridWidth(1);
  # For the frame:
  ROOT.gStyle.SetFrameBorderMode(0);
  ROOT.gStyle.SetFrameBorderSize(1);
  ROOT.gStyle.SetFrameFillColor(0);
  ROOT.gStyle.SetFrameFillStyle(0);
  ROOT.gStyle.SetFrameLineColor(1);
  ROOT.gStyle.SetFrameLineStyle(1);
  ROOT.gStyle.SetFrameLineWidth(1);
  # For the histo:
  ROOT.gStyle.SetHistLineColor(1);
  ROOT.gStyle.SetHistLineStyle(0);
  ROOT.gStyle.SetHistLineWidth(1);
  ROOT.gStyle.SetEndErrorSize(2);
  ROOT.gStyle.SetMarkerStyle(20);
  #For the date:
  ROOT.gStyle.SetOptDate(0);
  # For the statistics box:
  ROOT.gStyle.SetOptFile(0);
  #ROOT.gStyle.SetOptStat(0);
  ROOT.gStyle.SetOptStat("eou");
  ROOT.gStyle.SetStatColor(0);
  ROOT.gStyle.SetStatFont(42);
  ROOT.gStyle.SetStatFontSize(0.04);#/---> ROOT.gStyle.SetStatFontSize(0.025);
  ROOT.gStyle.SetStatTextColor(1);
  ROOT.gStyle.SetStatFormat("6.4g");
  ROOT.gStyle.SetStatBorderSize(1);
  ROOT.gStyle.SetStatH(0.15);
  ROOT.gStyle.SetStatW(0.3);#/---> ROOT.gStyle.SetStatW(0.15);
  # Margins:
  ROOT.gStyle.SetPadTopMargin(0.05);
  ROOT.gStyle.SetPadBottomMargin(0.13);
  ROOT.gStyle.SetPadLeftMargin(0.16);
  ROOT.gStyle.SetPadRightMargin(0.04);
  # For the Global title:
  ROOT.gStyle.SetOptTitle(0);
  # For the axis titles:
  ROOT.gStyle.SetTitleColor(1, "XYZ");
  ROOT.gStyle.SetTitleFont(42, "XYZ");
  ROOT.gStyle.SetTitleSize(0.06, "XYZ");
  ROOT.gStyle.SetTitleXOffset(0.9);
  ROOT.gStyle.SetTitleYOffset(1.25);
  # For the axis labels:
  ROOT.gStyle.SetLabelColor(1, "XYZ");
  ROOT.gStyle.SetLabelFont(42, "XYZ");
  ROOT.gStyle.SetLabelOffset(0.007, "XYZ");
  ROOT.gStyle.SetLabelSize(0.05, "XYZ");
  # For the axis:
  ROOT.gStyle.SetAxisColor(1, "XYZ");
  ROOT.gStyle.SetStripDecimals(True);
  ROOT.gStyle.SetTickLength(0.03, "XYZ");
  ROOT.gStyle.SetNdivisions(510, "XYZ");
  ROOT.gStyle.SetPadTickX(1);  # To get tick marks on the opposite side of the frame
  ROOT.gStyle.SetPadTickY(1);
  # Postscript options:
  ROOT.gStyle.SetPaperSize(20.,20.);
  ## OVERRIDES
  ROOT.gStyle.SetHistMinimumZero(1)
  ROOT.gStyle.SetErrorX(0.5);
  ROOT.gStyle.SetOptStat(False)
  # Done
  ROOT.gROOT.ForceStyle();

# ###############################
def testLumiEnv():
  lumi_cmd='lumiCalc2.py -h'
  status_cmd,lumiout=commands.getstatusoutput(lumi_cmd)
  if status_cmd != 0 :
         print "\n ERROR with Environment for lumiCalc.py ==> %s \n"%lumiout
         sys.exit(1)

# ###############################
def ReBinRunHisto(histo,histoname):

  ## find the new binning based on runs with muons
  nbinX = histo.GetNbinsX();
  nbinOutX =0
  first = True
  bin_first = int(round(histo.GetBinCenter(1),0));
  bin_last = int(round(histo.GetBinCenter(nbinX),0));
  for ix in range(nbinX):
    run = round(histo.GetBinCenter(ix),0)
    irun = int(run)
    if histo.GetBinContent(ix)>0 :  # skip runs without JPsi events
       first = False
       if (first):
           bin_first = int(round(histo.GetBinCenter(ix),0));
       bin_last = int(round(histo.GetBinCenter(ix),0));
       nbinOutX = nbinOutX + 1
       binmap[ix]=histo.GetBinContent(ix)

  print 'First Run= %s  Last Run= %s'%(bin_first,bin_last);
  histoOut = ROOT.TH1F(histoname, 'Nb. of JPsi vs run', nbinOutX,bin_first,bin_last )
  ixnew = 1
  ## fill the histogram with the new binning
  #for ix in range(nbinX):
  for ix, val in binmap.iteritems():
    if val>0 : # skip runs without JPsi events
       run = round(histo.GetBinCenter(ix),0)
       irun = int(run)
       #    print "====> run %i is now in bin %s"%(irun,ixnew)
       runcontent = val #histo.GetBinContent(ix)
       histoOut.SetBinContent(ixnew,runcontent)
       histoOut.SetBinError(ixnew,histo.GetBinError(ix))
       histoOut.GetXaxis().SetBinLabel(ixnew,str(irun))
       ixnew = ixnew + 1
  return histoOut

# ###############################
def GetMuPerLumiHisto(histo,hlt):
  nbinX = histo.GetNbinsX();
  histo.GetNbinsX();
#
  histo.Sumw2();
#
  histoOut = histo.Clone("lumiJPsirate" + hlt);


  #lumi_cmd='lumiCalc2.py -i '+jsonfile1+' overview > all_lumi_nobkg.txt'
  lumi_cmd='brilcalc lumi -i '+jsonfile1+' > all_lumi_nobkg.txt'
  lumi_cmd="echo"
  print(lumi_cmd)
  status_cmd,lumiout=commands.getstatusoutput(lumi_cmd)
  if status_cmd != 0 :
     print "ERROR %s"%lumiout
     sys.exit(2)


  #for ix in range(1,nbinX):
  for ix, val in binmap.iteritems():
   nmuperlumi = 0.
   lumi=999999.
   #run = round(histo.GetBinCenter(ix),0)
   run = round(histo.GetBinCenter(ix),0)
   irun = int(run)
   #if histo.GetBinContent(ix)>0 : # for runs with JPsi
   if val>0:
      #runcontent = histo.GetBinContent(ix)
      runcontent = val
      lumi_cmd='cat brilCalcOutput.txt | grep '+str(irun)+': | grep -v Warning | grep -v WARNING | cut -d"|" -f7'
      status_cmd,lumiout=commands.getstatusoutput(lumi_cmd)
      #print(irun)
      #print(lumiout)
      lumiclean = string.strip(lumiout)
      lumiclean = lumiout.strip(" ");
      lumiclean = lumiclean.strip("\n");
      lumiclean = lumiclean.strip("\b");
      lumiclean = lumiclean.strip("\t");
      #
      #print(lumiclean)
      #print(status_cmd)
      if status_cmd != 0:
         print "ERROR %s for run %i"%(lumiout,irun)
      else:
         #print(lumiout)
         #lumiout=
         #lumituple=lumiout.partition("(/")
         #lumiunit=lumituple[2]
         #lumi = float(string.strip(lumituple[0]))
         lumi =float(string.strip(lumiclean)) / 1000.
         #if (lumiunit=="nb)"):
           # lumi=1000.*lumi
         #if (lumiunit=="pb)"):
           # lumi=1000000.*lumi
         nmuperlumi = runcontent/lumi
         err=sqrt(runcontent)/lumi
        #  print(lumi)
        #  print(runcontent)
        #  print(nmuperlumi)
         if nmuperlumi < 1.5 or nmuperlumi > 2.5 and hlt=="":
             strangeruns[hlt].append(irun)
         if nmuperlumi < 0.19 or nmuperlumi > 0.41 and hlt=="_HLT4":
             strangeruns[hlt].append(irun)
         if nmuperlumi < 0.4 or nmuperlumi > 0.65and hlt=="_HLT8":
             strangeruns[hlt].append(irun)
         print("run=%i #JPsi=%i lumi=%s  #JPsi per lumi=%s +/- %s"%(irun, runcontent,lumi,nmuperlumi,err))
        # except:
        #  try:
        #
        #   print "ERROR for run=%i "%(irun)
        #   nmuperlumi = 0
   #histoOut.GetXaxis().SetBinLabel(ix,str(irun))
   histoOut.SetBinContent(ix,nmuperlumi)
   histoOut.SetBinError(ix,histo.GetBinError(ix)/lumi)
  return histoOut
# ###############################
# ###############################
def printRunHisto(histo,title):
    histo.SetXTitle("Run")
    histo.SetYTitle(title)
    histo.SetMarkerStyle(20)
    histo.Draw("E1")


################################
#
################################
if __name__ == '__main__':
  ROOT.gBenchmark.Start('run info')
  tdrStyle()
  hfile = TFile( 'run_%s_%s_%s.root'%(selection,period,trigger), 'RECREATE', 'ROOT file with histograms' )

  hlts = ["","_HLT4","_HLT8"]

  totCanvas = TCanvas("totCanvas","totCanvas",1200, 900)

  for hlt in hlts:
      strangeruns[hlt] = []
      fileIn.cd()
      histo = fileIn.Get("JPsi_vs_run" + hlt)
      ## #mu vs run
      #c1 = TCanvas("c1","c1",1200, 900)
      #c1.SetLogy(1)
      #printRunHisto(histo,'Nb. JPsi ')
      #c1.SaveAs("JPsivsrun%s.eps"%(hlt))


      #c2 = TCanvas("c2","c2",1200, 900)
      #c2.SetLogy(1)
      histo_rebin=ReBinRunHisto(histo,'JPsi_rebin' + hlt)
      #printRunHisto(histo_rebin,'Nb. JPsi ')
      #c2.SaveAs("JPsivsrunrebin%s.eps"%(trigger))

      #testLumiEnv()
      histoout=GetMuPerLumiHisto(histo,hlt)
      #cl = TCanvas("cl","cl",1200, 900)
      #cl.SetLogy(0)
      #printRunHisto(histoout,'# JPsi (from X)/lumi (mb)' + hlt)
      #cl.SaveAs("JPsivsrun_lumi%s.eps"%(hlt))

      #clr = TCanvas("clr","clr",1200, 900)
      #clr.SetLogy(0)
      histoout_rebin=ReBinRunHisto(histoout,'JPsilumi_rebin' + hlt)
      #printRunHisto(histoout_rebin,'# JPsi (from X)/lumi (mb)')
      #clr.SaveAs("JPsivsrunrebin_lumi%s.eps"%(hlt))

      hfile.cd()
      histo.Write()
      histo_rebin.Write()
      histoout.Write()
      histoout_rebin.Write()

      totCanvas.cd()
      if hlt=="":
        histoout_rebin.Draw()
      else:
        histoout_rebin.Draw("same")

  totCanvas.SaveAs("JPsivsrunrebin_lumi_tot.eps")
  ROOT.gBenchmark.Show( 'run info' )
  hfile.Write()

  print("Strange runs :")
  print(strangeruns)
