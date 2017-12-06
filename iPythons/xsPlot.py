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


# In[4]:


rootfile = "/Users/adrianodiflorio/Desktop/X4140_roots/mumukk_tree.root" 
inputfile = TFile(rootfile,"READ") 
xTuple = (inputfile.Get("outuple")) 


# In[6]:


myReader = TTreeReader("outuple", inputfile)
nentries = xTuple.GetEntries()
print nentries


# In[8]:


massbins = (6.0 - 4.0)/0.005
mass = RooRealVar("xM","M(#mu#muKK)[GeV]",4.0,6.0)
mass.setBins(400)
lxy = RooRealVar("xL","l(xy)",0.0,10000.)
hlt = RooRealVar("xHlt","xHlt",0.0,20.0)
masskk = RooRealVar("kkM","kkM",0.5,1.5)
masskk.setBins(int(200))
massmumu = RooRealVar("mumuM","mumuM",2.5,3.5)


# In[9]:


alldata = RooDataSet("alldata","alldata",xTuple,RooArgSet(masskk,mass,lxy,hlt,massmumu))#,cutFormula)


# In[10]:


alldata.numEntries()


# In[ ]:


xdataPrompt = (alldata.reduce('xM<4.8')).reduce('xM>4.0').reduce("xL<2.0").reduce("kkM<1.020+0.03").reduce("kkM>1.020-0.03")
xdataPrompt.numEntries()


# In[ ]:


a0 = RooRealVar("a0","a0",0.001,-1.,1.)
a1 = RooRealVar("a1","a1",0.001,-0.5,0.5)
a2 = RooRealVar("a2","a2",-0.00001,-2.,2.)
a3 = RooRealVar("a3","a3",0.0)#
a4 = RooRealVar("a4","a4",0.0,-0.1,0.1)
a5 = RooRealVar("a5","a5",0.0,-0.025,0.05)
a6 = RooRealVar("a6","a6",0.0,-0.001,0.001)

aset = RooArgList(a0,a1,a2)#,a3,a4,a5)
phimean = 1.020
sigma = RooRealVar("sigma","width of gaussian",0.0013)
gamma = RooRealVar("gamma","gamma of bw",0.004253)#,0.001,0.01)
mean = RooRealVar("mean","mean of gaussian",phimean,phimean-0.005,phimean+0.005);

nSig = RooRealVar("nSig","nSig",1E6,0.,5.0E6)
nBkg = RooRealVar("nBkg","nBkg",5E5,0.,5.0E6)
cheb = RooChebychev("cheb","Background",masskk,aset)
#gauss = RooGaussian("gauss","gaussian PDF ",mass,mean,sigma)
signal = RooVoigtian("signal","signal",masskk,mean,gamma,sigma)

tot = RooAddPdf("tot","g+cheb",RooArgList(signal,cheb),RooArgList(nSig,nBkg))


# In[ ]:


rPhifit = tot.fitTo(xdataPrompt,Range(1.020-0.03,1.020+0.03))


# In[ ]:


c = TCanvas("canvas","canvas",1200,800) 
phiFrame = masskk.frame(Range(1.020-0.03,1.020+0.03))
xdataPrompt.plotOn(phiFrame,RooLinkedList())
tot.plotOn(phiFrame)

phiFrame.Draw()
c.SaveAs("phiMassSPlot.png")
c.Clear()


# In[ ]:


cD=TCanvas("cD","cD",750,600);cD.cd()
splot   = RooStats.SPlot ( "sPlot","sPlot", b0dataNonPrompt, tot, RooArgList(nSig,nBkg))
dstree  = b0dataNonPrompt.store().tree()
dstree.GetEntryNumber(88)


# In[19]:


shist   = TH1F('shist','shist', 100, 1.00, 1.05)


# In[20]:


shist.Sumw2()
shist.SetLineColor(2)    
shist.SetMarkerColor(2); shist.SetMinimum(0.)
dstree.Project('shist','kkM','nSig_sw');  


# In[21]:


shist.Draw('e0');
cD.SaveAs('OtherPlotX.gif')

