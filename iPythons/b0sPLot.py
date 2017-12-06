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


# In[5]:


massmin = 5.2
massmax = 5.55


# In[6]:


myReader = TTreeReader("outuple", inputfile)
nentries = xTuple.GetEntries()
print nentries


# In[8]:


massbins = (6.0 - 4.0)/0.005
mass = RooRealVar("xM","M(#mu#muKK)[GeV]",4.0,6.0)
mass.setBins(int(massbins))
lxy = RooRealVar("xL","l(xy)",0.0,10000.)
hlt = RooRealVar("xHlt","xHlt",0.0,20.0)
masskk = RooRealVar("kkM","kkM",0.5,1.5)
massbins = 200
masskk.setBins(int(massbins))
massmumu = RooRealVar("mumuM","mumuM",2.5,3.5)
cutFormula = RooFormulaVar("cutFormula","cutFormula","xHlt!=8.0",RooArgList(hlt))


# In[9]:


alldata = RooDataSet("alldata","alldata",xTuple,RooArgSet(masskk,mass,lxy,hlt,massmumu))#,cutFormula)
datasetfile = TFile("xMassDataset.root","RECREATE") 
datasetfile.cd()
alldata.Write()


# In[10]:


alldata.numEntries()


# In[ ]:


#xb->setRange("alt","x_coarse_bin1,x_coarse_bin3,x_coarse_bin5,x_coarse_bin7,x_coarse_bin9") ;
b0dataNonPrompt = ((alldata.reduce('xHlt!=8')).reduce('xM>5.2')).reduce("xL>3.0")


# In[ ]:


b0dataNonPromptMass = b0dataNonPrompt.reduce(SelectVars(RooArgSet(mass)))
b0dataNonPrompt.numEntries()


# In[ ]:


c = TCanvas("canvas","canvas",1200,800) 
mass.setRange("fitRange",massmin,massmax)
mass.setBins(200)
massFrame = mass.frame(Range(massmin,massmax))
b0dataNonPrompt.plotOn(massFrame)
massFrame.Draw()
c.SaveAs("testmass.png")


# In[13]:


mean = RooRealVar("mean","mean of gaussian",5.38,5.31,5.41);
sigma1 = RooRealVar("sigma1","width of gaussian1",0.002,0.0005,0.05);
sigma2 = RooRealVar("sigma2","width of gaussian2",0.004,0.004,0.01);

a0 = RooRealVar("a0","a0",0.001,-1.,1.)
a1 = RooRealVar("a1","a1",0.001,-0.5,0.5)
a2 = RooRealVar("a2","a2",-0.00001,-2.,2.)
a3 = RooRealVar("a3","a3",-0.000001,-0.1,0.1)
a4 = RooRealVar("a4","a4",-0.000001,-2.,2.)
a5 = RooRealVar("a5","a5",-0.000001)
a6 = RooRealVar("a6","a6",-0.000001,-0.01,0.01)

aset = RooArgList(a0,a1,a2)#,a3)
gaussFrac = RooRealVar("sig1frac","fraction of component 1 in signal",0.3,0.0,1.0)
nSig = RooRealVar("nSig","nSig",100000,0.,10E6)
nBkg = RooRealVar("nBkg","nBkg",55000,0.,10E6)


# In[14]:


cheb = RooChebychev("cheb","Background",mass,aset)
gauss1 = RooGaussian("gauss1","gaussian PDF 1",mass,mean,sigma1)
gauss2 = RooGaussian("gauss2","gaussian PDF 2",mass,mean,sigma2)

model  = RooAddPdf("model","g1+g2",gauss1,gauss2,gaussFrac)
tot = RooAddPdf("tot","g1+g2+cheb",RooArgList(model,cheb),RooArgList(nSig,nBkg))


# In[ ]:





# In[15]:


rfit = tot.fitTo(b0dataNonPrompt,Range(massmin,massmax))


# In[16]:


massFrame = mass.frame(Range(massmin,massmax))
b0dataNonPrompt.plotOn(massFrame,RooLinkedList())
tot.plotOn(massFrame)

massFrame.Draw()
c.SaveAs("testmassFit.png")


# In[ ]:





# In[17]:


cD=TCanvas("cD","cD",750,600);cD.cd()
splot   = RooStats.SPlot ( "sPlot","sPlot", b0dataNonPrompt, tot, RooArgList(nSig,nBkg))


# In[18]:


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
cD.SaveAs('OtherPlot.gif')


# In[39]:


sys.exit()



