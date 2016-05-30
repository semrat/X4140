//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Nov  3 15:18:57 2014 by ROOT version 5.34/10
// from TTree Z_data/Z(4430) Data
// found on file: MuOniaRun2012D_13Jul_MuMuPiPiPAT_ntpl.root
//////////////////////////////////////////////////////////

#ifndef psiPrimePiK_MC_h
#define psiPrimePiK_MC_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

#include <TSystem.h>
#include <TTree.h>
#include <TBranch.h>
#include <TCint.h>
#include <TRandom.h>
#include <TMath.h>
#include <TDirectory.h>
#include "TEnv.h"
#include <TString.h>
#include <TProof.h>
#include <TProofOutputFile.h>
#include <TLorentzVector.h>
#include "TPoint.h"
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
//#include <tuple> //  must be enabled with the -std=c++11 or -std=gnu++11 compiler options
/*
////////// RooFit //////////
//gInterpreter->AddIncludePath("/cvmfs/cms.cern.ch/slc6_amd64_gcc481/lcg/roofit/5.34.09-cms5/include") ;
//gROOT->ProcessLine(".include /cvmfs/cms.cern.ch/slc6_amd64_gcc481/lcg/roofit/5.34.09-cms5/include") ;
// one of the 2 above lines need to be executed interactively before following includes:
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooGaussian.h>
#include <RooConstVar.h>
#include <RooChebychev.h>
#include <RooAddPdf.h>
#include <RooAbsReal.h>
#include <RooPlot.h>
*/
// Header file for the classes stored in the TTree if any.



// Fixed size dimensions of array or collections stored in the TTree if any.

class psiPrimePiK_MC : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   ////////// Declaration of leaf types //////////
   vector<unsigned int> *TrigRes;
   vector<string>  *TrigNames;
   vector<string>  *MatchTriggerNames;
   vector<unsigned int> *L1TrigRes;
   UInt_t          evtNum;
   UInt_t          runNum;
   UInt_t          lumiNum;
   Float_t         priVtx_X;
   Float_t         priVtx_Y;
   Float_t         priVtx_Z;
   Float_t         priVtx_XE;
   Float_t         priVtx_YE;
   Float_t         priVtx_ZE;
   Float_t         priVtx_NormChi2;
   Float_t         priVtx_Chi2;
   Float_t         priVtx_CL;
   UInt_t          priVtx_tracks;
   // MC Analysis
   UInt_t          nMCAll ;
   vector<int>     *MCPdgIdAll;
   UInt_t          nMCB0 ;
   vector<float>   *PriVtxGen_X ;
   vector<float>   *PriVtxGen_Y ;
   vector<float>   *PriVtxGen_Z ;
   vector<float>   *PriVtxGen_EX ;
   vector<float>   *PriVtxGen_EY ;
   vector<float>   *PriVtxGen_EZ ;
   vector<float>   *PriVtxGen_Chi2 ;
   vector<float>   *PriVtxGen_CL ;
   vector<int>     *PriVtxGen_Ndof ;
   vector<int>     *PriVtxGen_tracks ;
   //
   vector<float>   *genEvtVtx_X ;
   vector<float>   *genEvtVtx_Y ;
   vector<float>   *genEvtVtx_Z ;
   vector<int>     *genEvtVtx_particlesOut ;
   //
   vector<float>   *MCpsi2SPx ;
   vector<float>   *MCpsi2SPy ;
   vector<float>   *MCpsi2SPz ;
   vector<float>   *MCmupPx ;
   vector<float>   *MCmupPy ;
   vector<float>   *MCmupPz ;
   vector<float>   *MCmumPx ;
   vector<float>   *MCmumPy ;
   vector<float>   *MCmumPz ;
   vector<float>   *MCpionPx ;
   vector<float>   *MCpionPy ;
   vector<float>   *MCpionPz ;
   vector<float>   *MCkaonPx ;
   vector<float>   *MCkaonPy ;
   vector<float>   *MCkaonPz ;
   vector<float>   *MCPx;
   vector<float>   *MCPy;
   vector<float>   *MCPz;
   vector<int>     *MCpionCh;
   vector<int>     *MCkaonCh;

   // Generic Muons
   UInt_t          nMu;
   vector<float>   *muPx;
   vector<float>   *muPy;
   vector<float>   *muPz;
   vector<float>   *muD0;
   vector<float>   *muDz;
   vector<float>   *muChi2;
   vector<int>     *muNDF;
   vector<int>     *muPhits;
   vector<int>     *muShits;
   vector<int>     *muLayersTr;
   vector<int>     *muLayersPix;
   vector<float>   *muD0E;
   vector<float>   *muDzVtxErr;
   vector<unsigned int> *muKey;
   vector<int>     *muGlMuHits;
   vector<float>   *muGlChi2;
   vector<int>     *muGlNDF;
   vector<int>     *muGlMatchedStation;
   vector<float>   *muGlDzVtx;
   vector<float>   *muGlDxyVtx;
   vector<int>     *muType;
   vector<int>     *muQual;
   vector<int>     *muTrack;
   vector<float>   *muCharge;
   vector<float>   *mufHits;
   vector<bool>    *muFirstBarrel;
   vector<bool>    *muFirstEndCap;
   vector<float>   *muDzVtx;
   vector<float>   *muDxyVtx;
   vector<float>   *trNotRef;
   vector<float>   *trRef;
   vector<float>   *trackPx;
   vector<float>   *trackPy;
   vector<float>   *trackPz;
   vector<float>   *trackEnergy;
   vector<int>     *trackNDF;
   vector<int>     *trackPhits;
   vector<int>     *trackShits;
   vector<float>   *trackChi2;
   vector<float>   *trackD0;
   vector<float>   *trackD0Err;
   vector<float>   *trackCharge;
   vector<int>     *TrackHighPurity;
   vector<int>     *TrackTight;
   vector<float>   *trackfHits;
   vector<bool>    *trackFirstBarrel;
   vector<bool>    *trackFirstEndCap;
   vector<float>   *trackDzVtx;
   vector<float>   *trackDxyVtx;
   vector<double>  *tr_nsigdedx;
   vector<float>   *tr_dedx;
   vector<float>   *tr_dedxMass;
   vector<float>   *tr_theo;
   vector<float>   *tr_sigma;
   vector<float>   *tr_dedx_byHits;
   vector<float>   *tr_dedxErr_byHits;
   vector<int>     *tr_saturMeas_byHits;
   vector<int>     *tr_Meas_byHits;
   UInt_t          nMuMu;
   vector<float>   *MuMuMass;
   vector<float>   *MuMuPx;
   vector<float>   *MuMuPy;
   vector<float>   *MuMuPz;
   vector<float>   *MuMuVtx_CL;
   vector<float>   *MuMuVtx_Chi2;
   vector<float>   *MuMuDecayVtx_X;
   vector<float>   *MuMuDecayVtx_Y;
   vector<float>   *MuMuDecayVtx_Z;
   vector<float>   *MuMuDecayVtx_XE;
   vector<float>   *MuMuDecayVtx_YE;
   vector<float>   *MuMuDecayVtx_ZE;
   vector<int>     *mu1Idx;
   vector<int>     *mu2Idx;
   vector<float>   *mu1Px_MuMu;
   vector<float>   *mu1Py_MuMu;
   vector<float>   *mu1Pz_MuMu;
   vector<float>   *mu1Chi_MuMu2;
   vector<int>     *mu1NDF_MuMu;
   vector<float>   *mu2Px_MuMu;
   vector<float>   *mu2Py_MuMu;
   vector<float>   *mu2Pz_MuMu;
   vector<float>   *mu2Chi2_MuMu;
   vector<int>     *mu2NDF_MuMu;
   vector<int>     *MuMuType;
   vector<bool>    *MuMuMuonTrigMatch;
   vector<float>   *PriVtxMuMuCorr_X;
   vector<float>   *PriVtxMuMuCorr_Y;
   vector<float>   *PriVtxMuMuCorr_Z;
   vector<double>  *PriVtxMuMuCorr_EX;
   vector<double>  *PriVtxMuMuCorr_EY;
   vector<double>  *PriVtxMuMuCorr_EZ;
   vector<float>   *PriVtxMuMuCorr_Chi2;
   vector<float>   *PriVtxMuMuCorr_CL;
   vector<int>     *PriVtxMuMuCorr_tracks;
   UInt_t          nB0;
   vector<float>   *B0Mass;
   vector<float>   *B0Px;
   vector<float>   *B0Py;
   vector<float>   *B0Pz;
   vector<double>  *B0PxE;
   vector<double>  *B0PyE;
   vector<double>  *B0PzE;
   vector<float>   *B0Vtx_CL;
   vector<float>   *B0Vtx_Chi2;
   vector<float>   *B0DecayVtx_X;
   vector<float>   *B0DecayVtx_Y;
   vector<float>   *B0DecayVtx_Z;
   vector<double>  *B0DecayVtx_XE;
   vector<double>  *B0DecayVtx_YE;
   vector<double>  *B0DecayVtx_ZE;
   // Primary Vertex with largest B0_cos(alpha)
   vector<float>   *PriVtx_B0CosAlpha_X, *PriVtx_B0CosAlpha_Y, *PriVtx_B0CosAlpha_Z ; 
   vector<double>  *PriVtx_B0CosAlpha_EX, *PriVtx_B0CosAlpha_EY, *PriVtx_B0CosAlpha_EZ ;
   vector<float>   *PriVtx_B0CosAlpha_Chi2, *PriVtx_B0CosAlpha_CL;
   vector<int>     *PriVtx_B0CosAlpha_tracks ;
   vector<float>   *B0LessPV_tracksPtSq, *B0LessPV_4tracksPtSq ;
   vector<float>   *PriVtxB0Less_X, *PriVtxB0Less_Y, *PriVtxB0Less_Z ; 
   vector<double>  *PriVtxB0Less_EX, *PriVtxB0Less_EY, *PriVtxB0Less_EZ ;
   vector<float>   *PriVtxB0Less_Chi2, *PriVtxB0Less_CL;
   vector<int>     *PriVtxB0Less_tracks ;
   vector<float>   *PriVtxB0Less_B0CosAlpha_X, *PriVtxB0Less_B0CosAlpha_Y, *PriVtxB0Less_B0CosAlpha_Z ; 
   vector<double>  *PriVtxB0Less_B0CosAlpha_EX, *PriVtxB0Less_B0CosAlpha_EY, *PriVtxB0Less_B0CosAlpha_EZ ;
   vector<float>   *PriVtxB0Less_B0CosAlpha_Chi2, *PriVtxB0Less_B0CosAlpha_CL;
   vector<int>     *PriVtxB0Less_B0CosAlpha_tracks ;
   // Primary Vertex with largest B0_cos(alpha) 3D
   vector<float>   *PriVtx_B0CosAlpha3D_X, *PriVtx_B0CosAlpha3D_Y, *PriVtx_B0CosAlpha3D_Z ; 
   vector<double>  *PriVtx_B0CosAlpha3D_EX, *PriVtx_B0CosAlpha3D_EY, *PriVtx_B0CosAlpha3D_EZ ;
   vector<float>   *PriVtx_B0CosAlpha3D_Chi2, *PriVtx_B0CosAlpha3D_CL;
   vector<int>     *PriVtx_B0CosAlpha3D_tracks ;
   vector<float>   *PriVtxB0Less_B0CosAlpha3D_X, *PriVtxB0Less_B0CosAlpha3D_Y, *PriVtxB0Less_B0CosAlpha3D_Z ; 
   vector<double>  *PriVtxB0Less_B0CosAlpha3D_EX, *PriVtxB0Less_B0CosAlpha3D_EY, *PriVtxB0Less_B0CosAlpha3D_EZ ;
   vector<float>   *PriVtxB0Less_B0CosAlpha3D_Chi2, *PriVtxB0Less_B0CosAlpha3D_CL;
   vector<int>     *PriVtxB0Less_B0CosAlpha3D_tracks ;
   // Primary Vertex with "B0 correction" // check name
   vector<float>   *PriVtxB0Corr_X, *PriVtxB0Corr_Y, *PriVtxB0Corr_Z;
   vector<double>  *PriVtxB0Corr_EX, *PriVtxB0Corr_EY, *PriVtxB0Corr_EZ;
   vector<float>   *PriVtxB0Corr_Chi2, *PriVtxB0Corr_CL;
   vector<int>     *PriVtxB0Corr_tracks;
   // Lifetimes
   vector<double>  *B0CosAlphaBS, *B0CosAlpha3DBS, *B0CTauBS, *B0CTauBSE, *B0LxyBS ;
   vector<double>  *B0CosAlphaPV, *B0CosAlpha3DPV, *B0CTauPV, *B0CTauPVE, *B0LxyPV ;
   vector<double>  *B0CosAlphaPVCosAlpha, *B0CTauPVCosAlpha, *B0CTauPVCosAlphaE, *B0LxyPVCosAlpha ;
   vector<double>  *B0CosAlphaPVCosAlpha3D, *B0CTauPVCosAlpha3D, *B0CTauPVCosAlpha3DE, *B0LxyPVCosAlpha3D ;
   vector<double>  *B0CosAlphaB0LessPV, *B0CosAlpha3DB0LessPV, *B0CTauB0LessPV, *B0CTauB0LessPVE, *B0LxyB0LessPV ;
   vector<double>  *B0CosAlphaB0LessPVCosAlpha, *B0CTauB0LessPVCosAlpha, *B0CTauB0LessPVCosAlphaE, *B0LxyB0LessPVCosAlpha ;
   vector<double>  *B0CosAlphaB0LessPVCosAlpha3D, *B0CTauB0LessPVCosAlpha3D, *B0CTauB0LessPVCosAlpha3DE, *B0LxyB0LessPVCosAlpha3D ;
   vector<double>  *B0CosAlphaPVX, *B0CTauPVX, *B0CTauPVXE, *B0LxyPVX ;
   vector<float>   *B0CTauPVX_3D, *B0CTauPVX_3D_err ;
   vector<int>     *B0MuMuIdx, *B0PionIdx, *B0KaonIdx ;
   vector<float>   *PiPiMass_err;
   vector<float>   *Muon1Px_MuMuPiK, *Muon1Py_MuMuPiK, *Muon1Pz_MuMuPiK, *Muon1E_MuMuPiK;
   vector<float>   *Muon2Px_MuMuPiK, *Muon2Py_MuMuPiK, *Muon2Pz_MuMuPiK, *Muon2E_MuMuPiK;
   vector<float>   *PionPx_MuMuPiK, *PionPy_MuMuPiK, *PionPz_MuMuPiK, *PionE_MuMuPiK;
   vector<double>  *pion_nsigdedx;
   vector<float>   *pion_dedx, *pion_dedxMass, *pion_theo, *pion_sigma;
   vector<float>   *pion_dedx_byHits, *pion_dedxErr_byHits;
   vector<int>     *pion_saturMeas_byHits, *pion_Meas_byHits;
   vector<float>   *KaonPx_MuMuPiK, *KaonPy_MuMuPiK, *KaonPz_MuMuPiK, *KaonE_MuMuPiK;
   vector<double>  *kaon_nsigdedx;
   vector<float>   *kaon_dedx, *kaon_dedxMass, *kaon_theo, *kaon_sigma;
   vector<float>   *kaon_dedx_byHits, *kaon_dedxErr_byHits;
   vector<int>     *kaon_saturMeas_byHits, *kaon_Meas_byHits;

   // List of branches
   TBranch        *b_TrigRes;   //!
   TBranch        *b_TrigNames;   //!
   TBranch        *b_MatchTriggerNames;   //!
   TBranch        *b_L1TrigRes;   //!
   TBranch        *b_evtNum;   //!
   TBranch        *b_runNum;   //!
   TBranch        *b_lumiNum;   //!
   TBranch        *b_priVtx_X;   //!
   TBranch        *b_priVtx_Y;   //!
   TBranch        *b_priVtx_Z;   //!
   TBranch        *b_priVtx_XE;   //!
   TBranch        *b_priVtx_YE;   //!
   TBranch        *b_priVtx_ZE;   //!
   TBranch        *b_priVtx_NormChi2;   //!
   TBranch        *b_priVtx_Chi2;   //!
   TBranch        *b_priVtx_CL;   //!
   TBranch        *b_priVtx_tracks;   //!
   // MC Analysis
   TBranch        *b_nMCAll ;
   TBranch        *b_MCPdgIdAll ;
   TBranch        *b_nMCB0 ;
   TBranch        *b_PriVtxGen_X ;
   TBranch        *b_PriVtxGen_Y ;
   TBranch        *b_PriVtxGen_Z ;
   TBranch        *b_PriVtxGen_EX ;
   TBranch        *b_PriVtxGen_EY ;
   TBranch        *b_PriVtxGen_EZ ;
   TBranch        *b_PriVtxGen_Chi2 ;
   TBranch        *b_PriVtxGen_CL ;
   TBranch        *b_PriVtxGen_Ndof ;
   TBranch        *b_PriVtxGen_tracks ;
   TBranch        *b_genEvtVtx_X ;
   TBranch        *b_genEvtVtx_Y ;
   TBranch        *b_genEvtVtx_Z ;
   TBranch        *b_genEvtVtx_particlesOut ;
   TBranch        *b_MCpsi2SPx ;
   TBranch        *b_MCpsi2SPy ;
   TBranch        *b_MCpsi2SPz ;
   TBranch        *b_MCmupPx ;
   TBranch        *b_MCmupPy ;
   TBranch        *b_MCmupPz ;
   TBranch        *b_MCmumPx ;
   TBranch        *b_MCmumPy ;
   TBranch        *b_MCmumPz ;
   TBranch        *b_MCpionPx ;
   TBranch        *b_MCpionPy ;
   TBranch        *b_MCpionPz ;
   TBranch        *b_MCkaonPx ;
   TBranch        *b_MCkaonPy ;
   TBranch        *b_MCkaonPz ;
   TBranch        *b_MCPx;   
   TBranch        *b_MCPy;   
   TBranch        *b_MCPz;   
   TBranch        *b_MCpionCh;
   TBranch        *b_MCkaonCh;

   // Generic Muons
   TBranch        *b_nMu;   //!
   TBranch        *b_muPx;   //!
   TBranch        *b_muPy;   //!
   TBranch        *b_muPz;   //!
   TBranch        *b_muD0;   //!
   TBranch        *b_muDz;   //!
   TBranch        *b_muChi2;   //!
   TBranch        *b_muNDF;   //!
   TBranch        *b_muPhits;   //!
   TBranch        *b_muShits;   //!
   TBranch        *b_muLayersTr;   //!
   TBranch        *b_muLayersPix;   //!
   TBranch        *b_muD0E;   //!
   TBranch        *b_muDzVtxErr;   //!
   TBranch        *b_muKey;   //!
   TBranch        *b_muGlMuHits;   //!
   TBranch        *b_muGlChi2;   //!
   TBranch        *b_muGlNDF;   //!
   TBranch        *b_muGlMatchedStation;   //!
   TBranch        *b_muGlDzVtx;   //!
   TBranch        *b_muGlDxyVtx;   //!
   TBranch        *b_muType;   //!
   TBranch        *b_muQual;   //!
   TBranch        *b_muTrack;   //!
   TBranch        *b_muCharge;   //!
   TBranch        *b_mufHits;   //!
   TBranch        *b_muFirstBarrel;   //!
   TBranch        *b_muFirstEndCap;   //!
   TBranch        *b_muDzVtx;   //!
   TBranch        *b_muDxyVtx;   //!
   TBranch        *b_trNotRef;   //!
   TBranch        *b_trRef;   //!
   TBranch        *b_trackPx;   //!
   TBranch        *b_trackPy;   //!
   TBranch        *b_trackPz;   //!
   TBranch        *b_trackEnergy;   //!
   TBranch        *b_trackNDF;   //!
   TBranch        *b_trackPhits;   //!
   TBranch        *b_trackShits;   //!
   TBranch        *b_trackChi2;   //!
   TBranch        *b_trackD0;   //!
   TBranch        *b_trackD0Err;   //!
   TBranch        *b_trackCharge;   //!
   TBranch        *b_TrackHighPurity;   //!
   TBranch        *b_TrackTight;   //!
   TBranch        *b_trackfHits;   //!
   TBranch        *b_trackFirstBarrel;   //!
   TBranch        *b_trackFirstEndCap;   //!
   TBranch        *b_trackDzVtx;   //!
   TBranch        *b_trackDxyVtx;   //!
   TBranch        *b_tr_nsigdedx;   //!
   TBranch        *b_tr_dedx;   //!
   TBranch        *b_tr_dedxMass;   //!
   TBranch        *b_tr_theo;   //!
   TBranch        *b_tr_sigma;   //!
   TBranch        *b_tr_dedx_byHits;   //!
   TBranch        *b_tr_dedxErr_byHits;   //!
   TBranch        *b_tr_saturMeas_byHits;   //!
   TBranch        *b_tr_Meas_byHits;   //!
   TBranch        *b_nMuMu;   //!
   TBranch        *b_MuMuMass;   //!
   TBranch        *b_MuMuPx;   //!
   TBranch        *b_MuMuPy;   //!
   TBranch        *b_MuMuPz;   //!
   TBranch        *b_MuMuVtx_CL;   //!
   TBranch        *b_MuMuVtx_Chi2;   //!
   TBranch        *b_MuMuDecayVtx_X;   //!
   TBranch        *b_MuMuDecayVtx_Y;   //!
   TBranch        *b_MuMuDecayVtx_Z;   //!
   TBranch        *b_MuMuDecayVtx_XE;   //!
   TBranch        *b_MuMuDecayVtx_YE;   //!
   TBranch        *b_MuMuDecayVtx_ZE;   //!
   TBranch        *b_mu1Idx;   //!
   TBranch        *b_mu2Idx;   //!
   TBranch        *b_mu1Px_MuMu;   //!
   TBranch        *b_mu1Py_MuMu;   //!
   TBranch        *b_mu1Pz_MuMu;   //!
   TBranch        *b_mu1Chi_MuMu2;   //!
   TBranch        *b_mu1NDF_MuMu;   //!
   TBranch        *b_mu2Px_MuMu;   //!
   TBranch        *b_mu2Py_MuMu;   //!
   TBranch        *b_mu2Pz_MuMu;   //!
   TBranch        *b_mu2Chi2_MuMu;   //!
   TBranch        *b_mu2NDF_MuMu;   //!
   TBranch        *b_MuMuType;   //!
   TBranch        *b_MuMuMuonTrigMatch;   //!
   TBranch        *b_PriVtxMuMuCorr_X;   //!
   TBranch        *b_PriVtxMuMuCorr_Y;   //!
   TBranch        *b_PriVtxMuMuCorr_Z;   //!
   TBranch        *b_PriVtxMuMuCorr_EX;   //!
   TBranch        *b_PriVtxMuMuCorr_EY;   //!
   TBranch        *b_PriVtxMuMuCorr_EZ;   //!
   TBranch        *b_PriVtxMuMuCorr_Chi2;   //!
   TBranch        *b_PriVtxMuMuCorr_CL;   //!
   TBranch        *b_PriVtxMuMuCorr_tracks;   //!
   TBranch        *b_nB0;   //!
   TBranch        *b_B0Mass;   //!
   TBranch        *b_B0Px;   //!
   TBranch        *b_B0Py;   //!
   TBranch        *b_B0Pz;   //!
   TBranch        *b_B0PxE;   //!
   TBranch        *b_B0PyE;   //!
   TBranch        *b_B0PzE;   //!
   TBranch        *b_B0Vtx_CL;   //!
   TBranch        *b_B0Vtx_Chi2;   //!
   TBranch        *b_B0DecayVtx_X;   //!
   TBranch        *b_B0DecayVtx_Y;   //!
   TBranch        *b_B0DecayVtx_Z;   //!
   TBranch        *b_B0DecayVtx_XE;   //!
   TBranch        *b_B0DecayVtx_YE;   //!
   TBranch        *b_B0DecayVtx_ZE;   //!
   // Primary Vertex with largest B0_cos(alpha)
   TBranch        *b_PriVtx_B0CosAlpha_X, *b_PriVtx_B0CosAlpha_Y, *b_PriVtx_B0CosAlpha_Z ; 
   TBranch        *b_PriVtx_B0CosAlpha_EX, *b_PriVtx_B0CosAlpha_EY, *b_PriVtx_B0CosAlpha_EZ ;
   TBranch	  *b_PriVtx_B0CosAlpha_Chi2, *b_PriVtx_B0CosAlpha_CL;
   TBranch        *b_PriVtx_B0CosAlpha_tracks ;
   TBranch        *b_B0LessPV_tracksPtSq, *b_B0LessPV_4tracksPtSq ;
   TBranch        *b_PriVtxB0Less_X, *b_PriVtxB0Less_Y, *b_PriVtxB0Less_Z ; 
   TBranch        *b_PriVtxB0Less_EX, *b_PriVtxB0Less_EY, *b_PriVtxB0Less_EZ ;
   TBranch	  *b_PriVtxB0Less_Chi2, *b_PriVtxB0Less_CL;
   TBranch        *b_PriVtxB0Less_tracks ;
   TBranch        *b_PriVtxB0Less_B0CosAlpha_X, *b_PriVtxB0Less_B0CosAlpha_Y, *b_PriVtxB0Less_B0CosAlpha_Z ; 
   TBranch        *b_PriVtxB0Less_B0CosAlpha_EX, *b_PriVtxB0Less_B0CosAlpha_EY, *b_PriVtxB0Less_B0CosAlpha_EZ ;
   TBranch	  *b_PriVtxB0Less_B0CosAlpha_Chi2, *b_PriVtxB0Less_B0CosAlpha_CL;
   TBranch        *b_PriVtxB0Less_B0CosAlpha_tracks ;
   // Primary Vertex with largest B0_cos(alpha) 3D
   TBranch        *b_PriVtx_B0CosAlpha3D_X, *b_PriVtx_B0CosAlpha3D_Y, *b_PriVtx_B0CosAlpha3D_Z ; 
   TBranch        *b_PriVtx_B0CosAlpha3D_EX, *b_PriVtx_B0CosAlpha3D_EY, *b_PriVtx_B0CosAlpha3D_EZ ;
   TBranch	  *b_PriVtx_B0CosAlpha3D_Chi2, *b_PriVtx_B0CosAlpha3D_CL;
   TBranch        *b_PriVtx_B0CosAlpha3D_tracks ;
   TBranch        *b_PriVtxB0Less_B0CosAlpha3D_X, *b_PriVtxB0Less_B0CosAlpha3D_Y, *b_PriVtxB0Less_B0CosAlpha3D_Z ; 
   TBranch        *b_PriVtxB0Less_B0CosAlpha3D_EX, *b_PriVtxB0Less_B0CosAlpha3D_EY, *b_PriVtxB0Less_B0CosAlpha3D_EZ ;
   TBranch	  *b_PriVtxB0Less_B0CosAlpha3D_Chi2, *b_PriVtxB0Less_B0CosAlpha3D_CL;
   TBranch        *b_PriVtxB0Less_B0CosAlpha3D_tracks ;
   // Primary Vertex with "B0 correction" // check name
   TBranch        *b_PriVtxB0Corr_X, *b_PriVtxB0Corr_Y, *b_PriVtxB0Corr_Z;
   TBranch        *b_PriVtxB0Corr_EX, *b_PriVtxB0Corr_EY, *b_PriVtxB0Corr_EZ;
   TBranch	  *b_PriVtxB0Corr_Chi2, *b_PriVtxB0Corr_CL;
   TBranch        *b_PriVtxB0Corr_tracks;
   // Lifetimes
   TBranch        *b_B0CosAlphaBS, *b_B0CosAlpha3DBS, *b_B0CTauBS, *b_B0CTauBSE, *b_B0LxyBS ;
   TBranch        *b_B0CosAlphaPV, *b_B0CosAlpha3DPV, *b_B0CTauPV, *b_B0CTauPVE, *b_B0LxyPV ;
   TBranch        *b_B0CosAlphaPVCosAlpha, *b_B0CTauPVCosAlpha, *b_B0CTauPVCosAlphaE, *b_B0LxyPVCosAlpha ;
   TBranch        *b_B0CosAlphaPVCosAlpha3D, *b_B0CTauPVCosAlpha3D, *b_B0CTauPVCosAlpha3DE, *b_B0LxyPVCosAlpha3D ;
   TBranch        *b_B0CosAlphaB0LessPV, *b_B0CosAlpha3DB0LessPV, *b_B0CTauB0LessPV, *b_B0CTauB0LessPVE, *b_B0LxyB0LessPV ;
   TBranch        *b_B0CosAlphaB0LessPVCosAlpha, *b_B0CTauB0LessPVCosAlpha, *b_B0CTauB0LessPVCosAlphaE, *b_B0LxyB0LessPVCosAlpha ;
   TBranch        *b_B0CosAlphaB0LessPVCosAlpha3D, *b_B0CTauB0LessPVCosAlpha3D, *b_B0CTauB0LessPVCosAlpha3DE, *b_B0LxyB0LessPVCosAlpha3D ;
   TBranch        *b_B0CosAlphaPVX, *b_B0CTauPVX, *b_B0CTauPVXE, *b_B0LxyPVX ;
   TBranch	  *b_B0CTauPVX_3D, *b_B0CTauPVX_3D_err;
   TBranch        *b_B0MuMuIdx, *b_B0PionIdx, *b_B0KaonIdx;   //!
   TBranch        *b_PiPiMass_err;   //!
   TBranch        *b_Muon1Px_MuMuPiK, *b_Muon1Py_MuMuPiK, *b_Muon1Pz_MuMuPiK, *b_Muon1E_MuMuPiK;   //!
   TBranch        *b_Muon2Px_MuMuPiK, *b_Muon2Py_MuMuPiK, *b_Muon2Pz_MuMuPiK, *b_Muon2E_MuMuPiK;   //!
   TBranch        *b_PionPx_MuMuPiK, *b_PionPy_MuMuPiK, *b_PionPz_MuMuPiK, *b_PionE_MuMuPiK;   //!
   TBranch        *b_pion_nsigdedx, *b_pion_dedx, *b_pion_dedxMass, *b_pion_theo, *b_pion_sigma;   //!
   TBranch        *b_pion_dedx_byHits, *b_pion_dedxErr_byHits, *b_pion_saturMeas_byHits, *b_pion_Meas_byHits;   //!
   TBranch        *b_KaonPx_MuMuPiK, *b_KaonPy_MuMuPiK, *b_KaonPz_MuMuPiK, *b_KaonE_MuMuPiK;   //!
   TBranch        *b_kaon_nsigdedx, *b_kaon_dedx, *b_kaon_dedxMass, *b_kaon_theo, *b_kaon_sigma;   //!
   TBranch        *b_kaon_dedx_byHits, *b_kaon_dedxErr_byHits, *b_kaon_saturMeas_byHits, *b_kaon_Meas_byHits;   //!

   psiPrimePiK_MC(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~psiPrimePiK_MC() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   int getXVarType(const double jpsivar, vector<double> vectbin);

   //virtual bool     isMuonInAccept(float eta, float pt, float p);
   virtual bool     isMuonInAccept(float eta, float pt);
   virtual float    GetWeightPtShape(float pt, TH1F* histo );
   virtual float    ComputedEdxMass(int itrack, float TrackPt, float K, float C);
   virtual double   getMumuHelicityAngle(TLorentzVector B0, TLorentzVector Kstar, TLorentzVector mumu, TLorentzVector muM);
   virtual double   getBelleMumuHelicityAngle(TLorentzVector B0, TLorentzVector Kstar, TLorentzVector mumu, TLorentzVector muM);
   virtual double   getAleMumuHelicityAngle(TLorentzVector mumu, TLorentzVector muM);
   virtual double   getPlanesAngle(TLorentzVector B0, TLorentzVector K, TLorentzVector Pi, TLorentzVector muM, TLorentzVector muP);
   virtual double   getBellePlanesAngle(TLorentzVector B0, TLorentzVector K, TLorentzVector Pi, TLorentzVector muM, TLorentzVector muP);
   ///////////////////////////////////////// adde here my stuff

   ///////////// Input file //////////////////
   //TFile            *f1;

   ///////////// Output files //////////////////
   TProofOutputFile *OutFile, *SaraFile, *AlexisFile, *tightCutsFile ;
   TFile            *fOut, *SaraOut, *AlexisOut, *tightCutsOut ;
   Bool_t newSoftMuonID, only1psiPrime;
   Bool_t MC, officialMC, oldNtuple, priVtxsInfo, Kai_cuts, HLT_Displaced, fixedMassWind, sideband ;

   //TProofOutputFile *myOutFile, *myOutFileA, *myOutFileB;
   //ofstream *myoutfile, *myoutfileA, *myoutfileB; // works!

   /*
   ofstream myoutfile("DatasetAandB_KaonTrackRefit.txt"); // in.C
   ofstream myoutfileA("DatasetA_KaonTrackRefit.txt"); // in.C
   ofstream myoutfileB("DatasetB_KaonTrackRefit.txt"); // in.C
   */

   ///////////// Histograms //////////////////
   // 2011
   vector< vector<TH1F*> > myInclusiveMuMuMass_2011_HLT ;
   //
   // 2012
   TH1F* myInclusiveMuMuMass, *myInclusiveMuMuMass_noTrigg;
   //
   // HLT0_v*
   TH1F* myInclusiveMuMuMass_HLT_0_v3, *myInclusiveMuMuMass_HLT_0_v4, *myInclusiveMuMuMass_HLT_0_v5, *myInclusiveMuMuMass_HLT_0_v6 ;
   TH1F* myInclusiveMuMuMass_HLT_0_v3A, *myInclusiveMuMuMass_HLT_0_v3B ;
   TH1F* myInclusiveMuMuMass_HLT_0_v6C, *myInclusiveMuMuMass_HLT_0_v6D ;
   TH1F* myInclusiveMuMuMass_HLT_0_v3Bv4, *myInclusiveMuMuMass_HLT_0_v5v6C ;
   // HLT5_v*
   TH1F* myInclusiveMuMuMass_HLT_5_v3, *myInclusiveMuMuMass_HLT_5_v4, *myInclusiveMuMuMass_HLT_5_v5, *myInclusiveMuMuMass_HLT_5_v6 ;
   TH1F* myInclusiveMuMuMass_HLT_5_v3A, *myInclusiveMuMuMass_HLT_5_v3B ;
   TH1F* myInclusiveMuMuMass_HLT_5_v6C, *myInclusiveMuMuMass_HLT_5_v6D ;
   TH1F* myInclusiveMuMuMass_HLT_5_v3Bv4, *myInclusiveMuMuMass_HLT_5_v5v6C ;
   // HLT7_v*
   TH1F* myInclusiveMuMuMass_HLT_7_v1, *myInclusiveMuMuMass_HLT_7_v2, *myInclusiveMuMuMass_HLT_7_v3 ;
   TH1F* myInclusiveMuMuMass_HLT_7_v3C, *myInclusiveMuMuMass_HLT_7_v2v3C, *myInclusiveMuMuMass_HLT_7_v3D ;
   // datasets
   TH1F* myInclusiveMuMuMass_12A;
   TH1F* myInclusiveMuMuMass_12A_HLT_0_v3;
   TH1F* myInclusiveMuMuMass_12A_HLT_5_v3;
   TH1F* myInclusiveMuMuMass_12A_HLT_9_v9;
   //
   TH1F* myInclusiveMuMuMass_12B;
   TH1F* myInclusiveMuMuMass_12B_HLT_0_v3v4;
   TH1F* myInclusiveMuMuMass_12B_HLT_5_v3v4;
   TH1F* myInclusiveMuMuMass_12B_HLT_7_v1_9_v9;
   //
   TH1F* myInclusiveMuMuMass_12C;
   TH1F* myInclusiveMuMuMass_12C_HLT_0_v5;
   TH1F* myInclusiveMuMuMass_12C_HLT_5_v5;
   TH1F* myInclusiveMuMuMass_12C_HLT_7_v2;
   TH1F* myInclusiveMuMuMass_12C_HLT_5_v5_or_7_v2;
   //
   TH1F* myInclusiveMuMuMass_12D;
   TH1F* myInclusiveMuMuMass_12D_HLT_0_v6;
   TH1F* myInclusiveMuMuMass_12D_HLT_5_v6;
   TH1F* myInclusiveMuMuMass_12D_HLT_7_v3;
   TH1F* myInclusiveMuMuMass_out;
   //
   TH1F* myMuMuMass, *myMuMuMass_noMuID, *myMuMuMass_MuID, *myMuMuMass_Mu1ID, *myMuMuMass_Mu2ID ;
   TH1F *myMuMuMass_HLT, *myMuMuMass_HLT7, *myMuMuMass_HLT5_7comp, *myMuMuMass_HLT3p5 ; //, *myMuMuMass_HLT_yMu, *myMuMuMass_HLT_yMu_ypsi,
   //TH1F* myMuMuMass_test ;
   TH1F* myMuMuMass_HLT_final, *myMuMuMass_HLT_DR1, *myMuMuMass_HLT_VtxCLDiff_cut ;
   TH1F* myMuMuMass_HLT_max2psiPrime, *myMuMuMass_HLT_2psiPrime, *myMuMuMass_HLT_1psiPrime ;
   TH1F* myMuMuKMass ;
   TH1F* myMuMuPiMass , *myMuMuPiMass_nonpr, *myMuMuPiMass_pr, *myMuMuPiMass_quasipr, *myMuMuPiMass_Kai, *myMuMuPiMass_Kai_HLT7 ;
   TH1F* myMuMuPiMass_DeltaR0p1_nonpr5, *myMuMuPiMass_DeltaR0p2_nonpr5, *myMuMuPiMass_DeltaR0p1_nonpr8, *myMuMuPiMass_DeltaR0p2_nonpr8 ;
   TH2F* psi2S_DeltaR_pi_vs_mumupiMass_nonpr ;
   TH1F* myMuMuPiMass_HLT5to3p5, *myMuMuPiMass_HLT3p5_allCuts, *myMuMuKMass_HLT5to3p5, *myMuMuKMass_HLT3p5_allCuts ;
   TH1F* myMuMuKMass_Kai, *myMuMuMass_KaiBF, *myMuMuKMass_KaiBF, *myMuMuMass_noMuID_KaiBF, *myMuMuMass_MuID_KaiBF ;
   TH2F *dEdx_vs_p_Bplus ;
   TH1F* myMuMuPiMass_pr_DzCut, *myMuMuPiMass_pr_Dz01, *myMuMuPiMass_nonpr_Dz1inf ;
   TH1F* myMuMuPiMass_barrel, *myMuMuPiMass_barrel_pr, *myMuMuPiMass_nonBarrel_pr, *myMuMuPiMass_barrel_nonpr ;
   // Sara plots
   vector<TH1F*> MuMuKMass_ZCTauBSOverE, MuMuKMass_mu1Pt, MuMuKMass_mu2Pt, MuMuKMass_trk1Pt, MuMuKMass_trk2Pt, MuMuKMass_ZCL ;
   // Nairit plots
   TH1F* MuMuKMass_Nairit, *MuMuKMass_sb, *trk_pT_test_h ;
   // K*
   TH1F* myExclusiveMuMuKPiMass ;
   TH1F* trk_dxy, *trk_dz, *hmyPsiPMassSel, *hPsiPrimefromBMassSel, *hmyPsiPKPiMassSelAlt, *hmyPsiPKPiMassSelAltZoom, *hmyPsiPKPiMassSelZoom ;
   TH2F* psi2SPiSqMass_vs_KPiSqMass_AC ;
   // tight cuts
   vector <Float_t> piPt_B0, KPt_B0, vtxCL_B0, cTau_B0 ;
   vector <TH1F*> myMuMuKPiMass_piPt, myMuMuKPiMass_KPt, myMuMuKPiMass_vtxCL, myMuMuKPiMass_cTau ;
   //vector < tuple<Float_t, Float_t, Float_t> > trksPt_B0VtxCL_B0CTau ;
   struct tightCuts {
     Float_t trksPt; Float_t B0VtxCL; Float_t B0CTau;
   };
   vector < tightCuts > trksPt_B0VtxCL_B0CTau ;
   vector <TH1F*> myMuMuKPiMass_trksAndB0vtxCLAndB0Ctau ;
   //
   // PV matching
   TH1F* priVtx_deltaZ ;
   TH2F* B0_pointingAngle_PVvsB0LessPV[2], *z_PVvsB0LessPV[2] ;
   TH1F* B0_PVCTau_h, *B0_CosAlphaCTau_h, *B0_PVXCTau_h, *B0_PVX3DCTau_h ;
   TH1I* priVtx_matched ;
   //
   TH2F* kaonP_vs_pionP, *kaonPt_vs_pionPt ;
   TH1F* hmyPsiPKPiMass_B0pTtest, *hPionFromKstar_p, *hKaonFromKstar_p ;
   TH1F* hmyPsiPKPiMass_deltaRtest, *hmyPsiPKPiMass_B0cTautest ;
   TH1F* hmyPsiPKPiMassBaseSelAlt ;
   TH1F* hmyPsiPKPiMass, *hmyPsiPKPiMass_KPiExchanged, *hmyPsiPKPiMass_bothComb ;
   TH1F* hB0VtxCL, *hB0_Pt ;
   TH1F* hB0CTauPVSignif, *hB0CTauPVXSignif ;
   TH1F* hPsiPrimefromBMass, *hPsiPrime_fromB_Pt, *hmyPsiP_Pt ;
   TH1F* hPionfromBpt, *hKaonfromBpt, *hPionFromKstar_pT, *hKaonFromKstar_pT, *hPionFromKstar_ex_pT, *hKaonFromKstar_ex_pT, *hPionFromKstar_sb_pT, *hKaonFromKstar_sb_pT ;
   TH1F* myExclMuMuKPiMassSel, *hmyPsiPKPiMass_bothCombSel ;
   TH2F* psi2SPiSqMass_vs_KPiSqMass_ABaseC ;
   TH1F* hmyPsiPKPiMass_fromKstar, *hmyPsiPKPiMass_fromKstar_ex_noRange, *hmyPsiPKPiMass_fromKstar_diff, *hmyPsiPKPiMass_fromKstar_ex, *hmyPsiPKPiMass_fromKstar_sb ;
   TH1F* hmyKPiMass, *hmyKPiMass_ex, *hmyKPiMass_KpTcut, *hmyKPiMass_ex_KpTcut, *hmyKPiMass_sb;
   TH1F* hmyKPiMass_KpTcut_tightTr, *hmyKPiMass_ex_KpTcut_tightTr, *hmyKPiMass_KpTcut_HP, *hmyKPiMass_ex_KpTcut_HP ;
   TH2F* hpsiPK_vs_KPi, *hpsiPK_vs_KPi_ex, *hpsiPK_vs_KPi_sb ;
   TH2F* dEdx_vs_pT_Kstar, *dEdx_vs_p_Kstar ;
   TH2F* dEdx_vs_pT_piFromKstar, *dEdx_vs_pT_KFromKstar, *dEdx_vs_p_piFromKstar, *dEdx_vs_p_KFromKstar, *dEdx_vs_p_piFromKstar_KpCut, *dEdx_vs_p_KFromKstar_KpCut, *dEdxK_vs_dEdxPi ;
   TH1I* piFromKstar_stripHits_h, *KFromKstar_stripHits_h ;
   TH2F* piFromKstar_stripHits_vs_y_h, *KFromKstar_stripHits_vs_y_h ;
   TH2F* dEdx_vs_pT_KhqFromKstar, *dEdx_vs_p_KhqFromKstar ;
   TH1F* dEdxPi_pCut, *dEdxK_pCut, *dEdxPi_pCut_stripCut, *dEdxK_pCut_stripCut ;
   TH2F* Armenteros_Kstar_orig, *Armenteros_Kstar, *Armenteros_Kstar_orig_ex, *Armenteros_Kstar_orig_sb, *alpha_vs_alphaEx ;
   TH1F* pT_diff, *pT_diff_ex, *pT_diff_sb ;
   //   
   // my scan
   vector< vector< vector<TH1F*> > > myMuMuPiMass_piPt, myMuMuPiMass_ZPt, myMuMuPiMass_piAndZPt ;
   //
   //TH1F* psiPrime_rapidBin[24] ;
   vector<TH1F*> psiPrime_rapidBin, psiPrime_2011_rapidBin ;
   //
   // selection varables
   // Z- variables
   TH1I* nZ_h, *nB0_h, *nB0AC_noMassWin_h, *nB0AC_h, *nB0AC_signalWin_h, *nB0AC_purityWin_h, *nB0ACInMC_h, *nMCB0_h, *nMCB0_BT_h, *nTwins_h[2] ;
   TH1F* hB0Mass_1B0, *hB0Mass_1B0matched[2][2], *hB0Mass_noTwins_noSignalWinNotTwins ;
   TH2F* psi2SPiSqMass_vs_KPiSqMass_1B0[2] ;
   TH1F* psi2S_helicityAngle, *cos_psi2S_helicityAngle, *planesAngle ;
   TH1F* hB0CTau_1B0, *hB0CTau_1B0matched ;
   TH1F* B0_gen_p_h, *B0_gen_pT_h, *psi2S_gen_p_h, *psi2S_gen_pT_h, *K_gen_p_h, *K_gen_pT_h, *pi_gen_p_h, *pi_gen_pT_h ;
   TH1F* hMCDeltaRPi_2B0, *hMCDeltaRK_2B0 ;
   TH1F* B0_gen_mass_h ;
   TH2F* psi2SPiSqMass_vs_KPiSqMass_gen, *psi2SPiSqMass_vs_KPiSqMass_BT_gen, *KPiSqMass_vs_psi2SKSqMass_BT_gen, *piCh_vs_KCh_gen ;
   TH1F* psi2S_helicityAngle_BT_gen, *cos_psi2S_helicityAngle_BT_gen, *planesAngle_BT_gen ;
   TH2F* B0CosAlpha_2Dvs3D[4] ;
   TH1F* priVtx_z[7] ;
   TH1F* priVtx_delta[2][3], *priVtx_deltaX_test, *priVtx_deltaY_test, *priVtx_B0_delta[2][3], *priVtx_B0CosAlpha_delta[2][2][3], *priVtxB0Less_delta[2][3], *priVtxB0Less_B0CosAlpha_delta[2][2][3] ;
   TH1F* dz_pri_priB0CosAlpha[2][2], *dz_pri_priB0Less[2], *dz_pri_priB0LessB0CosAlpha[2][2], *dz_priB0CosAlpha_priB0LessB0CosAlpha[2][2], *dz_priB0Less_priB0LessB0CosAlpha[2][2] ;
   TH1F* hMCDeltaR_1B0, *hMCDeltaPt_1B0 ;
   TH1I* hMCTruthMatching_1B0[2][3], *hMCTruthMatching_2B0notTwin[3], *hMCTruthMatching_2B0twin[3], *hMCTruthMatching_2B0signalWin_notTwin ; 
   TH1F* hB0Mass_2B0notTwin, *hB0Mass_2B0notTwin_matched[3], *hB0Mass_2B0twin, *hB0Mass_2B0twin_matched[3] ;
   TH2F* hDeltaB0Mass_2B0twin[2], *hAlpha_2B0twin[2] ;
   TH1F* hB0Mass_2B0twin_hourglass[2][2] ;
   TH2F* hDeltaB0Mass_2B0notTwin[2] ;
   vector<TH1F*> hMCDeltaR_2B0, hMCDeltaPt_2B0 ;
   TH2F* hMCDeltaPt_1vs2_2B0, *hMCDeltaR_1vs2_2B0 ;
   TH2F* hMCDeltaR_1vs2_2B0twin, *hMCDeltaPt_1vs2_2B0twin, *hB0VtxCL_1vs2_2B0twin, *hB0Mass_1vs2_2B0twin, *hB0CTau_1vs2_2B0twin, *hAlpha_1vs2_2B0twin ;
   TH1F* hMCDeltaR_matchedTwin[2], *hMCDeltaPt_matchedTwin[2], *hB0VtxCL_matchedTwin[2], *hB0Mass_matchedTwin[2], *hB0CTau_matchedTwin[2], *hAlpha_matchedTwin[2] ;
   //TH1I*  hMCnGenX ;
   TH1F* SpectrumFromBPlus, *SpectrumFromBPlusCtaued ;  
   // pion matching
   TH1F* hMCDeltaRMuMuMatch, *hMCDeltaRPiMatch, *hMCDeltaRK_fromTwins ;
   TH2F* hMCPtRecoVsPtTrueMuMuMatch, *hMCPtRecoVsPtTruePiMatch, *hMCPtRecoVsPtTrueK_fromTwins ;
   TH1F* hMCpTTrueMuMuMatch, *hMCPtTruePiMatch, *hMCpTTrueK_fromTwins ;
   TH1F* hMCDeltaPtK_fromTwins, *hMCDeltaPtK_fromTwinsPiMatched ;
   TH2F* kaonDr_vs_pionDr_matched, *kaonP_vs_pionP_matched, *kaonPt_vs_pionPt_matched ;
   TH1F* hMCStripsHitsTruePiMatch, *hMCPixelHitsTruePiMatch, *hMCTotalHitsTruePiMatch, *hMCChiTruePiMatch, *hMCStripsHitsTrackNotMached, *hMCPixelHitsTrackNotMached, *hMCTotalHitsTrackNotMached, *hMCChiTrackNotMached ;


   TH1F* dEdxMass_h[2] ;
   TH2F* dEdxMass_vs_pT, *dEdx_vs_p, *dEdx_vs_pT[3] ;
   TH1F* B0_Mass_h, *myB0_Mass_h ;
   TH1F* Z_VtxCL_h ;
   TH1F* Z_CTauBS_h, *Z_CTauBSE_h, *Z_CTauBS_over_Z_CTauBSE_h ;
   TH1F* Z_pT_h ;
   TH1F* Z_pointingAngle[3] ;
   vector<TH1F*> Z_cosPA ;
   //vector< pair<TH1F*, TH1F*> > Z_cosPA ;
   TH2F* Z_vs_pi_pT_h ;
   // psi(2S) variables
   TH1F* psi2S_VtxCL_h, *psi2S_VtxCL_HLT5_h, *psi2S_VtxCL_1_h, *psi2S_VtxCL_2_h, *psi2S_VtxCL_diff_h ;
   TH1I* psi2SMuonTrigMatch_h ;
   TH1F* psi2S_pT_h ;
   vector<TH1F*> psi2S_2011_pT_h ;
   vector< vector<TH1F*> > psi2S_2011v_pT_h ;
   TH1F* psi2S_y_h, *psi2S_HLT5_yMu_h ;
   //TH1F* psi2S_DeltaR_pi_h, *psi2S_DeltaR_pi_pr_h, *psi2S_DeltaR_pi_quasipr_h, *psi2S_DeltaR_pi_nonpr_h ;
   vector<TH1F*> psi2S_DeltaR_pi_h ;
   vector<TH2F*> Z_DeltaR_vs_cosPA ;
   //TH1F* psi2S_pointingAngle, *psi2S_pointingAngle_sign, *psi2S_pointingAngle_sb ; // empty variable
   TH2F* pT_vs_y_psi2S_h ;
   TH1I* nPsi2S_h ;
   // mu variables
   TH1F* mu_px_psi2S, *mu_px, *fitmu_px ;
   TH1I* mu1_pixelHits_h, *mu2_pixelHits_h ;
   TH1I* mu1_stripHits_h , *mu2_stripHits_h ;
   TH1F* mu1_charge_h, *mu2_charge_h ; 
   TH1F* mu1_Chi2N_h , *mu2_Chi2N_h, *mu1Tr_Chi2N_h, *mu2Tr_Chi2N_h, *mu1GlTr_Chi2N_h, *mu2GlTr_Chi2N_h ;
   TH1F* mu1_Dz_h, *mu2_Dz_h ;
   TH1F* mu1_Dxy_h, *mu2_Dxy_h ;
   TH1F* mu1_y_h, *mu2_y_h ;
   TH1F* mu1_pT_h, *mu2_pT_h, *mu1_pT_test_h, *mu2_pT_test_h ;
   TH1F* mu1_eta_test_h, *mu2_eta_test_h ;
   TH1F* mu1_DeltaR_mu2_h ; //, *mu2_DeltaR_mu1_h ;
   //TH2F* mu1_vs_mu2_pT_h ;
   TH2F* pT_vs_y_mu1_h, *pT_vs_y_mu2_h ;
   TH2F* pT_vs_eta_mu_h, *pT10_vs_eta_mu_h ;
   // pi+ variables
   TH1I* pi_pixelHits_h ;
   TH1I* pi_stripHits_h ;
   TH1F* pi_Chi2N_h , *pi_Chi2N_test_h ;
   TH1F* pi_pT_h , *pi_pT_test_h ;
   TH1F* pi_y_h , *pi_y_test_h ;
   TH2F* pT_vs_pixelHits_pi_h, *pT_vs_stripHits_pi_h ; 
   TH2F* pT_vs_y_pi_h ;
   TH1F* pi_Dz_h, *pi_Dxy_h, *pi_Dz_pr_h, *pi_Dxy_pr_h, *pi_Dz_nonpr_h, *pi_Dxy_nonpr_h ;
   
   /////////////////////////////////////////
   
   //////// CUT THRESHOLDS & CONSTANTS /////////////
   Int_t B0_Id, psi2S_Id, Kaon_Id, pion_Id, muon_Id ;
   Double_t muon_mass, phi_mass, jpsi_mass, psi2S_mass, pionCh_mass, kaonCh_mass, B0_mass;
   Float_t B0_massFit, B0_sigmaFit ;
   //Double_t Bplus_mass;

   Float_t rapid_min, rapid_max, rapid_binSize ;  Int_t rapid_bins;
   Float_t psiPrimeMass_min, psiPrimeMass_max ;  Int_t psiPrimeMass_bins ; 
   
   Float_t cosAlpha_cut, cosAlpha_cut2 ;
   Float_t Z_cosPA_min , Z_cosPA_max ;
   Float_t scale_factor ;

   vector<Float_t> psi2S_2011_pT ;
   vector< vector<TString> > HLT_v ;

   // Sara plots
   vector<Float_t> Z_CTauBSOverE, mu_pT, trk_pT, Z_CL ;
   
   // scan
   vector<TString> mode ;
   vector<Float_t> pi_pT, Z_pT ;
   //vector< *Float_t > piAndZ_pT ; // never use vector <array>
   vector< pair<Float_t, Float_t> > piAndZ_pT ;
   vector<Float_t> deltaR_piPsi2S ;
   
   /////////////////////////////////////

   ClassDef(psiPrimePiK_MC,0);
};

#endif

#ifdef psiPrimePiK_MC_cxx
void psiPrimePiK_MC::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

  // all these varables cannot be read from SlaveBegin method because of this structure: https://root.cern.ch/drupal/content/developing-tselector
  ////////////////////////////////
  //MC = kFALSE ; oldNtuple = kTRUE ; // for June 2014 ntuples
  //MC = kTRUE ; oldNtuple = kFALSE ; // for late 2014 MC
  TString option = GetOption();
  if (option == "MC" || option == "officialMC") {
    MC = kTRUE ;
    if (option == "officialMC") 
      officialMC = kTRUE ;
  } else {
    MC = kFALSE ;  
    officialMC = kFALSE ;  
  }

  if ( MC || officialMC )
    oldNtuple = kFALSE ;
  else 
    oldNtuple = kTRUE ;

  priVtxsInfo = kFALSE ; priVtxsInfo = kTRUE ;

  B0_Id = 511 ; psi2S_Id = 100443 ; Kaon_Id = 321; pion_Id = 211; muon_Id = 13;
  
  muon_mass = 0.1056583715;
  // phi_mass = 1.019455;
  // jpsi_mass = 3.0967;
  psi2S_mass = 3.686109 ;
  pionCh_mass = 0.13957018;
  kaonCh_mass = 0.493677 ; 
  // Bplus_mass = 5.27925;
  B0_mass = 5.27958;
  B0_massFit = 5.27967; B0_sigmaFit = 0.01167;
  /*
  rapid_min = -2.4;
  rapid_max = 2.4;
  rapid_bins = 16 ;
  psiPrimeMass_min = 3.4 ;
  psiPrimeMass_max = 4.0 ;
  psiPrimeMass_bins = 120 ;
  */
  ///////////////////////////////

   // Set object pointer
   TrigRes = 0;
   TrigNames = 0;
   MatchTriggerNames = 0;
   L1TrigRes = 0;
   MCPdgIdAll = 0;
   PriVtxGen_X = 0; PriVtxGen_Y = 0; PriVtxGen_Z = 0;
   PriVtxGen_EX = 0; PriVtxGen_EY = 0; PriVtxGen_EZ = 0;
   PriVtxGen_Chi2 = 0; PriVtxGen_CL = 0; PriVtxGen_Ndof = 0;
   PriVtxGen_tracks = 0;
   genEvtVtx_X = 0; genEvtVtx_Y = 0; genEvtVtx_Z = 0;
   genEvtVtx_particlesOut = 0;
   MCpsi2SPx = 0; MCpsi2SPy = 0; MCpsi2SPz = 0;
   MCmupPx = 0; MCmupPy = 0; MCmupPz = 0;
   MCmumPx = 0; MCmumPy = 0; MCmumPz = 0;
   MCpionPx = 0; MCpionPy = 0; MCpionPz = 0;
   MCkaonPx = 0; MCkaonPy = 0; MCkaonPz = 0;
   MCPx = 0; MCPy = 0; MCPz = 0;
   MCpionCh = 0; MCkaonCh = 0;

   // Generic Muons
   muPx = 0; muPy = 0; muPz = 0;
   muD0 = 0; muDz = 0;
   muChi2 = 0; muNDF = 0;
   muPhits = 0; muShits = 0;
   muLayersTr = 0; muLayersPix = 0;
   muD0E = 0; muDzVtxErr = 0;
   muKey = 0;
   muGlMuHits = 0;
   muGlChi2 = 0; muGlNDF = 0;
   muGlMatchedStation = 0;
   muGlDzVtx = 0; muGlDxyVtx = 0;
   muType = 0;
   muQual = 0;
   muTrack = 0;
   muCharge = 0;
   mufHits = 0;
   muFirstBarrel = 0; muFirstEndCap = 0;
   muDzVtx = 0; muDxyVtx = 0;
   trNotRef = 0; trRef = 0;
   trackPx = 0; trackPy = 0; trackPz = 0;
   trackEnergy = 0;
   trackNDF = 0;
   trackPhits = 0; trackShits = 0;
   trackChi2 = 0;
   trackD0 = 0; trackD0Err = 0;
   trackCharge = 0;
   TrackHighPurity = 0;
   TrackTight = 0;
   trackfHits = 0;
   trackFirstBarrel = 0; trackFirstEndCap = 0;
   trackDzVtx = 0; trackDxyVtx = 0;
   tr_nsigdedx = 0; tr_dedx = 0; tr_dedxMass = 0; tr_theo = 0; tr_sigma = 0;
   tr_dedx_byHits = 0; tr_dedxErr_byHits = 0; tr_saturMeas_byHits = 0; tr_Meas_byHits = 0;
   MuMuMass = 0;
   MuMuPx = 0; MuMuPy = 0; MuMuPz = 0;
   MuMuVtx_CL = 0; MuMuVtx_Chi2 = 0;
   MuMuDecayVtx_X = 0; MuMuDecayVtx_Y = 0; MuMuDecayVtx_Z = 0;
   MuMuDecayVtx_XE = 0; MuMuDecayVtx_YE = 0; MuMuDecayVtx_ZE = 0;
   mu1Idx = 0; mu2Idx = 0;
   mu1Px_MuMu = 0; mu1Py_MuMu = 0; mu1Pz_MuMu = 0;
   mu1Chi_MuMu2 = 0; mu1NDF_MuMu = 0;
   mu2Px_MuMu = 0; mu2Py_MuMu = 0; mu2Pz_MuMu = 0;
   mu2Chi2_MuMu = 0; mu2NDF_MuMu = 0;
   MuMuType = 0;
   MuMuMuonTrigMatch = 0;
   PriVtxMuMuCorr_X = 0; PriVtxMuMuCorr_Y = 0; PriVtxMuMuCorr_Z = 0;
   PriVtxMuMuCorr_EX = 0; PriVtxMuMuCorr_EY = 0; PriVtxMuMuCorr_EZ = 0;
   PriVtxMuMuCorr_Chi2 = 0; PriVtxMuMuCorr_CL = 0;
   PriVtxMuMuCorr_tracks = 0;
   B0Mass = 0;
   B0Px = 0; B0Py = 0; B0Pz = 0;
   B0PxE = 0; B0PyE = 0; B0PzE = 0;
   B0Vtx_CL = 0; B0Vtx_Chi2 = 0;
   B0DecayVtx_X = 0; B0DecayVtx_Y = 0; B0DecayVtx_Z = 0;
   B0DecayVtx_XE = 0; B0DecayVtx_YE = 0; B0DecayVtx_ZE = 0;
   // Primary Vertex with largest B0_cos(alpha)
   PriVtx_B0CosAlpha_X = 0; PriVtx_B0CosAlpha_Y = 0; PriVtx_B0CosAlpha_Z = 0;
   PriVtx_B0CosAlpha_EX = 0; PriVtx_B0CosAlpha_EY = 0; PriVtx_B0CosAlpha_EZ = 0;
   PriVtx_B0CosAlpha_Chi2 = 0; PriVtx_B0CosAlpha_CL = 0;
   PriVtx_B0CosAlpha_tracks = 0;
   B0LessPV_tracksPtSq = 0; B0LessPV_4tracksPtSq = 0;
   PriVtxB0Less_X = 0; PriVtxB0Less_Y = 0; PriVtxB0Less_Z = 0;
   PriVtxB0Less_EX = 0; PriVtxB0Less_EY = 0; PriVtxB0Less_EZ = 0;
   PriVtxB0Less_Chi2 = 0; PriVtxB0Less_CL = 0;
   PriVtxB0Less_tracks = 0;
   PriVtxB0Less_B0CosAlpha_X = 0; PriVtxB0Less_B0CosAlpha_Y = 0; PriVtxB0Less_B0CosAlpha_Z = 0;
   PriVtxB0Less_B0CosAlpha_EX = 0; PriVtxB0Less_B0CosAlpha_EY = 0; PriVtxB0Less_B0CosAlpha_EZ = 0;
   PriVtxB0Less_B0CosAlpha_Chi2 = 0; PriVtxB0Less_B0CosAlpha_CL = 0;
   PriVtxB0Less_B0CosAlpha_tracks = 0;
   // Primary Vertex with largest B0_cos(alpha) 3D
   PriVtx_B0CosAlpha3D_X = 0; PriVtx_B0CosAlpha3D_Y = 0; PriVtx_B0CosAlpha3D_Z = 0;
   PriVtx_B0CosAlpha3D_EX = 0; PriVtx_B0CosAlpha3D_EY = 0; PriVtx_B0CosAlpha3D_EZ = 0;
   PriVtx_B0CosAlpha3D_Chi2 = 0; PriVtx_B0CosAlpha3D_CL = 0;
   PriVtx_B0CosAlpha3D_tracks = 0;
   PriVtxB0Less_B0CosAlpha3D_X = 0; PriVtxB0Less_B0CosAlpha3D_Y = 0; PriVtxB0Less_B0CosAlpha3D_Z = 0;
   PriVtxB0Less_B0CosAlpha3D_EX = 0; PriVtxB0Less_B0CosAlpha3D_EY = 0; PriVtxB0Less_B0CosAlpha3D_EZ = 0;
   PriVtxB0Less_B0CosAlpha3D_Chi2 = 0; PriVtxB0Less_B0CosAlpha3D_CL = 0;
   PriVtxB0Less_B0CosAlpha3D_tracks = 0;
   // Primary Vertex with "B0 correction"
   PriVtxB0Corr_X = 0; PriVtxB0Corr_Y = 0; PriVtxB0Corr_Z = 0;
   PriVtxB0Corr_EX = 0; PriVtxB0Corr_EY = 0; PriVtxB0Corr_EZ = 0;
   PriVtxB0Corr_Chi2 = 0; PriVtxB0Corr_CL = 0;
   PriVtxB0Corr_tracks = 0;
   // Lifetime variables
   B0CosAlphaBS = 0; B0CosAlpha3DBS = 0; B0CTauBS = 0; B0CTauBSE = 0; B0LxyBS = 0; 
   B0CosAlphaPV = 0; B0CosAlpha3DPV = 0; B0CTauPV = 0; B0CTauPVE = 0; B0LxyPV = 0; 
   B0CosAlphaPVCosAlpha = 0; B0CTauPVCosAlpha = 0; B0CTauPVCosAlphaE = 0; B0LxyPVCosAlpha = 0; 
   B0CosAlphaPVCosAlpha3D = 0; B0CTauPVCosAlpha3D = 0; B0CTauPVCosAlpha3DE = 0; B0LxyPVCosAlpha3D = 0; 
   B0CosAlphaB0LessPV = 0; B0CosAlpha3DB0LessPV = 0; B0CTauB0LessPV = 0; B0CTauB0LessPVE = 0; B0LxyB0LessPV = 0; 
   B0CosAlphaB0LessPVCosAlpha = 0; B0CTauB0LessPVCosAlpha = 0; B0CTauB0LessPVCosAlphaE = 0; B0LxyB0LessPVCosAlpha = 0; 
   B0CosAlphaB0LessPVCosAlpha3D = 0; B0CTauB0LessPVCosAlpha3D = 0; B0CTauB0LessPVCosAlpha3DE = 0; B0LxyB0LessPVCosAlpha3D = 0; 
   B0CosAlphaPVX = 0; B0CTauPVX = 0; B0CTauPVXE = 0; B0LxyPVX = 0; 
   B0CTauPVX_3D = 0; B0CTauPVX_3D_err = 0;
   B0MuMuIdx = 0; B0PionIdx = 0; B0KaonIdx = 0;
   PiPiMass_err = 0;
   Muon1Px_MuMuPiK = 0; Muon1Py_MuMuPiK = 0; Muon1Pz_MuMuPiK = 0; Muon1E_MuMuPiK = 0;
   Muon2Px_MuMuPiK = 0; Muon2Py_MuMuPiK = 0; Muon2Pz_MuMuPiK = 0; Muon2E_MuMuPiK = 0;
   PionPx_MuMuPiK = 0; PionPy_MuMuPiK = 0; PionPz_MuMuPiK = 0; PionE_MuMuPiK = 0;
   pion_nsigdedx = 0; pion_dedx = 0; pion_dedxMass = 0; pion_theo = 0; pion_sigma = 0;
   pion_dedx_byHits = 0; pion_dedxErr_byHits = 0; pion_saturMeas_byHits = 0; pion_Meas_byHits = 0;
   KaonPx_MuMuPiK = 0; KaonPy_MuMuPiK = 0; KaonPz_MuMuPiK = 0; KaonE_MuMuPiK = 0;
   kaon_nsigdedx = 0; kaon_dedx = 0; kaon_dedxMass = 0; kaon_theo = 0; kaon_sigma = 0;
   kaon_dedx_byHits = 0; kaon_dedxErr_byHits = 0; kaon_saturMeas_byHits = 0; kaon_Meas_byHits = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("TrigRes", &TrigRes, &b_TrigRes);
   fChain->SetBranchAddress("TrigNames", &TrigNames, &b_TrigNames);
   fChain->SetBranchAddress("MatchTriggerNames", &MatchTriggerNames, &b_MatchTriggerNames);
   fChain->SetBranchAddress("L1TrigRes", &L1TrigRes, &b_L1TrigRes);
   fChain->SetBranchAddress("evtNum", &evtNum, &b_evtNum);
   fChain->SetBranchAddress("runNum", &runNum, &b_runNum);
   fChain->SetBranchAddress("lumiNum", &lumiNum, &b_lumiNum);
   fChain->SetBranchAddress("priVtx_X", &priVtx_X, &b_priVtx_X);
   fChain->SetBranchAddress("priVtx_Y", &priVtx_Y, &b_priVtx_Y);
   fChain->SetBranchAddress("priVtx_Z", &priVtx_Z, &b_priVtx_Z);
   fChain->SetBranchAddress("priVtx_XE", &priVtx_XE, &b_priVtx_XE);
   fChain->SetBranchAddress("priVtx_YE", &priVtx_YE, &b_priVtx_YE);
   fChain->SetBranchAddress("priVtx_ZE", &priVtx_ZE, &b_priVtx_ZE);
   fChain->SetBranchAddress("priVtx_NormChi2", &priVtx_NormChi2, &b_priVtx_NormChi2);
   fChain->SetBranchAddress("priVtx_Chi2", &priVtx_Chi2, &b_priVtx_Chi2);
   fChain->SetBranchAddress("priVtx_CL", &priVtx_CL, &b_priVtx_CL);
   fChain->SetBranchAddress("priVtx_tracks", &priVtx_tracks, &b_priVtx_tracks);
   if ( MC || officialMC ) { // MC Analysis
     fChain->SetBranchAddress("nMCAll", &nMCAll, &b_nMCAll);
     fChain->SetBranchAddress("MCPdgIdAll", &MCPdgIdAll, &b_MCPdgIdAll);
     fChain->SetBranchAddress("nMCB0", &nMCB0, &b_nMCB0);
     //
     fChain->SetBranchAddress("PriVtxGen_X", &PriVtxGen_X, &b_PriVtxGen_X);
     fChain->SetBranchAddress("PriVtxGen_Y", &PriVtxGen_Y, &b_PriVtxGen_Y);
     fChain->SetBranchAddress("PriVtxGen_Z", &PriVtxGen_Z, &b_PriVtxGen_Z);
     fChain->SetBranchAddress("PriVtxGen_EX", &PriVtxGen_EX, &b_PriVtxGen_EX);
     fChain->SetBranchAddress("PriVtxGen_EY", &PriVtxGen_EY, &b_PriVtxGen_EY);
     fChain->SetBranchAddress("PriVtxGen_EZ", &PriVtxGen_EZ, &b_PriVtxGen_EZ);
     fChain->SetBranchAddress("PriVtxGen_Chi2", &PriVtxGen_Chi2, &b_PriVtxGen_Chi2);
     fChain->SetBranchAddress("PriVtxGen_CL", &PriVtxGen_CL, &b_PriVtxGen_CL);
     fChain->SetBranchAddress("PriVtxGen_Ndof", &PriVtxGen_Ndof, &b_PriVtxGen_Ndof);
     fChain->SetBranchAddress("PriVtxGen_tracks", &PriVtxGen_tracks, &b_PriVtxGen_tracks);
     //
     fChain->SetBranchAddress("genEvtVtx_X", &genEvtVtx_X, &b_genEvtVtx_X);
     fChain->SetBranchAddress("genEvtVtx_Y", &genEvtVtx_Y, &b_genEvtVtx_Y);
     fChain->SetBranchAddress("genEvtVtx_Z", &genEvtVtx_Z, &b_genEvtVtx_Z);
     fChain->SetBranchAddress("genEvtVtx_particles", &genEvtVtx_particlesOut, &b_genEvtVtx_particlesOut);
     //
     fChain->SetBranchAddress("MCpsi2SPx", &MCpsi2SPx, &b_MCpsi2SPx);
     fChain->SetBranchAddress("MCpsi2SPy", &MCpsi2SPy, &b_MCpsi2SPy);
     fChain->SetBranchAddress("MCpsi2SPz", &MCpsi2SPz, &b_MCpsi2SPz);
     fChain->SetBranchAddress("MCmupPx", &MCmupPx, &b_MCmupPx);
     fChain->SetBranchAddress("MCmupPy", &MCmupPy, &b_MCmupPy);
     fChain->SetBranchAddress("MCmupPz", &MCmupPz, &b_MCmupPz);
     fChain->SetBranchAddress("MCmumPx", &MCmumPx, &b_MCmumPx);
     fChain->SetBranchAddress("MCmumPy", &MCmumPy, &b_MCmumPy);
     fChain->SetBranchAddress("MCmumPz", &MCmumPz, &b_MCmumPz);
     fChain->SetBranchAddress("MCpionPx", &MCpionPx, &b_MCpionPx);
     fChain->SetBranchAddress("MCpionPy", &MCpionPy, &b_MCpionPy);
     fChain->SetBranchAddress("MCpionPz", &MCpionPz, &b_MCpionPz);
     fChain->SetBranchAddress("MCkaonPx", &MCkaonPx, &b_MCkaonPx);
     fChain->SetBranchAddress("MCkaonPy", &MCkaonPy, &b_MCkaonPy);
     fChain->SetBranchAddress("MCkaonPz", &MCkaonPz, &b_MCkaonPz);
     fChain->SetBranchAddress("MCPx", &MCPx, &b_MCPx);
     fChain->SetBranchAddress("MCPy", &MCPy, &b_MCPy);
     fChain->SetBranchAddress("MCPz", &MCPz, &b_MCPz);
     fChain->SetBranchAddress("MCpionCh", &MCpionCh, &b_MCpionCh);
     fChain->SetBranchAddress("MCkaonCh", &MCkaonCh, &b_MCkaonCh);
   }
   // Generic Muons
   fChain->SetBranchAddress("nMu", &nMu, &b_nMu);
   fChain->SetBranchAddress("muPx", &muPx, &b_muPx);
   fChain->SetBranchAddress("muPy", &muPy, &b_muPy);
   fChain->SetBranchAddress("muPz", &muPz, &b_muPz);
   fChain->SetBranchAddress("muD0", &muD0, &b_muD0);
   fChain->SetBranchAddress("muDz", &muDz, &b_muDz);
   fChain->SetBranchAddress("muChi2", &muChi2, &b_muChi2);
   fChain->SetBranchAddress("muNDF", &muNDF, &b_muNDF);
   fChain->SetBranchAddress("muPhits", &muPhits, &b_muPhits);
   fChain->SetBranchAddress("muShits", &muShits, &b_muShits);
   fChain->SetBranchAddress("muLayersTr", &muLayersTr, &b_muLayersTr);
   fChain->SetBranchAddress("muLayersPix", &muLayersPix, &b_muLayersPix);
   fChain->SetBranchAddress("muD0E", &muD0E, &b_muD0E);
   fChain->SetBranchAddress("muDzVtxErr", &muDzVtxErr, &b_muDzVtxErr);
   fChain->SetBranchAddress("muKey", &muKey, &b_muKey);
   fChain->SetBranchAddress("muGlMuHits", &muGlMuHits, &b_muGlMuHits);
   fChain->SetBranchAddress("muGlChi2", &muGlChi2, &b_muGlChi2);
   fChain->SetBranchAddress("muGlNDF", &muGlNDF, &b_muGlNDF);
   fChain->SetBranchAddress("muGlMatchedStation", &muGlMatchedStation, &b_muGlMatchedStation);
   fChain->SetBranchAddress("muGlDzVtx", &muGlDzVtx, &b_muGlDzVtx);
   fChain->SetBranchAddress("muGlDxyVtx", &muGlDxyVtx, &b_muGlDxyVtx);
   fChain->SetBranchAddress("muType", &muType, &b_muType);
   fChain->SetBranchAddress("muQual", &muQual, &b_muQual);
   fChain->SetBranchAddress("muTrack", &muTrack, &b_muTrack);
   fChain->SetBranchAddress("muCharge", &muCharge, &b_muCharge);
   fChain->SetBranchAddress("mufHits", &mufHits, &b_mufHits);
   fChain->SetBranchAddress("muFirstBarrel", &muFirstBarrel, &b_muFirstBarrel);
   fChain->SetBranchAddress("muFirstEndCap", &muFirstEndCap, &b_muFirstEndCap);
   fChain->SetBranchAddress("muDzVtx", &muDzVtx, &b_muDzVtx);
   fChain->SetBranchAddress("muDxyVtx", &muDxyVtx, &b_muDxyVtx);
   if (!oldNtuple) {
     fChain->SetBranchAddress("trNotRef", &trNotRef, &b_trNotRef);
     fChain->SetBranchAddress("trRef", &trRef, &b_trRef);
   }
   fChain->SetBranchAddress("trackPx", &trackPx, &b_trackPx);
   fChain->SetBranchAddress("trackPy", &trackPy, &b_trackPy);
   fChain->SetBranchAddress("trackPz", &trackPz, &b_trackPz);
   fChain->SetBranchAddress("trackEnergy", &trackEnergy, &b_trackEnergy);
   fChain->SetBranchAddress("trackNDF", &trackNDF, &b_trackNDF);
   fChain->SetBranchAddress("trackPhits", &trackPhits, &b_trackPhits);
   fChain->SetBranchAddress("trackShits", &trackShits, &b_trackShits);
   fChain->SetBranchAddress("trackChi2", &trackChi2, &b_trackChi2);
   fChain->SetBranchAddress("trackD0", &trackD0, &b_trackD0);
   fChain->SetBranchAddress("trackD0Err", &trackD0Err, &b_trackD0Err);
   fChain->SetBranchAddress("trackCharge", &trackCharge, &b_trackCharge);
   fChain->SetBranchAddress("TrackHighPurity", &TrackHighPurity, &b_TrackHighPurity);
   fChain->SetBranchAddress("TrackTight", &TrackTight, &b_TrackTight);
   fChain->SetBranchAddress("trackfHits", &trackfHits, &b_trackfHits);
   fChain->SetBranchAddress("trackFirstBarrel", &trackFirstBarrel, &b_trackFirstBarrel);
   fChain->SetBranchAddress("trackFirstEndCap", &trackFirstEndCap, &b_trackFirstEndCap);
   fChain->SetBranchAddress("trackDzVtx", &trackDzVtx, &b_trackDzVtx);
   fChain->SetBranchAddress("trackDxyVtx", &trackDxyVtx, &b_trackDxyVtx);
   fChain->SetBranchAddress("tr_nsigdedx", &tr_nsigdedx, &b_tr_nsigdedx);
   fChain->SetBranchAddress("tr_dedx", &tr_dedx, &b_tr_dedx);
   fChain->SetBranchAddress("tr_dedxMass", &tr_dedxMass, &b_tr_dedxMass);
   fChain->SetBranchAddress("tr_theo", &tr_theo, &b_tr_theo);
   fChain->SetBranchAddress("tr_sigma", &tr_sigma, &b_tr_sigma);
   if (!oldNtuple) {
     fChain->SetBranchAddress("tr_dedx_byHits", &tr_dedx_byHits, &b_tr_dedx_byHits);
     fChain->SetBranchAddress("tr_dedxErr_byHits", &tr_dedxErr_byHits, &b_tr_dedxErr_byHits);
     fChain->SetBranchAddress("tr_saturMeas_byHits", &tr_saturMeas_byHits, &b_tr_saturMeas_byHits);
   fChain->SetBranchAddress("tr_Meas_byHits", &tr_Meas_byHits, &b_tr_Meas_byHits);
   }
   if (!oldNtuple) {
     fChain->SetBranchAddress("nMuMu", &nMuMu, &b_nMuMu);
     fChain->SetBranchAddress("MuMuMass", &MuMuMass, &b_MuMuMass);
     fChain->SetBranchAddress("MuMuPx", &MuMuPx, &b_MuMuPx);
     fChain->SetBranchAddress("MuMuPy", &MuMuPy, &b_MuMuPy);
     fChain->SetBranchAddress("MuMuPz", &MuMuPz, &b_MuMuPz);
     fChain->SetBranchAddress("MuMuVtx_CL", &MuMuVtx_CL, &b_MuMuVtx_CL);
     fChain->SetBranchAddress("MuMuVtx_Chi2", &MuMuVtx_Chi2, &b_MuMuVtx_Chi2);
     fChain->SetBranchAddress("MuMuDecayVtx_X", &MuMuDecayVtx_X, &b_MuMuDecayVtx_X);
     fChain->SetBranchAddress("MuMuDecayVtx_Y", &MuMuDecayVtx_Y, &b_MuMuDecayVtx_Y);
     fChain->SetBranchAddress("MuMuDecayVtx_Z", &MuMuDecayVtx_Z, &b_MuMuDecayVtx_Z);
     fChain->SetBranchAddress("MuMuDecayVtx_XE", &MuMuDecayVtx_XE, &b_MuMuDecayVtx_XE);
     fChain->SetBranchAddress("MuMuDecayVtx_YE", &MuMuDecayVtx_YE, &b_MuMuDecayVtx_YE);
     fChain->SetBranchAddress("MuMuDecayVtx_ZE", &MuMuDecayVtx_ZE, &b_MuMuDecayVtx_ZE);
   } else {
     fChain->SetBranchAddress("nPsiPrime", &nMuMu, &b_nMuMu);
     fChain->SetBranchAddress("PsiPrimeMass", &MuMuMass, &b_MuMuMass);
     fChain->SetBranchAddress("PsiPrimePx", &MuMuPx, &b_MuMuPx);
     fChain->SetBranchAddress("PsiPrimePy", &MuMuPy, &b_MuMuPy);
     fChain->SetBranchAddress("PsiPrimePz", &MuMuPz, &b_MuMuPz);
     fChain->SetBranchAddress("PsiPrimeVtx_CL", &MuMuVtx_CL, &b_MuMuVtx_CL);
     fChain->SetBranchAddress("PsiPrimeVtx_Chi2", &MuMuVtx_Chi2, &b_MuMuVtx_Chi2);
     fChain->SetBranchAddress("PsiPrimeDecayVtx_X", &MuMuDecayVtx_X, &b_MuMuDecayVtx_X);
     fChain->SetBranchAddress("PsiPrimeDecayVtx_Y", &MuMuDecayVtx_Y, &b_MuMuDecayVtx_Y);
     fChain->SetBranchAddress("PsiPrimeDecayVtx_Z", &MuMuDecayVtx_Z, &b_MuMuDecayVtx_Z);
     fChain->SetBranchAddress("PsiPrimeDecayVtx_XE", &MuMuDecayVtx_XE, &b_MuMuDecayVtx_XE);
     fChain->SetBranchAddress("PsiPrimeDecayVtx_YE", &MuMuDecayVtx_YE, &b_MuMuDecayVtx_YE);
     fChain->SetBranchAddress("PsiPrimeDecayVtx_ZE", &MuMuDecayVtx_ZE, &b_MuMuDecayVtx_ZE);
   }
   fChain->SetBranchAddress("mu1Idx", &mu1Idx, &b_mu1Idx);
   fChain->SetBranchAddress("mu2Idx", &mu2Idx, &b_mu2Idx);
   fChain->SetBranchAddress("mu1Px_MuMu", &mu1Px_MuMu, &b_mu1Px_MuMu);
   fChain->SetBranchAddress("mu1Py_MuMu", &mu1Py_MuMu, &b_mu1Py_MuMu);
   fChain->SetBranchAddress("mu1Pz_MuMu", &mu1Pz_MuMu, &b_mu1Pz_MuMu);
   fChain->SetBranchAddress("mu1Chi_MuMu2", &mu1Chi_MuMu2, &b_mu1Chi_MuMu2);
   fChain->SetBranchAddress("mu1NDF_MuMu", &mu1NDF_MuMu, &b_mu1NDF_MuMu);
   fChain->SetBranchAddress("mu2Px_MuMu", &mu2Px_MuMu, &b_mu2Px_MuMu);
   fChain->SetBranchAddress("mu2Py_MuMu", &mu2Py_MuMu, &b_mu2Py_MuMu);
   fChain->SetBranchAddress("mu2Pz_MuMu", &mu2Pz_MuMu, &b_mu2Pz_MuMu);
   fChain->SetBranchAddress("mu2Chi2_MuMu", &mu2Chi2_MuMu, &b_mu2Chi2_MuMu);
   fChain->SetBranchAddress("mu2NDF_MuMu", &mu2NDF_MuMu, &b_mu2NDF_MuMu);
   if (!oldNtuple) {
     fChain->SetBranchAddress("MuMuType", &MuMuType, &b_MuMuType);
     fChain->SetBranchAddress("MuMuMuonTrigMatch", &MuMuMuonTrigMatch, &b_MuMuMuonTrigMatch);
   } else {
     fChain->SetBranchAddress("jtype", &MuMuType, &b_MuMuType);
     fChain->SetBranchAddress("PsiPrimeMuonTrigMatch", &MuMuMuonTrigMatch, &b_MuMuMuonTrigMatch); 
   }
   fChain->SetBranchAddress("PriVtxMuMuCorr_X", &PriVtxMuMuCorr_X, &b_PriVtxMuMuCorr_X);
   fChain->SetBranchAddress("PriVtxMuMuCorr_Y", &PriVtxMuMuCorr_Y, &b_PriVtxMuMuCorr_Y);
   fChain->SetBranchAddress("PriVtxMuMuCorr_Z", &PriVtxMuMuCorr_Z, &b_PriVtxMuMuCorr_Z);
   fChain->SetBranchAddress("PriVtxMuMuCorr_EX", &PriVtxMuMuCorr_EX, &b_PriVtxMuMuCorr_EX);
   fChain->SetBranchAddress("PriVtxMuMuCorr_EY", &PriVtxMuMuCorr_EY, &b_PriVtxMuMuCorr_EY);
   fChain->SetBranchAddress("PriVtxMuMuCorr_EZ", &PriVtxMuMuCorr_EZ, &b_PriVtxMuMuCorr_EZ);
   fChain->SetBranchAddress("PriVtxMuMuCorr_Chi2", &PriVtxMuMuCorr_Chi2, &b_PriVtxMuMuCorr_Chi2);
   fChain->SetBranchAddress("PriVtxMuMuCorr_CL", &PriVtxMuMuCorr_CL, &b_PriVtxMuMuCorr_CL);
   fChain->SetBranchAddress("PriVtxMuMuCorr_tracks", &PriVtxMuMuCorr_tracks, &b_PriVtxMuMuCorr_tracks);
   fChain->SetBranchAddress("nB0", &nB0, &b_nB0);
   fChain->SetBranchAddress("B0Mass", &B0Mass, &b_B0Mass);
   fChain->SetBranchAddress("B0Px", &B0Px, &b_B0Px);
   fChain->SetBranchAddress("B0Py", &B0Py, &b_B0Py);
   fChain->SetBranchAddress("B0Pz", &B0Pz, &b_B0Pz);
   fChain->SetBranchAddress("B0PxE", &B0PxE, &b_B0PxE);
   fChain->SetBranchAddress("B0PyE", &B0PyE, &b_B0PyE);
   fChain->SetBranchAddress("B0PzE", &B0PzE, &b_B0PzE);
   fChain->SetBranchAddress("B0Vtx_CL", &B0Vtx_CL, &b_B0Vtx_CL);
   fChain->SetBranchAddress("B0Vtx_Chi2", &B0Vtx_Chi2, &b_B0Vtx_Chi2);
   fChain->SetBranchAddress("B0DecayVtx_X", &B0DecayVtx_X, &b_B0DecayVtx_X);
   fChain->SetBranchAddress("B0DecayVtx_Y", &B0DecayVtx_Y, &b_B0DecayVtx_Y);
   fChain->SetBranchAddress("B0DecayVtx_Z", &B0DecayVtx_Z, &b_B0DecayVtx_Z);
   fChain->SetBranchAddress("B0DecayVtx_XE", &B0DecayVtx_XE, &b_B0DecayVtx_XE);
   fChain->SetBranchAddress("B0DecayVtx_YE", &B0DecayVtx_YE, &b_B0DecayVtx_YE);
   fChain->SetBranchAddress("B0DecayVtx_ZE", &B0DecayVtx_ZE, &b_B0DecayVtx_ZE);
   //
   fChain->SetBranchAddress("B0CosAlphaBS", &B0CosAlphaBS, &b_B0CosAlphaBS);
   fChain->SetBranchAddress("B0CosAlpha3DBS", &B0CosAlpha3DBS, &b_B0CosAlpha3DBS);
   fChain->SetBranchAddress("B0CTauBS", &B0CTauBS, &b_B0CTauBS);
   fChain->SetBranchAddress("B0CTauBSE", &B0CTauBSE, &b_B0CTauBSE);
   fChain->SetBranchAddress("B0LxyBS", &B0LxyBS, &b_B0LxyBS);
   //
   fChain->SetBranchAddress("B0CosAlphaPV", &B0CosAlphaPV, &b_B0CosAlphaPV);
   fChain->SetBranchAddress("B0CosAlpha3DPV", &B0CosAlpha3DPV, &b_B0CosAlpha3DPV);
   fChain->SetBranchAddress("B0CTauPV", &B0CTauPV, &b_B0CTauPV);
   fChain->SetBranchAddress("B0CTauPVE", &B0CTauPVE, &b_B0CTauPVE);
   fChain->SetBranchAddress("B0LxyPV", &B0LxyPV, &b_B0LxyPV);
   //
   fChain->SetBranchAddress("PriVtx_B0CosAlpha_X", &PriVtx_B0CosAlpha_X, &b_PriVtx_B0CosAlpha_X);
   fChain->SetBranchAddress("PriVtx_B0CosAlpha_Y", &PriVtx_B0CosAlpha_Y, &b_PriVtx_B0CosAlpha_Y);
   fChain->SetBranchAddress("PriVtx_B0CosAlpha_Z", &PriVtx_B0CosAlpha_Z, &b_PriVtx_B0CosAlpha_Z);
   fChain->SetBranchAddress("PriVtx_B0CosAlpha_EX", &PriVtx_B0CosAlpha_EX, &b_PriVtx_B0CosAlpha_EX);
   fChain->SetBranchAddress("PriVtx_B0CosAlpha_EY", &PriVtx_B0CosAlpha_EY, &b_PriVtx_B0CosAlpha_EY);
   fChain->SetBranchAddress("PriVtx_B0CosAlpha_EZ", &PriVtx_B0CosAlpha_EZ, &b_PriVtx_B0CosAlpha_EZ);
   fChain->SetBranchAddress("PriVtx_B0CosAlpha_Chi2", &PriVtx_B0CosAlpha_Chi2, &b_PriVtx_B0CosAlpha_Chi2);
   fChain->SetBranchAddress("PriVtx_B0CosAlpha_CL", &PriVtx_B0CosAlpha_CL, &b_PriVtx_B0CosAlpha_CL);
   fChain->SetBranchAddress("PriVtx_B0CosAlpha_tracks", &PriVtx_B0CosAlpha_tracks, &b_PriVtx_B0CosAlpha_tracks);
   fChain->SetBranchAddress("B0CosAlphaPVCosAlpha", &B0CosAlphaPVCosAlpha, &b_B0CosAlphaPVCosAlpha);
   fChain->SetBranchAddress("B0CTauPVCosAlpha", &B0CTauPVCosAlpha, &b_B0CTauPVCosAlpha);
   fChain->SetBranchAddress("B0CTauPVCosAlphaE", &B0CTauPVCosAlphaE, &b_B0CTauPVCosAlphaE);
   fChain->SetBranchAddress("B0LxyPVCosAlpha", &B0LxyPVCosAlpha, &b_B0LxyPVCosAlpha);
   // 3D
   fChain->SetBranchAddress("PriVtx_B0CosAlpha3D_X", &PriVtx_B0CosAlpha3D_X, &b_PriVtx_B0CosAlpha3D_X);
   fChain->SetBranchAddress("PriVtx_B0CosAlpha3D_Y", &PriVtx_B0CosAlpha3D_Y, &b_PriVtx_B0CosAlpha3D_Y);
   fChain->SetBranchAddress("PriVtx_B0CosAlpha3D_Z", &PriVtx_B0CosAlpha3D_Z, &b_PriVtx_B0CosAlpha3D_Z);
   fChain->SetBranchAddress("PriVtx_B0CosAlpha3D_EX", &PriVtx_B0CosAlpha3D_EX, &b_PriVtx_B0CosAlpha3D_EX);
   fChain->SetBranchAddress("PriVtx_B0CosAlpha3D_EY", &PriVtx_B0CosAlpha3D_EY, &b_PriVtx_B0CosAlpha3D_EY);
   fChain->SetBranchAddress("PriVtx_B0CosAlpha3D_EZ", &PriVtx_B0CosAlpha3D_EZ, &b_PriVtx_B0CosAlpha3D_EZ);
   fChain->SetBranchAddress("PriVtx_B0CosAlpha3D_Chi2", &PriVtx_B0CosAlpha3D_Chi2, &b_PriVtx_B0CosAlpha3D_Chi2);
   fChain->SetBranchAddress("PriVtx_B0CosAlpha3D_CL", &PriVtx_B0CosAlpha3D_CL, &b_PriVtx_B0CosAlpha3D_CL);
   fChain->SetBranchAddress("PriVtx_B0CosAlpha3D_tracks", &PriVtx_B0CosAlpha3D_tracks, &b_PriVtx_B0CosAlpha3D_tracks);
   fChain->SetBranchAddress("B0CosAlphaPVCosAlpha3D", &B0CosAlphaPVCosAlpha3D, &b_B0CosAlphaPVCosAlpha3D);
   fChain->SetBranchAddress("B0CTauPVCosAlpha3D", &B0CTauPVCosAlpha3D, &b_B0CTauPVCosAlpha3D);
   fChain->SetBranchAddress("B0CTauPVCosAlpha3DE", &B0CTauPVCosAlpha3DE, &b_B0CTauPVCosAlpha3DE);
   fChain->SetBranchAddress("B0LxyPVCosAlpha3D", &B0LxyPVCosAlpha3D, &b_B0LxyPVCosAlpha3D);
   //
   if (priVtxsInfo) {
     fChain->SetBranchAddress("B0LessPV_tracksPtSq", &B0LessPV_tracksPtSq, &b_B0LessPV_tracksPtSq);
     fChain->SetBranchAddress("B0LessPV_4tracksPtSq", &B0LessPV_4tracksPtSq, &b_B0LessPV_4tracksPtSq);
     fChain->SetBranchAddress("PriVtxB0Less_X", &PriVtxB0Less_X, &b_PriVtxB0Less_X);
     fChain->SetBranchAddress("PriVtxB0Less_Y", &PriVtxB0Less_Y, &b_PriVtxB0Less_Y);
     fChain->SetBranchAddress("PriVtxB0Less_Z", &PriVtxB0Less_Z, &b_PriVtxB0Less_Z);
     fChain->SetBranchAddress("PriVtxB0Less_EX", &PriVtxB0Less_EX, &b_PriVtxB0Less_EX);
     fChain->SetBranchAddress("PriVtxB0Less_EY", &PriVtxB0Less_EY, &b_PriVtxB0Less_EY);
     fChain->SetBranchAddress("PriVtxB0Less_EZ", &PriVtxB0Less_EZ, &b_PriVtxB0Less_EZ);
     fChain->SetBranchAddress("PriVtxB0Less_Chi2", &PriVtxB0Less_Chi2, &b_PriVtxB0Less_Chi2);
     fChain->SetBranchAddress("PriVtxB0Less_CL", &PriVtxB0Less_CL, &b_PriVtxB0Less_CL);
     fChain->SetBranchAddress("PriVtxB0Less_tracks", &PriVtxB0Less_tracks, &b_PriVtxB0Less_tracks);
     fChain->SetBranchAddress("B0CosAlphaB0LessPV", &B0CosAlphaB0LessPV, &b_B0CosAlphaB0LessPV);
     fChain->SetBranchAddress("B0CosAlpha3DB0LessPV", &B0CosAlpha3DB0LessPV, &b_B0CosAlpha3DB0LessPV);
     fChain->SetBranchAddress("B0CTauB0LessPV", &B0CTauB0LessPV, &b_B0CTauB0LessPV);
     fChain->SetBranchAddress("B0CTauB0LessPVE", &B0CTauB0LessPVE, &b_B0CTauB0LessPVE);
     fChain->SetBranchAddress("B0LxyB0LessPV", &B0LxyB0LessPV, &b_B0LxyB0LessPV);
     //	
     fChain->SetBranchAddress("PriVtxB0Less_B0CosAlpha_X", &PriVtxB0Less_B0CosAlpha_X, &b_PriVtxB0Less_B0CosAlpha_X);
     fChain->SetBranchAddress("PriVtxB0Less_B0CosAlpha_Y", &PriVtxB0Less_B0CosAlpha_Y, &b_PriVtxB0Less_B0CosAlpha_Y);
     fChain->SetBranchAddress("PriVtxB0Less_B0CosAlpha_Z", &PriVtxB0Less_B0CosAlpha_Z, &b_PriVtxB0Less_B0CosAlpha_Z);
     fChain->SetBranchAddress("PriVtxB0Less_B0CosAlpha_EX", &PriVtxB0Less_B0CosAlpha_EX, &b_PriVtxB0Less_B0CosAlpha_EX);
     fChain->SetBranchAddress("PriVtxB0Less_B0CosAlpha_EY", &PriVtxB0Less_B0CosAlpha_EY, &b_PriVtxB0Less_B0CosAlpha_EY);
     fChain->SetBranchAddress("PriVtxB0Less_B0CosAlpha_EZ", &PriVtxB0Less_B0CosAlpha_EZ, &b_PriVtxB0Less_B0CosAlpha_EZ);
     fChain->SetBranchAddress("PriVtxB0Less_B0CosAlpha_Chi2", &PriVtxB0Less_B0CosAlpha_Chi2, &b_PriVtxB0Less_B0CosAlpha_Chi2);
     fChain->SetBranchAddress("PriVtxB0Less_B0CosAlpha_CL", &PriVtxB0Less_B0CosAlpha_CL, &b_PriVtxB0Less_B0CosAlpha_CL);
     fChain->SetBranchAddress("PriVtxB0Less_B0CosAlpha_tracks", &PriVtxB0Less_B0CosAlpha_tracks, &b_PriVtxB0Less_B0CosAlpha_tracks);
     fChain->SetBranchAddress("B0LxyB0LessPVCosAlpha", &B0LxyB0LessPVCosAlpha, &b_B0LxyB0LessPVCosAlpha);
     fChain->SetBranchAddress("B0CosAlphaB0LessPVCosAlpha", &B0CosAlphaB0LessPVCosAlpha, &b_B0CosAlphaB0LessPVCosAlpha);
     fChain->SetBranchAddress("B0CTauB0LessPVCosAlpha", &B0CTauB0LessPVCosAlpha, &b_B0CTauB0LessPVCosAlpha);
     fChain->SetBranchAddress("B0CTauB0LessPVCosAlphaE", &B0CTauB0LessPVCosAlphaE, &b_B0CTauB0LessPVCosAlphaE);
     // 3D	
     fChain->SetBranchAddress("PriVtxB0Less_B0CosAlpha3D_X", &PriVtxB0Less_B0CosAlpha3D_X, &b_PriVtxB0Less_B0CosAlpha3D_X);
     fChain->SetBranchAddress("PriVtxB0Less_B0CosAlpha3D_Y", &PriVtxB0Less_B0CosAlpha3D_Y, &b_PriVtxB0Less_B0CosAlpha3D_Y);
     fChain->SetBranchAddress("PriVtxB0Less_B0CosAlpha3D_Z", &PriVtxB0Less_B0CosAlpha3D_Z, &b_PriVtxB0Less_B0CosAlpha3D_Z);
     fChain->SetBranchAddress("PriVtxB0Less_B0CosAlpha3D_EX", &PriVtxB0Less_B0CosAlpha3D_EX, &b_PriVtxB0Less_B0CosAlpha3D_EX);
     fChain->SetBranchAddress("PriVtxB0Less_B0CosAlpha3D_EY", &PriVtxB0Less_B0CosAlpha3D_EY, &b_PriVtxB0Less_B0CosAlpha3D_EY);
     fChain->SetBranchAddress("PriVtxB0Less_B0CosAlpha3D_EZ", &PriVtxB0Less_B0CosAlpha3D_EZ, &b_PriVtxB0Less_B0CosAlpha3D_EZ);
     fChain->SetBranchAddress("PriVtxB0Less_B0CosAlpha3D_Chi2", &PriVtxB0Less_B0CosAlpha3D_Chi2, &b_PriVtxB0Less_B0CosAlpha3D_Chi2);
     fChain->SetBranchAddress("PriVtxB0Less_B0CosAlpha3D_CL", &PriVtxB0Less_B0CosAlpha3D_CL, &b_PriVtxB0Less_B0CosAlpha3D_CL);
     fChain->SetBranchAddress("PriVtxB0Less_B0CosAlpha3D_tracks", &PriVtxB0Less_B0CosAlpha3D_tracks, &b_PriVtxB0Less_B0CosAlpha3D_tracks);
     fChain->SetBranchAddress("B0CosAlphaB0LessPVCosAlpha3D", &B0CosAlphaB0LessPVCosAlpha3D, &b_B0CosAlphaB0LessPVCosAlpha3D);
     fChain->SetBranchAddress("B0CTauB0LessPVCosAlpha3D", &B0CTauB0LessPVCosAlpha3D, &b_B0CTauB0LessPVCosAlpha3D);
     fChain->SetBranchAddress("B0CTauB0LessPVCosAlpha3DE", &B0CTauB0LessPVCosAlpha3DE, &b_B0CTauB0LessPVCosAlpha3DE);
     fChain->SetBranchAddress("B0LxyB0LessPVCosAlpha3D", &B0LxyB0LessPVCosAlpha3D, &b_B0LxyB0LessPVCosAlpha3D);
   }
   //
   fChain->SetBranchAddress("PriVtxB0Corr_X", &PriVtxB0Corr_X, &b_PriVtxB0Corr_X);
   fChain->SetBranchAddress("PriVtxB0Corr_Y", &PriVtxB0Corr_Y, &b_PriVtxB0Corr_Y);
   fChain->SetBranchAddress("PriVtxB0Corr_Z", &PriVtxB0Corr_Z, &b_PriVtxB0Corr_Z);
   fChain->SetBranchAddress("PriVtxB0Corr_EX", &PriVtxB0Corr_EX, &b_PriVtxB0Corr_EX);
   fChain->SetBranchAddress("PriVtxB0Corr_EY", &PriVtxB0Corr_EY, &b_PriVtxB0Corr_EY);
   fChain->SetBranchAddress("PriVtxB0Corr_EZ", &PriVtxB0Corr_EZ, &b_PriVtxB0Corr_EZ);
   fChain->SetBranchAddress("PriVtxB0Corr_Chi2", &PriVtxB0Corr_Chi2, &b_PriVtxB0Corr_Chi2);
   fChain->SetBranchAddress("PriVtxB0Corr_CL", &PriVtxB0Corr_CL, &b_PriVtxB0Corr_CL);
   fChain->SetBranchAddress("PriVtxB0Corr_tracks", &PriVtxB0Corr_tracks, &b_PriVtxB0Corr_tracks);
   //
   fChain->SetBranchAddress("B0CosAlphaPVX", &B0CosAlphaPVX, &b_B0CosAlphaPVX);
   fChain->SetBranchAddress("B0CTauPVX", &B0CTauPVX, &b_B0CTauPVX);
   fChain->SetBranchAddress("B0CTauPVXE", &B0CTauPVXE, &b_B0CTauPVXE);
   fChain->SetBranchAddress("B0LxyPVX", &B0LxyPVX, &b_B0LxyPVX);
   fChain->SetBranchAddress("B0CTauPVX_3D", &B0CTauPVX_3D, &b_B0CTauPVX_3D);
   fChain->SetBranchAddress("B0CTauPVX_3D_err", &B0CTauPVX_3D_err, &b_B0CTauPVX_3D_err);
   if (!oldNtuple) {
     fChain->SetBranchAddress("B0MuMuIdx", &B0MuMuIdx, &b_B0MuMuIdx);
   } else {
     fChain->SetBranchAddress("B0PsiPrimeIdx", &B0MuMuIdx, &b_B0MuMuIdx);
   }
   fChain->SetBranchAddress("B0PionIdx", &B0PionIdx, &b_B0PionIdx);
   fChain->SetBranchAddress("B0KaonIdx", &B0KaonIdx, &b_B0KaonIdx);
   fChain->SetBranchAddress("PiPiMass_err", &PiPiMass_err, &b_PiPiMass_err);
   fChain->SetBranchAddress("Muon1Px_MuMuPiK", &Muon1Px_MuMuPiK, &b_Muon1Px_MuMuPiK);
   fChain->SetBranchAddress("Muon1Py_MuMuPiK", &Muon1Py_MuMuPiK, &b_Muon1Py_MuMuPiK);
   fChain->SetBranchAddress("Muon1Pz_MuMuPiK", &Muon1Pz_MuMuPiK, &b_Muon1Pz_MuMuPiK);
   fChain->SetBranchAddress("Muon1E_MuMuPiK", &Muon1E_MuMuPiK, &b_Muon1E_MuMuPiK);
   fChain->SetBranchAddress("Muon2Px_MuMuPiK", &Muon2Px_MuMuPiK, &b_Muon2Px_MuMuPiK);
   fChain->SetBranchAddress("Muon2Py_MuMuPiK", &Muon2Py_MuMuPiK, &b_Muon2Py_MuMuPiK);
   fChain->SetBranchAddress("Muon2Pz_MuMuPiK", &Muon2Pz_MuMuPiK, &b_Muon2Pz_MuMuPiK);
   fChain->SetBranchAddress("Muon2E_MuMuPiK", &Muon2E_MuMuPiK, &b_Muon2E_MuMuPiK);
   fChain->SetBranchAddress("PionPx_MuMuPiK", &PionPx_MuMuPiK, &b_PionPx_MuMuPiK);
   fChain->SetBranchAddress("PionPy_MuMuPiK", &PionPy_MuMuPiK, &b_PionPy_MuMuPiK);
   fChain->SetBranchAddress("PionPz_MuMuPiK", &PionPz_MuMuPiK, &b_PionPz_MuMuPiK);
   fChain->SetBranchAddress("PionE_MuMuPiK", &PionE_MuMuPiK, &b_PionE_MuMuPiK);
   if (!oldNtuple) { 
     fChain->SetBranchAddress("pion_nsigdedx", &pion_nsigdedx, &b_pion_nsigdedx);
     fChain->SetBranchAddress("pion_dedx", &pion_dedx, &b_pion_dedx);
     fChain->SetBranchAddress("pion_dedxMass", &pion_dedxMass, &b_pion_dedxMass);
     fChain->SetBranchAddress("pion_theo", &pion_theo, &b_pion_theo);
     fChain->SetBranchAddress("pion_sigma", &pion_sigma, &b_pion_sigma);
     fChain->SetBranchAddress("pion_dedx_byHits", &pion_dedx_byHits, &b_pion_dedx_byHits);
     fChain->SetBranchAddress("pion_dedxErr_byHits", &pion_dedxErr_byHits, &b_pion_dedxErr_byHits);
     fChain->SetBranchAddress("pion_saturMeas_byHits", &pion_saturMeas_byHits, &b_pion_saturMeas_byHits);
     fChain->SetBranchAddress("pion_Meas_byHits", &pion_Meas_byHits, &b_pion_Meas_byHits);
   }
   fChain->SetBranchAddress("KaonPx_MuMuPiK", &KaonPx_MuMuPiK, &b_KaonPx_MuMuPiK);
   fChain->SetBranchAddress("KaonPy_MuMuPiK", &KaonPy_MuMuPiK, &b_KaonPy_MuMuPiK);
   fChain->SetBranchAddress("KaonPz_MuMuPiK", &KaonPz_MuMuPiK, &b_KaonPz_MuMuPiK);
   fChain->SetBranchAddress("KaonE_MuMuPiK", &KaonE_MuMuPiK, &b_KaonE_MuMuPiK);
   if (!oldNtuple) {
     fChain->SetBranchAddress("kaon_nsigdedx", &kaon_nsigdedx, &b_kaon_nsigdedx);
     fChain->SetBranchAddress("kaon_dedx", &kaon_dedx, &b_kaon_dedx);
     fChain->SetBranchAddress("kaon_dedxMass", &kaon_dedxMass, &b_kaon_dedxMass);
     fChain->SetBranchAddress("kaon_theo", &kaon_theo, &b_kaon_theo);
     fChain->SetBranchAddress("kaon_sigma", &kaon_sigma, &b_kaon_sigma);
     fChain->SetBranchAddress("kaon_dedx_byHits", &kaon_dedx_byHits, &b_kaon_dedx_byHits);
     fChain->SetBranchAddress("kaon_dedxErr_byHits", &kaon_dedxErr_byHits, &b_kaon_dedxErr_byHits);
     fChain->SetBranchAddress("kaon_saturMeas_byHits", &kaon_saturMeas_byHits, &b_kaon_saturMeas_byHits);
     fChain->SetBranchAddress("kaon_Meas_byHits", &kaon_Meas_byHits, &b_kaon_Meas_byHits);
   }
}

Bool_t psiPrimePiK_MC::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef psiPrimePiK_MC_cxx

//rsync -vut --existing psiPrimePiK_MC.h cristella@cmssusy.ba.infn.it:/cmshome/cristella/work/Z_analysis/exclusive/clean_14ott/CMSSW_5_3_22/src/UserCode/MuMuPiKPAT/test/sanjay/selector/psiPrimePiK_MC.h
