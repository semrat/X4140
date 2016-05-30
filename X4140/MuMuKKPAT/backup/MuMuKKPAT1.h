// -*- C++ -*-
//
// Package:    MuMuKKPAT
// Class:      MuMuKKPAT
// 
/**\class MuMuKKPAT MuMuKKPAT.cc myAnalyzers/MuMuKKPAT/src/MuMuKKPAT.cc

   Description: <one line class summary>
   Make rootTuple for JPsiKK reconstruction

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  Keith Ulmer   keith.ulmer@colorado.edu
//
//

#ifndef _MuMuKKPAT_h
#define _MuMuKKPAT_h

// system include files
#include <memory>

/// user include files
#include "../interface/VertexReProducer.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackFromFTSFactory.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/V0Candidate/interface/V0Candidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "RecoVertex/V0Producer/interface/V0Producer.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"

#include "TFile.h"
#include "TTree.h"

#include <vector>
#include <string>
#include "RecoVertex/VertexTools/interface/InvariantMassFromVertex.h"
#include "RecoVertex/PrimaryVertexProducer/interface/VertexHigherPtSquared.h"


///
/// class decleration
///

using std::vector;
using namespace edm;
using namespace reco;
using namespace std;

class MuMuKKPAT : public edm::EDAnalyzer {
public:
  explicit MuMuKKPAT(const edm::ParameterSet&);
  ~MuMuKKPAT();
  
private:
  virtual void beginJob() ;
  virtual void beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  InvariantMassFromVertex massCalculator;
  const reco::DeDxDataValueMap *energyLoss;
  Int_t iexception_dedx;
  /// dE/dx hits
  edm::ValueMap<reco::DeDxData> dEdxTrack, dEdxTrack_Kaon;
  bool isAbHadron(int pdgID);
  bool isAMixedbHadron(int pdgID, int momPdgID);
  std::pair<int, float> findCandMCInfo(reco::GenParticleRef genCand);
  virtual double getSigmaOfLogdEdx(double logde);
  virtual float  getEnergyLoss(const reco::TrackRef & track);
  virtual double nsigmaofdedx(const reco::TrackRef & track, double & theo, double & sigma);
  virtual double getLogdEdx(double bg);
  virtual double GetMass(const reco::TrackRef & track);
  bool isSameMuon(const reco::Muon &mu1, const reco::Muon &mu2) const ;
  template<typename T> bool isBetterMuon(const T &mu1, const T &mu2) const ;

  /// ----------member data ---------------------------
  std::string proccessName_;
  HLTConfigProvider hltConfig_;

  edm::InputTag hlTriggerResults_;
  std::map<std::string,int> *HLTTrig; /// HLT trigger prescale for accepted paths

  edm::InputTag inputGEN_;
  std::string vtxSample;
  bool doData, doMC;
  int  MCParticle;
  bool MCExclusiveDecay;
  int  MCMother, MCDaughtersN;
  bool doMuMuMassConst;
  bool skipJPsi, skipPsi2S;
  int MuMinPixHits, MuMinSiHits;
  double MuMaxNormChi;
  double MuMaxD0;
  bool sharedFraction;
  int    TrMinSiHits;
  double TrMinPt;
  double TrMaxNormChi2;
  vector<string> TriggersForMatching_, FiltersForMatching_;
  bool resolveAmbiguity_; 
  int  MatchingTriggerResult[50];
  bool   addMuMulessPrimaryVertex_;
  double MuMuMinMass, MuMuMaxMass, JPsiMinMass, JPsiMaxMass; 
  //double KKMinMass, KKMaxMass, PhiMinMass, PhiMaxMass; 
  double JPsiPhiMaxXMass, JPsiPhiMinBs0Mass, JPsiPhiMaxBs0Mass; 
  double MuMuTrackMaxDR, Bs0TrackMaxDR; 
  bool   UseBs0DR ; 
  //bool   UseXDR ; 
  //double XTrackMaxDR; 
  double MuMuKKMinMass, MuMuKKMaxMass;
  bool   addBs0lessPrimaryVertex_; 
  bool   Debug_;
  std::string DeDxEstimator_, m_dEdxDiscrimTag, m_dEdxDiscrimTag_kaon ;
  TTree* Bs0_One_Tree_;
  //TTree* X_One_Tree_;
  unsigned int        runNum, evtNum, lumiNum;
  vector<unsigned int>* trigRes;
  vector<std::string>*  trigNames;
  vector<unsigned int>* L1TT;
  vector<std::string>*  MatchTriggerNames;

  /// counters for Bs0 & X(4140)
  unsigned int          nMu, nMuMu, nBs0; 
  unsigned int          nBs0_pre0, nBs0_pre1, nBs0_pre2, nBs0_pre3, nBs0_pre4, nBs0_pre5, nBs0_pre6, nBs0_pre7, nBs0_pre8, nBs0_pre9, nBs0_pre10, nBs0_pre11, nBs0_pre12, nBs0_pre13, nBs0_pre14; 
  //unsigned int          nX; 
  //unsigned int          nX_pre0, nX_pre1, nX_pre2, nX_pre3, nX_pre4, nX_pre5, nX_pre6, nX_pre7, nX_pre8, nX_pre9, nX_pre10, nX_pre11, nX_pre12, nX_pre13, nX_pre14;
  int                   priVtx_n;
  float                 priVtx_X, priVtx_Y, priVtx_Z, priVtx_XE, priVtx_YE, priVtx_ZE, priVtx_NormChi2, priVtx_Chi2, priVtx_CL;
  int                   priVtx_tracks;
  float                 priVtx_tracksPtSq;
  /// Indices
  vector<int>           *mu1Idx, *mu2Idx;
  vector<int>           *MuMuType;
  vector<int>           *Bs0_MuMuIdx, *Bs0_k1Idx, *Bs0_k2Idx; 
  //vector<int>           *X_MuMuIdx, *X_k1Idx, *X_k2Idx; 
  /// MC Analysis
  // Gen Primary Vertex
  unsigned int          n_genEvtVtx;
  vector<float>         *genEvtVtx_X, *genEvtVtx_Y, *genEvtVtx_Z; 
  vector<int>           *genEvtVtx_particles;
  vector<int>           *n_Bs0Ancestors; 
  unsigned int          nMCAll, nMCBs0, nMCBs0Vtx; 
  vector<int>           *MCPdgIdAll, *MCDanNumAll;
  // Gen Primary Vertex 
  vector<float>       *PriVtxGen_X, *PriVtxGen_Y, *PriVtxGen_Z ; 
  vector<double>      *PriVtxGen_EX, *PriVtxGen_EY, *PriVtxGen_EZ ;
  vector<float>	      *PriVtxGen_Chi2, *PriVtxGen_CL, *PriVtxGen_Ndof;
  vector<int>         *PriVtxGen_tracks ;
  //vector<float>       *MCpsi2SPx, *MCpsi2SPy, *MCpsi2SPz;
  vector<float>       *MCmupPx, *MCmupPy, *MCmupPz;
  vector<float>       *MCmumPx, *MCmumPy, *MCmumPz;
  vector<float>       *MCpionPx, *MCpionPy, *MCpionPz;
  vector<float>       *MCkaonPx, *MCkaonPy, *MCkaonPz;
  vector<int>         *MCpionCh, *MCkaonCh;
  vector<float>       *MCPx, *MCPy, *MCPz;
  /// Generic Muons
  vector<float>         *muPx, *muPy, *muPz, *muCharge;
  vector<int>           *muPhits, *muShits, *muLayersTr, *muLayersPix;
  vector<float>	        *muD0, *muD0E, *muDz, *muChi2 ;
  vector<int>           *muNDF;
  vector<float>         *mufHits;
  vector<bool>          *muFirstBarrel, *muFirstEndCap;
  vector<float>	        *muDzVtx, *muDxyVtx, *muDzVtxErr ;   
  vector<unsigned int>	*muKey;
  vector<bool> 	        *muIsGlobal, *muIsPF ;
  vector<int>           *muGlMuHits;
  vector<float>         *muGlChi2;
  vector<int>           *muGlNDF, *muGlMatchedStation;
  vector<float>         *muGlDzVtx, *muGlDxyVtx;
  vector<int>           *nMatchedStations;
  vector<int>           *muType, *muQual, *muTrack, *muNOverlap, *muNSharingSegWith;
  /// Generic tracks
  vector<float>         *trNotRef, *trRef;
  vector<float>         *trPx, *trPy, *trPz, *trE;
  vector<int>           *trNDF, *trPhits, *trShits;
  vector<float>         *trChi2;
  vector<float>         *trD0, *trD0E, *trCharge;
  vector<float>         *trfHits;
  vector<bool>          *trFirstBarrel, *trFirstEndCap;
  vector<float>         *trDzVtx, *trDxyVtx;
  vector<int>           *trQualityHighPurity, *trQualityTight;
  vector<double>        *tr_nsigdedx;
  vector<float>         *tr_dedx, *tr_dedxMass, *tr_theo, *tr_sigma;
  vector<float>         *tr_dedx_byHits, *tr_dedxErr_byHits ;
  vector<int>           *tr_saturMeas_byHits, *tr_Meas_byHits ;
  /// MuMu
  vector<float>         *MuMuMass, *MuMuPx, *MuMuPy, *MuMuPz;
  vector<float>         *MuMuVtx_CL, *MuMuVtx_Chi2;
  vector<float>         *MuMuDecayVtx_X, *MuMuDecayVtx_Y, *MuMuDecayVtx_Z, *MuMuDecayVtx_XE, *MuMuDecayVtx_YE, *MuMuDecayVtx_ZE;
  vector<bool>          *MuMuMuonTrigMatch;
  /// Muons after JPsi (MuMu) fit 
  vector<float>         *mu1_MuMu_Px, *mu1_MuMu_Py, *mu1_MuMu_Pz ;
  vector<float>         *mu1_MuMu_Chi2 ;
  vector<int>           *mu1_MuMu_NDF ;
  vector<float>         *mu2_MuMu_Px, *mu2_MuMu_Py, *mu2_MuMu_Pz ;
  vector<float>         *mu2_MuMu_Chi2 ;
  vector<int>           *mu2_MuMu_NDF ;
  /// Primary Vertex with "MuMu correction"
  vector<int>           *PriVtxMuMuCorr_n;
  vector<float>         *PriVtxMuMuCorr_X, *PriVtxMuMuCorr_Y, *PriVtxMuMuCorr_Z ; 
  vector<double>        *PriVtxMuMuCorr_EX, *PriVtxMuMuCorr_EY, *PriVtxMuMuCorr_EZ ;
  vector<float>	        *PriVtxMuMuCorr_Chi2, *PriVtxMuMuCorr_CL;
  vector<int>           *PriVtxMuMuCorr_tracks ;
  vector<int>           *nTrk ;
  /// Bs0 cand & X(4140) cand 
  vector<float>         *bs0Mass, *bs0Vtx_CL, *bs0Vtx_Chi2; 
  vector<float>         *bs0Px, *bs0Py, *bs0Pz ; 
  vector<double>        *bs0PxE, *bs0PyE, *bs0PzE ; 
  vector<float>         *bs0DecayVtx_X, *bs0DecayVtx_Y, *bs0DecayVtx_Z ; 
  vector<double>        *bs0DecayVtx_XE, *bs0DecayVtx_YE, *bs0DecayVtx_ZE ;
  //vector<float>         *xMass, *xVtx_CL, *xVtx_Chi2;
  //vector<float>         *xPx, *xPy, *xPz ;
  //vector<double>        *xPxE, *xPyE, *xPzE ;
  //vector<float>         *xDecayVtx_X, *xDecayVtx_Y, *xDecayVtx_Z ;
  //vector<double>        *xDecayVtx_XE, *xDecayVtx_YE, *xDecayVtx_ZE ; 
  /// Muons and tracks after Bs0 cand fit & X(4140) cand fit 
  vector<float>         *mu1Px_MuMuKK, *mu1Py_MuMuKK, *mu1Pz_MuMuKK, *mu1E_MuMuKK ;
  vector<float>         *mu2Px_MuMuKK, *mu2Py_MuMuKK, *mu2Pz_MuMuKK, *mu2E_MuMuKK ;
  vector<float>         *k1Px_MuMuKK, *k1Py_MuMuKK, *k1Pz_MuMuKK, *k1E_MuMuKK ;  
  vector<double>        *kaon1_nsigdedx; 
  vector<float>         *kaon1_dedx, *kaon1_dedxMass, *kaon1_theo, *kaon1_sigma ;
  vector<float>         *kaon1_dedx_byHits, *kaon1_dedxErr_byHits ; 
  vector<int>           *kaon1_saturMeas_byHits, *kaon1_Meas_byHits ; 
  vector<float>         *k2Px_MuMuKK, *k2Py_MuMuKK, *k2Pz_MuMuKK, *k2E_MuMuKK ;
  vector<double>        *kaon2_nsigdedx; 
  vector<float>         *kaon2_dedx, *kaon2_dedxMass, *kaon2_theo, *kaon2_sigma ; 
  vector<float>         *kaon2_dedx_byHits, *kaon2_dedxErr_byHits ; 
  vector<int>           *kaon2_saturMeas_byHits, *kaon2_Meas_byHits ; 
  //vector<float>         *X_mu1Px_MuMuKK, *X_mu1Py_MuMuKK, *X_mu1Pz_MuMuKK, *X_mu1E_MuMuKK ; 
  //vector<float>         *X_mu2Px_MuMuKK, *X_mu2Py_MuMuKK, *X_mu2Pz_MuMuKK, *X_mu2E_MuMuKK ;
  //vector<float>         *X_k1Px_MuMuKK, *X_k1Py_MuMuKK, *X_k1Pz_MuMuKK, *X_k1E_MuMuKK ; 
  //vector<double>        *X_kaon1_nsigdedx; 
  //vector<float>         *X_kaon1_dedx, *X_kaon1_dedxMass, *X_kaon1_theo, *X_kaon1_sigma; 
  //vector<float>         *X_kaon1_dedx_byHits, *X_kaon1_dedxErr_byHits ; 
  //vector<int>           *X_kaon1_saturMeas_byHits, *X_kaon1_Meas_byHits ; 
  //vector<float>         *X_k2Px_MuMuKK, *X_k2Py_MuMuKK, *X_k2Pz_MuMuKK, *X_k2E_MuMuKK ; 
  //vector<double>        *X_kaon2_nsigdedx; 
  //vector<float>         *X_kaon2_dedx, *X_kaon2_dedxMass, *X_kaon2_theo, *X_kaon2_sigma ;
  //vector<float>         *X_kaon2_dedx_byHits, *X_kaon2_dedxErr_byHits ; 
  //vector<int>           *X_kaon2_saturMeas_byHits, *X_kaon2_Meas_byHits ; 
  /// Primary Vertex with largest Bs0_cos(alpha) & largest X(4140)_cos(alpha)
  vector<int>           *PriVtx_Bs0CosAlpha_n; 
  vector<float>         *PriVtx_Bs0CosAlpha_X, *PriVtx_Bs0CosAlpha_Y, *PriVtx_Bs0CosAlpha_Z ; 
  vector<double>        *PriVtx_Bs0CosAlpha_EX, *PriVtx_Bs0CosAlpha_EY, *PriVtx_Bs0CosAlpha_EZ ; 
  vector<float>	        *PriVtx_Bs0CosAlpha_Chi2, *PriVtx_Bs0CosAlpha_CL; 
  vector<int>           *PriVtx_Bs0CosAlpha_tracks ; 
  vector<int>           *PriVtx_Bs0CosAlpha3D_n; 
  vector<float>         *PriVtx_Bs0CosAlpha3D_X, *PriVtx_Bs0CosAlpha3D_Y, *PriVtx_Bs0CosAlpha3D_Z ; 
  vector<double>        *PriVtx_Bs0CosAlpha3D_EX, *PriVtx_Bs0CosAlpha3D_EY, *PriVtx_Bs0CosAlpha3D_EZ ; 
  vector<float>	        *PriVtx_Bs0CosAlpha3D_Chi2, *PriVtx_Bs0CosAlpha3D_CL;
  vector<int>           *PriVtx_Bs0CosAlpha3D_tracks ;
  vector<float>         *Bs0LessPV_tracksPtSq, *Bs0LessPV_4tracksPtSq ;
  vector<int>           *PriVtxBs0Less_n;
  vector<float>         *PriVtxBs0Less_X, *PriVtxBs0Less_Y, *PriVtxBs0Less_Z ; 
  vector<double>        *PriVtxBs0Less_EX, *PriVtxBs0Less_EY, *PriVtxBs0Less_EZ ;
  vector<float>	        *PriVtxBs0Less_Chi2, *PriVtxBs0Less_CL;
  vector<int>           *PriVtxBs0Less_tracks ;
  vector<int>           *PriVtxBs0Less_Bs0CosAlpha_n;
  vector<float>         *PriVtxBs0Less_Bs0CosAlpha_X, *PriVtxBs0Less_Bs0CosAlpha_Y, *PriVtxBs0Less_Bs0CosAlpha_Z ; 
  vector<double>        *PriVtxBs0Less_Bs0CosAlpha_EX, *PriVtxBs0Less_Bs0CosAlpha_EY, *PriVtxBs0Less_Bs0CosAlpha_EZ ;
  vector<float>	        *PriVtxBs0Less_Bs0CosAlpha_Chi2, *PriVtxBs0Less_Bs0CosAlpha_CL;
  vector<int>           *PriVtxBs0Less_Bs0CosAlpha_tracks ;
  vector<int>           *PriVtxBs0Less_Bs0CosAlpha3D_n;
  vector<float>         *PriVtxBs0Less_Bs0CosAlpha3D_X, *PriVtxBs0Less_Bs0CosAlpha3D_Y, *PriVtxBs0Less_Bs0CosAlpha3D_Z ; 
  vector<double>        *PriVtxBs0Less_Bs0CosAlpha3D_EX, *PriVtxBs0Less_Bs0CosAlpha3D_EY, *PriVtxBs0Less_Bs0CosAlpha3D_EZ ;
  vector<float>	        *PriVtxBs0Less_Bs0CosAlpha3D_Chi2, *PriVtxBs0Less_Bs0CosAlpha3D_CL;
  vector<int>           *PriVtxBs0Less_Bs0CosAlpha3D_tracks ;
  //vector<int>           *PriVtx_XCosAlpha_n;
  //vector<float>         *PriVtx_XCosAlpha_X, *PriVtx_XCosAlpha_Y, *PriVtx_XCosAlpha_Z ;
  //vector<double>        *PriVtx_XCosAlpha_EX, *PriVtx_XCosAlpha_EY, *PriVtx_XCosAlpha_EZ ;
  //vector<float>         *PriVtx_XCosAlpha_Chi2, *PriVtx_XCosAlpha_CL;
  //vector<int>           *PriVtx_XCosAlpha_tracks ;
  //vector<int>           *PriVtx_XCosAlpha3D_n;
  //vector<float>         *PriVtx_XCosAlpha3D_X, *PriVtx_XCosAlpha3D_Y, *PriVtx_XCosAlpha3D_Z ;
  //vector<double>        *PriVtx_XCosAlpha3D_EX, *PriVtx_XCosAlpha3D_EY, *PriVtx_XCosAlpha3D_EZ ;
  //vector<float>         *PriVtx_XCosAlpha3D_Chi2, *PriVtx_XCosAlpha3D_CL;
  //vector<int>           *PriVtx_XCosAlpha3D_tracks ;
  /// Primary Vertex with "Bs0 correction" & "X(4140) correction"
  vector<int>           *PriVtxBs0Corr_n; 
  vector<float>         *PriVtxBs0Corr_X, *PriVtxBs0Corr_Y, *PriVtxBs0Corr_Z; 
  vector<double>        *PriVtxBs0Corr_EX, *PriVtxBs0Corr_EY, *PriVtxBs0Corr_EZ; 
  vector<float>	        *PriVtxBs0Corr_Chi2, *PriVtxBs0Corr_CL; 
  vector<int>           *PriVtxBs0Corr_tracks; 
  //vector<int>           *PriVtxXCorr_n;
  //vector<float>         *PriVtxXCorr_X, *PriVtxXCorr_Y, *PriVtxXCorr_Z;
  //vector<double>        *PriVtxXCorr_EX, *PriVtxXCorr_EY, *PriVtxXCorr_EZ;
  //vector<float>         *PriVtxXCorr_Chi2, *PriVtxXCorr_CL;
  //vector<int>           *PriVtxXCorr_tracks;
  /// Lifetimes variables for Bs0 & X(4140)
  vector<double>        *bs0CosAlphaBS, *bs0CosAlpha3DBS, *bs0CTauBS, *bs0CTauBSE, *bs0LxyBS, *bs0LxyBSE, *bs0LxyzBS, *bs0LxyzBSE ;
  vector<double>        *bs0CosAlphaPV, *bs0CosAlpha3DPV, *bs0CTauPV, *bs0CTauPVE, *bs0LxyPV, *bs0LxyPVE, *bs0LxyzPV, *bs0LxyzPVE ;
  vector<double>        *bs0CosAlphaPVCosAlpha, *bs0CosAlpha3DPVCosAlpha, *bs0CTauPVCosAlpha, *bs0CTauPVCosAlphaE, *bs0LxyPVCosAlpha, *bs0LxyPVCosAlphaE, *bs0LxyzPVCosAlpha, *bs0LxyzPVCosAlphaE ;
  vector<double>        *bs0CosAlphaPVCosAlpha3D, *bs0CosAlpha3DPVCosAlpha3D, *bs0CTauPVCosAlpha3D, *bs0CTauPVCosAlpha3DE, *bs0LxyPVCosAlpha3D, *bs0LxyPVCosAlpha3DE, *bs0LxyzPVCosAlpha3D, *bs0LxyzPVCosAlpha3DE ;
  vector<double>        *bs0CosAlphaBs0LessPV, *bs0CosAlpha3DBs0LessPV, *bs0CTauBs0LessPV, *bs0CTauBs0LessPVE, *bs0LxyBs0LessPV, *bs0LxyBs0LessPVE, *bs0LxyzBs0LessPV, *bs0LxyzBs0LessPVE ;
  vector<double>        *bs0CosAlphaBs0LessPVCosAlpha, *bs0CosAlpha3DBs0LessPVCosAlpha, *bs0CTauBs0LessPVCosAlpha, *bs0CTauBs0LessPVCosAlphaE, *bs0LxyBs0LessPVCosAlpha, *bs0LxyBs0LessPVCosAlphaE, *bs0LxyzBs0LessPVCosAlpha, *bs0LxyzBs0LessPVCosAlphaE ;
  vector<double>        *bs0CosAlphaBs0LessPVCosAlpha3D, *bs0CosAlpha3DBs0LessPVCosAlpha3D, *bs0CTauBs0LessPVCosAlpha3D, *bs0CTauBs0LessPVCosAlpha3DE, *bs0LxyBs0LessPVCosAlpha3D, *bs0LxyBs0LessPVCosAlpha3DE, *bs0LxyzBs0LessPVCosAlpha3D, *bs0LxyzBs0LessPVCosAlpha3DE ;
  vector<double>        *bs0CosAlphaPVX, *bs0CTauPVX, *bs0CTauPVXE, *bs0LxyPVX, *bs0LxyPVXE, *bs0LxyzPVX, *bs0LxyzPVXE ;
  vector<float>	        *bs0CTauPVX_3D, *bs0CTauPVX_3D_err;
  //vector<double>        *xCosAlphaBS, *xCosAlpha3DBS, *xCTauBS, *xCTauBSE, *xLxyBS, *xLxyBSE, *xLxyzBS, *xLxyzBSE ;
  //vector<double>        *xCosAlphaPV, *xCosAlpha3DPV, *xCTauPV, *xCTauPVE, *xLxyPV, *xLxyPVE, *xLxyzPV, *xLxyzPVE ;
  //vector<double>        *xCosAlphaPVCosAlpha, *xCosAlpha3DPVCosAlpha, *xCTauPVCosAlpha, *xCTauPVCosAlphaE, *xLxyPVCosAlpha, *xLxyPVCosAlphaE, *xLxyzPVCosAlpha, *xLxyzPVCosAlphaE ;  
  //vector<double>        *xCosAlphaPVCosAlpha3D, *xCosAlpha3DPVCosAlpha3D, *xCTauPVCosAlpha3D, *xCTauPVCosAlpha3DE, *xLxyPVCosAlpha3D, *xLxyPVCosAlpha3DE, *xLxyzPVCosAlpha3D, *xLxyzPVCosAlpha3DE ;  
  //vector<double>        *xCosAlphaPVX, *xCTauPVX, *xCTauPVXE, *xLxyPVX, *xLxyPVXE, *xLxyzPVX, *xLxyzPVXE ;
  //vector<float>         *xCTauPVX_3D, *xCTauPVX_3D_err;

  //vector<float>         *PiPiMass_err;

};

#endif

// rsync -vut --existing interface/MuMuPiKPAT.h semrat@lxplus.cern.ch:/afs/cern.ch/work/s/semrat/private/TetraQuark/CMSSW_5_3_22/src/X4140/MuMuKKPAT/interface
