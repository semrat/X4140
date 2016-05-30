// -*- C++ -*-
//
// Package:    MuMuPiKPAT
// Class:      MuMuPiKPAT
// 
/**\class MuMuPiKPAT MuMuPiKPAT.cc myAnalyzers/MuMuPiKPAT/src/MuMuPiKPAT.cc

   Description: <one line class summary>
   Make rootTuple for JPsiPiPi reconstruction

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  Keith Ulmer   keith.ulmer@colorado.edu
//
//

#ifndef _MuMuPiKPAT_h
#define _MuMuPiKPAT_h

// system include files
#include <memory>

// user include files
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
//#include "DataFormats/PatCandidates/interface/GenericParticle.h" // for namespace pat

//
// class decleration
//

using std::vector;
using namespace edm;
using namespace reco;
using namespace std;

class MuMuPiKPAT : public edm::EDAnalyzer {
public:
  explicit MuMuPiKPAT(const edm::ParameterSet&);
  ~MuMuPiKPAT();
  
private:
  virtual void beginJob() ;
  virtual void beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
	
  InvariantMassFromVertex massCalculator;

  const reco::DeDxDataValueMap *energyLoss;
  Int_t iexception_dedx;
  // dE/dx hits
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

  // ----------member data ---------------------------
  std::string proccessName_;
  HLTConfigProvider hltConfig_;

  edm::InputTag hlTriggerResults_;
  std::map<std::string,int> *HLTTrig; // HLT trigger prescale for accepted paths

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

  bool addMuMulessPrimaryVertex_;
  double MuMuMinMass, MuMuMaxMass, JPsiMinMass, JPsiMaxMass, PsiPrimeMinMass, PsiPrimeMaxMass ;
  double MuMuTrackMaxDR, B0TrackMaxDR ;
  bool   UseB0DR ;
  double MuMuPiKMinMass, MuMuPiKMaxMass;
  bool addB0lessPrimaryVertex_;
  bool Debug_;

  std::string DeDxEstimator_, m_dEdxDiscrimTag, m_dEdxDiscrimTag_kaon ;

  TTree* Z_One_Tree_;

  unsigned int        runNum, evtNum, lumiNum;
  //
  vector<unsigned int>* trigRes;
  vector<std::string>*  trigNames;
  //
  vector<unsigned int>* L1TT;
  vector<std::string>*  MatchTriggerNames;

  unsigned int          nMu, nMuMu, nB0;
  unsigned int          nB0_pre0, nB0_pre1, nB0_pre2, nB0_pre3, nB0_pre4, nB0_pre5, nB0_pre6, nB0_pre7, nB0_pre8, nB0_pre9, nB0_pre10, nB0_pre11, nB0_pre12, nB0_pre13, nB0_pre14;

  int                   priVtx_n;
  float                 priVtx_X, priVtx_Y, priVtx_Z, priVtx_XE, priVtx_YE, priVtx_ZE, priVtx_NormChi2, priVtx_Chi2, priVtx_CL;
  int                   priVtx_tracks;
  float                 priVtx_tracksPtSq;
  //// Indices
  vector<int>           *mu1Idx, *mu2Idx;
  vector<int>           *MuMuType;
  vector<int>           *B0_MuMuIdx, *B0_piIdx, *B0_kIdx;

  //// MC Analysis
  // Gen Primary Vertex
  unsigned int          n_genEvtVtx;
  vector<float>         *genEvtVtx_X, *genEvtVtx_Y, *genEvtVtx_Z; 
  vector<int>           *genEvtVtx_particles;
  vector<int>           *n_B0Ancestors;

  unsigned int        nMCAll, nMCB0, nMCB0Vtx;
  vector<int>         *MCPdgIdAll, *MCDanNumAll;
  // Gen Primary Vertex 
  vector<float>       *PriVtxGen_X, *PriVtxGen_Y, *PriVtxGen_Z ; 
  vector<double>      *PriVtxGen_EX, *PriVtxGen_EY, *PriVtxGen_EZ ;
  vector<float>	      *PriVtxGen_Chi2, *PriVtxGen_CL, *PriVtxGen_Ndof;
  vector<int>         *PriVtxGen_tracks ;
  vector<float>       *MCpsi2SPx, *MCpsi2SPy, *MCpsi2SPz;
  vector<float>       *MCmupPx, *MCmupPy, *MCmupPz;
  vector<float>       *MCmumPx, *MCmumPy, *MCmumPz;
  vector<float>       *MCpionPx, *MCpionPy, *MCpionPz;
  vector<float>       *MCkaonPx, *MCkaonPy, *MCkaonPz;
  vector<int>         *MCpionCh, *MCkaonCh;
  vector<float>       *MCPx, *MCPy, *MCPz;

  //// Generic Muons
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

  //// Generic tracks
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
  // dE/dx hits
  vector<float>         *tr_dedx_byHits, *tr_dedxErr_byHits ;
  vector<int>           *tr_saturMeas_byHits, *tr_Meas_byHits ;

  //// MuMu
  vector<float>         *MuMuMass, *MuMuPx, *MuMuPy, *MuMuPz;
  vector<float>         *MuMuVtx_CL, *MuMuVtx_Chi2;
  vector<float>         *MuMuDecayVtx_X, *MuMuDecayVtx_Y, *MuMuDecayVtx_Z, *MuMuDecayVtx_XE, *MuMuDecayVtx_YE, *MuMuDecayVtx_ZE;
  vector<bool>          *MuMuMuonTrigMatch;
  //// Muons after MuMu fit
  vector<float>         *mu1_MuMu_Px, *mu1_MuMu_Py, *mu1_MuMu_Pz ;
  vector<float>         *mu1_MuMu_Chi2 ;
  vector<int>           *mu1_MuMu_NDF ;
  vector<float>         *mu2_MuMu_Px, *mu2_MuMu_Py, *mu2_MuMu_Pz ;
  vector<float>         *mu2_MuMu_Chi2 ;
  vector<int>           *mu2_MuMu_NDF ;
  // Primary Vertex with "MuMu correction"
  vector<int>           *PriVtxMuMuCorr_n;
  vector<float>         *PriVtxMuMuCorr_X, *PriVtxMuMuCorr_Y, *PriVtxMuMuCorr_Z ; 
  vector<double>        *PriVtxMuMuCorr_EX, *PriVtxMuMuCorr_EY, *PriVtxMuMuCorr_EZ ;
  vector<float>	        *PriVtxMuMuCorr_Chi2, *PriVtxMuMuCorr_CL;
  vector<int>           *PriVtxMuMuCorr_tracks ;
  vector<int>           *nTrk ;
  // B0 cand 
  vector<float>         *b0Mass, *b0Vtx_CL, *b0Vtx_Chi2;
  vector<float>         *b0Px, *b0Py, *b0Pz ;
  vector<double>        *b0PxE, *b0PyE, *b0PzE ;
  vector<float>         *b0DecayVtx_X, *b0DecayVtx_Y, *b0DecayVtx_Z ;
  vector<double>        *b0DecayVtx_XE, *b0DecayVtx_YE, *b0DecayVtx_ZE ;
  // Muons and tracks after B0 cand fit
  vector<float>         *mu1Px_MuMuPiK, *mu1Py_MuMuPiK, *mu1Pz_MuMuPiK, *mu1E_MuMuPiK ;
  vector<float>         *mu2Px_MuMuPiK, *mu2Py_MuMuPiK, *mu2Pz_MuMuPiK, *mu2E_MuMuPiK ;
  vector<float>         *piPx_MuMuPiK, *piPy_MuMuPiK, *piPz_MuMuPiK, *piE_MuMuPiK ;
  vector<double>        *pion_nsigdedx;
  vector<float>         *pion_dedx, *pion_dedxMass, *pion_theo, *pion_sigma;
  vector<float>         *pion_dedx_byHits, *pion_dedxErr_byHits ;
  vector<int>           *pion_saturMeas_byHits, *pion_Meas_byHits ;
  vector<float>         *kPx_MuMuPiK, *kPy_MuMuPiK, *kPz_MuMuPiK, *kE_MuMuPiK ;
  vector<double>        *kaon_nsigdedx;
  vector<float>         *kaon_dedx, *kaon_dedxMass, *kaon_theo, *kaon_sigma;
  vector<float>         *kaon_dedx_byHits, *kaon_dedxErr_byHits ;
  vector<int>           *kaon_saturMeas_byHits, *kaon_Meas_byHits ;

  // Primary Vertex with largest B0_cos(alpha)
  vector<int>           *PriVtx_B0CosAlpha_n;
  vector<float>         *PriVtx_B0CosAlpha_X, *PriVtx_B0CosAlpha_Y, *PriVtx_B0CosAlpha_Z ; 
  vector<double>        *PriVtx_B0CosAlpha_EX, *PriVtx_B0CosAlpha_EY, *PriVtx_B0CosAlpha_EZ ;
  vector<float>	        *PriVtx_B0CosAlpha_Chi2, *PriVtx_B0CosAlpha_CL;
  vector<int>           *PriVtx_B0CosAlpha_tracks ;
  vector<int>           *PriVtx_B0CosAlpha3D_n;
  vector<float>         *PriVtx_B0CosAlpha3D_X, *PriVtx_B0CosAlpha3D_Y, *PriVtx_B0CosAlpha3D_Z ; 
  vector<double>        *PriVtx_B0CosAlpha3D_EX, *PriVtx_B0CosAlpha3D_EY, *PriVtx_B0CosAlpha3D_EZ ;
  vector<float>	        *PriVtx_B0CosAlpha3D_Chi2, *PriVtx_B0CosAlpha3D_CL;
  vector<int>           *PriVtx_B0CosAlpha3D_tracks ;
  vector<float>         *B0LessPV_tracksPtSq, *B0LessPV_4tracksPtSq ;
  vector<int>           *PriVtxB0Less_n;
  vector<float>         *PriVtxB0Less_X, *PriVtxB0Less_Y, *PriVtxB0Less_Z ; 
  vector<double>        *PriVtxB0Less_EX, *PriVtxB0Less_EY, *PriVtxB0Less_EZ ;
  vector<float>	        *PriVtxB0Less_Chi2, *PriVtxB0Less_CL;
  vector<int>           *PriVtxB0Less_tracks ;
  vector<int>           *PriVtxB0Less_B0CosAlpha_n;
  vector<float>         *PriVtxB0Less_B0CosAlpha_X, *PriVtxB0Less_B0CosAlpha_Y, *PriVtxB0Less_B0CosAlpha_Z ; 
  vector<double>        *PriVtxB0Less_B0CosAlpha_EX, *PriVtxB0Less_B0CosAlpha_EY, *PriVtxB0Less_B0CosAlpha_EZ ;
  vector<float>	        *PriVtxB0Less_B0CosAlpha_Chi2, *PriVtxB0Less_B0CosAlpha_CL;
  vector<int>           *PriVtxB0Less_B0CosAlpha_tracks ;
  vector<int>           *PriVtxB0Less_B0CosAlpha3D_n;
  vector<float>         *PriVtxB0Less_B0CosAlpha3D_X, *PriVtxB0Less_B0CosAlpha3D_Y, *PriVtxB0Less_B0CosAlpha3D_Z ; 
  vector<double>        *PriVtxB0Less_B0CosAlpha3D_EX, *PriVtxB0Less_B0CosAlpha3D_EY, *PriVtxB0Less_B0CosAlpha3D_EZ ;
  vector<float>	        *PriVtxB0Less_B0CosAlpha3D_Chi2, *PriVtxB0Less_B0CosAlpha3D_CL;
  vector<int>           *PriVtxB0Less_B0CosAlpha3D_tracks ;
  // Primary Vertex with "B0 correction" // check name
  vector<int>           *PriVtxB0Corr_n;
  vector<float>         *PriVtxB0Corr_X, *PriVtxB0Corr_Y, *PriVtxB0Corr_Z;
  vector<double>        *PriVtxB0Corr_EX, *PriVtxB0Corr_EY, *PriVtxB0Corr_EZ;
  vector<float>	        *PriVtxB0Corr_Chi2, *PriVtxB0Corr_CL;
  vector<int>           *PriVtxB0Corr_tracks;
  // Lifetimes
  vector<double>        *b0CosAlphaBS, *b0CosAlpha3DBS, *b0CTauBS, *b0CTauBSE, *b0LxyBS, *b0LxyBSE, *b0LxyzBS, *b0LxyzBSE ;
  vector<double>        *b0CosAlphaPV, *b0CosAlpha3DPV, *b0CTauPV, *b0CTauPVE, *b0LxyPV, *b0LxyPVE, *b0LxyzPV, *b0LxyzPVE ;
  vector<double>        *b0CosAlphaPVCosAlpha, *b0CosAlpha3DPVCosAlpha, *b0CTauPVCosAlpha, *b0CTauPVCosAlphaE, *b0LxyPVCosAlpha, *b0LxyPVCosAlphaE, *b0LxyzPVCosAlpha, *b0LxyzPVCosAlphaE ;
  vector<double>        *b0CosAlphaPVCosAlpha3D, *b0CosAlpha3DPVCosAlpha3D, *b0CTauPVCosAlpha3D, *b0CTauPVCosAlpha3DE, *b0LxyPVCosAlpha3D, *b0LxyPVCosAlpha3DE, *b0LxyzPVCosAlpha3D, *b0LxyzPVCosAlpha3DE ;
  vector<double>        *b0CosAlphaB0LessPV, *b0CosAlpha3DB0LessPV, *b0CTauB0LessPV, *b0CTauB0LessPVE, *b0LxyB0LessPV, *b0LxyB0LessPVE, *b0LxyzB0LessPV, *b0LxyzB0LessPVE ;
  vector<double>        *b0CosAlphaB0LessPVCosAlpha, *b0CosAlpha3DB0LessPVCosAlpha, *b0CTauB0LessPVCosAlpha, *b0CTauB0LessPVCosAlphaE, *b0LxyB0LessPVCosAlpha, *b0LxyB0LessPVCosAlphaE, *b0LxyzB0LessPVCosAlpha, *b0LxyzB0LessPVCosAlphaE ;
  vector<double>        *b0CosAlphaB0LessPVCosAlpha3D, *b0CosAlpha3DB0LessPVCosAlpha3D, *b0CTauB0LessPVCosAlpha3D, *b0CTauB0LessPVCosAlpha3DE, *b0LxyB0LessPVCosAlpha3D, *b0LxyB0LessPVCosAlpha3DE, *b0LxyzB0LessPVCosAlpha3D, *b0LxyzB0LessPVCosAlpha3DE ;
  vector<double>        *b0CosAlphaPVX, *b0CTauPVX, *b0CTauPVXE, *b0LxyPVX, *b0LxyPVXE, *b0LxyzPVX, *b0LxyzPVXE ;
  vector<float>	        *b0CTauPVX_3D, *b0CTauPVX_3D_err;

  vector<float>         *PiPiMass_err;

};

#endif

// rsync -vut --existing interface/MuMuPiKPAT.h cristella@cmssusy.ba.infn.it:/cmshome/cristella/work/Z_analysis/exclusive/clean_14ott/CMSSW_5_3_22/src/UserCode/MuMuPiKPAT/interface/MuMuPiKPAT.h
