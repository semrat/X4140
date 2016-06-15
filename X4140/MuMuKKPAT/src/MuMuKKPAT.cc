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


/// system include files
#include <memory>

/// user include files
#include "../interface/MuMuKKPAT.h"
#include "../interface/VertexReProducer.h"
//#include "DataFormats/Candidate/interface/OverlapChecker.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GtFdlWord.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"

#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

//#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Math/interface/Error.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/Vector3D.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"

/// for 53x 
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"

#include <vector>
#include <utility>

#include "TMath.h"
#include "Math/VectorUtil.h"

/// useless so far
//#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
//#include "HepMC/GenVertex.h"
//#include <HepMC/GenVertex.h>
//#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 

///
/// constants, enums and typedefs
///

typedef math::Error<3>::type CovarianceMatrix;

const ParticleMass muon_mass = 0.10565837; //pdg mass
const ParticleMass kaon_mass = 0.493667; //pdg mass
ParticleMass JPsi_mass = 3.096916;
const ParticleMass Phi_mass = 1.0194; 


/// Setting insignificant mass sigma to avoid singularities in the covariance matrix.
float small_sigma = muon_mass*1.e-6;				
//float small_sigma = kaon_mass*1.e-6; /// SEMRA

///
/// static data member definitions
///

///
/// constructors and destructor
///
MuMuKKPAT::MuMuKKPAT(const edm::ParameterSet& iConfig) :
  hlTriggerResults_(iConfig.getUntrackedParameter<edm::InputTag>("HLTriggerResults",edm::InputTag("TriggerResults::HLT")) ),
  inputGEN_(iConfig.getUntrackedParameter<edm::InputTag>("inputGEN",edm::InputTag("genParticles"))),
  vtxSample(iConfig.getUntrackedParameter<std::string>("VtxSample",std::string("offlinePrimaryVertices"))), 
  
  doData( iConfig.getUntrackedParameter<bool>("DoDataAnalysis", true) ),
  doMC( iConfig.getUntrackedParameter<bool>("DoMonteCarloTree", true) ),
  MCParticle( iConfig.getUntrackedParameter<int>("MonteCarloParticleId", 20443) ), /// 20443 X, 100443 Psi(2S), 9120443 X from B / decide later for X(4140)
  MCExclusiveDecay( iConfig.getUntrackedParameter<bool>("MonteCarloExclusiveDecay", true) ),
  MCMother( iConfig.getUntrackedParameter<int>("MonteCarloMotherId", 511) ), /// 511 B0 (=anti-B0), 531 Bs0 / decide later MCMotherId for X(4140)
  MCDaughtersN( iConfig.getUntrackedParameter<int>(" MonteCarloDaughtersN", 3) ), /// will be same  
  doMuMuMassConst( iConfig.getUntrackedParameter<bool>("DoMuMuMassConstraint", true) ),
  skipJPsi(iConfig.getUntrackedParameter<bool>("SkipJPsi", false)), 

  MuMinPixHits(iConfig.getUntrackedParameter<int>("MinNumMuPixHits", 0)),
  MuMinSiHits(iConfig.getUntrackedParameter<int>("MinNumMuSiHits", 0)),
  MuMaxNormChi(iConfig.getUntrackedParameter<double>("MaxMuNormChi2", 1000)),
  MuMaxD0(iConfig.getUntrackedParameter<double>("MaxMuD0", 1000)),
  sharedFraction(iConfig.getUntrackedParameter<double>("sharedFraction", 0.5)),

  TrMinSiHits(iConfig.getUntrackedParameter<int>("MinNumTrSiHits", 0)),
  TrMinPt(iConfig.getUntrackedParameter<double>("MinTrPt", 0)),
  TrMaxNormChi2(iConfig.getUntrackedParameter<double>("MaxTrChi2NDF", 10)),
  TriggersForMatching_(iConfig.getUntrackedParameter<std::vector<std::string> >("TriggersForMatching")),
  FiltersForMatching_(iConfig.getUntrackedParameter<std::vector<std::string> >("FiltersForMatching")),
  resolveAmbiguity_(iConfig.getUntrackedParameter<bool>("resolvePileUpAmbiguity",true)),
  addMuMulessPrimaryVertex_(iConfig.getUntrackedParameter<bool>("addMuMulessPrimaryVertex", true)),

  JPsiMinMass(iConfig.getUntrackedParameter<double>("MinJPsiMass", 2.8)), 
  JPsiMaxMass(iConfig.getUntrackedParameter<double>("MaxJPsiMass", 3.4)), 
  PhiMinMass(iConfig.getUntrackedParameter<double>("MinPhiMass", 0.97)), 
  PhiMaxMass(iConfig.getUntrackedParameter<double>("MaxPhiMass", 1.07)), 
  JPsiPhiMaxXMass(iConfig.getUntrackedParameter<double>("MaxJPsiPhiXMass", 4.35)), 
  JPsiPhiMinBs0Mass(iConfig.getUntrackedParameter<double>("MinJPsiPhiBs0Mass", 5.1)), 
  JPsiPhiMaxBs0Mass(iConfig.getUntrackedParameter<double>("MaxJPsiPhiBs0Mass", 5.6)), 
  MuMuTrackMaxDR(iConfig.getUntrackedParameter<double>("MaxMuMuTrackDR", 1)), 

  Bs0TrackMaxDR(iConfig.getUntrackedParameter<double>("MaxBs0CandTrackDR", 1.1)),  
  UseBs0DR(iConfig.getUntrackedParameter<bool>("UseBs0DR", false)),  
  UseXDR(iConfig.getUntrackedParameter<bool>("UseXDR", false)), 
  XTrackMaxDR(iConfig.getUntrackedParameter<double>("MaxXCandTrackDR", 1.1)),
  MuMuKKMinBs0Mass(iConfig.getUntrackedParameter<double>("MinMuMuKKBs0Mass", 0)),
  MuMuKKMaxBs0Mass(iConfig.getUntrackedParameter<double>("MaxMuMuKKBs0Mass", 10)), 
  MuMuKKMaxXMass(iConfig.getUntrackedParameter<double>("MaxMuMuKKXMass", 10)), 
  addBs0lessPrimaryVertex_(iConfig.getUntrackedParameter<bool>("addBs0lessPrimaryVertex", true)), 

  Debug_(iConfig.getUntrackedParameter<bool>("Debug_Output",true)), 
  DeDxEstimator_(iConfig.getUntrackedParameter<std::string>("DeDxEstimator", std::string("dedxHarmonic2"))),
  m_dEdxDiscrimTag(iConfig.getUntrackedParameter<std::string>("DeDxDiscriminator", std::string("dedxHarmonic2"))),
  m_dEdxDiscrimTag_kaon(iConfig.getUntrackedParameter<std::string>("DeDxDiscriminator", std::string("dedxHarmonic2"))),
  
  X_One_Tree_(0),
  runNum(0), evtNum(0), lumiNum(0),
  trigRes(0), trigNames(0), L1TT(0), MatchTriggerNames(0), 
  /// counters for Bs0 & X(4140)
  nMu(0), nMuMu(0), nBs0(0), nKK(0),  
  nBs0_pre0(0), nBs0_pre1(0), nBs0_pre2(0), nBs0_pre3(0), nBs0_pre4(0), nBs0_pre5(0), nBs0_pre6(0), nBs0_pre7(0), nBs0_pre8(0), nBs0_pre9(0), nBs0_pre10(0), nBs0_pre11(0), nBs0_pre12(0), nBs0_pre13(0), nBs0_pre14(0), nBs0_pre15(0), nBs0_pre16(0), 
  nX(0), 
  nX_pre0(0), nX_pre1(0), nX_pre2(0), nX_pre3(0), nX_pre4(0), nX_pre5(0), nX_pre6(0), nX_pre7(0), nX_pre8(0), nX_pre9(0), nX_pre10(0), nX_pre11(0), nX_pre12(0), nX_pre13(0), nX_pre14(0), nX_pre15(0), nX_pre16(0), 
  priVtx_n(0), priVtx_X(0), priVtx_Y(0), priVtx_Z(0), priVtx_XE(0), priVtx_YE(0), priVtx_ZE(0), priVtx_NormChi2(0), priVtx_Chi2(0), priVtx_CL(0), priVtx_tracks(0), priVtx_tracksPtSq(0),   
  /// indices
  mu1Idx(0), mu2Idx(0), MuMuType(0), ka1Idx(0), ka2Idx(0),
  Bs0_MuMuIdx(0), Bs0_ka1Idx(0), Bs0_ka2Idx(0), 
  X_MuMuIdx(0), X_ka1Idx(0), X_ka2Idx(0), 
  /// MC Analysis /// n_Bs0Ancestors & no for X(4140)!! 
  n_genEvtVtx(0), genEvtVtx_X(0), genEvtVtx_Y(0), genEvtVtx_Z(0), genEvtVtx_particles(0), n_Bs0Ancestors(0),
  nMCAll(0), nMCBs0(0), /*nMCBs0Vtx(0),*/ MCPdgIdAll(0), MCDanNumAll(0),
  // Gen Primary Vertex
  PriVtxGen_X(0), PriVtxGen_Y(0), PriVtxGen_Z(0), PriVtxGen_EX(0), PriVtxGen_EY(0), PriVtxGen_EZ(0),
  PriVtxGen_Chi2(0), PriVtxGen_CL(0), PriVtxGen_Ndof(0), PriVtxGen_tracks(0),
  MCJPsiPx(0), MCJPsiPy(0), MCJPsiPz(0),
  MCmupPx(0), MCmupPy(0), MCmupPz(0), 
  MCmumPx(0), MCmumPy(0), MCmumPz(0),
  MCPhiPx(0), MCPhiPy(0), MCPhiPz(0),
  MCkpPx(0), MCkpPy(0), MCkpPz(0),
  MCkmPx(0), MCkmPy(0), MCkmPz(0), 
  //MCpionPx(0), MCpionPy(0), MCpionPz(0), 
  //MCkaonPx(0), MCkaonPy(0), MCkaonPz(0),
  //MCpionCh(0), MCkaonCh(0),
  MCPx(0), MCPy(0), MCPz(0), 
  /// generic muons
  muPx(0), muPy(0), muPz(0), muCharge(0),
  muPhits(0), muShits(0), muLayersTr(0), muLayersPix(0),  
  muD0(0),  muD0E(0), muDz(0), muChi2(0), muNDF(0),
  mufHits(0), muFirstBarrel(0), muFirstEndCap(0),
  muDzVtx(0), muDxyVtx(0), muDzVtxErr(0), muKey(0),
  muIsGlobal(0), muIsPF(0),
  muGlMuHits(0), muGlChi2(0), muGlNDF(0), muGlMatchedStation(0), 
  muGlDzVtx(0), muGlDxyVtx(0),
  nMatchedStations(0), 
  muType(0), muQual(0), muTrack(0), muNOverlap(0), muNSharingSegWith(0), 
  /// generic tracks
  trNotRef(0), trRef(0),
  trPx(0), trPy(0), trPz(0), trE(0),
  trNDF(0), trPhits(0), trShits(0), trChi2(0),
  trD0(0), trD0E(0), trCharge(0),
  trfHits(0), trFirstBarrel(0), trFirstEndCap(0),
  trDzVtx(0), trDxyVtx(0),
  trQualityHighPurity(0), trQualityTight(0),
  tr_nsigdedx(0), tr_dedx(0), tr_dedxMass(0), tr_theo(0), tr_sigma(0),
  tr_dedx_byHits(0), tr_dedxErr_byHits(0), tr_saturMeas_byHits(0), tr_Meas_byHits(0),
  /// MuMu cand & KK cand
  MuMuMass(0), MuMuPx(0), MuMuPy(0), MuMuPz(0),
  MuMuVtx_CL(0), MuMuVtx_Chi2(0), 
  MuMuDecayVtx_X(0), MuMuDecayVtx_Y(0), MuMuDecayVtx_Z(0),
  MuMuDecayVtx_XE(0), MuMuDecayVtx_YE(0), MuMuDecayVtx_ZE(0),
  MuMuMuonTrigMatch(0),
  KKMass(0), KKPx(0), KKPy(0), KKPz(0),
  KKVtx_CL(0), KKVtx_Chi2(0),
  KKDecayVtx_X(0), KKDecayVtx_Y(0), KKDecayVtx_Z(0),
  KKDecayVtx_XE(0), KKDecayVtx_YE(0), KKDecayVtx_ZE(0),
  /// muons after JPsi (MuMu) fit &kaons after Phi (KK) fit 
  mu1_MuMu_Px(0), mu1_MuMu_Py(0), mu1_MuMu_Pz(0), mu1_MuMu_Chi2(0), mu1_MuMu_NDF(0),
  mu2_MuMu_Px(0), mu2_MuMu_Py(0), mu2_MuMu_Pz(0), mu2_MuMu_Chi2(0), mu2_MuMu_NDF(0),
  ka1_KK_Px(0), ka1_KK_Py(0), ka1_KK_Pz(0), ka1_KK_Chi2(0), ka1_KK_NDF(0),
  ka2_KK_Px(0), ka2_KK_Py(0), ka2_KK_Pz(0), ka2_KK_Chi2(0), ka2_KK_NDF(0),
  DRMuMuK1(0), DRMuMuK2(0), DRbs0K1(0), DRbs0K2(0), DRxK1(0), DRxK2(0),
  /// Primary Vertex with "MuMu correction"
  PriVtxMuMuCorr_n(0),
  PriVtxMuMuCorr_X(0), PriVtxMuMuCorr_Y(0), PriVtxMuMuCorr_Z(0), PriVtxMuMuCorr_EX(0), PriVtxMuMuCorr_EY(0), PriVtxMuMuCorr_EZ(0),
  PriVtxMuMuCorr_Chi2(0), PriVtxMuMuCorr_CL(0), PriVtxMuMuCorr_tracks(0),
  nTrk(0),
  /// Bs0 cand & X(4140) cand 
  bs0Mass(0), bs0Vtx_CL(0), bs0Vtx_Chi2(0),
  bs0Px(0), bs0Py(0), bs0Pz(0), bs0PxE(0), bs0PyE(0), bs0PzE(0), 
  bs0DecayVtx_X(0), bs0DecayVtx_Y(0), bs0DecayVtx_Z(0), bs0DecayVtx_XE(0), bs0DecayVtx_YE(0), bs0DecayVtx_ZE(0),
  xMass(0), xVtx_CL(0), xVtx_Chi2(0),
  xPx(0), xPy(0), xPz(0), xPxE(0), xPyE(0), xPzE(0),
  xDecayVtx_X(0), xDecayVtx_Y(0), xDecayVtx_Z(0), xDecayVtx_XE(0), xDecayVtx_YE(0), xDecayVtx_ZE(0),
  /// Muons and tracks after Bs0 cand fit & X(4140) cand fit 
  Bs0_mu1Px_MuMuKK(0), Bs0_mu1Py_MuMuKK(0), Bs0_mu1Pz_MuMuKK(0), Bs0_mu1E_MuMuKK(0),
  Bs0_mu2Px_MuMuKK(0), Bs0_mu2Py_MuMuKK(0), Bs0_mu2Pz_MuMuKK(0), Bs0_mu2E_MuMuKK(0),
  Bs0_k1Px_MuMuKK(0), Bs0_k1Py_MuMuKK(0), Bs0_k1Pz_MuMuKK(0), Bs0_k1E_MuMuKK(0), 
  Bs0_kaon1_nsigdedx(0), Bs0_kaon1_dedx(0), Bs0_kaon1_dedxMass(0), Bs0_kaon1_theo(0), Bs0_kaon1_sigma(0), 
  Bs0_kaon1_dedx_byHits(0), Bs0_kaon1_dedxErr_byHits(0), Bs0_kaon1_saturMeas_byHits(0), Bs0_kaon1_Meas_byHits(0), 
  Bs0_k2Px_MuMuKK(0), Bs0_k2Py_MuMuKK(0), Bs0_k2Pz_MuMuKK(0), Bs0_k2E_MuMuKK(0),
  Bs0_kaon2_nsigdedx(0), Bs0_kaon2_dedx(0), Bs0_kaon2_dedxMass(0), Bs0_kaon2_theo(0), Bs0_kaon2_sigma(0), 
  Bs0_kaon2_dedx_byHits(0), Bs0_kaon2_dedxErr_byHits(0), Bs0_kaon2_saturMeas_byHits(0), Bs0_kaon2_Meas_byHits(0), 
  X_mu1Px_MuMuKK(0), X_mu1Py_MuMuKK(0), X_mu1Pz_MuMuKK(0), X_mu1E_MuMuKK(0),
  X_mu2Px_MuMuKK(0), X_mu2Py_MuMuKK(0), X_mu2Pz_MuMuKK(0), X_mu2E_MuMuKK(0),
  X_k1Px_MuMuKK(0), X_k1Py_MuMuKK(0), X_k1Pz_MuMuKK(0), X_k1E_MuMuKK(0), 
  X_kaon1_nsigdedx(0), X_kaon1_dedx(0), X_kaon1_dedxMass(0), X_kaon1_theo(0), X_kaon1_sigma(0), 
  X_kaon1_dedx_byHits(0), X_kaon1_dedxErr_byHits(0), X_kaon1_saturMeas_byHits(0), X_kaon1_Meas_byHits(0), 
  X_k2Px_MuMuKK(0), X_k2Py_MuMuKK(0), X_k2Pz_MuMuKK(0), X_k2E_MuMuKK(0),
  X_kaon2_nsigdedx(0), X_kaon2_dedx(0), X_kaon2_dedxMass(0), X_kaon2_theo(0), X_kaon2_sigma(0), 
  X_kaon2_dedx_byHits(0), X_kaon2_dedxErr_byHits(0), X_kaon2_saturMeas_byHits(0), X_kaon2_Meas_byHits(0), 
  /// Primary Vertex with largest Bs0_cos(alpha) & largest X(4140)_cos(alpha) & no less values for X(4140)  
  PriVtx_Bs0CosAlpha_n(0),
  PriVtx_Bs0CosAlpha_X(0), PriVtx_Bs0CosAlpha_Y(0), PriVtx_Bs0CosAlpha_Z(0), PriVtx_Bs0CosAlpha_EX(0), PriVtx_Bs0CosAlpha_EY(0), PriVtx_Bs0CosAlpha_EZ(0),
  PriVtx_Bs0CosAlpha_Chi2(0), PriVtx_Bs0CosAlpha_CL(0), PriVtx_Bs0CosAlpha_tracks(0),
  PriVtx_Bs0CosAlpha3D_n(0),
  PriVtx_Bs0CosAlpha3D_X(0), PriVtx_Bs0CosAlpha3D_Y(0), PriVtx_Bs0CosAlpha3D_Z(0), PriVtx_Bs0CosAlpha3D_EX(0), PriVtx_Bs0CosAlpha3D_EY(0), PriVtx_Bs0CosAlpha3D_EZ(0),
  PriVtx_Bs0CosAlpha3D_Chi2(0), PriVtx_Bs0CosAlpha3D_CL(0), PriVtx_Bs0CosAlpha3D_tracks(0),
  Bs0LessPV_tracksPtSq(0), Bs0LessPV_4tracksPtSq(0),
  PriVtxBs0Less_n(0),
  PriVtxBs0Less_X(0), PriVtxBs0Less_Y(0), PriVtxBs0Less_Z(0), PriVtxBs0Less_EX(0), PriVtxBs0Less_EY(0), PriVtxBs0Less_EZ(0),
  PriVtxBs0Less_Chi2(0), PriVtxBs0Less_CL(0), PriVtxBs0Less_tracks(0),
  PriVtxBs0Less_Bs0CosAlpha_n(0),
  PriVtxBs0Less_Bs0CosAlpha_X(0), PriVtxBs0Less_Bs0CosAlpha_Y(0), PriVtxBs0Less_Bs0CosAlpha_Z(0), PriVtxBs0Less_Bs0CosAlpha_EX(0), PriVtxBs0Less_Bs0CosAlpha_EY(0), PriVtxBs0Less_Bs0CosAlpha_EZ(0),
  PriVtxBs0Less_Bs0CosAlpha_Chi2(0), PriVtxBs0Less_Bs0CosAlpha_CL(0), PriVtxBs0Less_Bs0CosAlpha_tracks(0),
  PriVtxBs0Less_Bs0CosAlpha3D_n(0),
  PriVtxBs0Less_Bs0CosAlpha3D_X(0), PriVtxBs0Less_Bs0CosAlpha3D_Y(0), PriVtxBs0Less_Bs0CosAlpha3D_Z(0), PriVtxBs0Less_Bs0CosAlpha3D_EX(0), PriVtxBs0Less_Bs0CosAlpha3D_EY(0), PriVtxBs0Less_Bs0CosAlpha3D_EZ(0),
  PriVtxBs0Less_Bs0CosAlpha3D_Chi2(0), PriVtxBs0Less_Bs0CosAlpha3D_CL(0), PriVtxBs0Less_Bs0CosAlpha3D_tracks(0),
  PriVtx_XCosAlpha_n(0),
  PriVtx_XCosAlpha_X(0), PriVtx_XCosAlpha_Y(0), PriVtx_XCosAlpha_Z(0), PriVtx_XCosAlpha_EX(0), PriVtx_XCosAlpha_EY(0), PriVtx_XCosAlpha_EZ(0),
  PriVtx_XCosAlpha_Chi2(0), PriVtx_XCosAlpha_CL(0), PriVtx_XCosAlpha_tracks(0),
  PriVtx_XCosAlpha3D_n(0),
  PriVtx_XCosAlpha3D_X(0), PriVtx_XCosAlpha3D_Y(0), PriVtx_XCosAlpha3D_Z(0), PriVtx_XCosAlpha3D_EX(0), PriVtx_XCosAlpha3D_EY(0), PriVtx_XCosAlpha3D_EZ(0),
  PriVtx_XCosAlpha3D_Chi2(0), PriVtx_XCosAlpha3D_CL(0), PriVtx_XCosAlpha3D_tracks(0),
  /// Primary Vertex with "Bs0 correction" & "X(4140) correction"
  PriVtxBs0Corr_n(0),
  PriVtxBs0Corr_X(0), PriVtxBs0Corr_Y(0), PriVtxBs0Corr_Z(0), PriVtxBs0Corr_EX(0), PriVtxBs0Corr_EY(0), PriVtxBs0Corr_EZ(0),
  PriVtxBs0Corr_Chi2(0), PriVtxBs0Corr_CL(0), PriVtxBs0Corr_tracks(0),
  PriVtxXCorr_n(0),  
  PriVtxXCorr_X(0), PriVtxXCorr_Y(0), PriVtxXCorr_Z(0), PriVtxXCorr_EX(0), PriVtxXCorr_EY(0), PriVtxXCorr_EZ(0),
  PriVtxXCorr_Chi2(0), PriVtxXCorr_CL(0), PriVtxXCorr_tracks(0),
  /// Lifetime variables for Bs0 & X(4140)
  bs0CosAlphaBS(0), bs0CosAlpha3DBS(0), bs0CTauBS(0), bs0CTauBSE(0), bs0LxyBS(0), bs0LxyBSE(0), bs0LxyzBS(0), bs0LxyzBSE(0),
  bs0CosAlphaPV(0), bs0CosAlpha3DPV(0), bs0CTauPV(0), bs0CTauPVE(0), bs0LxyPV(0), bs0LxyPVE(0), bs0LxyzPV(0), bs0LxyzPVE(0), 
  bs0CosAlphaPVCosAlpha(0), bs0CosAlpha3DPVCosAlpha(0), bs0CTauPVCosAlpha(0), bs0CTauPVCosAlphaE(0), bs0LxyPVCosAlpha(0), bs0LxyPVCosAlphaE(0), bs0LxyzPVCosAlpha(0), bs0LxyzPVCosAlphaE(0), 
  bs0CosAlphaPVCosAlpha3D(0), bs0CosAlpha3DPVCosAlpha3D(0), bs0CTauPVCosAlpha3D(0), bs0CTauPVCosAlpha3DE(0), bs0LxyPVCosAlpha3D(0), bs0LxyPVCosAlpha3DE(0), bs0LxyzPVCosAlpha3D(0), bs0LxyzPVCosAlpha3DE(0), 
  bs0CosAlphaBs0LessPV(0), bs0CosAlpha3DBs0LessPV(0),bs0CTauBs0LessPV(0), 
  bs0CTauBs0LessPVE(0), bs0LxyBs0LessPV(0), bs0LxyBs0LessPVE(0), bs0LxyzBs0LessPV(0), bs0LxyzBs0LessPVE(0), 
  bs0CosAlphaBs0LessPVCosAlpha(0), bs0CosAlpha3DBs0LessPVCosAlpha(0), bs0CTauBs0LessPVCosAlpha(0), bs0CTauBs0LessPVCosAlphaE(0), bs0LxyBs0LessPVCosAlpha(0), bs0LxyBs0LessPVCosAlphaE(0), bs0LxyzBs0LessPVCosAlpha(0), bs0LxyzBs0LessPVCosAlphaE(0), 
  bs0CosAlphaBs0LessPVCosAlpha3D(0), bs0CosAlpha3DBs0LessPVCosAlpha3D(0), bs0CTauBs0LessPVCosAlpha3D(0), bs0CTauBs0LessPVCosAlpha3DE(0), bs0LxyBs0LessPVCosAlpha3D(0), bs0LxyBs0LessPVCosAlpha3DE(0), bs0LxyzBs0LessPVCosAlpha3D(0), bs0LxyzBs0LessPVCosAlpha3DE(0), 
  bs0CosAlphaPVX(0), bs0CTauPVX(0), bs0CTauPVXE(0), bs0LxyPVX(0), bs0LxyzPVX(0), bs0LxyzPVXE(0), 
  bs0CTauPVX_3D(0), bs0CTauPVX_3D_err(0),
  xCosAlphaBS(0), xCosAlpha3DBS(0), xCTauBS(0), xCTauBSE(0), xLxyBS(0), xLxyBSE(0), xLxyzBS(0), xLxyzBSE(0),
  xCosAlphaPV(0), xCosAlpha3DPV(0), xCTauPV(0), xCTauPVE(0), xLxyPV(0), xLxyPVE(0), xLxyzPV(0), xLxyzPVE(0),
  xCosAlphaPVCosAlpha3D(0), xCosAlpha3DPVCosAlpha3D(0), xCTauPVCosAlpha3D(0), xCTauPVCosAlpha3DE(0), xLxyPVCosAlpha3D(0), xLxyPVCosAlpha3DE(0), xLxyzPVCosAlpha3D(0), xLxyzPVCosAlpha3DE(0),
  xCosAlphaPVX(0), xCTauPVX(0), xCTauPVXE(0), xLxyPVX(0), xLxyzPVX(0), xLxyzPVXE(0),
  xCTauPVX_3D(0), xCTauPVX_3D_err(0),
  /// dxy, dz, dxyE, dzE for kaons from PV, BS, Bs0LessPV
  kaon1_dxy_PV(0), kaon1_dz_PV(0), kaon2_dxy_PV(0), kaon2_dz_PV(0),
  kaon1_dxy_BS(0), kaon1_dz_BS(0), kaon2_dxy_BS(0), kaon2_dz_BS(0),
  kaon1_dxy_Bs0LessPV(0), kaon1_dz_Bs0LessPV(0), kaon2_dxy_Bs0LessPV(0), kaon2_dz_Bs0LessPV(0), 
  kaon1_dxyE(0), kaon1_dzE(0), kaon2_dxyE(0), kaon2_dzE(0),

  KKMass_err(0), Kaon1FromPV(0), Kaon2FromPV(0)

{
  /// now do what ever initialization is needed
  MuMuMinMass = JPsiMinMass;
  MuMuMaxMass = JPsiMaxMass;
  KKMinMass = PhiMinMass; 
  KKMaxMass = PhiMaxMass; 
  MuMuKKMinBs0Mass = JPsiPhiMinBs0Mass; 
  MuMuKKMaxBs0Mass = JPsiPhiMaxBs0Mass;
  MuMuKKMaxXMass = JPsiPhiMaxXMass;
}

MuMuKKPAT::~MuMuKKPAT()
{
  /// do anything here that needs to be done at desctruction time
  /// (e.g. close files, deallocate resources etc.)
	
}


///
/// member functions
///

/// ------------ method called to for each event  ------------
void MuMuKKPAT::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{	
  /// get event content information
  bool decayChainOK = false;	
  runNum = iEvent.id().run();
  evtNum = iEvent.id().event();
  lumiNum = iEvent.id().luminosityBlock();

  bool hasRequestedTrigger = false;
  ESHandle<MagneticField> bFieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);
	
  /// first get HLT results
  map<string,int> HLTPreScaleMap;
  edm::Handle<edm::TriggerResults> hltresults;
  try {
    iEvent.getByLabel(hlTriggerResults_, hltresults);
  }
  catch ( ... ) {
    cout << "Couldn't get handle on HLT Trigger!" << endl;
  }
  if (!hltresults.isValid()) {
    cout << "No Trigger Results!" << endl;
  } 
  else {
    int ntrigs = hltresults->size();
    if (ntrigs==0){
      cout << "No trigger name given in TriggerResults of the input " << endl;
    } 
		
    /// get hold of trigger names - based on TriggerResults object!
    edm::TriggerNames triggerNames_;
    triggerNames_ = iEvent.triggerNames(*hltresults);
    int ntriggers = TriggersForMatching_.size();
    for (int MatchTrig = 0; MatchTrig < ntriggers; MatchTrig++) { // initialize MatchingTriggerResult array
      MatchingTriggerResult[MatchTrig] = 0;
    }

    for (int itrig = 0; itrig < ntrigs; itrig++) {
      string trigName = triggerNames_.triggerName(itrig);
      int hltflag = (*hltresults)[itrig].accept();
      if (Debug_) if (hltflag) cout << trigName << " " <<hltflag <<endl;
      trigRes->push_back(hltflag);
      trigNames->push_back(trigName);
      
      int ntriggers = TriggersForMatching_.size();
      for (int MatchTrig = 0; MatchTrig < ntriggers; MatchTrig++) {
	if (TriggersForMatching_[MatchTrig] == triggerNames_.triggerName(itrig)){
	  MatchingTriggerResult[MatchTrig] = hltflag;
	  if (hltflag==1) hasRequestedTrigger = true;
	  break;
	}
      }
    }
    for (int MatchTrig = 0; MatchTrig<ntriggers; MatchTrig++){
            if (Debug_) cout << TriggersForMatching_[MatchTrig]<<endl;
      MatchTriggerNames->push_back(TriggersForMatching_[MatchTrig]);
    }
    
    ///
    /// Get HLT map : triggername associated with its prescale, saved only for accepted trigger
    ///
    for (unsigned int i=0; i<triggerNames_.size(); i++){
      if ( hltresults->accept(i) ) { //  save trigger info only for accepted paths
	/// get the prescale from the HLTConfiguration, initialized at beginRun
	int prescale = hltConfig_.prescaleValue(iEvent,iSetup,triggerNames_.triggerNames().at(i));
	if (Debug_) std::cout<<" HLT===> "<<triggerNames_.triggerNames().at(i)<<" prescale ="<<prescale<<std::endl;
	HLTPreScaleMap[triggerNames_.triggerNames().at(i)] = prescale;
      }
    }
    HLTTrig = &HLTPreScaleMap; // store in the branch

  } /// end valid trigger
	
  /// get L1 trigger info
  edm::ESHandle<L1GtTriggerMenu> menuRcd;
  iSetup.get<L1GtTriggerMenuRcd>().get(menuRcd) ;	
  edm::Handle< L1GlobalTriggerReadoutRecord > gtRecord;
  iEvent.getByLabel( edm::InputTag("gtDigis"), gtRecord);
  const DecisionWord dWord = gtRecord->decisionWord();  
  const TechnicalTriggerWord ttWord = gtRecord->technicalTriggerWord();
  for(unsigned int l1i = 0; l1i != ttWord.size(); ++l1i){
    L1TT->push_back(ttWord.at(l1i));
  }
	
  Vertex thePrimaryVtx, theBeamSpotVtx;
  math::XYZPoint RefVtx;
  Int_t thePrimaryVtx_multiplicity = -1 ;

  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);
  if ( beamSpotHandle.isValid() ) {
    beamSpot = *beamSpotHandle; 
    theBeamSpotVtx = Vertex(beamSpot.position(), beamSpot.covariance3D());
  }
  else cout << "No beam spot available from EventSetup" << endl;
		
  Handle<VertexCollection> recVtxs;
  iEvent.getByLabel(vtxSample, recVtxs);
  unsigned int nVtxTrks = 0;
  if ( recVtxs->begin() != recVtxs->end() ) {
    thePrimaryVtx_multiplicity = recVtxs->size() ;

    if (addMuMulessPrimaryVertex_ || addBs0lessPrimaryVertex_ || resolveAmbiguity_) { 
      //thePrimaryVtx = Vertex(*(recVtxs->begin()));
      cout <<"here" <<endl;
      thePrimaryVtx = *(recVtxs->begin());
    }
    else {
      for ( reco::VertexCollection::const_iterator vtx = recVtxs->begin(); vtx != recVtxs->end(); ++vtx) {
	if (nVtxTrks < vtx->tracksSize() ) {
	  nVtxTrks = vtx->tracksSize();
	  thePrimaryVtx = Vertex(*vtx);
	}				
      }
    }
  } else {
    thePrimaryVtx = Vertex(beamSpot.position(), beamSpot.covariance3D());
    thePrimaryVtx_multiplicity = 1 ;
  }
	
  edm::ESHandle<TransientTrackBuilder> theTTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);
  KalmanVertexFitter vtxFitter(true);

  RefVtx = thePrimaryVtx.position(); /// reference primary vertex choosen
  priVtx_n = thePrimaryVtx_multiplicity ;
  priVtx_X = (thePrimaryVtx.position().x()) ;
  priVtx_Y = (thePrimaryVtx.position().y()) ;
  priVtx_Z = (thePrimaryVtx.position().z()) ;
  priVtx_XE = (thePrimaryVtx.xError()) ;	
  priVtx_YE = (thePrimaryVtx.yError()) ;
  priVtx_ZE = (thePrimaryVtx.zError()) ;
  priVtx_NormChi2 = (thePrimaryVtx.normalizedChi2()) ;
  priVtx_Chi2 = thePrimaryVtx.chi2() ;
  priVtx_CL = ChiSquaredProbability( (double)(thePrimaryVtx.chi2()), (double)(thePrimaryVtx.ndof())) ;
  priVtx_tracks = thePrimaryVtx.tracksSize() ;
  VertexHigherPtSquared vertexHigherPtSquared ;
  priVtx_tracksPtSq = vertexHigherPtSquared.sumPtSquared(thePrimaryVtx) ;
	
  /// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /// try reconstruction without fitting modules
  /// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Handle< vector<pat::GenericParticle> > thePATTrackHandle; 
  iEvent.getByLabel("cleanPatTrackCands", thePATTrackHandle); /// container of tracks with pion mass hypothesis
  Handle< vector<pat::GenericParticle> > theKaonRefittedPATTrackHandle;
  iEvent.getByLabel("cleanPatTrackKaonCands", theKaonRefittedPATTrackHandle); /// container of tracks with kaon mass hypothesis

  for ( vector<pat::GenericParticle>::const_iterator TrackNotRefitted = thePATTrackHandle->begin(); TrackNotRefitted != thePATTrackHandle->end(); ++TrackNotRefitted ) { 
    for ( vector<pat::GenericParticle>::const_iterator TrackRefitted = theKaonRefittedPATTrackHandle->begin(); TrackRefitted != theKaonRefittedPATTrackHandle->end(); ++TrackRefitted ) { 
      if ( TrackNotRefitted->track().key() == TrackRefitted->track().key() ) {
	trNotRef->push_back( TrackNotRefitted->p() ) ;
	trRef->push_back( TrackRefitted->p() ) ;
	break ;
      }
    }
    break ;
  }

  Handle< vector<pat::Muon> > thePATMuonHandle;
  iEvent.getByLabel("patMuonsWithTrigger", thePATMuonHandle);

  Handle<reco::DeDxDataValueMap> elossCollection;
  energyLoss = 0;
  iexception_dedx = 0;
  try {
    iEvent.getByLabel(DeDxEstimator_, elossCollection);
    energyLoss = elossCollection.product();
  } catch ( cms::Exception& ex ) {
    if (evtNum < 100) edm::LogError("Analyzer") <<"Warning can't get collection with label: elossCollection";
    iexception_dedx = 1;
  }

  /// dE/dx hits
  Handle<edm::ValueMap<reco::DeDxData> > dEdxTrackHandle;
  try {
    iEvent.getByLabel(m_dEdxDiscrimTag, dEdxTrackHandle);
    dEdxTrack = *dEdxTrackHandle.product();
  } catch ( cms::Exception& ex ) {
    if (evtNum < 100) edm::LogError("Analyzer") <<"Warning can't get collection with label: dEdxTrackHandle";
    iexception_dedx = 1;
  }

  Handle<edm::ValueMap<reco::DeDxData> > dEdxTrackHandle_Kaon;
  try {
    iEvent.getByLabel(m_dEdxDiscrimTag_kaon, dEdxTrackHandle_Kaon);
    dEdxTrack_Kaon = *dEdxTrackHandle_Kaon.product();
  } catch ( cms::Exception& ex ) {
    if (evtNum < 100) edm::LogError("Analyzer") <<"Warning can't get collection with label: dEdxTrackHandle_Kaon";
    iexception_dedx = 1;
  }


  ////////////////// check MC truth //////////////////
  if (doMC) {
    /*
    // Get generated event
    //Handle<edm::HepMCProduct> hepEv;
    //iEvent.getByLabel("generator", hepEv);
    Handle<GenEventInfoProduct> genEvtInfo;
    iEvent.getByLabel("generator", genEvtInfo);

    //const HepMC::GenEvent *myGenEvent = hepEv->GetEvent();
    const HepMC::GenEvent *myGenEvent = genEvtInfo->GetEvent();
    n_genEvtVtx = myGenEvent->vertices_size() ;

    HepMC::GenVertex* primaryGenVtx = *(myGenEvent->vertices_begin()) ;

    genEvtVtx_X->push_back( primaryGenVtx->point3d().x() );
    genEvtVtx_Y->push_back( primaryGenVtx->point3d().y() );
    genEvtVtx_Z->push_back( primaryGenVtx->point3d().z() );
    //genEvtVtx_XE = (primaryGenVtx->xError()) ;	
    //genEvtVtx_YE = (primaryGenVtx->yError()) ;
    //genEvtVtx_ZE = (primaryGenVtx->zError()) ;
    //genEvtVtx_NormChi2 = (primaryGenVtx->normalizedChi2()) ;
    //genEvtVtx_Chi2 = primaryGenVtx->chi2() ;
    //genEvtVtx_CL = ChiSquaredProbability( (double)(primaryGenVtx.chi2()), (double)(primaryGenVtx.ndof())) ;
    genEvtVtx_particles->push_back( primaryGenVtx->particles_out_size() );
    */

    Handle< vector< PileupSummaryInfo > >  PupInfo;
    iEvent.getByLabel("addPileupInfo", PupInfo);
    vector<PileupSummaryInfo>::const_iterator PVI;
    if (Debug_) cout <<"\nBunchXing multiplicity = " <<PupInfo->size() <<endl ;
    for (PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI)
      if (Debug_) cout <<"Pileup Information: bunchXing, nvtx: " <<PVI->getBunchCrossing() <<" " <<PVI->getPU_NumInteractions() <<endl;

    Handle<GenParticleCollection> genParticles;
    iEvent.getByLabel("genParticles", genParticles);
    if (Debug_) cout << "############### GenParticles Analysis ###############" << endl;
    float jpsiPx=0., jpsiPy=0., jpsiPz=0.;
    float  mupPx=0., mupPy=0., mupPz=0., mumPx=0., mumPy=0., mumPz=0.; 
    float phiPx=0., phiPy=0., phiPz=0.;
    float  kpPx=0., kpPy=0., kpPz=0., kmPx=0., kmPy=0., kmPz=0.;
    //float pionPx=0., pionPy=0., pionPz=0., kaonPx=0., kaonPy=0., kaonPz=0.;
    //int pionCh=0, kaonCh=0 ;

    for (size_t i = 0; i < genParticles->size(); ++ i) {
      nMCAll++;
      const reco::GenParticle &p = (*genParticles)[i];
      int pdgid = p.pdgId() ; 
      int dauNum = p.numberOfDaughters();
      MCPdgIdAll->push_back( pdgid );
      MCDanNumAll->push_back( dauNum );

      if ( MCExclusiveDecay ) {
	/// check if there is a MCMother which has MCDaughtersN daughters
	if ( abs(pdgid) == MCMother  &&  dauNum == MCDaughtersN ) {
	  bool mumuOK = false; 
          bool kkOK = false;
          //bool pionOK = false, kaonOK = false;

	  for (int j=0; j<dauNum; ++j) {
	    const Candidate *dau = p.daughter(j);
	    if (Debug_) cout << "dauPdgId = " << dau->pdgId() << endl;
	  
	    /// check if one of B0 daughters is a psi(nS) whitch has 2 muons as daughters /// SEMRA ask again !!!
	    int mumuId = 0 ;
	    if (skipJPsi) /// SEMRA cleaned skipPsi2S
	      if (Debug_) cout <<"Skipping J/psi!" <<endl ; /// SEMRA cleaned skipPsi2S
	    //else if (skipPsi2S) /// SEMRA
	    //  mumuId = 443 ; /// SEMRA (JPsi ID)
 
	    if ( ((skipJPsi) && (dau->pdgId() == mumuId)) || 
		 ((!skipJPsi) && (dau->pdgId()%1000 == 443)) ) { 
	      jpsiPx = dau->px(); jpsiPy = dau->py(); jpsiPz = dau->pz(); 
	      int jpsiDauNum = dau->numberOfDaughters();
	      if (Debug_) cout << "jpsiDauNum = " << jpsiDauNum << endl;
	      int muNum = 0;	
	      for (int k=0; k<jpsiDauNum; ++k) {
		const Candidate *grandDau = dau->daughter(k); 
		if (Debug_)  cout << "grandDauPdgId = " << grandDau->pdgId() << endl;
		if ( abs(grandDau->pdgId()) == 13 ) {
		  muNum++;
		  if (grandDau->pdgId() < 0) {
		    mupPx = grandDau->px(); mupPy = grandDau->py(); mupPz = grandDau->pz();
		  } else {
		    mumPx = grandDau->px(); mumPy = grandDau->py(); mumPz = grandDau->pz();
		  }
		}
	      } 
	      if ( muNum == 2 ) mumuOK = true ;
	      
	    } /// end check if one of the MCMother daughters is a J/Psi or psi'
     
	     /// for Phi 
              phiPx = dau->px(); phiPy = dau->py(); phiPz = dau->pz();
              int phiDauNum = dau->numberOfDaughters();
              if (Debug_) cout << "phiDauNum = " << phiDauNum << endl;
              int kNum = 0;
              for (int n=0; n<phiDauNum; ++n) {
                const Candidate *grandDau = dau->daughter(n);
                if (Debug_)  cout << "grandDauPdgId = " << grandDau->pdgId() << endl;
                if ( abs(grandDau->pdgId()) == 321 ) {
                  kNum++;
                  if (grandDau->pdgId() < 0) {
                    kpPx = grandDau->px(); kpPy = grandDau->py(); kpPz = grandDau->pz();
                  } else {
                    kmPx = grandDau->px(); kmPy = grandDau->py(); kmPz = grandDau->pz();
                  }
                }
              }
	      if ( kNum == 2 ) kkOK = true ;

	    
	    /*else if ( abs(dau->pdgId()) == 211 ) { // check if one of B0 daughters is a pion /// SEMRA ask again !!!
	      pionPx = dau->px(); pionPy = dau->py(); pionPz = dau->pz();
	      pionCh = (dau->pdgId() == 211)? 1 : -1;
	      pionOK = true; /// SEMRA pions change with kaons for Bs0 ?
	    } else if ( abs(dau->pdgId()) == 321 ) { // check if one of B0 daughters is a kaon /// SEMRA ask again !!!
	      kaonPx = dau->px(); kaonPy=dau->py(); kaonPz=dau->pz();
	      kaonCh = (dau->pdgId() == 321)? 1 : -1;
	      kaonOK = true;
	    }*/
	    
	  } /// end loop on MCMother daughters

	  if (Debug_) cout << "mumuOK = " << mumuOK << ", kkOK = " << kkOK << endl;
	  if ( mumuOK && kkOK ) {
	    if (Debug_) {
	      cout <<"\nnumber of Bs0 mothers = " <<p.numberOfMothers() <<endl ; 
	      cout <<"Bs0 mother pdgID = " <<p.mother(0)->pdgId() <<endl ; 
	    }
	    ++nMCBs0 ;  
	      PriVtxGen_X->push_back( p.vx() ) ;
	      PriVtxGen_Y->push_back( p.vy() ) ; 
	      PriVtxGen_Z->push_back( p.vz() ) ;
	      PriVtxGen_CL->push_back( p.vertexNormalizedChi2() ) ;
	      PriVtxGen_Chi2->push_back( p.vertexChi2() ) ;
	      PriVtxGen_Ndof->push_back( p.vertexNdof() ) ;
	      
	      Bool_t status = kTRUE ;
	      const Candidate *bs0_ancestor = p.mother(0) ; /// a particle can have several mothers
	      Int_t n_ancestors = 1 ;
	      while ( status ) { 
		if ( abs(bs0_ancestor->pdgId()) <= 8 || bs0_ancestor->pdgId() == 21 || bs0_ancestor->status() == 3 ) {
		  status = kFALSE ;
		  if (Debug_) cout <<"Bs0 ancestor ID = " <<bs0_ancestor->pdgId() <<endl ; 
		  genEvtVtx_X->push_back( bs0_ancestor->daughter(0)->vx() ) ;
		  genEvtVtx_Y->push_back( bs0_ancestor->daughter(0)->vy() ) ;
		  genEvtVtx_Z->push_back( bs0_ancestor->daughter(0)->vz() ) ;
		  genEvtVtx_particles->push_back( bs0_ancestor->numberOfDaughters() ) ;
		  n_Bs0Ancestors->push_back( n_ancestors ) ; 
		}
		else {
		  bs0_ancestor = bs0_ancestor->mother(0) ;
		  n_ancestors++ ;
		}
	      }

	    MCJPsiPx->push_back(jpsiPx); MCJPsiPy->push_back(jpsiPy); MCJPsiPz->push_back(jpsiPz);
	    MCmupPx->push_back(mupPx); MCmupPy->push_back(mupPy); MCmupPz->push_back(mupPz);
	    MCmumPx->push_back(mumPx); MCmumPy->push_back(mumPy); MCmumPz->push_back(mumPz);
            MCPhiPx->push_back(phiPx); MCPhiPy->push_back(phiPy); MCPhiPz->push_back(phiPz);
            MCkpPx->push_back(kpPx); MCkpPy->push_back(kpPy); MCkpPz->push_back(kpPz);
            MCkmPx->push_back(kmPx); MCkmPy->push_back(kmPy); MCkmPz->push_back(kmPz);
	    //MCpionPx->push_back(pionPx); MCpionPy->push_back(pionPy); MCpionPz->push_back(pionPz);
	    //MCkaonPx->push_back(kaonPx); MCkaonPy->push_back(kaonPy); MCkaonPz->push_back(kaonPz);
	    //MCpionCh->push_back(pionCh) ; MCkaonCh->push_back(kaonCh) ;
	    decayChainOK = true;
	    MCPx->push_back( p.px() );
	    MCPy->push_back( p.py() );
	    MCPz->push_back( p.pz() );
	  }
	  if (Debug_) cout << "decayChainOK = " << decayChainOK << endl;
	} // if ( abs(pdgid) == MCMother  &&  dauNum == 3 )
      } // if ( !MCExclusiveDecay )

    } // for (size_t i = 0; i < genParticles->size(); ++ i)
  } // if (doMC)
  
  /// reconstruction only for events with B decaying in psi(nS)+Pi+K /// SEMRA JPsiPhi !!!
  if ( (doMC && !MCExclusiveDecay) || (doMC && (MCExclusiveDecay && decayChainOK)) || doData ) {
	
    bool isEventWithInvalidMu = false;
	
    if (Debug_) cout << "starting event with " << thePATTrackHandle->size() << " tracks, and " << thePATMuonHandle->size() << " muons" << endl;
    
    if ((thePATMuonHandle->size()) * (thePATTrackHandle->size()) > 20000) {
      cout << "Too many Muons: " << thePATMuonHandle->size() << ", and Tracks: " << thePATTrackHandle->size() << endl;
    } else //if (thePATMuonHandle->size() >= 2) { // check
      if (thePATMuonHandle->size() >= 2  && hasRequestedTrigger) {
	if (Debug_) cout <<"============================  evt: " <<evtNum <<" Accept event with 2 mu and TRIGGER ==============================================" <<endl;
	
	////////////////// filling track tree //////////////////
	for ( vector<pat::GenericParticle>::const_iterator iTr = thePATTrackHandle->begin(); iTr != thePATTrackHandle->end(); ++iTr ) { 
	  pat::GenericParticle tr = *iTr;
	  trPx->push_back(tr.px());
	  trPy->push_back(tr.py());
	  trPz->push_back(tr.pz());
	  trE->push_back(tr.energy());
	  trPhits->push_back(tr.track()->hitPattern().numberOfValidPixelHits());
	  trShits->push_back(tr.track()->hitPattern().numberOfValidStripHits());
	  trChi2->push_back(tr.track()->chi2());
	  trNDF->push_back(tr.track()->ndof());
	  trD0->push_back(tr.track()->d0());
	  trD0E->push_back(tr.track()->d0Error());
	  trCharge->push_back(tr.charge());
	  float hits = (1.0*tr.track()->found() )/ (tr.track()->found()+ tr.track()->lost() + tr.track()->trackerExpectedHitsInner().numberOfHits() + tr.track()->trackerExpectedHitsOuter().numberOfHits());
	  trfHits->push_back(hits);
	  trFirstBarrel->push_back(tr.track()->hitPattern().hasValidHitInFirstPixelBarrel());
	  trFirstEndCap->push_back(tr.track()->hitPattern().hasValidHitInFirstPixelEndcap());
	  trDzVtx->push_back(tr.track()->dz(RefVtx));
	  trDxyVtx->push_back(tr.track()->dxy(RefVtx));
	  double theo = 0., sigma = 0. ;
	  tr_nsigdedx->push_back(nsigmaofdedx(tr.track(),theo,sigma));
	  tr_dedx->push_back(getEnergyLoss(tr.track()));
	  tr_dedxMass->push_back(GetMass(tr.track()));
	  tr_theo->push_back(theo);
	  tr_sigma->push_back(sigma);
	  tr_dedx_byHits->push_back( (dEdxTrack)[tr.track()].dEdx() );
	  tr_dedxErr_byHits->push_back( (dEdxTrack)[tr.track()].dEdxError() );
	  tr_saturMeas_byHits->push_back( (dEdxTrack)[tr.track()].numberOfSaturatedMeasurements() );
	  tr_Meas_byHits->push_back( (dEdxTrack)[tr.track()].numberOfMeasurements() );
	  /// Track quality:            
	  /// loose=0, tight=1, highPurity=2, confirmed=3, goodIterative=4, looseSetWithPV=5, highPuritySetWithPV=6
	  bool ishighPurity = tr.track()->quality(reco::TrackBase::highPurity);
	  trQualityHighPurity->push_back(ishighPurity);
	  trQualityTight->push_back(tr.track()->quality(reco::TrackBase::tight));
	}
	
	/// get MuMu cands
	for ( std::vector<pat::Muon>::const_iterator Muon1 = thePATMuonHandle->begin(); Muon1 != thePATMuonHandle->end(); ++Muon1 ) {

	  /// push back all muon information
	  ++nMu;
	  const reco::Muon* rmu1 = dynamic_cast<const reco::Muon * >(Muon1->originalObject());
	  muPx->push_back(rmu1->px());
	  muPy->push_back(rmu1->py());
	  muPz->push_back(rmu1->pz());
	  muCharge->push_back(rmu1->charge());
	
	  if (rmu1->track().isNull()) { // rmu->track() returns innerTrack();
	    cout << "no track for " << std::distance(thePATMuonHandle->begin(), Muon1) << " filling defaults" << endl; 
	    /// AF
	    muD0->push_back(0);
	    muDz->push_back(0);
	    muChi2->push_back(0);
	    muNDF->push_back(-1);
	    muPhits->push_back(0);
	    muShits->push_back(0);
	    muLayersTr->push_back(0);
	    muLayersPix->push_back(0);
	    muDzVtx->push_back(0);
	    muDxyVtx->push_back(0);
	    mufHits->push_back(0);
	    muFirstBarrel->push_back(0);
	    muFirstEndCap->push_back(0);
	    muD0E->push_back(0);
	    muDzVtxErr->push_back(0);
	    muKey->push_back(0);
	    muGlChi2->push_back(0); 
	    muGlNDF->push_back(-1);
	    muGlMuHits->push_back(0);
	    muGlMatchedStation->push_back(0);
	    muGlDzVtx->push_back(0);
	    muGlDxyVtx->push_back(0); 
	    nMatchedStations->push_back(0) ;

	    if (Debug_) cout <<"evt:" <<evtNum << "no track for PAT muon " <<std::distance(thePATMuonHandle->begin(), Muon1) <<" skipping muon... should skip event instead" <<endl;
	    isEventWithInvalidMu = true;
	    continue;
	  }
	  else {
	    muD0->push_back(rmu1->track()->d0());
	    muDz->push_back(rmu1->track()->dz());
	    muChi2->push_back(rmu1->track()->chi2());
	    muNDF->push_back(rmu1->track()->ndof());
	    muPhits->push_back(rmu1->track()->hitPattern().numberOfValidPixelHits());
	    muShits->push_back(rmu1->track()->hitPattern().numberOfValidStripHits());
	    if (Debug_) cout <<"evt:" <<evtNum <<" trackerLayersWithMeasurement=" <<rmu1->track()->hitPattern().trackerLayersWithMeasurement() <<endl;
	    if ( !(rmu1->track()->hitPattern().trackerLayersWithMeasurement()) ) { 
	      isEventWithInvalidMu = true;
	      if (Debug_) cout <<"evt:" <<evtNum <<" problem with trackerLayersWithMeasurement" <<endl;
	      continue ;
	    }
	    if ( !(rmu1->track()->hitPattern().pixelLayersWithMeasurement()) ) {
	      isEventWithInvalidMu = true;
	      continue ;
	    }
	    muLayersTr->push_back(rmu1->track()->hitPattern().trackerLayersWithMeasurement());
	    muLayersPix->push_back(rmu1->track()->hitPattern().pixelLayersWithMeasurement());	
	    muDzVtx->push_back(rmu1->track()->dz(RefVtx));
	    muDxyVtx->push_back(rmu1->track()->dxy(RefVtx));
	    mufHits->push_back((1.0*rmu1->track()->found())/ (rmu1->track()->found()+ rmu1->track()->lost() + rmu1->track()->trackerExpectedHitsInner().numberOfHits() + rmu1->track()->trackerExpectedHitsOuter().numberOfHits() ) );
	    if (Debug_) cout <<"mu found " <<rmu1->track()->found() <<" fHits=" <<(1.0*rmu1->track()->found())/ (rmu1->track()->found()+ rmu1->track()->lost() + rmu1->track()->trackerExpectedHitsInner().numberOfHits() + rmu1->track()->trackerExpectedHitsOuter().numberOfHits() ) <<endl;
	    muFirstBarrel->push_back(rmu1->track()->hitPattern().hasValidHitInFirstPixelBarrel());
	    muFirstEndCap->push_back(rmu1->track()->hitPattern().hasValidHitInFirstPixelEndcap());
	    muD0E->push_back(rmu1->track()->d0Error());
	    muDzVtxErr->push_back(rmu1->track()->dzError());
	    muKey->push_back(rmu1->track().key());
	  }
	  muIsGlobal->push_back( rmu1->isGlobalMuon() ) ;
	  muIsPF->push_back( rmu1->isPFMuon() ) ;
	  if ( rmu1->globalTrack().isNull() ) { 
	    muGlMuHits->push_back(0);
	    muGlChi2->push_back(0);
	    muGlNDF->push_back(-1);
	    muGlMatchedStation->push_back(0);
	    muGlDzVtx->push_back(-1);
	    muGlDxyVtx->push_back(-1);
	  }
	  else {
	    muGlMuHits->push_back(rmu1->globalTrack()->hitPattern().numberOfValidMuonHits());
	    muGlChi2->push_back(rmu1->globalTrack()->chi2());
	    muGlNDF->push_back(rmu1->globalTrack()->ndof());
	    muGlMatchedStation->push_back(rmu1->numberOfMatchedStations());
	    muGlDzVtx->push_back(rmu1->globalTrack()->dz(RefVtx));
	    muGlDxyVtx->push_back(rmu1->globalTrack()->dxy(RefVtx));
	  }
	  nMatchedStations->push_back(rmu1->numberOfMatchedStations()) ;
	  muType->push_back(rmu1->type());
	  int qm = 0;
	  for (int qi=1; qi!= 24; ++qi) {
	    if (muon::isGoodMuon(*rmu1, muon::SelectionType(qi)))
	      qm += 1<<qi;
	  }
	  muQual->push_back(qm);
	  muTrack->push_back(-1);// not implemented yet
	
	  ////////////////// muon cleaning //////////////////
	  int nOverlapMus = 0, nSharingSegWith = -1;
	  int nSegments1 = rmu1->numberOfMatches(reco::Muon::SegmentArbitration);
	  for ( std::vector<pat::Muon>::const_iterator Muon2 = Muon1+1; Muon2 != thePATMuonHandle->end(); ++Muon2) {
	    const reco::Muon* rmu2 = dynamic_cast<const reco::Muon*>(Muon2->originalObject());
	    if ( isSameMuon(*rmu1, *rmu2)) continue;
	    if ( !muon::isGoodMuon(*rmu2, muon::TMOneStationTight) ) continue;
	    /// geometric overlap
	    if ( muon::overlap( *rmu1, *rmu2 ) )
	      nOverlapMus++ ;
	    /// shared segments 
	    int nSegments2 = rmu2->numberOfMatches(reco::Muon::SegmentArbitration);
	      if (nSegments2 == 0 || nSegments1 == 0) continue;
	      double sf = muon::sharedSegments(*rmu1, *rmu2) / std::min<double>(nSegments1, nSegments2);
	      if (sf > sharedFraction) {
		nSharingSegWith = 0;
		if ( !isBetterMuon(*rmu1, *rmu2) ) 
		  nSharingSegWith++ ;
	      }
	  }
	  muNOverlap->push_back( nOverlapMus ) ;
	  muNSharingSegWith->push_back( nSharingSegWith ) ;


	  ////////////////// check for muon1 //////////////////
	  TrackRef muTrack1 = Muon1->track();
	  if ( muTrack1.isNull() )
	    continue;      
	  /// cuts on muon1
	  if (rmu1->track()->hitPattern().numberOfValidPixelHits() < MuMinPixHits
	      || rmu1->track()->hitPattern().numberOfValidStripHits() < MuMinSiHits
	      || rmu1->track()->chi2()/rmu1->track()->ndof() > MuMaxNormChi
	      || fabs(rmu1->track()->dxy(RefVtx)) > MuMaxD0) {
	    continue ;
	  }
			
	  ////////////////// check for muon2 //////////////////
	  for ( std::vector<pat::Muon>::const_iterator Muon2 = Muon1+1; Muon2 != thePATMuonHandle->end(); ++Muon2) {
	    if(Muon2->charge() * Muon1->charge() > 0)
	      continue ; 
	    const reco::Muon* rmu2 = dynamic_cast<const reco::Muon *>(Muon2->originalObject()) ;	 
	    if (muon::overlap(*rmu1, *rmu2) )
	      continue ;	
	    TrackRef muTrack2 = Muon2->track() ;
	    if ( muTrack2.isNull() )
	      continue ;	
	    /// cuts on muon2
	    if (rmu2->track()->hitPattern().numberOfValidPixelHits() < MuMinPixHits
		|| rmu2->track()->hitPattern().numberOfValidStripHits() < MuMinSiHits
		|| rmu2->track()->chi2()/rmu1->track()->ndof() > MuMaxNormChi
		|| fabs(rmu2->track()->dxy(RefVtx)) > MuMaxD0) {
	      continue ;
	    }


	    ////////////////// get the MuMu information //////////////////
	    TransientTrack muon1TT( muTrack1, &(*bFieldHandle) );
	    TransientTrack muon2TT( muTrack2, &(*bFieldHandle) );			
	    KinematicParticleFactoryFromTransientTrack pFactory;
			
	    /// initial chi2 and ndf before kinematic fits
	    float chi = 0., ndf = 0.;
	    vector<RefCountedKinematicParticle> muons; /// the final state muons produced by the KinematicParticleFactory 
	    muons.push_back( pFactory.particle( muon1TT, muon_mass, chi, ndf, small_sigma));
	    muons.push_back( pFactory.particle( muon2TT, muon_mass, chi, ndf, small_sigma));
	    KinematicParticleVertexFitter MuMuFitter; /// creating the vertex fitter for JPsi  
	    RefCountedKinematicTree MuMuVertexFitTree;
	    MuMuVertexFitTree = MuMuFitter.fit(muons); 
	    if (!MuMuVertexFitTree->isValid())  
	      continue ; 
	    MuMuVertexFitTree->movePointerToTheTop(); 				
	    RefCountedKinematicParticle MuMuCand_fromFit = MuMuVertexFitTree->currentParticle();
	    RefCountedKinematicVertex MuMuCand_vertex_fromFit = MuMuVertexFitTree->currentDecayVertex();
	    MuMuVertexFitTree->movePointerToTheFirstChild();
	    RefCountedKinematicParticle Mu1Cand_fromFit = MuMuVertexFitTree->currentParticle(); 
	    MuMuVertexFitTree->movePointerToTheNextChild();
	    RefCountedKinematicParticle Mu2Cand_fromFit = MuMuVertexFitTree->currentParticle();
	    KinematicParameters Mu1Cand_KP = Mu1Cand_fromFit->currentState().kinematicParameters();
	    KinematicParameters Mu2Cand_KP = Mu2Cand_fromFit->currentState().kinematicParameters();
	    
            ////////////////// fill the MuMu vectors //////////////////
	    if (MuMuCand_fromFit->currentState().mass() < MuMuMinMass  ||  MuMuCand_fromFit->currentState().mass() > MuMuMaxMass) 
	      continue ;
	    MuMuMass->push_back( MuMuCand_fromFit->currentState().mass() );	
	    MuMuDecayVtx_X->push_back( MuMuCand_vertex_fromFit->position().x() );
	    MuMuDecayVtx_Y->push_back( MuMuCand_vertex_fromFit->position().y() );
	    MuMuDecayVtx_Z->push_back( MuMuCand_vertex_fromFit->position().z() );
	    MuMuDecayVtx_XE->push_back( sqrt( MuMuCand_vertex_fromFit->error().cxx()) );
	    MuMuDecayVtx_YE->push_back( sqrt( MuMuCand_vertex_fromFit->error().cyy()) );
	    MuMuDecayVtx_ZE->push_back( sqrt( MuMuCand_vertex_fromFit->error().czz()) );
	    MuMuVtx_CL->push_back( ChiSquaredProbability((double)( MuMuCand_vertex_fromFit->chiSquared()),(double)( MuMuCand_vertex_fromFit->degreesOfFreedom())) );
	    MuMuVtx_Chi2->push_back( MuMuCand_vertex_fromFit->chiSquared() ) ;
	    MuMuPx->push_back( Mu1Cand_KP.momentum().x() + Mu2Cand_KP.momentum().x() );
	    MuMuPy->push_back( Mu1Cand_KP.momentum().y() + Mu2Cand_KP.momentum().y() );
	    MuMuPz->push_back( Mu1Cand_KP.momentum().z() + Mu2Cand_KP.momentum().z() );
	    mu1Idx->push_back(std::distance(thePATMuonHandle->begin(), Muon1)); 
	    mu2Idx->push_back(std::distance(thePATMuonHandle->begin(), Muon2));

            ////////////////// JPsi (MuMu) fit //////////////////
            mu1_MuMu_Px->push_back( Mu1Cand_KP.momentum().x()); /// SEMRA for JPsi
	    mu1_MuMu_Py->push_back( Mu1Cand_KP.momentum().y());
	    mu1_MuMu_Pz->push_back( Mu1Cand_KP.momentum().z());
	    mu1_MuMu_Chi2->push_back( Mu1Cand_fromFit->chiSquared());
	    mu1_MuMu_NDF->push_back( Mu1Cand_fromFit->degreesOfFreedom());
	    mu2_MuMu_Px->push_back( Mu2Cand_KP.momentum().x());
	    mu2_MuMu_Py->push_back( Mu2Cand_KP.momentum().y());
	    mu2_MuMu_Pz->push_back( Mu2Cand_KP.momentum().z());
	    mu2_MuMu_Chi2->push_back( Mu2Cand_fromFit->chiSquared());
	    mu2_MuMu_NDF->push_back( Mu2Cand_fromFit->degreesOfFreedom());

	    Int_t dimuonType = 0;   //0 nothing,  1 J/psi  , 2 psi(2S)   
	    if ( MuMuCand_fromFit->currentState().mass() > JPsiMinMass  &&  MuMuCand_fromFit->currentState().mass() < JPsiMaxMass ) {
	      dimuonType = 1 ;
	    } 
	    if (Debug_) cout <<dimuonType <<endl;
     
	    if (Debug_) cout <<"evt:" <<evtNum <<" MuMu with diMuonType = " <<dimuonType <<endl;
	    //cout << "POINT 0" << endl; 
	    MuMuType->push_back(dimuonType);
	    //cout << "POINT 2" << endl;	

	    int ntriggers = TriggersForMatching_.size();
	    if (Debug_) cout << "ntriggers: " << ntriggers << endl; 
	    for (int MatchTrig = 0; MatchTrig < ntriggers; MatchTrig++) 
	      {
	        if (Debug_) cout << "MatchingTriggerResult[" << MatchTrig << "]: " << MatchingTriggerResult[MatchTrig] << endl; 
		if ( MatchingTriggerResult[MatchTrig]!=0 ) 
		  {
		    //cout << "POINT 3" << endl; 
		    if (Debug_) cout << "CHECKING FiltersForMatching_[" << MatchTrig << "]: " << FiltersForMatching_[MatchTrig] << endl; 
		    //cout << "POINT 4" << endl; 
		    pat::TriggerObjectStandAloneCollection mu1HLTMatches = Muon1->triggerObjectMatchesByFilter( FiltersForMatching_[MatchTrig] );
		    //cout << "POINT 5" << endl;
		    pat::TriggerObjectStandAloneCollection mu2HLTMatches = Muon2->triggerObjectMatchesByFilter( FiltersForMatching_[MatchTrig] );
	            //cout << "POINT 6" << endl; 
		    bool pass1 = mu1HLTMatches.size() > 0;
		    bool pass2 = mu2HLTMatches.size() > 0;
		    //cout << "POINT 7" << endl; 	  
		    if ((pass1) && (pass2))
		      {
			//cout << "POINT 8" << endl; 
			MuMuMuonTrigMatch->push_back(true);
			if (Debug_) cout <<"Matched MuMu" <<endl ;
		      } else
	              //cout << "POINT 9" << endl; 
		      MuMuMuonTrigMatch->push_back(false);
		  }
		else
		  //cout << "POINT 10" << endl; 
		  MuMuMuonTrigMatch->push_back(false);
	      }
  
	    /// vertex without matched muons 
	    vector<TransientVertex> pvs ;
	    Vertex MuMuLessPV = thePrimaryVtx ;

	    if (addMuMulessPrimaryVertex_)
	      {
		VertexReProducer revertex(recVtxs, iEvent);
		Handle<TrackCollection> pvtracks;   
		iEvent.getByLabel(revertex.inputTracks(),   pvtracks);
		Handle<BeamSpot>        pvbeamspot;
		iEvent.getByLabel(revertex.inputBeamSpot(), pvbeamspot);
		
		if ( pvbeamspot.isValid() < 0 ) 
		  continue ; 
		if (pvbeamspot.id() != beamSpotHandle.id()) {
		  edm::LogWarning("Inconsistency") << "The BeamSpot used for PV reco is not the same used in this analyzer.";
		}
		const reco::Muon *rmu_1 = dynamic_cast<const reco::Muon*>( Muon1->originalObject() ) ;
		const reco::Muon *rmu_2 = dynamic_cast<const reco::Muon*>( Muon2->originalObject() ) ;
		
		if (rmu_1 != 0  &&  rmu_2 != 0  &&  rmu_1->track().id() == pvtracks.id()  &&  rmu_2->track().id() == pvtracks.id() ) { 
		  TrackCollection MuMuLess;
		  MuMuLess.reserve(pvtracks->size());
		  for (size_t i = 0, n = pvtracks->size(); i < n; ++i) {
		    if (i == rmu_1->track().key()) continue;
		    if (i == rmu_2->track().key()) continue; 
		    MuMuLess.push_back((*pvtracks)[i]);
		  }
		  if (Debug_) cout <<"pvbeamspot.isValid() = " <<pvbeamspot.isValid() <<endl ;
		  pvs = revertex.makeVertices(MuMuLess, *pvbeamspot, iSetup) ;
		  if (!pvs.empty()) {
		    MuMuLessPV = Vertex(pvs.front());
		  }
		}
	      }							
	    PriVtxMuMuCorr_n->push_back( pvs.size() ) ;
	    PriVtxMuMuCorr_X->push_back( MuMuLessPV.position().x() ) ;
	    PriVtxMuMuCorr_Y->push_back( MuMuLessPV.position().y() ) ; 
	    PriVtxMuMuCorr_Z->push_back( MuMuLessPV.position().z() ) ;
	    PriVtxMuMuCorr_EX->push_back( MuMuLessPV.xError() ) ;
	    PriVtxMuMuCorr_EY->push_back( MuMuLessPV.yError() ) ;
	    PriVtxMuMuCorr_EZ->push_back( MuMuLessPV.zError() ) ;
	    PriVtxMuMuCorr_CL->push_back( ChiSquaredProbability( (double)(MuMuLessPV.chi2()), (double)(MuMuLessPV.ndof())) ) ;
	    PriVtxMuMuCorr_Chi2->push_back( MuMuLessPV.chi2() ) ;
	    PriVtxMuMuCorr_tracks->push_back( MuMuLessPV.tracksSize() ) ;
		
            ++nMuMu;
	    muons.clear();

	    //////////////////////////////////////////////////////////////////////
            /// for Bs0 & X(4140)
	    if (Debug_) cout <<"evt:"<<evtNum<< " is Invalid Muon ?  " <<isEventWithInvalidMu << endl;
	    //if (skipJPsi && ( dimuonType == 1 ));
            //cout<< "POINT 11" <<endl;
	    nTrk->push_back( thePATTrackHandle->size() ) ;
	    //cout<< "POINT 12" <<endl;
            if (thePATTrackHandle->size() < 2) {
            //cout<< "POINT 13" <<endl;
            continue ; nBs0_pre0++; nX_pre0++ ; 
	    }  
            if (Debug_) cout<<"nmumu : "<<nMuMu<<endl;


	    ////////////////// cuts on MuMu mass window for Bs0 ////////////////////////////
	    if (MuMuMass->at(nMuMu-1) < MuMuMinMass  ||  MuMuMass->at(nMuMu-1) > MuMuMaxMass){
	       nBs0_pre1++ ; nX_pre1++ ;
            }
	    //cout<< "POINT 14" <<endl;


	    ////////////////// check tracks for kaon1 for Bs0 & X(4140)//////////////////
	    for ( vector<pat::GenericParticle>::const_iterator Track1 = theKaonRefittedPATTrackHandle->begin(); Track1 != theKaonRefittedPATTrackHandle->end(); ++Track1 ) {
	      //cout<< "POINT 15" <<endl;
	      /// check track doesn't overlap with the MuMu candidate tracks
	      if (Track1->track().key() == rmu1->track().key()  ||  Track1->track().key() == rmu2->track().key())
		 continue ; nBs0_pre2++ ; nX_pre2++ ;
      
	      //cout<< "POINT 16" <<endl;
	      /// cuts on charged tracks	
	      if (( Track1->track()->chi2()/Track1->track()->ndof() > TrMaxNormChi2 )  ||  Track1->pt() < TrMinPt)
		continue ; nBs0_pre3++ ; nX_pre3++ ;
              
	      //cout<< "POINT 17" <<endl;

	   ////////////////// check tracks for kaon2 for Bs0 & X(4140) //////////////////
	   for ( vector<pat::GenericParticle>::const_iterator Track2 = Track1+1; Track2 != theKaonRefittedPATTrackHandle->end(); ++Track2 ){
	     /// check that this second track doesn't overlap with the the first track candidate
	     if (Track2->track().key() == Track1->track().key()) 
	        continue ; nBs0_pre4++ ; nX_pre4++ ;
            
            /// check track doesn't overlap with the MuMu candidate tracks
            if (Track2->track().key() == rmu1->track().key()  ||  Track2->track().key() == rmu2->track().key())
              continue ; nBs0_pre5++ ; nX_pre5++ ;
            if (Track1->charge() * Track2->charge() > 0)
	      continue ; nBs0_pre6++ ; nX_pre6++ ; 
	    /// cuts on charged tracks	
            if ((Track2->track()->chi2() / Track2->track()->ndof() > TrMaxNormChi2)  ||  Track2->pt() < TrMinPt)
	      continue; nBs0_pre7++ ; nX_pre7++ ;
		
		
	   ////////////////// get the KK information //////////////////
	   TransientTrack kaon1TT( Track1->track(), &(*bFieldHandle) );
           TransientTrack kaon2TT( Track2->track(), &(*bFieldHandle) );
           KinematicParticleFactoryFromTransientTrack pFactory;

           /// initial chi2 and ndf before kinematic fits
	   float chi = 0., ndf = 0.;
           vector<RefCountedKinematicParticle> kaons;
           kaons.push_back( pFactory.particle( kaon1TT, kaon_mass, chi, ndf, small_sigma));
           kaons.push_back( pFactory.particle( kaon2TT, kaon_mass, chi, ndf, small_sigma));
           KinematicParticleVertexFitter KKFitter;
           RefCountedKinematicTree KKVertexFitTree;
           KKVertexFitTree = KKFitter.fit(kaons);
            if (!KKVertexFitTree->isValid())
              continue ;
           KKVertexFitTree->movePointerToTheTop();
           RefCountedKinematicParticle KKCand_fromFit = KKVertexFitTree->currentParticle();
           RefCountedKinematicVertex KKCand_vertex_fromFit = KKVertexFitTree->currentDecayVertex();
           KKVertexFitTree->movePointerToTheFirstChild();
           RefCountedKinematicParticle Ka1Cand_fromFit = KKVertexFitTree->currentParticle();
           KKVertexFitTree->movePointerToTheNextChild();
           RefCountedKinematicParticle Ka2Cand_fromFit = KKVertexFitTree->currentParticle();
           KinematicParameters Ka1Cand_KP = Ka1Cand_fromFit->currentState().kinematicParameters();
           KinematicParameters Ka2Cand_KP = Ka2Cand_fromFit->currentState().kinematicParameters();

           ////////////////// fill the KK vectors //////////////////   
	   if (KKCand_fromFit->currentState().mass() < KKMinMass  ||  KKCand_fromFit->currentState().mass() > KKMaxMass)
              continue ;
           KKMass->push_back( KKCand_fromFit->currentState().mass() );
           KKDecayVtx_X->push_back( KKCand_vertex_fromFit->position().x() );
           KKDecayVtx_Y->push_back( KKCand_vertex_fromFit->position().y() );
           KKDecayVtx_Z->push_back( KKCand_vertex_fromFit->position().z() );
           KKDecayVtx_XE->push_back( sqrt( KKCand_vertex_fromFit->error().cxx()) );
           KKDecayVtx_YE->push_back( sqrt( KKCand_vertex_fromFit->error().cyy()) );
           KKDecayVtx_ZE->push_back( sqrt( KKCand_vertex_fromFit->error().czz()) );
           KKVtx_CL->push_back( ChiSquaredProbability((double)( KKCand_vertex_fromFit->chiSquared()),(double)( KKCand_vertex_fromFit->degreesOfFreedom())) );
           KKVtx_Chi2->push_back( MuMuCand_vertex_fromFit->chiSquared() ) ;
           KKPx->push_back( Ka1Cand_KP.momentum().x() + Ka2Cand_KP.momentum().x() );
           KKPy->push_back( Ka1Cand_KP.momentum().y() + Ka2Cand_KP.momentum().y() );
           KKPz->push_back( Ka1Cand_KP.momentum().z() + Ka2Cand_KP.momentum().z() );
           ka1Idx->push_back(std::distance(theKaonRefittedPATTrackHandle->begin(), Track1));
           ka2Idx->push_back(std::distance(theKaonRefittedPATTrackHandle->begin(), Track2));

           ////////////////// Phi (KK) fit //////////////////
	   ka1_KK_Px->push_back( Ka1Cand_KP.momentum().x()); 
           ka1_KK_Py->push_back( Ka1Cand_KP.momentum().y());
           ka1_KK_Pz->push_back( Ka1Cand_KP.momentum().z());
           ka1_KK_Chi2->push_back( Ka1Cand_fromFit->chiSquared());
           ka1_KK_NDF->push_back( Ka1Cand_fromFit->degreesOfFreedom());
           ka2_KK_Px->push_back( Ka2Cand_KP.momentum().x());
           ka2_KK_Py->push_back( Ka2Cand_KP.momentum().y());
           ka2_KK_Pz->push_back( Ka2Cand_KP.momentum().z());
           ka2_KK_Chi2->push_back( Ka2Cand_fromFit->chiSquared());
           ka2_KK_NDF->push_back( Ka2Cand_fromFit->degreesOfFreedom());	

	   ++nKK;
           kaons.clear();

            ////////////////// cuts on tracks' delta R for Bs0 & X(4140) //////////////////	
            math::XYZTLorentzVector MuMu = (rmu1->p4() + rmu2->p4()); 
	    math::XYZTLorentzVector bs0 = (rmu1->p4() + rmu2->p4() + Track1->p4() + Track2->p4());   
	    math::XYZTLorentzVector x = (rmu1->p4() + rmu2->p4() + Track1->p4() + Track2->p4());
	    float MuMuK1DR = sqrt( pow(MuMu.eta() - Track1->p4().eta(),2) + pow(MuMu.phi() - Track1->p4().phi(), 2) );
	    float MuMuK2DR = sqrt( pow(MuMu.eta() - Track2->p4().eta(),2) + pow(MuMu.phi() - Track2->p4().phi(), 2) );
            float bs0K1DR = sqrt( pow(bs0.eta() - Track1->p4().eta(),2) + pow(bs0.phi() - Track1->p4().phi(), 2) );
            float bs0K2DR = sqrt( pow(MuMu.eta() - Track2->p4().eta(),2) + pow(MuMu.phi() - Track2->p4().phi(), 2) );
            float xK1DR = sqrt( pow(x.eta() - Track1->p4().eta(),2) + pow(x.phi() - Track1->p4().phi(), 2) );
            float xK2DR = sqrt( pow(MuMu.eta() - Track2->p4().eta(),2) + pow(MuMu.phi() - Track2->p4().phi(), 2) );		 

	    DRMuMuK1->push_back(MuMuK1DR);
	    DRMuMuK2->push_back(MuMuK2DR);
            DRbs0K1->push_back(bs0K1DR);
            DRbs0K2->push_back(bs0K2DR);
            DRxK1->push_back(xK1DR);
            DRxK2->push_back(xK2DR);
	    

		  if (UseBs0DR) {
		    if (bs0K1DR > Bs0TrackMaxDR || bs0K2DR > Bs0TrackMaxDR)
		      Bs0TrackMaxDR = 2; 
		  } else {
		    if (MuMuK1DR > MuMuTrackMaxDR || MuMuK2DR > MuMuTrackMaxDR)
		      MuMuTrackMaxDR = 3.5;
		  }
		  nBs0_pre8++ ; 
 
		  if (UseXDR) {
		    if (xK1DR > XTrackMaxDR || xK2DR > XTrackMaxDR)
		   XTrackMaxDR = 2;   
                  } else {
		    if (MuMuK1DR > MuMuTrackMaxDR || MuMuK2DR > MuMuTrackMaxDR)
		     MuMuTrackMaxDR = 3.5;
		  }
		  nX_pre8++ ;

            ////////////////// cuts on MuMuKK mass window for Bs0 & X(4140) //////////////////
            if ((Track1->p4() + Track2->p4() + MuMu).M() < MuMuKKMaxBs0Mass  ||  (Track1->p4() + Track2->p4() + MuMu).M() > MuMuKKMinBs0Mass)
	      nBs0_pre9++ ; 
	    if ((Track1->p4() + Track2->p4() + MuMu).M() > MuMuKKMaxXMass) /// < changed with >
	      nX_pre9++ ;   
	    //cout << "mass window" << endl;



                  /// having two oppositely charged muons, and two oppositely charged tracks: try to vertex them
                 //TransientTrack kaon1TT( Track1->track(), &(*bFieldHandle) ); 
                 //TransientTrack kaon2TT( Track2->track(), &(*bFieldHandle) );
		
		  TransientTrack kaon2TT_notRefit ; 
		  Bool_t notRefittedPartner = false ;
		  for ( vector<pat::GenericParticle>::const_iterator Track2_notRefit = thePATTrackHandle->begin(); Track2_notRefit != thePATTrackHandle->end(); ++Track2_notRefit )
		    if ( Track2_notRefit->track().key() == Track2->track().key() ) {
		      notRefittedPartner = true ;
		      kaon2TT_notRefit = TransientTrack( Track2_notRefit->track(), &(*bFieldHandle) ) ; 
		      break ;
		    }
			      
		  /// do mass constraint for MuMu cand and do mass constrained vertex fit for Bs0
		  vector<RefCountedKinematicParticle> bs0Daughters;
		  bs0Daughters.push_back(pFactory.particle( muon1TT, muon_mass, chi, ndf, small_sigma));
		  bs0Daughters.push_back(pFactory.particle( muon2TT, muon_mass, chi, ndf, small_sigma));
		  bs0Daughters.push_back(pFactory.particle( kaon1TT, kaon_mass, chi, ndf, small_sigma)); 
		  bs0Daughters.push_back(pFactory.particle( kaon2TT, kaon_mass, chi, ndf, small_sigma)); 
		
		  RefCountedKinematicTree Bs0VertexFitTree, Bs0VertexFitTree_noKrefit ;
		  KinematicConstrainedVertexFitter Bs0Fitter ; 
		
		  if (doMuMuMassConst) { // MassConst = 'MC' in the following
		    MultiTrackKinematicConstraint *MuMu = 0;
		    if (dimuonType == 1) { // constrain to JPsi mass
		      MuMu = new TwoTrackMassKinematicConstraint(JPsi_mass);
		    } //else if (dimuonType == 2) { // constrain to Psi(2S) mass /// SEMRA will we use this or not ?
		      //MuMu = new TwoTrackMassKinematicConstraint(psi2S_mass);
		    //} // already asked for: if (dimuonType == 0) continue ;
		  
		    Bs0VertexFitTree = Bs0Fitter.fit( bs0Daughters, MuMu ); 
		    if (notRefittedPartner) { // use not refitted kaons
		      bs0Daughters.pop_back() ;
		      bs0Daughters.push_back(pFactory.particle( kaon2TT_notRefit, kaon_mass, chi, ndf, small_sigma)); 
		      Bs0VertexFitTree_noKrefit = Bs0Fitter.fit( bs0Daughters, MuMu ); 
		    }
		  } 
		  else {
		    Bs0VertexFitTree = Bs0Fitter.fit( bs0Daughters ); 
		    if (notRefittedPartner) { // use not refitted kaons
		      bs0Daughters.pop_back() ;
		      bs0Daughters.push_back(pFactory.particle( kaon2TT_notRefit, kaon_mass, chi, ndf, small_sigma)); 
		      Bs0VertexFitTree_noKrefit = Bs0Fitter.fit( bs0Daughters ); 
		    }
		  }
		  
		  if ( !Bs0VertexFitTree->isValid() ) 
	             continue ; nBs0_pre10++ ;     
 
		  Bs0VertexFitTree->movePointerToTheTop(); 
		  RefCountedKinematicParticle Bs0Cand_fromMCFit = Bs0VertexFitTree->currentParticle(); 
		  RefCountedKinematicVertex Bs0Cand_vertex_fromMCFit = Bs0VertexFitTree->currentDecayVertex(); 
		  		  
		  if ( !Bs0Cand_vertex_fromMCFit->vertexIsValid() )  
	             continue ; nBs0_pre11++ ; 
		  		  
		  if ( Bs0Cand_vertex_fromMCFit->chiSquared() < 0  ||  Bs0Cand_vertex_fromMCFit->chiSquared() > 10000 ) nBs0_pre12++ ; 

	          if (Bs0Cand_vertex_fromMCFit->chiSquared() / Bs0Cand_vertex_fromMCFit->degreesOfFreedom() > 7 ) nBs0_pre13++; 
 
		  if ( Bs0Cand_fromMCFit->currentState().mass() > 100 ) nBs0_pre14++ ; 
		  
		  double bs0VtxProb = ChiSquaredProbability((double)(Bs0Cand_vertex_fromMCFit->chiSquared()), (double)(Bs0Cand_vertex_fromMCFit->degreesOfFreedom())); 
		  if ( bs0VtxProb < 0.005 ) //0.0001 ) 
		     nBs0_pre15++ ;
					
		 /// do mass constraint for MuMu cand and do mass constrained vertex fit for X(4140)
		 vector<RefCountedKinematicParticle> xDaughters;
                 xDaughters.push_back(pFactory.particle( muon1TT, muon_mass, chi, ndf, small_sigma));
                 xDaughters.push_back(pFactory.particle( muon2TT, muon_mass, chi, ndf, small_sigma));
                 xDaughters.push_back(pFactory.particle( kaon1TT, kaon_mass, chi, ndf, small_sigma)); 
                 xDaughters.push_back(pFactory.particle( kaon2TT, kaon_mass, chi, ndf, small_sigma)); 

                 RefCountedKinematicTree XVertexFitTree, XVertexFitTree_noKrefit ;  
                 KinematicConstrainedVertexFitter XFitter ;
		 if (doMuMuMassConst) { // MassConst = 'MC' in the following
                    MultiTrackKinematicConstraint *MuMu = 0;
                    if (dimuonType == 1) { // constrain to JPsi mass
                      MuMu = new TwoTrackMassKinematicConstraint(JPsi_mass);
                    } //else if (dimuonType == 2) { // constrain to Psi(2S) mass
                     // MuMu = new TwoTrackMassKinematicConstraint(psi2S_mass);
                    //} 
                    XVertexFitTree = XFitter.fit( xDaughters, MuMu ); 
                  if (notRefittedPartner) { // use not refitted kaons
                      xDaughters.pop_back() ;
                      xDaughters.push_back(pFactory.particle( kaon2TT_notRefit, kaon_mass, chi, ndf, small_sigma)); 
                      XVertexFitTree_noKrefit = XFitter.fit( bs0Daughters, MuMu ); 
                    }
                  }
                  else {
                    XVertexFitTree = XFitter.fit( xDaughters );  
                    if (notRefittedPartner) { // use not refitted kaons
                      xDaughters.pop_back() ;
                      xDaughters.push_back(pFactory.particle( kaon2TT_notRefit, kaon_mass, chi, ndf, small_sigma)); 
                      XVertexFitTree_noKrefit = XFitter.fit( xDaughters ); 
                    }
                  }
                  if ( !XVertexFitTree->isValid() ) 
		     continue ; nX_pre10++ ;        
                  
		  XVertexFitTree->movePointerToTheTop(); 
		  RefCountedKinematicParticle XCand_fromMCFit = XVertexFitTree->currentParticle();  
                  RefCountedKinematicVertex XCand_vertex_fromMCFit = XVertexFitTree->currentDecayVertex();
  
                  if ( !XCand_vertex_fromMCFit->vertexIsValid() ) 
	             continue ; nX_pre11++ ;
 
                  if ( XCand_vertex_fromMCFit->chiSquared() < 0  ||  XCand_vertex_fromMCFit->chiSquared() > 10000 ) nX_pre12++ ;

		 if (XCand_vertex_fromMCFit->chiSquared() / XCand_vertex_fromMCFit->degreesOfFreedom() > 7) nX_pre13++;		
                  //cout << "chi2 cut" << endl;
 
                  if ( XCand_fromMCFit->currentState().mass() > 100 ) nX_pre14++ ; 
                  double xVtxProb = ChiSquaredProbability((double)(XCand_vertex_fromMCFit->chiSquared()), (double)(XCand_vertex_fromMCFit->degreesOfFreedom())); 
 
                  if ( xVtxProb < 0.005 ) nX_pre15++ ; 

											
		  //////////////////// Lifetimes calculations for Bs0 & X(4140) //////////////////// 
		  TVector3 Bs0_vtx((*Bs0Cand_vertex_fromMCFit).position().x(), (*Bs0Cand_vertex_fromMCFit).position().y(), 0) ; /// Bs0		
		  TVector3 Bs0_pperp(Bs0Cand_fromMCFit->currentState().globalMomentum().x(), Bs0Cand_fromMCFit->currentState().globalMomentum().y(), 0); 
		  TVector3 Bs0_vtx3D((*Bs0Cand_vertex_fromMCFit).position().x(), (*Bs0Cand_vertex_fromMCFit).position().y(), (*Bs0Cand_vertex_fromMCFit).position().z()) ;
		  TVector3 Bs0_pperp3D(Bs0Cand_fromMCFit->currentState().globalMomentum().x(), Bs0Cand_fromMCFit->currentState().globalMomentum().y(), Bs0Cand_fromMCFit->currentState().globalMomentum().z()); 
		  AlgebraicVector3 Bs0_v3pperp ; 
		  Bs0_v3pperp[0] = Bs0_pperp.x(); Bs0_v3pperp[1] = Bs0_pperp.y(); Bs0_v3pperp[2] = 0.; 
		  TVector3 Bs0_pvtx, Bs0_pvtx3D, Bs0_vdiff, Bs0_vdiff3D ; 
		  double Bs0_cosAlpha, Bs0_cosAlpha3D, Bs0_ctau ; 
		  VertexDistanceXY Bs0_vdistXY ; 
		  Measurement1D Bs0_distXY ; 
		  GlobalError Bs0_v1e = (Vertex(*Bs0Cand_vertex_fromMCFit)).error(); 
		  GlobalError Bs0_v2e ; 
		  AlgebraicSymMatrix33 Bs0_vXYe ;
		  double Bs0_ctauErr ; 
		  float Bs0_lxy, Bs0_lxyErr, Bs0_lxyz, Bs0_lxyzErr ; 
		  ROOT::Math::SVector<double, 3> Bs0_vDiff, Bs0_vDiff3D ; // needed by Similarity method 
		  
		  TVector3 X_vtx((*XCand_vertex_fromMCFit).position().x(), (*XCand_vertex_fromMCFit).position().y(), 0) ; /// X(4140)
                  TVector3 X_pperp(XCand_fromMCFit->currentState().globalMomentum().x(), XCand_fromMCFit->currentState().globalMomentum().y(), 0);
                  TVector3 X_vtx3D((*XCand_vertex_fromMCFit).position().x(), (*XCand_vertex_fromMCFit).position().y(), (*XCand_vertex_fromMCFit).position().z()) ;
                  TVector3 X_pperp3D(XCand_fromMCFit->currentState().globalMomentum().x(), XCand_fromMCFit->currentState().globalMomentum().y(), XCand_fromMCFit->currentState().globalMomentum().z());
                  AlgebraicVector3 X_v3pperp ;
                  X_v3pperp[0] = X_pperp.x(); X_v3pperp[1] = X_pperp.y(); X_v3pperp[2] = 0.;
                  TVector3 X_pvtx, X_pvtx3D, X_vdiff, X_vdiff3D ;
                  double X_cosAlpha, X_cosAlpha3D, X_ctau ;
                  VertexDistanceXY X_vdistXY ;
                  Measurement1D X_distXY ;
                  GlobalError X_v1e = (Vertex(*XCand_vertex_fromMCFit)).error();
                  GlobalError X_v2e ;
                  AlgebraicSymMatrix33 X_vXYe ;
                  double X_ctauErr ;
                  float X_lxy, X_lxyErr, X_lxyz, X_lxyzErr ;
                  ROOT::Math::SVector<double, 3> X_vDiff, X_vDiff3D ;


		  ////////////////// Lifetime wrt PV for Bs0 & X(4140) //////////////////
		  Bs0_v2e = thePrimaryVtx.error(); /// Bs0
		  Bs0_vXYe = Bs0_v1e.matrix() + Bs0_v2e.matrix() ; 
		  /// 2D
		  Bs0_pvtx.SetXYZ(thePrimaryVtx.position().x(), thePrimaryVtx.position().y(), 0) ; 
		  Bs0_vdiff = Bs0_vtx - Bs0_pvtx ; 
		  Bs0_cosAlpha = Bs0_vdiff.Dot(Bs0_pperp) / (Bs0_vdiff.Perp()*Bs0_pperp.Perp()) ; 
		  Bs0_lxy = Bs0_vdiff.Perp(); 
		  Bs0_vDiff[0] = Bs0_vdiff.x(); Bs0_vDiff[1] = Bs0_vdiff.y(); Bs0_vDiff[2] = 0 ; // needed by Similarity method 
		  Bs0_lxyErr = sqrt(ROOT::Math::Similarity(Bs0_vDiff,Bs0_vXYe)) / Bs0_vdiff.Perp(); 
		  Bs0_distXY = Bs0_vdistXY.distance(Vertex(*Bs0Cand_vertex_fromMCFit), Vertex(thePrimaryVtx)); 
		  Bs0_ctau = Bs0_distXY.value() * Bs0_cosAlpha * Bs0Cand_fromMCFit->currentState().mass() / Bs0_pperp.Perp(); 
		  Bs0_ctauErr = sqrt(ROOT::Math::Similarity(Bs0_v3pperp,Bs0_vXYe)) * Bs0Cand_fromMCFit->currentState().mass() / (Bs0_pperp.Perp2()) ; 
		  /// 3D
		  Bs0_pvtx3D.SetXYZ(thePrimaryVtx.position().x(), thePrimaryVtx.position().y(), thePrimaryVtx.position().z()); 
		  Bs0_vdiff3D = Bs0_vtx3D - Bs0_pvtx3D; 
		  Bs0_cosAlpha3D = Bs0_vdiff3D.Dot(Bs0_pperp3D)/(Bs0_vdiff3D.Mag()*Bs0_pperp3D.Mag()); 
		  Bs0_lxyz = Bs0_vdiff3D.Mag(); 
		  Bs0_vDiff3D[0] = Bs0_vdiff3D.x(); Bs0_vDiff3D[1] = Bs0_vdiff3D.y(); Bs0_vDiff3D[2] = Bs0_vdiff3D.z() ;
		  Bs0_lxyzErr = sqrt(ROOT::Math::Similarity(Bs0_vDiff3D,Bs0_vXYe)) / Bs0_vdiff3D.Mag(); 

		  X_v2e = thePrimaryVtx.error(); /// X(4140)
                  X_vXYe = X_v1e.matrix() + X_v2e.matrix() ;  		  
		  /// 2D
		  X_pvtx.SetXYZ(thePrimaryVtx.position().x(), thePrimaryVtx.position().y(), 0) ;
		  X_vdiff = X_vtx - X_pvtx ;
		  X_cosAlpha = X_vdiff.Dot(X_pperp) / (X_vdiff.Perp()*X_pperp.Perp()) ;  
		  X_lxy = X_vdiff.Perp();
		  X_vDiff[0] = X_vdiff.x(); X_vDiff[1] = X_vdiff.y(); X_vDiff[2] = 0 ;  
                  X_lxyErr = sqrt(ROOT::Math::Similarity(X_vDiff,X_vXYe)) / X_vdiff.Perp();
		  X_distXY = X_vdistXY.distance(Vertex(*XCand_vertex_fromMCFit), Vertex(thePrimaryVtx));
                  X_ctau = X_distXY.value() * X_cosAlpha * XCand_fromMCFit->currentState().mass() / X_pperp.Perp();
                  X_ctauErr = sqrt(ROOT::Math::Similarity(X_v3pperp,X_vXYe)) * XCand_fromMCFit->currentState().mass() / (X_pperp.Perp2()) ;
		  /// 3D
		  X_pvtx3D.SetXYZ(thePrimaryVtx.position().x(), thePrimaryVtx.position().y(), thePrimaryVtx.position().z());
                  X_vdiff3D = X_vtx3D - X_pvtx3D;
                  X_cosAlpha3D = X_vdiff3D.Dot(X_pperp3D)/(X_vdiff3D.Mag()*X_pperp3D.Mag());
                  X_lxyz = X_vdiff3D.Mag();
                  X_vDiff3D[0] = X_vdiff3D.x(); X_vDiff3D[1] = X_vdiff3D.y(); X_vDiff3D[2] = X_vdiff3D.z() ;
                  X_lxyzErr = sqrt(ROOT::Math::Similarity(X_vDiff3D,X_vXYe)) / X_vdiff3D.Mag(); 


		  ////////////////// Last cuts for Bs0 & X(4140 )//////////////////
		  if ( !(Bs0_ctau/Bs0_ctauErr > 2.8) || !(Bs0_cosAlpha > 0.8) )
		    continue ;
		  if ( !(X_ctau/X_ctauErr > 2.8) || !(X_cosAlpha > 0.8) )
                    continue ;


		  ////////////////// fill Bs0 & X(4140) candidate variables //////////////////
		  bs0Mass->push_back( Bs0Cand_fromMCFit->currentState().mass()) ; /// Bs0
		  bs0Px->push_back( Bs0Cand_fromMCFit->currentState().globalMomentum().x()) ;
		  bs0Py->push_back( Bs0Cand_fromMCFit->currentState().globalMomentum().y()) ;
		  bs0Pz->push_back( Bs0Cand_fromMCFit->currentState().globalMomentum().z()) ;			
		  bs0PxE->push_back( sqrt( Bs0Cand_fromMCFit->currentState().kinematicParametersError().matrix()(3,3) ) ) ;
		  bs0PyE->push_back( sqrt( Bs0Cand_fromMCFit->currentState().kinematicParametersError().matrix()(4,4) ) ) ;
		  bs0PzE->push_back( sqrt( Bs0Cand_fromMCFit->currentState().kinematicParametersError().matrix()(5,5) ) ) ;
		  bs0Vtx_CL->push_back( bs0VtxProb );
		  bs0Vtx_Chi2->push_back( Bs0Cand_vertex_fromMCFit->chiSquared() ) ;
		  bs0DecayVtx_X->push_back((*Bs0Cand_vertex_fromMCFit).position().x());
		  bs0DecayVtx_Y->push_back((*Bs0Cand_vertex_fromMCFit).position().y());
		  bs0DecayVtx_Z->push_back((*Bs0Cand_vertex_fromMCFit).position().z());
		  bs0DecayVtx_XE->push_back(sqrt((*Bs0Cand_vertex_fromMCFit).error().cxx()));
		  bs0DecayVtx_YE->push_back(sqrt((*Bs0Cand_vertex_fromMCFit).error().cyy()));
		  bs0DecayVtx_ZE->push_back(sqrt((*Bs0Cand_vertex_fromMCFit).error().czz()));
		  Bs0VertexFitTree->movePointerToTheFirstChild(); 
		  RefCountedKinematicParticle Bs0_mu1_MuMuKK = Bs0VertexFitTree->currentParticle();
		  Bs0VertexFitTree->movePointerToTheNextChild();
		  RefCountedKinematicParticle Bs0_mu2_MuMuKK = Bs0VertexFitTree->currentParticle();
		  Bs0VertexFitTree->movePointerToTheNextChild();
		  RefCountedKinematicParticle Bs0_k1_MuMuKK = Bs0VertexFitTree->currentParticle();
		  Bs0VertexFitTree->movePointerToTheNextChild();
		  RefCountedKinematicParticle Bs0_k2_MuMuKK = Bs0VertexFitTree->currentParticle();
		  /// muon1 & muon2			
		  Bs0_mu1Px_MuMuKK->push_back( Bs0_mu1_MuMuKK->currentState().globalMomentum().x() );
		  Bs0_mu1Py_MuMuKK->push_back( Bs0_mu1_MuMuKK->currentState().globalMomentum().y() );
		  Bs0_mu1Pz_MuMuKK->push_back( Bs0_mu1_MuMuKK->currentState().globalMomentum().z() );
		  Bs0_mu1E_MuMuKK->push_back( Bs0_mu1_MuMuKK->currentState().kinematicParameters().energy() );
		  Bs0_mu2Px_MuMuKK->push_back( Bs0_mu2_MuMuKK->currentState().globalMomentum().x() );
		  Bs0_mu2Py_MuMuKK->push_back( Bs0_mu2_MuMuKK->currentState().globalMomentum().y() );
		  Bs0_mu2Pz_MuMuKK->push_back( Bs0_mu2_MuMuKK->currentState().globalMomentum().z() );
		  Bs0_mu2E_MuMuKK->push_back( Bs0_mu2_MuMuKK->currentState().kinematicParameters().energy() );
		  /// kaon1 & kaon2
		  Bs0_k1Px_MuMuKK->push_back( Bs0_k1_MuMuKK->currentState().globalMomentum().x() );
		  Bs0_k1Py_MuMuKK->push_back( Bs0_k1_MuMuKK->currentState().globalMomentum().y() );
		  Bs0_k1Pz_MuMuKK->push_back( Bs0_k1_MuMuKK->currentState().globalMomentum().z() );
		  Bs0_k1E_MuMuKK->push_back( Bs0_k1_MuMuKK->currentState().kinematicParameters().energy() );
		  Double_t theo = 0., sigma = 0. ;
		  Bs0_kaon1_nsigdedx->push_back( nsigmaofdedx(Track1->track(),theo,sigma) );
		  Bs0_kaon1_dedx->push_back( getEnergyLoss(Track1->track()) );
		  Bs0_kaon1_dedxMass->push_back( GetMass(Track1->track()) );
		  Bs0_kaon1_theo->push_back( theo );
		  Bs0_kaon1_sigma->push_back( sigma );
		  Bs0_kaon1_dedx_byHits->push_back( (dEdxTrack)[Track1->track()].dEdx() );
		  Bs0_kaon1_dedxErr_byHits->push_back( (dEdxTrack)[Track1->track()].dEdxError() );
		  Bs0_kaon1_saturMeas_byHits->push_back( (dEdxTrack)[Track1->track()].numberOfSaturatedMeasurements() );
		  Bs0_kaon1_Meas_byHits->push_back( (dEdxTrack)[Track1->track()].numberOfMeasurements() );
		  Bs0_k2Px_MuMuKK->push_back( Bs0_k2_MuMuKK->currentState().globalMomentum().x() );
		  Bs0_k2Py_MuMuKK->push_back( Bs0_k2_MuMuKK->currentState().globalMomentum().y() );
		  Bs0_k2Pz_MuMuKK->push_back( Bs0_k2_MuMuKK->currentState().globalMomentum().z() );
		  Bs0_k2E_MuMuKK->push_back( Bs0_k2_MuMuKK->currentState().kinematicParameters().energy() );
		  theo = 0.; sigma = 0. ;
		  Bs0_kaon2_nsigdedx->push_back(nsigmaofdedx(Track2->track(),theo,sigma));
		  Bs0_kaon2_dedx->push_back(getEnergyLoss(Track2->track()));
		  Bs0_kaon2_dedxMass->push_back(GetMass(Track2->track()));
		  Bs0_kaon2_theo->push_back(theo);
		  Bs0_kaon2_sigma->push_back(sigma);
		  Bs0_kaon2_dedx_byHits->push_back( (dEdxTrack_Kaon)[Track2->track()].dEdx() );
		  Bs0_kaon2_dedxErr_byHits->push_back( (dEdxTrack_Kaon)[Track2->track()].dEdxError() );
		  Bs0_kaon2_saturMeas_byHits->push_back( (dEdxTrack_Kaon)[Track2->track()].numberOfSaturatedMeasurements() );
		  Bs0_kaon2_Meas_byHits->push_back( (dEdxTrack_Kaon)[Track2->track()].numberOfMeasurements() );
		  /// PV 
		  bs0CosAlphaPV->push_back( Bs0_cosAlpha ); bs0CosAlpha3DPV->push_back( Bs0_cosAlpha3D );
		  bs0CTauPV->push_back( Bs0_ctau ); bs0CTauPVE->push_back( Bs0_ctauErr );
		  bs0LxyPV->push_back( Bs0_lxy ); bs0LxyPVE->push_back( Bs0_lxyErr );
		  bs0LxyzPV->push_back( Bs0_lxyz ); bs0LxyzPVE->push_back( Bs0_lxyzErr );

		  xMass->push_back( XCand_fromMCFit->currentState().mass()) ; /// X(4140)
                  xPx->push_back( XCand_fromMCFit->currentState().globalMomentum().x()) ;
                  xPy->push_back( XCand_fromMCFit->currentState().globalMomentum().y()) ;
                  xPz->push_back( XCand_fromMCFit->currentState().globalMomentum().z()) ;
                  xPxE->push_back( sqrt( XCand_fromMCFit->currentState().kinematicParametersError().matrix()(3,3) ) ) ;
                  xPyE->push_back( sqrt( XCand_fromMCFit->currentState().kinematicParametersError().matrix()(4,4) ) ) ;
                  xPzE->push_back( sqrt( XCand_fromMCFit->currentState().kinematicParametersError().matrix()(5,5) ) ) ;
                  xVtx_CL->push_back( xVtxProb );
                  xVtx_Chi2->push_back( XCand_vertex_fromMCFit->chiSquared() ) ;
                  xDecayVtx_X->push_back((*XCand_vertex_fromMCFit).position().x());
                  xDecayVtx_Y->push_back((*XCand_vertex_fromMCFit).position().y());
                  xDecayVtx_Z->push_back((*XCand_vertex_fromMCFit).position().z());
                  xDecayVtx_XE->push_back(sqrt((*XCand_vertex_fromMCFit).error().cxx()));
                  xDecayVtx_YE->push_back(sqrt((*XCand_vertex_fromMCFit).error().cyy()));
                  xDecayVtx_ZE->push_back(sqrt((*XCand_vertex_fromMCFit).error().czz()));
                  XVertexFitTree->movePointerToTheFirstChild();
                  RefCountedKinematicParticle X_mu1_MuMuKK = XVertexFitTree->currentParticle();
                  XVertexFitTree->movePointerToTheNextChild();
                  RefCountedKinematicParticle X_mu2_MuMuKK = XVertexFitTree->currentParticle();
                  XVertexFitTree->movePointerToTheNextChild();
                  RefCountedKinematicParticle X_k1_MuMuKK = XVertexFitTree->currentParticle();
                  XVertexFitTree->movePointerToTheNextChild();
                  RefCountedKinematicParticle X_k2_MuMuKK = XVertexFitTree->currentParticle();		  
		  /// muon1 & muon2
		  X_mu1Px_MuMuKK->push_back( X_mu1_MuMuKK->currentState().globalMomentum().x() );
                  X_mu1Py_MuMuKK->push_back( X_mu1_MuMuKK->currentState().globalMomentum().y() );
                  X_mu1Pz_MuMuKK->push_back( X_mu1_MuMuKK->currentState().globalMomentum().z() );
                  X_mu1E_MuMuKK->push_back( X_mu1_MuMuKK->currentState().kinematicParameters().energy() );
		  X_mu2Px_MuMuKK->push_back( X_mu2_MuMuKK->currentState().globalMomentum().x() );
                  X_mu2Py_MuMuKK->push_back( X_mu2_MuMuKK->currentState().globalMomentum().y() );
                  X_mu2Pz_MuMuKK->push_back( X_mu2_MuMuKK->currentState().globalMomentum().z() );
                  X_mu2E_MuMuKK->push_back( X_mu2_MuMuKK->currentState().kinematicParameters().energy() );
		  /// kaon1 & kaon2
		  X_k1Px_MuMuKK->push_back( X_k1_MuMuKK->currentState().globalMomentum().x() );
                  X_k1Py_MuMuKK->push_back( X_k1_MuMuKK->currentState().globalMomentum().y() );
                  X_k1Pz_MuMuKK->push_back( X_k1_MuMuKK->currentState().globalMomentum().z() );
                  X_k1E_MuMuKK->push_back( X_k1_MuMuKK->currentState().kinematicParameters().energy() );
                  X_kaon1_nsigdedx->push_back( nsigmaofdedx(Track1->track(),theo,sigma) );
                  X_kaon1_dedx->push_back( getEnergyLoss(Track1->track()) );
                  X_kaon1_dedxMass->push_back( GetMass(Track1->track()) );
                  X_kaon1_theo->push_back( theo );
                  X_kaon1_sigma->push_back( sigma );
                  X_kaon1_dedx_byHits->push_back( (dEdxTrack)[Track1->track()].dEdx() );
                  X_kaon1_dedxErr_byHits->push_back( (dEdxTrack)[Track1->track()].dEdxError() );
                  X_kaon1_saturMeas_byHits->push_back( (dEdxTrack)[Track1->track()].numberOfSaturatedMeasurements() );
                  X_kaon1_Meas_byHits->push_back( (dEdxTrack)[Track1->track()].numberOfMeasurements() );
		  X_k2Px_MuMuKK->push_back( X_k2_MuMuKK->currentState().globalMomentum().x() );
                  X_k2Py_MuMuKK->push_back( X_k2_MuMuKK->currentState().globalMomentum().y() );
                  X_k2Pz_MuMuKK->push_back( X_k2_MuMuKK->currentState().globalMomentum().z() );
                  X_k2E_MuMuKK->push_back( X_k2_MuMuKK->currentState().kinematicParameters().energy() );
                  X_kaon2_nsigdedx->push_back(nsigmaofdedx(Track2->track(),theo,sigma));
                  X_kaon2_dedx->push_back(getEnergyLoss(Track2->track()));
                  X_kaon2_dedxMass->push_back(GetMass(Track2->track()));
                  X_kaon2_theo->push_back(theo);
                  X_kaon2_sigma->push_back(sigma);
                  X_kaon2_dedx_byHits->push_back( (dEdxTrack_Kaon)[Track2->track()].dEdx() );
                  X_kaon2_dedxErr_byHits->push_back( (dEdxTrack_Kaon)[Track2->track()].dEdxError() );
                  X_kaon2_saturMeas_byHits->push_back( (dEdxTrack_Kaon)[Track2->track()].numberOfSaturatedMeasurements() );
		  X_kaon2_Meas_byHits->push_back( (dEdxTrack_Kaon)[Track2->track()].numberOfMeasurements() );
		  /// PV
		  xCosAlphaPV->push_back( X_cosAlpha ); xCosAlpha3DPV->push_back( X_cosAlpha3D );
                  xCTauPV->push_back( X_ctau ); xCTauPVE->push_back( X_ctauErr );
                  xLxyPV->push_back( X_lxy ); xLxyPVE->push_back( X_lxyErr );
                  xLxyzPV->push_back( X_lxyz ); xLxyzPVE->push_back( X_lxyzErr );
	          /// dxy, dz, dxyE, dzE for kaons from PV 
                  kaon1_dxy_PV->push_back( Track1->track()->dxy(RefVtx) ); 
		  kaon1_dz_PV->push_back( Track1->track()->dz(RefVtx) );
                  kaon2_dxy_PV->push_back( Track2->track()->dxy(RefVtx) );
                  kaon2_dz_PV->push_back( Track2->track()->dz(RefVtx) );


		  ////////////////// Lifetime wrt BS for Bs0 & X(4140) //////////////////
		  Bs0_v2e = theBeamSpotVtx.error(); /// Bs0
		  Bs0_vXYe = Bs0_v1e.matrix() + Bs0_v2e.matrix();
		  /// 2D
		  Bs0_pvtx.SetXYZ(theBeamSpotVtx.position().x(), theBeamSpotVtx.position().y(), 0);
		  Bs0_vdiff = Bs0_vtx - Bs0_pvtx;
		  Bs0_cosAlpha = Bs0_vdiff.Dot(Bs0_pperp)/(Bs0_vdiff.Perp()*Bs0_pperp.Perp());
		  Bs0_lxy = Bs0_vdiff.Perp();
		  Bs0_vDiff[0] = Bs0_vdiff.x(); Bs0_vDiff[1] = Bs0_vdiff.y(); Bs0_vDiff[2] = 0 ; // needed by Similarity method
		  Bs0_lxyErr = sqrt(ROOT::Math::Similarity(Bs0_vDiff,Bs0_vXYe)) / Bs0_vdiff.Perp();
		  Bs0_distXY = Bs0_vdistXY.distance(Vertex(*Bs0Cand_vertex_fromMCFit), Vertex(theBeamSpotVtx));
		  Bs0_ctau = Bs0_distXY.value() * Bs0_cosAlpha * (Bs0Cand_fromMCFit->currentState().mass() / Bs0_pperp.Perp()) ;
		  Bs0_ctauErr = sqrt(ROOT::Math::Similarity(Bs0_v3pperp,Bs0_vXYe)) * Bs0Cand_fromMCFit->currentState().mass()/Bs0_pperp.Perp2();
		  /// 3D
		  Bs0_pvtx3D.SetXYZ(theBeamSpotVtx.position().x(), theBeamSpotVtx.position().y(), theBeamSpotVtx.position().z());
		  Bs0_vdiff3D = Bs0_vtx3D - Bs0_pvtx3D;
		  Bs0_cosAlpha3D = Bs0_vdiff3D.Dot(Bs0_pperp3D)/(Bs0_vdiff3D.Mag()*Bs0_pperp3D.Mag());
		  Bs0_lxyz = Bs0_vdiff3D.Mag();
		  Bs0_vDiff3D[0] = Bs0_vdiff3D.x(); Bs0_vDiff3D[1] = Bs0_vdiff3D.y(); Bs0_vDiff3D[2] = Bs0_vdiff3D.z() ;
		  Bs0_lxyzErr = sqrt(ROOT::Math::Similarity(Bs0_vDiff3D,Bs0_vXYe)) / Bs0_vdiff3D.Mag();
		  		  				
		  X_v2e = theBeamSpotVtx.error(); /// X(4140)
                  X_vXYe = X_v1e.matrix() + X_v2e.matrix();
		  /// 2D
		  X_pvtx.SetXYZ(theBeamSpotVtx.position().x(), theBeamSpotVtx.position().y(), 0);
                  X_vdiff = X_vtx - X_pvtx;
                  X_cosAlpha = X_vdiff.Dot(X_pperp)/(X_vdiff.Perp()*X_pperp.Perp());
                  X_lxy = X_vdiff.Perp();
                  X_vDiff[0] = X_vdiff.x(); X_vDiff[1] = X_vdiff.y(); X_vDiff[2] = 0 ; 
                  X_lxyErr = sqrt(ROOT::Math::Similarity(X_vDiff,X_vXYe)) / X_vdiff.Perp();
                  X_distXY = X_vdistXY.distance(Vertex(*XCand_vertex_fromMCFit), Vertex(theBeamSpotVtx));
                  X_ctau = X_distXY.value() * X_cosAlpha * (XCand_fromMCFit->currentState().mass() / X_pperp.Perp()) ;
                  X_ctauErr = sqrt(ROOT::Math::Similarity(X_v3pperp,X_vXYe)) * XCand_fromMCFit->currentState().mass()/X_pperp.Perp2();
		  /// 3D
		  X_pvtx3D.SetXYZ(theBeamSpotVtx.position().x(), theBeamSpotVtx.position().y(), theBeamSpotVtx.position().z());
                  X_vdiff3D = X_vtx3D - X_pvtx3D;
                  X_cosAlpha3D = X_vdiff3D.Dot(X_pperp3D)/(X_vdiff3D.Mag()*X_pperp3D.Mag());
                  X_lxyz = X_vdiff3D.Mag();
                  X_vDiff3D[0] = X_vdiff3D.x(); X_vDiff3D[1] = X_vdiff3D.y(); X_vDiff3D[2] = X_vdiff3D.z() ;
                  X_lxyzErr = sqrt(ROOT::Math::Similarity(X_vDiff3D,X_vXYe)) / X_vdiff3D.Mag();


		  ////////////////// BS (beam spot) for Bs0 ////////////////// 	
		  bs0CosAlphaBS->push_back( Bs0_cosAlpha ); bs0CosAlpha3DBS->push_back( Bs0_cosAlpha3D );
		  bs0CTauBS->push_back( Bs0_ctau ); bs0CTauBSE->push_back( Bs0_ctauErr );
		  bs0LxyBS->push_back( Bs0_lxy ); bs0LxyBSE->push_back( Bs0_lxyErr );
		  bs0LxyzBS->push_back( Bs0_lxyz ); bs0LxyzBSE->push_back( Bs0_lxyzErr );

		  vector<TransientVertex> Bs0_pvs ;  
		  Vertex Bs0LessPV = thePrimaryVtx ;

		  if (addBs0lessPrimaryVertex_) 
		    {
		      VertexReProducer revertex(recVtxs, iEvent);
		      Handle<TrackCollection> pvtracks;   
		      iEvent.getByLabel(revertex.inputTracks(), pvtracks);
		      Handle<BeamSpot>        pvbeamspot;
		      iEvent.getByLabel(revertex.inputBeamSpot(), pvbeamspot);
		      
		      if (pvbeamspot.id() != beamSpotHandle.id() ) 
			edm::LogWarning("Inconsistency") << "The BeamSpot used for PV reco is not the same used in this analyzer.";
		      
		      const reco::Muon *Bs0rmu_1 = dynamic_cast<const reco::Muon *>(Muon1->originalObject());
		      const reco::Muon *Bs0rmu_2 = dynamic_cast<const reco::Muon *>(Muon2->originalObject());
		      
		      if (Bs0rmu_1 != 0  &&  Bs0rmu_2 != 0  &&  Bs0rmu_1->track().id() == pvtracks.id()  &&  Bs0rmu_2->track().id() == pvtracks.id() 
			  &&  Track1->track().id() == pvtracks.id()  &&  Track2->track().id() ==  pvtracks.id()) { 
			vector<TransientTrack> Bs0Less; // need TransientTrack to keep the TrackRef
			Bs0Less.reserve( pvtracks->size() );
			Double_t removedTrksPtSq = 0. ;
			for (size_t i = 0, n = pvtracks->size(); i < n; ++i) { 
			  if (i == Bs0rmu_1->track().key()) { removedTrksPtSq += (Bs0rmu_1->track()->pt())*(Bs0rmu_1->track()->pt()) ;
			    continue; }
			  if (i == Bs0rmu_2->track().key()) { removedTrksPtSq += (Bs0rmu_2->track()->pt())*(Bs0rmu_2->track()->pt()) ;
			    continue; }
			  if (i == Track1->track().key()) { removedTrksPtSq += (Track1->track()->pt())*(Track1->track()->pt()) ;
			    continue; }
			  if (i == Track2->track().key()) { removedTrksPtSq += (Track2->track()->pt())*(Track2->track()->pt()) ;
			    continue; } 

			  reco::TrackRef trk_now(pvtracks, i) ;
			  TransientTrack transientTrack = theTTBuilder->build( trk_now ); 
			  transientTrack.setBeamSpot( beamSpot );
			  Bs0Less.push_back( transientTrack );
			}
			if ( removedTrksPtSq > 0. ) {
			  Bs0_pvs = revertex.makeVertices(Bs0Less, *pvbeamspot, iSetup) ; // list of PV
			} else
			  if (Debug_) cout <<"\n\\\\\\\\\\\\\\\\\\\\ excluded tracks pT^2 = 0 \\\\\\\\\\\\\\\\\\\\\n" <<endl ;
			if ( !Bs0_pvs.empty() ) {
			  Bs0LessPV = Vertex(Bs0_pvs.front());
			  Bs0LessPV_tracksPtSq->push_back( vertexHigherPtSquared.sumPtSquared(Bs0LessPV) ) ;
			  Bs0LessPV_4tracksPtSq->push_back( removedTrksPtSq ) ;
			  if (Debug_) {
			    cout <<"\nBs0LessPV_z = " <<Bs0LessPV.position().z() <<endl ;
			    cout <<"Bs0LessPV_tracks = " <<Bs0LessPV.tracksSize() <<endl ;
			    cout <<"Bs0LessPV_tracksPtSq = " <<vertexHigherPtSquared.sumPtSquared(Bs0LessPV) <<endl ;
			    cout <<"Bs0LessPV_removedTracksPtSq = " <<removedTrksPtSq <<endl ;
			    cout <<"Bs0_pvs->size() = " <<Bs0_pvs.size() <<endl ;
			    cout <<"priVtx_z = " << priVtx_Z <<endl ;
			    cout <<"priVtx_tracks = " <<priVtx_tracks <<endl ;
			    cout <<"priVtx_tracksPtSq = " <<priVtx_tracksPtSq <<endl ;
			    cout <<"recVtxs->size() = " <<recVtxs->size() <<endl ;
			  }
			}
		      }
		    }

		  PriVtxBs0Less_n->push_back( Bs0_pvs.size() ) ;
		  PriVtxBs0Less_X->push_back( Bs0LessPV.position().x() ) ;
		  PriVtxBs0Less_Y->push_back( Bs0LessPV.position().y() ) ;
		  PriVtxBs0Less_Z->push_back( Bs0LessPV.position().z() ) ; 
		  PriVtxBs0Less_EX->push_back( Bs0LessPV.xError() ) ;
		  PriVtxBs0Less_EY->push_back( Bs0LessPV.yError() ) ;
		  PriVtxBs0Less_EZ->push_back( Bs0LessPV.zError() ) ;
		  PriVtxBs0Less_CL->push_back( ChiSquaredProbability( (double)(Bs0LessPV.chi2()), (double)(Bs0LessPV.ndof())) );
		  PriVtxBs0Less_Chi2->push_back( Bs0LessPV.chi2() ) ;
		  PriVtxBs0Less_tracks->push_back( Bs0LessPV.tracksSize() ) ;
		 
		  /// dxy, dz, dxyE, dzE for kaons from BS
    		  math::XYZPoint BSVtx;
		  BSVtx = theBeamSpotVtx.position();		
	          kaon1_dxy_BS->push_back( Track1->track()->dxy(BSVtx) ); 
                  kaon1_dz_BS->push_back( Track1->track()->dz(BSVtx) );
                  kaon2_dxy_BS->push_back( Track2->track()->dxy(BSVtx) );
                  kaon2_dz_BS->push_back( Track2->track()->dz(BSVtx) );


		  ////////////////// Lifetime wrt Bs0LessPV for Bs0 ////////////////// 
		  Bs0_v2e = Bs0LessPV.error();
		  Bs0_vXYe = Bs0_v1e.matrix() + Bs0_v2e.matrix();
		  /// 2D
		  Bs0_pvtx.SetXYZ(Bs0LessPV.position().x(), Bs0LessPV.position().y(), 0) ;	
		  Bs0_vdiff = Bs0_vtx - Bs0_pvtx ;
		  Bs0_cosAlpha = Bs0_vdiff.Dot(Bs0_pperp)/(Bs0_vdiff.Perp()*Bs0_pperp.Perp());
		  Bs0_lxy = Bs0_vdiff.Perp();
		  Bs0_vDiff[0] = Bs0_vdiff.x(); Bs0_vDiff[1] = Bs0_vdiff.y(); Bs0_vDiff[2] = 0 ; // needed by Similarity method
		  Bs0_lxyErr = sqrt(ROOT::Math::Similarity(Bs0_vDiff,Bs0_vXYe)) / Bs0_vdiff.Perp();
		  Bs0_distXY = Bs0_vdistXY.distance( Vertex(*Bs0Cand_vertex_fromMCFit), Vertex(Bs0LessPV) ) ;
		  Bs0_ctau = Bs0_distXY.value() * Bs0_cosAlpha * Bs0Cand_fromMCFit->currentState().mass() / Bs0_pperp.Perp();
		  Bs0_ctauErr = sqrt(ROOT::Math::Similarity(Bs0_v3pperp,Bs0_vXYe)) * Bs0Cand_fromMCFit->currentState().mass() / (Bs0_pperp.Perp2());
		  /// 3D
		  Bs0_pvtx3D.SetXYZ(Bs0LessPV.position().x(), Bs0LessPV.position().y(), Bs0LessPV.position().z());
		  Bs0_vdiff3D = Bs0_vtx3D - Bs0_pvtx3D;
		  Bs0_cosAlpha3D = Bs0_vdiff3D.Dot(Bs0_pperp3D)/( Bs0_vdiff3D.Mag()*Bs0_pperp3D.Mag() );
		  Bs0_lxyz = Bs0_vdiff3D.Mag();
		  Bs0_vDiff3D[0] = Bs0_vdiff3D.x(); Bs0_vDiff3D[1] = Bs0_vdiff3D.y(); Bs0_vDiff3D[2] = Bs0_vdiff3D.z() ;
		  Bs0_lxyzErr = sqrt(ROOT::Math::Similarity(Bs0_vDiff3D,Bs0_vXYe)) / Bs0_vdiff3D.Mag();

		  bs0CosAlphaBs0LessPV->push_back( Bs0_cosAlpha ) ; bs0CosAlpha3DBs0LessPV->push_back( Bs0_cosAlpha3D ) ;
		  bs0CTauBs0LessPV->push_back( Bs0_ctau ) ; bs0CTauBs0LessPVE->push_back( Bs0_ctauErr ) ;
		  bs0LxyBs0LessPV->push_back( Bs0_lxy ) ; bs0LxyBs0LessPVE->push_back( Bs0_lxyErr ) ;
		  bs0LxyzBs0LessPV->push_back( Bs0_lxyz ) ; bs0LxyzBs0LessPVE->push_back( Bs0_lxyzErr ) ;

		  /// dxy, dz, dxyE, dzE for kaons from Bs0LessPV
		  math::XYZPoint Bs0LessPVvtx; 
     		  Bs0LessPVvtx = Bs0LessPV.position();
                  kaon1_dxy_Bs0LessPV->push_back( Track1->track()->dxy(Bs0LessPVvtx) );
                  kaon1_dz_Bs0LessPV->push_back( Track1->track()->dz(Bs0LessPVvtx) );
                  kaon2_dxy_Bs0LessPV->push_back( Track2->track()->dxy(Bs0LessPVvtx) );
                  kaon2_dz_Bs0LessPV->push_back( Track2->track()->dz(Bs0LessPVvtx) );
 
                  kaon1_dxyE->push_back( Track1->track()->dxyError() );
                  kaon1_dzE->push_back( Track1->track()->dzError() );
		  kaon2_dxyE->push_back( Track2->track()->dxyError() );
		  kaon2_dzE->push_back( Track2->track()->dzError() );

 
		  /// Find the PV among the original offlinePV with the largest Bs0_cos(alpha)
		  Vertex theCosAlphaV = thePrimaryVtx ; 
		  float maxCosAlpha = -1. ;
		  
		  for (VertexCollection::const_iterator itv = recVtxs->begin(), itvend = recVtxs->end(); itv != itvend; ++itv) {
		    Bs0_pvtx.SetXYZ(itv->position().x(), itv->position().y(), 0) ;	
		    Bs0_vdiff = Bs0_vtx - Bs0_pvtx ;
		    float cosAlpha_temp = Bs0_vdiff.Dot(Bs0_pperp) / (Bs0_vdiff.Perp()*Bs0_pperp.Perp()) ; // Perp() == Mag() when z = 0
		  
		    if ( cosAlpha_temp > maxCosAlpha ) {
		      maxCosAlpha = cosAlpha_temp ;    
		      theCosAlphaV = Vertex(*itv) ;
		    }
		  }
		  
		  PriVtx_Bs0CosAlpha_n->push_back( recVtxs->size() ) ;
		  PriVtx_Bs0CosAlpha_X->push_back( theCosAlphaV.position().x() ) ;
		  PriVtx_Bs0CosAlpha_Y->push_back( theCosAlphaV.position().y() ) ;
		  PriVtx_Bs0CosAlpha_Z->push_back( theCosAlphaV.position().z() ) ;
		  PriVtx_Bs0CosAlpha_EX->push_back( theCosAlphaV.xError() ) ;
		  PriVtx_Bs0CosAlpha_EY->push_back( theCosAlphaV.yError() ) ;
		  PriVtx_Bs0CosAlpha_EZ->push_back( theCosAlphaV.zError() ) ;
		  PriVtx_Bs0CosAlpha_CL->push_back( ChiSquaredProbability((double)(theCosAlphaV.chi2()), (double)(theCosAlphaV.ndof())) ) ;
		  PriVtx_Bs0CosAlpha_Chi2->push_back( theCosAlphaV.chi2() ) ;
		  PriVtx_Bs0CosAlpha_tracks->push_back( theCosAlphaV.tracksSize() ) ;


		  ////////////////// Lifetime wrt PV with largest Bs0_cos(alpha) & X_cos(alpha) candidate //////////////////  
		  Bs0_v2e = theCosAlphaV.error(); /// Bs0
		  Bs0_vXYe = Bs0_v1e.matrix() + Bs0_v2e.matrix();
		  /// 2D
		  Bs0_pvtx.SetXYZ(theCosAlphaV.position().x(), theCosAlphaV.position().y(), 0) ;	
		  Bs0_vdiff = Bs0_vtx - Bs0_pvtx ;
		  Bs0_cosAlpha =  maxCosAlpha ;
		  Bs0_lxy = Bs0_vdiff.Perp();
		  Bs0_vDiff[0] = Bs0_vdiff.x(); Bs0_vDiff[1] = Bs0_vdiff.y(); Bs0_vDiff[2] = 0 ; // needed by Similarity method
		  Bs0_lxyErr = sqrt(ROOT::Math::Similarity(Bs0_vDiff,Bs0_vXYe)) / Bs0_vdiff.Perp();
		  Bs0_distXY = Bs0_vdistXY.distance( Vertex(*Bs0Cand_vertex_fromMCFit), Vertex(theCosAlphaV) ) ;
		  Bs0_ctau = Bs0_distXY.value() * Bs0_cosAlpha * Bs0Cand_fromMCFit->currentState().mass() / Bs0_pperp.Perp();
		  Bs0_ctauErr = sqrt(ROOT::Math::Similarity(Bs0_v3pperp,Bs0_vXYe)) * Bs0Cand_fromMCFit->currentState().mass() / (Bs0_pperp.Perp2());
		  Bs0_lxy = Bs0_vdiff.Dot(Bs0_pperp) / Bs0_pperp.Mag() ;
		  /// 3D
		  Bs0_pvtx3D.SetXYZ(theCosAlphaV.position().x(), theCosAlphaV.position().y(), theCosAlphaV.position().z());
		  Bs0_vdiff3D = Bs0_vtx3D - Bs0_pvtx3D;
		  Bs0_cosAlpha3D = Bs0_vdiff3D.Dot(Bs0_pperp3D)/( Bs0_vdiff3D.Mag()*Bs0_pperp3D.Mag() );
		  Bs0_lxyz = Bs0_vdiff3D.Mag();
		  Bs0_vDiff3D[0] = Bs0_vdiff3D.x(); Bs0_vDiff3D[1] = Bs0_vdiff3D.y(); Bs0_vDiff3D[2] = Bs0_vdiff3D.z() ;
		  Bs0_lxyzErr = sqrt(ROOT::Math::Similarity(Bs0_vDiff3D,Bs0_vXYe)) / Bs0_vdiff3D.Mag();
		
		  bs0CosAlphaPVCosAlpha->push_back( Bs0_cosAlpha ) ; bs0CosAlpha3DPVCosAlpha->push_back( Bs0_cosAlpha3D ) ;
		  bs0CTauPVCosAlpha->push_back( Bs0_ctau ) ; bs0CTauPVCosAlphaE->push_back( Bs0_ctauErr ) ;
		  bs0LxyPVCosAlpha->push_back( Bs0_lxy ) ; bs0LxyPVCosAlphaE->push_back( Bs0_lxyErr ) ;
		  bs0LxyzPVCosAlpha->push_back( Bs0_lxyz ) ; bs0LxyzPVCosAlphaE->push_back( Bs0_lxyzErr ) ;

		  /// Find the PV among the original offlinePV with the largest Bs0_cos(alpha) 3D
		  Vertex theCosAlpha3DV = thePrimaryVtx ; 
		  float maxCosAlpha3D = -1. ;
		  
		  for (VertexCollection::const_iterator itv = recVtxs->begin(), itvend = recVtxs->end(); itv != itvend; ++itv) {
		    Bs0_pvtx3D.SetXYZ(itv->position().x(), itv->position().y(), itv->position().z()) ;	
		    Bs0_vdiff3D = Bs0_vtx3D - Bs0_pvtx3D ;
		    float cosAlpha_temp3D = Bs0_vdiff3D.Dot(Bs0_pperp3D) / (Bs0_vdiff3D.Mag()*Bs0_pperp3D.Mag()) ;
		  
		    if ( cosAlpha_temp3D > maxCosAlpha3D ) {
		      maxCosAlpha3D = cosAlpha_temp3D ;    
		      theCosAlpha3DV = Vertex(*itv) ;
		    }
		  }
		  PriVtx_Bs0CosAlpha3D_n->push_back( recVtxs->size() ) ;
		  PriVtx_Bs0CosAlpha3D_X->push_back( theCosAlpha3DV.position().x() ) ;
		  PriVtx_Bs0CosAlpha3D_Y->push_back( theCosAlpha3DV.position().y() ) ;
		  PriVtx_Bs0CosAlpha3D_Z->push_back( theCosAlpha3DV.position().z() ) ;
		  PriVtx_Bs0CosAlpha3D_EX->push_back( theCosAlpha3DV.xError() ) ;
		  PriVtx_Bs0CosAlpha3D_EY->push_back( theCosAlpha3DV.yError() ) ;
		  PriVtx_Bs0CosAlpha3D_EZ->push_back( theCosAlpha3DV.zError() ) ;
		  PriVtx_Bs0CosAlpha3D_CL->push_back( ChiSquaredProbability((double)(theCosAlpha3DV.chi2()), (double)(theCosAlpha3DV.ndof())) ) ;
		  PriVtx_Bs0CosAlpha3D_Chi2->push_back( theCosAlpha3DV.chi2() ) ;
		  PriVtx_Bs0CosAlpha3D_tracks->push_back( theCosAlpha3DV.tracksSize() ) ;


		  X_v2e = theCosAlphaV.error(); /// X(4140)
                  X_vXYe = X_v1e.matrix() + X_v2e.matrix();
		  /// 2D
		  X_pvtx.SetXYZ(theCosAlphaV.position().x(), theCosAlphaV.position().y(), 0) ;
                  X_vdiff = X_vtx - X_pvtx ;
                  X_cosAlpha =  maxCosAlpha ;
                  X_lxy = X_vdiff.Perp();
                  X_vDiff[0] = X_vdiff.x(); X_vDiff[1] = X_vdiff.y(); X_vDiff[2] = 0 ; 
                  X_lxyErr = sqrt(ROOT::Math::Similarity(X_vDiff,X_vXYe)) / X_vdiff.Perp();
                  X_distXY = X_vdistXY.distance( Vertex(*XCand_vertex_fromMCFit), Vertex(theCosAlphaV) ) ;
                  X_ctau = X_distXY.value() * X_cosAlpha * XCand_fromMCFit->currentState().mass() / X_pperp.Perp();
                  X_ctauErr = sqrt(ROOT::Math::Similarity(X_v3pperp, X_vXYe)) * XCand_fromMCFit->currentState().mass() / (X_pperp.Perp2());
                  X_lxy = X_vdiff.Dot(X_pperp) / X_pperp.Mag() ;
		  /// 3D
		  X_pvtx3D.SetXYZ(theCosAlphaV.position().x(), theCosAlphaV.position().y(), theCosAlphaV.position().z());
                  X_vdiff3D = X_vtx3D - X_pvtx3D;
                  X_cosAlpha3D = X_vdiff3D.Dot(X_pperp3D)/( X_vdiff3D.Mag()*X_pperp3D.Mag() );
                  X_lxyz = X_vdiff3D.Mag();
                  X_vDiff3D[0] = X_vdiff3D.x(); X_vDiff3D[1] = X_vdiff3D.y(); X_vDiff3D[2] = X_vdiff3D.z() ;
                  X_lxyzErr = sqrt(ROOT::Math::Similarity(X_vDiff3D,X_vXYe)) / X_vdiff3D.Mag();

                  xCosAlphaPVCosAlpha->push_back( X_cosAlpha ) ; xCosAlpha3DPVCosAlpha->push_back( X_cosAlpha3D ) ;
                  xCTauPVCosAlpha->push_back( X_ctau ) ; xCTauPVCosAlphaE->push_back( X_ctauErr ) ;
                  xLxyPVCosAlpha->push_back( X_lxy ) ; xLxyPVCosAlphaE->push_back( X_lxyErr ) ;
                  xLxyzPVCosAlpha->push_back( X_lxyz ) ; xLxyzPVCosAlphaE->push_back( X_lxyzErr ) ;

		  /// Find the PV among the original offlinePV with the largest X_cos(alpha) 3D
		  //Vertex theCosAlpha3DV = thePrimaryVtx ;
                  //float maxCosAlpha3D = -1. ;

                  for (VertexCollection::const_iterator itv = recVtxs->begin(), itvend = recVtxs->end(); itv != itvend; ++itv) {
                    X_pvtx3D.SetXYZ(itv->position().x(), itv->position().y(), itv->position().z()) ;
                    X_vdiff3D = X_vtx3D - X_pvtx3D ;
                    float cosAlpha_temp3D = X_vdiff3D.Dot(X_pperp3D) / (X_vdiff3D.Mag()*X_pperp3D.Mag()) ;

                    if ( cosAlpha_temp3D > maxCosAlpha3D ) {
                      maxCosAlpha3D = cosAlpha_temp3D ;
                      theCosAlpha3DV = Vertex(*itv) ;
                    }
                  }

                  PriVtx_XCosAlpha3D_n->push_back( recVtxs->size() ) ;
                  PriVtx_XCosAlpha3D_X->push_back( theCosAlpha3DV.position().x() ) ;
                  PriVtx_XCosAlpha3D_Y->push_back( theCosAlpha3DV.position().y() ) ;
                  PriVtx_XCosAlpha3D_Z->push_back( theCosAlpha3DV.position().z() ) ;
                  PriVtx_XCosAlpha3D_EX->push_back( theCosAlpha3DV.xError() ) ;
                  PriVtx_XCosAlpha3D_EY->push_back( theCosAlpha3DV.yError() ) ;
                  PriVtx_XCosAlpha3D_EZ->push_back( theCosAlpha3DV.zError() ) ;
                  PriVtx_XCosAlpha3D_CL->push_back( ChiSquaredProbability((double)(theCosAlpha3DV.chi2()), (double)(theCosAlpha3DV.ndof())) ) ;
                  PriVtx_XCosAlpha3D_Chi2->push_back( theCosAlpha3DV.chi2() ) ;
                  PriVtx_XCosAlpha3D_tracks->push_back( theCosAlpha3DV.tracksSize() ) ;



		  ////////////////// Lifetime wrt PV with largest Bs0_cos(alpha) & X_cos(alpha) 3D candidate //////////////////  
		  Bs0_v2e = theCosAlpha3DV.error(); /// Bs0
		  Bs0_vXYe = Bs0_v1e.matrix() + Bs0_v2e.matrix();
		  /// 2D
		  Bs0_pvtx.SetXYZ(theCosAlpha3DV.position().x(), theCosAlpha3DV.position().y(), 0) ;	
		  Bs0_vdiff = Bs0_vtx - Bs0_pvtx ;
		  Bs0_cosAlpha = Bs0_vdiff.Dot(Bs0_pperp)/(Bs0_vdiff.Perp()*Bs0_pperp.Perp()); ;
		  Bs0_lxy = Bs0_vdiff.Perp();
		  Bs0_vDiff[0] = Bs0_vdiff.x(); Bs0_vDiff[1] = Bs0_vdiff.y(); Bs0_vDiff[2] = 0 ; // needed by Similarity method
		  Bs0_lxyErr = sqrt(ROOT::Math::Similarity(Bs0_vDiff,Bs0_vXYe)) / Bs0_vdiff.Perp();
		  Bs0_distXY = Bs0_vdistXY.distance( Vertex(*Bs0Cand_vertex_fromMCFit), Vertex(theCosAlphaV) ) ;
		  Bs0_ctau = Bs0_distXY.value() * Bs0_cosAlpha * Bs0Cand_fromMCFit->currentState().mass() / Bs0_pperp.Perp();
		  Bs0_ctauErr = sqrt(ROOT::Math::Similarity(Bs0_v3pperp,Bs0_vXYe)) * Bs0Cand_fromMCFit->currentState().mass() / (Bs0_pperp.Perp2());
		  Bs0_lxy = Bs0_vdiff.Dot(Bs0_pperp) / Bs0_pperp.Mag() ;
		  /// 3D
		  Bs0_pvtx3D.SetXYZ(theCosAlpha3DV.position().x(), theCosAlpha3DV.position().y(), theCosAlpha3DV.position().z()) ;	
		  Bs0_vdiff3D = Bs0_vtx3D - Bs0_pvtx3D ;
		  Bs0_cosAlpha3D =  maxCosAlpha3D ;
		  Bs0_lxyz = Bs0_vdiff3D.Mag();
		  Bs0_vDiff3D[0] = Bs0_vdiff3D.x(); Bs0_vDiff3D[1] = Bs0_vdiff3D.y(); Bs0_vDiff3D[2] = Bs0_vdiff3D.z() ;
		  Bs0_lxyzErr = sqrt(ROOT::Math::Similarity(Bs0_vDiff3D,Bs0_vXYe)) / Bs0_vdiff3D.Mag();

		  bs0CosAlphaPVCosAlpha3D->push_back( Bs0_cosAlpha ) ; bs0CosAlpha3DPVCosAlpha3D->push_back( Bs0_cosAlpha3D ) ;
		  bs0CTauPVCosAlpha3D->push_back( Bs0_ctau ) ; bs0CTauPVCosAlpha3DE->push_back( Bs0_ctauErr ) ;
		  bs0LxyPVCosAlpha3D->push_back( Bs0_lxy ) ; bs0LxyPVCosAlpha3DE->push_back( Bs0_lxyErr ) ;
		  bs0LxyzPVCosAlpha3D->push_back( Bs0_lxyz ) ; bs0LxyzPVCosAlpha3DE->push_back( Bs0_lxyzErr ) ;

		  /// Find the PV among the B0lessPV with the largest Bs0_cos(alpha)
		  Vertex theBs0LessCosAlphaV = thePrimaryVtx ;
		  maxCosAlpha = -1. ; 
		  
		  for (vector<TransientVertex>::iterator itv = Bs0_pvs.begin(), itvend = Bs0_pvs.end(); itv != itvend; ++itv) {
		    Bs0_pvtx.SetXYZ(itv->position().x(), itv->position().y(), 0) ;	
		    Bs0_vdiff = Bs0_vtx - Bs0_pvtx ;
		    float cosAlpha_temp = Bs0_vdiff.Dot(Bs0_pperp) / (Bs0_vdiff.Perp()*Bs0_pperp.Perp()) ; // Perp() == Mag() when z = 0
		  
 		    if ( cosAlpha_temp > maxCosAlpha ) {
		      maxCosAlpha = cosAlpha_temp ;    
		      theBs0LessCosAlphaV = Vertex(*itv) ;
		    }
		  }
	
		  PriVtxBs0Less_Bs0CosAlpha_n->push_back( Bs0_pvs.size() ) ;
		  PriVtxBs0Less_Bs0CosAlpha_X->push_back( theBs0LessCosAlphaV.position().x() ) ;
		  PriVtxBs0Less_Bs0CosAlpha_Y->push_back( theBs0LessCosAlphaV.position().y() ) ;
		  PriVtxBs0Less_Bs0CosAlpha_Z->push_back( theBs0LessCosAlphaV.position().z() ) ;
		  PriVtxBs0Less_Bs0CosAlpha_EX->push_back( theBs0LessCosAlphaV.xError() ) ;
		  PriVtxBs0Less_Bs0CosAlpha_EY->push_back( theBs0LessCosAlphaV.yError() ) ;
		  PriVtxBs0Less_Bs0CosAlpha_EZ->push_back( theBs0LessCosAlphaV.zError() ) ;
		  PriVtxBs0Less_Bs0CosAlpha_CL->push_back( ChiSquaredProbability((double)(theBs0LessCosAlphaV.chi2()), (double)(theBs0LessCosAlphaV.ndof())) ) ;
		  PriVtxBs0Less_Bs0CosAlpha_Chi2->push_back( theBs0LessCosAlphaV.chi2() ) ;
		  PriVtxBs0Less_Bs0CosAlpha_tracks->push_back( theBs0LessCosAlphaV.tracksSize() ) ;

		  X_v2e = theCosAlpha3DV.error(); /// X(4140)
                  X_vXYe = X_v1e.matrix() + X_v2e.matrix();
		  /// 2D
		  X_pvtx.SetXYZ(theCosAlpha3DV.position().x(), theCosAlpha3DV.position().y(), 0) ;
                  X_vdiff = X_vtx - X_pvtx ;
                  X_cosAlpha = X_vdiff.Dot(X_pperp)/(X_vdiff.Perp()*X_pperp.Perp()); ;
                  X_lxy = X_vdiff.Perp();
                  X_vDiff[0] = X_vdiff.x(); X_vDiff[1] = X_vdiff.y(); X_vDiff[2] = 0 ; 
                  X_lxyErr = sqrt(ROOT::Math::Similarity(X_vDiff,X_vXYe)) / X_vdiff.Perp();
                  X_distXY = X_vdistXY.distance( Vertex(*XCand_vertex_fromMCFit), Vertex(theCosAlphaV) ) ;
                  X_ctau = X_distXY.value() * X_cosAlpha * XCand_fromMCFit->currentState().mass() / X_pperp.Perp();
                  X_ctauErr = sqrt(ROOT::Math::Similarity(X_v3pperp,X_vXYe)) * XCand_fromMCFit->currentState().mass() / (X_pperp.Perp2());
                  X_lxy = X_vdiff.Dot(X_pperp) / X_pperp.Mag() ;
		  /// 3D
		  X_pvtx3D.SetXYZ(theCosAlpha3DV.position().x(), theCosAlpha3DV.position().y(), theCosAlpha3DV.position().z()) ;
                  X_vdiff3D = X_vtx3D - X_pvtx3D ;
                  X_cosAlpha3D =  maxCosAlpha3D ;
                  X_lxyz = X_vdiff3D.Mag();
                  X_vDiff3D[0] = X_vdiff3D.x(); X_vDiff3D[1] = X_vdiff3D.y(); X_vDiff3D[2] = X_vdiff3D.z() ;
                  X_lxyzErr = sqrt(ROOT::Math::Similarity(X_vDiff3D,Bs0_vXYe)) / X_vdiff3D.Mag();

                  xCosAlphaPVCosAlpha3D->push_back( X_cosAlpha ) ; xCosAlpha3DPVCosAlpha3D->push_back( X_cosAlpha3D ) ;
                  xCTauPVCosAlpha3D->push_back( X_ctau ) ; xCTauPVCosAlpha3DE->push_back( X_ctauErr ) ;
                  xLxyPVCosAlpha3D->push_back( X_lxy ) ; xLxyPVCosAlpha3DE->push_back( X_lxyErr ) ;
                  xLxyzPVCosAlpha3D->push_back( X_lxyz ) ; xLxyzPVCosAlpha3DE->push_back( X_lxyzErr ) ;
	


		  ////////////////// Lifetime wrt B0LessPV with largest Bs0_cos(alpha) candidate 
		  Bs0_v2e = theBs0LessCosAlphaV.error();
		  Bs0_vXYe = Bs0_v1e.matrix() + Bs0_v2e.matrix();
		  /// 2D
		  Bs0_pvtx.SetXYZ(theBs0LessCosAlphaV.position().x(), theBs0LessCosAlphaV.position().y(), 0) ;	
		  Bs0_vdiff = Bs0_vtx - Bs0_pvtx ;
		  Bs0_cosAlpha =  maxCosAlpha ;
		  Bs0_lxy = Bs0_vdiff.Perp();
		  Bs0_vDiff[0] = Bs0_vdiff.x(); Bs0_vDiff[1] = Bs0_vdiff.y(); Bs0_vDiff[2] = 0 ; // needed by Similarity method
		  Bs0_lxyErr = sqrt(ROOT::Math::Similarity(Bs0_vDiff,Bs0_vXYe)) / Bs0_vdiff.Perp() ;
		  Bs0_distXY = Bs0_vdistXY.distance( Vertex(*Bs0Cand_vertex_fromMCFit), Vertex(theBs0LessCosAlphaV) ) ;
		  Bs0_ctau = Bs0_distXY.value() * Bs0_cosAlpha * Bs0Cand_fromMCFit->currentState().mass() / Bs0_pperp.Perp();
		  Bs0_ctauErr = sqrt(ROOT::Math::Similarity(Bs0_v3pperp,Bs0_vXYe)) * Bs0Cand_fromMCFit->currentState().mass() / (Bs0_pperp.Perp2());
		  Bs0_lxy = Bs0_vdiff.Dot(Bs0_pperp) / Bs0_pperp.Mag() ;
		  /// 3D
		  Bs0_pvtx3D.SetXYZ(theBs0LessCosAlphaV.position().x(), theBs0LessCosAlphaV.position().y(), theBs0LessCosAlphaV.position().z());
		  Bs0_vdiff3D = Bs0_vtx3D - Bs0_pvtx3D;
		  Bs0_cosAlpha3D = Bs0_vdiff3D.Dot(Bs0_pperp3D)/( Bs0_vdiff3D.Mag()*Bs0_pperp3D.Mag() );
		  Bs0_lxyz = Bs0_vdiff3D.Mag();
		  Bs0_vDiff3D[0] = Bs0_vdiff3D.x(); Bs0_vDiff3D[1] = Bs0_vdiff3D.y(); Bs0_vDiff3D[2] = Bs0_vdiff3D.z() ;
		  Bs0_lxyzErr = sqrt(ROOT::Math::Similarity(Bs0_vDiff3D,Bs0_vXYe)) / Bs0_vdiff3D.Mag();
		
		  bs0CosAlphaBs0LessPVCosAlpha->push_back( Bs0_cosAlpha ) ; bs0CosAlpha3DBs0LessPVCosAlpha->push_back( Bs0_cosAlpha3D ) ;
		  bs0CTauBs0LessPVCosAlpha->push_back( Bs0_ctau ) ; bs0CTauBs0LessPVCosAlphaE->push_back( Bs0_ctauErr ) ;
		  bs0LxyBs0LessPVCosAlpha->push_back( Bs0_lxy ) ; bs0LxyBs0LessPVCosAlphaE->push_back( Bs0_lxyErr ) ;
		  bs0LxyzBs0LessPVCosAlpha->push_back( Bs0_lxyz ) ; bs0LxyzBs0LessPVCosAlphaE->push_back( Bs0_lxyzErr ) ;

		  /// Find the PV among the Bs0lessPV with the largest Bs0_cos(alpha) 3D
		  Vertex theBs0LessCosAlpha3DV = thePrimaryVtx ;
		  maxCosAlpha3D = -1. ; 
		  
		  for (vector<TransientVertex>::iterator itv = Bs0_pvs.begin(), itvend = Bs0_pvs.end(); itv != itvend; ++itv) {
		    Bs0_pvtx3D.SetXYZ(itv->position().x(), itv->position().y(), itv->position().z()) ;	
		    Bs0_vdiff3D = Bs0_vtx3D - Bs0_pvtx3D ;
		    float cosAlpha_temp3D = Bs0_vdiff3D.Dot(Bs0_pperp3D) / (Bs0_vdiff3D.Mag()*Bs0_pperp3D.Mag()) ;
		  
		    if ( cosAlpha_temp3D > maxCosAlpha3D ) {
		      maxCosAlpha3D = cosAlpha_temp3D ;   
		      theBs0LessCosAlpha3DV = Vertex(*itv) ;
		    }
		  }
		  	
		  PriVtxBs0Less_Bs0CosAlpha3D_n->push_back( Bs0_pvs.size() ) ;
		  PriVtxBs0Less_Bs0CosAlpha3D_X->push_back( theBs0LessCosAlpha3DV.position().x() ) ;
		  PriVtxBs0Less_Bs0CosAlpha3D_Y->push_back( theBs0LessCosAlpha3DV.position().y() ) ;
		  PriVtxBs0Less_Bs0CosAlpha3D_Z->push_back( theBs0LessCosAlpha3DV.position().z() ) ;
		  PriVtxBs0Less_Bs0CosAlpha3D_EX->push_back( theBs0LessCosAlpha3DV.xError() ) ;
		  PriVtxBs0Less_Bs0CosAlpha3D_EY->push_back( theBs0LessCosAlpha3DV.yError() ) ;
		  PriVtxBs0Less_Bs0CosAlpha3D_EZ->push_back( theBs0LessCosAlpha3DV.zError() ) ;
		  PriVtxBs0Less_Bs0CosAlpha3D_CL->push_back( ChiSquaredProbability((double)(theBs0LessCosAlpha3DV.chi2()), (double)(theBs0LessCosAlpha3DV.ndof())) ) ;
		  PriVtxBs0Less_Bs0CosAlpha3D_Chi2->push_back( theBs0LessCosAlpha3DV.chi2() ) ;
		  PriVtxBs0Less_Bs0CosAlpha3D_tracks->push_back( theBs0LessCosAlpha3DV.tracksSize() ) ;

		  ////////////////// Lifetime wrt Bs0LessPV with largest Bs0_cos(alpha) 3D candidate  
		  Bs0_v2e = theBs0LessCosAlpha3DV.error();
		  Bs0_vXYe = Bs0_v1e.matrix() + Bs0_v2e.matrix();
		  /// 2D
		  Bs0_pvtx.SetXYZ(theBs0LessCosAlpha3DV.position().x(), theBs0LessCosAlpha3DV.position().y(), 0) ;	
		  Bs0_vdiff = Bs0_vtx - Bs0_pvtx ;
		  Bs0_cosAlpha = Bs0_vdiff.Dot(Bs0_pperp)/(Bs0_vdiff.Perp()*Bs0_pperp.Perp());
		  Bs0_lxy = Bs0_vdiff.Perp();
		  Bs0_vDiff[0] = Bs0_vdiff.x(); Bs0_vDiff[1] = Bs0_vdiff.y(); Bs0_vDiff[2] = 0 ; // needed by Similarity method
		  Bs0_lxyErr = sqrt(ROOT::Math::Similarity(Bs0_vDiff,Bs0_vXYe)) / Bs0_vdiff.Perp();
		  Bs0_distXY = Bs0_vdistXY.distance( Vertex(*Bs0Cand_vertex_fromMCFit), Vertex(theCosAlphaV) ) ;
		  Bs0_ctau = Bs0_distXY.value() * Bs0_cosAlpha * Bs0Cand_fromMCFit->currentState().mass() / Bs0_pperp.Perp();
		  Bs0_ctauErr = sqrt(ROOT::Math::Similarity(Bs0_v3pperp,Bs0_vXYe)) * Bs0Cand_fromMCFit->currentState().mass() / (Bs0_pperp.Perp2());
		  Bs0_lxy = Bs0_vdiff.Dot(Bs0_pperp) / Bs0_pperp.Mag() ;
		  /// 3D
		  Bs0_pvtx3D.SetXYZ(theBs0LessCosAlpha3DV.position().x(), theBs0LessCosAlpha3DV.position().y(), theBs0LessCosAlpha3DV.position().z()) ;	
		  Bs0_vdiff3D = Bs0_vtx3D - Bs0_pvtx3D ;
		  Bs0_cosAlpha3D =  maxCosAlpha3D ;
		  Bs0_lxyz = Bs0_vdiff3D.Mag();
		  Bs0_vDiff3D[0] = Bs0_vdiff3D.x(); Bs0_vDiff3D[1] = Bs0_vdiff3D.y(); Bs0_vDiff3D[2] = Bs0_vdiff3D.z() ;
		  Bs0_lxyzErr = sqrt(ROOT::Math::Similarity(Bs0_vDiff3D,Bs0_vXYe)) / Bs0_vdiff3D.Mag();
		  Bs0_lxy = Bs0_vdiff3D.Dot(Bs0_pperp) / Bs0_pperp.Mag() ;
		
		  bs0CosAlphaBs0LessPVCosAlpha3D->push_back( Bs0_cosAlpha ) ;  bs0CosAlpha3DBs0LessPVCosAlpha3D->push_back( Bs0_cosAlpha3D ) ;
		  bs0CTauBs0LessPVCosAlpha3D->push_back( Bs0_ctau ) ; bs0CTauBs0LessPVCosAlpha3DE->push_back( Bs0_ctauErr ) ;
		  bs0LxyBs0LessPVCosAlpha3D->push_back( Bs0_lxy ) ; bs0LxyBs0LessPVCosAlpha3DE->push_back( Bs0_lxyErr ) ;
		  bs0LxyzBs0LessPVCosAlpha3D->push_back( Bs0_lxyz ) ; bs0LxyzBs0LessPVCosAlpha3DE->push_back( Bs0_lxyzErr ) ;

		
		  Vertex theOtherV = thePrimaryVtx; 
		  	  
		  if (resolveAmbiguity_) {
		    float minDz = 999999. ;
		    if (!addBs0lessPrimaryVertex_) {
		      for (VertexCollection::const_iterator itv = recVtxs->begin(), itvend = recVtxs->end(); itv != itvend; ++itv)
			{
			  float deltaZ = fabs((*Bs0Cand_vertex_fromMCFit).position().z() - itv->position().z()) ;
			  if ( deltaZ < minDz ) {
			    minDz = deltaZ;    
			    thePrimaryVtx = Vertex(*itv);
			    theOtherV = thePrimaryVtx;
			  }
			}
		    } else {
		      for (vector<TransientVertex>::iterator itv2 = Bs0_pvs.begin(), itvend2 = Bs0_pvs.end(); itv2 != itvend2; ++itv2)
			{
			  float deltaZ = fabs((*Bs0Cand_vertex_fromMCFit).position().z() - itv2->position().z()) ;
			  if ( deltaZ < minDz ) {
			    minDz = deltaZ;    
			    Vertex Bs0LessPV = Vertex(*itv2); 
			    thePrimaryVtx = Bs0LessPV;
			    theOtherV = Bs0LessPV;
			  }
			}
		    }
		  } 
		  
		  Vertex TheOtherVertex3D = thePrimaryVtx; 
		  if (Debug_) cout<<" choose PV ="<< endl;
		  Int_t theBs0CorrPV_multiplicity = -1 ; 
		  if (resolveAmbiguity_) {
		    float minDz = 999999.;
		    if (!addBs0lessPrimaryVertex_) {
		      theBs0CorrPV_multiplicity = recVtxs->size() ;
		      for (VertexCollection::const_iterator itv = recVtxs->begin(), itvend = recVtxs->end(); itv != itvend; ++itv) {
			float deltaZ = fabs((*Bs0Cand_vertex_fromMCFit).position().z() - itv->position().z()) ;
			if ( deltaZ < minDz ) {
			  minDz = deltaZ;    
			  TheOtherVertex3D = Vertex(*itv);
			}
		      }
		    } else {
		      theBs0CorrPV_multiplicity = Bs0_pvs.size() ;
		      for (vector<TransientVertex>::iterator itv2 = Bs0_pvs.begin(), itvend2 = Bs0_pvs.end(); itv2 != itvend2; ++itv2) {
			VertexDistance3D a3d;
			float deltaZ   = a3d.distance(Vertex(*itv2), Vertex(*Bs0Cand_vertex_fromMCFit)).value();
			if ( deltaZ < minDz ) {
			  minDz = deltaZ;    
			  Vertex XLessPV = Vertex(*itv2); 
			  TheOtherVertex3D = XLessPV;
			  //cout<<" z(X) - z(vtx) min="<<minDz<<endl; 
			}
		      
		      }
		    }
		  } 
		 
		  PriVtxBs0Corr_n->push_back( theBs0CorrPV_multiplicity ) ; 
		  PriVtxBs0Corr_X->push_back( thePrimaryVtx.position().x() ) ; 
		  PriVtxBs0Corr_Y->push_back( thePrimaryVtx.position().y() ) ; 
		  PriVtxBs0Corr_Z->push_back( thePrimaryVtx.position().z() ) ; 
		  PriVtxBs0Corr_EX->push_back( thePrimaryVtx.xError() ) ; 
		  PriVtxBs0Corr_EY->push_back( thePrimaryVtx.yError() ) ; 
		  PriVtxBs0Corr_EZ->push_back( thePrimaryVtx.zError() ) ; 
		  PriVtxBs0Corr_CL->push_back( ChiSquaredProbability( (double)(thePrimaryVtx.chi2()), (double)(thePrimaryVtx.ndof())) ); 
		  PriVtxBs0Corr_Chi2->push_back( thePrimaryVtx.chi2() ) ; 
		  PriVtxBs0Corr_tracks->push_back( thePrimaryVtx.tracksSize() ) ; 
		  			
		
				
		  ////////////////// Lifetime wrt PV with smaller longitudinal X impact parameter for Bs0 & X(4140) //////////////////  
		  Bs0_pvtx.SetXYZ(theOtherV.position().x(), theOtherV.position().y(), 0); /// Bs0
		  Bs0_vdiff = Bs0_vtx - Bs0_pvtx;
		  Bs0_cosAlpha = Bs0_vdiff.Dot(Bs0_pperp) / (Bs0_vdiff.Perp()*Bs0_pperp.Perp());
		  Bs0_distXY = Bs0_vdistXY.distance(Vertex(*Bs0Cand_vertex_fromMCFit), Vertex(theOtherV));
		  double Bs0_ctauPVX = Bs0_distXY.value() * Bs0_cosAlpha * Bs0Cand_fromMCFit->currentState().mass() / Bs0_pperp.Perp();
		  GlobalError Bs0_v1eX = (Vertex(*Bs0Cand_vertex_fromMCFit)).error();
		  GlobalError Bs0_v2eX = theOtherV.error();
		  AlgebraicSymMatrix33 Bs0_vXYeX = Bs0_v1eX.matrix() + Bs0_v2eX.matrix();
		  double ctauErrPVX = sqrt(ROOT::Math::Similarity(Bs0_v3pperp,Bs0_vXYeX)) * Bs0Cand_fromMCFit->currentState().mass() / (Bs0_pperp.Perp2());
		  float lxyPVX = Bs0_vdiff.Dot(Bs0_pperp) / Bs0_pperp.Mag() ;
		  float lxyzPVX = Bs0_vdiff3D.Dot(Bs0_pperp3D) / Bs0_pperp3D.Mag() ;
		  bs0CosAlphaPVX->push_back(Bs0_cosAlpha);
		  bs0CTauPVX->push_back(Bs0_ctauPVX); bs0CTauPVXE->push_back(ctauErrPVX);
		  bs0LxyPVX->push_back(lxyPVX);
		  bs0LxyzPVX->push_back(lxyzPVX);
		  VertexDistance3D a3d; 
		  float Dist3DPV     = a3d.distance(TheOtherVertex3D, Vertex(*Bs0Cand_vertex_fromMCFit)).value() ;
		  float Dist3DPV_err = a3d.distance(TheOtherVertex3D, Vertex(*Bs0Cand_vertex_fromMCFit)).error() ;
		  bs0CTauPVX_3D->push_back(Dist3DPV);
		  bs0CTauPVX_3D_err->push_back(Dist3DPV_err);
		  //cout << Dist3DPV << " " << Dist3DPV_err << endl; 
		  Bs0_MuMuIdx->push_back(nMuMu-1);
		  Bs0_ka1Idx->push_back(std::distance(theKaonRefittedPATTrackHandle->begin(), Track1)); 
		  Bs0_ka2Idx->push_back(std::distance(theKaonRefittedPATTrackHandle->begin(), Track2)); 
		  nBs0++;
		  bs0Daughters.clear();
		  
		  X_pvtx.SetXYZ(theOtherV.position().x(), theOtherV.position().y(), 0); /// X(4140)
                  X_vdiff = X_vtx - X_pvtx;
                  X_cosAlpha = X_vdiff.Dot(X_pperp) / (X_vdiff.Perp()*X_pperp.Perp());
                  X_distXY = X_vdistXY.distance(Vertex(*XCand_vertex_fromMCFit), Vertex(theOtherV));
                  double X_ctauPVX = X_distXY.value() * X_cosAlpha * XCand_fromMCFit->currentState().mass() / X_pperp.Perp();
                  GlobalError X_v1eX = (Vertex(*XCand_vertex_fromMCFit)).error();
                  GlobalError X_v2eX = theOtherV.error();
                  AlgebraicSymMatrix33 X_vXYeX = X_v1eX.matrix() + X_v2eX.matrix();
                  double X_ctauErrPVX = sqrt(ROOT::Math::Similarity(X_v3pperp,X_vXYeX)) * XCand_fromMCFit->currentState().mass() / (X_pperp.Perp2());
                  float X_lxyPVX = X_vdiff.Dot(X_pperp) / X_pperp.Mag() ;
                  float X_lxyzPVX = X_vdiff3D.Dot(X_pperp3D) / X_pperp3D.Mag() ;
                  xCosAlphaPVX->push_back(X_cosAlpha);
                  xCTauPVX->push_back(X_ctauPVX); xCTauPVXE->push_back(X_ctauErrPVX);
                  xLxyPVX->push_back(X_lxyPVX);
                  xLxyzPVX->push_back(X_lxyzPVX);
                  //VertexDistance3D a3d;
                  //float Dist3DPV     = a3d.distance(TheOtherVertex3D, Vertex(*XCand_vertex_fromMCFit)).value() ;
                  //float Dist3DPV_err = a3d.distance(TheOtherVertex3D, Vertex(*XCand_vertex_fromMCFit)).error() ;
                  xCTauPVX_3D->push_back(Dist3DPV);
                  xCTauPVX_3D_err->push_back(Dist3DPV_err);
		  //cout << Dist3DPV << " " << Dist3DPV_err << endl;
		  X_MuMuIdx->push_back(nMuMu-1);
                  X_ka1Idx->push_back(std::distance(theKaonRefittedPATTrackHandle->begin(), Track1)); 
                  X_ka2Idx->push_back(std::distance(theKaonRefittedPATTrackHandle->begin(), Track2)); 
                  nX++;
                  xDaughters.clear();

		  ////////////////// flag for checking the Kaons from PV or not PV //////////////////
		  /// flag for kaon1
		  vector<TransientTrack> vertexTracksKaon1;
                  pat::GenericParticle patKaon1 = *Track1;
                  //cout << "\nthePrimaryVtx.tracksSize() = " << thePrimaryVtx.tracksSize() << endl;
                  //cout << "thePrimaryVtx.nTracks() = " << thePrimaryVtx.nTracks() << endl;
		  for ( std::vector<TrackBaseRef >::const_iterator iTrack = thePrimaryVtx.tracks_begin(); iTrack != thePrimaryVtx.tracks_end(); ++iTrack) {

		    TrackRef trackRefkaon1 = iTrack->castTo<TrackRef>();
		    //cout << "\ntrackRefkaon1 = " << trackRefkaon1 << endl;
		    //cout <<"before match" ;
		    //if ( (patKaon1.track() == trackRefkaon1) ) {
		    if ( (Track1->track().key() == trackRefkaon1.key()) ) {
		      cout << "\ninside match" << endl;                 
		      TransientTrack kaon1TT(trackRefkaon1, &(*bFieldHandle) );
		      vertexTracksKaon1.push_back(kaon1TT);
		    }      
		  }
                  //cout << "\nvertexTracksKaon1.size() = " << vertexTracksKaon1.size() << endl; 
		  if (vertexTracksKaon1.size()==0)
		    Kaon1FromPV->push_back(false);
		  else
		    Kaon1FromPV->push_back(true);
		   
		  /// flag for kaon2
		  vector<TransientTrack> vertexTracksKaon2;
	          pat::GenericParticle patKaon2 = *Track2;
                  //cout << "\nthePrimaryVtx.tracksSize() = " << thePrimaryVtx.tracksSize() << endl;
                  //cout << "thePrimaryVtx.nTracks() = " << thePrimaryVtx.nTracks() << endl;
		  for ( std::vector<TrackBaseRef >::const_iterator iTrack = thePrimaryVtx.tracks_begin(); iTrack != thePrimaryVtx.tracks_end(); ++iTrack) {
	
		    TrackRef trackRefkaon2 = iTrack->castTo<TrackRef>();
                    //cout << "\ntrackRefkaon2 = " << trackRefkaon2 << endl;
                    //cout <<"before match" ;
                    //if ( (patKaon2.track() == trackRefkaon2) ) {
		    if (  (Track2->track().key() == trackRefkaon2.key()) ) {
		      TransientTrack kaon2TT(trackRefkaon2, &(*bFieldHandle) );
		      vertexTracksKaon2.push_back(kaon2TT);
		    }
		  }
		  //cout << "\nvertexTracksKaon1.size() = " << vertexTracksKaon1.size() << endl; 
		  if (vertexTracksKaon2.size()==0)
                    Kaon2FromPV->push_back(false);
		  else 
		    Kaon2FromPV->push_back(true);
		  
		} // 2nd loop over track (look for k2)
	    } // 1st loop over track (look for k1)	 
	  } // 2nd loop over muons (look for mu-)
	} //first loop over muons (look for mu+)
      } // if (thePATMuonHandle->size() >= 2  && hasRequestedTrigger) {
  } // if (doMC || doData) 
//} 
  // AT THE END OF THE EVENT fill the tree and clear the vectors
  // ===========================================================
 
  if (nX > 0) 
    X_One_Tree_->Fill() ;

  /// trigger stuff
  trigRes->clear(); trigNames->clear(); L1TT->clear(); MatchTriggerNames->clear(); 
  /// event numbers
  runNum = 0; evtNum = 0; lumiNum = 0;
  /// counters for Bs0 & X(4140)
  nMu = 0; nMuMu = 0; nBs0 = 0; nKK = 0; 
  nBs0_pre0 = 0; nBs0_pre1 = 0; nBs0_pre2 = 0; nBs0_pre3 = 0; nBs0_pre4 = 0; nBs0_pre5 = 0; nBs0_pre6 = 0; nBs0_pre7 = 0; nBs0_pre8 = 0; nBs0_pre9 = 0; nBs0_pre10 = 0; nBs0_pre11 = 0; nBs0_pre12 = 0; nBs0_pre13 = 0; nBs0_pre14 = 0; nBs0_pre15 = 0; nBs0_pre16 = 0;
  nX = 0;
  nX_pre0 = 0; nX_pre1 = 0; nX_pre2 = 0; nX_pre3 = 0; nX_pre4 = 0; nX_pre5 = 0; nX_pre6 = 0; nX_pre7 = 0; nX_pre8 = 0; nX_pre9 = 0; nX_pre10 = 0; nX_pre11 = 0; nX_pre12 = 0; nX_pre13 = 0; nX_pre14 = 0; nX_pre15 = 0; nX_pre16 = 0;
  /// indices
  mu1Idx->clear(); mu2Idx->clear(); 
  ka1Idx->clear(); ka2Idx->clear();
  Bs0_MuMuIdx->clear(); Bs0_ka1Idx->clear(); Bs0_ka2Idx->clear(); 
  X_MuMuIdx->clear(); X_ka1Idx->clear(); X_ka2Idx->clear();

  /// MC Analysis
  if (doMC) {
    // Gen Primary Vertex
    n_genEvtVtx = 0;
    genEvtVtx_X->clear(); genEvtVtx_Y->clear(); genEvtVtx_Z->clear(); 
    genEvtVtx_particles->clear();
    n_Bs0Ancestors->clear(); 
    nMCAll = 0, nMCBs0 = 0; //nMCB0Vtx = 0; 
    // Gen Primary Vertex
    PriVtxGen_X->clear(); PriVtxGen_Y->clear(); PriVtxGen_Z->clear(); 
    PriVtxGen_EX->clear(); PriVtxGen_EY->clear(); PriVtxGen_EZ->clear();  
    PriVtxGen_Chi2->clear(); PriVtxGen_CL->clear(); PriVtxGen_Ndof->clear();
    PriVtxGen_tracks->clear();
    
    MCPdgIdAll->clear(); MCDanNumAll->clear();
    MCJPsiPx->clear(); MCJPsiPy->clear(); MCJPsiPz->clear(); 
    MCmupPx->clear(); MCmupPy->clear(); MCmupPz->clear(); 
    MCmumPx->clear(); MCmumPy->clear(); MCmumPz->clear();
    MCPhiPx->clear(); MCPhiPy->clear(); MCPhiPz->clear();
    MCkpPx->clear(); MCkpPy->clear(); MCkpPz->clear();
    MCkmPx->clear(); MCkmPy->clear(); MCkmPz->clear(); 
    //MCpionPx->clear(); MCpionPy->clear(); MCpionPz->clear(); 
    //MCkaonPx->clear(); MCkaonPy->clear(); MCkaonPz->clear(); 
    //MCpionCh->clear(); MCkaonCh->clear();
    MCPx->clear(); MCPy->clear(); MCPz->clear();
  }
  if (Debug_) cout <<"after MC stuff clear" <<endl ;
  /// Primary Vertex	
  priVtx_n = 0;
  priVtx_X = 0; priVtx_Y = 0; priVtx_Z = 0 ; 
  priVtx_XE = 0; priVtx_YE = 0; priVtx_ZE = 0 ; 
  priVtx_NormChi2 = 0; priVtx_Chi2 = 0; priVtx_CL = 0; priVtx_tracks = 0; priVtx_tracksPtSq = 0 ;
  /// MuMu cand & KK cand
  MuMuMass->clear(); MuMuVtx_CL->clear(); MuMuVtx_Chi2->clear(); 
  MuMuPx->clear(); MuMuPy->clear(); MuMuPz->clear();
  MuMuDecayVtx_X->clear(); MuMuDecayVtx_Y->clear(); MuMuDecayVtx_Z->clear();
  MuMuDecayVtx_XE->clear(); MuMuDecayVtx_YE->clear(); MuMuDecayVtx_ZE->clear();
  MuMuMuonTrigMatch->clear();
  KKMass->clear(); KKPx->clear(); KKPy->clear(); KKPz->clear();
  KKVtx_CL->clear(); KKVtx_Chi2->clear();
  KKDecayVtx_X->clear(); KKDecayVtx_Y->clear(); KKDecayVtx_Z->clear();
  KKDecayVtx_XE->clear(); KKDecayVtx_YE->clear(); KKDecayVtx_ZE->clear();
  /// muons from JPsi (MuMu) fit & kaons from Phi (KK) fit
  mu1_MuMu_Px->clear(); mu1_MuMu_Py->clear(); mu1_MuMu_Pz->clear(); mu1_MuMu_Chi2->clear(); mu1_MuMu_NDF->clear();
  mu2_MuMu_Px->clear(); mu2_MuMu_Py->clear(); mu2_MuMu_Pz->clear(); mu2_MuMu_Chi2->clear(); mu2_MuMu_NDF->clear();
  MuMuType->clear();
  ka1_KK_Px->clear(); ka1_KK_Py->clear(); ka1_KK_Pz->clear(); ka1_KK_Chi2->clear(); ka1_KK_NDF->clear();
  ka2_KK_Px->clear(); ka2_KK_Py->clear();  ka2_KK_Pz->clear(); ka2_KK_Chi2->clear(); ka2_KK_NDF->clear();
  DRMuMuK1->clear(); DRMuMuK2->clear(); DRbs0K1->clear(); DRbs0K2->clear(); DRxK1->clear(); DRxK2->clear();
  /// Primary Vertex with "MuMu correction"
  PriVtxMuMuCorr_n->clear();
  PriVtxMuMuCorr_X->clear(); PriVtxMuMuCorr_Y->clear(); PriVtxMuMuCorr_Z->clear(); 
  PriVtxMuMuCorr_EX->clear(); PriVtxMuMuCorr_EY->clear(); PriVtxMuMuCorr_EZ->clear();  
  PriVtxMuMuCorr_Chi2->clear(); PriVtxMuMuCorr_CL->clear(); PriVtxMuMuCorr_tracks->clear();
  nTrk->clear();
  /// Bs0 cand & X(4140) cand	
  bs0Mass->clear(); bs0Vtx_CL->clear(); bs0Vtx_Chi2->clear(); 
  bs0Px->clear(); bs0Py->clear(); bs0Pz->clear(); 
  bs0PxE->clear(); bs0PyE->clear(); bs0PzE->clear();
  bs0DecayVtx_X->clear(); bs0DecayVtx_Y->clear(); bs0DecayVtx_Z->clear(); 
  bs0DecayVtx_XE->clear(); bs0DecayVtx_YE->clear(); bs0DecayVtx_ZE->clear(); 
  xMass->clear(); xVtx_CL->clear(); xVtx_Chi2->clear();
  xPx->clear(); xPy->clear(); xPz->clear();
  xPxE->clear(); xPyE->clear(); xPzE->clear();
  xDecayVtx_X->clear(); xDecayVtx_Y->clear(); xDecayVtx_Z->clear();
  xDecayVtx_XE->clear(); xDecayVtx_YE->clear(); xDecayVtx_ZE->clear();
  /// muons and tracks after Bs0 cand fit & X(4140) cand fit 
  Bs0_mu1Px_MuMuKK->clear(); Bs0_mu1Py_MuMuKK->clear(); Bs0_mu1Pz_MuMuKK->clear(); Bs0_mu1E_MuMuKK->clear();
  Bs0_mu2Px_MuMuKK->clear(); Bs0_mu2Py_MuMuKK->clear(); Bs0_mu2Pz_MuMuKK->clear(); Bs0_mu2E_MuMuKK->clear();
  Bs0_k1Px_MuMuKK->clear(); Bs0_k1Py_MuMuKK->clear(); Bs0_k1Pz_MuMuKK->clear(); Bs0_k1E_MuMuKK->clear(); 
  Bs0_kaon1_nsigdedx->clear(); Bs0_kaon1_dedx->clear(); Bs0_kaon1_dedxMass->clear(); Bs0_kaon1_theo->clear(); Bs0_kaon1_sigma->clear(); 
  Bs0_kaon1_dedx_byHits->clear(); Bs0_kaon1_dedxErr_byHits->clear(); Bs0_kaon1_saturMeas_byHits->clear(); Bs0_kaon1_Meas_byHits->clear(); 
  Bs0_k2Px_MuMuKK->clear(); Bs0_k2Py_MuMuKK->clear(); Bs0_k2Pz_MuMuKK->clear(); Bs0_k2E_MuMuKK->clear(); 
  Bs0_kaon2_nsigdedx->clear(); Bs0_kaon2_dedx->clear(); Bs0_kaon2_dedxMass->clear(); Bs0_kaon2_theo->clear(); Bs0_kaon2_sigma->clear(); 
  Bs0_kaon2_dedx_byHits->clear(); Bs0_kaon2_dedxErr_byHits->clear(); Bs0_kaon2_saturMeas_byHits->clear(); Bs0_kaon2_Meas_byHits->clear(); 
  X_mu1Px_MuMuKK->clear(); X_mu1Py_MuMuKK->clear(); X_mu1Pz_MuMuKK->clear(); X_mu1E_MuMuKK->clear();
  X_mu2Px_MuMuKK->clear(); X_mu2Py_MuMuKK->clear(); X_mu2Pz_MuMuKK->clear(); X_mu2E_MuMuKK->clear();
  X_k1Px_MuMuKK->clear(); X_k1Py_MuMuKK->clear(); X_k1Pz_MuMuKK->clear(); X_k1E_MuMuKK->clear(); 
  X_kaon1_nsigdedx->clear(); X_kaon1_dedx->clear(); X_kaon1_dedxMass->clear(); X_kaon1_theo->clear(); X_kaon1_sigma->clear();
  X_kaon1_dedx_byHits->clear(); X_kaon1_dedxErr_byHits->clear(); X_kaon1_saturMeas_byHits->clear(); X_kaon1_Meas_byHits->clear(); 
  X_k2Px_MuMuKK->clear(); X_k2Py_MuMuKK->clear(); X_k2Pz_MuMuKK->clear(); X_k2E_MuMuKK->clear(); 
  X_kaon2_nsigdedx->clear(); X_kaon2_dedx->clear(); X_kaon2_dedxMass->clear(); X_kaon2_theo->clear(); X_kaon2_sigma->clear(); 
  X_kaon2_dedx_byHits->clear(); X_kaon2_dedxErr_byHits->clear(); X_kaon2_saturMeas_byHits->clear(); X_kaon2_Meas_byHits->clear();
  /// Primary Vertex with "Bs0 correction" 
  PriVtxBs0Less_n->clear();
  PriVtxBs0Less_X->clear(); PriVtxBs0Less_Y->clear(); PriVtxBs0Less_Z->clear(); 
  PriVtxBs0Less_EX->clear(); PriVtxBs0Less_EY->clear(); PriVtxBs0Less_EZ->clear();  
  PriVtxBs0Less_Chi2->clear(); PriVtxBs0Less_CL->clear(); PriVtxBs0Less_tracks->clear();
  /// Primary Vertex with largest Bs0_cos(alpha) & largest X(4140)_cos(alpha)
  Bs0LessPV_tracksPtSq->clear(); Bs0LessPV_4tracksPtSq->clear();
  PriVtx_Bs0CosAlpha_n->clear();
  PriVtx_Bs0CosAlpha_X->clear(); PriVtx_Bs0CosAlpha_Y->clear(); PriVtx_Bs0CosAlpha_Z->clear(); 
  PriVtx_Bs0CosAlpha_EX->clear(); PriVtx_Bs0CosAlpha_EY->clear(); PriVtx_Bs0CosAlpha_EZ->clear();  
  PriVtx_Bs0CosAlpha_Chi2->clear(); PriVtx_Bs0CosAlpha_CL->clear(); PriVtx_Bs0CosAlpha_tracks->clear();
  PriVtxBs0Less_Bs0CosAlpha_n->clear();
  PriVtxBs0Less_Bs0CosAlpha_X->clear(); PriVtxBs0Less_Bs0CosAlpha_Y->clear(); PriVtxBs0Less_Bs0CosAlpha_Z->clear(); 
  PriVtxBs0Less_Bs0CosAlpha_EX->clear(); PriVtxBs0Less_Bs0CosAlpha_EY->clear(); PriVtxBs0Less_Bs0CosAlpha_EZ->clear();  
  PriVtxBs0Less_Bs0CosAlpha_Chi2->clear(); PriVtxBs0Less_Bs0CosAlpha_CL->clear(); PriVtxBs0Less_Bs0CosAlpha_tracks->clear();
  PriVtx_XCosAlpha_n->clear();
  PriVtx_XCosAlpha_X->clear(); PriVtx_XCosAlpha_Y->clear(); PriVtx_XCosAlpha_Z->clear();
  PriVtx_XCosAlpha_EX->clear(); PriVtx_XCosAlpha_EY->clear(); PriVtx_XCosAlpha_EZ->clear();
  PriVtx_XCosAlpha_Chi2->clear(); PriVtx_XCosAlpha_CL->clear(); PriVtx_XCosAlpha_tracks->clear();
  /// Primary Vertex with "Bs0 correction" & "X(4140) correction" 
  PriVtxBs0Corr_n->clear();
  PriVtxBs0Corr_X->clear(); PriVtxBs0Corr_Y->clear(); PriVtxBs0Corr_Z->clear(); 
  PriVtxBs0Corr_EX->clear(); PriVtxBs0Corr_EY->clear(); PriVtxBs0Corr_EZ->clear();  
  PriVtxBs0Corr_Chi2->clear(); PriVtxBs0Corr_CL->clear(); PriVtxBs0Corr_tracks->clear();
  PriVtxXCorr_n->clear();
  PriVtxXCorr_X->clear(); PriVtxXCorr_Y->clear(); PriVtxXCorr_Z->clear();
  PriVtxXCorr_EX->clear(); PriVtxXCorr_EY->clear(); PriVtxXCorr_EZ->clear();
  PriVtxXCorr_Chi2->clear(); PriVtxXCorr_CL->clear(); PriVtxXCorr_tracks->clear();
  /// Lifetime variables for Bs0 & X(4140)
  bs0CosAlphaBS->clear(); bs0CosAlpha3DBS->clear(); bs0CTauBS->clear(); bs0CTauBSE->clear(); bs0LxyBS->clear(); bs0LxyBSE->clear(); bs0LxyzBS->clear(); bs0LxyzBSE->clear(); 
  bs0CosAlphaPV->clear(); bs0CosAlpha3DPV->clear(); bs0CTauPV->clear(); bs0CTauPVE->clear(); bs0LxyPV->clear(); bs0LxyPVE->clear(); bs0LxyzPV->clear(); bs0LxyzPVE->clear();
  bs0CosAlphaPVCosAlpha->clear(); bs0CosAlpha3DPVCosAlpha->clear(); bs0CTauPVCosAlpha->clear(); bs0CTauPVCosAlphaE->clear(); bs0LxyPVCosAlpha->clear(); bs0LxyPVCosAlphaE->clear(); bs0LxyzPVCosAlpha->clear(); bs0LxyzPVCosAlphaE->clear(); 
  bs0CosAlphaPVCosAlpha3D->clear(); bs0CosAlpha3DPVCosAlpha3D->clear(); bs0CTauPVCosAlpha3D->clear(); bs0CTauPVCosAlpha3DE->clear(); bs0LxyPVCosAlpha3D->clear(); bs0LxyPVCosAlpha3DE->clear(); bs0LxyzPVCosAlpha3D->clear(); bs0LxyzPVCosAlpha3DE->clear();
  bs0CosAlphaBs0LessPV->clear(); bs0CosAlpha3DBs0LessPV->clear(); bs0CTauBs0LessPV->clear() ; bs0CTauBs0LessPVE->clear() ; bs0LxyBs0LessPV->clear() ; bs0LxyBs0LessPVE->clear() ; bs0LxyzBs0LessPV->clear() ; bs0LxyzBs0LessPVE->clear() ;
  bs0CosAlphaBs0LessPVCosAlpha->clear(); bs0CosAlpha3DBs0LessPVCosAlpha->clear(); bs0CTauBs0LessPVCosAlpha->clear() ; bs0CTauBs0LessPVCosAlphaE->clear() ; bs0LxyBs0LessPVCosAlpha->clear() ; bs0LxyBs0LessPVCosAlphaE->clear() ; bs0LxyzBs0LessPVCosAlpha->clear() ; bs0LxyzBs0LessPVCosAlphaE->clear() ;
  bs0CosAlphaBs0LessPVCosAlpha3D->clear(); bs0CosAlpha3DBs0LessPVCosAlpha3D->clear(); bs0CTauBs0LessPVCosAlpha3D->clear() ; bs0CTauBs0LessPVCosAlpha3DE->clear() ; bs0LxyBs0LessPVCosAlpha3D->clear() ; bs0LxyBs0LessPVCosAlpha3DE->clear() ; bs0LxyzBs0LessPVCosAlpha3D->clear() ; bs0LxyzBs0LessPVCosAlpha3DE->clear() ;
  bs0CosAlphaPVX->clear(); bs0CTauPVX->clear(); bs0CTauPVXE->clear(); bs0LxyPVX->clear(); bs0LxyzPVX->clear(); 
  bs0CTauPVX_3D->clear(); bs0CTauPVX_3D_err->clear();
  xCosAlphaBS->clear(); xCosAlpha3DBS->clear(); xCTauBS->clear(); xCTauBSE->clear(); xLxyBS->clear(); xLxyBSE->clear(); xLxyzBS->clear(); xLxyzBSE->clear();
  xCosAlphaPV->clear(); xCosAlpha3DPV->clear(); xCTauPV->clear(); xCTauPVE->clear(); xLxyPV->clear(); xLxyPVE->clear(); xLxyzPV->clear(); xLxyzPVE->clear(); xCosAlphaPVCosAlpha->clear(); xCosAlpha3DPVCosAlpha->clear(); xCTauPVCosAlpha->clear(); xCTauPVCosAlphaE->clear(); xLxyPVCosAlpha->clear(); xLxyPVCosAlphaE->clear(); xLxyzPVCosAlpha->clear(); xLxyzPVCosAlphaE->clear();
  xCosAlphaPVCosAlpha3D->clear(); xCosAlpha3DPVCosAlpha3D->clear(); xCTauPVCosAlpha3D->clear(); xCTauPVCosAlpha3DE->clear(); xLxyPVCosAlpha3D->clear(); xLxyPVCosAlpha3DE->clear(); xLxyzPVCosAlpha3D->clear(); xLxyzPVCosAlpha3DE->clear();  
  xCosAlphaPVX->clear(); xCTauPVX->clear(); xCTauPVXE->clear(); xLxyPVX->clear(); xLxyzPVX->clear();
  xCTauPVX_3D->clear(); xCTauPVX_3D_err->clear();
  /// dxy, dz, dxyE, dzE for kaons from PV, BS, Bs0LessPV
  kaon1_dxy_PV->clear(); kaon1_dz_PV->clear(); kaon2_dxy_PV->clear(); kaon2_dz_PV->clear(); 
  kaon1_dxy_BS->clear(); kaon1_dz_BS->clear(); kaon2_dxy_BS->clear(); kaon2_dz_BS->clear(); 
  kaon1_dxy_Bs0LessPV->clear(); kaon1_dz_Bs0LessPV->clear(); kaon2_dxy_Bs0LessPV->clear(); kaon2_dz_Bs0LessPV->clear(); 
  kaon1_dxyE->clear(); kaon1_dzE->clear(); kaon2_dxyE->clear(); kaon2_dzE->clear();

  KKMass_err->clear(); Kaon1FromPV->clear(); Kaon2FromPV->clear();

  if (Debug_) cout <<"before muon stuff clear" <<endl ;
  /// muons
  muPx->clear(); muPy->clear(); muPz->clear(); muCharge->clear();
  muD0->clear(); muDz->clear(); muChi2->clear(); muGlChi2->clear();
  mufHits->clear(); muFirstBarrel->clear(); muFirstEndCap->clear(); muD0E->clear() ;  muDzVtxErr->clear() ; muKey->clear() ;
  muIsGlobal->clear(); muIsPF->clear();
  muDzVtx->clear(); muDxyVtx->clear(); muGlMatchedStation->clear(); muGlDzVtx->clear(); muGlDxyVtx->clear(); 
  nMatchedStations->clear();
  muNDF->clear(); muGlNDF->clear(); muPhits->clear(); muShits->clear(); muGlMuHits->clear(); muType->clear(); 
  muQual->clear(); muTrack->clear(); muNOverlap->clear(); muNSharingSegWith->clear(); 
  
  if (Debug_) cout <<"after muon stuff clear" <<endl ;
  /// tracks
  trNotRef->clear(); trRef->clear();
  trPx->clear(); trPy->clear(); trPz->clear(); trE->clear();	
  trNDF->clear(); trPhits->clear(); trShits->clear(); trChi2->clear();
  trD0->clear(); trD0E->clear(); trCharge->clear();
  trQualityHighPurity->clear(); trQualityTight->clear();
  trfHits->clear(); trFirstBarrel->clear(); trFirstEndCap->clear();
  trDzVtx->clear(); trDxyVtx->clear();
  tr_nsigdedx->clear(); tr_dedx->clear(); tr_dedxMass->clear(); tr_theo->clear(); tr_sigma->clear();
  tr_dedx_byHits->clear(); tr_dedxErr_byHits->clear(); tr_saturMeas_byHits->clear(); tr_Meas_byHits->clear();
  
  if (Debug_) cout <<"end of branches clear" <<endl ;
 }
//}/// analyze
/// ------------ method called once each job just before starting event loop  ------------
void MuMuKKPAT::beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup)
{
}
void MuMuKKPAT::beginJob()
{
  edm::Service<TFileService> fs;
	
  X_One_Tree_ = fs->make<TTree>("X_data", "X(4140) Data");
  
  X_One_Tree_->Branch("TrigRes", &trigRes);
  X_One_Tree_->Branch("TrigNames", &trigNames);
  X_One_Tree_->Branch("MatchTriggerNames", &MatchTriggerNames);
  X_One_Tree_->Branch("L1TrigRes", &L1TT);  
  X_One_Tree_->Branch("evtNum", &evtNum,"evtNum/i");
  X_One_Tree_->Branch("runNum", &runNum,"runNum/i");
  X_One_Tree_->Branch("lumiNum", &lumiNum, "lumiNum/i");
  X_One_Tree_->Branch("priVtx_n", &priVtx_n, "priVtx_n/i");
  X_One_Tree_->Branch("priVtx_X", &priVtx_X, "priVtx_X/f");
  X_One_Tree_->Branch("priVtx_Y", &priVtx_Y, "priVtx_Y/f");
  X_One_Tree_->Branch("priVtx_Z", &priVtx_Z, "priVtx_Z/f");
  X_One_Tree_->Branch("priVtx_XE", &priVtx_XE, "priVtx_XE/f");
  X_One_Tree_->Branch("priVtx_YE", &priVtx_YE, "priVtx_YE/f");
  X_One_Tree_->Branch("priVtx_ZE", &priVtx_ZE, "priVtx_ZE/f");
  X_One_Tree_->Branch("priVtx_NormChi2",&priVtx_NormChi2, "priVtx_NormChi2/f");
  X_One_Tree_->Branch("priVtx_Chi2",&priVtx_Chi2, "priVtx_Chi2/f");
  X_One_Tree_->Branch("priVtx_CL",&priVtx_CL, "priVtx_CL/f");
  X_One_Tree_->Branch("priVtx_tracks", &priVtx_tracks, "priVtx_tracks/i");
  X_One_Tree_->Branch("priVtx_tracksPtSq", &priVtx_tracksPtSq, "priVtx_tracksPtSq/f");
  /// MC Analysis
  if (doMC) {
    // Gen Primary Vertex
    X_One_Tree_->Branch("genEvtVtx_X", &genEvtVtx_X); 
    X_One_Tree_->Branch("genEvtVtx_Y", &genEvtVtx_Y);
    X_One_Tree_->Branch("genEvtVtx_Z", &genEvtVtx_Z);
    X_One_Tree_->Branch("genEvtVtx_particles", &genEvtVtx_particles);
    X_One_Tree_->Branch("n_Bs0Ancestors", &n_Bs0Ancestors);  
    X_One_Tree_->Branch("nMCAll", &nMCAll, "nMCAll/i");
    X_One_Tree_->Branch("MCPdgIdAll", &MCPdgIdAll);
    X_One_Tree_->Branch("MCDanNumAll", &MCDanNumAll);
    X_One_Tree_->Branch("nMCBs0",&nMCBs0,"nMCBs0/i"); 
    // Gen Primary Vertex
    X_One_Tree_->Branch("PriVtxGen_X",&PriVtxGen_X);
    X_One_Tree_->Branch("PriVtxGen_Y",&PriVtxGen_Y);
    X_One_Tree_->Branch("PriVtxGen_Z",&PriVtxGen_Z);
    X_One_Tree_->Branch("PriVtxGen_EX",&PriVtxGen_EX);
    X_One_Tree_->Branch("PriVtxGen_EY",&PriVtxGen_EY);
    X_One_Tree_->Branch("PriVtxGen_EZ",&PriVtxGen_EZ);
    X_One_Tree_->Branch("PriVtxGen_Chi2",&PriVtxGen_Chi2);
    X_One_Tree_->Branch("PriVtxGen_CL",&PriVtxGen_CL);
    X_One_Tree_->Branch("PriVtxGen_Ndof",&PriVtxGen_Ndof);
    X_One_Tree_->Branch("PriVtxGen_tracks",&PriVtxGen_tracks);
    X_One_Tree_->Branch("MCJPsiPx",&MCJPsiPx);
    X_One_Tree_->Branch("MCJPsiPy",&MCJPsiPy);
    X_One_Tree_->Branch("MCJPsiPz",&MCJPsiPz);
    X_One_Tree_->Branch("MCmupPx",&MCmupPx);
    X_One_Tree_->Branch("MCmupPy",&MCmupPy);
    X_One_Tree_->Branch("MCmupPz",&MCmupPz);
    X_One_Tree_->Branch("MCmumPx",&MCmumPx);
    X_One_Tree_->Branch("MCmumPy",&MCmumPy);
    X_One_Tree_->Branch("MCmumPz",&MCmumPz);
    X_One_Tree_->Branch("MCPhiPx",&MCPhiPx);
    X_One_Tree_->Branch("MCPhiPy",&MCPhiPy);
    X_One_Tree_->Branch("MCPhiPz",&MCPhiPz);
    X_One_Tree_->Branch("MCkpPx",&MCkpPx);
    X_One_Tree_->Branch("MCkpPy",&MCkpPy);
    X_One_Tree_->Branch("MCkpPz",&MCkpPz);
    X_One_Tree_->Branch("MCkmPx",&MCkmPx);
    X_One_Tree_->Branch("MCkmPy",&MCkmPy);
    X_One_Tree_->Branch("MCkmPz",&MCkmPz);
    //X_One_Tree_->Branch("MCpionPx",&MCpionPx);
    //X_One_Tree_->Branch("MCpionPy",&MCpionPy);
    //X_One_Tree_->Branch("MCpionPz",&MCpionPz);
    //X_One_Tree_->Branch("MCpionCh",&MCpionCh);
    //X_One_Tree_->Branch("MCkaonPx",&MCkaonPx);
    //X_One_Tree_->Branch("MCkaonPy",&MCkaonPy);
    //X_One_Tree_->Branch("MCkaonPz",&MCkaonPz);
    //X_One_Tree_->Branch("MCkaonCh",&MCkaonCh);
    X_One_Tree_->Branch("MCPx", &MCPx);
    X_One_Tree_->Branch("MCPy", &MCPy);
    X_One_Tree_->Branch("MCPz", &MCPz);
  }
  /// generic tracks
  X_One_Tree_->Branch("trNotRef", &trNotRef);
  X_One_Tree_->Branch("trRef", &trRef);
  X_One_Tree_->Branch("trackPx", &trPx);
  X_One_Tree_->Branch("trackPy", &trPy);
  X_One_Tree_->Branch("trackPz", &trPz);
  X_One_Tree_->Branch("trackEnergy", &trE);
  X_One_Tree_->Branch("trackNDF", &trNDF);
  X_One_Tree_->Branch("trackPhits", &trPhits);
  X_One_Tree_->Branch("trackShits", &trShits);
  X_One_Tree_->Branch("trackChi2", &trChi2);
  X_One_Tree_->Branch("trackD0", &trD0);
  X_One_Tree_->Branch("trackD0Err", &trD0E);
  X_One_Tree_->Branch("trackCharge", &trCharge);
  X_One_Tree_->Branch("TrackHighPurity", &trQualityHighPurity);
  X_One_Tree_->Branch("TrackTight", &trQualityTight);
  X_One_Tree_->Branch("trackfHits", &trfHits);
  X_One_Tree_->Branch("trackFirstBarrel", &trFirstBarrel);
  X_One_Tree_->Branch("trackFirstEndCap", &trFirstEndCap);
  X_One_Tree_->Branch("trackDzVtx", &trDzVtx);
  X_One_Tree_->Branch("trackDxyVtx", &trDxyVtx);
  X_One_Tree_->Branch("tr_nsigdedx", &tr_nsigdedx);
  X_One_Tree_->Branch("tr_dedx", &tr_dedx);
  X_One_Tree_->Branch("tr_dedxMass", &tr_dedxMass);
  X_One_Tree_->Branch("tr_theo", &tr_theo);
  X_One_Tree_->Branch("tr_sigma", &tr_sigma);
  X_One_Tree_->Branch("tr_dedx_byHits", &tr_dedx_byHits );
  X_One_Tree_->Branch("tr_dedxErr_byHits", &tr_dedxErr_byHits );
  X_One_Tree_->Branch("tr_saturMeas_byHits", &tr_saturMeas_byHits );
  X_One_Tree_->Branch("tr_Meas_byHits", &tr_Meas_byHits );
  /// Generic muons
  X_One_Tree_->Branch("nMu", &nMu, "nMu/i");
  X_One_Tree_->Branch("muPx",&muPx);
  X_One_Tree_->Branch("muPy",&muPy);
  X_One_Tree_->Branch("muPz",&muPz);
  X_One_Tree_->Branch("muCharge", &muCharge);
  X_One_Tree_->Branch("muD0",&muD0);
  X_One_Tree_->Branch("muDz",&muDz);
  X_One_Tree_->Branch("muChi2",&muChi2);
  X_One_Tree_->Branch("muNDF",&muNDF);
  X_One_Tree_->Branch("muPhits",&muPhits);
  X_One_Tree_->Branch("muShits",&muShits);
  X_One_Tree_->Branch("muLayersTr",&muLayersTr);
  X_One_Tree_->Branch("muLayersPix",&muLayersPix);
  X_One_Tree_->Branch("muD0E",&muD0E);
  X_One_Tree_->Branch("muDzVtxErr",&muDzVtxErr);
  X_One_Tree_->Branch("muKey",&muKey);
  X_One_Tree_->Branch("muIsGlobal",&muIsGlobal);
  X_One_Tree_->Branch("muIsPF",&muIsPF);
  X_One_Tree_->Branch("muGlMuHits",&muGlMuHits);
  X_One_Tree_->Branch("muGlChi2",&muGlChi2);
  X_One_Tree_->Branch("muGlNDF",&muGlNDF);
  X_One_Tree_->Branch("muGlMatchedStation",&muGlMatchedStation);
  X_One_Tree_->Branch("muGlDzVtx", &muGlDzVtx);
  X_One_Tree_->Branch("muGlDxyVtx", &muGlDxyVtx);
  X_One_Tree_->Branch("nMatchedStations", &nMatchedStations);
  X_One_Tree_->Branch("muType",&muType);
  X_One_Tree_->Branch("muQual",&muQual);
  X_One_Tree_->Branch("muTrack",&muTrack);
  X_One_Tree_->Branch("muNOverlap",&muNOverlap);
  X_One_Tree_->Branch("muNSharingSegWith",&muNSharingSegWith);
  X_One_Tree_->Branch("mufHits", &mufHits);
  X_One_Tree_->Branch("muFirstBarrel", &muFirstBarrel);
  X_One_Tree_->Branch("muFirstEndCap", &muFirstEndCap);
  X_One_Tree_->Branch("muDzVtx", &muDzVtx);
  X_One_Tree_->Branch("muDxyVtx", &muDxyVtx); 
  /// MuMu cand  
  X_One_Tree_->Branch("nMuMu",&nMuMu,"nMuMu/i");
  X_One_Tree_->Branch("MuMuMass",&MuMuMass);
  X_One_Tree_->Branch("MuMuPx",&MuMuPx);
  X_One_Tree_->Branch("MuMuPy",&MuMuPy);
  X_One_Tree_->Branch("MuMuPz",&MuMuPz);
  X_One_Tree_->Branch("MuMuVtx_CL",&MuMuVtx_CL);
  X_One_Tree_->Branch("MuMuVtx_Chi2",&MuMuVtx_Chi2);
  X_One_Tree_->Branch("MuMuDecayVtx_X",&MuMuDecayVtx_X);
  X_One_Tree_->Branch("MuMuDecayVtx_Y",&MuMuDecayVtx_Y);
  X_One_Tree_->Branch("MuMuDecayVtx_Z",&MuMuDecayVtx_Z);
  X_One_Tree_->Branch("MuMuDecayVtx_XE",&MuMuDecayVtx_XE);
  X_One_Tree_->Branch("MuMuDecayVtx_YE",&MuMuDecayVtx_YE);
  X_One_Tree_->Branch("MuMuDecayVtx_ZE",&MuMuDecayVtx_ZE);
  /// muons from JPsi (MuMu) fit  
  X_One_Tree_->Branch("mu1Idx",&mu1Idx);
  X_One_Tree_->Branch("mu2Idx",&mu2Idx);
  X_One_Tree_->Branch("mu1Px_MuMu",&mu1_MuMu_Px);
  X_One_Tree_->Branch("mu1Py_MuMu",&mu1_MuMu_Py);
  X_One_Tree_->Branch("mu1Pz_MuMu",&mu1_MuMu_Pz);
  X_One_Tree_->Branch("mu1Chi2_MuMu",&mu1_MuMu_Chi2);
  X_One_Tree_->Branch("mu1NDF_MuMu",&mu1_MuMu_NDF);
  X_One_Tree_->Branch("mu2Px_MuMu",&mu2_MuMu_Px);
  X_One_Tree_->Branch("mu2Py_MuMu",&mu2_MuMu_Py);
  X_One_Tree_->Branch("mu2Pz_MuMu",&mu2_MuMu_Pz);
  X_One_Tree_->Branch("mu2Chi2_MuMu",&mu2_MuMu_Chi2);
  X_One_Tree_->Branch("mu2NDF_MuMu",&mu2_MuMu_NDF);
  X_One_Tree_->Branch("MuMuType",&MuMuType);
  X_One_Tree_->Branch("MuMuMuonTrigMatch",&MuMuMuonTrigMatch);
  /// Primary Vertex with "MuMu correction"
  X_One_Tree_->Branch("PriVtxMuMuCorr_n", &PriVtxMuMuCorr_n);
  X_One_Tree_->Branch("PriVtxMuMuCorr_X", &PriVtxMuMuCorr_X);
  X_One_Tree_->Branch("PriVtxMuMuCorr_Y", &PriVtxMuMuCorr_Y);
  X_One_Tree_->Branch("PriVtxMuMuCorr_Z", &PriVtxMuMuCorr_Z);
  X_One_Tree_->Branch("PriVtxMuMuCorr_EX", &PriVtxMuMuCorr_EX);
  X_One_Tree_->Branch("PriVtxMuMuCorr_EY", &PriVtxMuMuCorr_EY);
  X_One_Tree_->Branch("PriVtxMuMuCorr_EZ", &PriVtxMuMuCorr_EZ);
  X_One_Tree_->Branch("PriVtxMuMuCorr_Chi2", &PriVtxMuMuCorr_Chi2);
  X_One_Tree_->Branch("PriVtxMuMuCorr_CL", &PriVtxMuMuCorr_CL);
  X_One_Tree_->Branch("PriVtxMuMuCorr_tracks", &PriVtxMuMuCorr_tracks);
  X_One_Tree_->Branch("nTrk_afterMuMu", &nTrk);
  /// KK cand
  X_One_Tree_->Branch("nKK",&nKK,"nKK/i");
  X_One_Tree_->Branch("KKMass",&KKMass);
  X_One_Tree_->Branch("KKPx",&KKPx);
  X_One_Tree_->Branch("KKPy",&KKPy);
  X_One_Tree_->Branch("KKPz",&KKPz);
  X_One_Tree_->Branch("KKVtx_CL",&KKVtx_CL);
  X_One_Tree_->Branch("KKVtx_Chi2",&KKVtx_Chi2);
  X_One_Tree_->Branch("KKDecayVtx_X",&KKDecayVtx_X);
  X_One_Tree_->Branch("KKDecayVtx_Y",&KKDecayVtx_Y);
  X_One_Tree_->Branch("KKDecayVtx_Z",&KKDecayVtx_Z);
  X_One_Tree_->Branch("KKDecayVtx_XE",&KKDecayVtx_XE);
  X_One_Tree_->Branch("KKDecayVtx_YE",&KKDecayVtx_YE);
  X_One_Tree_->Branch("KKDecayVtx_ZE",&KKDecayVtx_ZE);
  /// kaons from Phi (KK) fit
  X_One_Tree_->Branch("ka1Idx",&ka1Idx);
  X_One_Tree_->Branch("ka2Idx",&ka2Idx);
  X_One_Tree_->Branch("ka1Px_KK",&ka1_KK_Px);
  X_One_Tree_->Branch("ka1Py_KK",&ka1_KK_Py);
  X_One_Tree_->Branch("ka1Pz_KK",&ka1_KK_Pz);
  X_One_Tree_->Branch("ka1Chi2_KK",&ka1_KK_Chi2);
  X_One_Tree_->Branch("ka1NDF_KK",&ka1_KK_NDF);
  X_One_Tree_->Branch("ka2Px_KK",&ka2_KK_Px);
  X_One_Tree_->Branch("ka2Py_KK",&ka2_KK_Py);
  X_One_Tree_->Branch("ka2Pz_KK",&ka2_KK_Pz);
  X_One_Tree_->Branch("ka2Chi2_KK",&ka2_KK_Chi2);
  X_One_Tree_->Branch("ka2NDF_KK",&ka2_KK_NDF);
  X_One_Tree_->Branch("DRMuMuK1",&DRMuMuK1);
  X_One_Tree_->Branch("DRMuMuK2",&DRMuMuK2);
  X_One_Tree_->Branch("DRbs0K1",&DRbs0K1);
  X_One_Tree_->Branch("DRbs0K2",&DRbs0K2);
  X_One_Tree_->Branch("DRxK1",&DRxK1);
  X_One_Tree_->Branch("DRxK2",&DRxK2);  
  /// counters for Bs0 & X(4140)
  X_One_Tree_->Branch("nBs0",&nBs0,"nBs0/i");
  X_One_Tree_->Branch("nBs0_pre0",&nBs0_pre0,"nBs0_pre0/i");
  X_One_Tree_->Branch("nBs0_pre1",&nBs0_pre1,"nBs0_pre1/i");
  X_One_Tree_->Branch("nBs0_pre2",&nBs0_pre2,"nBs0_pre2/i");
  X_One_Tree_->Branch("nBs0_pre3",&nBs0_pre3,"nBs0_pre3/i");
  X_One_Tree_->Branch("nBs0_pre4",&nBs0_pre4,"nBs0_pre4/i");
  X_One_Tree_->Branch("nBs0_pre5",&nBs0_pre5,"nBs0_pre5/i");
  X_One_Tree_->Branch("nBs0_pre6",&nBs0_pre6,"nBs0_pre6/i");
  X_One_Tree_->Branch("nBs0_pre7",&nBs0_pre7,"nBs0_pre7/i");
  X_One_Tree_->Branch("nBs0_pre8",&nBs0_pre8,"nBs0_pre8/i");
  X_One_Tree_->Branch("nBs0_pre9",&nBs0_pre9,"nBs0_pre9/i");
  X_One_Tree_->Branch("nBs0_pre10",&nBs0_pre10,"nBs0_pre10/i");
  X_One_Tree_->Branch("nBs0_pre11",&nBs0_pre11,"nBs0_pre11/i");
  X_One_Tree_->Branch("nBs0_pre12",&nBs0_pre12,"nBs0_pre12/i");
  X_One_Tree_->Branch("nBs0_pre13",&nBs0_pre13,"nBs0_pre13/i");
  X_One_Tree_->Branch("nBs0_pre14",&nBs0_pre14,"nBs0_pre14/i");
  X_One_Tree_->Branch("nBs0_pre15",&nBs0_pre15,"nBs0_pre15/i");
  X_One_Tree_->Branch("nBs0_pre16",&nBs0_pre16,"nBs0_pre16/i");
  X_One_Tree_->Branch("nX",&nX,"nX/i");
  X_One_Tree_->Branch("nX_pre0",&nX_pre0,"nX_pre0/i");
  X_One_Tree_->Branch("nX_pre1",&nX_pre1,"nX_pre1/i");
  X_One_Tree_->Branch("nX_pre2",&nX_pre2,"nX_pre2/i");
  X_One_Tree_->Branch("nX_pre3",&nX_pre3,"nX_pre3/i");
  X_One_Tree_->Branch("nX_pre4",&nX_pre4,"nX_pre4/i");
  X_One_Tree_->Branch("nX_pre5",&nX_pre5,"nX_pre5/i");
  X_One_Tree_->Branch("nX_pre6",&nX_pre6,"nX_pre6/i");
  X_One_Tree_->Branch("nX_pre7",&nX_pre7,"nX_pre7/i");
  X_One_Tree_->Branch("nX_pre8",&nX_pre8,"nX_pre8/i");
  X_One_Tree_->Branch("nX_pre9",&nX_pre9,"nX_pre9/i");
  X_One_Tree_->Branch("nX_pre10",&nX_pre10,"nX_pre10/i");
  X_One_Tree_->Branch("nX_pre11",&nX_pre11,"nX_pre11/i");
  X_One_Tree_->Branch("nX_pre12",&nX_pre12,"nX_pre12/i");
  X_One_Tree_->Branch("nX_pre13",&nX_pre13,"nX_pre13/i");
  X_One_Tree_->Branch("nX_pre14",&nX_pre14,"nX_pre14/i");
  X_One_Tree_->Branch("nX_pre15",&nX_pre15,"nX_pre15/i");
  X_One_Tree_->Branch("nX_pre16",&nX_pre16,"nX_pre16/i");
  /// Bs0 cand & X(4140) cand 
  X_One_Tree_->Branch("Bs0Mass",&bs0Mass);
  X_One_Tree_->Branch("Bs0Px",&bs0Px);
  X_One_Tree_->Branch("Bs0Py",&bs0Py);
  X_One_Tree_->Branch("Bs0Pz",&bs0Pz);
  X_One_Tree_->Branch("Bs0PxE",&bs0PxE);
  X_One_Tree_->Branch("Bs0PyE",&bs0PyE);
  X_One_Tree_->Branch("Bs0PzE",&bs0PzE);
  X_One_Tree_->Branch("Bs0Vtx_CL",&bs0Vtx_CL);
  X_One_Tree_->Branch("Bs0Vtx_Chi2",&bs0Vtx_Chi2);
  X_One_Tree_->Branch("Bs0DecayVtx_X",&bs0DecayVtx_X);
  X_One_Tree_->Branch("Bs0DecayVtx_Y",&bs0DecayVtx_Y);
  X_One_Tree_->Branch("Bs0DecayVtx_Z",&bs0DecayVtx_Z);
  X_One_Tree_->Branch("Bs0DecayVtx_XE",&bs0DecayVtx_XE);
  X_One_Tree_->Branch("Bs0DecayVtx_YE",&bs0DecayVtx_YE);
  X_One_Tree_->Branch("Bs0DecayVtx_ZE",&bs0DecayVtx_ZE);
  X_One_Tree_->Branch("Bs0CosAlphaBS", &bs0CosAlphaBS);
  X_One_Tree_->Branch("Bs0CosAlpha3DBS", &bs0CosAlpha3DBS);
  X_One_Tree_->Branch("Bs0CTauBS", &bs0CTauBS);
  X_One_Tree_->Branch("Bs0CTauBSE", &bs0CTauBSE);
  X_One_Tree_->Branch("Bs0LxyBS", &bs0LxyBS);
  X_One_Tree_->Branch("Bs0LxyBSE", &bs0LxyBSE);
  X_One_Tree_->Branch("Bs0LxyzBS", &bs0LxyzBS);
  X_One_Tree_->Branch("Bs0LxyzBSE", &bs0LxyzBSE);
  X_One_Tree_->Branch("Bs0CosAlphaPV", &bs0CosAlphaPV);
  X_One_Tree_->Branch("Bs0CosAlpha3DPV", &bs0CosAlpha3DPV);
  X_One_Tree_->Branch("Bs0CTauPV", &bs0CTauPV);
  X_One_Tree_->Branch("Bs0CTauPVE", &bs0CTauPVE);
  X_One_Tree_->Branch("Bs0LxyPV", &bs0LxyPV);
  X_One_Tree_->Branch("Bs0LxyPVE", &bs0LxyPVE);
  X_One_Tree_->Branch("Bs0LxyzPV", &bs0LxyzPV);
  X_One_Tree_->Branch("Bs0LxyzPVE", &bs0LxyzPVE);
  X_One_Tree_->Branch("XMass",&xMass);
  X_One_Tree_->Branch("XPx",&xPx);
  X_One_Tree_->Branch("XPy",&xPy);
  X_One_Tree_->Branch("XPz",&xPz);
  X_One_Tree_->Branch("XPxE",&xPxE);
  X_One_Tree_->Branch("XPyE",&xPyE);
  X_One_Tree_->Branch("XPzE",&xPzE);
  X_One_Tree_->Branch("XVtx_CL",&xVtx_CL);
  X_One_Tree_->Branch("XVtx_Chi2",&xVtx_Chi2);
  X_One_Tree_->Branch("XDecayVtx_X",&xDecayVtx_X);
  X_One_Tree_->Branch("XDecayVtx_Y",&xDecayVtx_Y);
  X_One_Tree_->Branch("XDecayVtx_Z",&xDecayVtx_Z);
  X_One_Tree_->Branch("XDecayVtx_XE",&xDecayVtx_XE);
  X_One_Tree_->Branch("XDecayVtx_YE",&xDecayVtx_YE);
  X_One_Tree_->Branch("XDecayVtx_ZE",&xDecayVtx_ZE);
  X_One_Tree_->Branch("XCosAlphaBS", &xCosAlphaBS);
  X_One_Tree_->Branch("XCosAlpha3DBS", &xCosAlpha3DBS);
  X_One_Tree_->Branch("XCTauBS", &xCTauBS);
  X_One_Tree_->Branch("XCTauBSE", &xCTauBSE);
  X_One_Tree_->Branch("XLxyBS", &xLxyBS);
  X_One_Tree_->Branch("XLxyBSE", &xLxyBSE);
  X_One_Tree_->Branch("XLxyzBS", &xLxyzBS);
  X_One_Tree_->Branch("XLxyzBSE", &xLxyzBSE);
  X_One_Tree_->Branch("XCosAlphaPV", &xCosAlphaPV);
  X_One_Tree_->Branch("XCosAlpha3DPV", &xCosAlpha3DPV);
  X_One_Tree_->Branch("XCTauPV", &xCTauPV);
  X_One_Tree_->Branch("XCTauPVE", &xCTauPVE);
  X_One_Tree_->Branch("XLxyPV", &xLxyPV);
  X_One_Tree_->Branch("XLxyPVE", &xLxyPVE);
  X_One_Tree_->Branch("XLxyzPV", &xLxyzPV);
  X_One_Tree_->Branch("XLxyzPVE", &xLxyzPVE); 
  /// Primary Vertex with largest Bs0_cos(alpha) & largest X(4140)_cos(alpha)
  X_One_Tree_->Branch("PriVtx_Bs0CosAlpha_n",&PriVtx_Bs0CosAlpha_n);
  X_One_Tree_->Branch("PriVtx_Bs0CosAlpha_X",&PriVtx_Bs0CosAlpha_X);
  X_One_Tree_->Branch("PriVtx_Bs0CosAlpha_Y",&PriVtx_Bs0CosAlpha_Y);
  X_One_Tree_->Branch("PriVtx_Bs0CosAlpha_Z",&PriVtx_Bs0CosAlpha_Z);
  X_One_Tree_->Branch("PriVtx_Bs0CosAlpha_EX",&PriVtx_Bs0CosAlpha_EX);
  X_One_Tree_->Branch("PriVtx_Bs0CosAlpha_EY",&PriVtx_Bs0CosAlpha_EY);
  X_One_Tree_->Branch("PriVtx_Bs0CosAlpha_EZ",&PriVtx_Bs0CosAlpha_EZ);
  X_One_Tree_->Branch("PriVtx_Bs0CosAlpha_Chi2",&PriVtx_Bs0CosAlpha_Chi2);
  X_One_Tree_->Branch("PriVtx_Bs0CosAlpha_CL",&PriVtx_Bs0CosAlpha_CL);
  X_One_Tree_->Branch("PriVtx_Bs0CosAlpha_tracks",&PriVtx_Bs0CosAlpha_tracks);
  X_One_Tree_->Branch("Bs0CosAlphaPVCosAlpha", &bs0CosAlphaPVCosAlpha);
  X_One_Tree_->Branch("Bs0CosAlpha3DPVCosAlpha", &bs0CosAlpha3DPVCosAlpha);
  X_One_Tree_->Branch("Bs0CTauPVCosAlpha", &bs0CTauPVCosAlpha);
  X_One_Tree_->Branch("Bs0CTauPVCosAlphaE", &bs0CTauPVCosAlphaE);
  X_One_Tree_->Branch("Bs0LxyPVCosAlpha", &bs0LxyPVCosAlpha);
  X_One_Tree_->Branch("Bs0LxyPVCosAlphaE", &bs0LxyPVCosAlphaE);
  X_One_Tree_->Branch("Bs0LxyzPVCosAlpha", &bs0LxyzPVCosAlpha);
  X_One_Tree_->Branch("Bs0LxyzPVCosAlphaE", &bs0LxyzPVCosAlphaE);
  X_One_Tree_->Branch("PriVtx_XCosAlpha_n",&PriVtx_XCosAlpha_n);
  X_One_Tree_->Branch("PriVtx_XCosAlpha_X",&PriVtx_XCosAlpha_X);
  X_One_Tree_->Branch("PriVtx_XCosAlpha_Y",&PriVtx_XCosAlpha_Y);
  X_One_Tree_->Branch("PriVtx_XCosAlpha_Z",&PriVtx_XCosAlpha_Z);
  X_One_Tree_->Branch("PriVtx_XCosAlpha_EX",&PriVtx_XCosAlpha_EX);
  X_One_Tree_->Branch("PriVtx_XCosAlpha_EY",&PriVtx_XCosAlpha_EY);
  X_One_Tree_->Branch("PriVtx_XCosAlpha_EZ",&PriVtx_XCosAlpha_EZ);
  X_One_Tree_->Branch("PriVtx_XCosAlpha_Chi2",&PriVtx_XCosAlpha_Chi2);
  X_One_Tree_->Branch("PriVtx_XCosAlpha_CL",&PriVtx_XCosAlpha_CL);
  X_One_Tree_->Branch("PriVtx_XCosAlpha_tracks",&PriVtx_XCosAlpha_tracks);
  X_One_Tree_->Branch("XCosAlphaPVCosAlpha", &xCosAlphaPVCosAlpha);
  X_One_Tree_->Branch("XCosAlpha3DPVCosAlpha", &xCosAlpha3DPVCosAlpha);
  X_One_Tree_->Branch("XCTauPVCosAlpha", &xCTauPVCosAlpha);
  X_One_Tree_->Branch("XCTauPVCosAlphaE", &xCTauPVCosAlphaE);
  X_One_Tree_->Branch("XLxyPVCosAlpha", &xLxyPVCosAlpha);
  X_One_Tree_->Branch("XLxyPVCosAlphaE", &xLxyPVCosAlphaE);
  X_One_Tree_->Branch("XLxyzPVCosAlpha", &xLxyzPVCosAlpha);
  X_One_Tree_->Branch("XLxyzPVCosAlphaE", &xLxyzPVCosAlphaE);
  X_One_Tree_->Branch("PriVtx_Bs0CosAlpha3D_n",&PriVtx_Bs0CosAlpha3D_n);
  X_One_Tree_->Branch("PriVtx_Bs0CosAlpha3D_X",&PriVtx_Bs0CosAlpha3D_X);
  X_One_Tree_->Branch("PriVtx_Bs0CosAlpha3D_Y",&PriVtx_Bs0CosAlpha3D_Y);
  X_One_Tree_->Branch("PriVtx_Bs0CosAlpha3D_Z",&PriVtx_Bs0CosAlpha3D_Z);
  X_One_Tree_->Branch("PriVtx_Bs0CosAlpha3D_EX",&PriVtx_Bs0CosAlpha3D_EX);
  X_One_Tree_->Branch("PriVtx_Bs0CosAlpha3D_EY",&PriVtx_Bs0CosAlpha3D_EY);
  X_One_Tree_->Branch("PriVtx_Bs0CosAlpha3D_EZ",&PriVtx_Bs0CosAlpha3D_EZ);
  X_One_Tree_->Branch("PriVtx_Bs0CosAlpha3D_Chi2",&PriVtx_Bs0CosAlpha3D_Chi2);
  X_One_Tree_->Branch("PriVtx_Bs0CosAlpha3D_CL",&PriVtx_Bs0CosAlpha3D_CL);
  X_One_Tree_->Branch("PriVtx_Bs0CosAlpha3D_tracks",&PriVtx_Bs0CosAlpha3D_tracks);
  X_One_Tree_->Branch("Bs0CosAlphaPVCosAlpha3D", &bs0CosAlphaPVCosAlpha3D);
  X_One_Tree_->Branch("Bs0CosAlpha3DPVCosAlpha3D", &bs0CosAlpha3DPVCosAlpha3D);
  X_One_Tree_->Branch("Bs0CTauPVCosAlpha3D", &bs0CTauPVCosAlpha3D);
  X_One_Tree_->Branch("Bs0CTauPVCosAlpha3DE", &bs0CTauPVCosAlpha3DE);
  X_One_Tree_->Branch("Bs0LxyPVCosAlpha3D", &bs0LxyPVCosAlpha3D);
  X_One_Tree_->Branch("Bs0LxyPVCosAlpha3DE", &bs0LxyPVCosAlpha3DE);
  X_One_Tree_->Branch("Bs0LxyzPVCosAlpha3D", &bs0LxyzPVCosAlpha3D);
  X_One_Tree_->Branch("Bs0LxyzPVCosAlpha3DE", &bs0LxyzPVCosAlpha3DE);
  X_One_Tree_->Branch("PriVtx_XCosAlpha3D_n",&PriVtx_XCosAlpha3D_n);
  X_One_Tree_->Branch("PriVtx_XCosAlpha3D_X",&PriVtx_XCosAlpha3D_X);
  X_One_Tree_->Branch("PriVtx_XCosAlpha3D_Y",&PriVtx_XCosAlpha3D_Y);
  X_One_Tree_->Branch("PriVtx_XCosAlpha3D_Z",&PriVtx_XCosAlpha3D_Z);
  X_One_Tree_->Branch("PriVtx_XCosAlpha3D_EX",&PriVtx_XCosAlpha3D_EX);
  X_One_Tree_->Branch("PriVtx_XCosAlpha3D_EY",&PriVtx_XCosAlpha3D_EY);
  X_One_Tree_->Branch("PriVtx_XCosAlpha3D_EZ",&PriVtx_XCosAlpha3D_EZ);
  X_One_Tree_->Branch("PriVtx_XCosAlpha3D_Chi2",&PriVtx_XCosAlpha3D_Chi2);
  X_One_Tree_->Branch("PriVtx_XCosAlpha3D_CL",&PriVtx_XCosAlpha3D_CL);
  X_One_Tree_->Branch("PriVtx_XCosAlpha3D_tracks",&PriVtx_XCosAlpha3D_tracks);
  X_One_Tree_->Branch("XCosAlphaPVCosAlpha3D", &xCosAlphaPVCosAlpha3D);
  X_One_Tree_->Branch("XCosAlpha3DPVCosAlpha3D", &xCosAlpha3DPVCosAlpha3D);
  X_One_Tree_->Branch("XCTauPVCosAlpha3D", &xCTauPVCosAlpha3D);
  X_One_Tree_->Branch("XCTauPVCosAlpha3DE", &xCTauPVCosAlpha3DE);
  X_One_Tree_->Branch("XLxyPVCosAlpha3D", &xLxyPVCosAlpha3D);
  X_One_Tree_->Branch("XLxyPVCosAlpha3DE", &xLxyPVCosAlpha3DE);
  X_One_Tree_->Branch("XLxyzPVCosAlpha3D", &xLxyzPVCosAlpha3D);
  X_One_Tree_->Branch("XLxyzPVCosAlpha3DE", &xLxyzPVCosAlpha3DE);
 
  X_One_Tree_->Branch("Bs0LessPV_tracksPtSq",&Bs0LessPV_tracksPtSq);
  X_One_Tree_->Branch("Bs0LessPV_4tracksPtSq",&Bs0LessPV_4tracksPtSq);
  X_One_Tree_->Branch("PriVtxBs0Less_n",&PriVtxBs0Less_n);
  X_One_Tree_->Branch("PriVtxBs0Less_X",&PriVtxBs0Less_X);
  X_One_Tree_->Branch("PriVtxBs0Less_Y",&PriVtxBs0Less_Y);
  X_One_Tree_->Branch("PriVtxBs0Less_Z",&PriVtxBs0Less_Z);
  X_One_Tree_->Branch("PriVtxBs0Less_EX",&PriVtxBs0Less_EX);
  X_One_Tree_->Branch("PriVtxBs0Less_EY",&PriVtxBs0Less_EY);
  X_One_Tree_->Branch("PriVtxBs0Less_EZ",&PriVtxBs0Less_EZ);
  X_One_Tree_->Branch("PriVtxBs0Less_Chi2",&PriVtxBs0Less_Chi2);
  X_One_Tree_->Branch("PriVtxBs0Less_CL",&PriVtxBs0Less_CL);
  X_One_Tree_->Branch("PriVtxBs0Less_tracks",&PriVtxBs0Less_tracks);
  X_One_Tree_->Branch("Bs0CosAlphaBs0LessPV", &bs0CosAlphaBs0LessPV);
  X_One_Tree_->Branch("Bs0CosAlpha3DBs0LessPV", &bs0CosAlpha3DBs0LessPV);
  X_One_Tree_->Branch("Bs0CTauBs0LessPV", &bs0CTauBs0LessPV);
  X_One_Tree_->Branch("Bs0CTauBs0LessPVE", &bs0CTauBs0LessPVE);
  X_One_Tree_->Branch("Bs0LxyBs0LessPV", &bs0LxyBs0LessPV);
  X_One_Tree_->Branch("Bs0LxyBs0LessPVE", &bs0LxyBs0LessPVE);
  X_One_Tree_->Branch("Bs0LxyzBs0LessPV", &bs0LxyzBs0LessPV);
  X_One_Tree_->Branch("Bs0LxyzBs0LessPVE", &bs0LxyzBs0LessPVE);  	
  X_One_Tree_->Branch("PriVtxBs0Less_Bs0CosAlpha_n",&PriVtxBs0Less_Bs0CosAlpha_n);
  X_One_Tree_->Branch("PriVtxBs0Less_Bs0CosAlpha_X",&PriVtxBs0Less_Bs0CosAlpha_X);
  X_One_Tree_->Branch("PriVtxBs0Less_Bs0CosAlpha_Y",&PriVtxBs0Less_Bs0CosAlpha_Y);
  X_One_Tree_->Branch("PriVtxBs0Less_Bs0CosAlpha_Z",&PriVtxBs0Less_Bs0CosAlpha_Z);
  X_One_Tree_->Branch("PriVtxBs0Less_Bs0CosAlpha_EX",&PriVtxBs0Less_Bs0CosAlpha_EX);
  X_One_Tree_->Branch("PriVtxBs0Less_Bs0CosAlpha_EY",&PriVtxBs0Less_Bs0CosAlpha_EY);
  X_One_Tree_->Branch("PriVtxBs0Less_Bs0CosAlpha_EZ",&PriVtxBs0Less_Bs0CosAlpha_EZ);
  X_One_Tree_->Branch("PriVtxBs0Less_Bs0CosAlpha_Chi2",&PriVtxBs0Less_Bs0CosAlpha_Chi2);
  X_One_Tree_->Branch("PriVtxBs0Less_Bs0CosAlpha_CL",&PriVtxBs0Less_Bs0CosAlpha_CL);
  X_One_Tree_->Branch("PriVtxBs0Less_Bs0CosAlpha_tracks",&PriVtxBs0Less_Bs0CosAlpha_tracks);
  X_One_Tree_->Branch("Bs0CosAlphaBs0LessPVCosAlpha", &bs0CosAlphaBs0LessPVCosAlpha);
  X_One_Tree_->Branch("Bs0CosAlpha3DBs0LessPVCosAlpha", &bs0CosAlpha3DBs0LessPVCosAlpha);
  X_One_Tree_->Branch("Bs0CTauBs0LessPVCosAlpha", &bs0CTauBs0LessPVCosAlpha);
  X_One_Tree_->Branch("Bs0CTauBs0LessPVCosAlphaE", &bs0CTauBs0LessPVCosAlphaE);
  X_One_Tree_->Branch("Bs0LxyBs0LessPVCosAlpha", &bs0LxyBs0LessPVCosAlpha);
  X_One_Tree_->Branch("Bs0LxyBs0LessPVCosAlphaE", &bs0LxyBs0LessPVCosAlphaE);
  X_One_Tree_->Branch("Bs0LxyzBs0LessPVCosAlpha", &bs0LxyzBs0LessPVCosAlpha);
  X_One_Tree_->Branch("Bs0LxyzBs0LessPVCosAlphaE", &bs0LxyzBs0LessPVCosAlphaE);
  X_One_Tree_->Branch("PriVtxBs0Less_Bs0CosAlpha3D_n",&PriVtxBs0Less_Bs0CosAlpha3D_n);
  X_One_Tree_->Branch("PriVtxBs0Less_Bs0CosAlpha3D_X",&PriVtxBs0Less_Bs0CosAlpha3D_X);
  X_One_Tree_->Branch("PriVtxBs0Less_Bs0CosAlpha3D_Y",&PriVtxBs0Less_Bs0CosAlpha3D_Y);
  X_One_Tree_->Branch("PriVtxBs0Less_Bs0CosAlpha3D_Z",&PriVtxBs0Less_Bs0CosAlpha3D_Z);
  X_One_Tree_->Branch("PriVtxBs0Less_Bs0CosAlpha3D_EX",&PriVtxBs0Less_Bs0CosAlpha3D_EX);
  X_One_Tree_->Branch("PriVtxBs0Less_Bs0CosAlpha3D_EY",&PriVtxBs0Less_Bs0CosAlpha3D_EY);
  X_One_Tree_->Branch("PriVtxBs0Less_Bs0CosAlpha3D_EZ",&PriVtxBs0Less_Bs0CosAlpha3D_EZ);
  X_One_Tree_->Branch("PriVtxBs0Less_Bs0CosAlpha3D_Chi2",&PriVtxBs0Less_Bs0CosAlpha3D_Chi2);
  X_One_Tree_->Branch("PriVtxBs0Less_Bs0CosAlpha3D_CL",&PriVtxBs0Less_Bs0CosAlpha3D_CL);
  X_One_Tree_->Branch("PriVtxBs0Less_Bs0CosAlpha3D_tracks",&PriVtxBs0Less_Bs0CosAlpha3D_tracks);
  X_One_Tree_->Branch("Bs0CosAlphaBs0LessPVCosAlpha3D", &bs0CosAlphaBs0LessPVCosAlpha3D);
  X_One_Tree_->Branch("Bs0CosAlpha3DBs0LessPVCosAlpha3D", &bs0CosAlpha3DBs0LessPVCosAlpha3D);
  X_One_Tree_->Branch("Bs0CTauBs0LessPVCosAlpha3D", &bs0CTauBs0LessPVCosAlpha3D);
  X_One_Tree_->Branch("Bs0CTauBs0LessPVCosAlpha3DE", &bs0CTauBs0LessPVCosAlpha3DE);
  X_One_Tree_->Branch("Bs0LxyBs0LessPVCosAlpha3D", &bs0LxyBs0LessPVCosAlpha3D);
  X_One_Tree_->Branch("Bs0LxyBs0LessPVCosAlpha3DE", &bs0LxyBs0LessPVCosAlpha3DE);
  X_One_Tree_->Branch("Bs0LxyzBs0LessPVCosAlpha3D", &bs0LxyzBs0LessPVCosAlpha3D);
  X_One_Tree_->Branch("Bs0LxyzBs0LessPVCosAlpha3DE", &bs0LxyzBs0LessPVCosAlpha3DE);
  /// Primary Vertex with "Bs0 correction" & "X(4140) correction" 
  X_One_Tree_->Branch("PriVtxBs0Corr_n",&PriVtxBs0Corr_n);
  X_One_Tree_->Branch("PriVtxBs0Corr_X",&PriVtxBs0Corr_X);
  X_One_Tree_->Branch("PriVtxBs0Corr_Y",&PriVtxBs0Corr_Y);
  X_One_Tree_->Branch("PriVtxBs0Corr_Z",&PriVtxBs0Corr_Z);
  X_One_Tree_->Branch("PriVtxBs0Corr_EX",&PriVtxBs0Corr_EX);
  X_One_Tree_->Branch("PriVtxBs0Corr_EY",&PriVtxBs0Corr_EY);
  X_One_Tree_->Branch("PriVtxBs0Corr_EZ",&PriVtxBs0Corr_EZ);
  X_One_Tree_->Branch("PriVtxBs0Corr_Chi2",&PriVtxBs0Corr_Chi2);
  X_One_Tree_->Branch("PriVtxBs0Corr_CL",&PriVtxBs0Corr_CL);
  X_One_Tree_->Branch("PriVtxBs0Corr_tracks",&PriVtxBs0Corr_tracks);
  X_One_Tree_->Branch("PriVtxXCorr_n",&PriVtxXCorr_n);
  X_One_Tree_->Branch("PriVtxXCorr_X",&PriVtxXCorr_X);
  X_One_Tree_->Branch("PriVtxXCorr_Y",&PriVtxXCorr_Y);
  X_One_Tree_->Branch("PriVtxXCorr_Z",&PriVtxXCorr_Z);
  X_One_Tree_->Branch("PriVtxXCorr_EX",&PriVtxXCorr_EX);
  X_One_Tree_->Branch("PriVtxXCorr_EY",&PriVtxXCorr_EY);
  X_One_Tree_->Branch("PriVtxXCorr_EZ",&PriVtxXCorr_EZ);
  X_One_Tree_->Branch("PriVtxXCorr_Chi2",&PriVtxXCorr_Chi2);
  X_One_Tree_->Branch("PriVtxXCorr_CL",&PriVtxXCorr_CL);
  X_One_Tree_->Branch("PriVtxXCorr_tracks",&PriVtxXCorr_tracks);
  /// Lifetime variables for Bs0 & X(4140)
  X_One_Tree_->Branch("Bs0CosAlphaPVX", &bs0CosAlphaPVX);
  X_One_Tree_->Branch("Bs0CTauPVX", &bs0CTauPVX);
  X_One_Tree_->Branch("Bs0CTauPVXE", &bs0CTauPVXE);
  X_One_Tree_->Branch("Bs0LxyPVX", &bs0LxyPVX);
  X_One_Tree_->Branch("Bs0LxyzPVX", &bs0LxyzPVX);
  X_One_Tree_->Branch("Bs0CTauPVX_3D", &bs0CTauPVX_3D);
  X_One_Tree_->Branch("Bs0CTauPVX_3D_err", &bs0CTauPVX_3D_err);
  X_One_Tree_->Branch("XCosAlphaPVX", &xCosAlphaPVX);
  X_One_Tree_->Branch("XCTauPVX", &xCTauPVX);
  X_One_Tree_->Branch("XCTauPVXE", &xCTauPVXE);
  X_One_Tree_->Branch("XLxyPVX", &xLxyPVX);
  X_One_Tree_->Branch("XLxyzPVX", &xLxyzPVX);
  X_One_Tree_->Branch("XCTauPVX_3D", &xCTauPVX_3D);
  X_One_Tree_->Branch("XCTauPVX_3D_err", &xCTauPVX_3D_err);
  /// dxy, dz, dxyE, dzE for kaons from PV, BS, Bs0LessPV 
  X_One_Tree_->Branch("kaon1_dxy_PV", &kaon1_dxy_PV);
  X_One_Tree_->Branch("kaon1_dz_PV", &kaon1_dz_PV);
  X_One_Tree_->Branch("kaon2_dxy_PV", &kaon2_dxy_PV);
  X_One_Tree_->Branch("kaon2_dz_PV", &kaon2_dz_PV);
  X_One_Tree_->Branch("kaon1_dxy_BS", &kaon1_dxy_BS);
  X_One_Tree_->Branch("kaon1_dz_BS", &kaon1_dz_BS);
  X_One_Tree_->Branch("kaon2_dxy_BS", &kaon2_dxy_BS);
  X_One_Tree_->Branch("kaon2_dz_BS", &kaon2_dz_BS);
  X_One_Tree_->Branch("kaon1_dxy_Bs0LessPV", &kaon1_dxy_Bs0LessPV);
  X_One_Tree_->Branch("kaon1_dz_Bs0LessPV", &kaon1_dz_Bs0LessPV);
  X_One_Tree_->Branch("kaon2_dxy_Bs0LessPV", &kaon2_dxy_Bs0LessPV);
  X_One_Tree_->Branch("kaon2_dz_Bs0LessPV", &kaon2_dz_Bs0LessPV);
  X_One_Tree_->Branch("kaon1_dxyE", &kaon1_dxyE);
  X_One_Tree_->Branch("kaon1_dzE", &kaon1_dzE);
  X_One_Tree_->Branch("kaon2_dxyE", &kaon2_dxyE);
  X_One_Tree_->Branch("kaon2_dzE", &kaon2_dzE);
 
  X_One_Tree_->Branch("Bs0MuMuIdx", &Bs0_MuMuIdx);
  X_One_Tree_->Branch("Bs0Kaon1Idx", &Bs0_ka1Idx);
  X_One_Tree_->Branch("Bs0Kaon2Idx", &Bs0_ka2Idx); 
  X_One_Tree_->Branch("XMuMuIdx", &X_MuMuIdx);
  X_One_Tree_->Branch("XKaon1Idx", &X_ka1Idx); 
  X_One_Tree_->Branch("XKaon2Idx", &X_ka2Idx);

  X_One_Tree_->Branch("KKMass_err",&KKMass_err);
  X_One_Tree_->Branch("Kaon1FromPV",&Kaon1FromPV);
  X_One_Tree_->Branch("Kaon2FromPV",&Kaon2FromPV );

  /// Muons and tracks after Bs0 cand fit & X(4140) cand fit
  X_One_Tree_->Branch("Bs0_Muon1Px_MuMuKK", &Bs0_mu1Px_MuMuKK);
  X_One_Tree_->Branch("Bs0_Muon1Py_MuMuKK", &Bs0_mu1Py_MuMuKK);
  X_One_Tree_->Branch("Bs0_Muon1Pz_MuMuKK", &Bs0_mu1Pz_MuMuKK);
  X_One_Tree_->Branch("Bs0_Muon1E_MuMuKK", &Bs0_mu1E_MuMuKK);
  X_One_Tree_->Branch("X_Muon1Px_MuMuKK", &X_mu1Px_MuMuKK);
  X_One_Tree_->Branch("X_Muon1Py_MuMuKK", &X_mu1Py_MuMuKK);
  X_One_Tree_->Branch("X_Muon1Pz_MuMuKK", &X_mu1Pz_MuMuKK);
  X_One_Tree_->Branch("X_Muon1E_MuMuKK", &X_mu1E_MuMuKK); 
  X_One_Tree_->Branch("Bs0_Muon2Px_MuMuKK", &Bs0_mu2Px_MuMuKK);
  X_One_Tree_->Branch("Bs0_Muon2Py_MuMuKK", &Bs0_mu2Py_MuMuKK);
  X_One_Tree_->Branch("Bs0_Muon2Pz_MuMuKK", &Bs0_mu2Pz_MuMuKK);
  X_One_Tree_->Branch("Bs0_Muon2E_MuMuKK", &Bs0_mu2E_MuMuKK);
  X_One_Tree_->Branch("X_Muon2Px_MuMuKK", &X_mu2Px_MuMuKK);
  X_One_Tree_->Branch("X_Muon2Py_MuMuKK", &X_mu2Py_MuMuKK);
  X_One_Tree_->Branch("X_Muon2Pz_MuMuKK", &X_mu2Pz_MuMuKK);
  X_One_Tree_->Branch("X_Muon2E_MuMuKK", &X_mu2E_MuMuKK);
  X_One_Tree_->Branch("Bs0_Kaon1Px_MuMuKK", &Bs0_k1Px_MuMuKK); 
  X_One_Tree_->Branch("Bs0_Kaon1Py_MuMuKK", &Bs0_k1Py_MuMuKK); 
  X_One_Tree_->Branch("Bs0_Kaon1Pz_MuMuKK", &Bs0_k1Pz_MuMuKK); 
  X_One_Tree_->Branch("Bs0_Kion1E_MuMuKK", &Bs0_k1E_MuMuKK); 
  X_One_Tree_->Branch("Bs0_kaon1_nsigdedx", &Bs0_kaon1_nsigdedx); 
  X_One_Tree_->Branch("Bs0_kaon1_dedx", &Bs0_kaon1_dedx); 
  X_One_Tree_->Branch("Bs0_kaon1_dedxMass", &Bs0_kaon1_dedxMass); 
  X_One_Tree_->Branch("Bs0_kaon1_theo", &Bs0_kaon1_theo); 
  X_One_Tree_->Branch("Bs0_kaon1_sigma", &Bs0_kaon1_sigma); 
  X_One_Tree_->Branch("Bs0_kaon1_dedx_byHits", &Bs0_kaon1_dedx_byHits); 
  X_One_Tree_->Branch("Bs0_kaon1_dedxErr_byHits", &Bs0_kaon1_dedxErr_byHits); 
  X_One_Tree_->Branch("Bs0_kaon1_saturMeas_byHits", &Bs0_kaon1_saturMeas_byHits); 
  X_One_Tree_->Branch("Bs0_kaon1_Meas_byHits", &Bs0_kaon1_Meas_byHits); 
  X_One_Tree_->Branch("X_Kaon1Px_MuMuKK", &X_k1Px_MuMuKK); 
  X_One_Tree_->Branch("X_Kaon1Py_MuMuKK", &X_k1Py_MuMuKK); 
  X_One_Tree_->Branch("X_Kaon1Pz_MuMuKK", &X_k1Pz_MuMuKK); 
  X_One_Tree_->Branch("X_Kion1E_MuMuKK", &X_k1E_MuMuKK); 
  X_One_Tree_->Branch("X_kaon1_nsigdedx", &X_kaon1_nsigdedx); 
  X_One_Tree_->Branch("X_kaon1_dedx", &X_kaon1_dedx); 
  X_One_Tree_->Branch("X_kaon1_dedxMass", &X_kaon1_dedxMass); 
  X_One_Tree_->Branch("X_kaon1_theo", &X_kaon1_theo); 
  X_One_Tree_->Branch("X_kaon1_sigma", &X_kaon1_sigma); 
  X_One_Tree_->Branch("X_kaon1_dedx_byHits", &X_kaon1_dedx_byHits);
  X_One_Tree_->Branch("X_kaon1_dedxErr_byHits", &X_kaon1_dedxErr_byHits); 
  X_One_Tree_->Branch("X_kaon1_saturMeas_byHits", &X_kaon1_saturMeas_byHits);
  X_One_Tree_->Branch("X_kaon1_Meas_byHits", &X_kaon1_Meas_byHits); 
  X_One_Tree_->Branch("Bs0_Kaon2Px_MuMuKK", &Bs0_k2Px_MuMuKK); 
  X_One_Tree_->Branch("Bs0_Kaon2Py_MuMuKK", &Bs0_k2Py_MuMuKK); 
  X_One_Tree_->Branch("Bs0_Kaon2Pz_MuMuKK", &Bs0_k2Pz_MuMuKK); 
  X_One_Tree_->Branch("Bs0_Kaon2E_MuMuKK", &Bs0_k2E_MuMuKK); 
  X_One_Tree_->Branch("Bs0_kaon2_nsigdedx", &Bs0_kaon2_nsigdedx); 
  X_One_Tree_->Branch("Bs0_kaon2_dedx", &Bs0_kaon2_dedx); 
  X_One_Tree_->Branch("Bs0_kaon2_dedxMass", &Bs0_kaon2_dedxMass); 
  X_One_Tree_->Branch("Bs0_kaon2_theo", &Bs0_kaon2_theo); 
  X_One_Tree_->Branch("Bs0_kaon2_sigma", &Bs0_kaon2_sigma); 
  X_One_Tree_->Branch("Bs0_kaon2_dedx_byHits", &Bs0_kaon2_dedx_byHits); 
  X_One_Tree_->Branch("Bs0_kaon2_dedxErr_byHits", &Bs0_kaon2_dedxErr_byHits); 
  X_One_Tree_->Branch("Bs0_kaon2_saturMeas_byHits", &Bs0_kaon2_saturMeas_byHits); 
  X_One_Tree_->Branch("Bs0_kaon2_Meas_byHits", &Bs0_kaon2_Meas_byHits); 
  X_One_Tree_->Branch("X_Kaon2Px_MuMuKK", &X_k2Px_MuMuKK); 
  X_One_Tree_->Branch("X_Kaon2Py_MuMuKK", &X_k2Py_MuMuKK); 
  X_One_Tree_->Branch("X_Kaon2Pz_MuMuKK", &X_k2Pz_MuMuKK); 
  X_One_Tree_->Branch("X_Kaon2E_MuMuKK", &X_k2E_MuMuKK); 
  X_One_Tree_->Branch("X_kaon2_nsigdedx", &X_kaon2_nsigdedx); 
  X_One_Tree_->Branch("X_kaon2_dedx", &X_kaon2_dedx); 
  X_One_Tree_->Branch("X_kaon2_dedxMass", &X_kaon2_dedxMass); 
  X_One_Tree_->Branch("X_kaon2_theo", &X_kaon2_theo); 
  X_One_Tree_->Branch("X_kaon2_sigma", &X_kaon2_sigma);
  X_One_Tree_->Branch("X_kaon2_dedx_byHits", &X_kaon2_dedx_byHits); 
  X_One_Tree_->Branch("X_kaon2_dedxErr_byHits", &X_kaon2_dedxErr_byHits); 
  X_One_Tree_->Branch("X_kaon2_saturMeas_byHits", &X_kaon2_saturMeas_byHits); 
  X_One_Tree_->Branch("X_kaon2_Meas_byHits", &X_kaon2_Meas_byHits); 

}/// begin Job

/// ------------ method called once each job just after ending the event loop  ------------
void MuMuKKPAT::endJob() {
  X_One_Tree_->GetDirectory()->cd();
  X_One_Tree_->Write();
}/// endjob


bool MuMuKKPAT::isAbHadron(int pdgID) {
    
    if (abs(pdgID) == 511 || abs(pdgID) == 521 || abs(pdgID) == 531 || abs(pdgID) == 5122) return true;
    return false;
    
}

bool MuMuKKPAT::isAMixedbHadron(int pdgID, int momPdgID) {
    
    if ((abs(pdgID) == 511 && abs(momPdgID) == 511 && pdgID*momPdgID < 0) || 
        (abs(pdgID) == 531 && abs(momPdgID) == 531 && pdgID*momPdgID < 0)) 
        return true;
    return false;
    
}          

std::pair<int, float> MuMuKKPAT::findCandMCInfo(reco::GenParticleRef genCand) {
    
    int momJpsiID = 0;
    float trueLife = -99.;
    //cout <<"externalmodule"<<endl;
    
    if (genCand->numberOfMothers()>0) {
        
        TVector3 trueVtx(0.0,0.0,0.0);
        TVector3 trueP(0.0,0.0,0.0);
        TVector3 trueVtxMom(0.0,0.0,0.0);
        
        trueVtx.SetXYZ(genCand->vertex().x(),genCand->vertex().y(),genCand->vertex().z());
        trueP.SetXYZ(genCand->momentum().x(),genCand->momentum().y(),genCand->momentum().z());
        
        bool aBhadron = false;
        reco::GenParticleRef Candmom = genCand->motherRef();       // find mothers
        if (Candmom.isNull()) {
            std::pair<int, float> result = std::make_pair(momJpsiID, trueLife);
            return result;
        } else {
            reco::GenParticleRef CandGrandMom = Candmom->motherRef();
            if (isAbHadron(Candmom->pdgId())) {
                if (CandGrandMom.isNonnull() && isAMixedbHadron(Candmom->pdgId(),CandGrandMom->pdgId())) {
                    momJpsiID = CandGrandMom->pdgId();
                    trueVtxMom.SetXYZ(CandGrandMom->vertex().x(),CandGrandMom->vertex().y(),CandGrandMom->vertex().z());
                } else {
                    momJpsiID = Candmom->pdgId();
                    trueVtxMom.SetXYZ(Candmom->vertex().x(),Candmom->vertex().y(),Candmom->vertex().z());
                }
                aBhadron = true;
            } else {
                if (CandGrandMom.isNonnull() && isAbHadron(CandGrandMom->pdgId())) {
                    reco::GenParticleRef JpsiGrandgrandmom = CandGrandMom->motherRef();
                    if (JpsiGrandgrandmom.isNonnull() && isAMixedbHadron(CandGrandMom->pdgId(),JpsiGrandgrandmom->pdgId())) {
                        momJpsiID = JpsiGrandgrandmom->pdgId();
                        trueVtxMom.SetXYZ(JpsiGrandgrandmom->vertex().x(),JpsiGrandgrandmom->vertex().y(),JpsiGrandgrandmom->vertex().z());
                    } else {
                        momJpsiID = CandGrandMom->pdgId();
                        trueVtxMom.SetXYZ(CandGrandMom->vertex().x(),CandGrandMom->vertex().y(),CandGrandMom->vertex().z());
                    }
                    aBhadron = true;
                }
            }
            if (!aBhadron) {
                momJpsiID = Candmom->pdgId();
                trueVtxMom.SetXYZ(Candmom->vertex().x(),Candmom->vertex().y(),Candmom->vertex().z()); 
            }
        } 
        
        TVector3 vdiff = trueVtx - trueVtxMom;
        trueLife = vdiff.Perp()*genCand->mass()/trueP.Perp();
 }
    std::pair<int, float> result = std::make_pair(momJpsiID, trueLife);
    return result;
    
}

double MuMuKKPAT::getSigmaOfLogdEdx(double logde)
{
  return 0.3;
}

float MuMuKKPAT::getEnergyLoss(const reco::TrackRef & track)
{
  if (iexception_dedx==1) return 9999.;
  const reco::DeDxDataValueMap & eloss = *energyLoss;
  return eloss[track].dEdx();
}

double MuMuKKPAT::nsigmaofdedx(const reco::TrackRef & track, double & theo, double & sigma)
{  
    
  // no usable dE/dx if p > 2	
  double nsigma = 99 ;
  if (iexception_dedx==1) return nsigma ;
	
  double m  = 0.13957; 
  double bg = track->p() / m;
	
  theo = getLogdEdx(bg);
	
 
  int nhitr = track->numberOfValidHits();	
  double meas = log(getEnergyLoss(track));
  sigma = getSigmaOfLogdEdx(theo) * pow(nhitr,-0.65);
  if (sigma>0)
    nsigma = (meas-theo) / sigma ;
  return nsigma;
}


double MuMuKKPAT::getLogdEdx(double bg)
{
  const double a =  3.25 ;
  const double b =  0.288;
  const double c = -0.852;
	
  double beta = bg/sqrt(bg*bg + 1);
  double dedx = log( a/(beta*beta) + b * log(bg) + c );
	
  return dedx;
	
}


double MuMuKKPAT::GetMass(const reco::TrackRef & track){
  double P = track->p();
  double C = 2.625;
  double K = 2.495;
  double I = getEnergyLoss(track);
  return sqrt((I-C)/K)*P;
}


template<typename T>
bool MuMuKKPAT::isBetterMuon(const T &mu1, const T &mu2) const {
  if (mu2.track().isNull()) return true;
  if (mu1.track().isNull()) return false;
  if (mu1.isPFMuon() != mu2.isPFMuon()) return mu1.isPFMuon();
  if (mu1.isGlobalMuon() != mu2.isGlobalMuon()) return mu1.isGlobalMuon();
  if (mu1.charge() == mu2.charge() && deltaR2(mu1,mu2) < 0.0009) {
    return mu1.track()->ptError()/mu1.track()->pt() < mu2.track()->ptError()/mu2.track()->pt();
  } else {
    int nm1 = mu1.numberOfMatches(reco::Muon::SegmentArbitration);
    int nm2 = mu2.numberOfMatches(reco::Muon::SegmentArbitration);
    return (nm1 != nm2 ? nm1 > nm2 : mu1.pt() > mu2.pt());
  }
}

bool MuMuKKPAT::isSameMuon(const reco::Muon &mu1, const reco::Muon &mu2) const {
  return (& mu1 == & mu2) ||
    //(mu1.originalObjectRef() == mu2.originalObjectRef()) ||
    (mu1.reco::Muon::innerTrack().isNonnull() ?
     mu1.reco::Muon::innerTrack() == mu2.reco::Muon::innerTrack() :
     mu1.reco::Muon::outerTrack() == mu2.reco::Muon::outerTrack());
}


/// define this as a plug-in
DEFINE_FWK_MODULE(MuMuKKPAT);

// rsync -vut --existing src/MuMuPiKPAT.cc semrat@lxplus.cern.ch:/afs/cern.ch/user/s/semrat/scratch0/CMSSW_5_3_22/src/X4140/MuMuKKPAT/src/MuMuPiKPAT.cc

