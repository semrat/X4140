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
  MCMother( iConfig.getUntrackedParameter<int>("MonteCarloMotherId", 511) ), /// 511 B0 (=anti-B0), 531 B0 / decide later MCMotherId for X(4140)
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
  JPsiPhiMinB0Mass(iConfig.getUntrackedParameter<double>("MinJPsiPhiB0Mass", 5.1)),
  JPsiPhiMaxB0Mass(iConfig.getUntrackedParameter<double>("MaxJPsiPhiB0Mass", 5.6)),
  MuMuTrackMaxDR(iConfig.getUntrackedParameter<double>("MaxMuMuTrackDR", 1)),

  B0TrackMaxDR(iConfig.getUntrackedParameter<double>("MaxB0CandTrackDR", 1.1)),
  UseB0DR(iConfig.getUntrackedParameter<bool>("UseB0DR", false)),
  MuMuKKMinB0Mass(iConfig.getUntrackedParameter<double>("MinMuMuKKB0Mass", 0)),
  MuMuKKMaxB0Mass(iConfig.getUntrackedParameter<double>("MaxMuMuKKB0Mass", 10)),
  MuMuKKMaxXMass(iConfig.getUntrackedParameter<double>("MaxMuMuKKXMass", 10)),
  addB0lessPrimaryVertex_(iConfig.getUntrackedParameter<bool>("addB0lessPrimaryVertex", true)),

  Debug_(iConfig.getUntrackedParameter<bool>("Debug_Output",true)),
  DeDxEstimator_(iConfig.getUntrackedParameter<std::string>("DeDxEstimator", std::string("dedxHarmonic2"))),
  m_dEdxDiscrimTag(iConfig.getUntrackedParameter<std::string>("DeDxDiscriminator", std::string("dedxHarmonic2"))),
  m_dEdxDiscrimTag_kaon(iConfig.getUntrackedParameter<std::string>("DeDxDiscriminator", std::string("dedxHarmonic2"))),

  X_One_Tree_(0),
  runNum(0), evtNum(0), lumiNum(0),
  trigRes(0), trigNames(0), L1TT(0), MatchTriggerNames(0),
  /// counters for B0 
  nMu(0), nMuMu(0), nB0(0), nKK(0),
  nB0_pre0(0), nB0_pre1(0), nB0_pre2(0), nB0_pre3(0), nB0_pre4(0), nB0_pre5(0), nB0_pre6(0), nB0_pre7(0), nB0_pre8(0), nB0_pre9(0), nB0_pre10(0), nB0_pre11(0), nB0_pre12(0), nB0_pre13(0), nB0_pre14(0), nB0_pre15(0), 
  //nX(0),
  priVtx_n(0), priVtx_X(0), priVtx_Y(0), priVtx_Z(0), priVtx_XE(0), priVtx_YE(0), priVtx_ZE(0), priVtx_NormChi2(0), priVtx_Chi2(0), priVtx_CL(0), priVtx_tracks(0), priVtx_tracksPtSq(0),
  /// indices
  mu1Idx(0), mu2Idx(0), MuMuType(0), ka1Idx(0), ka2Idx(0),
  B0_MuMuIdx(0), B0_ka1Idx(0), B0_ka2Idx(0),
  /// MC Analysis /// n_B0Ancestors & no for X(4140)!!
  n_genEvtVtx(0), genEvtVtx_X(0), genEvtVtx_Y(0), genEvtVtx_Z(0), genEvtVtx_particles(0), n_B0Ancestors(0),
  nMCAll(0), nMCB0(0), /*nMCB0Vtx(0),*/ MCPdgIdAll(0), MCDanNumAll(0),
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
  DRMuMuK1(0), DRMuMuK2(0), DRb0K1(0), DRb0K2(0), 
  /// Primary Vertex with "MuMu correction"
  PriVtxMuMuCorr_n(0),
  PriVtxMuMuCorr_X(0), PriVtxMuMuCorr_Y(0), PriVtxMuMuCorr_Z(0), PriVtxMuMuCorr_EX(0), PriVtxMuMuCorr_EY(0), PriVtxMuMuCorr_EZ(0),
  PriVtxMuMuCorr_Chi2(0), PriVtxMuMuCorr_CL(0), PriVtxMuMuCorr_tracks(0),
  nTrk(0),
  /// B0 cand
  b0Mass(0), b0Vtx_CL(0), b0Vtx_Chi2(0),
  b0Px(0), b0Py(0), b0Pz(0), b0PxE(0), b0PyE(0), b0PzE(0),
  b0DecayVtx_X(0), b0DecayVtx_Y(0), b0DecayVtx_Z(0), b0DecayVtx_XE(0), b0DecayVtx_YE(0), b0DecayVtx_ZE(0),
  /// Muons and tracks after B0 cand fit 
  mu1Px_MuMuKK(0), mu1Py_MuMuKK(0), mu1Pz_MuMuKK(0), mu1E_MuMuKK(0),
  mu2Px_MuMuKK(0), mu2Py_MuMuKK(0), mu2Pz_MuMuKK(0), mu2E_MuMuKK(0),
  k1Px_MuMuKK(0), k1Py_MuMuKK(0), k1Pz_MuMuKK(0), k1E_MuMuKK(0),
  kaon1_nsigdedx(0), kaon1_dedx(0), kaon1_dedxMass(0), kaon1_theo(0), kaon1_sigma(0),
  kaon1_dedx_byHits(0), kaon1_dedxErr_byHits(0), kaon1_saturMeas_byHits(0), kaon1_Meas_byHits(0),
  k2Px_MuMuKK(0), k2Py_MuMuKK(0), k2Pz_MuMuKK(0), k2E_MuMuKK(0),
  kaon2_nsigdedx(0), kaon2_dedx(0), kaon2_dedxMass(0), kaon2_theo(0), kaon2_sigma(0),
  kaon2_dedx_byHits(0), kaon2_dedxErr_byHits(0), kaon2_saturMeas_byHits(0), kaon2_Meas_byHits(0),
  /// Primary Vertex with largest B0_cos(alpha) no less values for X(4140)
  PriVtx_B0CosAlpha_n(0),
  PriVtx_B0CosAlpha_X(0), PriVtx_B0CosAlpha_Y(0), PriVtx_B0CosAlpha_Z(0), PriVtx_B0CosAlpha_EX(0), PriVtx_B0CosAlpha_EY(0), PriVtx_B0CosAlpha_EZ(0),
  PriVtx_B0CosAlpha_Chi2(0), PriVtx_B0CosAlpha_CL(0), PriVtx_B0CosAlpha_tracks(0),
  PriVtx_B0CosAlpha3D_n(0),
  PriVtx_B0CosAlpha3D_X(0), PriVtx_B0CosAlpha3D_Y(0), PriVtx_B0CosAlpha3D_Z(0), PriVtx_B0CosAlpha3D_EX(0), PriVtx_B0CosAlpha3D_EY(0), PriVtx_B0CosAlpha3D_EZ(0),
  PriVtx_B0CosAlpha3D_Chi2(0), PriVtx_B0CosAlpha3D_CL(0), PriVtx_B0CosAlpha3D_tracks(0),
  B0LessPV_tracksPtSq(0), B0LessPV_4tracksPtSq(0),
  PriVtxB0Less_n(0),
  PriVtxB0Less_X(0), PriVtxB0Less_Y(0), PriVtxB0Less_Z(0), PriVtxB0Less_EX(0), PriVtxB0Less_EY(0), PriVtxB0Less_EZ(0),
  PriVtxB0Less_Chi2(0), PriVtxB0Less_CL(0), PriVtxB0Less_tracks(0),
  PriVtxB0Less_B0CosAlpha_n(0),
  PriVtxB0Less_B0CosAlpha_X(0), PriVtxB0Less_B0CosAlpha_Y(0), PriVtxB0Less_B0CosAlpha_Z(0), PriVtxB0Less_B0CosAlpha_EX(0), PriVtxB0Less_B0CosAlpha_EY(0), PriVtxB0Less_B0CosAlpha_EZ(0),
  PriVtxB0Less_B0CosAlpha_Chi2(0), PriVtxB0Less_B0CosAlpha_CL(0), PriVtxB0Less_B0CosAlpha_tracks(0),
  PriVtxB0Less_B0CosAlpha3D_n(0),
  PriVtxB0Less_B0CosAlpha3D_X(0), PriVtxB0Less_B0CosAlpha3D_Y(0), PriVtxB0Less_B0CosAlpha3D_Z(0), PriVtxB0Less_B0CosAlpha3D_EX(0), PriVtxB0Less_B0CosAlpha3D_EY(0), PriVtxB0Less_B0CosAlpha3D_EZ(0),
  PriVtxB0Less_B0CosAlpha3D_Chi2(0), PriVtxB0Less_B0CosAlpha3D_CL(0), PriVtxB0Less_B0CosAlpha3D_tracks(0),
  /// Primary Vertex with "B0 correction" 
  PriVtxB0Corr_n(0),
  PriVtxB0Corr_X(0), PriVtxB0Corr_Y(0), PriVtxB0Corr_Z(0), PriVtxB0Corr_EX(0), PriVtxB0Corr_EY(0), PriVtxB0Corr_EZ(0),
  PriVtxB0Corr_Chi2(0), PriVtxB0Corr_CL(0), PriVtxB0Corr_tracks(0),
  /// Lifetime variables for B0 
  b0CosAlphaBS(0), b0CosAlpha3DBS(0), b0CTauBS(0), b0CTauBSE(0), b0LxyBS(0), b0LxyBSE(0), b0LxyzBS(0), b0LxyzBSE(0),
  b0CosAlphaPV(0), b0CosAlpha3DPV(0), b0CTauPV(0), b0CTauPVE(0), b0LxyPV(0), b0LxyPVE(0), b0LxyzPV(0), b0LxyzPVE(0),
  b0CosAlphaPVCosAlpha(0), b0CosAlpha3DPVCosAlpha(0), b0CTauPVCosAlpha(0), b0CTauPVCosAlphaE(0), b0LxyPVCosAlpha(0), b0LxyPVCosAlphaE(0), b0LxyzPVCosAlpha(0), b0LxyzPVCosAlphaE(0),
  b0CosAlphaPVCosAlpha3D(0), b0CosAlpha3DPVCosAlpha3D(0), b0CTauPVCosAlpha3D(0), b0CTauPVCosAlpha3DE(0), b0LxyPVCosAlpha3D(0), b0LxyPVCosAlpha3DE(0), b0LxyzPVCosAlpha3D(0), b0LxyzPVCosAlpha3DE(0),
  b0CosAlphaB0LessPV(0), b0CosAlpha3DB0LessPV(0),b0CTauB0LessPV(0),
  b0CTauB0LessPVE(0), b0LxyB0LessPV(0), b0LxyB0LessPVE(0), b0LxyzB0LessPV(0), b0LxyzB0LessPVE(0),
  b0CosAlphaB0LessPVCosAlpha(0), b0CosAlpha3DB0LessPVCosAlpha(0), b0CTauB0LessPVCosAlpha(0), b0CTauB0LessPVCosAlphaE(0), b0LxyB0LessPVCosAlpha(0), b0LxyB0LessPVCosAlphaE(0), b0LxyzB0LessPVCosAlpha(0), b0LxyzB0LessPVCosAlphaE(0),
  b0CosAlphaB0LessPVCosAlpha3D(0), b0CosAlpha3DB0LessPVCosAlpha3D(0), b0CTauB0LessPVCosAlpha3D(0), b0CTauB0LessPVCosAlpha3DE(0), b0LxyB0LessPVCosAlpha3D(0), b0LxyB0LessPVCosAlpha3DE(0), b0LxyzB0LessPVCosAlpha3D(0), b0LxyzB0LessPVCosAlpha3DE(0),
  b0CosAlphaPVX(0), b0CTauPVX(0), b0CTauPVXE(0), b0LxyPVX(0), b0LxyzPVX(0), b0LxyzPVXE(0),
  b0CTauPVX_3D(0), b0CTauPVX_3D_err(0),
  /// dxy, dz, dxyE, dzE for kaons from PV, BS, B0LessPV
  kaon1_dxy_PV(0), kaon1_dz_PV(0), kaon2_dxy_PV(0), kaon2_dz_PV(0),
  kaon1_dxy_BS(0), kaon1_dz_BS(0), kaon2_dxy_BS(0), kaon2_dz_BS(0),
  kaon1_dxy_B0LessPV(0), kaon1_dz_B0LessPV(0), kaon2_dxy_B0LessPV(0), kaon2_dz_B0LessPV(0),
  kaon1_dxyE(0), kaon1_dzE(0), kaon2_dxyE(0), kaon2_dzE(0),

  KKMass_err(0), Kaon1FromPV(0), Kaon2FromPV(0)

{
  /// now do what ever initialization is needed
  MuMuMinMass = JPsiMinMass;
  MuMuMaxMass = JPsiMaxMass;
  KKMinMass = PhiMinMass;
  KKMaxMass = PhiMaxMass;
  MuMuKKMinB0Mass = JPsiPhiMinB0Mass;
  MuMuKKMaxB0Mass = JPsiPhiMaxB0Mass;
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

    if (addMuMulessPrimaryVertex_ || addB0lessPrimaryVertex_ || resolveAmbiguity_) {
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
	      pionOK = true; /// SEMRA pions change with kaons for B0 ?
	    } else if ( abs(dau->pdgId()) == 321 ) { // check if one of B0 daughters is a kaon /// SEMRA ask again !!!
	      kaonPx = dau->px(); kaonPy=dau->py(); kaonPz=dau->pz();
	      kaonCh = (dau->pdgId() == 321)? 1 : -1;
	      kaonOK = true;
	    }*/

	  } /// end loop on MCMother daughters

	  if (Debug_) cout << "mumuOK = " << mumuOK << ", kkOK = " << kkOK << endl;
	  if ( mumuOK && kkOK ) {
	    if (Debug_) {
	      cout <<"\nnumber of B0 mothers = " <<p.numberOfMothers() <<endl ;
	      cout <<"B0 mother pdgID = " <<p.mother(0)->pdgId() <<endl ;
	    }
	    ++nMCB0 ;
	      PriVtxGen_X->push_back( p.vx() ) ;
	      PriVtxGen_Y->push_back( p.vy() ) ;
	      PriVtxGen_Z->push_back( p.vz() ) ;
	      PriVtxGen_CL->push_back( p.vertexNormalizedChi2() ) ;
	      PriVtxGen_Chi2->push_back( p.vertexChi2() ) ;
	      PriVtxGen_Ndof->push_back( p.vertexNdof() ) ;

	      Bool_t status = kTRUE ;
	      const Candidate *b0_ancestor = p.mother(0) ; /// a particle can have several mothers
	      Int_t n_ancestors = 1 ;
	      while ( status ) {
		if ( abs(b0_ancestor->pdgId()) <= 8 || b0_ancestor->pdgId() == 21 || b0_ancestor->status() == 3 ) {
		  status = kFALSE ;
		  if (Debug_) cout <<"B0 ancestor ID = " <<b0_ancestor->pdgId() <<endl ;
		  genEvtVtx_X->push_back( b0_ancestor->daughter(0)->vx() ) ;
		  genEvtVtx_Y->push_back( b0_ancestor->daughter(0)->vy() ) ;
		  genEvtVtx_Z->push_back( b0_ancestor->daughter(0)->vz() ) ;
		  genEvtVtx_particles->push_back( b0_ancestor->numberOfDaughters() ) ;
		  n_B0Ancestors->push_back( n_ancestors ) ;
		}
		else {
		  b0_ancestor = b0_ancestor->mother(0) ;
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
            /// for B0 
	    if (Debug_) cout <<"evt:"<<evtNum<< " is Invalid Muon ?  " <<isEventWithInvalidMu << endl;
	    //if (skipJPsi && ( dimuonType == 1 ));
            //cout<< "POINT 11" <<endl;
	    nTrk->push_back( thePATTrackHandle->size() ) ;
	    //cout<< "POINT 12" <<endl;
            if (thePATTrackHandle->size() < 2) {
            //cout<< "POINT 13" <<endl;
             nB0_pre0++;
	    }
            if (Debug_) cout<<"nmumu : "<<nMuMu<<endl;


	    ////////////////// cuts on MuMu mass window for B0 ////////////////////////////
	    if (MuMuMass->at(nMuMu-1) < MuMuMinMass  ||  MuMuMass->at(nMuMu-1) > MuMuMaxMass){
	       continue ; nB0_pre1++ ; 
            }
	    //cout<< "POINT 14" <<endl;


	    ////////////////// check tracks for kaon1 for B0 //////////////////
	    for ( vector<pat::GenericParticle>::const_iterator Track1 = theKaonRefittedPATTrackHandle->begin(); Track1 != theKaonRefittedPATTrackHandle->end(); ++Track1 ) {
	      //cout<< "POINT 15" <<endl;
	      /// check track doesn't overlap with the MuMu candidate tracks
	      if (Track1->track().key() == rmu1->track().key()  ||  Track1->track().key() == rmu2->track().key())
		 continue ; nB0_pre2++ ; 

	      //cout<< "POINT 16" <<endl;
	      /// cuts on charged tracks
	      if (( Track1->track()->chi2()/Track1->track()->ndof() > TrMaxNormChi2 )  ||  Track1->pt() < TrMinPt)
		continue ; nB0_pre3++ ; 

	      //cout<< "POINT 17" <<endl;

	   ////////////////// check tracks for kaon2 for B0 //////////////////
	   for ( vector<pat::GenericParticle>::const_iterator Track2 = Track1+1; Track2 != theKaonRefittedPATTrackHandle->end(); ++Track2 ){
	     /// check that this second track doesn't overlap with the the first track candidate
	     if (Track2->track().key() == Track1->track().key())
	        continue ; nB0_pre4++ ; 

            /// check track doesn't overlap with the MuMu candidate tracks
            if (Track2->track().key() == rmu1->track().key()  ||  Track2->track().key() == rmu2->track().key())
              continue ; nB0_pre5++ ; 
            if (Track1->charge() * Track2->charge() > 0)
	      continue ; nB0_pre6++ ; 
	    /// cuts on charged tracks
            if ((Track2->track()->chi2() / Track2->track()->ndof() > TrMaxNormChi2)  ||  Track2->pt() < TrMinPt)
	      continue; nB0_pre7++ ; 


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

            ////////////////// cuts on tracks' delta R for B0 //////////////////
            math::XYZTLorentzVector MuMu = (rmu1->p4() + rmu2->p4());
	    math::XYZTLorentzVector b0 = (rmu1->p4() + rmu2->p4() + Track1->p4() + Track2->p4());
	    float MuMuK1DR = sqrt( pow(MuMu.eta() - Track1->p4().eta(),2) + pow(MuMu.phi() - Track1->p4().phi(), 2) );
	    float MuMuK2DR = sqrt( pow(MuMu.eta() - Track2->p4().eta(),2) + pow(MuMu.phi() - Track2->p4().phi(), 2) );
            float b0K1DR = sqrt( pow(b0.eta() - Track1->p4().eta(),2) + pow(b0.phi() - Track1->p4().phi(), 2) );
            float b0K2DR = sqrt( pow(b0.eta() - Track2->p4().eta(),2) + pow(b0.phi() - Track2->p4().phi(), 2) );

	    DRMuMuK1->push_back(MuMuK1DR);
	    DRMuMuK2->push_back(MuMuK2DR);
            DRb0K1->push_back(b0K1DR);
            DRb0K2->push_back(b0K2DR);


		  if (UseB0DR) {
		    if (b0K1DR > B0TrackMaxDR || b0K2DR > B0TrackMaxDR)
		      B0TrackMaxDR = 2;
		  } else {
		    if (MuMuK1DR > MuMuTrackMaxDR || MuMuK2DR > MuMuTrackMaxDR)
		      MuMuTrackMaxDR = 3.5;
		  }
		  nB0_pre8++ ;


	    ////////////////// cuts on MuMuKK mass window for B0 //////////////////
            if (((Track1->p4() + Track2->p4() + MuMu).M() > MuMuKKMaxB0Mass  ||  (Track1->p4() + Track2->p4() + MuMu).M() < MuMuKKMinB0Mass) && ((Track1->p4() + Track2->p4() + MuMu).M() >  MuMuKKMaxXMass))
	       continue ; nB0_pre9++ ; 



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

		  /// do mass constraint for MuMu cand and do mass constrained vertex fit for B0
		  vector<RefCountedKinematicParticle> b0Daughters;
		  b0Daughters.push_back(pFactory.particle( muon1TT, muon_mass, chi, ndf, small_sigma));
		  b0Daughters.push_back(pFactory.particle( muon2TT, muon_mass, chi, ndf, small_sigma));
		  b0Daughters.push_back(pFactory.particle( kaon1TT, kaon_mass, chi, ndf, small_sigma));
		  b0Daughters.push_back(pFactory.particle( kaon2TT, kaon_mass, chi, ndf, small_sigma));

		  RefCountedKinematicTree B0VertexFitTree, B0VertexFitTree_noKrefit ;
		  KinematicConstrainedVertexFitter B0Fitter ;

		  if (doMuMuMassConst) { // MassConst = 'MC' in the following
		    MultiTrackKinematicConstraint *MuMu = 0;
		    if (dimuonType == 1) { // constrain to JPsi mass
		      MuMu = new TwoTrackMassKinematicConstraint(JPsi_mass);
		    } //else if (dimuonType == 2) { // constrain to Psi(2S) mass /// SEMRA will we use this or not ?
		      //MuMu = new TwoTrackMassKinematicConstraint(psi2S_mass);
		    //} // already asked for: if (dimuonType == 0) continue ;

		    B0VertexFitTree = B0Fitter.fit( b0Daughters, MuMu );
		    if (notRefittedPartner) { // use not refitted kaons
		      b0Daughters.pop_back() ;
		      b0Daughters.push_back(pFactory.particle( kaon2TT_notRefit, kaon_mass, chi, ndf, small_sigma));
		      B0VertexFitTree_noKrefit = B0Fitter.fit( b0Daughters, MuMu );
		    }
		  }
		  else {
		    B0VertexFitTree = B0Fitter.fit( b0Daughters );
		    if (notRefittedPartner) { // use not refitted kaons
		      b0Daughters.pop_back() ;
		      b0Daughters.push_back(pFactory.particle( kaon2TT_notRefit, kaon_mass, chi, ndf, small_sigma));
		      B0VertexFitTree_noKrefit = B0Fitter.fit( b0Daughters );
		    }
		  }

		  if ( !B0VertexFitTree->isValid() ) /// B0 variables started 
                    continue ; nB0_pre10++ ;

		  B0VertexFitTree->movePointerToTheTop();
		  RefCountedKinematicParticle B0Cand_fromMCFit = B0VertexFitTree->currentParticle();
		  RefCountedKinematicVertex B0Cand_vertex_fromMCFit = B0VertexFitTree->currentDecayVertex();

		  if ( !B0Cand_vertex_fromMCFit->vertexIsValid() )
	             continue ; nB0_pre11++ ;

		  if ( B0Cand_vertex_fromMCFit->chiSquared() < 0  ||  B0Cand_vertex_fromMCFit->chiSquared() > 10000 )
		     continue ; nB0_pre12++ ;

	          if (B0Cand_vertex_fromMCFit->chiSquared() / B0Cand_vertex_fromMCFit->degreesOfFreedom() > 7 )
                     continue ; nB0_pre13++;

		  if ( B0Cand_fromMCFit->currentState().mass() > 100 ) 
                     continue ; nB0_pre14++ ;

		  double b0VtxProb = ChiSquaredProbability((double)(B0Cand_vertex_fromMCFit->chiSquared()), (double)(B0Cand_vertex_fromMCFit->degreesOfFreedom()));
		  if ( b0VtxProb < 0.005 ) //0.0001 )
		     continue ; nB0_pre15++ ;

		  
                  //////////////////// Lifetimes calculations for B0 ////////////////////
		  TVector3 B0_vtx((*B0Cand_vertex_fromMCFit).position().x(), (*B0Cand_vertex_fromMCFit).position().y(), 0) ; 
		  TVector3 B0_pperp(B0Cand_fromMCFit->currentState().globalMomentum().x(), B0Cand_fromMCFit->currentState().globalMomentum().y(), 0);
		  TVector3 B0_vtx3D((*B0Cand_vertex_fromMCFit).position().x(), (*B0Cand_vertex_fromMCFit).position().y(), (*B0Cand_vertex_fromMCFit).position().z()) ;
		  TVector3 B0_pperp3D(B0Cand_fromMCFit->currentState().globalMomentum().x(), B0Cand_fromMCFit->currentState().globalMomentum().y(), B0Cand_fromMCFit->currentState().globalMomentum().z());

                  AlgebraicVector3 B0_v3pperp ;
		  B0_v3pperp[0] = B0_pperp.x(); B0_v3pperp[1] = B0_pperp.y(); B0_v3pperp[2] = 0.;
		  TVector3 B0_pvtx, B0_pvtx3D, B0_vdiff, B0_vdiff3D ;
		  double B0_cosAlpha, B0_cosAlpha3D, B0_ctau ;
		  VertexDistanceXY B0_vdistXY ;
		  Measurement1D B0_distXY ;
		  GlobalError B0_v1e = (Vertex(*B0Cand_vertex_fromMCFit)).error();
		  GlobalError B0_v2e ;
		  AlgebraicSymMatrix33 B0_vXYe ;
		  double B0_ctauErr ;
		  float B0_lxy, B0_lxyErr, B0_lxyz, B0_lxyzErr ;
		  ROOT::Math::SVector<double, 3> B0_vDiff, B0_vDiff3D ; // needed by Similarity method


		  ////////////////// Lifetime wrt PV for B0 //////////////////
		  B0_v2e = thePrimaryVtx.error(); 
		  B0_vXYe = B0_v1e.matrix() + B0_v2e.matrix() ;
		  /// 2D
		  B0_pvtx.SetXYZ(thePrimaryVtx.position().x(), thePrimaryVtx.position().y(), 0) ;
		  B0_vdiff = B0_vtx - B0_pvtx ;
		  B0_cosAlpha = B0_vdiff.Dot(B0_pperp) / (B0_vdiff.Perp()*B0_pperp.Perp()) ;
		  B0_lxy = B0_vdiff.Perp();
		  B0_vDiff[0] = B0_vdiff.x(); B0_vDiff[1] = B0_vdiff.y(); B0_vDiff[2] = 0 ; // needed by Similarity method
		  B0_lxyErr = sqrt(ROOT::Math::Similarity(B0_vDiff,B0_vXYe)) / B0_vdiff.Perp();
		  B0_distXY = B0_vdistXY.distance(Vertex(*B0Cand_vertex_fromMCFit), Vertex(thePrimaryVtx));
		  B0_ctau = B0_distXY.value() * B0_cosAlpha * B0Cand_fromMCFit->currentState().mass() / B0_pperp.Perp();
		  B0_ctauErr = sqrt(ROOT::Math::Similarity(B0_v3pperp,B0_vXYe)) * B0Cand_fromMCFit->currentState().mass() / (B0_pperp.Perp2()) ;
		  /// 3D
		  B0_pvtx3D.SetXYZ(thePrimaryVtx.position().x(), thePrimaryVtx.position().y(), thePrimaryVtx.position().z());
		  B0_vdiff3D = B0_vtx3D - B0_pvtx3D;
		  B0_cosAlpha3D = B0_vdiff3D.Dot(B0_pperp3D)/(B0_vdiff3D.Mag()*B0_pperp3D.Mag());
		  B0_lxyz = B0_vdiff3D.Mag();
		  B0_vDiff3D[0] = B0_vdiff3D.x(); B0_vDiff3D[1] = B0_vdiff3D.y(); B0_vDiff3D[2] = B0_vdiff3D.z() ;
		  B0_lxyzErr = sqrt(ROOT::Math::Similarity(B0_vDiff3D,B0_vXYe)) / B0_vdiff3D.Mag();
		    
                 
                  ////////////////// Last cuts for B0 //////////////////
		  if ( !(B0_ctau/B0_ctauErr > 2.8) || !(B0_cosAlpha > 0.8) )
		    continue ;
		  		  

                  ////////////////// fill B0 candidate variables //////////////////
		  b0Mass->push_back( B0Cand_fromMCFit->currentState().mass()) ; 
		  b0Px->push_back( B0Cand_fromMCFit->currentState().globalMomentum().x()) ;
		  b0Py->push_back( B0Cand_fromMCFit->currentState().globalMomentum().y()) ;
		  b0Pz->push_back( B0Cand_fromMCFit->currentState().globalMomentum().z()) ;
		  b0PxE->push_back( sqrt( B0Cand_fromMCFit->currentState().kinematicParametersError().matrix()(3,3) ) ) ;
		  b0PyE->push_back( sqrt( B0Cand_fromMCFit->currentState().kinematicParametersError().matrix()(4,4) ) ) ;
		  b0PzE->push_back( sqrt( B0Cand_fromMCFit->currentState().kinematicParametersError().matrix()(5,5) ) ) ;
		  b0Vtx_CL->push_back( b0VtxProb );
		  b0Vtx_Chi2->push_back( B0Cand_vertex_fromMCFit->chiSquared() ) ;
		  b0DecayVtx_X->push_back((*B0Cand_vertex_fromMCFit).position().x());
		  b0DecayVtx_Y->push_back((*B0Cand_vertex_fromMCFit).position().y());
		  b0DecayVtx_Z->push_back((*B0Cand_vertex_fromMCFit).position().z());
		  b0DecayVtx_XE->push_back(sqrt((*B0Cand_vertex_fromMCFit).error().cxx()));
		  b0DecayVtx_YE->push_back(sqrt((*B0Cand_vertex_fromMCFit).error().cyy()));
		  b0DecayVtx_ZE->push_back(sqrt((*B0Cand_vertex_fromMCFit).error().czz()));
		  B0VertexFitTree->movePointerToTheFirstChild();
		  RefCountedKinematicParticle mu1_MuMuKK = B0VertexFitTree->currentParticle();
		  B0VertexFitTree->movePointerToTheNextChild();
		  RefCountedKinematicParticle mu2_MuMuKK = B0VertexFitTree->currentParticle();
		  B0VertexFitTree->movePointerToTheNextChild();
		  RefCountedKinematicParticle k1_MuMuKK = B0VertexFitTree->currentParticle();
		  B0VertexFitTree->movePointerToTheNextChild();
		  RefCountedKinematicParticle k2_MuMuKK = B0VertexFitTree->currentParticle();
		  /// muon1 & muon2
		  mu1Px_MuMuKK->push_back( mu1_MuMuKK->currentState().globalMomentum().x() );
		  mu1Py_MuMuKK->push_back( mu1_MuMuKK->currentState().globalMomentum().y() );
		  mu1Pz_MuMuKK->push_back( mu1_MuMuKK->currentState().globalMomentum().z() );
		  mu1E_MuMuKK->push_back( mu1_MuMuKK->currentState().kinematicParameters().energy() );
		  mu2Px_MuMuKK->push_back( mu2_MuMuKK->currentState().globalMomentum().x() );
		  mu2Py_MuMuKK->push_back( mu2_MuMuKK->currentState().globalMomentum().y() );
		  mu2Pz_MuMuKK->push_back( mu2_MuMuKK->currentState().globalMomentum().z() );
		  mu2E_MuMuKK->push_back( mu2_MuMuKK->currentState().kinematicParameters().energy() );
		  /// kaon1 & kaon2
		  k1Px_MuMuKK->push_back( k1_MuMuKK->currentState().globalMomentum().x() );
		  k1Py_MuMuKK->push_back( k1_MuMuKK->currentState().globalMomentum().y() );
		  k1Pz_MuMuKK->push_back( k1_MuMuKK->currentState().globalMomentum().z() );
		  k1E_MuMuKK->push_back( k1_MuMuKK->currentState().kinematicParameters().energy() );
		  Double_t theo = 0., sigma = 0. ;
		  kaon1_nsigdedx->push_back( nsigmaofdedx(Track1->track(),theo,sigma) );
		  kaon1_dedx->push_back( getEnergyLoss(Track1->track()) );
		  kaon1_dedxMass->push_back( GetMass(Track1->track()) );
		  kaon1_theo->push_back( theo );
		  kaon1_sigma->push_back( sigma );
		  kaon1_dedx_byHits->push_back( (dEdxTrack)[Track1->track()].dEdx() );
		  kaon1_dedxErr_byHits->push_back( (dEdxTrack)[Track1->track()].dEdxError() );
		  kaon1_saturMeas_byHits->push_back( (dEdxTrack)[Track1->track()].numberOfSaturatedMeasurements() );
		  kaon1_Meas_byHits->push_back( (dEdxTrack)[Track1->track()].numberOfMeasurements() );
		  k2Px_MuMuKK->push_back( k2_MuMuKK->currentState().globalMomentum().x() );
		  k2Py_MuMuKK->push_back( k2_MuMuKK->currentState().globalMomentum().y() );
		  k2Pz_MuMuKK->push_back( k2_MuMuKK->currentState().globalMomentum().z() );
		  k2E_MuMuKK->push_back( k2_MuMuKK->currentState().kinematicParameters().energy() );
		  theo = 0.; sigma = 0. ;
		  kaon2_nsigdedx->push_back(nsigmaofdedx(Track2->track(),theo,sigma));
		  kaon2_dedx->push_back(getEnergyLoss(Track2->track()));
		  kaon2_dedxMass->push_back(GetMass(Track2->track()));
		  kaon2_theo->push_back(theo);
		  kaon2_sigma->push_back(sigma);
		  kaon2_dedx_byHits->push_back( (dEdxTrack_Kaon)[Track2->track()].dEdx() );
		  kaon2_dedxErr_byHits->push_back( (dEdxTrack_Kaon)[Track2->track()].dEdxError() );
		  kaon2_saturMeas_byHits->push_back( (dEdxTrack_Kaon)[Track2->track()].numberOfSaturatedMeasurements() );
		  kaon2_Meas_byHits->push_back( (dEdxTrack_Kaon)[Track2->track()].numberOfMeasurements() );
		  /// PV
		  b0CosAlphaPV->push_back( B0_cosAlpha ); b0CosAlpha3DPV->push_back( B0_cosAlpha3D );
		  b0CTauPV->push_back( B0_ctau ); b0CTauPVE->push_back( B0_ctauErr );
		  b0LxyPV->push_back( B0_lxy ); b0LxyPVE->push_back( B0_lxyErr );
		  b0LxyzPV->push_back( B0_lxyz ); b0LxyzPVE->push_back( B0_lxyzErr );
		  /// dxy, dz, dxyE, dzE for kaons from PV 
		  kaon1_dxy_PV->push_back( Track1->track()->dxy(RefVtx) );
                  kaon1_dz_PV->push_back( Track1->track()->dz(RefVtx) );
                  kaon2_dxy_PV->push_back( Track2->track()->dxy(RefVtx) );
                  kaon2_dz_PV->push_back( Track2->track()->dz(RefVtx) );		    


                  ////////////////// Lifetime wrt BS for B0 //////////////////
		  B0_v2e = theBeamSpotVtx.error(); 
		  B0_vXYe = B0_v1e.matrix() + B0_v2e.matrix();
		  /// 2D
		  B0_pvtx.SetXYZ(theBeamSpotVtx.position().x(), theBeamSpotVtx.position().y(), 0);
		  B0_vdiff = B0_vtx - B0_pvtx;
		  B0_cosAlpha = B0_vdiff.Dot(B0_pperp)/(B0_vdiff.Perp()*B0_pperp.Perp());
		  B0_lxy = B0_vdiff.Perp();
		  B0_vDiff[0] = B0_vdiff.x(); B0_vDiff[1] = B0_vdiff.y(); B0_vDiff[2] = 0 ; // needed by Similarity method
		  B0_lxyErr = sqrt(ROOT::Math::Similarity(B0_vDiff,B0_vXYe)) / B0_vdiff.Perp();
		  B0_distXY = B0_vdistXY.distance(Vertex(*B0Cand_vertex_fromMCFit), Vertex(theBeamSpotVtx));
		  B0_ctau = B0_distXY.value() * B0_cosAlpha * (B0Cand_fromMCFit->currentState().mass() / B0_pperp.Perp()) ;
		  B0_ctauErr = sqrt(ROOT::Math::Similarity(B0_v3pperp,B0_vXYe)) * B0Cand_fromMCFit->currentState().mass()/B0_pperp.Perp2();
		  /// 3D
		  B0_pvtx3D.SetXYZ(theBeamSpotVtx.position().x(), theBeamSpotVtx.position().y(), theBeamSpotVtx.position().z());
		  B0_vdiff3D = B0_vtx3D - B0_pvtx3D;
		  B0_cosAlpha3D = B0_vdiff3D.Dot(B0_pperp3D)/(B0_vdiff3D.Mag()*B0_pperp3D.Mag());
		  B0_lxyz = B0_vdiff3D.Mag();
		  B0_vDiff3D[0] = B0_vdiff3D.x(); B0_vDiff3D[1] = B0_vdiff3D.y(); B0_vDiff3D[2] = B0_vdiff3D.z() ;
		  B0_lxyzErr = sqrt(ROOT::Math::Similarity(B0_vDiff3D,B0_vXYe)) / B0_vdiff3D.Mag();
		  	  

                  ////////////////// BS (beam spot) for B0 //////////////////
		  b0CosAlphaBS->push_back( B0_cosAlpha ); b0CosAlpha3DBS->push_back( B0_cosAlpha3D );
		  b0CTauBS->push_back( B0_ctau ); b0CTauBSE->push_back( B0_ctauErr );
		  b0LxyBS->push_back( B0_lxy ); b0LxyBSE->push_back( B0_lxyErr );
		  b0LxyzBS->push_back( B0_lxyz ); b0LxyzBSE->push_back( B0_lxyzErr );

		  vector<TransientVertex> B0_pvs ;
		  Vertex B0LessPV = thePrimaryVtx ;

		  if (addB0lessPrimaryVertex_)
		    {
		      VertexReProducer revertex(recVtxs, iEvent);
		      Handle<TrackCollection> pvtracks;
		      iEvent.getByLabel(revertex.inputTracks(), pvtracks);
		      Handle<BeamSpot>        pvbeamspot;
		      iEvent.getByLabel(revertex.inputBeamSpot(), pvbeamspot);

		      if (pvbeamspot.id() != beamSpotHandle.id() )
			edm::LogWarning("Inconsistency") << "The BeamSpot used for PV reco is not the same used in this analyzer.";

		      const reco::Muon *B0rmu_1 = dynamic_cast<const reco::Muon *>(Muon1->originalObject());
		      const reco::Muon *B0rmu_2 = dynamic_cast<const reco::Muon *>(Muon2->originalObject());

		      if (B0rmu_1 != 0  &&  B0rmu_2 != 0  &&  B0rmu_1->track().id() == pvtracks.id()  &&  B0rmu_2->track().id() == pvtracks.id()
			  &&  Track1->track().id() == pvtracks.id()  &&  Track2->track().id() ==  pvtracks.id()) {
			vector<TransientTrack> B0Less; // need TransientTrack to keep the TrackRef
			B0Less.reserve( pvtracks->size() );
			Double_t removedTrksPtSq = 0. ;
			for (size_t i = 0, n = pvtracks->size(); i < n; ++i) {
			  if (i == B0rmu_1->track().key()) { removedTrksPtSq += (B0rmu_1->track()->pt())*(B0rmu_1->track()->pt()) ;
			    continue; }
			  if (i == B0rmu_2->track().key()) { removedTrksPtSq += (B0rmu_2->track()->pt())*(B0rmu_2->track()->pt()) ;
			    continue; }
			  if (i == Track1->track().key()) { removedTrksPtSq += (Track1->track()->pt())*(Track1->track()->pt()) ;
			    continue; }
			  if (i == Track2->track().key()) { removedTrksPtSq += (Track2->track()->pt())*(Track2->track()->pt()) ;
			    continue; }

			  reco::TrackRef trk_now(pvtracks, i) ;
			  TransientTrack transientTrack = theTTBuilder->build( trk_now );
			  transientTrack.setBeamSpot( beamSpot );
			  B0Less.push_back( transientTrack );
			}
			if ( removedTrksPtSq > 0. ) {
			  B0_pvs = revertex.makeVertices(B0Less, *pvbeamspot, iSetup) ; // list of PV
			} else
			  if (Debug_) cout <<"\n\\\\\\\\\\\\\\\\\\\\ excluded tracks pT^2 = 0 \\\\\\\\\\\\\\\\\\\\\n" <<endl ;
			if ( !B0_pvs.empty() ) {
			  B0LessPV = Vertex(B0_pvs.front());
			  B0LessPV_tracksPtSq->push_back( vertexHigherPtSquared.sumPtSquared(B0LessPV) ) ;
			  B0LessPV_4tracksPtSq->push_back( removedTrksPtSq ) ;
			  if (Debug_) {
			    cout <<"\nB0LessPV_z = " <<B0LessPV.position().z() <<endl ;
			    cout <<"B0LessPV_tracks = " <<B0LessPV.tracksSize() <<endl ;
			    cout <<"B0LessPV_tracksPtSq = " <<vertexHigherPtSquared.sumPtSquared(B0LessPV) <<endl ;
			    cout <<"B0LessPV_removedTracksPtSq = " <<removedTrksPtSq <<endl ;
			    cout <<"B0_pvs->size() = " <<B0_pvs.size() <<endl ;
			    cout <<"priVtx_z = " << priVtx_Z <<endl ;
			    cout <<"priVtx_tracks = " <<priVtx_tracks <<endl ;
			    cout <<"priVtx_tracksPtSq = " <<priVtx_tracksPtSq <<endl ;
			    cout <<"recVtxs->size() = " <<recVtxs->size() <<endl ;
			  }
			}
		      }
		    }

		  PriVtxB0Less_n->push_back( B0_pvs.size() ) ;
		  PriVtxB0Less_X->push_back( B0LessPV.position().x() ) ;
		  PriVtxB0Less_Y->push_back( B0LessPV.position().y() ) ;
		  PriVtxB0Less_Z->push_back( B0LessPV.position().z() ) ;
		  PriVtxB0Less_EX->push_back( B0LessPV.xError() ) ;
		  PriVtxB0Less_EY->push_back( B0LessPV.yError() ) ;
		  PriVtxB0Less_EZ->push_back( B0LessPV.zError() ) ;
		  PriVtxB0Less_CL->push_back( ChiSquaredProbability( (double)(B0LessPV.chi2()), (double)(B0LessPV.ndof())) );
		  PriVtxB0Less_Chi2->push_back( B0LessPV.chi2() ) ;
		  PriVtxB0Less_tracks->push_back( B0LessPV.tracksSize() ) ;

		  /// dxy, dz, dxyE, dzE for kaons from BS
    		  math::XYZPoint BSVtx;
		  BSVtx = theBeamSpotVtx.position();
	          kaon1_dxy_BS->push_back( Track1->track()->dxy(BSVtx) );
                  kaon1_dz_BS->push_back( Track1->track()->dz(BSVtx) );
                  kaon2_dxy_BS->push_back( Track2->track()->dxy(BSVtx) );
                  kaon2_dz_BS->push_back( Track2->track()->dz(BSVtx) );


		  ////////////////// Lifetime wrt B0LessPV for B0 //////////////////
		  B0_v2e = B0LessPV.error();
		  B0_vXYe = B0_v1e.matrix() + B0_v2e.matrix();
		  /// 2D
		  B0_pvtx.SetXYZ(B0LessPV.position().x(), B0LessPV.position().y(), 0) ;
		  B0_vdiff = B0_vtx - B0_pvtx ;
		  B0_cosAlpha = B0_vdiff.Dot(B0_pperp)/(B0_vdiff.Perp()*B0_pperp.Perp());
		  B0_lxy = B0_vdiff.Perp();
		  B0_vDiff[0] = B0_vdiff.x(); B0_vDiff[1] = B0_vdiff.y(); B0_vDiff[2] = 0 ; // needed by Similarity method
		  B0_lxyErr = sqrt(ROOT::Math::Similarity(B0_vDiff,B0_vXYe)) / B0_vdiff.Perp();
		  B0_distXY = B0_vdistXY.distance( Vertex(*B0Cand_vertex_fromMCFit), Vertex(B0LessPV) ) ;
		  B0_ctau = B0_distXY.value() * B0_cosAlpha * B0Cand_fromMCFit->currentState().mass() / B0_pperp.Perp();
		  B0_ctauErr = sqrt(ROOT::Math::Similarity(B0_v3pperp,B0_vXYe)) * B0Cand_fromMCFit->currentState().mass() / (B0_pperp.Perp2());
		  /// 3D
		  B0_pvtx3D.SetXYZ(B0LessPV.position().x(), B0LessPV.position().y(), B0LessPV.position().z());
		  B0_vdiff3D = B0_vtx3D - B0_pvtx3D;
		  B0_cosAlpha3D = B0_vdiff3D.Dot(B0_pperp3D)/( B0_vdiff3D.Mag()*B0_pperp3D.Mag() );
		  B0_lxyz = B0_vdiff3D.Mag();
		  B0_vDiff3D[0] = B0_vdiff3D.x(); B0_vDiff3D[1] = B0_vdiff3D.y(); B0_vDiff3D[2] = B0_vdiff3D.z() ;
		  B0_lxyzErr = sqrt(ROOT::Math::Similarity(B0_vDiff3D,B0_vXYe)) / B0_vdiff3D.Mag();

		  b0CosAlphaB0LessPV->push_back( B0_cosAlpha ) ; b0CosAlpha3DB0LessPV->push_back( B0_cosAlpha3D ) ;
		  b0CTauB0LessPV->push_back( B0_ctau ) ; b0CTauB0LessPVE->push_back( B0_ctauErr ) ;
		  b0LxyB0LessPV->push_back( B0_lxy ) ; b0LxyB0LessPVE->push_back( B0_lxyErr ) ;
		  b0LxyzB0LessPV->push_back( B0_lxyz ) ; b0LxyzB0LessPVE->push_back( B0_lxyzErr ) ;

		  /// dxy, dz, dxyE, dzE for kaons from B0LessPV
		  math::XYZPoint B0LessPVvtx;
     		  B0LessPVvtx = B0LessPV.position();
                  kaon1_dxy_B0LessPV->push_back( Track1->track()->dxy(B0LessPVvtx) );
                  kaon1_dz_B0LessPV->push_back( Track1->track()->dz(B0LessPVvtx) );
                  kaon2_dxy_B0LessPV->push_back( Track2->track()->dxy(B0LessPVvtx) );
                  kaon2_dz_B0LessPV->push_back( Track2->track()->dz(B0LessPVvtx) );

                  kaon1_dxyE->push_back( Track1->track()->dxyError() );
                  kaon1_dzE->push_back( Track1->track()->dzError() );
		  kaon2_dxyE->push_back( Track2->track()->dxyError() );
		  kaon2_dzE->push_back( Track2->track()->dzError() );


		  /// Find the PV among the original offlinePV with the largest B0_cos(alpha)
		  Vertex theCosAlphaV = thePrimaryVtx ;
		  float maxCosAlpha = -1. ;

		  for (VertexCollection::const_iterator itv = recVtxs->begin(), itvend = recVtxs->end(); itv != itvend; ++itv) {
		    B0_pvtx.SetXYZ(itv->position().x(), itv->position().y(), 0) ;
		    B0_vdiff = B0_vtx - B0_pvtx ;
		    float cosAlpha_temp = B0_vdiff.Dot(B0_pperp) / (B0_vdiff.Perp()*B0_pperp.Perp()) ; // Perp() == Mag() when z = 0

		    if ( cosAlpha_temp > maxCosAlpha ) {
		      maxCosAlpha = cosAlpha_temp ;
		      theCosAlphaV = Vertex(*itv) ;
		    }
		  }

		  PriVtx_B0CosAlpha_n->push_back( recVtxs->size() ) ;
		  PriVtx_B0CosAlpha_X->push_back( theCosAlphaV.position().x() ) ;
		  PriVtx_B0CosAlpha_Y->push_back( theCosAlphaV.position().y() ) ;
		  PriVtx_B0CosAlpha_Z->push_back( theCosAlphaV.position().z() ) ;
		  PriVtx_B0CosAlpha_EX->push_back( theCosAlphaV.xError() ) ;
		  PriVtx_B0CosAlpha_EY->push_back( theCosAlphaV.yError() ) ;
		  PriVtx_B0CosAlpha_EZ->push_back( theCosAlphaV.zError() ) ;
		  PriVtx_B0CosAlpha_CL->push_back( ChiSquaredProbability((double)(theCosAlphaV.chi2()), (double)(theCosAlphaV.ndof())) ) ;
		  PriVtx_B0CosAlpha_Chi2->push_back( theCosAlphaV.chi2() ) ;
		  PriVtx_B0CosAlpha_tracks->push_back( theCosAlphaV.tracksSize() ) ;
  

                  /// Find the PV among the original offlinePV with the largest B0_cos(alpha) 3D
		  Vertex theCosAlpha3DV = thePrimaryVtx ;
		  float maxCosAlpha3D = -1. ;

		  for (VertexCollection::const_iterator itv = recVtxs->begin(), itvend = recVtxs->end(); itv != itvend; ++itv) {
		    B0_pvtx3D.SetXYZ(itv->position().x(), itv->position().y(), itv->position().z()) ;
		    B0_vdiff3D = B0_vtx3D - B0_pvtx3D ;
		    float cosAlpha_temp3D = B0_vdiff3D.Dot(B0_pperp3D) / (B0_vdiff3D.Mag()*B0_pperp3D.Mag()) ;

		    if ( cosAlpha_temp3D > maxCosAlpha3D ) {
		      maxCosAlpha3D = cosAlpha_temp3D ;
		      theCosAlpha3DV = Vertex(*itv) ;
		    }
		  }


		  ////////////////// Lifetime wrt PV with largest B0_cos(alpha) 3D candidate //////////////////
		  B0_v2e = theCosAlpha3DV.error(); 
		  B0_vXYe = B0_v1e.matrix() + B0_v2e.matrix();
		  /// 2D
		  B0_pvtx.SetXYZ(theCosAlpha3DV.position().x(), theCosAlpha3DV.position().y(), 0) ;
		  B0_vdiff = B0_vtx - B0_pvtx ;
		  B0_cosAlpha = B0_vdiff.Dot(B0_pperp)/(B0_vdiff.Perp()*B0_pperp.Perp()); ;
		  B0_lxy = B0_vdiff.Perp();
		  B0_vDiff[0] = B0_vdiff.x(); B0_vDiff[1] = B0_vdiff.y(); B0_vDiff[2] = 0 ; // needed by Similarity method
		  B0_lxyErr = sqrt(ROOT::Math::Similarity(B0_vDiff,B0_vXYe)) / B0_vdiff.Perp();
		  B0_distXY = B0_vdistXY.distance( Vertex(*B0Cand_vertex_fromMCFit), Vertex(theCosAlphaV) ) ;
		  B0_ctau = B0_distXY.value() * B0_cosAlpha * B0Cand_fromMCFit->currentState().mass() / B0_pperp.Perp();
		  B0_ctauErr = sqrt(ROOT::Math::Similarity(B0_v3pperp,B0_vXYe)) * B0Cand_fromMCFit->currentState().mass() / (B0_pperp.Perp2());
		  B0_lxy = B0_vdiff.Dot(B0_pperp) / B0_pperp.Mag() ;
		  /// 3D
		  B0_pvtx3D.SetXYZ(theCosAlpha3DV.position().x(), theCosAlpha3DV.position().y(), theCosAlpha3DV.position().z()) ;
		  B0_vdiff3D = B0_vtx3D - B0_pvtx3D ;
		  B0_cosAlpha3D =  maxCosAlpha3D ;
		  B0_lxyz = B0_vdiff3D.Mag();
		  B0_vDiff3D[0] = B0_vdiff3D.x(); B0_vDiff3D[1] = B0_vdiff3D.y(); B0_vDiff3D[2] = B0_vdiff3D.z() ;
		  B0_lxyzErr = sqrt(ROOT::Math::Similarity(B0_vDiff3D,B0_vXYe)) / B0_vdiff3D.Mag();

		  b0CosAlphaPVCosAlpha3D->push_back( B0_cosAlpha ) ; b0CosAlpha3DPVCosAlpha3D->push_back( B0_cosAlpha3D ) ;
		  b0CTauPVCosAlpha3D->push_back( B0_ctau ) ; b0CTauPVCosAlpha3DE->push_back( B0_ctauErr ) ;
		  b0LxyPVCosAlpha3D->push_back( B0_lxy ) ; b0LxyPVCosAlpha3DE->push_back( B0_lxyErr ) ;
		  b0LxyzPVCosAlpha3D->push_back( B0_lxyz ) ; b0LxyzPVCosAlpha3DE->push_back( B0_lxyzErr ) ;

		  
                  /// Find the PV among the B0lessPV with the largest B0_cos(alpha)
		  Vertex theB0LessCosAlphaV = thePrimaryVtx ;
		  maxCosAlpha = -1. ;

		  for (vector<TransientVertex>::iterator itv = B0_pvs.begin(), itvend = B0_pvs.end(); itv != itvend; ++itv) {
		    B0_pvtx.SetXYZ(itv->position().x(), itv->position().y(), 0) ;
		    B0_vdiff = B0_vtx - B0_pvtx ;
		    float cosAlpha_temp = B0_vdiff.Dot(B0_pperp) / (B0_vdiff.Perp()*B0_pperp.Perp()) ; // Perp() == Mag() when z = 0

 		    if ( cosAlpha_temp > maxCosAlpha ) {
		      maxCosAlpha = cosAlpha_temp ;
		      theB0LessCosAlphaV = Vertex(*itv) ;
		    }
		  }

		  PriVtxB0Less_B0CosAlpha_n->push_back( B0_pvs.size() ) ;
		  PriVtxB0Less_B0CosAlpha_X->push_back( theB0LessCosAlphaV.position().x() ) ;
		  PriVtxB0Less_B0CosAlpha_Y->push_back( theB0LessCosAlphaV.position().y() ) ;
		  PriVtxB0Less_B0CosAlpha_Z->push_back( theB0LessCosAlphaV.position().z() ) ;
		  PriVtxB0Less_B0CosAlpha_EX->push_back( theB0LessCosAlphaV.xError() ) ;
		  PriVtxB0Less_B0CosAlpha_EY->push_back( theB0LessCosAlphaV.yError() ) ;
		  PriVtxB0Less_B0CosAlpha_EZ->push_back( theB0LessCosAlphaV.zError() ) ;
		  PriVtxB0Less_B0CosAlpha_CL->push_back( ChiSquaredProbability((double)(theB0LessCosAlphaV.chi2()), (double)(theB0LessCosAlphaV.ndof())) ) ;
		  PriVtxB0Less_B0CosAlpha_Chi2->push_back( theB0LessCosAlphaV.chi2() ) ;
		  PriVtxB0Less_B0CosAlpha_tracks->push_back( theB0LessCosAlphaV.tracksSize() ) ;


		  ////////////////// Lifetime wrt PV with largest B0_cos(alpha) candidate //////////////////
		  B0_v2e = theCosAlphaV.error(); 
		  B0_vXYe = B0_v1e.matrix() + B0_v2e.matrix();
		  /// 2D
		  B0_pvtx.SetXYZ(theCosAlphaV.position().x(), theCosAlphaV.position().y(), 0) ;
		  B0_vdiff = B0_vtx - B0_pvtx ;
		  B0_cosAlpha =  maxCosAlpha ;
		  B0_lxy = B0_vdiff.Perp();
		  B0_vDiff[0] = B0_vdiff.x(); B0_vDiff[1] = B0_vdiff.y(); B0_vDiff[2] = 0 ; // needed by Similarity method
		  B0_lxyErr = sqrt(ROOT::Math::Similarity(B0_vDiff,B0_vXYe)) / B0_vdiff.Perp();
		  B0_distXY = B0_vdistXY.distance( Vertex(*B0Cand_vertex_fromMCFit), Vertex(theCosAlphaV) ) ;
		  B0_ctau = B0_distXY.value() * B0_cosAlpha * B0Cand_fromMCFit->currentState().mass() / B0_pperp.Perp();
		  B0_ctauErr = sqrt(ROOT::Math::Similarity(B0_v3pperp,B0_vXYe)) * B0Cand_fromMCFit->currentState().mass() / (B0_pperp.Perp2());
		  B0_lxy = B0_vdiff.Dot(B0_pperp) / B0_pperp.Mag() ;
		  /// 3D
		  B0_pvtx3D.SetXYZ(theCosAlphaV.position().x(), theCosAlphaV.position().y(), theCosAlphaV.position().z());
		  B0_vdiff3D = B0_vtx3D - B0_pvtx3D;
		  B0_cosAlpha3D = B0_vdiff3D.Dot(B0_pperp3D)/( B0_vdiff3D.Mag()*B0_pperp3D.Mag() );
		  B0_lxyz = B0_vdiff3D.Mag();
		  B0_vDiff3D[0] = B0_vdiff3D.x(); B0_vDiff3D[1] = B0_vdiff3D.y(); B0_vDiff3D[2] = B0_vdiff3D.z() ;
		  B0_lxyzErr = sqrt(ROOT::Math::Similarity(B0_vDiff3D,B0_vXYe)) / B0_vdiff3D.Mag();

		  b0CosAlphaPVCosAlpha->push_back( B0_cosAlpha ) ; b0CosAlpha3DPVCosAlpha->push_back( B0_cosAlpha3D ) ;
		  b0CTauPVCosAlpha->push_back( B0_ctau ) ; b0CTauPVCosAlphaE->push_back( B0_ctauErr ) ;
		  b0LxyPVCosAlpha->push_back( B0_lxy ) ; b0LxyPVCosAlphaE->push_back( B0_lxyErr ) ;
		  b0LxyzPVCosAlpha->push_back( B0_lxyz ) ; b0LxyzPVCosAlphaE->push_back( B0_lxyzErr ) ;

		  PriVtx_B0CosAlpha3D_n->push_back( recVtxs->size() ) ;
		  PriVtx_B0CosAlpha3D_X->push_back( theCosAlpha3DV.position().x() ) ;
		  PriVtx_B0CosAlpha3D_Y->push_back( theCosAlpha3DV.position().y() ) ;
		  PriVtx_B0CosAlpha3D_Z->push_back( theCosAlpha3DV.position().z() ) ;
		  PriVtx_B0CosAlpha3D_EX->push_back( theCosAlpha3DV.xError() ) ;
		  PriVtx_B0CosAlpha3D_EY->push_back( theCosAlpha3DV.yError() ) ;
		  PriVtx_B0CosAlpha3D_EZ->push_back( theCosAlpha3DV.zError() ) ;
		  PriVtx_B0CosAlpha3D_CL->push_back( ChiSquaredProbability((double)(theCosAlpha3DV.chi2()), (double)(theCosAlpha3DV.ndof())) ) ;
		  PriVtx_B0CosAlpha3D_Chi2->push_back( theCosAlpha3DV.chi2() ) ;
		  PriVtx_B0CosAlpha3D_tracks->push_back( theCosAlpha3DV.tracksSize() ) ;
  

                  ////////////////// Lifetime wrt B0LessPV with largest B0_cos(alpha) candidate
		  B0_v2e = theB0LessCosAlphaV.error();
		  B0_vXYe = B0_v1e.matrix() + B0_v2e.matrix();
		  /// 2D
		  B0_pvtx.SetXYZ(theB0LessCosAlphaV.position().x(), theB0LessCosAlphaV.position().y(), 0) ;
		  B0_vdiff = B0_vtx - B0_pvtx ;
		  B0_cosAlpha =  maxCosAlpha ;
		  B0_lxy = B0_vdiff.Perp();
		  B0_vDiff[0] = B0_vdiff.x(); B0_vDiff[1] = B0_vdiff.y(); B0_vDiff[2] = 0 ; // needed by Similarity method
		  B0_lxyErr = sqrt(ROOT::Math::Similarity(B0_vDiff,B0_vXYe)) / B0_vdiff.Perp() ;
		  B0_distXY = B0_vdistXY.distance( Vertex(*B0Cand_vertex_fromMCFit), Vertex(theB0LessCosAlphaV) ) ;
		  B0_ctau = B0_distXY.value() * B0_cosAlpha * B0Cand_fromMCFit->currentState().mass() / B0_pperp.Perp();
		  B0_ctauErr = sqrt(ROOT::Math::Similarity(B0_v3pperp,B0_vXYe)) * B0Cand_fromMCFit->currentState().mass() / (B0_pperp.Perp2());
		  B0_lxy = B0_vdiff.Dot(B0_pperp) / B0_pperp.Mag() ;
		  /// 3D
		  B0_pvtx3D.SetXYZ(theB0LessCosAlphaV.position().x(), theB0LessCosAlphaV.position().y(), theB0LessCosAlphaV.position().z());
		  B0_vdiff3D = B0_vtx3D - B0_pvtx3D;
		  B0_cosAlpha3D = B0_vdiff3D.Dot(B0_pperp3D)/( B0_vdiff3D.Mag()*B0_pperp3D.Mag() );
		  B0_lxyz = B0_vdiff3D.Mag();
		  B0_vDiff3D[0] = B0_vdiff3D.x(); B0_vDiff3D[1] = B0_vdiff3D.y(); B0_vDiff3D[2] = B0_vdiff3D.z() ;
		  B0_lxyzErr = sqrt(ROOT::Math::Similarity(B0_vDiff3D,B0_vXYe)) / B0_vdiff3D.Mag();

		  b0CosAlphaB0LessPVCosAlpha->push_back( B0_cosAlpha ) ; b0CosAlpha3DB0LessPVCosAlpha->push_back( B0_cosAlpha3D ) ;
		  b0CTauB0LessPVCosAlpha->push_back( B0_ctau ) ; b0CTauB0LessPVCosAlphaE->push_back( B0_ctauErr ) ;
		  b0LxyB0LessPVCosAlpha->push_back( B0_lxy ) ; b0LxyB0LessPVCosAlphaE->push_back( B0_lxyErr ) ;
		  b0LxyzB0LessPVCosAlpha->push_back( B0_lxyz ) ; b0LxyzB0LessPVCosAlphaE->push_back( B0_lxyzErr ) ;


		  /// Find the PV among the B0lessPV with the largest B0_cos(alpha) 3D
		  Vertex theB0LessCosAlpha3DV = thePrimaryVtx ;
		  maxCosAlpha3D = -1. ;

		  for (vector<TransientVertex>::iterator itv = B0_pvs.begin(), itvend = B0_pvs.end(); itv != itvend; ++itv) {
		    B0_pvtx3D.SetXYZ(itv->position().x(), itv->position().y(), itv->position().z()) ;
		    B0_vdiff3D = B0_vtx3D - B0_pvtx3D ;
		    float cosAlpha_temp3D = B0_vdiff3D.Dot(B0_pperp3D) / (B0_vdiff3D.Mag()*B0_pperp3D.Mag()) ;

		    if ( cosAlpha_temp3D > maxCosAlpha3D ) {
		      maxCosAlpha3D = cosAlpha_temp3D ;
		      theB0LessCosAlpha3DV = Vertex(*itv) ;
		    }
		  }

		  PriVtxB0Less_B0CosAlpha3D_n->push_back( B0_pvs.size() ) ;
		  PriVtxB0Less_B0CosAlpha3D_X->push_back( theB0LessCosAlpha3DV.position().x() ) ;
		  PriVtxB0Less_B0CosAlpha3D_Y->push_back( theB0LessCosAlpha3DV.position().y() ) ;
		  PriVtxB0Less_B0CosAlpha3D_Z->push_back( theB0LessCosAlpha3DV.position().z() ) ;
		  PriVtxB0Less_B0CosAlpha3D_EX->push_back( theB0LessCosAlpha3DV.xError() ) ;
		  PriVtxB0Less_B0CosAlpha3D_EY->push_back( theB0LessCosAlpha3DV.yError() ) ;
		  PriVtxB0Less_B0CosAlpha3D_EZ->push_back( theB0LessCosAlpha3DV.zError() ) ;
		  PriVtxB0Less_B0CosAlpha3D_CL->push_back( ChiSquaredProbability((double)(theB0LessCosAlpha3DV.chi2()), (double)(theB0LessCosAlpha3DV.ndof())) ) ;
		  PriVtxB0Less_B0CosAlpha3D_Chi2->push_back( theB0LessCosAlpha3DV.chi2() ) ;
		  PriVtxB0Less_B0CosAlpha3D_tracks->push_back( theB0LessCosAlpha3DV.tracksSize() ) ;


		  ////////////////// Lifetime wrt B0LessPV with largest B0_cos(alpha) 3D candidate
		  B0_v2e = theB0LessCosAlpha3DV.error();
		  B0_vXYe = B0_v1e.matrix() + B0_v2e.matrix();
		  /// 2D
		  B0_pvtx.SetXYZ(theB0LessCosAlpha3DV.position().x(), theB0LessCosAlpha3DV.position().y(), 0) ;
		  B0_vdiff = B0_vtx - B0_pvtx ;
		  B0_cosAlpha = B0_vdiff.Dot(B0_pperp)/(B0_vdiff.Perp()*B0_pperp.Perp());
		  B0_lxy = B0_vdiff.Perp();
		  B0_vDiff[0] = B0_vdiff.x(); B0_vDiff[1] = B0_vdiff.y(); B0_vDiff[2] = 0 ; // needed by Similarity method
		  B0_lxyErr = sqrt(ROOT::Math::Similarity(B0_vDiff,B0_vXYe)) / B0_vdiff.Perp();
		  B0_distXY = B0_vdistXY.distance( Vertex(*B0Cand_vertex_fromMCFit), Vertex(theCosAlphaV) ) ;
		  B0_ctau = B0_distXY.value() * B0_cosAlpha * B0Cand_fromMCFit->currentState().mass() / B0_pperp.Perp();
		  B0_ctauErr = sqrt(ROOT::Math::Similarity(B0_v3pperp,B0_vXYe)) * B0Cand_fromMCFit->currentState().mass() / (B0_pperp.Perp2());
		  B0_lxy = B0_vdiff.Dot(B0_pperp) / B0_pperp.Mag() ;
		  /// 3D
		  B0_pvtx3D.SetXYZ(theB0LessCosAlpha3DV.position().x(), theB0LessCosAlpha3DV.position().y(), theB0LessCosAlpha3DV.position().z()) ;
		  B0_vdiff3D = B0_vtx3D - B0_pvtx3D ;
		  B0_cosAlpha3D =  maxCosAlpha3D ;
		  B0_lxyz = B0_vdiff3D.Mag();
		  B0_vDiff3D[0] = B0_vdiff3D.x(); B0_vDiff3D[1] = B0_vdiff3D.y(); B0_vDiff3D[2] = B0_vdiff3D.z() ;
		  B0_lxyzErr = sqrt(ROOT::Math::Similarity(B0_vDiff3D,B0_vXYe)) / B0_vdiff3D.Mag();
		  B0_lxy = B0_vdiff3D.Dot(B0_pperp) / B0_pperp.Mag() ;

		  b0CosAlphaB0LessPVCosAlpha3D->push_back( B0_cosAlpha ) ; b0CosAlpha3DB0LessPVCosAlpha3D->push_back( B0_cosAlpha3D ) ;
		  b0CTauB0LessPVCosAlpha3D->push_back( B0_ctau ) ; b0CTauB0LessPVCosAlpha3DE->push_back( B0_ctauErr ) ;
		  b0LxyB0LessPVCosAlpha3D->push_back( B0_lxy ) ; b0LxyB0LessPVCosAlpha3DE->push_back( B0_lxyErr ) ;
		  b0LxyzB0LessPVCosAlpha3D->push_back( B0_lxyz ) ; b0LxyzB0LessPVCosAlpha3DE->push_back( B0_lxyzErr ) ;

		  Vertex theOtherV = thePrimaryVtx;
		  if (resolveAmbiguity_) {
		    float minDz = 999999. ;
		    if (!addB0lessPrimaryVertex_) {
		      for (VertexCollection::const_iterator itv = recVtxs->begin(), itvend = recVtxs->end(); itv != itvend; ++itv)
			{
			  float deltaZ = fabs((*B0Cand_vertex_fromMCFit).position().z() - itv->position().z()) ;
			  if ( deltaZ < minDz ) {
			    minDz = deltaZ;
			    thePrimaryVtx = Vertex(*itv);
			    theOtherV = thePrimaryVtx;
			  }
			}
		    } else {
		      for (vector<TransientVertex>::iterator itv2 = B0_pvs.begin(), itvend2 = B0_pvs.end(); itv2 != itvend2; ++itv2)
			{
			  float deltaZ = fabs((*B0Cand_vertex_fromMCFit).position().z() - itv2->position().z()) ;
			  if ( deltaZ < minDz ) {
			    minDz = deltaZ;
			    Vertex B0LessPV = Vertex(*itv2);
			    thePrimaryVtx = B0LessPV;
			    theOtherV = B0LessPV;
			  }
			}
		    }
		  }

		  Vertex TheOtherVertex3D = thePrimaryVtx;
		  if (Debug_) cout<<" choose PV ="<< endl;
		  Int_t theB0CorrPV_multiplicity = -1 ;
		  if (resolveAmbiguity_) {
		    float minDz = 999999.;
		    if (!addB0lessPrimaryVertex_) {
		      theB0CorrPV_multiplicity = recVtxs->size() ;
		      for (VertexCollection::const_iterator itv = recVtxs->begin(), itvend = recVtxs->end(); itv != itvend; ++itv) {
			float deltaZ = fabs((*B0Cand_vertex_fromMCFit).position().z() - itv->position().z()) ;
			if ( deltaZ < minDz ) {
			  minDz = deltaZ;
			  TheOtherVertex3D = Vertex(*itv);
			}
		      }
		    } else {
		      theB0CorrPV_multiplicity = B0_pvs.size() ;
		      for (vector<TransientVertex>::iterator itv2 = B0_pvs.begin(), itvend2 = B0_pvs.end(); itv2 != itvend2; ++itv2) {
			VertexDistance3D a3d;
			float deltaZ   = a3d.distance(Vertex(*itv2), Vertex(*B0Cand_vertex_fromMCFit)).value();
			if ( deltaZ < minDz ) {
			  minDz = deltaZ;
			  Vertex XLessPV = Vertex(*itv2);
			  TheOtherVertex3D = XLessPV;
			  //cout<<" z(X) - z(vtx) min="<<minDz<<endl;
			}

		      }
		    }
		  }

		  PriVtxB0Corr_n->push_back( theB0CorrPV_multiplicity ) ;
		  PriVtxB0Corr_X->push_back( thePrimaryVtx.position().x() ) ;
		  PriVtxB0Corr_Y->push_back( thePrimaryVtx.position().y() ) ;
		  PriVtxB0Corr_Z->push_back( thePrimaryVtx.position().z() ) ;
		  PriVtxB0Corr_EX->push_back( thePrimaryVtx.xError() ) ;
		  PriVtxB0Corr_EY->push_back( thePrimaryVtx.yError() ) ;
		  PriVtxB0Corr_EZ->push_back( thePrimaryVtx.zError() ) ;
		  PriVtxB0Corr_CL->push_back( ChiSquaredProbability( (double)(thePrimaryVtx.chi2()), (double)(thePrimaryVtx.ndof())) );
		  PriVtxB0Corr_Chi2->push_back( thePrimaryVtx.chi2() ) ;
		  PriVtxB0Corr_tracks->push_back( thePrimaryVtx.tracksSize() ) ;


		  ////////////////// Lifetime wrt PV with smaller longitudinal X impact parameter for B0  //////////////////
		  B0_pvtx.SetXYZ(theOtherV.position().x(), theOtherV.position().y(), 0); 
		  B0_vdiff = B0_vtx - B0_pvtx;
		  B0_cosAlpha = B0_vdiff.Dot(B0_pperp) / (B0_vdiff.Perp()*B0_pperp.Perp());
		  B0_distXY = B0_vdistXY.distance(Vertex(*B0Cand_vertex_fromMCFit), Vertex(theOtherV));
		  double B0_ctauPVX = B0_distXY.value() * B0_cosAlpha * B0Cand_fromMCFit->currentState().mass() / B0_pperp.Perp();
		  GlobalError B0_v1eX = (Vertex(*B0Cand_vertex_fromMCFit)).error();
		  GlobalError B0_v2eX = theOtherV.error();
		  AlgebraicSymMatrix33 B0_vXYeX = B0_v1eX.matrix() + B0_v2eX.matrix();
		  double ctauErrPVX = sqrt(ROOT::Math::Similarity(B0_v3pperp,B0_vXYeX)) * B0Cand_fromMCFit->currentState().mass() / (B0_pperp.Perp2());
		  float lxyPVX = B0_vdiff.Dot(B0_pperp) / B0_pperp.Mag() ;
		  float lxyzPVX = B0_vdiff3D.Dot(B0_pperp3D) / B0_pperp3D.Mag() ;
		  b0CosAlphaPVX->push_back(B0_cosAlpha);
		  b0CTauPVX->push_back(B0_ctauPVX); b0CTauPVXE->push_back(ctauErrPVX);
		  b0LxyPVX->push_back(lxyPVX);
		  b0LxyzPVX->push_back(lxyzPVX);
		  VertexDistance3D a3d;
		  float Dist3DPV     = a3d.distance(TheOtherVertex3D, Vertex(*B0Cand_vertex_fromMCFit)).value() ;
		  float Dist3DPV_err = a3d.distance(TheOtherVertex3D, Vertex(*B0Cand_vertex_fromMCFit)).error() ;
		  b0CTauPVX_3D->push_back(Dist3DPV);
		  b0CTauPVX_3D_err->push_back(Dist3DPV_err);
		  //cout << Dist3DPV << " " << Dist3DPV_err << endl;
		  B0_MuMuIdx->push_back(nMuMu-1);
		  B0_ka1Idx->push_back(std::distance(theKaonRefittedPATTrackHandle->begin(), Track1));
		  B0_ka2Idx->push_back(std::distance(theKaonRefittedPATTrackHandle->begin(), Track2));
		  nB0++;
		  b0Daughters.clear();



		  ////////////////// flag for checking the Kaons from PV or not PV //////////////////
		  /// flag for kaon1
		  vector<TransientTrack> vertexTracksKaon1;
                  //cout << "\nthePrimaryVtx.tracksSize() = " << thePrimaryVtx.tracksSize() << endl;
                  //cout << "thePrimaryVtx.nTracks() = " << thePrimaryVtx.nTracks() << endl;
		  for ( std::vector<TrackBaseRef >::const_iterator iTrack = thePrimaryVtx.tracks_begin(); iTrack != thePrimaryVtx.tracks_end(); ++iTrack) {

		    TrackRef trackRefkaon1 = iTrack->castTo<TrackRef>();
		    //cout << "\ntrackRefkaon1 = " << trackRefkaon1 << endl;
		    //cout <<"before match" ;
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
                  //cout << "\nthePrimaryVtx.tracksSize() = " << thePrimaryVtx.tracksSize() << endl;
                  //cout << "thePrimaryVtx.nTracks() = " << thePrimaryVtx.nTracks() << endl;
		  for ( std::vector<TrackBaseRef >::const_iterator iTrack = thePrimaryVtx.tracks_begin(); iTrack != thePrimaryVtx.tracks_end(); ++iTrack) {

		    TrackRef trackRefkaon2 = iTrack->castTo<TrackRef>();
                    //cout << "\ntrackRefkaon2 = " << trackRefkaon2 << endl;
                    //cout <<"before match" ;
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
  // AT THE END OF THE EVENT fill the tree and clear the vectors
  // ===========================================================

  if (nB0 > 0)
    X_One_Tree_->Fill() ;

  /// trigger stuff
  trigRes->clear(); trigNames->clear(); L1TT->clear(); MatchTriggerNames->clear();
  /// event numbers
  runNum = 0; evtNum = 0; lumiNum = 0;
  /// counters for B0 
  nMu = 0; nMuMu = 0; nB0 = 0; nKK = 0;
  nB0_pre0 = 0; nB0_pre1 = 0; nB0_pre2 = 0; nB0_pre3 = 0; nB0_pre4 = 0; nB0_pre5 = 0; nB0_pre6 = 0; nB0_pre7 = 0; nB0_pre8 = 0; nB0_pre9 = 0; nB0_pre10 = 0; nB0_pre11 = 0; nB0_pre12 = 0; nB0_pre13 = 0; nB0_pre14 = 0; nB0_pre15 = 0; 
  //nX = 0;
  /// indices
  mu1Idx->clear(); mu2Idx->clear();
  ka1Idx->clear(); ka2Idx->clear();
  B0_MuMuIdx->clear(); B0_ka1Idx->clear(); B0_ka2Idx->clear();

  /// MC Analysis
  if (doMC) {
    // Gen Primary Vertex
    n_genEvtVtx = 0;
    genEvtVtx_X->clear(); genEvtVtx_Y->clear(); genEvtVtx_Z->clear();
    genEvtVtx_particles->clear();
    n_B0Ancestors->clear();
    nMCAll = 0, nMCB0 = 0; //nMCB0Vtx = 0;
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
  DRMuMuK1->clear(); DRMuMuK2->clear(); DRb0K1->clear(); DRb0K2->clear(); 
  /// Primary Vertex with "MuMu correction"
  PriVtxMuMuCorr_n->clear();
  PriVtxMuMuCorr_X->clear(); PriVtxMuMuCorr_Y->clear(); PriVtxMuMuCorr_Z->clear();
  PriVtxMuMuCorr_EX->clear(); PriVtxMuMuCorr_EY->clear(); PriVtxMuMuCorr_EZ->clear();
  PriVtxMuMuCorr_Chi2->clear(); PriVtxMuMuCorr_CL->clear(); PriVtxMuMuCorr_tracks->clear();
  nTrk->clear();
  /// B0 cand 
  b0Mass->clear(); b0Vtx_CL->clear(); b0Vtx_Chi2->clear();
  b0Px->clear(); b0Py->clear(); b0Pz->clear();
  b0PxE->clear(); b0PyE->clear(); b0PzE->clear();
  b0DecayVtx_X->clear(); b0DecayVtx_Y->clear(); b0DecayVtx_Z->clear();
  b0DecayVtx_XE->clear(); b0DecayVtx_YE->clear(); b0DecayVtx_ZE->clear();
  /// muons and tracks after B0 cand fit
  mu1Px_MuMuKK->clear(); mu1Py_MuMuKK->clear(); mu1Pz_MuMuKK->clear(); mu1E_MuMuKK->clear();
  mu2Px_MuMuKK->clear(); mu2Py_MuMuKK->clear(); mu2Pz_MuMuKK->clear(); mu2E_MuMuKK->clear();
  k1Px_MuMuKK->clear(); k1Py_MuMuKK->clear(); k1Pz_MuMuKK->clear(); k1E_MuMuKK->clear();
  kaon1_nsigdedx->clear(); kaon1_dedx->clear(); kaon1_dedxMass->clear(); kaon1_theo->clear(); kaon1_sigma->clear();
  kaon1_dedx_byHits->clear(); kaon1_dedxErr_byHits->clear(); kaon1_saturMeas_byHits->clear(); kaon1_Meas_byHits->clear();
  k2Px_MuMuKK->clear(); k2Py_MuMuKK->clear(); k2Pz_MuMuKK->clear(); k2E_MuMuKK->clear();
  kaon2_nsigdedx->clear(); kaon2_dedx->clear(); kaon2_dedxMass->clear(); kaon2_theo->clear(); kaon2_sigma->clear();
  kaon2_dedx_byHits->clear(); kaon2_dedxErr_byHits->clear(); kaon2_saturMeas_byHits->clear(); kaon2_Meas_byHits->clear();
  /// Primary Vertex with "B0 correction"
  PriVtxB0Less_n->clear();
  PriVtxB0Less_X->clear(); PriVtxB0Less_Y->clear(); PriVtxB0Less_Z->clear();
  PriVtxB0Less_EX->clear(); PriVtxB0Less_EY->clear(); PriVtxB0Less_EZ->clear();
  PriVtxB0Less_Chi2->clear(); PriVtxB0Less_CL->clear(); PriVtxB0Less_tracks->clear();
  /// Primary Vertex with largest B0_cos(alpha)
  B0LessPV_tracksPtSq->clear(); B0LessPV_4tracksPtSq->clear();
  PriVtx_B0CosAlpha_n->clear();
  PriVtx_B0CosAlpha_X->clear(); PriVtx_B0CosAlpha_Y->clear(); PriVtx_B0CosAlpha_Z->clear();
  PriVtx_B0CosAlpha_EX->clear(); PriVtx_B0CosAlpha_EY->clear(); PriVtx_B0CosAlpha_EZ->clear();
  PriVtx_B0CosAlpha_Chi2->clear(); PriVtx_B0CosAlpha_CL->clear(); PriVtx_B0CosAlpha_tracks->clear();
  PriVtxB0Less_B0CosAlpha_n->clear();
  PriVtxB0Less_B0CosAlpha_X->clear(); PriVtxB0Less_B0CosAlpha_Y->clear(); PriVtxB0Less_B0CosAlpha_Z->clear();
  PriVtxB0Less_B0CosAlpha_EX->clear(); PriVtxB0Less_B0CosAlpha_EY->clear(); PriVtxB0Less_B0CosAlpha_EZ->clear();
  PriVtxB0Less_B0CosAlpha_Chi2->clear(); PriVtxB0Less_B0CosAlpha_CL->clear(); PriVtxB0Less_B0CosAlpha_tracks->clear();
  /// Primary Vertex with "B0 correction" 
  PriVtxB0Corr_n->clear();
  PriVtxB0Corr_X->clear(); PriVtxB0Corr_Y->clear(); PriVtxB0Corr_Z->clear();
  PriVtxB0Corr_EX->clear(); PriVtxB0Corr_EY->clear(); PriVtxB0Corr_EZ->clear();
  PriVtxB0Corr_Chi2->clear(); PriVtxB0Corr_CL->clear(); PriVtxB0Corr_tracks->clear();
  /// Lifetime variables for B0 
  b0CosAlphaBS->clear(); b0CosAlpha3DBS->clear(); b0CTauBS->clear(); b0CTauBSE->clear(); b0LxyBS->clear(); b0LxyBSE->clear(); b0LxyzBS->clear(); b0LxyzBSE->clear();
  b0CosAlphaPV->clear(); b0CosAlpha3DPV->clear(); b0CTauPV->clear(); b0CTauPVE->clear(); b0LxyPV->clear(); b0LxyPVE->clear(); b0LxyzPV->clear(); b0LxyzPVE->clear();
  b0CosAlphaPVCosAlpha->clear(); b0CosAlpha3DPVCosAlpha->clear(); b0CTauPVCosAlpha->clear(); b0CTauPVCosAlphaE->clear(); b0LxyPVCosAlpha->clear(); b0LxyPVCosAlphaE->clear(); b0LxyzPVCosAlpha->clear(); b0LxyzPVCosAlphaE->clear();
  b0CosAlphaPVCosAlpha3D->clear(); b0CosAlpha3DPVCosAlpha3D->clear(); b0CTauPVCosAlpha3D->clear(); b0CTauPVCosAlpha3DE->clear(); b0LxyPVCosAlpha3D->clear(); b0LxyPVCosAlpha3DE->clear(); b0LxyzPVCosAlpha3D->clear(); b0LxyzPVCosAlpha3DE->clear();
  b0CosAlphaB0LessPV->clear(); b0CosAlpha3DB0LessPV->clear(); b0CTauB0LessPV->clear() ; b0CTauB0LessPVE->clear() ; b0LxyB0LessPV->clear() ; b0LxyB0LessPVE->clear() ; b0LxyzB0LessPV->clear() ; b0LxyzB0LessPVE->clear() ;
  b0CosAlphaB0LessPVCosAlpha->clear(); b0CosAlpha3DB0LessPVCosAlpha->clear(); b0CTauB0LessPVCosAlpha->clear() ; b0CTauB0LessPVCosAlphaE->clear() ; b0LxyB0LessPVCosAlpha->clear() ; b0LxyB0LessPVCosAlphaE->clear() ; b0LxyzB0LessPVCosAlpha->clear() ; b0LxyzB0LessPVCosAlphaE->clear() ;
  b0CosAlphaB0LessPVCosAlpha3D->clear(); b0CosAlpha3DB0LessPVCosAlpha3D->clear(); b0CTauB0LessPVCosAlpha3D->clear() ; b0CTauB0LessPVCosAlpha3DE->clear() ; b0LxyB0LessPVCosAlpha3D->clear() ; b0LxyB0LessPVCosAlpha3DE->clear() ; b0LxyzB0LessPVCosAlpha3D->clear() ; b0LxyzB0LessPVCosAlpha3DE->clear() ;
  b0CosAlphaPVX->clear(); b0CTauPVX->clear(); b0CTauPVXE->clear(); b0LxyPVX->clear(); b0LxyzPVX->clear();
  b0CTauPVX_3D->clear(); b0CTauPVX_3D_err->clear();
  /// dxy, dz, dxyE, dzE for kaons from PV, BS, B0LessPV
  kaon1_dxy_PV->clear(); kaon1_dz_PV->clear(); kaon2_dxy_PV->clear(); kaon2_dz_PV->clear();
  kaon1_dxy_BS->clear(); kaon1_dz_BS->clear(); kaon2_dxy_BS->clear(); kaon2_dz_BS->clear();
  kaon1_dxy_B0LessPV->clear(); kaon1_dz_B0LessPV->clear(); kaon2_dxy_B0LessPV->clear(); kaon2_dz_B0LessPV->clear();
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
    X_One_Tree_->Branch("n_B0Ancestors", &n_B0Ancestors);
    X_One_Tree_->Branch("nMCAll", &nMCAll, "nMCAll/i");
    X_One_Tree_->Branch("MCPdgIdAll", &MCPdgIdAll);
    X_One_Tree_->Branch("MCDanNumAll", &MCDanNumAll);
    X_One_Tree_->Branch("nMCB0",&nMCB0,"nMCB0/i");
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
  X_One_Tree_->Branch("DRb0K1",&DRb0K1);
  X_One_Tree_->Branch("DRb0K2",&DRb0K2);
  /// counters for B0
  X_One_Tree_->Branch("nB0",&nB0,"nB0/i");
  X_One_Tree_->Branch("nB0_pre0",&nB0_pre0,"nB0_pre0/i");
  X_One_Tree_->Branch("nB0_pre1",&nB0_pre1,"nB0_pre1/i");
  X_One_Tree_->Branch("nB0_pre2",&nB0_pre2,"nB0_pre2/i");
  X_One_Tree_->Branch("nB0_pre3",&nB0_pre3,"nB0_pre3/i");
  X_One_Tree_->Branch("nB0_pre4",&nB0_pre4,"nB0_pre4/i");
  X_One_Tree_->Branch("nB0_pre5",&nB0_pre5,"nB0_pre5/i");
  X_One_Tree_->Branch("nB0_pre6",&nB0_pre6,"nB0_pre6/i");
  X_One_Tree_->Branch("nB0_pre7",&nB0_pre7,"nB0_pre7/i");
  X_One_Tree_->Branch("nB0_pre8",&nB0_pre8,"nB0_pre8/i");
  X_One_Tree_->Branch("nB0_pre9",&nB0_pre9,"nB0_pre9/i");
  X_One_Tree_->Branch("nB0_pre10",&nB0_pre10,"nB0_pre10/i");
  X_One_Tree_->Branch("nB0_pre11",&nB0_pre11,"nB0_pre11/i");
  X_One_Tree_->Branch("nB0_pre12",&nB0_pre12,"nB0_pre12/i");
  X_One_Tree_->Branch("nB0_pre13",&nB0_pre13,"nB0_pre13/i");
  X_One_Tree_->Branch("nB0_pre14",&nB0_pre14,"nB0_pre14/i");
  X_One_Tree_->Branch("nB0_pre15",&nB0_pre15,"nB0_pre15/i");
  //X_One_Tree_->Branch("nX",&nX,"nX/i");
  /// B0 cand 
  X_One_Tree_->Branch("B0Mass",&b0Mass);
  X_One_Tree_->Branch("B0Px",&b0Px);
  X_One_Tree_->Branch("B0Py",&b0Py);
  X_One_Tree_->Branch("B0Pz",&b0Pz);
  X_One_Tree_->Branch("B0PxE",&b0PxE);
  X_One_Tree_->Branch("B0PyE",&b0PyE);
  X_One_Tree_->Branch("B0PzE",&b0PzE);
  X_One_Tree_->Branch("B0Vtx_CL",&b0Vtx_CL);
  X_One_Tree_->Branch("B0Vtx_Chi2",&b0Vtx_Chi2);
  X_One_Tree_->Branch("B0DecayVtx_X",&b0DecayVtx_X);
  X_One_Tree_->Branch("B0DecayVtx_Y",&b0DecayVtx_Y);
  X_One_Tree_->Branch("B0DecayVtx_Z",&b0DecayVtx_Z);
  X_One_Tree_->Branch("B0DecayVtx_XE",&b0DecayVtx_XE);
  X_One_Tree_->Branch("B0DecayVtx_YE",&b0DecayVtx_YE);
  X_One_Tree_->Branch("B0DecayVtx_ZE",&b0DecayVtx_ZE);
  X_One_Tree_->Branch("B0CosAlphaBS", &b0CosAlphaBS);
  X_One_Tree_->Branch("B0CosAlpha3DBS", &b0CosAlpha3DBS);
  X_One_Tree_->Branch("B0CTauBS", &b0CTauBS);
  X_One_Tree_->Branch("B0CTauBSE", &b0CTauBSE);
  X_One_Tree_->Branch("B0LxyBS", &b0LxyBS);
  X_One_Tree_->Branch("B0LxyBSE", &b0LxyBSE);
  X_One_Tree_->Branch("B0LxyzBS", &b0LxyzBS);
  X_One_Tree_->Branch("B0LxyzBSE", &b0LxyzBSE);
  X_One_Tree_->Branch("B0CosAlphaPV", &b0CosAlphaPV);
  X_One_Tree_->Branch("B0CosAlpha3DPV", &b0CosAlpha3DPV);
  X_One_Tree_->Branch("B0CTauPV", &b0CTauPV);
  X_One_Tree_->Branch("B0CTauPVE", &b0CTauPVE);
  X_One_Tree_->Branch("B0LxyPV", &b0LxyPV);
  X_One_Tree_->Branch("B0LxyPVE", &b0LxyPVE);
  X_One_Tree_->Branch("B0LxyzPV", &b0LxyzPV);
  X_One_Tree_->Branch("B0LxyzPVE", &b0LxyzPVE);
  /// Primary Vertex with largest B0_cos(alpha) 
  X_One_Tree_->Branch("PriVtx_B0CosAlpha_n",&PriVtx_B0CosAlpha_n);
  X_One_Tree_->Branch("PriVtx_B0CosAlpha_X",&PriVtx_B0CosAlpha_X);
  X_One_Tree_->Branch("PriVtx_B0CosAlpha_Y",&PriVtx_B0CosAlpha_Y);
  X_One_Tree_->Branch("PriVtx_B0CosAlpha_Z",&PriVtx_B0CosAlpha_Z);
  X_One_Tree_->Branch("PriVtx_B0CosAlpha_EX",&PriVtx_B0CosAlpha_EX);
  X_One_Tree_->Branch("PriVtx_B0CosAlpha_EY",&PriVtx_B0CosAlpha_EY);
  X_One_Tree_->Branch("PriVtx_B0CosAlpha_EZ",&PriVtx_B0CosAlpha_EZ);
  X_One_Tree_->Branch("PriVtx_B0CosAlpha_Chi2",&PriVtx_B0CosAlpha_Chi2);
  X_One_Tree_->Branch("PriVtx_B0CosAlpha_CL",&PriVtx_B0CosAlpha_CL);
  X_One_Tree_->Branch("PriVtx_B0CosAlpha_tracks",&PriVtx_B0CosAlpha_tracks);
  X_One_Tree_->Branch("B0CosAlphaPVCosAlpha", &b0CosAlphaPVCosAlpha);
  X_One_Tree_->Branch("B0CosAlpha3DPVCosAlpha", &b0CosAlpha3DPVCosAlpha);
  X_One_Tree_->Branch("B0CTauPVCosAlpha", &b0CTauPVCosAlpha);
  X_One_Tree_->Branch("B0CTauPVCosAlphaE", &b0CTauPVCosAlphaE);
  X_One_Tree_->Branch("B0LxyPVCosAlpha", &b0LxyPVCosAlpha);
  X_One_Tree_->Branch("B0LxyPVCosAlphaE", &b0LxyPVCosAlphaE);
  X_One_Tree_->Branch("B0LxyzPVCosAlpha", &b0LxyzPVCosAlpha);
  X_One_Tree_->Branch("B0LxyzPVCosAlphaE", &b0LxyzPVCosAlphaE);
  X_One_Tree_->Branch("PriVtx_B0CosAlpha3D_n",&PriVtx_B0CosAlpha3D_n);
  X_One_Tree_->Branch("PriVtx_B0CosAlpha3D_X",&PriVtx_B0CosAlpha3D_X);
  X_One_Tree_->Branch("PriVtx_B0CosAlpha3D_Y",&PriVtx_B0CosAlpha3D_Y);
  X_One_Tree_->Branch("PriVtx_B0CosAlpha3D_Z",&PriVtx_B0CosAlpha3D_Z);
  X_One_Tree_->Branch("PriVtx_B0CosAlpha3D_EX",&PriVtx_B0CosAlpha3D_EX);
  X_One_Tree_->Branch("PriVtx_B0CosAlpha3D_EY",&PriVtx_B0CosAlpha3D_EY);
  X_One_Tree_->Branch("PriVtx_B0CosAlpha3D_EZ",&PriVtx_B0CosAlpha3D_EZ);
  X_One_Tree_->Branch("PriVtx_B0CosAlpha3D_Chi2",&PriVtx_B0CosAlpha3D_Chi2);
  X_One_Tree_->Branch("PriVtx_B0CosAlpha3D_CL",&PriVtx_B0CosAlpha3D_CL);
  X_One_Tree_->Branch("PriVtx_B0CosAlpha3D_tracks",&PriVtx_B0CosAlpha3D_tracks);
  X_One_Tree_->Branch("B0CosAlphaPVCosAlpha3D", &b0CosAlphaPVCosAlpha3D);
  X_One_Tree_->Branch("B0CosAlpha3DPVCosAlpha3D", &b0CosAlpha3DPVCosAlpha3D);
  X_One_Tree_->Branch("B0CTauPVCosAlpha3D", &b0CTauPVCosAlpha3D);
  X_One_Tree_->Branch("B0CTauPVCosAlpha3DE", &b0CTauPVCosAlpha3DE);
  X_One_Tree_->Branch("B0LxyPVCosAlpha3D", &b0LxyPVCosAlpha3D);
  X_One_Tree_->Branch("B0LxyPVCosAlpha3DE", &b0LxyPVCosAlpha3DE);
  X_One_Tree_->Branch("B0LxyzPVCosAlpha3D", &b0LxyzPVCosAlpha3D);
  X_One_Tree_->Branch("B0LxyzPVCosAlpha3DE", &b0LxyzPVCosAlpha3DE);

  X_One_Tree_->Branch("B0LessPV_tracksPtSq",&B0LessPV_tracksPtSq);
  X_One_Tree_->Branch("B0LessPV_4tracksPtSq",&B0LessPV_4tracksPtSq);
  X_One_Tree_->Branch("PriVtxB0Less_n",&PriVtxB0Less_n);
  X_One_Tree_->Branch("PriVtxB0Less_X",&PriVtxB0Less_X);
  X_One_Tree_->Branch("PriVtxB0Less_Y",&PriVtxB0Less_Y);
  X_One_Tree_->Branch("PriVtxB0Less_Z",&PriVtxB0Less_Z);
  X_One_Tree_->Branch("PriVtxB0Less_EX",&PriVtxB0Less_EX);
  X_One_Tree_->Branch("PriVtxB0Less_EY",&PriVtxB0Less_EY);
  X_One_Tree_->Branch("PriVtxB0Less_EZ",&PriVtxB0Less_EZ);
  X_One_Tree_->Branch("PriVtxB0Less_Chi2",&PriVtxB0Less_Chi2);
  X_One_Tree_->Branch("PriVtxB0Less_CL",&PriVtxB0Less_CL);
  X_One_Tree_->Branch("PriVtxB0Less_tracks",&PriVtxB0Less_tracks);
  X_One_Tree_->Branch("B0CosAlphaB0LessPV", &b0CosAlphaB0LessPV);
  X_One_Tree_->Branch("B0CosAlpha3DB0LessPV", &b0CosAlpha3DB0LessPV);
  X_One_Tree_->Branch("B0CTauB0LessPV", &b0CTauB0LessPV);
  X_One_Tree_->Branch("B0CTauB0LessPVE", &b0CTauB0LessPVE);
  X_One_Tree_->Branch("B0LxyB0LessPV", &b0LxyB0LessPV);
  X_One_Tree_->Branch("B0LxyB0LessPVE", &b0LxyB0LessPVE);
  X_One_Tree_->Branch("B0LxyzB0LessPV", &b0LxyzB0LessPV);
  X_One_Tree_->Branch("B0LxyzB0LessPVE", &b0LxyzB0LessPVE);
  X_One_Tree_->Branch("PriVtxB0Less_B0CosAlpha_n",&PriVtxB0Less_B0CosAlpha_n);
  X_One_Tree_->Branch("PriVtxB0Less_B0CosAlpha_X",&PriVtxB0Less_B0CosAlpha_X);
  X_One_Tree_->Branch("PriVtxB0Less_B0CosAlpha_Y",&PriVtxB0Less_B0CosAlpha_Y);
  X_One_Tree_->Branch("PriVtxB0Less_B0CosAlpha_Z",&PriVtxB0Less_B0CosAlpha_Z);
  X_One_Tree_->Branch("PriVtxB0Less_B0CosAlpha_EX",&PriVtxB0Less_B0CosAlpha_EX);
  X_One_Tree_->Branch("PriVtxB0Less_B0CosAlpha_EY",&PriVtxB0Less_B0CosAlpha_EY);
  X_One_Tree_->Branch("PriVtxB0Less_B0CosAlpha_EZ",&PriVtxB0Less_B0CosAlpha_EZ);
  X_One_Tree_->Branch("PriVtxB0Less_B0CosAlpha_Chi2",&PriVtxB0Less_B0CosAlpha_Chi2);
  X_One_Tree_->Branch("PriVtxB0Less_B0CosAlpha_CL",&PriVtxB0Less_B0CosAlpha_CL);
  X_One_Tree_->Branch("PriVtxB0Less_B0CosAlpha_tracks",&PriVtxB0Less_B0CosAlpha_tracks);
  X_One_Tree_->Branch("B0CosAlphaB0LessPVCosAlpha", &b0CosAlphaB0LessPVCosAlpha);
  X_One_Tree_->Branch("B0CosAlpha3DB0LessPVCosAlpha", &b0CosAlpha3DB0LessPVCosAlpha);
  X_One_Tree_->Branch("B0CTauB0LessPVCosAlpha", &b0CTauB0LessPVCosAlpha);
  X_One_Tree_->Branch("B0CTauB0LessPVCosAlphaE", &b0CTauB0LessPVCosAlphaE);
  X_One_Tree_->Branch("B0LxyB0LessPVCosAlpha", &b0LxyB0LessPVCosAlpha);
  X_One_Tree_->Branch("B0LxyB0LessPVCosAlphaE", &b0LxyB0LessPVCosAlphaE);
  X_One_Tree_->Branch("B0LxyzB0LessPVCosAlpha", &b0LxyzB0LessPVCosAlpha);
  X_One_Tree_->Branch("B0LxyzB0LessPVCosAlphaE", &b0LxyzB0LessPVCosAlphaE);
  X_One_Tree_->Branch("PriVtxB0Less_B0CosAlpha3D_n",&PriVtxB0Less_B0CosAlpha3D_n);
  X_One_Tree_->Branch("PriVtxB0Less_B0CosAlpha3D_X",&PriVtxB0Less_B0CosAlpha3D_X);
  X_One_Tree_->Branch("PriVtxB0Less_B0CosAlpha3D_Y",&PriVtxB0Less_B0CosAlpha3D_Y);
  X_One_Tree_->Branch("PriVtxB0Less_B0CosAlpha3D_Z",&PriVtxB0Less_B0CosAlpha3D_Z);
  X_One_Tree_->Branch("PriVtxB0Less_B0CosAlpha3D_EX",&PriVtxB0Less_B0CosAlpha3D_EX);
  X_One_Tree_->Branch("PriVtxB0Less_B0CosAlpha3D_EY",&PriVtxB0Less_B0CosAlpha3D_EY);
  X_One_Tree_->Branch("PriVtxB0Less_B0CosAlpha3D_EZ",&PriVtxB0Less_B0CosAlpha3D_EZ);
  X_One_Tree_->Branch("PriVtxB0Less_B0CosAlpha3D_Chi2",&PriVtxB0Less_B0CosAlpha3D_Chi2);
  X_One_Tree_->Branch("PriVtxB0Less_B0CosAlpha3D_CL",&PriVtxB0Less_B0CosAlpha3D_CL);
  X_One_Tree_->Branch("PriVtxB0Less_B0CosAlpha3D_tracks",&PriVtxB0Less_B0CosAlpha3D_tracks);
  X_One_Tree_->Branch("B0CosAlphaB0LessPVCosAlpha3D", &b0CosAlphaB0LessPVCosAlpha3D);
  X_One_Tree_->Branch("B0CosAlpha3DB0LessPVCosAlpha3D", &b0CosAlpha3DB0LessPVCosAlpha3D);
  X_One_Tree_->Branch("B0CTauB0LessPVCosAlpha3D", &b0CTauB0LessPVCosAlpha3D);
  X_One_Tree_->Branch("B0CTauB0LessPVCosAlpha3DE", &b0CTauB0LessPVCosAlpha3DE);
  X_One_Tree_->Branch("B0LxyB0LessPVCosAlpha3D", &b0LxyB0LessPVCosAlpha3D);
  X_One_Tree_->Branch("B0LxyB0LessPVCosAlpha3DE", &b0LxyB0LessPVCosAlpha3DE);
  X_One_Tree_->Branch("B0LxyzB0LessPVCosAlpha3D", &b0LxyzB0LessPVCosAlpha3D);
  X_One_Tree_->Branch("B0LxyzB0LessPVCosAlpha3DE", &b0LxyzB0LessPVCosAlpha3DE);
  /// Primary Vertex with "B0 correction" 
  X_One_Tree_->Branch("PriVtxB0Corr_n",&PriVtxB0Corr_n);
  X_One_Tree_->Branch("PriVtxB0Corr_X",&PriVtxB0Corr_X);
  X_One_Tree_->Branch("PriVtxB0Corr_Y",&PriVtxB0Corr_Y);
  X_One_Tree_->Branch("PriVtxB0Corr_Z",&PriVtxB0Corr_Z);
  X_One_Tree_->Branch("PriVtxB0Corr_EX",&PriVtxB0Corr_EX);
  X_One_Tree_->Branch("PriVtxB0Corr_EY",&PriVtxB0Corr_EY);
  X_One_Tree_->Branch("PriVtxB0Corr_EZ",&PriVtxB0Corr_EZ);
  X_One_Tree_->Branch("PriVtxB0Corr_Chi2",&PriVtxB0Corr_Chi2);
  X_One_Tree_->Branch("PriVtxB0Corr_CL",&PriVtxB0Corr_CL);
  X_One_Tree_->Branch("PriVtxB0Corr_tracks",&PriVtxB0Corr_tracks);
  /// Lifetime variables for B0 
  X_One_Tree_->Branch("B0CosAlphaPVX", &b0CosAlphaPVX);
  X_One_Tree_->Branch("B0CTauPVX", &b0CTauPVX);
  X_One_Tree_->Branch("B0CTauPVXE", &b0CTauPVXE);
  X_One_Tree_->Branch("B0LxyPVX", &b0LxyPVX);
  X_One_Tree_->Branch("B0LxyzPVX", &b0LxyzPVX);
  X_One_Tree_->Branch("B0CTauPVX_3D", &b0CTauPVX_3D);
  X_One_Tree_->Branch("B0CTauPVX_3D_err", &b0CTauPVX_3D_err);
  /// dxy, dz, dxyE, dzE for kaons from PV, BS, B0LessPV
  X_One_Tree_->Branch("kaon1_dxy_PV", &kaon1_dxy_PV);
  X_One_Tree_->Branch("kaon1_dz_PV", &kaon1_dz_PV);
  X_One_Tree_->Branch("kaon2_dxy_PV", &kaon2_dxy_PV);
  X_One_Tree_->Branch("kaon2_dz_PV", &kaon2_dz_PV);
  X_One_Tree_->Branch("kaon1_dxy_BS", &kaon1_dxy_BS);
  X_One_Tree_->Branch("kaon1_dz_BS", &kaon1_dz_BS);
  X_One_Tree_->Branch("kaon2_dxy_BS", &kaon2_dxy_BS);
  X_One_Tree_->Branch("kaon2_dz_BS", &kaon2_dz_BS);
  X_One_Tree_->Branch("kaon1_dxy_B0LessPV", &kaon1_dxy_B0LessPV);
  X_One_Tree_->Branch("kaon1_dz_B0LessPV", &kaon1_dz_B0LessPV);
  X_One_Tree_->Branch("kaon2_dxy_B0LessPV", &kaon2_dxy_B0LessPV);
  X_One_Tree_->Branch("kaon2_dz_B0LessPV", &kaon2_dz_B0LessPV);
  X_One_Tree_->Branch("kaon1_dxyE", &kaon1_dxyE);
  X_One_Tree_->Branch("kaon1_dzE", &kaon1_dzE);
  X_One_Tree_->Branch("kaon2_dxyE", &kaon2_dxyE);
  X_One_Tree_->Branch("kaon2_dzE", &kaon2_dzE);

  X_One_Tree_->Branch("B0MuMuIdx", &B0_MuMuIdx);
  X_One_Tree_->Branch("B0Kaon1Idx", &B0_ka1Idx);
  X_One_Tree_->Branch("B0Kaon2Idx", &B0_ka2Idx);

  X_One_Tree_->Branch("KKMass_err",&KKMass_err);
  X_One_Tree_->Branch("Kaon1FromPV",&Kaon1FromPV);
  X_One_Tree_->Branch("Kaon2FromPV",&Kaon2FromPV );

  /// Muons and tracks after B0 cand fit 
  X_One_Tree_->Branch("Muon1Px_MuMuKK", &mu1Px_MuMuKK);
  X_One_Tree_->Branch("Muon1Py_MuMuKK", &mu1Py_MuMuKK);
  X_One_Tree_->Branch("Muon1Pz_MuMuKK", &mu1Pz_MuMuKK);
  X_One_Tree_->Branch("Muon1E_MuMuKK", &mu1E_MuMuKK);
  X_One_Tree_->Branch("Muon2Px_MuMuKK", &mu2Px_MuMuKK);
  X_One_Tree_->Branch("Muon2Py_MuMuKK", &mu2Py_MuMuKK);
  X_One_Tree_->Branch("Muon2Pz_MuMuKK", &mu2Pz_MuMuKK);
  X_One_Tree_->Branch("Muon2E_MuMuKK", &mu2E_MuMuKK);
  X_One_Tree_->Branch("Kaon1Px_MuMuKK", &k1Px_MuMuKK);
  X_One_Tree_->Branch("Kaon1Py_MuMuKK", &k1Py_MuMuKK);
  X_One_Tree_->Branch("Kaon1Pz_MuMuKK", &k1Pz_MuMuKK);
  X_One_Tree_->Branch("Kion1E_MuMuKK", &k1E_MuMuKK);
  X_One_Tree_->Branch("kaon1_nsigdedx", &kaon1_nsigdedx);
  X_One_Tree_->Branch("kaon1_dedx", &kaon1_dedx);
  X_One_Tree_->Branch("kaon1_dedxMass", &kaon1_dedxMass);
  X_One_Tree_->Branch("kaon1_theo", &kaon1_theo);
  X_One_Tree_->Branch("kaon1_sigma", &kaon1_sigma);
  X_One_Tree_->Branch("kaon1_dedx_byHits", &kaon1_dedx_byHits);
  X_One_Tree_->Branch("kaon1_dedxErr_byHits", &kaon1_dedxErr_byHits);
  X_One_Tree_->Branch("kaon1_saturMeas_byHits", &kaon1_saturMeas_byHits);
  X_One_Tree_->Branch("kaon1_Meas_byHits", &kaon1_Meas_byHits);
  X_One_Tree_->Branch("Kaon2Px_MuMuKK", &k2Px_MuMuKK);
  X_One_Tree_->Branch("Kaon2Py_MuMuKK", &k2Py_MuMuKK);
  X_One_Tree_->Branch("Kaon2Pz_MuMuKK", &k2Pz_MuMuKK);
  X_One_Tree_->Branch("Kaon2E_MuMuKK", &k2E_MuMuKK);
  X_One_Tree_->Branch("kaon2_nsigdedx", &kaon2_nsigdedx);
  X_One_Tree_->Branch("kaon2_dedx", &kaon2_dedx);
  X_One_Tree_->Branch("kaon2_dedxMass", &kaon2_dedxMass);
  X_One_Tree_->Branch("kaon2_theo", &kaon2_theo);
  X_One_Tree_->Branch("kaon2_sigma", &kaon2_sigma);
  X_One_Tree_->Branch("kaon2_dedx_byHits", &kaon2_dedx_byHits);
  X_One_Tree_->Branch("kaon2_dedxErr_byHits", &kaon2_dedxErr_byHits);
  X_One_Tree_->Branch("kaon2_saturMeas_byHits", &kaon2_saturMeas_byHits);
  X_One_Tree_->Branch("kaon2_Meas_byHits", &kaon2_Meas_byHits);

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

// rsync -vut --existing src/MuMuKKPAT.cc semrat@lxplus.cern.ch:/afs/cern.ch/user/s/semrat/scratch0/CMSSW_5_3_22/src/X4140/MuMuKKPAT/src/MuMuKKPAT.cc

