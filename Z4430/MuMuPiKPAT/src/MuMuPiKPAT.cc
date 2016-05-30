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


// system include files
#include <memory>

// user include files
#include "../interface/MuMuPiKPAT.h"
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

//for 53x 
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"

#include <vector>
#include <utility>

#include "TMath.h"
#include "Math/VectorUtil.h"

// useless so far
//#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
//#include "HepMC/GenVertex.h"
//#include <HepMC/GenVertex.h>
//#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 

//
// constants, enums and typedefs
//

typedef math::Error<3>::type CovarianceMatrix;

const ParticleMass muon_mass = 0.10565837; //pdg mass
const ParticleMass pion_mass = 0.13957018; //pdg mass
const ParticleMass kaon_mass = 0.493667;
ParticleMass Jpsi_mass = 3.096916;
ParticleMass psi2S_mass = 3.686093;
const ParticleMass proton_mass = 0.938272046; //pdg mass

//Setting insignificant mass sigma to avoid singularities in the covariance matrix.
float small_sigma = muon_mass*1.e-6;				

//
// static data member definitions
//

//
// constructors and destructor
//
MuMuPiKPAT::MuMuPiKPAT(const edm::ParameterSet& iConfig) :
  hlTriggerResults_(iConfig.getUntrackedParameter<edm::InputTag>("HLTriggerResults",edm::InputTag("TriggerResults::HLT")) ),
  inputGEN_(iConfig.getUntrackedParameter<edm::InputTag>("inputGEN",edm::InputTag("genParticles"))),
  vtxSample(iConfig.getUntrackedParameter<std::string>("VtxSample",std::string("offlinePrimaryVertices"))), 
  ///// 
  doData( iConfig.getUntrackedParameter<bool>("DoDataAnalysis", true) ),
  doMC( iConfig.getUntrackedParameter<bool>("DoMonteCarloTree", true) ),
  MCParticle( iConfig.getUntrackedParameter<int>("MonteCarloParticleId", 20443) ), // 20443 X, 100443 Psi(2S), 9120443 X from B
  MCExclusiveDecay( iConfig.getUntrackedParameter<bool>("MonteCarloExclusiveDecay", true) ),
  MCMother( iConfig.getUntrackedParameter<int>("MonteCarloMotherId", 511) ), // 511 B0 (=anti-B0), 531 Bs0
  MCDaughtersN( iConfig.getUntrackedParameter<int>(" MonteCarloDaughtersN", 3) ), // 3 for exclusive B0->psi'KPi
  ////////
  doMuMuMassConst( iConfig.getUntrackedParameter<bool>("DoMuMuMassConstraint", true) ),

  skipJPsi(iConfig.getUntrackedParameter<bool>("SkipJPsi", false)),
  skipPsi2S(iConfig.getUntrackedParameter<bool>("SkipPsi2S", false)),

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

  JPsiMinMass(iConfig.getUntrackedParameter<double>("MinJPsiMass", 2.7)),
  JPsiMaxMass(iConfig.getUntrackedParameter<double>("MaxJPsiMass", 3.2)),
  PsiPrimeMinMass(iConfig.getUntrackedParameter<double>("MinPsiPrimeMass", 3.)),
  PsiPrimeMaxMass(iConfig.getUntrackedParameter<double>("MaxPsiPrimeMass", 4.)),

  MuMuTrackMaxDR(iConfig.getUntrackedParameter<double>("MaxMuMuTrackDR", 1)),
  // B0
  B0TrackMaxDR(iConfig.getUntrackedParameter<double>("MaxB0CandTrackDR", 1.1)),
  UseB0DR(iConfig.getUntrackedParameter<bool>("UseB0DR", false)),
  MuMuPiKMinMass(iConfig.getUntrackedParameter<double>("MinMuMuPiKMass", 0)),
  MuMuPiKMaxMass(iConfig.getUntrackedParameter<double>("MaxMuMuPiKMass", 10)),
  addB0lessPrimaryVertex_(iConfig.getUntrackedParameter<bool>("addB0lessPrimaryVertex", true)),

  Debug_(iConfig.getUntrackedParameter<bool>("Debug_Output",false)),
  DeDxEstimator_(iConfig.getUntrackedParameter<std::string>("DeDxEstimator", std::string("dedxHarmonic2"))),
  m_dEdxDiscrimTag(iConfig.getUntrackedParameter<std::string>("DeDxDiscriminator", std::string("dedxHarmonic2"))),
//m_dEdxDiscrimTag_kaon(m_dEdxDiscrimTag, std::string("dedxHarmonic2Kaon"))), // not working
  m_dEdxDiscrimTag_kaon(iConfig.getUntrackedParameter<std::string>("DeDxDiscriminator", std::string("dedxHarmonic2"))),
  
  Z_One_Tree_(0),

  runNum(0), evtNum(0), lumiNum(0),

  trigRes(0), trigNames(0), L1TT(0), MatchTriggerNames(0), 

  nMu(0), nMuMu(0), nB0(0), 
  nB0_pre0(0), nB0_pre1(0), nB0_pre2(0), nB0_pre3(0), nB0_pre4(0), nB0_pre5(0), nB0_pre6(0), nB0_pre7(0), nB0_pre8(0), nB0_pre9(0), nB0_pre10(0), nB0_pre11(0), nB0_pre12(0), nB0_pre13(0), nB0_pre14(0),
  priVtx_n(0), priVtx_X(0), priVtx_Y(0), priVtx_Z(0), priVtx_XE(0), priVtx_YE(0), priVtx_ZE(0), priVtx_NormChi2(0), priVtx_Chi2(0), priVtx_CL(0), priVtx_tracks(0), priVtx_tracksPtSq(0), 

  // indices
  mu1Idx(0), mu2Idx(0), MuMuType(0),
  // B0
  B0_MuMuIdx(0), B0_piIdx(0), B0_kIdx(0),   

  // MC Analysis
  n_genEvtVtx(0), genEvtVtx_X(0), genEvtVtx_Y(0), genEvtVtx_Z(0), genEvtVtx_particles(0), n_B0Ancestors(0),
  nMCAll(0), nMCB0(0), /*nMCB0Vtx(0),*/ MCPdgIdAll(0), MCDanNumAll(0),
  // Gen Primary Vertex
  PriVtxGen_X(0), PriVtxGen_Y(0), PriVtxGen_Z(0), PriVtxGen_EX(0), PriVtxGen_EY(0), PriVtxGen_EZ(0),
  PriVtxGen_Chi2(0), PriVtxGen_CL(0), PriVtxGen_Ndof(0), PriVtxGen_tracks(0),
  MCpsi2SPx(0), MCpsi2SPy(0), MCpsi2SPz(0),
  MCmupPx(0), MCmupPy(0), MCmupPz(0), 
  MCmumPx(0), MCmumPy(0), MCmumPz(0), 
  MCpionPx(0), MCpionPy(0), MCpionPz(0), 
  MCkaonPx(0), MCkaonPy(0), MCkaonPz(0),
  MCpionCh(0), MCkaonCh(0),
  MCPx(0), MCPy(0), MCPz(0), 

  // generic muons
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

  // generic tracks
  trNotRef(0), trRef(0),
  //
  trPx(0), trPy(0), trPz(0), trE(0),
  trNDF(0), trPhits(0), trShits(0), trChi2(0),

  trD0(0), trD0E(0), trCharge(0),
  trfHits(0), trFirstBarrel(0), trFirstEndCap(0),
  trDzVtx(0), trDxyVtx(0),
  trQualityHighPurity(0), trQualityTight(0),

  tr_nsigdedx(0), tr_dedx(0), tr_dedxMass(0), tr_theo(0), tr_sigma(0),
  // by hits
  tr_dedx_byHits(0), tr_dedxErr_byHits(0), tr_saturMeas_byHits(0), tr_Meas_byHits(0),

  // MuMu
  MuMuMass(0), MuMuPx(0), MuMuPy(0), MuMuPz(0),
  MuMuVtx_CL(0), MuMuVtx_Chi2(0), 
  MuMuDecayVtx_X(0), MuMuDecayVtx_Y(0), MuMuDecayVtx_Z(0),
  MuMuDecayVtx_XE(0), MuMuDecayVtx_YE(0), MuMuDecayVtx_ZE(0),
  MuMuMuonTrigMatch(0),
  // muons after Psi' fit
  mu1_MuMu_Px(0), mu1_MuMu_Py(0), mu1_MuMu_Pz(0), mu1_MuMu_Chi2(0), mu1_MuMu_NDF(0),
  mu2_MuMu_Px(0), mu2_MuMu_Py(0), mu2_MuMu_Pz(0), mu2_MuMu_Chi2(0), mu2_MuMu_NDF(0),
  // Primary Vertex with "MuMu correction"
  PriVtxMuMuCorr_n(0),
  PriVtxMuMuCorr_X(0), PriVtxMuMuCorr_Y(0), PriVtxMuMuCorr_Z(0), PriVtxMuMuCorr_EX(0), PriVtxMuMuCorr_EY(0), PriVtxMuMuCorr_EZ(0),
  PriVtxMuMuCorr_Chi2(0), PriVtxMuMuCorr_CL(0), PriVtxMuMuCorr_tracks(0),
  nTrk(0),
  // B0  
  b0Mass(0), b0Vtx_CL(0), b0Vtx_Chi2(0),
  b0Px(0), b0Py(0), b0Pz(0), b0PxE(0), b0PyE(0), b0PzE(0), 
  b0DecayVtx_X(0), b0DecayVtx_Y(0), b0DecayVtx_Z(0), b0DecayVtx_XE(0), b0DecayVtx_YE(0), b0DecayVtx_ZE(0),
  // Muons and pions/tracks after B0 Cand fit
  mu1Px_MuMuPiK(0), mu1Py_MuMuPiK(0), mu1Pz_MuMuPiK(0), mu1E_MuMuPiK(0),
  mu2Px_MuMuPiK(0), mu2Py_MuMuPiK(0), mu2Pz_MuMuPiK(0), mu2E_MuMuPiK(0),
  piPx_MuMuPiK(0), piPy_MuMuPiK(0), piPz_MuMuPiK(0), piE_MuMuPiK(0),
  pion_nsigdedx(0), pion_dedx(0), pion_dedxMass(0), pion_theo(0), pion_sigma(0), 
  pion_dedx_byHits(0), pion_dedxErr_byHits(0), pion_saturMeas_byHits(0), pion_Meas_byHits(0), 
  kPx_MuMuPiK(0), kPy_MuMuPiK(0), kPz_MuMuPiK(0), kE_MuMuPiK(0), 
  kaon_nsigdedx(0), kaon_dedx(0), kaon_dedxMass(0), kaon_theo(0), kaon_sigma(0), 
  kaon_dedx_byHits(0), kaon_dedxErr_byHits(0), kaon_saturMeas_byHits(0), kaon_Meas_byHits(0), 
  // Primary Vertex with largest B0_cos(alpha)
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
  // Primary Vertex with "B0 correction"
  PriVtxB0Corr_n(0),
  PriVtxB0Corr_X(0), PriVtxB0Corr_Y(0), PriVtxB0Corr_Z(0), PriVtxB0Corr_EX(0), PriVtxB0Corr_EY(0), PriVtxB0Corr_EZ(0),
  PriVtxB0Corr_Chi2(0), PriVtxB0Corr_CL(0), PriVtxB0Corr_tracks(0),
  // Lifetime variables
  b0CosAlphaBS(0), b0CosAlpha3DBS(0), b0CTauBS(0), b0CTauBSE(0), b0LxyBS(0), b0LxyBSE(0), b0LxyzBS(0), b0LxyzBSE(0),
  b0CosAlphaPV(0), b0CosAlpha3DPV(0), b0CTauPV(0), b0CTauPVE(0), b0LxyPV(0), b0LxyPVE(0), b0LxyzPV(0), b0LxyzPVE(0), 
  b0CosAlphaPVCosAlpha(0), b0CosAlpha3DPVCosAlpha(0), b0CTauPVCosAlpha(0), b0CTauPVCosAlphaE(0), b0LxyPVCosAlpha(0), b0LxyPVCosAlphaE(0), b0LxyzPVCosAlpha(0), b0LxyzPVCosAlphaE(0), 
  b0CosAlphaPVCosAlpha3D(0), b0CosAlpha3DPVCosAlpha3D(0), b0CTauPVCosAlpha3D(0), b0CTauPVCosAlpha3DE(0), b0LxyPVCosAlpha3D(0), b0LxyPVCosAlpha3DE(0), b0LxyzPVCosAlpha3D(0), b0LxyzPVCosAlpha3DE(0), 
  b0CosAlphaB0LessPV(0), b0CosAlpha3DB0LessPV(0), b0CTauB0LessPV(0), b0CTauB0LessPVE(0), b0LxyB0LessPV(0), b0LxyB0LessPVE(0), b0LxyzB0LessPV(0), b0LxyzB0LessPVE(0), 
  b0CosAlphaB0LessPVCosAlpha(0), b0CosAlpha3DB0LessPVCosAlpha(0), b0CTauB0LessPVCosAlpha(0), b0CTauB0LessPVCosAlphaE(0), b0LxyB0LessPVCosAlpha(0), b0LxyB0LessPVCosAlphaE(0), b0LxyzB0LessPVCosAlpha(0), b0LxyzB0LessPVCosAlphaE(0), 
  b0CosAlphaB0LessPVCosAlpha3D(0), b0CosAlpha3DB0LessPVCosAlpha3D(0), b0CTauB0LessPVCosAlpha3D(0), b0CTauB0LessPVCosAlpha3DE(0), b0LxyB0LessPVCosAlpha3D(0), b0LxyB0LessPVCosAlpha3DE(0), b0LxyzB0LessPVCosAlpha3D(0), b0LxyzB0LessPVCosAlpha3DE(0), 
  b0CosAlphaPVX(0), b0CTauPVX(0), b0CTauPVXE(0), b0LxyPVX(0), b0LxyzPVX(0), b0LxyzPVXE(0), 
  b0CTauPVX_3D(0), b0CTauPVX_3D_err(0),

  PiPiMass_err(0)

{
  //string DoubleMu0 = "HLT_DoubleMu0"; 
  //TriggersForMatching.push_back(DoubleMu0);
  //string DoubleMu0_Quark = "HLT_DoubleMuO_Quarkonium_V1";
  //TriggersForMatching.push_back(DoubleMu0_Quark);

  //now do what ever initialization is needed
  MuMuMinMass = JPsiMinMass;
  MuMuMaxMass = PsiPrimeMaxMass;

}


MuMuPiKPAT::~MuMuPiKPAT()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
	
}


//
// member functions
//

// ------------ method called to for each event  ------------
void MuMuPiKPAT::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{	
  // get event content information
  bool decayChainOK = false;
	
  runNum = iEvent.id().run();
  evtNum = iEvent.id().event();
  lumiNum = iEvent.id().luminosityBlock();

  if (evtNum != 221385081) {
    return ;
    cout <<"\nafter" <<endl ;
  }
  bool hasRequestedTrigger = false;
  ESHandle<MagneticField> bFieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);
	
  // first get HLT results
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
		
    // get hold of trigger names - based on TriggerResults object!
    edm::TriggerNames triggerNames_;
    //     triggerNames_.init(*hltresults); 
    triggerNames_ = iEvent.triggerNames(*hltresults);
	
    int ntriggers = TriggersForMatching_.size();
    for (int MatchTrig = 0; MatchTrig < ntriggers; MatchTrig++) { // initialize MatchingTriggerResult array
      MatchingTriggerResult[MatchTrig] = 0;
    }

    for (int itrig = 0; itrig < ntrigs; itrig++) {
      string trigName = triggerNames_.triggerName(itrig);
      int hltflag = (*hltresults)[itrig].accept();
      //if (hltflag && Debug_) cout << trigName << " " <<hltflag <<endl;
      if (hltflag) cout << trigName << " " <<hltflag <<endl;
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
    //int ntriggers = TriggersForMatching_.size();
    for (int MatchTrig = 0; MatchTrig<ntriggers; MatchTrig++){
      //      cout << TriggersForMatching_[MatchTrig]<<endl;
      MatchTriggerNames->push_back(TriggersForMatching_[MatchTrig]);
    }
    
    //
    // Get HLT map : triggername associated with its prescale, saved only for accepted trigger
    //
    for (unsigned int i=0; i<triggerNames_.size(); i++){
      if ( hltresults->accept(i) ) { //  save trigger info only for accepted paths
	// get the prescale from the HLTConfiguration, initialized at beginRun
	int prescale = hltConfig_.prescaleValue(iEvent,iSetup,triggerNames_.triggerNames().at(i));
	//std::cout<<" HLT===> "<<triggerNames_.triggerNames().at(i)<<" prescale ="<<prescale<<std::endl;
	HLTPreScaleMap[triggerNames_.triggerNames().at(i)] = prescale;
      }
    }
    HLTTrig = &HLTPreScaleMap; // store in the branch

  } // end valid trigger

	
  // get L1 trigger info
	
  edm::ESHandle<L1GtTriggerMenu> menuRcd;
  iSetup.get<L1GtTriggerMenuRcd>().get(menuRcd) ;

  //const L1GtTriggerMenu* menu = menuRcd.product();
	
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
      thePrimaryVtx = Vertex(*(recVtxs->begin()));
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

  RefVtx = thePrimaryVtx.position(); //reference primary vertex choosen

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
	
  //  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // try reconstruction without fitting modules
  //  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Handle< vector<pat::GenericParticle> > thePATTrackHandle;
  iEvent.getByLabel("cleanPatTrackCands", thePATTrackHandle);
  Handle< vector<pat::GenericParticle> > theKaonRefittedPATTrackHandle;
  iEvent.getByLabel("cleanPatTrackKaonCands", theKaonRefittedPATTrackHandle);

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
    //iEvent.getByLabel("dedxHarmonic2", elossCollection);
    iEvent.getByLabel(DeDxEstimator_, elossCollection);
    energyLoss = elossCollection.product();
  } catch ( cms::Exception& ex ) {
    if (evtNum < 100) edm::LogError("Analyzer") <<"Warning can't get collection with label: elossCollection";
    iexception_dedx = 1;
  }

  // dE/dx hits
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


  /////// check MC truth
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
    // (then, for example, you can do)
    if (Debug_) cout <<"\nBunchXing multiplicity = " <<PupInfo->size() <<endl ;
    for (PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI)
      if (Debug_) cout <<"Pileup Information: bunchXing, nvtx: " <<PVI->getBunchCrossing() <<" " <<PVI->getPU_NumInteractions() <<endl;

    Handle<GenParticleCollection> genParticles;
    iEvent.getByLabel("genParticles", genParticles);
    
    if (Debug_) cout << "############### GenParticles Analysis ###############" << endl;
    float psi2SPx=0., psi2SPy=0., psi2SPz=0., mupPx=0., mupPy=0., mupPz=0., mumPx=0., mumPy=0., mumPz=0., pionPx=0., pionPy=0., pionPz=0., kaonPx=0., kaonPy=0., kaonPz=0.;
    int pionCh=0, kaonCh=0 ;

    for (size_t i = 0; i < genParticles->size(); ++ i) {
      nMCAll++;
      const reco::GenParticle &p = (*genParticles)[i];
      int pdgid = p.pdgId() ; 
      int dauNum = p.numberOfDaughters();
      MCPdgIdAll->push_back( pdgid );
      MCDanNumAll->push_back( dauNum );

      if ( MCExclusiveDecay ) {
	// check if there is a MCMother which has MCDaughtersN daughters
	if ( abs(pdgid) == MCMother  &&  dauNum == MCDaughtersN ) {
	  bool mumuOK = false, pionOK = false, kaonOK = false;

	  for (int j=0; j<dauNum; ++j) {
	    const Candidate *dau = p.daughter(j);
	    if (Debug_) cout << "dauPdgId = " << dau->pdgId() << endl;
	  
	    // check if one of B0 daughters is a psi(nS) whitch has 2 muons as daughters
	    int mumuId = 0 ;
	    if (skipJPsi && skipPsi2S) 
	      cout <<"Skipping both J/psi and psi(2S)!" <<endl ;
	    else if (skipJPsi)
	      mumuId = 100443 ;
	    else if (skipPsi2S)
	      mumuId = 443 ;
 
	    if ( ((skipJPsi || skipPsi2S) && (dau->pdgId() == mumuId)) ||
		 ((!skipJPsi && !skipPsi2S) && (dau->pdgId()%1000 == 443)) ) {
	      psi2SPx = dau->px(); psi2SPy = dau->py(); psi2SPz = dau->pz();
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
	      
	    } // end check if one of the MCMother daughters is a J/Psi or psi'
	    
	    else if ( abs(dau->pdgId()) == 211 ) { // check if one of B0 daughters is a pion
	      pionPx = dau->px(); pionPy = dau->py(); pionPz = dau->pz();
	      pionCh = (dau->pdgId() == 211)? 1 : -1;
	      pionOK = true;
	    } else if ( abs(dau->pdgId()) == 321 ) { // check if one of B0 daughters is a kaon
	      kaonPx = dau->px(); kaonPy=dau->py(); kaonPz=dau->pz();
	      kaonCh = (dau->pdgId() == 321)? 1 : -1;
	      kaonOK = true;
	    }
	    
	  } // end loop on MCMother daughters

	  if (Debug_) cout << "mumuOK = " << mumuOK << ", pionOK = " << pionOK << ", kaonOK = " << kaonOK << endl;
	  
	  if ( mumuOK && pionOK && kaonOK ) {
	    if (Debug_) {
	      cout <<"\nnumber of B0 mothers = " <<p.numberOfMothers() <<endl ;
	      cout <<"B0 mother pdgID = " <<p.mother(0)->pdgId() <<endl ;
	    }
	    ++nMCB0 ;		  
	    //const Point genVertex = p.vertex();
	    //if (genVertex != 0) {
	      //++nMCB0Vtx ;
	      PriVtxGen_X->push_back( p.vx() ) ;
	      PriVtxGen_Y->push_back( p.vy() ) ; 
	      PriVtxGen_Z->push_back( p.vz() ) ;
	      //PriVtxGen_EX->push_back( no method ) ;
	      //PriVtxGen_EY->push_back( no method ) ;
	      //PriVtxGen_EZ->push_back( no method ) ;
	      PriVtxGen_CL->push_back( p.vertexNormalizedChi2() ) ;
	      PriVtxGen_Chi2->push_back( p.vertexChi2() ) ;
	      PriVtxGen_Ndof->push_back( p.vertexNdof() ) ;
	      //PriVtxGen_tracks->push_back( no method ) ;
	      //}
	      
	      Bool_t status = kTRUE ;
	      const Candidate *b0_ancestor = p.mother(0) ; // a particle can have several mothers
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


	    MCpsi2SPx->push_back(psi2SPx); MCpsi2SPy->push_back(psi2SPy); MCpsi2SPz->push_back(psi2SPz);
	    MCmupPx->push_back(mupPx); MCmupPy->push_back(mupPy); MCmupPz->push_back(mupPz);
	    MCmumPx->push_back(mumPx); MCmumPy->push_back(mumPy); MCmumPz->push_back(mumPz);
	    MCpionPx->push_back(pionPx); MCpionPy->push_back(pionPy); MCpionPz->push_back(pionPz);
	    MCkaonPx->push_back(kaonPx); MCkaonPy->push_back(kaonPy); MCkaonPz->push_back(kaonPz);
	    MCpionCh->push_back(pionCh) ; MCkaonCh->push_back(kaonCh) ;
	    decayChainOK = true;
	    //
	    MCPx->push_back( p.px() );
	    MCPy->push_back( p.py() );
	    MCPz->push_back( p.pz() );
	  }
	  if (Debug_) cout << "decayChainOK = " << decayChainOK << endl;
	} // if ( abs(pdgid) == MCMother  &&  dauNum == 3 )
      } // if ( !MCExclusiveDecay )

    } // for (size_t i = 0; i < genParticles->size(); ++ i)
  } // if (doMC)
  //
  // reconstruction only for events with B decaying in psi(nS)+Pi+K
  if ( (doMC && !MCExclusiveDecay) || (doMC && (MCExclusiveDecay && decayChainOK)) || doData ) {
	
    bool isEventWithInvalidMu = false;
	
    if (Debug_) cout << "starting event with " << thePATTrackHandle->size() << " tracks, and " << thePATMuonHandle->size() << " muons" << endl;
    
    if ((thePATMuonHandle->size()) * (thePATTrackHandle->size()) > 20000) {
      cout << "Too many Muons: " << thePATMuonHandle->size() << ", and Tracks: " << thePATTrackHandle->size() << endl;
    } else //if (thePATMuonHandle->size() >= 2) { // check
      if (thePATMuonHandle->size() >= 2  && hasRequestedTrigger) {
	if (Debug_) cout <<"============================  evt: " <<evtNum <<" Accept event with 2 mu and TRIGGER ==============================================" <<endl;
	
	//filling track tree
      	
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
	  // dE/dx hits
	  tr_dedx_byHits->push_back( (dEdxTrack)[tr.track()].dEdx() );
	  tr_dedxErr_byHits->push_back( (dEdxTrack)[tr.track()].dEdxError() );
	  tr_saturMeas_byHits->push_back( (dEdxTrack)[tr.track()].numberOfSaturatedMeasurements() );
	  tr_Meas_byHits->push_back( (dEdxTrack)[tr.track()].numberOfMeasurements() );
	  //  Track quality:            
	  //// loose=0, tight=1, highPurity=2, confirmed=3, goodIterative=4, looseSetWithPV=5, highPuritySetWithPV=6
	  bool ishighPurity = tr.track()->quality(reco::TrackBase::highPurity);
	  trQualityHighPurity->push_back(ishighPurity);
	  trQualityTight->push_back(tr.track()->quality(reco::TrackBase::tight));
	}
	
	//get MuMu cands
	for ( std::vector<pat::Muon>::const_iterator Muon1 = thePATMuonHandle->begin(); Muon1 != thePATMuonHandle->end(); ++Muon1 ) {

	  // push back all muon information
	  ++nMu;
	  const reco::Muon* rmu1 = dynamic_cast<const reco::Muon * >(Muon1->originalObject());
	  muPx->push_back(rmu1->px());
	  muPy->push_back(rmu1->py());
	  muPz->push_back(rmu1->pz());
	  muCharge->push_back(rmu1->charge());
	
	  if (rmu1->track().isNull()) { // rmu->track() returns innerTrack();
	    //	 cout << "no track for " << std::distance(thePATMuonHandle->begin(), Muon1) << " filling defaults" << endl;
	    ///AF
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
	  //
	  muType->push_back(rmu1->type());
	  int qm = 0;
	  for (int qi=1; qi!= 24; ++qi) {
	    if (muon::isGoodMuon(*rmu1, muon::SelectionType(qi)))
	      qm += 1<<qi;
	  }
	  muQual->push_back(qm);
	  muTrack->push_back(-1);// not implemented yet
	
	  // muon cleaning
	  int nOverlapMus = 0, nSharingSegWith = -1;
	  int nSegments1 = rmu1->numberOfMatches(reco::Muon::SegmentArbitration);
	  for ( std::vector<pat::Muon>::const_iterator Muon2 = Muon1+1; Muon2 != thePATMuonHandle->end(); ++Muon2) {
	    const reco::Muon* rmu2 = dynamic_cast<const reco::Muon*>(Muon2->originalObject());
	    if ( isSameMuon(*rmu1, *rmu2)) continue;
	    if ( !muon::isGoodMuon(*rmu2, muon::TMOneStationTight) ) continue;
	    // geometric overlap
	    if ( muon::overlap( *rmu1, *rmu2 ) )
	      nOverlapMus++ ;
	    // shared segments 
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

	  //check for muon1
	  TrackRef muTrack1 = Muon1->track();
	  if ( muTrack1.isNull() )
	    continue;      
	  // cuts on muon1
	  if (rmu1->track()->hitPattern().numberOfValidPixelHits() < MuMinPixHits
	      || rmu1->track()->hitPattern().numberOfValidStripHits() < MuMinSiHits
	      || rmu1->track()->chi2()/rmu1->track()->ndof() > MuMaxNormChi
	      || fabs(rmu1->track()->dxy(RefVtx)) > MuMaxD0) {
	    continue ;
	  }
			
	  //next check for muon2
	  for ( std::vector<pat::Muon>::const_iterator Muon2 = Muon1+1; Muon2 != thePATMuonHandle->end(); ++Muon2) {
	    //checks for muon2
	    if(Muon2->charge() * Muon1->charge() > 0)
	      continue ;
	    //
	    const reco::Muon* rmu2 = dynamic_cast<const reco::Muon *>(Muon2->originalObject()) ;	 
	    if (muon::overlap(*rmu1, *rmu2) )
	      continue ;	
	    //
	    TrackRef muTrack2 = Muon2->track() ;
	    if ( muTrack2.isNull() )
	      continue ;	
	    // cuts on muon2
	    if (rmu2->track()->hitPattern().numberOfValidPixelHits() < MuMinPixHits
		|| rmu2->track()->hitPattern().numberOfValidStripHits() < MuMinSiHits
		|| rmu2->track()->chi2()/rmu1->track()->ndof() > MuMaxNormChi
		|| fabs(rmu2->track()->dxy(RefVtx)) > MuMaxD0) {
	      continue ;
	    }

	    //Get the MuMu information				
	    TransientTrack muon1TT( muTrack1, &(*bFieldHandle) );
	    TransientTrack muon2TT( muTrack2, &(*bFieldHandle) );
				
	    KinematicParticleFactoryFromTransientTrack pFactory;
				
	    //initial chi2 and ndf before kinematic fits.
	    float chi = 0., ndf = 0. ;

	    vector<RefCountedKinematicParticle> muons;
	    muons.push_back( pFactory.particle( muon1TT, muon_mass, chi, ndf, small_sigma));
	    muons.push_back( pFactory.particle( muon2TT, muon_mass, chi, ndf, small_sigma));
				
	    KinematicParticleVertexFitter MuMuFitter;   
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
				
	    //Fill the MuMu vectors
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
				
	    mu1_MuMu_Px->push_back( Mu1Cand_KP.momentum().x());
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
	    } else if ( MuMuCand_fromFit->currentState().mass() > PsiPrimeMinMass  &&  MuMuCand_fromFit->currentState().mass() < PsiPrimeMaxMass ) {
	      dimuonType = 2 ;
	    }
	    //cout <<dimuonType <<endl;

	    if (Debug_) cout <<"evt:" <<evtNum <<" MuMu with diMuonType = " <<dimuonType <<endl;

	    MuMuType->push_back(dimuonType);
	

	    int ntriggers = TriggersForMatching_.size();
	    for (int MatchTrig = 0; MatchTrig < ntriggers; MatchTrig++) 
	      {
		if ( MatchingTriggerResult[MatchTrig]!=0 ) 
		  {
		    pat::TriggerObjectStandAloneCollection mu1HLTMatches = Muon1->triggerObjectMatchesByFilter( FiltersForMatching_[MatchTrig] );
		    pat::TriggerObjectStandAloneCollection mu2HLTMatches = Muon2->triggerObjectMatchesByFilter( FiltersForMatching_[MatchTrig] );
		    bool pass1 = mu1HLTMatches.size() > 0;
		    bool pass2 = mu2HLTMatches.size() > 0;
		  
		    if ((pass1) && (pass2))
		      {
			MuMuMuonTrigMatch->push_back(true);
			if (Debug_) cout <<"Matched MuMu" <<endl ;
		      } else
		      MuMuMuonTrigMatch->push_back(false);
		  }
		else
		  MuMuMuonTrigMatch->push_back(false);
	      }
	  
	    // vertex without matched muons 
	    vector<TransientVertex> pvs ;
	    Vertex MuMuLessPV = thePrimaryVtx ;

	    if (addMuMulessPrimaryVertex_)
	      {
		VertexReProducer revertex(recVtxs, iEvent);
		Handle<TrackCollection> pvtracks;   
		iEvent.getByLabel(revertex.inputTracks(),   pvtracks);
		Handle<BeamSpot>        pvbeamspot;
		iEvent.getByLabel(revertex.inputBeamSpot(), pvbeamspot);
		//
		if ( pvbeamspot.isValid() < 0 ) 
		  continue ; 
		if (pvbeamspot.id() != beamSpotHandle.id()) {
		  edm::LogWarning("Inconsistency") << "The BeamSpot used for PV reco is not the same used in this analyzer.";
		}
		//
		const reco::Muon *rmu_1 = dynamic_cast<const reco::Muon*>( Muon1->originalObject() ) ;
		const reco::Muon *rmu_2 = dynamic_cast<const reco::Muon*>( Muon2->originalObject() ) ;
		//
		if (rmu_1 != 0  &&  rmu_2 != 0  &&  rmu_1->track().id() == pvtracks.id()  &&  rmu_2->track().id() == pvtracks.id() ) { 
		  TrackCollection MuMuLess;
		  MuMuLess.reserve(pvtracks->size());
		  for (size_t i = 0, n = pvtracks->size(); i < n; ++i) {
		    if (i == rmu_1->track().key()) continue;
		    if (i == rmu_2->track().key()) continue;
		    //
		    MuMuLess.push_back((*pvtracks)[i]);
		  }
		  if (Debug_) cout <<"pvbeamspot.isValid() = " <<pvbeamspot.isValid() <<endl ;
		  pvs = revertex.makeVertices(MuMuLess, *pvbeamspot, iSetup) ;
		  //
		  if (!pvs.empty()) {
		    MuMuLessPV = Vertex(pvs.front());
		    //thePrimaryVtx = MuMuLessPV;
		  }
		}
	      }
	  
	    //							
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

	    if (dimuonType == 0) 
	      continue ;
	  	
	    if (Debug_) cout <<"evt:"<<evtNum<< " is Invalid Muon ?  " <<isEventWithInvalidMu << endl;
	    /*
	      if (isEventWithInvalidMu ) { 
	      if (Debug_) cout <<"evt:"<<evtNum<< "invalid muon ... skip event" << endl;
	      continue;
	      }
	    */
	    if ( skipJPsi  &&  (dimuonType == 1) ) // skip J/psi events
	      continue ;           
	    if ( skipPsi2S  &&  (dimuonType == 2) ) // skip Psi(2S) events
	      continue ;           

	    nTrk->push_back( thePATTrackHandle->size() ) ;
	    if (thePATTrackHandle->size() < 2)
	      continue ; 
	    nB0_pre0++ ;
	    /////// cuts on MuMu mass window
	    if (MuMuMass->at(nMuMu-1) < MuMuMinMass  ||  MuMuMass->at(nMuMu-1) > MuMuMaxMass)
	      continue ; nB0_pre1++ ;
	
	    // next check tracks for pion
	    //
	    for ( vector<pat::GenericParticle>::const_iterator Track1 = thePATTrackHandle->begin(); Track1 != thePATTrackHandle->end(); ++Track1 ) {
					
	      //check track doesn't overlap with the MuMu candidate tracks
	      if (Track1->track().key() == rmu1->track().key()  ||  Track1->track().key() == rmu2->track().key())
		continue; nB0_pre2++ ;
	    
	      // cuts on charged tracks	
	      if (( Track1->track()->chi2()/Track1->track()->ndof() > TrMaxNormChi2 )  ||  Track1->pt() < TrMinPt)
		continue ; nB0_pre3++ ;

	      //next check tracks for kaon
	      //for ( vector<pat::GenericParticle>::const_iterator Track2 = Track1+1; Track2 != thePATTrackHandle->end(); ++Track2 )
	      for ( vector<pat::GenericParticle>::const_iterator Track2 = theKaonRefittedPATTrackHandle->begin(); Track2 != theKaonRefittedPATTrackHandle->end(); ++Track2 )
		{
		  //check that this second track doesn't overlap with the the first track candidate
		  if (Track2->track().key() == Track1->track().key())
		    continue ; nB0_pre4++ ;	      
		  //check track doesn't overlap with the MuMu candidate tracks
		  if (Track2->track().key() == rmu1->track().key()  ||  Track2->track().key() == rmu2->track().key())
		    continue ; nB0_pre5++ ;      
		  
		  if (Track1->charge() * Track2->charge() > 0)
		    continue ; nB0_pre6++ ;
	      				  
		  ///// cuts on charged tracks	
		  if ((Track2->track()->chi2() / Track2->track()->ndof() > TrMaxNormChi2)  ||  Track2->pt() < TrMinPt)
		    continue; nB0_pre7++ ;
		
		  // cuts on tracks' delta R								
		  math::XYZTLorentzVector MuMu = (rmu1->p4() + rmu2->p4()); 
		  math::XYZTLorentzVector b0 = (rmu1->p4() + rmu2->p4() + Track1->p4() + Track2->p4()); 
		  //
		  float MuMuPiDR = sqrt( pow(MuMu.eta() - Track1->p4().eta(),2) + pow(MuMu.phi() - Track1->p4().phi(), 2) );
		  float MuMuKDR = sqrt( pow(MuMu.eta() - Track2->p4().eta(),2) + pow(MuMu.phi() - Track2->p4().phi(), 2) );
		
		  float b0PiDR = sqrt( pow(b0.eta() - Track1->p4().eta(),2) + pow(b0.phi() - Track1->p4().phi(), 2));
		  float b0KDR = sqrt( pow(b0.eta() - Track2->p4().eta(),2) + pow(b0.phi() - Track2->p4().phi(), 2));

		  
		  // B0
		  if (UseB0DR) {
		    if (b0PiDR > B0TrackMaxDR  ||  b0KDR > B0TrackMaxDR)
		      continue ; // B0TrackMaxDR = 2
		  } else {
		    if (MuMuPiDR > MuMuTrackMaxDR  ||  MuMuKDR > MuMuTrackMaxDR)
		      continue ; // MuMuTrackMaxDR = 3.5
		  }
		  nB0_pre8++ ;
		  ////////// cuts on MuMu+pi+k mass window
		  if ((Track1->p4() + Track2->p4() + MuMu).M() > MuMuPiKMaxMass  ||  (Track1->p4() + Track2->p4() + MuMu).M() < MuMuPiKMinMass)
		    continue ; nB0_pre9++ ;
		
		  // having two oppositely charged muons, and two oppositely charged tracks: try to vertex them
		  TransientTrack pionTT( Track1->track(), &(*bFieldHandle) ); 
		  TransientTrack kaonTT( Track2->track(), &(*bFieldHandle) );
		
		  TransientTrack kaonTT_notRefit ;
		  Bool_t notRefittedPartner = false ;
		  for ( vector<pat::GenericParticle>::const_iterator Track2_notRefit = thePATTrackHandle->begin(); Track2_notRefit != thePATTrackHandle->end(); ++Track2_notRefit )
		    if ( Track2_notRefit->track().key() == Track2->track().key() ) {
		      notRefittedPartner = true ;
		      kaonTT_notRefit = TransientTrack( Track2_notRefit->track(), &(*bFieldHandle) ) ;
		      break ;
		    }
			      
		  // Do mass constraint for MuMu cand and do mass constrained vertex fit
		  vector<RefCountedKinematicParticle> b0Daughters;
		  b0Daughters.push_back(pFactory.particle( muon1TT, muon_mass, chi, ndf, small_sigma));
		  b0Daughters.push_back(pFactory.particle( muon2TT, muon_mass, chi, ndf, small_sigma));
		  //b0Daughters.push_back(pFactory.particle( pionTT, pion_mass, chi, ndf, small_sigma));
		  //b0Daughters.push_back(pFactory.particle( kaonTT, kaon_mass, chi, ndf, small_sigma));
		  b0Daughters.push_back(pFactory.particle( pionTT, proton_mass, chi, ndf, small_sigma));
		  b0Daughters.push_back(pFactory.particle( kaonTT, proton_mass, chi, ndf, small_sigma));
		
		  RefCountedKinematicTree B0VertexFitTree, B0VertexFitTree_noKrefit ;
		  KinematicConstrainedVertexFitter B0Fitter ;
		
		  if (doMuMuMassConst) { // MassConst = 'MC' in the following
		    //MultiTrackKinematicConstraint *MuMu = new  TwoTrackMassKinematicConstraint(psi2S_mass);
		    MultiTrackKinematicConstraint *MuMu = 0;
		    if (dimuonType == 1) { // constrain to JPsi mass
		      MuMu = new TwoTrackMassKinematicConstraint(Jpsi_mass);
		    } else if (dimuonType == 2) { // constrain to Psi(2S) mass
		      MuMu = new TwoTrackMassKinematicConstraint(psi2S_mass);
		    } // already asked for: if (dimuonType == 0) continue ;
		  
		    B0VertexFitTree = B0Fitter.fit( b0Daughters, MuMu );
		    if (notRefittedPartner) { // use not refitted kaons
		      b0Daughters.pop_back() ;
		      b0Daughters.push_back(pFactory.particle( kaonTT_notRefit, kaon_mass, chi, ndf, small_sigma));
		      B0VertexFitTree_noKrefit = B0Fitter.fit( b0Daughters, MuMu );
		    }
		  } 
		  else {
		    B0VertexFitTree = B0Fitter.fit( b0Daughters );
		    if (notRefittedPartner) { // use not refitted kaons
		      b0Daughters.pop_back() ;
		      b0Daughters.push_back(pFactory.particle( kaonTT_notRefit, kaon_mass, chi, ndf, small_sigma));
		      B0VertexFitTree_noKrefit = B0Fitter.fit( b0Daughters );
		    }
		  }
		  /*
		  if ( !B0VertexFitTree->isValid() )
		    if (Debug_) cout <<"B0VertexFitTree is NOT valid!" <<endl ;
		  if ( notRefittedPartner  &&  !B0VertexFitTree_noKrefit->isValid() )
		    if (Debug_) cout <<"B0VertexFitTree_noKrefit is NOT valid!" <<endl ;
		  */
		  if ( !B0VertexFitTree->isValid() )
		    //if ( !B0VertexFitTree->isValid() || !B0VertexFitTree_noKrefit->isValid() )
		    continue ; nB0_pre10++ ;      
		  B0VertexFitTree->movePointerToTheTop();
		  RefCountedKinematicParticle B0Cand_fromMCFit = B0VertexFitTree->currentParticle();
		  //B0VertexFitTree_noKrefit->movePointerToTheTop();
		  //RefCountedKinematicParticle B0Cand_noKrefit_fromMCFit = B0VertexFitTree_noKrefit->currentParticle();

		  RefCountedKinematicVertex B0Cand_vertex_fromMCFit = B0VertexFitTree->currentDecayVertex();	      
		  //
		  
		  if ( !B0Cand_vertex_fromMCFit->vertexIsValid() )
		    continue ; nB0_pre11++ ;
		  		  
		  if ( B0Cand_vertex_fromMCFit->chiSquared() < 0  ||  B0Cand_vertex_fromMCFit->chiSquared() > 10000 )
		    continue ; nB0_pre12++ ;
		  
		  if ( B0Cand_fromMCFit->currentState().mass() > 100 )
		    continue ; nB0_pre13++ ;
		  
		  double b0VtxProb = ChiSquaredProbability((double)(B0Cand_vertex_fromMCFit->chiSquared()), (double)(B0Cand_vertex_fromMCFit->degreesOfFreedom()));
		  if ( b0VtxProb < 0.005 ) //0.0001 ) 
		    continue ; nB0_pre14++ ;
		
		  ///////////////////////////////////						
													
		  //////////////////// Lifetimes calculations ////////////////////
		  TVector3 B0_vtx((*B0Cand_vertex_fromMCFit).position().x(), (*B0Cand_vertex_fromMCFit).position().y(), 0) ;			
		  TVector3 B0_pperp(B0Cand_fromMCFit->currentState().globalMomentum().x(), B0Cand_fromMCFit->currentState().globalMomentum().y(), 0);
		  TVector3 B0_vtx3D((*B0Cand_vertex_fromMCFit).position().x(), (*B0Cand_vertex_fromMCFit).position().y(), (*B0Cand_vertex_fromMCFit).position().z()) ;
		  TVector3 B0_pperp3D(B0Cand_fromMCFit->currentState().globalMomentum().x(), B0Cand_fromMCFit->currentState().globalMomentum().y(), B0Cand_fromMCFit->currentState().globalMomentum().z());
		  // needed by Similarity method
		  //AlgebraicVector B0_vpperp(3); 
		  //B0_vpperp[0] = B0_pperp.x(); B0_vpperp[1] = B0_pperp.y(); B0_vpperp[2] = 0.;
		  AlgebraicVector3 B0_v3pperp ;
		  B0_v3pperp[0] = B0_pperp.x(); B0_v3pperp[1] = B0_pperp.y(); B0_v3pperp[2] = 0.;

		  TVector3 B0_pvtx, B0_pvtx3D, B0_vdiff, B0_vdiff3D ;
		  double B0_cosAlpha, B0_cosAlpha3D, B0_ctau ;
		  VertexDistanceXY B0_vdistXY ;
		  Measurement1D B0_distXY ;
		  GlobalError B0_v1e = (Vertex(*B0Cand_vertex_fromMCFit)).error();
		  GlobalError B0_v2e ;
		  //AlgebraicSymMatrix B0_vXYe ;
		  AlgebraicSymMatrix33 B0_vXYe ;
		  double B0_ctauErr ;
		  float B0_lxy, B0_lxyErr, B0_lxyz, B0_lxyzErr ;
		  ROOT::Math::SVector<double, 3> B0_vDiff, B0_vDiff3D ; // needed by Similarity method


		  // Lifetime wrt PV
		  B0_v2e = thePrimaryVtx.error();
		  B0_vXYe = B0_v1e.matrix() + B0_v2e.matrix() ;
		  ///// 2D
		  B0_pvtx.SetXYZ(thePrimaryVtx.position().x(), thePrimaryVtx.position().y(), 0) ;
		  B0_vdiff = B0_vtx - B0_pvtx ;
		  B0_cosAlpha = B0_vdiff.Dot(B0_pperp) / (B0_vdiff.Perp()*B0_pperp.Perp()) ;
		  //B0_lxy = B0_vdiff.Dot(B0_pperp)/B0_pperp.Mag();
		  B0_lxy = B0_vdiff.Perp();
		  B0_vDiff[0] = B0_vdiff.x(); B0_vDiff[1] = B0_vdiff.y(); B0_vDiff[2] = 0 ; // needed by Similarity method
		  B0_lxyErr = sqrt(ROOT::Math::Similarity(B0_vDiff,B0_vXYe)) / B0_vdiff.Perp();
		  //
		  B0_distXY = B0_vdistXY.distance(Vertex(*B0Cand_vertex_fromMCFit), Vertex(thePrimaryVtx));
		  B0_ctau = B0_distXY.value() * B0_cosAlpha * B0Cand_fromMCFit->currentState().mass() / B0_pperp.Perp();
		  //B0_ctauErrPV = sqrt(B0_vXYe.similarity(B0_vpperp)) * B0Cand_fromMCFit->currentState().mass() / (B0_pperp.Perp2()) ;
		  B0_ctauErr = sqrt(ROOT::Math::Similarity(B0_v3pperp,B0_vXYe)) * B0Cand_fromMCFit->currentState().mass() / (B0_pperp.Perp2()) ;
		  ///// 3D
		  B0_pvtx3D.SetXYZ(thePrimaryVtx.position().x(), thePrimaryVtx.position().y(), thePrimaryVtx.position().z());
		  B0_vdiff3D = B0_vtx3D - B0_pvtx3D;
		  B0_cosAlpha3D = B0_vdiff3D.Dot(B0_pperp3D)/(B0_vdiff3D.Mag()*B0_pperp3D.Mag());
		  //B0_lxyz = B0_vdiff3D.Dot(B0_pperp3D)/B0_pperp3D.Mag();
		  B0_lxyz = B0_vdiff3D.Mag();
		  B0_vDiff3D[0] = B0_vdiff3D.x(); B0_vDiff3D[1] = B0_vdiff3D.y(); B0_vDiff3D[2] = B0_vdiff3D.z() ;
		  B0_lxyzErr = sqrt(ROOT::Math::Similarity(B0_vDiff3D,B0_vXYe)) / B0_vdiff3D.Mag();


		  //////////////////// last cuts ////////////////////
		  if ( !(B0_ctau/B0_ctauErr > 2.8) || !(B0_cosAlpha > 0.8) )
		    continue ;


		  // fill B0 candidate variables now
		  //if (Debug_) cout <<"B0 mass with refitted Kaon = " <<B0Cand_fromMCFit->currentState().mass() <<endl ;
		  //if (Debug_) cout <<"B0 mass with NOT refitted Kaon = " <<B0Cand_noKrefit_fromMCFit->currentState().mass() <<endl ;
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
		  RefCountedKinematicParticle mu1_MuMuPiK = B0VertexFitTree->currentParticle();
		  B0VertexFitTree->movePointerToTheNextChild();
		  RefCountedKinematicParticle mu2_MuMuPiK = B0VertexFitTree->currentParticle();
		  B0VertexFitTree->movePointerToTheNextChild();
		  RefCountedKinematicParticle pi_MuMuPiK = B0VertexFitTree->currentParticle();
		  B0VertexFitTree->movePointerToTheNextChild();
		  RefCountedKinematicParticle k_MuMuPiK = B0VertexFitTree->currentParticle();
		  // muon1			
		  mu1Px_MuMuPiK->push_back( mu1_MuMuPiK->currentState().globalMomentum().x() );
		  mu1Py_MuMuPiK->push_back( mu1_MuMuPiK->currentState().globalMomentum().y() );
		  mu1Pz_MuMuPiK->push_back( mu1_MuMuPiK->currentState().globalMomentum().z() );
		  mu1E_MuMuPiK->push_back( mu1_MuMuPiK->currentState().kinematicParameters().energy() );
		  // muon2
		  mu2Px_MuMuPiK->push_back( mu2_MuMuPiK->currentState().globalMomentum().x() );
		  mu2Py_MuMuPiK->push_back( mu2_MuMuPiK->currentState().globalMomentum().y() );
		  mu2Pz_MuMuPiK->push_back( mu2_MuMuPiK->currentState().globalMomentum().z() );
		  mu2E_MuMuPiK->push_back( mu2_MuMuPiK->currentState().kinematicParameters().energy() );
		  // pion
		  piPx_MuMuPiK->push_back( pi_MuMuPiK->currentState().globalMomentum().x() );
		  piPy_MuMuPiK->push_back( pi_MuMuPiK->currentState().globalMomentum().y() );
		  piPz_MuMuPiK->push_back( pi_MuMuPiK->currentState().globalMomentum().z() );
		  piE_MuMuPiK->push_back( pi_MuMuPiK->currentState().kinematicParameters().energy() );
		
		  // dE/dx	
		  Double_t theo = 0., sigma = 0. ;
		  pion_nsigdedx->push_back( nsigmaofdedx(Track1->track(),theo,sigma) );
		  pion_dedx->push_back( getEnergyLoss(Track1->track()) );
		  pion_dedxMass->push_back( GetMass(Track1->track()) );
		  pion_theo->push_back( theo );
		  pion_sigma->push_back( sigma );
		  // dE/dx hits
		  pion_dedx_byHits->push_back( (dEdxTrack)[Track1->track()].dEdx() );
		  pion_dedxErr_byHits->push_back( (dEdxTrack)[Track1->track()].dEdxError() );
		  pion_saturMeas_byHits->push_back( (dEdxTrack)[Track1->track()].numberOfSaturatedMeasurements() );
		  pion_Meas_byHits->push_back( (dEdxTrack)[Track1->track()].numberOfMeasurements() );
		  // kaon
		  kPx_MuMuPiK->push_back( k_MuMuPiK->currentState().globalMomentum().x() );
		  kPy_MuMuPiK->push_back( k_MuMuPiK->currentState().globalMomentum().y() );
		  kPz_MuMuPiK->push_back( k_MuMuPiK->currentState().globalMomentum().z() );
		  kE_MuMuPiK->push_back( k_MuMuPiK->currentState().kinematicParameters().energy() );
		  // dE/dx
		  /*theo = 0.; sigma = 0. ;
		    kaon_nsigdedx->push_back(nsigmaofdedx(Track2->track(),theo,sigma));
		    kaon_dedx->push_back(getEnergyLoss(Track2->track()));
		    kaon_dedxMass->push_back(GetMass(Track2->track()));
		    kaon_theo->push_back(theo);
		    kaon_sigma->push_back(sigma);*/
		  // dE/dx hits
		  kaon_dedx_byHits->push_back( (dEdxTrack_Kaon)[Track2->track()].dEdx() );
		  kaon_dedxErr_byHits->push_back( (dEdxTrack_Kaon)[Track2->track()].dEdxError() );
		  kaon_saturMeas_byHits->push_back( (dEdxTrack_Kaon)[Track2->track()].numberOfSaturatedMeasurements() );
		  kaon_Meas_byHits->push_back( (dEdxTrack_Kaon)[Track2->track()].numberOfMeasurements() );

		  // PV
		  b0CosAlphaPV->push_back( B0_cosAlpha ); b0CosAlpha3DPV->push_back( B0_cosAlpha3D );
		  b0CTauPV->push_back( B0_ctau ); b0CTauPVE->push_back( B0_ctauErr );
		  b0LxyPV->push_back( B0_lxy ); b0LxyPVE->push_back( B0_lxyErr );
		  b0LxyzPV->push_back( B0_lxyz ); b0LxyzPVE->push_back( B0_lxyzErr );


		  ////////// Lifetime wrt BS
		  B0_v2e = theBeamSpotVtx.error();
		  B0_vXYe = B0_v1e.matrix() + B0_v2e.matrix();
		  ///// 2D
		  B0_pvtx.SetXYZ(theBeamSpotVtx.position().x(), theBeamSpotVtx.position().y(), 0);
		  B0_vdiff = B0_vtx - B0_pvtx;
		  B0_cosAlpha = B0_vdiff.Dot(B0_pperp)/(B0_vdiff.Perp()*B0_pperp.Perp());
		  //B0_lxy = B0_vdiff.Dot(B0_pperp)/B0_pperp.Mag();
		  B0_lxy = B0_vdiff.Perp();
		  B0_vDiff[0] = B0_vdiff.x(); B0_vDiff[1] = B0_vdiff.y(); B0_vDiff[2] = 0 ; // needed by Similarity method
		  B0_lxyErr = sqrt(ROOT::Math::Similarity(B0_vDiff,B0_vXYe)) / B0_vdiff.Perp();
		  //
		  B0_distXY = B0_vdistXY.distance(Vertex(*B0Cand_vertex_fromMCFit), Vertex(theBeamSpotVtx));
		  B0_ctau = B0_distXY.value() * B0_cosAlpha * (B0Cand_fromMCFit->currentState().mass() / B0_pperp.Perp()) ;
		  //B0_ctauErr = sqrt(B0_vXYe.similarity(B0_vpperp))*B0Cand_fromMCFit->currentState().mass()/(B0_pperp.Perp2());
		  B0_ctauErr = sqrt(ROOT::Math::Similarity(B0_v3pperp,B0_vXYe)) * B0Cand_fromMCFit->currentState().mass()/B0_pperp.Perp2();
		  ///// 3D
		  B0_pvtx3D.SetXYZ(theBeamSpotVtx.position().x(), theBeamSpotVtx.position().y(), theBeamSpotVtx.position().z());
		  B0_vdiff3D = B0_vtx3D - B0_pvtx3D;
		  B0_cosAlpha3D = B0_vdiff3D.Dot(B0_pperp3D)/(B0_vdiff3D.Mag()*B0_pperp3D.Mag());
		  //B0_lxyz = B0_vdiff3D.Dot(B0_pperp3D)/B0_pperp3D.Mag();
		  B0_lxyz = B0_vdiff3D.Mag();
		  B0_vDiff3D[0] = B0_vdiff3D.x(); B0_vDiff3D[1] = B0_vdiff3D.y(); B0_vDiff3D[2] = B0_vdiff3D.z() ;
		  B0_lxyzErr = sqrt(ROOT::Math::Similarity(B0_vDiff3D,B0_vXYe)) / B0_vdiff3D.Mag();
		  		  				
		  // BS			
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
		      //
		      if (pvbeamspot.id() != beamSpotHandle.id() ) 
			edm::LogWarning("Inconsistency") << "The BeamSpot used for PV reco is not the same used in this analyzer.";
		      //
		      const reco::Muon *B0rmu_1 = dynamic_cast<const reco::Muon *>(Muon1->originalObject());
		      const reco::Muon *B0rmu_2 = dynamic_cast<const reco::Muon *>(Muon2->originalObject());
		      //
		      if (B0rmu_1 != 0  &&  B0rmu_2 != 0  &&  B0rmu_1->track().id() == pvtracks.id()  &&  B0rmu_2->track().id() == pvtracks.id() 
			  &&  Track1->track().id() == pvtracks.id()  &&  Track2->track().id() ==  pvtracks.id()) { 
			//TrackCollection B0Less;
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
			  //B0Less.push_back((*pvtracks)[i]); // need TransientTrack to keep the TrackRef
			  //TransientTrack transientTrack = theTTBuilder->build((*pvtracks)[i]); // need to biuld from a TrackRef to keep the TrackRef
			  reco::TrackRef trk_now(pvtracks, i) ;
			  //B0Less.push_back( *trk_now ); // need to build from a TrackRef to keep the TrackRef
			  //
			  TransientTrack transientTrack = theTTBuilder->build( trk_now ); 
			  transientTrack.setBeamSpot( beamSpot );
			  B0Less.push_back( transientTrack );
			  //
			  //B0Less.push_back( transientTrack.track() ); // same error
			}
			if ( removedTrksPtSq > 0. ) {
			  B0_pvs = revertex.makeVertices(B0Less, *pvbeamspot, iSetup) ; // list of PV
			  //
			  //AdaptiveVertexFitter* theFitter = new AdaptiveVertexFitter();
			  //TransientVertex myVertex = theFitter->vertex(B0Less, beamSpot); // just 1 PV
			  ////TransientVertex myVertex = theFitter->vertex(B0Less);
			  //B0_pvs.push_back( myVertex ) ;
			} else
			  cout <<"\n\\\\\\\\\\\\\\\\\\\\ excluded tracks pT^2 = 0 \\\\\\\\\\\\\\\\\\\\\n" <<endl ;
			if ( !B0_pvs.empty() ) {
			  //sort(B0_pvs.begin(), B0_pvs.end(), VertexHigherPtSquared()); // already called from revertex.makeVertices()
			  B0LessPV = Vertex(B0_pvs.front());
			  //thePrimaryVtx = B0LessPV;
			  
			  B0LessPV_tracksPtSq->push_back( vertexHigherPtSquared.sumPtSquared(B0LessPV) ) ;
			  B0LessPV_4tracksPtSq->push_back( removedTrksPtSq ) ;
			  if (Debug_) {
			    cout <<"\nB0LessPV_z = " <<B0LessPV.position().z() <<endl ;
			    cout <<"B0LessPV_tracks = " <<B0LessPV.tracksSize() <<endl ;
			    cout <<"B0LessPV_tracksPtSq = " <<vertexHigherPtSquared.sumPtSquared(B0LessPV) <<endl ;
			    cout <<"B0LessPV_removedTracksPtSq = " <<removedTrksPtSq <<endl ;
			    cout <<"B0_pvs->size() = " <<B0_pvs.size() <<endl ;
			    //
			    cout <<"priVtx_z = " << priVtx_Z <<endl ;
			    cout <<"priVtx_tracks = " <<priVtx_tracks <<endl ;
			    cout <<"priVtx_tracksPtSq = " <<priVtx_tracksPtSq <<endl ;
			    cout <<"recVtxs->size() = " <<recVtxs->size() <<endl ;
			  }
			}
		      }
		    }
		  //
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
		
		  // Lifetime wrt B0LessPV 
		  B0_v2e = B0LessPV.error();
		  B0_vXYe = B0_v1e.matrix() + B0_v2e.matrix();
		  ///// 2D
		  B0_pvtx.SetXYZ(B0LessPV.position().x(), B0LessPV.position().y(), 0) ;	
		  B0_vdiff = B0_vtx - B0_pvtx ;
		  B0_cosAlpha = B0_vdiff.Dot(B0_pperp)/(B0_vdiff.Perp()*B0_pperp.Perp());
		  //B0_lxy = B0_vdiff.Dot(B0_pperp)/B0_pperp.Mag();
		  B0_lxy = B0_vdiff.Perp();
		  B0_vDiff[0] = B0_vdiff.x(); B0_vDiff[1] = B0_vdiff.y(); B0_vDiff[2] = 0 ; // needed by Similarity method
		  B0_lxyErr = sqrt(ROOT::Math::Similarity(B0_vDiff,B0_vXYe)) / B0_vdiff.Perp();
		  //
		  B0_distXY = B0_vdistXY.distance( Vertex(*B0Cand_vertex_fromMCFit), Vertex(B0LessPV) ) ;
		  B0_ctau = B0_distXY.value() * B0_cosAlpha * B0Cand_fromMCFit->currentState().mass() / B0_pperp.Perp();
		  //B0_ctauErr = sqrt(B0_vXYe.similarity(B0_vpperp))*B0Cand_fromMCFit->currentState().mass()/(B0_pperp.Perp2());
		  B0_ctauErr = sqrt(ROOT::Math::Similarity(B0_v3pperp,B0_vXYe)) * B0Cand_fromMCFit->currentState().mass() / (B0_pperp.Perp2());
		  ///// 3D
		  B0_pvtx3D.SetXYZ(B0LessPV.position().x(), B0LessPV.position().y(), B0LessPV.position().z());
		  B0_vdiff3D = B0_vtx3D - B0_pvtx3D;
		  B0_cosAlpha3D = B0_vdiff3D.Dot(B0_pperp3D)/( B0_vdiff3D.Mag()*B0_pperp3D.Mag() );
		  //B0_lxyz = B0_vdiff3D.Dot(B0_pperp3D)/B0_pperp3D.Mag();
		  B0_lxyz = B0_vdiff3D.Mag();
		  B0_vDiff3D[0] = B0_vdiff3D.x(); B0_vDiff3D[1] = B0_vdiff3D.y(); B0_vDiff3D[2] = B0_vdiff3D.z() ;
		  B0_lxyzErr = sqrt(ROOT::Math::Similarity(B0_vDiff3D,B0_vXYe)) / B0_vdiff3D.Mag();

		  b0CosAlphaB0LessPV->push_back( B0_cosAlpha ) ; b0CosAlpha3DB0LessPV->push_back( B0_cosAlpha3D ) ;
		  b0CTauB0LessPV->push_back( B0_ctau ) ; b0CTauB0LessPVE->push_back( B0_ctauErr ) ;
		  b0LxyB0LessPV->push_back( B0_lxy ) ; b0LxyB0LessPVE->push_back( B0_lxyErr ) ;
		  b0LxyzB0LessPV->push_back( B0_lxyz ) ; b0LxyzB0LessPVE->push_back( B0_lxyzErr ) ;
 

		  // Find the PV among the original offlinePV with the largest B0_cos(alpha)
		  Vertex theCosAlphaV = thePrimaryVtx ; 
		  float maxCosAlpha = -1. ;
		  //
		  for (VertexCollection::const_iterator itv = recVtxs->begin(), itvend = recVtxs->end(); itv != itvend; ++itv) {
		    B0_pvtx.SetXYZ(itv->position().x(), itv->position().y(), 0) ;	
		    B0_vdiff = B0_vtx - B0_pvtx ;
		    float cosAlpha_temp = B0_vdiff.Dot(B0_pperp) / (B0_vdiff.Perp()*B0_pperp.Perp()) ; // Perp() == Mag() when z = 0
		  
		    if ( cosAlpha_temp > maxCosAlpha ) {
		      maxCosAlpha = cosAlpha_temp ;    
		      //thePrimaryVtx = Vertex(*itv);
		      //theCosAlphaV = thePrimaryVtx;
		      theCosAlphaV = Vertex(*itv) ;
		    }
		  }
		  //
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

		  // Lifetime wrt PV with largest B0_cos(alpha) candidate 
		  B0_v2e = theCosAlphaV.error();
		  B0_vXYe = B0_v1e.matrix() + B0_v2e.matrix();
		  ///// 2D
		  B0_pvtx.SetXYZ(theCosAlphaV.position().x(), theCosAlphaV.position().y(), 0) ;	
		  B0_vdiff = B0_vtx - B0_pvtx ;
		  B0_cosAlpha =  maxCosAlpha ;
		  //B0_lxy = B0_vdiff.Dot(B0_pperp)/B0_pperp.Mag();
		  B0_lxy = B0_vdiff.Perp();
		  B0_vDiff[0] = B0_vdiff.x(); B0_vDiff[1] = B0_vdiff.y(); B0_vDiff[2] = 0 ; // needed by Similarity method
		  B0_lxyErr = sqrt(ROOT::Math::Similarity(B0_vDiff,B0_vXYe)) / B0_vdiff.Perp();
		  //
		  B0_distXY = B0_vdistXY.distance( Vertex(*B0Cand_vertex_fromMCFit), Vertex(theCosAlphaV) ) ;
		  B0_ctau = B0_distXY.value() * B0_cosAlpha * B0Cand_fromMCFit->currentState().mass() / B0_pperp.Perp();
		  //B0_ctauErr = sqrt(B0_vXYe.similarity(B0_vpperp))*B0Cand_fromMCFit->currentState().mass()/(B0_pperp.Perp2());
		  B0_ctauErr = sqrt(ROOT::Math::Similarity(B0_v3pperp,B0_vXYe)) * B0Cand_fromMCFit->currentState().mass() / (B0_pperp.Perp2());
		  B0_lxy = B0_vdiff.Dot(B0_pperp) / B0_pperp.Mag() ;
		  ///// 3D
		  B0_pvtx3D.SetXYZ(theCosAlphaV.position().x(), theCosAlphaV.position().y(), theCosAlphaV.position().z());
		  B0_vdiff3D = B0_vtx3D - B0_pvtx3D;
		  B0_cosAlpha3D = B0_vdiff3D.Dot(B0_pperp3D)/( B0_vdiff3D.Mag()*B0_pperp3D.Mag() );
		  //B0_lxyz = B0_vdiff3D.Dot(B0_pperp3D)/B0_pperp3D.Mag();
		  B0_lxyz = B0_vdiff3D.Mag();
		  B0_vDiff3D[0] = B0_vdiff3D.x(); B0_vDiff3D[1] = B0_vdiff3D.y(); B0_vDiff3D[2] = B0_vdiff3D.z() ;
		  B0_lxyzErr = sqrt(ROOT::Math::Similarity(B0_vDiff3D,B0_vXYe)) / B0_vdiff3D.Mag();
		
		  b0CosAlphaPVCosAlpha->push_back( B0_cosAlpha ) ; b0CosAlpha3DPVCosAlpha->push_back( B0_cosAlpha3D ) ;
		  b0CTauPVCosAlpha->push_back( B0_ctau ) ; b0CTauPVCosAlphaE->push_back( B0_ctauErr ) ;
		  b0LxyPVCosAlpha->push_back( B0_lxy ) ; b0LxyPVCosAlphaE->push_back( B0_lxyErr ) ;
		  b0LxyzPVCosAlpha->push_back( B0_lxyz ) ; b0LxyzPVCosAlphaE->push_back( B0_lxyzErr ) ;


		  // Find the PV among the original offlinePV with the largest B0_cos(alpha) 3D
		  Vertex theCosAlpha3DV = thePrimaryVtx ; 
		  float maxCosAlpha3D = -1. ;
		  //
		  for (VertexCollection::const_iterator itv = recVtxs->begin(), itvend = recVtxs->end(); itv != itvend; ++itv) {
		    B0_pvtx3D.SetXYZ(itv->position().x(), itv->position().y(), itv->position().z()) ;	
		    B0_vdiff3D = B0_vtx3D - B0_pvtx3D ;
		    float cosAlpha_temp3D = B0_vdiff3D.Dot(B0_pperp3D) / (B0_vdiff3D.Mag()*B0_pperp3D.Mag()) ;
		  
		    if ( cosAlpha_temp3D > maxCosAlpha3D ) {
		      maxCosAlpha3D = cosAlpha_temp3D ;    
		      //thePrimaryVtx = Vertex(*itv);
		      //theCosAlpha3DV = thePrimaryVtx;
		      theCosAlpha3DV = Vertex(*itv) ;
		    }
		  }
		  //
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

		  // Lifetime wrt PV with largest B0_cos(alpha) 3D candidate 
		  B0_v2e = theCosAlpha3DV.error();
		  B0_vXYe = B0_v1e.matrix() + B0_v2e.matrix();
		  ///// 2D
		  B0_pvtx.SetXYZ(theCosAlpha3DV.position().x(), theCosAlpha3DV.position().y(), 0) ;	
		  B0_vdiff = B0_vtx - B0_pvtx ;
		  B0_cosAlpha = B0_vdiff.Dot(B0_pperp)/(B0_vdiff.Perp()*B0_pperp.Perp()); ;
		  //B0_lxy = B0_vdiff.Dot(B0_pperp)/B0_pperp.Mag();
		  B0_lxy = B0_vdiff.Perp();
		  B0_vDiff[0] = B0_vdiff.x(); B0_vDiff[1] = B0_vdiff.y(); B0_vDiff[2] = 0 ; // needed by Similarity method
		  B0_lxyErr = sqrt(ROOT::Math::Similarity(B0_vDiff,B0_vXYe)) / B0_vdiff.Perp();
		  //
		  B0_distXY = B0_vdistXY.distance( Vertex(*B0Cand_vertex_fromMCFit), Vertex(theCosAlphaV) ) ;
		  B0_ctau = B0_distXY.value() * B0_cosAlpha * B0Cand_fromMCFit->currentState().mass() / B0_pperp.Perp();
		  //B0_ctauErr = sqrt(B0_vXYe.similarity(B0_vpperp))*B0Cand_fromMCFit->currentState().mass()/(B0_pperp.Perp2());
		  B0_ctauErr = sqrt(ROOT::Math::Similarity(B0_v3pperp,B0_vXYe)) * B0Cand_fromMCFit->currentState().mass() / (B0_pperp.Perp2());
		  B0_lxy = B0_vdiff.Dot(B0_pperp) / B0_pperp.Mag() ;
		  ///// 3D
		  B0_pvtx3D.SetXYZ(theCosAlpha3DV.position().x(), theCosAlpha3DV.position().y(), theCosAlpha3DV.position().z()) ;	
		  B0_vdiff3D = B0_vtx3D - B0_pvtx3D ;
		  B0_cosAlpha3D =  maxCosAlpha3D ;
		  //B0_lxyz = B0_vdiff3D.Dot(B0_pperp3D)/B0_pperp3D.Mag();
		  B0_lxyz = B0_vdiff3D.Mag();
		  B0_vDiff3D[0] = B0_vdiff3D.x(); B0_vDiff3D[1] = B0_vdiff3D.y(); B0_vDiff3D[2] = B0_vdiff3D.z() ;
		  B0_lxyzErr = sqrt(ROOT::Math::Similarity(B0_vDiff3D,B0_vXYe)) / B0_vdiff3D.Mag();

		  b0CosAlphaPVCosAlpha3D->push_back( B0_cosAlpha ) ; b0CosAlpha3DPVCosAlpha3D->push_back( B0_cosAlpha3D ) ;
		  b0CTauPVCosAlpha3D->push_back( B0_ctau ) ; b0CTauPVCosAlpha3DE->push_back( B0_ctauErr ) ;
		  b0LxyPVCosAlpha3D->push_back( B0_lxy ) ; b0LxyPVCosAlpha3DE->push_back( B0_lxyErr ) ;
		  b0LxyzPVCosAlpha3D->push_back( B0_lxyz ) ; b0LxyzPVCosAlpha3DE->push_back( B0_lxyzErr ) ;


		  // Find the PV among the B0lessPV with the largest B0_cos(alpha)
		  Vertex theB0LessCosAlphaV = thePrimaryVtx ;
		  maxCosAlpha = -1. ; 
		  //
		  for (vector<TransientVertex>::iterator itv = B0_pvs.begin(), itvend = B0_pvs.end(); itv != itvend; ++itv) {
		    B0_pvtx.SetXYZ(itv->position().x(), itv->position().y(), 0) ;	
		    B0_vdiff = B0_vtx - B0_pvtx ;
		    float cosAlpha_temp = B0_vdiff.Dot(B0_pperp) / (B0_vdiff.Perp()*B0_pperp.Perp()) ; // Perp() == Mag() when z = 0
		  
		    if ( cosAlpha_temp > maxCosAlpha ) {
		      maxCosAlpha = cosAlpha_temp ;    
		      //thePrimaryVtx = Vertex(*itv);
		      //theB0LessCosAlphaV = thePrimaryVtx;
		      theB0LessCosAlphaV = Vertex(*itv) ;
		    }
		  }
		  //	
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

		  // Lifetime wrt B0LessPV with largest B0_cos(alpha) candidate 
		  B0_v2e = theB0LessCosAlphaV.error();
		  B0_vXYe = B0_v1e.matrix() + B0_v2e.matrix();
		  ///// 2D
		  B0_pvtx.SetXYZ(theB0LessCosAlphaV.position().x(), theB0LessCosAlphaV.position().y(), 0) ;	
		  B0_vdiff = B0_vtx - B0_pvtx ;
		  B0_cosAlpha =  maxCosAlpha ;
		  //B0_lxy = B0_vdiff.Dot(B0_pperp)/B0_pperp.Mag();
		  B0_lxy = B0_vdiff.Perp();
		  B0_vDiff[0] = B0_vdiff.x(); B0_vDiff[1] = B0_vdiff.y(); B0_vDiff[2] = 0 ; // needed by Similarity method
		  B0_lxyErr = sqrt(ROOT::Math::Similarity(B0_vDiff,B0_vXYe)) / B0_vdiff.Perp();
		  //
		  B0_distXY = B0_vdistXY.distance( Vertex(*B0Cand_vertex_fromMCFit), Vertex(theB0LessCosAlphaV) ) ;
		  B0_ctau = B0_distXY.value() * B0_cosAlpha * B0Cand_fromMCFit->currentState().mass() / B0_pperp.Perp();
		  //B0_ctauErr = sqrt(B0_vXYe.similarity(B0_vpperp))*B0Cand_fromMCFit->currentState().mass()/(B0_pperp.Perp2());
		  B0_ctauErr = sqrt(ROOT::Math::Similarity(B0_v3pperp,B0_vXYe)) * B0Cand_fromMCFit->currentState().mass() / (B0_pperp.Perp2());
		  B0_lxy = B0_vdiff.Dot(B0_pperp) / B0_pperp.Mag() ;
		  ///// 3D
		  B0_pvtx3D.SetXYZ(theB0LessCosAlphaV.position().x(), theB0LessCosAlphaV.position().y(), theB0LessCosAlphaV.position().z());
		  B0_vdiff3D = B0_vtx3D - B0_pvtx3D;
		  B0_cosAlpha3D = B0_vdiff3D.Dot(B0_pperp3D)/( B0_vdiff3D.Mag()*B0_pperp3D.Mag() );
		  //B0_lxyz = B0_vdiff3D.Dot(B0_pperp3D)/B0_pperp3D.Mag();
		  B0_lxyz = B0_vdiff3D.Mag();
		  B0_vDiff3D[0] = B0_vdiff3D.x(); B0_vDiff3D[1] = B0_vdiff3D.y(); B0_vDiff3D[2] = B0_vdiff3D.z() ;
		  B0_lxyzErr = sqrt(ROOT::Math::Similarity(B0_vDiff3D,B0_vXYe)) / B0_vdiff3D.Mag();
		
		  b0CosAlphaB0LessPVCosAlpha->push_back( B0_cosAlpha ) ; b0CosAlpha3DB0LessPVCosAlpha->push_back( B0_cosAlpha3D ) ;
		  b0CTauB0LessPVCosAlpha->push_back( B0_ctau ) ; b0CTauB0LessPVCosAlphaE->push_back( B0_ctauErr ) ;
		  b0LxyB0LessPVCosAlpha->push_back( B0_lxy ) ; b0LxyB0LessPVCosAlphaE->push_back( B0_lxyErr ) ;
		  b0LxyzB0LessPVCosAlpha->push_back( B0_lxyz ) ; b0LxyzB0LessPVCosAlphaE->push_back( B0_lxyzErr ) ;

		
		  // Find the PV among the B0lessPV with the largest B0_cos(alpha) 3D
		  Vertex theB0LessCosAlpha3DV = thePrimaryVtx ;
		  maxCosAlpha3D = -1. ; 
		  //
		  for (vector<TransientVertex>::iterator itv = B0_pvs.begin(), itvend = B0_pvs.end(); itv != itvend; ++itv) {
		    B0_pvtx3D.SetXYZ(itv->position().x(), itv->position().y(), itv->position().z()) ;	
		    B0_vdiff3D = B0_vtx3D - B0_pvtx3D ;
		    float cosAlpha_temp3D = B0_vdiff3D.Dot(B0_pperp3D) / (B0_vdiff3D.Mag()*B0_pperp3D.Mag()) ;
		  
		    if ( cosAlpha_temp3D > maxCosAlpha3D ) {
		      maxCosAlpha3D = cosAlpha_temp3D ;    
		      //thePrimaryVtx = Vertex(*itv);
		      //theB0LessCosAlpha3DV = thePrimaryVtx;
		      theB0LessCosAlpha3DV = Vertex(*itv) ;
		    }
		  }
		  //	
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

		  // Lifetime wrt B0LessPV with largest B0_cos(alpha) 3D candidate 
		  B0_v2e = theB0LessCosAlpha3DV.error();
		  B0_vXYe = B0_v1e.matrix() + B0_v2e.matrix();
		  ///// 2D
		  B0_pvtx.SetXYZ(theB0LessCosAlpha3DV.position().x(), theB0LessCosAlpha3DV.position().y(), 0) ;	
		  B0_vdiff = B0_vtx - B0_pvtx ;
		  B0_cosAlpha = B0_vdiff.Dot(B0_pperp)/(B0_vdiff.Perp()*B0_pperp.Perp()); ;
		  //B0_lxy = B0_vdiff.Dot(B0_pperp)/B0_pperp.Mag();
		  B0_lxy = B0_vdiff.Perp();
		  B0_vDiff[0] = B0_vdiff.x(); B0_vDiff[1] = B0_vdiff.y(); B0_vDiff[2] = 0 ; // needed by Similarity method
		  B0_lxyErr = sqrt(ROOT::Math::Similarity(B0_vDiff,B0_vXYe)) / B0_vdiff.Perp();
		  //
		  B0_distXY = B0_vdistXY.distance( Vertex(*B0Cand_vertex_fromMCFit), Vertex(theCosAlphaV) ) ;
		  B0_ctau = B0_distXY.value() * B0_cosAlpha * B0Cand_fromMCFit->currentState().mass() / B0_pperp.Perp();
		  //B0_ctauErr = sqrt(B0_vXYe.similarity(B0_vpperp))*B0Cand_fromMCFit->currentState().mass()/(B0_pperp.Perp2());
		  B0_ctauErr = sqrt(ROOT::Math::Similarity(B0_v3pperp,B0_vXYe)) * B0Cand_fromMCFit->currentState().mass() / (B0_pperp.Perp2());
		  B0_lxy = B0_vdiff.Dot(B0_pperp) / B0_pperp.Mag() ;
		  ///// 3D
		  B0_pvtx3D.SetXYZ(theB0LessCosAlpha3DV.position().x(), theB0LessCosAlpha3DV.position().y(), theB0LessCosAlpha3DV.position().z()) ;	
		  B0_vdiff3D = B0_vtx3D - B0_pvtx3D ;
		  B0_cosAlpha3D =  maxCosAlpha3D ;
		  //B0_lxyz = B0_vdiff3D.Dot(B0_pperp3D)/B0_pperp3D.Mag();
		  B0_lxyz = B0_vdiff3D.Mag();
		  B0_vDiff3D[0] = B0_vdiff3D.x(); B0_vDiff3D[1] = B0_vdiff3D.y(); B0_vDiff3D[2] = B0_vdiff3D.z() ;
		  B0_lxyzErr = sqrt(ROOT::Math::Similarity(B0_vDiff3D,B0_vXYe)) / B0_vdiff3D.Mag();
		  B0_lxy = B0_vdiff3D.Dot(B0_pperp) / B0_pperp.Mag() ;
		
		  b0CosAlphaB0LessPVCosAlpha3D->push_back( B0_cosAlpha ) ;  b0CosAlpha3DB0LessPVCosAlpha3D->push_back( B0_cosAlpha3D ) ;
		  b0CTauB0LessPVCosAlpha3D->push_back( B0_ctau ) ; b0CTauB0LessPVCosAlpha3DE->push_back( B0_ctauErr ) ;
		  b0LxyB0LessPVCosAlpha3D->push_back( B0_lxy ) ; b0LxyB0LessPVCosAlpha3DE->push_back( B0_lxyErr ) ;
		  b0LxyzB0LessPVCosAlpha3D->push_back( B0_lxyz ) ; b0LxyzB0LessPVCosAlpha3DE->push_back( B0_lxyzErr ) ;

		
		  Vertex theOtherV = thePrimaryVtx; 
		  //			
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
		  //cout<<" choose PV ="<<endl;
		  Int_t theB0CorrPV_multiplicity = -1 ;
		  if (resolveAmbiguity_) {
		    float minDz = 999999.;
		    //float minDzTrack = 999999.;
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
		  
		  //
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
		  //
		  			
		  ///////////////////////////////////////////////////////////////////////////////
						
						
		  // Lifetime wrt PV with smaller longitudinal X impact parameter 
		  B0_pvtx.SetXYZ(theOtherV.position().x(), theOtherV.position().y(), 0);
		  B0_vdiff = B0_vtx - B0_pvtx;
		  B0_cosAlpha = B0_vdiff.Dot(B0_pperp) / (B0_vdiff.Perp()*B0_pperp.Perp());
		  B0_distXY = B0_vdistXY.distance(Vertex(*B0Cand_vertex_fromMCFit), Vertex(theOtherV));
		  double B0_ctauPVX = B0_distXY.value() * B0_cosAlpha * B0Cand_fromMCFit->currentState().mass() / B0_pperp.Perp();

		  GlobalError B0_v1eX = (Vertex(*B0Cand_vertex_fromMCFit)).error();
		  GlobalError B0_v2eX = theOtherV.error();
		  AlgebraicSymMatrix33 B0_vXYeX = B0_v1eX.matrix() + B0_v2eX.matrix();
		  //double ctauErrPVX = sqrt(B0_vXYeX.similarity(B0_vpperp))*B0Cand_fromMCFit->currentState().mass()/(B0_pperp.Perp2());
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
		  B0_piIdx->push_back(std::distance(thePATTrackHandle->begin(), Track1));
		  //B0_kIdx->push_back(std::distance(thePATTrackHandle->begin(), Track2));
		  B0_kIdx->push_back(std::distance(theKaonRefittedPATTrackHandle->begin(), Track2));
		
		  nB0++;
						
		  b0Daughters.clear();
		} // 2nd loop over track (look for k)
	    } // 1st loop over track (look for pi)
	    
	  } // 2nd loop over muons (look for mu-)
	} //first loop over muons (look for mu+)
      } // if (thePATMuonHandle->size() >= 2  && hasRequestedTrigger) {
  } // if (doMC || doData) 
 	
  // AT THE END OF THE EVENT fill the tree and clear the vectors
  // ===========================================================
	
  //if (nMuMu > 0)
  if (nB0 > 0)
    Z_One_Tree_->Fill() ;
  //
  // trigger stuff
  trigRes->clear(); trigNames->clear(); L1TT->clear(); MatchTriggerNames->clear(); 
  //
  // event numbers
  runNum = 0; evtNum = 0; lumiNum = 0;
  //
  // counters
  nMu = 0; nMuMu = 0; nB0 = 0; 
  nB0_pre0 = 0; nB0_pre1 = 0; nB0_pre2 = 0; nB0_pre3 = 0; nB0_pre4 = 0; nB0_pre5 = 0; nB0_pre6 = 0; nB0_pre7 = 0; nB0_pre8 = 0; nB0_pre9 = 0; nB0_pre10 = 0; nB0_pre11 = 0; nB0_pre12 = 0; nB0_pre13 = 0; nB0_pre14 = 0;
  //
  // indices
  mu1Idx->clear(); mu2Idx->clear();
  // B0
  B0_MuMuIdx->clear(); B0_piIdx->clear(); B0_kIdx->clear();
  //
  // MC Analysis
  if (doMC) {
    // Gen Primary Vertex
    n_genEvtVtx = 0;
    genEvtVtx_X->clear(); genEvtVtx_Y->clear(); genEvtVtx_Z->clear(); 
    genEvtVtx_particles->clear();
    n_B0Ancestors->clear();
    // 
    nMCAll = 0, nMCB0 = 0; //nMCB0Vtx = 0;
    // Gen Primary Vertex
    PriVtxGen_X->clear(); PriVtxGen_Y->clear(); PriVtxGen_Z->clear(); 
    PriVtxGen_EX->clear(); PriVtxGen_EY->clear(); PriVtxGen_EZ->clear();  
    PriVtxGen_Chi2->clear(); PriVtxGen_CL->clear(); PriVtxGen_Ndof->clear();
    PriVtxGen_tracks->clear();
    //
    MCPdgIdAll->clear(); MCDanNumAll->clear();
    MCpsi2SPx->clear(); MCpsi2SPy->clear(); MCpsi2SPz->clear(); 
    MCmupPx->clear(); MCmupPy->clear(); MCmupPz->clear(); 
    MCmumPx->clear(); MCmumPy->clear(); MCmumPz->clear(); 
    MCpionPx->clear(); MCpionPy->clear(); MCpionPz->clear(); 
    MCkaonPx->clear(); MCkaonPy->clear(); MCkaonPz->clear(); 
    MCpionCh->clear(); MCkaonCh->clear();
    MCPx->clear(); MCPy->clear(); MCPz->clear();
  }
  if (Debug_) cout <<"after MC stuff clear" <<endl ;
  // Primary Vertex	
  priVtx_n = 0;
  priVtx_X = 0; priVtx_Y = 0; priVtx_Z = 0 ; 
  priVtx_XE = 0; priVtx_YE = 0; priVtx_ZE = 0 ; 
  priVtx_NormChi2 = 0; priVtx_Chi2 = 0; priVtx_CL = 0; priVtx_tracks = 0; priVtx_tracksPtSq = 0 ;
  //
  // MuMu cand
  MuMuMass->clear(); MuMuVtx_CL->clear(); MuMuVtx_Chi2->clear(); 
  MuMuPx->clear(); MuMuPy->clear(); MuMuPz->clear();
  MuMuDecayVtx_X->clear(); MuMuDecayVtx_Y->clear(); MuMuDecayVtx_Z->clear();
  MuMuDecayVtx_XE->clear(); MuMuDecayVtx_YE->clear(); MuMuDecayVtx_ZE->clear();
  MuMuMuonTrigMatch->clear();
  //
  // muons from MuMu fit
  mu1_MuMu_Px->clear(); mu1_MuMu_Py->clear(); mu1_MuMu_Pz->clear(); mu1_MuMu_Chi2->clear(); mu1_MuMu_NDF->clear();
  mu2_MuMu_Px->clear(); mu2_MuMu_Py->clear(); mu2_MuMu_Pz->clear(); mu2_MuMu_Chi2->clear(); mu2_MuMu_NDF->clear();
  MuMuType->clear();
  // Primary Vertex with "MuMu correction"
  PriVtxMuMuCorr_n->clear();
  PriVtxMuMuCorr_X->clear(); PriVtxMuMuCorr_Y->clear(); PriVtxMuMuCorr_Z->clear(); 
  PriVtxMuMuCorr_EX->clear(); PriVtxMuMuCorr_EY->clear(); PriVtxMuMuCorr_EZ->clear();  
  PriVtxMuMuCorr_Chi2->clear(); PriVtxMuMuCorr_CL->clear(); PriVtxMuMuCorr_tracks->clear();
  nTrk->clear();
  //
  // B0 cand	
  b0Mass->clear(); b0Vtx_CL->clear(); b0Vtx_Chi2->clear(); 
  b0Px->clear(); b0Py->clear(); b0Pz->clear(); 
  b0PxE->clear(); b0PyE->clear(); b0PzE->clear();
  b0DecayVtx_X->clear(); b0DecayVtx_Y->clear(); b0DecayVtx_Z->clear(); 
  b0DecayVtx_XE->clear(); b0DecayVtx_YE->clear(); b0DecayVtx_ZE->clear(); 
  // muons and tracks after B0 Cand fit
  mu1Px_MuMuPiK->clear(); mu1Py_MuMuPiK->clear(); mu1Pz_MuMuPiK->clear(); mu1E_MuMuPiK->clear();
  mu2Px_MuMuPiK->clear(); mu2Py_MuMuPiK->clear(); mu2Pz_MuMuPiK->clear(); mu2E_MuMuPiK->clear();
  piPx_MuMuPiK->clear(); piPy_MuMuPiK->clear(); piPz_MuMuPiK->clear(); piE_MuMuPiK->clear();
  pion_nsigdedx->clear(); pion_dedx->clear(); pion_dedxMass->clear(); pion_theo->clear(); pion_sigma->clear();
  pion_dedx_byHits->clear(); pion_dedxErr_byHits->clear(); pion_saturMeas_byHits->clear(); pion_Meas_byHits->clear();
  kPx_MuMuPiK->clear(); kPy_MuMuPiK->clear(); kPz_MuMuPiK->clear(); kE_MuMuPiK->clear();
  kaon_nsigdedx->clear(); kaon_dedx->clear(); kaon_dedxMass->clear(); kaon_theo->clear(); kaon_sigma->clear();
  kaon_dedx_byHits->clear(); kaon_dedxErr_byHits->clear(); kaon_saturMeas_byHits->clear(); kaon_Meas_byHits->clear();
  // Primary Vertex with "B0 correction"
  PriVtxB0Less_n->clear();
  PriVtxB0Less_X->clear(); PriVtxB0Less_Y->clear(); PriVtxB0Less_Z->clear(); 
  PriVtxB0Less_EX->clear(); PriVtxB0Less_EY->clear(); PriVtxB0Less_EZ->clear();  
  PriVtxB0Less_Chi2->clear(); PriVtxB0Less_CL->clear(); PriVtxB0Less_tracks->clear();
  // Primary Vertex with largest B0_cos(alpha)
  B0LessPV_tracksPtSq->clear(); B0LessPV_4tracksPtSq->clear();
  PriVtx_B0CosAlpha_n->clear();
  PriVtx_B0CosAlpha_X->clear(); PriVtx_B0CosAlpha_Y->clear(); PriVtx_B0CosAlpha_Z->clear(); 
  PriVtx_B0CosAlpha_EX->clear(); PriVtx_B0CosAlpha_EY->clear(); PriVtx_B0CosAlpha_EZ->clear();  
  PriVtx_B0CosAlpha_Chi2->clear(); PriVtx_B0CosAlpha_CL->clear(); PriVtx_B0CosAlpha_tracks->clear();
  PriVtxB0Less_B0CosAlpha_n->clear();
  PriVtxB0Less_B0CosAlpha_X->clear(); PriVtxB0Less_B0CosAlpha_Y->clear(); PriVtxB0Less_B0CosAlpha_Z->clear(); 
  PriVtxB0Less_B0CosAlpha_EX->clear(); PriVtxB0Less_B0CosAlpha_EY->clear(); PriVtxB0Less_B0CosAlpha_EZ->clear();  
  PriVtxB0Less_B0CosAlpha_Chi2->clear(); PriVtxB0Less_B0CosAlpha_CL->clear(); PriVtxB0Less_B0CosAlpha_tracks->clear();
  //
  // Primary Vertex with "B0 correction"
  PriVtxB0Corr_n->clear();
  PriVtxB0Corr_X->clear(); PriVtxB0Corr_Y->clear(); PriVtxB0Corr_Z->clear(); 
  PriVtxB0Corr_EX->clear(); PriVtxB0Corr_EY->clear(); PriVtxB0Corr_EZ->clear();  
  PriVtxB0Corr_Chi2->clear(); PriVtxB0Corr_CL->clear(); PriVtxB0Corr_tracks->clear();
  // B0 lifetime
  b0CosAlphaBS->clear(); b0CosAlpha3DBS->clear(); b0CTauBS->clear(); b0CTauBSE->clear(); b0LxyBS->clear(); b0LxyBSE->clear(); b0LxyzBS->clear(); b0LxyzBSE->clear(); 
  b0CosAlphaPV->clear(); b0CosAlpha3DPV->clear(); b0CTauPV->clear(); b0CTauPVE->clear(); b0LxyPV->clear(); b0LxyPVE->clear(); b0LxyzPV->clear(); b0LxyzPVE->clear();
  b0CosAlphaPVCosAlpha->clear(); b0CosAlpha3DPVCosAlpha->clear(); b0CTauPVCosAlpha->clear(); b0CTauPVCosAlphaE->clear(); b0LxyPVCosAlpha->clear(); b0LxyPVCosAlphaE->clear(); b0LxyzPVCosAlpha->clear(); b0LxyzPVCosAlphaE->clear(); 
  b0CosAlphaPVCosAlpha3D->clear(); b0CosAlpha3DPVCosAlpha3D->clear(); b0CTauPVCosAlpha3D->clear(); b0CTauPVCosAlpha3DE->clear(); b0LxyPVCosAlpha3D->clear(); b0LxyPVCosAlpha3DE->clear(); b0LxyzPVCosAlpha3D->clear(); b0LxyzPVCosAlpha3DE->clear();
  b0CosAlphaB0LessPV->clear(); b0CosAlpha3DB0LessPV->clear(); b0CTauB0LessPV->clear() ; b0CTauB0LessPVE->clear() ; b0LxyB0LessPV->clear() ; b0LxyB0LessPVE->clear() ; b0LxyzB0LessPV->clear() ; b0LxyzB0LessPVE->clear() ;
  b0CosAlphaB0LessPVCosAlpha->clear(); b0CosAlpha3DB0LessPVCosAlpha->clear(); b0CTauB0LessPVCosAlpha->clear() ; b0CTauB0LessPVCosAlphaE->clear() ; b0LxyB0LessPVCosAlpha->clear() ; b0LxyB0LessPVCosAlphaE->clear() ; b0LxyzB0LessPVCosAlpha->clear() ; b0LxyzB0LessPVCosAlphaE->clear() ;
  b0CosAlphaB0LessPVCosAlpha3D->clear(); b0CosAlpha3DB0LessPVCosAlpha3D->clear(); b0CTauB0LessPVCosAlpha3D->clear() ; b0CTauB0LessPVCosAlpha3DE->clear() ; b0LxyB0LessPVCosAlpha3D->clear() ; b0LxyB0LessPVCosAlpha3DE->clear() ; b0LxyzB0LessPVCosAlpha3D->clear() ; b0LxyzB0LessPVCosAlpha3DE->clear() ;
  b0CosAlphaPVX->clear(); b0CTauPVX->clear(); b0CTauPVXE->clear(); b0LxyPVX->clear(); b0LxyzPVX->clear(); 
  b0CTauPVX_3D->clear(); b0CTauPVX_3D_err->clear();
  //
  if (Debug_) cout <<"before muon stuff clear" <<endl ;
  // muons
  muPx->clear(); muPy->clear(); muPz->clear(); muCharge->clear();
  muD0->clear(); muDz->clear(); muChi2->clear(); muGlChi2->clear();
  mufHits->clear(); muFirstBarrel->clear(); muFirstEndCap->clear(); muD0E->clear() ;  muDzVtxErr->clear() ; muKey->clear() ;
  muIsGlobal->clear(); muIsPF->clear();
  muDzVtx->clear(); muDxyVtx->clear(); muGlMatchedStation->clear(); muGlDzVtx->clear(); muGlDxyVtx->clear(); 
  nMatchedStations->clear();
  muNDF->clear(); muGlNDF->clear(); muPhits->clear(); muShits->clear(); muGlMuHits->clear(); muType->clear(); 
  muQual->clear(); muTrack->clear(); muNOverlap->clear(); muNSharingSegWith->clear(); 
  //
  if (Debug_) cout <<"after muon stuff clear" <<endl ;
  // tracks
  trNotRef->clear(); trRef->clear();
  //
  trPx->clear(); trPy->clear(); trPz->clear(); trE->clear();	
  trNDF->clear(); trPhits->clear(); trShits->clear(); trChi2->clear();
  trD0->clear(); trD0E->clear(); trCharge->clear();
  trQualityHighPurity->clear(); trQualityTight->clear();
  trfHits->clear(); trFirstBarrel->clear(); trFirstEndCap->clear();
  trDzVtx->clear(); trDxyVtx->clear();
  tr_nsigdedx->clear(); tr_dedx->clear(); tr_dedxMass->clear(); tr_theo->clear(); tr_sigma->clear();
  tr_dedx_byHits->clear(); tr_dedxErr_byHits->clear(); tr_saturMeas_byHits->clear(); tr_Meas_byHits->clear();
  //
  if (Debug_) cout <<"end of branches clear" <<endl ;

}// analyze


// ------------ method called once each job just before starting event loop  ------------
void MuMuPiKPAT::beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup){
  //   bool changed = true;
  //   proccessName_="HLT";
  //   hltConfig_.init(iRun,iSetup,proccessName_,changed);
}


void MuMuPiKPAT::beginJob()
{
  edm::Service<TFileService> fs;
	
  Z_One_Tree_ = fs->make<TTree>("Z_data", "Z(4430) Data");
  
  Z_One_Tree_->Branch("TrigRes", &trigRes);
  Z_One_Tree_->Branch("TrigNames", &trigNames);
  Z_One_Tree_->Branch("MatchTriggerNames", &MatchTriggerNames);
  Z_One_Tree_->Branch("L1TrigRes", &L1TT);
  //
  Z_One_Tree_->Branch("evtNum", &evtNum,"evtNum/i");
  Z_One_Tree_->Branch("runNum", &runNum,"runNum/i");
  Z_One_Tree_->Branch("lumiNum", &lumiNum, "lumiNum/i");
  //
  Z_One_Tree_->Branch("priVtx_n", &priVtx_n, "priVtx_n/i");
  Z_One_Tree_->Branch("priVtx_X", &priVtx_X, "priVtx_X/f");
  Z_One_Tree_->Branch("priVtx_Y", &priVtx_Y, "priVtx_Y/f");
  Z_One_Tree_->Branch("priVtx_Z", &priVtx_Z, "priVtx_Z/f");
  Z_One_Tree_->Branch("priVtx_XE", &priVtx_XE, "priVtx_XE/f");
  Z_One_Tree_->Branch("priVtx_YE", &priVtx_YE, "priVtx_YE/f");
  Z_One_Tree_->Branch("priVtx_ZE", &priVtx_ZE, "priVtx_ZE/f");
  Z_One_Tree_->Branch("priVtx_NormChi2",&priVtx_NormChi2, "priVtx_NormChi2/f");
  Z_One_Tree_->Branch("priVtx_Chi2",&priVtx_Chi2, "priVtx_Chi2/f");
  Z_One_Tree_->Branch("priVtx_CL",&priVtx_CL, "priVtx_CL/f");
  Z_One_Tree_->Branch("priVtx_tracks", &priVtx_tracks, "priVtx_tracks/i");
  Z_One_Tree_->Branch("priVtx_tracksPtSq", &priVtx_tracksPtSq, "priVtx_tracksPtSq/f");
  // MC Analysis
  if (doMC) {
    // Gen Primary Vertex
    //Z_One_Tree_->Branch("n_genEvtVtx", &n_genEvtVtx, "n_genEvtVtx/i");
    Z_One_Tree_->Branch("genEvtVtx_X", &genEvtVtx_X); 
    Z_One_Tree_->Branch("genEvtVtx_Y", &genEvtVtx_Y);
    Z_One_Tree_->Branch("genEvtVtx_Z", &genEvtVtx_Z);
    Z_One_Tree_->Branch("genEvtVtx_particles", &genEvtVtx_particles);
    Z_One_Tree_->Branch("n_B0Ancestors", &n_B0Ancestors);
    //	 
    Z_One_Tree_->Branch("nMCAll", &nMCAll, "nMCAll/i");
    Z_One_Tree_->Branch("MCPdgIdAll", &MCPdgIdAll);
    Z_One_Tree_->Branch("MCDanNumAll", &MCDanNumAll);
    Z_One_Tree_->Branch("nMCB0",&nMCB0,"nMCB0/i");
    // Gen Primary Vertex
    //Z_One_Tree_->Branch("nMCB0Vtx",&nMCB0Vtx,"nMCB0Vtx/i");
    Z_One_Tree_->Branch("PriVtxGen_X",&PriVtxGen_X);
    Z_One_Tree_->Branch("PriVtxGen_Y",&PriVtxGen_Y);
    Z_One_Tree_->Branch("PriVtxGen_Z",&PriVtxGen_Z);
    Z_One_Tree_->Branch("PriVtxGen_EX",&PriVtxGen_EX);
    Z_One_Tree_->Branch("PriVtxGen_EY",&PriVtxGen_EY);
    Z_One_Tree_->Branch("PriVtxGen_EZ",&PriVtxGen_EZ);
    Z_One_Tree_->Branch("PriVtxGen_Chi2",&PriVtxGen_Chi2);
    Z_One_Tree_->Branch("PriVtxGen_CL",&PriVtxGen_CL);
    Z_One_Tree_->Branch("PriVtxGen_Ndof",&PriVtxGen_Ndof);
    Z_One_Tree_->Branch("PriVtxGen_tracks",&PriVtxGen_tracks);
    Z_One_Tree_->Branch("MCpsi2SPx",&MCpsi2SPx);
    Z_One_Tree_->Branch("MCpsi2SPy",&MCpsi2SPy);
    Z_One_Tree_->Branch("MCpsi2SPz",&MCpsi2SPz);
    Z_One_Tree_->Branch("MCmupPx",&MCmupPx);
    Z_One_Tree_->Branch("MCmupPy",&MCmupPy);
    Z_One_Tree_->Branch("MCmupPz",&MCmupPz);
    Z_One_Tree_->Branch("MCmumPx",&MCmumPx);
    Z_One_Tree_->Branch("MCmumPy",&MCmumPy);
    Z_One_Tree_->Branch("MCmumPz",&MCmumPz);
    Z_One_Tree_->Branch("MCpionPx",&MCpionPx);
    Z_One_Tree_->Branch("MCpionPy",&MCpionPy);
    Z_One_Tree_->Branch("MCpionPz",&MCpionPz);
    Z_One_Tree_->Branch("MCpionCh",&MCpionCh);
    Z_One_Tree_->Branch("MCkaonPx",&MCkaonPx);
    Z_One_Tree_->Branch("MCkaonPy",&MCkaonPy);
    Z_One_Tree_->Branch("MCkaonPz",&MCkaonPz);
    Z_One_Tree_->Branch("MCkaonCh",&MCkaonCh);
    Z_One_Tree_->Branch("MCPx", &MCPx);
    Z_One_Tree_->Branch("MCPy", &MCPy);
    Z_One_Tree_->Branch("MCPz", &MCPz);
  }
  // Generic muons
  Z_One_Tree_->Branch("nMu", &nMu, "nMu/i");
  Z_One_Tree_->Branch("muPx",&muPx);
  Z_One_Tree_->Branch("muPy",&muPy);
  Z_One_Tree_->Branch("muPz",&muPz);
  Z_One_Tree_->Branch("muCharge", &muCharge);
  Z_One_Tree_->Branch("muD0",&muD0);
  Z_One_Tree_->Branch("muDz",&muDz);
  Z_One_Tree_->Branch("muChi2",&muChi2);
  Z_One_Tree_->Branch("muNDF",&muNDF);
  Z_One_Tree_->Branch("muPhits",&muPhits);
  Z_One_Tree_->Branch("muShits",&muShits);

  Z_One_Tree_->Branch("muLayersTr",&muLayersTr);
  Z_One_Tree_->Branch("muLayersPix",&muLayersPix);
    
  Z_One_Tree_->Branch("muD0E",&muD0E);
  Z_One_Tree_->Branch("muDzVtxErr",&muDzVtxErr);
  Z_One_Tree_->Branch("muKey",&muKey);

  Z_One_Tree_->Branch("muIsGlobal",&muIsGlobal);
  Z_One_Tree_->Branch("muIsPF",&muIsPF);

  Z_One_Tree_->Branch("muGlMuHits",&muGlMuHits);
  Z_One_Tree_->Branch("muGlChi2",&muGlChi2);
  Z_One_Tree_->Branch("muGlNDF",&muGlNDF);
  Z_One_Tree_->Branch("muGlMatchedStation",&muGlMatchedStation);
  Z_One_Tree_->Branch("muGlDzVtx", &muGlDzVtx);
  Z_One_Tree_->Branch("muGlDxyVtx", &muGlDxyVtx);

  Z_One_Tree_->Branch("nMatchedStations", &nMatchedStations);

  Z_One_Tree_->Branch("muType",&muType);
  Z_One_Tree_->Branch("muQual",&muQual);
  Z_One_Tree_->Branch("muTrack",&muTrack);
  Z_One_Tree_->Branch("muNOverlap",&muNOverlap);
  Z_One_Tree_->Branch("muNSharingSegWith",&muNSharingSegWith);
  //
  Z_One_Tree_->Branch("mufHits", &mufHits);
  Z_One_Tree_->Branch("muFirstBarrel", &muFirstBarrel);
  Z_One_Tree_->Branch("muFirstEndCap", &muFirstEndCap);
  Z_One_Tree_->Branch("muDzVtx", &muDzVtx);
  Z_One_Tree_->Branch("muDxyVtx", &muDxyVtx);
  //
  Z_One_Tree_->Branch("trNotRef", &trNotRef);
  Z_One_Tree_->Branch("trRef", &trRef);
  //
  Z_One_Tree_->Branch("trackPx", &trPx);
  Z_One_Tree_->Branch("trackPy", &trPy);
  Z_One_Tree_->Branch("trackPz", &trPz);
  Z_One_Tree_->Branch("trackEnergy", &trE);
  Z_One_Tree_->Branch("trackNDF", &trNDF);
  Z_One_Tree_->Branch("trackPhits", &trPhits);
  Z_One_Tree_->Branch("trackShits", &trShits);
  Z_One_Tree_->Branch("trackChi2", &trChi2);
  Z_One_Tree_->Branch("trackD0", &trD0);
  Z_One_Tree_->Branch("trackD0Err", &trD0E);
  Z_One_Tree_->Branch("trackCharge", &trCharge);
  Z_One_Tree_->Branch("TrackHighPurity", &trQualityHighPurity);
  Z_One_Tree_->Branch("TrackTight", &trQualityTight);
  Z_One_Tree_->Branch("trackfHits", &trfHits);
  Z_One_Tree_->Branch("trackFirstBarrel", &trFirstBarrel);
  Z_One_Tree_->Branch("trackFirstEndCap", &trFirstEndCap);
  Z_One_Tree_->Branch("trackDzVtx", &trDzVtx);
  Z_One_Tree_->Branch("trackDxyVtx", &trDxyVtx);
  //
  Z_One_Tree_->Branch("tr_nsigdedx", &tr_nsigdedx);
  Z_One_Tree_->Branch("tr_dedx", &tr_dedx);
  Z_One_Tree_->Branch("tr_dedxMass", &tr_dedxMass);
  Z_One_Tree_->Branch("tr_theo", &tr_theo);
  Z_One_Tree_->Branch("tr_sigma", &tr_sigma);
  // by hits
  Z_One_Tree_->Branch("tr_dedx_byHits", &tr_dedx_byHits );
  Z_One_Tree_->Branch("tr_dedxErr_byHits", &tr_dedxErr_byHits );
  Z_One_Tree_->Branch("tr_saturMeas_byHits", &tr_saturMeas_byHits );
  Z_One_Tree_->Branch("tr_Meas_byHits", &tr_Meas_byHits );
  // MuMu cand
  Z_One_Tree_->Branch("nMuMu",&nMuMu,"nMuMu/i");
  Z_One_Tree_->Branch("MuMuMass",&MuMuMass);
  Z_One_Tree_->Branch("MuMuPx",&MuMuPx);
  Z_One_Tree_->Branch("MuMuPy",&MuMuPy);
  Z_One_Tree_->Branch("MuMuPz",&MuMuPz);
  Z_One_Tree_->Branch("MuMuVtx_CL",&MuMuVtx_CL);
  Z_One_Tree_->Branch("MuMuVtx_Chi2",&MuMuVtx_Chi2);

  Z_One_Tree_->Branch("MuMuDecayVtx_X",&MuMuDecayVtx_X);
  Z_One_Tree_->Branch("MuMuDecayVtx_Y",&MuMuDecayVtx_Y);
  Z_One_Tree_->Branch("MuMuDecayVtx_Z",&MuMuDecayVtx_Z);
  Z_One_Tree_->Branch("MuMuDecayVtx_XE",&MuMuDecayVtx_XE);
  Z_One_Tree_->Branch("MuMuDecayVtx_YE",&MuMuDecayVtx_YE);
  Z_One_Tree_->Branch("MuMuDecayVtx_ZE",&MuMuDecayVtx_ZE);
  // muons from psi' fit
  Z_One_Tree_->Branch("mu1Idx",&mu1Idx);
  Z_One_Tree_->Branch("mu2Idx",&mu2Idx);
  Z_One_Tree_->Branch("mu1Px_MuMu",&mu1_MuMu_Px);
  Z_One_Tree_->Branch("mu1Py_MuMu",&mu1_MuMu_Py);
  Z_One_Tree_->Branch("mu1Pz_MuMu",&mu1_MuMu_Pz);
  Z_One_Tree_->Branch("mu1Chi_MuMu2",&mu1_MuMu_Chi2);
  Z_One_Tree_->Branch("mu1NDF_MuMu",&mu1_MuMu_NDF);
  Z_One_Tree_->Branch("mu2Px_MuMu",&mu2_MuMu_Px);
  Z_One_Tree_->Branch("mu2Py_MuMu",&mu2_MuMu_Py);
  Z_One_Tree_->Branch("mu2Pz_MuMu",&mu2_MuMu_Pz);
  Z_One_Tree_->Branch("mu2Chi2_MuMu",&mu2_MuMu_Chi2);
  Z_One_Tree_->Branch("mu2NDF_MuMu",&mu2_MuMu_NDF);
  Z_One_Tree_->Branch("MuMuType",&MuMuType);
  Z_One_Tree_->Branch("MuMuMuonTrigMatch",&MuMuMuonTrigMatch);
  //
  Z_One_Tree_->Branch("PriVtxMuMuCorr_n", &PriVtxMuMuCorr_n);
  Z_One_Tree_->Branch("PriVtxMuMuCorr_X", &PriVtxMuMuCorr_X);
  Z_One_Tree_->Branch("PriVtxMuMuCorr_Y", &PriVtxMuMuCorr_Y);
  Z_One_Tree_->Branch("PriVtxMuMuCorr_Z", &PriVtxMuMuCorr_Z);
  Z_One_Tree_->Branch("PriVtxMuMuCorr_EX", &PriVtxMuMuCorr_EX);
  Z_One_Tree_->Branch("PriVtxMuMuCorr_EY", &PriVtxMuMuCorr_EY);
  Z_One_Tree_->Branch("PriVtxMuMuCorr_EZ", &PriVtxMuMuCorr_EZ);
  Z_One_Tree_->Branch("PriVtxMuMuCorr_Chi2", &PriVtxMuMuCorr_Chi2);
  Z_One_Tree_->Branch("PriVtxMuMuCorr_CL", &PriVtxMuMuCorr_CL);
  Z_One_Tree_->Branch("PriVtxMuMuCorr_tracks", &PriVtxMuMuCorr_tracks);
  Z_One_Tree_->Branch("nTrk_afterMuMu", &nTrk);
  // B0 cand
  Z_One_Tree_->Branch("nB0",&nB0,"nB0/i");
  Z_One_Tree_->Branch("nB0_pre0",&nB0_pre0,"nB0_pre0/i");
  Z_One_Tree_->Branch("nB0_pre1",&nB0_pre1,"nB0_pre1/i");
  Z_One_Tree_->Branch("nB0_pre2",&nB0_pre2,"nB0_pre2/i");
  Z_One_Tree_->Branch("nB0_pre3",&nB0_pre3,"nB0_pre3/i");
  Z_One_Tree_->Branch("nB0_pre4",&nB0_pre4,"nB0_pre4/i");
  Z_One_Tree_->Branch("nB0_pre5",&nB0_pre5,"nB0_pre5/i");
  Z_One_Tree_->Branch("nB0_pre6",&nB0_pre6,"nB0_pre6/i");
  Z_One_Tree_->Branch("nB0_pre7",&nB0_pre7,"nB0_pre7/i");
  Z_One_Tree_->Branch("nB0_pre8",&nB0_pre8,"nB0_pre8/i");
  Z_One_Tree_->Branch("nB0_pre9",&nB0_pre9,"nB0_pre9/i");
  Z_One_Tree_->Branch("nB0_pre10",&nB0_pre10,"nB0_pre10/i");
  Z_One_Tree_->Branch("nB0_pre11",&nB0_pre11,"nB0_pre11/i");
  Z_One_Tree_->Branch("nB0_pre12",&nB0_pre12,"nB0_pre12/i");
  Z_One_Tree_->Branch("nB0_pre13",&nB0_pre13,"nB0_pre13/i");
  Z_One_Tree_->Branch("nB0_pre14",&nB0_pre14,"nB0_pre14/i");
  Z_One_Tree_->Branch("B0Mass",&b0Mass);
  Z_One_Tree_->Branch("B0Px",&b0Px);
  Z_One_Tree_->Branch("B0Py",&b0Py);
  Z_One_Tree_->Branch("B0Pz",&b0Pz);
  Z_One_Tree_->Branch("B0PxE",&b0PxE);
  Z_One_Tree_->Branch("B0PyE",&b0PyE);
  Z_One_Tree_->Branch("B0PzE",&b0PzE);
  Z_One_Tree_->Branch("B0Vtx_CL",&b0Vtx_CL);
  Z_One_Tree_->Branch("B0Vtx_Chi2",&b0Vtx_Chi2);
  //
  Z_One_Tree_->Branch("B0DecayVtx_X",&b0DecayVtx_X);
  Z_One_Tree_->Branch("B0DecayVtx_Y",&b0DecayVtx_Y);
  Z_One_Tree_->Branch("B0DecayVtx_Z",&b0DecayVtx_Z);
  Z_One_Tree_->Branch("B0DecayVtx_XE",&b0DecayVtx_XE);
  Z_One_Tree_->Branch("B0DecayVtx_YE",&b0DecayVtx_YE);
  Z_One_Tree_->Branch("B0DecayVtx_ZE",&b0DecayVtx_ZE);
  //
  Z_One_Tree_->Branch("B0CosAlphaBS", &b0CosAlphaBS);
  Z_One_Tree_->Branch("B0CosAlpha3DBS", &b0CosAlpha3DBS);
  Z_One_Tree_->Branch("B0CTauBS", &b0CTauBS);
  Z_One_Tree_->Branch("B0CTauBSE", &b0CTauBSE);
  Z_One_Tree_->Branch("B0LxyBS", &b0LxyBS);
  Z_One_Tree_->Branch("B0LxyBSE", &b0LxyBSE);
  Z_One_Tree_->Branch("B0LxyzBS", &b0LxyzBS);
  Z_One_Tree_->Branch("B0LxyzBSE", &b0LxyzBSE);
  //
  Z_One_Tree_->Branch("B0CosAlphaPV", &b0CosAlphaPV);
  Z_One_Tree_->Branch("B0CosAlpha3DPV", &b0CosAlpha3DPV);
  Z_One_Tree_->Branch("B0CTauPV", &b0CTauPV);
  Z_One_Tree_->Branch("B0CTauPVE", &b0CTauPVE);
  Z_One_Tree_->Branch("B0LxyPV", &b0LxyPV);
  Z_One_Tree_->Branch("B0LxyPVE", &b0LxyPVE);
  Z_One_Tree_->Branch("B0LxyzPV", &b0LxyzPV);
  Z_One_Tree_->Branch("B0LxyzPVE", &b0LxyzPVE);
  //
  Z_One_Tree_->Branch("PriVtx_B0CosAlpha_n",&PriVtx_B0CosAlpha_n);
  Z_One_Tree_->Branch("PriVtx_B0CosAlpha_X",&PriVtx_B0CosAlpha_X);
  Z_One_Tree_->Branch("PriVtx_B0CosAlpha_Y",&PriVtx_B0CosAlpha_Y);
  Z_One_Tree_->Branch("PriVtx_B0CosAlpha_Z",&PriVtx_B0CosAlpha_Z);
  Z_One_Tree_->Branch("PriVtx_B0CosAlpha_EX",&PriVtx_B0CosAlpha_EX);
  Z_One_Tree_->Branch("PriVtx_B0CosAlpha_EY",&PriVtx_B0CosAlpha_EY);
  Z_One_Tree_->Branch("PriVtx_B0CosAlpha_EZ",&PriVtx_B0CosAlpha_EZ);
  Z_One_Tree_->Branch("PriVtx_B0CosAlpha_Chi2",&PriVtx_B0CosAlpha_Chi2);
  Z_One_Tree_->Branch("PriVtx_B0CosAlpha_CL",&PriVtx_B0CosAlpha_CL);
  Z_One_Tree_->Branch("PriVtx_B0CosAlpha_tracks",&PriVtx_B0CosAlpha_tracks);
  Z_One_Tree_->Branch("B0CosAlphaPVCosAlpha", &b0CosAlphaPVCosAlpha);
  Z_One_Tree_->Branch("B0CosAlpha3DPVCosAlpha", &b0CosAlpha3DPVCosAlpha);
  Z_One_Tree_->Branch("B0CTauPVCosAlpha", &b0CTauPVCosAlpha);
  Z_One_Tree_->Branch("B0CTauPVCosAlphaE", &b0CTauPVCosAlphaE);
  Z_One_Tree_->Branch("B0LxyPVCosAlpha", &b0LxyPVCosAlpha);
  Z_One_Tree_->Branch("B0LxyPVCosAlphaE", &b0LxyPVCosAlphaE);
  Z_One_Tree_->Branch("B0LxyzPVCosAlpha", &b0LxyzPVCosAlpha);
  Z_One_Tree_->Branch("B0LxyzPVCosAlphaE", &b0LxyzPVCosAlphaE);
  //
  Z_One_Tree_->Branch("PriVtx_B0CosAlpha3D_n",&PriVtx_B0CosAlpha3D_n);
  Z_One_Tree_->Branch("PriVtx_B0CosAlpha3D_X",&PriVtx_B0CosAlpha3D_X);
  Z_One_Tree_->Branch("PriVtx_B0CosAlpha3D_Y",&PriVtx_B0CosAlpha3D_Y);
  Z_One_Tree_->Branch("PriVtx_B0CosAlpha3D_Z",&PriVtx_B0CosAlpha3D_Z);
  Z_One_Tree_->Branch("PriVtx_B0CosAlpha3D_EX",&PriVtx_B0CosAlpha3D_EX);
  Z_One_Tree_->Branch("PriVtx_B0CosAlpha3D_EY",&PriVtx_B0CosAlpha3D_EY);
  Z_One_Tree_->Branch("PriVtx_B0CosAlpha3D_EZ",&PriVtx_B0CosAlpha3D_EZ);
  Z_One_Tree_->Branch("PriVtx_B0CosAlpha3D_Chi2",&PriVtx_B0CosAlpha3D_Chi2);
  Z_One_Tree_->Branch("PriVtx_B0CosAlpha3D_CL",&PriVtx_B0CosAlpha3D_CL);
  Z_One_Tree_->Branch("PriVtx_B0CosAlpha3D_tracks",&PriVtx_B0CosAlpha3D_tracks);
  Z_One_Tree_->Branch("B0CosAlphaPVCosAlpha3D", &b0CosAlphaPVCosAlpha3D);
  Z_One_Tree_->Branch("B0CosAlpha3DPVCosAlpha3D", &b0CosAlpha3DPVCosAlpha3D);
  Z_One_Tree_->Branch("B0CTauPVCosAlpha3D", &b0CTauPVCosAlpha3D);
  Z_One_Tree_->Branch("B0CTauPVCosAlpha3DE", &b0CTauPVCosAlpha3DE);
  Z_One_Tree_->Branch("B0LxyPVCosAlpha3D", &b0LxyPVCosAlpha3D);
  Z_One_Tree_->Branch("B0LxyPVCosAlpha3DE", &b0LxyPVCosAlpha3DE);
  Z_One_Tree_->Branch("B0LxyzPVCosAlpha3D", &b0LxyzPVCosAlpha3D);
  Z_One_Tree_->Branch("B0LxyzPVCosAlpha3DE", &b0LxyzPVCosAlpha3DE);
  //	
  Z_One_Tree_->Branch("B0LessPV_tracksPtSq",&B0LessPV_tracksPtSq);
  Z_One_Tree_->Branch("B0LessPV_4tracksPtSq",&B0LessPV_4tracksPtSq);
  Z_One_Tree_->Branch("PriVtxB0Less_n",&PriVtxB0Less_n);
  Z_One_Tree_->Branch("PriVtxB0Less_X",&PriVtxB0Less_X);
  Z_One_Tree_->Branch("PriVtxB0Less_Y",&PriVtxB0Less_Y);
  Z_One_Tree_->Branch("PriVtxB0Less_Z",&PriVtxB0Less_Z);
  Z_One_Tree_->Branch("PriVtxB0Less_EX",&PriVtxB0Less_EX);
  Z_One_Tree_->Branch("PriVtxB0Less_EY",&PriVtxB0Less_EY);
  Z_One_Tree_->Branch("PriVtxB0Less_EZ",&PriVtxB0Less_EZ);
  Z_One_Tree_->Branch("PriVtxB0Less_Chi2",&PriVtxB0Less_Chi2);
  Z_One_Tree_->Branch("PriVtxB0Less_CL",&PriVtxB0Less_CL);
  Z_One_Tree_->Branch("PriVtxB0Less_tracks",&PriVtxB0Less_tracks);
  Z_One_Tree_->Branch("B0CosAlphaB0LessPV", &b0CosAlphaB0LessPV);
  Z_One_Tree_->Branch("B0CosAlpha3DB0LessPV", &b0CosAlpha3DB0LessPV);
  Z_One_Tree_->Branch("B0CTauB0LessPV", &b0CTauB0LessPV);
  Z_One_Tree_->Branch("B0CTauB0LessPVE", &b0CTauB0LessPVE);
  Z_One_Tree_->Branch("B0LxyB0LessPV", &b0LxyB0LessPV);
  Z_One_Tree_->Branch("B0LxyB0LessPVE", &b0LxyB0LessPVE);
  Z_One_Tree_->Branch("B0LxyzB0LessPV", &b0LxyzB0LessPV);
  Z_One_Tree_->Branch("B0LxyzB0LessPVE", &b0LxyzB0LessPVE);
  //	
  Z_One_Tree_->Branch("PriVtxB0Less_B0CosAlpha_n",&PriVtxB0Less_B0CosAlpha_n);
  Z_One_Tree_->Branch("PriVtxB0Less_B0CosAlpha_X",&PriVtxB0Less_B0CosAlpha_X);
  Z_One_Tree_->Branch("PriVtxB0Less_B0CosAlpha_Y",&PriVtxB0Less_B0CosAlpha_Y);
  Z_One_Tree_->Branch("PriVtxB0Less_B0CosAlpha_Z",&PriVtxB0Less_B0CosAlpha_Z);
  Z_One_Tree_->Branch("PriVtxB0Less_B0CosAlpha_EX",&PriVtxB0Less_B0CosAlpha_EX);
  Z_One_Tree_->Branch("PriVtxB0Less_B0CosAlpha_EY",&PriVtxB0Less_B0CosAlpha_EY);
  Z_One_Tree_->Branch("PriVtxB0Less_B0CosAlpha_EZ",&PriVtxB0Less_B0CosAlpha_EZ);
  Z_One_Tree_->Branch("PriVtxB0Less_B0CosAlpha_Chi2",&PriVtxB0Less_B0CosAlpha_Chi2);
  Z_One_Tree_->Branch("PriVtxB0Less_B0CosAlpha_CL",&PriVtxB0Less_B0CosAlpha_CL);
  Z_One_Tree_->Branch("PriVtxB0Less_B0CosAlpha_tracks",&PriVtxB0Less_B0CosAlpha_tracks);
  Z_One_Tree_->Branch("B0CosAlphaB0LessPVCosAlpha", &b0CosAlphaB0LessPVCosAlpha);
  Z_One_Tree_->Branch("B0CosAlpha3DB0LessPVCosAlpha", &b0CosAlpha3DB0LessPVCosAlpha);
  Z_One_Tree_->Branch("B0CTauB0LessPVCosAlpha", &b0CTauB0LessPVCosAlpha);
  Z_One_Tree_->Branch("B0CTauB0LessPVCosAlphaE", &b0CTauB0LessPVCosAlphaE);
  Z_One_Tree_->Branch("B0LxyB0LessPVCosAlpha", &b0LxyB0LessPVCosAlpha);
  Z_One_Tree_->Branch("B0LxyB0LessPVCosAlphaE", &b0LxyB0LessPVCosAlphaE);
  Z_One_Tree_->Branch("B0LxyzB0LessPVCosAlpha", &b0LxyzB0LessPVCosAlpha);
  Z_One_Tree_->Branch("B0LxyzB0LessPVCosAlphaE", &b0LxyzB0LessPVCosAlphaE);
  //
  Z_One_Tree_->Branch("PriVtxB0Less_B0CosAlpha3D_n",&PriVtxB0Less_B0CosAlpha3D_n);
  Z_One_Tree_->Branch("PriVtxB0Less_B0CosAlpha3D_X",&PriVtxB0Less_B0CosAlpha3D_X);
  Z_One_Tree_->Branch("PriVtxB0Less_B0CosAlpha3D_Y",&PriVtxB0Less_B0CosAlpha3D_Y);
  Z_One_Tree_->Branch("PriVtxB0Less_B0CosAlpha3D_Z",&PriVtxB0Less_B0CosAlpha3D_Z);
  Z_One_Tree_->Branch("PriVtxB0Less_B0CosAlpha3D_EX",&PriVtxB0Less_B0CosAlpha3D_EX);
  Z_One_Tree_->Branch("PriVtxB0Less_B0CosAlpha3D_EY",&PriVtxB0Less_B0CosAlpha3D_EY);
  Z_One_Tree_->Branch("PriVtxB0Less_B0CosAlpha3D_EZ",&PriVtxB0Less_B0CosAlpha3D_EZ);
  Z_One_Tree_->Branch("PriVtxB0Less_B0CosAlpha3D_Chi2",&PriVtxB0Less_B0CosAlpha3D_Chi2);
  Z_One_Tree_->Branch("PriVtxB0Less_B0CosAlpha3D_CL",&PriVtxB0Less_B0CosAlpha3D_CL);
  Z_One_Tree_->Branch("PriVtxB0Less_B0CosAlpha3D_tracks",&PriVtxB0Less_B0CosAlpha3D_tracks);
  Z_One_Tree_->Branch("B0CosAlphaB0LessPVCosAlpha3D", &b0CosAlphaB0LessPVCosAlpha3D);
  Z_One_Tree_->Branch("B0CosAlpha3DB0LessPVCosAlpha3D", &b0CosAlpha3DB0LessPVCosAlpha3D);
  Z_One_Tree_->Branch("B0CTauB0LessPVCosAlpha3D", &b0CTauB0LessPVCosAlpha3D);
  Z_One_Tree_->Branch("B0CTauB0LessPVCosAlpha3DE", &b0CTauB0LessPVCosAlpha3DE);
  Z_One_Tree_->Branch("B0LxyB0LessPVCosAlpha3D", &b0LxyB0LessPVCosAlpha3D);
  Z_One_Tree_->Branch("B0LxyB0LessPVCosAlpha3DE", &b0LxyB0LessPVCosAlpha3DE);
  Z_One_Tree_->Branch("B0LxyzB0LessPVCosAlpha3D", &b0LxyzB0LessPVCosAlpha3D);
  Z_One_Tree_->Branch("B0LxyzB0LessPVCosAlpha3DE", &b0LxyzB0LessPVCosAlpha3DE);
  //
  Z_One_Tree_->Branch("PriVtxB0Corr_n",&PriVtxB0Corr_n);
  Z_One_Tree_->Branch("PriVtxB0Corr_X",&PriVtxB0Corr_X);
  Z_One_Tree_->Branch("PriVtxB0Corr_Y",&PriVtxB0Corr_Y);
  Z_One_Tree_->Branch("PriVtxB0Corr_Z",&PriVtxB0Corr_Z);
  Z_One_Tree_->Branch("PriVtxB0Corr_EX",&PriVtxB0Corr_EX);
  Z_One_Tree_->Branch("PriVtxB0Corr_EY",&PriVtxB0Corr_EY);
  Z_One_Tree_->Branch("PriVtxB0Corr_EZ",&PriVtxB0Corr_EZ);
  Z_One_Tree_->Branch("PriVtxB0Corr_Chi2",&PriVtxB0Corr_Chi2);
  Z_One_Tree_->Branch("PriVtxB0Corr_CL",&PriVtxB0Corr_CL);
  Z_One_Tree_->Branch("PriVtxB0Corr_tracks",&PriVtxB0Corr_tracks);
  //
  Z_One_Tree_->Branch("B0CosAlphaPVX", &b0CosAlphaPVX);
  Z_One_Tree_->Branch("B0CTauPVX", &b0CTauPVX);
  Z_One_Tree_->Branch("B0CTauPVXE", &b0CTauPVXE);
  Z_One_Tree_->Branch("B0LxyPVX", &b0LxyPVX);
  Z_One_Tree_->Branch("B0LxyzPVX", &b0LxyzPVX);
  Z_One_Tree_->Branch("B0CTauPVX_3D", &b0CTauPVX_3D);
  Z_One_Tree_->Branch("B0CTauPVX_3D_err", &b0CTauPVX_3D_err);
  //
  Z_One_Tree_->Branch("B0MuMuIdx", &B0_MuMuIdx);
  Z_One_Tree_->Branch("B0PionIdx", &B0_piIdx);
  Z_One_Tree_->Branch("B0KaonIdx", &B0_kIdx);
  //
  Z_One_Tree_->Branch("PiPiMass_err",& PiPiMass_err);
  //
  Z_One_Tree_->Branch("Muon1Px_MuMuPiK", &mu1Px_MuMuPiK);
  Z_One_Tree_->Branch("Muon1Py_MuMuPiK", &mu1Py_MuMuPiK);
  Z_One_Tree_->Branch("Muon1Pz_MuMuPiK", &mu1Pz_MuMuPiK);
  Z_One_Tree_->Branch("Muon1E_MuMuPiK", &mu1E_MuMuPiK);
  //
  Z_One_Tree_->Branch("Muon2Px_MuMuPiK", &mu2Px_MuMuPiK);
  Z_One_Tree_->Branch("Muon2Py_MuMuPiK", &mu2Py_MuMuPiK);
  Z_One_Tree_->Branch("Muon2Pz_MuMuPiK", &mu2Pz_MuMuPiK);
  Z_One_Tree_->Branch("Muon2E_MuMuPiK", &mu2E_MuMuPiK);
  //
  Z_One_Tree_->Branch("PionPx_MuMuPiK", &piPx_MuMuPiK);
  Z_One_Tree_->Branch("PionPy_MuMuPiK", &piPy_MuMuPiK);
  Z_One_Tree_->Branch("PionPz_MuMuPiK", &piPz_MuMuPiK);
  Z_One_Tree_->Branch("PionE_MuMuPiK", &piE_MuMuPiK);
  Z_One_Tree_->Branch("pion_nsigdedx", &pion_nsigdedx);
  Z_One_Tree_->Branch("pion_dedx", &pion_dedx);
  Z_One_Tree_->Branch("pion_dedxMass", &pion_dedxMass);
  Z_One_Tree_->Branch("pion_theo", &pion_theo);
  Z_One_Tree_->Branch("pion_sigma", &pion_sigma);
  Z_One_Tree_->Branch("pion_dedx_byHits", &pion_dedx_byHits);
  Z_One_Tree_->Branch("pion_dedxErr_byHits", &pion_dedxErr_byHits);
  Z_One_Tree_->Branch("pion_saturMeas_byHits", &pion_saturMeas_byHits);
  Z_One_Tree_->Branch("pion_Meas_byHits", &pion_Meas_byHits);
  //
  Z_One_Tree_->Branch("KaonPx_MuMuPiK", &kPx_MuMuPiK);
  Z_One_Tree_->Branch("KaonPy_MuMuPiK", &kPy_MuMuPiK);
  Z_One_Tree_->Branch("KaonPz_MuMuPiK", &kPz_MuMuPiK);
  Z_One_Tree_->Branch("KaonE_MuMuPiK", &kE_MuMuPiK);
  Z_One_Tree_->Branch("kaon_nsigdedx", &kaon_nsigdedx);
  Z_One_Tree_->Branch("kaon_dedx", &kaon_dedx);
  Z_One_Tree_->Branch("kaon_dedxMass", &kaon_dedxMass);
  Z_One_Tree_->Branch("kaon_theo", &kaon_theo);
  Z_One_Tree_->Branch("kaon_sigma", &kaon_sigma);
  Z_One_Tree_->Branch("kaon_dedx_byHits", &kaon_dedx_byHits);
  Z_One_Tree_->Branch("kaon_dedxErr_byHits", &kaon_dedxErr_byHits);
  Z_One_Tree_->Branch("kaon_saturMeas_byHits", &kaon_saturMeas_byHits);
  Z_One_Tree_->Branch("kaon_Meas_byHits", &kaon_Meas_byHits);

}// begin Job

// ------------ method called once each job just after ending the event loop  ------------
void 
MuMuPiKPAT::endJob() {
  Z_One_Tree_->GetDirectory()->cd();
  Z_One_Tree_->Write();
}//endjob


bool MuMuPiKPAT::isAbHadron(int pdgID) {
    
    if (abs(pdgID) == 511 || abs(pdgID) == 521 || abs(pdgID) == 531 || abs(pdgID) == 5122) return true;
    return false;
    
}

bool MuMuPiKPAT::isAMixedbHadron(int pdgID, int momPdgID) {
    
    if ((abs(pdgID) == 511 && abs(momPdgID) == 511 && pdgID*momPdgID < 0) || 
        (abs(pdgID) == 531 && abs(momPdgID) == 531 && pdgID*momPdgID < 0)) 
        return true;
    return false;
    
}          

std::pair<int, float> MuMuPiKPAT::findCandMCInfo(reco::GenParticleRef genCand) {
    
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
        //trueLife = vdiff.Perp()*3.09688/trueP.Perp();
        trueLife = vdiff.Perp()*genCand->mass()/trueP.Perp();
    }
    std::pair<int, float> result = std::make_pair(momJpsiID, trueLife);
    return result;
    
}

double MuMuPiKPAT::getSigmaOfLogdEdx(double logde)
{
  return 0.3;
}

float MuMuPiKPAT::getEnergyLoss(const reco::TrackRef & track)
{
  if (iexception_dedx==1) return 9999.;
  //   const DeDxDataValueMap & eloss = *energyLoss;
  const reco::DeDxDataValueMap & eloss = *energyLoss;
  return eloss[track].dEdx();
}

double MuMuPiKPAT::nsigmaofdedx(const reco::TrackRef & track, double & theo, double & sigma)
{  
	
  //   ch = (track->charge() > 0 ? 0 : 1);
    
  // no usable dE/dx if p > 2
	
  double nsigma = 99 ;
  if (iexception_dedx==1) return nsigma ;
	
  //   if(track->p() > 2) return nsigma;
	
  double m  = 0.13957; 
  double bg = track->p() / m;
	
  theo = getLogdEdx(bg);
	
  // !!!!!!
  int nhitr = track->numberOfValidHits();
	
  double meas = log(getEnergyLoss(track));
  sigma = getSigmaOfLogdEdx(theo) * pow(nhitr,-0.65);
  //   double errdedxTrk = eloss[trk1Ref].dEdxError();
  if (sigma>0)
    nsigma = (meas-theo) / sigma ;
  return nsigma;
}


double MuMuPiKPAT::getLogdEdx(double bg)
{
  const double a =  3.25 ;
  const double b =  0.288;
  const double c = -0.852;
	
  double beta = bg/sqrt(bg*bg + 1);
  double dedx = log( a/(beta*beta) + b * log(bg) + c );
	
  return dedx;
	
}


double MuMuPiKPAT::GetMass(const reco::TrackRef & track){
  double P = track->p();
  double C = 2.625;
  double K = 2.495;
  double I = getEnergyLoss(track);
  return sqrt((I-C)/K)*P;
}


template<typename T>
bool MuMuPiKPAT::isBetterMuon(const T &mu1, const T &mu2) const {
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

bool MuMuPiKPAT::isSameMuon(const reco::Muon &mu1, const reco::Muon &mu2) const {
  return (& mu1 == & mu2) ||
    //(mu1.originalObjectRef() == mu2.originalObjectRef()) ||
    (mu1.reco::Muon::innerTrack().isNonnull() ?
     mu1.reco::Muon::innerTrack() == mu2.reco::Muon::innerTrack() :
     mu1.reco::Muon::outerTrack() == mu2.reco::Muon::outerTrack());
}


//define this as a plug-in
DEFINE_FWK_MODULE(MuMuPiKPAT);

// rsync -vut --existing src/MuMuPiKPAT.cc cristella@cmssusy.ba.infn.it:/cmshome/cristella/work/Z_analysis/exclusive/clean_14ott/CMSSW_5_3_22/src/UserCode/MuMuPiKPAT/src/MuMuPiKPAT.cc

