#include "../interface/FourOnia2MuMu.h"

//Headers for the data items
#include <DataFormats/TrackReco/interface/TrackFwd.h>
#include <DataFormats/TrackReco/interface/Track.h>
#include <DataFormats/MuonReco/interface/MuonFwd.h>
#include <DataFormats/MuonReco/interface/Muon.h>
#include <DataFormats/Common/interface/View.h>
#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/VertexReco/interface/VertexFwd.h>

//Headers for services and tools
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "TMath.h"
#include "Math/VectorUtil.h"
#include "TVector3.h"
#include "../interface/FourOniaVtxReProducer.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

FourOnia2KKPAT::FourOnia2KKPAT(const edm::ParameterSet& iConfig):
muons_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
km_tracks_(consumes<edm::View<reco::TrackCollection>>(iConfig.getParameter<edm::InputTag>("kmtracks"))),
km_tracks_(consumes<edm::View<reco::TrackCollection>>(iConfig.getParameter<edm::InputTag>("kptracks"))),
thebeamspot_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotTag"))),
thePVs_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertexTag"))),
higherPuritySelection_(iConfig.getParameter<std::string>("higherPuritySelection")),
lowerPuritySelection_(iConfig.getParameter<std::string>("lowerPuritySelection")),
dimuonSelection_(iConfig.existsAs<std::string>("dimuonSelection") ? iConfig.getParameter<std::string>("dimuonSelection") : ""),
addCommonVertex_(iConfig.getParameter<bool>("addCommonVertex")),
addMuonlessPrimaryVertex_(iConfig.getParameter<bool>("addMuonlessPrimaryVertex")),
resolveAmbiguity_(iConfig.getParameter<bool>("resolvePileUpAmbiguity")),
addMCTruth_(iConfig.getParameter<bool>("addMCTruth")),
HLTFilters_(iConfig.getParameter<std::vector<std::string>>("HLTFilters"))
{
  revtxtrks_ = consumes<reco::TrackCollection>((edm::InputTag)"generalTracks"); //if that is not true, we will raise an exception
  revtxbs_ = consumes<reco::BeamSpot>((edm::InputTag)"offlineBeamSpot");
  produces<pat::CompositeCandidateCollection>();
}


FourOnia2KKPAT::~FourOnia2KKPAT()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

UInt_t FourOnia2KKPAT::isTriggerMatched(pat::CompositeCandidate *diMuon_cand) {
  const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(diMuon_cand->daughter("muon1"));
  const pat::Muon* muon2 = dynamic_cast<const pat::Muon*>(diMuon_cand->daughter("muon2"));
  UInt_t matched = 0;  // if no list is given, is not matched

  // if matched a given trigger, set the bit, in the same order as listed
  for (unsigned int iTr = 0; iTr<HLTFilters_.size(); iTr++ ) {
    // std::cout << HLTFilters_[iTr] << std::endl;
    const pat::TriggerObjectStandAloneCollection mu1HLTMatches = muon1->triggerObjectMatchesByFilter(HLTFilters_[iTr]);
    const pat::TriggerObjectStandAloneCollection mu2HLTMatches = muon2->triggerObjectMatchesByFilter(HLTFilters_[iTr]);
    if (!mu1HLTMatches.empty() && !mu2HLTMatches.empty()) matched += (1<<iTr);
  }
  // const pat::TriggerObjectStandAloneCollection muon1Collection = muon1->triggerObjectMatches();
  // const pat::TriggerObjectStandAloneCollection muon2Collection = muon2->triggerObjectMatches();
  //
  //
  // for ( size_t i = 0; i < muon1Collection.size(); ++i )
  // {
  //   std::cout << "Muon1 collection : " << muon1->triggerObjectMatch(i)->collection() << std::endl;
  //
  //   std::cout << "Filters" << std::endl;
  //   for ( size_t j = 0; j < muon1->triggerObjectMatch(i)->filterLabels().size(); ++j )
  //     std::cout << (muon1->triggerObjectMatch(i)->filterLabels())[j] << std::endl;
  //
  //   std::cout << "Paths" << std::endl;
  //   for ( size_t j = 0; j < muon1->triggerObjectMatch(i)->pathNames().size(); ++j )
  //     std::cout << (muon1->triggerObjectMatch(i)->pathNames())[j] << std::endl;
  //
  //   std::cout << "Algos" << std::endl;
  //   for ( size_t j = 0; j < muon1->triggerObjectMatch(i)->algorithmNames().size(); ++j )
  //     std::cout << (muon1->triggerObjectMatch(i)->algorithmNames())[j] << std::endl;
  //
  //   std::cout << "Conditions" << std::endl;
  //   for ( size_t j = 0; j < muon1->triggerObjectMatch(i)->conditionNames().size(); ++j )
  //     std::cout << (muon1->triggerObjectMatch(i)->conditionNames())[j] << std::endl;
  //
  //   //std::cout << (muon1->triggerObjectMatch(i)->hasL3Filter())<< std::endl;
  //
  // }
  //
  // for ( size_t i = 0; i < muon2Collection.size(); ++i )
  // {
  //   std::cout << "Muon2 collection : " << muon2->triggerObjectMatch(i)->collection() << std::endl;
  //
  //   std::cout << "Filters" << std::endl;
  //   for ( size_t j = 0; j < muon2->triggerObjectMatch(i)->filterLabels().size(); ++j )
  //     std::cout << (muon2->triggerObjectMatch(i)->filterLabels())[j] << std::endl;
  //   std::cout << "Paths" << std::endl;
  //   for ( size_t j = 0; j < muon2->triggerObjectMatch(i)->pathNames().size(); ++j )
  //     std::cout << (muon2->triggerObjectMatch(i)->pathNames())[j] << std::endl;
  // }
  //
  //
  //
  // std::cout << "Triggers matched : " << matched << std::endl;
  // std::cout << "Sizes : " << muon1Collection.size() << " - " << muon2Collection.size() << std::endl;

  return matched;
}

// ------------ method called to produce the data  ------------

void
FourOnia2KKPAT::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;
  typedef Candidate::LorentzVector LorentzVector;

  vector<double> kMasses;
  kMasses.push_back( 0.493677 );
  kMasses.push_back( 0.493677 );

  std::unique_ptr<pat::CompositeCandidateCollection> phiOutput(new pat::CompositeCandidateCollection);
  //std::cout<<"Four muonia producer"<<std::endl;
  Vertex thePrimaryV;
  Vertex theBeamSpotV;

  ESHandle<MagneticField> magneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
  const MagneticField* field = magneticField.product();

  Handle<BeamSpot> theBeamSpot;
  iEvent.getByToken(thebeamspot_,theBeamSpot);
  BeamSpot bs = *theBeamSpot;
  theBeamSpotV = Vertex(bs.position(), bs.covariance3D());

  Handle<VertexCollection> priVtxs;
  iEvent.getByToken(thePVs_, priVtxs);
  if ( priVtxs->begin() != priVtxs->end() ) {
    thePrimaryV = Vertex(*(priVtxs->begin()));
  }
  else {
    thePrimaryV = Vertex(bs.position(), bs.covariance3D());
  }

  Handle< View<reco::TrackCollection> > kaonMTracks;
  iEvent.getByToken(km_tracks_,kaonMTracks);

  Handle< View<reco::TrackCollection> > kaonPTracks;
  iEvent.getByToken(kp_tracks_,kaonPTracks);

  edm::ESHandle<TransientTrackBuilder> theTTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);
  KalmanVertexFitter vtxFitter(true);
  TrackCollection kaonLess;

  // JPsi candidates only from muons

  std::cout<<"M"<< std:endl;
  for(View<reco::TrackCollection>::const_iterator it = kaonMTracks->begin(), itend = kaonMTracks->end(); it != itend; ++it)
  {
    std::cout << it->pt() << std::endl;
  }
  std::cout<<"P"<< std:endl;
  for(View<reco::TrackCollection>::const_iterator it = kaonPTracks->begin(), itend = kaonPTracks->end(); it != itend; ++it){
  {
    std::cout << it->pt() << std::endl;
  }

      iEvent.put(std::move(phiOutput));

    }


    bool
    FourOnia2KKPAT::isAbHadron(int pdgID) {

      if (abs(pdgID) == 511 || abs(pdgID) == 521 || abs(pdgID) == 531 || abs(pdgID) == 5122) return true;
      return false;

    }

    bool
    FourOnia2KKPAT::isAMixedbHadron(int pdgID, int momPdgID) {

      if ((abs(pdgID) == 511 && abs(momPdgID) == 511 && pdgID*momPdgID < 0) ||
      (abs(pdgID) == 531 && abs(momPdgID) == 531 && pdgID*momPdgID < 0))
      return true;
      return false;

    }

    std::pair<int, float>
    FourOnia2KKPAT::findJpsiMCInfo(reco::GenParticleRef genJpsi) {

      int momJpsiID = 0;
      float trueLife = -99.;

      if (genJpsi->numberOfMothers()>0) {

        TVector3 trueVtx(0.0,0.0,0.0);
        TVector3 trueP(0.0,0.0,0.0);
        TVector3 trueVtxMom(0.0,0.0,0.0);

        trueVtx.SetXYZ(genJpsi->vertex().x(),genJpsi->vertex().y(),genJpsi->vertex().z());
        trueP.SetXYZ(genJpsi->momentum().x(),genJpsi->momentum().y(),genJpsi->momentum().z());

        bool aBhadron = false;
        reco::GenParticleRef Jpsimom = genJpsi->motherRef();       // find mothers
        if (Jpsimom.isNull()) {
          std::pair<int, float> result = std::make_pair(momJpsiID, trueLife);
          return result;
        } else {
          reco::GenParticleRef Jpsigrandmom = Jpsimom->motherRef();
          if (isAbHadron(Jpsimom->pdgId())) {
            if (Jpsigrandmom.isNonnull() && isAMixedbHadron(Jpsimom->pdgId(),Jpsigrandmom->pdgId())) {
              momJpsiID = Jpsigrandmom->pdgId();
              trueVtxMom.SetXYZ(Jpsigrandmom->vertex().x(),Jpsigrandmom->vertex().y(),Jpsigrandmom->vertex().z());
            } else {
              momJpsiID = Jpsimom->pdgId();
              trueVtxMom.SetXYZ(Jpsimom->vertex().x(),Jpsimom->vertex().y(),Jpsimom->vertex().z());
            }
            aBhadron = true;
          } else {
            if (Jpsigrandmom.isNonnull() && isAbHadron(Jpsigrandmom->pdgId())) {
              reco::GenParticleRef JpsiGrandgrandmom = Jpsigrandmom->motherRef();
              if (JpsiGrandgrandmom.isNonnull() && isAMixedbHadron(Jpsigrandmom->pdgId(),JpsiGrandgrandmom->pdgId())) {
                momJpsiID = JpsiGrandgrandmom->pdgId();
                trueVtxMom.SetXYZ(JpsiGrandgrandmom->vertex().x(),JpsiGrandgrandmom->vertex().y(),JpsiGrandgrandmom->vertex().z());
              } else {
                momJpsiID = Jpsigrandmom->pdgId();
                trueVtxMom.SetXYZ(Jpsigrandmom->vertex().x(),Jpsigrandmom->vertex().y(),Jpsigrandmom->vertex().z());
              }
              aBhadron = true;
            }
          }
          if (!aBhadron) {
            momJpsiID = Jpsimom->pdgId();
            trueVtxMom.SetXYZ(Jpsimom->vertex().x(),Jpsimom->vertex().y(),Jpsimom->vertex().z());
          }
        }

        TVector3 vdiff = trueVtx - trueVtxMom;
        //trueLife = vdiff.Perp()*3.09688/trueP.Perp();
        trueLife = vdiff.Perp()*genJpsi->mass()/trueP.Perp();
      }
      std::pair<int, float> result = std::make_pair(momJpsiID, trueLife);
      return result;

    }

    // ------------ method called once each job just before starting event loop  ------------
    void
    FourOnia2KKPAT::beginJob()
    {
    }

    // ------------ method called once each job just after ending the event loop  ------------
    void
    FourOnia2KKPAT::endJob() {
    }

    //define this as a plug-in
    DEFINE_FWK_MODULE(FourOnia2KKPAT);
