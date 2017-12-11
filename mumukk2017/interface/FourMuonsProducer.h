#ifndef HeavyFlavorAnalysis_Onia2MuMuMuMu_FourMuonsProducer_h
#define HeavyFlavorAnalysis_Onia2MuMuMuMu_FourMuonsProducer_h


// system include files
#include <memory>

// FW include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/Utils/interface/PtComparator.h"

// DataFormat includes
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/PatCandidates/interface/Muon.h>

#include <CommonTools/UtilAlgos/interface/StringCutObjectSelector.h>
#include "RecoVertex/VertexTools/interface/InvariantMassFromVertex.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"

template<typename T>
struct GreaterByVProb {
  bool operator()( const T & t1, const T & t2 ) const {
    return t1.userFloat("vProb") > t2.userFloat("vProb");
  }
};


//
// class decleration
//

class oniaMuMuMuMuPAT : public edm::EDProducer {
 public:
  explicit oniaMuMuMuMuPAT(const edm::ParameterSet&);
  ~oniaMuMuMuMuPAT() override;

 private:
  void beginJob() override ;
  void produce(edm::Event&, const edm::EventSetup&) override;
  void endJob() override ;
  bool isAbHadron(int pdgID);
  bool isAMixedbHadron(int pdgID, int momPdgID);
  std::pair<int, float> findJpsiMCInfo(reco::GenParticleRef genJpsi);
  bool uniqueMuons(const std::vector<pat::Muon> fourMuons) const;

  // ----------member data ---------------------------
 private:

  edm::EDGetTokenT<edm::View<pat::Muon>> muons_;
  edm::EDGetTokenT<reco::BeamSpot> thebeamspot_;
  edm::EDGetTokenT<reco::VertexCollection> thePVs_;
  edm::EDGetTokenT<reco::TrackCollection> revtxtrks_;
  edm::EDGetTokenT<reco::BeamSpot> revtxbs_;
  StringCutObjectSelector<pat::Muon> higherPuritySelection_;
  StringCutObjectSelector<pat::Muon> lowerPuritySelection_;
  StringCutObjectSelector<reco::Candidate, true> quadmuonSelection_;
  bool addCommonVertex_, addMuonlessPrimaryVertex_;
  bool resolveAmbiguity_;
  bool addMCTruth_;
  GreaterByVProb<pat::CompositeCandidate> vPComparator_;

  InvariantMassFromVertex massCalculator;

};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//

#endif
