//#include "HeavyFlavorAnalysis/Onia2MuMu/interface/oniaMuMuMuMuPAT.h"

#include "../interface/FourMuonsProducer.h"
#include <iostream>
#include <string>

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
#include "HeavyFlavorAnalysis/Onia2MuMu/interface/OniaVtxReProducer.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

oniaMuMuMuMuPAT::oniaMuMuMuMuPAT(const edm::ParameterSet& iConfig):
muons_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
thebeamspot_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotTag"))),
thePVs_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertexTag"))),
higherPuritySelection_(iConfig.getParameter<std::string>("higherPuritySelection")),
lowerPuritySelection_(iConfig.getParameter<std::string>("lowerPuritySelection")),
quadmuonSelection_(iConfig.existsAs<std::string>("quadmuonSelection") ? iConfig.getParameter<std::string>("quadmuonSelection") : ""),
addCommonVertex_(iConfig.getParameter<bool>("addCommonVertex")),
addMuonlessPrimaryVertex_(iConfig.getParameter<bool>("addMuonlessPrimaryVertex")),
resolveAmbiguity_(iConfig.getParameter<bool>("resolvePileUpAmbiguity")),
addMCTruth_(iConfig.getParameter<bool>("addMCTruth"))
{
  revtxtrks_ = consumes<reco::TrackCollection>((edm::InputTag)"generalTracks"); //if that is not true, we will raise an exception
  revtxbs_ = consumes<reco::BeamSpot>((edm::InputTag)"offlineBeamSpot");
  produces<pat::CompositeCandidateCollection>();
}


oniaMuMuMuMuPAT::~oniaMuMuMuMuPAT()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
oniaMuMuMuMuPAT::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;
  typedef Candidate::LorentzVector LorentzVector;

  vector<double> muMasses;
  muMasses.push_back( 0.1056583715);
  muMasses.push_back( 0.1056583715);
  muMasses.push_back( 0.1056583715);
  muMasses.push_back( 0.1056583715);

  std::unique_ptr<pat::CompositeCandidateCollection> oniaOutput(new pat::CompositeCandidateCollection);

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

  Handle< View<pat::Muon> > muons;
  iEvent.getByToken(muons_,muons);

  edm::ESHandle<TransientTrackBuilder> theTTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);
  KalmanVertexFitter vtxFitter(true);
  TrackCollection muonLess;
  std::string bufferstring;
  // xcand candidates only from muons
  for(View<pat::Muon>::const_iterator it = muons->begin(), itend = muons->end(); it != itend; ++it){
    // both must pass low quality
    if(!lowerPuritySelection_(*it)) continue;
    bufferstring = "Good 1 muon!";
    // std::cout<<bufferstring<<std::endl;
    for(View<pat::Muon>::const_iterator it2 = muons->begin(), itend = muons->end(); it2 != itend; ++it2){
      // both must pass low quality
      if(!lowerPuritySelection_(*it2)) continue;
      for(View<pat::Muon>::const_iterator it3 = muons->begin(), itend = muons->end(); it3 != itend; ++it3){
        if(!lowerPuritySelection_(*it3)) continue;
        if((it3->charge() + it2->charge() + it->charge())>1)
          continue;

        if((it3->charge() + it2->charge() + it->charge())<-1)
          continue;
        for(View<pat::Muon>::const_iterator it4 = muons->begin(), itend = muons->end(); it4 != itend; ++it4){
          if((it3->charge() + it2->charge() + it4->charge() + it->charge())!=0)
            continue;
          if(!lowerPuritySelection_(*it4)) continue;
          // one must pass tight quality
          if (!(higherPuritySelection_(*it) || higherPuritySelection_(*it2) || higherPuritySelection_(*it3) || higherPuritySelection_(*it4))) continue;
          bufferstring = "Good 4 muons!";
          std::cout<<bufferstring<<std::endl;

          pat::CompositeCandidate myCand;
          vector<TransientVertex> pvs;
          std::vector<pat::Muon> fourMuons;
          fourMuons.push_back(*it);
          fourMuons.push_back(*it2);
          fourMuons.push_back(*it3);
          fourMuons.push_back(*it4);


          if (oniaMuMuMuMuPAT::uniqueMuons(fourMuons))
            continue;

          // ---- no explicit order defined ----
          myCand.addDaughter(*it, "muon1");
          myCand.addDaughter(*it2,"muon2");
          myCand.addDaughter(*it3, "muon3");
          myCand.addDaughter(*it4,"muon4");

          // ---- define and set candidate's 4momentum  ----
          LorentzVector xcand = it->p4() + it2->p4() + it3->p4() + it4->p4();
          myCand.setP4(xcand);
          myCand.setCharge(it->charge() + it2->charge() + it3->charge() + it4->charge());

          if(myCand.charge()!=0)
            continue;

          if(!quadmuonSelection_(myCand)) continue;

          oniaOutput->push_back(myCand);

          }

        }

      }
    }

      std::sort(oniaOutput->begin(),oniaOutput->end(),vPComparator_);

      iEvent.put(std::move(oniaOutput));

    }


    bool
    oniaMuMuMuMuPAT::isAbHadron(int pdgID) {

      if (abs(pdgID) == 511 || abs(pdgID) == 521 || abs(pdgID) == 531 || abs(pdgID) == 5122) return true;
      return false;

    }

    bool
    oniaMuMuMuMuPAT::isAMixedbHadron(int pdgID, int momPdgID) {

      if ((abs(pdgID) == 511 && abs(momPdgID) == 511 && pdgID*momPdgID < 0) ||
      (abs(pdgID) == 531 && abs(momPdgID) == 531 && pdgID*momPdgID < 0))
      return true;
      return false;

    }

    bool oniaMuMuMuMuPAT::uniqueMuons(const std::vector<pat::Muon> fourMuons) const{

      bool same = false;
      // for (size_t i = 0; i < fourMuons.size(); i++)
      // for (size_t j = i+1; j < fourMuons.size(); j++)
      // same = same || (fourMuons[i] == fourMuons[j])

      for (size_t i = 0; i < fourMuons.size(); i++)
      for (size_t j = i+1; j < fourMuons.size(); j++)
      same = same || (fourMuons[i].track() == fourMuons[j].track());

      return same;
    }

    std::pair<int, float>
    oniaMuMuMuMuPAT::findJpsiMCInfo(reco::GenParticleRef genJpsi) {

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
    oniaMuMuMuMuPAT::beginJob()
    {
    }

    // ------------ method called once each job just after ending the event loop  ------------
    void
    oniaMuMuMuMuPAT::endJob() {
    }

    //define this as a plug-in
    DEFINE_FWK_MODULE(oniaMuMuMuMuPAT);
