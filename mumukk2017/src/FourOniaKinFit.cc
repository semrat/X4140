// -*- C++ -*-
//
// Package:    FourOniaKinFit
// Class:      FourOniaKinFit
//
/**

Description: Kinematic Fit for 4Mu-Onia

Implementation:
Original work from Stefano Argiro, and the Torino group
Adapter for MINIAOD by Alberto Sanchez-Hernandez
Adapter for JPsiPhi->4Mu search by Adriano Di Florio
*/


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

///For kinematic fit:
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/PatCandidates/interface/PackedCandidate.h>

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"

#include <boost/foreach.hpp>
#include <string>

//
// class declaration
//

class FourOniaKinFit : public edm::EDProducer {
public:
  explicit FourOniaKinFit(const edm::ParameterSet&);
  ~FourOniaKinFit() override {};
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  double x_mass_;
  std::string product_name_;
  int pdgID_;
  edm::EDGetTokenT<pat::CompositeCandidateCollection> x_Label;

  template<typename T>
  struct GreaterByVProb {
    typedef T first_argument_type;
    typedef T second_argument_type;
    bool operator()( const T & t1, const T & t2 ) const {
      return t1.userFloat("vProb") > t2.userFloat("vProb");
    }
  };
};

FourOniaKinFit::FourOniaKinFit(const edm::ParameterSet& iConfig) {
  x_Label     = consumes<pat::CompositeCandidateCollection>(iConfig.getParameter< edm::InputTag>("x_cand"));
  x_mass_ = iConfig.getParameter<double>("x_mass");
  pdgID_ = iConfig.getParameter<int>("pdgID");
  product_name_ = iConfig.getParameter<std::string>("product_name");
  produces<pat::CompositeCandidateCollection>(product_name_);
}

// ------------ method called to produce the data  ------------
void FourOniaKinFit::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // Grab paramenters
  edm::Handle<pat::CompositeCandidateCollection> xCandHandle;
  iEvent.getByToken(x_Label, xCandHandle);

  //Kinemati refit collection
  std::unique_ptr<pat::CompositeCandidateCollection> xCompCandRefitColl(new pat::CompositeCandidateCollection);

  // Kinematic fit
  edm::ESHandle<TransientTrackBuilder> theB;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

  int index=-1;

  for (pat::CompositeCandidateCollection::const_iterator xCand=xCandHandle->begin();xCand!=xCandHandle->end();++xCand) {

    index++;

    reco::TrackRef JpsiTk[2]={
      ( dynamic_cast<const pat::Muon*>(xCand->daughter("phi")->daughter("muon1") ) )->innerTrack(),
      ( dynamic_cast<const pat::Muon*>(xCand->daughter("phi")->daughter("muon2") ) )->innerTrack(),
      ( dynamic_cast<const pat::Muon*>(xCand->daughter("jpsi")->daughter("muon3") ) )->innerTrack(),
      ( dynamic_cast<const pat::Muon*>(xCand->daughter("jpsi")->daughter("muon4") ) )->innerTrack()
    };


    std::vector<reco::TransientTrack> PhiMuMuTT;
    PhiMuMuTT.push_back((*theB).build(&JpsiTk[0]));
    PhiMuMuTT.push_back((*theB).build(&JpsiTk[1]));

    std::vector<reco::TransientTrack> JPsiMuMuTT;
    JPsiMuMuTT.push_back((*theB).build(&JpsiTk[0]));
    JPsiMuMuTT.push_back((*theB).build(&JpsiTk[1]));

    const ParticleMass muonMass(0.1056584);
    float muonSigma = muonMass*1E-6;

    std::vector<RefCountedKinematicParticle> allXDaughters;
    allXDaughters.push_back(pFactory.particle (PhiMuMuTT[0], muonMass, float(0), float(0), muonSigma));
    allXDaughters.push_back(pFactory.particle (PhiMuMuTT[1], muonMass, float(0), float(0), muonSigma));
    allXDaughters.push_back(pFactory.particle (JPsiMuMuTT[0], muonMass, float(0), float(0), muonSigma));
    allXDaughters.push_back(pFactory.particle (JPsiMuMuTT[1], muonMass, float(0), float(0), muonSigma));

    KinematicConstrainedVertexFitter constVertexFitter;

    MultiTrackKinematicConstraint *xcand_mtc = new  TwoTrackMassKinematicConstraint(x_mass_);
    RefCountedKinematicTree XTree = constVertexFitter.fit(allXDaughters,xcand_mtc);

    if (!XTree->isEmpty()) {

      XTree->movePointerToTheTop();
      RefCountedKinematicParticle fitX = XTree->currentParticle();
      RefCountedKinematicVertex XDecayVertex = XTree->currentDecayVertex();

      if (fitX->currentState().isValid()) { //Get get x

        float XM_fit  = fitX->currentState().mass();
        float XPx_fit = fitX->currentState().kinematicParameters().momentum().x();
        float XPy_fit = fitX->currentState().kinematicParameters().momentum().y();
        float XPz_fit = fitX->currentState().kinematicParameters().momentum().z();
        float XVtxX_fit = XDecayVertex->position().x();
        float XVtxY_fit = XDecayVertex->position().y();
        float XVtxZ_fit = XDecayVertex->position().z();
        float XVtxP_fit = ChiSquaredProbability((double)(XDecayVertex->chiSquared()),
        (double)(XDecayVertex->degreesOfFreedom()));

        reco::CompositeCandidate recoX(0, math::XYZTLorentzVector(XPx_fit, XPy_fit, XPz_fit,
        sqrt(XM_fit*XM_fit + XPx_fit*XPx_fit + XPy_fit*XPy_fit + XPz_fit*XPz_fit)), math::XYZPoint(XVtxX_fit,
        XVtxY_fit, XVtxZ_fit), pdgID_);

        pat::CompositeCandidate patX(recoX);
        patX.addUserFloat("vProb",XVtxP_fit);
        patX.addUserInt("Index",index);  // this also holds the index of the current xCand

        //get first muon
        bool child = XTree->movePointerToTheFirstChild();
        RefCountedKinematicParticle fitMu1 = XTree->currentParticle();
        if (!child) break;

        float mu1M_fit  = fitMu1->currentState().mass();
        float mu1Q_fit  = fitMu1->currentState().particleCharge();
        float mu1Px_fit = fitMu1->currentState().kinematicParameters().momentum().x();
        float mu1Py_fit = fitMu1->currentState().kinematicParameters().momentum().y();
        float mu1Pz_fit = fitMu1->currentState().kinematicParameters().momentum().z();

        reco::CompositeCandidate recoMu1(mu1Q_fit, math::XYZTLorentzVector(mu1Px_fit, mu1Py_fit, mu1Pz_fit,
        sqrt(mu1M_fit*mu1M_fit + mu1Px_fit*mu1Px_fit + mu1Py_fit*mu1Py_fit +
        mu1Pz_fit*mu1Pz_fit)), math::XYZPoint(XVtxX_fit, XVtxY_fit, XVtxZ_fit), 13);

        pat::CompositeCandidate patMu1(recoMu1);

        //get second muon
        child = XTree->movePointerToTheNextChild();
        RefCountedKinematicParticle fitMu2 = XTree->currentParticle();

        if (!child) break;

        float mu2M_fit  = fitMu2->currentState().mass();
        float mu2Q_fit  = fitMu2->currentState().particleCharge();
        float mu2Px_fit = fitMu2->currentState().kinematicParameters().momentum().x();
        float mu2Py_fit = fitMu2->currentState().kinematicParameters().momentum().y();
        float mu2Pz_fit = fitMu2->currentState().kinematicParameters().momentum().z();
        reco::CompositeCandidate recoMu2(mu2Q_fit, math::XYZTLorentzVector(mu2Px_fit, mu2Py_fit, mu2Pz_fit,
        sqrt(mu2M_fit*mu2M_fit + mu2Px_fit*mu2Px_fit + mu2Py_fit*mu2Py_fit +
        mu2Pz_fit*mu2Pz_fit)), math::XYZPoint(XVtxX_fit, XVtxY_fit, XVtxZ_fit), 13);

        pat::CompositeCandidate patMu2(recoMu2);

        //get third muon
        child = XTree->movePointerToTheFirstChild();
        RefCountedKinematicParticle fitMu3 = XTree->currentParticle();
        if (!child) break;

        float mu3M_fit  = fitMu3->currentState().mass();
        float mu3Q_fit  = fitMu3->currentState().particleCharge();
        float mu3Px_fit = fitMu3->currentState().kinematicParameters().momentum().x();
        float mu3Py_fit = fitMu3->currentState().kinematicParameters().momentum().y();
        float mu3Pz_fit = fitMu3->currentState().kinematicParameters().momentum().z();

        reco::CompositeCandidate recoMu3(mu3Q_fit, math::XYZTLorentzVector(mu3Px_fit, mu3Py_fit, mu3Pz_fit,
        sqrt(mu3M_fit*mu3M_fit + mu3Px_fit*mu3Px_fit + mu3Py_fit*mu3Py_fit +
        mu3Pz_fit*mu3Pz_fit)), math::XYZPoint(XVtxX_fit, XVtxY_fit, XVtxZ_fit), 13);

        pat::CompositeCandidate patMu3(recoMu3);


        //get fourth muon
        child = XTree->movePointerToTheFirstChild();
        RefCountedKinematicParticle fitMu4 = XTree->currentParticle();
        if (!child) break;

        float mu4M_fit  = fitMu4->currentState().mass();
        float mu4Q_fit  = fitMu4->currentState().particleCharge();
        float mu4Px_fit = fitMu4->currentState().kinematicParameters().momentum().x();
        float mu4Py_fit = fitMu4->currentState().kinematicParameters().momentum().y();
        float mu4Pz_fit = fitMu4->currentState().kinematicParameters().momentum().z();

        reco::CompositeCandidate recoMu4(mu4Q_fit, math::XYZTLorentzVector(mu4Px_fit, mu4Py_fit, mu4Pz_fit,
        sqrt(mu4M_fit*mu4M_fit + mu4Px_fit*mu4Px_fit + mu4Py_fit*mu4Py_fit +
        mu4Pz_fit*mu4Pz_fit)), math::XYZPoint(XVtxX_fit, XVtxY_fit, XVtxZ_fit), 13);

        pat::CompositeCandidate patMu4(recoMu4);

        //Define Onia from two muons
        pat::CompositeCandidate phi;
        phi.addDaughter(patMu1,"muon1");
        phi.addDaughter(patMu2,"muon2");
        phi.setP4(patMu1.p4()+patMu2.p4());

        pat::CompositeCandidate jps;
        jps.addDaughter(patMu3,"muon1");
        jps.addDaughter(patMu4,"muon2");
        jps.setP4(patMu3.p4()+patMu4.p4());

        patX.addDaughter(phi,"phi");
        patX.addDaughter(jps,"jpsi");

        xCompCandRefitColl->push_back(patX);

      }
    }
  }
                  // End kinematic fit

                  // ...ash not sorted, since we will use the best un-refitted candidate
                  // now sort by vProb

        FourOniaKinFit::GreaterByVProb<pat::CompositeCandidate> vPComparator;
        std::sort(xCompCandRefitColl->begin(),xCompCandRefitColl->end(), vPComparator);

        iEvent.put(std::move(xCompCandRefitColl),product_name_);

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void FourOniaKinFit::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(FourOniaKinFit);
