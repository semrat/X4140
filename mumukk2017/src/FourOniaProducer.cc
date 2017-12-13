#include "../interface/FourOniaProducer.h"

FourOniaProducer::FourOniaProducer(const edm::ParameterSet& ps):
phi_dimuon_Label(consumes<pat::CompositeCandidateCollection>(ps.getParameter< edm::InputTag>("phidimuons"))),
jpsi_dimuon_Label(consumes<pat::CompositeCandidateCollection>(ps.getParameter< edm::InputTag>("jpsidimuons"))),
pi0OnlineSwitch_(ps.getParameter<bool>("pi0OnlineSwitch")),
deltaMass_(ps.getParameter<std::vector<double> >("deltaMass")),
dzMax_(ps.getParameter<double>("dzmax")),
triggerMatch_(ps.getParameter<bool>("triggerMatch"))
{
  produces<pat::CompositeCandidateCollection>();
  candidates = 0;
  delta_mass_fail = 0;
  dz_phi_cut_fail = 0;
  dz_jps_cut_fail = 0;
}


void FourOniaProducer::produce(edm::Event& event, const edm::EventSetup& esetup){
  std::unique_ptr<pat::CompositeCandidateCollection> xCandColl(new pat::CompositeCandidateCollection);

  edm::Handle<pat::CompositeCandidateCollection> dimuonsPhi;
  event.getByToken(phi_dimuon_Label,dimuonsPhi);

  edm::Handle<pat::CompositeCandidateCollection> dimuonsJPsi;
  event.getByToken(jpsi_dimuon_Label,dimuonsJPsi);

  // Note: since Dimuon cand are sorted by decreasing vertex probability then the first chi cand is the one associated with the "best" dimuon
  for (pat::CompositeCandidateCollection::const_iterator  phiCand = dimuonsPhi->begin(); phiCand!= dimuonsPhi->end(); ++phiCand){

    // use only trigger-matched Jpsi or Upsilon if so requested
    if (triggerMatch_)
    if (!phiCand->userInt("isTriggerMatched"))
    continue;

    for (pat::CompositeCandidateCollection::const_iterator  jpsiCand = dimuonsJPsi->begin(); jpsiCand!= dimuonsJPsi->end(); ++jpsiCand){

      if (triggerMatch_)
      if (!jpsiCand->userInt("isTriggerMatched"))
      continue;

      if(isOverlappedMuons(*phiCand,*jpsiCand))
      continue;

      pat::CompositeCandidate xCand = makeCandidate(*phiCand, *jpsiCand);

      const reco::Vertex *ipv = phiCand->userData<reco::Vertex>("commonVertex");
      float dzPhi = fabs(Getdz(*phiCand,ipv->position()));              // onia2mumu stores vertex as userData
      xCand.addUserFloat("dzPhi",dzPhi);

      if (!cutdz(dzPhi)){
        dz_phi_cut_fail++;
        continue;
      }

      ipv = jpsiCand->userData<reco::Vertex>("commonVertex");
      float dzJpsi = fabs(Getdz(*jpsiCand,ipv->position()));              // onia2mumu stores vertex as userData
      xCand.addUserFloat("dzJpsi",dzJpsi);

      if (!cutdz(dzJpsi)){
        dz_jps_cut_fail++;
        continue;
      }


      xCandColl->push_back(xCand);
      candidates++;
    }
  }
  event.put(std::move(xCandColl));
}

float FourOniaProducer::Getdz(const pat::CompositeCandidate& c, const reco::Candidate::Point &p) {

  const reco::Candidate::LorentzVector& mom = c.p4();
  const reco::Candidate::Point& vtx = c.vertex();

  double dz = (vtx.Z()-p.Z()) - ((vtx.X()-p.X())*mom.X()+(vtx.Y()-p.Y())*mom.Y())/mom.Rho() * mom.Z()/mom.Rho();
  return (float) dz;

}

void FourOniaProducer::endJob(){
  std::cout << "###########################" << std::endl;
  std::cout << "Chi Candidate producer report:" << std::endl;
  std::cout << "###########################" << std::endl;
  // std::cout << "Delta mass fail: " << delta_mass_fail << std::endl;
  std::cout << "Dz phi fail:      " << dz_phi_cut_fail << std::endl;
  std::cout << "Dz jps fail:      " << dz_jps_cut_fail << std::endl;
  std::cout << "###########################" << std::endl;
  std::cout << "Found " << candidates << " Chi candidates." << std::endl;
  std::cout << "###########################" << std::endl;
}
//
// UInt_t FourOniaProducer::isTriggerMatched(const pat::CompositeCandidate *diMuon_cand) {
//   const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(diMuon_cand->daughter("muon1"));
//   const pat::Muon* muon2 = dynamic_cast<const pat::Muon*>(diMuon_cand->daughter("muon2"));
//   UInt_t matched = 0;  // if no list is given, is not matched
//
//   // if matched a given trigger, set the bit, in the same order as listed
//   for (unsigned int iTr = 0; iTr<HLTFilters_.size(); iTr++ ) {
//     const pat::TriggerObjectStandAloneCollection mu1HLTMatches = muon1->triggerObjectMatchesByFilter(HLTFilters_[iTr]);
//     const pat::TriggerObjectStandAloneCollection mu2HLTMatches = muon2->triggerObjectMatchesByFilter(HLTFilters_[iTr]);
//     if (!mu1HLTMatches.empty() && !mu2HLTMatches.empty()) matched += (1<<iTr);
//   }
//   return matched;
// }

bool FourOniaProducer::isOverlappedMuons(const pat::CompositeCandidate *phi,const pat::CompositeCandidate *jpsi) {

  bool same = false;

  std::vector<pat::Muon*> fourMuons;

  fourMuons.push_back(dynamic_cast<const pat::Muon*>(phi->daughter("muon1")));
  fourMuons.push_back(dynamic_cast<const pat::Muon*>(phi->daughter("muon2")));

  fourMuons.push_back(dynamic_cast<const pat::Muon*>(jpsi->daughter("muon1")));
  fourMuons.push_back(dynamic_cast<const pat::Muon*>(jpsi->daughter("muon2")));

  for (size_t i = 0; i < fourMuons.size(); i++)
  for (size_t j = i+1; j < fourMuons.size(); j++)
  same = same || (fourMuons[i].track() == fourMuons[j].track());

  // const pat::Muon* muon2 = dynamic_cast<const pat::Muon*>(phi->daughter("muon2"));
  // const pat::Muon* muon3 = dynamic_cast<const pat::Muon*>(jpsi->daughter("muon1"));
  // const pat::Muon* muon4 = dynamic_cast<const pat::Muon*>(jpsi->daughter("muon2"));

  return same;

}

const pat::CompositeCandidate FourOniaProducer::makeCandidate(const pat::CompositeCandidate& phi,
  const pat::CompositeCandidate& jpsi){
    pat::CompositeCandidate xCand;
    xCand.addDaughter(phi,"phi");
    xCand.addDaughter(jpsi,"jpsi");
    reco::Candidate::LorentzVector vX = phi.p4() + jpsi.p4();
    xCand.setP4(vX);
    return xCand;
  }

  // check if the mass difference is in desired range
  bool FourOniaProducer::cutDeltaMass(const pat::CompositeCandidate& xCand,
    const pat::CompositeCandidate& dimuonCand){
      float deltam = xCand.p4().M() - dimuonCand.p4().M();
      float m1     = deltaMass_[0];
      float m2     = deltaMass_[1];
      return (deltam > m1 && deltam < m2);
    }

    //define this as a plug-in
    DEFINE_FWK_MODULE(FourOniaProducer);
