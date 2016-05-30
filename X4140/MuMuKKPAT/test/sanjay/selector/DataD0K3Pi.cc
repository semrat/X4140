// -*- C++ -*-
//
// Package:    DataD0K3Pi
// Class:      DataD0K3Pi
// 
/**\class DataD0K3Pi DataD0K3Pi.cc UserCode/DataD0K3Pi/src/DataD0K3Pi.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author: A.P.- INFN Bari 
//         Created:  Mon May  9 12:29:46 CEST 2011
// $Id$
//
//


//--- system include files
#include <memory>
#include <string>
#include <vector>
#include <cmath>

//--- basic C++ 
#include <iostream>


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
//#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"
//
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Provenance/interface/ParameterSetID.h"
//
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
//
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h" //for the get<TrackRef>() Call
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
//#include "DataFormats/RecoCandidate/interface/FitQuality.h"
#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexSorter.h"
//#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
//
//magnetic Field
//#include "MagneticField/Engine/interface/MagneticField.h"
//#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
//#include "MagneticField/VolumeBasedEngine/interface/VolumeBasedMagneticField.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "DetectorDescription/Core/interface/DDLogicalPart.h"
#include "DetectorDescription/Core/interface/DDSolid.h"
//
//#include "PhysicsTools/CandUtils/interface/helicityAngle.h"
#include "DataFormats/Candidate/interface/OverlapChecker.h"
//
#include "CommonTools/CandUtils/interface/AddFourMomenta.h" 
//
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
//
class TFile;
class TH1D;
class TH2D;
class TH1I;
class TH1F;
#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>
#include <typeinfo>
//
#include "TMath.h"
#include <Math/VectorUtil.h>
//
#include "FWCore/MessageLogger/interface/MessageLogger.h"
//
#include<algorithm>
//
using namespace edm;
using namespace reco;
using namespace std;
//
using std::vector;
using std::cout;
using std::endl;
//
using reco::Track;
using reco::TransientTrack; 
using reco::TrackRef;
//
//
// class declaration
//
class DataD0K3Pi : public edm::EDAnalyzer {
   public:
      explicit DataD0K3Pi(const edm::ParameterSet&);
      ~DataD0K3Pi();
   //
   private:
      virtual void beginJob(void) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      std::string outputFileName_;
      //
      edm::InputTag muons_;
      edm::InputTag tracks_;
      edm::InputTag pions_;
      edm::InputTag softpions_;
      edm::InputTag kaons_;
      //edm::InputTag k3pi_;
      //
      // output TTree
      TTree* Tree_;
      //
   // output branches
   //
   // general variables
   //
   int nK3Pi;
   //int isBSValid;
   int runId;
   int evtId;
   int lsId;
   double PVx;
   double PVy;
   double PVz;
   double PVxErr;
   double PVyErr;
   double PVzErr;
   int ntracks;
   int nGlbMuons;
   int npions;
   int nsoftpions;
   int nkaons;
   //bool hasBq; //-MC-only
   //bool hasGenD0; //-MC-only
   //
   //////////////////////////// D0 candidate variables
   //
   double K3Pi_Mass[10000]; 
   double K3Pi_pseudoP[10000]; 
   double K3Pi_pseudoPt[10000]; 
   double K3Pi_pseudoEta[10000]; 
   double K3Pi_pseudoPhi[10000]; 
   int K3Pi_pseudoCharge[10000]; 
   //
   double K3Pi_deltaR_kp1[10000];
   double K3Pi_deltaR_kp2[10000];
   double K3Pi_deltaR_kp3[10000];
   double K3Pi_deltaR_p1p3[10000];
   double K3Pi_deltaR_p2p3[10000];
   double K3Pi_deltaR_p1p2[10000];
   //
   double K3Pi_pointingAngle[10000];
   double K3Pi_2DpointingAngle[10000];
   //
   double K3Pi_VtxChi2[10000]; 
   double K3Pi_VtxChi2Prob[10000]; 
   double K3Pi_VtxCalcDistXY[10000]; 
   double K3Pi_VtxDistXY[10000]; 
   double K3Pi_VtxDistXYError[10000]; 
   double K3Pi_VtxDistXYSignif[10000];
   //
   //////////////////////////D0 daughters variables (kaon-pion)
   //
   double KDau_P[10000];
   double P1Dau_P[10000];
   double P2Dau_P[10000];
   double P3Dau_P[10000];
   //
   double KDau_Pt[10000];
   double P1Dau_Pt[10000];
   double P2Dau_Pt[10000];
   double P3Dau_Pt[10000];
   //
   double KDau_Eta[10000];
   double P1Dau_Eta[10000];
   double P2Dau_Eta[10000];
   double P3Dau_Eta[10000];
   //
   double KDau_Phi[10000];
   double P1Dau_Phi[10000];
   double P2Dau_Phi[10000];
   double P3Dau_Phi[10000];
   //
   int KDau_Charge[10000];
   int P1Dau_Charge[10000];
   int P2Dau_Charge[10000];
   int P3Dau_Charge[10000];
   //
   double KDau_NormChi2[10000];
   double P1Dau_NormChi2[10000];
   double P2Dau_NormChi2[10000];
   double P3Dau_NormChi2[10000];
   //
   double KDau_Chi2[10000];
   double P1Dau_Chi2[10000];
   double P2Dau_Chi2[10000];
   double P3Dau_Chi2[10000];
   //
   double KDau_Prob[10000];
   double P1Dau_Prob[10000];
   double P2Dau_Prob[10000];
   double P3Dau_Prob[10000];
   //
   //int KDau_NHits[100000];
   //int P1Dau_NHits[100000];
   //int P2Dau_NHits[100000];
   //int P3Dau_NHits[100000];
   //
   int KDau_NValHits[10000];
   int P1Dau_NValHits[10000];
   int P2Dau_NValHits[10000];
   int P3Dau_NValHits[10000];
   //
   //int KDau_NLostHits[100000];
   //int P1Dau_NLostHits[100000];
   //int P2Dau_NLostHits[100000];
   //int P3Dau_NLostHits[100000];
   // 
   int KDau_NTrackerHits[10000];
   int P1Dau_NTrackerHits[10000];
   int P2Dau_NTrackerHits[10000];
   int P3Dau_NTrackerHits[10000];
   //
   int KDau_NPixelHits[10000];
   int P1Dau_NPixelHits[10000];
   int P2Dau_NPixelHits[10000];
   int P3Dau_NPixelHits[10000];
   //
   int KDau_NStripHits[10000]; 
   int P1Dau_NStripHits[10000];
   int P2Dau_NStripHits[10000];
   int P3Dau_NStripHits[10000];
   //
   int KDau_NTrackerLayers[10000];
   int P1Dau_NTrackerLayers[10000];
   int P2Dau_NTrackerLayers[10000];
   int P3Dau_NTrackerLayers[10000];
   //
   int KDau_NPixelLayers[10000];
   int P1Dau_NPixelLayers[10000];
   int P2Dau_NPixelLayers[10000];
   int P3Dau_NPixelLayers[10000];
   //
   int KDau_NStripLayers[10000];
   int P1Dau_NStripLayers[10000];
   int P2Dau_NStripLayers[10000];
   int P3Dau_NStripLayers[10000];
   //
   int KDau_hasValHitPxB[10000];
   int P1Dau_hasValHitPxB[10000];
   int P2Dau_hasValHitPxB[10000];
   int P3Dau_hasValHitPxB[10000];
   //
   int KDau_hasValHitPxE[10000];
   int P1Dau_hasValHitPxE[10000];
   int P2Dau_hasValHitPxE[10000];
   int P3Dau_hasValHitPxE[10000];
   //
   double KDau_dxyPV[10000];
   double P1Dau_dxyPV[10000];
   double P2Dau_dxyPV[10000];
   double P3Dau_dxyPV[10000];
   //
   double KDau_dxyPVErr[10000];
   double P1Dau_dxyPVErr[10000];
   double P2Dau_dxyPVErr[10000];
   double P3Dau_dxyPVErr[10000];
   //
   double KDau_dzPV[10000];
   double P1Dau_dzPV[10000];
   double P2Dau_dzPV[10000];
   double P3Dau_dzPV[10000];
   //
   double KDau_dzPVErr[10000];
   double P1Dau_dzPVErr[10000];
   double P2Dau_dzPVErr[10000];
   double P3Dau_dzPVErr[10000];
   //
   double KDau_dszPV[10000];
   double P1Dau_dszPV[10000];
   double P2Dau_dszPV[10000];
   double P3Dau_dszPV[10000];
   //
   double KDau_dszPVErr[10000];
   double P1Dau_dszPVErr[10000];
   double P2Dau_dszPVErr[10000];
   double P3Dau_dszPVErr[10000];
   //
   double KDau_dxyPVSignif[10000];
   double P1Dau_dxyPVSignif[10000];
   double P2Dau_dxyPVSignif[10000];
   double P3Dau_dxyPVSignif[10000];
   //
   double KDau_dzPVSignif[10000];
   double P1Dau_dzPVSignif[10000];
   double P2Dau_dzPVSignif[10000];
   double P3Dau_dzPVSignif[10000];
   //
   double KDau_dszPVSignif[10000];
   double P1Dau_dszPVSignif[10000];
   double P2Dau_dszPVSignif[10000];
   double P3Dau_dszPVSignif[10000];
   //
   double K3PiPiSoft_Mass[10000];
   double K3PiPiSoft_pseudoPt[10000];
   double K3PiPiSoft_pseudoEta[10000];
   int K3PiPiSoft_pseudoCharge[10000];
   //
   double K3PiPiSoft_deltaMass[10000];
   //
   double PiSoftDau_Pt[10000];
   double PiSoftDau_Eta[10000];
   double PiSoftDau_Phi[10000];
   double PiSoftDau_NormChi2[10000];
   int PiSoftDau_Charge[10000];
   int PiSoftDau_hasValHitPxB[10000]; 
   int PiSoftDau_hasValHitPxE[10000]; 
   int PiSoftDau_NTrackerHits[10000];
   double PiSoftDau_dxyPVSignif[10000];
  // 
   //- output ROOT file
   TFile* theFile_ ;
   // 
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
DataD0K3Pi::DataD0K3Pi(const edm::ParameterSet& cfg) :
   outputFileName_ (cfg.getUntrackedParameter<string>("histOutputFile") ),
   muons_          (cfg.getParameter<InputTag>( "muons" ) ),
   tracks_         (cfg.getParameter<InputTag>( "tracks" ) ),
   pions_          (cfg.getParameter<InputTag>( "pions" ) ),
   softpions_      (cfg.getParameter<InputTag>( "softpions" ) ),
   kaons_          (cfg.getParameter<InputTag>( "kaons" ) ) {
   //   kaons_          (cfg.getParameter<InputTag>( "kaons" ) ),
   //k3pi_           (cfg.getParameter<InputTag>( "k3pi" ) ) {
}


DataD0K3Pi::~DataD0K3Pi()
{ 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  theFile_->Write() ;
  theFile_->Close() ;
  Tree_=0;
  delete theFile_;
  delete Tree_;
  //
}

//
// member functions
//

// ------------ method called to for each event  ------------
void
DataD0K3Pi::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  //using namespace edm;
  //
  //std::cout << "Processing NEW EVENT: event # " << event.id() <<", time: " << event.time().value() << std::endl;
  //std::cout << "================================================" << endl;
  //
  //
  //////////////////// MC STUFF ////////////

  //Handle<SimTrackContainer> simTracks;
  //event.getByLabel( "g4SimHits", simTracks );
  //hasD0 = false;
  //
  ////Handle<GenParticleCollection> genParticles;
  ////event.getByLabel("genParticles", genParticles);
  ////hasBq = false;
  ////hasGenD0 = false;
  //
  //for (size_t k=0;  k<simTracks->size(); k++ )  
  //   {
  //    const SimTrack & simTr = (*simTracks)[k];
  //    //
      //double Pz = simTr.momentum().pz();
      //double P = sqrt((simTr.momentum().px()*simTr.momentum().px())+
      //                (simTr.momentum().py()*simTr.momentum().py())+
      //                (simTr.momentum().pz()*simTr.momentum().pz()));
      //double Eta = .5*log((1+Pz/P)/(1-Pz/P));
      //if (abs(simTr.type())==13 && simTr.momentum().pt()>=15. && Eta>-2.4 && Eta<2.4)
      //cout << "simtrack id = " << simTr.type() << endl;
  //    if ( abs(simTr.type())==421 ) hasD0 = true;
  //   } 
  // 
  // Get gen Particles and figure out if the event has a B
  //
  ////for (size_t i = 0; i < genParticles->size(); ++ i) 
  ////   {
  ////    const GenParticle & p = (*genParticles)[i];
  ////    int pdgid = p.pdgId();
  ////    if ( abs(pdgid)==421 ) hasGenD0 = true;
  ////    if ( flavour(pdgid)==5 ) hasBq = true;
  ////   }
  //
  //////////////////////////////////////////
  //
  // Get the geometry
  ESHandle<GlobalTrackingGeometry> theTrackingGeometry;
  eventSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);
  //
  // get the track builder
  ESHandle<TransientTrackBuilder> trackBuilder;
  eventSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",trackBuilder);
  //
  ////////////////////// Get primary vertex collection
  //
  Handle<VertexCollection> vertexHandle;
  event.getByLabel("offlinePrimaryVertices",vertexHandle);
  if ( !vertexHandle.isValid() ) return;
  VertexCollection vertexCollection = *(vertexHandle.product());
  int vertex_number = vertexCollection.size();
  //cout << "Number of Primary Vertices = " << vertex_number << endl;
  //
  //////////////////////////////////////
  typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
  typedef ROOT::Math::SVector<double, 3> SVector3;
  //////////////////////////////////////
  reco::Vertex primVertex;
  bool pvfound = (vertex_number != 0);
  if (pvfound)
    {
     PrimaryVertexSorter pvs;
     vector<reco::Vertex> sortedList = pvs.sortedList(*(vertexHandle.product()) );
     primVertex = (sortedList.front());
    } else {
            //creating a dummy PV
            reco::Vertex::Point p(0,0,0);          
            reco::Vertex::Error e;
            e(0,0) = 0.0015*0.0015;
            e(1,1) = 0.0015*0.0015;
            e(2,2) = 15.*15.;
            primVertex = reco::Vertex(p,e,1,1,1);
    }
  //
  Vertex::CovarianceMatrix primVertexCov = primVertex.covariance();
  //
  PVx = primVertex.position().x();
  PVy = primVertex.position().y();
  PVz = primVertex.position().z();
  PVxErr = primVertex.xError();
  PVyErr = primVertex.yError();
  PVzErr = primVertex.zError();
  //
  //GlobalPoint pVertex(primVertex.position().x(),primVertex.position().y(),primVertex.position().z());
  //Basic3DVector<double> pVertexErr(primVertex.xError(),primVertex.yError(),primVertex.zError());
  //
  //
  /////////////////////// not used in this round
  //
//   Handle<BeamSpot> theBeamSpot;
//   event.getByLabel("offlineBeamSpot",theBeamSpot);
//   bool beamspotfound = (theBeamSpot.isValid());
//   //
//   if (theBeamSpot.isValid())
//     {
//       isBSValid = 1;
//     } else {
//             isBSValid = 0;
//     }

  //
  ///////////////////////////////
  //
  // da mettere sotto 
  //se beamspot e' definita
  //if(haveBeamSpot)
  //{ PionBs_dxy[] = track->dxy(theBeamSpot->position());
  // PionBs_dz[] = track->dz(theBeamSpot->position());
  // PionBs_dsz[] = track->dsz(theBeamSpot->position());
  // PionBs_dxy_E[] = track->dxyError());
  // PionBs_dz_E[] = track->dzError());
  // PionBs_dsz_E[] = track->dszError();}
  //
  //else{ PionBs_dxy[nDiPion] =-9999.;
  // PionBs_dz[nDiPion] =-9999. ;
  // PionBs_dsz[nDiPion] =-9999.;
  // PionBs_dxy_E[nDiPion] =-9999.;
  // PionBs_dz_E[nDiPion] =-9999.;
  // PionBs_dsz_E[nDiPion] =-9999.;
  //variabile[];
  /////////////////////// Get run number
  //
  runId = event.id().run();
  evtId = event.id().event();
  lsId = event.luminosityBlock();
  //
  //cout << "Dumping the information of run = " << runId << " , event = " << evtId << " and lumisection = " << lsId << std::endl;
  //
  //OverlapChecker overlap; // for D*'s pi soft checks
  //
  ////////////////////////////////////////////////
  //
  // TRACKS 
  //
  Handle<CandidateView> tracks;
  event.getByLabel(tracks_, tracks);
  ntracks = tracks->size();
  //cout << "ntracks" << ntracks << endl;
  //
  // PIONS
  //
  Handle<CandidateView> kaons;
  event.getByLabel(kaons_, kaons);
  nkaons = kaons->size();
  //
  // KAONS
  //
  Handle<CandidateView> pions;
  event.getByLabel(pions_, pions);
  npions = pions->size();
  //
  //cout << "ntracks = " << ntracks << " & nkaons = " << nkaons << " & npions = " << npions << endl;
  //
  // SOFT PIONS
  //
  Handle<CandidateView> softpions;
  event.getByLabel(softpions_, softpions);
  nsoftpions = softpions->size();
  //
  //////////////////////////////////////////////////
  //
  // MUONS 
  //
  //Handle<CandidateView> glbmuons;
  //event.getByLabel(muons_, glbmuons);
  //nGlbMuons = glbmuons->size();
  //
  Handle<MuonCollection> glbmuons;
  event.getByLabel(muons_, glbmuons);
  nGlbMuons = glbmuons->size();
  //cout << "number of globalmuons = " << glbmuons->size() << endl;
  //
  //////////////////////////////////////////////////
  //
  //Handle<CandidateView> k3pi;
  //event.getByLabel(k3pi_, k3pi);
  //int numk3pi = k3pi->size();
  //cout << "number of K3Pi candidates = " << numk3pi << endl;
  //
  ////////////////////////////////////////////
  const double dstarNominalDeltaMass = 0.145421;
  ///////////////////////////////////////////
  //
  KalmanVertexFitter theKalmanFitter;
  vector<TransientTrack> transTracks; //This vector will hold the tracks to be vertexed  
  //
  nK3Pi = 0; // multeplicity
  //
  ////////////////////////////// preliminary loop for tagging the best candidates (by pt value)
  //
//   double maxtemp_pt = 0.; double max2ndtemp_pt = 0.;
//   double maxtemp_eta = 0.; 
//   double maxtemp_phi = 0.;
//   double best_pt = -1.; double best_eta = -9.; double best_phi = -9.;
//   double best2nd_pt = -1.; double best2nd_eta = -9.; double best2nd_phi= -9.;
//   //
//   for ( int k = 0; k < nk3pi; k++ ) 
//      {
//       const Candidate & kpiCandScan = (*kpi)[ k ];
//       //
//       if (kpiCandScan.pt() > maxtemp_pt) 
//         {
//          best2nd_pt = maxtemp_pt; max2ndtemp_pt = best2nd_pt;
//          best2nd_eta = maxtemp_eta;
//          best2nd_phi = maxtemp_phi;
//          //
//          maxtemp_pt = kpiCandScan.pt(); best_pt = maxtemp_pt;
//          maxtemp_eta = kpiCandScan.eta(); best_eta = maxtemp_eta;
//          maxtemp_phi = kpiCandScan.phi(); best_phi = maxtemp_phi;
//          //
//         }
//       else if(kpiCandScan.pt() > max2ndtemp_pt && kpiCandScan.pt() < maxtemp_pt)
//         {
//          max2ndtemp_pt = kpiCandScan.pt(); best2nd_pt = max2ndtemp_pt;
//          //
//          best2nd_eta = kpiCandScan.eta();
//          best2nd_phi = kpiCandScan.phi();
//         }
//       else {}
  //}
  //
  //cout << "Best PT is = " << best_pt << " while 2nd best is = " << best2nd_pt << endl;
  //cout << "Best ETA is = " << best_eta << " while 2nd best is = " << best2nd_eta << endl;
  //cout << "Best PHI is = " << best_phi << " while 2nd best is = " << best2nd_phi << endl;
  //
  ///////////////////////////////////////////////////////////////////
  //
  ////-- DAUGHTERS
  //
  KDau_P[nK3Pi] = -1.;   P1Dau_P[nK3Pi] = -1.;   P2Dau_P[nK3Pi] = -1.;   P3Dau_P[nK3Pi] = -1.;
  KDau_Pt[nK3Pi] = -1.;  P1Dau_Pt[nK3Pi] = -1.;  P2Dau_Pt[nK3Pi] = -1.;  P3Dau_Pt[nK3Pi] = -1.;
  KDau_Eta[nK3Pi] = -9.; P1Dau_Eta[nK3Pi] = -9.; P2Dau_Eta[nK3Pi] = -9.; P3Dau_Eta[nK3Pi] = -9.; 
  KDau_Phi[nK3Pi] = -9.; P1Dau_Phi[nK3Pi] = -9.; P2Dau_Phi[nK3Pi] = -9.; P3Dau_Phi[nK3Pi] = -9.;
  //
  KDau_NormChi2[nK3Pi] = -1.; P1Dau_NormChi2[nK3Pi] = -1.; P2Dau_NormChi2[nK3Pi] = -1.; P3Dau_NormChi2[nK3Pi] = -1.;
  KDau_Chi2[nK3Pi] = -1.;     P1Dau_Chi2[nK3Pi] = -1.;     P2Dau_Chi2[nK3Pi] = -1.;     P3Dau_Chi2[nK3Pi] = -1.;
  KDau_Prob[nK3Pi] = -1.;     P1Dau_Prob[nK3Pi] = -1.;     P2Dau_Prob[nK3Pi] = -1.;     P3Dau_Prob[nK3Pi] = -1.;
  //
  KDau_Charge[nK3Pi] = 0.; P1Dau_Charge[nK3Pi] = 0.; P2Dau_Charge[nK3Pi] = 0.; P3Dau_Charge[nK3Pi] = 0.; 
  //
  //KDau_NHits[nK3Pi] = -1;     P1Dau_NHits[nK3Pi] = -1;     P2Dau_NHits[nK3Pi] = -1;     P3Dau_NHits[nK3Pi] = -1;
  KDau_NValHits[nK3Pi] = -1;  P1Dau_NValHits[nK3Pi] = -1;  P2Dau_NValHits[nK3Pi] = -1;  P3Dau_NValHits[nK3Pi] = -1;
  //KDau_NLostHits[nK3Pi] = -1; P1Dau_NLostHits[nK3Pi] = -1; P2Dau_NLostHits[nK3Pi] = -1; P3Dau_NLostHits[nK3Pi] = -1;
  //
  KDau_NTrackerHits[nK3Pi] = -1; P1Dau_NTrackerHits[nK3Pi] = -1; P2Dau_NTrackerHits[nK3Pi] = -1; P3Dau_NTrackerHits[nK3Pi] = -1;
  KDau_NPixelHits[nK3Pi] = -1;   P1Dau_NPixelHits[nK3Pi] = -1;   P2Dau_NPixelHits[nK3Pi] = -1;   P3Dau_NPixelHits[nK3Pi] = -1;
  KDau_NStripHits[nK3Pi] = -1;   P1Dau_NStripHits[nK3Pi] = -1;   P2Dau_NStripHits[nK3Pi] = -1;   P3Dau_NStripHits[nK3Pi] = -1;
  //
  KDau_NTrackerLayers[nK3Pi] = -1; P1Dau_NTrackerLayers[nK3Pi] = -1; P2Dau_NTrackerLayers[nK3Pi] = -1; P3Dau_NTrackerLayers[nK3Pi] = -1;
  KDau_NPixelLayers[nK3Pi] = -1; P1Dau_NPixelLayers[nK3Pi] = -1; P2Dau_NPixelLayers[nK3Pi] = -1;       P3Dau_NPixelLayers[nK3Pi] = -1;
  KDau_NStripLayers[nK3Pi] = -1; P1Dau_NStripLayers[nK3Pi] = -1; P2Dau_NStripLayers[nK3Pi] = -1;       P3Dau_NStripLayers[nK3Pi] = -1;
  //
  KDau_hasValHitPxB[nK3Pi] = 0; P1Dau_hasValHitPxB[nK3Pi] = 0; P2Dau_hasValHitPxB[nK3Pi] = 0; P3Dau_hasValHitPxB[nK3Pi] = 0;
  KDau_hasValHitPxE[nK3Pi] = 0; P1Dau_hasValHitPxE[nK3Pi] = 0; P2Dau_hasValHitPxE[nK3Pi] = 0; P3Dau_hasValHitPxE[nK3Pi] = 0;
  //
  KDau_dxyPV[nK3Pi] = -999.;  P1Dau_dxyPV[nK3Pi] = -999.;  P2Dau_dxyPV[nK3Pi] = -999.;  P3Dau_dxyPV[nK3Pi] = -999.;
  KDau_dxyPVErr[nK3Pi] = -1.; P1Dau_dxyPVErr[nK3Pi] = -1.; P2Dau_dxyPVErr[nK3Pi] = -1.; P3Dau_dxyPVErr[nK3Pi] = -1.;
  KDau_dzPV[nK3Pi] = -999.;   P1Dau_dzPV[nK3Pi] = -999.;   P2Dau_dzPV[nK3Pi] = -999.;   P3Dau_dzPV[nK3Pi] = -999.; 
  KDau_dzPVErr[nK3Pi] = -1.;  P1Dau_dzPVErr[nK3Pi] = -1.;  P2Dau_dzPVErr[nK3Pi] = -1.;  P3Dau_dzPVErr[nK3Pi] = -1.; 
  KDau_dszPV[nK3Pi] = -999.;  P1Dau_dszPV[nK3Pi] = -999.;  P2Dau_dszPV[nK3Pi] = -999.;  P3Dau_dszPV[nK3Pi] = -999.; 
  KDau_dszPVErr[nK3Pi] = -1.; P1Dau_dszPVErr[nK3Pi] = -1.; P2Dau_dszPVErr[nK3Pi] = -1.; P3Dau_dszPVErr[nK3Pi] = -1.; 
  //
  KDau_dxyPVSignif[nK3Pi] = -9.; P1Dau_dxyPVSignif[nK3Pi] = -9.; P2Dau_dxyPVSignif[nK3Pi] = -9.; P3Dau_dxyPVSignif[nK3Pi] = -9.;
  KDau_dzPVSignif[nK3Pi] = -9.;  P1Dau_dzPVSignif[nK3Pi] = -9.;  P2Dau_dzPVSignif[nK3Pi] = -9.;  P3Dau_dzPVSignif[nK3Pi] = -9.; 
  KDau_dszPVSignif[nK3Pi] = -9.; P1Dau_dszPVSignif[nK3Pi] = -9.; P2Dau_dszPVSignif[nK3Pi] = -9.; P3Dau_dszPVSignif[nK3Pi] = -9.;
  //
  /////////////////////-- MOTHER
  //
  K3Pi_Mass[nK3Pi] = 0.;             // D0 candidates (invariant mass)
  K3Pi_pseudoP[nK3Pi] = -1.;         // D0 pseudo-momentum
  K3Pi_pseudoPt[nK3Pi] = -1.;        // D0 pseudo-transverse-momentum
  K3Pi_pseudoEta[nK3Pi] = -9.;       // D0 pseudo-eta
  K3Pi_pseudoPhi[nK3Pi] = -9.;       // D0 pseudo-phi
  K3Pi_pseudoCharge[nK3Pi] = 3.;    // D0 pseudo-charge (check: should be always zero)
  //
  K3Pi_deltaR_kp1[nK3Pi] = -99;        // deltaR different angles
  K3Pi_deltaR_kp2[nK3Pi] = -99;        
  K3Pi_deltaR_kp3[nK3Pi] = -99;        
  K3Pi_deltaR_p1p3[nK3Pi] = -99;        
  K3Pi_deltaR_p2p3[nK3Pi]  = -99;        
  K3Pi_deltaR_p1p2[nK3Pi] = -99;        
  //
  K3Pi_pointingAngle[nK3Pi] = -9;   // pointing angle 3D
  K3Pi_2DpointingAngle[nK3Pi] = -9; // pointing angle 2D
  //
  K3Pi_VtxChi2[nK3Pi] = -1.;          // normalized chi2 of K-Pi vertex fit
  K3Pi_VtxChi2Prob[nK3Pi] = -1.;      // chi2 probability of K-Pi vertex fit
  K3Pi_VtxCalcDistXY[nK3Pi] = -1.;    // D0 vertex transverse distance calculated (w.r.t PV)
  K3Pi_VtxDistXY[nK3Pi] = -9.;        // D0 vertex transverse distance (w.r.t PV)
  K3Pi_VtxDistXYError[nK3Pi] = -9.;   // D0 vertex transverse distance error
  K3Pi_VtxDistXYSignif[nK3Pi] = -9.;  // D0 vertex transverse distance significance (distance/error)
  //
  /////// ASSOCIATED BEST D*
  //
  K3PiPiSoft_Mass[nK3Pi] = 0.;
  K3PiPiSoft_pseudoPt[nK3Pi] = -1.;
  K3PiPiSoft_pseudoEta[nK3Pi] = -9.;
  K3PiPiSoft_pseudoCharge[nK3Pi] = 3;
  //
  K3PiPiSoft_deltaMass[nK3Pi] = 0.;
  //
  PiSoftDau_Pt[nK3Pi] = -1.;
  PiSoftDau_Eta[nK3Pi] = -9.;
  PiSoftDau_Phi[nK3Pi] = -9.;
  PiSoftDau_NormChi2[nK3Pi] = -1;
  PiSoftDau_Charge[nK3Pi] = 3;
  PiSoftDau_hasValHitPxB[nK3Pi] = 0; 
  PiSoftDau_hasValHitPxE[nK3Pi] = 0; 
  PiSoftDau_NTrackerHits[nK3Pi] = -1;
  PiSoftDau_dxyPVSignif[nK3Pi] = -9.;
  //
  /////////////////////////////// CALCOLO DEL COMBINATORIO ////////////
  //
  //CompositeCandidate k3piCand;
  std::list<CompositeCandidate> myD0cands;
  std::list<CompositeCandidate> myDstarCands;
  //
  for ( int i=0; i<nkaons; i++)
     { 
       const Candidate & kCand = (*kaons)[i];
       TrackRef kTrack = kCand.get<TrackRef>();
       //
       for (int j=0; j<npions; j++)
	  {
	   const Candidate & p1Cand = (*pions)[j];
           TrackRef p1Track = p1Cand.get<TrackRef>();
	   if( p1Track->pt() == kTrack->pt() && p1Track->eta() == kTrack->eta() && p1Track->phi() == kTrack->phi() ) continue;
	   //
	   for ( int jj=j+1; jj<npions; jj++)
	      {
	       const Candidate & p2Cand = (*pions)[jj];
	       TrackRef p2Track = p2Cand.get<TrackRef>();
	       //
	       if ( ((kTrack->charge())+(p1Track->charge())+(p2Track->charge()))==3 ) continue;
	       if( p2Track->pt() == kTrack->pt() && p2Track->eta() == kTrack->eta() && p2Track->phi() == kTrack->phi() ) continue;
	       //
	       for ( int jjj=jj+1; jjj<npions; jjj++)
		  {
		   const Candidate & p3Cand = (*pions)[jjj];
		   TrackRef p3Track = p3Cand.get<TrackRef>();
		   //
		   if( p3Track->pt() == kTrack->pt() && p3Track->eta() == kTrack->eta() && p3Track->phi() == kTrack->phi() ) continue;
		   if ( ((kTrack->charge())+(p1Track->charge())+(p2Track->charge())+(p3Track->charge()) ) != 0 ) continue;
		   if( kTrack->pt() < 1.0 && p1Track->pt() < 1.0 && p2Track->pt() < 1.0 && p3Track->pt() < 1.0 ) continue; // at least 1 track with pt>1.0
		   if( kTrack->normalizedChi2()>2.5 || p1Track->normalizedChi2()>2.5 || p2Track->normalizedChi2()>2.5 || p3Track->normalizedChi2()>2.5 ) continue;
                   if( kTrack->hitPattern().numberOfValidTrackerHits()<7 || p1Track->hitPattern().numberOfValidTrackerHits()<7 || 
		       p2Track->hitPattern().numberOfValidTrackerHits()<7 || p3Track->hitPattern().numberOfValidTrackerHits()<7) continue;
		   //
		   bool tobevetoed = false;
		   double ptveto; double etaveto; 
		   //
		   //cout << "The 4 tracks have pt = " << kTrack->pt() <<" , "<< p1Track->pt() <<" , "<< p2Track->pt() <<" , "<< p3Track->pt() << endl;
		   //
		   for (unsigned int iMuon = 0; iMuon < glbmuons->size(); ++iMuon)
		      {
			ptveto = glbmuons->at(iMuon).innerTrack()->pt();
			etaveto = glbmuons->at(iMuon).innerTrack()->eta();
			//
			if ( ( (kTrack->pt() == ptveto) && (kTrack->eta() == etaveto)  ) ||
			     ( (p1Track->pt() == ptveto) && (p1Track->eta() == etaveto)  ) ||
			     ( (p2Track->pt() == ptveto) && (p2Track->eta() == etaveto) ) ||
			     ( (p3Track->pt() == ptveto) && (p3Track->eta() == etaveto) ) )
			  {
			    tobevetoed = true;
			  } else {
			          tobevetoed = false;
			         } 
		      }
		   //
		   if( tobevetoed == true) continue;
		   //
		   CompositeCandidate k3piCand;
		   k3piCand.addDaughter(kCand);
		   k3piCand.addDaughter(p1Cand);
		   k3piCand.addDaughter(p2Cand);
		   k3piCand.addDaughter(p3Cand);
		   AddFourMomenta addP4;
		   //k3piCand.setPdgId(421); /// cosa farebbe?
		   addP4.set( k3piCand );
		   if ( (k3piCand.p4().mass() < 1.73 ) || ( k3piCand.p4().mass() > 2.0 ) ) continue;
		   //
		   myD0cands.push_back( k3piCand );
		   //
		   /////////////////////////////
		   //
		   vector<double> deltaMdist(nsoftpions);
		   vector<double>::iterator pos;
		   double dstarMgoodUpForCheck = -998.;
		   double dstarMgoodDownForCheck = -998.;
		   //
		   AddFourMomenta sumP4;
		   CompositeCandidate dstarCand;
		   dstarCand.addDaughter(k3piCand);
		   //
		   std::list<CompositeCandidate> allDstarCands;
		   //
		   for ( int l=0; l<nsoftpions; l++ )
		      { 
                       deltaMdist[l] = 9999.;
		       //
		       const Candidate & softpCand = (*softpions)[l]; 
		       TrackRef softTrack = softpCand.get<TrackRef>();
		       //
		       //if( softTrack->hitPattern().numberOfValidTrackerHits() < 3 ) continue;
		       //if( softTrack->normalizedChi2() > 2.5) continue;
		       //
		       if ( ( (kTrack->pt() == softTrack->pt()) && (kTrack->eta() == softTrack->eta()) ) ||
			    ( (p1Track->pt() == softTrack->pt()) && (p1Track->eta() == softTrack->eta()) ) ||
			    ( (p2Track->pt() == softTrack->pt()) && (p2Track->eta() == softTrack->eta()) ) ||
			    ( (p3Track->pt() == softTrack->pt()) && (p3Track->eta() == softTrack->eta()) ) ) continue;
		       //
		       dstarCand.addDaughter(softpCand);
		       sumP4.set(dstarCand);
		       //
		       double deltaMass = 9999.;
		       deltaMass = (dstarCand.p4().mass()- k3piCand.p4().mass());
		       deltaMdist[l] = abs(deltaMass - dstarNominalDeltaMass);
		       //cout << "index l = " << l << " and deltaMdist = " << deltaMdist[l] << endl;
		       //
		       //if ( dstarCand.p4().mass()<2.15 && dstarCand.p4().mass()>1.87 ) && deltaM < 0.18 )
		       //{
		       allDstarCands.push_back( dstarCand );
		       //}
		      } // end loop on softpions
		   //
		   //////////
		   pos = min_element(deltaMdist.begin(),deltaMdist.end()); // position in the array with value being the minimum distance
		   double deltaMdistMin = *pos;
		   //cout << "MINIMUM DELTA_M = " << deltaMdistMin << endl;
		   double deltaMgoodUp = deltaMdistMin + dstarNominalDeltaMass;
		   double deltaMgoodDown = dstarNominalDeltaMass - deltaMdistMin;    
		   //
		   double dstarMgoodUp = k3piCand.p4().mass() + deltaMgoodUp;
		   double dstarMgoodDown = k3piCand.p4().mass() + deltaMgoodDown; 
		   //
		   dstarMgoodUpForCheck = dstarMgoodUp;
		   dstarMgoodDownForCheck = dstarMgoodDown;
		   //
		   //cout << "2 D* masses possible are = " << dstarMgoodUpForCheck << " & " << dstarMgoodDownForCheck << endl;
		   //cout << "ALL D* cands list size = " << allDstarCands.size() << endl;
		   //
		   for (std::list<CompositeCandidate>::iterator allDstar = allDstarCands.begin(); allDstar != allDstarCands.end(); ++ allDstar )  
		      {
			const Candidate & scanDstar = * allDstar ;
			if ( scanDstar.p4().mass() == dstarMgoodDownForCheck || scanDstar.p4().mass() == dstarMgoodUpForCheck)
			  {
			    //cout << "FOUND MY BEST DSTAR!!!!!! with mass = " << scanDstar.p4().mass() << endl;
			   myDstarCands.push_back( *allDstar );
			  }
		      }

		   ///////////
		   //
		   //
		  } // end loop on jjj
	       //
	      } // end loop on jj
	   //
	  } // end loop on j
     } // end loop on i
  //
  ////////////////////////////////  
  //
  //cout << "D0 cands list size = " << myD0cands.size() << " that should be equal (or <=)to the D* cands list size which is = " << myDstarCands.size() << endl;
  //
  for (std::list<CompositeCandidate>::iterator myD0 = myD0cands.begin(); myD0 != myD0cands.end(); ++ myD0 )  
     {
       transTracks.clear();
       //
       const Candidate & myk3pi = * myD0 ;
       //
       /////////////////////////////////////
       //
       const Candidate * kdau = myD0->daughter(0);
       const Candidate * pi1dau = myD0->daughter(1);
       const Candidate * pi2dau = myD0->daughter(2);
       const Candidate * pi3dau = myD0->daughter(3); 
       //
       const Candidate & kdauCand = *(myD0->daughter(0));
       const Candidate & pi1dauCand = *(myD0->daughter(1));
       const Candidate & pi2dauCand = *(myD0->daughter(2));
       const Candidate & pi3dauCand = *(myD0->daughter(3)); 
       //
       TrackRef kDauTrack = kdau->get<TrackRef>();
       TrackRef pi1DauTrack = pi1dau->get<TrackRef>();
       TrackRef pi2DauTrack = pi2dau->get<TrackRef>();
       TrackRef pi3DauTrack = pi3dau->get<TrackRef>();
       // 
       ////////////////////////////////////
       //
       //cout << "READING K3PI cand with mass= " << myD0->p4().mass() << endl;
       //cout << "The 4 daughters have pt = " << kdau->pt() << " , " << pi1dau->pt() << " , " << pi2dau->pt() << " , " << pi3dau->pt() << endl;
       //
       K3Pi_Mass[nK3Pi] = myD0->p4().mass();
       K3Pi_pseudoP[nK3Pi] = myD0->p();
       K3Pi_pseudoPt[nK3Pi] = myD0->p4().pt();
       K3Pi_pseudoEta[nK3Pi] = myD0->p4().eta();
       K3Pi_pseudoPhi[nK3Pi] = myD0->p4().phi();
       K3Pi_pseudoCharge[nK3Pi] = myD0->charge();
       //
       K3Pi_deltaR_kp1[nK3Pi] = sqrt( ( (kdau->eta() - pi1dau->eta()) * (kdau->eta() - pi1dau->eta()) ) + ( (kdau->phi() - pi1dau->phi()) * (kdau->phi() - pi1dau->phi()) ) );
       K3Pi_deltaR_kp2[nK3Pi] = sqrt( ( (kdau->eta() - pi2dau->eta()) * (kdau->eta() - pi2dau->eta()) ) + ( (kdau->phi() - pi2dau->phi()) * (kdau->phi() - pi2dau->phi()) ) );
       K3Pi_deltaR_kp3[nK3Pi] = sqrt(( (kdau->eta() - pi3dau->eta())*(kdau->eta() - pi3dau->eta()) ) + ( (kdau->phi() - pi3dau->phi())*(kdau->phi() - pi3dau->phi()) ));
       K3Pi_deltaR_p1p3[nK3Pi] = sqrt( ( (pi1dau->eta() - pi3dau->eta()) * (pi1dau->eta() - pi3dau->eta()) ) + 
                                       ( (pi1dau->phi() - pi3dau->phi()) * (pi1dau->phi() - pi3dau->phi()) ) );
       K3Pi_deltaR_p2p3[nK3Pi] = sqrt( ( (pi2dau->eta() - pi3dau->eta()) * (pi2dau->eta() - pi3dau->eta()) ) + 
                                       ( (pi2dau->phi() - pi3dau->phi()) * (pi2dau->phi() - pi3dau->phi()) ) );
       K3Pi_deltaR_p1p2[nK3Pi] = sqrt( ( (pi1dau->eta() - pi2dau->eta()) * (pi1dau->eta() - pi2dau->eta()) ) + 
                                       ( (pi1dau->phi() - pi2dau->phi()) * (pi1dau->phi() - pi2dau->phi()) ) );
       //
       //---- variables related to candidate daughters
       //
       // Fill the vector of TransientTracks 
       transTracks.push_back( trackBuilder->build(kDauTrack) );
       transTracks.push_back( trackBuilder->build(pi1DauTrack) );
       transTracks.push_back( trackBuilder->build(pi2DauTrack) );
       transTracks.push_back( trackBuilder->build(pi3DauTrack) );
       //
       // K daughter
       //
       KDau_P[nK3Pi] = kDauTrack->p();
       KDau_Pt[nK3Pi] = kDauTrack->pt();   // KDau_PtErr[nK3Pi] = kDauTrack->ptError();
       KDau_Eta[nK3Pi] = kDauTrack->eta(); // KDau_EtaErr[nK3Pi] = kDauTrack->etaError();
       KDau_Phi[nK3Pi] = kDauTrack->phi(); // KDau_PhiErr[nK3Pi] = kDauTrack->phiError();
       //
       KDau_NormChi2[nK3Pi] = kDauTrack->normalizedChi2();
       KDau_Chi2[nK3Pi] = kDauTrack->chi2();               // KDau_Ndof[nK3Pi] = kDauTrack->ndof();
       KDau_Prob[nK3Pi] = TMath::Prob(kDauTrack->chi2(),(int)kDauTrack->ndof());
       //
       KDau_Charge[nK3Pi] = kDauTrack->charge();           // KDau_QoverP[nK3Pi] = kDauTrack_ptHi->qoverp();
       //
       //KDau_NHits[nK3Pi] = kDauTrack->hitPattern().numberOfHits();
       KDau_NValHits[nK3Pi] = kDauTrack->hitPattern().numberOfValidHits();
       //KDau_NLostHits[nK3Pi] = ( kDauTrack->hitPattern().numberOfHits() - kDauTrack->hitPattern().numberOfValidHits() );
       //
       KDau_NTrackerHits[nK3Pi] = kDauTrack->hitPattern().numberOfValidTrackerHits(); // tracker (pixel+strip) valid hits
       KDau_NPixelHits[nK3Pi] = kDauTrack->hitPattern().numberOfValidPixelHits();     // pixel valid hits
       KDau_NStripHits[nK3Pi] = kDauTrack->hitPattern().numberOfValidStripHits();     // strip valid hits
       //
       KDau_NTrackerLayers[nK3Pi] = kDauTrack->hitPattern().trackerLayersWithMeasurement(); // tracker (pixel+strip) layers
       KDau_NPixelLayers[nK3Pi] = kDauTrack->hitPattern().pixelLayersWithMeasurement();     // pixel layers
       KDau_NStripLayers[nK3Pi] = kDauTrack->hitPattern().stripLayersWithMeasurement();     // strip layers
       //
       KDau_hasValHitPxB[nK3Pi] = kDauTrack->hitPattern().hasValidHitInFirstPixelBarrel();   // has valid hit in pixel barrel
       KDau_hasValHitPxE[nK3Pi] = kDauTrack->hitPattern().hasValidHitInFirstPixelEndcap();   // has valid hit in pixel endcap (on newer releases)
       //
       KDau_dxyPV[nK3Pi] = kDauTrack->dxy(primVertex.position());
       KDau_dxyPVErr[nK3Pi] = kDauTrack->dxyError();
       KDau_dzPV[nK3Pi] = kDauTrack->dz(primVertex.position());
       KDau_dzPVErr[nK3Pi] = kDauTrack->dzError();
       KDau_dszPV[nK3Pi] = kDauTrack->dsz(primVertex.position());
       KDau_dszPVErr[nK3Pi] = kDauTrack->dszError();
       //
       KDau_dxyPVSignif[nK3Pi] = abs(kDauTrack->dxy(primVertex.position())) / kDauTrack->dxyError();
       KDau_dzPVSignif[nK3Pi] = abs(kDauTrack->dz(primVertex.position())) / kDauTrack->dzError();
       KDau_dszPVSignif[nK3Pi] = abs(kDauTrack->dsz(primVertex.position())) / kDauTrack->dszError();
       //
       //- PION 1
       //
       P1Dau_P[nK3Pi] = pi1DauTrack->p();
       P1Dau_Pt[nK3Pi] = pi1DauTrack->pt();   
       P1Dau_Eta[nK3Pi] = pi1DauTrack->eta();
       P1Dau_Phi[nK3Pi] = pi1DauTrack->phi(); 
       //
       P1Dau_NormChi2[nK3Pi] = pi1DauTrack->normalizedChi2();
       P1Dau_Chi2[nK3Pi] = pi1DauTrack->chi2();               
       P1Dau_Prob[nK3Pi] = TMath::Prob(pi1DauTrack->chi2(),(int)pi1DauTrack->ndof());
       //
       P1Dau_Charge[nK3Pi] = pi1DauTrack->charge();         
       //
       //P1Dau_NHits[nK3Pi] = pi1DauTrack->hitPattern().numberOfHits();
       P1Dau_NValHits[nK3Pi] = pi1DauTrack->hitPattern().numberOfValidHits();
       //P1Dau_NLostHits[nK3Pi] = ( pi1DauTrack->hitPattern().numberOfHits() - pi1DauTrack->hitPattern().numberOfValidHits() );
       //
       P1Dau_NTrackerHits[nK3Pi] = pi1DauTrack->hitPattern().numberOfValidTrackerHits(); // tracker (pixel+strip) valid hits
       P1Dau_NPixelHits[nK3Pi] = pi1DauTrack->hitPattern().numberOfValidPixelHits();     // pixel valid hits
       P1Dau_NStripHits[nK3Pi] = pi1DauTrack->hitPattern().numberOfValidStripHits();     // strip valid hits
       //
       P1Dau_NTrackerLayers[nK3Pi] = pi1DauTrack->hitPattern().trackerLayersWithMeasurement(); // tracker (pixel+strip) layers
       P1Dau_NPixelLayers[nK3Pi] = pi1DauTrack->hitPattern().pixelLayersWithMeasurement();     // pixel layers
       P1Dau_NStripLayers[nK3Pi] = pi1DauTrack->hitPattern().stripLayersWithMeasurement();     // strip layers
       //
       P1Dau_hasValHitPxB[nK3Pi] = pi1DauTrack->hitPattern().hasValidHitInFirstPixelBarrel();   // has valid hit in pixel barrel
       P1Dau_hasValHitPxE[nK3Pi] = pi1DauTrack->hitPattern().hasValidHitInFirstPixelEndcap();   // has valid hit in pixel endcap (on newer releases)
       //
       P1Dau_dxyPV[nK3Pi] = pi1DauTrack->dxy(primVertex.position());
       P1Dau_dxyPVErr[nK3Pi] = pi1DauTrack->dxyError();
       P1Dau_dzPV[nK3Pi] = pi1DauTrack->dz(primVertex.position());
       P1Dau_dzPVErr[nK3Pi] = pi1DauTrack->dzError();
       P1Dau_dszPV[nK3Pi] = pi1DauTrack->dsz(primVertex.position());
       P1Dau_dszPVErr[nK3Pi] = pi1DauTrack->dszError();
       //
       P1Dau_dxyPVSignif[nK3Pi] = abs(pi1DauTrack->dxy(primVertex.position())) / pi1DauTrack->dxyError();
       P1Dau_dzPVSignif[nK3Pi] = abs(pi1DauTrack->dz(primVertex.position())) / pi1DauTrack->dzError();
       P1Dau_dszPVSignif[nK3Pi] = abs(pi1DauTrack->dsz(primVertex.position())) / pi1DauTrack->dszError();
       //
       //- PION 2
       //
       P2Dau_P[nK3Pi] = pi2DauTrack->p();
       P2Dau_Pt[nK3Pi] = pi2DauTrack->pt();   
       P2Dau_Eta[nK3Pi] = pi2DauTrack->eta();
       P2Dau_Phi[nK3Pi] = pi2DauTrack->phi(); 
       //
       P2Dau_NormChi2[nK3Pi] = pi2DauTrack->normalizedChi2();
       P2Dau_Chi2[nK3Pi] = pi2DauTrack->chi2();               
       P2Dau_Prob[nK3Pi] = TMath::Prob(pi2DauTrack->chi2(),(int)pi2DauTrack->ndof());
       //
       P2Dau_Charge[nK3Pi] = pi2DauTrack->charge();         
       //
       //P2Dau_NHits[nK3Pi] = pi2DauTrack->hitPattern().numberOfHits();
       P2Dau_NValHits[nK3Pi] = pi2DauTrack->hitPattern().numberOfValidHits();
       //P2Dau_NLostHits[nK3Pi] = ( pi2DauTrack->hitPattern().numberOfHits() - pi2DauTrack->hitPattern().numberOfValidHits() );
       //
       P2Dau_NTrackerHits[nK3Pi] = pi2DauTrack->hitPattern().numberOfValidTrackerHits(); // tracker (pixel+strip) valid hits
       P2Dau_NPixelHits[nK3Pi] = pi2DauTrack->hitPattern().numberOfValidPixelHits();     // pixel valid hits
       P2Dau_NStripHits[nK3Pi] = pi2DauTrack->hitPattern().numberOfValidStripHits();     // strip valid hits
       //
       P2Dau_NTrackerLayers[nK3Pi] = pi2DauTrack->hitPattern().trackerLayersWithMeasurement(); // tracker (pixel+strip) layers
       P2Dau_NPixelLayers[nK3Pi] = pi2DauTrack->hitPattern().pixelLayersWithMeasurement();     // pixel layers
       P2Dau_NStripLayers[nK3Pi] = pi2DauTrack->hitPattern().stripLayersWithMeasurement();     // strip layers
       //
       P2Dau_hasValHitPxB[nK3Pi] = pi2DauTrack->hitPattern().hasValidHitInFirstPixelBarrel();   // has valid hit in pixel barrel
       P2Dau_hasValHitPxE[nK3Pi] = pi2DauTrack->hitPattern().hasValidHitInFirstPixelEndcap();   // has valid hit in pixel endcap (on newer releases)
       //
       P2Dau_dxyPV[nK3Pi] = pi2DauTrack->dxy(primVertex.position());
       P2Dau_dxyPVErr[nK3Pi] = pi2DauTrack->dxyError();
       P2Dau_dzPV[nK3Pi] = pi2DauTrack->dz(primVertex.position());
       P2Dau_dzPVErr[nK3Pi] = pi2DauTrack->dzError();
       P2Dau_dszPV[nK3Pi] = pi2DauTrack->dsz(primVertex.position());
       P2Dau_dszPVErr[nK3Pi] = pi2DauTrack->dszError();
       //
       P2Dau_dxyPVSignif[nK3Pi] = abs(pi2DauTrack->dxy(primVertex.position())) / pi2DauTrack->dxyError();
       P2Dau_dzPVSignif[nK3Pi] = abs(pi2DauTrack->dz(primVertex.position())) / pi2DauTrack->dzError();
       P2Dau_dszPVSignif[nK3Pi] = abs(pi2DauTrack->dsz(primVertex.position())) / pi2DauTrack->dszError();
       //
       //- PION 3
       //
       P3Dau_P[nK3Pi] = pi3DauTrack->p();
       P3Dau_Pt[nK3Pi] = pi3DauTrack->pt();   
       P3Dau_Eta[nK3Pi] = pi3DauTrack->eta();
       P3Dau_Phi[nK3Pi] = pi3DauTrack->phi(); 
       //
       P3Dau_NormChi2[nK3Pi] = pi3DauTrack->normalizedChi2();
       P3Dau_Chi2[nK3Pi] = pi3DauTrack->chi2();               
       P3Dau_Prob[nK3Pi] = TMath::Prob(pi3DauTrack->chi2(),(int)pi3DauTrack->ndof());
       //
       P3Dau_Charge[nK3Pi] = pi3DauTrack->charge();         
       //
       //P3Dau_NHits[nK3Pi] = pi3DauTrack->hitPattern().numberOfHits();
       P3Dau_NValHits[nK3Pi] = pi3DauTrack->hitPattern().numberOfValidHits();
       //P3Dau_NLostHits[nK3Pi] = ( pi3DauTrack->hitPattern().numberOfHits() - pi3DauTrack->hitPattern().numberOfValidHits() );
       //
       P3Dau_NTrackerHits[nK3Pi] = pi3DauTrack->hitPattern().numberOfValidTrackerHits(); // tracker (pixel+strip) valid hits
       P3Dau_NPixelHits[nK3Pi] = pi3DauTrack->hitPattern().numberOfValidPixelHits();     // pixel valid hits
       P3Dau_NStripHits[nK3Pi] = pi3DauTrack->hitPattern().numberOfValidStripHits();     // strip valid hits
       //
       P3Dau_NTrackerLayers[nK3Pi] = pi3DauTrack->hitPattern().trackerLayersWithMeasurement(); // tracker (pixel+strip) layers
       P3Dau_NPixelLayers[nK3Pi] = pi3DauTrack->hitPattern().pixelLayersWithMeasurement();     // pixel layers
       P3Dau_NStripLayers[nK3Pi] = pi3DauTrack->hitPattern().stripLayersWithMeasurement();     // strip layers
       //
       P3Dau_hasValHitPxB[nK3Pi] = pi3DauTrack->hitPattern().hasValidHitInFirstPixelBarrel();   // has valid hit in pixel barrel
       P3Dau_hasValHitPxE[nK3Pi] = pi3DauTrack->hitPattern().hasValidHitInFirstPixelEndcap();   // has valid hit in pixel endcap (on newer releases)
       //
       P3Dau_dxyPV[nK3Pi] = pi3DauTrack->dxy(primVertex.position());
       P3Dau_dxyPVErr[nK3Pi] = pi3DauTrack->dxyError();
       P3Dau_dzPV[nK3Pi] = pi3DauTrack->dz(primVertex.position());
       P3Dau_dzPVErr[nK3Pi] = pi3DauTrack->dzError();
       P3Dau_dszPV[nK3Pi] = pi3DauTrack->dsz(primVertex.position());
       P3Dau_dszPVErr[nK3Pi] = pi3DauTrack->dszError();
       //
       P3Dau_dxyPVSignif[nK3Pi] = abs(pi3DauTrack->dxy(primVertex.position())) / pi3DauTrack->dxyError();
       P3Dau_dzPVSignif[nK3Pi] = abs(pi3DauTrack->dz(primVertex.position())) / pi3DauTrack->dzError();
       P3Dau_dszPVSignif[nK3Pi] = abs(pi3DauTrack->dsz(primVertex.position())) / pi3DauTrack->dszError();
       //
       /////////////////////////////////////////////////////////////////////
       //
       //double mydeltaMbest = 0.;
       bool dstar_exist = false;
       //
       for (std::list<CompositeCandidate>::iterator myDstar = myDstarCands.begin(); myDstar != myDstarCands.end(); ++ myDstar )  
          { 
	   const Candidate & myBestDstar = * myDstar;
	   //
	   const Candidate * d0dauCand = myDstar->daughter(0);
	   const Candidate * psoftDauCand = myDstar->daughter(1); 
	   TrackRef psoftDauTrack = psoftDauCand->get<TrackRef>();
	   //
	   if ( d0dauCand->p4().mass() == myD0->p4().mass() )
	     {
	       dstar_exist = true;
	       // I enter here (hopefully) only once per D0 candidate!
	       // cout << "dstar daughter mass is = " << d0dauCand->p4().mass() << endl;
	       // cout << "k3pi mass is = " << myD0->p4().mass() << endl; // cout << "psoft pt = " << psoftDauTrack->pt() << endl;
	       //
	       //cout << "My best D* mass = " << myBestDstar.p4().mass() << endl;
	       // mydeltaMbest = myBestDstar.p4().mass() - d0dauCand->p4().mass() ;
	       //cout << "My deltaM best = " << mydeltaMbest << endl;
	       //
	       K3PiPiSoft_Mass[nK3Pi] = myBestDstar.p4().mass();
	       K3PiPiSoft_pseudoPt[nK3Pi] = myBestDstar.p4().pt();
	       K3PiPiSoft_pseudoEta[nK3Pi] = myBestDstar.p4().eta();
	       K3PiPiSoft_pseudoCharge[nK3Pi] = myBestDstar.charge();
	       //
	       K3PiPiSoft_deltaMass[nK3Pi] = myBestDstar.p4().mass() - d0dauCand->p4().mass() ;
	       //
	       PiSoftDau_Pt[nK3Pi] = psoftDauTrack->pt();
	       PiSoftDau_Eta[nK3Pi] = psoftDauTrack->eta();
	       PiSoftDau_Phi[nK3Pi] = psoftDauTrack->phi();
	       PiSoftDau_NormChi2[nK3Pi] = psoftDauTrack->normalizedChi2();
	       PiSoftDau_Charge[nK3Pi] = psoftDauTrack->charge();
	       PiSoftDau_hasValHitPxB[nK3Pi] = psoftDauTrack->hitPattern().hasValidHitInFirstPixelBarrel(); 
	       PiSoftDau_hasValHitPxE[nK3Pi] = psoftDauTrack->hitPattern().hasValidHitInFirstPixelEndcap(); 
	       PiSoftDau_NTrackerHits[nK3Pi] = psoftDauTrack->hitPattern().numberOfValidTrackerHits();
	       PiSoftDau_dxyPVSignif[nK3Pi] = abs(psoftDauTrack->dxy(primVertex.position())) / psoftDauTrack->dxyError();
	     }
	   //
	  }
       //
       // if there is no D* associated to the D0:
       if ( dstar_exist == false)
	 {
	   //cout << "NO ASSOCIATED D* for this D0" << endl;
	   K3PiPiSoft_pseudoCharge[nK3Pi] = -3;
	   K3PiPiSoft_Mass[nK3Pi] = 0.;
	   K3PiPiSoft_pseudoPt[nK3Pi] = 0.;
           K3PiPiSoft_pseudoEta[nK3Pi] = -9.;
	   K3PiPiSoft_deltaMass[nK3Pi] = 0.;
	   //
	   PiSoftDau_Pt[nK3Pi] = 0.;
	   PiSoftDau_Eta[nK3Pi] = -9.;
	   PiSoftDau_Phi[nK3Pi] = -9.;
	   PiSoftDau_NormChi2[nK3Pi] = -1.;
	   PiSoftDau_Charge[nK3Pi] = -3;
	   PiSoftDau_hasValHitPxB[nK3Pi] = 0;
	   PiSoftDau_hasValHitPxE[nK3Pi] = 0;
	   PiSoftDau_NTrackerHits[nK3Pi] = 0;
	   PiSoftDau_dxyPVSignif[nK3Pi] = -9.;
	 }
       /////////////////////////////////////////////////////////////////////
       //
       TransientVertex theRecoVertex = theKalmanFitter.vertex(transTracks);
       //
       if ( theRecoVertex.isValid() && theRecoVertex.totalChiSquared() >= 0. ) 
         {
         // Create reco::Vertex object for use in creating the Candidate
         reco::Vertex theVtx = theRecoVertex;
         //
         //GlobalPoint vtxPos(theVtx.x(), theVtx.y(), theVtx.z());
         //cout << "x-vtx = " << theVtx.x() << " , y-vtx = " << theVtx.y() << " , z-vtx = " << theVtx.z() << endl;
         //cout << "Normalized chi2 = " << theVtx.normalizedChi2() << endl;
         //const Vertex::CovarianceMatrix vtxCov(theVtx.covariance());
         //
         double vtxChi2Prob = TMath::Prob(theVtx.chi2(),(int)theVtx.ndof());
         //
         K3Pi_VtxChi2[nK3Pi] = theVtx.normalizedChi2();
         K3Pi_VtxChi2Prob[nK3Pi] = vtxChi2Prob;
         //
         //////////////////////////////////////////////////////
         // 
         double xydist = sqrt(((PVx - theVtx.x())*(PVx - theVtx.x())) + ((PVy - theVtx.y())*(PVy - theVtx.y())));
         // this is just for check wrt vtxDistXY 
         K3Pi_VtxCalcDistXY[nK3Pi] = xydist;
         //
         ////////////////////////////////////////////////////
         //
         //somma scorrelata delle matrici di covarianza
         SMatrixSym3D totalCov = primVertexCov + theVtx.covariance();
         //
         //calculate transverse distance wrt PV:
         SVector3 distanceVector( theVtx.x() - primVertex.position().x(), theVtx.y() - primVertex.position().y(), 0. );
         //
         double vtxDistXY = ROOT::Math::Mag(distanceVector);
         //
         double vtxDistXYError = sqrt(ROOT::Math::Similarity(totalCov, distanceVector)) / vtxDistXY ;
         //
         K3Pi_VtxDistXY[nK3Pi] = vtxDistXY;
         K3Pi_VtxDistXYError[nK3Pi] = vtxDistXYError;
         K3Pi_VtxDistXYSignif[nK3Pi] = vtxDistXY / vtxDistXYError;
         //
         //-- pointing angle
         //
         // calculate 3D distance
         math::XYZVector pseudoMomVector( myD0->p4().Px(), myD0->p4().Py(), myD0->p4().Pz() );
         math::XYZVector distance3DVector( theVtx.x() - primVertex.position().x(), theVtx.y() - primVertex.position().y(), theVtx.z() - primVertex.position().z() );
         //
         K3Pi_pointingAngle[nK3Pi] = ROOT::Math::VectorUtil::Angle(pseudoMomVector,distance3DVector);
         //
         // calculate 2D distance
         math::XYZVector pseudoTransvMomVector( myD0->p4().Px(), myD0->p4().Py(), 0. );
         math::XYZVector distance2DVector( theVtx.x() - primVertex.position().x(), theVtx.y() - primVertex.position().y(), 0. );
         //
         K3Pi_2DpointingAngle[nK3Pi] = ROOT::Math::VectorUtil::Angle(pseudoTransvMomVector,distance2DVector);
         //
         } else {
                 K3Pi_VtxChi2[nK3Pi] = -9.;
                 K3Pi_VtxChi2Prob[nK3Pi] = -9.;
                 K3Pi_VtxCalcDistXY[nK3Pi] = -9.;
                 //
                 K3Pi_VtxDistXY[nK3Pi] = -9.;
                 K3Pi_VtxDistXYError[nK3Pi] = -9.;
                 K3Pi_VtxDistXYSignif[nK3Pi] = -9.;
                 //
                 K3Pi_pointingAngle[nK3Pi] = -9.;
                 K3Pi_2DpointingAngle[nK3Pi] = -9.;
         }
       // 
       /////////////////////////////
       //
       nK3Pi++;
       //
     }
  //
  //cout << "out of the loop nK3Pi = " << nK3Pi << endl;
  //
  if( myD0cands.size()< 1000) Tree_->Fill();
  //
}


// ------------ method called once each job just before starting event loop  ------------
void 
DataD0K3Pi::beginJob( void )
{
   theFile_    = new TFile( outputFileName_.c_str(), "RECREATE" ) ;
  //
  /////////////////////////
  // initialize the tree
  Tree_ = new TTree("D0K3PiTree","D0K3Pi Tree");
  //
  Tree_->Branch("nK3Pi",&nK3Pi,"nK3Pi/i");
  //Tree_->Branch("isBSValid",&isBSValid,"isBSValid/i");
  Tree_->Branch("runId",&runId,"runId/i");
  Tree_->Branch("evtId",&evtId,"evtId/i");
  Tree_->Branch("lsId",&lsId,"lsId/i");  
  Tree_->Branch("ntracks",&ntracks,"ntracks/i");  
  Tree_->Branch("nGlbMuons",&nGlbMuons,"nGlbMuons/i");  
  //
  // MC-only
  ////Tree_->Branch("hasGenD0",&hasGenD0,"hasGenD0/b");
  ////Tree_->Branch("hasBq",&hasBq,"hasBq/b");
  //
  Tree_->Branch("PVx",&PVx,"PVx/D");   
  Tree_->Branch("PVy",&PVy,"PVy/D");   
  Tree_->Branch("PVz",&PVz,"PVz/D");   
  Tree_->Branch("PVxErr",&PVxErr,"PVxErr/D");   
  Tree_->Branch("PVyErr",&PVyErr,"PVyErr/D");   
  Tree_->Branch("PVzErr",&PVzErr,"PVzErr/D"); 
  //
  ////////// Mother
  //
  Tree_->Branch("K3Pi_Mass",&K3Pi_Mass,"K3Pi_Mass[nK3Pi]/D");
  Tree_->Branch("K3Pi_pseudoP",&K3Pi_pseudoP,"K3Pi_pseudoP[nK3Pi]/D");
  Tree_->Branch("K3Pi_pseudoPt",&K3Pi_pseudoPt,"K3Pi_pseudoPt[nK3Pi]/D");
  Tree_->Branch("K3Pi_pseudoEta",&K3Pi_pseudoEta,"K3Pi_pseudoEta[nK3Pi]/D");
  Tree_->Branch("K3Pi_pseudoPhi",&K3Pi_pseudoPhi,"K3Pi_pseudoPhi[nK3Pi]/D");
  Tree_->Branch("K3Pi_pseudoCharge",&K3Pi_pseudoCharge,"K3Pi_pseudoCharge[nK3Pi]/I");
  //
  //Tree_->Branch("K3Pi_helicityAngle",&K3Pi_helicityAngle,"K3Pi_helicityAngle[nK3Pi]/D");
  //Tree_->Branch("ptHigh_helicityAngle",&ptHigh_helicityAngle,"ptHigh_helicityAngle[nK3Pi]/D");
  //
  Tree_->Branch("K3Pi_deltaR_kp1",&K3Pi_deltaR_kp1,"K3Pi_deltaR_kp1[nK3Pi]/D");
  Tree_->Branch("K3Pi_deltaR_kp2",&K3Pi_deltaR_kp2,"K3Pi_deltaR_kp2[nK3Pi]/D");
  Tree_->Branch("K3Pi_deltaR_kp3",&K3Pi_deltaR_kp3,"K3Pi_deltaR_kp3[nK3Pi]/D");
  Tree_->Branch("K3Pi_deltaR_p1p3",&K3Pi_deltaR_p1p3,"K3Pi_deltaR_p1p3[nK3Pi]/D");
  Tree_->Branch("K3Pi_deltaR_p2p3",&K3Pi_deltaR_p2p3,"K3Pi_deltaR_p2p3[nK3Pi]/D");
  Tree_->Branch("K3Pi_deltaR_p1p2",&K3Pi_deltaR_p1p2,"K3Pi_deltaR_p1p2[nK3Pi]/D");
  //
  Tree_->Branch("K3Pi_pointingAngle",&K3Pi_pointingAngle,"K3Pi_pointingAngle[nK3Pi]/D");
  Tree_->Branch("K3Pi_2DpointingAngle",&K3Pi_2DpointingAngle,"K3Pi_2DpointingAngle[nK3Pi]/D");
  //
  //
  Tree_->Branch("K3Pi_VtxChi2",&K3Pi_VtxChi2,"K3Pi_VtxChi2[nK3Pi]/D");
  Tree_->Branch("K3Pi_VtxChi2Prob",&K3Pi_VtxChi2Prob,"K3Pi_VtxChi2Prob[nK3Pi]/D");
  Tree_->Branch("K3Pi_VtxCalcDistXY",&K3Pi_VtxCalcDistXY,"K3Pi_VtxCalcDistXY[nK3Pi]/D");
  //
  Tree_->Branch("K3Pi_VtxDistXY",&K3Pi_VtxDistXY,"K3Pi_VtxDistXY[nK3Pi]/D");
  Tree_->Branch("K3Pi_VtxDistXYError",&K3Pi_VtxDistXYError,"K3Pi_VtxDistXYError[nK3Pi]/D");
  Tree_->Branch("K3Pi_VtxDistXYSignif",&K3Pi_VtxDistXYSignif,"K3Pi_VtxDistXYSignif[nK3Pi]/D");
  //
  ///////// Daughters
  //
  Tree_->Branch("KDau_P",&KDau_P,"KDau_P[nK3Pi]/D");   
  Tree_->Branch("P1Dau_P",&P1Dau_P,"P1Dau_P[nK3Pi]/D");    
  Tree_->Branch("P2Dau_P",&P2Dau_P,"P2Dau_P[nK3Pi]/D");   
  Tree_->Branch("P3Dau_P",&P3Dau_P,"P3Dau_P[nK3Pi]/D"); 
  //
  Tree_->Branch("KDau_Pt",&KDau_Pt,"KDau_Pt[nK3Pi]/D");   
  Tree_->Branch("P1Dau_Pt",&P1Dau_Pt,"P1Dau_Pt[nK3Pi]/D");   
  Tree_->Branch("P2Dau_Pt",&P2Dau_Pt,"P2Dau_Pt[nK3Pi]/D");   
  Tree_->Branch("P3Dau_Pt",&P3Dau_Pt,"P3Dau_Pt[nK3Pi]/D"); 
  //
  Tree_->Branch("KDau_Eta",&KDau_Eta,"KDau_Eta[nK3Pi]/D");     
  Tree_->Branch("P1Dau_Eta",&P1Dau_Eta,"P1Dau_Eta[nK3Pi]/D");    
  Tree_->Branch("P2Dau_Eta",&P2Dau_Eta,"P2Dau_Eta[nK3Pi]/D");     
  Tree_->Branch("P3Dau_Eta",&P3Dau_Eta,"P3Dau_Eta[nK3Pi]/D");   
  //
  Tree_->Branch("KDau_Phi",&KDau_Phi,"KDau_Phi[nK3Pi]/D");     
  Tree_->Branch("P1Dau_Phi",&P1Dau_Phi,"P1Dau_Phi[nK3Pi]/D");     
  Tree_->Branch("P2Dau_Phi",&P2Dau_Phi,"P2Dau_Phi[nK3Pi]/D");     
  Tree_->Branch("P3Dau_Phi",&P3Dau_Phi,"P3Dau_Phi[nK3Pi]/D");    
  //
  Tree_->Branch("KDau_Charge",&KDau_Charge,"KDau_Charge[nK3Pi]/I"); 
  Tree_->Branch("P1Dau_Charge",&P1Dau_Charge,"P1Dau_Charge[nK3Pi]/I");
  Tree_->Branch("P2Dau_Charge",&P2Dau_Charge,"P2Dau_Charge[nK3Pi]/I");
  Tree_->Branch("P3Dau_Charge",&P3Dau_Charge,"P3Dau_Charge[nK3Pi]/I"); 
  //
  Tree_->Branch("KDau_NormChi2",&KDau_NormChi2,"KDau_NormChi2[nK3Pi]/D"); 
  Tree_->Branch("P1Dau_NormChi2",&P1Dau_NormChi2,"P1Dau_NormChi2[nK3Pi]/D"); 
  Tree_->Branch("P2Dau_NormChi2",&P2Dau_NormChi2,"P2Dau_NormChi2[nK3Pi]/D"); 
  Tree_->Branch("P3Dau_NormChi2",&P3Dau_NormChi2,"P3Dau_NormChi2[nK3Pi]/D"); 
  //
  Tree_->Branch("KDau_Chi2",&KDau_Chi2,"KDau_Chi2[nK3Pi]/D"); 
  Tree_->Branch("P1Dau_Chi2",&P1Dau_Chi2,"P1Dau_Chi2[nK3Pi]/D");
  Tree_->Branch("P2Dau_Chi2",&P2Dau_Chi2,"P2Dau_Chi2[nK3Pi]/D"); 
  Tree_->Branch("P3Dau_Chi2",&P3Dau_Chi2,"P3Dau_Chi2[nK3Pi]/D"); 
  //
  Tree_->Branch("KDau_Prob",&KDau_Prob,"KDau_Prob[nK3Pi]/D"); 
  Tree_->Branch("P1Dau_Prob",&P1Dau_Prob,"P1Dau_Prob[nK3Pi]/D"); 
  Tree_->Branch("P2Dau_Prob",&P2Dau_Prob,"P2Dau_Prob[nK3Pi]/D"); 
  Tree_->Branch("P3Dau_Prob",&P3Dau_Prob,"P3Dau_Prob[nK3Pi]/D"); 
  //
  //Tree_->Branch("KDau_NHits",&KDau_NHits,"KDau_NHits[nK3Pi]/I");     
  //Tree_->Branch("P1Dau_NHits",&P1Dau_NHits,"P1Dau_NHits[nK3Pi]/I");    
  //Tree_->Branch("P2Dau_NHits",&P2Dau_NHits,"P2Dau_NHits[nK3Pi]/I");      
  //Tree_->Branch("P3Dau_NHits",&P3Dau_NHits,"P3Dau_NHits[nK3Pi]/I"); 
  //
  Tree_->Branch("KDau_NValHits",&KDau_NValHits,"KDau_NValHits[nK3Pi]/I");   
  Tree_->Branch("P1Dau_NValHits",&P1Dau_NValHits,"P1Dau_NValHits[nK3Pi]/I");   
  Tree_->Branch("P2Dau_NValHits",&P2Dau_NValHits,"P2Dau_NValHits[nK3Pi]/I");  
  Tree_->Branch("P3Dau_NValHits",&P3Dau_NValHits,"P3Dau_NValHits[nK3Pi]/I"); 
  //
  //Tree_->Branch("KDau_NLostHits",&KDau_NLostHits,"KDau_NLostHits[nK3Pi]/I");  
  //Tree_->Branch("P1Dau_NLostHits",&P1Dau_NLostHits,"P1Dau_NLostHits[nK3Pi]/I");  
  //Tree_->Branch("P2Dau_NLostHits",&P2Dau_NLostHits,"P2Dau_NLostHits[nK3Pi]/I");  
  //Tree_->Branch("P3Dau_NLostHits"&P3Dau_NLostHits,"P3Dau_NLostHits[nK3Pi]/I"); 
  //
  Tree_->Branch("KDau_NTrackerHits",&KDau_NTrackerHits,"KDau_NTrackerHits[nK3Pi]/I");  
  Tree_->Branch("P1Dau_NTrackerHits",&P1Dau_NTrackerHits,"P1Dau_NTrackerHits[nK3Pi]/I");  
  Tree_->Branch("P2Dau_NTrackerHits",&P2Dau_NTrackerHits,"P2Dau_NTrackerHits[nK3Pi]/I");  
  Tree_->Branch("P3Dau_NTrackerHits",&P3Dau_NTrackerHits,"P3Dau_NTrackerHits[nK3Pi]/I");
  // 
  Tree_->Branch("KDau_NPixelHits",&KDau_NPixelHits,"KDau_NPixelHits[nK3Pi]/I");    
  Tree_->Branch("P1Dau_NPixelHits",&P1Dau_NPixelHits,"P1Dau_NPixelHits[nK3Pi]/I");    
  Tree_->Branch("P2Dau_NPixelHits",&P2Dau_NPixelHits,"P2Dau_NPixelHits[nK3Pi]/I");    
  Tree_->Branch("P3Dau_NPixelHits",&P3Dau_NPixelHits,"P3Dau_NPixelHits[nK3Pi]/I"); 
  //
  Tree_->Branch("KDau_NStripHits",&KDau_NStripHits,"KDau_NStripHits[nK3Pi]/I");  
  Tree_->Branch("P1Dau_NStripHits",&P1Dau_NStripHits,"P1Dau_NStripHits[nK3Pi]/I");    
  Tree_->Branch("P2Dau_NStripHits",&P2Dau_NStripHits,"P2Dau_NStripHits[nK3Pi]/I");    
  Tree_->Branch("P3Dau_NStripHits",&P3Dau_NStripHits,"P3Dau_NStripHits[nK3Pi]/I"); 
  //
  Tree_->Branch("KDau_NTrackerLayers",&KDau_NTrackerLayers,"KDau_NTrackerLayers[nK3Pi]/I"); 
  Tree_->Branch("P1Dau_NTrackerLayers",&P1Dau_NTrackerLayers,"P1Dau_NTrackerLayers[nK3Pi]/I"); 
  Tree_->Branch("P2Dau_NTrackerLayers",&P2Dau_NTrackerLayers,"P2Dau_NTrackerLayers[nK3Pi]/I"); 
  Tree_->Branch("P3Dau_NTrackerLayers",&P3Dau_NTrackerLayers,"P3Dau_NTrackerLayers[nK3Pi]/I"); 
  //
  Tree_->Branch("KDau_NPixelLayers",&KDau_NPixelLayers,"KDau_NPixelLayers[nK3Pi]/I"); 
  Tree_->Branch("P1Dau_NPixelLayers",&P1Dau_NPixelLayers,"P1Dau_NPixelLayers[nK3Pi]/I"); 
  Tree_->Branch("P2Dau_NPixelLayers",&P2Dau_NPixelLayers,"P2Dau_NPixelLayers[nK3Pi]/I");
  Tree_->Branch("P3Dau_NPixelLayers",&P3Dau_NPixelLayers,"P3Dau_NPixelLayers[nK3Pi]/I"); 
  //
  Tree_->Branch("KDau_NStripLayers",&KDau_NStripLayers,"KDau_NStripLayers[nK3Pi]/I"); 
  Tree_->Branch("P1Dau_NStripLayers",&P1Dau_NStripLayers,"P1Dau_NStripLayers[nK3Pi]/I"); 
  Tree_->Branch("P2Dau_NStripLayers",&P2Dau_NStripLayers,"P2Dau_NStripLayers[nK3Pi]/I"); 
  Tree_->Branch("P3Dau_NStripLayers",&P3Dau_NStripLayers,"P3Dau_NStripLayers[nK3Pi]/I"); 
  //
  Tree_->Branch("KDau_hasValHitPxB",&KDau_hasValHitPxB,"KDau_hasValHitPxB[nK3Pi]/b"); 
  Tree_->Branch("P1Dau_hasValHitPxB",&P1Dau_hasValHitPxB,"P1Dau_hasValHitPxB[nK3Pi]/b"); 
  Tree_->Branch("P2Dau_hasValHitPxB",&P2Dau_hasValHitPxB,"P2Dau_hasValHitPxB[nK3Pi]/b");
  Tree_->Branch("P3Dau_hasValHitPxB",&P3Dau_hasValHitPxB,"P3Dau_hasValHitPxB[nK3Pi]/b"); 
  //
  Tree_->Branch("KDau_hasValHitPxE",&KDau_hasValHitPxE,"KDau_hasValHitPxE[nK3Pi]/b"); 
  Tree_->Branch("P1Dau_hasValHitPxE",&P1Dau_hasValHitPxE,"P1Dau_hasValHitPxE[nK3Pi]/b"); 
  Tree_->Branch("P2Dau_hasValHitPxE",&P2Dau_hasValHitPxE,"P2Dau_hasValHitPxE[nK3Pi]/b");
  Tree_->Branch("P3Dau_hasValHitPxE",&P3Dau_hasValHitPxE,"P3Dau_hasValHitPxE[nK3Pi]/b");
  //
  Tree_->Branch("KDau_dxyPV",&KDau_dxyPV,"KDau_dxyPV[nK3Pi]/D"); 
  Tree_->Branch("P1Dau_dxyPV",&P1Dau_dxyPV,"P1Dau_dxyPV[nK3Pi]/D");
  Tree_->Branch("P2Dau_dxyPV",&P2Dau_dxyPV,"P2Dau_dxyPV[nK3Pi]/D");
  Tree_->Branch("P3Dau_dxyPV",&P3Dau_dxyPV,"P3Dau_dxyPV[nK3Pi]/D");
  //
  Tree_->Branch("KDau_dxyPVErr",&KDau_dxyPVErr,"KDau_dxyPVErr[nK3Pi]/D"); 
  Tree_->Branch("P1Dau_dxyPVErr",&P1Dau_dxyPVErr,"P1Dau_dxyPVErr[nK3Pi]/D"); 
  Tree_->Branch("P2Dau_dxyPVErr",&P2Dau_dxyPVErr,"P2Dau_dxyPVErr[nK3Pi]/D"); 
  Tree_->Branch("P3Dau_dxyPVErr",&P3Dau_dxyPVErr,"P3Dau_dxyPVErr[nK3Pi]/D"); 
  //
  Tree_->Branch("KDau_dzPV",&KDau_dzPV,"KDau_dzPV[nK3Pi]/D");    
  Tree_->Branch("P1Dau_dzPV",&P1Dau_dzPV,"P1Dau_dzPV[nK3Pi]/D");     
  Tree_->Branch("P2Dau_dzPV",&P2Dau_dzPV,"P2Dau_dzPV[nK3Pi]/D");    
  Tree_->Branch("P3Dau_dzPV",&P3Dau_dzPV,"P3Dau_dzPV[nK3Pi]/D");  
  //
  Tree_->Branch("KDau_dzPVErr",&KDau_dzPVErr,"KDau_dzPVErr[nK3Pi]/D"); 
  Tree_->Branch("P1Dau_dzPVErr",&P1Dau_dzPVErr,"P1Dau_dzPVErr[nK3Pi]/D");  
  Tree_->Branch("P2Dau_dzPVErr",&P2Dau_dzPVErr,"P2Dau_dzPVErr[nK3Pi]/D");  
  Tree_->Branch("P3Dau_dzPVErr",&P3Dau_dzPVErr,"P3Dau_dzPVErr[nK3Pi]/D"); 
  //
  Tree_->Branch("KDau_dszPV",&KDau_dszPV,"KDau_dszPV[nK3Pi]/D"); 
  Tree_->Branch("P1Dau_dszPV",&P1Dau_dszPV,"P1Dau_dszPV[nK3Pi]/D"); 
  Tree_->Branch("P2Dau_dszPV",&P2Dau_dszPV,"P2Dau_dszPV[nK3Pi]/D"); 
  Tree_->Branch("P3Dau_dszPV",&P3Dau_dszPV,"P3Dau_dszPV[nK3Pi]/D"); 
  // 
  Tree_->Branch("KDau_dszPVErr",&KDau_dszPVErr,"KDau_dszPVErr[nK3Pi]/D"); 
  Tree_->Branch("P1Dau_dszPVErr",&P1Dau_dszPVErr,"P1Dau_dszPVErr[nK3Pi]/D"); 
  Tree_->Branch("P2Dau_dszPVErr",&P2Dau_dszPVErr,"P2Dau_dszPVErr[nK3Pi]/D"); 
  Tree_->Branch("P3Dau_dszPVErr",&P3Dau_dszPVErr,"P3Dau_dszPVErr[nK3Pi]/D"); 
  //
  Tree_->Branch("KDau_dxyPVSignif",&KDau_dxyPVSignif,"KDau_dxyPVSignif[nK3Pi]/D"); 
  Tree_->Branch("P1Dau_dxyPVSignif",&P1Dau_dxyPVSignif,"P1Dau_dxyPVSignif[nK3Pi]/D");
  Tree_->Branch("P2Dau_dxyPVSignif",&P2Dau_dxyPVSignif,"P2Dau_dxyPVSignif[nK3Pi]/D"); 
  Tree_->Branch("P3Dau_dxyPVSignif",&P3Dau_dxyPVSignif,"P3Dau_dxyPVSignif[nK3Pi]/D"); 
  //
  Tree_->Branch("KDau_dzPVSignif",&KDau_dzPVSignif,"KDau_dzPVSignif[nK3Pi]/D"); 
  Tree_->Branch("P1Dau_dzPVSignif",&P1Dau_dzPVSignif,"P1Dau_dzPVSignif[nK3Pi]/D"); 
  Tree_->Branch("P2Dau_dzPVSignif",&P2Dau_dzPVSignif,"P2Dau_dzPVSignif[nK3Pi]/D");
  Tree_->Branch("P3Dau_dzPVSignif",&P3Dau_dzPVSignif,"P3Dau_dzPVSignif[nK3Pi]/D"); 
  //
  Tree_->Branch("KDau_dszPVSignif",&KDau_dszPVSignif,"KDau_dszPVSignif[nK3Pi]/D"); 
  Tree_->Branch("P1Dau_dszPVSignif",&P1Dau_dszPVSignif,"P1Dau_dszPVSignif[nK3Pi]/D");
  Tree_->Branch("P2Dau_dszPVSignif",&P2Dau_dszPVSignif,"P2Dau_dszPVSignif[nK3Pi]/D");
  Tree_->Branch("P3Dau_dszPVSignif",&P3Dau_dszPVSignif,"P3Dau_dszPVSignif[nK3Pi]/D"); 
  //
  ///////
  //
  Tree_->Branch("K3PiPiSoft_Mass",&K3PiPiSoft_Mass,"K3PiPiSoft_Mass[nK3Pi]/D");
  Tree_->Branch("K3PiPiSoft_pseudoPt",&K3PiPiSoft_pseudoPt,"K3PiPiSoft_pseudoPt[nK3Pi]/D");
  Tree_->Branch("K3PiPiSoft_pseudoEta",&K3PiPiSoft_pseudoEta,"K3PiPiSoft_pseudoEta[nK3Pi]/D");
  Tree_->Branch("K3PiPiSoft_pseudoCharge",&K3PiPiSoft_pseudoCharge,"K3PiPiSoft_pseudoCharge[nK3Pi]/I");
  //
  Tree_->Branch("K3PiPiSoft_deltaMass",&K3PiPiSoft_deltaMass,"K3PiPiSoft_deltaMass[nK3Pi]/D");
  //
  Tree_->Branch("PiSoftDau_Pt",&PiSoftDau_Pt,"PiSoftDau_Pt[nK3Pi]/D");
  Tree_->Branch("PiSoftDau_Eta",&PiSoftDau_Eta,"PiSoftDau_Eta[nK3Pi]/D");
  Tree_->Branch("PiSoftDau_Phi",&PiSoftDau_Phi,"PiSoftDau_Phi[nK3Pi]/D");
  Tree_->Branch("PiSoftDau_NormChi2",&PiSoftDau_NormChi2,"PiSoftDau_NormChi2[nK3Pi]/D");
  Tree_->Branch("PiSoftDau_Charge",&PiSoftDau_Charge,"PiSoftDau_Charge[nK3Pi]/I");
  Tree_->Branch("PiSoftDau_hasValHitPxB",&PiSoftDau_hasValHitPxB,"PiSoftDau_hasValHitPxB[nK3Pi]/b");
  Tree_->Branch("PiSoftDau_hasValHitPxE",&PiSoftDau_hasValHitPxE,"PiSoftDau_hasValHitPxE[nK3Pi]/b");
  Tree_->Branch("PiSoftDau_NTrackerHits",&PiSoftDau_NTrackerHits,"PiSoftDau_NTrackerHits[nK3Pi]/I");
  Tree_->Branch("PiSoftDau_dxyPVSignif",&PiSoftDau_dxyPVSignif,"PiSoftDau_dxyPVSignif[nK3Pi]/D");
  //
  /////////////////////////
  //
  //theFile_ ->cd();
  //
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DataD0K3Pi::endJob() {
  //
  //theFile_->Close();
}

#include "FWCore/Framework/interface/MakerMacros.h"
//define this as a plug-in
DEFINE_FWK_MODULE(DataD0K3Pi);
