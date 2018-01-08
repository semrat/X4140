#input_filename = '/store/data/Run2017B/MuOnia/MINIAOD/PromptReco-v1/000/297/723/00000/9040368C-DE5E-E711-ACFF-02163E0134FF.root'
ouput_filename = 'rootuple.root'
input_filename = 'file:FABC2662-9AC8-E711-BF94-02163E019BB9.root'

import FWCore.ParameterSet.Config as cms
process = cms.Process("Rootuple")

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '92X_dataRun2_Prompt_v4', '')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 20000

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(input_filename))
process.TFileService = cms.Service("TFileService",fileName = cms.string(ouput_filename))
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False))

process.load("mmkk.mmkk.slimmedMuonsTriggerMatcher2017_cfi")

process.oniaSelectedMuons = cms.EDFilter('PATMuonSelector',
   src = cms.InputTag('slimmedMuonsWithTrigger'),
   cut = cms.string('muonID(\"TMOneStationTight\")'
                    ' && abs(innerTrack.dxy) < 0.3'
                    ' && abs(innerTrack.dz)  < 20.'
                    ' && innerTrack.hitPattern.trackerLayersWithMeasurement > 5'
                    ' && innerTrack.hitPattern.pixelLayersWithMeasurement > 0'
                    ' && innerTrack.quality(\"highPurity\")'
                    ' && (abs(eta) <= 2.5 && pt > 0.2)'
   ),
   filter = cms.bool(True)
)

# process.FourOnia2MuMuJPsi = cms.EDProducer('FourOnia2MuMuPAT',
#         muons                       = cms.InputTag('oniaSelectedMuons'),
#         primaryVertexTag            = cms.InputTag('offlineSlimmedPrimaryVertices'),
#         beamSpotTag                 = cms.InputTag('offlineBeamSpot'),
#         higherPuritySelection       = cms.string(""),
#         lowerPuritySelection        = cms.string(""),
#         dimuonSelection             = cms.string("2.9 < mass && mass < 3.3 && charge==0 "),
#         addCommonVertex             = cms.bool(True),
#         addMuonlessPrimaryVertex    = cms.bool(False),
#         addMCTruth                  = cms.bool(False),
#         resolvePileUpAmbiguity      = cms.bool(True),
#         HLTFilters                  = cms.vstring('hltDiMuonGlbOrTrk0zFiltered0p2v2','hltMumuFilterDoubleMu2Jpsi')
# )

process.load("HeavyFlavorAnalysis.Onia2MuMu.onia2MuMuPAT_cfi")
process.onia2MuMuPAT.muons=cms.InputTag('oniaSelectedMuons')
process.onia2MuMuPAT.primaryVertexTag=cms.InputTag('offlineSlimmedPrimaryVertices')
process.onia2MuMuPAT.beamSpotTag=cms.InputTag('offlineBeamSpot')
process.onia2MuMuPAT.higherPuritySelection=cms.string("")
process.onia2MuMuPAT.lowerPuritySelection=cms.string("")
process.onia2MuMuPAT.dimuonSelection=cms.string("2.9 < mass && mass < 3.15")
process.onia2MuMuPAT.addMCTruth = cms.bool(False)

process.Onia2MuMuFilteredJpsi = cms.EDProducer('DiMuonFilter',
      OniaTag             = cms.InputTag("onia2MuMuPAT"),
      singlemuonSelection = cms.string(""),
      dimuonSelection     = cms.string("2.95 < mass && mass < 3.15 && userFloat('vProb') > 0.0 && charge == 0"),
      do_trigger_match    = cms.bool(False),
      HLTFilters          = cms.vstring('hltDiMuonGlbOrTrk0zFiltered0p2v2','hltMumuFilterDoubleMu2Jpsi')
)

process.DiMuonCounterJPsi = cms.EDFilter('CandViewCountFilter',
    src       = cms.InputTag("Onia2MuMuFilteredJpsi"),
    minNumber = cms.uint32(1),
    filter    = cms.bool(True)
)

process.DiMuonCounterPhi = cms.EDFilter('CandViewCountFilter',
    src       = cms.InputTag("Onia2MuMuFilteredPhi"),
    minNumber = cms.uint32(1),
    filter    = cms.bool(True)
)

# process.xFitter = cms.EDProducer('FourOniaKinFit',
#                           x_cand = cms.InputTag("jProducer"),
#                           x_mass = cms.double(5.36679), # GeV   1S = 9.46030   2S = 10.02326    3S = 10.35520  J/psi=3.0969
#                           product_name = cms.string("xCand"),
#                           pdgID = cms.int32(531)
#                          )


process.triggerSelection = cms.EDFilter("TriggerResultsFilter",
                                        triggerConditions = cms.vstring(
                                                                        'HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi_v*',
                                                                        "HLT_Dimuon0_Jpsi_v*","HLT_Dimuon0_Jpsi_NoVertexing_v*",
                                                                  "HLT_Dimuon25_Jpsi_v*","HLT_Dimuon20_Jpsi_Barrel_Seagulls_v*",
                                                            	      "HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi_v*","HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi_v*"
                                                                       ),
                                        hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
                                        l1tResults = cms.InputTag( "" ),
                                        throw = cms.bool(False)
                                        )

process.jCandSequence = cms.Sequence(
                   process.triggerSelection *
                   process.slimmedMuonsWithTriggerSequence *
				   process.oniaSelectedMuons *
				   process.onia2MuMuPAT *
                   process.Onia2MuMuFilteredJpsi
				   )

process.rootuple = cms.EDAnalyzer('jpsiRootupler',
			              j_cand = cms.InputTag("Onia2MuMuFilteredJpsi"),
                          primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                          TriggerResults  = cms.InputTag("TriggerResults", "", "HLT"),
                          isMC = cms.bool(False)
                         )
process.p = cms.Path(process.jCandSequence * process.rootuple)
