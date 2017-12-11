#input_filename = '/store/data/Run2017B/MuOnia/MINIAOD/PromptReco-v1/000/297/723/00000/9040368C-DE5E-E711-ACFF-02163E0134FF.root'
ouput_filename = 'rootuple.root'
input_filename = '/store/data/Run2017F/MuOnia/MINIAOD/PromptReco-v1/000/306/460/00000/FABC2662-9AC8-E711-BF94-02163E019BB9.root'

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

process.load("mmkk.mmkk.slimmedMuons")

# In MiniAOD, the PATMuons are already present. We just need to run Onia2MuMu, with a selection of muons.
process.oniaSelectedMuons = cms.EDFilter('PATMuonSelector',
   src = cms.InputTag('slimmedMuonsWithTrigger'),
   cut = cms.string('muonID(\"TMOneStationTight\")'
                    ' && abs(innerTrack.dxy) < 0.3'
                    ' && abs(innerTrack.dz)  < 20.'
                    ' && innerTrack.hitPattern.trackerLayersWithMeasurement > 5'
                    ' && innerTrack.hitPattern.pixelLayersWithMeasurement > 0'
                    ' && innerTrack.quality(\"highPurity\")'
                    ' && (abs(eta) <= 1.4 && pt > 1.)'
   ),
   filter = cms.bool(True)
)



process.triggerSelection = cms.EDFilter("TriggerResultsFilter",
                                        triggerConditions = cms.vstring('HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi_v*',
                                                                        'HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi_v*'
                                                                       ),
                                        hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
                                        l1tResults = cms.InputTag( "" ),
                                        throw = cms.bool(False)
                                        )

process.Onia2MuMuMuMu = cms.EDProducer('oniaMuMuMuMuPAT',
      muons=cms.InputTag('oniaSelectedMuons'),
      primaryVertexTag=cms.InputTag('offlineSlimmedPrimaryVertices'),
      beamSpotTag=cms.InputTag('offlineBeamSpot'),
      quadmuonSelection=cms.string("4.0 < mass && mass < 6.0"),
      addMCTruth = cms.bool(False),
      higherPuritySelection=cms.string(""),
      lowerPuritySelection=cms.string(""),
      addCommonVertex=cms.bool(False),
      resolvePileUpAmbiguity=cms.bool(False),
      addMuonlessPrimaryVertex=cms.bool(False),
)

process.xCandSequence = cms.Sequence(
                   process.triggerSelection *
				   process.slimmedMuonsWithTriggerSequence *
				   process.oniaSelectedMuons *
				   process.Onia2MuMuMuMu
				   )

process.p = cms.Path(process.xCandSequence)
