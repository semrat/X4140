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

process.FourOnia2MuMuPhi = cms.EDProducer('FourOnia2MuMuPAT',
        muons                       = cms.InputTag('oniaSelectedMuons'),
        primaryVertexTag            = cms.InputTag('offlineSlimmedPrimaryVertices'),
        beamSpotTag                 = cms.InputTag('offlineBeamSpot'),
        higherPuritySelection       = cms.string(""),
        lowerPuritySelection        = cms.string(""),
        dimuonSelection             = cms.string("0.9 < mass && mass < 1.2"),
        addCommonVertex             = cms.bool(True),
        addMuonlessPrimaryVertex    = cms.bool(False),
        addMCTruth                  = cms.bool(False),
        resolvePileUpAmbiguity      = cms.bool(True)
)

process.FourOnia2MuMuJPsi = cms.EDProducer('FourOnia2MuMuPAT',
        muons                       = cms.InputTag('oniaSelectedMuons'),
        primaryVertexTag            = cms.InputTag('offlineSlimmedPrimaryVertices'),
        beamSpotTag                 = cms.InputTag('offlineBeamSpot'),
        higherPuritySelection       = cms.string(""),
        lowerPuritySelection        = cms.string(""),
        dimuonSelection             = cms.string("2.9 < mass && mass < 3.3"),
        addCommonVertex             = cms.bool(True),
        addMuonlessPrimaryVertex    = cms.bool(False),
        addMCTruth                  = cms.bool(False),
        resolvePileUpAmbiguity      = cms.bool(True)
)

process.Onia2MuMuFilteredJpsi = cms.EDProducer('DiMuonFilter',
      OniaTag             = cms.InputTag("FourOnia2MuMuJPsi"),
      singlemuonSelection = cms.string(""),
      dimuonSelection     = cms.string("2.95 < mass && mass < 3.15"),
      do_trigger_match    = cms.bool(True),
      HLTFilters          = cms.vstring('hltDiMuonGlbOrTrk0zFiltered0p2v2'),
)

process.Onia2MuMuFilteredPhi = cms.EDProducer('DiMuonFilter',
      OniaTag             = cms.InputTag("FourOnia2MuMuPhi"),
      singlemuonSelection = cms.string(""),
      dimuonSelection     = cms.string("0.92 < mass && mass < 1.12"),
      do_trigger_match    = cms.bool(True),
      HLTFilters          = cms.vstring('hltDiMuonGlbOrTrk0zFiltered0p2v2'),

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

process.xProducer = cms.EDProducer('FourOniaProducer',
    phidimuons                  = cms.InputTag("Onia2MuMuFilteredPhi"),
    jpsidimuons                 = cms.InputTag("Onia2MuMuFilteredJpsi"),
    primaryVertexTag            = cms.InputTag('offlineSlimmedPrimaryVertices'),
    beamSpotTag                 = cms.InputTag('offlineBeamSpot'),
    dzmax                       = cms.double(20.0),
    triggerMatch                = cms.bool(True),
    addCommonVertex             = cms.bool(True),
    addMuonlessPrimaryVertex    = cms.bool(True),
    resolvePileUpAmbiguity      = cms.bool(True),
    quadmuonSelection           = cms.string("4.0 < mass && mass < 6.0"),
    deltaMass                   = cms.vdouble(0.0,2.0)  # trigger match is performed in Onia2MuMuFiltered
)

# process.xFitter = cms.EDProducer('FourOniaKinFit',
#                           x_cand = cms.InputTag("xProducer"),
#                           x_mass = cms.double(5.36679), # GeV   1S = 9.46030   2S = 10.02326    3S = 10.35520  J/psi=3.0969
#                           product_name = cms.string("xCand"),
#                           pdgID = cms.int32(531)
#                          )


process.triggerSelection = cms.EDFilter("TriggerResultsFilter",
                                        triggerConditions = cms.vstring('HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi_v*',
                                                                        'HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi_v*'
                                                                       ),
                                        hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
                                        l1tResults = cms.InputTag( "" ),
                                        throw = cms.bool(False)
                                        )

process.xCandSequence = cms.Sequence(
                   process.triggerSelection *
                   process.slimmedMuonsWithTriggerSequence *
				   process.oniaSelectedMuons *
                   process.FourOnia2MuMuPhi *
                   process.Onia2MuMuFilteredPhi *
                   #process.DiMuonCounterPhi *
				   process.FourOnia2MuMuJPsi *
                   process.Onia2MuMuFilteredJpsi *
                   #process.DiMuonCounterJPsi *
                   process.xProducer
                   #process.xFitter
				   )

process.rootuple = cms.EDAnalyzer('x4MuRootupler',
                          #chi_cand = cms.InputTag("chiProducer"),
			              x_cand = cms.InputTag("xProducer"),
                          #xrefit = cms.InputTag("xFitter","xCand"),
			              # refit2S  = cms.InputTag("chiFitter2S","y2S"),
			              # refit3S  = cms.InputTag("chiFitter3S","y3S"),
                          primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                          TriggerResults  = cms.InputTag("TriggerResults", "", "HLT"),
                          isMC = cms.bool(False)
                         )
process.p = cms.Path(process.xCandSequence * process.rootuple)
