import FWCore.ParameterSet.Config as cms

process = cms.Process('NTUPLE')

process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True)
	#,SkipEvent = cms.untracked.vstring('ProductNotFound')
)
# import of standard configurations
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.suppressInfo = cms.untracked.vstring( "mkcands" )
process.MessageLogger.suppressWarning = cms.untracked.vstring( "mkcands" )
process.MessageLogger.cerr.FwkReport.reportEvery = 1
MC = False

# Input source
process.source = cms.Source("PoolSource",
                            skipEvents = cms.untracked.uint32( 0 ), #with 11976 Processing run: 201707 lumi: 281 event: 383901681
                            fileNames = cms.untracked.vstring(
	#'file:/gpfs_data/local/cms/store/data/Run2012A/MuOnia/AOD/13Jul2012-v1/00000/FCF2DA6C-52CF-E111-99D1-003048678FE0.root'
	#'file:/gpfs_data/local/cms/store/data/Run2012D/MuOnia/AOD/16Jan2013-v1/10000/FCD0323A-B76B-E211-9EB6-003048FFD76E.root'
        #'/store/data/Run2012A/MuOnia/AOD/22Jan2013-v1/30000/000D2FF5-EE82-E211-BEBA-0026189438A5.root',
        #'/store/data/Run2012A/MuOnia/AOD/22Jan2013-v1/30000/000F9808-D483-E211-975C-003048FFD7BE.root',
        #'/store/data/Run2012A/MuOnia/AOD/22Jan2013-v1/30000/004800B2-CF83-E211-A917-003048FFCB74.root',
        #'/store/data/Run2012A/MuOnia/AOD/22Jan2013-v1/30000/00B2FD6B-CC83-E211-B74E-0025905938A8.root',
        #'/store/data/Run2012A/MuOnia/AOD/22Jan2013-v1/30000/00B38A17-D983-E211-BA46-003048FFD79C.root'
	#
	#'/store/data/Run2012C/MuOniaParked/AOD/22Jan2013-v1/30000/1E71D761-D870-E211-9343-00215E25A5E2.root'
        #'/store/data/Run2011A/MuOnia/RECO/05Aug2011-v1/0001/989142ED-ADC0-E011-A937-002618943810.root'
        #'/store/data/Run2011A/MuOnia/AOD/PromptReco-v4/000/165/098/D8E6DF38-DA7F-E011-9A5B-003048F118DE.root'
        #
        #'file:/lustre/cms/store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/194/075/8C70102D-8A9E-E111-BE2F-BCAEC532970D.root',
        #'file:/lustre/cms/store/data/Run2011A/MuOnia/RECO/05Aug2011-v1/0000/0039AA65-CDC0-E011-8FDC-0026189438E7.root'
        # RECO
        # MuOnia
        #'file:/lustre/cms/store/data/Run2012D/MuOnia/RECO/16Jan2013-v1/10000/0A3A1C33-276B-E211-B137-00248C55CC40.root',
        #'file:/lustre/cms/store/data/Run2012D/MuOnia/RECO/16Jan2013-v1/10000/0CF96EB0-276B-E211-8376-002618FDA21D.root',
        #'file:/lustre/cms/store/data/Run2012D/MuOnia/RECO/16Jan2013-v1/10000/2A2AF16E-516B-E211-AC81-0025905964A6.root',
        # DoubleMu
        #'file:/lustre/cms/store/data/Run2012D/DoubleMu/RECO/16Jan2013-v2/10000/24DC3AB8-486C-E211-8AE5-E0CB4E1A1188.root',
        #'file:/lustre/cms/store/data/Run2012D/DoubleMu/RECO/16Jan2013-v2/10000/40D891C0-AE6B-E211-A1C2-485B39897219.root',
        #'file:/lustre/cms/store/data/Run2012D/DoubleMu/RECO/16Jan2013-v2/10000/6C84B066-C66B-E211-A31F-485B39800C34.root',
        # DoubleMuParked
        #'file:/lustre/cms/store/data/Run2012D/DoubleMuParked/RECO/22Jan2013-v1/20001/2AA3F254-0282-E211-AC6F-001E4F3F1E86.root',
        #'file:/lustre/cms/store/data/Run2012D/DoubleMuParked/RECO/22Jan2013-v1/20001/2AAC14FC-0782-E211-9425-00259073E50A.root',
        #'file:/lustre/cms/store/data/Run2012D/DoubleMuParked/RECO/22Jan2013-v1/20001/2EF32F14-1B82-E211-882C-00259074AE6A.root'
        # Sanjay
        #'file:PYTHIA6_Bd2Psi2SKpi_TuneZ2star_8TeV_cff_py_RAW2DIGI_L1Reco_RECO.root'
        # CRAB
        #'file:/lustre/cms/store/group/cristella/Bd2Psi2SKpi-PHSP/MC_generation/141028_153606/merge/MC_Bd2Psi2SKpi_first111.root'
        #'file:/lustre/cms/store/group/cristella/Bd2Psi2SKpi-PHSP/MC_generation/141028_153606/merge/MC_Bd2Psi2SKpi_1of2.root'
        'file:/lustre/cms/store/group/cristella/Bd2Psi2SKpi-PHSP/MC_generation/141028_153606/merge/MC_Bd2Psi2SKpi.root'
        ) 
)

process.source.inputCommands = cms.untracked.vstring(
        "keep *",
        "drop L1GlobalTriggerObjectMapRecord_hltL1GtObjectMap__RECO",
        "drop *_MEtoEDMConverter_*_*"
	)

process.maxEvents = cms.untracked.PSet(
	#input = cms.untracked.int32( 100000 ) # = ?Mb for MC_Bd2Psi2SKpi.root in ?' #
	input = cms.untracked.int32( -1 ) # = 356Mb for 5M events MC_Bd2Psi2SKpi.root in 5h30' #
	#input = cms.untracked.int32( 24000 ) # = 870Kb # timeout after 24500 for Run2012A/MuOnia
	#input = cms.untracked.int32( -1 ) # = 5718Kb # timeout after 3700 for Run2012A/MuOnia
	)

#Output size of CRAB jobs ~200MB usually works well. (max 300-500 Mb according to Cesare) 

process.load('Configuration.Geometry.GeometryIdeal_cff') # 53x

process.load("Configuration.StandardSequences.GeometryExtended_cff") # from Lucia
process.load("Configuration.StandardSequences.Reconstruction_cff") # from Lucia

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = 'FT_53_V6_AN3::All'
process.GlobalTag.globaltag = 'START53_V19F::All'
#process.GlobalTag.globaltag = 'START53_V7C::All'


process.load('Configuration/EventContent/EventContent_cff')
#
#  Load common sequences
#
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskAlgoTrigConfig_cff')
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')

####################################################################################
##################################good collisions############################################
    
#### 44x
#process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
#                                                      vertexCollection = cms.InputTag('offlinePrimaryVertices'),
#                                                      minimumNDOF = cms.uint32(4) ,
#                                                      maxAbsZ = cms.double(24),
#                                                      maxd0 = cms.double(2)
#                                           )

## 53x                                    
pvSelection = cms.PSet(
        minNdof = cms.double( 4. )
        , maxZ    = cms.double( 24. )
        , maxRho  = cms.double( 2. )
)

process.goodOfflinePrimaryVertices = cms.EDFilter("PrimaryVertexObjectFilter", # checks for fake PVs automatically
                                                  filterParams = pvSelection,
                                                  filter       = cms.bool( False ), # use only as producer
                                                  src          = cms.InputTag( 'offlinePrimaryVertices' )
                                          )

process.primaryVertexFilter = process.goodOfflinePrimaryVertices.clone( filter = True )


process.noscraping = cms.EDFilter("FilterOutScraping",
                                  applyfilter = cms.untracked.bool(True),
                                  debugOn = cms.untracked.bool(False),
                                  #debugOn = cms.untracked.bool(True),
                                  numtrack = cms.untracked.uint32(10),
                                  thresh = cms.untracked.double(0.25)
                          )


# PAT Layer 0+1
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.load("PhysicsTools.PatAlgos.cleaningLayer1.genericTrackCleaner_cfi")
process.cleanPatTracks.checkOverlaps.muons.requireNoOverlaps = cms.bool(False)
process.cleanPatTracks.checkOverlaps.electrons.requireNoOverlaps = cms.bool(False)
from PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi import *
patMuons.embedTrack = cms.bool(True)
patMuons.embedPickyMuon = cms.bool(False)
patMuons.embedTpfmsMuon = cms.bool(False)

# Prune generated particles to muons and their parents
process.genMuons = cms.EDProducer("GenParticlePruner",
                                  src = cms.InputTag("genParticles"),
                                  select = cms.vstring(
                                          "drop  *  ",                     # this is the default
                                          "++keep abs(pdgId) = 13",        # keep muons and their parents
                                          "drop pdgId == 21 && status = 2" # remove intermediate qcd spam carrying no flavour info
                                  )
 )



process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")
from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import  addMCinfo, useExistingPATMuons, useL1MatchingWindowForSinglets, changeTriggerProcessName, switchOffAmbiguityResolution, addDiMuonTriggers
    # with some customization
if MC:
        addMCinfo(process)
        # since we match inner tracks, keep the matching tight and make it one-to-one
        process.muonMatch.maxDeltaR = 0.05
        process.muonMatch.resolveByMatchQuality = True

addDiMuonTriggers(process)
useExistingPATMuons(process,'cleanPatMuons',addL1Info=False)
changeTriggerProcessName(process, 'HLT')
switchOffAmbiguityResolution(process) # Switch off ambiguity resolution: allow multiple reco muons to match to the same trigger muon
useL1MatchingWindowForSinglets(process)

process.muonL1Info.maxDeltaR     = 0.3
process.muonL1Info.fallbackToME1 = True
process.muonMatchHLTL1.maxDeltaR     = 0.3
process.muonMatchHLTL1.fallbackToME1 = True
process.muonMatchHLTL2.maxDeltaR = 0.3
process.muonMatchHLTL2.maxDPtRel = 10.0
process.muonMatchHLTL3.maxDeltaR = 0.1
process.muonMatchHLTL3.maxDPtRel = 10.0
process.muonMatchHLTCtfTrack.maxDeltaR = 0.1
process.muonMatchHLTCtfTrack.maxDPtRel = 10.0
process.muonMatchHLTTrackMu.maxDeltaR = 0.1
process.muonMatchHLTTrackMu.maxDPtRel = 10.0

from PhysicsTools.PatAlgos.tools.trackTools import *
######## adding tracks refitted with different mass
from RecoTracker.TrackProducer.TrackRefitters_cff import *
from TrackingTools.MaterialEffects.RungeKuttaTrackerPropagator_cfi import *
#process.RungeKuttaTrackerPropagatorForMuons = TrackingTools.MaterialEffects.RungeKuttaTrackerPropagator_cfi.RungeKuttaTrackerPropagator.clone( Mass = cms.double(0.10565837), ComponentName = cms.string('RungeKuttaTrackerPropagatorForMuons') )
#process.refittedGeneralTracksMuon = RecoTracker.TrackProducer.TrackRefitter_cfi.TrackRefitter.clone( Propagator = "RungeKuttaTrackerPropagatorForMuons" )
process.RungeKuttaTrackerPropagatorForPions = TrackingTools.MaterialEffects.RungeKuttaTrackerPropagator_cfi.RungeKuttaTrackerPropagator.clone( Mass = cms.double(0.13957), ComponentName = cms.string('RungeKuttaTrackerPropagatorForPions') )
process.refittedGeneralTracksPion = RecoTracker.TrackProducer.TrackRefitter_cfi.TrackRefitter.clone( Propagator = "RungeKuttaTrackerPropagatorForPions" )
makeTrackCandidates( process,                                # patAODTrackCands
                     label = 'TrackCands',                   # output collection will be 'allLayer0TrackCands', 'allLayer1TrackCands', 'selectedLayer1TrackCands'
                     tracks = cms.InputTag('generalTracks'), # input track collection
                     #tracks = cms.InputTag('refittedGeneralTracksMuon'), # input track collection               // AP changed from generalTracks
                     #tracks = cms.InputTag('refittedGeneralTracksPion'), # input track collection               // AP changed from generalTracks
                     #particleType = 'mu+',                   # particle type (for assigning a mass)
                     particleType = 'pi+',                   # particle type (for assigning a mass)
                     preselection = 'pt > 0.35',              # preselection cut on candidates. Only methods of 'reco::Candidate' are available
                     selection = 'pt > 0.35',                 # Selection on PAT Layer 1 objects ('selectedLayer1TrackCands')
                     #selection = 'p < 0.7',                 # Selection on PAT Layer 1 objects ('selectedLayer1TrackCands')
                     isolation = {},                         # Isolations to use ('source':deltaR; set to {} for None)
                     isoDeposits = [],
                     mcAs = None                           # Replicate MC match as the one used for Muons
             );                                    # you can specify more than one collection for this

l1cands = getattr(process, 'patTrackCands')
l1cands.addGenMatch = False

######## adding tracks refitted with Kaon mass
#process.RungeKuttaTrackerPropagator.Mass = cms.double(0.493677)
process.RungeKuttaTrackerPropagatorForKaons = TrackingTools.MaterialEffects.RungeKuttaTrackerPropagator_cfi.RungeKuttaTrackerPropagator.clone(
        Mass = cms.double(0.493677), ComponentName = cms.string('RungeKuttaTrackerPropagatorForKaons') )
process.refittedGeneralTracksKaon = RecoTracker.TrackProducer.TrackRefitter_cfi.TrackRefitter.clone( Propagator = "RungeKuttaTrackerPropagatorForKaons" )
###################################################
makeTrackCandidates( process,                                        # patAODTrackCands
                     label = 'TrackKaonCands',                       # output collection will be 'allLayer0TrackCands', 'allLayer1TrackCands', 'selectedLayer1TrackCands'
                     #tracks = cms.InputTag('refittedGeneralTracksKaon'), # input track collection               // AP changed from generalTracks
                     tracks = cms.InputTag('generalTracks'), # input track collection               // AP changed from generalTracks
                     particleType = 'K+',                            # particle type (for assigning a mass)  // AP changed from pi to K
                     #particleType = 'pi+',                            # particle type (for assigning a mass)  // AP changed from pi to K
                     #particleType = 'mu+',                            # particle type (for assigning a mass)  // AP changed from pi to K
                     preselection = 'pt > 0.35',                      # preselection cut on candidates. Only methods of 'reco::Candidate' are available
                     selection = 'pt > 0.35',                         # Selection on PAT Layer 1 objects ('selectedLayer1TrackCands')
                     #selection = 'p < 0.7',                         # Selection on PAT Layer 1 objects ('selectedLayer1TrackCands')
                     isolation = {},                                 # Isolations to use ('source':deltaR; set to {} for None)
                     isoDeposits = [],
                     #mcAs = 'muon'                                   # Replicate MC match as the one used for Muons # AP "=None"  ??
                     mcAs = None                                    # Replicate MC match as the one used for Muons
             );                                                      # you can specify more than one collection for this

l1cands = getattr(process, 'patTrackKaonCands')
l1cands.addGenMatch = False

process.load("RecoTracker.DeDx.dedxHarmonic2_cfi")
process.dedxHarmonic2Kaon = RecoTracker.DeDx.dedxHarmonic2_cfi.dedxHarmonic2.clone (
        tracks = 'refittedGeneralTracksKaon',
        trajectoryTrackAssociation = 'refittedGeneralTracksKaon'
)
# dE/dx hits
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cff")
#process.load("RecoTracker.TrackProducer.TrackRefitters_cff") #already imported above
#process.TrackRefitter.src = 'generalTracks'
#process.TrackRefitter.src = 'refittedGeneralTracksPion'

#process.dedxHitInfo = cms.EDProducer("HSCPDeDxInfoProducer",
#                                     #tracks = cms.InputTag("TrackRefitter"),
#                                     #trajectoryTrackAssociation = cms.InputTag("TrackRefitter"),
#                                     tracks = cms.InputTag("refittedGeneralTracksPion"),
#                                     trajectoryTrackAssociation = cms.InputTag("refittedGeneralTracksPion"),
#                                     
#                                     UseStrip  = cms.bool(True),
#                                     UsePixel  = cms.bool(True),
#                                     MeVperADCStrip = cms.double(3.61e-06*265),
#                                     MeVperADCPixel = cms.double(3.61e-06),
#                                     
#                                     UseCalibration = cms.bool(False),
#                                     calibrationPath = cms.string("/afs/cern.ch/user/q/querten/workspace/public/dEdx/CMSSW_5_2_4/src/dEdx/ppGridProject/Gains.root"),
#                                     ShapeTest = cms.bool(True),
#                             )
#
#process.dedxHitInfoKaon = cms.EDProducer("HSCPDeDxInfoProducer",
#                                         tracks = cms.InputTag("refittedGeneralTracksKaon"),
#                                         trajectoryTrackAssociation = cms.InputTag("refittedGeneralTracksKaon"),
#                                         
#                                         UseStrip  = cms.bool(True),
#                                         UsePixel  = cms.bool(True),
#                                         MeVperADCStrip = cms.double(3.61e-06*265),
#                                         MeVperADCPixel = cms.double(3.61e-06),
#                                         
#                                         UseCalibration = cms.bool(False),
#                                         calibrationPath = cms.string("/afs/cern.ch/user/q/querten/workspace/public/dEdx/CMSSW_5_2_4/src/dEdx/ppGridProject/Gains.root"),
#                                         ShapeTest = cms.bool(True),
#                                 )


#process.PATfilter = cms.EDFilter("X3872FilterPAT")
process.PATfilter = cms.EDFilter("Z4430FilterPAT")

process.mkcands = cms.EDAnalyzer("MuMuPiKPAT",
                                 HLTriggerResults = cms.untracked.InputTag("TriggerResults","","HLT"),
                                 inputGEN  = cms.untracked.InputTag("genParticles"),
                                 VtxSample   = cms.untracked.string('offlinePrimaryVertices'),
                                 DoMuMuMassConstraint = cms.untracked.bool(True),
                                 SkipJPsi = cms.untracked.bool(True),
                                 SkipPsi2S = cms.untracked.bool(False),
                                 SameSign = cms.untracked.bool(False),
                                 DoMonteCarloTree = cms.untracked.bool(True),
                                 MonteCarloParticleId = cms.untracked.int32(20443),
                                 MinNumMuPixHits = cms.untracked.int32(1),
                                 MinNumMuSiHits = cms.untracked.int32(8),
                                 MaxMuNormChi2 = cms.untracked.double(7),
                                 MaxMuD0 = cms.untracked.double(10.0),

                                 MinJPsiMass = cms.untracked.double(2.9), 
                                 MaxJPsiMass = cms.untracked.double(3.3), 
                                 MinPsiPrimeMass = cms.untracked.double(3.55), 
                                 MaxPsiPrimeMass = cms.untracked.double(3.8), 

                                 MinNumTrSiHits = cms.untracked.int32(4),
                                 MinTrPt = cms.untracked.double(0.350),
                                 Chi2NDF_Track =  cms.untracked.double(7.0),
				 # Delta R
				 MaxMuMuTrackDR = cms.untracked.double(1.5), 
                                 MaxB0CandTrackDR = cms.untracked.double(1.5),   
                                 UseB0Dr = cms.untracked.bool(True),            

                                 MinMuMuPiKMass = cms.untracked.double(4.8),
                                 MaxMuMuPiKMass = cms.untracked.double(5.6),

                                 resolvePileUpAmbiguity = cms.untracked.bool(True),
                                 #addMuMulessPrimaryVertex = cms.untracked.bool(True),
                                 addMuMulessPrimaryVertex = cms.untracked.bool(False),
                                 addB0lessPrimaryVertex = cms.untracked.bool(True),
                                 #Debug_Output = cms.untracked.bool(False),
                                 Debug_Output = cms.untracked.bool(True),
                                 ##
                                 ##  use the correct trigger path
                                 ##
                                 TriggersForMatching = cms.untracked.vstring(
                                         # 2012 displaced J/psi = Alessandra
                                         #"HLT_DoubleMu4_Jpsi_Displaced_v9", "HLT_DoubleMu4_Jpsi_Displaced_v10", "HLT_DoubleMu4_Jpsi_Displaced_v11", "HLT_DoubleMu4_Jpsi_Displaced_v12",
                                         # Lucia
                                         # 2010
                                         #"HLT_DoubleMu3_Quarkonium_v1", "HLT_DoubleMu3_Quarkonium_v2",
                                         #"HLT_Dimuon6p5_Barrel_PsiPrime_v1",
                                         # 2011
                                         #"HLT_Dimuon7_PsiPrime_v1", "HLT_Dimuon7_PsiPrime_v2", "HLT_Dimuon7_PsiPrime_v3", "HLT_Dimuon7_PsiPrime_v4", "HLT_Dimuon7_PsiPrime_v5",
                                         #"HLT_Dimuon9_PsiPrime_v1", "HLT_Dimuon9_PsiPrime_v4", "HLT_Dimuon9_PsiPrime_v5",
                                         #"HLT_Dimuon11_PsiPrime_v1", "HLT_Dimuon11_PsiPrime_v4", "HLT_Dimuon11_PsiPrime_v5",
                                         # inclusive psi(2S)
                                         #"HLT_Dimuon0_PsiPrime_v3", "HLT_Dimuon0_PsiPrime_v4", "HLT_Dimuon0_PsiPrime_v5", "HLT_Dimuon0_PsiPrime_v6",
                                         "HLT_Dimuon5_PsiPrime_v3", "HLT_Dimuon5_PsiPrime_v4", "HLT_Dimuon5_PsiPrime_v5", "HLT_Dimuon5_PsiPrime_v6",
                                         #"HLT_Dimuon7_PsiPrime_v1", "HLT_Dimuon7_PsiPrime_v2", "HLT_Dimuon7_PsiPrime_v3", "HLT_Dimuon9_PsiPrime_v9",
                                         #"HLT_DoubleMu3p5_LowMass_Displaced_v3", "HLT_DoubleMu3p5_LowMass_Displaced_v4", "HLT_DoubleMu3p5_LowMass_Displaced_v5", "HLT_DoubleMu3p5_LowMass_Displaced_v6"
                                 ),
                                FiltersForMatching = cms.untracked.vstring(
                                        # Alessandra
                                        #"hltDisplacedmumuFilterDoubleMu4Jpsi", "hltDisplacedmumuFilterDoubleMu4Jpsi", "hltDisplacedmumuFilterDoubleMu4Jpsi", "hltDisplacedmumuFilterDoubleMu4Jpsi"
                                        # Kay
                                        "hltVertexmumuFilterDimuon5PsiPrime", "hltVertexmumuFilterDimuon5PsiPrime", "hltVertexmumuFilterDimuon5PsiPrime", "hltVertexmumuFilterDimuon5PsiPrime"#, "hltVertexmumuFilterDimuon7PsiPrime", "hltVertexmumuFilterDimuon7PsiPrime", "hltVertexmumuFilterDimuon7PsiPrime", "hltVertexmumuFilterDimuon7PsiPrime"                               
                                        #hltDoubleMu4JpsiDisplacedL3Filtered         
                                )
                                 
                                 
                         )


process.TFileService = cms.Service("TFileService",
    fileName = cms.string('MuOniaRun2012D_13Jul_MuMuPiPiPAT_ntpl.root')
)


# turn off MC matching for the process
from PhysicsTools.PatAlgos.tools.coreTools import *
# old: removeMCMatching(process, ['All'], outputInProcess = False)
removeMCMatching(process,['All'],"",None,[])

process.patDefaultSequence.remove(process.patJetCorrFactors)
process.patDefaultSequence.remove(process.patJetCharge)
process.patDefaultSequence.remove(process.patJetPartonMatch)
process.patDefaultSequence.remove(process.patJetGenJetMatch)
process.patDefaultSequence.remove(process.patJetPartons)
## error in 5_3_22, so removing it
#process.patDefaultSequence.remove(process.patJetPartonAssociation)
process.patDefaultSequence.remove(process.patJetFlavourAssociation)
process.patDefaultSequence.remove(process.patJets)
## error in 53x, so removing it
#process.patDefaultSequence.remove(process.metJESCorAK5CaloJet)
#process.patDefaultSequence.remove(process.metJESCorAK5CaloJetMuons)
process.patDefaultSequence.remove(process.patMETs)
process.patDefaultSequence.remove(process.selectedPatJets)
process.patDefaultSequence.remove(process.cleanPatJets)
process.patDefaultSequence.remove(process.countPatJets)

process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('onia2MuMuPAT.root'),
                               outputCommands = cms.untracked.vstring('drop *',
                                                #'keep *_genMuons_*_Onia2MuMuPAT',                      # generated muons and parents
                                                'keep patMuons_patMuonsWithTrigger_*_NTUPLE',    # All PAT muons including general tracks and matches to triggers
                                                              )
                       )

process.filter = cms.Sequence(
        process.goodOfflinePrimaryVertices
        + process.primaryVertexFilter
        + process.noscraping
)
#44x process.filter = cms.Sequence(process.primaryVertexFilter+process.noscraping)

process.ntup = cms.Path(
        #process.refittedGeneralTracksPion *
        #process.refittedGeneralTracksMuon *
        #process.refittedGeneralTracksKaon *
        #process.offlineBeamSpot * process.TrackRefitter * process.dedxHitInfo
        #process.dedxHarmonic2Kaon * 
        process.offlineBeamSpot #* process.dedxHitInfo
        * process.filter
        * process.patDefaultSequence
        * process.patMuonsWithTriggerSequence
        * process.PATfilter
        * process.mkcands
)

process.schedule = cms.Schedule(process.ntup)


# rsync -vut --existing test/sanjay/runMuMuPiPiPAT_data_Run2012_13Jul.py cristella@cmssusy.ba.infn.it:/cmshome/cristella/work/Z_analysis/exclusive/clean_14ott/CMSSW_5_3_22/src/UserCode/MuMuPiKPAT/test/sanjay/runMuMuPiPiPAT_data_Run2012_13Jul.py
