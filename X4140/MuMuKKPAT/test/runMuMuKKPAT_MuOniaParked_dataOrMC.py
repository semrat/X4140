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

#MC = False
MC = True
if MC :
        #official = False
        official = True
MCMotherId = 100443 # Y(4140)
MCMotherId = 531 # Bs0
if MCMotherId == 100443 : # Y(4140)
    MCExclusiveDecay = True
elif MCMotherId == 531 : # Bs0
    MCExclusiveDecay = False

# Input source
process.source = cms.Source("PoolSource",
                            skipEvents = cms.untracked.uint32( 0 ), #with 11976 Processing run: 201707 lumi: 281 event: 383901681
                            fileNames = cms.untracked.vstring()
)

if (not MC) :
    sourceFiles = cms.untracked.vstring( # 'root://cms-xrd-global.cern.ch/' prefix could help sometimes
	    #'root://cms-xrd-global.cern.ch//store/data/Run2012C/MuOniaParked/AOD/22Jan2013-v1/30000/1E71D761-D870-E211-9343-00215E25A5E2.root'
	    #'root://xrootd.unl.edu//store/data/Run2012C/MuOniaParked/AOD/22Jan2013-v1/20000/00109B2A-2E77-E211-893B-E41F1318165C.root'
            #'/store/data/Run2012B/MuOniaParked/AOD/22Jan2013-v1/20002/74EA8A9D-2D6C-E211-A15A-00215E21D588.root'
            #'root://cmsxrootd.hep.wisc.edu//store/data/Run2012B/MuOniaParked/AOD/22Jan2013-v1/20003/CE2C94D8-036F-E211-A034-00215E21D8EE.root'
            '/store/data/Run2012B/MuOniaParked/AOD/22Jan2013-v1/20001/EE35A843-3E69-E211-9292-00215E2226AC.root'
    )
elif MC :
        if MCMotherId == 100443 : # Y(4140)
                if (not official) :
                        sourceFiles = cms.untracked.vstring(
                                # CRAB
                                # private
                                #'file:/lustre/cms/store/group/cristella/Bd2Psi2SKpi-PHSP/MC_generation/141028_153606/merge/MC_Bd2Psi2SKpi_first111.root'
                                #'file:/lustre/cms/store/group/cristella/Bd2Psi2SKpi-PHSP/MC_generation/141028_153606/merge/MC_Bd2Psi2SKpi_1of2.root'
                                #'file:/lustre/cms/store/group/cristella/Bd2Psi2SKpi-PHSP/MC_generation/141028_153606/merge/MC_Bd2Psi2SKpi.root'
                        )
		else :
                        sourceFiles = cms.untracked.vstring(
                                # offcial MC
                                #'file:/lustre/cms/store/mc/Summer12DR53X/BdToPsi2SKPi_MSEL5_TuneZ2star_8TeV-pythia6/AODSIM/PU_RD2_START53_V19F-v1/10000/5A7A64C1-C9C3-E411-99C8-C4346BC808B8.root',
                                #'file:/lustre/cms/store/mc/Summer12DR53X/BdToPsi2SKPi_MSEL5_TuneZ2star_8TeV-pythia6/AODSIM/PU_RD2_START53_V19F-v1/10000/A8203BEC-BFC3-E411-BCE7-00266CFFC948.root'

				#'/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/50000/00AC758F-4558-E611-AF39-02163E014BCA.root',
				#'/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/50000/188EC304-1F58-E611-B1E5-FA163E7918CE.root'
				#'/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/50000/4018D9C2-2058-E611-A698-02163E00F483.root'
				#'/store/mc/Summer12DR53X/Y4700ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v1/40000/0A5CD34A-51FF-E711-AC26-FA163EB89C43.root'
				#'/store/mc/Summer12DR53X/Y4700ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v1/40000/0AEBCE6F-C5FE-E711-A0D5-FA163EF56069.root'				
				'root://cms-xrd-global.cern.ch//store/mc/Summer12DR53X/Y4700ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v1/40000/0E064110-4804-E811-939E-02163E00C230.root'
				#'/store/mc/Summer12DR53X/Y4700ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v1/40000/0EB2CFFC-6F00-E811-B945-FA163E795CFE.root'
				#'/store/mc/Summer12DR53X/Y4300ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v1/40000/02A29428-60FF-E711-890A-FA163E4207B8.root'
				#'/store/mc/Summer12DR53X/Y4300ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v1/40000/046C632E-5CFE-E711-9EA5-FA163ECA7B40.root'
				#'/store/mc/Summer12DR53X/Y4300ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v1/40000/04748CFE-6700-E811-8A6E-FA163EAE0AE5.root'
				#'/store/mc/Summer12DR53X/Y4300ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v1/40000/0499F198-81FE-E711-84A0-FA163EBB9798.root'
				#'/store/mc/Summer12DR53X/Y4300ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v1/40000/04F1380B-84FE-E711-B278-FA163E42B283.root'
                                #'root://cms-xrd-global.cern.ch//store/mc/Summer12DR53X/Y4300ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v1/40000/0855E514-51FF-E711-BFAD-FA163ED679E0.root'
                                #'root://cms-xrd-global.cern.ch//store/mc/Summer12DR53X/Y4300ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v1/40000/0A072C76-6700-E811-B6DB-FA163E3084E5.root'
				#'root://cms-xrd-global.cern.ch//store/mc/Summer12DR53X/Y4300ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v1/40000/0E4E7C64-6A00-E811-A39E-FA163E056376.root'
                                #'root://cms-xrd-global.cern.ch//store/mc/Summer12DR53X/Y4300ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v1/40000/0EEDA64B-7FFE-E711-A972-FA163E9E291A.root'
				#'root://cms-xrd-global.cern.ch//store/mc/Summer12DR53X/Y4300ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v1/40000/120BADCF-78FE-E711-95B8-FA163E2BE6D1.root'
				#'root://cms-xrd-global.cern.ch//store/mc/Summer12DR53X/Y4300ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v1/40000/12DEB4D3-6CFE-E711-B5D8-FA163E94291A.root'
				#'root://cms-xrd-global.cern.ch//store/mc/Summer12DR53X/Y4300ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v1/40000/14A88CFB-B9FE-E711-89EF-FA163E24A420.root'
           )
	elif MCMotherId == 531 : # Bs0
                sourceFiles = cms.untracked.vstring(
		
		#Bs0
		#'/store/mc/Summer12DR53X/Y4300ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v1/40000/00D14AB2-6A00-E811-B9AB-FA163ED7B479.root'            
		#'/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/50000/3C06A515-2258-E611-98CA-FA163E9AE9D0.root',
		'/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/50000/4018D9C2-2058-E611-A698-02163E00F483.root'

		)

process.PoolSource.fileNames = sourceFiles ;


process.source.inputCommands = cms.untracked.vstring(
        "keep *",
        "drop L1GlobalTriggerObjectMapRecord_hltL1GtObjectMap__RECO",
        "drop *_MEtoEDMConverter_*_*"
	)

process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32( 50 ) # 256Kb in 2' for 100 events, 1Mb in 7' for 1k events, 6Mb in 50' for 8650 events, 11Mb in 66' for 10k events, 100Mb in 14h for 150k events, 1.4Gb in 4 days for 1.2M events of official MC
        #input = cms.untracked.int32( 1000 ) # 310Kb in 3' for 1k events of private MC
        #input = cms.untracked.int32( 100 ) # = 20Mb in 2h for 15k events, 2Mb in 10' for 1k events of Run2012C/MuOniaParked/AOD/22Jan2013-v1
	#input = cms.untracked.int32( 1000 ) # = 3Mb for 6546 events, 85Kb for 100, 800kb for 1k events of BsToPsiMuMu
	#input = cms.untracked.int32( 24000 ) # = 870Kb # timeout after 24500 for Run2012A/MuOnia
	#input = cms.untracked.int32( -1 ) # = 5718Kb # timeout after 3700 for Run2012A/MuOnia
	)

#Output size of CRAB jobs ~200MB usually works well. (max 300-500 Mb according to Cesare)

process.load('Configuration.Geometry.GeometryIdeal_cff') # 53x

process.load("Configuration.StandardSequences.GeometryExtended_cff") # from Lucia
process.load("Configuration.StandardSequences.Reconstruction_cff") # from Lucia

process.load("Configuration.StandardSequences.MagneticField_cff") # for using TransientTrackBuilder
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff") # for using TransientTrackBuilder
#process.GlobalTag.globaltag = 'FT_53_V6_AN3::All' # for using TransientTrackBuilder
#process.GlobalTag.globaltag = 'START53_V19F::All' # from CMSdas my MC dataset information
#process.GlobalTag.globaltag = 'START53_V7C::All'
process.GlobalTag.globaltag = 'FT_R_53_V18::All' #Global tag for 2012B data
#process.GlobalTag.globaltag = 'FT53_V21A_AN6::All' #Global tag for 2012C data
#process.GlobalTag.globaltag = 'START53_V27::All",' # from Leonardo
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
                     #particleType = 'mu+',                   # particle type (for assigning a mass) # not working, everything is a pion
                     particleType = 'pi+',                   # particle type (for assigning a mass) # not working, everything is a pion
                     preselection = 'pt > 0.4',              # preselection cut on candidates. Only methods of 'reco::Candidate' are available
                     #selection = 'pt > 0.35',                 # Selection on PAT Layer 1 objects ('selectedLayer1TrackCands')
                     #selection = 'p > 0.5',                 # Selection on PAT Layer 1 objects ('selectedLayer1TrackCands')
                     selection = 'pt > 0.4 && p > 0.5',     # Selection on PAT Layer 1 objects ('selectedLayer1TrackCands')
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
                     particleType = 'K+',                            # particle type (for assigning a mass)  // AP changed from pi to K # not working, everything is a pion
                     #particleType = 'pi+',                            # particle type (for assigning a mass)  // AP changed from pi to K # not working, everything is a pion
                     #particleType = 'mu+',                            # particle type (for assigning a mass)  // AP changed from pi to K # not working, everything is a pion
                     preselection = 'pt > 0.4',                      # preselection cut on candidates. Only methods of 'reco::Candidate' are available
                     #selection = 'pt > 0.35',                 # Selection on PAT Layer 1 objects ('selectedLayer1TrackCands')
                     #selection = 'p > 0.5',                 # Selection on PAT Layer 1 objects ('selectedLayer1TrackCands')
                     selection = 'pt > 0.4 && p > 0.5',     # Selection on PAT Layer 1 objects ('selectedLayer1TrackCands')
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
process.PATfilter = cms.EDFilter("X4140FilterPAT")

process.mkcands = cms.EDAnalyzer("MuMuKKPAT",
                                 HLTriggerResults = cms.untracked.InputTag("TriggerResults","","HLT"),
                                 inputGEN  = cms.untracked.InputTag("genParticles"),
                                 VtxSample   = cms.untracked.string('offlinePrimaryVertices'),
                                 SameSign = cms.untracked.bool(False),
                                 DoMonteCarloTree = cms.untracked.bool( MC ),
                                 MonteCarloParticleId = cms.untracked.int32(100443), #20443
                                 MonteCarloExclusiveDecay = cms.untracked.bool( MCExclusiveDecay ),
                                 MonteCarloMotherId = cms.untracked.int32( MCMotherId ),
                                 MonteCarloDaughtersN = cms.untracked.int32( 2 ), # my daughter JPsi & Phi
				 #MonteCarloDaughtersN = cms.untracked.int32( MCDaughtersN ),
				 #MonteCarloDaughterID = cms.untracked.vuint32( MCDaughterID ),
                                 #
                                 DoMuMuMassConstraint = cms.untracked.bool(True),
                                 #SkipJPsi = cms.untracked.bool(True),
                                 SkipJPsi = cms.untracked.bool(False),
                                 SkipPsi2S = cms.untracked.bool(False), # SEMRA ask closing or not
                                 MinNumMuPixHits = cms.untracked.int32(1),
                                 MinNumMuSiHits = cms.untracked.int32(8),
                                 MaxMuNormChi2 = cms.untracked.double(7),
                                 MaxMuD0 = cms.untracked.double(10.0),
                                 sharedFraction = cms.untracked.double(0.5),

                                 MinJPsiMass = cms.untracked.double(2.8), # SEMRA changed
                                 MaxJPsiMass = cms.untracked.double(3.4), # SEMRA changed
				 MinPhiMass = cms.untracked.double (0.97), # SEMRA added
 				 MaxPhiMass = cms.untracked.double (1.07), # SEMRA added
				 MaxJPsiPhiXMass = cms.untracked.double (4.35), # SEMRA added
				 MinJPsiPhiB0Mass = cms.untracked.double (5.1), # SEMRA added
				 MaxJPsiPhiB0Mass = cms.untracked.double (5.6), # SEMRA added

                                 MinNumTrSiHits = cms.untracked.int32(4),
                                 MinTrPt = cms.untracked.double(0.350),
                                 Chi2NDF_Track =  cms.untracked.double(7.0),
				 # Delta R
				 MaxMuMuTrackDR = cms.untracked.double(1.5),
                                 MaxXCandTrackDR = cms.untracked.double(1.5),
                                 UseXDr = cms.untracked.bool(True),

                                 resolvePileUpAmbiguity = cms.untracked.bool(False),
                                 addMuMulessPrimaryVertex = cms.untracked.bool(True),
                                 #addMuMulessPrimaryVertex = cms.untracked.bool(False),
                                 addXlessPrimaryVertex = cms.untracked.bool(True),
                                 Debug_Output = cms.untracked.bool(True), # False
                                 ##
                                 ##  use the correct trigger path
                                 ##
                                 TriggersForMatching = cms.untracked.vstring(
                                         #2012 displaced J/psi = Alessandra
                                         "HLT_DoubleMu4_Jpsi_Displaced_v9",
					 "HLT_DoubleMu4_Jpsi_Displaced_v10",
					 "HLT_DoubleMu4_Jpsi_Displaced_v11", "HLT_DoubleMu4_Jpsi_Displaced_v12",
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
                                         #"HLT_Dimuon5_PsiPrime_v3", "HLT_Dimuon5_PsiPrime_v4", "HLT_Dimuon5_PsiPrime_v5", "HLT_Dimuon5_PsiPrime_v6",
                                         #"HLT_Dimuon7_PsiPrime_v1", "HLT_Dimuon7_PsiPrime_v2", "HLT_Dimuon7_PsiPrime_v3", "HLT_Dimuon9_PsiPrime_v9",
                                         #"HLT_DoubleMu3p5_LowMass_Displaced_v3", "HLT_DoubleMu3p5_LowMass_Displaced_v4", "HLT_DoubleMu3p5_LowMass_Displaced_v5", "HLT_DoubleMu3p5_LowMass_Displaced_v6"
					 # inclusive J/psi
					 "HLT_Dimuon8_Jpsi_v3",
					 "HLT_Dimuon8_Jpsi_v4",
					 "HLT_Dimuon8_Jpsi_v5", "HLT_Dimuon8_Jpsi_v6", "HLT_Dimuon8_Jpsi_v7",
                                 ),
				 FiltersForMatching = cms.untracked.vstring(
                                         # Alessandra
                                         "hltDisplacedmumuFilterDoubleMu4Jpsi", "hltDisplacedmumuFilterDoubleMu4Jpsi", "hltDisplacedmumuFilterDoubleMu4Jpsi", "hltDisplacedmumuFilterDoubleMu4Jpsi"
                                         # Kay
                                         #"hltVertexmumuFilterDimuon5PsiPrime", "hltVertexmumuFilterDimuon5PsiPrime", "hltVertexmumuFilterDimuon5PsiPrime", "hltVertexmumuFilterDimuon5PsiPrime", #"hltVertexmumuFilterDimuon7PsiPrime", "hltVertexmumuFilterDimuon7PsiPrime", "hltVertexmumuFilterDimuon7PsiPrime", "hltVertexmumuFilterDimuon7PsiPrime"
                                         #hltDoubleMu4JpsiDisplacedL3Filtered
					 # inclusive J/psi (https://espace.cern.ch/cms-quarkonia/trigger-bph/SitePages/2012-InclusiveJPsi.aspx)
					 "hltVertexmumuFilterDimuon8Jpsi", "hltVertexmumuFilterDimuon8Jpsi", "hltVertexmumuFilterDimuon8Jpsi", "hltVertexmumuFilterDimuon8Jpsi", "hltVertexmumuFilterDimuon8Jpsi",
                                )


                         )


process.TFileService = cms.Service("TFileService",
    fileName = cms.string('set_below.root')
)
if (not MC) :
    process.TFileService.fileName = cms.string('MuOniaParked_Run2012B_MuMuKKPAT_ntpl.root')
elif MC :
    if MCMotherId == 100443 :
            if (not official) :
                    process.TFileService.fileName = cms.string('BdToPsiKpi_18Mar_MuMuPiKPAT_ntpl.root')
            else :
                    process.TFileService.fileName = cms.string('official_Y4140_MuMuKKPAT_ntpl.root')
    elif MCMotherId == 531 :
        process.TFileService.fileName = cms.string('official_Bs0_MuMuKKPAT_ntpl.root')


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
                                             #'keep *_genMuons_*_Onia2MuMuPAT', # generated muons and parents
                                             'keep patMuons_patMuonsWithTrigger_*_NTUPLE', # All PAT muons including general tracks and matches to triggers
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

# rsync -vut --existing test/runMuMuKKPAT_dataOrMC.py semrat@lxplus.cern.ch:/afs/cern.ch/user/s/semrat/scratch0/CMSSW_5_3_22/src/X4140/MuMuKKPAT/test/runMuMuKKPAT_dataOrMC.py
