from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferOutputs = True

# Data
#config.General.requestName = 'dataset_C'
#config.General.requestName = 'dataset_C_missingLumis' # for resubmission
#config.General.requestName = 'dataset_C_missingLumis_backup'
#config.General.requestName = 'psi1and2S_dataset_D' # wrong J/psi triggers spelling
#config.General.requestName = 'psi1SokAnd2S_dataset_D'
#config.General.requestName = 'psi1SokAnd2S_dataset_C'
#config.General.requestName = 'psi1SokAnd2S_dataset_C_backup'
#config.General.requestName = 'psi1SokAnd2S_dataset_C_Pisa'
#config.General.requestName = 'psi1SokAnd2S_dataset_B_backup'
config.General.requestName = 'B0toMuMuPP_datasetB'

#config.General.requestName = 'official_MC'
#config.General.requestName = 'official_MC_withGenPart'
#config.General.requestName = 'official_MC_withGenVert'
#config.General.requestName = 'official_MC_withCosAlpha3D'
#config.General.requestName = 'official_MC_ext1_Legnaro'
#config.General.requestName = 'official_MC_withGenMuons'

config.General.transferLogs = True

config.section_('JobType')
#config.JobType.psetName = 'runMuMuPiKPAT_dataOrMC_03Mar.py'
#config.JobType.psetName = 'runMuMuPiKPAT_dataOrMC_09Apr.py' # for resubmission
config.JobType.psetName = 'runMuMuPiKPAT_dataOrMC_09Apr.py'

config.JobType.pluginName = 'Analysis'

config.section_('Data')
config.Data.inputDataset = '/MuOniaParked/Run2012B-22Jan2013-v1/AOD'
#config.Data.inputDataset = '/MuOniaParked/Run2012C-22Jan2013-v1/AOD'
#config.Data.inputDataset = '/MuOniaParked/Run2012D-22Jan2013-v1/AOD'
#config.Data.inputDataset = '/BdToPsi2SKPi_MSEL5_TuneZ2star_8TeV-pythia6/Summer12DR53X-PU_RD2_START53_V19F-v1/AODSIM'
#config.Data.inputDataset = '/BdToPsi2SKPi_MSEL5_TuneZ2star_8TeV-pythia6/Summer12DR53X-PU_RD2_START53_V19F_ext1-v1/AODSIM'

# Data
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 6 # 10: Task consists of 7687 jobs to process 76868 lumis for datasetD
#config.Data.unitsPerJob = 1 # for resubmission

# MC
#config.Data.splitting = 'FileBased'
#config.Data.unitsPerJob = 1 # 365 (= number of files) jobs created... forced choice since next one doesn't work
##config.Data.splitting = 'EventAwareLumiBased'
##config.Data.unitsPerJob = 10 # not working

config.Data.totalUnits = -1

config.Data.lumiMask = 'Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON_MuonPhys.txt'
#config.Data.lumiMask = 'crab_dataset_C/results/missingLumiSummary.json' # for resubmission

config.Data.publication = False

#config.Data.outLFNDirBase = '/store/user/lecriste/april2015/'
#config.Data.outLFNDirBase = '/store/user/lecriste/pentaQuark/'
#config.Data.outLFNDirBase = '/store/group/phys_muon/lecriste/TnP/' # for CERN
#/eos/cms/store/group/phys_muon/lecriste/TnP

#config.Data.ignoreLocality = True
#config.Data.ignoreLocality = False # for resubmission
#config.Data.ignoreLocality = True
config.Data.ignoreLocality = False

config.section_('User')

config.section_('Site')
#config.Site.blacklist = ['T0', 'T1'] # T0 blacklisted by default
config.Site.blacklist = ['T1*']
#config.Site.blacklist = ['T1*','T2_US_UCSD','T3_UK_London_QMUL','T2_US_Vanderbilt','T2_HU_Budapest','T3_US_TTU']

config.Site.storageSite = 'T2_IT_Bari'
#config.Site.storageSite = 'T2_IT_Legnaro'
#config.Site.storageSite = 'T2_CH_CERN'
#config.Site.storageSite = 'T2_IT_Pisa'


