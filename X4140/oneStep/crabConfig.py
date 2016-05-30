from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.workArea = 'crab_projects'

config.section_("JobType")
config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = 'PYTHIA6_Bd2Psi2SKpi_TuneZ2star_8TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_RAW2DIGI_L1Reco_RECO.py'

config.section_("Data")
config.Data.primaryDataset = 'Bd2Psi2SKpi-PHSP'
config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = 5000
NJOBS = 1000  # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.publishDbsUrl = 'phys03'
config.Data.publishDataName = 'MC_generation'
config.Data.outlfn = '/store/group/cristella'
config.Data.publication = True

config.section_("Site")
config.Site.storageSite = 'T2_IT_Bari'
