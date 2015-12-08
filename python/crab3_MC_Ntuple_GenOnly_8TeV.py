from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'Ntuple_MC_GenOnly_8TeV_v2'
config.General.workArea = 'crab3_MC_Ntuple_GenOnly_v2'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'btopimumu_MC_GenOnly.py'
config.JobType.allowUndistributedCMSSW = True
config.JobType.disableAutomaticOutputCollection = False
config.JobType.outputFiles = ['BToPiMuMu.root']

config.Data.inputDataset = '/PYTHIA6_Bu2MuMuPi_GenOnly_8TeV/nsahoo-crab_edm_BToPiMuMu_MC_GenOnly_8TeV_v2-a3efb7e846d3ea8e877cf211ad265e84/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5


config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.outputDatasetTag = ''

config.Site.storageSite = 'T2_IN_TIFR'

