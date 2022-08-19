
from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'DsMuMu_analysis_v2'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../DsMuMuRootupler.py'

config.Data.inputDataset = '/Charmonium/Run2018C-PromptReco-v2/MINIAOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 20
#config.Data.lumiMask = 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions16/13TeV/Cert_271036-275783_13TeV_PromptReco_Collisions16_JSON.txt'
#config.Data.runRange = '275776-275782'
config.Data.publication = False
config.Data.outputDatasetTag = 'DsMuMu_test2'

config.Data.outLFNDirBase = '/store/user/skamble/'
config.Site.storageSite = 'T3_CH_CERNBOX'
