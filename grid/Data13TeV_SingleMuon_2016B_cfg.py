from WMCore.Configuration import Configuration
import os
config = Configuration()

config.section_("General")
config.General.requestName = "Data13TeV_SingleMuon_2016B"
config.General.workArea = "grid"
config.General.transferOutputs=True

config.section_("JobType")
config.JobType.pluginName = "Analysis"
config.JobType.psetName = "/afs/cern.ch/work/m/mgul/public/Hto_ttbar/8_0_11_tuples/CMSSW_8_0_11/src/Tupel/Tupel/scripts/simple_run_80X_cfg.py"
config.JobType.disableAutomaticOutputCollection = False
config.JobType.pyCfgParams = ['runOnData=True']
config.JobType.inputFiles = ['Spring16_25nsV6_DATA.db']

config.section_("Data")
config.Data.inputDataset = "/SingleMuon/Run2016B-PromptReco-v2/MINIAOD"
config.Data.inputDBS = "global"
config.Data.splitting = "LumiBased"
config.Data.unitsPerJob = 1
config.Data.totalUnits = 1
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-277148_13TeV_PromptReco_Collisions16_JSON.txt'
config.Data.publication = False
config.Data.ignoreLocality = False
config.Data.outLFNDirBase = '/store/user/mgul/test/0e0f217/'

config.section_("Site")
config.Site.storageSite = "T2_BE_IIHE"
