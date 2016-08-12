from WMCore.Configuration import Configuration
import os
config = Configuration()

config.section_("General")
config.General.requestName = "MC13TeV_TTJets_powheg"
config.General.workArea = "grid"
config.General.transferOutputs=True

config.section_("JobType")
config.JobType.pluginName = "Analysis"
config.JobType.psetName = "/afs/cern.ch/work/m/mgul/public/Hto_ttbar/8_0_11_tuples/CMSSW_8_0_11/src/Tupel/Tupel/scripts/simple_run_80X_cfg.py"
config.JobType.disableAutomaticOutputCollection = False
config.JobType.pyCfgParams = ['runOnData=False']
config.JobType.inputFiles = ['Spring16_25nsV6_MC.db']

config.section_("Data")
config.Data.inputDataset = "/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISpring16MiniAODv1-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_v3_ext3-v1/MINIAODSIM"
config.Data.inputDBS = "global"
config.Data.splitting = "LumiBased"
config.Data.unitsPerJob = 1
config.Data.totalUnits = 1
config.Data.publication = False
config.Data.ignoreLocality = False
config.Data.outLFNDirBase = '/store/user/mgul/test/0e0f217/'

config.section_("Site")
config.Site.storageSite = "T2_BE_IIHE"
