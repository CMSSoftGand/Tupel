import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register('runOnData', True,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Run this on real data"
                 )
options.parseArguments()
process = cms.Process("S2")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load('RecoEgamma.ElectronIdentification.Identification.cutBasedElectronHLTPreselecition_Summer16_V1_cff')
from JetMETCorrections.Configuration.DefaultJEC_cff import *
from JetMETCorrections.Configuration.JetCorrectionServices_cff import *

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_2016SeptRepro_v7' if options.runOnData else '80X_mcRun2_asymptotic_2016_TrancheIV_v8')
#dataFile='file:pickevents.root'
dataFile='/store/mc/RunIISpring16MiniAODv2/ST_t-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/MINIAODSIM/premix_withHLT_80X_mcRun2_asymptotic_v14_ext1-v1/80000/02236588-E871-E611-BDA6-D8D385AE85C0.root'
#dataFile='/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/0693E0E7-97BE-E611-B32F-0CC47A78A3D8.root'
#dataFile='/store/mc/RunIISummer16MiniAODv2/QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/02F9F731-90D5-E611-9551-B083FED04D68.root'
jecLevels = ['L1FastJet', 'L2Relative', 'L3Absolute']
#jecFile='sqlite:Summer16_23Sep2016V3_MC.db'
jerFile='sqlite:jerSpring16_25nsV10_MC.db'
#jecTag='JetCorrectorParametersCollection_Summer16_23Sep2016V3_MC_AK4PFchs'
if options.runOnData :
	jecLevels = ['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual']
#	jecFile='sqlite:Summer16_23Sep2016AllV3_DATA.db'
#	jecTag='JetCorrectorParametersCollection_Summer16_23Sep2016AllV3_DATA_AK4PFchs'	
#	dataFile='/store/data/Run2016G/SingleMuon/MINIAOD/23Sep2016-v1/90000/02148252-C198-E611-9790-0CC47A6C0682.root'
	dataFile='/store/data/Run2016H/SingleMuon/MINIAOD/03Feb2017_ver2-v1/110000/0A4364ED-83EA-E611-B619-002481CFE888.root'

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(dataFile)

)
process.load('Configuration.StandardSequences.Services_cff')
process.load("JetMETCorrections.Modules.JetResolutionESProducer_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.load("CondCore.CondDB.CondDB_cfi")
'''
process.jec = cms.ESSource("PoolDBESSource",
        DBParameters = cms.PSet(
        messageLevel = cms.untracked.int32(0)
        ),
        toGet = cms.VPSet(
                cms.PSet(
                record = cms.string('JetCorrectionsRecord'),
                tag    = cms.string(jecTag),
                label  = cms.untracked.string('AK4PFchs')
                ),
        ),
        connect = cms.string(jecFile)
        )
process.es_prefer_jec = cms.ESPrefer('PoolDBESSource', 'jec')
'''
#----JER  https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyResolution
from CondCore.DBCommon.CondDBSetup_cfi import *
process.jer = cms.ESSource("PoolDBESSource",
        CondDBSetup,
        toGet = cms.VPSet(
            # Resolution
            cms.PSet(
                record = cms.string('JetResolutionRcd'),
                tag    = cms.string('JR_Spring16_25nsV10_MC_PtResolution_AK4PFchs'),
                label  = cms.untracked.string('AK4PFchs_pt')
                ),
            cms.PSet(
                record = cms.string('JetResolutionRcd'),
                tag    = cms.string('JR_Spring16_25nsV10_MC_PtResolution_AK4PFchs'),
                label  = cms.untracked.string('AK4PFchs_phi')
                ),

            # Scale factors
            cms.PSet(
                record = cms.string('JetResolutionScaleFactorRcd'),
                tag    = cms.string('JR_Spring16_25nsV10_MC_SF_AK4PFchs'),
                label  = cms.untracked.string('AK4PFchs')
                ),
            ),
        connect = cms.string(jerFile)
        )

process.es_prefer_jer = cms.ESPrefer('PoolDBESSource', 'jer')


from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection
updateJetCollection(
   process,
   jetSource = cms.InputTag('slimmedJets'),
   labelName = 'UpdatedJEC',
   jetCorrections = ('AK4PFchs', cms.vstring(jecLevels), 'None'),  # Do not forget 'L2L3Residual' on data!
   btagDiscriminators = [
    'pfJetBProbabilityBJetTags',
    'deepFlavourJetTags:probudsg'        ,
    'deepFlavourJetTags:probb'           ,
    'deepFlavourJetTags:probc'           ,
    'deepFlavourJetTags:probbb'          ,
    'deepFlavourJetTags:probcc'          ,
        ]
     )
#process.updatedPatJetsTransientCorrectedUpdatedJECBTag.addTagInfos = cms.bool(True)
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
#from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import slimmedMETsMuEGClean
runMetCorAndUncFromMiniAOD(process,
                           isData=options.runOnData,
                           metType = 'PF',
                           postfix=''
                           )

# Check interactively the good lumis
#if options.runOnData :
#	import FWCore.PythonUtilities.LumiList as LumiList
#	process.source.lumisToProcess = LumiList.LumiList(filename='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Final/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt').getVLuminosityBlockRange()	

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)
my_id_modules = [
    'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronHLTPreselecition_Summer16_V1_cff',
    'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff'
    ]
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('ntuple.root' )
)

jetsrcc="selectedUpdatedPatJetsUpdatedJEC"

trigPath="HLT"
trigFiltPath="PAT"
eventFilter="PAT"
metTag="slimmedMETs"
if options.runOnData :
	trigFiltPath="RECO"
	eventFilter="PAT"
	metTag="slimmedMETsMuEGClean"
process.tupel = cms.EDAnalyzer("Tupel",
  triggerfilters = cms.InputTag("TriggerResults","",trigFiltPath),
  EvtFilterSrc   = cms.InputTag("TriggerResults","",eventFilter ),
  HLTSrc         = cms.InputTag("TriggerResults","",trigPath), 
  photonSrc      = cms.InputTag("slimmedPhotons"),
  electronSrc    = cms.InputTag("slimmedElectrons"),
  muonSrc        = cms.InputTag("slimmedMuons"),
  #tauSrc        = cms.untracked.InputTag("slimmedPatTaus"),
  pfcandSrc	 = cms.InputTag("packedPFCandidates"),
  jetSrc         = cms.InputTag(jetsrcc),
  metSrc         = cms.InputTag("patMETsPF"),
  genSrc         = cms.InputTag("prunedGenParticles"),
  pgenSrc        =cms.InputTag("packedGenParticles"),
  gjetSrc        = cms.InputTag('slimmedGenJets'),
  muonMatch    = cms.string( 'muonTriggerMatchHLTMuons' ),
  muonMatch2    = cms.string( 'muonTriggerMatchHLTMuons2' ),
  elecMatch    = cms.string( 'elecTriggerMatchHLTElecs' ),
  cutBasedElectronID_Summer16_80X_V1_veto = cms.InputTag('egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto'),
  cutBasedElectronID_Summer16_80X_V1_loose = cms.InputTag('egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose'),
  cutBasedElectronID_Summer16_80X_V1_medium = cms.InputTag('egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium'),
  cutBasedElectronID_Summer16_80X_V1_tight = cms.InputTag('egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight'),
  cutBasedElectronHLTPreselection_Summer16_V1 = cms.InputTag('egmGsfElectronIDs:cutBasedElectronHLTPreselection-Summer16-V1'),
  channel    = cms.string( 'noselection' ),
  keepparticlecoll    = cms.bool(False),
  mSrcRho      = cms.InputTag('fixedGridRhoFastjetAll'),#arbitrary rho now
  CalojetLabel = cms.InputTag('slimmedJets'), #same collection now BB 
  metSource = cms.VInputTag(metTag,metTag),
  lheSource=cms.InputTag("source")
)

from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
#process.goodOfflinePrimaryVertices = cms.EDFilter(
#    "PrimaryVertexObjectFilter",
#    filterParams = pvSelector.clone( minNdof = cms.double(4.0), maxZ = cms.double(24.0),maxd0 = cms.double(2.0) ),
#    src=cms.InputTag('offlineSlimmedPrimaryVertices')
#    )

process.goodOfflinePrimaryVertices = cms.EDFilter(
    "PrimaryVertexObjectFilter",
    filterParams = pvSelector.clone( NPV     = cms.int32(1),minNdof = cms.double(4.0), maxZ = cms.double(24.0),maxRho = cms.double(2.0) ),
   src=cms.InputTag('offlineSlimmedPrimaryVertices')
    )

process.p = cms.Path(
#    process.patJets
# + process.patMETs
# + process.inclusiveSecondaryVertexFinderTagInfos
#+process.selectedMuons
#+process.selectedElectrons 
# +process.goodOfflinePrimaryVertices 
#+process.kt6PFJets
#process.pseudoTop
#+process.JetResolutionESProducer_AK4PFchs
#+process.JetResolutionESProducer_SF_AK4PFchs 
 process.tupel 
)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.suppressWarning = cms.untracked.vstring('ecalLaserCorrFilter','manystripclus53X','toomanystripclus53X')
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)


process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('test.root'),
#    outputCommands = cms.untracked.vstring(['drop *','keep patJets_patJets_*_*','keep *_*_*_PAT','keep recoTracks_unp*_*_*','keep recoVertexs_unp*_*_*'])
#    outputCommands = cms.untracked.vstring(['drop *'])
outputCommands = cms.untracked.vstring(['drop *', 'keep *_slimmedMETs*_*_*'])
)
#process.endpath= cms.EndPath(process.out)


#from PhysicsTools.PatAlgos.tools.trigTools import *
#switchOnTrigger( process ) # This is optional and can be omitted.

# Switch to selected PAT objects in the trigger matching
#removeCleaningFromTriggerMatching( process )
##############################
iFileName = "fileNameDump_cfg.py"
file = open(iFileName,'w')
file.write(str(process.dumpPython()))
file.close()


