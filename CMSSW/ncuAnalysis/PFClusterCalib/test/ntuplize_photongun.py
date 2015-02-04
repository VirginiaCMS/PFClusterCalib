import FWCore.ParameterSet.Config as cms

process = cms.Process('ncuAnalysis')

process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.GlobalTag.globaltag = cms.string('MCRUN2_72_V3')

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source('PoolSource',
    fileNames = cms.untracked.vstring('/store/mc/Fall14DR/DoublePhotonNoMaterial_FlatPt-0p01To100/GEN-SIM-RAW-RECO/Flat20to50bx50_NoMaterial_MCRUN2_72_V3-v1/00000/0057A7EE-008B-E411-BF09-0025905A60B8.root')
)

process.ntuplizer = cms.EDAnalyzer('PFClusterNtuplizer',
    pileupSummaryLabel = cms.InputTag('addPileupInfo'),
    genParticlesLabel  = cms.InputTag('genParticles'),
    pfClustersLabel    = cms.InputTag('particleFlowClusterECAL'),
    recHitsEBLabel     = cms.InputTag('ecalRecHit', 'EcalRecHitsEB'),
    recHitsEELabel     = cms.InputTag('ecalRecHit', 'EcalRecHitsEE'),
    verticesLabel      = cms.InputTag('offlinePrimaryVertices')
)

process.TFileService = cms.Service('TFileService',
    fileName = cms.string('ntuple_photongun.root')
)

process.load('RecoParticleFlow.PFClusterProducer.particleFlowCluster_cff')

process.reclustering_step = cms.Path(process.pfClusteringPS * process.pfClusteringECAL)
process.ntuplizer_step = cms.Path(process.ntuplizer)
process.endjob_step    = cms.EndPath(process.endOfProcess)

process.schedule = cms.Schedule(
    process.reclustering_step,
    process.ntuplizer_step,
    process.endjob_step
)
