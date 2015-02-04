# This file was initially auto generated with the following command:
#   cmsDriver.py step2 --filein file:EGM-Fall14DR-00001_step1.root --fileout file:EGM-Fall14DR-00001.root \
#       --mc --eventcontent FEVTSIM --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1 \
#       --datatier GEN-SIM-RAW-RECO --conditions MCRUN2_72_V3 --step RAW2DIGI,L1Reco,RECO --magField 38T_PostLS1 --no_exec

import FWCore.ParameterSet.Config as cms

process = cms.Process('RERECO')

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(5))

process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring('/store/mc/Fall14DR/DoublePhotonNoMaterial_FlatPt-0p01To100/GEN-SIM-RAW-RECO/Flat20to50bx50_NoMaterial_MCRUN2_72_V3-v1/00000/0057A7EE-008B-E411-BF09-0025905A60B8.root')
)

process.options = cms.untracked.PSet()

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.19 $'),
    annotation = cms.untracked.string('step2 nevts:1'),
    name = cms.untracked.string('Applications')
)

process.ntuplizer = cms.EDAnalyzer('PFClusterNtuplizer',
    pileupSummaryLabel = cms.InputTag('addPileupInfo'),
    genParticlesLabel  = cms.InputTag('genParticles'),
    pfClustersLabel    = cms.InputTag('particleFlowClusterECAL', '', 'RERECO'),
    recHitsEBLabel     = cms.InputTag('ecalRecHit', 'EcalRecHitsEB', 'RERECO'),
    recHitsEELabel     = cms.InputTag('ecalRecHit', 'EcalRecHitsEE', 'RERECO'),
    verticesLabel      = cms.InputTag('offlinePrimaryVertices',  '', 'RERECO')
)

process.TFileService = cms.Service('TFileService',
    fileName = cms.string('ntuple_photongun.root')
)

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'MCRUN2_72_V3', '')

process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.ntuplizer_step = cms.Path(process.ntuplizer)
process.endjob_step = cms.EndPath(process.endOfProcess)

process.schedule = cms.Schedule(
    process.raw2digi_step,
    process.L1Reco_step,
    process.reconstruction_step,
    process.ntuplizer_step,
    process.endjob_step,
)

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.postLS1Customs
from SLHCUpgradeSimulations.Configuration.postLS1Customs import customisePostLS1 

#call to customisation function customisePostLS1 imported from SLHCUpgradeSimulations.Configuration.postLS1Customs
process = customisePostLS1(process)
