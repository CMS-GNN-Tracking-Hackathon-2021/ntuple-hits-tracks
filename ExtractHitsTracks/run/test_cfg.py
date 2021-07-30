import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C9_cff import Phase2C9
process = cms.Process("TEST",Phase2C9)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(10)
)

process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T15', '')

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring('file:/eos/user/b/bainbrid/WORK/14-GNNTracking/RelValZEE_14/RECO.root'),
    skipEvents = cms.untracked.uint32(0), #6 gives a muon, 71-77 gives central pion
)

process.load("RecoLocalTracker.SiPixelRecHits.SiPixelRecHits_cfi")
process.load("RecoLocalTracker.Phase2TrackerRecHits.Phase2TrackerRecHits_cfi")
process.load("SimTracker.TrackerHitAssociation.tpClusterProducerDefault_cfi")
# default: https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/PatAlgos/python/slimming/prunedGenParticles_cfi.py

process.load("PhysicsTools.PatAlgos.slimming.prunedGenParticles_cfi")
process.prunedGenParticles.select = ["drop *","keep abs(pdgId) == 211 && pt > 1. && abs(eta) < 0.5",]

process.load("SimGeneral.TrackingAnalysis.simHitTPAssociation_cfi")
process.load("SimTracker.TrackAssociatorProducers.quickTrackAssociatorByHits_cfi")
process.load("SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi")
process.trackingParticleRecoTrackAsssociation.associator = 'quickTrackAssociatorByHits'
process.load("test.ExtractHitsTracks.TrackingParticleSelector_cfi")
process.load("test.ExtractHitsTracks.Ntuplizer_cfi")
process.ntuplizer.verbose = 0 # switch verbose mode on if >0
process.ntuplizer.usePrunedGenParticles = False
process.ntuplizer.activeTrackingRegions = [1,2,3] # IT only

#process.load("RecoTracker.TkTrackingRegions.candidatePointSeededTrackingRegionsFromBeamSpot_cfi")

process.p = cms.Path(process.siPixelRecHits*
                     process.siPhase2RecHits*
                     process.tpClusterProducer*
                     process.prunedGenParticles*
                     process.simHitTPAssocProducer*
                     process.quickTrackAssociatorByHits*
                     process.trackingParticleRecoTrackAsssociation*
                     process.trackingParticleSelector*
                     #process.candidatePointSeededTrackingRegionsFromBeamSpot*
                     process.ntuplizer
)

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("ntuple.root")
)
