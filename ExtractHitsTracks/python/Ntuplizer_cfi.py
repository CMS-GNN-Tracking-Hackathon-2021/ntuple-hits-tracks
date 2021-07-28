import FWCore.ParameterSet.Config as cms 

ntuplizer = cms.EDFilter(
    "Ntuplizer",
    verbose = cms.int32(0),
    pixelRecHits = cms.InputTag("siPixelRecHits"),
    trackerRecHits = cms.InputTag("siPhase2RecHits"),
    clustersToTP = cms.InputTag("tpClusterProducer"),
    genParticles = cms.InputTag("genParticles"),
    prunedGenParticles = cms.InputTag("prunedGenParticles"),
    tpToSimHits = cms.InputTag("simHitTPAssocProducer"),
    ctfTracks = cms.InputTag("generalTracks"),
    trackingParticles = cms.InputTag("mix","MergedTrackTruth"),
    prunedTrackingParticles = cms.InputTag("trackingParticleSelector"),
    tpToTracks = cms.string('quickTrackAssociatorByHits'),
)
