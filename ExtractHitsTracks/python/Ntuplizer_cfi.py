import FWCore.ParameterSet.Config as cms 

ntuplizer = cms.EDFilter(
    "Ntuplizer",
    verbose = cms.int32(0),
    pixelRecHits = cms.InputTag("siPixelRecHits"),
    trackerRecHits = cms.InputTag("siPhase2RecHits"),
    activeTrackingRegions = cms.vint32(), # e.g. cms.vint32(2,5),
    clustersToTP = cms.InputTag("tpClusterProducer"),
    genParticles = cms.InputTag("genParticles"),
    usePrunedGenParticles = cms.bool(True),
    prunedGenParticles = cms.InputTag("prunedGenParticles"),
    tpToSimHits = cms.InputTag("simHitTPAssocProducer"),
    ctfTracks = cms.InputTag("generalTracks"),
    trackingParticles = cms.InputTag("mix","MergedTrackTruth"),
    #prunedTrackingParticles = cms.InputTag("trackingParticleSelector"),
    tpToTracks = cms.string('quickTrackAssociatorByHits'),
)
