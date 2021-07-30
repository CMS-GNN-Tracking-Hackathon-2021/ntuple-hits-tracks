import FWCore.ParameterSet.Config as cms 

trackingParticleSelector = cms.EDProducer(
    "TrackingParticleSelector",
    trackingParticles = cms.InputTag("mix","MergedTrackTruth"),
    prunedGenParticles = cms.InputTag("prunedGenParticles"),
)
