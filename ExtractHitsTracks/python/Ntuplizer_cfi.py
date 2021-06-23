import FWCore.ParameterSet.Config as cms 

ntuplizer = cms.EDFilter(
    "Ntuplizer", 
    verbose = cms.int32(0), 
    ctfTracks = cms.InputTag("generalTracks"),
    siPixelClusters = cms.InputTag("siPixelClusters"),
    siPhase2Clusters = cms.InputTag("siPhase2Clusters"),
)
