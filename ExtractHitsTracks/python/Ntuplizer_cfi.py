import FWCore.ParameterSet.Config as cms 

ntuplizer = cms.EDFilter(
	"Ntuplizer", 
	verbose = cms.int32(0), 
	tracks = cms.InputTag("generalTracks")
	)
