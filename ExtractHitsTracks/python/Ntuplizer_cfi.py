import FWCore.ParameterSet.Config as cms 

ntuplizer = cms.EDFilter("Ntuplizer", 
	verbose = cms.int32(0), 
	ctfTracks = cms.InputTag("generalTracks"),
    	beamspot = cms.InputTag("offlineBeamSpot")
	)
