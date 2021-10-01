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
    input = cms.untracked.int32(3)
)

process.source = cms.Source(
    "PoolSource",
     #fileNames = cms.untracked.vstring('file:/afs/cern.ch/user/l/livage/CMSSW_11_2_4/src/test/ExtractHitsTracks/run/ntuple_PU200_numEvent1000.root'),
    #fileNames = cms.untracked.vstring('file:/eos/home-l/livage/RelValTTbar_noPU_CMSSW_11_2_4_GEN-SIM-RECO_numEvent1000.root'),
    fileNames = cms.untracked.vstring('/store/relval/CMSSW_11_2_4/RelValTTbar_14TeV/GEN-SIM-RECO/PU_112X_mcRun4_realistic_v7_2026D49PU200_rsb-v1/00000/06aa3f74-9954-42a4-ab0f-3c03c2a17904.root', 
'/store/relval/CMSSW_11_2_4/RelValTTbar_14TeV/GEN-SIM-RECO/PU_112X_mcRun4_realistic_v7_2026D49PU200_rsb-v1/00000/07230d61-7c8e-4241-b9ba-825c64d4c4a2.root', 
#'/store/relval/CMSSW_11_2_4/RelValTTbar_14TeV/GEN-SIM-RECO/PU_112X_mcRun4_realistic_v7_2026D49PU200_rsb-v1/00000/07ab0af6-0b1f-492e-9575-3c93bce23c59.root', 
#'/store/relval/CMSSW_11_2_4/RelValTTbar_14TeV/GEN-SIM-RECO/PU_112X_mcRun4_realistic_v7_2026D49PU200_rsb-v1/00000/07b945f1-b4ef-4d3f-878d-cbd5f484aae6.root', 
#'/store/relval/CMSSW_11_2_4/RelValTTbar_14TeV/GEN-SIM-RECO/PU_112X_mcRun4_realistic_v7_2026D49PU200_rsb-v1/00000/096cfba5-bbee-4f53-88bf-bac61e0d2d02.root', 
#'/store/relval/CMSSW_11_2_4/RelValTTbar_14TeV/GEN-SIM-RECO/PU_112X_mcRun4_realistic_v7_2026D49PU200_rsb-v1/00000/10a9bdb8-f547-439b-bd6d-e4d76251f67b.root', 
#'/store/relval/CMSSW_11_2_4/RelValTTbar_14TeV/GEN-SIM-RECO/PU_112X_mcRun4_realistic_v7_2026D49PU200_rsb-v1/00000/1120e680-e2d7-4d08-ae10-025cfaeab3a6.root', 
#'/store/relval/CMSSW_11_2_4/RelValTTbar_14TeV/GEN-SIM-RECO/PU_112X_mcRun4_realistic_v7_2026D49PU200_rsb-v1/00000/12cdcc62-062c-4382-a3dc-67d8db76096a.root', 
#'/store/relval/CMSSW_11_2_4/RelValTTbar_14TeV/GEN-SIM-RECO/PU_112X_mcRun4_realistic_v7_2026D49PU200_rsb-v1/00000/13d6e17e-3b76-455d-bc7b-61b5a66f8ed0.root', 
'/store/relval/CMSSW_11_2_4/RelValTTbar_14TeV/GEN-SIM-RECO/PU_112X_mcRun4_realistic_v7_2026D49PU200_rsb-v1/00000/13f757ea-087d-48c2-871f-f963e6ebc2f8.root'), 
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
process.ntuplizer.verbose = 3 # switch verbose mode on if >0
# added flag to say whether or not to build the graph
process.ntuplizer.buildGraph = True 
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
