import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(10)
)

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
        #'file:/eos/user/b/bainbrid/WORK/14-GNNTracking/RelValZEE_14/PREPROCESSED.root'
        'file:./output.root'
    )
)

process.load("test.ExtractHitsTracks.Ntuplizer_cfi")
process.ntuplizer.verbose = 0 # switch verbose mode on if >0

process.p = cms.Path(process.ntuplizer)

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("ntuple.root")
)
