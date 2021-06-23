import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.setDefault('inputFiles','/store/relval/CMSSW_11_2_4/RelValTTbar_14TeV/GEN-SIM-RECO/112X_mcRun4_realistic_v7_2026D49noPU-v1/00000/b45f20a7-3655-4e29-8dda-a272c9d00c84.root')
options.setDefault('outputFile','output.root')

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(-1)
)

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
        #options.inputFiles
        'file:/afs/cern.ch/user/b/bainbrid/work/public/14-GNNTracking/RECO.root'
    )
)

process.load("test.ExtractHitsTracks.Ntuplizer_cfi")

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string(options.outputFile)
)

process.output = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('output_filtered.root'),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('ntuplizer_path'))
)

process.p = cms.Path(process.ntuplizer)
