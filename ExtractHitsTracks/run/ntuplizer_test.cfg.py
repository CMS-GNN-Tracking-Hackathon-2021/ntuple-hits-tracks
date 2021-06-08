import FWCore.ParameterSet.Config as cms
#import ExtractHitsTracks.Ntuplizer_cff


process = cms.Process("TEST")

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
#options.setDefault('inputFiles','file:input.root')
options.setDefault('outputFile','output.root')


#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(-1)
)
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.load("test.ExtractHitsTracks.Ntuplizer_cfi")

#process.ntuples = cms.EDFilter(
#	"Ntuplizer"
#)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
            #'file:/afs/cern.ch/cms/Tutorials/TWIKI_DATA/TTJets_8TeV_53X.root'
	    'file:/eos/home-l/livage/RelValTTbar_noPU_CMSSW_11_2_4_GEN-SIM-RECO_numEvent1000.root'             
   )
                            )

#process.demo = cms.EDAnalyzer('Ntuplizer')                              

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string(options.outputFile)
    )

process.output = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('output_filtered.root'),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('ntuplizer_path'))
)


#process.filterSeq = cms.Sequence(process.Ntuplizer)
#process.ntuplizer_seq = cms.Sequence()
#process.load('ExtractHitsTracks.Ntuplizer_cfi')

process.p = cms.Path(process.ntuplizer)
