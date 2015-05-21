import FWCore.ParameterSet.Config as cms

process= cms.Process("TEST")
# process.GlobalTag.globaltag = "GR_R_52_V7::All"

process.load('FWCore.MessageService.MessageLogger_cfi')

process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring('/store/relval/CMSSW_7_4_0/RelValTTbarLepton_13/MINIAODSIM/MCRUN2_74_V7_gensim_740pre7-v1/00000/603E7935-4EDD-E411-B16E-0025905A612E.root')
)

process.maxEvents = cms.untracked.PSet(
   input= cms.untracked.int32(100)
)

process.load('NtupleTools.InfoProducers.infoProducerSequence_cff')
process.load('NtupleTools.Filters.DarkMatterFilter_cfi')
process.output = cms.EDAnalyzer('MonoJetFlatTreeMaker',
    ntupleFileName = cms.string('monojet.root')
)

process.dmPath = cms.Path( process.darkMatterFilter * process.infoProducerSequence * process.output )


process.schedule = cms.Schedule(
    process.dmPath
)

# Spit out filter efficiency at the end.
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
