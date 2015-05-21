import FWCore.ParameterSet.Config as cms

from NtupleTools.InfoProducers.jetInfoProducer_cfi import *

infoProducerSequence     = cms.Sequence(
    jetInfoProducer
)
