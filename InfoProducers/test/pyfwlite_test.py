#! /usr/bin/env python

# import ROOT in batch mode
import sys
oldargv = sys.argv[:]
sys.argv = [ '-b-' ]
import ROOT
ROOT.gROOT.SetBatch(True)
sys.argv = oldargv

# load FWLite C++ libraries
ROOT.gSystem.Load("libFWCoreFWLite.so");
ROOT.gSystem.Load("libDataFormatsFWLite.so");
ROOT.AutoLibraryLoader.enable()

#cms python data types
import FWCore.ParameterSet.Config as cms

# load FWlite python libraries
from DataFormats.FWLite import Handle, Events

jets, jetsLabel = Handle("std::vector<pat::Jet>"), "slimmedJets"

events = Events("output.root")
for iev,event in enumerate(events):
    if iev > 100: break
    event.getByLabel(jetsLabel, jets)
    print "Event %d: run %6d, lumi %4d, event %12d" % (iev,event.eventAuxiliary().run(), 
                                                       event.eventAuxiliary().luminosityBlock(),
                                                       event.eventAuxiliary().event())
    for i,jet in enumerate(jets.product()): 
        print "jet %2d: pt %4.1f" % (i, jet.pt())
