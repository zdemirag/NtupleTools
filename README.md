Tools in this package are designed to be used within CMSSW framework
to produce compact "flat" ntuples and/or skim data. The package
consists of the following major parts:

* InfoProducers - code needed to add information to edm::Event that is
  not present in the input data

* Filters - tools needed to skim events. The selection may depend on
  information prepared by InfoProducers

* NtupleMakers - final "flat" ntuple producers for specific analysis
  needs. These are typically small ntuples targeted to a single
  analysis or a few analyses with similar signatures.

To simplify code we assume that the input data are in MiniAOD
format. If for some reason we need to access AOD/RECO/Full level of
information we should produce MiniAOD on the fly.

Core/ is a special area where we keep framework independent code that
can be used in other frameworks. Only dependences on ROOT and standard
C++ libraries are allowed.

-----------------
 E X A M P L E S
-----------------

- Produce extra information and skim data. Output is in EDM format
  keeping all objects.

cmsRun NtupleTools/InfoProducers/test/test_produce_and_skim.py 

- Produce MonoJet ntuples
cmsRun NtupleTools/NtupleMakers/test/test_monojet.py

Other examples

-------
 Q & A 
-------


