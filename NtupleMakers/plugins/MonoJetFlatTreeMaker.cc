// -*- C++ -*-
//
// Package:    NtupleTools/MonoJetFlatTreeMaker
// Class:      MonoJetFlatTreeMaker
// 
/**\class MonoJetFlatTreeMaker MonoJetFlatTreeMaker.cc NtupleTools/MonoJetFlatTreeMaker/plugins/MonoJetFlatTreeMaker.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "NtupleTools/Core/MonoJetFlatTree.h"
#include "TFile.h"

class MonoJetFlatTreeMaker : public edm::EDAnalyzer {
   public:
      explicit MonoJetFlatTreeMaker(const edm::ParameterSet&);
      ~MonoJetFlatTreeMaker();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      MonoJetFlatTree* m_tree;
      TFile* m_file;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MonoJetFlatTreeMaker::MonoJetFlatTreeMaker(const edm::ParameterSet& iConfig)
{
  m_file = TFile::Open(iConfig.getParameter<std::string>("ntupleFileName").c_str(), "RECREATE");
  m_tree = new MonoJetFlatTree();
  m_tree->CreateTree();
  m_tree->tree_->SetDirectory(m_file);
}


MonoJetFlatTreeMaker::~MonoJetFlatTreeMaker()
{
  m_file->cd();
  m_tree->tree_->Write();
  m_file->Close(); 
}


//
// member functions
//

// ------------ method called for each event  ------------
void
MonoJetFlatTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  m_tree->InitVariables();
  m_tree->run_       = iEvent.id().run();
  m_tree->event_     = iEvent.id().event();
  m_tree->lumi_      = iEvent.luminosityBlock();

  m_tree->tree_->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
MonoJetFlatTreeMaker::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MonoJetFlatTreeMaker::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MonoJetFlatTreeMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MonoJetFlatTreeMaker);
