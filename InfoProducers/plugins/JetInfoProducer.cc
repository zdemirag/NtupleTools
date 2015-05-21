// -*- C++ -*-
//
// Package:    NtupleTools/JetInfoProducer
// Class:      JetInfoProducer
// 
/**\class JetInfoProducer JetInfoProducer.cc NtupleTools/JetInfoProducer/plugins/JetInfoProducer.cc

 Description: Additional information about jets

 Implementation:
     Add information missing in the input edm files.
*/


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/JetReco/interface/Jet.h"

class JetInfoProducer : public edm::EDProducer {
   public:
      explicit JetInfoProducer(const edm::ParameterSet&);
      ~JetInfoProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      // ----------member data ---------------------------
};
JetInfoProducer::JetInfoProducer(const edm::ParameterSet& iConfig)
{
  //register your products
  produces<edm::ValueMap<bool> >("monojetSelection");
}


JetInfoProducer::~JetInfoProducer(){}

void
JetInfoProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // define input
  edm::Handle<edm::View<reco::Jet> > ak4Jets;
  iEvent.getByLabel("slimmedJets", ak4Jets);
  
  // prepare output
  std::vector<float> monojet_selection_values;
  
  // fill information
  for ( const auto& jet : ak4Jets->ptrs() ){
    bool selected = false;
    if (jet->pt()>30) selected = true;
    monojet_selection_values.push_back(selected);
  }

  // finalize
  std::auto_ptr<edm::ValueMap<bool> > monojet_selection(new edm::ValueMap<bool>());
  edm::ValueMap<bool>::Filler filler(*monojet_selection);
  filler.insert(ak4Jets, monojet_selection_values.begin(), monojet_selection_values.end());
  filler.fill();
  iEvent.put(monojet_selection,"monojetSelection");
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
JetInfoProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(JetInfoProducer);
