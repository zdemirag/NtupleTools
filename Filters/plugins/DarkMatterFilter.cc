// -*- C++ -*-
//
// Package:    NtupleTools/DarkMatterFilter
// Class:      DarkMatterFilter
// 
/**\class DarkMatterFilter DarkMatterFilter.cc NtupleTools/DarkMatterFilter/plugins/DarkMatterFilter.cc

 Description: skim for DarkMatter searches

 Implementation:
     Need to keep both signal and control regions
*/


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/METReco/interface/MET.h"

class DarkMatterFilter : public edm::EDFilter {
   public:
      explicit DarkMatterFilter(const edm::ParameterSet&);
      ~DarkMatterFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      // ----------member data ---------------------------
};

DarkMatterFilter::DarkMatterFilter(const edm::ParameterSet& iConfig)
{
}


DarkMatterFilter::~DarkMatterFilter(){}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
DarkMatterFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<edm::View<reco::MET> > pfmet;
  iEvent.getByLabel("slimmedMETs", pfmet);
  if (pfmet->front().pt()>200 ) return true;
  return false;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DarkMatterFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(DarkMatterFilter);
