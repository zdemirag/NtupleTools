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

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "TLorentzVector.h"


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

      edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
      edm::EDGetTokenT<pat::MuonCollection> muonToken_;
      edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
      edm::EDGetTokenT<pat::TauCollection> tauToken_;
      edm::EDGetTokenT<pat::PhotonCollection> photonToken_;
      edm::EDGetTokenT<pat::JetCollection> jetToken_;
      edm::EDGetTokenT<pat::JetCollection> fatjetToken_;
      edm::EDGetTokenT<pat::METCollection> metToken_;
      edm::EDGetTokenT<pat::JetCollection> fatterjetToken_;
      edm::EDGetTokenT<pat::METCollection> skimmetToken_;
      edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
      bool verbose_, mc_;


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
MonoJetFlatTreeMaker::MonoJetFlatTreeMaker(const edm::ParameterSet& iConfig):
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
  tauToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
  photonToken_(consumes<pat::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
  jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
  fatjetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("fatjets"))),
  metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
  fatterjetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("fatterjets"))),
  skimmetToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("skimmets"))),
  pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands")))
{
  m_file = TFile::Open(iConfig.getParameter<std::string>("ntupleFileName").c_str(), "RECREATE");
  m_tree = new MonoJetFlatTree();
  m_tree->CreateTree();
  m_tree->tree_->SetDirectory(m_file);
  verbose_ = iConfig.getParameter<bool>("verbose_");
  mc_      = iConfig.getParameter<bool>("mc_");
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

  // Met
  edm::Handle<pat::METCollection> mets;
  iEvent.getByToken(metToken_, mets);
  const pat::MET &met = mets->front();
  if(verbose_){
      printf("MET: pt %5.1f, phi %+4.2f, sumEt (%.1f). genMET %.1f. MET with JES up/down: %.1f/%.1f\n",
             met.pt(), met.phi(), met.sumEt(),
             met.genMET()->pt(),
             met.shiftedPt(pat::MET::JetEnUp), met.shiftedPt(pat::MET::JetEnDown));
  }
  m_tree->pfmet_.SetPxPyPzE(met.px(), met.py(), met.pz(), met.energy());
  if(mc_) m_tree->genmet_.SetPxPyPzE(met.genMET()->px(), met.genMET()->py(), met.genMET()->pz(), met.genMET()->energy());
  else    m_tree->genmet_.SetPxPyPzE(0.0,0.0,0.0,0.0);
  
  // Jets
  edm::Handle<pat::JetCollection> jets;
  iEvent.getByToken(jetToken_, jets);
  int njets = 0;
  if (m_tree->ak4jet_->GetLast()!=-1) m_tree->ak4jet_->Clear();
  m_tree->ak4jet_btag_.clear(); m_tree->ak4jet_puid_.clear();
  for (const pat::Jet &j : *jets) {
      if (j.pt() < 20 ) continue;
      if(verbose_) 
          printf("jet  with pt %5.1f (raw pt %5.1f), eta %+4.2f, btag CSV %.3f, pileup mva disc %+.2f\n",
                 j.pt(), j.pt()*j.jecFactor("Uncorrected"), j.eta(), std::max(0.f,j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")), 
                 j.userFloat("pileupJetId:fullDiscriminant"));
      new ( (*(m_tree->ak4jet_))[njets] ) TLorentzVector(j.px(), j.py(), j.pz(), j.energy());
      m_tree->ak4jet_btag_.push_back(std::max(0.f,j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")));
      m_tree->ak4jet_puid_.push_back(j.userFloat("pileupJetId:fullDiscriminant"));
      njets++;
  }
     
  // Fat Jets
  edm::Handle<pat::JetCollection> fatjets;
  iEvent.getByToken(fatjetToken_, fatjets);
  int nfatjets = 0;
  m_tree->ak8jet_tau1_.clear(); m_tree->ak8jet_tau2_.clear(); m_tree->ak8jet_tau3_.clear();
  m_tree->ak8jet_hasSubjet_.clear(),  m_tree->ak8jet_softdropmass_.clear(); 
  m_tree->ak8jet_trimmedmass_.clear(); m_tree->ak8jet_prunedmass_.clear(); m_tree->ak8jet_filteredmass_.clear(); m_tree->ak8subjet_btag_.clear();
  if (m_tree->ak8jet_p4_ ->GetLast()!=-1) m_tree->ak8jet_p4_ ->Clear();
  if (m_tree->ak8_subjet_->GetLast()!=-1) m_tree->ak8_subjet_->Clear();
  for (const pat::Jet &j : *fatjets) {
      if (j.pt() < 200 ) continue;
      if(verbose_) printf("AK8j with pt %5.1f (raw pt %5.1f), eta %+4.2f, mass %5.1f ungroomed, %5.1f pruned, %5.1f trimmed, %5.1f filtered. CMS TopTagger %.1f\n",j.pt(), j.pt()*j.jecFactor("Uncorrected"), j.eta(), j.mass(), j.userFloat("ak8PFJetsCHSPrunedMass"), j.userFloat("ak8PFJetsCHSTrimmedMass"), j.userFloat("ak8PFJetsCHSFilteredMass"), j.userFloat("cmsTopTagPFJetsCHSMassAK8"));
      new ( (* (m_tree->ak8jet_p4_))[nfatjets] ) TLorentzVector(j.px(), j.py(), j.pz(), j.energy());
      
      m_tree->ak8jet_tau1_.push_back(j.userFloat("NjettinessAK8:tau1"));
      m_tree->ak8jet_tau2_.push_back(j.userFloat("NjettinessAK8:tau2"));
      m_tree->ak8jet_tau3_.push_back(j.userFloat("NjettinessAK8:tau3"));
      m_tree->ak8jet_softdropmass_.push_back(j.userFloat("ak8PFJetsCHSSoftDropMass"));
      m_tree->ak8jet_trimmedmass_.push_back(j.userFloat("ak8PFJetsCHSTrimmedMass"));
      m_tree->ak8jet_prunedmass_.push_back(j.userFloat("ak8PFJetsCHSPrunedMass"));
      m_tree->ak8jet_filteredmass_.push_back(j.userFloat("ak8PFJetsCHSFilteredMass"));
      m_tree->ak8jet_hasSubjet_.push_back(j.hasSubjets("SoftDrop"));
      
      int nsubjet = 0;
      auto Subjets = j.subjets("SoftDrop");
      for ( auto const & i : Subjets ) {
          new ( (*(m_tree->ak8_subjet_))[nsubjet] ) TLorentzVector(i->px(), i->py(), i->pz(), i->energy());
          m_tree->ak8subjet_btag_.push_back(i->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
          if(verbose_) printf(" %i  w subjet with pt %5.1f (raw pt %5.1f), eta %+4.2f, mass %5.1f ungroomed, btagging %5.1f\n",
                              nsubjet, i->pt(), i->pt()*i->jecFactor("Uncorrected"), i->eta(), i->mass(), 
                              i->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") );
          nsubjet++;
      }
  }

  //Photons
  edm::Handle<pat::PhotonCollection> photons;
  iEvent.getByToken(photonToken_, photons);
  if (m_tree->photon_->GetLast()!=-1) m_tree->photon_->Clear();
  int nphotons = 0;
  for (const pat::Photon &pho : *photons) {
      if (pho.pt() < 20 or pho.chargedHadronIso()/pho.pt() > 0.3) continue;
      if (verbose_)  printf("phot with pt %4.1f, supercluster eta %+5.3f, sigmaIetaIeta %.3f (%.3f with full5x5 shower shapes)\n",
                            pho.pt(), pho.superCluster()->eta(), pho.sigmaIetaIeta(), pho.full5x5_sigmaIetaIeta());
      new ( (*(m_tree->photon_))[nphotons] ) TLorentzVector(pho.px(), pho.py(), pho.pz(), pho.energy());
      nphotons++;
  }
  
  //Taus
  edm::Handle<pat::TauCollection> taus;
  iEvent.getByToken(tauToken_, taus);
  if (m_tree->tau_->GetLast()!=-1) m_tree->tau_->Clear();
  int ntaus = 0;
  for (const pat::Tau &tau : *taus) {
      if (tau.pt() < 20) continue;
      if (verbose_) printf("tau  with pt %4.1f, dxy signif %.1f, ID(byMediumCombinedIsolationDeltaBetaCorr3Hits) %.1f, lead candidate pt %.1f, pdgId %d \n",tau.pt(), tau.dxy_Sig(), tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits"), tau.leadCand()->pt(), tau.leadCand()->pdgId());
      new ( (*(m_tree->tau_))[ntaus] ) TLorentzVector(tau.px(), tau.py(), tau.pz(), tau.energy());
      ntaus++;
  }  
  
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
