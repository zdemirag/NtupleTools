#ifndef MonoJetFlatTree_H
#define MonoJetFlatTree_H

#include <set>
#include <vector>
#include <string>
#include <utility>

#include "TFile.h"
#include "TTree.h"
#include "TError.h"
#include "Math/VectorUtil.h"
#include "Math/LorentzVector.h"

//
// Mono-jet Dark Matter analysis
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/Monojet
//

class MonoJetFlatTree {
 public:
  /// float doesn't have dictionary before ROOT6, so use double
  typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

  /// bit map
  /// DON'T CHANGE ORDER
  enum Selection {
    BaseLine           = 1UL<<0,  // synchronization selection (all needed requirements are applied)
    Trigger            = 1UL<<30  // passed a set of triggers
  };

  /// DON'T CHANGE ORDER
  enum DataType {
    data,
    ttbar,
    tw,
    dyee,
    dymm,
    dytt,
    dynn,
    wjets,
    wz,
    zz,
    wgamma,
    qcd,
    other
  };

  /// Core variables
  unsigned int   run_;        // Run number
  unsigned int   lumi_;       // Lumi number
  unsigned int   event_;      // Event number
  DataType       dstype_;     // Sample type
  unsigned int   cuts_;       // Precomputed selection
  unsigned int   nEvents_;    // Number of generated events processed 
  float          scale1fb_;   // Event weight to scale event yields to 1/fb of data

  /// variables
  unsigned int   nvtx_;       // Number of primary vertices
  float          rho_;        // rho of the event
  unsigned int   nMuons_;     // Number of muons passing the "loose" selection
  unsigned int   nElectrons_; // Number of electrons passing the "loose" selection
  unsigned int   nTaus_;      // Number of taus passing the "loose" selection
  unsigned int   nPhotons_;   // Number of photons passing the "loose" selection
  LorentzVector  jet1_;       // Leading jet p4
  float          jet1RawPt_;  // Raw pt of the leading jet 
  float          jet1CHFrac_; // Charged hadron fraction of the leading jet
  float          jet1NHFrac_; // Neutral hadron fraction of the leading jet
  float          jet1EmFrac_; // Neutral electromagnetic fraction of the leading jet
  LorentzVector  jet2_;       // 2nd jet p4
  float          jet2RawPt_;  // Raw pt of the 2nd jet 
  float          jet2CHFrac_; // Charged hadron fraction of the 2nd jet
  float          jet2NHFrac_; // Neutral hadron fraction of the 2nd jet
  float          jet2EmFrac_; // Neutral electromagnetic fraction of the 2nd jet
  LorentzVector  jet3_;       // 3d jet p4
  float          jet3RawPt_;  // Raw pt of the 3d jet 
  float          jet3CHFrac_; // Charged hadron fraction of the 3d jet
  float          jet3NHFrac_; // Neutral hadron fraction of the 3d jet
  float          jet3EmFrac_; // Neutral electromagnetic fraction of the 3d jet
  float          met_;        // MET in the event
  float          metPhi_;     // MET in the event
  TNamed         info_;       // Information about ntuple

  private:
  // need a pointer for each LorentzVector
  LorentzVector* jet1_ptr_;
  LorentzVector* jet2_ptr_;
  LorentzVector* jet3_ptr_;
  
 public:
  /// this is the main element
  TTree *tree_;
  TFile *f_;
  
  /// default constructor  
  MonoJetFlatTree():info_("info","MonoJet Dark Matter analysis"),
    jet1_ptr_(&jet1_), jet2_ptr_(&jet2_), jet3_ptr_(&jet3_){}

  /// default destructor
  ~MonoJetFlatTree(){ 
    if (f_) f_->Close();  
  };

  /// initialize varibles and fill list of available variables
  void InitVariables();

  /// load a MonoJetFlatTree
  void LoadTree(const char* file, int type = -1){
    f_ = TFile::Open(file);
    assert(f_);
    tree_ = dynamic_cast<TTree*>(f_->Get("tree"));
    assert(tree_);
  }

  /// create a MonoJetFlatTree
  void CreateTree(){
    tree_ = new TTree("tree","MonoJet Dark Matter analysis");
    f_ = 0;
    InitVariables();
    //book the branches
    tree_->Branch("run"          , &run_          ,   "run/i");
    tree_->Branch("lumi"         , &lumi_         ,   "lumi/i");
    tree_->Branch("event"        , &event_        ,   "event/i");
    tree_->Branch("dstype"       , &dstype_       ,   "dstype/I");
    tree_->Branch("cuts"         , &cuts_         ,   "cuts/i");
    tree_->Branch("nEvents"      , &nEvents_      ,   "dstype/I");
    tree_->Branch("scale1fb"     , &scale1fb_     ,   "scale1fb/F");

    tree_->Branch("nvtx"         , &nvtx_         ,   "nvtx/i");
    tree_->Branch("rho"          , &rho_          ,   "rho/F");
    tree_->Branch("nMuons"       , &nMuons_       ,   "nMuons/I");
    tree_->Branch("nElectrons"   , &nElectrons_   ,   "nElectrons/I");
    tree_->Branch("nTaus"        , &nTaus_        ,   "nTaus/I");
    tree_->Branch("nPhotons"     , &nPhotons_     ,   "nPhotons/I");
    tree_->Branch("jet1"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jet1_ptr_);
    tree_->Branch("jet1RawPt"    , &jet1RawPt_    ,   "jet1RawPt/F");
    tree_->Branch("jet1CHFrac"   , &jet1CHFrac_   ,   "jet1CHFrac/F");
    tree_->Branch("jet1NHFrac"   , &jet1NHFrac_   ,   "jet1NHFrac/F");
    tree_->Branch("jet1EmFrac"   , &jet1EmFrac_   ,   "jet1EmFrac/F");
    tree_->Branch("jet2"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jet2_ptr_);
    tree_->Branch("jet2RawPt"    , &jet2RawPt_    ,   "jet2RawPt/F");
    tree_->Branch("jet2CHFrac"   , &jet2CHFrac_   ,   "jet2CHFrac/F");
    tree_->Branch("jet2NHFrac"   , &jet2NHFrac_   ,   "jet2NHFrac/F");
    tree_->Branch("jet2EmFrac"   , &jet2EmFrac_   ,   "jet2EmFrac/F");
    tree_->Branch("jet3"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jet3_ptr_);
    tree_->Branch("jet3RawPt"    , &jet3RawPt_    ,   "jet3RawPt/F");
    tree_->Branch("jet3CHFrac"   , &jet3CHFrac_   ,   "jet3CHFrac/F");
    tree_->Branch("jet3NHFrac"   , &jet3NHFrac_   ,   "jet3NHFrac/F");
    tree_->Branch("jet3EmFrac"   , &jet3EmFrac_   ,   "jet3EmFrac/F");
    tree_->Branch("met"          , &met_          ,   "met/F");
    tree_->Branch("metPhi"       , &metPhi_       ,   "metPhi/F");
  }

  // initialze a MonoJetFlatTree
  void InitTree(){
    assert(tree_);
    // don't forget to set pointers to zero before you set address
    // or you will fully appreciate that "ROOT sucks" :)
    InitVariables();
    //Set branch address
    Int_t currentState = gErrorIgnoreLevel;
    // gErrorIgnoreLevel = kError;
    gErrorIgnoreLevel = kBreak;
    tree_->SetBranchAddress("run",           &run_);
    tree_->SetBranchAddress("lumi",          &lumi_);
    tree_->SetBranchAddress("event",         &event_);
    tree_->SetBranchAddress("dstype",        &dstype_);
    tree_->SetBranchAddress("cuts",          &cuts_);
    tree_->SetBranchAddress("nEvents",       &nEvents_);
    tree_->SetBranchAddress("scale1fb",      &scale1fb_);

    tree_->SetBranchAddress("nvtx",          &nvtx_);
    tree_->SetBranchAddress("rho"          , &rho_);
    tree_->SetBranchAddress("nMuons"       , &nMuons_);
    tree_->SetBranchAddress("nElectrons"   , &nElectrons_);
    tree_->SetBranchAddress("nTaus"        , &nTaus_);
    tree_->SetBranchAddress("nPhotons"     , &nPhotons_);
    tree_->SetBranchAddress("jet1"         , &jet1_ptr_);
    tree_->SetBranchAddress("jet1RawPt"    , &jet1RawPt_);
    tree_->SetBranchAddress("jet1CHFrac"   , &jet1CHFrac_);
    tree_->SetBranchAddress("jet1NHFrac"   , &jet1NHFrac_);
    tree_->SetBranchAddress("jet1EmFrac"   , &jet1EmFrac_);
    tree_->SetBranchAddress("jet2"         , &jet2_ptr_);
    tree_->SetBranchAddress("jet2RawPt"    , &jet2RawPt_);
    tree_->SetBranchAddress("jet2CHFrac"   , &jet2CHFrac_);
    tree_->SetBranchAddress("jet2NHFrac"   , &jet2NHFrac_);
    tree_->SetBranchAddress("jet2EmFrac"   , &jet2EmFrac_);
    tree_->SetBranchAddress("jet3"         , &jet3_ptr_);
    tree_->SetBranchAddress("jet3RawPt"    , &jet3RawPt_);
    tree_->SetBranchAddress("jet3CHFrac"   , &jet3CHFrac_);
    tree_->SetBranchAddress("jet3NHFrac"   , &jet3NHFrac_);
    tree_->SetBranchAddress("jet3EmFrac"   , &jet3EmFrac_);
    tree_->SetBranchAddress("met"          , &met_); 
    tree_->SetBranchAddress("metPhi"       , &metPhi_);
    gErrorIgnoreLevel = currentState;
  }

  /// transform DateType to string
  static std::string name(DataType type){
    switch (type){
    case data:   return "data";
    case ttbar:  return "ttbar";
    case tw:     return "tw";
    case dyee:   return "dyee";
    case dymm:   return "dymm";
    case dytt:   return "dytt";
    case dynn:   return "dynn";
    case wjets:  return "wjets";
    case wz:     return "wz";
    case zz:     return "zz";
    case wgamma: return "wgamma";
    case qcd:    return "qcd";
    case other:  return "other";
    default:     return "unknown";
    }
  };

}; 

inline void 
MonoJetFlatTree::InitVariables(){
  // inizialize variables
  run_           = 0;
  lumi_          = 0;
  event_         = 0;
  dstype_        = data;
  cuts_          = 0;
  nEvents_       = 0;
  scale1fb_      = 0;

  nvtx_          = 0;
  rho_           = 0;
  nMuons_        = 0;
  nElectrons_    = 0;
  nTaus_         = 0;
  nPhotons_      = 0;
  jet1RawPt_     = -999.;
  jet1CHFrac_    = -999.;
  jet1NHFrac_    = -999.;
  jet1EmFrac_    = -999.;
  jet2RawPt_     = -999.;
  jet2CHFrac_    = -999.;
  jet2NHFrac_    = -999.;
  jet2EmFrac_    = -999.;
  jet3RawPt_     = -999.;
  jet3CHFrac_    = -999.;
  jet3NHFrac_    = -999.;
  jet3EmFrac_    = -999.;
  met_           = -999.;
  metPhi_        = -999.;
}

#endif
