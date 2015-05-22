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
#include "TClonesArray.h"
#include "TLorentzVector.h"

#include "TMath.h"
#include "TH1F.h"

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
  TNamed         info_;       // Information about ntuple
 

  //AK4 Jets
  TClonesArray       *ak4jet_;
  std::vector<float> ak4jet_btag_;
  std::vector<int>   ak4jet_puid_;

  //AK8 Jets
  TClonesArray *ak8jet_p4_, *ak8_subjet_;
  std::vector<float> ak8jet_tau1_, ak8jet_tau2_, ak8jet_tau3_;
  std::vector<float> ak8jet_softdropmass_, ak8jet_trimmedmass_;
  std::vector<float> ak8jet_prunedmass_, ak8jet_filteredmass_;
  std::vector<int>   ak8jet_hasSubjet_, ak8subjet_btag_;

  //MET
  TLorentzVector pfmet_, genmet_;
  //Photon
  TClonesArray *photon_;
  //Tau
  TClonesArray *tau_;
  
  private:
  // need a pointer for each LorentzVector
  
 public:
  /// this is the main element
  TTree *tree_;
  TFile *f_;
  
  /// default constructor  
 MonoJetFlatTree():info_("info","MonoJet Dark Matter analysis"){}

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
    //tree_->Branch("dstype"       , &dstype_       ,   "dstype/I");
    //tree_->Branch("cuts"         , &cuts_         ,   "cuts/i");
    //tree_->Branch("nEvents"      , &nEvents_      ,   "dstype/I");
    //tree_->Branch("scale1fb"     , &scale1fb_     ,   "scale1fb/F");
    //tree_->Branch("nvtx"         , &nvtx_         ,   "nvtx/I");
    //tree_->Branch("rho"          , &rho_          ,   "rho/F");

    //MET
    tree_->Branch("pfmet","TLorentzVector",&pfmet_);
    tree_->Branch("genmet","TLorentzVector",&genmet_);

    // AK4 Jets
    ak4jet_ = new TClonesArray("TLorentzVector", 20);
    tree_->Branch("ak4jet"       , "TClonesArray" ,   &ak4jet_, 128000, 0);
    tree_->Branch("ak4jet_btag"  , &ak4jet_btag_);
    tree_->Branch("ak4jet_puid"  , &ak4jet_puid_);
    
    // Fat Jets
    ak8jet_p4_  = new TClonesArray("TLorentzVector", 20);
    ak8_subjet_ = new TClonesArray("TLorentzVector", 20);
    tree_->Branch("ak8jet"        ,"TClonesArray", &ak8jet_p4_, 128000, 0);
    tree_->Branch("ak8_subjet"    ,"TClonesArray", &ak8_subjet_, 128000, 0);
    tree_->Branch("ak8jet_tau1"        , &ak8jet_tau1_);
    tree_->Branch("ak8jet_tau2"        , &ak8jet_tau2_);
    tree_->Branch("ak8jet_tau3"        , &ak8jet_tau3_);
    tree_->Branch("ak8jet_softdropmass", &ak8jet_softdropmass_);
    tree_->Branch("ak8jet_trimmedmass" , &ak8jet_trimmedmass_);
    tree_->Branch("ak8jet_prunedmass"  , &ak8jet_prunedmass_);
    tree_->Branch("ak8jet_filteredmass", &ak8jet_filteredmass_);
    tree_->Branch("ak8subjet_btag"     , &ak8subjet_btag_);

    photon_ = new TClonesArray("TLorentzVector", 20);
    tree_->Branch("photon"        ,"TClonesArray", &photon_, 128000, 0);

    tau_ = new TClonesArray("TLorentzVector", 20);
    tree_->Branch("tau"        ,"TClonesArray", &tau_, 128000, 0);

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

    tree_->SetBranchAddress("ak4jet_",       &ak4jet_);
    tree_->SetBranchAddress("ak4jet_btag_",  &ak4jet_btag_);
    tree_->SetBranchAddress("ak4jet_puid_",  &ak4jet_puid_);

    tree_->SetBranchAddress("ak8jet_p4_"          , &ak8jet_p4_);
    tree_->SetBranchAddress("ak8_subjet_"         , &ak8_subjet_);
    tree_->SetBranchAddress("ak8jet_tau1_"        , &ak8jet_tau1_);
    tree_->SetBranchAddress("ak8jet_tau2_"        , &ak8jet_tau2_);
    tree_->SetBranchAddress("ak8jet_tau3_"        , &ak8jet_tau3_);
    tree_->SetBranchAddress("ak8jet_softdropmass_", &ak8jet_softdropmass_);
    tree_->SetBranchAddress("ak8jet_trimmedmass_" , &ak8jet_trimmedmass_);
    tree_->SetBranchAddress("ak8jet_prunedmass_"  , &ak8jet_prunedmass_);
    tree_->SetBranchAddress("ak8jet_filteredmass_", &ak8jet_filteredmass_);
    tree_->SetBranchAddress("ak8subjet_btag_"     , &ak8subjet_btag_);

    tree_->SetBranchAddress("pfmet_",&pfmet_);
    tree_->SetBranchAddress("genmet_",&genmet_);

    tree_->SetBranchAddress("photon_",&photon_);
    tree_->SetBranchAddress("tau_",&tau_);

 
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

}

#endif
