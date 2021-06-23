#include "test/ExtractHitsTracks/interface/Ntuple.h"
#include "TTree.h"

////////////////////////////////////////////////////////////////////////////////
//
void Ntuple::link_tree( TTree *tree ) {
  tree->Branch("run",  &run_ , "run/i");
  tree->Branch("lumi", &lumi_, "lumi/i");
  tree->Branch("evt",  &evt_ , "evt/i");
}

////////////////////////////////////////////////////////////////////////////////
//
void Ntuple::fill_evt( const edm::EventID& id ) {
  run_  = id.run();
  lumi_ = id.luminosityBlock();
  evt_  = id.event();
}
