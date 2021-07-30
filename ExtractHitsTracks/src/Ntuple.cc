#include "test/ExtractHitsTracks/interface/Ntuple.h"
#include "TTree.h"

////////////////////////////////////////////////////////////////////////////////
//
void Ntuple::reset() {
  // create a new object 
  Ntuple dummy;
  // delete this object
  // clear vectors???
  // use assignment to reset
  *this = dummy; 
}

////////////////////////////////////////////////////////////////////////////////
//
void Ntuple::link_tree(TTree *tree) {

  // Event-level scalars

  tree->Branch("run", &run_ );
  tree->Branch("lumi", &lumi_);
  tree->Branch("evt", &evt_);

  tree->Branch("nhit", &nhit_); // Number of RecHits in an event
  tree->Branch("hit_n", &hit_n_); // Number of RecHits *considered* (up to ARRAY_SIZE_MAX)

  // Event-level arrays

  // RecHit 
  tree->Branch("hit_id", &hit_id_);
  tree->Branch("x", &x_);
  tree->Branch("y", &y_);
  tree->Branch("z", &z_);
  
  // GEN
  tree->Branch("particle_id", &particle_id_);
  tree->Branch("pdg_id", &pdg_id_);
//  tree->Branch("px", &px_);
//  tree->Branch("py", &py_);
  tree->Branch("gen_pt", &gen_pt_);
  tree->Branch("gen_eta", &gen_eta_);
  tree->Branch("gen_phi", &gen_phi_);

  // SimHit 
  tree->Branch("sim_type", &sim_type_);
  tree->Branch("sim_id", &sim_id_);
  tree->Branch("sim_dxy_sig", &sim_dxy_sig_);
  tree->Branch("sim_pt", &sim_pt_);
  tree->Branch("sim_eta", &sim_eta_);
  tree->Branch("sim_phi", &sim_phi_);
  
  // Geometry
  tree->Branch("volume_id", &volume_id_);
  tree->Branch("layer_id", &layer_id_);
  tree->Branch("module_id", &module_id_);
  
}

////////////////////////////////////////////////////////////////////////////////
//
void Ntuple::fill_evt(const edm::EventID& id) {
  run_  = id.run();
  lumi_ = id.luminosityBlock();
  evt_  = id.event();
}

////////////////////////////////////////////////////////////////////////////////
//
void Ntuple::fill_data(const std::vector<ntuple::Data>& vdata, size_t& index ) {

  for ( auto const& data : vdata ) {
    if ( index >= ARRAY_SIZE_MAX ) { continue; }

    // RecHit 
    hit_id_.push_back(int(data.hitId_));
    
    x_.push_back(float(data.recHit_ == nullptr ? 0. : data.recHit_->globalPosition().x()));
    y_.push_back(float(data.recHit_ == nullptr ? 0. : data.recHit_->globalPosition().y()));
    z_.push_back(float(data.recHit_ == nullptr ? 0. : data.recHit_->globalPosition().z()));

    // GEN
    particle_id_.push_back(int(data.gen_.isNull() ? -1 : data.gen_.key()));
    pdg_id_.push_back(int(data.gen_.isNull() ? 0 : data.gen_->pdgId()));
//    px_.push_back(float(data.gen_.isNull() ? 0 : data.gen_->px()));
//    py_.push_back(float(data.gen_.isNull() ? 0 : data.gen_->py()));
    gen_pt_.push_back(float(data.gen_.isNull() ? 0 : data.gen_->pt()));
    gen_eta_.push_back(float(data.gen_.isNull() ? 0 : data.gen_->eta()));
    gen_phi_.push_back(float(data.gen_.isNull() ? 0 : data.gen_->phi()));

    // SimHit 
    sim_type_.push_back(int(data.type_));
    sim_id_.push_back(int(data.simHit_.isNull() ? -1 : data.simHit_.key()));
    sim_dxy_sig_.push_back(float(data.dxy_sig_));
    sim_pt_.push_back(float(data.simTrack_ == nullptr ? 0 : data.simTrack_->momentum().pt()));
    sim_eta_.push_back(float(data.simTrack_ == nullptr ? 0 : data.simTrack_->momentum().eta()));
    sim_phi_.push_back(float(data.simTrack_ == nullptr ? 0 : data.simTrack_->momentum().phi()));
    
    // Geometry
    volume_id_.push_back(int(data.topo_.volume_));
    layer_id_.push_back(int(data.topo_.layer_));
    module_id_.push_back(int(data.topo_.module_));
    
    // Counter
    index++;
  }

  hit_n_ = index; // ... index reflects signal and bkgd, so final call will set correctly
  nhit_ += vdata.size(); // Increment each time this method is called (signal+bkgd)

}
