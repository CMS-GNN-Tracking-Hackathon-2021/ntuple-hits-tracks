#include "test/ExtractHitsTracks/interface/Ntuple.h"
#include "TTree.h"

////////////////////////////////////////////////////////////////////////////////
//
void Ntuple::link_tree( TTree *tree ) {

  // Event scalars
  tree->Branch("run",  &run_ , "run/I");
  tree->Branch("lumi", &lumi_, "lumi/I");
  tree->Branch("evt",  &evt_ , "evt/I");

  // RecHit 
  tree->Branch("nhit", &nhit_, "nhit/I"); // Number of RecHits in an event
  tree->Branch("hit_n", &hit_n_, "hit_n/I"); // Number of RecHits *considered* (up to ARRAY_SIZE_MAX)
  tree->Branch("hit_id", &hit_id_, "hit_id[hit_n]/I");
  tree->Branch("x", &x_, "x[hit_n]/F");
  tree->Branch("y", &y_, "y[hit_n]/F");
  tree->Branch("z", &z_, "z[hit_n]/F");
  
  // GEN
  tree->Branch("particle_id", &particle_id_, "particle_id[hit_n]/I");
  tree->Branch("pdg_id", &pdg_id_, "pdg_id[hit_n]/I");
  tree->Branch("px", &px_, "px[hit_n]/F");
  tree->Branch("py", &py_, "py[hit_n]/F");

  // SimHit 
  tree->Branch("sim_id", &sim_id_, "sim_id[hit_n]/I");
  tree->Branch("sim_dxy_sig", &sim_dxy_sig_, "sim_dxy_sig[hit_n]/F");
  
  // Geometry
  tree->Branch("volume_id", &volume_id_, "volume_id[hit_n]/I");
  tree->Branch("layer_id", &layer_id_, "layer_id[hit_n]/I");
  tree->Branch("module_id", &module_id_, "module_id[hit_n]/I");
  
}

////////////////////////////////////////////////////////////////////////////////
//
void Ntuple::fill_evt( const edm::EventID& id ) {
  run_  = id.run();
  lumi_ = id.luminosityBlock();
  evt_  = id.event();
}

////////////////////////////////////////////////////////////////////////////////
//
void Ntuple::fill_data( const std::vector<ntuple::Data>& vdata ) {

  //std::cout << "vdata.size(): " << vdata.size() << std::endl;

  unsigned int idx = 0;
  for ( auto const& data : vdata ) {
    if ( idx >= ARRAY_SIZE_MAX ) { break; }

    // RecHit 
    hit_id_[idx] = idx;
    x_[idx] = (data.recHit_ == nullptr ? 0. : data.recHit_->globalPosition().x());
    y_[idx] = (data.recHit_ == nullptr ? 0. : data.recHit_->globalPosition().y());
    z_[idx] = (data.recHit_ == nullptr ? 0. : data.recHit_->globalPosition().z());

    // GEN
    particle_id_[idx] = (data.gen_.isNull() ? -1 : data.gen_.key());
    pdg_id_[idx] = (data.gen_.isNull() ? 0 : data.gen_->pdgId());
    px_[idx] = (data.gen_.isNull() ? 0 : data.gen_->px());
    py_[idx] = (data.gen_.isNull() ? 0 : data.gen_->py());

    // SimHit 
    sim_id_[idx] = (data.simHit_.isNull() ? -1 : data.simHit_.key());
    sim_dxy_sig_[idx] = data.dxy_sig_;

    // Geometry
    volume_id_[idx] = int(data.det_);
    layer_id_[idx] = data.topo_.layer_;
    module_id_[idx] = data.topo_.module_;
    
    // Counter
    idx++;
  }
  hit_n_ = idx;
  nhit_ = vdata.size();

}

//////////////////////////////////////////////////////////////////////////////////
////
//void Ntuple::fill_clu( edmNew::DetSetVector<SiPixelCluster> const* siPixelClusters,
//		       edmNew::DetSetVector<Phase2TrackerCluster1D> const* siPhase2Clusters ) {
//
//  // siPixelClusters
//  //std::cout << "siPixelClusters: ";
//  unsigned int pixel_clu_n = 0;
//  for ( auto const& detset : *siPixelClusters ) {
//    //std::cout << " detset.detId(): " << detset.detId();
//    for ( auto const& clu : detset ) {
//      if ( pixel_clu_n >= ARRAY_SIZE_MAX ) { break; }
//      //std::cout << " clu.size(): " << clu.size();
//      pixel_clu_size_[pixel_clu_n] = clu.size();
//      pixel_clu_n++;
//    }
//    //std::cout << std::endl;
//  }
//  pixel_clu_n_ = pixel_clu_n;
//
//  // siPhase2Clusters
//  //std::cout << "siPhase2Clusters: ";
//  unsigned int strip_clu_n = 0;
//  for ( auto const& detset : *siPhase2Clusters ) {
//    //std::cout << " detset.detId(): " << detset.detId();
//    for ( auto const& clu : detset ) {
//      if ( strip_clu_n >= ARRAY_SIZE_MAX ) { break; }
//      //std::cout << " clu.size(): " << clu.size();
//      strip_clu_size_[strip_clu_n] = clu.size();
//      strip_clu_n++;
//    }
//    //std::cout << std::endl;
//  }
//  strip_clu_n_ = strip_clu_n;
//
//}
