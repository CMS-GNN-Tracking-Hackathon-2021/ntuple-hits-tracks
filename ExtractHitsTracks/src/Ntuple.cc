#include "test/ExtractHitsTracks/interface/Ntuple.h"
#include "TTree.h"

////////////////////////////////////////////////////////////////////////////////
//
void Ntuple::link_tree( TTree *tree ) {

  // Event scalars
  tree->Branch("run",  &run_ , "run/i");
  tree->Branch("lumi", &lumi_, "lumi/i");
  tree->Branch("evt",  &evt_ , "evt/i");

  // Cluster features
  tree->Branch("pixel_clu_n", &pixel_clu_n_, "pixel_clu_n/I");
  tree->Branch("pixel_clu_size", &pixel_clu_size_, "pixel_clu_size[pixel_clu_n]/I");
  tree->Branch("strip_clu_n", &strip_clu_n_, "strip_clu_n/I");
  tree->Branch("strip_clu_size", &strip_clu_size_, "strip_clu_size[strip_clu_n]/I");
  
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
void Ntuple::fill_clu( edmNew::DetSetVector<SiPixelCluster> const* siPixelClusters,
		       edmNew::DetSetVector<Phase2TrackerCluster1D> const* siPhase2Clusters ) {

  // siPixelClusters
  //std::cout << "siPixelClusters: ";
  unsigned int pixel_clu_n = 0;
  for ( auto const& detset : *siPixelClusters ) {
    //std::cout << " detset.detId(): " << detset.detId();
    for ( auto const& clu : detset ) {
      if ( pixel_clu_n >= ARRAY_SIZE_MAX ) { break; }
      //std::cout << " clu.size(): " << clu.size();
      pixel_clu_size_[pixel_clu_n] = clu.size();
      pixel_clu_n++;
    }
    //std::cout << std::endl;
  }
  pixel_clu_n_ = pixel_clu_n;

  // siPhase2Clusters
  //std::cout << "siPhase2Clusters: ";
  unsigned int strip_clu_n = 0;
  for ( auto const& detset : *siPhase2Clusters ) {
    //std::cout << " detset.detId(): " << detset.detId();
    for ( auto const& clu : detset ) {
      if ( strip_clu_n >= ARRAY_SIZE_MAX ) { break; }
      //std::cout << " clu.size(): " << clu.size();
      strip_clu_size_[strip_clu_n] = clu.size();
      strip_clu_n++;
    }
    //std::cout << std::endl;
  }
  strip_clu_n_ = strip_clu_n;

}
