#ifndef test_ExtractHitsTracks_Ntuple
#define test_ExtractHitsTracks_Ntuple

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "test/ExtractHitsTracks/interface/Data.h"
#include <vector>

class TTree;

class Ntuple {
  
 public:
  
  static constexpr size_t ARRAY_SIZE_MAX = 100000;
  static constexpr int NEG_INT = -10;
  static constexpr float NEG_FLOAT = -10.;
  static constexpr float NEG_FLOATSQ = -1.*NEG_FLOAT*NEG_FLOAT;
  
  Ntuple() {}
  
  void reset() {
    Ntuple dummy;  // create a new object 
    *this = dummy; // use assignment to reset
  }
  
  void link_tree( TTree* tree );

  void fill_evt( const edm::EventID& id );
  void fill_data( const std::vector<ntuple::Data>& data );
  
 public:

  // Event scalars
  Int_t run_ = 0;
  Int_t lumi_ = 0;
  Int_t evt_ = 0;

  // RecHit
  Int_t nhit_ = 0;
  Int_t hit_n_ = 0;
  Int_t hit_id_[ARRAY_SIZE_MAX] = {};
  Float_t x_[ARRAY_SIZE_MAX] = {};
  Float_t y_[ARRAY_SIZE_MAX] = {};
  Float_t z_[ARRAY_SIZE_MAX] = {};
  
  // GEN
  Int_t particle_id_[ARRAY_SIZE_MAX] = {};
  Int_t pdg_id_[ARRAY_SIZE_MAX] = {};
  Float_t px_[ARRAY_SIZE_MAX] = {};
  Float_t py_[ARRAY_SIZE_MAX] = {};
  
  // SimHit
  Int_t sim_id_[ARRAY_SIZE_MAX] = {};
  Float_t sim_dxy_sig_[ARRAY_SIZE_MAX] = {};
  
  // Geometry
  Int_t volume_id_[ARRAY_SIZE_MAX] = {};
  Int_t layer_id_[ARRAY_SIZE_MAX] = {};
  Int_t module_id_[ARRAY_SIZE_MAX] = {};

}; 

#endif // test_ExtractHitsTracks_Ntuple
