#ifndef test_ExtractHitsTracks_Ntuple
#define test_ExtractHitsTracks_Ntuple

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "test/ExtractHitsTracks/interface/Data.h"
#include <vector>

class TTree;

class Ntuple {
  
 public:
  
  static constexpr size_t ARRAY_SIZE_MAX = 500000;
  
  Ntuple() {}
  void reset();
  void link_tree(TTree* tree);

  void fill_evt(const edm::EventID& id);
  void fill_data(const std::vector<ntuple::Data>& data, size_t& index);
  
 public:

  // Event-level scalars

  Int_t run_ = 0;
  Int_t lumi_ = 0;
  Int_t evt_ = 0;

  Int_t nhit_ = 0;
  Int_t hit_n_ = 0;

  // Event-level arrays

  // RecHit
  std::vector<int> hit_id_;

  std::vector<float> x_;
  std::vector<float> y_;
  std::vector<float> z_;

  // GEN
  std::vector<int> particle_id_;
  std::vector<int> pdg_id_;
  std::vector<float> gen_pt_;
  std::vector<float> gen_eta_;
  std::vector<float> gen_phi_;

  // SimHit
  std::vector<int> sim_type_;
  std::vector<int> sim_id_;
  std::vector<float> sim_dxy_sig_;
  std::vector<float> sim_pt_;
  std::vector<float> sim_eta_;
  std::vector<float> sim_phi_;

  // Geometry
  std::vector<int> volume_id_;
  std::vector<int> layer_id_;
  std::vector<int> module_id_;

}; 

#endif // test_ExtractHitsTracks_Ntuple
