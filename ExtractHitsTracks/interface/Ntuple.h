#ifndef test_ExtractHitsTracks_Ntuple
#define test_ExtractHitsTrakcs_Ntuple

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/ParticleFlowReco/interface/PreIdFwd.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Phase2TrackerCluster/interface/Phase2TrackerCluster1D.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include <vector>

class TTree;

constexpr size_t ARRAY_SIZE_MAX = 10000;

class Ntuple {
  
 public:
  
  static constexpr size_t NHITS_MAX = 30;
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
  void fill_clu( edmNew::DetSetVector<SiPixelCluster> const* siPixelClusters,
		 edmNew::DetSetVector<Phase2TrackerCluster1D> const* siPhase2Clusters );
  
 public:
  
  unsigned int run_ = 0;
  unsigned int lumi_ = 0;
  unsigned long long evt_ = 0;
  
  int pixel_clu_n_ = 0;
  int pixel_clu_size_[ARRAY_SIZE_MAX] = {};
  int strip_clu_n_ = 0;
  int strip_clu_size_[ARRAY_SIZE_MAX] = {};

}; 

#endif // test_ExtractHitsTracks_Ntuple
