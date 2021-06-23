#ifndef test_ExtractHitsTracks_Ntuple
#define test_ExtractHitsTrakcs_Ntuple

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/ParticleFlowReco/interface/PreIdFwd.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include <vector>

class TTree;

namespace reco { typedef edm::Ptr<Track> TrackPtr; }

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
  void fill_trk( const reco::TrackPtr& trk,
		 const reco::BeamSpot& spot );
  
  unsigned int run_ = 0;
  unsigned int lumi_ = 0;
  unsigned long long evt_ = 0;

}; 

#endif // test_ExtractHitsTracks_Ntuple
