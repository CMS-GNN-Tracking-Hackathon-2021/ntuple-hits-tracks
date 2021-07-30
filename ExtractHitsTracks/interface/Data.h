#ifndef test_ExtractHitsTracks_Data
#define test_ExtractHitsTracks_Data

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/TrackerRecHit2D/interface/OmniClusterRef.h"
#include "DataFormats/TrackerRecHit2D/interface/TrackerSingleRecHit.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "test/ExtractHitsTracks/interface/Topo.h"

class SiPixelRecHit;
class Phase2TrackerRecHit1D;

namespace ntuple {
  
  ////////////////////////////////////////////////////////////////////////////////
  //
  enum class Type { Signal = 0, ITPileup = 1, OOTPileup = 2, Noise = 3, Unknown = 99 };
  enum class Det { Inner = 0, Outer = 1, Unknown = 99 };
  typedef edm::Ref<std::vector<PSimHit> > SimHitRef;
  
  ////////////////////////////////////////////////////////////////////////////////
  //
  class Data {

  public:

    Data() {;}
    Data(size_t hit_id) { hitId_ = hit_id; }
    ~Data() {;}

    bool operator< (const Data& data) const {
      float data_R = data.recHit_ == nullptr ? -1 : 
	sqrt( data.recHit_->globalPosition().x()*data.recHit_->globalPosition().x() + 
	      data.recHit_->globalPosition().y()*data.recHit_->globalPosition().y() + 
	      data.recHit_->globalPosition().z()*data.recHit_->globalPosition().z() );
      float this_R = this->recHit_ == nullptr ? -1 : 
	sqrt( this->recHit_->globalPosition().x()*this->recHit_->globalPosition().x() + 
	      this->recHit_->globalPosition().y()*this->recHit_->globalPosition().y() + 
	      this->recHit_->globalPosition().z()*this->recHit_->globalPosition().z() );
      return ( (int(this->gen_.key()) < int(data.gen_.key())) // Order by GEN idx
	       || ( (int(this->gen_.key()) == int(data.gen_.key())) && (this_R < data_R) ) ); // Then R (3D)
      if (0) { // TURNED OFF - Order by GEN idx then DetId info 
	return ( (int(this->gen_.key()) < int(data.gen_.key()))
		 || ( (int(this->gen_.key()) == int(data.gen_.key())) 
		      && (int(this->det_) < int(data.det_)) )
		 || ( (int(this->gen_.key()) == int(data.gen_.key())) 
		      && (int(this->det_) == int(data.det_))
		      && (int(this->topo_.sub_) < int(data.topo_.sub_)) )
		 || ( (int(this->gen_.key()) == int(data.gen_.key())) 
		      && (int(this->det_) == int(data.det_))
		      && (int(this->topo_.sub_) == int(data.topo_.sub_))
		      && (int(this->topo_.layer_) < int(data.topo_.layer_)) )
		 || ( (int(this->gen_.key()) == int(data.gen_.key())) 
		      && (int(this->det_) == int(data.det_))
		      && (int(this->topo_.sub_) == int(data.topo_.sub_))
		      && (int(this->topo_.layer_) == int(data.topo_.layer_))
		      && (int(this->topo_.module_) < int(data.topo_.module_)) )
		 );
      }
    }

  public:

    // RecHit and cluster info
    size_t hitId_ = 0;
    const TrackerSingleRecHit* recHit_ = nullptr;
    const OmniClusterRef* clu_ = nullptr;

    // DET id
    Det det_ = Det::Unknown;
    Topo topo_;

    // GEN (for target) and TP info
    uint32_t idx_ = -1;
    reco::GenParticleRef gen_;
    TrackingParticleRef tp_;

    //  SimHit info
    SimHitRef simHit_;
    Type type_ = Type::Unknown;
    float dxy_sig_ = 0.;
    float dxy_val_ = 0.;
    float dxy_err_ = 0.;
    const SimTrack* simTrack_ = nullptr;

    // RECO track
    edm::RefToBase<reco::Track> track_;

  };
  
  ////////////////////////////////////////////////////////////////////////////////
  //
  template <typename T> Det innerOrOuter( const T* hit ) { return Det::Unknown; }
  inline Det innerOrOuter( const SiPixelRecHit* hit ) { return Det::Inner; }
  inline Det innerOrOuter( const Phase2TrackerRecHit1D* hit ) { return Det::Outer; }
  
  ////////////////////////////////////////////////////////////////////////////////
  //
  void x_y_z( float x, float y, float z, std::stringstream& ss, bool global = true );
  void x_y_z( const GlobalPoint& p, std::stringstream& ss );
  void x_y_z( const LocalPoint& p, std::stringstream& ss );
  void r_phi_z( float r, float phi, float z, std::stringstream& ss );
  void r_phi_z( GlobalPoint p, std::stringstream& ss );
  void pt_eta_phi( double pt, double eta, double phi, std::stringstream& ss );
  template <typename T> void pt_eta_phi( const T* ptr, std::stringstream& ss );
  std::ostream& operator<< ( std::ostream& out, const Data& data );
  std::ostream& header( std::ostream& out );
  std::ostream& terse( std::ostream& out, const Data& data );
  std::ostream& print( std::ostream& out, const std::vector<Data>& vdata, int n = -1 );

}; 

#endif // test_ExtractHitsTracks_Data
