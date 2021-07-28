#ifndef test_ExtractHitsTracks_Data
#define test_ExtractHitsTracks_Data

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/TrackerRecHit2D/interface/OmniClusterRef.h"
#include "DataFormats/TrackerRecHit2D/interface/TrackerSingleRecHit.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include <boost/core/demangle.hpp>

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
  class Topo {
  public:
    Topo() {;}
    Topo(const TrackerTopology* topo, const DetId& id) { init(topo,id); }
    void init(const TrackerTopology* topo, const DetId& id) {
      det_ = id.det();
      sub_ = id.subdetId();
      layer_ = topo->layer(id);
      module_ = topo->module(id);
      side_ = topo->side(id);
      raw_ = id.rawId();
    };
  public:
    int det_ = 0; // Tracker = 1, https://github.com/cms-sw/cmssw/blob/master/DataFormats/DetId/interface/DetId.h
    int sub_ = 0; // PixelBarrel = 1, PixelEndcap = 2, https://github.com/cms-sw/cmssw/blob/master/DataFormats/SiPixelDetId/interface/PixelSubdetector.h
    int layer_ = 0;
    int module_ = 0;
    int side_ = 0; // Barrel = 0, NegEndcap = 1, PosEndcap = 2, https://github.com/cms-sw/cmssw/blob/master/DataFormats/TrackerCommon/interface/TrackerDetSide.h
    uint32_t raw_ = 0;
  };

  ////////////////////////////////////////////////////////////////////////////////
  //
  std::ostream& operator<< (std::ostream& out, const Topo& topo) {
    std::stringstream ss;
    ss << std::setprecision(0) 
       << " DET: " << std::setw(2) << topo.det_
       << " SUB: " << std::setw(2) << topo.sub_
       << " LAYER: " << std::setw(2) << topo.layer_
       << " MODULE: " << std::setw(2) << topo.module_
       << " SIDE: " << std::setw(2) << topo.side_
       << " RAW: " << std::setw(-1) << topo.raw_;
    out << ss.str();
    return out;
  };
  
  ////////////////////////////////////////////////////////////////////////////////
  //
  class Data {
  public:
    Data() {;}
    ~Data() { init(); }
    void init() {;}
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
    uint32_t idx_ = -1;
    Det det_ = Det::Unknown;
    Topo topo_;
    TrackingParticleRef tp_;
    reco::GenParticleRef gen_;
    Type type_ = Type::Unknown;
    const SimTrack* simTrack_ = nullptr;
    edm::RefToBase<reco::Track> track_;
    const TrackerSingleRecHit* recHit_ = nullptr;
    const OmniClusterRef* clu_ = nullptr;
    SimHitRef simHit_;
    float dxy_sig_ = 0.;
    float dxy_val_ = 0.;
    float dxy_err_ = 0.;
  };
  
  ////////////////////////////////////////////////////////////////////////////////
  //
  template <typename T> Det innerOrOuter( const T* hit ) { return Det::Unknown; }
  Det innerOrOuter( const SiPixelRecHit* hit ) { return Det::Inner; }
  Det innerOrOuter( const Phase2TrackerRecHit1D* hit ) { return Det::Outer; }
  
  ////////////////////////////////////////////////////////////////////////////////
  //
  void x_y_z( float x, float y, float z, std::stringstream& ss, bool global = true ) {
    if (global)
      ss << " x,y,z: " 
	 << std::fixed
	 << std::setprecision(1) 
	 << std::setw(6) << x << ", " 
	 << std::setw(6) << y << ", " 
	 << std::setw(6) << z;
    else 
      ss << " x,y,z: " 
	 << std::fixed
	 << std::setprecision(2) 
	 << std::setw(6) << x << ", " 
	 << std::setw(6) << y << ", " 
	 << std::setw(6) << z;
      
  };  

  ////////////////////////////////////////////////////////////////////////////////
  //
  void x_y_z( const GlobalPoint& p, std::stringstream& ss ) {
    x_y_z( p.x(), p.y(), p.z(), ss, true );
  };  

  ////////////////////////////////////////////////////////////////////////////////
  //
  void x_y_z( const LocalPoint& p, std::stringstream& ss ) {
    x_y_z( p.x(), p.y(), p.z(), ss, false );
  };  

  ////////////////////////////////////////////////////////////////////////////////
  //
  void r_phi_z( float r, float phi, float z, std::stringstream& ss ) {
    ss << " r,phi,z: " 
       << std::fixed
       << std::setprecision(2) 
       << std::setw(5) << r << ", " 
       << std::setprecision(1) 
       << std::setw(5) << phi << ", " 
       << std::setprecision(2) 
       << std::setw(6) << z;
  };  

  ////////////////////////////////////////////////////////////////////////////////
  //
  void r_phi_z( GlobalPoint p, std::stringstream& ss ) {
    r_phi_z( p.perp(), p.phi(), p.z(), ss );
  };  

  ////////////////////////////////////////////////////////////////////////////////
  //
  void pt_eta_phi( double pt, double eta, double phi, std::stringstream& ss ) {
    ss << " pt,eta,phi: " 
       << std::fixed
       << std::setprecision(2) 
       << std::setw(5) << pt << ", " 
       << std::setw(4) << eta << ", " 
       << std::setw(4) << phi;
  };  

  ////////////////////////////////////////////////////////////////////////////////
  //
  template <typename T> 
    void pt_eta_phi( const T* ptr, std::stringstream& ss ) {
    if ( ptr == nullptr ) { ss << " Null pointer..."; }
    else { pt_eta_phi( ptr->pt(), ptr->eta(), ptr->phi(), ss ); }
  };  
  
  ////////////////////////////////////////////////////////////////////////////////
  //
  std::ostream& operator<< ( std::ostream& out, const Data& data ) {
    std::stringstream ss;
    ss << "Data";
    ss << "\n Trajectories:";
    ss << "\n  GEN:  "; pt_eta_phi(data.gen_.get(),ss); 
    ss << " idx: " << data.gen_.key() << " pdgId: " << data.gen_->pdgId() << " ID: " << data.gen_.id();
    ss << "\n  TP:   "; pt_eta_phi(data.tp_.get(),ss);
    ss << "\n  SIM:  "; pt_eta_phi(data.simTrack_->momentum().pt(),
				   data.simTrack_->momentum().eta(),
				   data.simTrack_->momentum().phi(),
				   ss);
    ss << "\n  RECO: "; pt_eta_phi(data.track_.get(),ss);
    ss << "\n DetId: "
       << " IT/OT: " << int(data.det_)
       << " " << data.topo_;
    ss << "\n Hit:";
    ss << "\n  REC:  ";
    r_phi_z(data.recHit_->globalPosition(), ss);
    //x_y_z(data.recHit_->globalPosition(), ss);
    ss << " (global)";
    ss << "\n  REC:  ";
    x_y_z(data.recHit_->localPosition(), ss);
    ss << " (local)";
    //ss << "\n  CLU:  ";
    //x_y_z(data.clu_->localPosition(), ss);
    //ss << " (local)";
    ss << "\n  SIM:  ";
    x_y_z(data.simHit_->localPosition(), ss);
    ss << " (local)";
    out << ss.str();
    return out;
  };
  
  ////////////////////////////////////////////////////////////////////////////////
  //
  std::ostream& header( std::ostream& out ) {
    std::stringstream ss;
    ss << " " << std::setw(5) << "IDX"
       << " " << std::setw(5) << "PDGid"
       << " " << std::setw(6) << "GEN_pt"
       << " " << std::setw(5) << "eta"
       << " " << std::setw(5) << "phi"
       << " " << std::setw(6) << "REC_pt"
       << " " << std::setw(5) << "eta"
       << " " << std::setw(5) << "phi"
       << " " << std::setw(3) << "OT?"
       << " " << std::setw(3) << "DET"
       << " " << std::setw(3) << "SUB"
       << " " << std::setw(3) << "LAY"
       << " " << std::setw(3) << "MOD"
       << " " << std::setw(3) << "SID"
       << " " << std::setw(9) << "DetId"
       << " " << std::setw(6) << "R_2d"
       << " " << std::setw(5) << "phi"
       << " " << std::setw(7) << "Z"
       << " " << std::setw(6) << "R_3d"
       << " " << std::setw(5) << "SIM"
      //<< " " << std::setw(5) << "PDGid"
       << " " << std::setw(7) << "dxy*1E4"
       << " " << std::setw(7) << "err*1E4"
       << " " << std::setw(7) << "sig";
    out << ss.str();
    return out;
  }
  
  ////////////////////////////////////////////////////////////////////////////////
  //
  std::ostream& terse( std::ostream& out, const Data& data ) {
    std::stringstream ss;
    float R = data.recHit_ == nullptr ? -1 : sqrt( data.recHit_->globalPosition().x()*data.recHit_->globalPosition().x() + 
						  data.recHit_->globalPosition().y()*data.recHit_->globalPosition().y() + 
						  data.recHit_->globalPosition().z()*data.recHit_->globalPosition().z() );
    ss << std::fixed
       << std::setprecision(0) << " " << std::setw(5) << (data.gen_.isNull() ? -1 : data.gen_.key())
       << std::setprecision(0) << " " << std::setw(5) << (data.gen_.isNull() ? -1 : data.gen_->pdgId())
       << std::setprecision(2) << " " << std::setw(6) << (data.gen_.isNull() ? -1. : data.gen_->pt())
       << std::setprecision(2) << " " << std::setw(5) << (data.gen_.isNull() ? 0. : data.gen_->eta())
       << std::setprecision(2) << " " << std::setw(5) << (data.gen_.isNull() ? 0. : data.gen_->phi())
       << std::setprecision(2) << " " << std::setw(6) << (data.track_.isNull() ? -1. : data.track_->pt())
       << std::setprecision(2) << " " << std::setw(5) << (data.track_.isNull() ? 0. : data.track_->eta())
       << std::setprecision(2) << " " << std::setw(5) << (data.track_.isNull() ? 0. : data.track_->phi())
       << std::setprecision(2) << " " << std::setw(3) << int(data.det_)
       << std::setprecision(0) << " " << std::setw(3) << data.topo_.det_
       << std::setprecision(0) << " " << std::setw(3) << data.topo_.sub_
       << std::setprecision(0) << " " << std::setw(3) << data.topo_.layer_
       << std::setprecision(0) << " " << std::setw(3) << data.topo_.module_
       << std::setprecision(0) << " " << std::setw(3) << data.topo_.side_
       << std::setprecision(0) << " " << std::setw(9) << data.topo_.raw_
       << std::setprecision(2) << " " << std::setw(6) << (data.recHit_ == nullptr ? -1. : data.recHit_->globalPosition().perp())
       << std::setprecision(2) << " " << std::setw(5) << (data.recHit_ == nullptr ? 0. : float(data.recHit_->globalPosition().phi()))
       << std::setprecision(2) << " " << std::setw(7) << (data.recHit_ == nullptr ? 0. : data.recHit_->globalPosition().z())
       << std::setprecision(2) << " " << std::setw(6) << R
       << std::setprecision(0) << " " << std::setw(5) << (data.simHit_.isNull() ? -1 : data.simHit_.key())
      //<< std::setprecision(0) << " " << std::setw(5) << (data.simHit_.isNull() ? -1 : data.simHit_->particleType())
       << std::setprecision(2) << " " << std::setw(7) << (data.dxy_val_*1.e4 < 10000. ? data.dxy_val_*1.e4 : 9999.99)
       << std::setprecision(4) << " " << std::setw(7) << (data.dxy_err_*1.e4 < 100. ? data.dxy_err_*1.e4 : 99.9999)
       << std::setprecision(1) << " " << std::setw(7) << data.dxy_sig_;
    out << ss.str();
    return out;
  };
  
  ////////////////////////////////////////////////////////////////////////////////
  //
  std::ostream& print( std::ostream& out, const std::vector<Data>& vdata ) {
    std::stringstream ss;
    ss << "Print summary of std::vector<Data> content, size: "
       << vdata.size()
       << std::endl;
    unsigned int idx = std::numeric_limits<unsigned int>::max();
    for ( auto data : vdata ) { 
      if ( data.gen_.key() != idx ) {
	ntuple::header(ss);
	ss << std::endl;
	idx = data.gen_.key();
      } 
      ntuple::terse(ss,data);
      ss << std::endl;
    }
    out << ss.str();
    return out;
  };
  
}; 

#endif // test_ExtractHitsTracks_Data
