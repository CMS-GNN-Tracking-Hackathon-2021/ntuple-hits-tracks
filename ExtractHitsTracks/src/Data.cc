#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/TrackerRecHit2D/interface/OmniClusterRef.h"
//#include "DataFormats/TrackerRecHit2D/interface/Phase2TrackerRecHit1D.h"
//#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/TrackerSingleRecHit.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "test/ExtractHitsTracks/interface/Data.h"
#include <sstream>
#include <iomanip>

namespace ntuple {
    
  ////////////////////////////////////////////////////////////////////////////////
  //
  void x_y_z( float x, float y, float z, std::stringstream& ss, bool global ) {
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
    ss << " " << std::setw(6) << "HITid"
       << " " << std::setw(6) << "R_2d"
       << " " << std::setw(5) << "phi"
       << " " << std::setw(7) << "Z"
       << " " << std::setw(6) << "R_3d"
       << " " << std::setw(5) << "SIMid"
      //<< " " << std::setw(5) << "PDGid"
       << " " << std::setw(7) << "dxy*1E4"
       << " " << std::setw(7) << "err*1E4"
       << " " << std::setw(7) << "dxy_sig"
       << " " << std::setw(3) << "DET"
       << " " << std::setw(3) << "SID"
       << " " << std::setw(3) << "OT?"
       << " " << std::setw(3) << "SUB"
       << " " << std::setw(3) << "VOL"
       << " " << std::setw(3) << "LAY"
       << " " << std::setw(3) << "MOD"
       << " " << std::setw(9) << "DetId"
       << " " << std::setw(5) << "GENid"
       << " " << std::setw(5) << "PDGid"
       << " " << std::setw(6) << "GEN_pt"
       << " " << std::setw(5) << "eta"
       << " " << std::setw(5) << "phi"
       << " " << std::setw(6) << "TRK_pt"
       << " " << std::setw(5) << "eta"
       << " " << std::setw(5) << "phi";
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
       << std::setprecision(0) << " " << std::setw(6) << data.hitId_
       << std::setprecision(2) << " " << std::setw(6) << float(data.recHit_ == nullptr ? -1. : data.recHit_->globalPosition().perp())
       << std::setprecision(2) << " " << std::setw(5) << float(data.recHit_ == nullptr ? 0. : float(data.recHit_->globalPosition().phi()))
       << std::setprecision(2) << " " << std::setw(7) << float(data.recHit_ == nullptr ? 0. : data.recHit_->globalPosition().z())
       << std::setprecision(2) << " " << std::setw(6) << R
       << std::setprecision(0) << " " << std::setw(5) << int(data.simHit_.isNull() ? -1 : data.simHit_.key())
      //<< std::setprecision(0) << " " << std::setw(5) << int(data.simHit_.isNull() ? -1 : data.simHit_->particleType())
       << std::setprecision(2) << " " << std::setw(7) << float(data.dxy_val_*1.e4 < 10000. ? data.dxy_val_*1.e4 : 9999.99)
       << std::setprecision(4) << " " << std::setw(7) << float(data.dxy_err_*1.e4 < 100. ? data.dxy_err_*1.e4 : 99.9999)
       << std::setprecision(1) << " " << std::setw(7) << data.dxy_sig_
       << std::setprecision(0) << " " << std::setw(3) << int(data.topo_.det_)
       << std::setprecision(0) << " " << std::setw(3) << int(data.topo_.side_)
       << std::setprecision(0) << " " << std::setw(3) << int(data.det_)
       << std::setprecision(0) << " " << std::setw(3) << int(data.topo_.sub_)
       << std::setprecision(0) << " " << std::setw(3) << int(data.topo_.volume_)
       << std::setprecision(0) << " " << std::setw(3) << int(data.topo_.layer_)
       << std::setprecision(0) << " " << std::setw(3) << int(data.topo_.module_)
       << std::setprecision(0) << " " << std::setw(9) << int(data.topo_.raw_)
       << std::setprecision(0) << " " << std::setw(5) << int(data.gen_.isNull() ? -1 : data.gen_.key())
       << std::setprecision(0) << " " << std::setw(5) << int(data.gen_.isNull() ? -1 : data.gen_->pdgId())
       << std::setprecision(2) << " " << std::setw(6) << float(data.gen_.isNull() ? -1. : data.gen_->pt())
       << std::setprecision(2) << " " << std::setw(5) << float(data.gen_.isNull() ? 0. : data.gen_->eta())
       << std::setprecision(2) << " " << std::setw(5) << float(data.gen_.isNull() ? 0. : data.gen_->phi())
       << std::setprecision(2) << " " << std::setw(6) << float(data.track_.isNull() ? -1. : data.track_->pt())
       << std::setprecision(2) << " " << std::setw(5) << float(data.track_.isNull() ? 0. : data.track_->eta())
       << std::setprecision(2) << " " << std::setw(5) << float(data.track_.isNull() ? 0. : data.track_->phi());
    out << ss.str();
    return out;
  };
  
  ////////////////////////////////////////////////////////////////////////////////
  //
  std::ostream& print( std::ostream& out, const std::vector<Data>& vdata, int n ) {
    std::stringstream ss;
    ss << "Print summary of hits in std::vector<Data> (size: "
       << vdata.size()
       << ") from first " 
       << n
       << " particles: "
       << std::endl;
    unsigned int idx = std::numeric_limits<unsigned int>::max();
    int cntr = 0;
    for ( auto data : vdata ) { 
      if ( n > 0 && cntr > n ) { continue; } 
      if ( data.gen_.key() != idx ) {
  	header(ss);
  	ss << std::endl;
  	idx = data.gen_.key();
  	++cntr;
      } 
      terse(ss,data);
      ss << std::endl;
    }
    out << ss.str();
    return out;
  };

};
