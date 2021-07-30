#include "test/ExtractHitsTracks/interface/Topo.h"
#include <sstream>
#include <iomanip>

namespace ntuple {
    
  void Topo::init(const TrackerTopology* topo, const DetId& id) {
    det_ = id.det();
    sub_ = id.subdetId();
    layer_ = topo->layer(id);
    module_ = topo->module(id);
    side_ = topo->side(id);
    raw_ = id.rawId();
    volume_ = Topo::volume(topo,id);
  };
  
  uint32_t Topo::volume(const TrackerTopology* topo, const DetId& id) {
    if      ( topo == nullptr )                       return 0;
    else if ( id.subdetId()==2 && topo->side(id)==1 ) return 1; // IT endcap -ve
    else if ( id.subdetId()==1 && topo->side(id)==0 ) return 2; // IT barrel
    else if ( id.subdetId()==2 && topo->side(id)==2 ) return 3; // IT endcap +ve
    else if ( id.subdetId()==4 && topo->side(id)==1 ) return 4; // OT endcap -ve
    else if ( id.subdetId()==5 && topo->side(id)==0 ) return 5; // OT barrel
    else if ( id.subdetId()==4 && topo->side(id)==2 ) return 6; // OT endcap +ve
    else                                              return 0;
  };
  
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

};
