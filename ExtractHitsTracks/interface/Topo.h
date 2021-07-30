#ifndef test_ExtractHitsTracks_Topo
#define test_ExtractHitsTracks_Topo

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include <boost/core/demangle.hpp>
#include<iostream>

namespace ntuple {
  
  ////////////////////////////////////////////////////////////////////////////////
  //
  class Topo {

  public:

    Topo() {;}
    Topo(const TrackerTopology* topo, const DetId& id) { init(topo,id); }

    void init(const TrackerTopology* topo, const DetId& id);
    static uint32_t volume(const TrackerTopology* topo, const DetId& id);

  public:

    // https://github.com/cms-sw/cmssw/blob/master/DataFormats/DetId/interface/DetId.h
    // Tracker = 1, etc
    int det_ = 0; 

    // https://github.com/cms-sw/cmssw/blob/master/DataFormats/SiPixelDetId/interface/PixelSubdetector.h
    // https://github.com/cms-sw/cmssw/blob/master/DataFormats/SiStripDetId/interface/StripSubdetector.h
    // https://github.com/cms-sw/cmssw/tree/master/Geometry/TrackerGeometryBuilder
    // SiStrip interpretation: Skype chat with Marco Munisch, TRK POG for GNN, 29 Jul 2021 
    // (DetId::subDetId()==1) or (DetId::subdetId()==PixelSubdetector::PixelBarrel) means Phase2 Inner Tracker Barrel
    // (DetId::subDetId()==2) or (DetId::subdetId()==PixelSubdetector::PixelBarrel) means Phase2 Inner Tracker Endcap
    // (DetId::subDetId()==5) or (DetId::subdetId()==StripSubdetector::TOB) means Phase2 Outer Tracker Barrel
    // (DetId::subDetId()==4) or (DetId::subdetId()==StripSubdetector::TID) means Phase2 Outer Tracker Endcap
    int sub_ = 0; 

    int layer_ = 0;

    int module_ = 0;
    
    // https://github.com/cms-sw/cmssw/blob/master/DataFormats/TrackerCommon/interface/TrackerDetSide.h
    // Barrel = 0, NegEndcap = 1, PosEndcap = 2
    int side_ = 0; 
    
    // See Topo::volume for definition
    uint32_t volume_ = 0; 
    
    // Raw DetId
    uint32_t raw_ = 0;

  };

  ////////////////////////////////////////////////////////////////////////////////
  // Print detector information
  std::ostream& operator<< (std::ostream& out, const Topo& topo);  

}; 

#endif // test_ExtractHitsTracks_Topo
