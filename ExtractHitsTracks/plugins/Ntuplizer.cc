#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Phase2TrackerCluster/interface/Phase2TrackerCluster1D.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/Phase2TrackerRecHit1D.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimGeneral/TrackingAnalysis/interface/SimHitTPAssociationProducer.h"
#include "SimTracker/TrackerHitAssociation/interface/ClusterTPAssociation.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "test/ExtractHitsTracks/interface/Ntuple.h"
#include "TTree.h"
#include<iostream> 
#include<string> 

////////////////////////////////////////////////////////////////////////////////
//
class Ntuplizer : public edm::EDFilter {
  
public:
  
  explicit Ntuplizer(const edm::ParameterSet&);
  ~Ntuplizer();
  
  virtual bool filter(edm::Event&, const edm::EventSetup&) override;

  template <typename T> void match(const edmNew::DetSetVector<T>* hits, size_t& hit_id);
  void match2(std::vector<ntuple::Data>&);
  void readCollections(edm::Event&, const edm::EventSetup&);
  void debug();

  // LIKELY TO BE DEPRECATED, USED BY match2()
  using value_type = ClusterTPAssociation::value_type;
  static bool compare(const value_type& i, const value_type& j) { 
    return i.second.key() > j.second.key(); 
  }
  static bool compareSort(const value_type& i, const value_type& j) {
    return compare(i,j) || (i.second.key() == j.second.key() && i.first.rawIndex() > j.first.rawIndex());
  }
  
private: 
  
  edm::Service<TFileService> fs_;
  TTree* tree_;	
  Ntuple ntuple_;
  int verbose_;
  std::vector<ntuple::Data> signal_;
  std::vector<ntuple::Data> bkgd_;

  std::vector<int> activeTrackingRegions_;

  edm::EDGetTokenT<SiPixelRecHitCollection> pixelRecHitsToken_;
  SiPixelRecHitCollection const* pixelRecHits_;

  edm::EDGetTokenT<Phase2TrackerRecHit1DCollectionNew> trackerRecHitsToken_;
  Phase2TrackerRecHit1DCollectionNew const* trackerRecHits_;

  edm::EDGetTokenT<ClusterTPAssociation> clustersToTPToken_;
  ClusterTPAssociation const* clustersToTP_;

  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
  edm::Handle<reco::GenParticleCollection> genParticles_;
  //reco::GenParticleCollection const* genParticles_;

  bool usePrunedGenParticles_;
  edm::EDGetTokenT<edm::Association<reco::GenParticleCollection> > prunedGenParticlesToken_;
  //edm::Handle<edm::Association<reco::GenParticleCollection> > prunedGenParticles_;
  edm::Association<reco::GenParticleCollection> const* prunedGenParticles_;

  edm::EDGetTokenT<SimHitTPAssociationProducer::SimHitTPAssociationList> tpToSimHitsMapToken_;
  SimHitTPAssociationProducer::SimHitTPAssociationList const* tpToSimHitsMap_;

  const edm::EDGetTokenT<edm::View<reco::Track> > ctfTracksToken_;
  edm::Handle<edm::View<reco::Track> > ctfTracks_;

  const edm::EDGetTokenT<TrackingParticleCollection> trackingParticlesToken_;
  edm::Handle<TrackingParticleCollection> trackingParticles_;

  //const edm::EDGetTokenT<TrackingParticleCollection> prunedTrackingParticlesToken_;
  //TrackingParticleCollection const* prunedTrackingParticles_;

  edm::EDGetTokenT<reco::TrackToTrackingParticleAssociator> tpToTracksToken_;
  const reco::TrackToTrackingParticleAssociator* tpToTracks_;

  const edm::ESGetToken<TrackerTopology,TrackerTopologyRcd> topoToken_;
  const TrackerTopology* topo_;


}; 

////////////////////////////////////////////////////////////////////////////////
//
Ntuplizer::~Ntuplizer(){}

////////////////////////////////////////////////////////////////////////////////
//
Ntuplizer::Ntuplizer( const edm::ParameterSet& cfg ) : 
  tree_(nullptr),
  ntuple_(),
  verbose_(cfg.getParameter<int>("verbose")),
  signal_(),
  bkgd_(),
  activeTrackingRegions_(cfg.getParameter<std::vector<int> >("activeTrackingRegions")),
  pixelRecHitsToken_(consumes<SiPixelRecHitCollection>(cfg.getParameter<edm::InputTag>("pixelRecHits"))),
  pixelRecHits_(),
  trackerRecHitsToken_(consumes<Phase2TrackerRecHit1DCollectionNew>(cfg.getParameter<edm::InputTag>("trackerRecHits"))),
  trackerRecHits_(),
  clustersToTPToken_(consumes<ClusterTPAssociation>(cfg.getParameter<edm::InputTag>("clustersToTP"))),
  clustersToTP_(),
  genParticlesToken_(consumes<reco::GenParticleCollection>(cfg.getParameter<edm::InputTag>("genParticles"))),
  genParticles_(),
  usePrunedGenParticles_(cfg.getParameter<bool>("usePrunedGenParticles")),
  prunedGenParticlesToken_(consumes<edm::Association<reco::GenParticleCollection> >(cfg.getParameter<edm::InputTag>("prunedGenParticles"))),
  prunedGenParticles_(),
  tpToSimHitsMapToken_(consumes<SimHitTPAssociationProducer::SimHitTPAssociationList>(cfg.getParameter<edm::InputTag>("tpToSimHits"))),
  tpToSimHitsMap_(),
  ctfTracksToken_(consumes<edm::View<reco::Track> >(cfg.getParameter<edm::InputTag>("ctfTracks"))),
  ctfTracks_(),
  trackingParticlesToken_(consumes<TrackingParticleCollection>(cfg.getParameter<edm::InputTag>("trackingParticles"))),
  trackingParticles_(),
//prunedTrackingParticlesToken_(consumes<TrackingParticleCollection>(cfg.getParameter<edm::InputTag>("prunedTrackingParticles"))),
//prunedTrackingParticles_(),
  tpToTracksToken_(consumes<reco::TrackToTrackingParticleAssociator>(cfg.getParameter<std::string>("tpToTracks"))),
  tpToTracks_(),
  topoToken_(esConsumes())
{
  tree_ = fs_->make<TTree>("tree","tree");
  ntuple_.link_tree(tree_);
  if (verbose_>2) std::cout << "[Ntuplizer::Ntuplizer]" << std::endl;
  std::cout << "[Ntuplizer::Ntuplizer]" << std::endl
	    << " Verbosity level: "<< verbose_ << std::endl
	    << " Number of activeTrackingRegions: "<< activeTrackingRegions_.size() << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
//
bool Ntuplizer::filter(edm::Event& event, const edm::EventSetup& setup ) {
  if (verbose_>2) std::cout << "[Ntuplizer::filter]" << std::endl;
  
  // Init
  ntuple_.reset();
  readCollections(event,setup);
  signal_.clear();
  bkgd_.clear();

  // Populate
  size_t hit_id = 0;
  match<SiPixelRecHit>(pixelRecHits_,hit_id);
  match<Phase2TrackerRecHit1D>(trackerRecHits_,hit_id);

  std::sort(signal_.begin(), signal_.end());
  if (verbose_>1 && !signal_.empty()) {
    std::cout << "SIGNAL HITS" << std::endl;
    std::stringstream ss; 
    ntuple::print(ss,signal_,10); // Print hits for first 10 particles only
    std::cout << ss.str() << std::endl;
  }

  if (verbose_>2 && !bkgd_.empty()) {
    int n = 20;
    std::cout << "BKGD HITS" << std::endl;
    std::stringstream ss;
    ss << "Print summary of first "
       << n
       << " hits in std::vector<Data> (size: "
       << bkgd_.size()
       << ")" 
       << std::endl;
    int cntr = 0;
    ntuple::header(ss);
    ss << std::endl;
    for ( auto data : bkgd_ ) { 
      if ( n > 0 && cntr > n ) { continue; } 
      ntuple::terse(ss,data);
      ss << std::endl;
      ++cntr;
    }
    std::cout << ss.str() << std::endl;
  }

  //////////
  // LIKELY TO BE DEPRECATED
  //  std::vector<ntuple::Data> vdata;
  //  match2(vdata);
  //  std::sort(vdata.begin(), vdata.end());
  //  if (verbose_>1) {
  //    std::stringstream ss; 
  //    ntuple::print(ss,vdata);
  //    std::cout << ss.str() << std::endl;
  //  }
  //////////

  // Fill Ntuple
  ntuple_.fill_evt(event.id());
  size_t index = 0;
  ntuple_.fill_data(signal_,index);
  ntuple_.fill_data(bkgd_,index);
  tree_->Fill(); 
  return true; 

} 

////////////////////////////////////////////////////////////////////////////////
//
void Ntuplizer::readCollections( edm::Event& event, const edm::EventSetup& setup ) { 
  if (verbose_>2) std::cout << "[Ntuplizer::readCollections]" << std::endl;
  pixelRecHits_ = &event.get(pixelRecHitsToken_);
  trackerRecHits_ = &event.get(trackerRecHitsToken_);
  clustersToTP_ = &event.get(clustersToTPToken_);
  //genParticles_ = &event.get(genParticlesToken_);
  event.getByToken(genParticlesToken_,genParticles_);
  prunedGenParticles_ = &event.get(prunedGenParticlesToken_);
  tpToSimHitsMap_ = &event.get(tpToSimHitsMapToken_);
  event.getByToken(trackingParticlesToken_,trackingParticles_);
  //prunedTrackingParticles_ = &event.get(prunedTrackingParticlesToken_);
  event.getByToken(ctfTracksToken_,ctfTracks_);
  tpToTracks_ = &event.get(tpToTracksToken_);
  topo_ = &setup.getData(topoToken_);
  debug();
}

////////////////////////////////////////////////////////////////////////////////
//
void Ntuplizer::debug() {
  if (verbose_>2) std::cout << "[Ntuplizer::debug]" << std::endl;
  if (verbose_<1) { return; }

  int pixelRecHits = 0;
  for ( auto const& detset : *pixelRecHits_ ) { pixelRecHits += detset.size(); }
  int trackerRecHits = 0;
  for ( auto const& detset : *trackerRecHits_ ) { trackerRecHits += detset.size(); }
  int tpToSimHitsMap = 0; // see .dev
  int prunedGenParticles = 0; 
  for ( unsigned int idx = 0; idx < genParticles_->size(); ++idx ) {
    if ( prunedGenParticles_->get(idx).isNonnull() ) { prunedGenParticles++; }
  }
  //int tpToTracks = 0; // see .dev
  std::cout << "[Ntuplizer::readCollections]" << std::endl
	    << " pixelRecHits_->size(): " << pixelRecHits
	    << " (dsv: " << pixelRecHits_->size() << ")" << std::endl
	    << " trackerRecHits_->size(): " << trackerRecHits
	    << " (dsv: " << trackerRecHits_->size() << ")" << std::endl
	    << " clustersToTP_->size(): " << clustersToTP_->size() << std::endl
	    << " genParticles_->size(): " << genParticles_->size() << std::endl
	    << " prunedGenParticles: " << prunedGenParticles
	    << " (all: " << prunedGenParticles_->size() << ")" << std::endl
	    << " tpToSimHitsMap_->size(): " << tpToSimHitsMap
	    << " (tp: " << tpToSimHitsMap_->size() << ")" << std::endl
	    << " trackingParticles_->size(): " << trackingParticles_->size() << std::endl
    //<< " prunedTrackingParticles_->size(): " << prunedTrackingParticles_->size() << std::endl
	    << " ctfTracks_->size(): " << ctfTracks_->size() << std::endl;

//  for ( unsigned int idx = 0; idx < genParticles_->size(); ++idx ) {
//    reco::GenParticleRef gen(genParticles_,idx);
//    if ( (*prunedGenParticles_)[gen].isNonnull() ) {
//      std::cout << "pruned:"
//		<< " idx: " << gen.key() 
//		<< " pdg: " << gen->pdgId()
//		<< std::endl;
//    }
//  }

}

////////////////////////////////////////////////////////////////////////////////
//
template <typename T> 
void Ntuplizer::match( const edmNew::DetSetVector<T>* hits, size_t& hit_id ) {
  if (verbose_>2) std::cout << "[Ntuplizer::match]" << std::endl;

  reco::SimToRecoCollection associations = tpToTracks_->associateSimToReco(ctfTracks_,trackingParticles_);

  // Iterate through DetSetVector of RecHits
  for ( auto detSet = hits->begin(); detSet != hits->end(); ++detSet ) {
    for ( auto recHit = detSet->begin(); recHit != detSet->end(); ++recHit ) {

      // Initialise container (with unique id)
      ntuple::Data data(hit_id);
      
      // Retrieve cluster
      const OmniClusterRef& clu = recHit->firstClusterRef();

      // Store RecHit, Cluster, and DetId info (via Topo class)
      data.recHit_ = &*recHit;
      data.clu_ = &clu;
      data.topo_ = ntuple::Topo(topo_,detSet->detId());
      data.det_ = ntuple::innerOrOuter(&*recHit);
      
      // Check if RecHit is in one of the active tracking regions (if set)
      auto iter = std::find(activeTrackingRegions_.begin(),
			    activeTrackingRegions_.end(), 
			    int(data.topo_.volume_));
      if (!activeTrackingRegions_.empty() && 
	  iter == activeTrackingRegions_.end()) { continue; }

      // Find associated TrackingParticles (can be more than one per cluster!)
      auto range_tp = clustersToTP_->equal_range(clu);

      // Skip if none
      if ( range_tp.first == range_tp.second ) { continue; }

      // Loop through TPs
      for ( auto tp = range_tp.first; tp != range_tp.second; ++tp ) {

	// Set TP
	data.tp_ = tp->second;

	// If match already found, skip remaining TPs
	//if ( data.gen_.isNonnull() ) { continue; }

	// Match one of GenParticles from TP against the pruned collection
	const reco::GenParticleRefVector& genParticles = tp->second->genParticles();
	for ( unsigned int idx = 0; idx < genParticles.size(); ++idx ) {
	  const reco::GenParticleRef& gen = genParticles[idx];
	  if ( !usePrunedGenParticles_ || (*prunedGenParticles_)[gen].isNonnull() ) {
	    data.gen_ = gen;
	    data.idx_ = gen.key();
	    break; // register first GEN match
	  }
	}
	
	// If no match found, move onto next TP
	//if ( data.gen_.isNull() ) { continue; }
	
	// Record local position of RecHit (to find closest SimHit)
	float x_rh = 0., y_rh = 0., xx_rh = 0., xy_rh = 0., yy_rh = 0.;
	x_rh = data.recHit_->localPosition().x();
	y_rh = data.recHit_->localPosition().y();
	xx_rh = data.recHit_->localPositionError().xx();
	xy_rh = data.recHit_->localPositionError().xy();
	yy_rh = data.recHit_->localPositionError().yy();

	// Store reco::Track
	reco::SimToRecoCollection::const_iterator iter = associations.find(tp->second);
	if (iter != associations.end()) {
	  edm::RefToBase<reco::Track> trk = iter->val.front().first;
	  data.track_ = iter->val.front().first;
	}

	// Is the TrackingParticle from signal, or in-time PU, or OOT PU, or noise (or UNKNOWN by default)
	data.type_ = ntuple::Type::Noise;
	const auto event = data.tp_->eventId().event();
	const auto bx = data.tp_->eventId().bunchCrossing();
	ntuple::Type type = ntuple::Type::OOTPileup;
	if ( bx == 0 ) { type = (event == 0 ? ntuple::Type::Signal : ntuple::Type::ITPileup); }
	data.type_ = static_cast<ntuple::Type>( std::min(static_cast<int>(data.type_), static_cast<int>(type)) );

	// Store SimTrack
	if ( !(data.tp_->g4Tracks().empty()) ) { data.simTrack_ = &(data.tp_->g4Tracks()[0]); }
	
	// Find matched SimHits
	std::pair<TrackingParticleRef, TrackPSimHitRef> dummy(tp->second, TrackPSimHitRef());
	auto range_simhits = std::equal_range(tpToSimHitsMap_->begin(),
					      tpToSimHitsMap_->end(),
					      dummy, // The match is based only on the TP
					      SimHitTPAssociationProducer::simHitTPAssociationListGreater);
	
	// Store SimHits
	float sig_min = -1.;
	float dxy_val = 0.;
	float dxy_err = 0.;
	auto iter_min = range_simhits.second; 
	for ( auto iter = range_simhits.first; iter != range_simhits.second; ++iter ) {
	  ntuple::SimHitRef simHit = iter->second;
	  DetId simHitDetId = DetId(simHit->detUnitId());
	  if (simHitDetId.rawId() == data.topo_.raw_) {
	    // Skip electron SimHits for non-electron TPs (why? doesn't ever seem to occur...)
	    if (std::abs(simHit->particleType()) == 11 && std::abs(data.tp_->pdgId()) != 11) { continue; }
	    // Significance of RecHit-SimHit residual w.r.t. RecHit uncertainties
	    // https://github.com/cms-sw/cmssw/blob/master/Validation/RecoTrack/plugins/TrackingNtuple.cc#L2470
            float dx = simHit->localPosition().x() - x_rh;
            float dy = simHit->localPosition().y() - y_rh;
            float err = xx_rh*yy_rh - xy_rh*xy_rh;
	    float dxy = dx*dx*yy_rh - 2.*dx*dy*xy_rh + dy*dy*xx_rh;
            float sig = ( err > 0. ? dxy / err : 1.e4 );
            if (sig_min < 0. || sig < sig_min) {
              sig_min = sig;
	      dxy_val = dxy;
	      dxy_err = err;
	      iter_min = iter;
            }
	  } // If same DetId
	} // Iterate through SimHits

	// Store SimHit with smallest residual significance w.r.t. RecHit
	if ( iter_min != range_simhits.second ) { 
	  data.simHit_ = iter_min->second;//.get();
	  data.dxy_sig_ = sig_min;
	  data.dxy_val_ = dxy_val;
	  data.dxy_err_ = dxy_err;
	}
	
      } // Loop through TPs

      // Add container to one of the list, depending if GEN match is found
      if ( data.gen_.isNonnull() ) { signal_.push_back(data); }
      else                         { bkgd_.push_back(data); }
      
      // Increment (unique) hit_id counter, important!
      ++hit_id;

    } // Loop through RecHits in DS
  } // Loop through DSs in DSV
  
}

////////////////////////////////////////////////////////////////////////////////
// LIKELY TO BE DEPRECATED, DO NOT USE
void Ntuplizer::match2(std::vector<ntuple::Data>& vdata) {
  if (verbose_>2) std::cout << "[Ntuplizer::match]" << std::endl;

  reco::SimToRecoCollection associations = tpToTracks_->associateSimToReco(ctfTracks_,trackingParticles_);
  
  // Copy entire ClustersToTP map and sort by TP (then clusters)
  ClusterTPAssociation::map_type tpToClusters(clustersToTP_->map());
  std::sort(tpToClusters.begin(), tpToClusters.end(), compareSort);

  for ( unsigned int idx = 0; idx < trackingParticles_->size(); ++idx ) {
    TrackingParticleRef tpref(trackingParticles_,idx);
//    std::cout << " tp idx " << idx 
//	      << std::endl;

    // Filter TP by GenParticles from the pruned collection
    reco::GenParticleRef gen;
    const reco::GenParticleRefVector& genParticles = tpref->genParticles();
    //for ( unsigned int idx = 0; idx < genParticles.size(); ++idx ) {
    for ( auto& tmp : genParticles ) {
      if ( (*prunedGenParticles_)[tmp].isNonnull() ) { 
//	std::cout << " tp idx " << idx
//		  << " out of " << trackingParticles_->size()
//		  << ", gen idx " << tmp.key() 
//		  << ", gen pt " << tmp->pt()
//		  << std::endl;
	gen = tmp;
	break; // register first GEN match
      }
    }
    
    // If no match found, move onto next TP
    if ( gen.isNull() ) { continue; }
    //std::cout << " HERE " << std::endl;
    
    // Find all clusters matching the TP
    auto range_clusters = std::equal_range(tpToClusters.begin(),
					   tpToClusters.end(),
					   value_type(OmniClusterRef(),tpref), // The match is based only on the TP
					   compare);
//    std::cout << " range size " << std::distance(range_clusters.first,range_clusters.second)
//	      << std::endl;
    
    for ( auto iter = range_clusters.first; iter != range_clusters.second; ++iter ) {
//      std::cout << " tp idx " << idx 
//		<< " tp pdg " << tpref->pdgId() 
//		<< " tp key " << iter->second.key()
//		<< " clu key " << iter->first.key()
//		<< std::endl;
      
      // Initialise container
      vdata.emplace_back(ntuple::Data());
      const OmniClusterRef* clu = &(iter->first);
      vdata.back().clu_ = clu;

//      uint32_t detId = clu->isPixel() ? 
//	clu->cluster_pixel().detUnit().geographicalId().det() :
//	clu->isPhase2() ?
//	clu->phase2OTCluster().detUnit().geographicalId().det() :
//	0;
      
      //vdata.back().data.topo_ = ntuple::Topo(topo_,detSet->detId());
      //vdata.back().recHit_ = ???;
      //vdata.back().data.det_ = ntuple::innerOrOuter(&*recHit);
      vdata.back().tp_ = iter->second;//.get();
      vdata.back().gen_ = gen;
      vdata.back().idx_ = gen.key();
    }

//    for ( const auto& tp : tpToClusters ) {
//      if ( tpref.key() == tp.second.key() ) {
//	std::cout << " clu key " << tp.first.key() 
//		  << " tp clu " << tp.second.key() 
//		  << std::endl;
//      }
//    }
    
  }
  
}
      
////////////////////////////////////////////////////////////////////////////////
//
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(Ntuplizer);
