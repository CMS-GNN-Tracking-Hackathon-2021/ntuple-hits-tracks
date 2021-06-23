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
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "test/ExtractHitsTracks/interface/Ntuple.h"
#include "TTree.h"
#include<iostream> 

////////////////////////////////////////////////////////////////////////////////
//
class Ntuplizer : public edm::EDFilter {
  
public:
  
  explicit Ntuplizer( const edm::ParameterSet& );
  ~Ntuplizer();
  
  virtual void beginRun( const edm::Run&, const edm::EventSetup& ) override;
  virtual bool filter( edm::Event&, const edm::EventSetup& ) override;
  
  void fill( const edm::Event& event, const edm::EventSetup& setup );
  void readCollections( const edm::Event&, const edm::EventSetup& );
  
private: 
  
  edm::Service<TFileService> fs_;
  TTree* tree_;	
  Ntuple ntuple_;
  int verbose_;
  
  const edm::EDGetTokenT<reco::TrackCollection> ctfTracksToken_;
  reco::TrackCollection const* ctfTracks_;

  const edm::EDGetTokenT< edmNew::DetSetVector<SiPixelCluster> > siPixelClustersToken_;
  edmNew::DetSetVector<SiPixelCluster> const* siPixelClusters_;

  const edm::EDGetTokenT< edmNew::DetSetVector<Phase2TrackerCluster1D> > siPhase2ClustersToken_;
  edmNew::DetSetVector<Phase2TrackerCluster1D> const* siPhase2Clusters_;

}; 

////////////////////////////////////////////////////////////////////////////////
//
Ntuplizer::~Ntuplizer(){}

////////////////////////////////////////////////////////////////////////////////
//
Ntuplizer::Ntuplizer( const edm::ParameterSet& cfg ) 
  : tree_(nullptr),
    ntuple_(),
    verbose_(cfg.getParameter<int>("verbose")),
    ctfTracksToken_(consumes<reco::TrackCollection>(cfg.getParameter<edm::InputTag>("ctfTracks"))),
    ctfTracks_(),
    siPixelClustersToken_(consumes< edmNew::DetSetVector<SiPixelCluster> >(cfg.getParameter<edm::InputTag>("siPixelClusters"))),
    siPixelClusters_(),
    siPhase2ClustersToken_(consumes< edmNew::DetSetVector<Phase2TrackerCluster1D> >(cfg.getParameter<edm::InputTag>("siPhase2Clusters"))),
    siPhase2Clusters_()
  {
    tree_ = fs_->make<TTree>("tree","tree");
    ntuple_.link_tree(tree_);
    std::cout << "[Ntuplizer::Ntuplizer] Verbosity level: "<< verbose_ << std::endl;
  }

////////////////////////////////////////////////////////////////////////////////
//
void Ntuplizer::beginRun( const edm::Run& run, const edm::EventSetup& setup ) {
}

////////////////////////////////////////////////////////////////////////////////
//
bool Ntuplizer::filter(edm::Event& event, const edm::EventSetup& setup ) {
  readCollections(event,setup);
  fill(event,setup);
  return false; 
} 

////////////////////////////////////////////////////////////////////////////////
//
void Ntuplizer::readCollections( const edm::Event& event, const edm::EventSetup& setup ) {
  ctfTracks_ = &event.get(ctfTracksToken_);
  siPixelClusters_ = &event.get(siPixelClustersToken_);
  siPhase2Clusters_ = &event.get(siPhase2ClustersToken_);
  if (verbose_>0) {
    std::cout << "[Ntuplizer::readCollections]"
	      << " ctfTracks_->size(): " << ctfTracks_->size()
	      << " siPixelClusters_->size(): " << siPixelClusters_->size()
	      << " siPhase2Clusters_->size(): " << siPhase2Clusters_->size()
	      << std::endl;
  }
}

////////////////////////////////////////////////////////////////////////////////
//
void Ntuplizer::fill( const edm::Event& event,
		      const edm::EventSetup& setup ) {
  ntuple_.reset();
  ntuple_.fill_evt(event.id());
  ntuple_.fill_clu(siPixelClusters_,siPhase2Clusters_);
  tree_->Fill(); 
}

////////////////////////////////////////////////////////////////////////////////
//
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(Ntuplizer);
