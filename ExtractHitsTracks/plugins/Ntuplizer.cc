#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/Common/interface/View.h"
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
  
  const edm::EDGetTokenT< edm::View<reco::Track> > ctfTracks_;
  edm::Handle< edm::View<reco::Track> > ctfTracksH_;
  
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
    ctfTracks_(consumes< edm::View<reco::Track> >(cfg.getParameter<edm::InputTag>("ctfTracks"))),
    ctfTracksH_()
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
  event.getByToken(ctfTracks_, ctfTracksH_);
}

////////////////////////////////////////////////////////////////////////////////
//
void Ntuplizer::fill( const edm::Event& event,
		      const edm::EventSetup& setup ) {
  ntuple_.reset();
  ntuple_.fill_evt(event.id());
  tree_->Fill(); 
}

////////////////////////////////////////////////////////////////////////////////
//
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(Ntuplizer);
