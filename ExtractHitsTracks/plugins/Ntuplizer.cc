#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/Common/interface/View.h"
#include "FWCore/Framework/interface/EDFilter.h" // EDAnalyzer.h
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "test/ExtractHitsTracks/interface/Ntuple.h"
#include "TTree.h"
#include<iostream> 

class Ntuplizer : public edm::EDFilter { // edm::EDAnalyzer
  
public:
  
  explicit Ntuplizer( const edm::ParameterSet& );
  ~Ntuplizer();
  
  virtual void beginRun( const edm::Run&, const edm::EventSetup& ) override;
  virtual bool filter( edm::Event&, const edm::EventSetup& ) override; // analyze(const,const)

  // Fills tree per ElectronChain object
  void fill( const edm::Event& event, const edm::EventSetup& setup );
  void readCollections( const edm::Event&, const edm::EventSetup& );
  

private: 

  edm::Service<TFileService> fs_;
  TTree* tree_;	
  Ntuple ntuple_;
  int verbose_;


  const edm::EDGetTokenT< edm::View<reco::Track> > ctfTracks_; // AOD
  edm::Handle< edm::View<reco::Track> > ctfTracksH_;


  const edm::EDGetTokenT<reco::BeamSpot> beamspot_;
  edm::Handle<reco::BeamSpot> beamspotH_;

}; 

Ntuplizer::~Ntuplizer(){}

reco::TrackPtr trk_;

Ntuplizer::Ntuplizer( const edm::ParameterSet& cfg ) 
  : tree_(nullptr),
    ntuple_(),
    verbose_(cfg.getParameter<int>("verbose")),
    ctfTracks_(consumes< edm::View<reco::Track> >(cfg.getParameter<edm::InputTag>("ctfTracks"))),
    ctfTracksH_(),
    beamspot_(consumes<reco::BeamSpot>(cfg.getParameter<edm::InputTag>("beamspot"))),
    beamspotH_()
  {
    tree_ = fs_->make<TTree>("tree","tree");
    ntuple_.link_tree(tree_);
    std::cout << "[NTuplizer::NTuplizer] Verbosity level: "<< verbose_ << std::endl;
  }
// Initialise the weights LUT to filter fake tracks 
 void Ntuplizer::beginRun( const edm::Run& run, const edm::EventSetup& es ) {
//   //@@ ?
 }


bool Ntuplizer::filter(edm::Event& event, const edm::EventSetup& setup ) {
	readCollections(event,setup);
	fill(event,setup);
	return false; 
} 


void Ntuplizer::readCollections( const edm::Event& event, const edm::EventSetup& setup ) {
	event.getByToken(ctfTracks_, ctfTracksH_);
  	event.getByToken(beamspot_, beamspotH_);
}

void Ntuplizer::fill( const edm::Event& event,
			const edm::EventSetup& setup ) {

    ntuple_.reset();
    ntuple_.fill_evt(event.id());
    //ntuple_.fill_trk( trk_, *beamspotH_ );
    tree_->Fill(); 

}
#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(Ntuplizer);

