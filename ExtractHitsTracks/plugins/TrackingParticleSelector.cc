#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

////////////////////////////////////////////////////////////////////////////////
//
class TrackingParticleSelector : public edm::stream::EDProducer<> {
  
public:
  
  explicit TrackingParticleSelector( const edm::ParameterSet& );
  ~TrackingParticleSelector();
  void produce(edm::Event&, const edm::EventSetup&) override;
  //static void fillDescriptions(edm::ConfigurationDescriptions&);
  
private: 

  const edm::EDGetTokenT<TrackingParticleCollection> tpToken_;
  const edm::EDGetTokenT<edm::Association<reco::GenParticleCollection> > prunedToken_;

}; 

////////////////////////////////////////////////////////////////////////////////
//
TrackingParticleSelector::~TrackingParticleSelector(){}

////////////////////////////////////////////////////////////////////////////////
//
TrackingParticleSelector::TrackingParticleSelector( const edm::ParameterSet& cfg ) : 
  tpToken_(consumes<TrackingParticleCollection>(cfg.getParameter<edm::InputTag>("trackingParticles"))),
  prunedToken_(consumes<edm::Association<reco::GenParticleCollection> >(cfg.getParameter<edm::InputTag>("prunedGenParticles")))
{
  produces<TrackingParticleCollection>();
}

////////////////////////////////////////////////////////////////////////////////
//
void TrackingParticleSelector::produce(edm::Event& event, const edm::EventSetup& setup) {
  const auto& tpcoll = event.get(tpToken_);
  const auto& pruned = event.get(prunedToken_);
  auto out = std::make_unique<TrackingParticleCollection>();
  out->reserve(tpcoll.size());
  for ( auto&& tp : tpcoll ) {
    //auto& genParticles = tp.genParticles();
    for ( const auto& gen : tp.genParticles() ) {
      if ( pruned[gen].isNonnull() ) {
	out->emplace_back(tp);
	break;
      }
    }
  }
  event.put(std::move(out));
}    

//auto out = std::make_unique<GenParticleCollection>();
//GenParticleRefProd outRef = evt.getRefBeforePut<GenParticleCollection>();
//out->reserve(counter);
//
//edm::OrphanHandle<reco::GenParticleCollection> oh = evt.put(std::move(out));
//auto orig2new = std::make_unique<edm::Association<reco::GenParticleCollection>>(oh);
//edm::Association<reco::GenParticleCollection>::Filler orig2newFiller(*orig2new);
//orig2newFiller.insert(src, flags_.begin(), flags_.end());
//orig2newFiller.fill();
//evt.put(std::move(orig2new));


//////////////////////////////////////////////////////////////////////////////////
////
//void TrackingParticleSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
//  edm::ParameterSetDescription desc;
//  desc.add<edm::InputTag>("trackingParticles", edm::InputTag("mix","MergedTrackTruth"));
//  desc.add<edm::InputTag>("prunedGenParticles", edm::InputTag("prunedGenParticles"));
//  descriptions.add("TrackingParticleSelectorDefault", desc);
//}

////////////////////////////////////////////////////////////////////////////////
//
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TrackingParticleSelector);
