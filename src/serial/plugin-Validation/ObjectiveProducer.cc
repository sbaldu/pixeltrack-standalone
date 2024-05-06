#include "CUDADataFormats/PixelTrackHeterogeneous.h"
#include "CUDADataFormats/SiPixelClustersSoA.h"
#include "CUDADataFormats/SiPixelDigisSoA.h"
#include "CUDADataFormats/TrackingRecHit2DHeterogeneous.h"
#include "CUDADataFormats/ZVertexHeterogeneous.h"
#include "Framework/EventSetup.h"
#include "Framework/Event.h"
#include "Framework/PluginFactory.h"
#include "Framework/EDProducer.h"

#include "SimpleAtomicHisto.h"

#include <map>
#include <fstream>

class ObjectiveProducer : public edm::EDProducer {
public:
  explicit ObjectiveProducer(edm::ProductRegistry& reg);

private:
  void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override;
  void endJob() override;

  edm::EDGetTokenT<TrackingRecHit2DCPU> hitToken_;
  edm::EDGetTokenT<PixelTrackHeterogeneous> trackToken_;

  static std::map<std::string, float_t> result;
  uint16_t simulatedTracks_;
};

std::map<std::string, float_t> ObjectiveProducer::result = {
  {"reconstructed", 0.f},
  {"simulated", 0.f},
  {"fakes", 0.f},
  {"matching", 0.f},
  {"efficiency", 0.f},
  {"fake_rate", 0.f}
};

ObjectiveProducer::ObjectiveProducer(edm::ProductRegistry& reg)
    : hitToken_(reg.consumes<TrackingRecHit2DCPU>()),
      trackToken_(reg.consumes<PixelTrackHeterogeneous>()) {}

void ObjectiveProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // TODO: find a way to do this once per job
  auto const* hits = iEvent.get(hitToken_).view();
  auto const nHits = hits->nHits(); 
  std::set<uint16_t> uniques;
  for (uint32_t i = 0; i < nHits; ++i) {
    uniques.insert(hits->particleIndex(i));
  }
  simulatedTracks_ = uniques.size();

  // TODO: loop through reconstructed
  auto const& tracks = iEvent.get(trackToken_);
    std::cout<<tracks->hitIndices.off<<'\n';
  
    // TODO: count tracks with >= 0.75 hits with same p_id  

  // auto const* hits = iEvent.get(hitToken_).view();

  // auto const nHits = hits->nHits();
  // histos["hit_n"].fill(nHits);
  // for (uint32_t i = 0; i < nHits; ++i) {
  //   histos["hit_lx"].fill(hits->xLocal(i));
  //   histos["hit_ly"].fill(hits->yLocal(i));
  //   histos["hit_lex"].fill(hits->xerrLocal(i));
  //   histos["hit_ley"].fill(hits->yerrLocal(i));
  //   histos["hit_gx"].fill(hits->xGlobal(i));
  //   histos["hit_gy"].fill(hits->yGlobal(i));
  //   histos["hit_gz"].fill(hits->zGlobal(i));
  //   histos["hit_gr"].fill(hits->rGlobal(i));
  //   histos["hit_charge"].fill(hits->charge(i));
  //   histos["hit_sizex"].fill(hits->clusterSizeX(i));
  //   histos["hit_sizey"].fill(hits->clusterSizeY(i));
  // }

  // {
  //   auto const& tracks = iEvent.get(trackToken_);

  //   int nTracks = 0;
  //   for (int i = 0; i < tracks->stride(); ++i) {
  //     if (tracks->nHits(i) > 0 and tracks->quality(i) >= trackQuality::loose) {
  //       ++nTracks;
  //       histos["track_nhits"].fill(tracks->nHits(i));
  //       histos["track_chi2"].fill(tracks->chi2(i));
  //       histos["track_pt"].fill(tracks->pt(i));
  //       histos["track_eta"].fill(tracks->eta(i));
  //       histos["track_phi"].fill(tracks->phi(i));
  //       histos["track_tip"].fill(tracks->tip(i));
  //       histos["track_tip_zoom"].fill(tracks->tip(i));
  //       histos["track_zip"].fill(tracks->zip(i));
  //       histos["track_zip_zoom"].fill(tracks->zip(i));
  //       histos["track_quality"].fill(tracks->quality(i));
  //     }
  //   }

  //   histos["track_n"].fill(nTracks);
  // }

}

void ObjectiveProducer::endJob() {
  result["simulated"] = simulatedTracks_;
  result["efficiency"] = result["matching"]/result["simulated"];
  result["fakes"] = result["reconstructed"]-result["matching"];
  result["fake_rate"] = result["fakes"]/result["reconstructed"];
  std::ofstream out("objectives.txt");
  for (auto const& elem : result) {
    out << elem.first << " " << elem.second << "\n";
  }
}

DEFINE_FWK_MODULE(ObjectiveProducer);
