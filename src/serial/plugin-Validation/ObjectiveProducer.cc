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
  
  auto const* hits = iEvent.get(hitToken_).view();
  auto const nHits = hits->nHits(); 
  std::set<int64_t> uniques;
  for (size_t i = 0; i < nHits; ++i) {
    if (hits->particlePT(i) > 0.899999976158)
      uniques.insert(hits->particleIndex(i));
  }
  result["simulated"] = uniques.size();

  auto const& tracks = iEvent.get(trackToken_);
  auto* soa = tracks.get();

  std::vector<uint16_t> indeces;
  for (auto &e : soa->hitIndices){
    indeces.push_back(e);
  }

  // assert(indeces.size() > 0);

  for (auto & d : soa->hitIndices.off){
    bool same = true;
    std::cout << d << '\r';
    if (d >= indeces.size()) break;
    auto particle = hits->particleIndex(indeces.at(d));
    auto pt = hits->particlePT(indeces.at(d));
    for (int i = 1; i < 3; ++i){
      if (particle != hits->particleIndex(indeces.at(d+i))){
        same = false;
        break;
      }
    }
    if (same && particle != 0){
      ++result["matching"];
    }
    ++result["reconstructed"];
  }
  std::cout<<'\n';

}

void ObjectiveProducer::endJob() {
  result["efficiency"] = result["matching"]/result["simulated"];
  result["fakes"] = result["reconstructed"]-result["matching"];
  result["fake_rate"] = result["fakes"]/result["reconstructed"];
  std::ofstream out("objectives.txt");
  for (auto const& elem : result) {
    out << elem.first << " " << elem.second << "\n";
  }
}

DEFINE_FWK_MODULE(ObjectiveProducer);
