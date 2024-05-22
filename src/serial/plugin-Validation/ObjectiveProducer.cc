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
  {"recotosim", 0.f},
  {"simtoreco", 0.f},
  {"efficiency", 0.f},
  {"fake_rate", 0.f}
};

ObjectiveProducer::ObjectiveProducer(edm::ProductRegistry& reg)
    : hitToken_(reg.consumes<TrackingRecHit2DCPU>()),
      trackToken_(reg.consumes<PixelTrackHeterogeneous>()) {}

void ObjectiveProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  auto const* hits = iEvent.get(hitToken_).view();
  auto const nHits = hits->nHits(); 

  auto const& tracks = iEvent.get(trackToken_);
  auto* soa = tracks.get();

  std::vector<uint16_t> indeces;
  for (auto &e : soa->hitIndices){
    indeces.push_back(e);
  }

  // assert(indeces.size() > 0);
  std::map<int, int> hitsPerLayer;



  std::set<int64_t> uniques;
  for (size_t i = 0; i < nHits; ++i) {
    if (hits->particlePT(i) > 0.899999976158)
      uniques.insert(hits->particleIndex(i));
  }
  result["simulated"] = uniques.size();
  
  std::map<int64_t, int> recos;

  for (int i = 0; i < soa->stride(); ++i){
    auto nHits = soa->nHits(i);
    if (nHits == 0)
      break;
    float same = 1.;
    auto offset = soa->hitIndices.off[i];
    if (offset >= indeces.size()) break;




    auto quality = soa->quality(i);
    if (quality == trackQuality::bad)
      continue;
    auto particle = hits->particleIndex(indeces.at(offset));
    auto pt = hits->particlePT(indeces.at(offset));
    for (int j = 1; j < nHits; ++j){
      hitsPerLayer[hits->detectorIndex(indeces.at(offset))]++;
      if (particle != hits->particleIndex(indeces.at(offset+j))){
        continue;
      }
      same++;
    }
    float threshold = same/nHits;
    if (threshold>=0.75 && particle != 0){
      if (recos.find(particle) == recos.end())
        recos[particle] = 1;
      recos[particle]++;      
      ++result["recotosim"];
    }
    ++result["reconstructed"];
  }

  for (auto & sim: uniques){
    if (recos.find(sim) != recos.end()){
      ++result["simtoreco"];
    }
  }
  
  // print the hitsPerLayer map
  for (auto const& elem : hitsPerLayer) {
    std::cout << elem.first << " " << elem.second << "\n";
  }

  // for (auto const& elem : recos) {
  //   if (elem.second > 1){
  //     std::cout << elem.first << " " << elem.second << "\n";
  //   }
  // }

}


void ObjectiveProducer::endJob() {
  result["efficiency"] = result["simtoreco"]/result["simulated"];
  result["fakes"] = result["reconstructed"]-result["recotosim"];
  result["fake_rate"] = result["fakes"]/result["reconstructed"];
  std::ofstream out("objectives.txt");
  for (auto const& elem : result) {
    out << elem.first << " " << elem.second << "\n";
  }
}

DEFINE_FWK_MODULE(ObjectiveProducer);
