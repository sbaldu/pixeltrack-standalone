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
  {"duplicates", 0.f},
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
  // std::map<int, int> hitsPerLayer;



  std::map<int64_t, std::tuple<float, float, float, int>> uniques;
  for (size_t i = 0; i < nHits; ++i) {
    if (hits->particlePT(i) > 0.9 && hits->particleNHits(i) > 3){
      auto pT = hits->particlePT(i);
      auto dR = hits->particledR(i);
      auto vz = hits->particleVz(i);
      auto pnHits = hits->particleNHits(i);
      auto index = hits->particleIndex(i);

      auto values = std::make_tuple(pT, dR, vz, pnHits);
      uniques.insert(std::make_pair(index, values));
      }
  }
  result["simulated"] += uniques.size();
  
  std::map<int64_t, int> recos;

    auto duplicates = 0;
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
    // auto chi2 = soa->chi2(i);
    // std::cout << "chi2: " << chi2 << std::endl;
    auto majorityIndex{indeces.at(offset)};
    int count{1};
    for (int j = 1; j < nHits; ++j){
      // std::cout << hits->detectorIndex(indeces.at(offset+j)) << '\t';
      // hitsPerLayer[hits->detectorIndex(indeces.at(offset+j))]++;
      if (hits->particleIndex(indeces.at(offset+j)) != hits->particleIndex(majorityIndex)){
        if (--count == 0){
          majorityIndex = indeces.at(offset+j);
          count = 1;
        }
        continue;
      }
      ++count;
    }
    for (int j = 0; j < nHits; ++j)
      if (hits->particleIndex(indeces.at(offset+j)) == hits->particleIndex(majorityIndex))
        ++same ;
    auto particle = hits->particleIndex(majorityIndex);
    float threshold = same/nHits;
    auto fake = true;
    if (threshold>0.75 && particle != 0){
      if (recos.find(particle) == recos.end())
        recos[particle] = 0;
      ++recos[particle];
      if (recos[particle] > 1){
        ++duplicates;
        std::cout << "Duplicate: ";
        for (int j = 0; j < nHits; ++j){
          // std::cout << hits->particleIndex(indeces.at(offset+j)) << '\t';
          std::cout << hits->detectorIndex(indeces.at(offset+j)) << '\t';
          // std::cout << hits->xGlobal(indeces.at(offset+j)) << '\t';
        }
        std::cout << std::endl;
      }
      else{
        std::cout << "Reco: ";
        for (int j = 0; j < nHits; ++j){
          std::cout << hits->particleIndex(indeces.at(offset+j)) << '\t';
        }
        std::cout << std::endl;
      }
      ++result["recotosim"];
      fake = false;
    }
    if (fake){
      std::cout << "Fake: ";
      for (int j = 0; j < nHits; ++j){
        std::cout << hits->particleIndex(indeces.at(offset+j)) << '\t';
      }
    std::cout << std::endl;
    }
    ++result["reconstructed"];
  }
  std::cout << "Duplicates: " << duplicates << std::endl;
  // std::ofstream out("simulated.csv");
  // out << "index, pT, dR, vz, nHits\n";
  // for (auto & sim: uniques){
  //     auto pInd = sim.first;
  //     auto pT = std::get<0>(sim.second);
  //     auto dR = std::get<1>(sim.second);
  //     auto vz = std::get<2>(sim.second);
  //     auto nHits = std::get<3>(sim.second);
  //     out << pInd << ","<< pT << "," << dR << "," << vz << "," << nHits << "\n";
  //  }
  // out.close();

  // out.open("simtoreco.csv");
  // out << "index, pT, dR, vz, nHits\n";
  for (auto & sim: uniques){
    // auto pInd = sim.first;
    // auto pT = std::get<0>(sim.second);
    // auto dR = std::get<1>(sim.second);
    // auto vz = std::get<2>(sim.second);
    // auto nHits = std::get<3>(sim.second);
    bool found_simtoreco = recos.find(sim.first) != recos.end();
    if (found_simtoreco) {
      // out << pInd << ","<< pT << "," << dR << "," << vz << "," << nHits << "\n";
      ++result["simtoreco"];
    }
  }
}


void ObjectiveProducer::endJob() {
  result["efficiency"] = result["simtoreco"]/result["simulated"];
  result["fakes"] = result["reconstructed"]-result["recotosim"];
  result["duplicates"] = result["recotosim"] - result["simtoreco"];
  result["fake_rate"] = (result["fakes"]+result["duplicates"])/result["reconstructed"];
  std::ofstream out("objectives.txt");
  for (auto const& elem : result) {
    out << elem.first << " " << elem.second << "\n";
  }
}

DEFINE_FWK_MODULE(ObjectiveProducer);
