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

std::map<std::string, float_t> ObjectiveProducer::result = {{"reconstructed", 0.f},
                                                            {"simulated", 0.f},
                                                            {"fakes", 0.f},
                                                            {"duplicates", 0.f},
                                                            {"recotosim", 0.f},
                                                            {"simtoreco", 0.f},
                                                            {"efficiency", 0.f},
                                                            {"fake_rate", 0.f}};

ObjectiveProducer::ObjectiveProducer(edm::ProductRegistry& reg)
    : hitToken_(reg.consumes<TrackingRecHit2DCPU>()),
      trackToken_(reg.consumes<PixelTrackHeterogeneous>()) {}

void ObjectiveProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  auto const* hits = iEvent.get(hitToken_).view();
  auto const nHits = hits->nHits();

  auto const& tracks = iEvent.get(trackToken_);
  auto* soa = tracks.get();

  // Cache frequently accessed static counters to avoid string lookups
  float reconstructed_count = 0.f;
  float recotosim_count = 0.f;
  float simtoreco_count = 0.f;

  // Create indices vector but reserve exact size to avoid reallocations
  std::vector<uint16_t> indeces;
  indeces.reserve(soa->hitIndices.size());
  for (const auto& e : soa->hitIndices) {
    indeces.push_back(e);
  }

  std::unordered_map<int64_t, std::tuple<float, float, float, int>> uniques;
  uniques.reserve(nHits);
  
  // Cache particle property lookups
  for (size_t i = 0; i < nHits; ++i) {
    auto pT = hits->particlePT(i);
    auto pnHits = hits->particleNHits(i);
    if (pT > 0.9f && pnHits > 3) {
      auto index = hits->particleIndex(i);
      // Use emplace with hint for better performance
      uniques.emplace(index, std::make_tuple(pT, hits->particledR(i), hits->particleVz(i), pnHits));
    }
  }
  result["simulated"] += uniques.size();

  std::unordered_map<int64_t, int> recos;
  recos.reserve(uniques.size());
  std::string dir = "/eos/user/s/srossiti/track-ml/";

  // auto dump_file = std::ofstream(dir + "reco_dump.csv");
  // dump_file << "type, particle_id,"
  //           << "hit1_id, hit1_x, hit1_y, hit1_z, hit1_r, hit1_module,"
  //           << "hit2_id, hit2_x, hit2_y, hit2_z, hit2_r, hit2_module,"
  //           << "hit3_id, hit3_x, hit3_y, hit3_z, hit3_r, hit3_module,"
  //           << "hit4_id, hit4_x, hit4_y, hit4_z, hit4_r, hit4_module,"
  //           << "hit5_id, hit5_x, hit5_y, hit5_z, hit5_r, hit5_module,"
  //           << "hit6_id, hit6_x, hit6_y, hit6_z, hit6_r, hit6_module,"
  //           << "hit7_id, hit7_x, hit7_y, hit7_z, hit7_r, hit7_module,"
  //           <<'\n';

  auto duplicates = 0;
  const auto stride = soa->stride();
  
  for (int i = 0; i < stride; ++i) {
    auto nHits = soa->nHits(i);
    if (nHits == 0)
      break;
    auto offset = soa->hitIndices.off[i];
    if (offset >= indeces.size())
      break;
    auto quality = soa->quality(i);
    if (quality == trackQuality::bad || quality == trackQuality::dup)
      continue;
      
    // Find majority particle index and count matches in one pass
    std::unordered_map<int64_t, int> particle_counts;
    particle_counts.reserve(nHits);
    for (int j = 0; j < nHits; ++j) {
      // Use bounds-checked access but avoid repeated .at() calls
      int64_t pidx = hits->particleIndex(indeces[offset + j]);
      ++particle_counts[pidx];
    }
    
    // Find majority more efficiently
    int64_t majorityIndex = 0;
    int max_count = 0;
    for (const auto& [particle_id, count] : particle_counts) {
      if (count > max_count) {
        majorityIndex = particle_id;
        max_count = count;
      }
    }
    
    float threshold = static_cast<float>(max_count) / nHits;
    if (threshold > 0.75f && majorityIndex != 0) {
      // Use insert for more efficient map operations
      auto [it, inserted] = recos.insert({majorityIndex, 0});
      ++it->second;
      if (it->second > 1) {
        ++duplicates;
      }
      ++recotosim_count;
    }
    ++reconstructed_count;
  }
  // std::cout << "Duplicates: " << duplicates << std::endl;
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

  // Count simtoreco more efficiently
  for (const auto& [particle_id, _] : uniques) {
    if (recos.count(particle_id)) {
      ++simtoreco_count;
    }
  }
  
  // Update global counters at the end
  result["reconstructed"] += reconstructed_count;
  result["recotosim"] += recotosim_count;
  result["simtoreco"] += simtoreco_count;
  result["duplicates"] += duplicates;
}

void ObjectiveProducer::endJob() {
  result["efficiency"] = result["simtoreco"] / result["simulated"];
  result["fakes"] = result["reconstructed"] - result["recotosim"];
  result["fake+duplicates_rate"] = (result["fakes"] + result["duplicates"]) / result["reconstructed"];
  result["fake_rate"] = result["fakes"] / result["reconstructed"];
  std::cout << "Writing objectives.txt" << '\n';
  std::ofstream out("objectives.txt");
  for (auto const& elem : result) {
    out << elem.first << " " << elem.second << "\n";
  }
}

DEFINE_FWK_MODULE(ObjectiveProducer);

