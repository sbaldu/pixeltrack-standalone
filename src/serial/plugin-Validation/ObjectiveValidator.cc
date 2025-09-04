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

class ObjectiveValidator : public edm::EDProducer {
public:
  explicit ObjectiveValidator(edm::ProductRegistry& reg);

private:
  void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override;
  void endJob() override;

  edm::EDGetTokenT<TrackingRecHit2DCPU> hitToken_;
  edm::EDGetTokenT<PixelTrackHeterogeneous> trackToken_;

  static std::map<std::string, std::unordered_map<uint16_t,std::vector<float_t>>> result;
};

std::map<std::string, std::unordered_map<uint16_t,std::vector<float_t>>> ObjectiveValidator::result = {{"pt_sim", {}},
                                                                                          {"pt_sim2reco", {}},
                                                                                          {"pt_reco", {}},
                                                                                          {"pt_reco2sim", {}},
                                                                                          {"eta_sim", {}},
                                                                                          {"eta_sim2reco", {}},
                                                                                          {"eta_reco", {}},
                                                                                          {"eta_reco2sim", {}},
                                                                                          {"phi_sim", {}},
                                                                                          {"phi_sim2reco", {}},
                                                                                          {"phi_reco", {}},
                                                                                          {"phi_reco2sim", {}}};

ObjectiveValidator::ObjectiveValidator(edm::ProductRegistry& reg)
    : hitToken_(reg.consumes<TrackingRecHit2DCPU>()),
      trackToken_(reg.consumes<PixelTrackHeterogeneous>()) {}

void ObjectiveValidator::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  auto const* hits = iEvent.get(hitToken_).view();
  auto const nHits = hits->nHits();

  auto const& tracks = iEvent.get(trackToken_);
  auto* soa = tracks.get();

  // Create indices vector but reserve exact size to avoid reallocations
  std::vector<uint16_t> indeces;
  indeces.reserve(soa->hitIndices.size());
  for (const auto& e : soa->hitIndices) {
    indeces.push_back(e);
  }
  // Create pt vectors
  std::vector<float> pt_sim;
  std::vector<float> pt_sim2reco;
  std::vector<float> pt_reco;
  std::vector<float> pt_reco2sim;
  
  // Create eta vectors
  std::vector<float> eta_sim;
  std::vector<float> eta_sim2reco;
  std::vector<float> eta_reco;
  std::vector<float> eta_reco2sim;

  // Create phi vectors
  std::vector<float> phi_sim;
  std::vector<float> phi_sim2reco;
  std::vector<float> phi_reco;
  std::vector<float> phi_reco2sim;
  

  std::unordered_map<int64_t, std::tuple<float, float, float, int>> uniques;
  uniques.reserve(nHits);
  
  // Cache particle property lookups
  for (size_t i = 0; i < nHits; ++i) {
    auto pT = hits->particlePT(i);
    auto pnHits = hits->particleNHits(i);
    if (pT > 0.9f && pnHits > 3) {
      auto index = hits->particleIndex(i);
      // Use emplace with hint for better performance
      uniques.emplace(index, std::make_tuple(pT, hits->particleEta(i), hits->particlePhi(i), pnHits));
    }
  }

  std::unordered_map<int64_t, int> recos;
  recos.reserve(uniques.size());

  const auto stride = soa->stride();
  
  for (int i = 0; i < stride; ++i) {
    auto nHits = soa->nHits(i);
    float pt = soa->pt(i);
    float eta = soa->eta(i);
    float phi = soa->phi(i);
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
      pt_reco2sim.push_back(pt);
      eta_reco2sim.push_back(eta);
      phi_reco2sim.push_back(phi);
    }
    pt_reco.push_back(pt);
    eta_reco.push_back(eta);
    phi_reco.push_back(phi);
  }

  // Count simtoreco more efficiently
  for (const auto& [particle_id, _] : uniques) {
    if (recos.count(particle_id)) {
      pt_sim2reco.push_back(std::get<0>(uniques[particle_id]));
      eta_sim2reco.push_back(std::get<1>(uniques[particle_id]));
      phi_sim2reco.push_back(std::get<2>(uniques[particle_id]));
    }
    pt_sim.push_back(std::get<0>(uniques[particle_id]));
    eta_sim.push_back(std::get<1>(uniques[particle_id]));
    phi_sim.push_back(std::get<2>(uniques[particle_id]));
  }

  // Update global vectors at the end
  auto event_id = iEvent.eventID();
  result["pt_sim"][event_id] = pt_sim;
  result["pt_sim2reco"][event_id] = pt_sim2reco;
  result["pt_reco"][event_id] = pt_reco;
  result["pt_reco2sim"][event_id] = pt_reco2sim;
  result["eta_sim"][event_id] = eta_sim;
  result["eta_sim2reco"][event_id] = eta_sim2reco;
  result["eta_reco"][event_id] = eta_reco;
  result["eta_reco2sim"][event_id] = eta_reco2sim;
  result["phi_sim"][event_id] = phi_sim;
  result["phi_sim2reco"][event_id] = phi_sim2reco;
  result["phi_reco"][event_id] = phi_reco;
  result["phi_reco2sim"][event_id] = phi_reco2sim;

}

void ObjectiveValidator::endJob() {
  auto dir = std::string("optimization/");
  for (const auto& [key, events] : result) {
    std::cout << "Writing " << key << ".csv" << '\n';
    std::ofstream out(dir + key + ".csv");
    for (const auto& [id, values] : events) {
      for (const auto& val : values) {
        out << val << "\n";
      }
    }
  }
}

DEFINE_FWK_MODULE(ObjectiveValidator);

