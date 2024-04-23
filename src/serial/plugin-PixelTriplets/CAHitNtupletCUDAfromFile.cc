#include "Framework/EventSetup.h"
#include "Framework/Event.h"
#include "Framework/PluginFactory.h"
#include "Framework/EDProducer.h"
#include "Framework/RunningAverage.h"

#include "CAHitNtupletGeneratorOnGPU.h"
#include "CUDADataFormats/PixelTrackHeterogeneous.h"
#include "CUDADataFormats/TrackingRecHit2DHeterogeneous.h"

class CAHitNtupletCUDAfromFile : public edm::EDProducer {
public:
  explicit CAHitNtupletCUDAfromFile(edm::ProductRegistry& reg, const std::string& filename);
  ~CAHitNtupletCUDAfromFile() override = default;

private:
  void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override;

  edm::EDGetTokenT<TrackingRecHit2DCPU> tokenHitCPU_;
  edm::EDPutTokenT<PixelTrackHeterogeneous> tokenTrackCPU_;

  CAHitNtupletGeneratorOnGPU gpuAlgo_;
};

CAHitNtupletCUDAfromFile::CAHitNtupletCUDAfromFile(edm::ProductRegistry& reg, const std::string& filename)
    : tokenHitCPU_{reg.consumes<TrackingRecHit2DCPU>()},
      tokenTrackCPU_{reg.produces<PixelTrackHeterogeneous>()},
      gpuAlgo_(reg, filename) {}

void CAHitNtupletCUDAfromFile::produce(edm::Event& iEvent, const edm::EventSetup& es) {
  auto bf = 0.0114256972711507;  // 1/fieldInGeV

  auto const& hits = iEvent.get(tokenHitCPU_);

  iEvent.emplace(tokenTrackCPU_, gpuAlgo_.makeTuples(hits, bf));
}

DEFINE_FWK_MODULE(CAHitNtupletCUDAfromFile);
