#include "Framework/EventSetup.h"
#include "Framework/Event.h"
#include "Framework/PluginFactory.h"
#include "Framework/EDProducer.h"
#include "Framework/RunningAverage.h"

#include "CAHitNtupletGeneratorOnGPU.h"
#include "CUDADataFormats/PixelTrackHeterogeneous.h"
#include "CUDADataFormats/TrackingRecHit2DHeterogeneous.h"
#include "plugin-SiPixelRecHits/PixelRecHits.h"

class TrackMLProducer : public edm::EDProducer {
public:
  explicit TrackMLProducer(edm::ProductRegistry& reg);
  ~TrackMLProducer() override = default;

private:
  void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override;

  pixelgpudetails::PixelRecHitGPUKernel algo_;
  CAHitNtupletGeneratorOnGPU gpuAlgo_;
  edm::EDPutTokenT<TrackingRecHit2DCPU> tokenHitCPU_;
};

TrackMLProducer::TrackMLProducer(edm::ProductRegistry& reg)
    : algo_(), gpuAlgo_(reg), tokenHitCPU_(reg.produces<TrackingRecHit2DCPU>()) {}

void TrackMLProducer::produce(edm::Event& iEvent, const edm::EventSetup& es) {
  iEvent.emplace(tokenHitCPU_, algo_.makeHits(iEvent.eventID()));
}

DEFINE_FWK_MODULE(TrackMLProducer);
