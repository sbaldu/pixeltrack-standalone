#include "AlpakaCore/ScopedContext.h"
#include "AlpakaCore/config.h"
#include "AlpakaDataFormats/alpaka/PixelTrackAlpaka.h"
#include "AlpakaDataFormats/alpaka/ZVertexAlpaka.h"
#include "Framework/EDProducer.h"
#include "Framework/Event.h"
#include "Framework/EventSetup.h"
#include "Framework/PluginFactory.h"
#include "Framework/RunningAverage.h"

#include "clueVertexFinder.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class PixelVertexProducerAlpakaCLUE : public edm::EDProducer {
  public:
    explicit PixelVertexProducerAlpakaCLUE(edm::ProductRegistry& reg);
    ~PixelVertexProducerAlpakaCLUE() override = default;

  private:
    void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override;

    edm::EDGetTokenT<cms::alpakatools::Product<Queue, PixelTrackAlpaka>> tokenTrack_;
    edm::EDPutTokenT<cms::alpakatools::Product<Queue, ZVertexAlpaka>> tokenVertex_;

    // Tracking cuts before sending tracks to vertex algo
    const float m_ptMin;
    const float m_dc;
    const float m_rhoc;
    const float m_dm;
    const float m_seed_dc;
  };

  PixelVertexProducerAlpakaCLUE::PixelVertexProducerAlpakaCLUE(edm::ProductRegistry& reg)
      : tokenTrack_(reg.consumes<cms::alpakatools::Product<Queue, PixelTrackAlpaka>>()),
        tokenVertex_(reg.produces<cms::alpakatools::Product<Queue, ZVertexAlpaka>>()),
        m_ptMin{0.5f},  // 0.5 GeV
        m_dc{0.04f},
        m_rhoc{0.04f},
        m_dm{0.04f},
        m_seed_dc{0.04f} {}

  void PixelVertexProducerAlpakaCLUE::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    cms::alpakatools::Product<Queue, PixelTrackAlpaka> const& tracksWrapped = iEvent.get(tokenTrack_);
    cms::alpakatools::ScopedContextProduce<Queue> ctx{tracksWrapped};
    auto const& tracks = ctx.get(tracksWrapped);
    gpuVertexFinder::CLUEVertexProducer algo(m_dc, m_rhoc, m_dm, m_seed_dc);
    ctx.emplace(iEvent, tokenVertex_, algo.makeAsync(tracks.data(), m_ptMin, ctx.stream()));
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

DEFINE_FWK_ALPAKA_MODULE(PixelVertexProducerAlpakaCLUE);
