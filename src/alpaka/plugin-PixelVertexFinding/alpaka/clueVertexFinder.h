
#pragma once

#include "AlpakaCore/config.h"
#include "AlpakaDataFormats/alpaka/PixelTrackAlpaka.h"
#include "AlpakaDataFormats/alpaka/ZVertexAlpaka.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  namespace gpuVertexFinder {

    class CLUEVertexProducer {
    public:
      using TkSoA = pixelTrack::TrackSoA;

      CLUEVertexProducer(float dc, float rhoc, float dm, float seed_dc)
          : m_dc{dc}, m_rhoc{rhoc}, m_dm{dm}, m_seed_dc{seed_dc} {}

      ~CLUEVertexProducer() = default;

      ZVertexAlpaka makeAsync(TkSoA const* tksoa, float ptMin, Queue& queue) const;

    private:
      float m_dc;
      float m_rhoc;
      float m_dm;
      float m_seed_dc;
    };

  }  // namespace gpuVertexFinder
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
