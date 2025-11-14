
#pragma once

#include "AlpakaCore/config.h"
#include "AlpakaDataFormats/alpaka/PixelTrackAlpaka.h"
#include "AlpakaDataFormats/alpaka/ZVertexAlpaka.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  namespace clueVertexFinder {

    class Producer {
    public:
      Producer(float dc, float rhoc, float dm, float seed_dc) : m_dc{dc}, m_rhoc{rhoc}, m_dm{dm}, m_seed_dc{seed_dc} {}

      ~Producer() = default;

      ZVertexAlpaka makeAsync(TkSoA const* tksoa, float ptMin, Queue& queue) const;

    private:
      float m_dc;
      float m_rhoc;
      float m_dm;
      float m_seed_dc;
    };

  }  // namespace clueVertexFinder
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
