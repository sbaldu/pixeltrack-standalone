
#pragma once

#include "AlpakaDataFormats/alpaka/PixelTrackAlpaka.h"
#include "AlpakaDataFormats/alpaka/ZVertexAlpaka.h"
#include "plugin-PixelVertexFinding/alpaka/VertexWorkspace.h"
#include <alpaka/alpaka.hpp>

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  namespace gpuVertexFinder {

    struct loadTracks {
      template <typename TAcc>
      ALPAKA_FN_ACC void operator()(
          const TAcc& acc, TkSoA const* ptracks, ZVertexSoA* soa, WorkSpace* pws, float ptMin) const {
        ALPAKA_ASSERT_ACC(ptracks);
        ALPAKA_ASSERT_ACC(soa);
        auto const& tracks = *ptracks;
        auto const& fit = tracks.stateAtBS;
        auto const* quality = tracks.qualityData();

        for (auto idx : alpaka::uniformElement(acc, ptracks->m_ntracks)) {
          auto nHits = tracks.nHits(idx);
          if (nHits == 0)
            return;  // this is a guard: maybe we need to move to nTracks...

          // initialize soa...
          soa->idv[idx] = -1;

          if (nHits < 4)
            return;  // no triplets
          if (quality[idx] != trackQuality::loose)
            return;

          auto pt = tracks.pt(idx);

          if (pt < ptMin)
            return;

          auto& data = *pws;
          auto it = alpaka::atomicAdd(acc, &data.ntrks, 1u, alpaka::hierarchy::Blocks{});
          data.itrk[it] = idx;
          data.zt[it] = tracks.zip(idx);
          data.ezt2[it] = fit.covariance(idx)(14);
          data.ptt2[it] = pt * pt;
        }
      }
    };

  }  // namespace gpuVertexFinder
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
