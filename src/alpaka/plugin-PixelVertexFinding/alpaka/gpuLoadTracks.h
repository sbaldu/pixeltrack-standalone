#ifndef plugin_PixelVertexFinding_alpaka_LoadTracks_h
#define plugin_PixelVertexFinding_alpaka_LoadTracks_h

#include "AlpakaDataFormats/alpaka/PixelTrackAlpaka.h"
#include "AlpakaDataFormats/alpaka/ZVertexAlpaka.h"
#include "VertexWorkspace.h"
#include <alpaka/alpaka.hpp>

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  namespace gpuVertexFinder {

    struct LoadTracks {
      template <typename TAcc>
      ALPAKA_FN_ACC void operator()(
          const TAcc &acc, pixelTrack::TrackSoA const *ptracks, ZVertexSoA *soa, WorkSpaceView pws, float ptMin) const {
        ALPAKA_ASSERT_ACC(ptracks);
        ALPAKA_ASSERT_ACC(soa);
        auto const &tracks = *ptracks;
        auto const &fit = tracks.stateAtBS;
        auto const *quality = tracks.qualityData();

        cms::alpakatools::for_each_element_in_grid_strided(acc, pixelTrack::TrackSoA::stride(), [&](uint32_t idx) {
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

          auto it = alpaka::atomicAdd(acc, pws.ntrks, 1u, alpaka::hierarchy::Blocks{});
          pws.itrk[it] = idx;
          pws.zt[it] = tracks.zip(idx);
          pws.ezt2[it] = fit.covariance(idx)(14);
          pws.ptt2[it] = pt * pt;
        });
      }
    };

  }  // namespace gpuVertexFinder
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif
