#ifndef RecoLocalTracker_SiPixelRecHits_plugins_gpuPixelDoublets_h
#define RecoLocalTracker_SiPixelRecHits_plugins_gpuPixelDoublets_h

#include <cstdint>
#include "gpuPixelDoubletsAlgos.h"

namespace gpuPixelDoublets {

  constexpr int nPairs = 17;
  static_assert(nPairs <= CAConstants::maxNumberOfLayerPairs());

  // start constants
  // clang-format off

  // the pairs of layers of the TrackML detector
  constexpr uint8_t layerPairs[2 * nPairs] = {
    0, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10,				// barrel to left disk layers
    0, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15, 16, 16, 17,	// barrel to right disk layers
    0, 1, 1, 2, 2, 3										// barrel to barrel                                
  };

  constexpr int16_t phi0p05 = 522;  // round(521.52189...) = phi2short(0.05);
  constexpr int16_t phi0p06 = 626;  // round(625.82270...) = phi2short(0.06);
  constexpr int16_t phi0p07 = 730;  // round(730.12648...) = phi2short(0.07);

  constexpr int16_t phicuts[nPairs]{phi0p05, phi0p05, phi0p05, phi0p05};

  // clang-format on

  using CellNeighbors = CAConstants::CellNeighbors;
  using CellTracks = CAConstants::CellTracks;
  using CellNeighborsVector = CAConstants::CellNeighborsVector;
  using CellTracksVector = CAConstants::CellTracksVector;

  void initDoublets(GPUCACell::OuterHitOfCell* isOuterHitOfCell,
                    int nHits,
                    CellNeighborsVector* cellNeighbors,
                    CellNeighbors* cellNeighborsContainer,
                    CellTracksVector* cellTracks,
                    CellTracks* cellTracksContainer) {
    assert(isOuterHitOfCell);
    int first = 0;
    for (int i = first; i < nHits; i++)
      isOuterHitOfCell[i].reset();

    if (0 == first) {
      cellNeighbors->construct(CAConstants::maxNumOfActiveDoublets(), cellNeighborsContainer);
      cellTracks->construct(CAConstants::maxNumOfActiveDoublets(), cellTracksContainer);
      auto i = cellNeighbors->extend();
      assert(0 == i);
      (*cellNeighbors)[0].reset();
      i = cellTracks->extend();
      assert(0 == i);
      (*cellTracks)[0].reset();
    }
  }

  constexpr auto getDoubletsFromHistoMaxBlockSize = 64;  // for both x and y
  constexpr auto getDoubletsFromHistoMinBlocksPerMP = 16;

#ifdef __CUDACC__
  __launch_bounds__(getDoubletsFromHistoMaxBlockSize, getDoubletsFromHistoMinBlocksPerMP)
#endif
      void getDoubletsFromHisto(GPUCACell* cells,
                                uint32_t* nCells,
                                CellNeighborsVector* cellNeighbors,
                                CellTracksVector* cellTracks,
                                TrackingRecHit2DSOAView const* __restrict__ hhp,
                                GPUCACell::OuterHitOfCell* isOuterHitOfCell,
                                int nActualPairs,
				bool ideal_cond,
                                bool doClusterCut,
                                bool doZ0Cut,
                                bool doPtCut,
                                uint32_t maxNumOfDoublets,
                                int16_t* phiCuts=nullptr
				) {
    auto const& __restrict__ hh = *hhp;
    if (phiCuts==nullptr) {
        phiCuts=const_cast<int16_t*>(phicuts);
    }
    doubletsFromHisto(layerPairs,
                      nActualPairs,
                      cells,
                      nCells,
                      cellNeighbors,
                      cellTracks,
                      hh,
                      isOuterHitOfCell,
                      phicuts,
                      ideal_cond,
                      doClusterCut,
                      doZ0Cut,
                      doPtCut,
                      maxNumOfDoublets);
  }

}  // namespace gpuPixelDoublets

#endif  // RecoLocalTracker_SiPixelRecHits_plugins_gpuPixelDouplets_h
