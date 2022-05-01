#ifndef RecoLocalTracker_SiPixelRecHits_plugins_gpuPixelDoublets_h
#define RecoLocalTracker_SiPixelRecHits_plugins_gpuPixelDoublets_h

#include "gpuPixelDoubletsAlgos.h"

namespace gpuPixelDoublets {

  constexpr int nPairs = 13 + 2 + 4;    // Devo cambiare npairs
  static_assert(nPairs <= CAConstants::maxNumbedddrOfLayerPairs());   // devo cambiare pure questo (20)

  // start constants
  // clang-format off

  /*
  constexpr uint8_t layerPairs[2 * nPairs] = {
      0, 1, 0, 4, 0, 7,              // BPIX1 (3)
      1, 2, 1, 4, 1, 7,              // BPIX2 (5)
      4, 5, 7, 8,                    // FPIX1 (8)
      2, 3, 2, 4, 2, 7, 5, 6, 8, 9,  // BPIX3 & FPIX2 (13)
      0, 2, 1, 3,                    // Jumping Barrel (15)
      0, 5, 0, 8,                    // Jumping Forward (BPIX1,FPIX2)
      4, 6, 7, 9                     // Jumping Forward (19)
  };
  */
  constexpr uint8_t layerPairs[2 * nPairs] = {
    0, 1, 1, 2, 2, 3                                // BV8
    18, 19, 19, 20, 20, 21                          // BV13
    34, 35                                          // BV17
    4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10             // DV7
    11, 12, 12, 13, 13, 14, 14, 15, 15, 16, 16, 17  // DV9
    22, 23, 23, 24, 24, 25, 25, 26, 26, 27          // DV12
    28, 29, 29, 30, 30, 31, 31, 32, 32, 33          // DV14
    36, 37, 37, 38, 38, 39, 39, 40, 40, 41          // DV16
                                                    // DV18
    0, 4                                            // BV8DV7
    4, 27                                           // DV7DV12
    0, 11, 1, 11, 2, 11                             // BV8DV9
    1, 27, 3, 22, 0, 27, 3, 24, 2, 25               // BV8DV12
    3, 18                                           // BV8BV13
    17, 33, 16, 31, 14, 29, 12, 28, 15, 30          // DV9DV14
                                                    // BV13DV12
                                                    // BV13DV14
                                                    // BV13DV16
                                                    // BV13BV17
                                                    // BV13DV18
                                                    // DV14DV18
  };

  constexpr int16_t phi0p05 = 522;  // round(521.52189...) = phi2short(0.05);
  constexpr int16_t phi0p06 = 626;  // round(625.82270...) = phi2short(0.06);
  constexpr int16_t phi0p07 = 730;  // round(730.12648...) = phi2short(0.07);

  constexpr int16_t phicuts[nPairs]{phi0p05,
                                             phi0p07,
                                             phi0p07,
                                             phi0p05,
                                             phi0p06,
                                             phi0p06,
                                             phi0p05,
                                             phi0p05,
                                             phi0p06,
                                             phi0p06,
                                             phi0p06,
                                             phi0p05,
                                             phi0p05,
                                             phi0p05,
                                             phi0p05,
                                             phi0p05,
                                             phi0p05,
                                             phi0p05,
                                             phi0p05};
  //   phi0p07, phi0p07, phi0p06,phi0p06, phi0p06,phi0p06};  // relaxed cuts

  constexpr float minz[nPairs] = {
      -20., 0., -30., -22., 10., -30., -70., -70., -22., 15., -30, -70., -70., -20., -22., 0, -30., -70., -70.};
  constexpr float maxz[nPairs] = {
      20., 30., 0., 22., 30., -10., 70., 70., 22., 30., -15., 70., 70., 20., 22., 30., 0., 70., 70.};
  constexpr float maxr[nPairs] = {
      20., 9., 9., 20., 7., 7., 5., 5., 20., 6., 6., 5., 5., 20., 20., 9., 9., 9., 9.};

  // end constants
  // clang-format on

  using CellNeighbors = CAConstants::CellNeighbors;
  using CellTracks = CAConstants::CellTracks;
  using CellNeighborsVector = CAConstants::CellNeighborsVector;
  using CellTracksVector = CAConstants::CellTracksVector;

  __global__ void initDoublets(GPUCACell::OuterHitOfCell* isOuterHitOfCell,
                               int nHits,
                               CellNeighborsVector* cellNeighbors,
                               CellNeighbors* cellNeighborsContainer,
                               CellTracksVector* cellTracks,
                               CellTracks* cellTracksContainer) {
    std::cout << "The max number of active doublets is " << CAConstants::maxNumOfActiveDoublets() << '\n';
    std::cout << "isOuterHitOfCell = " << isOuterHitOfCell << '\n';

    assert(isOuterHitOfCell);
    int first = blockIdx.x * blockDim.x + threadIdx.x;
    for (int i = first; i < nHits; i += gridDim.x * blockDim.x)
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

  __global__
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
                                uint32_t maxNumOfDoublets) {
    auto const& __restrict__ hh = *hhp;
    doubletsFromHisto(layerPairs,
                      nActualPairs,
                      cells,
                      nCells,
                      cellNeighbors,
                      cellTracks,
                      hh,
                      isOuterHitOfCell,
                      phicuts,
                      minz,
                      maxz,
                      maxr,
                      ideal_cond,
                      doClusterCut,
                      doZ0Cut,
                      doPtCut,
                      maxNumOfDoublets);
  }

}  // namespace gpuPixelDoublets

#endif  // RecoLocalTracker_SiPixelRecHits_plugins_gpuPixelDouplets_h
