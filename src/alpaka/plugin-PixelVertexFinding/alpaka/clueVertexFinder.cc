
#include "AlpakaCore/config.h"
#include "AlpakaCore/memory.h"
#include "AlpakaCore/workdivision.h"

#include "gpuVertexFinder.h"
#include "gpuClusterTracksByDensity.h"
#include "gpuClusterTracksDBSCAN.h"
#include "gpuClusterTracksIterative.h"
#include "gpuFitVertices.h"
#include "gpuSortByPt2.h"
#include "gpuSplitVertices.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  namespace clueVertexFinder {

    struct loadTracks {
      template <typename TAcc>
      ALPAKA_FN_ACC void operator()(
          const TAcc& acc, gpuVertexFinder::TkSoA const* ptracks, ZVertexSoA* soa, WorkSpace* pws, float ptMin) const {
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

    ZVertexAlpaka Producer::makeAsync(::ALPAKA_ACCELERATOR_NAMESPACE::gpuVertexFinder::TkSoA const* tksoa,
                                      float ptMin,
                                      Queue& queue) const {
      ALPAKA_ASSERT_ACC(tksoa);
      const auto maxTracks = tksoa::stride();
      std::cout << "max tracks = " << maxTracks << std::endl;
      auto vertices = cms::alpakatools::make_device_buffer<ZVertexSoA>(queue);
      auto verticesView = vertices.view();
      auto vertexTrackDataView = vertices.view<::reco::ZVertexTracksSoA>();

      // Initialize the workspace
      auto workspace = cms::alpakatools::make_device_buffer<WorkSpace>(queue);
      auto workspaceView = workspace.data();
      // const auto initWorkDiv = cms::alpakatools::make_workdiv<Acc1D>(1, 1);
      // alpaka::exec<Acc1D>(
      //     queue, initWorkDiv, ALPAKA_ACCELERATOR_NAMESPACE::vertexFinder::Init{}, verticesView, workspaceView);

      //Load Tracks
      const uint32_t blockSize = 128;
      const uint32_t gridSize = cms::alpakatools::divide_up_by(maxTracks, blockSize);
      const auto loadTracksWorkDiv = cms::alpakatools::make_workdiv<Acc1D>(gridSize, blockSize);
      alpaka::exec<Acc1D>(queue,
                          loadTracksWorkDiv,
                          LoadTracks{},
                          tracks_view,
                          verticesView,
                          vertexTrackDataView,
                          workspaceView,
                          ptMin,
                          ptMax);

      // Copy number of tracks to host
      auto nTracksBuf = cms::alpakatools::make_host_buffer<uint32_t>(queue);
      alpaka::memcpy(queue, nTracksBuf, cms::alpakatools::make_device_view<uint32_t>(queue, workspaceView.ntrks()));
      alpaka::wait(queue);
      const auto nTracks = *nTracksBuf;

      // Run CLUEstering
      if (nTracks > 0) {
        clue::Clusterer<1> clusterer(queue, m_dc, m_rhoc, m_dm);
        clue::PointsDevice<1, Device> d_points(
            queue, nTracks, workspaceView.zt(), workspaceView.ptt2(), workspaceView.iv());
        clusterer.make_clusters(queue, d_points);
        auto h_points = clue::copyToHost(queue, h_points, d_points);
        alpaka::wait(queue);
        const auto nVertices = h_points.n_clusters();
        alpaka::memcpy(queue,
                       cms::alpakatools::make_device_view<uint32_t>(queue, verticesView.nvFinal()),
                       cms::alpakatools::make_host_view<uint32_t>(nVertices));
        alpaka::memcpy(queue,
                       cms::alpakatools::make_device_view<uint32_t>(queue, workspaceView.nvIntermediate()),
                       cms::alpakatools::make_host_view<uint32_t>(nVertices));
      }
      const auto finderSorterWorkDiv = cms::alpakatools::make_workdiv<Acc1D>(1, 1024 - 128);
      alpaka::exec<Acc1D>(queue,
                          finderSorterWorkDiv,
                          ALPAKA_ACCELERATOR_NAMESPACE::gpuVertexFinder::fitVerticesKernel{},
                          verticesView,
                          vertexTrackDataView,
                          workspaceView,
                          maxChi2ForFirstFit);
      const auto splitterFitterWorkDiv = cms::alpakatools::make_workdiv<Acc1D>(1024, 128);
      alpaka::exec<Acc1D>(queue,
                          splitterFitterWorkDiv,
                          ALPAKA_ACCELERATOR_NAMESPACE::gpuVertexFinder::splitVerticesKernel{},
                          verticesView,
                          vertexTrackDataView,
                          workspaceView,
                          maxChi2ForSplit);
      alpaka::exec<Acc1D>(queue,
                          finderSorterWorkDiv,
                          ALPAKA_ACCELERATOR_NAMESPACE::gpuVertexFinder::fitVerticesKernel{},
                          verticesView,
                          vertexTrackDataView,
                          workspaceView,
                          maxChi2ForFinalFit);
      alpaka::exec<Acc1D>(queue,
                          finderSorterWorkDiv,
                          ALPAKA_ACCELERATOR_NAMESPACE::gpuVertexFinder::sortByPt2Kernel{},
                          verticesView,
                          vertexTrackDataView,
                          workspaceView);

      return vertices;
    }

  }  // namespace clueVertexFinder
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
