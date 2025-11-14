
#include "AlpakaCore/config.h"
#include "AlpakaCore/memory.h"
#include "AlpakaCore/workdivision.h"

#include "gpuVertexFinder.h"
#include "clueVertexFinder.h"
#include "gpuClusterTracksByDensity.h"
#include "gpuClusterTracksDBSCAN.h"
#include "gpuClusterTracksIterative.h"
#include "gpuFitVertices.h"
#include "gpuSortByPt2.h"
#include "gpuSplitVertices.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  namespace gpuVertexFinder {

    ZVertexAlpaka CLUEVertexProducer::makeAsync(TkSoA const* tksoa, float ptMin, Queue& queue) const {
      ALPAKA_ASSERT_ACC(tksoa);
      const auto maxTracks = TkSoA::stride();
      std::cout << "max tracks = " << maxTracks << std::endl;
      auto vertices = cms::alpakatools::make_device_buffer<ZVertexSoA>(queue);
      auto verticesView = vertices.view();
      auto vertexTrackDataView = vertices.view<::reco::ZVertexTracksSoA>();

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
                          fitVerticesKernel{},
                          verticesView,
                          vertexTrackDataView,
                          workspaceView,
                          maxChi2ForFirstFit);
      const auto splitterFitterWorkDiv = cms::alpakatools::make_workdiv<Acc1D>(1024, 128);
      alpaka::exec<Acc1D>(queue,
                          splitterFitterWorkDiv,
                          splitVerticesKernel{},
                          verticesView,
                          vertexTrackDataView,
                          workspaceView,
                          maxChi2ForSplit);
      alpaka::exec<Acc1D>(queue,
                          finderSorterWorkDiv,
                          fitVerticesKernel{},
                          verticesView,
                          vertexTrackDataView,
                          workspaceView,
                          maxChi2ForFinalFit);
      alpaka::exec<Acc1D>(
          queue, finderSorterWorkDiv, sortByPt2Kernel{}, verticesView, vertexTrackDataView, workspaceView);

      return vertices;
    }

  }  // namespace gpuVertexFinder
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
