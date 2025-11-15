#ifndef plugin_PixelVertexFinding_alpaka_gpuVertexFinder_h
#define plugin_PixelVertexFinding_alpaka_gpuVertexFinder_h

#include "AlpakaCore/config.h"
#include "AlpakaDataFormats/alpaka/PixelTrackAlpaka.h"
#include "AlpakaDataFormats/alpaka/ZVertexAlpaka.h"
#include "VertexWorkspace.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  namespace gpuVertexFinder {

    using ZVertices = ZVertexSoA;
    using TkSoA = pixelTrack::TrackSoA;

    class Producer {
    public:
      using ZVertices = ZVertexSoA;
      using WorkSpace = gpuVertexFinder::WorkSpace;
      using TkSoA = pixelTrack::TrackSoA;

      Producer(bool oneKernel,
               bool useDensity,
               bool useDBSCAN,
               bool useIterative,
               int iminT,      // min number of neighbours to be "core"
               float ieps,     // max absolute distance to cluster
               float ierrmax,  // max error to be "seed"
               float ichi2max  // max normalized distance to cluster
               )
          : oneKernel_(oneKernel && !(useDBSCAN || useIterative)),
            useDensity_(useDensity),
            useDBSCAN_(useDBSCAN),
            useIterative_(useIterative),
            minT(iminT),
            eps(ieps),
            errmax(ierrmax),
            chi2max(ichi2max) {}

      ~Producer() = default;

      ZVertexAlpaka makeAsync(TkSoA const* tksoa, float ptMin, Queue& queue) const;

    private:
      const bool oneKernel_;
      const bool useDensity_;
      const bool useDBSCAN_;
      const bool useIterative_;

      int minT;       // min number of neighbours to be "core"
      float eps;      // max absolute distance to cluster
      float errmax;   // max error to be "seed"
      float chi2max;  // max normalized distance to cluster
    };

  }  // namespace gpuVertexFinder
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // plugin_PixelVertexFinding_alpaka_gpuVertexFinder_h
