
#ifndef alpaka_vertex_workspace_h
#define alpaka_vertex_workspace_h

#include "AlpakaDataFormats/ZVertexSoA.h"
#include <cstdint>

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  namespace gpuVertexFinder {

    struct WorkSpace {
      static constexpr uint32_t MAXTRACKS = ZVertexSoA::MAXTRACKS;
      static constexpr uint32_t MAXVTX = ZVertexSoA::MAXVTX;

      uint32_t ntrks;            // number of "selected tracks"
      uint16_t itrk[MAXTRACKS];  // index of original track
      float zt[MAXTRACKS];       // input track z at bs
      float ezt2[MAXTRACKS];     // input error^2 on the above
      float ptt2[MAXTRACKS];     // input pt^2 on the above
      uint8_t izt[MAXTRACKS];    // interized z-position of input tracks
      int32_t iv[MAXTRACKS];     // vertex index for each associated track

      uint32_t nvIntermediate;  // the number of vertices after splitting pruning etc.
    };

  }  // namespace gpuVertexFinder
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // alpaka_vertex_workspace_h
