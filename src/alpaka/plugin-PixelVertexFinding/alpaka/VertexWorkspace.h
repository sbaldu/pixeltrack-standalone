
#ifndef alpaka_vertex_workspace_h
#define alpaka_vertex_workspace_h

#include "AlpakaDataFormats/ZVertexSoA.h"
#include "AlpakaCore/memory.h"
#include <cstddef>
#include <cstdint>

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  namespace gpuVertexFinder {

    // struct WorkSpace {
    //   static constexpr uint32_t MAXTRACKS = ZVertexSoA::MAXTRACKS;
    //   static constexpr uint32_t MAXVTX = ZVertexSoA::MAXVTX;

    //   uint32_t ntrks;            // number of "selected tracks"
    //   uint16_t itrk[MAXTRACKS];  // index of original track
    //   float zt[MAXTRACKS];       // input track z at bs
    //   float ezt2[MAXTRACKS];     // input error^2 on the above
    //   float ptt2[MAXTRACKS];     // input pt^2 on the above
    //   uint8_t izt[MAXTRACKS];    // interized z-position of input tracks
    //   int32_t iv[MAXTRACKS];     // vertex index for each associated track

    //   uint32_t nvIntermediate;  // the number of vertices after splitting pruning etc.
    // };

  }  // namespace gpuVertexFinder

  struct WorkSpaceView {
    static constexpr auto size_bytes() {
      return sizeof(int32_t) + sizeof(uint16_t) + sizeof(uint8_t) + 3 * sizeof(float);
    }
    static constexpr uint32_t MAXTRACKS = ZVertexSoA::MAXTRACKS;
    static constexpr uint32_t MAXVTX = ZVertexSoA::MAXVTX;

    uint16_t* itrk;  // index of original track
    float* zt;       // input track z at bs
    float* ezt2;     // input error^2 on the above
    float* ptt2;     // input pt^2 on the above
    uint8_t* izt;    // interized z-position of input tracks
    int32_t* iv;     // vertex index for each associated track

    uint32_t* ntrks;              // number of "selected tracks"
    uint32_t nvIntermediate = 0;  // the number of vertices after splitting pruning etc.

    WorkSpaceView(std::byte* buffer, uint32_t n_tracks) {
      itrk = reinterpret_cast<uint16_t*>(buffer);
      zt = reinterpret_cast<float*>(buffer + n_tracks * (sizeof(uint16_t)));
      ezt2 = reinterpret_cast<float*>(buffer + n_tracks * (sizeof(uint16_t) + sizeof(float)));
      ptt2 = reinterpret_cast<float*>(buffer + n_tracks * (sizeof(uint16_t) + 2 * sizeof(float)));
      izt = reinterpret_cast<uint8_t*>(buffer + n_tracks * (sizeof(uint16_t) + 3 * sizeof(float)));
      iv = reinterpret_cast<int32_t*>(buffer + n_tracks * (sizeof(uint16_t) + 3 * sizeof(float) + sizeof(uint8_t)));
      ntrks = reinterpret_cast<uint32_t*>(
          buffer + n_tracks * (sizeof(uint16_t) + 3 * sizeof(float) + sizeof(uint8_t) + sizeof(int32_t)));
    }
  };

  template <typename TDev>
  class WorkSpace {
  private:
    cms::alpakatools::device_buffer<TDev, std::byte[]> m_buffer;
    WorkSpaceView m_view;

  public:
    template <typename TQueue>
    WorkSpace(TQueue& queue, uint32_t n_tracks)
        : m_buffer{cms::alpakatools::make_device_buffer<std::byte[]>(
              queue, n_tracks * WorkSpaceView::size_bytes() + sizeof(uint32_t))},
          m_view{m_buffer.data(), n_tracks} {}

    const auto& view() const { return m_view; }
    auto& view() { return m_view; }
  };
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // alpaka_vertex_workspace_h
