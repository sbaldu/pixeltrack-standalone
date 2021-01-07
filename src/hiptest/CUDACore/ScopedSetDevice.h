#ifndef HeterogeneousCore_CUDAUtilities_ScopedSetDevice_h
#define HeterogeneousCore_CUDAUtilities_ScopedSetDevice_h

#include "CUDACore/cudaCheck.h"

#include <hip/hip_runtime.h>

namespace cms {
  namespace cuda {
    class ScopedSetDevice {
    public:
      explicit ScopedSetDevice(int newDevice) {
        cudaCheck(hipGetDevice(&prevDevice_));
        cudaCheck(hipSetDevice(newDevice));
      }

      ~ScopedSetDevice() {
        // Intentionally don't check the return value to avoid
        // exceptions to be thrown. If this call fails, the process is
        // doomed anyway.
        hipSetDevice(prevDevice_);
      }

    private:
      int prevDevice_;
    };
  }  // namespace cuda
}  // namespace cms

#endif
