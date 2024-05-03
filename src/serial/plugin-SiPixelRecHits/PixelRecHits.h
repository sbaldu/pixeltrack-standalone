#ifndef RecoLocalTracker_SiPixelRecHits_plugins_PixelRecHits_h
#define RecoLocalTracker_SiPixelRecHits_plugins_PixelRecHits_h

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>

#include "DataFormats/BeamSpotPOD.h"
#include "DataFormats/approx_atan2.h"
#include "CUDADataFormats/SiPixelClustersSoA.h"
#include "CUDADataFormats/SiPixelDigisSoA.h"
#include "CUDADataFormats/TrackingRecHit2DHeterogeneous.h"

#include "DataFormats/HitsCoordsSoA.h"

namespace pixelgpudetails {

  class PixelRecHitGPUKernel {
  public:
    PixelRecHitGPUKernel() = default;
    ~PixelRecHitGPUKernel() = default;

    PixelRecHitGPUKernel(const PixelRecHitGPUKernel&) = delete;
    PixelRecHitGPUKernel(PixelRecHitGPUKernel&&) = delete;
    PixelRecHitGPUKernel& operator=(const PixelRecHitGPUKernel&) = delete;
    PixelRecHitGPUKernel& operator=(PixelRecHitGPUKernel&&) = delete;

    TrackingRecHit2DCPU makeHits(SiPixelDigisSoA const& digis_d,
                                 SiPixelClustersSoA const& clusters_d,
                                 BeamSpotPOD const& bs_d,
                                 pixelCPEforGPU::ParamsOnGPU const* cpeParams) const;
    TrackingRecHit2DCPU makeHits() {
      HitsCoordsSoA hits;
      uint32_t nHits{};

      std::ifstream iFile("data/track-ml/hits_1000.csv");
      if (!iFile.is_open()) {
        std::cerr << "Error opening file" << std::endl;
      }

      std::string line;
      // TODO: maybe add a bit of error handling in the header
      getline(iFile, line);
      while (getline(iFile, line)) {
        std::stringstream fileStream(line);
        std::string temp;

        getline(fileStream, temp, ',');
        float x{std::stof(temp)};
        hits.x.push_back(x);
        getline(fileStream, temp, ',');
        float y{std::stof(temp)};
        hits.y.push_back(y);
        getline(fileStream, temp, ',');
        hits.z.push_back(std::stof(temp));

        hits.r.push_back(std::sqrt(x * x + y * y));

        int16_t phi{phi2short(std::atan(y / x))};
        hits.phi.push_back(phi);

        getline(fileStream, temp, ',');
        hits.global_indexes.push_back(std::stoi(temp));

        ++nHits;
      }

      std::ifstream iFileTruth("data/track-ml/hits_1000.csv");
      if (!iFileTruth.is_open()) {
        std::cerr << "Error opening file" << std::endl;
      }

      // TODO: maybe add a bit of error handling in the header
      getline(iFileTruth, line);
      while (getline(iFileTruth, line)) {
        std::stringstream fileStream(line);
        std::string temp;

        getline(fileStream, temp, ',');
        hits.particle_indexes.push_back(std::stoi(temp));
      }

      

      std::vector<uint32_t> layerStart_ = {0};
      for (size_t j{1}; j < hits.global_indexes.size() - 1; ++j) {
        if (hits.global_indexes[j + 1] != hits.global_indexes[j]) {
          layerStart_.push_back(j + 1);
        }
      }

      hits.reset();  // reset the view
      TrackingRecHit2DCPU hits_d(nHits, std::move(hits), std::move(layerStart_), nullptr);
      return hits_d;
    }
  };
}  // namespace pixelgpudetails

#endif  // RecoLocalTracker_SiPixelRecHits_plugins_PixelRecHits_h
