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
#include <iomanip>

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
    TrackingRecHit2DCPU makeHits(int eventID){
      HitsCoordsSoA hits;
      uint32_t nHits{};
      std::string dir = "/eos/user/s/srossiti/track-ml/preprocessed_train_100_events/";
      //filename is trackml_049.csv with 49 being the eventID
      std::stringstream eventid_ss;
      eventid_ss << std::setw(3) << std::setfill('0') << (eventID-1);
      std::string filename = "trackml_" + eventid_ss.str() + ".csv";
      std::ifstream iFile(dir + filename);
      // std::ifstream iFile("data/track-ml/hits_high_eta.csv");
      // std::ifstream iFile("data/track-ml/hits_1000.csv");
      std::cout << "Event ID: " << eventID << " Filename: " << filename << std::endl;
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
        hits.particle_indexes.push_back(std::stol(temp));
        getline(fileStream, temp, ',');
        hits.particle_Vzs.push_back(std::stof(temp));
        getline(fileStream, temp, ',');
        getline(fileStream, temp, ',');
        hits.global_indexes.push_back(std::stoi(temp));
        getline(fileStream, temp, ',');
        hits.particle_pTs.push_back(std::stof(temp));
        getline(fileStream, temp, ',');
        hits.particle_dRs.push_back(std::stof(temp));
        getline(fileStream, temp, ',');
        hits.particle_nHits.push_back(std::stoi(temp));
        ++nHits;
      }

      

      std::vector<std::size_t> indexes(hits.global_indexes.size());
      std::iota(indexes.begin(), indexes.end(), 0);
      std::sort(indexes.begin(), indexes.end(), [&hits](std::size_t i1, std::size_t i2) {
        return hits.global_indexes[i1] < hits.global_indexes[i2];
      });
      //sort hits.x, hits.y, hits.z, hits.r, hits.phi, hits.global_indexes, hits.particle_indexes, hits.particle_pTs by hits.global_indexes
      HitsCoordsSoA hits_sorted;
      for (size_t i{0}; i < hits.global_indexes.size(); ++i) {
        hits_sorted.x.push_back(hits.x[indexes[i]]);
        hits_sorted.y.push_back(hits.y[indexes[i]]);
        hits_sorted.z.push_back(hits.z[indexes[i]]);
        hits_sorted.r.push_back(hits.r[indexes[i]]);
        hits_sorted.phi.push_back(hits.phi[indexes[i]]);
        hits_sorted.global_indexes.push_back(hits.global_indexes[indexes[i]]);
        hits_sorted.particle_indexes.push_back(hits.particle_indexes[indexes[i]]);
        hits_sorted.particle_pTs.push_back(hits.particle_pTs[indexes[i]]);
        hits_sorted.particle_dRs.push_back(hits.particle_dRs[indexes[i]]);
        hits_sorted.particle_Vzs.push_back(hits.particle_Vzs[indexes[i]]);
        hits_sorted.particle_nHits.push_back(hits.particle_nHits[indexes[i]]);
      }
      auto currentLayerIndex = 0;
      auto lastLayerIndex = 0;
      std::vector<uint32_t> layerStart_ = {0};
      for (size_t j{0}; j < hits_sorted.global_indexes.size() - 1; ++j) {
        if (hits_sorted.global_indexes[j + 1] != hits_sorted.global_indexes[j]) {
          currentLayerIndex=hits_sorted.global_indexes[j]+1;
          layerStart_.push_back(j + 1);
          while (currentLayerIndex < hits_sorted.global_indexes[j + 1]) {
            currentLayerIndex++;
            layerStart_.push_back(j + 1);
            lastLayerIndex = j + 1;
          }
        }
      }
      
      // std::cout << "nHits"<< nHits << std::endl;
      

      // for (auto & l_off: layerStart_){
      //   std::cout << "Layer offset: " << l_off << std::endl;
      // }
      // for(auto & r: hits_sorted.r)
      // {
      //   std::cout << "Hit r: " << r << std::endl;
      // }
      // hits.reset();  // reset the view
      hits_sorted.reset();  // reset the view
      TrackingRecHit2DCPU hits_d(nHits, std::move(hits_sorted), std::move(layerStart_), nullptr);
      return hits_d;
    }
  };
}  // namespace pixelgpudetails

#endif  // RecoLocalTracker_SiPixelRecHits_plugins_PixelRecHits_h
