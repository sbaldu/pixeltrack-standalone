//
// Original Author: Felice Pantaleo, CERN
//

#include <array>
#include <cassert>
#include <fstream>
#include <iostream>
#include <functional>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>

#include "Framework/Event.h"

#include "CAHitNtupletGeneratorOnGPU.h"

namespace {

  template <typename T>
  T sqr(T x) {
    return x * x;
  }

  inline int16_t stos(std::string str) {
    int value = std::stoi(str);

    if (value < std::numeric_limits<short>::min() || value > std::numeric_limits<short>::max()) {
      throw std::range_error("String value out of range for short");
    }

    return static_cast<short>(value);
  }

  int16_t* makePhiCuts() {
    constexpr int nPairs = 13 + 2 + 4;

    constexpr int16_t phi0p05 = 522;  // round(521.52189...) = phi2short(0.05);
    constexpr int16_t phi0p06 = 626;  // round(625.82270...) = phi2short(0.06);
    constexpr int16_t phi0p07 = 730;  // round(730.12648...) = phi2short(0.07);

    return new int16_t[nPairs]{phi0p05,
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
  }

  cAHitNtupletGenerator::QualityCuts makeQualityCuts() {
    auto coeff =
        std::vector<double>{0.68177776, 0.74609577, -0.08035491, 0.00315399};  // chi2Coeff
    return cAHitNtupletGenerator::QualityCuts{
        // polynomial coefficients for the pT-dependent chi2 cut
        {(float)coeff[0], (float)coeff[1], (float)coeff[2], (float)coeff[3]},
        // max pT used to determine the chi2 cut
        10.f,  // chi2MaxPt
               // chi2 scale factor: 30 for broken line fit, 45 for Riemann fit
        30.f,  // chi2Scale
               // regional cuts for triplets
        {
            0.3f,  //tripletMaxTip
            0.5f,  // tripletMinPt
            12.f   // tripletMaxZip
        },
        // regional cuts for quadruplets
        {
            0.5f,  // quadrupletMaxTip
            0.3f,  // quadrupletMinPt
            12.f   // quadrupletMaxZip
        }};
  }
}  // namespace

using namespace std;
CAHitNtupletGeneratorOnGPU::CAHitNtupletGeneratorOnGPU(edm::ProductRegistry& reg)
    : m_params(false,              // onGPU
               3,                  // minHitsPerNtuplet,
               458752,             // maxNumberOfDoublets
               false,              //useRiemannFit
               true,               // fit5as4,
               true,               //includeJumpingForwardDoublets
               true,               // earlyFishbone
               false,              // lateFishbone
               true,               // idealConditions
               false,              //fillStatistics
               true,               // doClusterCut
               true,               // doZ0Cut
               true,               // doPtCut
               0.899999976158,     // ptmin
               0.00200000009499,   // CAThetaCutBarrel
               0.00300000002608,   // CAThetaCutForward
               0.0328407224959,    // hardCurvCut
               0.15000000596,      // dcaCutInnerTriplet
               0.25,               // dcaCutOuterTriplet
               makeQualityCuts(),  // QualityCuts
               makePhiCuts())      // phiCuts
{
#ifdef DUMP_GPU_TK_TUPLES
  printf("TK: %s %s % %s %s %s %s %s %s %s %s %s %s %s %s %s\n",
         "tid",
         "qual",
         "nh",
         "charge",
         "pt",
         "eta",
         "phi",
         "tip",
         "zip",
         "chi2",
         "h1",
         "h2",
         "h3",
         "h4",
         "h5");
#endif

  m_counters = new Counters();
  memset(m_counters, 0, sizeof(Counters));
}

CAHitNtupletGeneratorOnGPU::CAHitNtupletGeneratorOnGPU(edm::ProductRegistry& reg, std::string filename)
    : m_params(false,              // onGPU
               3,                  // minHitsPerNtuplet,
               458752,             // maxNumberOfDoublets
               false,              //useRiemannFit
               true,               // fit5as4,
               true,               //includeJumpingForwardDoublets
               true,               // earlyFishbone
               false,              // lateFishbone
               true,               // idealConditions
               false,              //fillStatistics
               true,               // doClusterCut
               true,               // doZ0Cut
               true,               // doPtCut
               0.899999976158,     // ptmin
               0.00200000009499,   // CAThetaCutBarrel
               0.00300000002608,   // CAThetaCutForward
               0.0328407224959,    // hardCurvCut
               0.15000000596,      // dcaCutInnerTriplet
               0.25,               // dcaCutOuterTriplet
               makeQualityCuts(),  // QualityCuts
               makePhiCuts())      // phiCuts
{
  std::ifstream file(filename);
  if (!file.is_open()) {
    printf("No file opened\n");
    return;
  }

  std::vector<std::string> params_from_file;
  std::string line;
  while (std::getline(file, line)) {
    std::string value;
    std::vector<std::string> row;
    std::stringstream sstr(line);
    while (std::getline(sstr, value, ','))
      row.push_back(value);
    params_from_file.insert(params_from_file.end(), row.begin(), row.end());
  }
  std::cout << "File opened\n";
  float CAThetaCutBarrel_ = std::stof(params_from_file[0]);
  float CAThetaCutForward_ = std::stof(params_from_file[1]);
  float dcaCutInnerTriplet_ = std::stof(params_from_file[2]);
  float dcaCutOuterTriplet_ = std::stof(params_from_file[3]);
  float hardCurvCut_ = std::stof(params_from_file[4]);
  bool doZ0Cut_ = std::stof(params_from_file[5]);
  int16_t phiCuts[47] = {
      stos(params_from_file[6]), stos(params_from_file[7]), stos(params_from_file[8]), stos(params_from_file[9]),
      stos(params_from_file[0]), stos(params_from_file[0]), stos(params_from_file[0]), stos(params_from_file[0]),
      stos(params_from_file[0]), stos(params_from_file[0]), stos(params_from_file[0]), stos(params_from_file[0]),
      stos(params_from_file[0]), stos(params_from_file[0]), stos(params_from_file[0]), stos(params_from_file[0]),
      stos(params_from_file[0]), stos(params_from_file[0]), stos(params_from_file[0]), stos(params_from_file[0]),
      stos(params_from_file[0]), stos(params_from_file[0]), stos(params_from_file[0]), stos(params_from_file[0]),
      stos(params_from_file[0]), stos(params_from_file[0]), stos(params_from_file[0]), stos(params_from_file[0]),
      stos(params_from_file[0]), stos(params_from_file[0]), stos(params_from_file[0]), stos(params_from_file[0]),
      stos(params_from_file[0]), stos(params_from_file[0]), stos(params_from_file[0]), stos(params_from_file[0]),
      stos(params_from_file[0]), stos(params_from_file[0]), stos(params_from_file[0]), stos(params_from_file[0]),
      stos(params_from_file[0]), stos(params_from_file[0]), stos(params_from_file[0]), stos(params_from_file[0]),
      stos(params_from_file[0]), stos(params_from_file[0]), stos(params_from_file[0])};

  m_params.doZ0Cut_=doZ0Cut_;
  m_params.CAThetaCutBarrel_=CAThetaCutBarrel_;
  m_params.CAThetaCutForward_=CAThetaCutForward_;
  m_params.hardCurvCut_=hardCurvCut_;
  m_params.dcaCutInnerTriplet_=dcaCutInnerTriplet_;
  m_params.dcaCutOuterTriplet_=dcaCutOuterTriplet_;
  m_params.phiCuts_=phiCuts;


#ifdef DUMP_GPU_TK_TUPLES
  printf("TK: %s %s % %s %s %s %s %s %s %s %s %s %s %s %s %s\n",
         "tid",
         "qual",
         "nh",
         "charge",
         "pt",
         "eta",
         "phi",
         "tip",
         "zip",
         "chi2",
         "h1",
         "h2",
         "h3",
         "h4",
         "h5");
#endif

  m_counters = new Counters();
  memset(m_counters, 0, sizeof(Counters));
}

CAHitNtupletGeneratorOnGPU::~CAHitNtupletGeneratorOnGPU() {
  if (m_params.doStats_) {
    CAHitNtupletGeneratorKernelsCPU::printCounters(m_counters);
  }
  delete m_counters;
}

PixelTrackHeterogeneous CAHitNtupletGeneratorOnGPU::makeTuples(
    TrackingRecHit2DCPU const& hits_d, float bfield) const {
  PixelTrackHeterogeneous tracks(std::make_unique<pixelTrack::TrackSoA>());

  auto* soa = tracks.get();
  assert(soa);

  CAHitNtupletGeneratorKernelsCPU kernels(m_params);
  kernels.counters_ = m_counters;
  kernels.allocateOnGPU(nullptr);

  kernels.buildDoublets(hits_d, nullptr);
  kernels.launchKernels(hits_d, soa, nullptr);
  kernels.fillHitDetIndices(
      hits_d.view(), soa, nullptr);  // in principle needed only if Hits not "available"

  return tracks;
}
