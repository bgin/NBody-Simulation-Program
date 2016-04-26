
#include "NBody_Helpers.h"


auto nbody_implementation::NBodyHelpers::check_particles(_In_ const std::size_t nParticles)->void {
	constexpr std::size_t lo_bound{ 128U };
	if (nParticles <= lo_bound)
		throw std::runtime_error(std::string("***Fatal Error***: nParticle out_of_boundary\n"));
}

auto nbody_implementation::NBodyHelpers::check_threads(_In_ const int nThreads)->void {
	constexpr int lo_bound{ 2U };
	if (nThreads <= lo_bound)
		throw std::runtime_error(std::string("***Fatal Error***: nThreads out-of-boundary\n"));
}

auto  nbody_implementation::NBodyHelpers::check_step(_In_ const double nSteps)->void {
	constexpr double zero{ 0.0 };
	if (nSteps == zero)
		throw std::runtime_error(std::string("***Fatal Error***: nSteps out-of-boundary"));
}