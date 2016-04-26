
#include "NBodyVecAOSImpl.h"

nbody_implementation::NBodyVecAos::NBodyVecAos(_In_ const std::size_t nParticles, _In_ const int nThreads, _In_ const double nStep) :
m_uinParticles( (nbody_implementation::NBodyHelpers::check_particles(nParticles),nParticles) ),
m_inThreads( (nbody_implementation::NBodyHelpers::check_threads(nThreads),nThreads) ),
m_dnStep( (nbody_implementation::NBodyHelpers::check_step(nStep),nStep) )
{
	// Default initialization with zeroes.
	this->m_vNBodies = std::vector<NBodyObject>(this->m_uinParticles);
	
	
	
}

nbody_implementation::NBodyVecAos::NBodyVecAos(_In_ const NBodyVecAos &rhs) :
m_uinParticles{ rhs.m_uinParticles },
m_inThreads{ rhs.m_inThreads },
m_dnStep{ rhs.m_dnStep },
m_vNBodies{ rhs.m_vNBodies }
{

}

nbody_implementation::NBodyVecAos::NBodyVecAos(_In_ NBodyVecAos &&rhs) :
m_uinParticles{ std::move(rhs.m_uinParticles) },
m_inThreads{ std::move(rhs.m_inThreads) },
m_dnStep{ std::move(rhs.m_dnStep) },
m_vNBodies{ std::move(rhs.m_vNBodies) }
{

}

void    nbody_implementation::NBodyVecAos::run_simulation(_In_ double lo, _In_ double hi) {

	if ((std::fabs(lo - hi)) <= std::numeric_limits<double>::epsilon())
		throw std::runtime_error(std::string("Invalid input in: NBodyDArray::run_simualtion\n"));

#if defined VERBOSE
	std::cout << "Allocated: " << static_cast<double>((6 * this->m_uinParticles * sizeof(double))) / 1024.0 << " KiB of heap memory" << std::endl;

#endif
	clock_t seed{ clock() };
	auto rand = std::bind(std::uniform_real_distribution < double > {lo, hi}, std::default_random_engine(seed));
	const std::string RNG_TYPE{ typeid(rand).name() };
#if defined VERBOSE
	std::cout << "Created Random number generator of type: " << RNG_TYPE.c_str() << "with the interval of randomness" << "[" << -1.0 << " -- " << 1.0 << "]" << std::endl;
#endif
#if defined VERBOSE
	std::cout << "Initialization of particles vector position started... ";
#endif
	for (std::size_t i{ 0 }; i != this->m_uinParticles; ++i){
		this->m_vNBodies.operator[](i).m_dX = rand();
		this->m_vNBodies.operator[](i).m_dY = rand();
		this->m_vNBodies.operator[](i).m_dZ = rand();
	}
#if defined VERBOSE
	std::cout << "Done!!" << std::endl;
	std::cout << "Utilizing single CPU core " << std::endl;
#endif
#if (EXECUTE_NTIMES) == 0x1
	for (int k{ 0 }; k != NTIMES; ++k) {
#endif
		for (std::size_t i{ 0 }; i != this->m_uinParticles; ++i) {

			double Fx{ 0.0 }, Fy{ 0.0 }, Fz{ 0.0 };

			for (std::size_t j{ 0 }; j != this->m_uinParticles; ++j) {

				if (j != i) {

					// Prefetching no really possible here because
					// of AOS access:
					// [X0,Y0,Z0,VX0,VY0,VZ0.....Xn,Yn,Zn,VXn,VYn,VZn]
					const double dx{ this->m_vNBodies[j].m_dX - this->m_vNBodies[i].m_dX };
					const double dy{ this->m_vNBodies[j].m_dY - this->m_vNBodies[i].m_dY };
					const double dz{ this->m_vNBodies[j].m_dZ - this->m_vNBodies[i].m_dZ };

					const double distSquared{ dx * dx + dy * dy + dz * dz };
					const double invdistCube{ 1.0 / (distSquared * ::sqrt(distSquared)) };

					Fx += dx * invdistCube;
					Fy += dy * invdistCube;
					Fz += dz * invdistCube;
				}
			}
			this->m_vNBodies[i].m_dVX = Fx * this->m_dnStep;
			this->m_vNBodies[i].m_dVY = Fy * this->m_dnStep;
			this->m_vNBodies[i].m_dVZ = Fz * this->m_dnStep;
		}
		for (std::size_t i{ 0 }; i != this->m_uinParticles; ++i) {
			this->m_vNBodies[i].m_dX = this->m_vNBodies[i].m_dVX + this->m_dnStep;
			this->m_vNBodies[i].m_dY = this->m_vNBodies[i].m_dVY + this->m_dnStep;
			this->m_vNBodies[i].m_dZ = this->m_vNBodies[i].m_dVZ + this->m_dnStep;
		}
#if (EXECUTE_NTIMES) == 0x1
	}
#endif

}

void    nbody_implementation::NBodyVecAos::run_simulation_omp(_In_ const double lo, _In_ const double hi) {

	if ((std::fabs(lo - hi)) <= std::numeric_limits<double>::epsilon())
		throw std::runtime_error(std::string("Invalid input in: NBodyDArray::run_simualtion\n"));

#if defined VERBOSE
	std::cout << "Allocated: " << static_cast<double>((6 * this->m_uinParticles * sizeof(double))) / 1024.0 << " KiB of heap memory" << std::endl;

#endif
	clock_t seed{ clock() };
	auto rand = std::bind(std::uniform_real_distribution < double > {lo, hi}, std::default_random_engine(seed));
	const std::string RNG_TYPE{ typeid(rand).name() };
#if defined VERBOSE
	std::cout << "Created Random number generator of type: " << RNG_TYPE.c_str() << "with the interval of randomness" << "[" << -1.0 << " -- " << 1.0 << "]" << std::endl;
#endif
#if defined VERBOSE
	std::cout << "Initialization of particles vector position started... ";
#endif

	/* Copy object state to automatic variables for OpenMP MT */
	std::vector<nbody_implementation::NBodyObject> a_vNbodies(this->m_vNBodies);
	const std::size_t nParticles{ this->m_uinParticles };
	const int nThreads{ this->m_inThreads };
	const double nStep{ this->m_dnStep };
	auto const chunk_size = nParticles / nThreads;
	/* Initializing position vectors using unaligned store. */
	for (std::size_t i{ 0 }; i != nParticles; ++i) {
		a_vNbodies.operator[](i).m_dX = rand();
		a_vNbodies.operator[](i).m_dY = rand();
		a_vNbodies.operator[](i).m_dZ = rand();
	}
#if defined VERBOSE
	std::cout << "Done!!" << std::endl;
	std::cout << "OpenMP in use with: " << this->m_inThreads << " Threads." << std::endl;
#endif
	::omp_set_num_threads(nThreads);
	int k;
#if (EXECUTE_NTIMES) == 0x1
	for ( k = 0; k != NTIMES; ++k) {
	const double dstart_time{ ::omp_get_wtime() };
#endif
#pragma omp parallel for schedule(static,chunk_size)
	for (std::size_t i = 0; i < nParticles; ++i) {

		double Fx{ 0.0 }, Fy{ 0.0 }, Fz{ 0.0 };

		for (std::size_t j{ 0 }; j < nParticles; ++j) {

			if (j != i) {

				const double dx{ a_vNbodies.operator[](j).m_dX - a_vNbodies.operator[](i).m_dX };
				const double dy{ a_vNbodies.operator[](j).m_dY - a_vNbodies.operator[](i).m_dY };
				const double dz{ a_vNbodies.operator[](j).m_dZ - a_vNbodies.operator[](i).m_dZ };

				const double distSquared{ dx * dx + dy * dy + dz * dz };
				const double invdistCube{ 1.0 / (distSquared * ::sqrt(distSquared)) };

				Fx += dx * invdistCube;
				Fy += dy * invdistCube;
				Fz += dz * invdistCube;
			}
		}
		a_vNbodies.operator[](i).m_dVX = Fx * nStep;
		a_vNbodies.operator[](i).m_dVY = Fy * nStep;
		a_vNbodies.operator[](i).m_dVZ = Fz * nStep;
	}
	for (std::size_t i{ 0 }; i != nParticles; ++i) {
		a_vNbodies.operator[](i).m_dX = a_vNbodies.operator[](i).m_dVX + nStep;
		a_vNbodies.operator[](i).m_dY = a_vNbodies.operator[](i).m_dVY + nStep;
		a_vNbodies.operator[](i).m_dZ = a_vNbodies.operator[](i).m_dVZ + nStep;
	}

#if (EXECUTE_NTIMES) == 0x1
	const double dstop_time{ ::omp_get_wtime() };
	std::cout << "Iteration #: " << k << " executed in: " << std::fixed << std::setprecision(9) << dstop_time - dstart_time << " sec, " <<
		(dstop_time - dstart_time) * 1.0E+9 << " nanesconds." << std::endl;
}
#endif
	/* Copy automatic variables to object state */
	this->m_vNBodies.operator=(a_vNbodies);
}

/* Not implemented because the run_simulation implementes scalar version.*/
void    nbody_implementation::NBodyVecAos::run_simulation_scalar(_In_ const double lo, const double hi) {

}

void    nbody_implementation::NBodyVecAos::display_state() const {

	std::cout << "****Dumping object state!!****" << std::endl;
	std::cout << "Object of type: " << typeid(NBodyVecAos).name() << std::endl;
	std::cout << "this at address: " << std::hex <<  this << std::endl;
	std::cout << "Started dump of Position Vectors values" << std::endl;
	std::cout << "X:              Y:              Z:" << std::endl;
	for (std::size_t i{ 0 }; i != this->m_uinParticles; ++i){

		std::cout << std::fixed << std::setprecision(10) << std::setw(4) << this->m_vNBodies[i].m_dX <<
			std::setw(17) << this->m_vNBodies[i].m_dY << std::setw(26) << this->m_vNBodies[i].m_dZ << std::endl;
	}
	std::cout << "****End of Position vectors values dump****" << std::endl;
	std::cout << "****End of object dump state****" << std::endl;
}

void    nbody_implementation::NBodyVecAos::dump_address_range()const {

	std::cout << "****Begining dump of Position vectors address range.****\n\n";
	std::cout << "&X             &Y            &Z" << std::endl;
	for (std::size_t i{ 0 }; i != this->m_uinParticles; ++i) {
		std::cout << std::setw(4) << std::hex << &this->m_vNBodies[i].m_dX <<
			std::setw(15) << std::hex << &this->m_vNBodies[i].m_dY <<
			std::setw(23) << std::hex << &this->m_vNBodies[i].m_dZ << std::endl;
	}
	std::cout << "****End of Position vectors address range dump.****" << std::endl;
}