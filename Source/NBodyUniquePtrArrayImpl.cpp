
#include "NbodyUniquePtrArrayImpl.h"

nbody_implementation::NBodyUptrArray::NBodyUptrArray(_In_ const std::size_t nParticles, _In_ const int nThreads, _In_ const double nStep) : 

m_uinParticles( (checkParticles(nParticles), nParticles) ),
m_inThreads( (checkThreads(nThreads),nThreads) ),
m_dStep( (checkStep(nStep),nStep)),
m_spX(new double[this->m_uinParticles], [](double* ptr)->void {delete[] ptr; }),
m_spY(new double[this->m_uinParticles], [](double* ptr)->void {delete[] ptr; }),
m_spZ(new double[this->m_uinParticles], [](double* ptr)->void {delete[] ptr; }),
m_spVX(new double[this->m_uinParticles], [](double* ptr)->void {delete[] ptr; }),
m_spVY(new double[this->m_uinParticles], [](double* ptr)->void {delete[] ptr; }),
m_spVZ(new double[this->m_uinParticles], [](double* ptr)->void {delete[] ptr; })
{

}

void  nbody_implementation::NBodyUptrArray::checkParticles(_In_ const std::size_t nParticles) {

	if (nParticles <= 128)
		throw std::runtime_error(std::string("Invalid input in: Ctor").append(typeid(NBodyUptrArray).name()));
}

void  nbody_implementation::NBodyUptrArray::checkThreads(_In_ const int nThreads) {

	if (nThreads <= 2)
		throw std::runtime_error(std::string("Invalid input in: Ctor").append(typeid(NBodyUptrArray).name()));
}

void  nbody_implementation::NBodyUptrArray::checkStep(_In_ const double nStep) {

	if (nStep == 0.0)
		throw std::runtime_error(std::string("Invalid input in: Ctor").append(typeid(NBodyUptrArray).name()));
}

/***********************************************************************
                        !! Warning
              Assumes that arrays size is divisable modulo 32.
************************************************************************/

void   nbody_implementation::NBodyUptrArray::run_simulation(_In_  double lo, _In_  double hi) {

	
	if (std::abs(lo - hi) <= std::numeric_limits<double>::epsilon()) {
		lo = 0.01;
		hi = 1.0;
	}
#if defined VERBOSE
	std::cout << "Allocated: " << static_cast<double>((6 * this->m_uinParticles * sizeof(double))) / 1024.0 << " KiB of heap memory" << std::endl;

#endif
	clock_t seed{ clock() };
	auto rand = std::bind(std::uniform_real_distribution < double > {lo, hi}, std::default_random_engine(seed));
	const std::string RNG_TYPE{ typeid(rand).name() };
#if defined VERBOSE
	std::cout << "Created Random number generator of type: " << RNG_TYPE.c_str() << "with the interval of randomness" << "[" << -1.0 << " -- " << 1.0 << "]" << std::endl;
#endif
	const __m256d v_dt = _mm256_set1_pd(this->m_dStep);
	const __m256d v_Zero = _mm256_setzero_pd();
#if defined VERBOSE
	std::cout << "Initialization of particles vector position started... ";
#endif
	
	for (std::size_t i{ 0 }; i != this->m_uinParticles; i += 4) {
		_mm256_storeu_pd(&this->m_spX.operator->()[i], _mm256_set_pd(rand(), rand(), rand(), rand()));
		_mm256_storeu_pd(&this->m_spY.operator->()[i], _mm256_set_pd(rand(), rand(), rand(), rand()));
		_mm256_storeu_pd(&this->m_spZ.operator->()[i], _mm256_set_pd(rand(), rand(), rand(), rand()));
	}
#if defined VERBOSE
	std::cout << "Done!!" << std::endl;
	std::cout << "Utilizing single CPU core " << std::endl;
#endif
	const double d_start{ ::omp_get_wtime() };

	for (std::size_t i{ 0 }; i != this->m_uinParticles; i += 4) {
		__m256d Fx(_mm256_setzero_pd());
		__m256d Fy(_mm256_setzero_pd());
		__m256d Fz(_mm256_setzero_pd());
		__m256d tx(_mm256_setzero_pd());
		__m256d ty(_mm256_setzero_pd());
		__m256d tz(_mm256_setzero_pd());

		for (std::size_t j{ 0 }; j != this->m_uinParticles; j += 4) {

			if (j != i) {

				/* Prfetch to L1D */
				_mm_prefetch(reinterpret_cast<const char*>(&this->m_spX.operator->()[j + 4]), _MM_HINT_T0);
				__m256d dx ( _mm256_sub_pd(_mm256_loadu_pd(&this->m_spX.operator->()[j]), _mm256_loadu_pd(&this->m_spX.operator->()[i])));
				/* Prefetch to L1D */
				_mm_prefetch(reinterpret_cast<const char*>(&this->m_spY.operator->()[j + 4]), _MM_HINT_T0);
				__m256d dy( _mm256_sub_pd(_mm256_loadu_pd(&this->m_spY.operator->()[j]), _mm256_loadu_pd(&this->m_spY.operator->()[i])));
				/* Prefetch to L1D */
				_mm_prefetch(reinterpret_cast<const char*>(&this->m_spZ.operator->()[j + 4]), _MM_HINT_T0);
				__m256d dz( _mm256_sub_pd(_mm256_loadu_pd(&this->m_spZ.operator->()[j]), _mm256_loadu_pd(&this->m_spZ.operator->()[i])));

				const __m256d dstSquared(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(dx, dx), _mm256_mul_pd(dy, dy)), _mm256_mul_pd(dz, dz)));
				__m256d inv_dstSquared(_mm256_castps_pd(_mm256_rsqrt_ps(_mm256_castpd_ps(dstSquared))));
				__m256d dstCubed(_mm256_mul_pd(_mm256_mul_pd(inv_dstSquared, inv_dstSquared), inv_dstSquared));

				Fx = _mm256_fmadd_pd(Fx, dstCubed, dx);
				Fy = _mm256_fmadd_pd(Fy, dstCubed, dy);
				Fz = _mm256_fmadd_pd(Fz, dstCubed, dz);

				tx = _mm256_add_pd(tx, Fx);
				ty = _mm256_add_pd(ty, Fy);
				tz = _mm256_add_pd(tz, Fz);

			}
		}
		_mm256_storeu_pd(&this->m_spVX.operator->()[i], _mm256_mul_pd(_mm256_hadd_pd(tx, tx), v_dt));
		_mm256_storeu_pd(&this->m_spVY.operator->()[i], _mm256_mul_pd(_mm256_hadd_pd(ty, ty), v_dt));
		_mm256_storeu_pd(&this->m_spVZ.operator->()[i], _mm256_mul_pd(_mm256_hadd_pd(tz, tz), v_dt));
	}

	for (std::size_t i{ 0 }; i != this->m_uinParticles; i += 8) {
		_mm256_storeu_pd(&this->m_spX.operator->()[i], _mm256_add_pd(_mm256_load_pd(&this->m_spVX.operator->()[i]), v_dt));
		_mm256_storeu_pd(&this->m_spX.operator->()[i + 4], _mm256_add_pd(_mm256_load_pd(&this->m_spVX.operator->()[i + 4]), v_dt));
		_mm256_storeu_pd(&this->m_spY.operator->()[i], _mm256_add_pd(_mm256_load_pd(&this->m_spVY.operator->()[i]), v_dt));
		_mm256_storeu_pd(&this->m_spY.operator->()[i + 4], _mm256_add_pd(_mm256_load_pd(&this->m_spVY.operator->()[i + 4]), v_dt));
		_mm256_storeu_pd(&this->m_spZ.operator->()[i], _mm256_add_pd(_mm256_load_pd(&this->m_spVZ.operator->()[i]), v_dt));
		_mm256_storeu_pd(&this->m_spZ.operator->()[i + 4], _mm256_add_pd(_mm256_load_pd(&this->m_spVZ.operator->()[i + 4]), v_dt));
	}
	const double d_end{ ::omp_get_wtime() };
#if defined VERBOSE
	std::cout << "Crude timing result: " << d_end - d_start << " sec, " << (d_end - d_start) * 1.0E+9 << " ns" << std::endl;
#endif

}

void  nbody_implementation::NBodyUptrArray::display_state() const {

	std::cout << " Dumping state of postion vectors\n\n";
	std::cout << "Pos. X" << "Pos. Y" << "Pos. Z" << std::endl;
	for (int i{ 0 }; i != this->m_uinParticles; ++i) {
		std::cout << std::fixed << std::showpoint << std::setprecision(10) << std::setw(4) <<
			this->m_spX.operator->()[i] << std::setw(18) << this->m_spY.operator->()[i] << std::setw(25) << this->m_spZ.operator->()[i] << std::endl;
	}
	std::cout << "Reached end of dump" << std::endl;
}

void  nbody_implementation::NBodyUptrArray::run_simulation_omp(_In_ const double lo, _In_ const double hi) {

	if ((std::fabs(lo - hi)) <= std::numeric_limits<double>::epsilon())
		throw std::runtime_error(std::string("Invalid input in: NBodyDArray::run_simualtion_omp\n"));

#if defined VERBOSE
	std::cout << "Allocated: " << static_cast<double>((6 * this->m_uinParticles * sizeof(double))) / 1024.0 << " KiB of heap memory" << std::endl;

#endif
	clock_t seed{ clock() };
	auto rand = std::bind(std::uniform_real_distribution < double > {lo, hi}, std::default_random_engine(seed));
	const std::string RNG_TYPE{ typeid(rand).name() };
#if defined VERBOSE
	std::cout << "Created Random number generator of type: " << RNG_TYPE.c_str() << "with the interval of randomness" << "[" << -1.0 << " -- " << 1.0 << "]" << std::endl;
#endif
	const __m256d v_dt = _mm256_set1_pd(this->m_dStep);
	const __m256d v_Zero = _mm256_setzero_pd();
	const std::size_t nParticles{ this->m_uinParticles };
	const int nThreads{ this->m_inThreads };
	auto const  chunk_size = this->m_uinParticles / this->m_inThreads;
	/* Allocating automatic smart pointers for OpenMP threading */
	std::unique_ptr<double, void(*)(double*)> a_spX(new double[this->m_uinParticles], [](double* ptr)->void {delete[] ptr; });
	std::unique_ptr<double, void(*)(double*)> a_spY(new double[this->m_uinParticles], [](double* ptr)->void {delete[] ptr; });
	std::unique_ptr<double, void(*)(double*)> a_spZ(new double[this->m_uinParticles], [](double* ptr)->void {delete[] ptr; });
	std::unique_ptr<double, void(*)(double*)> a_spVX(new double[this->m_uinParticles], [](double* ptr)->void {delete[] ptr; });
	std::unique_ptr<double, void(*)(double*)> a_spVY(new double[this->m_uinParticles], [](double* ptr)->void {delete[] ptr; });
	std::unique_ptr<double, void(*)(double*)> a_spVZ(new double[this->m_uinParticles], [](double* ptr)->void {delete[] ptr; });
	/* Initialization of Position vectors 
	   using unaligned store*/
	for (std::size_t i{ 0 }; i != this->m_uinParticles; i += 4) {
		// Position vector X
		_mm256_storeu_pd(&a_spX.operator->()[i], _mm256_set_pd(rand(), rand(), rand(), rand()));
		// Position vector Y
		_mm256_storeu_pd(&a_spY.operator->()[i], _mm256_set_pd(rand(), rand(), rand(), rand()));
		// Position vector Z
		_mm256_storeu_pd(&a_spZ.operator->()[i], _mm256_set_pd(rand(), rand(), rand(), rand()));
	}
#if defined VERBOSE
	std::cout << "Done!!" << std::endl;
	std::cout << "OpenMP in use with: " << this->m_inThreads << " threads" << std::endl;
#endif
	::omp_set_num_threads(nThreads);
	const double d_start_time{ ::omp_get_wtime() };
#pragma omp parallel for schedule(static,chunk_size)
	for (std::size_t i = 0; i < nParticles; i += 4) {
		__m256d Fx(_mm256_setzero_pd());
		__m256d Fy(_mm256_setzero_pd());
		__m256d Fz(_mm256_setzero_pd());
		__m256d tx(_mm256_setzero_pd());
		__m256d ty(_mm256_setzero_pd());
		__m256d tz(_mm256_setzero_pd());

		for (std::size_t j = 0; j < nParticles; j += 4) {

			if (j != i) {

				// Begin prefetching to L1D only.
				_mm_prefetch(reinterpret_cast<const char*>(&a_spX.operator->()[j + 4]), _MM_HINT_T0);
				__m256d dx(_mm256_sub_pd(_mm256_loadu_pd(&a_spX.operator->()[j]), _mm256_loadu_pd(&a_spX.operator->()[i])));

				_mm_prefetch(reinterpret_cast<const char*>(&a_spY.operator->()[j + 4]), _MM_HINT_T0);
				__m256d dy(_mm256_sub_pd(_mm256_loadu_pd(&a_spY.operator->()[j]), _mm256_loadu_pd(&a_spY.operator->()[i])));

				_mm_prefetch(reinterpret_cast<const char*>(&a_spZ.operator->()[j + 4]), _MM_HINT_T0);
				__m256d dz(_mm256_sub_pd(_mm256_loadu_pd(&a_spZ.operator->()[j]), _mm256_loadu_pd(&a_spZ.operator->()[i])));

				const __m256d dstSquared(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(dx, dx), _mm256_mul_pd(dy, dy)), _mm256_mul_pd(dz, dz)));
				__m256d inv_dstSquared(_mm256_castps_pd(_mm256_rsqrt_ps(_mm256_castpd_ps(dstSquared))));
				__m256d dstCubed(_mm256_mul_pd(_mm256_mul_pd(inv_dstSquared, inv_dstSquared), inv_dstSquared));

				Fx = _mm256_fmadd_pd(Fx, dstCubed, dx);
				Fy = _mm256_fmadd_pd(Fy, dstCubed, dy);
				Fz = _mm256_fmadd_pd(Fz, dstCubed, dz);

				tx = _mm256_add_pd(tx, Fx);
				ty = _mm256_add_pd(ty, Fy);
				tz = _mm256_add_pd(tz, Fz);
			}
		}
		_mm256_storeu_pd(&a_spVX.operator->()[i], _mm256_mul_pd(_mm256_hadd_pd(tx, tx), v_dt));
		_mm256_storeu_pd(&a_spVY.operator->()[i], _mm256_mul_pd(_mm256_hadd_pd(ty, ty), v_dt));
		_mm256_storeu_pd(&a_spVZ.operator->()[i], _mm256_mul_pd(_mm256_hadd_pd(tz, tz), v_dt));
	}
	for (int i{ 0 }; i != this->m_uinParticles; i += 8) {

		_mm256_storeu_pd(&a_spX.operator->()[i], _mm256_add_pd(_mm256_loadu_pd(&a_spVX.operator->()[i]), v_dt));
		_mm256_storeu_pd(&a_spX.operator->()[i + 4], _mm256_add_pd(_mm256_loadu_pd(&a_spVX.operator->()[i + 4]), v_dt));
		_mm256_storeu_pd(&a_spY.operator->()[i], _mm256_add_pd(_mm256_loadu_pd(&a_spVY.operator->()[i]), v_dt));
		_mm256_storeu_pd(&a_spY.operator->()[i + 4], _mm256_add_pd(_mm256_loadu_pd(&a_spVY.operator->()[i + 4]), v_dt));
		_mm256_storeu_pd(&a_spZ.operator->()[i], _mm256_add_pd(_mm256_loadu_pd(&a_spVZ.operator->()[i]), v_dt));
		_mm256_storeu_pd(&a_spZ.operator->()[i + 4], _mm256_add_pd(_mm256_loadu_pd(&a_spVZ.operator->()[i + 4]), v_dt));
	}
	const double d_stop_time{ ::omp_get_wtime() };
#if defined VERBOSE
	std::cout << "Crude timing result: " << std::showpoint << std::fixed << d_stop_time - d_start_time << " sec, " << (d_stop_time - d_start_time) * 1.0E+9 << " ns" << std::endl;
#endif
	// Copy to object state.
	// Using fast vectorised copy instead of scalar copy loop as implemented by the c
	// compiler while calling std::copy
	// Unrolling 2x.
	for (std::size_t i{ 0 }; i != this->m_uinParticles; i += 8) {
		_mm256_storeu_pd(&this->m_spVX.operator->()[i], _mm256_loadu_pd(&a_spVX.operator->()[i]));
		_mm256_storeu_pd(&this->m_spVX.operator->()[i + 4], _mm256_loadu_pd(&a_spVX.operator->()[i + 4]));
		_mm256_storeu_pd(&this->m_spVY.operator->()[i], _mm256_loadu_pd(&a_spVY.operator->()[i]));
		_mm256_storeu_pd(&this->m_spVY.operator->()[i + 4], _mm256_loadu_pd(&a_spVY.operator->()[i + 4]));
		_mm256_storeu_pd(&this->m_spVZ.operator->()[i], _mm256_loadu_pd(&a_spVZ.operator->()[i]));
		_mm256_storeu_pd(&this->m_spVZ.operator->()[i + 4], _mm256_loadu_pd(&a_spVZ.operator->()[i + 4]));
	}
}

void  nbody_implementation::NBodyUptrArray::run_simulation_scalar(_In_ const double lo, _In_ const double hi) {
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
	for (std::size_t i{ 0 }; i != this->m_uinParticles; ++i) {

		this->m_spX.operator->()[i] = rand();
		this->m_spY.operator->()[i] = rand();
		this->m_spZ.operator->()[i] = rand();
	}
	std::cout << "Finished initialization of " << this->m_uinParticles << "Particles with:" << typeid(rand).name() << std::endl;
#ifdef EXECUTE_NTIMES
	std::cout << "SStarting excuting NBody system: " << NTIMES << std::endl;
	for (int k{ 0 }; k != NTIMES; ++k) {
#endif
		for (std::size_t i{ 0 }; i != this->m_uinParticles; ++i) {
			double Fx{ 0.0 };
			double Fy{ 0.0 };
			double Fz{ 0.0 };

			for (std::size_t j{ 0 }; j != this->m_uinParticles; ++j) {

				if (j != i) {

					_mm_prefetch(reinterpret_cast<const char*>(&this->m_spX.operator->()[j + 7]), _MM_HINT_T0);
					const double dx{ this->m_spX.operator->()[j] - this->m_spX.operator->()[i] };

					_mm_prefetch(reinterpret_cast<const char*>(&this->m_spY.operator->()[j + 7]) , _MM_HINT_T0);
					const double dy{ this->m_spVY.operator->()[j] - this->m_spY.operator->()[i] };

					_mm_prefetch(reinterpret_cast<const char*>(&this->m_spZ.operator->()[j + 7]), _MM_HINT_T0);
					const double dz{ this->m_spZ.operator->()[j] - this->m_spZ.operator->()[i] };

					const double distSquared{ dx * dx + dy * dy + dz * dz };
					const double invdistCube{ 1.0 / (distSquared * ::sqrt(distSquared)) };

					Fx += dx * invdistCube;
					Fy += dy * invdistCube;
					Fz += dz * invdistCube;
				}
			}
			this->m_spVX.operator->()[i] = Fx * this->m_dStep;
			this->m_spVY.operator->()[i] = Fy * this->m_dStep;
			this->m_spVZ.operator->()[i] = Fz * this->m_dStep;
		}
		for (std::size_t i{ 0 }; i != this->m_uinParticles; ++i) {
			this->m_spX.operator->()[i] = this->m_spVX.operator->()[i] + this->m_dStep;
			this->m_spY.operator->()[i] = this->m_spVY.operator->()[i] + this->m_dStep;
			this->m_spZ.operator->()[i] = this->m_spVZ.operator->()[i] + this->m_dStep;
		}

#ifdef EXECUTE_NTIMES
	}
		std::cout << "Finished executing NBody system: " << NTIMES << std::endl;
#endif
}

void   nbody_implementation::NBodyUptrArray::dump_address_range() const {

	std::cout << "****Begining dump of Position vectors address range.****\n\n";
	std::cout << "&X             &Y            &Z" << std::endl;
	for (std::size_t i{ 0 }; i != this->m_uinParticles; ++i) {
		std::cout << std::setw(4) << std::hex << &this->m_spX.operator->()[i] <<
			std::setw(15) << std::hex << &this->m_spY.operator->()[i] <<
			std::setw(23) << std::hex << &this->m_spZ.operator->()[i] << std::endl;
	}
	std::cout << "****End of Position vectors address range dump.****" << std::endl;
}