#include "NBodyVectorImpl.h"

nbody_implementation::NBodyVec::NBodyVec(_In_ const std::size_t nParticles, _In_ const int nThreads,
	_In_ const double Step) : 
	m_uinParticles{nParticles},
	m_inThreads{nThreads},
	m_dStep{Step}
	
	
{ 
	if ((this->m_uinParticles <= 128) || (this->m_inThreads <= 2) || 
		(this->m_dStep == 0.0) || ((this->m_uinParticles % 32) != 0))
		throw std::runtime_error(std::string("Invalid input in: ").append( typeid(NBodyVec).name()));

	this->m_vX = std::vector<double>(this->m_uinParticles);
	this->m_vY = std::vector<double>(this->m_uinParticles);
	this->m_vZ = std::vector<double>(this->m_uinParticles);
	this->m_vVX = std::vector<double>(this->m_uinParticles);
	this->m_vVY = std::vector<double>(this->m_uinParticles);
	this->m_vVZ = std::vector<double>(this->m_uinParticles);

#if VERBOSE
	nbody_vec_internal_state();
#endif 
}

nbody_implementation::NBodyVec::NBodyVec(_In_ const NBodyVec &rhs) :
m_uinParticles{ rhs.m_uinParticles },
m_inThreads{ rhs.m_inThreads },
m_dStep{ rhs.m_dStep },
m_vX{ rhs.m_vX },
m_vY{ rhs.m_vY },
m_vZ{ rhs.m_vZ },
m_vVX{ rhs.m_vVX },
m_vVY{ rhs.m_vY },
m_vVZ{ rhs.m_vZ }
{

}

nbody_implementation::NBodyVec::NBodyVec(_In_ NBodyVec &&rhs) :
m_uinParticles{ std::move(rhs.m_uinParticles) },
m_inThreads{ std::move(rhs.m_inThreads) },
m_dStep{ std::move(rhs.m_dStep) },
m_vX{ std::move(rhs.m_vX) },
m_vY{ std::move(rhs.m_vY) },
m_vZ{ std::move(rhs.m_vZ) },
m_vVX{ std::move(rhs.m_vVX) },
m_vVY{ std::move(rhs.m_vVY) },
m_vVZ{ std::move(rhs.m_vVZ) }
{

}

void   nbody_implementation::NBodyVec::run_simulation(_In_ const double lo, _In_ const double hi) {

	if ((std::fabs(lo - hi)) <= std::numeric_limits<double>::epsilon())
		throw std::runtime_error(std::string("Invalid input in: NBodyVec::run_simualtion\n"));

#if defined VERBOSE
	std::cout << "Allocated: " << static_cast<double>((6 * this->m_uinParticles * sizeof(double))) / 1024.0 << " KB's of heap memory" << std::endl;

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
	for (int i{ 0 }; i != this->m_uinParticles; i += 4) {

		_mm256_store_pd(&this->m_vX[i], _mm256_set_pd(rand(), rand(), rand(), rand()));
		_mm256_store_pd(&this->m_vY[i], _mm256_set_pd(rand(), rand(), rand(), rand()));
		_mm256_store_pd(&this->m_vZ[i], _mm256_set_pd(rand(), rand(), rand(), rand()));
	}
#if defined VERBOSE
	std::cout << "Done!!" << std::endl;
	std::cout << "Utilizing single CPU core " << std::endl;
#endif
	const double d_start{ ::omp_get_wtime() };

	for (int i{ 0 }; i != this->m_uinParticles; i += 4) {
		__m256d Fx = _mm256_setzero_pd();
		__m256d Fy = _mm256_setzero_pd();
		__m256d Fz = _mm256_setzero_pd();
		__m256d tx = _mm256_setzero_pd();
		__m256d ty = _mm256_setzero_pd();
		__m256d tz = _mm256_setzero_pd();

		for (int j{ 0 }; j != this->m_uinParticles; j += 4) {

			if (j != i) {

				_mm_prefetch(reinterpret_cast<const char*>(&this->m_vX[j + 2]), _MM_HINT_T0);
				__m256d dx(_mm256_sub_pd(_mm256_load_pd(&this->m_vX[j]), _mm256_load_pd(&this->m_vX[i])));
				
				_mm_prefetch(reinterpret_cast<const char*>(&this->m_vY[j + 2]), _MM_HINT_T0);
				__m256d dy(_mm256_sub_pd(_mm256_load_pd(&this->m_vY[j]), _mm256_load_pd(&this->m_vY[i])));

				_mm_prefetch(reinterpret_cast<const char*>(&this->m_vZ[j + 2]), _MM_HINT_T0);
				__m256d dz(_mm256_sub_pd(_mm256_load_pd(&this->m_vZ[j]), _mm256_load_pd(&this->m_vZ[i])));

				const __m256d dstSquared( _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(dx, dx), _mm256_mul_pd(dy, dy)), _mm256_mul_pd(dz, dz)));
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
		_mm256_store_pd(&this->m_vVX[i], _mm256_mul_pd(_mm256_hadd_pd(tx, tx), v_dt));
		_mm256_store_pd(&this->m_vVY[i], _mm256_mul_pd(_mm256_hadd_pd(ty, ty), v_dt));
		_mm256_store_pd(&this->m_vVZ[i], _mm256_mul_pd(_mm256_hadd_pd(tz, tz), v_dt));
		
	}
	for (int i{ 0 }; i != this->m_uinParticles; i += 8) {

		_mm256_store_pd(&this->m_vX[i], _mm256_add_pd(_mm256_load_pd(&this->m_vVX[i]), v_dt));
		_mm256_store_pd(&this->m_vX[i + 4], _mm256_add_pd(_mm256_load_pd(&this->m_vVX[i + 4]), v_dt));
		_mm256_store_pd(&this->m_vY[i], _mm256_add_pd(_mm256_load_pd(&this->m_vVY[i]), v_dt));
		_mm256_store_pd(&this->m_vY[i + 4], _mm256_add_pd(_mm256_load_pd(&this->m_vVY[i + 4]), v_dt));
		_mm256_store_pd(&this->m_vZ[i], _mm256_add_pd(_mm256_load_pd(&this->m_vVZ[i]), v_dt));
		_mm256_store_pd(&this->m_vZ[i + 4], _mm256_add_pd(_mm256_load_pd(&this->m_vVZ[i + 4]), v_dt));
	}
	const double d_end{ ::omp_get_wtime() };
#if defined VERBOSE
	std::cout << "Crude timing result: " << d_end - d_start << " sec, " << (d_end - d_start) * 1.0E+9 << " ns" << std::endl;
#endif

}

void      nbody_implementation::NBodyVec::display_state() const {

	std::cout << " Dumping state of postion vectors\n\n";
	std::cout << "Pos. X"          << "Pos. Y"        << "Pos. Z" << std::endl;
	for (int i{ 0 }; i != this->m_uinParticles; ++i) {
		std::cout << std::fixed << std::showpoint << std::setprecision(10) << std::setw(4) <<
			this->m_vX[i] << std::setw(18) << this->m_vY[i] << std::setw(25) << this->m_vZ[i] << std::endl;
	}
	std::cout << "Reached end of dump" << std::endl;
}

void      nbody_implementation::NBodyVec::run_simulation_omp(_In_ const double lo, _In_ const double hi) {
	if ((std::fabs(lo - hi)) <= std::numeric_limits<double>::epsilon())
		throw std::runtime_error(std::string("Invalid input in NBodyVec::run_simulation_omp\n"));

#if defined VERBOSE
	std::cout << "Allocated: " << static_cast<double>((6 * this->m_uinParticles * sizeof(double))) / 1024 << " KB's of heap memory" << std::endl;

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
	std::vector<double> a_vX( this->m_vX.size() );
	std::vector<double> a_vY(this->m_vY.size() );
	std::vector<double> a_vZ( this->m_vZ.size() );
	std::vector<double> a_vXV( this->m_vVX.size()) ;
	std::vector<double> a_vYV( this->m_vVY.size());
	std::vector<double> a_vZV( this->m_vVZ.size() );
	auto const chunk_size = nParticles / nThreads;
#if defined VERBOSE
	std::cout << "Initialization of position vectors started... ";
#endif
	for (int i{ 0 }; i != this->m_uinParticles; i += 4) {

		_mm256_store_pd(&a_vX[i], _mm256_set_pd(rand(), rand(), rand(), rand()));
		_mm256_store_pd(&a_vY[i], _mm256_set_pd(rand(), rand(), rand(), rand()));
		_mm256_store_pd(&a_vZ[i], _mm256_set_pd(rand(), rand(), rand(), rand()));
	}
#if defined VERBOSE
	std::cout << "Done!!" << std::endl;
	std::cout << "OpenMP in use with: " << this->m_inThreads << " threads" << std::endl;
#endif
	::omp_set_num_threads(nThreads);
	const double d_time_start{ ::omp_get_wtime() };
#pragma omp parallel for schedule(static, chunk_size)
	for (int i = 0; i < nParticles; i += 4) {
		__m256d Fx = _mm256_setzero_pd();
		__m256d Fy = _mm256_setzero_pd();
		__m256d Fz = _mm256_setzero_pd();
		__m256d tx = _mm256_setzero_pd();
		__m256d ty = _mm256_setzero_pd();
		__m256d tz = _mm256_setzero_pd();

		for (int j = 0; j < nParticles; j += 4) {

			if (j != i) {

				_mm_prefetch(reinterpret_cast<const char*>(&a_vX[j + 2]), _MM_HINT_T0);
				__m256d dx(_mm256_sub_pd(_mm256_load_pd(&a_vX[j]), _mm256_load_pd(&a_vX[i])));

				_mm_prefetch(reinterpret_cast<const char*>(&a_vY[j + 2]), _MM_HINT_T0);
				__m256d dy(_mm256_sub_pd(_mm256_load_pd(&a_vY[j]), _mm256_load_pd(&a_vY[i])));

				_mm_prefetch(reinterpret_cast<const char*>(&a_vZ[j + 2]), _MM_HINT_T0);
				__m256d dz(_mm256_sub_pd(_mm256_load_pd(&a_vZ[j + 2]), _mm256_load_pd(&a_vZ[i])));

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
		_mm256_store_pd(&a_vXV[i], _mm256_mul_pd(_mm256_hadd_pd(tx, tx), v_dt));
		_mm256_store_pd(&a_vYV[i], _mm256_mul_pd(_mm256_hadd_pd(ty, ty), v_dt));
		_mm256_store_pd(&a_vZV[i], _mm256_mul_pd(_mm256_hadd_pd(tz, tz), v_dt));
	}
	for (int i{ 0 }; i != this->m_uinParticles; i += 8) {
		
		_mm256_store_pd(&a_vX[i], _mm256_add_pd(_mm256_load_pd(&a_vXV[i]), v_dt));
		_mm256_store_pd(&a_vX[i + 4], _mm256_add_pd(_mm256_load_pd(&a_vX[i + 4]), v_dt));
		_mm256_store_pd(&a_vY[i], _mm256_add_pd(_mm256_load_pd(&a_vYV[i]), v_dt));
		_mm256_store_pd(&a_vY[i+4], _mm256_add_pd(_mm256_load_pd(&a_vYV[i+4]), v_dt));
		_mm256_store_pd(&a_vZ[i], _mm256_add_pd(_mm256_load_pd(&a_vZV[i]), v_dt));
		_mm256_store_pd(&a_vZ[i+4], _mm256_add_pd(_mm256_load_pd(&a_vZV[i+4]), v_dt));
	}
	const double d_stop_time{ ::omp_get_wtime() };
#if defined VERBOSE
	std::cout << "Crude timing result: " << std::showpoint << std::fixed << d_stop_time - d_time_start << " sec, " << (d_stop_time - d_time_start) * 1.0E+9 << " ns" << std::endl;
#endif
	this->m_vX.operator=(a_vX);
	this->m_vY.operator=(a_vY);
	this->m_vZ.operator=(a_vZ);
}

/* Implement it later */
auto   nbody_implementation::NBodyVec::run_simulation_scalar(_In_ const double lo, _In_ const double hi)->void {

}

#if defined VERBOSE
auto  nbody_implementation::NBodyVec::nbody_vec_internal_state()const ->void {
	std::cout << typeid(NBodyVec).name() << " executing now at loc:" << __LINE__ << "at address: " << this << std::endl;
	std::cout << typeid(NBodyVec).name() << " members at: " << std::endl;
	std::cout << "m_vX at: " << "0x" << &this->m_vX << std::endl;
	std::cout << "m_vY at: " << "0x" << &this->m_vY << std::endl;
	std::cout << "m_vZ at: " << "0x" << &this->m_vZ << std::endl;
	std::cout << "m_vVX at:" << "0x" << &this->m_vVX << std::endl;
	std::cout << "m_vVY at:" << "0x" << &this->m_vVY << std::endl;
	std::cout << "m_vVZ at:" << "0x" << &this->m_vVZ << std::endl;
} 
#endif
