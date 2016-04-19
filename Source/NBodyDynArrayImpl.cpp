
#include "NBodyDynArrayImpl.h"

nbody_implementation::NBodyDArray::NBodyDArray(_In_ const std::size_t nParticles, _In_ const int nThreads, _In_ const double nSteps) :
m_uinParticles{ nParticles },
m_inThreads{ nThreads },
m_dStep{ nSteps }
{
	// Exit on any failure deallocate already allocated memory
	// in order to prevent leak.
	// Object state consistent only partially before executing catch block.
	try {
		this->m_pX = new double[this->m_uinParticles];
		this->m_pY = new double[this->m_uinParticles];
		this->m_pZ = new double[this->m_uinParticles];
		this->m_pVX = new double[this->m_uinParticles];
		this->m_pVY = new double[this->m_uinParticles];
		this->m_pVZ = new double[this->m_uinParticles];
	}
	catch (std::bad_alloc& ba){
		delete[] this->m_pX; this->m_pX = nullptr;
		
		delete[] this->m_pY; this->m_pY = nullptr;
		delete[] this->m_pZ; this->m_pZ = nullptr;
		delete[] this->m_pVX; this->m_pVX = nullptr;
		delete[] this->m_pVY; this->m_pVY = nullptr;
		delete[] this->m_pVZ; this->m_pVZ = nullptr;
		std::cout << "General memory allocation failure at loc #: " << __LINE__ << std::endl;
		std::cout << " In:" << typeid(NBodyDArray).name() << "NBodyDArray" << std::endl;
		std::cerr << ba.what() << std::endl;
		std::cout << "!!Critical failure calling... exit(-1) now!!" << std::endl;
		std::exit(EXIT_FAILURE);

	}
#ifdef VERBOSE
	std::cout << "Executing in NBodyDArray::NBodyDArray(std::size_t,int,double)" << std::endl;
#endif
}

nbody_implementation::NBodyDArray::NBodyDArray(_In_ const NBodyDArray &rhs) :
m_uinParticles{ rhs.m_uinParticles },
m_inThreads{ rhs.m_inThreads },
m_dStep{ rhs.m_dStep }

{
	std::copy(rhs.m_pX, rhs.m_pX + rhs.m_uinParticles, this->m_pX);
	std::copy(rhs.m_pY, rhs.m_pY + rhs.m_uinParticles, this->m_pY);
	std::copy(rhs.m_pZ, rhs.m_pZ + rhs.m_uinParticles, this->m_pZ);
	std::copy(rhs.m_pVX, rhs.m_pVX + rhs.m_uinParticles, this->m_pVX);
	std::copy(rhs.m_pVY, rhs.m_pVY + rhs.m_uinParticles, this->m_pVY);
	std::copy(rhs.m_pVZ, rhs.m_pVZ + rhs.m_uinParticles, this->m_pVZ);
}

nbody_implementation::NBodyDArray::NBodyDArray(_In_ NBodyDArray &&rhs) :
m_uinParticles{ std::move(rhs.m_uinParticles) },
m_inThreads{ std::move(rhs.m_inThreads) },
m_dStep{ std::move(rhs.m_dStep) }
{
	// Nullify the state.
	rhs.m_dStep = 0.0;
	rhs.m_inThreads = 0;
	rhs.m_uinParticles = 0;
	std::move(rhs.m_pX, rhs.m_pX + rhs.m_uinParticles, this->m_pX);
	delete[] rhs.m_pX;
	rhs.m_pX = nullptr;
	std::move(rhs.m_pY, rhs.m_pY + rhs.m_uinParticles, this->m_pY);
	delete[] rhs.m_pY;
	rhs.m_pY = nullptr;
	std::move(rhs.m_pZ, rhs.m_pZ + rhs.m_uinParticles, this->m_pZ);
	delete[] rhs.m_pZ;
	rhs.m_pZ = nullptr;
	std::move(rhs.m_pVX, rhs.m_pVX + rhs.m_uinParticles, this->m_pVX);
	delete[] rhs.m_pVX;
	rhs.m_pVX = nullptr;
	std::move(rhs.m_pVY, rhs.m_pVY + rhs.m_uinParticles, this->m_pVY);
	delete[] rhs.m_pVY;
	rhs.m_pVY = nullptr;
	std::move(rhs.m_pVZ, rhs.m_pVZ + rhs.m_uinParticles, this->m_pVZ);
	delete[] rhs.m_pVZ;
	rhs.m_pVZ = nullptr;

}

nbody_implementation::NBodyDArray::~NBodyDArray() {
#ifdef VERBOSE
	std::cout << "Executing in: NBodyDArray::~NBodyDArray()" << std::endl;
#endif
	
	 if(this->m_pX) delete[] this->m_pX; this->m_pX = nullptr;
	 if(this->m_pY) delete[] this->m_pY; this->m_pY = nullptr;
	 if(this->m_pZ) delete[] this->m_pZ; this->m_pZ = nullptr;
	 if(this->m_pVX) delete[] this->m_pVX; this->m_pVX = nullptr;
	 if(this->m_pVY) delete[] this->m_pVY; this->m_pVY = nullptr;
	 if(this->m_pVZ) delete[] this->m_pVZ; this->m_pVZ = nullptr;



}
/***********************************************************************
                       !! Warning
          Assumes that arrays size is divisable modulo 32.
************************************************************************/
void  nbody_implementation::NBodyDArray::run_simulation(_In_ const double lo, _In_ const double hi) {

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
	const __m256d v_dt = _mm256_set1_pd(this->m_dStep);
	const __m256d v_Zero = _mm256_setzero_pd();
#if defined VERBOSE
	std::cout << "Initialization of particles vector position started... ";
#endif
	for (std::size_t i{ 0 }; i != this->m_uinParticles; i += 4) {

		_mm256_storeu_pd(&this->m_pX[i], _mm256_set_pd(rand(), rand(), rand(), rand()));
		_mm256_storeu_pd(&this->m_pY[i], _mm256_set_pd(rand(), rand(), rand(), rand()));
		_mm256_storeu_pd(&this->m_pZ[i], _mm256_set_pd(rand(), rand(), rand(), rand()));
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

				_mm_prefetch(reinterpret_cast<const char*>(&this->m_pX[j + 2]), _MM_HINT_T0);
				__m256d dx(_mm256_sub_pd(_mm256_load_pd(&this->m_pX[j]), _mm256_load_pd(&this->m_pX[i])));

				_mm_prefetch(reinterpret_cast<const char*>(&this->m_pY[j + 2]), _MM_HINT_T0);
				__m256d dy(_mm256_sub_pd(_mm256_load_pd(&this->m_pY[j]), _mm256_load_pd(&this->m_pY[i])));

				_mm_prefetch(reinterpret_cast<const char*>(&this->m_pZ[j + 2]), _MM_HINT_T0);
				__m256d dz(_mm256_sub_pd(_mm256_load_pd(&this->m_pZ[j]), _mm256_load_pd(&this->m_pZ[i])));

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
		_mm256_store_pd(&this->m_pVX[i], _mm256_mul_pd(_mm256_hadd_pd(tx, tx), v_dt));
		_mm256_store_pd(&this->m_pVY[i], _mm256_mul_pd(_mm256_hadd_pd(ty, ty), v_dt));
		_mm256_store_pd(&this->m_pVZ[i], _mm256_mul_pd(_mm256_hadd_pd(tz, tz), v_dt));

	}
	for (int i{ 0 }; i != this->m_uinParticles; i += 8) {

		_mm256_store_pd(&this->m_pX[i], _mm256_add_pd(_mm256_load_pd(&this->m_pVX[i]), v_dt));
		_mm256_store_pd(&this->m_pX[i + 4], _mm256_add_pd(_mm256_load_pd(&this->m_pVX[i + 4]), v_dt));
		_mm256_store_pd(&this->m_pY[i], _mm256_add_pd(_mm256_load_pd(&this->m_pVY[i]), v_dt));
		_mm256_store_pd(&this->m_pY[i + 4], _mm256_add_pd(_mm256_load_pd(&this->m_pVY[i + 4]), v_dt));
		_mm256_store_pd(&this->m_pZ[i], _mm256_add_pd(_mm256_load_pd(&this->m_pVZ[i]), v_dt));
		_mm256_store_pd(&this->m_pZ[i + 4], _mm256_add_pd(_mm256_load_pd(&this->m_pVZ[i + 4]), v_dt));
	}
	const double d_end{ ::omp_get_wtime() };
#if defined VERBOSE
	std::cout << "Crude timing result: " << d_end - d_start << " sec, " << (d_end - d_start) * 1.0E+9 << " ns" << std::endl;
#endif

}

void   nbody_implementation::NBodyDArray::display_state() const {

	std::cout << " Dumping state of postion vectors\n\n";
	std::cout << "Pos. X" << "Pos. Y" << "Pos. Z" << std::endl;
	for (int i{ 0 }; i != this->m_uinParticles; ++i) {
		std::cout << std::fixed << std::showpoint << std::setprecision(10) << std::setw(4) <<
			this->m_pX[i] << std::setw(18) << this->m_pY[i] << std::setw(25) << this->m_pZ[i] << std::endl;
	}
	std::cout << "Reached end of dump" << std::endl;
}

void    nbody_implementation::NBodyDArray::dump_address_range() const {

	std::cout << "Dumping address range of position vectora\n\n" << std::endl;
	std::cout << " &X " << " &Y " << " &Z " << std::endl;
	for (int i{ 0 }; i != this->m_uinParticles; ++i) {
		std::cout << std::setw(4) << std::hex << &this->m_pX[i] << std::setw(18) << std::hex << &this->m_pY[i] <<
			std::setw(24) << std::hex << &this->m_pZ[i] << std::endl;
	}
	std::cout << " End of memory address range dump " << std::endl;
			
	
}

/*double  *nbody_implementation::NBodyDArray::allocate_array1D(_In_ const std::size_t dataSize) {

	if (dataSize <= 0)
		throw std::runtime_error(std::string{ "Invalid input in: " }.append(typeid(NBodyDArray).name()));

	double* pArray1D = nullptr;
	try{
		pArray1D = new double[dataSize];
	}
	catch (std::bad_alloc& ba) {
		std::cerr << "Fatal Error in: " << typeid(NBodyDArray).name() << "::allocate_array1D, calling exit..." << std::endl;
		std::cerr << ba.what() << std::endl;
		delete[] pArray1D;
		pArray1D = nullptr;
		std::exit(EXIT_FAILURE);
	}
	return (pArray1D);
}*/

std::tuple<double*, double*, double*, double*, double*, double*> nbody_implementation::NBodyDArray::allocate_array1D(_In_ const std::size_t nParticles) {
	if (nParticles <= 0)
		throw std::runtime_error(std::string{ "Invalid input in: " }.append(typeid(NBodyDArray).name()));
	
	double* pAr1, *pAr2, *pAr3, *pAr4, *pAr5, *pAr6;
	
	
	try{
		pAr1 = new double[nParticles];
		pAr2 = new double[nParticles];
		pAr3 = new double[nParticles];
		pAr4 = new double[nParticles];
		pAr5 = new double[nParticles];
		pAr6 = new double[nParticles];
		
	}
	catch (std::bad_alloc& ba) {
		std::cerr << "Fatal Error in: " << typeid(NBodyDArray).name() << "::allocate_array1D, calling exit..." << std::endl;
		std::cerr << ba.what() << std::endl;
		delete[] pAr1; pAr1 = nullptr;
		delete[] pAr2; pAr2 = nullptr;
		delete[] pAr3; pAr3 = nullptr;
		delete[] pAr4; pAr4 = nullptr;
		delete[] pAr5; pAr5 = nullptr;
		delete[] pAr6; pAr6 = nullptr;
		std::exit(EXIT_FAILURE);
	}
	return std::make_tuple(pAr1, pAr2, pAr3, pAr4, pAr5, pAr6);
}

void   nbody_implementation::NBodyDArray::run_simulation_omp(_In_ const double lo, _In_ const double hi) {

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
	auto ptrs = allocate_array1D(nParticles);
	std::fill(&std::get<3>(ptrs)[0], &std::get<3>(ptrs)[0] + this->m_uinParticles, 0.0);
	std::fill(&std::get<4>(ptrs)[0], &std::get<4>(ptrs)[0] + this->m_uinParticles, 0.0);
	std::fill(&std::get<5>(ptrs)[0], &std::get<5>(ptrs)[0] + this->m_uinParticles, 0.0);
	/* allocating automatic arrays for OpenMP threading.*/
	for (std::size_t i{ 0 }; i != this->m_uinParticles; i += 4) {
		// Pos X
		_mm256_storeu_pd(&std::get<0>(ptrs)[i], _mm256_set_pd(rand(), rand(), rand(), rand()));
		// Pos Y
		_mm256_storeu_pd(&std::get<1>(ptrs)[i], _mm256_set_pd(rand(), rand(), rand(), rand()));
		// Pos Z
		_mm256_storeu_pd(&std::get<2>(ptrs)[i], _mm256_set_pd(rand(), rand(), rand(), rand()));
	}
#if defined VERBOSE
	std::cout << "Done!!" << std::endl;
	std::cout << "OpenMP in use with: " << this->m_inThreads << " threads" << std::endl;
#endif
	::omp_set_num_threads(nThreads);
	const double d_time_start{ ::omp_get_wtime() };
#pragma omp parallel for schedule(runtime)
	for (int i = 0; i < nParticles; i += 4) {
		__m256d Fx = _mm256_setzero_pd();
		__m256d Fy = _mm256_setzero_pd();
		__m256d Fz = _mm256_setzero_pd();
		__m256d tx = _mm256_setzero_pd();
		__m256d ty = _mm256_setzero_pd();
		__m256d tz = _mm256_setzero_pd();

		for (int j = 0; j < nParticles; j += 4) {

			if (j != i) {

				_mm_prefetch(reinterpret_cast<const char*>(&std::get<0>(ptrs)[j + 2]), _MM_HINT_T0);
				__m256d dx(_mm256_sub_pd(_mm256_load_pd(&std::get<0>(ptrs)[j]), _mm256_load_pd(&std::get<0>(ptrs)[i])));

				_mm_prefetch(reinterpret_cast<const char*>(&std::get<1>(ptrs)[j + 2]), _MM_HINT_T0);
				__m256d dy(_mm256_sub_pd(_mm256_load_pd(&std::get<1>(ptrs)[j]), _mm256_load_pd(&std::get<1>(ptrs)[i])));

				_mm_prefetch(reinterpret_cast<const char*>(&std::get<2>(ptrs)[j + 2]), _MM_HINT_T0);
				__m256d dz(_mm256_sub_pd(_mm256_load_pd(&std::get<2>(ptrs)[j + 2]), _mm256_load_pd(&std::get<2>(ptrs)[i])));

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
		_mm256_store_pd(&std::get<3>(ptrs)[i], _mm256_mul_pd(_mm256_hadd_pd(tx, tx), v_dt));
		_mm256_store_pd(&std::get<4>(ptrs)[i], _mm256_mul_pd(_mm256_hadd_pd(ty, ty), v_dt));
		_mm256_store_pd(&std::get<5>(ptrs)[i], _mm256_mul_pd(_mm256_hadd_pd(tz, tz), v_dt));
	}

	for (int i{ 0 }; i != this->m_uinParticles; i += 8) {

		_mm256_store_pd(&std::get<0>(ptrs)[i], _mm256_add_pd(_mm256_load_pd(&std::get<3>(ptrs)[i]), v_dt));
		_mm256_store_pd(&std::get<0>(ptrs)[i + 4], _mm256_add_pd(_mm256_load_pd(&std::get<3>(ptrs)[i + 4]), v_dt));
		_mm256_store_pd(&std::get<1>(ptrs)[i], _mm256_add_pd(_mm256_load_pd(&std::get<4>(ptrs)[i]), v_dt));
		_mm256_store_pd(&std::get<1>(ptrs)[i + 4], _mm256_add_pd(_mm256_load_pd(&std::get<4>(ptrs)[i + 4]), v_dt));
		_mm256_store_pd(&std::get<2>(ptrs)[i], _mm256_add_pd(_mm256_load_pd(&std::get<5>(ptrs)[i]), v_dt));
		_mm256_store_pd(&std::get<2>(ptrs)[i + 4], _mm256_add_pd(_mm256_load_pd(&std::get<5>(ptrs)[i + 4]), v_dt));
	}
	const double d_stop_time{ ::omp_get_wtime() };
#if defined VERBOSE
	std::cout << "Crude timing result: " << std::showpoint << std::fixed << d_stop_time - d_time_start << " sec, " << (d_stop_time - d_time_start) * 1.0E+9 << " ns" << std::endl;
#endif
	// Copy to member arrays.
	std::copy(&std::get<0>(ptrs)[0], &std::get<0>(ptrs)[0] + this->m_uinParticles, &this->m_pX[0]);
	std::copy(&std::get<1>(ptrs)[0], &std::get<1>(ptrs)[0] + this->m_uinParticles, &this->m_pY[0]);
	std::copy(&std::get<2>(ptrs)[0], &std::get<2>(ptrs)[0] + this->m_uinParticles, &this->m_pZ[0]);
	// Deallocate automatic arrays
	delete[] std::get<0>(ptrs);
	delete[] std::get<1>(ptrs);
	delete[] std::get<2>(ptrs);
	delete[] std::get<3>(ptrs);
	delete[] std::get<4>(ptrs);
	delete[] std::get<5>(ptrs);
}

void   nbody_implementation::NBodyDArray::run_simulation_scalar(_In_ const double lo, _In_ const double hi) {


}
	




