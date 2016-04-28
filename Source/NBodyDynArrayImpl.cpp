
#include "NBodyDynArrayImpl.h"

nbody_implementation::NBodyDArray::NBodyDArray(_In_ const std::size_t nParticles, _In_ const int nThreads, _In_ const double nSteps) :
m_uinParticles{ nParticles },
m_inThreads{ nThreads },
m_dStep{ nSteps }
{
	// Exit on any failure deallocate already allocated memory
	// in order to prevent leak.
	// Object state consistent only partially before executing catch block.
	
	try{
		this->m_pX = new double[this->m_uinParticles];
	}
	catch (std::bad_alloc& ba) {
#ifdef DEBUG_VERBOSE
		std::cout << "General memory allocation failure at loc #: " << __LINE__ << std::endl;
		std::cout << " In:" << typeid(NBodyDArray).name() << "::NBodyDArray" << std::endl;
		std::cout << "!!Critical failure calling... exit(1) now!!" << std::endl;
#else
		std::cout << ba.what() << std::endl;
#endif
		delete[] this->m_pX;
		std::exit(EXIT_FAILURE);
	}
	
	try{
		this->m_pY = new double[this->m_uinParticles];
	}
	catch (std::bad_alloc& ba) {
#ifdef DEBUG_VERBOSE
		std::cout << "General memory allocation failure at loc #: " << __LINE__ << std::endl;
		std::cout << " In:" << typeid(NBodyDArray).name() << "::NBodyDArray" << std::endl;
		std::cout << "!!Critical failure calling... exit(1) now!!" << std::endl;
#else
		std::cout << ba.what() << std::endl;
#endif
		delete[] this->m_pX; delete[] this->m_pY;
		std::exit(EXIT_FAILURE);
	}
	
	
	try{
		this->m_pZ = new double[this->m_uinParticles];
	}
	catch (std::bad_alloc& ba) {
#ifdef DEBUG_VERBOSE
		std::cout << "General memory allocation failure at loc #: " << __LINE__ << std::endl;
		std::cout << " In:" << typeid(NBodyDArray).name() << "::NBodyDArray" << std::endl;
		std::cout << "!!Critical failure calling... exit(1) now!!" << std::endl;
#else
		std::cout << ba.what() << std::endl;
#endif
		delete[] this->m_pX; delete[] this->m_pY; delete[] this->m_pZ;
		std::exit(EXIT_FAILURE);
	}
		
	

	try{
		this->m_pVX = new double[this->m_uinParticles];
	}
	catch (std::bad_alloc& ba) {
#ifdef DEBUG_VERBOSE
		std::cout << "General memory allocation failure at loc #: " << __LINE__ << std::endl;
		std::cout << " In:" << typeid(NBodyDArray).name() << "::NBodyDArray" << std::endl;
		std::cout << "!!Critical failure calling... exit(1) now!!" << std::endl;
#else
		std::cout << ba.what() << std::endl;
#endif
		delete[] this->m_pX; delete[] this->m_pY; delete[] this->m_pZ;
		delete[] this->m_pVX;
		std::exit(EXIT_FAILURE);
	}
	try{
		this->m_pVY = new double[this->m_uinParticles];
	}
	catch (std::bad_alloc& ba) {
#ifdef DEBUG_VERBOSE
		std::cout << "General memory allocation failure at loc #: " << __LINE__ << std::endl;
		std::cout << " In:" << typeid(NBodyDArray).name() << "::NBodyDArray" << std::endl;
		std::cout << "!!Critical failure calling... exit(1) now!!" << std::endl;
#else
		std::cout << ba.what() << std::endl;
#endif
		delete[] this->m_pX; delete[] this->m_pY; delete[] this->m_pZ;
		delete[] this->m_pVX; delete[] this->m_pVY;
		std::exit(EXIT_FAILURE);
	}
	try{
		this->m_pVZ = new double[this->m_uinParticles];
	}
	catch (std::bad_alloc& ba) {
#ifdef DEBUG_VERBOSE
		std::cout << "General memory allocation failure at loc #: " << __LINE__ << std::endl;
		std::cout << " In:" << typeid(NBodyDArray).name() << "::NBodyDArray" << std::endl;
		std::cout << "!!Critical failure calling... exit(1) now!!" << std::endl;
#else
		std::cout << ba.what() << std::endl;
#endif
		delete[] this->m_pX; delete[] this->m_pY; delete[] this->m_pZ;
		delete[] this->m_pVX; delete[] this->m_pVY; delete[] this->m_pZ;
		std::exit(EXIT_FAILURE);
	}
#ifdef ICC_HIGHEST_OPTIM_LEVEL 

	std::fill(&this->m_pX[0], &this->m_pX[0] + this->m_uinParticles, 0.0);
	std::fill(&this->m_pY[0], &this->m_pY[0] + this->m_uinParticles, 0.0);
	std::fill(&this->m_pZ[0], &this->m_pZ[0] + this->m_uinParticles, 0.0);
	std::fill(&this->m_pVX[0], &this->m_pVX[0] + this->m_uinParticles, 0.0);
	std::fill(&this->m_pY[0], &this->m_pY[0] + this->m_uinParticles, 0.0);
	std::fill(&this->m_pVZ[0], &this->m_pVZ[0] + this->m_uinParticles, 0.0);
#else
	/* Using vector memory load/store 
	   Assumes that Core i7 Sandy Bridge is present
	   Unrolling by 2 in use.Unaligned store by default */
	for (std::size_t i{ 0 }; i != this->m_uinParticles; i += 8) {
		_mm256_storeu_pd(&this->m_pX[i], _mm256_setzero_pd());
		_mm256_storeu_pd(&this->m_pX[i + 4], _mm256_setzero_pd());
		_mm256_storeu_pd(&this->m_pY[i], _mm256_setzero_pd());
		_mm256_storeu_pd(&this->m_pVY[i + 4], _mm256_setzero_pd());
		_mm256_storeu_pd(&this->m_pZ[i], _mm256_setzero_pd());
		_mm256_storeu_pd(&this->m_pZ[i + 4], _mm256_setzero_pd());
		_mm256_storeu_pd(&this->m_pVX[i], _mm256_setzero_pd());
		_mm256_storeu_pd(&this->m_pVX[i + 4], _mm256_setzero_pd());
		_mm256_storeu_pd(&this->m_pVY[i], _mm256_setzero_pd());
		_mm256_storeu_pd(&this->m_pVY[i + 4], _mm256_setzero_pd());
		_mm256_storeu_pd(&this->m_pVZ[i], _mm256_setzero_pd());
		_mm256_storeu_pd(&this->m_pVZ[i + 4], _mm256_setzero_pd());
	}
#endif

#ifdef VERBOSE
	std::cout << "Executing in NBodyDArray::NBodyDArray(std::size_t,int,double)" << std::endl;
#endif
}

nbody_implementation::NBodyDArray::NBodyDArray(_In_ const NBodyDArray &rhs) :
m_uinParticles{ rhs.m_uinParticles },
m_inThreads{ rhs.m_inThreads },
m_dStep{ rhs.m_dStep }

{
	try{
		this->m_pX = new double[this->m_uinParticles];
	}
	catch (std::bad_alloc& ba) {
#ifdef DEBUG_VERBOSE
		std::cout << "General memory allocation failure at loc #: " << __LINE__ << std::endl;
		std::cout << " In:" << typeid(NBodyDArray).name() << "::NBodyDArray" << std::endl;
		std::cout << "!!Critical failure calling... exit(1) now!!" << std::endl;
#else
		std::cout << ba.what() << std::endl;
#endif
		delete[] this->m_pX;
		std::exit(EXIT_FAILURE);
	}

	try{
		this->m_pY = new double[this->m_uinParticles];
	}
	catch (std::bad_alloc& ba) {
#ifdef DEBUG_VERBOSE
		std::cout << "General memory allocation failure at loc #: " << __LINE__ << std::endl;
		std::cout << " In:" << typeid(NBodyDArray).name() << "::NBodyDArray" << std::endl;
		std::cout << "!!Critical failure calling... exit(1) now!!" << std::endl;
#else
		std::cout << ba.what() << std::endl;
#endif
		delete[] this->m_pX; delete[] this->m_pY;
		std::exit(EXIT_FAILURE);
	}


	try{
		this->m_pZ = new double[this->m_uinParticles];
	}
	catch (std::bad_alloc& ba) {
#ifdef DEBUG_VERBOSE
		std::cout << "General memory allocation failure at loc #: " << __LINE__ << std::endl;
		std::cout << " In:" << typeid(NBodyDArray).name() << "::NBodyDArray" << std::endl;
		std::cout << "!!Critical failure calling... exit(1) now!!" << std::endl;
#else
		std::cout << ba.what() << std::endl;
#endif
		delete[] this->m_pX; delete[] this->m_pY; delete[] this->m_pZ;
		std::exit(EXIT_FAILURE);
	}



	try{
		this->m_pVX = new double[this->m_uinParticles];
	}
	catch (std::bad_alloc& ba) {
#ifdef DEBUG_VERBOSE
		std::cout << "General memory allocation failure at loc #: " << __LINE__ << std::endl;
		std::cout << " In:" << typeid(NBodyDArray).name() << "::NBodyDArray" << std::endl;
		std::cout << "!!Critical failure calling... exit(1) now!!" << std::endl;
#else
		std::cout << ba.what() << std::endl;
#endif
		delete[] this->m_pX; delete[] this->m_pY; delete[] this->m_pZ;
		delete[] this->m_pVX;
		std::exit(EXIT_FAILURE);
	}
	try{
		this->m_pVY = new double[this->m_uinParticles];
	}
	catch (std::bad_alloc& ba) {
#ifdef DEBUG_VERBOSE
		std::cout << "General memory allocation failure at loc #: " << __LINE__ << std::endl;
		std::cout << " In:" << typeid(NBodyDArray).name() << "::NBodyDArray" << std::endl;
		std::cout << "!!Critical failure calling... exit(1) now!!" << std::endl;
#else
		std::cout << ba.what() << std::endl;
#endif
		delete[] this->m_pX; delete[] this->m_pY; delete[] this->m_pZ;
		delete[] this->m_pVX; delete[] this->m_pVY;
		std::exit(EXIT_FAILURE);
	}
	try{
		this->m_pVZ = new double[this->m_uinParticles];
	}
	catch (std::bad_alloc& ba) {
#ifdef DEBUG_VERBOSE
		std::cout << "General memory allocation failure at loc #: " << __LINE__ << std::endl;
		std::cout << " In:" << typeid(NBodyDArray).name() << "::NBodyDArray" << std::endl;
		std::cout << "!!Critical failure calling... exit(1) now!!" << std::endl;
#else
		std::cout << ba.what() << std::endl;
#endif
		delete[] this->m_pX; delete[] this->m_pY; delete[] this->m_pZ;
		delete[] this->m_pVX; delete[] this->m_pVY; delete[] this->m_pZ;
		std::exit(EXIT_FAILURE);
	}
	/* At highest optimization level ICC will inline and vectorise call to
	   std::copy so undef it when high or highest level is set.*/
#ifdef ICC_HIGHEST_OPTIM_LEVEL
	std::copy(&rhs.m_pX[0], &rhs.m_pX[0] + rhs.m_uinParticles, &this->m_pX[0]);
	std::copy(&rhs.m_pY[0], &rhs.m_pY[0] + rhs.m_uinParticles, &this->m_pY[0]);
	std::copy(&rhs.m_pZ[0], &rhs.m_pZ[0] + rhs.m_uinParticles, &this->m_pZ[0]);
	std::copy(&rhs.m_pVX[0], &rhs.m_pVX[0] + rhs.m_uinParticles, &this->m_pVX[0]);
	std::copy(&rhs.m_pVY[0], &rhs.m_pVY[0] + rhs.m_uinParticles, &this->m_pVY[0]);
	std::copy(&rhs.m_pVZ[0], &rhs.m_pVZ[0] + rhs.m_uinParticles, &this->m_pVZ[0]);
#else
	for (std::size_t i{ 0 }; i != this->m_uinParticles; i += 8) {
		_mm256_storeu_pd(&this->m_pX[i], _mm256_loadu_pd(&rhs.m_pX[i]));
		_mm256_storeu_pd(&this->m_pX[i + 4], _mm256_loadu_pd(&rhs.m_pX[i + 4]));
		_mm256_storeu_pd(&this->m_pY[i], _mm256_loadu_pd(&rhs.m_pY[i]));
		_mm256_storeu_pd(&this->m_pY[i + 4], _mm256_loadu_pd(&rhs.m_pY[i + 4]));
		_mm256_storeu_pd(&this->m_pZ[i], _mm256_loadu_pd(&rhs.m_pZ[i]));
		_mm256_storeu_pd(&this->m_pZ[i + 4], _mm256_loadu_pd(&rhs.m_pZ[i + 4]));
		_mm256_storeu_pd(&this->m_pVX[i], _mm256_loadu_pd(&rhs.m_pVX[i]));
		_mm256_storeu_pd(&this->m_pVX[i + 4], _mm256_loadu_pd(&rhs.m_pVX[i + 4]));
		_mm256_storeu_pd(&this->m_pY[i], _mm256_loadu_pd(&rhs.m_pVY[i]));
		_mm256_storeu_pd(&this->m_pY[i + 4], _mm256_loadu_pd(&rhs.m_pY[i + 4]));
		_mm256_storeu_pd(&this->m_pVZ[i], _mm256_loadu_pd(&rhs.m_pVZ[i]));
		_mm256_storeu_pd(&this->m_pVZ[i + 4], _mm256_loadu_pd(&rhs.m_pZ[i + 4]));
	}
#endif
#ifndef DEBUG_VERBOSE
	std::cout << "Executing in NBodyDArray::NBodyDArray(NBodyDArray& )" << std::endl;
#endif
}

nbody_implementation::NBodyDArray::NBodyDArray(_In_ NBodyDArray &&rhs) :
m_uinParticles{},
m_inThreads{},
m_dStep{},
m_pX( nullptr ),
m_pY( nullptr ),
m_pZ( nullptr ),
m_pVX( nullptr ),
m_pVY( nullptr ),
m_pVZ( nullptr )
{

	this->m_pX = rhs.m_pX;
	this->m_pY = rhs.m_pY;
	this->m_pZ = rhs.m_pZ;
	this->m_pVX = rhs.m_pVX;
	this->m_pVY = rhs.m_pVY;
	this->m_pVZ = rhs.m_pVZ;
	this->m_uinParticles = rhs.m_uinParticles;
	this->m_inThreads = rhs.m_inThreads;
	this->m_dStep = rhs.m_dStep;


	// Nullify the rhs state.
	rhs.m_pX = nullptr;
	rhs.m_pY = nullptr;
	rhs.m_pZ = nullptr;
	rhs.m_pVX = nullptr;
	rhs.m_pVY = nullptr;
	rhs.m_pVZ = nullptr;
	rhs.m_dStep = 0.0;
	rhs.m_inThreads = 0;
	rhs.m_uinParticles = 0;
#ifndef DEBUG_VERBOSE 
	std::cout << "Executing in NBodyDArray::NBodyDArray(NBodyDArray&& )" << std::endl;
#endif

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
		throw std::invalid_argument(std::string("Invalid input in: NBodyDArray::run_simualtion\n"));

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



std::tuple<double*, double*, double*, double*, double*, double*> nbody_implementation::NBodyDArray::allocate_array1D(_In_ const std::size_t nParticles) {
	if (nParticles <= 0)
		throw std::runtime_error(std::string{ "Invalid input in: " }.append(typeid(NBodyDArray).name()));
	
	double* pAr1, *pAr2, *pAr3, *pAr4, *pAr5, *pAr6;
	
	
	
	try{
		pAr1 = new double[nParticles];
	}
	catch (std::bad_alloc& ba) {
#ifdef DEBUG_VERBOSE
		std::cerr << "Memory allocation failure in: " << typeid(NBodyDArray).name() << "::allocate_array1D, calling exit..." << std::endl;
		std::cerr << "Error occurred at line #:" << __LINE__ << std::endl;
#else
		std::cerr << ba.what() << "at line #: " << __LINE__ << std::endl;
#endif
		delete[] pAr1; pAr1 = nullptr;
		::exit(EXIT_FAILURE);
	}
	try {
		pAr2 = new double[nParticles];
	}
	catch (std::bad_alloc& ba){
#ifdef DEBUG_VERBOSE
		std::cerr << "Memory allocation failure in: " << typeid(NBodyDArray).name() << "::allocate_array1D, calling exit..." << std::endl;
		std::cerr << "Error occurred at line #:" << __LINE__ << std::endl;
#else
		std::cerr << ba.what() << "at line #: " << __LINE__ << std::endl;
#endif
		delete[] pAr1; pAr1 = nullptr;
		delete[] pAr2; pAr2 = nullptr;
		::exit(EXIT_FAILURE);
	}
	try{
		pAr3 = new double[nParticles];
	}
	catch (std::bad_alloc& ba) {
#ifdef DEBUG_VERBOSE
		std::cerr << "Memory allocation failure in: " << typeid(NBodyDArray).name() << "::allocate_array1D, calling exit..." << std::endl;
		std::cerr << "Error occurred at line #:" << __LINE__ << std::endl;
#else
		std::cerr << ba.what() << "at line #: " << __LINE__ << std::endl;
#endif
		delete[] pAr1; pAr1 = nullptr;
		delete[] pAr2; pAr2 = nullptr;
		delete[] pAr3; pAr3 = nullptr;
		::exit(EXIT_FAILURE);
	}
	try{
		pAr4 = new double[nParticles];
	}
	catch (std::bad_alloc& ba) {
#ifdef DEBUG_VERBOSE
		std::cerr << "Memory allocation failure in: " << typeid(NBodyDArray).name() << "::allocate_array1D, calling exit..." << std::endl;
		std::cerr << "Error occurred at line #:" << __LINE__ << std::endl;
#else
		std::cerr << ba.what() << "at line #: " << __LINE__ << std::endl;
#endif
		delete[] pAr1; pAr1 = nullptr;
		delete[] pAr2; pAr2 = nullptr;
		delete[] pAr3; pAr3 = nullptr;
		delete[] pAr4; pAr4 = nullptr;
	}
	try{
		pAr5 = new double[nParticles];
	}
	catch (std::bad_alloc& ba) {
#ifdef DEBUG_VERBOSE
		std::cerr << "Memory allocation failure in: " << typeid(NBodyDArray).name() << "::allocate_array1D, calling exit..." << std::endl;
		std::cerr << "Error occurred at line #:" << __LINE__ << std::endl;
#else
		std::cerr << ba.what() << "at line #: " << __LINE__ << std::endl;
#endif
		delete[] pAr1; pAr1 = nullptr;
		delete[] pAr2; pAr2 = nullptr;
		delete[] pAr3; pAr3 = nullptr;
		delete[] pAr4; pAr4 = nullptr;
		delete[] pAr5; pAr5 = nullptr;
	}
	try{
		pAr6 = new double[nParticles];
	}
	catch (std::bad_alloc& ba) {
#ifdef DEBUG_VERBOSE
		std::cerr << "Memory allocation failure in: " << typeid(NBodyDArray).name() << "::allocate_array1D, calling exit..." << std::endl;
		std::cerr << "Error occurred at line #:" << __LINE__ << std::endl;
#else
		std::cerr << ba.what() << "at line #: " << __LINE__ << std::endl;
#endif
		delete[] pAr1; pAr1 = nullptr;
		delete[] pAr2; pAr2 = nullptr;
		delete[] pAr3; pAr3 = nullptr;
		delete[] pAr4; pAr4 = nullptr;
		delete[] pAr5; pAr5 = nullptr;
		delete[] pAr6; pAr6 = nullptr;
	}
	return (std::make_tuple(pAr1, pAr2, pAr3, pAr4, pAr5, pAr6));
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

				_mm_prefetch(reinterpret_cast<const char*>(&std::get<0>(ptrs)[j + 4]), _MM_HINT_T0);
				__m256d dx(_mm256_sub_pd(_mm256_load_pd(&std::get<0>(ptrs)[j]), _mm256_load_pd(&std::get<0>(ptrs)[i])));

				_mm_prefetch(reinterpret_cast<const char*>(&std::get<1>(ptrs)[j + 4]), _MM_HINT_T0);
				__m256d dy(_mm256_sub_pd(_mm256_load_pd(&std::get<1>(ptrs)[j]), _mm256_load_pd(&std::get<1>(ptrs)[i])));

				_mm_prefetch(reinterpret_cast<const char*>(&std::get<2>(ptrs)[j + 4]), _MM_HINT_T0);
				__m256d dz(_mm256_sub_pd(_mm256_load_pd(&std::get<2>(ptrs)[j]), _mm256_load_pd(&std::get<2>(ptrs)[i])));

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

		this->m_pX[i] = rand();
		this->m_pY[i] = rand();
		this->m_pZ[i] = rand();
	}
	std::cout << "Finished initialization of " << this->m_uinParticles << "Particles with:" << typeid(rand).name() << std::endl;
#ifdef EXECUTE_NTIMES
	std::cout << "Starting executing NBody system: " << NTIMES << std::endl;
	for (int k{ 0 }; k != NTIMES; ++k) {
#endif
		for (std::size_t i{ 0 }; i != this->m_uinParticles; ++i) {
			double Fx{ 0.0 }; double Fy{ 0.0 }; double Fz{ 0.0 };

			for (std::size_t j{ 0 }; j != this->m_uinParticles; ++j) {

				if (j != i) {

					_mm_prefetch(reinterpret_cast<const char*>(&m_pX[j + 7]), _MM_HINT_T0);
					const double dx{ this->m_pX[j] - this->m_pX[i] };
					_mm_prefetch(reinterpret_cast<const char*>(&m_pY[j + 7]), _MM_HINT_T0);
					const double dy{ this->m_pY[j] - this->m_pY[i] };
					_mm_prefetch(reinterpret_cast<const char*>(&m_pZ[j + 7]), _MM_HINT_T0);
					const double dz{ this->m_pZ[j] - this->m_pZ[i] };

					const double distSquared{ dx * dx + dy * dy + dz * dz };
					const double invdistCube{ 1.0 / (distSquared * ::sqrt(distSquared)) };

					Fx += dx * invdistCube;
					Fy += dy * invdistCube;
					Fz += dz * invdistCube;

				}
			}
			this->m_pVX[i] = Fx * this->m_dStep;
			this->m_pVY[i] = Fy * this->m_dStep;
			this->m_pVZ[i] = Fz * this->m_dStep;

		}
		for (std::size_t i{ 0 }; i != this->m_uinParticles; ++i) {
			this->m_pX[i] = this->m_pVX[i] + this->m_dStep;
			this->m_pY[i] = this->m_pVY[i] + this->m_dStep;
			this->m_pZ[i] = this->m_pVZ[i] + this->m_dStep;
		}
#ifdef EXECUTE_NTIMES
		std::cout << "Finished executing NBody system: " << NTIMES << std::endl;
	}
#endif

}
#ifdef DEBUG_VERBOSE
void    nbody_implementation::NBodyDArray::init_mem_dbg(_In_ const double val) {

	std::fill(&this->m_pX[0], &this->m_pX[0] + this->m_uinParticles, val);
	std::fill(&this->m_pY[0], &this->m_pY[0] + this->m_uinParticles, val);
	std::fill(&this->m_pZ[0], &this->m_pZ[0] + this->m_uinParticles, val);
	std::fill(&this->m_pVX[0], &this->m_pVX[0] + this->m_uinParticles, val);
	std::fill(&this->m_pY[0], &this->m_pY[0] + this->m_uinParticles, val);
	std::fill(&this->m_pVZ[0], &this->m_pVZ[0] + this->m_uinParticles, val);
}
#endif	




