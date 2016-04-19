
#include "ArraysTest.h"
#include <array>
auto  test::ArraysMT::allocateVec1D(_In_ const int datumSize)->std::valarray < double > {

	_ASSERT((datumSize != 0) && (datumSize % 32 != 0));
	std::valarray<double> d_Data = std::valarray<double>(datumSize);
	return d_Data;
}

auto  test::ArraysMT::initVec1D(_In_ std::valarray<double> & Data)->void {
	
	clock_t seed{ ::clock() };
	auto rand_gen = std::bind(std::uniform_real_distribution < double > {},
		std::default_random_engine(seed));
	int i{ 0 };
	omp_set_num_threads(8);
	long long start{ _rdtsc() };
#pragma omp parallel for private(i) default(shared)
	
	for (i; i < Data.size(); i += 8) {
		// Vectorize only
		_mm256_storeu_pd(&Data.operator[](i), _mm256_set_pd(rand_gen(), rand_gen(), rand_gen(), rand_gen()));
		_mm256_storeu_pd(&Data.operator[](i + 4), _mm256_set_pd(rand_gen(), rand_gen(), rand_gen(), rand_gen()));
	}
	 long long stop{ _rdtsc() };
	 if ((stop - start) > 0)
	   std::printf("initVec1D executed in:%llu TSC_CYCLES\n", stop - start);
	 else
	   std::printf("Fatal Error: Negative timing\n");
}

auto    test::ArraysMT::runComputation(_In_ std::valarray<double> &Data, _Out_ std::valarray<double> &Out)->void {
	if (Data.size() != Out.size())
		throw std::runtime_error("Invalid Array size in runComputation\n");
	constexpr int nThreads{ 4 };
	
	constexpr double L{ 2.5 }; /* length in meters*/
	constexpr double invL{ 0.4 };
	constexpr double r{ 0.1 }; /* radius of cylindrical coordinates */
	const int K = Data.size();
	__m256d r_over_L(_mm256_div_pd(_mm256_set1_pd(r), _mm256_set1_pd(L)));
	__m256d r_over_Ls(_mm256_mul_pd(r_over_L, r_over_L));
	omp_set_num_threads(nThreads);
	int i;
#pragma omp parallel for private(i) default(shared)
	for (i = 0; i < K; i += 8) {
		__m256d Qk = _mm256_sub_pd(_mm256_mul_pd(_mm256_loadu_pd(&Data[i]), _mm256_set1_pd(invL)), _mm256_mul_pd(_mm256_loadu_pd(&Data[i + 4]), _mm256_set1_pd(invL)));
		__m256d Qk2 = _mm256_mul_pd(Qk, Qk);
		__m256d Denom = _mm256_add_pd(Qk2, r_over_Ls);
		__m256d SqrtDenom = _mm256_pow_pd(Denom, _mm256_set1_pd(1.5));
		_mm256_storeu_pd(&Out[i], _mm256_div_pd(Qk, SqrtDenom));
		_mm256_storeu_pd(&Out[i + 4], _mm256_div_pd(Qk, SqrtDenom));
	}
	
}

auto    test::ArraysMT::printData(_In_ std::valarray<double> & Datum)->void {

	for (int i{ 0 }; i != Datum.size(); ++i)std::printf("Datum[%d]=%.16f\n",i, Datum[i]);
}

auto    test::ArraysMT::testArray2D()->void {
	constexpr int Rows{ 1024 };
	constexpr int Cols{ 1024 };
	clock_t seed{ ::clock() };
	auto rand_gen = std::bind(std::uniform_real_distribution < double > {},
		std::default_random_engine(seed));
	std::array<std::array<double, Rows>, Cols> dData{};
	
	long long start{ _rdtsc() };
	omp_set_num_threads(8);
	int i,j;
#pragma omp parallel for private(i) default(shared)
	for ( i = 0 ; i != Rows; ++i)
		for (int j{ 0 }; j != Cols; j += 8){
		_mm256_storeu_pd(&dData[i][j], _mm256_set_pd(rand_gen(), rand_gen(), rand_gen(), rand_gen()));
		_mm256_storeu_pd(&dData[i][j + 4], _mm256_set_pd(rand_gen(), rand_gen(), rand_gen(), rand_gen()));
		
		}
	long long stop{ _rdtsc() };
	if ((stop - start) > 0)
		std::printf("testArray2D executed in:%llu TSC_CYCLES\n", stop - start);
	else
		std::printf("Fatal Error: Negative timing\n");
}

auto      test::ArraysMT::simple_test()-> void{
	constexpr int arLength{ 1024 };
	__declspec(align(32))double ar1[arLength] = { 0.0 };
	std::array<double, arLength> ar2{};
	std::valarray<double> ar3(arLength);
	std::vector<double> ar4(arLength);
	std::vector<double>::iterator it;
	for (int i{ 0 }; i != arLength; ++i) ar1[i] = static_cast<double>(i);
	
	for (int i{ 0 }; i != arLength; ++i) ar2[i] = static_cast<double>(i);
	for (int i{ 0 }; i != arLength; ++i) ar3[i] = static_cast<double>(i);

	int i = 1;
	for (it = ar4.begin(); it != ar4.end(); ++it)
		*it = static_cast<double>(++i);

	}
		
//   *******************    Particle Simulation for different data types   ***************************  //	


/*
@brief         Built-in static arrays.
*/
	
   auto    test::ArraysMT::PartSimArray( _In_ const int nSteps, _In_ const int nThreads,
	_In_ const double dt)->void {

	//static_assert(N != 0, " Number of particles cannot be zero");
	   if ((nSteps <= 0) || (nThreads <= 0) || (dt == 0.0))
		   throw std::runtime_error("Wrong arguments in PartSimArray!!\n");

		std::cout << "Built-in arrays allocation started" << std::endl;
	constexpr int nParticles{ 1 << 13 };
	__declspec(align(32))double X[nParticles] = { 0.0 };
	__declspec(align(32))double Y[nParticles] = { 0.0 };
	__declspec(align(32))double Z[nParticles] = { 0.0 };
	__declspec(align(32))double VX[nParticles] = { 0.0 };
	__declspec(align(32))double VY[nParticles] = { 0.0 };
	__declspec(align(32))double VZ[nParticles] = { 0.0 };
	 
	 int i, ii, j;
	std::cout << "Allocated 6 static arrays of: " << 6 * nParticles * sizeof(double) << " bytes" << std::endl;
	std::cout << "Begin initialization of Particles position and velocity arrays\n";
	clock_t seed{ clock() };
	auto rand = std::bind(std::uniform_real_distribution < double > {-1.0, 1.0}, std::default_random_engine(seed));
	//omp_set_num_threads(nThreads);
//#pragma omp parallel for private(i) default(shared)
	for (i = 0; i != nParticles; i += 8) {
		_mm256_store_pd(&X[i], _mm256_set1_pd(rand()));
		_mm256_store_pd(&X[i + 4], _mm256_set1_pd(rand()));
		_mm256_store_pd(&Y[i], _mm256_set1_pd(rand()));
		_mm256_store_pd(&Y[i + 4], _mm256_set1_pd(rand()));
		_mm256_store_pd(&Z[i], _mm256_set1_pd(rand()));
		_mm256_store_pd(&Z[i + 4], _mm256_set1_pd(rand()));
		_mm256_store_pd(&VX[i], _mm256_set1_pd(rand()));
		_mm256_store_pd(&VX[i + 4], _mm256_set1_pd(rand()));
		_mm256_store_pd(&VY[i], _mm256_set1_pd(rand()));
		_mm256_store_pd(&VY[i + 4], _mm256_set1_pd(rand()));
		_mm256_store_pd(&VZ[i], _mm256_set1_pd(rand()));
		_mm256_store_pd(&VZ[i + 4], _mm256_set1_pd(rand()));
	}
	std::cout << "Finished initialization of " << nParticles << "Particles with:" << typeid(rand).name() << std::endl;
	double rate{ 0.0 };
	double d_rate{ 0.0 };
	const int skipSteps{ 1 };
#if 0
	int nIters = 0;
#endif
	
	std::cout << "Starting:" << nSteps << "iterations outer loop\n";
	omp_set_num_threads(nThreads);
	
	for (int i{ 1 }; i != nSteps; ++i) {
		
		const double start{ omp_get_wtime() };

//#pragma omp parallel for private(ii,j) default(none) shared(X,Y,Z,VX,VY,VZ)
#pragma  omp parallel for schedule(dynamic)	
		for (ii = 0;  ii < nParticles;  ++ii){
			double Fx{ 0.0 }; double Fy{ 0.0 }; double Fz{ 0.0 };
			//std::cout << "ii loop\n";
			for (j = 0;   j < nParticles;   ++j) {
				//std::cout << "j loop\n";
				if (j != ii){
#if 0
					++nIters;
#endif
					const double dx{ X[j] - X[ii] };
					const double dy{ Y[j] - Y[ii] };
					const double dz{ Z[j] - Z[ii] };
					const double rdSquared{ dx * dx + dy * dy + dz * dz };
					const double dCbrtSquared{1.0 / (rdSquared * sqrt(rdSquared)) };

					Fx += dx * dCbrtSquared;
					Fy += dy * dCbrtSquared;
					Fz += dz * dCbrtSquared;
				}
			}
			VX[ii] = Fx * dt;
			VY[ii] = Fy * dt;
			VZ[ii] = Fz * dt;
		}
		for (int i{ 0 }; i != nParticles; i += 8){
			//std::cout << "j loop\n";
			_mm256_store_pd(&X[i], _mm256_add_pd(_mm256_load_pd(&VX[i]), _mm256_set1_pd(dt)));
			_mm256_store_pd(&X[i + 4], _mm256_add_pd(_mm256_load_pd(&VX[i + 4]), _mm256_set1_pd(dt)));
			_mm256_store_pd(&Y[i], _mm256_add_pd(_mm256_load_pd(&VY[i]), _mm256_set1_pd(dt)));
			_mm256_store_pd(&Y[i + 4], _mm256_add_pd(_mm256_load_pd(&VY[i + 4]), _mm256_set1_pd(dt)));
			_mm256_store_pd(&Z[i], _mm256_add_pd(_mm256_load_pd(&VZ[i]), _mm256_set1_pd(dt)));
			_mm256_store_pd(&Z[i + 4], _mm256_add_pd(_mm256_load_pd(&VZ[i + 4]), _mm256_set1_pd(dt)));
		}
		const double end{ omp_get_wtime() };
		if (nSteps > skipSteps){
			rate += 1.0 / (end - start);
			d_rate += 1.0 / ((end - start) * (end - start));
		}
		std::cout << std::showpoint << std::scientific << "Time per step= " << end - start << std::endl;
		std::cout << std::showpoint << std::scientific << "Step no:" << nSteps << "Time Delta=" << (end - start) << std::endl;


	}
	rate /= static_cast<double>(nSteps - skipSteps);
	d_rate = sqrt(d_rate / static_cast<double>(nSteps - skipSteps) - rate * rate);
	std::cout << std::showpoint << std::scientific << "Rate= " << rate << "Delta Rate= " << d_rate << std::endl;
#if 0
	std::cout << "nIters = " << nIters << std::endl;
#endif
}

/*
              std::valarray<double> Simple based Particle SSimulation.
*/
auto    test::ArraysMT::PartSimVArray(_In_ const std::size_t nParticles, _In_ const int nSteps, _In_ const int nThreads,
	_In_ const double dt)->void {

	if ((nParticles <= 512) || (nSteps <= 1) || (nThreads <= 2) || (dt == 0.0))
		throw std::runtime_error("Invalid input argument in PartSimVArray\n");

	std::cout << "Construction of 6 " << typeid(std::valarray<double>).name() << "started...";
	std::valarray<double> X(nParticles);
	std::valarray<double> Y(nParticles);
	std::valarray<double> Z(nParticles);
	std::valarray<double> VX(nParticles);
	std::valarray<double> VY(nParticles);
	std::valarray<double> VZ(nParticles);
	std::cout << "Finished construction\n";
	std::cout << "Allocated: " << 6 * nParticles * sizeof(double) << "bytes of heap memory\n";
	clock_t seed{ clock() };
	auto rand = std::bind(std::uniform_real_distribution < double > {-1.0, 1.0}, 
		std::default_random_engine(seed));
	
	//   Unroll 2x times. 
	//   2 load ports/core.
	std::cout << "valarrays initialization started\n";
	std::cout << "2x unrolling in use\n";
	for (std::size_t i{ 0 }; i != X.size(); i += 8) {
		/*
		    Assume valarrays are unaligned on 32-byte boundary.
		*/
		_mm256_storeu_pd(&X[i], _mm256_set_pd(rand(), rand(), rand(), rand()));
		_mm256_storeu_pd(&X[i + 4], _mm256_set_pd(rand(), rand(), rand(), rand()));
		_mm256_storeu_pd(&Y[i], _mm256_set_pd(rand(), rand(), rand(), rand()));
		_mm256_storeu_pd(&Y[i + 4], _mm256_set_pd(rand(), rand(), rand(), rand()));
		_mm256_storeu_pd(&Z[i], _mm256_set_pd(rand(), rand(), rand(), rand()));
		_mm256_storeu_pd(&VX[i], _mm256_set_pd(rand(), rand(), rand(), rand()));
		_mm256_storeu_pd(&VX[i + 4], _mm256_set_pd(rand(), rand(), rand(), rand()));
		_mm256_storeu_pd(&VY[i], _mm256_set_pd(rand(), rand(), rand(), rand()));
		_mm256_storeu_pd(&VY[i + 4], _mm256_set_pd(rand(), rand(), rand(), rand()));
		_mm256_storeu_pd(&VZ[i], _mm256_set_pd(rand(), rand(), rand(), rand()));
		_mm256_storeu_pd(&VZ[i + 4], _mm256_set_pd(rand(), rand(), rand(), rand()));

	}
	std::cout << "Initialization completed\n";
	std::cout << "Made " << 48 * nParticles << "calls to initialize data sets\n";
	std::cout << "Data sets initialized with:" << typeid(rand).name() << std::endl;
	double rate{ 0.0 };
	double d_rate{ 0.0 };
	const int skipSteps{ 1 };
#if 0
	int nIters = 0;
#endif
	std::cout << "Starting " << nSteps << "Simple particle simulation" << std::endl;
	std::cout << "Parallelizing outer loop with " << nThreads << "OpenMP threads\n";
	int ii, j;
	omp_set_num_threads(nThreads);
	
	for (int i{ 1 }; i != nSteps; ++i){
		const double start{ omp_get_wtime() };

#pragma omp parallel for schedule(dynamic)
		for (ii = 0; ii < nParticles; ++ii){
			double Fx{ 0.0 }; double Fy{ 0.0 }; double Fz{ 0.0 };
			//std::cout << "ii loop\n";
			for (j = 0; j < nParticles; ++j) {
				//std::cout << "j loop\n";
				if (j != ii){
#if 0
					++nIters;
#endif
					const double dx{ X[j] - X[ii] };
					const double dy{ Y[j] - Y[ii] };
					const double dz{ Z[j] - Z[ii] };
					const double rdSquared{ dx * dx + dy * dy + dz * dz };
					/*
					     1.0 / (rdSquared * sqrt(rdSquared) created large hotspot
						 which consumed approximately ~45% of CPU time.
					*/
					const double dCbrtSquared{ 1.0 / (rdSquared * sqrt(rdSquared)) };

					Fx += dx * dCbrtSquared;
					Fy += dy * dCbrtSquared;
					Fz += dz * dCbrtSquared;
				}
			}
			VX[ii] = Fx * dt;
			VY[ii] = Fy * dt;
			VZ[ii] = Fz * dt;
		}
		for (int i{ 0 }; i != nParticles; i += 8){
			//std::cout << "j loop\n";
			_mm256_store_pd(&X[i], _mm256_add_pd(_mm256_load_pd(&VX[i]), _mm256_set1_pd(dt)));
			_mm256_store_pd(&X[i + 4], _mm256_add_pd(_mm256_load_pd(&VX[i + 4]), _mm256_set1_pd(dt)));
			_mm256_store_pd(&Y[i], _mm256_add_pd(_mm256_load_pd(&VY[i]), _mm256_set1_pd(dt)));
			_mm256_store_pd(&Y[i + 4], _mm256_add_pd(_mm256_load_pd(&VY[i + 4]), _mm256_set1_pd(dt)));
			_mm256_store_pd(&Z[i], _mm256_add_pd(_mm256_load_pd(&VZ[i]), _mm256_set1_pd(dt)));
			_mm256_store_pd(&Z[i + 4], _mm256_add_pd(_mm256_load_pd(&VZ[i + 4]), _mm256_set1_pd(dt)));
		}
		const double end{ omp_get_wtime() };
		if (nSteps > skipSteps){
			rate += 1.0 / (end - start);
			d_rate += 1.0 / ((end - start) * (end - start));
		}
		std::cout << std::setw(12) << "Step no:" << i << std::endl;
		std::cout << std::setprecision(10) << std::showpoint << std::fixed << "Time per step= " << end - start << std::endl;
		std::cout << std::setprecision(10) << std::showpoint << std::fixed  << "Time Delta=" << (end - start) << std::endl;


	}
	rate /= static_cast<double>(nSteps - skipSteps);
	d_rate = sqrt(d_rate / static_cast<double>(nSteps - skipSteps) - rate * rate);
	std::cout << std::setprecision(10) << std::showpoint << std::fixed << "Rate= " << rate << "Delta Rate= " << d_rate << std::endl;
#if 0
	std::cout << "nIters = " << nIters << std::endl;
#endif
	}

auto   test::ArraysMT::PartSimStdArray(_In_ const int nSteps, _In_ const int nThreads, _In_ const double dt)->void {

	if (nSteps <= 1 || nThreads <= 2 || dt == 0.0)
		throw std::runtime_error("Invalid input argument in PartSimStdArray\n");

	// Number of Particles.
	constexpr int datumSize{ 1 << 13 };
	std::cout << "Construction of 6 " << typeid(std::array<double, 1>).name() << "started...";
	__declspec(align(32))std::array<double, datumSize> X{};
	__declspec(align(32))std::array<double, datumSize> Y{};
	__declspec(align(32))std::array<double, datumSize> Z{};
	__declspec(align(32))std::array<double, datumSize> VX{};
	__declspec(align(32))std::array<double, datumSize> VY{};
	__declspec(align(32))std::array<double, datumSize> VZ{};
	std::cout << "Construction finished" << std::endl;
	std::cout << "Allocated: " << 6 * datumSize * sizeof(double) << "bytes of stack memory" << std::endl;
	clock_t seed{ clock() };
	auto rand = std::bind(std::uniform_real_distribution < double > {0.01, 1.0}, std::default_random_engine(seed));
	const std::string RNG_NAME{ typeid(rand).name() };
	std::cout << "Created RNG of type: " << RNG_NAME.c_str() <<
		"with Interval of randomness: [" << std::fixed << -1.0 << "--" << +1.0 << "]" << std::endl;
	std::cout << "std::array initialization started\n";

	
	for (int i{ 0 }; i != datumSize; i += 4) {
		_mm256_store_pd(&X[i], _mm256_set_pd(rand(), rand(), rand(), rand()));

		_mm256_store_pd(&Y[i], _mm256_set_pd(rand(), rand(), rand(), rand()));

		_mm256_store_pd(&Z[i], _mm256_set_pd(rand(), rand(), rand(), rand()));
	}
	
	
	std::cout << "Initialization completed" << std::endl;
	std::cout << "Made " << 12 * datumSize << "calls to " << RNG_NAME.c_str() << std::endl;
	double rate{ 0.0 };
	double d_rate{ 0.0 };
	const int skipSteps{ 1 };
	std::cout << "Starting " << nSteps << "Simple particle simulation" << std::endl;
	std::cout << "Parallelizing outer loop with " << nThreads << " OpenMP threads" << std::endl;
	int ii, j;
	omp_set_num_threads(nThreads);
	
		const double start{ omp_get_wtime() };

#pragma omp parallel for schedule(dynamic)
		for (ii = 0; ii < datumSize; ++ii){
			double Fx{ 0.0 }; double Fy{ 0.0 }; double Fz{ 0.0 };
			//std::cout << "ii loop\n";
			for (j = 0; j < datumSize; ++j) {
				//std::cout << "j loop\n";
				if (j != ii){
#if 0
					++nIters;
#endif
					const double dx{ X[j] - X[ii] };
					const double dy{ Y[j] - Y[ii] };
					const double dz{ Z[j] - Z[ii] };
					const double rdSquared{ dx * dx + dy * dy + dz * dz };
#if 0
					std::cout << "rdSquared: " << rdSquared << std::endl;
#endif
					/*
					1.0 / (rdSquared * sqrt(rdSquared) created large hotspot
					which consumed approximately ~45% of CPU time.
					*/
					//const double dCbrtSquared{ 1.0 / (rdSquared * sqrt(rdSquared)) };
					const double dCbrtSquared{ pow(rdSquared, 1.5) };
					Fx += dx * dCbrtSquared;
					Fy += dy * dCbrtSquared;
					Fz += dz * dCbrtSquared;
#if 0
					std::cout << "Fx: " << Fx << "Fy: " << Fy << "Fz: " << Fz << std::endl;
#endif
				}
			}
			VX[ii] = Fx * dt;
			VY[ii] = Fy * dt;
			VZ[ii] = Fz * dt;
			
		}
		for (int i{ 0 }; i != datumSize; i += 8){
			//std::cout << "j loop\n";
			_mm256_store_pd(&X[i], _mm256_add_pd(_mm256_load_pd(&VX[i]), _mm256_set1_pd(dt)));
			_mm256_store_pd(&X[i + 4], _mm256_add_pd(_mm256_load_pd(&VX[i + 4]), _mm256_set1_pd(dt)));
			_mm256_store_pd(&Y[i], _mm256_add_pd(_mm256_load_pd(&VY[i]), _mm256_set1_pd(dt)));
			_mm256_store_pd(&Y[i + 4], _mm256_add_pd(_mm256_load_pd(&VY[i + 4]), _mm256_set1_pd(dt)));
			_mm256_store_pd(&Z[i], _mm256_add_pd(_mm256_load_pd(&VZ[i]), _mm256_set1_pd(dt)));
			_mm256_store_pd(&Z[i + 4], _mm256_add_pd(_mm256_load_pd(&VZ[i + 4]), _mm256_set1_pd(dt)));
			std::cout << X[i] << std::endl;
		}
		const double end{ omp_get_wtime() };
		if (nSteps > skipSteps){
			rate += 1.0 / (end - start);
			d_rate += 1.0 / ((end - start) * (end - start));
		}
		
		std::cout << std::setprecision(10) << std::showpoint << std::fixed << "Time per step= " << end - start << std::endl;
		//std::cout << std::setprecision(10) << std::showpoint << std::fixed << "Time Delta=" << (end - start) << std::endl;
	
	rate /= static_cast<double>(nSteps - skipSteps);
	d_rate = sqrt(d_rate / static_cast<double>(nSteps - skipSteps) - (rate * rate));
	
	std::cout << std::setprecision(10) << std::showpoint << std::fixed << "Avrage of rate is: " << rate << "Average of delta rate: " << d_rate << std::endl;

#if 0
	std::cout << "Dumping the results \n\n";
	for (int i{ 0 }; i != X.size(); ++i)
		std::cout << std::setw(4) << X[i] << std::setw(19) << Y[i] << std::setw(25) << Z[i] << std::endl;
#endif
}
		

auto   test::ArraysMT::PartSimVector(_In_ const int nParticles, _In_ const int nSteps, _In_ const int nThreads,
	_In_ const double dt)->void {

	if ((nParticles <= 512) || (nSteps <= 1) || (nThreads <= 2) || (dt == 0.0))
		throw std::runtime_error("Invalid argument(s) in PartSimVector\n");
	std::cout << "Construction of 6 " << typeid(std::vector<double>).name() << "vectors...";
	std::vector<double> X(nParticles);
	std::vector<double> Y(nParticles);
	std::vector<double> Z(nParticles);
	std::vector<double> VX(nParticles,0.0);
	std::vector<double> VY(nParticles,0.0);
	std::vector<double> VZ(nParticles,0.0);
	std::cout << "Construction completed" << std::endl;
	std::cout <<"Allocated: " <<6 * nParticles * sizeof(double) << " heap memory bytes\n";
	clock_t seed{ clock() };
	auto rand = std::bind(std::uniform_real_distribution < double > {0.01, 1.0}, std::default_random_engine(seed));
	const std::string RNG_TYPE{ typeid(rand).name() };
	std::cout << "Created Random number generator of type: " << RNG_TYPE.c_str() << "with the interval of randomness" << "[" << -1.0 << " -- " << 1.0 << "]" << std::endl;
			

	std::cout << "std::array initialization started\n";
	//std::cout << "2x unrolling in use\n";
	const double start_init{ ::omp_get_wtime() };
	//omp_set_num_threads(3);
//#pragma omp parallel for schedule(runtime)
	for (int i{ 0 }; i != nParticles; i += 4) {
		_mm256_store_pd(&X[i], _mm256_set_pd(rand(), rand(), rand(), rand()));

		_mm256_store_pd(&Y[i], _mm256_set_pd(rand(), rand(), rand(), rand()));

		_mm256_store_pd(&Z[i], _mm256_set_pd(rand(), rand(), rand(), rand()));

	}
	const double end_init{ ::omp_get_wtime() };
	std::cout << "Initialization loop executed in: " << (end_init - start_init) << std::endl;
		
	
	std::cout << "Initialization completed" << std::endl;
	std::cout << "Made " << 12 * nParticles  << " calls to " << RNG_TYPE.c_str() << std::endl;
	double rate{ 0.0 };
	double d_rate{ 0.0 };
	const int skipSteps{ 1 };
	std::cout << "Starting " << nSteps << "Simple particle simulation" << std::endl;
	std::cout << "Parallelizing outer loop with " << nThreads << " OpenMP threads" << std::endl;
	int ii, j;
	
	omp_set_num_threads(nThreads);
	for (int i{ 1 }; i != nSteps; ++i){
		const double start{ omp_get_wtime() };

#pragma omp parallel for schedule(static,2046)
		for (ii = 0; ii < nParticles; ++ii){
			//double Fx{ 0.0 }; double Fy{ 0.0 }; double Fz{ 0.0 };
			__m256d Fx = _mm256_setzero_pd();
			__m256d Fy = _mm256_setzero_pd();
			__m256d Fz = _mm256_setzero_pd();
			//std::cout << "ii loop\n";
			for (j = 0; j < nParticles; j += 4) {
				//std::cout << "j loop\n";
				if (j != ii){
#if 0
					++nIters;
#endif
					
					_mm_prefetch(reinterpret_cast<const char*>(&X[ii + 8]), _MM_HINT_T0);
					//const double dx{ X[j] - X[ii] };
					const double t0 = X[ii];
					const __m256d dx = _mm256_sub_pd(_mm256_loadu_pd(&X[j]), _mm256_set1_pd(t0));
					_mm_prefetch(reinterpret_cast<const char*>(&Y[ii + 8]), _MM_HINT_T0);
					//const double dy{ Y[j] - Y[ii] };
					const double t1 = Y[ii];
					const __m256d dy = _mm256_sub_pd(_mm256_loadu_pd(&Y[j]), _mm256_set1_pd(t1));
					_mm_prefetch(reinterpret_cast<const char*>(&Z[ii + 8]), _MM_HINT_T0);
					//const double dz{ Z[j] - Z[ii] };
					const double t2 = Z[ii];
					const __m256d dz = _mm256_sub_pd(_mm256_loadu_pd(&Z[j]), _mm256_set1_pd(t2));
					//const double rdSquared{ dx * dx + dy * dy + dz * dz };
					const __m256d rdSquared = _mm256_mul_pd(_mm256_fmadd_pd(dx, dx, dy), _mm256_fmadd_pd(dz, dz, dy));
					/*
					1.0 / (rdSquared * sqrt(rdSquared) created large hotspot
					which consumed approximately ~45% of CPU time.
					*/
					//const double dCbrtSquared{ 1.0 / (rdSquared * sqrt(rdSquared)) };
					
					__m256d  dCbrtSqrt = _mm256_mul_pd(rdSquared, _mm256_castps_pd(_mm256_rsqrt_ps(_mm256_castpd_ps(rdSquared))));
					// __m256d dCbrtSqrt = _mm256_castps_pd( _mm256_rcp_ps( _mm256_castpd_ps( _mm256_mul_pd(rdSquared, _mm256_sqrt_pd(rdSquared)))));
					//const double dCbrtSquared{ pow(rdSquared, 1.5) };
					//Fx += dx * dCbrtSquared;
					 Fx = _mm256_fmadd_pd(dx, dCbrtSqrt, Fx);
					 Fy = _mm256_fmadd_pd(dy, dCbrtSqrt, Fy);
					 Fz = _mm256_fmadd_pd(dz, dCbrtSqrt, Fz);

				}
			}
			_mm256_storeu_pd(&VX[ii], _mm256_mul_pd( Fx , _mm256_set1_pd(dt)));
			_mm256_storeu_pd(&VY[ii], _mm256_mul_pd( Fy ,_mm256_set1_pd( dt)));
			_mm256_storeu_pd(&VZ[ii], _mm256_mul_pd(Fz , _mm256_set1_pd(dt)));
		}
		for (int i{ 0 }; i != nParticles; i += 8){
			//std::cout << "j loop\n";
			_mm256_store_pd(&X[i], _mm256_add_pd(_mm256_load_pd(&VX[i]), _mm256_set1_pd(dt)));
			_mm256_store_pd(&X[i + 4], _mm256_add_pd(_mm256_load_pd(&VX[i + 4]), _mm256_set1_pd(dt)));
			_mm256_store_pd(&Y[i], _mm256_add_pd(_mm256_load_pd(&VY[i]), _mm256_set1_pd(dt)));
			_mm256_store_pd(&Y[i + 4], _mm256_add_pd(_mm256_load_pd(&VY[i + 4]), _mm256_set1_pd(dt)));
			_mm256_store_pd(&Z[i], _mm256_add_pd(_mm256_load_pd(&VZ[i]), _mm256_set1_pd(dt)));
			_mm256_store_pd(&Z[i + 4], _mm256_add_pd(_mm256_load_pd(&VZ[i + 4]), _mm256_set1_pd(dt)));
		}
		const double end{ omp_get_wtime() };
		if (nSteps > skipSteps){
			rate += 1.0 / (end - start);
			d_rate += 1.0 / ((end - start) * (end - start));
		}
		std::cout << std::setw(12) << "Step no:" << i << std::endl;
		std::cout << std::setprecision(10) << std::showpoint << std::fixed << "Time per step= " << end - start << std::endl;
		//std::cout << std::setprecision(10) << std::showpoint << std::fixed << "Time Delta=" << (end - start) << std::endl;
	}
	rate /= static_cast<double>(nSteps - skipSteps);
	d_rate = sqrt(d_rate / static_cast<double>(nSteps - skipSteps) - (rate * rate));
	std::cout << std::setprecision(10) << std::showpoint << std::fixed << "Average rate is: " << rate << " Delta of Average rate is: " << d_rate << std::endl;

}
	
auto      test::ArraysMT::ParsSimDynArray(_In_ const int nParticles, _In_ const int nSteps, _In_ const int nThreads, _In_ const double dt)->void {
	if ((nParticles <= 512) || (nSteps <= 1) || (nThreads <= 2) || (dt == 0.0))
		throw std::runtime_error("Invalid argument(s) in PartSimVector\n");
	std::cout << "Allocation of 6 dynamic arrays started... ";
	double*  X =  reinterpret_cast<double*>(_mm_malloc(nParticles * sizeof(double),64));

	double*  Y =  reinterpret_cast<double*>(_mm_malloc(nParticles * sizeof(double),64));

	double*  Z = reinterpret_cast<double*>(_mm_malloc(nParticles * sizeof(double), 64));

	double*  VX = reinterpret_cast<double*>(_mm_malloc(nParticles * sizeof(double), 64));

	double*  VY = reinterpret_cast<double*>(_mm_malloc(nParticles * sizeof(double), 64));

	double*  VZ = reinterpret_cast<double*>(_mm_malloc(nParticles * sizeof(double), 64));
	
	if ((X == NULL) || (Y == NULL) || (Z == NULL) || (VX == NULL) || (VY == NULL) || (VZ == NULL)){
		std::cerr << "Failed to allocate heap memory" << std::endl;
		return;
	}

	std::cout << "Done!!" << std::endl;
	
	std::cout << "Allocated: " << 6 * nParticles * sizeof(double) << " bytes of heap memory" << std::endl;
	clock_t seed{ clock() };
	auto rand = std::bind(std::uniform_real_distribution < double > {0.01, 1.0}, std::default_random_engine(seed));


	const std::string RNG_TYPE{ typeid(rand).name() };
	std::cout << "Created Random number generator of type: " << RNG_TYPE.c_str() << "with the interval of randomness" << "[" << -1.0 << " -- " << 1.0 << "]" << std::endl;
	std::cout << "Initialization of particles vector position started... ";
	const __m256d  ZERO = _mm256_setzero_pd();
	const __m256d DT = _mm256_set1_pd(dt);
	for (int i{ 0 }; i != nParticles; i += 4){

		_mm256_store_pd(&X[i], _mm256_set_pd(rand(), rand(), rand(), rand()));
		_mm256_store_pd(&Y[i], _mm256_set_pd(rand(), rand(), rand(), rand()));
		_mm256_store_pd(&Z[i], _mm256_set_pd(rand(), rand(), rand(), rand()));


		_mm256_store_pd(&VX[i], ZERO);
		_mm256_store_pd(&VY[i], ZERO);
		_mm256_store_pd(&VZ[i], ZERO);
	}
	std::cout << "Done!!" << std::endl;
	double rate{ 0.0 };
	double d_rate{ 0.0 };
	const int skipSteps{ 1 };
	std::cout << "Starting " << nSteps << "Simple particle simulation" << std::endl;
	std::cout << "Parallelizing outer loop with " << nThreads << " OpenMP threads" << std::endl;
	int ii, j;
	omp_set_num_threads(nThreads);
	//for (int i{ 1 }; i != nSteps; ++i){
		const double start{ omp_get_wtime() };

#pragma omp parallel for schedule(static,2046)
		for (ii = 0; ii < nParticles; ii += 4){
			//double Fx{ 0.0 }; double Fy{ 0.0 }; double Fz{ 0.0 };
			__m256d Fx = _mm256_setzero_pd();
			__m256d Fy = _mm256_setzero_pd();
			__m256d Fz = _mm256_setzero_pd();
			__m256d tx = _mm256_setzero_pd();
			__m256d ty = _mm256_setzero_pd();
			__m256d tz = _mm256_setzero_pd();
			//std::cout << "ii loop\n";
			for (j = 0; j < nParticles; j += 4) {
				//std::cout << "j loop\n";
				if (j != ii){


					_mm_prefetch(reinterpret_cast<const char*>(&X[j + 2]), _MM_HINT_T0);
					//const double dx{ X[j] - X[ii] };
					//const double t0 = X[ii];

					 __m256d dx = _mm256_sub_pd(_mm256_load_pd(&X[j]), _mm256_load_pd(&X[ii]));
					_mm_prefetch(reinterpret_cast<const char*>(&Y[j + 2]), _MM_HINT_T0);
					//const double dy{ Y[j] - Y[ii] };
					//const double t1 = Y[ii];
					 __m256d dy = _mm256_sub_pd(_mm256_load_pd(&Y[j]), _mm256_load_pd(&Y[ii]));
					_mm_prefetch(reinterpret_cast<const char*>(&Z[j + 2]), _MM_HINT_T0);
					//const double dz{ Z[j] - Z[ii] };
					//const double t2 = Z[ii];
					 __m256d dz = _mm256_sub_pd(_mm256_load_pd(&Z[j]), _mm256_load_pd(&Z[ii]));
					//const double rdSquared{ dx * dx + dy * dy + dz * dz };
					//const __m256d rdSquared = _mm256_mul_pd(_mm256_fmadd_pd(dx, dx, dy), _mm256_fmadd_pd(dz, dz, dy));
					const __m256d rdSquared = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(dx, dx), _mm256_mul_pd(dy, dy)), _mm256_mul_pd(dz, dz));
#if 0
					std::printf("%.7f,%.7f,%.7f,%.7f\n", rdSquared.m256d_f64[0], rdSquared.m256d_f64[0],
						rdSquared.m256d_f64[2], rdSquared.m256d_f64[3]);
#endif
					/*
					1.0 / (rdSquared * sqrt(rdSquared) created large hotspot
					which consumed approximately ~45% of CPU time.
					*/
					//const double dCbrtSquared{ (rdSquared * 1.0/sqrt(rdSquared)) };
					// pow(rdSquared,1.5);
					//__m256d  dCbrtSqrt = _mm256_mul_pd(rdSquared, _mm256_castps_pd(_mm256_rsqrt_ps(_mm256_castpd_ps(rdSquared))));
					__m256d dCbrtSqrt = _mm256_castps_pd(_mm256_rcp_ps(_mm256_castpd_ps(_mm256_mul_pd(rdSquared, _mm256_sqrt_pd(rdSquared)))));
					//const double dCbrtSquared{ pow(rdSquared, 1.5) };
					//Fx += dx * dCbrtSquared;
#if 0
					std::printf("%.15f%,%.15f,%.15f,%.15f\n", dCbrtSqrt.m256d_f64[0],
						dCbrtSqrt.m256d_f64[1], dCbrtSqrt.m256d_f64[2], dCbrtSqrt.m256d_f64[3]);
#endif
					Fx = _mm256_fmadd_pd(dx, dCbrtSqrt, Fx);

					Fy = _mm256_fmadd_pd(dy, dCbrtSqrt, Fy);
					Fz = _mm256_fmadd_pd(dz, dCbrtSqrt, Fz);
					// reduction 
					tx = _mm256_add_pd(tx, Fx);

					ty = _mm256_add_pd(ty, Fy);

					tz = _mm256_add_pd(tz, Fz);

#if 0
					std::printf("%.15f%,%.15f,%.15f,%.15f\n", Fx.m256d_f64[0],
						Fx.m256d_f64[1], Fx.m256d_f64[2], Fx.m256d_f64[3]);
#endif
				}
			}
			_mm256_store_pd(&VX[ii], _mm256_mul_pd(_mm256_hadd_pd(tx, tx), DT));
			_mm256_store_pd(&VY[ii], _mm256_mul_pd(_mm256_hadd_pd(ty, ty), DT));
			_mm256_store_pd(&VZ[ii], _mm256_mul_pd(_mm256_hadd_pd(tz, tz), DT));
#if 0
			std::cout << "VX: " << VX[ii] << std::endl;
			std::cout << "VY: " << VY[ii] << std::endl;
			std::cout << "VZ: " << VZ[ii] << std::endl;
#endif 
		}

		for (int i{ 0 }; i != nParticles; i += 8){
			//std::cout << "j loop\n";
			_mm256_store_pd(&X[i], _mm256_add_pd(_mm256_load_pd(&VX[i]), DT));
			_mm256_store_pd(&X[i + 4], _mm256_add_pd(_mm256_load_pd(&VX[i + 4]), DT));
			_mm256_store_pd(&Y[i], _mm256_add_pd(_mm256_load_pd(&VY[i]), DT));
			_mm256_store_pd(&Y[i + 4], _mm256_add_pd(_mm256_load_pd(&VY[i + 4]), DT));
			_mm256_store_pd(&Z[i], _mm256_add_pd(_mm256_load_pd(&VZ[i]), DT));
			_mm256_store_pd(&Z[i + 4], _mm256_add_pd(_mm256_load_pd(&VZ[i + 4]), DT));
			

		}
		const double end{ omp_get_wtime() };

		/*if (nSteps > skipSteps){
			rate += 1.0 / (end - start);
			d_rate += 1.0 / ((end - start) * (end - start));
		}*/
		//std::cout << std::setw(12) << "Step no:" << i << std::endl;
		std::cout << std::setprecision(10) << std::showpoint << std::fixed << "Time per step= " << end - start << std::endl;
		//std::cout << std::setprecision(10) << std::showpoint << std::fixed << "Time Delta=" << (end - start) << std::endl;


	//}

	//rate /= static_cast<double>(nSteps - skipSteps);
	//d_rate = sqrt(d_rate / static_cast<double>(nSteps - skipSteps) - (rate * rate));
	//std::cout << std::setprecision(10) << std::showpoint << std::fixed << "Average rate is: " << rate << " Delta of Average rate is: " << d_rate << std::endl;

	//delete[]X; delete[]Y; delete[]Z;
	//delete[]VX; delete[]VY; delete[]VZ;
		_mm_free(X); _mm_free(Y); _mm_free(Z);
		_mm_free(VX); _mm_free(VY); _mm_free(VZ);
}

	

  
/*auto    test::ArraysMT::stdMthreading()->void {

	std::thread t1([]()->void { std::cout << "1\n"; });
	std::thread t2([]()->void {std::cout << "2\n"; });
	std::thread t3([]()->void {std::cout << "3\n"; });
	
	//t1.join();
	//t2.join();
	//t3.join();
	
}*/
	

		
	

