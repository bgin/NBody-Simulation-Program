#ifndef _ARRAYS_TEST_H_
#define _ARRAYS_TEST_H_


#include <omp.h>
#include <iostream>
#include <random>
#include <valarray>
#include <immintrin.h>
#include <ctime>
#include <functional>
#include <iomanip>


namespace   test {


	class  ArraysMT{

	public:


	static	auto  allocateVec1D(_In_ const int)->std::valarray<double>;

		

	static	auto  initVec1D(_In_ std::valarray<double> &)->void;

	static	auto  runComputation(_In_ std::valarray<double> &, _Out_ std::valarray<double> &)->void;

	static  auto  printData(_In_ std::valarray<double> &)->void;

	static  auto  testArray2D()->void;

	static  auto  simple_test()->void;

	/*
	     Simple particle simulation. 
		 Data type: static array.
	*/
	   static  auto  PartSimArray( _In_ const int, _In_ const int, _In_ const double)->void;

	/*
	     Simple particle simulation.
		 Data typ: std::valarray<double>
	*/
	   static  auto  PartSimVArray(_In_ const std::size_t, _In_ const int, _In_ const int, _In_ const double)->void;
	
	/*
	     Simple particle simulation.
		 Datat type: std::array.
	*/
	   static  auto  PartSimStdArray(_In_ const int, _In_ const int, _In_ const double)->void;

	/*
	     Simple particle simulation.
		 Data type: std::vector.
	*/
	static  auto  PartSimVector(_In_ const int, _In_ const int, _In_ const int, _In_ const double)->void;

	/*
	     Using C-style dynamic arrays allocated by the call to mm_malloc
	*/
	static  auto  ParsSimDynArray(_In_ const int, _In_ const int, _In_ const int, _In_ const double)->void;
	};
}
#endif  /*_ARRAYS_TEST_H_*/

