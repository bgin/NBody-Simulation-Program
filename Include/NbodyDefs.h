#ifndef _NBODY_DEFS_H_04_14_16
#define _NBODY_DEFS_H_04_14_16

/* Frequently used header files */
#include <iostream>
#include <iomanip>
#include <valarray>
#include <vector>
#include <functional>
#include <memory>
#include <random>
#include <time.h>
#include <omp.h>

#if defined (__INTEL_COMPILER) 
#include <immintrin.h>
#else
#include <intrin.h>
#endif
#ifndef VERBOSE
#define VERBOSE 0x1
#endif

#endif /*_NBODY_DEFS_H_04_14_16*/