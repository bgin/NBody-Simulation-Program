#ifndef _NBODY_HELPERS_H_04_25_16
#define _NBODY_HELPERS_H_04_25_16 0x100

/********************************************
  Class containing few static helpers for
  checking Ctors arguments in initializing
  list.
********************************************/
#include "NbodyDefs.h"
#include <cstddef>

namespace nbody_implementation {

	class NBodyHelpers {


	public:

		/* Checks the amount of particles passed as an argument
		   to various Ctors.
		   throws an exception of type std::runtime_error*/
		static auto check_particles(_In_ const std::size_t)noexcept(false)->void;

		/* Checks the number of OpenMP threads passed as an
		   argument to various Ctors.
		   throws an exception of type std::runtime_error*/
		static  auto check_threads(_In_ const int)noexcept(false)->void;

		/* Checks the value of increament Step passed to various
		   to various Ctors.
		   throws an exception of type std::runntime*/
		static  auto  check_step(_In_ const double)noexcept(false)->void;
	};
}
#endif  /* _NBODY_HELPERS_H_04_25_16*/