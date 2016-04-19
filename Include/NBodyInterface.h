#ifndef _NBODY_INTERFACE_H_04_14_16
#define _NBODY_INTERFACE_H_04_14_16

#include "NbodyDefs.h"
/*************************************
  Pure abstract class for NBody
  Simulation.
**************************************/
namespace nbody_interface {

	class NBody {

	public:

		/*****************************
		  Pure  Abstract	Methods
		******************************/

		/* This function runs the simualtion overriden by the subclass*/
		virtual	void  run_simulation(_In_ const double, _In_ const double) = 0;

		/* This function runs the multithreaded simulation overriden by the subclass */
		virtual void  run_simulation_omp(_In_ const double, _In_ const double) = 0;

		/* This function runs the single threaded scalar simulation, overriden by the sub
		 classes*/
		virtual void  run_simulation_scalar(_In_ const double, _In_ const double) = 0;

		/* This function displays the result overriden by the subclass */
		virtual void  display_state() const = 0;

		/* Virtual Dtor*/
		virtual ~NBody() {};
	};
}

#endif /*_NBODY_INTERFACE_H_04_14_16*/