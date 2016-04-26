#ifndef _NBODY_VEC_AOS_IMPL_H_04_22_16
#define _NBODY_VEC_AOS_IMPL_H_04_22_16 0x100

#include "NBodyInterface.h"
#include "NBodyAOS.h"
#include "NBody_Helpers.h"

namespace  nbody_implementation {

	/********************************************
	  This derived class implements abstract
	  functions of base interface class.
	*********************************************/

	class NBodyVecAos final : private nbody_interface::NBody {


	public:

		/************************************************
		         Construction and Destruction
		*************************************************/

		/* explicitly default Ctor */
		NBodyVecAos() = default;

		/* "Main" class Ctor initializes NBodies */
		NBodyVecAos( _In_ const std::size_t, _In_ const int, _In_ const double)noexcept(false);

		

		/* Copy - Ctor */
		NBodyVecAos(_In_ const NBodyVecAos &);

		/* Move - Ctor */
		NBodyVecAos(_In_ NBodyVecAos &&);

		/* Class Dtor = default */
		~NBodyVecAos() = default;

		/*************************************
		Overriden pure abstract functions
		*************************************/

		/* Single threaded implementation of NBody
		algorithm. Runtime complexity O(n^2)*/
		virtual    void   run_simulation(_In_  double, _In_  double)  override;

		/* Multithreaded implementation of NBody
		algorithm. Runtime complexity O(n^2)*/
		virtual    void   run_simulation_omp(_In_ const double, _In_ const double) override;

		/* Single threaded implementation of NBody
		algorithm. Runtime complexity O(n^2) */
		virtual    void   run_simulation_scalar(_In_ const double, _In_ const double) override;

		/* Display final particles velocity vectors.
		Runtime complexity O(n) */
		virtual   void   display_state() const override;

		/*Non virtual member function, dumps range address of Position vectors.*/
		void      dump_address_range() const;


	private:

		/* class member m_uinParticles, denotes 
		   number of particles in simualtion.*/
		std::size_t m_uinParticles;

		/* class member m_uiThreads, denotes
		   number of OpeMP threads created
		   to run the simulation.*/
		int m_inThreads;

		/* class member m_dnStep, denotes
		   simulation time step increment.*/
		double m_dnStep;

		/* class member m_vNBodies, denotes vector 
		   container of   NBodyObject types*/
		  
		std::vector<NBodyObject> m_vNBodies;
	};

	
}

#endif /*_NBODY_VEC_AOS_IMPL_H_04_22_16*/