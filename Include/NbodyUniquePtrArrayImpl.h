
#ifndef _NBODY_UNIQUE_PTR_ARRAY_IMPL_H_04_21_16
#define _NBODY_UNIQUE_PTR_ARRAY_IMPL_H_04_21_16

#include "NBodyInterface.h"

/**************************************************************
 Implementation of NBody Interface. This implementation is based
 on smart pointers wrapping naked new allocation.
***************************************************************/

namespace nbody_implementation {

	class NBodyUptrArray final : nbody_interface::NBody {

		/***************************************************
							Constructors
							Obviously there is no need for Destructor becuse
							of usage of smart pointers.
							****************************************************/


		/************************************************
				  Only one Ctor in this class
				  ************************************************/

	public:



		NBodyUptrArray(_In_ const std::size_t, _In_ const int, _In_ const double);

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

		/* class member m_uinParticles denotes number
		   of Particles */
		std::size_t  m_uinParticles;

		/* class member m_inThreads denotes a number
		   of OpenMP threads  scheduled to run the
		   simulation.
		   */
		int   m_inThreads;

		/* class member m_dnStep denotes number of step size.*/
		double m_dStep;

		/* class member m_spX denotes Position vector X
		   wrapped inside unique_ptr */
		std::unique_ptr<double, void(*)(double*)> m_spX;

		/* class member m_spY denotes Position vector Y
		   wrapped inside unique_ptr. */
		std::unique_ptr<double, void(*)(double*)> m_spY;

		/* class member m_spZ denotes Position vector Z
		   wrapped inside unique_ptr.*/
		std::unique_ptr<double, void(*)(double*)> m_spZ;

		/* class member m_spVX denotes Velocity vector
		   VX wrpped inside unique_ptr.*/
		std::unique_ptr<double, void(*)(double*)> m_spVX;

		/* class member m_spVY denotes Velocity vector
		   VY wrapped inside unique_ptr. */
		std::unique_ptr<double, void(*)(double*)> m_spVY;

		/* class member m_spVZ denotes Velocity vector
		   VZ wrapped inside unique_ptr.*/
		std::unique_ptr<double, void(*)(double*)> m_spVZ;

		/* simple static functions for input checking called from
		   initialization list*/

		/* check NParticles argument
		   throws std::runtime_error*/
		static void  checkParticles(_In_ const std::size_t) noexcept(false);

		/* check nThreads argument
		   throws std::runtime_error*/
		static void  checkThreads(_In_ const int) noexcept(false);

		/* check nStep argument 
		   throws std::runtime_error */
		static void   checkStep(_In_ const double) noexcept(false);
	};

}


#endif /*_NBODY_UNIQUE_PTR_ARRAY_IMPL_H_04_21_16*/