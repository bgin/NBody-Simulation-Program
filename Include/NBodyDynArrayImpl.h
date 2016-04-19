#ifndef _NBODY_DYN_ARRAY_IMPL_H_04_14_16
#define _NBODY_DYN_ARRAY_IMPL_H_04_14_16

#include "NBodyInterface.h"

namespace nbody_implementation {

	/***********************************************
	 Implemetation of NBody Interface based on 
	 dynamically allocated arrays by the operator new.
	*************************************************/
	class NBodyDArray final : private nbody_interface::NBody {

    

	public:

		/*  Default Ctor = deleted*/
		NBodyDArray() = delete;

		/* " Main" class Ctor initilizes class members and allocates resources */
		NBodyDArray(_In_ const std::size_t, _In_ const int, _In_ const double)noexcept(false);

		/* Copy Ctor */
		NBodyDArray(_In_ const NBodyDArray &);

		/* Move Ctor */
		NBodyDArray(_In_ NBodyDArray &&);

		/* Dtor*/
		virtual ~NBodyDArray();

		/*************************************
		Overriden pure abstract functions
		*************************************/

		/* Single threaded implementation of NBody
		algorithm. Runtime complexity O(n^2)*/
		virtual    void   run_simulation(_In_ const double, _In_ const double)  override;

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

		/*   Class member m_pX denotes a NBody
		     position vector X. */

		double*    m_pX;

		/*   Class member m_pY denotes a NBody
		     position vector Y.*/
		double*    m_pY;

		/*   Class member m_pZ denotes a NBody
		     position vector Z.*/
		double*    m_pZ;

		/*   Class member m_pVX denotes a NBody
		     velocity vector VX*/
		double*    m_pVX;

		/*   Class member m_pVY denotes a NBody
		     velocity vector VY.*/
		double*    m_pVY;

		/*   Class member m_pVZ denotes a NBody
		     velocity vector VZ.*/
		double*   m_pVZ;

		/*   Class member m_uinParticles */
		std::size_t  m_uinParticles;

		/*   Class member m_inThreads*/
		int     m_inThreads;

		/*   Class member m_dStep*/
		double  m_dStep;

		/* allocates dynamic array and returns the pointer to it
		   Caller will have to free memory.
		   Catches an exception of type std::bad_alloc 
		   Throws an exception of type std::runtime 
		   when dataSize is <= 0
		   */
		//double *allocate_array1D(_In_ const std::size_t) noexcept(false);
		
		std::tuple<double*, double*, double*, double*, double*, double*> allocate_array1D(_In_ const std::size_t) noexcept(false);

	};
}
#endif /*_NBODY_DYN_ARRAY_IMPL_H_04_14_16*/