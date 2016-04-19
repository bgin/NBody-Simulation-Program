#ifndef _NBODY_VECTOR_IMPL_H_04_14_16
#define _NBODY_VECTOR_IMPL_H_04_14_16

#include "NBodyInterface.h"
/*****************************************
  Implementation of NBody Interface
 Based on std::vector STL container.
*****************************************/
namespace nbody_implementation {



	class NBodyVec : private  nbody_interface::NBody {


	public:


		/*
		      Default CTor 
		*/
		NBodyVec() = default;

		/*
		   "Main" class Ctor - initializes class member variables
		   nParticles, nThreads, Step.
		*/
		NBodyVec(_In_ const std::size_t, _In_ const int, _In_ const double)noexcept(false);

		/*
		    Copy-Ctor
		*/
		NBodyVec(_In_ const NBodyVec &);

		/*
		    Move-Ctor
		*/
		NBodyVec(_In_ NBodyVec &&);

		/*
		    Dtor
		*/
		~NBodyVec() = default;

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

		

	private:

		/*
		   std::vector<double> m_X -- X position vector.
		*/
		std::vector<double> m_vX;

		/*
		    std::vector<double> m_Y -- Y position vector.
		*/
		std::vector<double> m_vY;

		/*
		    std::vector<double> m_Z --  Z postion vector.
		*/
		std::vector<double> m_vZ;
		/*
		     std::vector<double> m_VX -- VX velocity vector.
		*/
		std::vector<double> m_vVX;

		/*
		     std::vector<double> m_VY -- VY velocity vector.
		*/
		std::vector<double> m_vVY;

		/*
		     std::vector<double> m_VZ --  VZ velocity vector.
		*/
		std::vector<double> m_vVZ;

		/*
		   Number of Particles 
		*/
		std::size_t m_uinParticles;

		/*
		   Number of OpenMP threads
		*/
		int  m_inThreads;

		/*
		    Simulation time-step
		*/
		double  m_dStep;

		/* Dump object state */
#if defined VERBOSE
		auto nbody_vec_internal_state()const ->void;
#endif
	};

}
#endif /*_NBODY_VECTOR_IMPL_H_04_14_16*/