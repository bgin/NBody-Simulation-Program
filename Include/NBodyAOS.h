
#ifndef _NBODY_AOS_H_04_22_16
#define _NBODY_AOS_H_04_22_16


#include <initializer_list>
#include <algorithm>
namespace  nbody_implementation {

	/************************************************
	  This class represents single NBody object.
	  This object is a part of derived
	  class compositions. This data layout
	  is known as a AOS, i.e Array of Structures
	  and from the point of cache-friendly programming
	  is not as effective as SOA layout.

	************************************************/

	class NBodyObject {

	public:

		/******************************
		  Construction and Destruction
		*******************************/

		/* Default Ctor initiliazes NBodyObject
		   to null state. */
		NBodyObject()noexcept(true);

		/* "Main" class Ctor initializes
		   NBodyObject from the user passed
		   scalar arguments.*/
		NBodyObject(_In_ const double, _In_ const double, _In_ const double,
			_In_ const double, _In_ const  double, _In_ const double) noexcept(true);

		/* "Main" class Ctor initializes
		   NBodyObject from the user passed
		   std::initializer_list*/
		NBodyObject(_In_ const std::initializer_list<double> &);

		/* Copy-Ctor */
		NBodyObject(_In_ const NBodyObject &)noexcept(true);

		/* Move-Ctor */
		NBodyObject(_In_ NBodyObject &&)noexcept(true);

		/* Class Dtor explicitly default */
		~NBodyObject() = default;

		/* operator copy assignment */
		NBodyObject & operator=(_In_ const NBodyObject &);

		/* operator move assignment */
		NBodyObject & operator=(_In_ NBodyObject &&);

		// Public members for convenience of access.
		// Can be used setters/getters also 
	public:

		/* class member m_dX, denotes particle's 
		   position at time t0...tn. */
		double   m_dX;

		/*b class member m_dY, denotes particle's
		  position at time t0...tn.*/
		double   m_dY;

		/* class member m_dZ, denotes particle's
		 position at time t0...tn.*/
		double   m_dZ;

		/* class member m_dVX, denotes particle's
		 velocity at time t0...tn.*/
		double   m_dVX;

		/* class member m_dVY, denotes particle's
		 velocity at time t0...tn.*/
		double   m_dVY;

		/* class member m_dVZ, denotes particle's
		 velocity at time t0...tn*/
		double   m_dVZ;
	};
}
#endif /*_NBODY_AOS_H_04_22_16*/