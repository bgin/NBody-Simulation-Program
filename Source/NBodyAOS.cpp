
#include "NBodyAOS.h"

/* class NBodyObject implementation  */

/*         Default Ctor     */
nbody_implementation::NBodyObject::NBodyObject() :
  m_dX{ 0.0 },
  m_dY{ 0.0 },
  m_dZ{ 0.0 },
  m_dVX{ 0.0 },
  m_dVY{ 0.0 },
  m_dVZ{ 0.0 }
{

}

/* Scalar arguments Ctor   */
nbody_implementation::NBodyObject::NBodyObject(_In_ const double X, _In_ const double Y,
	_In_ const double Z, _In_ const double VX,
	_In_ const double VY, _In_ const double VZ) :
  m_dX{ X },
  m_dY{ Y },
  m_dZ{ Z },
  m_dVX{ VX },
  m_dVY{ VY },
  m_dVZ{ VZ }
{

}

/*  Construct from the std::initializer_list*/
nbody_implementation::NBodyObject::NBodyObject(_In_ const std::initializer_list<double> &ilist) :
 m_dX{ ilist.begin()[0] },
 m_dY{ ilist.begin()[1] },
 m_dZ{ ilist.begin()[2] },
 m_dVX{ ilist.begin()[3] },
 m_dVY{ ilist.begin()[4] },
 m_dVZ{ ilist.begin()[5] }
{

}

/*  Copy Construct */
nbody_implementation::NBodyObject::NBodyObject(_In_ const NBodyObject &rhs) :
 m_dX{ rhs.m_dX },
 m_dY{ rhs.m_dY },
 m_dZ{ rhs.m_dZ },
 m_dVX{ rhs.m_dVX },
 m_dVY{ rhs.m_dVY },
 m_dVZ{ rhs.m_dVZ }
{

}

/* Move Construct */
nbody_implementation::NBodyObject::NBodyObject(_In_  NBodyObject &&rhs) :
m_dX{ std::move(rhs.m_dX) },
m_dY{ std::move(rhs.m_dY) },
m_dZ{ std::move(rhs.m_dZ) },
m_dVX{ std::move(rhs.m_dVX) },
m_dVY{ std::move(rhs.m_dVY) },
m_dVZ{ std::move(rhs.m_dVZ) }
{

}

/* operator copy-assign */
nbody_implementation::NBodyObject &  nbody_implementation::NBodyObject::operator=(_In_ const NBodyObject &rhs) {
	if (this == &rhs) return (*this);
	
	NBodyObject temp{ rhs };
	std::swap(*this, temp);
	return (*this);
}

/* operator move-assign */
nbody_implementation::NBodyObject &  nbody_implementation::NBodyObject::operator=(_In_  NBodyObject &&rhs) {
	if (this == &rhs) return (*this);
	NBodyObject temp{ rhs };
	std::swap(*this, temp);
	return (*this);
}

