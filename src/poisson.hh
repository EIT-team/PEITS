#ifndef POISSON_PROBLEMS_HH
#define POISSON_PROBLEMS_HH

#include <cassert>
#include <cmath>

#include "probleminterface.hh"

// -laplace u = 0 is the problem we want to solve in EIT
template <class FunctionSpace> 
class EIT_solver : public ProblemInterface < FunctionSpace >
{
  typedef ProblemInterface < FunctionSpace >  BaseType;
public:
  typedef typename BaseType :: RangeType            RangeType;
  typedef typename BaseType :: DomainType           DomainType;
  typedef typename BaseType :: JacobianRangeType    JacobianRangeType;
  typedef typename BaseType :: DiffusionTensorType  DiffusionTensorType;

  enum { dimRange  = BaseType :: dimRange };
  enum { dimDomain = BaseType :: dimDomain };

  //! the right hand side data
  virtual void f(const DomainType& x,
                 RangeType& phi) const
  {
    phi = 1;
  }

  //! the exact solution
  virtual void u(const DomainType& x,
                 RangeType& phi) const
  {
	phi = 0;
	for( int i=0; i < dimDomain; ++i )
	    phi += 2*x[i] + 3;
  }

  //! the jacobian of the exact solution
  virtual void uJacobian(const DomainType& x,
                         JacobianRangeType& ret) const
  {
	for( int i = 0; i < dimDomain; ++i )
	{
	  ret[ 0 ][ i ] = 2;
	}
  }
  //! mass coefficient has to be 1 for this problem
  virtual void m(const DomainType& x, RangeType &m) const
  {
    m = 1;
  }
};

/*********************************************************/

#endif
