/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/

#pragma once

/**
 * @file AlgebraicCurveFromOrderedPoints.h
 * @brief Representation of an algebraic curve
 * @author Tristan Roussillon (\c tristan.roussillon@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 * @date 2011/09/22
 *
 * @brief Header file for module AlgebraicCurveFromOrderedPoints.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(AlgebraicCurveFromOrderedPoints_RECURSES)
#error Recursive header files inclusion detected in AlgebraicCurveFromOrderedPoints.h
#else // defined(AlgebraicCurveFromOrderedPoints_RECURSES)
/** Prevents recursive inclusion of headers. */
#define AlgebraicCurveFromOrderedPoints_RECURSES

#if !defined AlgebraicCurveFromOrderedPoints_h
/** Prevents repeated inclusion of headers. */
#define AlgebraicCurveFromOrderedPoints_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>

#include "DGtal/base/Common.h"
#include "DGtal/kernel/PointVector.h"
#include "DGtal/kernel/NumberTraits.h"
#include "DGtal/io/Color.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  namespace details
  {
    template<typename TInteger>
    struct CircleFromPoints 
    {
    public: 
      static const DGtal::Dimension N = 3; 

      typedef TInteger Integer;
      typedef DGtal::PointVector<3,Integer> Point; 
      typedef Integer Value;  //TODO to promote
     
    public: 
      Integer
      operator() (const boost::array<Point, N>& a, const Point& aP) const 
      {
 	typedef DGtal::PointVector<2,Integer> Vector; 

	Vector u( (a[0][0]-aP[0])*(a[1][1]-aP[1])-  (a[1][0]-aP[0])*(a[0][1]-aP[1]), 
		  (a[1][0]-aP[0])*(a[1][0]-a[0][0])+(a[1][1]-aP[1])*(a[1][1]-a[0][1]) );
	Vector v( (a[0][0]-aP[0])*(a[2][1]-aP[1])-  (a[2][0]-aP[0])*(a[0][1]-aP[1]), 
		  (a[2][0]-aP[0])*(a[2][0]-a[0][0])+(a[2][1]-aP[1])*(a[2][1]-a[0][1]) );
	return -( (u[0] * v[1]) - (u[1] * v[0]) ); 

      }
    }; 
  /////////////////////////////////////////////////////////////////////////////
  // template class AlgebraicCurveFromOrderedPointsBase
  /**
   * \brief Aim: Represents an algebraic curve, which 
   * is constrained to pass through @e F given points.
   *
   * @tparam TPolynome type of the polynome interpolating 
   * the @e F given points. It is a functor that returns 
   * the algebraic distance of a point to the curve.
   *
   * @tparam F number of free points (between 0 and N, 
   * the maximal number of points that uniquely define
   * a given curve)
   *
   */
    template <typename TPolynome, DGtal::Dimension F = TPolynome::N, typename TInteger = typename TPolynome::Integer>
  struct AlgebraicCurveFromOrderedPointsBase
  {

    BOOST_STATIC_ASSERT((F >= 0 ));
    BOOST_STATIC_ASSERT((F <= TPolynome::N ));

    /**
     * Total number of points required
     * to define uniquely the curve
    */
    static const DGtal::Dimension N = TPolynome::N; 
    /**
     * Total number of given points, 
     * the curve has to pass through
    */
    static const DGtal::Dimension G = N-F; 

    // ----------------------- associated types ------------------------------
    /**
     * Type of coordinate
    */
    typedef TInteger Coordinate;
    /**
     * Type of point
    */
    typedef DGtal::PointVector<2,Coordinate> Point; 
    /**
     * Type of point in homogeneous coordinates 
     * (one more coordinate, equal to 0 for point 
     * at infinity and different otherwise)
    */
    typedef DGtal::PointVector<3,Coordinate> HPoint;

    /**
     * Return type of the operator()
    */
    typedef Coordinate Value;//TODO to promote or to adapt (-1, 0, 1) ?

    typedef boost::array<HPoint,G> GArray;  //array of fixed points
    typedef boost::array<HPoint,F> FArray;  //array of variable points
    typedef boost::array<HPoint,N> Array;  //array of all points

    // ------------------------- Private Datas --------------------------------
  protected:
    /**
       Array of fixed points
    */
    GArray myGArray;
    /**
       Array of variable points
    */
    FArray myFArray;
    /**
     * Distance functor
    */
    TPolynome myPolynome;

    // ----------------------- Standard services ------------------------------
  public:


    /**
     * Constructor when no points are fixed.
     */
    AlgebraicCurveFromOrderedPointsBase() { BOOST_STATIC_ASSERT(( G == 0 )); };

    /**
     * Constructor.
     * @param aGArray an array of given points
     */
    AlgebraicCurveFromOrderedPointsBase(const GArray& aGArray): 
      myGArray( aGArray ), myFArray(), myPolynome() {};


    /**
     * Init.
     * @param itb begin iterator on points
     * @param ite end iterator on points
     */
    template <typename I>
    void init(const I& itb, const I& ite) 
      { 
	typename FArray::iterator ait = myFArray.begin(); 
	for (I it = itb; ( (it != ite)&&(ait != myFArray.end()) ); ++it, ++ait)
            *ait = toHPoint( *it );  
      };

    /**
     * Init.
     * @param aFArray an array of free points
     */
    void init(const FArray& aFArray) { myFArray = aFArray; };

    /**
     * Copy constructor.
     * @param other the object to clone.
     */
	template <typename Other>
	AlgebraicCurveFromOrderedPointsBase ( const Other & other ): 
	  myGArray( other.myGArray ), myFArray( other.myFArray ), myPolynome( other.myPolynome ) {};

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     */
    AlgebraicCurveFromOrderedPointsBase & operator= ( const AlgebraicCurveFromOrderedPointsBase & other )
      {
	if (this != &other)
          {
            myGArray = other.myGArray; 
            myFArray = other.myFArray;
	    myPolynome = other.myPolynome; 
          }
	return *this; 
      };

    /**
     * Destructor. Does nothing
     */
    ~AlgebraicCurveFromOrderedPointsBase() {};

    // ----------------------- Interface --------------------------------------
  public:

    /**
     * Fonction
     * @param aPoint any point.
     * @return value (<0, 0 on, >0 )
     */
    Value operator()(const Point& aP) const 
     {
	Array tmp; 
	//lexicographic order ?
	std::copy( myFArray.begin(), myFArray.end(), tmp.begin() ); 
	std::copy( myGArray.begin(), myGArray.end(), tmp.begin()+F ); 
	return myPolynome(tmp, toHPoint( aP ) ); 
     }

    //------------------ useful functions -------------------------------

    /**
     * Return a point to infinity (last homogeneous coordinate set to 0).
     * @return the point.
     */
    HPoint toInfinity() const
    {
      HPoint p = HPoint::diagonal(1); 
      p[2] = 0; 
      return p; 
    }


    //------------------ accessors -------------------------------
    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    bool isValid() const 
    { 
      //TODO must check if all points are not confounded (in O(N^2) )
      return true; 
    };
   

    //------------------ display -------------------------------

    /**
     * Writes/Displays the object on an output stream.
     * @param out the output stream where the object is written.
     */
    void selfDisplay ( std::ostream & out ) const
      {
         std::cout << "[AlgebraicCurveFromOrderedPointsBase] "
	           << F << " free / " << N << " points. "
                   << std::endl; 
	std::copy( myFArray.begin(), myFArray.end(), std::ostream_iterator<HPoint>(std::cout, ",") ); 
	std::copy( myGArray.begin(), myGArray.end(), std::ostream_iterator<HPoint>(std::cout, ",") ); 
      };
    

    // ------------------------- Internals ------------------------------------
  protected:

    /**
     * @return @e aPoint in homogeneous coordinates.
     */
    HPoint toHPoint(const Point& aPoint) const 
     { 
	HPoint tmp; 
	std::copy( aPoint.begin(), aPoint.end(), tmp.begin() ); 
	tmp[2] = 1;  
	return tmp; 
     };


  }; // end of class AlgebraicCurveFromOrderedPointsBase

  } //end namespace details

  /////////////////////////////////////////////////////////////////////////////
  // template class AlgebraicCurveFromOrderedPoints
  /**
   * \brief Aim: Represents an algebraic curve, which 
   * is constrained to pass through @e F given points.
   *
   * @tparam TPolynome type of the polynome interpolating 
   * the @e F given points. It is a functor that returns 
   * the algebraic distance of a point to the curve.
   *
   * @tparam T number of free points (between 0 and N, 
   * the maximal number of points that uniquely define
   * a given curve)
   *
   */
  template <typename TPolynome, DGtal::Dimension T = TPolynome::N, typename TInteger = typename TPolynome::Integer>
  class AlgebraicCurveFromOrderedPoints: 
    public details::AlgebraicCurveFromOrderedPointsBase<TPolynome, T, TInteger>
  {

    // ----------------------- Static constant ------------------------------
  public: 
    /*
     * Number of free points
     */
    static const DGtal::Dimension F = T;

    // ----------------------- Friend class ------------------------------
  public: 
    /*
     * (F+1)-class as friend for initFromUp methods
     */
    friend class AlgebraicCurveFromOrderedPoints<TPolynome, F+1, TInteger>;  

    // ----------------------- Nested type ------------------------------
  public: 
    /*
     * This F-class (curve constrained to pass through G given points)
     */
    typedef AlgebraicCurveFromOrderedPoints<TPolynome, F, TInteger> Self;  
    /*
     * Type of the parent class
     */
    typedef details::AlgebraicCurveFromOrderedPointsBase<TPolynome, F, TInteger> Super; 
    /*
     * Type of the (F-1)-class (curve contrained to pass through G+1 given points) 
     */
    typedef AlgebraicCurveFromOrderedPoints<TPolynome, F-1, TInteger> Up; 

    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Constructor when no points are fixed.
     */
    AlgebraicCurveFromOrderedPoints() { BOOST_STATIC_ASSERT(( Self::G == 0 )); };

    /**
     * Constructor.
     * @param aGArray an array of fixed points
     */
    AlgebraicCurveFromOrderedPoints(const typename Self::GArray& aGArray): Super( aGArray ) {};

    /**
     * Copy constructor.
     * @param other the object to clone.
     */
    AlgebraicCurveFromOrderedPoints ( const Self & other ): Super( other ) {};

    /**
     * Set the last free point as the first given point.
     */
    void shift() 
    {
      BOOST_STATIC_ASSERT(( Self::F > 0 ));
      this->myFArray[1] = this->myFArray[0]; 
    }

    /**
     * Init from an instance of Up
     * @param u an instance of Up
     */
    void initFromUp(const Up& u) 
      {
	/// variable points: 
	//all variable points of u
	std::copy( u.myFArray.begin(), u.myFArray.end(), this->myFArray.begin() ); 
	//+ the first (i.e. last added) fixed point of u
	ASSERT( Up::G > 0 ); 
	this->myFArray.back() = u.myGArray.front();
	/// fixed points: 
	std::copy( u.myGArray.begin()+1, u.myGArray.end(), this->myGArray.begin() ); 
      }

    /**
     * Build and return an instance of Up
     * @param aPoint a point 
     * @return an instance of Up ( @e aHPoint is taken as an extra fixed point)
     */
    Up getUp( const typename Self::Point& aPoint ) const
      {
	return getUp( this->toHPoint( aPoint ) ); 
      }

    /**
     * Build and return an instance of Up
     * @param aHPoint a point (in homogeneous coordinates)
     * @return an instance of Up ( @e aHPoint is taken as an extra fixed point)
     */
    Up getUp( const typename Self::HPoint& aHPoint ) const
      {
	boost::array<typename Self::HPoint,Self::G+1> tmp; 
	tmp.front() = aHPoint;  
	std::copy( this->myGArray.begin(), this->myGArray.end(), tmp.begin()+1 ); 
	return Up( tmp ); 
      }  

    /**
     * Destructor. Does nothing
     */
    ~AlgebraicCurveFromOrderedPoints() {};

   private: 

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     */
    Self & operator= ( const Self & other );

  }; // end of class AlgebraicCurveFromOrderedPoints

  /////////////////////////////////////////////////////////////////////////////
  // Specialization
  template <typename TPolynome, typename TInteger>
  class AlgebraicCurveFromOrderedPoints<TPolynome, 0, TInteger>: 
    public details::AlgebraicCurveFromOrderedPointsBase<TPolynome, 0, TInteger>
  {

    // ----------------------- Static constant ------------------------------
  public: 
    /*
     * number of free points (= 0)
     */
    static const DGtal::Dimension F = 0;

    // ----------------------- Friend class  ------------------------------
  public: 
    /*
     * (1)-class as friend for initFromUp methods
     */
    friend class AlgebraicCurveFromOrderedPoints<TPolynome, 1, TInteger>;  

    // ----------------------- Nested type ------------------------------
  public: 
    /*
     * Self type
     */
    typedef AlgebraicCurveFromOrderedPoints<TPolynome, 0, TInteger> Self; 
    /*
     * Type of the parent class
     */
    typedef details::AlgebraicCurveFromOrderedPointsBase<TPolynome, 0, TInteger> Super; 
    /*
     * Type of the (F-1)-class defined as an alias of Self (= 0) 
     */
    typedef Self Up;   

    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Constructor.
     * @param aGArray an array of fixed points
     */
    AlgebraicCurveFromOrderedPoints(const typename Self::GArray& aGArray): Super( aGArray ) {};

    /**
     * Copy constructor.
     * @param other the object to clone.
     */
    AlgebraicCurveFromOrderedPoints ( const Self & other ): Super( other ) {};

    /**
     * Destructor. Does nothing
     */
    ~AlgebraicCurveFromOrderedPoints() {};

    /**
     * Does nothing
     */
    void shift( )
      {
      }

    /**
     * Does nothing
     */
    void initFromUp( const Up& u ) const
      {
      }

    /**
     *
     * @param aPoint a point 
     * @return an instance of Up
     */
    Up getUp( const typename Self::Point& /*aPoint */) const
      {
	return *this; 
      }

    /**
     *
     * @param aHPoint a point (in homogeneous coordinates)
     * @return an instance of Up
     */
    Up getUp( const typename Self::HPoint& /*aHPoint*/ ) const
      {
	return *this; 
      }  

   private: 

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     */
    Self & operator= ( const Self & other );

  }; // end of class AlgebraicCurveFromOrderedPoints

  /**
   * Overloads 'operator<<' for displaying objects of class 'AlgebraicCurveFromOrderedPoints'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'AlgebraicCurveFromOrderedPoints' to write.
   * @return the output stream after the writing.
   */
  template <typename TPolynome, DGtal::Dimension T, typename TInteger>
  inline
  std::ostream&
  operator<< ( std::ostream & out, 
	       const details::AlgebraicCurveFromOrderedPointsBase<TPolynome, T, TInteger> & object )
  {
    object.selfDisplay( out );
    return out;
  }


} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
//#include "DGtal/shapes/fromPoints/AlgebraicCurveFromOrderedPoints.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined AlgebraicCurveFromOrderedPoints_h

#undef AlgebraicCurveFromOrderedPoints_RECURSES
#endif // else defined(AlgebraicCurveFromOrderedPoints_RECURSES)
