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
 * @brief Representation of a AlgebraicCurveFromOrderedPoints uniquely defined by two 2D points.
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
    struct AlgebraicDistanceToCircle 
    {
    public: 
      static const DGtal::Dimension N = 3; 

    public: 
      template<typename Point>
      typename Point::Coordinate //to promote
      operator() (const boost::array<Point, N>& a, const Point& aP) const 
      {
 	typedef DGtal::PointVector<2,typename Point::Coordinate> Vector; 

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
   * is contrained to pass through @e F given points.
   *
   * @tparam F number of free points (between 0 and N, 
   * the maximal number of points that uniquely define
   * a given curve)
   *
   * @tparam TDistance type of functor, which computes 
   * the algebraic distance to the curve
   */
    template <DGtal::Dimension F, typename TInteger = int, typename TDistance = AlgebraicDistanceToCircle>
  struct AlgebraicCurveFromOrderedPointsBase
  {

    BOOST_STATIC_ASSERT((F >= 0 ));
    BOOST_STATIC_ASSERT((F <= TDistance::N ));

	//total number of points
	//(defining uniquely a circle)
    static const DGtal::Dimension N = TDistance::N; 
 	//number of given points 
    static const DGtal::Dimension G = N-F; 

    // ----------------------- associated types ------------------------------
    typedef TInteger Coordinate;
    typedef DGtal::PointVector<2,Coordinate> Point; 
    typedef DGtal::PointVector<3,Coordinate> HPoint; //in homogeneous coordinates 

    typedef Coordinate Value;//to promote or to adapt (-1, 0, 1) ?

    typedef boost::array<HPoint,G> GArray;  //array of fixed points
    typedef boost::array<HPoint,F> FArray;  //array of variable points
    typedef boost::array<HPoint,N> Array;  //array of all points

    // ------------------------- Private Datas --------------------------------
  protected:
    /**
       Array of points (first V points, variable, last F points, fixed)
    */
    //    Array myArray;

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
    TDistance myDistance;

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
    AlgebraicCurveFromOrderedPointsBase(const GArray& aGArray): myGArray( aGArray ), myDistance() {};


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
	  myGArray( other.myGArray ), myFArray( other.myFArray ), myDistance( other.myDistance ) {};

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
	    myDistance = other.myDistance; 
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
     * @return value (<0 in, 0 on, >0 out)
     */
    Value operator()(const Point& aP) const 
     {
	Array tmp; 
	//lexicographic order ?
	std::copy( myFArray.begin(), myFArray.end(), tmp.begin() ); 
	std::copy( myGArray.begin(), myGArray.end(), tmp.begin()+F ); 
	return distance( tmp, toHPoint( aP ) ); 
     }

    //------------------ accessors -------------------------------
    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    bool isValid() const { return true; };
   

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

    /**
     * @return determinant sign (<0 in, 0 on, >0 out)
     */
    Value distance(const Array& a, const HPoint& aP) const 
     { 
       return myDistance(a, aP); 
     };

  }; // end of class AlgebraicCurveFromOrderedPointsBase

  } //end namespace details

  /////////////////////////////////////////////////////////////////////////////
  // template class AlgebraicCurveFromOrderedPoints
  /**
   * \brief Aim: Represents an algebraic curve 
   *
   */
  template <DGtal::Dimension T, typename TInteger = int, typename TDistance = details::AlgebraicDistanceToCircle>
  class AlgebraicCurveFromOrderedPoints: public details::AlgebraicCurveFromOrderedPointsBase<T,TInteger,TDistance>
  {

  public: 
    /*
     * (F+1)-class as friend for initFromUp methods
     */
    friend class AlgebraicCurveFromOrderedPoints<T+1, TInteger, TDistance>;  

    // ----------------------- Static constant ------------------------------
  public: 
    /*
     * Number of free points
     */
    static const DGtal::Dimension F = T;

    // ----------------------- Nested type ------------------------------
  public: 
    /*
     * This class (curve constrained to pass through G given points)
     */
    typedef AlgebraicCurveFromOrderedPoints<F, TInteger, TDistance> Self;  
    /*
     * Type of the parent class
     */
    typedef details::AlgebraicCurveFromOrderedPointsBase<F, TInteger, TDistance> Super; 
    /*
     * Type of the (F-1)-class (curve contrained to pass through G+1 given points) 
     */
    typedef AlgebraicCurveFromOrderedPoints<F-1, TInteger, TDistance> Up; 

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
  template <typename TInteger, typename TDistance>
  class AlgebraicCurveFromOrderedPoints<0,TInteger,TDistance>: 
    public details::AlgebraicCurveFromOrderedPointsBase<0, TInteger, TDistance>
  {

  public: 
    /*
     * (1)-class as friend for initFromUp methods
     */
    friend class AlgebraicCurveFromOrderedPoints<1, TInteger>;  

    // ----------------------- Static constant ------------------------------
  public: 
    /*
     * number of free points (= 0)
     */
    static const DGtal::Dimension F = 0;

    // ----------------------- Nested type ------------------------------
  public: 
    /*
     * Self type
     */
    typedef AlgebraicCurveFromOrderedPoints<0, TInteger, TDistance> Self; 
    /*
     * Type of the parent class
     */
    typedef details::AlgebraicCurveFromOrderedPointsBase<0, TInteger, TDistance> Super; 
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
  template <DGtal::Dimension T, typename TInteger, typename TDistance>
  inline
  std::ostream&
  operator<< ( std::ostream & out, 
	       const details::AlgebraicCurveFromOrderedPointsBase<T, TInteger, TDistance> & object )
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
