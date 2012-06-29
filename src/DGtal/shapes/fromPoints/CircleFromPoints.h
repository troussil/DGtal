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
 * @file CircleFromPoints.h
 * @brief Representation of a CircleFromPoints uniquely defined by two 2D points.
 * @author Tristan Roussillon (\c tristan.roussillon@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 * @date 2011/09/22
 *
 * @brief Header file for module CircleFromPoints.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(CircleFromPoints_RECURSES)
#error Recursive header files inclusion detected in CircleFromPoints.h
#else // defined(CircleFromPoints_RECURSES)
/** Prevents recursive inclusion of headers. */
#define CircleFromPoints_RECURSES

#if !defined CircleFromPoints_h
/** Prevents repeated inclusion of headers. */
#define CircleFromPoints_h

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

  /////////////////////////////////////////////////////////////////////////////
  // template class CircleFromPointsBase
  /**
   * \brief Aim: Represents a circle 
   *
   */
  template <DGtal::Dimension F, typename TInteger = int>
  struct CircleFromPointsBase
  {

    BOOST_STATIC_ASSERT((F >= 0 ));
    BOOST_STATIC_ASSERT((F <= 3 ));

	//total number of points
	//(defining uniquely a circle)
    static const DGtal::Dimension N = 3; 
 	//number of variable points 
    static const DGtal::Dimension V = N-F; 

    // ----------------------- associated types ------------------------------
    typedef TInteger Coordinate;
    typedef DGtal::PointVector<2,Coordinate> Point; 
    typedef DGtal::PointVector<3,Coordinate> HPoint; //in homogeneous coordinates 

    typedef Coordinate Value; //TODO: to promote

    typedef boost::array<HPoint,F> FArray;  //array of fixed points
    typedef boost::array<HPoint,V> VArray;  //array of variable points
    typedef boost::array<HPoint,N> Array;  //array of all points

    // ------------------------- Private Datas --------------------------------
  protected:
    /**
       Array of fixed points
    */
    FArray myFArray;
    /**
       Array of variable points
    */
    VArray myVArray;

    // ----------------------- Standard services ------------------------------
  public:


    /**
     * Constructor when no points are fixed.
     */
    CircleFromPointsBase() { BOOST_STATIC_ASSERT(( F == 0 )); };

    /**
     * Constructor.
     * @param aFArray an array of fixed points
     */
    CircleFromPointsBase(const FArray& aFArray): myFArray( aFArray ) {};


    /**
     * Init.
     * @param itb begin iterator on points
     * @param ite end iterator on points
     */
    template <typename I>
    void init(const I& itb, const I& ite) 
      { 
	typename VArray::iterator ait = myVArray.begin(); 
	for (I it = itb; ( (it != ite)&&(ait != myVArray.end()) ); ++it, ++ait)
            *ait = toHPoint( *it );  
      };

    /**
     * Init.
     * @param aVArray an array of variable points
     */
    void init(const VArray& aVArray) { myVArray = aVArray; };

    /**
     * Copy constructor.
     * @param other the object to clone.
     */
	template <typename Other>
    CircleFromPointsBase ( const Other & other ): myFArray( other.myFArray ), myVArray( other.myVArray ) {};

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     */
    CircleFromPointsBase & operator= ( const CircleFromPointsBase & other )
      {
	if (this != &other)
          {
            myFArray = other.myFArray; 
            myVArray = other.myVArray; 
          }
	return *this; 
      };

    /**
     * Destructor. Does nothing
     */
    ~CircleFromPointsBase() {};

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
	std::copy( myVArray.begin(), myVArray.end(), tmp.begin() ); 
	std::copy( myFArray.begin(), myFArray.end(), tmp.begin()+V ); 
	return determinant( tmp, toHPoint( aP ) ); 
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
         std::cout << "[CircleFromPointsBase] "
	           << F << " fixed / " << N << " points. "
                   << std::endl; 
	std::copy( myVArray.begin(), myVArray.end(), std::ostream_iterator<HPoint>(std::cout, ",") ); 
	std::copy( myFArray.begin(), myFArray.end(), std::ostream_iterator<HPoint>(std::cout, ",") ); 
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
    Value determinant(const Array& a, const HPoint& aP) const 
     { 

	//std::copy( a.begin(), a.end(), std::ostream_iterator<HPoint>(std::cout, ",") ); 
	//std::cout << aP << " " << a[2][1] << std::endl; 
	//TODO delegate this task. 

 	typedef DGtal::PointVector<2,Coordinate> Vector; 

	Vector u( (a[0][0]-aP[0])*(a[1][1]-aP[1])-  (a[1][0]-aP[0])*(a[0][1]-aP[1]), 
		  (a[1][0]-aP[0])*(a[1][0]-a[0][0])+(a[1][1]-aP[1])*(a[1][1]-a[0][1]) );
	Vector v( (a[0][0]-aP[0])*(a[2][1]-aP[1])-  (a[2][0]-aP[0])*(a[0][1]-aP[1]), 
		  (a[2][0]-aP[0])*(a[2][0]-a[0][0])+(a[2][1]-aP[1])*(a[2][1]-a[0][1]) );
	//std::cout << u << v << std::endl; 
	return -( (u[0] * v[1]) - (u[1] * v[0]) ); 
     };

  }; // end of class CircleFromPointsBase

  /////////////////////////////////////////////////////////////////////////////
  // template class CircleFromPoints
  /**
   * \brief Aim: Represents a circle 
   *
   */
  template <DGtal::Dimension T, typename TInteger = int>
  class CircleFromPoints: public CircleFromPointsBase<T,TInteger>
  {

  public: 
    //for down method
    friend class CircleFromPoints<T-1, TInteger>;  

    // ----------------------- Static constant ------------------------------
  public: 
	//number of fixed points
    static const DGtal::Dimension F = T;

    // ----------------------- Nested type ------------------------------
  public: 
    typedef CircleFromPoints<T,TInteger> Self; 
    typedef CircleFromPointsBase<T,TInteger> Super; 
    typedef CircleFromPoints<F+1, typename Self::Coordinate> Up; 

    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Constructor when no points are fixed.
     */
    CircleFromPoints() { BOOST_STATIC_ASSERT(( F == 0 )); };

    /**
     * Constructor.
     * @param aFArray an array of fixed points
     */
    CircleFromPoints(const typename Self::FArray& aFArray): Super( aFArray ) {};

    /**
     * Copy constructor.
     * @param other the object to clone.
     */
    CircleFromPoints ( const Self & other ): Super( other ) {};

    /**
     * Init from an instance of Up
     * @param u an instance of Up
     */
    void initFromUp(const Up& u) 
      {
	/// variable points: 
	//all variable points of u
	std::copy( u.myVArray.begin(), u.myVArray.end(), this->myVArray.begin() ); 
	//+ the first (i.e. last added) fixed point of u
	ASSERT( Up::F > 0 ); 
	this->myVArray.back() = u.myFArray.front();
	/// fixed points: 
	std::copy( u.myFArray.begin()+1, u.myFArray.end(), this->myFArray.begin() ); 
      }

    /**
     * Set
     * @param aPoint a point 
     * @return an instance of Up ( @e aHPoint is taken as an extra fixed point)
     */
    Up set( const typename Self::Point& aPoint ) const
      {
	return set( toHPoint( aPoint ) ); 
      }

    /**
     * Set
     * @param aHPoint a point (in homogeneous coordinates)
     * @return an instance of Up ( @e aHPoint is taken as an extra fixed point)
     */
    Up set( const typename Self::HPoint& aHPoint ) const
      {
	boost::array<typename Self::HPoint,F+1> tmp; 
	tmp.front() = aHPoint;  
	std::copy( this->myFArray.begin(), this->myFArray.end(), tmp.begin()+1 ); 
	return Up( tmp ); 
      }  

    /**
     * Destructor. Does nothing
     */
    ~CircleFromPoints() {};

   private: 

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     */
    Self & operator= ( const Self & other );

  }; // end of class CircleFromPoints

  /////////////////////////////////////////////////////////////////////////////
  // Specialization
  template <typename TInteger>
  class CircleFromPoints<3,TInteger>: public CircleFromPointsBase<3,TInteger>
  {

  public: 
    //for down method
    friend class CircleFromPoints<2, TInteger>;  

    // ----------------------- Static constant ------------------------------
  public: 
	//number of fixed points
    static const DGtal::Dimension F = 3;

    // ----------------------- Nested type ------------------------------
  public: 
    typedef CircleFromPoints<3,TInteger> Self; 
    typedef CircleFromPointsBase<3,TInteger> Super; 
    typedef Self Up;   

    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Constructor.
     * @param aFArray an array of fixed points
     */
    CircleFromPoints(const typename Self::FArray& aFArray): Super( aFArray ) {};


    /**
     * Copy constructor.
     * @param other the object to clone.
     */
    CircleFromPoints ( const Self & other ): Super( other ) {};

    /**
     * Destructor. Does nothing
     */
    ~CircleFromPoints() {};

    /**
     * Does nothing
     */
    void initFromUp( const Up& u ) const
      {
      }

    /**
     * Set
     * @param aPoint a point 
     * @return an instance of Up
     */
    Up set( const typename Self::Point& /*aPoint */) const
      {
	*this; 
      }

    /**
     * Set
     * @param aHPoint a point (in homogeneous coordinates)
     * @return an instance of Up
     */
    Up set( const typename Self::HPoint& /*aHPoint*/ ) const
      {
	return *this; 
      }  

    /**
     * @return the style name used for drawing this object.
     */
    std::string className() const { return "CircleFromPoints"; };

   private: 


    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     */
    Self & operator= ( const Self & other );

  }; // end of class CircleFromPoints

  /**
   * Overloads 'operator<<' for displaying objects of class 'CircleFromPoints'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'CircleFromPoints' to write.
   * @return the output stream after the writing.
   */
  template <DGtal::Dimension T, typename TInteger>
  inline
  std::ostream&
  operator<< ( std::ostream & out, 
        const CircleFromPointsBase<T, TInteger> & object )
  {
    object.selfDisplay( out );
    return out;
  }


} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
//#include "DGtal/shapes/fromPoints/CircleFromPoints.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined CircleFromPoints_h

#undef CircleFromPoints_RECURSES
#endif // else defined(CircleFromPoints_RECURSES)
