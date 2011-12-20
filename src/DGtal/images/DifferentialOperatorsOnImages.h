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
 * @file DifferentialOperatorsOnImages.h
 * @author Tristan Roussillon (\c tristan.roussillon@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 * @date 201/12/19
 *
 * Header file for module DifferentialOperatorsOnImages.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(DifferentialOperatorsOnImages_RECURSES)
#error Recursive header files inclusion detected in DifferentialOperatorsOnImages.h
#else // defined(DifferentialOperatorsOnImages_RECURSES)
/** Prevents recursive inclusion of headers. */
#define DifferentialOperatorsOnImages_RECURSES

#if !defined DifferentialOperatorsOnImages_h
/** Prevents repeated inclusion of headers. */
#define DifferentialOperatorsOnImages_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include "DGtal/images/CImageContainer.h"

//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // template class DifferentialOperatorsOnImages
  /**
   * Description of template class 'DifferentialOperatorsOnImages' <p>
   * \brief Aim: Computes some differential operators on an image. 
   *
   * @code 
   * @endcode
   *
   * @tparam TImage type of image 
   */
  template <typename TImage, typename TOutputValue = typename TImage::Value >
  class DifferentialOperatorsOnImages
  {

    //ASSERT on TImage
    BOOST_CONCEPT_ASSERT(( CImageContainer<TImage> )); 

    // ----------------------- Types ------------------------------
  public:

    typedef TImage Image;
    typedef TOutputValue OutputValue; 
    typedef typename Image::Value Value;
 
    typedef typename Image::Point Point;
    typedef typename Image::Vector Vector;  
    typedef typename Image::Domain Domain;

    typedef typename Image::Dimension Dimension;
    static const typename Image::Dimension dimension = Image::dimension;
   
    typedef DGtal::PointVector<dimension, OutputValue> Gradient; 

    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Constructor.
     * @param aStartingImage  any image of signed values
     */
    DifferentialOperatorsOnImages( Image& aStartingImage );

    /**
     * Constructor.
     *
     * @param aStartingImage  any image of signed values
     * @param aGridStep any length 
     */
    DifferentialOperatorsOnImages( Image& aStartingImage, const OutputValue& aGridStep );

    /**
     * Destructor. Does nothing.
     */
    ~DifferentialOperatorsOnImages();

    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    bool isValid() const;

    // ----------------------- First order ------------------------------
    /**
     * Forward difference.
     *
     * @param aPoint the point where the derivative is computed
     * @param aDim the axis along which the derivative is computed
     * @return first derivative along axis @a aDim at @ aPoint
     */
    OutputValue forwardDifference ( const Point& aPoint, const Dimension& aDim ) const; 

    /**
     * Backward difference.
     *
     * @param aPoint the point where the derivative is computed
     * @param aDim the axis along which the derivative is computed
     * @return first derivative along axis @a aDim at @ aPoint
     */
    OutputValue backwardDifference ( const Point& aPoint, const Dimension& aDim ) const; 

    /**
     * central difference.
     *
     * @param aPoint the point where the derivative is computed
     * @param aDim the axis along which the derivative is computed
     * @return first derivative along axis @a aDim at @ aPoint
     */
    OutputValue centralDifference ( const Point& aPoint, const Dimension& aDim ) const; 

    /**
     * Gradient
     *
     * @param aPoint the point where the gradient is computed
     * @return the gradient of @a myU at @a aPoint
     */
    Gradient centralGradient ( const Point& aPoint ) const; 

    /**
     * Gradient modulus
     *
     * @param aPoint the point where the gradient is computed
     * @return modulus of the gradient of @a myU at @a aPoint
     */
    double centralGradientModulus ( const Point& aPoint ) const; 

    /**
     * Gradient computed from the upwind scheme
     *
     * @param aPoint the point where the gradient is computed
     * @param aVector displacement vector 
     * @return the gradient of @a myU at @a aPoint
     */
    Gradient upwindGradient ( const Point& aPoint, const Vector& aVector ) const; 

    /**
     * Modulus of the gradient computed from the upwind scheme
     *
     * @param aPoint the point where the gradient is computed
     * @param aVector displacement vector 
     * @return modulus of the gradient of @a myU at @a aPoint
     */
    double upwindGradientModulus ( const Point& aPoint, const Vector& aVector ) const; 

    /**
     * Gradient computed from the Godunov scheme
     *
     * @param aPoint the point where the gradient is computed
     * @param isPositive boolean equal to 'true' if the displacement vector
     * has the same orientation as the normal to the interface
     * @return the gradient of @a myU at @a aPoint
     */
    DGtal::PointVector<TImage::dimension,TOutputValue>
    godunovGradient ( const Point& aPoint, bool isPositive ) const; 

    /**
     * Gradient modulus computed from the Godunov scheme
     *
     * @param aPoint the point where the gradient is computed
     * @param isPositive boolean equal to 'true' if the displacement vector
     * has the same orientation as the normal to the interface
     * @return modulus of the gradient of @a myU at @a aPoint
     */
    double godunovGradientModulus ( const Point& aPoint, bool isPositive ) const; 

    // ----------------------- second order ------------------------------
    /**
     * Forward/backward second order difference.
     *
     * @param aPoint the point where the derivative is computed
     * @param aDim the axis along which the derivative is computed
     * @return second derivative along axis @a aDim at @ aPoint
     */
    OutputValue forwardBackwardDifference2 ( const Point& aPoint, const Dimension& aDim ) const; 

    /**
     * Laplacian.
     *
     * @param aPoint the point where the laplacian is computed
     * @return laplacian of @myU at @ aPoint
     */
    OutputValue laplacian ( const Point& aPoint ) const; 
    

    /**
     * Weighted forward/backward second order difference.
     * The weights are the values of @w , divided by 
     * the gradients modulii.
     *
     * @param w image of weights
     * @param aPoint the point where the derivative is computed
     * @param aDim the axis along which the derivative is computed
     * @return second derivative along axis @a aDim at @ aPoint
     */
    OutputValue weightedDifference2 ( Image& w, const Point& aPoint, const Dimension& aDim ) const; 

    /**
     * Normalized forward/backward second order difference 
     * by the gradients modulii.
     *
     * @param aPoint the point where the derivative is computed
     * @param aDim the axis along which the derivative is computed
     * @return second derivative along axis @a aDim at @ aPoint
     */
    OutputValue normalizedDifference2 ( const Point& aPoint, const Dimension& aDim ) const; 

    /**
     * Weighted mean curvature
     *
     * @param aImg image of weights
     * @param aPoint the point where the curvature is computed
     * @return weighted mean curvature of @myU at @ aPoint
     */
    OutputValue weightedMeanCurvature ( Image& aImg, const Point& aPoint ) const; 

    /**
     * Mean curvature
     *
     * @param aPoint the point where the curvature is computed
     * @return mean curvature of @myU at @ aPoint
     */
    OutputValue meanCurvature ( const Point& aPoint ) const; 

    

    // ------------------------- Protected Datas ------------------------------
  protected:
    // ------------------------- Private Datas --------------------------------
  private:

    /**
     * Reference on an image
     */
    Image& myU; 

    /**
     * Grid step
     */
    OutputValue myH; 

    // ------------------------- Hidden services ------------------------------
  protected:

    /**
     * Constructor.
     * Forbidden by default (protected to avoid g++ warnings).
     */
    DifferentialOperatorsOnImages();

    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    DifferentialOperatorsOnImages ( const DifferentialOperatorsOnImages & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    DifferentialOperatorsOnImages & operator= ( const DifferentialOperatorsOnImages & other );

  private:


    // ------------------------- Internals ------------------------------------
  private:

    /**
     * Return the point following @a aPoint along axis @a aDim.
     *
     * @param aPoint any point
     * @param aDim the axis along which the @a aPoint is shifted
     * @return point following @a aPoint along axis @a aDim
     */
    Point getNext ( const Point& aPoint, const Dimension& aDim ) const; 

    /**
     * Return the point preceeding @a aPoint along axis @a aDim.
     *
     * @param aPoint any point
     * @param aDim the axis along which the @a aPoint is shifted
     * @return point preceeding @a aPoint along axis @a aDim
     */
    Point getPrevious ( const Point& aPoint, const Dimension& aDim ) const; 

    /**
     * Return the harmonic average of @a aV1 and @a aV2, 
     * @a aV1 and @a aV2 being normalized by @a aN1 and @a aN2.
     *
     * @param aV1 a first value
     * @param aN1 any value dividing @a aV1
     * @param aV2 a second value
     * @param aN2 any value dividing @a aV2
     * @return the normalized harmonic average of @a aV1 and @a aV2
     */
    double normalizedAverage ( const Value& aV1, const double& aN1, 
			       const Value& aV2, const double& aN2 ) const; 



  }; // end of class DifferentialOperatorsOnImages




} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "DifferentialOperatorsOnImages.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined DifferentialOperatorsOnImages_h

#undef DifferentialOperatorsOnImages_RECURSES
#endif // else defined(DifferentialOperatorsOnImages_RECURSES)
