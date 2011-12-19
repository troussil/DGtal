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
#include <list>
#include <cstdlib>
#include <cmath>
#include "DGtal/base/Common.h"
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
  template <typename TImage>
  class DifferentialOperatorsOnImages
  {

    //ASSERT on TImage
    BOOST_CONCEPT_ASSERT(( CImageContainer<TImage> )); 

    // ----------------------- Types ------------------------------
  public:

    typedef TImage Image;
    typedef typename Image::Value Value; 
    typedef typename Image::Point Point;
    typedef typename Image::Vector Vector;  
    typedef typename Image::Domain Domain;

    typedef typename Image::iterator Iterator;
    typedef typename Image::const_iterator ConstIterator;

    typedef typename Image::Dimension Dimension;
    static const typename Image::Dimension dimension = Image::dimension;
   


    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Constructor.
     * @param aStartingImage  any image of signed values
     */
    DifferentialOperatorsOnImages( Image& aStartingImage );

    /**
     * Destructor. Does nothing.
     */
    ~DifferentialOperatorsOnImages();

    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    bool isValid() const;


    /**
     * Writes/Displays the object on an output stream.
     * @param out the output stream where the object is written.
     */
    void selfDisplay ( std::ostream & out ) const;


    // ------------------------- Protected Datas ------------------------------
  protected:
    // ------------------------- Private Datas --------------------------------
  private:

    /**
     * Reference on an image
     */
    Image& myU; 

     /**
     * grid step
     */
    double myH; 


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

  }; // end of class DifferentialOperatorsOnImages


  /**
   * Overloads 'operator<<' for displaying objects of class 'DifferentialOperatorsOnImages'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'DifferentialOperatorsOnImages' to write.
   * @return the output stream after the writing.
   */
   template <typename TImage>
  std::ostream&
  operator<< ( std::ostream & out, const DifferentialOperatorsOnImages<TImage> & object );


} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "DifferentialOperatorsOnImages.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined DifferentialOperatorsOnImages_h

#undef DifferentialOperatorsOnImages_RECURSES
#endif // else defined(DifferentialOperatorsOnImages_RECURSES)
