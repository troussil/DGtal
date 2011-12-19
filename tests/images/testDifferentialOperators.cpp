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

/**
 * @file test_Image.cpp
 * @ingroup Tests
 * @author Tristan Roussillon (\c tristan.roussillon@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 * @date 2011/12/19
 *
 * @brief Implementation of inline methods defined in DifferentialOperatorsOnImages.h
 *
 * This file is part of the DGtal library.
 */



#include <iostream>

#include "DGtal/base/Common.h"
#include "DGtal/kernel/SpaceND.h"
#include "DGtal/kernel/domains/HyperRectDomain.h"

#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/images/DifferentialOperatorsOnImages.h"


//using namespace DGtal;
//using namespace std;

bool firstTest() 
{

  typedef int Label; 
  typedef DGtal::ImageContainerBySTLVector<HyperRectDomain<SpaceND<3,int > >, Label >  Image;
  typedef Image::Point Point; 
  Image img( Point(-5,-5,-5), Point(5,5,5) ); 
  DGtal::DifferentialOperatorsOnImages<Image> helper(img); 
  return true; 

}

int main(int argc, char** argv)
{

  trace.beginBlock ( "Testing class DifferentialOperatorsOnImages" );
  trace.info() << "Args:";
  for ( int i = 0; i < argc; ++i )
    trace.info() << " " << argv[ i ];
  trace.info() << endl;

  bool res = firstTest() 
;
  
  trace.emphase() << ( res ? "Passed." : "Error." ) << endl;
  trace.endBlock();
  return res ? 0 : 1;

}

