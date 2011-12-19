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


bool first2dTest() 
{

  //small img
  typedef DGtal::ImageContainerBySTLVector<HyperRectDomain<SpaceND<2,int > >, int >  Image;
  typedef Image::Point Point; 
  typedef Image::Vector Vector; 
  Image img( Point(-1,-1), Point(1,1) ); 
  DGtal::DifferentialOperatorsOnImages<Image> helper(img); 

  for (int i = 0; i < 3; ++i)
    {
      img.setValue( Point(i-1,-1), i*i ); 
      img.setValue( Point(i-1,0), i*i ); 
      img.setValue( Point(i-1,1), i*i ); 
    }

  Point c(0,0); 
  Point a(-1,0); 
  Point d(1,0); 

  //tests
  unsigned int nbok = 0;
  unsigned int nb = 0;  

  trace.beginBlock ( "simple test" );

  {
    bool flag = ( helper.forwardDifference(c,0) == 3 ) 
    && ( helper.forwardDifference(c,1) == 0 ) 
    && ( helper.forwardDifference(a,0) == 1 ) 
    && ( helper.forwardDifference(d,0) == 3 ); 
    trace.info() << ((flag)? "Passed": "Failed") << std::endl; 
    nbok += flag ? 1 : 0; 
    nb++;
  }

  {
    bool flag = ( helper.backwardDifference(c,0) == 1 ) 
    && ( helper.backwardDifference(c,1) == 0 ) 
    && ( helper.backwardDifference(a,0) == 1 ) 
    && ( helper.backwardDifference(d,0) == 3 ); 
    trace.info() << ((flag)? "Passed": "Failed") << std::endl; 
    nbok += flag ? 1 : 0; 
    nb++;
  }

  {
    bool flag = false; 
    flag = ( helper.centralDifference(c,0) == 2 ) 
    && ( helper.centralDifference(c,1) == 0 ) 
    && ( helper.centralDifference(a,0) == 1 ) 
    && ( helper.centralDifference(d,0) == 3 ); 
    trace.info() << ((flag)? "Passed": "Failed") << std::endl; 
    nbok += flag ? 1 : 0; 
    nb++;
  }

  {
    bool flag = ( helper.centralGradient(c) == Point(2,0) ) 
    && ( helper.centralGradientModulus(a) == 1 ); 
    trace.info() << ((flag)? "Passed": "Failed") << std::endl; 
    nbok += flag ? 1 : 0; 
    nb++;
  }

  {
    bool flag = ( helper.upwindGradientModulus(c,Vector(1,1)) == 1 )
      && ( helper.upwindGradient(c,Vector(-1,-1)) == Point(3,0) ); 
    trace.info() << ((flag)? "Passed": "Failed") << std::endl; 
    nbok += flag ? 1 : 0; 
    nb++;
  }

  {
    bool flag = ( helper.godunovGradientModulus(c,true) == 1 )
      && ( helper.godunovGradientModulus(c,false) == 3 ); 
    trace.info() << ((flag)? "Passed": "Failed") << std::endl; 
    nbok += flag ? 1 : 0; 
    nb++;
  }

  trace.info() << "(" << nbok << "/" << nb << ") " << std::endl;
  trace.endBlock();
  return (nbok == nb);
}

int main(int argc, char** argv)
{

  trace.beginBlock ( "Testing class DifferentialOperatorsOnImages" );
  trace.info() << "Args:";
  for ( int i = 0; i < argc; ++i )
    trace.info() << " " << argv[ i ];
  trace.info() << endl;

  bool res = first2dTest() 
 ;
  
  trace.emphase() << ( res ? "Passed." : "Error." ) << endl;
  trace.endBlock();
  return res ? 0 : 1;

}

