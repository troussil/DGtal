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

  trace.beginBlock ( "unit tests" );

  {
    DGtal::ForwardDifference<Image> fd(img); 

    bool flag = ( fd(c,0) == 3 ) 
    && ( fd(c,1) == 0 ) 
    && ( fd(a,0) == 1 ) 
    && ( fd(d,0) == 3 ); 
    trace.info() << ((flag)? "Passed": "Failed") << std::endl; 
    nbok += flag ? 1 : 0; 
    nb++;
  }

  {
    DGtal::BackwardDifference<Image> bd(img); 

    bool flag = ( bd(c,0) == 1 ) 
    && ( bd(c,1) == 0 ) 
    && ( bd(a,0) == 1 ) 
    && ( bd(d,0) == 3 ); 
    trace.info() << ((flag)? "Passed": "Failed") << std::endl; 
    nbok += flag ? 1 : 0; 
    nb++;
  }

  {
    DGtal::CentralDifference<Image> cd(img); 

    bool flag = ( cd(c,0) == 2 ) 
    && ( cd(c,1) == 0 ) 
    && ( cd(a,0) == 1 ) 
    && ( cd(d,0) == 3 ); 
    trace.info() << ((flag)? "Passed": "Failed") << std::endl; 
    nbok += flag ? 1 : 0; 
    nb++;
  }

  {
    DGtal::Gradient<DGtal::CentralDifference<Image> > g(img); 
    DGtal::GradientModulus<DGtal::Gradient<DGtal::CentralDifference<Image> > > m(img);  
    bool flag = ( g(c) == Point(2,0) ) 
    && ( m(a) == 1 ) 
    && ( m(c) == 2 )  
    && ( m(d) == 3 ); 
    trace.info() << ((flag)? "Passed": "Failed") << std::endl; 
    nbok += flag ? 1 : 0; 
    nb++;
  }

  {
    DGtal::UpwindGradient<Image> g(img); 
    bool flag = ( g(c,Vector(1,1)) == Point(1,0) )
      && ( g(c,Vector(-1,-1)) == Point(3,0) ); 
    trace.info() << ((flag)? "Passed": "Failed") << std::endl; 
    nbok += flag ? 1 : 0; 
    nb++;
  }

  {
    DGtal::GodunovGradient<Image> g(img); 
    DGtal::GradientModulus<DGtal::GodunovGradient<Image> > m(g); 
    DGtal::GodunovGradient<Image> g2(img,false); 
    DGtal::GradientModulus<DGtal::GodunovGradient<Image> > m2(g2); 

    bool flag = ( m(c) == 1 )
      && ( m2(c) == 3 ); 
    trace.info() << ((flag)? "Passed": "Failed") << std::endl; 
    nbok += flag ? 1 : 0; 
    nb++;
  }

  {
    DGtal::Difference2<Image> diff(img); 
    bool flag = ( diff(c,0) == 2 )
      && ( diff(a,0) == 0 ) 
      && ( diff(d,0) == 0 );
    trace.info() << ((flag)? "Passed": "Failed") << std::endl; 
    nbok += flag ? 1 : 0; 
    nb++;
  }

  {
    //weights
    Image w( Point(-1,-1), Point(1,1) ); 

    for (int i = 0; i < 3; ++i)
      {
	w.setValue( Point(i-1,-1), i+1 ); 
	w.setValue( Point(i-1,0), i+1 ); 
	w.setValue( Point(i-1,1), i+1 ); 
      }

    DGtal::WeightedDifference2<Image> wd(img,w); 
    DGtal::NormalizedDifference2<Image> nd(img); 
    DGtal::NormalizedDifference2<Image,double> nd2(img); 
    DGtal::Divergence<DGtal::WeightedDifference2<Image> > wmc(wd); 
    DGtal::Divergence<DGtal::NormalizedDifference2<Image,double> > mc(nd2); 
    DGtal::Divergence<DGtal::Difference2<Image> > laplacian(img); 

    bool flag = ( wd(c,0) == 2 )
      && ( nd(c,0) == 0 )
      && ( ( nd2(c,0) - (8/(double)15) ) < 0.0001 )
      && ( ( wmc(c) - 2 ) < 0.0001 )
      && ( ( mc(c) - (8/(double)15) ) < 0.0001 )
      && ( laplacian(c) == 2 )
      ;
    trace.info() << ((flag)? "Passed": "Failed") << std::endl; 
    nbok += flag ? 1 : 0; 
    nb++;
  }


  trace.info() << "(" << nbok << "/" << nb << ") " << std::endl;
  trace.endBlock();
  return (nbok == nb);
}


void applyOperatorsTest() 
{

  //small img
  typedef DGtal::ImageContainerBySTLVector<HyperRectDomain<SpaceND<2,int > >, int >  Image;
  typedef Image::Point Point; 
  Image img( Point(-1,-1), Point(1,1) ); 
  for (int i = 0; i < 3; ++i)
    {
      img.setValue( Point(i-1,-1), i*i ); 
      img.setValue( Point(i-1,0), i*i ); 
      img.setValue( Point(i-1,1), i*i ); 
    }
  
  trace.beginBlock ( "Apply operators" );

  //content of img
  std::copy(img.begin(), img.end(), std::ostream_iterator<int>(std::cout, ", ") ); 
  std::cout << std::endl; 

  //img2
  Image img2( Point(-1,-1), Point(1,1) ); 
  //operator 
  typedef DGtal::Divergence<DGtal::Difference2<Image> > Laplacian; 
  //fill img2 with the laplacian of img
  Image::Domain d = img.domain(); 
  std::transform(d.begin(),d.end(), img2.begin(), Laplacian(img) ); 

  //content of img2
  std::copy(img2.begin(), img2.end(), std::ostream_iterator<int>(std::cout, ", ") ); 
  std::cout << std::endl; 

  trace.endBlock();

}

int main(int argc, char** argv)
{

  trace.beginBlock ( "Testing class DifferentialOperatorsOnImages" );
  trace.info() << "Args:";
  for ( int i = 0; i < argc; ++i )
    trace.info() << " " << argv[ i ];
  trace.info() << endl;

  applyOperatorsTest(); 
  bool res = first2dTest() 
 ;
  
  trace.emphase() << ( res ? "Passed." : "Error." ) << endl;
  trace.endBlock();
  return res ? 0 : 1;

}

