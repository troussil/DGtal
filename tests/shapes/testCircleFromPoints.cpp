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
 * @file testCircleFromPoints.cpp
 * @author Tristan Roussillon (\c
 * tristan.roussillon@liris.cnrs.fr ) Laboratoire d'InfoRmatique en
 * Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS,
 * France
 *
 *
 * @date 2012/06/29
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <sstream>
#include <exception>
#include <boost/program_options.hpp>

#include "DGtal/base/Common.h"
#include "DGtal/kernel/PointVector.h"
#include "DGtal/shapes/fromPoints/CircleFromPoints.h"

//space / domain
#include "DGtal/kernel/SpaceND.h"
#include "DGtal/kernel/domains/HyperRectDomain.h"
#include "DGtal/topology/KhalimskySpaceND.h"

//shape and digitizer
#include "DGtal/shapes/ShapeFactory.h"
#include "DGtal/shapes/Shapes.h"
#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/shapes/GaussDigitizer.h"
#include "DGtal/geometry/curves/GridCurve.h"

#include "ConfigTest.h"

///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;




//////////////////////////////////////////////////////////////////////////////
// digital circle generator
template<typename TKSpace>
GridCurve<TKSpace>
ballGenerator(double aCx, double aCy, double aR, bool aFlagIsCW)
{

  // Types
  typedef TKSpace KSpace;  
  typedef typename KSpace::SCell SCell;
  typedef GridCurve<KSpace> GridCurve; 
  typedef typename KSpace::Space Space;  
  typedef Ball2D<Space> Shape;
  typedef typename Space::Point Point;
  typedef typename Space::RealPoint RealPoint;
  typedef HyperRectDomain<Space> Domain;

  //Forme
  Shape aShape(Point(aCx,aCy), aR);

  // Window for the estimation
  RealPoint xLow ( -aR-1, -aR-1 );
  RealPoint xUp( aR+1, aR+1 );
  GaussDigitizer<Space,Shape> dig;  
  dig.attach( aShape ); // attaches the shape.
  dig.init( xLow, xUp, 1 ); 
  Domain domain = dig.getDomain();
  // Create cellular space
  KSpace K;
  bool ok = K.init( dig.getLowerBound(), dig.getUpperBound(), true );
  if ( ! ok )
    {
      std::cerr << " "
		<< " error in creating KSpace." << std::endl;
      return GridCurve();
    }
  try 
    {

      // Extracts shape boundary
      SurfelAdjacency<KSpace::dimension> SAdj( true );
      SCell bel = Surfaces<KSpace>::findABel( K, dig, 10000 );
      // Getting the consecutive surfels of the 2D boundary
      std::vector<Point> points, points2;
      Surfaces<KSpace>::track2DBoundaryPoints( points, K, SAdj, dig, bel );
      //counter-clockwise oriented by default
      GridCurve c; 
      if (aFlagIsCW)
	{
	  points2.assign( points.rbegin(), points.rend() );
	  c.initFromVector(points2); 
	} 
      else 
	{
	  c.initFromVector(points); 
	}
      return c;
    }
  catch ( InputException e )
    {
      std::cerr << " "
		<< " error in finding a bel." << std::endl;
      return GridCurve();
    }
}

///////////////////////////////////////////////////////////////////////////////
// Seidel algo
///////////////////////////////////////////////////////////////////////////////

template <typename I, typename S>
bool algo(const I& itb, const I& ite, S& aShape); 

template <typename I, typename P, typename S>
bool algoUpdate(const I& itb, const I& ite, const P& p, S& aShape)
{
  if (S::F < S::N)
    {
      typename S::Up tmp = aShape.getUp( p ); 
      bool res = algo( itb, ite, tmp); 
      //trace.info()  << "tmp     : " << std::endl << tmp << std::endl; 
      aShape.initFromUp( tmp ); 
      //trace.info()  << "new init " << std::endl << aShape << std::endl; 
      return res; 
    }
  else 
    return false; 
}

template <typename I, typename S>
bool algo(const I& itb, const I& ite, S& aShape)
{

  //init
  typedef typename S::Point Point; 
  boost::array<Point,S::V> a;
  unsigned int counter = 0;  
  for (I it = itb; ( (it != ite)&&(counter < S::V) ); ++counter)
    {
      Point p = it->first; 
      while ( (it != ite)&&(p == it->first) ) ++it;  
      a[counter] = p; 
    }
  ASSERT( counter == S::V ); 
  aShape.init( a.begin(), a.end() ); 
  //trace.info() << "Init: " << std::endl << aShape << std::endl; 
  //TODO init pb when if the three first points have the wrong orientation

  //main loop
  bool res = true; 
  for (I it = itb; ( (it != ite)&&(res) ); ++it)
    {

      //trace.info() << " new pair " << it->first << it->second << std::endl; 

      if ( aShape( it->first ) < 0 )
	{//inner point outside
	  //trace.info() << "inner point " << it->first << " at wrong location " << std::endl; 
	  res = algoUpdate( itb, it , it->first, aShape ); 
	}
      if ( (res)&&( aShape( it->second ) > 0 ) )
	{//outer point inside
	  //trace.info() << "outer point " << it->second << " at wrong location " << std::endl; 
	  res = algoUpdate( itb, it, it->second, aShape ); 
	}
    } 
  return res; 

}


/**
 * Recogition of randomly generated digital circles
 */
bool testBallRecognition()
{

  typedef KhalimskySpaceND<2,int> KSpace; 
  GridCurve<KSpace> c; 
  
  unsigned int nbok = 0;
  unsigned int nb = 0;

  trace.beginBlock ( "Recognition" );
  
  for (unsigned int i = 0; i < 1; ++i)
    {
      //generate digital circle
      double cx = (rand()%100 ) / 100.0;
      double cy = (rand()%100 ) / 100.0;
      double radius = (rand()%100 )+0;
      // double cx = 0, cy = 0; 
      // double radius = 5.0; 
      c = ballGenerator<KSpace>( cx, cy, radius, ((i%2)==1) ); 
      trace.info() << " #ball #" << i << " c(" << cx << "," << cy << ") r=" << radius << endl; 
    
      //range
      typedef GridCurve<KSpace>::IncidentPointsRange Range; 
      Range r = c.getIncidentPointsRange();
    
      CircleFromPoints<0> circle; 
      bool flag = algo( r.begin(), r.end(), circle);

      trace.info() << std::endl << "Solution: " << circle << std::endl; 

      //conclusion
      nbok += flag ? 1 : 0; 
      nb++;
    }

  trace.endBlock();
  
  trace.info() << "(" << nbok << "/" << nb << ") " << endl;
  return nbok == nb;
}



///////////////////////////////////////////////////////////////////////////////
// testing Functions
///////////////////////////////////////////////////////////////////////////////

/**
 * Simple tests
 *
 */
template<typename T>
bool fun(const T& t) 
{
  typedef typename T::Point Point; 

  std::cout << t << std::endl; 

  if (T::F < T::N)
    return fun( t.getUp( Point(T::F, T::F*T::F+1) ) ); 
  else {
    typedef typename T::Value Det;
    Det d = t( Point::diagonal(0) );
    std::cout << d << std::endl; 
    return (d < 0); 
  }
}

bool testCircleFromPoints()
{

  trace.beginBlock("Simple test for CircleFromPoints"); 
  
  bool res = fun( CircleFromPoints<0,double>() ); 

  trace.endBlock(); 
  
  return res; 
}


bool testDeterminant()
{

  trace.beginBlock("Determinant test"); 

  CircleFromPoints<0> c;
  typedef CircleFromPoints<0>::Point Point; 
  std::vector<Point> v; 
  v.push_back( Point(0,1) ); 
  v.push_back( Point(150,18) ); 
  v.push_back( Point(250,-48) ); 
  c.init( v.begin(), v.end() ); 
  trace.info() << c << endl;

  typedef CircleFromPoints<0>::Value Det;
  Det d = c( Point(0,0) );   
  trace.info() << Point(0,0) << " is at distance " << d << endl;
  bool res = true; 
  if (d != 4026300) res = false;  

  trace.endBlock(); 
  
  return res; 
}

///////////////////////////////////////////////////////////////////////////////
// Standard services - public :

int main( int argc, char** argv )
{
  trace.beginBlock ( "Testing shapes from points" );
  trace.info() << "Args:";
  for ( int i = 0; i < argc; ++i )
    trace.info() << " " << argv[ i ];
  trace.info() << endl;


  bool res = testCircleFromPoints() && testDeterminant() && testBallRecognition();


  trace.emphase() << ( res ? "Passed." : "Error." ) << endl;
  trace.endBlock();
  
  return res ? 0 : 1;
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
