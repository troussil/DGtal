/**
 * @file test_Image.cpp
 * @author David Coeurjolly (\c david.coeurjolly@liris.cnrs.fr )
 *
 *
 * @date 2010/05/25
 *
 * This file is part of the DGtal library
 */

/**
 * Description of test_Image <p>
 * Aim: simple test of \ref Image
 */

#include <cstdio>
#include <cmath>
#include <iostream>

#include "DGtal/base/Common.h"
#include "DGtal/kernel/SpaceND.h"
#include "DGtal/kernel/domains/HyperRectDomain.h"
#include "DGtal/kernel/images/Image.h"
#include "DGtal/kernel/images/ImageContainerBySTLVector.h"
#include "DGtal/kernel/images/ImageContainerBySTLMap.h"

using namespace DGtal;
using namespace std;


/**
* Simple test of Image construction.
*
**/
bool testSimpleImage()
{

    typedef SpaceND<int,4> Space4Type;
    typedef Space4Type::Point Point;
    typedef HyperRectDomain<Space4Type> TDomain;
    typedef ImageContainerBySTLVector<Point, double> TContainer;

    const int t[ ] = { 1, 2, 3 ,4};
    const int t2[ ] = { 5, 5, 3 ,4};
    Point a ( t );
    Point b ( t2 );

    trace.beginBlock ( "Image init" );

    ///Domain characterized by points a and b
    Image<TDomain,double, TContainer > myImage ( a,b );
    trace.info() << myImage << std::endl;
    trace.endBlock();
    return myImage.isValid();
}


bool testImageContainer()
{
    typedef SpaceND<int,4> Space4Type;
    typedef Space4Type::Point Point;
    typedef HyperRectDomain<Space4Type> TDomain;
    typedef ImageContainerBySTLVector<Point, double> TContainerV;
    typedef ImageContainerBySTLMap<Point, double> TContainerM;

    bool res = true;

    const int t[ ] = { 1, 2, 3 ,4};
    const int t2[ ] = { 5, 5, 3 ,4};
    const int t3[ ] = { 5, 3, 3 ,4};
    Point a ( t );
    Point b ( t2 );

    Point c ( t3 );

    trace.beginBlock ( "Image Container" );

    ///Domain characterized by points a and b
    Image<TDomain,double, TContainerV> myImageV ( a,b );
    trace.info() << "Vector container, value at c="<<myImageV( c )<< std::endl;

    ///Domain characterized by points a and b
    Image<TDomain,double, TContainerM> myImageM ( a,b );

    res = res && myImageM.isValid() && myImageV.isValid();

    //We revert the bool flag to catch the exception
    res = !res;

    try {
        trace.info() << "Map container, value at c="<<myImageM( c )<< std::endl;
    }
    catch (std::bad_alloc e)
    {
        trace.warning() << "Exception bad_alloc catched.. this is normal for the map container"<<std::endl;
        res = !res;
    }

    trace.endBlock();

    return res;
}


bool testBuiltInIterators()
{
    typedef SpaceND<int,2> Space2Type;
    typedef Space2Type::Point Point;
    typedef HyperRectDomain<Space2Type> TDomain;
    typedef ImageContainerBySTLVector<Point, double> TContainerV;
    typedef ImageContainerBySTLMap<Point, double> TContainerM;

    const int t[ ] = { 1, 1};
    const int t2[ ] = { 5, 5};
    Point a ( t );
    Point b ( t2 );

    trace.beginBlock("Test of built-in iterators");

    Image<TDomain,double, TContainerV> myImageV ( a,b );
    Image<TDomain,double, TContainerM> myImageM ( a,b );

    trace.info() << " Vector iterator: ";
    for ( Image<TDomain,double, TContainerV>::Iterator it = myImageV.begin();
            it != myImageV.end();
            ++it)
        trace.info() << myImageV(it) <<" ";

    trace.info() << std::endl<<" Map iterator: ";


    for ( Image<TDomain,double, TContainerM>::Iterator it = myImageM.begin();
            it != myImageM.end();
            ++it)
        trace.info() << myImageM(it)<<" ";

    trace.info()<<std::endl;
    trace.endBlock();

    return true;
}

template <typename T>
struct Hop
{
    std::vector<double> t;
    typedef std::vector<double>::iterator Iterator;
		typedef std::vector<double>::iterator SpanIterator;
		typedef std::vector<double>::iterator ConstSpanIterator;
		typedef std::vector<double>::const_iterator ConstIterator;
    Hop(const T &a, const T &b) {};
    void  allocate(std::size_t hlop) {};
    Iterator begin() {
        return t.begin();
    };
};

bool testConcepts()
{
    typedef SpaceND<int,2> Space2Type;
    typedef Space2Type::Point Point;
    typedef HyperRectDomain<Space2Type> TDomain;
    typedef ImageContainerBySTLVector<Point, double> TContainerV;
    typedef ImageContainerBySTLMap<Point, double> TContainerM;

    const int t[ ] = { 1, 1};
    const int t2[ ] = { 5, 5};
    Point a ( t );
    Point b ( t2 );

    trace.beginBlock("Test of Concepts");

    //Image<TDomain,double, Hop<Point> > myImageV ( a,b );
    trace.endBlock();

    return true;
}




int main()
{

  if ( testSimpleImage() && testImageContainer() && testBuiltInIterators() && testConcepts() )
        return 0;
    else
        return 1;
}

