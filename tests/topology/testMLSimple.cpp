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
 * @file testMLSimple.cpp
 * @ingroup Tests
 * @author Tristan Roussillon (\c tristan.roussillon@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en LabelsImage et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 * @date 2011/12/19
 *
 * @brief Functions for testing class SimplePointHelper.
 *
 * This file is part of the DGtal library.
 */

#include "DGtal/topology/helpers/SimplePointHelper.h"
#include <iostream>




///////////////// main functions
template<DGtal::Dimension dimension, typename TIterator>
bool basicTest(const TIterator& itb, const TIterator& ite, const typename std::iterator_traits<TIterator>::value_type& aLabel, bool flag)
{
  typedef typename std::iterator_traits<TIterator>::value_type Label;
  typedef HyperRectDomain<SpaceND<dimension,int > > Domain;  
  typedef typename Domain::Point Point; 
  typedef DGtal::ImageContainerBySTLVector<Domain, Label >  Image;
  Image img( Domain( Point::diagonal(-1), Point::diagonal(1) ) ); 
  DGtal::SimplePointHelper<Image>::readConfiguration(img, itb, ite); 
  DGtal::SimplePointHelper<Image> h( img); 

  return ( h.isMLSimple( Point::diagonal(0), aLabel) == flag ); 
}


bool firstTest()
{

  trace.beginBlock ( "simple test" );
  std::string s = "aaaaabbbb";

  typedef std::iterator_traits<std::string::iterator>::value_type Label; 
  typedef HyperRectDomain<SpaceND<2,int > > Domain;  
  typedef Domain::Point Point; 
  typedef DGtal::ImageContainerBySTLVector<Domain, Label >  Image;
  Image img( Domain( Point::diagonal(-1), Point::diagonal(1) ) ); 

  DGtal::SimplePointHelper<Image>::readConfiguration(img, s.begin(), s.end()); 
  DGtal::SimplePointHelper<Image> h(img); 

  bool res = ( ( h.isValid() ) && ( h.isMLSimple( Point::diagonal(0),'b' ) == true ) );
  trace.info() << ((res)? "Passed": "Failed") << std::endl; 
  trace.endBlock();
  return (res == true);

}

bool known2dCasesTest()
{

  trace.beginBlock ( "known 2d Cases" );
  
  unsigned int nbok = 0;
  unsigned int nb = 0;  

  { //accepted config with only X and R (PRL2011 fig1a)
    unsigned int t[9] = {
      2,1,1,2,1,1,2,2,2 };

    bool res = (basicTest<2,unsigned int*>(t+0,t+9,2,true));
    trace.info() << ((res)? "Passed": "Failed") << std::endl; 
    nbok += res ? 1 : 0; 
    nb++;
  }

  { //accepted config with one region O distinct from X and R (PRL2011 fig1b)
    unsigned int t[9] = {
      3,1,1,3,1,1,2,2,2 };
  
    bool res = (basicTest<2,unsigned int*>(t+0,t+9,2,true));
    trace.info() << ((res)? "Passed": "Failed") << std::endl; 
    nbok += res ? 1 : 0; 
    nb++;
  }

  { //accepted config with one region O distinct from X and R
    unsigned int t[9] = {
      3,2,2,3,1,2,3,1,1 };
  
    bool res = (basicTest<2,unsigned int*>(t+0,t+9,2,true));
    trace.info() << ((res)? "Passed": "Failed") << std::endl; 
    nbok += res ? 1 : 0; 
    nb++;
  }
  
  { //rejected config because one region O, twice adjacent to X, 
    //is only once adjacent to X after the flip of x
    unsigned int t[9] = {
      4,3,2,3,1,2,3,1,1 };
  
    bool res = (basicTest<2,unsigned int*>(t+0,t+9,2,false));
    trace.info() << ((res)? "Passed": "Failed") << std::endl; 
    nbok += res ? 1 : 0; 
    nb++;
  }

  
  { //rejected config by cond 2 (PRL2011 fig2b)
    unsigned int t[9] = {
      1,3,3,1,1,1,4,2,2 };
  
    bool res = (basicTest<2,unsigned int*>(t+0,t+9,2,false));
    trace.info() << ((res)? "Passed": "Failed") << std::endl; 
    nbok += res ? 1 : 0; 
    nb++;
  }

  { //rejected config by cond 2 (PRL2011 fig2d)
    unsigned int t[9] = {
      3,1,4,3,1,1,2,2,2 };
  
    bool res = (basicTest<2,unsigned int*>(t+0,t+9,2,false));
    trace.info() << ((res)? "Passed": "Failed") << std::endl; 
    nbok += res ? 1 : 0; 
    nb++;
  }

  { //rejected config by cond 3 (PRL2011 fig2c)
    unsigned int t[9] = {
      3,1,1,3,1,1,4,2,2 };
  
    bool res = (basicTest<2,unsigned int*>(t+0,t+9,2,false));
    trace.info() << ((res)? "Passed": "Failed") << std::endl; 
    nbok += res ? 1 : 0; 
    nb++;
  }


  trace.info() << "(" << nbok << "/" << nb << ") " << endl;  
  trace.endBlock();
  return (nbok == nb);
}


bool known3dCasesTest()
{

  trace.beginBlock ( "known 3d Cases" );
  
  unsigned int nbok = 0;
  unsigned int nb = 0;  

  { //accepted configuration because there is no flip
    unsigned int t[27] = {
    1, 2, 2, 1, 1, 2, 1, 1, 1, 
    1, 2, 2, 1, 1, 2, 1, 1, 1, 
    1, 2, 2, 1, 1, 2, 1, 1, 1  };

    bool res = (basicTest<3,unsigned int*>(t+0,t+27,1,true));
    trace.info() << ((res)? "Passed": "Failed") << std::endl; 
    nbok += res ? 1 : 0; 
    nb++;
  }

  { //rejected configuration because a hole appears
    unsigned int t[27] = {
    1, 1, 1, 1, 1, 1, 1, 1, 1, 
    1, 1, 1, 1, 1, 1, 1, 1, 1, 
    1, 1, 1, 1, 1, 1, 1, 1, 1  };

    bool res = (basicTest<3,unsigned int*>(t+0,t+27,2,false));
    trace.info() << ((res)? "Passed": "Failed") << std::endl; 
    nbok += res ? 1 : 0; 
    nb++;
  }

  { //accepted configuration with only X (1) and R (2)
    unsigned int t[27] = {
    1, 2, 2, 1, 1, 2, 1, 1, 1, 
    1, 2, 2, 1, 1, 2, 1, 1, 1, 
    1, 2, 2, 1, 1, 2, 1, 1, 1  };

    bool res = (basicTest<3,unsigned int*>(t+0,t+27,2,true));
    trace.info() << ((res)? "Passed": "Failed") << std::endl; 
    nbok += res ? 1 : 0; 
    nb++;
  }
  
  { //accepted configuration with O (3) twice adjacent to x
    unsigned int t[27] = {
    3,3,3,3,3,3,3,3,3, 
  1, 2, 2, 1, 1, 2, 1, 1, 2, 
    3,3,3,3,3,3,3,3,3  };

    bool res = (basicTest<3,unsigned int*>(t+0,t+27,2,true));
    trace.info() << ((res)? "Passed": "Failed") << std::endl; 
    nbok += res ? 1 : 0; 
    nb++;
  }
  
  { //accepted configuration with a contact surface 
    //between O (3) and x that is homotopic to an annulus
  unsigned int t[27] = {
    3,3,3,3,1,3,3,3,3, 
    3,3,3,3,1,3,3,3,3, 
    3,3,3,3,2,3,3,3,3  };

    bool res = (basicTest<3,unsigned int*>(t+0,t+27,2,true));
    trace.info() << ((res)? "Passed": "Failed") << std::endl; 
    nbok += res ? 1 : 0; 
    nb++;
  }

  { //rejected configuration (PRL2011 fig2a)
  unsigned int t[27] = {
    2,2,2,2,2,2,2,2,2, 
    2,2,2,2,1,2,2,2,2, 
    3,3,3,3,3,3,3,3,3  };

    bool res = (basicTest<3,unsigned int*>(t+0,t+27,2,false));
    trace.info() << ((res)? "Passed": "Failed") << std::endl; 
    nbok += res ? 1 : 0; 
    nb++;
  }

  { //rejected config (PRL2011 fig3a where X=1, R=2, A=3, B=4)
  unsigned int t[27] = {
    2,2,2,2,2,2,2,2,2, 
    1,1,3,1,1,3,1,1,3, 
    1,1,4,1,1,4,1,1,4  };

    bool res = (basicTest<3,unsigned int*>(t+0,t+27,2,false));
    trace.info() << ((res)? "Passed": "Failed") << std::endl; 
    nbok += res ? 1 : 0; 
    nb++;
  }

  { //!! rejected config (PRL2011 fig3a where X=1, R=2, A=3, B=4, 
    //modified such that the middle point of B is flipped into A)
  unsigned int t[27] = {
    2,2,2,2,2,2,2,2,2, 
    1,1,3,1,1,3,1,1,3, 
    1,1,4,1,1,3,1,1,4  };

    bool res = (basicTest<3,unsigned int*>(t+0,t+27,2,false));
    trace.info() << ((res)? "Passed": "Failed") << std::endl; 
    nbok += res ? 1 : 0; 
    nb++;
  }

  { //accepted config (PRL2011 fig3a where X=1, R=2, A=3, 
    //modified such that all the points of B are flipped into A)
  unsigned int t[27] = {
    2,2,2,2,2,2,2,2,2, 
    1,1,3,1,1,3,1,1,3, 
    1,1,3,1,1,3,1,1,3  };

    bool res = (basicTest<3,unsigned int*>(t+0,t+27,2,true));
    trace.info() << ((res)? "Passed": "Failed") << std::endl; 
    nbok += res ? 1 : 0; 
    nb++;
  }

  { //rejected config (PRL2011 fig3b where X=1, R=2, A=3, B=4)
  unsigned int t[27] = {
    2,2,2,2,2,2,2,2,2, 
    3,1,1,3,1,1,3,3,1, 
    4,4,4,4,4,4,4,4,4  };

    bool res = (basicTest<3,unsigned int*>(t+0,t+27,2,false));
    trace.info() << ((res)? "Passed": "Failed") << std::endl; 
    nbok += res ? 1 : 0; 
    nb++;
  }
  
  trace.info() << "(" << nbok << "/" << nb << ") " << endl;
  trace.endBlock();
  return (nbok == nb);
}

template<DGtal::Dimension dimension>
bool randomTest(int size)
{

  trace.beginBlock ( "random test" );

  //labels
  typedef int Label;
  std::set<Label> labels; 
  for (int i = 1; i <= 5; ++i)
    labels.insert(i); 
  Label l = *(labels.begin()); 

  //image
  typedef DGtal::HyperRectDomain<DGtal::SpaceND<dimension,int > > Domain; 
  typedef typename Domain::Point Point; 
  typedef DGtal::ImageContainerBySTLVector<Domain, Label >  Image;
  Image img( Domain( Point::diagonal(-size), Point::diagonal(size) ) );
 
  bool flagIsValid = DGtal::SimplePointHelper<Image>::generateRandomConfiguration(img, labels, 0.5); 
  
  //helper
  DGtal::SimplePointHelper<Image> h( img ); 
  if ( (flagIsValid) && (!h.isValid()) )
    return false; 

  //testing all domain points
  unsigned int c = 0; 
  Domain d = img.domain(); 
  typename Domain::ConstIterator it = d.begin(); 
  typename Domain::ConstIterator itEnd = d.end(); 
  for ( ; it != itEnd; ++it)
    {
      if ( h.isMLSimple(*it, l) )
	++c; 
    }

  unsigned int tot = std::pow(2*size+1,dimension); 
  trace.info() << "# " << c << " points ML-simples / ";
  trace.info() << tot << " points" << endl;
  trace.endBlock();
  return (c <= tot); 
}

int main(int argc, char** argv)
{
  std::srand( std::time(NULL) ); 

  trace.beginBlock ( "Testing class SimplePointHelper" );
  trace.info() << "Args:";
  for ( int i = 0; i < argc; ++i )
    trace.info() << " " << argv[ i ];
  trace.info() << endl;

  bool res = firstTest() 
    && known2dCasesTest() && known3dCasesTest() 
    && randomTest<2>(25) && randomTest<3>(10); 
  
  trace.emphase() << ( res ? "Passed." : "Error." ) << endl;
  trace.endBlock();
  return res ? 0 : 1;

}

