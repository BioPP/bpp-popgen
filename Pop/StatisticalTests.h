/*
 * File StatisticalTests.h
 * Author : Eric Bazin <bazin@univ-montp2.fr>
 * Last modification : Friday December 5 2003
*/

//Secured inclusion of header's file
#ifndef _STATISTICALTESTS_H_
#define _STATISTICALTESTS_H_

//From the SeqLib library
#include <Seq/SiteIterator.h>
#include <Seq/SiteContainer.h>

//From the PolyLib library
#include "PolymorphismSequenceContainer.h"

using namespace std;

class StatisticalTests
{
    public: 
        // Class destructor:
        ~StatisticalTests();
	
        /*******************************************************************************/
    public:

	    static double polymorphicSiteNumber( SiteIterator & si );

	    static double polymorphicSiteNumber( const SiteContainer & v );

	    static double watterson75( const SiteContainer & v );
 
	    static double tajima83( const SiteContainer & v );
	    
	    static double watterson75( const PolymorphismSequenceContainer & v );
 
	    static double tajima83( const PolymorphismSequenceContainer & v );	
	        
        /*******************************************************************************/
};
#endif // _STATISTICALTESTS_H_
