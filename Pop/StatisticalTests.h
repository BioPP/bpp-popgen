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

        /**
         * @brief Compute number of polymorphic site in an alignment
	 *
         * Gaps are consider as mutations so if you want number of
	 * polymorphic site, you have to give a NonGapSiteIterator
	 *
         * @param si a SiteIterator
         */ 
	 static unsigned int polymorphicSiteNumber( SiteIterator & si );

        /**
         * @brief Compute number of polymorphic site in an alignment
         *
	 * Gaps are consider as mutation variant
	 *
         * @param v a SiteContainer
         */
	 static unsigned int polymorphicSiteNumber( const SiteContainer & v );

        /**
         * @brief Compute diversity estimator Theta of Watterson (1975)
         *
         * @param v a SiteContainer
         */
	 static double watterson75( const SiteContainer & v );

        /**
         * @brief Compute diversity estimator Theta of Tajima (1983)
         *
         * @param v a SiteContainer
         */ 
	 static double tajima83( const SiteContainer & v );

        /**
         * @brief Compute diversity estimator Theta of Watterson (1975)
         *
         * @param v a PolymorphismSequenceContainer
	 * @param gapflag flag set by default to true if you don't want to
	 * take gap into account
         */	    
	 static double watterson75( const PolymorphismSequenceContainer & v, bool gapflag = true );

        /**
         * @brief Compute diversity estimator Theta of Tajima (1983)
         *
         * @param v a PolymorphismSequenceContainer
	 * @param gapflag flag set by default to true if you don't want to
	 * take gap into account
         */ 
	 static double tajima83( const PolymorphismSequenceContainer & v, bool gapflag = true );

        /**
         * @brief Return the number of haplotype in the sample. 
	 * Depaulis and Veuille (1998)
         *
         * @param v a PolymorphismSequenceContainer
	 * @param gapflag flag set by default to true if you don't want to
	 * take gap into account
         */ 
	 static unsigned int DVK ( const PolymorphismSequenceContainer & v, bool gapflag = true );

        /**
         * @brief Return the haplotype diversity of a sample. 
	 * Depaulis and Veuille (1998)
         *
         * @param v a PolymorphismSequenceContainer
	 * @param gapflag flag set by default to true if you don't want to
	 * take gap into account
         */ 
	 static double DVH ( const PolymorphismSequenceContainer & v, bool gapflag = true );
              
        /*******************************************************************************/
};
#endif // _STATISTICALTESTS_H_
