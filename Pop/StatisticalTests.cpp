/*
 * File StatisticalTests.cpp
 * Author : Eric Bazin <bazin@univ-montp2.fr>
 * Last modification : Friday December 5 2003
*/

#include "StatisticalTests.h" // class's header file
#include "PolymorphismSequenceContainerTools.h" 

// From the STL:
#include <ctype.h>
#include <cmath>
#include <iostream>

// From SeqLib:
#include <Seq/Site.h>	
#include <Seq/SiteTools.h>

using namespace std;

StatisticalTests::~StatisticalTests() {}

// Methods to compute number of non gap polymorphic site in an alignment
// Arguments: a SiteIterator
// Return: Number of polymorphics sites	
double StatisticalTests::polymorphicSiteNumber( SiteIterator & si ) {
	int S=0;
	const Site *site;
	while ( si.hasMoreSites() ) {
		site=si.nextSite();
		if ( !SiteTools::isConstant(*site) ) {
			S++;
		}
	}
	return S;
}

// Methods to compute number of non gap polymorphic site in an ingroup alignment
// Arguments: a SiteContainer
// Return: Number of polymorphics sites	
double StatisticalTests::polymorphicSiteNumber( const SiteContainer & v ) {

	SiteIterator *si = new NoGapSiteIterator( v );
	double S = StatisticalTests::polymorphicSiteNumber( *si );
	
	return S;
}

// Methods to compute diversity estimator Theta of Watterson (1975)
// Arguments: a SiteContainer
// Return: theta of Watterson (1975)
double StatisticalTests::watterson75( const SiteContainer & v ) {
	double ThetaW;
	int n = v.getNumberOfSequences();
	double an = 0.0;
	SiteIterator *si = new NoGapSiteIterator( v );
	double S = StatisticalTests::polymorphicSiteNumber( *si );
	for ( int i = 1; i < n; i++ ) {
			an += (double) 1/i;
	}	
	ThetaW = S / an;
	return ThetaW;
}

// Methods to compute diversity estimator Theta of Tajima (1983)
// Arguments: a SiteContainer
// Return: theta of Tajima (1983)
double StatisticalTests::tajima83( const SiteContainer & v ) {
	double ThetaPi;
	double S = 0;
	const Site *site;
	int n = v.getNumberOfSequences();
	double etha[20];
	double somme = 0.0;
	SiteIterator *si = new NoGapSiteIterator( v );
	while ( si->hasMoreSites() ) {
		site = si->nextSite();
	if ( !SiteTools::isConstant(*site) ) {
			S++;
			for ( int i = 0; i < 4; i++ ) {
				etha[i] = 0;
			}
			for ( int j = 0; j < n; j++ ) {
				etha[site->getValue( j )]++;
			}			
			for ( int i = 0; i < 4; i++ ) {
				somme += (etha[i] * (etha[i] - 1)) / (n * (n  - 1));
			}
		}
	}
	ThetaPi = S - somme;
	return ThetaPi;
}

// Methods to compute diversity estimator Theta of Watterson (1975)
// Arguments: a PolymorphismSequenceContainer
// Return: theta of Watterson (1975)
double StatisticalTests::watterson75( const PolymorphismSequenceContainer & psc ) {
	double ThetaW;
	int n = psc.getNumberOfSequences();
	double an = 0.0;
	PolymorphismSequenceContainer *psci =
	PolymorphismSequenceContainerTools::extractIngroup( psc );
	SiteIterator *si = new NoGapSiteIterator( *psci );
	double S = StatisticalTests::polymorphicSiteNumber( *si );
	for ( int i = 1; i < n; i++ ) {
			an += (double) 1/i;
	}	
	ThetaW = S / an;
	delete psci;
	return ThetaW;
}

// Methods to compute diversity estimator Theta of Tajima (1983)
// Arguments: a PolymorphismSequenceContainer
// Return: theta of Tajima (1983)
double StatisticalTests::tajima83( const PolymorphismSequenceContainer & v ) {
	double ThetaPi;
	double S = 0;
	const Site *site;
	int n = v.getNumberOfSequences();
	double etha[20];
	double somme = 0.0;
	SiteIterator *si = new NoGapSiteIterator( v );
	while ( si->hasMoreSites() ) {
		site = si->nextSite();
	if ( !SiteTools::isConstant(*site) ) {
			S++;
			for ( int i = 0; i < 4; i++ ) {
				etha[i] = 0;
			}
			for ( int j = 0; j < n; j++ ) {
				etha[site->getValue( j )]++;
			}			
			for ( int i = 0; i < 4; i++ ) {
				somme += (etha[i] * (etha[i] - 1)) / (n * (n  - 1));
			}
		}
	}
	ThetaPi = S - somme;
	return ThetaPi;
}
