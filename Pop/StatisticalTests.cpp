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

// Method to compute number of polymorphic site in an alignment
// Arguments: a SiteIterator
// Return: Number of polymorphics sites	
int StatisticalTests::polymorphicSiteNumber( SiteIterator & si ) {
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

// Method to compute number of polymorphic site in an alignment
// Arguments: a SiteContainer
// Return: Number of polymorphics sites	
int StatisticalTests::polymorphicSiteNumber( const SiteContainer & v ) {

	SiteIterator *si = new NoGapSiteIterator( v );
	int S = StatisticalTests::polymorphicSiteNumber( *si );
	
	return S;
}

// Method to compute diversity estimator Theta of Watterson (1975)
// Arguments: a SiteContainer
// Return: theta of Watterson (1975)
double StatisticalTests::watterson75( const SiteContainer & v ) {
	double ThetaW;
	int n = v.getNumberOfSequences();
	double an = 0.0;
	SiteIterator *si = new NoGapSiteIterator( v );
	int S = StatisticalTests::polymorphicSiteNumber( *si );
	for ( int i = 1; i < n; i++ ) {
			an += (double) 1/i;
	}	
	ThetaW = (double) S / an;
	delete si;
	return ThetaW;
}

// Method to compute diversity estimator Theta of Tajima (1983)
// Arguments: a SiteContainer
// Return: theta of Tajima (1983)
double StatisticalTests::tajima83( const SiteContainer & v ) {
	double ThetaPi;
	int S = 0;
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
	delete si;
	return ThetaPi;
}

// Method to compute diversity estimator Theta of Watterson (1975)
// Arguments: a PolymorphismSequenceContainer
// Return: theta of Watterson (1975)
double StatisticalTests::watterson75( const PolymorphismSequenceContainer & psc ) {
	PolymorphismSequenceContainer *psci = PolymorphismSequenceContainerTools::extractIngroup(psc);
	double ThetaW;
	int n = psc.getNumberOfSequences();
	double an = 0.0;
	SiteIterator *si = new NoGapSiteIterator( *psci );
	int S = StatisticalTests::polymorphicSiteNumber( *si );
	for ( int i = 1; i < n; i++ ) {
			an += (double) 1/i;
	}	
	ThetaW = (double) S / an;	
	delete si;
	delete psci;
	return ThetaW;
}

// Method to compute diversity estimator Theta of Tajima (1983)
// Arguments: a PolymorphismSequenceContainer
// Return: theta of Tajima (1983)
double StatisticalTests::tajima83( const PolymorphismSequenceContainer & psc ) {
	PolymorphismSequenceContainer *psci = PolymorphismSequenceContainerTools::extractIngroup(psc);
	double ThetaPi;
	int S = 0;
	const Site *site;
	int n = psci -> getNumberOfSequences();
	double etha[20];
	double somme = 0.0;
	SiteIterator *si = new NoGapSiteIterator( *psci );
	while ( si -> hasMoreSites() ) {
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
	ThetaPi = (double) S - somme;
	delete si;
	delete psci;
	return ThetaPi;
}
