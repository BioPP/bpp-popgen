/*
 * File Statistics.cpp
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
unsigned int StatisticalTests::polymorphicSiteNumber( SiteIterator & si ) {
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
unsigned int StatisticalTests::polymorphicSiteNumber( const SiteContainer & v ) {

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
	unsigned int S = StatisticalTests::polymorphicSiteNumber( *si );
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
double StatisticalTests::watterson75( const PolymorphismSequenceContainer & psc, bool gapflag ) {
	PolymorphismSequenceContainer *psci = PolymorphismSequenceContainerTools::extractIngroup(psc);
	double ThetaW;
	int n = psc.getNumberOfSequences();
	double an = 0.0;
	SiteIterator *si;
	if (gapflag) 
		si = new NoGapSiteIterator( *psci );
	else
		si = new CompleteSiteIterator( *psci );
	unsigned int S = StatisticalTests::polymorphicSiteNumber( *si );
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
double StatisticalTests::tajima83( const PolymorphismSequenceContainer & psc, bool gapflag ) {
	PolymorphismSequenceContainer *psci = PolymorphismSequenceContainerTools::extractIngroup(psc);
	double ThetaPi;
	int S = 0;
	const Site *site;
	vector<double> etha;
	unsigned int alphasize; 
	double somme = 0.0;
	unsigned int samplesize = psci -> getNumberOfSequences();
	double denom = (double(samplesize)* (double(samplesize) - 1.0));
	SiteIterator *si;
	if (gapflag) {
		si = new NoGapSiteIterator( *psci );
		alphasize = (psc.getAlphabet())->getSize();
	} else {
		si = new CompleteSiteIterator( *psci );
		alphasize = (psc.getAlphabet())->getSize() + 1;
	}
	etha.resize(alphasize);
	while ( si -> hasMoreSites() ) {
		site = si->nextSite();
		if ( !SiteTools::isConstant(*site) ) {
			S++;
			for ( int i = 0; i < alphasize; i++ ) {
				etha[i] = 0;
			}
			for ( int j = 0; j < samplesize; j++ ) {
				if (site -> getValue( j ) < alphasize)
					etha[ site -> getValue( j ) ]+= psci -> getSequenceStrength( j );
			}			
			for ( int i = 0; i < alphasize; i++ ) {
				somme += (etha[i] * (etha[i] - 1)) / denom;
			}
		}
	}
	ThetaPi = (double) S - somme;
	delete si;
	delete psci;
	return ThetaPi;
}

// Return the number of haplotype in the sample. Depaulis and Veuille (1998)
// Arguments: a PolymorphismSequenceContainer
// Return: K (Depaulis and Veuille 1998)
unsigned int StatisticalTests::DVK( const PolymorphismSequenceContainer & psc, bool gapflag ) {
	PolymorphismSequenceContainer *psci = PolymorphismSequenceContainerTools::extractIngroup(psc);
	PolymorphismSequenceContainer *sc;
	if (gapflag)
		sc = PolymorphismSequenceContainerTools::getSitesWithoutGaps ( *psci );
	else
		sc = psci;
	int K = 0;
	vector<string> pscvector;
	pscvector.push_back(sc->toString(0));
	K++;
	for ( unsigned int i = 1; i < sc->getNumberOfSequences(); i++ ) {
		bool uniq = true;	
		string query = sc->toString(i);
		for ( vector<string>::iterator it = pscvector.begin(); it != pscvector.end(); it++ ) {
			if ( query.compare(*it) == 0 ) {
				uniq = false;
				break;
			}	
		}
		if (uniq) {
			K++;
			pscvector.push_back(query);
		}				
	}
	delete sc;
	if(gapflag)
		delete psci;	
	return K;
}
