/*
 * File Statistics.cpp
 * Author : Eric Bazin <bazin@univ-montp2.fr>
 *          Sylvain Gailard <yragael2001@yahoo.fr>
 * Last modification : Wednesday July 28 2004
 */

#include "SequenceStatistics.h" // class's header file
#include "PolymorphismSequenceContainerTools.h" 

// From the STL:
#include <ctype.h>
#include <cmath>
#include <iostream>

// From SeqLib:
#include <Seq/Site.h>	
#include <Seq/SiteTools.h>

using namespace std;

SequenceStatistics::~SequenceStatistics() {}

// Method to compute number of polymorphic site in an alignment
// Arguments: a SiteIterator
// Return: Number of polymorphics sites	
unsigned int SequenceStatistics::polymorphicSiteNumber( SiteIterator & si ) {
	unsigned int S=0;
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
unsigned int SequenceStatistics::polymorphicSiteNumber( const SiteContainer & v ) {

	SiteIterator *si = new NoGapSiteIterator( v );
	unsigned int S = SequenceStatistics::polymorphicSiteNumber( *si );
	
	return S;
}

// Method to compute number of singleton nucleotides in an alignment
// Arguments: a SiteIterator
// Return: Number of singleton nucleotides
unsigned int SequenceStatistics::countSingleton(SiteIterator & si) {
	unsigned int nus = 0;
	const Site * site;
	while (si.hasMoreSites()) {
		site = si.nextSite();
		map<int, unsigned int> states_count = SymbolListTools::getCounts(* site);
		for (map<int, unsigned int>::iterator it = states_count.begin() ; it != states_count.end() ; it++)
			if (it->second == 1)
				nus++;
	}
	return nus;
}

// Method to compute total number of mutation under an infinite site model in an alignment
// Arguments: a SiteIterator
// Return: Total number of mutations
unsigned int SequenceStatistics::totNumberMutations(SiteIterator & si) {
	unsigned int tnm = 0;
	const Site * site;
	while (si.hasMoreSites()) {
		site = si.nextSite();
		unsigned int tmp_count = 0;
		map<int, unsigned int> states_count = SymbolListTools::getCounts(* site);
		for (map<int, unsigned int>::iterator it = states_count.begin() ; it != states_count.end() ; it++)
			if (it->first >= 0)
				tmp_count++;
		if (tmp_count > 0)
			tmp_count--;
		tnm += tmp_count;
	}
	return tnm;
}

// Method to compute diversity estimator Theta of Watterson (1975)
// Arguments: a SiteContainer
// Return: theta of Watterson (1975)
double SequenceStatistics::watterson75( const SiteContainer & v ) {
	double ThetaW;
	int n = v.getNumberOfSequences();
	double an = 0.0;
	SiteIterator *si = new NoGapSiteIterator( v );
	unsigned int S = polymorphicSiteNumber( *si );
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
double SequenceStatistics::tajima83( const SiteContainer & v ) {
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
double SequenceStatistics::watterson75(const PolymorphismSequenceContainer & psc, bool gapflag) {
	double ThetaW;
	unsigned int n = psc.getNumberOfSequences();
	SiteIterator * si = NULL;
	if (gapflag) 
		si = new CompleteSiteIterator(psc);
	else
		si = new SimpleSiteIterator(psc);
	unsigned int S = polymorphicSiteNumber(* si);
	map<string, double> values = _getUsefullValues(n);
	ThetaW = (double) S / values["a1"];	
	delete si;
	return ThetaW;
}

// Method to compute diversity estimator Theta of Tajima (1983)
// Arguments: a PolymorphismSequenceContainer
// Return: theta of Tajima (1983)
double SequenceStatistics::tajima83(const PolymorphismSequenceContainer & psc, bool gapflag) {
	unsigned int alphabet_size = (psc.getAlphabet())->getSize();
	const Site * site;
	SiteIterator *si;
	double value = 0.;
	if (gapflag)
		si = new CompleteSiteIterator(psc);
	else
		si = new SimpleSiteIterator(psc);
	while (si->hasMoreSites()) {
		site = si->nextSite();
		if (! SiteTools::isConstant(* site)) {
			map<int, unsigned int> count = SymbolListTools::getCounts(* site);
			map<int, unsigned int> tmp_k;
			unsigned int tmp_n = 0;
			for (map<int, unsigned int>::iterator it = count.begin() ; it != count.end() ; it++)
				if (it->first >= 0 && it->first < (int) alphabet_size) {
					tmp_k[it->first] = it->second * (it->second - 1);
					tmp_n += it->second;
				}
			for(map<int, unsigned int>::iterator it = tmp_k.begin() ; it != tmp_k.end() ; it++)
				value += (double) it->second / (tmp_n * (tmp_n - 1));
		}
	}
	return 1 - value;
}

// Method to compute Tajima D test (1989)
// Arguments: a PolymorphismSequenceContainer
// Return: Tajima's D (1989)
double SequenceStatistics::tajimaD(const PolymorphismSequenceContainer & psc, bool gapflag) {
	SiteIterator * si;
	if (gapflag)
		si = new NoGapSiteIterator(psc);
	else
		si = new CompleteSiteIterator(psc);
	unsigned int S = polymorphicSiteNumber(* si);
	double tajima = tajima83(psc, gapflag);
	double watterson = watterson75(psc, gapflag);
	unsigned int n = psc.getNumberOfSequences();
	map<string, double> values = _getUsefullValues(n);
	return (tajima - watterson) / sqrt((values["e1"] * S) + (values["e2"] * S * (S - 1)));
}

// Return the number of haplotype in the sample. Depaulis and Veuille (1998)
// Arguments: a PolymorphismSequenceContainer
// Return: K (Depaulis and Veuille 1998)
unsigned int SequenceStatistics::DVK( const PolymorphismSequenceContainer & psc, bool gapflag ) {
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

// Return the haplotype diversity of a sample. Depaulis and Veuille (1998)
// Arguments: a PolymorphismSequenceContainer
// Return: H (Depaulis and Veuille 1998)
double SequenceStatistics::DVH( const PolymorphismSequenceContainer & psc, bool gapflag ) {
	PolymorphismSequenceContainer *psci = PolymorphismSequenceContainerTools::extractIngroup(psc);
	PolymorphismSequenceContainer *sc;
	if (gapflag)
		sc = PolymorphismSequenceContainerTools::getSitesWithoutGaps ( *psci );
	else
		sc = psci;
	double H = 0.0;
	unsigned int nbSeq;
	vector<string> pscvector;
	vector<int> effvector;
	pscvector.push_back(sc -> toString(0));
	effvector.push_back(sc -> getSequenceCount(0));
	nbSeq = sc -> getSequenceCount(0);
	for ( unsigned int i = 1; i < sc->getNumberOfSequences(); i++ ) {
		nbSeq += sc -> getSequenceCount(i);
		bool uniq = true;	
		string query = sc -> toString(i);
		for ( vector<string>::iterator it = pscvector.begin(); it != pscvector.end(); it++ ) {
			if ( query.compare(*it) == 0 ) {
				effvector[effvector.size() - 1] += sc -> getSequenceCount(i);
				uniq = false;
				break;
			}	
		}
		if (uniq) {
			pscvector.push_back(query);
			effvector.push_back(sc -> getSequenceCount(i));
		}				
	}
	for ( unsigned int i = 0; i < effvector.size(); i++ ) {
		H -= ( (double) effvector[i] / (double) nbSeq ) * ( (double) effvector[i] / (double) nbSeq );
	}
	H += 1.0;
	delete sc;
	if(gapflag)
		delete psci;	
	return H;
}

map<string, double> SequenceStatistics::_getUsefullValues(unsigned int n) {
	map<string, double> values;
	values["a1"] = 0.;
	values["a2"] = 0.;
	values["a1n"] = 0.;
	values["b1"] = 0.;
	values["b2"] = 0.;
	values["c1"] = 0.;
	values["c2"] = 0.;
	values["cn"] = 0.;
	values["dn"] = 0.;
	values["e1"] = 0.;
	values["e2"] = 0.;
	if (n > 1) {
		for (unsigned int i = 1 ; i < n ; i++) {
			values["a1"] += 1. / i;
			values["a2"] += 1. / (i * i);
		}
		values["a1n"] = values["a1"] + (1. / n);
		values["b1"] = (double) (n + 1) / (3 * (n - 1));
		values["b2"] = (double) 2 * ((n * n) + n + 3) / (9 * n * (n - 1));
		values["c1"] = values["b1"] - (1. / values["a1"]);
		values["c2"] = values["b2"] - ((n + 2) / (values["a1"] * n)) + (values["a2"] / (values["a1"] * values["a1"]));
		values["cn"] = 2. * ((n * values["a1"]) - (2 * (n - 1))) / ((n -1) * (n - 2));
		values["dn"] = values["cn"] + ((n - 2) / ((n - 1) * (n - 1))) + ((2 / (n - 1)) * ((3 / 2) - (((2 * values["a1n"]) - 3) / (n - 2)) - (1 / n)));
		values["e1"] = values["c1"] / values["a1"];
		values["e2"] = values["c2"] / ((values["a1"] * values["a1"]) + values["a2"]);
	}
	return values;
}
