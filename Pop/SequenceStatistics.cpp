/*
 * File SequenceStatistics.cpp
 * Author : Eric Bazin <bazin@univ-montp2.fr>
 *          Sylvain Gailard <yragael2001@yahoo.fr>
 * Last modification : Friday August 06 2004
 *
 * Copyright (C) 2004 Eric Bazin, Sylvain Gaillard and the
 *                    PopGenLib Development Core Team
 *
 * This file is part of PopGenLib.
 *
 * PopGenLib is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * PopGenLib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PopGenLib; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include "SequenceStatistics.h" // class's header file
#include "PolymorphismSequenceContainerTools.h"
#include "PolymorphismSequenceContainer.h"

// From the STL:
#include <ctype.h>
#include <cmath>
#include <iostream>
#include <math>


using namespace std;

// From SeqLib:
#include <Seq/Site.h>
#include <Seq/SiteTools.h>
#include <Seq/StringSequenceTools.h>
#include <Seq/CodonSiteTools.h>
#include <Seq/DNA.h>
#include <Seq/StandardCodonAlphabet.h>
#include <Seq/StandardGeneticCode.h>

// from NumCalc
#include <NumCalc/VectorTools.h>
using namespace VectorOperators;
using namespace VectorFunctions;
using namespace VectorStatTools;


SequenceStatistics::~SequenceStatistics() {}

//******************************************************************************************************************
//Basic statistics
//******************************************************************************************************************

// Method to compute number of polymorphic site in an alignment
// Return: Number of polymorphics sites
unsigned int SequenceStatistics::polymorphicSiteNumber(const PolymorphismSequenceContainer & psc, bool gapflag) {
	unsigned int S=0;
	const Site *site;
	SiteIterator * si = NULL;
	if (gapflag)
		si = new CompleteSiteIterator(psc);
	else
		si = new SimpleSiteIterator(psc);
	while ( si->hasMoreSites() ) {
		site=si->nextSite();
		if ( !SiteTools::isConstant(*site) ) {
			S++;
		}
	}
        delete si;
	return S;
}



// Method to compute number of parsimony informative sites in an alignment
// Arguments: a SiteIterator
// Return: Number of parsimony informative site
unsigned int SequenceStatistics::parsimonyInformativeSiteNumber( SiteIterator & si ) {
	unsigned int S=0;
	const Site *site;
	while ( si.hasMoreSites() ) {
		site=si.nextSite();
		if (SiteTools::isParsimonyInformativeSite(*site)) {
			S++;
		}
	}
	return S;
}

// Method to compute number of parsimony informative sites in an alignment
// Arguments: a SiteContainer   a boolean (gapflag: true: do not count site with gap or undetermined)
// Return: Number of parsimony informative site
unsigned int SequenceStatistics::parsimonyInformativeSiteNumber(const SiteContainer & v, bool gapflag) {
        SiteIterator *si;
        if(gapflag) si = new CompleteSiteIterator(v);
	else si = new SimpleSiteIterator(v);
	unsigned int S = SequenceStatistics::parsimonyInformativeSiteNumber(*si);
        delete si;
	return S;
}


// Method to compute number of singleton nucleotides in an alignment
// Return: Number of singleton nucleotides
unsigned int SequenceStatistics::countSingleton(const PolymorphismSequenceContainer & psc, bool gapflag) {
	unsigned int nus = 0;
	const Site * site;
	SiteIterator * si = NULL;
	if (gapflag)
		si = new CompleteSiteIterator(psc);
	else
		si = new SimpleSiteIterator(psc);
	while (si->hasMoreSites()) {
		site = si->nextSite();
		nus += _getSingletonNumber(* site);
	}
        delete si;
	return nus;
}


// Method to compute number of triplet sites in an alignment
// Arguments: a SiteIterator
// Return: Number of triplet sites
unsigned int SequenceStatistics::tripletNumber( SiteIterator & si ) {
	int S=0;
	const Site *site;
	while ( si.hasMoreSites() ) {
		site=si.nextSite();
		if ( SiteTools::isTriplet(*site) ) {
			S++;
		}
	}
	return S;
}


// Method to compute number of triplet sites in an alignment
// Arguments: a SiteContainer
// Return: Number of triplet sites
unsigned int SequenceStatistics::tripletNumber(const SiteContainer & v, bool gapflag) {
        SiteIterator *si;
        if(gapflag) si = new CompleteSiteIterator(v);
	else si = new SimpleSiteIterator(v);
	int S = SequenceStatistics::tripletNumber( *si );
        delete si;
	return S;
}




// Method to compute total number of mutation under an infinite site model in an alignment
// Return: Total number of mutations
unsigned int SequenceStatistics::totNumberMutations(const PolymorphismSequenceContainer & psc, bool gapflag) {
	unsigned int tnm = 0;
	const Site * site;
	SiteIterator * si = NULL;
	if (gapflag)
		si = new CompleteSiteIterator(psc);
	else
		si = new SimpleSiteIterator(psc);
	while (si->hasMoreSites()) {
		site = si->nextSite();
		tnm += _getMutationNumber(* site);
	}
        delete si;
	return tnm;
}

//******************************************************************************************************************
//GC statistics
//******************************************************************************************************************

// Method to compute mean GC content in an alignement
// Return: mean GC content
double SequenceStatistics::GCcontent(const PolymorphismSequenceContainer & psc) {
        SiteContainer* sc = new VectorSiteContainer(psc);
        map<int, double> freqs = SequenceContainerTools::getFrequencies(*sc);
        delete sc;
        return (freqs[1] + freqs[2]) / (freqs[0] + freqs[1] + freqs[2] + freqs[3]);
}



//Method that gives the number of GC alleles and the total number of allele at polymorphic sites
//G vs C and A vs T polymorphism are not taken into account
//Return: a vector with the total number of alleles and the number of GC alleles
vector<unsigned int> SequenceStatistics::GCpolymorphism(const PolymorphismSequenceContainer & psc, bool stopflag) {
	unsigned int nbMut = 0;
	unsigned int nbGC = 0;
	const unsigned int nbSeq = psc.getNumberOfSequences();
	vector<unsigned int> vect(2);
	const Site * site;
	SiteIterator * si = NULL;
    if(stopflag) si = new CompleteSiteIterator(psc);
    else si = new NoGapSiteIterator(psc);
	while (si->hasMoreSites()) {
		site = si->nextSite();
		if(!SiteTools::isConstant(*site)) {
        	double freqGC = StringSequenceTools::getGCcontent(site->toString(),0, nbSeq);
			if(freqGC > 0 && freqGC < 1) {
                                nbMut += (unsigned int) nbSeq;
                                double adGC = freqGC*nbSeq;
				nbGC += (unsigned int)adGC;
			}
		}
	}
	vect[0]=nbMut;
	vect[1]=nbGC;
    delete si;
    return vect;
}



//******************************************************************************************************************
//Diversity statistics
//******************************************************************************************************************

// Method to compute diversity estimator Theta of Watterson (1975)
// Arguments: a SiteContainer
// Return: theta of Watterson (1975)
double SequenceStatistics::watterson75( const SiteContainer & v ) {
	double ThetaW;
	int n = v.getNumberOfSequences();
	double an = 0.0;
	SiteIterator *si = new NoGapSiteIterator( v );
	unsigned int S = polymorphicSiteNumber( v );
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
	unsigned int S = polymorphicSiteNumber(psc, gapflag);
	map<string, double> values = _getUsefullValues(n);
	ThetaW = (double) S / values["a1"];
	return ThetaW;
}

// Method to compute diversity estimator Theta of Tajima (1983)
// Arguments: a PolymorphismSequenceContainer
// Return: theta of Tajima (1983)
double SequenceStatistics::tajima83(const PolymorphismSequenceContainer & psc, bool gapflag) {
	unsigned int alphabet_size = (psc.getAlphabet())->getSize();
	const Site * site;
	SiteIterator *si;
	double value2 = 0.;
	if (gapflag)
		si = new CompleteSiteIterator(psc);
	else
		si = new SimpleSiteIterator(psc);
	while (si->hasMoreSites()) {
		site = si->nextSite();
		if (! SiteTools::isConstant(* site)) {
			double value = 0.;
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
			value2 += 1. - value;
		}
	}
        delete si;
	return value2;
}


// Return the number of haplotype in the sample. Depaulis and Veuille (1998)
// Arguments: a PolymorphismSequenceContainer
// Return: K (Depaulis and Veuille 1998)
unsigned int SequenceStatistics::DVK( const PolymorphismSequenceContainer & psc, bool gapflag ) {
	PolymorphismSequenceContainer *sc;
	if (gapflag)
		sc = PolymorphismSequenceContainerTools::getSitesWithoutGaps(psc);
	else
		sc = new PolymorphismSequenceContainer(psc);
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
	return K;
}

// Return the haplotype diversity of a sample. Depaulis and Veuille (1998)
// Arguments: a PolymorphismSequenceContainer
// Return: H (Depaulis and Veuille 1998)
double SequenceStatistics::DVH( const PolymorphismSequenceContainer & psc, bool gapflag ) {
	PolymorphismSequenceContainer *sc;
	if (gapflag)
		sc = PolymorphismSequenceContainerTools::getSitesWithoutGaps(psc);
	else
		sc = new PolymorphismSequenceContainer(psc);
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
	return H;
}


//******************************************************************************************************************
//Synonymous and non-synonymous polymorphism
//******************************************************************************************************************


// Method to compute the number of codon sites with stop codon
// Arguments: a SiteIterator
// Return: Number of codon sites with stop codon
unsigned int SequenceStatistics::StopCodonSiteNumber(SiteIterator & si, const CodonAlphabet & alpha) {
	unsigned int S=0;
	const Site *site;
	while ( si.hasMoreSites() ) {
		site=si.nextSite();
		if (CodonSiteTools::hasStop(*site,alpha)) {
			S++;
		}
	}
	return S;
}

// Method to compute the number of codon sites with stop codon
// Arguments: a SiteContainer, a boolean
// Return: Number of codon sites with stop codon
unsigned int SequenceStatistics::StopCodonSiteNumber(const SiteContainer & v, bool gapflag) {
    SiteIterator *si = NULL;
    const NucleicAlphabet* na = new DNA();
    const CodonAlphabet* ca = new StandardCodonAlphabet(na);
    if(gapflag) si = new NoGapSiteIterator(v);
    else si = new SimpleSiteIterator(v);
    unsigned int S = SequenceStatistics::StopCodonSiteNumber(*si,*ca);
    delete si;
    delete na;
    delete ca;
    return S;
}


// Method to compute the number of polymorphic codon with only one mutated site
// Arguments: a SiteIterator, a NucleicAlphabet, a CodonAlphabet (Alphabets are given by the following method)
// Return: Number of monosite polymorphic codons
unsigned int SequenceStatistics::monoSitePolymorphicCodonNumber(SiteIterator & si, const NucleicAlphabet & na, const CodonAlphabet & ca) {
	unsigned int S=0;
	const Site *site;
	while ( si.hasMoreSites() ) {
		site=si.nextSite();
		if (CodonSiteTools::isMonoSitePolymorphic(*site,na,ca)) {
			S++;
		}
	}
	return S;
}

// Method to compute the number of polymorphic codon with only one mutated site
// Arguments: a SiteContainer, two booleans
// Return: Number of monosite polymorphic codons
unsigned int SequenceStatistics::monoSitePolymorphicCodonNumber(const SiteContainer & v, bool stopflag, bool gapflag) {
    SiteIterator *si = NULL;
    const NucleicAlphabet* na = new DNA();
    const CodonAlphabet* ca = new StandardCodonAlphabet(na);
    if(stopflag) si = new CompleteSiteIterator(v);
    else {
		if(gapflag) si = new NoGapSiteIterator(v);
		else si = new SimpleSiteIterator(v);
	}
    unsigned int S = SequenceStatistics::monoSitePolymorphicCodonNumber(*si,*na,*ca);
    delete si;
    delete na;
    delete ca;
    return S;
}


// Method to compute the number of synonymous polymorphic codon sites
// Arguments: a SiteIterator, a NucleicAlphabet, a CodonAlphabet (Alphabets are given by the following method)
// Return: Number of synonymous codon sites
unsigned int SequenceStatistics::synonymousPolymorphicCodonNumber(SiteIterator & si, const NucleicAlphabet & na, const CodonAlphabet & ca) {
	unsigned int S=0;
	const Site *site;
	while ( si.hasMoreSites() ) {
		site=si.nextSite();
		if (CodonSiteTools::isSynonymousPolymorphic(*site,na,ca)) {
			S++;
		}
	}
	return S;
}

// Method to compute the number of synonymous polymorphic codon sites
// Arguments: a SiteContainer, a boolean
// Return: Number of synonymous codon sites
unsigned int SequenceStatistics::synonymousPolymorphicCodonNumber(const SiteContainer & v, bool stopflag) {
    SiteIterator* si = NULL;
    const NucleicAlphabet* na = new DNA();
    const CodonAlphabet* ca = new StandardCodonAlphabet(na);
    if(stopflag) si = new CompleteSiteIterator(v);
    else new NoGapSiteIterator(v);
    unsigned int S = SequenceStatistics::synonymousPolymorphicCodonNumber(*si,*na,*ca);
    delete si;
    delete na;
    delete ca;
    return S;
}


// Method to compute the synonymous nucleotide diversity pi
// Arguments: a SiteIterator, a CodonAlphabet, a GeneticCode (given by the following method)
// Return: pi synonymous
double SequenceStatistics::piSynonymous(SiteIterator & si, const CodonAlphabet & ca, const GeneticCode & gc, bool minchange) {
	double S=0;
	const Site *site;
	while(si.hasMoreSites()) {
		site=si.nextSite();
                S += CodonSiteTools::piSynonymous(*site,ca,gc,minchange);
	}
	return S;
}


// Method to compute the synonymous nucleotide diversity pi
// Arguments: a SiteContainer, a boolean
// Return: pi synonymous
double SequenceStatistics::piSynonymous(const SiteContainer & v, bool stopflag, bool minchange) {
    SiteIterator *si = NULL;
    NucleicAlphabet* na = new DNA();
    CodonAlphabet* ca = new StandardCodonAlphabet(na);
    GeneticCode* gc = new StandardGeneticCode(na);
    if(stopflag) si = new CompleteSiteIterator(v);
    else si = new NoGapSiteIterator(v);
    double S = SequenceStatistics::piSynonymous(*si,*ca,*gc,minchange);
    delete si;
    delete na;
    delete ca;
    delete gc;
    return S;
}


// Method to compute the non-synonymous nucleotide diversity pi
// Arguments: a SiteIterator
// Return: pi synonymous
double SequenceStatistics::piNonSynonymous(SiteIterator & si, const NucleicAlphabet & na, const CodonAlphabet & ca, const GeneticCode & gc, bool minchange) {
	double S=0;
	const Site *site;
	while(si.hasMoreSites()) {
		site=si.nextSite();
                S += CodonSiteTools::piNonSynonymous(*site,na,ca,gc,minchange);
	}
	return S;
}


// Method to compute the non-synonymous nucleotide diversity pi
// Arguments: a SiteContainer, a boolean
// Return: pi synonymous
double SequenceStatistics::piNonSynonymous(const SiteContainer & v, bool stopflag, bool minchange) {
    SiteIterator *si = NULL;
    NucleicAlphabet* na = new DNA();
    CodonAlphabet* ca = new StandardCodonAlphabet(na);
    GeneticCode* gc = new StandardGeneticCode(na);
    if(stopflag) si = new CompleteSiteIterator(v);
    else si = new NoGapSiteIterator(v);
    double S = SequenceStatistics::piNonSynonymous(*si,*na,*ca,*gc,minchange);
    delete si;
    delete na;
    delete ca;
    delete gc;
    return S;
}

// Method to compute the mean number of synonymous site in an alignment
// Arguments: a SiteIterator
// Return: mean number of synonymous site
double SequenceStatistics::meanSynonymousSitesNumber(SiteIterator & si, const CodonAlphabet & ca, const GeneticCode &gc, double ratio) {
	double S=0;
	const Site *site;
	while(si.hasMoreSites()) {
		site=si.nextSite();
                S += CodonSiteTools::MeanNumberOfSynonymousPositions(*site,ca,gc,ratio);
	}
	return S;
}


// Method to compute the mean number of synonymous site in an alignment
// Arguments: a SiteIterator, a boolean
// Return: mean number of synonymous site
double SequenceStatistics::meanSynonymousSitesNumber(const SiteContainer & v, double ratio, bool stopflag) {
    SiteIterator *si = NULL;
    NucleicAlphabet* na = new DNA();
    CodonAlphabet* ca = new StandardCodonAlphabet(na);
    GeneticCode* gc = new StandardGeneticCode(na);
    if(stopflag) si = new CompleteSiteIterator(v);
    else si = new NoGapSiteIterator(v);
    double S = SequenceStatistics::meanSynonymousSitesNumber(*si,*ca,*gc,ratio);
    delete si;
    delete na;
    delete ca;
    delete gc;
    return S;
}


//******************************************************************************************************************
//Statistical tests
//******************************************************************************************************************


// Method to compute Tajima D test (1989)
// Arguments: a PolymorphismSequenceContainer
// Return: Tajima's D (1989)
double SequenceStatistics::tajimaDSS(const PolymorphismSequenceContainer & psc, bool gapflag) {
	unsigned int S = polymorphicSiteNumber(psc, gapflag);
	double tajima = tajima83(psc, gapflag);
	double watterson = watterson75(psc, gapflag);
	unsigned int n = psc.getNumberOfSequences();
	map<string, double> values = _getUsefullValues(n);
	return (tajima - watterson) / sqrt((values["e1"] * S) + (values["e2"] * S * (S - 1)));
}

// Method to compute Tajima D test (1989)
// Arguments: a PolymorphismSequenceContainer
// Return: Tajima's D (1989)
double SequenceStatistics::tajimaDTNM(const PolymorphismSequenceContainer & psc, bool gapflag) {
	unsigned int eta = totNumberMutations(psc, gapflag);
	double tajima = tajima83(psc, gapflag);
	unsigned int n = psc.getNumberOfSequences();
	map<string, double> values = _getUsefullValues(n);
	double eta_a1 = (double) eta / values["a1"];
	return (tajima - eta_a1) / sqrt((values["e1"] * eta) + (values["e2"] * eta * (eta - 1)));
}

double SequenceStatistics::fuliD(const PolymorphismSequenceContainer & ingroup, const PolymorphismSequenceContainer & outgroup) {
	unsigned int n = ingroup.getNumberOfSequences();
	double nn = (double) n;
	map<string, double> values = _getUsefullValues(n);
	double vD = 1. + (pow(values["a1"], 2) / (values["a2"] + pow(values["a1"], 2))) * (values["cn"] - ((nn + 1.) / (nn - 1.)));
	double uD = values["a1"] - 1. - vD;
	double eta = (double) totNumberMutations(ingroup);
	double etae = (double) countSingleton(outgroup);
	return (eta - values["a1"] * etae) / sqrt((uD * eta) + (vD * eta * eta));
}

double SequenceStatistics::fuliDstar(const PolymorphismSequenceContainer & group) {
	unsigned int n = group.getNumberOfSequences();
	double nn = (double) n;
	map<string, double> values = _getUsefullValues(n);

// Fu & Li 1993
	double _n = nn / (nn - 1.);
	double vDs = (
	               (_n * _n * values["a2"])
	             + (values["a1"] * values["a1"] * values["dn"])
	             - (2. * (nn * values["a1"] * (values["a1"] + 1.) / ((nn - 1.) * (nn - 1.))))
	             )
	             /
		           (pow(values["a1"], 2) + values["a2"]);
	double uDs = _n * (values["a1"] - _n) - vDs;

// Simonsen et al. 1995
/*	double vDs = (
	               (values["a2"] / pow(values["a1"], 2))
	             - (2./nn) * (1. + 1./values["a1"] - values["a1"] + values["a1"]/nn)
	             - 1./(nn*nn)
	             )
	             /
	             (pow(values["a1"], 2) + values["a2"]);
	double uDs = (((nn - 1.)/nn - 1./values["a1"]) / values["a1"]) - vDs;
*/
	double eta = (double) totNumberMutations(group);
	double etas = (double) countSingleton(group);

// Fu & Li 1993
	return ((_n * eta) - (values["a1"] * etas)) / sqrt(uDs * eta + vDs * eta * eta);

// Simonsen et al. 1995
//	return ((eta / values["a1"]) - (etas * ((n - 1) / n))) / sqrt(uDs * eta + vDs * eta * eta);
}

double SequenceStatistics::fuliF(const PolymorphismSequenceContainer & ingroup, const PolymorphismSequenceContainer & outgroup) {
	unsigned int n = ingroup.getNumberOfSequences();
	double nn = (double) n;
	map<string, double> values = _getUsefullValues(n);
	double pi = tajima83(ingroup, true);
	double vF = (values["cn"] + values["b2"] - 2. / (nn - 1.)) / (pow(values["a1"], 2) + values["a2"]);
	double uF = ((1. + values["b1"] - (4. * ((nn + 1.) / ((nn - 1.) * (nn - 1.)))) * (values["a1n"] - (2. * nn) / (nn + 1.))) / values["a1"]) - vF;
	double eta = (double) totNumberMutations(ingroup);
	double etae = (double) countSingleton(outgroup);
	return (pi - etae) / sqrt(uF * eta + vF * eta * eta);
}

double SequenceStatistics::fuliFstar(const PolymorphismSequenceContainer & group) {
	unsigned int n = group.getNumberOfSequences();
	double nn = (double) n;
	map<string, double> values = _getUsefullValues(n);
	double pi = tajima83(group, true);

// Fu & Li 1993
//	double vFs = (values["dn"] + values["b2"] - (2. / (nn - 1.)) * (4. * values["a2"] - 6. + 8. / nn)) / (pow(values["a1"], 2) + values["a2"]);
//	double uFs = (((nn / (nn - 1.)) + values["b1"] - (4. / (nn * (nn - 1.))) + 2. * ((nn + 1.) / (pow((nn - 1.), 2))) * (values["a1n"] - 2. * nn / (nn + 1.))) / values["a1"]) - vFs;

// Simonsen et al. 1995
	double vFs = (((2*nn*nn*nn + 110*nn*nn - 255*nn + 153) / (9*nn*nn*(nn-1))) + ((2*(n-1)*values["a1"]) / (nn*nn)) - 8*values["a2"]/nn) / (pow(values["a1"], 2) + values["a2"]);
	double uFs = (((4*nn*nn + 19*nn + 3 - 12*(nn+1.)*values["a1n"]) / (3*nn*(n-1))) / values["a1"]) - vFs;

	double eta = (double) totNumberMutations(group);
	double etas = (double) countSingleton(group);
// Fu & Li 1993
// Simonsen et al. 1995
	return (pi - ((nn - 1.) / nn * etas)) / sqrt(uFs * eta + vFs * eta * eta);
}


//******************************************************************************************************************
//Linkage disequilibrium statistics
//******************************************************************************************************************



	/**********************/
	/* Preliminary method */
	/**********************/

// Create a PolymorphismSequenceContainer with only polymorphic site and 0 (less frequent) and 1 (more frequent) alleles
// This psc is needed to compute Linkage Disequilibrium Statistics in the class SequenceStatistics
// Should be used before excluding gaps, but sites with gaps are not counted as polymorphic sites
// Singleton can be excluded
// Polymorphix site with the lowest frequency < threshold can be excluded
PolymorphismSequenceContainer * SequenceStatistics::GenerateLDContainer(const PolymorphismSequenceContainer & psc, bool keepsingleton, double freqmin) throw (Exception) {
	try {
                SiteSelection ss;
		// Extract polymorphic site with only two alleles
		for(unsigned int i=0; i<psc.getNumberOfSites(); i++){
			if(keepsingleton) {
				if(SiteTools::isComplete(*psc.getSite(i)) && !SiteTools::isConstant(*psc.getSite(i)) && !SiteTools::isTriplet(*psc.getSite(i))){
					ss.push_back(i);
				}
			}
			else{
				if(SiteTools::isComplete(*psc.getSite(i)) && !SiteTools::isConstant(*psc.getSite(i)) && !SiteTools::isTriplet(*psc.getSite(i)) && !SiteTools::hasSingleton(*psc.getSite(i))){
					ss.push_back(i);
				}
                        }
		}

		const SiteContainer* sc = SiteContainerTools::getSelectedSites(psc,ss);
                Alphabet* alpha = new DNA();
		PolymorphismSequenceContainer *ldpsc = new PolymorphismSequenceContainer(sc->getNumberOfSequences(),alpha);
		// Assign 1 to the more frequent and 0 to the less frequent alleles
		for(unsigned int i=0; i<sc->getNumberOfSites(); i++){
			const Site* site = sc->getSite(i);
			bool deletesite = false;
			map<int, double> freqs = SymbolListTools::getFrequencies(*site);
			for(unsigned int j=0; j<sc->getNumberOfSequences(); j++){
				if(freqs[site->getValue(j)]>=0.5){
					if(freqs[site->getValue(j)]<1-freqmin) site->setElement(j,1);
					else deletesite = true;
				}
				else site->setElement(j,0);
			}
                        if(!deletesite)	ldpsc->addSite(*site);
		}
                delete alpha;
		return ldpsc;
		}
	catch(...) {}

}


	/*************************************/
	/* Pairwise LD and distance measures */
	/*************************************/

// Return a vector with the pairwise distances between site positions corresponding to a LD PolymorphismSequenceContainer
// All sequences are supposed to have the same length
Vdouble SequenceStatistics::PairwiseDistances1(const PolymorphismSequenceContainer & psc,bool keepsingleton, double freqmin){
	//get Positions with sites of interest
	SiteSelection ss;
	for(unsigned int i=0; i<psc.getNumberOfSites(); i++){
		if(keepsingleton) {
			if(SiteTools::isComplete(*psc.getSite(i)) && !SiteTools::isConstant(*psc.getSite(i)) && !SiteTools::isTriplet(*psc.getSite(i))){
				const Site* site = psc.getSite(i);
				bool deletesite = false;
				map<int, double> freqs = SymbolListTools::getFrequencies(*site);
				for(unsigned int j=0; j<site->getAlphabet()->getSize(); j++){
					if(freqs[j]>=1-freqmin) deletesite = true;
				}
				if(!deletesite) ss.push_back(i);
			}
		}
		else{
			if(SiteTools::isComplete(*psc.getSite(i)) && !SiteTools::isConstant(*psc.getSite(i)) && !SiteTools::isTriplet(*psc.getSite(i)) && !SiteTools::hasSingleton(*psc.getSite(i))){
				ss.push_back(i);
				const Site* site = psc.getSite(i);
				bool deletesite = false;
				map<int, double> freqs = SymbolListTools::getFrequencies(*site);
				for(unsigned int j=0; j<site->getAlphabet()->getSize(); j++){
					if(freqs[j]>=1-freqmin) deletesite = true;
				}
				if(!deletesite) ss.push_back(i);
			}
                }
	}
	//compute pairwise distances
	Vdouble dist;
        if(ss.size()==0) return dist;
	for(unsigned int i=0; i<ss.size()-1; i++){
		for(unsigned int j=i+1; j<ss.size(); j++){
			dist.push_back(ss[j]-ss[i]);
		}
	}
	return dist;
}



// Return a vector with all the pairwise distances between two sites corresponding to a LD PolymorphismSequenceContainer
// This method take into account the fact that sequences may differ by their number of gaps
// Pairwise distance are computed for each sequence. The mean pairwise distance is then computed.
Vdouble SequenceStatistics::PairwiseDistances2(const PolymorphismSequenceContainer & psc, bool keepsingleton, double freqmin){
	SiteSelection ss;
	for(unsigned int i=0; i<psc.getNumberOfSites(); i++){
		if(keepsingleton) {
			if(SiteTools::isComplete(*psc.getSite(i)) && !SiteTools::isConstant(*psc.getSite(i)) && !SiteTools::isTriplet(*psc.getSite(i))){
				const Site* site = psc.getSite(i);
				bool deletesite = false;
				map<int, double> freqs = SymbolListTools::getFrequencies(*site);
				for(unsigned int j=0; j<site->getAlphabet()->getSize(); j++){
					if(freqs[j]>=1-freqmin) deletesite = true;
				}
				if(!deletesite) ss.push_back(i);
			}
		}
		else{
			if(SiteTools::isComplete(*psc.getSite(i)) && !SiteTools::isConstant(*psc.getSite(i)) && !SiteTools::isTriplet(*psc.getSite(i)) && !SiteTools::hasSingleton(*psc.getSite(i))){
				ss.push_back(i);
				const Site* site = psc.getSite(i);
				bool deletesite = false;
				map<int, double> freqs = SymbolListTools::getFrequencies(*site);
				for(unsigned int j=0; j<site->getAlphabet()->getSize(); j++){
					if(freqs[j]>=1-freqmin) deletesite = true;
				}
				if(!deletesite) ss.push_back(i);
			}
                }
	}
	unsigned int n = ss.size();
	Vdouble distance(n*(n-1)/2,0);
        if(n==0) return distance;
	unsigned int nbsite = psc.getNumberOfSites();
	for(unsigned int k=0; k<psc.getNumberOfSequences(); k++){
		const Sequence* seq = psc.getSequence(k);
		SiteSelection gap, newss = ss;
                Vdouble dist;
		for(unsigned int i=0; i<nbsite; i++){
			if(seq->getValue(i)==-1) gap.push_back(i);
		}
		//Site positions are re-numbered to take gaps into account
		for(unsigned int i=0; i<gap.size(); i++){
			for(unsigned int j=0; j<ss.size(); j++){
				if(ss[j]>gap[i]) newss[j]--;
			}
		}
		for(unsigned int i=0; i<n-1; i++){
			for(unsigned int j=i+1; j<n; j++){
				dist.push_back(newss[j]-newss[i]);
			}
		}
		distance += dist;
	}
	distance = distance/psc.getNumberOfSequences();
	return distance;
}

// Return a vector with all pairwise |D| measures between 2 sites (Lewontin & Kojima 1964)
Vdouble SequenceStatistics::PairwiseD(const PolymorphismSequenceContainer & psc, bool keepsingleton, double freqmin) {
	PolymorphismSequenceContainer* newpsc = SequenceStatistics::GenerateLDContainer(psc, keepsingleton,  freqmin);
	Vdouble D;
	unsigned int nbsite = newpsc->getNumberOfSites();
	unsigned int nbseq = newpsc->getNumberOfSequences();
        if(nbsite==0) return D;
	for(unsigned int i=0; i<nbsite-1; i++){
		for(unsigned int j=i+1; j<nbsite; j++){
			double haplo=0;
			const Site* site1 = newpsc->getSite(i);
			const Site* site2 = newpsc->getSite(j);
			map<int,double> freq1 = SymbolListTools::getFrequencies(*site1);
			map<int,double> freq2 = SymbolListTools::getFrequencies(*site2);
			for(unsigned int k=0; k<nbseq; k++){
				if(site1->getValue(k) + site2->getValue(k)==2) haplo++;
			}
			haplo = haplo/nbseq;
			D.push_back(std::abs(haplo-freq1[1]*freq2[1]));
		}
	}
	return D;
}



// Return a vector with all pairwise |D'| measures between 2 sites (Lewontin 1964)
Vdouble SequenceStatistics::PairwiseDprime(const PolymorphismSequenceContainer & psc, bool keepsingleton, double freqmin) {
	PolymorphismSequenceContainer* newpsc = SequenceStatistics::GenerateLDContainer(psc, keepsingleton, freqmin);
	Vdouble Dprime;
	unsigned int nbsite = newpsc->getNumberOfSites();
	unsigned int nbseq = newpsc->getNumberOfSequences();
        if(nbsite==0) return Dprime;
	for(unsigned int i=0; i<nbsite-1; i++){
		for(unsigned int j=i+1; j<nbsite; j++){
			double haplo=0;
			const Site* site1 = newpsc->getSite(i);
			const Site* site2 = newpsc->getSite(j);
			map<int,double> freq1 = SymbolListTools::getFrequencies(*site1);
			map<int,double> freq2 = SymbolListTools::getFrequencies(*site2);
			for(unsigned int k=0; k<nbseq; k++){
				if(site1->getValue(k) + site2->getValue(k)==2) haplo++;
			}
			haplo = haplo/nbseq;
			double d, D = (haplo-freq1[1]*freq2[1]);
			if(D>0){
				if(freq1[1]*freq2[0]<=freq1[0]*freq2[1]){
					d=std::abs(D)/(freq1[1]*freq2[0]);
				}
				else{
					d=std::abs(D)/(freq1[0]*freq2[1]);
				}
			}
			else{
				if(freq1[1]*freq2[1]<=freq1[0]*freq2[0]){
					d=std::abs(D)/(freq1[1]*freq2[1]);
				}
				else{
					d=std::abs(D)/(freq1[0]*freq2[0]);
				}
			}
			Dprime.push_back(d);
		}
	}
	return Dprime;
}


// Return a vector with all pairwise R� measures between 2 sites (Hill & Robertson 1968)
Vdouble SequenceStatistics::PairwiseR2(const PolymorphismSequenceContainer & psc, bool keepsingleton, double freqmin) {
	PolymorphismSequenceContainer* newpsc = SequenceStatistics::GenerateLDContainer(psc, keepsingleton, freqmin);
	Vdouble R2;
	unsigned int nbsite = newpsc->getNumberOfSites();
        if(nbsite==0) return R2;
	unsigned int nbseq = newpsc->getNumberOfSequences();
	for(unsigned int i=0; i<nbsite-1; i++){
		for(unsigned int j=i+1; j<nbsite; j++){
			double haplo=0;
			const Site* site1 = newpsc->getSite(i);
			const Site* site2 = newpsc->getSite(j);
			map<int,double> freq1 = SymbolListTools::getFrequencies(*site1);
			map<int,double> freq2 = SymbolListTools::getFrequencies(*site2);
			for(unsigned int k=0; k<nbseq; k++){
				if(site1->getValue(k) + site2->getValue(k)==2) haplo++;
			}
			haplo = haplo/nbseq;
			double r = ((haplo-freq1[1]*freq2[1])*(haplo-freq1[1]*freq2[1]))/(freq1[0]*freq1[1]*freq2[0]*freq2[1]);
			R2.push_back(r);
		}
	}
	return R2;
}




	/***********************************/
	/* Global LD and distance measures */
	/***********************************/


//Return the mean D over all pairwise comparisons
double SequenceStatistics::MeanD(const PolymorphismSequenceContainer & psc, bool keepsingleton, double freqmin){
	Vdouble D = SequenceStatistics::PairwiseD(psc,keepsingleton,freqmin);
	return mean(D);
}

//Return the mean D' over all pairwise comparisons
double SequenceStatistics::MeanDprime(const PolymorphismSequenceContainer & psc, bool keepsingleton, double freqmin){
	Vdouble Dprime = SequenceStatistics::PairwiseDprime(psc,keepsingleton,freqmin);
	return mean(Dprime);
}

//Return the mean R� over all pairwise comparisons
double SequenceStatistics::MeanR2(const PolymorphismSequenceContainer & psc, bool keepsingleton, double freqmin){
	Vdouble R2 = SequenceStatistics::PairwiseR2(psc,keepsingleton,freqmin);
	return mean(R2);
}

//Return the mean pairwise distances between sites / method 1: differences between sequences are not taken into account
double SequenceStatistics::MeanDistance1(const PolymorphismSequenceContainer & psc, bool keepsingleton, double freqmin){
	Vdouble dist = PairwiseDistances1(psc,keepsingleton,freqmin);
	return mean(dist);
}

//Return the mean pairwise distances between sites / method 2: differences between sequences are taken into account
double SequenceStatistics::MeanDistance2(const PolymorphismSequenceContainer & psc, bool keepsingleton, double freqmin){
	Vdouble dist = SequenceStatistics::PairwiseDistances2(psc,keepsingleton,freqmin);
	return mean(dist);
}

	/**********************/
	/* Regression methods */
	/**********************/


// Return the slope,a, of the regression |D| = 1+a*distance
// The slope is given in |D'| per kb
// Distance1 or distance2 are chose through the boolean distance1 (false by default)
double SequenceStatistics::OriginRegressionD(const PolymorphismSequenceContainer & psc, bool distance1, bool keepsingleton, double freqmin){
        Vdouble D = SequenceStatistics::PairwiseD(psc,keepsingleton,freqmin)-1;
        Vdouble dist;
        if(distance1) dist = PairwiseDistances1(psc,keepsingleton,freqmin)/1000;
        else  dist = PairwiseDistances2(psc,keepsingleton,freqmin)/1000;
        return sum(D*dist)/sum(dist*dist);
}


// Return the slope of the regression |D'| = 1+a*distance
// The slope is given in |D'| per kb
// Distance1 or distance2 are chose through the boolean distance1 (false by default)
double SequenceStatistics::OriginRegressionDprime(const PolymorphismSequenceContainer & psc, bool distance1, bool keepsingleton, double freqmin){
        Vdouble Dprime = SequenceStatistics::PairwiseDprime(psc,keepsingleton,freqmin)-1;
        Vdouble dist;
        if(distance1) dist = PairwiseDistances1(psc,keepsingleton,freqmin)/1000;
        else  dist = PairwiseDistances2(psc,keepsingleton,freqmin)/1000;
        return sum(Dprime*dist)/sum(dist*dist);
}

// Return the slope of the regression R� = 1+a*distance
// The slope is given in R� per kb
// Distance1 or distance2 are chose through the boolean distance1 (false by default)
double SequenceStatistics::OriginRegressionR2(const PolymorphismSequenceContainer & psc, bool distance1, bool keepsingleton, double freqmin){
        Vdouble R2 = SequenceStatistics::PairwiseR2(psc,keepsingleton,freqmin)-1;
        Vdouble dist;
        if(distance1) dist = PairwiseDistances1(psc,keepsingleton,freqmin)/1000;
        else  dist = PairwiseDistances2(psc,keepsingleton,freqmin)/1000;
        return sum(R2*dist)/sum(dist*dist);
}

// Return the slope and the origin of the regression |D| = a*distance + b
// The slope is given in |D| per kb
// Distance1 or distance2 are chose through the boolean distance1 (false by default)
Vdouble SequenceStatistics::LinearRegressionD(const PolymorphismSequenceContainer & psc, bool distance1, bool keepsingleton, double freqmin){
        Vdouble D = SequenceStatistics::PairwiseD(psc,keepsingleton,freqmin);
        Vdouble dist;
        Vdouble reg(2);
        if(distance1) dist = PairwiseDistances1(psc,keepsingleton,freqmin)/1000;
        else  dist = PairwiseDistances2(psc,keepsingleton,freqmin)/1000;
        reg[0]=cov(dist,D)/var(dist);
        reg[1]=mean(D)-reg[0]*mean(dist);
        return reg;
}

// Return the slope and the origin of the regression |D'| = a*distance + b
// The slope is given in |D'| per kb
// Distance1 or distance2 are chose through the boolean distance1 (false by default)
Vdouble SequenceStatistics::LinearRegressionDprime(const PolymorphismSequenceContainer & psc, bool distance1, bool keepsingleton, double freqmin){
        Vdouble Dprime = SequenceStatistics::PairwiseDprime(psc,keepsingleton,freqmin);
        Vdouble dist;
        Vdouble reg(2);
        if(distance1) dist = PairwiseDistances1(psc,keepsingleton,freqmin)/1000;
        else  dist = PairwiseDistances2(psc,keepsingleton,freqmin)/1000;
        reg[0]=cov(dist,Dprime)/var(dist);
        reg[1]=mean(Dprime)-reg[0]*mean(dist);
        return reg;
}

// Return the slope and the origin of the regression R� = a*distance + b
// The slope is given in R� per kb
// Distance1 or distance2 are chose through the boolean distance1 (false by default)
Vdouble SequenceStatistics::LinearRegressionR2(const PolymorphismSequenceContainer & psc, bool distance1, bool keepsingleton, double freqmin){
        Vdouble R2 = SequenceStatistics::PairwiseR2(psc,keepsingleton,freqmin);
        Vdouble dist;
        Vdouble reg(2);
        if(distance1) dist = PairwiseDistances1(psc,keepsingleton,freqmin)/1000;
        else  dist = PairwiseDistances2(psc,keepsingleton,freqmin)/1000;
        reg[0]=cov(dist,R2)/var(dist);
        reg[1]=mean(R2)-reg[0]*mean(dist);
        return reg;
}


// Return the slope the regression R� = 1/(1+a*distance)
// To fit the theoretical expectation R�=1/(1+4Nr)
// The slope is given in R� per kb
// Distance1 or distance2 are chose through the boolean distance1 (false by default)
double SequenceStatistics::InverseRegressionR2(const PolymorphismSequenceContainer & psc, bool distance1, bool keepsingleton, double freqmin){
        Vdouble R2 = SequenceStatistics::PairwiseR2(psc,keepsingleton,freqmin);
        Vdouble unit(R2.size(),1);
        Vdouble R2transformed = unit/R2 -1;
        Vdouble dist;
        if(distance1) dist = PairwiseDistances1(psc,keepsingleton,freqmin)/1000;
        else  dist = PairwiseDistances2(psc,keepsingleton,freqmin)/1000;
        return sum(R2transformed*dist)/sum(dist*dist);
}


//******************************************************************************************************************
//Private methods
//******************************************************************************************************************

unsigned int SequenceStatistics::_getMutationNumber(const Site & site) {
	unsigned int tmp_count = 0;
	map<int, unsigned int> states_count = SymbolListTools::getCounts(site);
	for (map<int, unsigned int>::iterator it = states_count.begin() ; it != states_count.end() ; it++)
		if (it->first >= 0)
			tmp_count++;
	if (tmp_count > 0)
		tmp_count--;
	return tmp_count;
}

unsigned int SequenceStatistics::_getSingletonNumber(const Site & site) {
	unsigned int nus = 0;
	map<int, unsigned int> states_count = SymbolListTools::getCounts(site);
	for (map<int, unsigned int>::iterator it = states_count.begin() ; it != states_count.end() ; it++)
		if (it->second == 1)
			nus++;
	return nus;
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
		double nn = (double) n;
		values["a1n"] = values["a1"] + (1. / nn);
		values["b1"] = (nn + 1.) / (3. * (nn - 1.));
		values["b2"] = 2. * ((nn * nn) + nn + 3.) / (9. * nn * (nn - 1.));
		values["c1"] = values["b1"] - (1. / values["a1"]);
		values["c2"] = values["b2"] - ((nn + 2.) / (values["a1"] * nn)) + (values["a2"] / (values["a1"] * values["a1"]));
		values["cn"] = 2. * ((nn * values["a1"]) - (2. * (nn - 1.))) / ((nn - 1.) * (nn - 2.));
		values["dn"] = values["cn"] + ((nn - 2.) / ((nn - 1.) * (nn - 1.))) + ((2. / (nn - 1.)) * ((3. / 2.) - (((2. * values["a1n"]) - 3.) / (nn - 2.)) - (1. / nn)));
		values["e1"] = values["c1"] / values["a1"];
		values["e2"] = values["c2"] / ((values["a1"] * values["a1"]) + values["a2"]);
	}
	return values;
}
