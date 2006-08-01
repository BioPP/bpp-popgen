/*
 * File MultilocusGenotypeStatistics.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Wednesday August 04 2004
 *
*/
/*
Copyright or © or Copr. CNRS, (November 17, 2004)


This software is a computer program whose purpose is to provide classes
for population genetics analysis.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/
#include "MultilocusGenotypeStatistics.h"
#include "PolymorphismMultiGContainerTools.h"

#include <iostream>
#include <math.h>

using namespace std;

vector<unsigned int> MultilocusGenotypeStatistics::getAllelesIdsForGroups(const PolymorphismMultiGContainer & pmgc, unsigned int locus_position, const set<unsigned int> & groups) throw (IndexOutOfBoundsException) {
	map<unsigned int, unsigned int> tmp_alleles;
	try {
		tmp_alleles = getAllelesMapForGroups(pmgc, locus_position, groups);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("MultilocusGenotypeStatistics::getAllelesIdsForGroups: locus_position out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
	return MapTools::getKeys(tmp_alleles);
}

unsigned int MultilocusGenotypeStatistics::countGametesForGroups(const PolymorphismMultiGContainer & pmgc, unsigned int locus_position, const set<unsigned int> & groups) throw (IndexOutOfBoundsException) {
	map<unsigned int, unsigned int> allele_count;
	unsigned int nb_tot_allele = 0;
	try {
		allele_count = getAllelesMapForGroups(pmgc, locus_position, groups);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("MultilocusGenotypeStatistics::countGametesForGroups: locus_position out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
	vector<unsigned int> counter = MapTools::getValues(allele_count);
	for (unsigned int i = 0 ; i < counter.size() ; i++)
		nb_tot_allele += counter[i];
	return nb_tot_allele;
}

map<unsigned int, unsigned int> MultilocusGenotypeStatistics::getAllelesMapForGroups(const PolymorphismMultiGContainer & pmgc, unsigned int locus_position, const set<unsigned int> & groups) throw (IndexOutOfBoundsException) {
	map<unsigned int, unsigned int> alleles_count;
	for (unsigned int i = 0 ; i < pmgc.size() ; i++) {
		try {

             if (! pmgc.getMultilocusGenotype(i)->isMonolocusGenotypeMissing(locus_position) &&  (groups.find(pmgc.getGroupId(i)) != groups.end()) ) {
			//if (! pmgc.getMultilocusGenotype(i)->isMonolocusGenotypeMissing(locus_position) &&  (find(groups.begin(), groups.end(), pmgc.getGroupId(i)) != groups.end()) ) {
				vector<unsigned int> tmp_alleles = pmgc.getMultilocusGenotype(i)->getMonolocusGenotype(locus_position)->getAlleleIndex();
				for (unsigned int j = 0 ; j < tmp_alleles.size() ; j++) alleles_count[tmp_alleles[j]]++;
			}
		}
		catch (IndexOutOfBoundsException & ioobe) {
			throw IndexOutOfBoundsException("MultilocusGenotypeStatistics::getAllelesMapForGroups: locus_position out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
		}
	}
	return alleles_count;
}

map<unsigned int, double> MultilocusGenotypeStatistics::getAllelesFrqForGroups(const PolymorphismMultiGContainer & pmgc, unsigned int locus_position, const set<unsigned int> & groups) throw (Exception) {
	map<unsigned int, double> alleles_frq;
	unsigned int nb_tot_allele = 0;
	map<unsigned int, unsigned int> tmp_alleles;
	try {
		tmp_alleles = getAllelesMapForGroups(pmgc, locus_position, groups);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("MultilocusGenotypeStatistics::getAllelesFrqForGroups: locus_position out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
	vector<unsigned int> counter = MapTools::getValues(tmp_alleles);
	for (unsigned int i = 0 ; i < counter.size() ; i++)
		nb_tot_allele += counter[i];
	if (nb_tot_allele == 0)
		throw ZeroDivisionException("MultilocusGenotypeStatistics::getAllelesFrqForGroups.");
	for (map<unsigned int, unsigned int>::iterator it = tmp_alleles.begin() ; it != tmp_alleles.end() ; it++)
		alleles_frq[it->first] = (double) it->second / (double) nb_tot_allele;
	return alleles_frq;
}

unsigned int MultilocusGenotypeStatistics::countNonMissingForGroups(const PolymorphismMultiGContainer & pmgc, unsigned int locus_position, const set<unsigned int> & groups) throw (IndexOutOfBoundsException) {
	unsigned int counter = 0;
	for (unsigned int i = 0 ; i < pmgc.size() ; i++) {
		try {
			//if (! pmgc.getMultilocusGenotype(i)->isMonolocusGenotypeMissing(locus_position) &&  (find(groups.begin(), groups.end(), pmgc.getGroupId(i)) != groups.end()) )
			 if (! pmgc.getMultilocusGenotype(i)->isMonolocusGenotypeMissing(locus_position) &&  (groups.find(pmgc.getGroupId(i) ) != groups.end()) )
				counter++;
		}
		catch (IndexOutOfBoundsException & ioobe) {
			throw IndexOutOfBoundsException("MultilocusGenotypeStatistics::countNonMissing: locus_position out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
		}
	}
	return counter;
}

unsigned int MultilocusGenotypeStatistics::countBiAllelicForGroups(const PolymorphismMultiGContainer & pmgc, unsigned int locus_position, const set<unsigned int> & groups) throw (IndexOutOfBoundsException) {
	unsigned int counter = 0;
	for (unsigned int i = 0 ; i < pmgc.size() ; i++) {
		try {
			//if (! pmgc.getMultilocusGenotype(i)->isMonolocusGenotypeMissing(locus_position) && (find(groups.begin(), groups.end(), pmgc.getGroupId(i)) != groups.end()) )
			 if (! pmgc.getMultilocusGenotype(i)->isMonolocusGenotypeMissing(locus_position) && (groups.find(pmgc.getGroupId(i)) != groups.end()) )
				if ((pmgc.getMultilocusGenotype(i)->getMonolocusGenotype(locus_position)->getAlleleIndex()).size() == 2)
					counter++;
		}
		catch (IndexOutOfBoundsException & ioobe) {
			throw IndexOutOfBoundsException("MultilocusGenotypeStatistics::countBiAllelic: locus_position out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
		}
	}
	return counter;
}

map<unsigned int, unsigned int> MultilocusGenotypeStatistics::countHeterozygousForGroups(const PolymorphismMultiGContainer & pmgc, unsigned int locus_position, const set<unsigned int> & groups) throw (IndexOutOfBoundsException) {
	map<unsigned int, unsigned int> counter;
	for (unsigned int i = 0 ; i < pmgc.size() ; i++) {
		try {
			if (! pmgc.getMultilocusGenotype(i)->isMonolocusGenotypeMissing(locus_position) && (groups.find(pmgc.getGroupId(i)) != groups.end() )) {
				const MonolocusGenotype * tmp_mg = pmgc.getMultilocusGenotype(i)->getMonolocusGenotype(locus_position);
				if ((tmp_mg->getAlleleIndex()).size() == 2) {
					if (! dynamic_cast<const BiAlleleMonolocusGenotype *>(tmp_mg)->isHomozygous()) {
						vector<unsigned int> tmp_alleles = tmp_mg->getAlleleIndex();
						for (unsigned int j = 0 ; j < tmp_alleles.size() ; j++) counter[tmp_alleles[j]]++;
					}
				}
			}
		}
		catch (IndexOutOfBoundsException & ioobe) {
			throw IndexOutOfBoundsException("MultilocusGenotypeStatistics::countHeterozygous: locus_position out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
		}
	}
	return counter;
}

map<unsigned int, double> MultilocusGenotypeStatistics::getHeterozygousFrqForGroups(const PolymorphismMultiGContainer & pmgc, unsigned int locus_position, const set<unsigned int> & groups) throw (Exception) {
	map<unsigned int, double> freq;
	unsigned int counter = 0;
	for (unsigned int i = 0 ; i < pmgc.size() ; i++) {
		try {
			if (! pmgc.getMultilocusGenotype(i)->isMonolocusGenotypeMissing(locus_position) && (groups.find(pmgc.getGroupId(i)) != groups.end()) ) {
				const MonolocusGenotype * tmp_mg = pmgc.getMultilocusGenotype(i)->getMonolocusGenotype(locus_position);
				if ((tmp_mg->getAlleleIndex()).size() == 2) {
					counter++;
					if (! dynamic_cast<const BiAlleleMonolocusGenotype *>(tmp_mg)->isHomozygous()) {
						vector<unsigned int> tmp_alleles = tmp_mg->getAlleleIndex();
						for (unsigned int j = 0 ; j < tmp_alleles.size() ; j++) freq[tmp_alleles[j]]++;
					}
				}
			}
		}
		catch (IndexOutOfBoundsException & ioobe) {
			throw IndexOutOfBoundsException("MultilocusGenotypeStatistics::getHeterozygousFrqForGroups: locus_position out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
		}
	}
	if (counter == 0)
		throw ZeroDivisionException("MultilocusGenotypeStatistics::getHeterozygousFrqForGroups.");
	for (map<unsigned int, double>::iterator i = freq.begin() ; i != freq.end() ; i++)
		i->second = (double) i->second / (double) counter;
	return freq;
}

double MultilocusGenotypeStatistics::getHobsForGroups(const PolymorphismMultiGContainer & pmgc, unsigned int locus_position, const set<unsigned int> & groups) throw (Exception) {
	map<unsigned int, double> heterozygous_frq;
	double frq = 0.;
	try {
		heterozygous_frq = getHeterozygousFrqForGroups(pmgc, locus_position, groups);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("MultilocusGenotypeStatistics::getHobsForGroups: locus_position out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
	catch (ZeroDivisionException & zde) {
		throw ZeroDivisionException("MultilocusGenotypeStatistics::getHobsForGroups.");
	}
	for (map<unsigned int, double>::iterator it = heterozygous_frq.begin() ; it != heterozygous_frq.end() ; it++)
		frq += it->second;
	return frq / heterozygous_frq.size();
}

double MultilocusGenotypeStatistics::getHexpForGroups(const PolymorphismMultiGContainer & pmgc, unsigned int locus_position, const set<unsigned int> & groups) throw (Exception) {
	map<unsigned int, double> allele_frq;
	double frqsqr = 0.;
	try {
		allele_frq = getAllelesFrqForGroups(pmgc, locus_position, groups);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("MultilocusGenotypeStatistics::getHexpForGroups: locus_position out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
	catch (ZeroDivisionException & zde) {
		throw ZeroDivisionException("MultilocusGenotypeStatistics::getHexpForGroups.");
	}
	for (map<unsigned int, double>::iterator it = allele_frq.begin() ; it != allele_frq.end() ; it++)
		frqsqr += it->second * it->second;
	return (1 - frqsqr);
}

double MultilocusGenotypeStatistics::getHnbForGroups(const PolymorphismMultiGContainer & pmgc, unsigned int locus_position, const set<unsigned int> & groups) throw (Exception) {
	unsigned int nb_alleles;
	double Hexp;
	try {
		nb_alleles = countGametesForGroups(pmgc, locus_position, groups);
		Hexp = getHexpForGroups(pmgc, locus_position, groups);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("MultilocusGenotypeStatistics::getHnbForGroups: locus_position out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
	catch (ZeroDivisionException & zde) {
		throw ZeroDivisionException("MultilocusGenotypeStatistics::getHnbForGroups.");
	}
	return (2 * nb_alleles * Hexp  / ((2 * nb_alleles) - 1));
}

double MultilocusGenotypeStatistics::getDnei72(const PolymorphismMultiGContainer & pmgc, vector<unsigned int> locus_positions, unsigned int grp1, unsigned int grp2) throw (Exception) {
	map<unsigned int, double> allele_frq1, allele_frq2;
	vector<unsigned int> allele_ids;
	set<unsigned int> group1_id;
	set<unsigned int> group2_id;
	set<unsigned int> groups_id;
	double Jx = 0.;
	double Jy = 0.;
	double Jxy = 0.;
	group1_id.insert(grp1);
	group2_id.insert(grp2);
	groups_id.insert(grp1);
	groups_id.insert(grp2);
	for (unsigned int i = 0 ; i < locus_positions.size() ; i++) {
		allele_ids.clear();
		allele_frq1.clear();
		allele_frq2.clear();
		try {
			allele_ids = getAllelesIdsForGroups(pmgc, locus_positions[i], groups_id);
			allele_frq1 = getAllelesFrqForGroups(pmgc, locus_positions[i], group1_id);
			allele_frq2 = getAllelesFrqForGroups(pmgc, locus_positions[i], group2_id);
		}
		catch (Exception & e) {
			throw e;
		}
		for (unsigned int j = 0 ; j < allele_ids.size() ; j++) {
			map<unsigned int, double>::iterator it1 = allele_frq1.find(allele_ids[j]);
			map<unsigned int, double>::iterator it2 = allele_frq2.find(allele_ids[j]);
			double tmp_frq1 = (it1 != allele_frq1.end()) ? it1->second : 0.;
			double tmp_frq2 = (it2 != allele_frq2.end()) ? it2->second : 0.;
			Jx += tmp_frq1 * tmp_frq1;
			Jy += tmp_frq2 * tmp_frq2;
			Jxy += tmp_frq1 * tmp_frq2;
		}
	}
	if (Jx * Jy == 0.)
		throw ZeroDivisionException("MultilocusGenotypeStatistics::getDnei72.");
	return -log(Jxy / sqrt(Jx * Jy));
}

double MultilocusGenotypeStatistics::getDnei78(const PolymorphismMultiGContainer & pmgc, vector<unsigned int> locus_positions, unsigned int grp1, unsigned int grp2) throw (Exception) {
	map<unsigned int, double> allele_frq1, allele_frq2;
	vector<unsigned int> allele_ids;
	set<unsigned int> group1_id;
	set<unsigned int> group2_id;
	set<unsigned int> groups_id;
	double Jx = 0.;
	double Jy = 0.;
	double Jxy = 0.;
	unsigned int nx = 0, ny = 0;
	group1_id.insert(grp1);
	group2_id.insert(grp2);
	groups_id.insert(grp1);
	groups_id.insert(grp2);
	for (unsigned int i = 0 ; i < locus_positions.size() ; i++) {
		allele_ids.clear();
		allele_frq1.clear();
		allele_frq2.clear();
		try {
			allele_ids = getAllelesIdsForGroups(pmgc, locus_positions[i], groups_id);
			allele_frq1 = getAllelesFrqForGroups(pmgc, locus_positions[i], group1_id);
			allele_frq2 = getAllelesFrqForGroups(pmgc, locus_positions[i], group2_id);
			nx = countBiAllelicForGroups(pmgc, locus_positions[i], group1_id);
			ny = countBiAllelicForGroups(pmgc, locus_positions[i], group2_id);
		}
		catch (Exception & e) {
			throw e;
		}
		double tmp_Jx = 0.;
		double tmp_Jy = 0.;
		for (unsigned int j = 0 ; j < allele_ids.size() ; j++) {
			map<unsigned int, double>::iterator it1 = allele_frq1.find(allele_ids[j]);
			map<unsigned int, double>::iterator it2 = allele_frq2.find(allele_ids[j]);
			double tmp_frq1 = (it1 != allele_frq1.end()) ? it1->second : 0.;
			double tmp_frq2 = (it2 != allele_frq2.end()) ? it2->second : 0.;
			tmp_Jx += tmp_frq1 * tmp_frq1;
			tmp_Jy += tmp_frq2 * tmp_frq2;
			Jxy += tmp_frq1 * tmp_frq2;
		}
		Jx += ((2. * (double) nx * tmp_Jx) - 1.) / ((2. * (double) nx) - 1.);
		Jy += ((2. * (double) ny * tmp_Jy) - 1.) / ((2. * (double) ny) - 1.);
	}
	double denom = Jx * Jy;
	if (denom == 0.)
		throw ZeroDivisionException("MultilocusGenotypeStatistics::getDnei78.");
	return -log(Jxy / sqrt(denom));
}

map<unsigned int, MultilocusGenotypeStatistics::Fstats> MultilocusGenotypeStatistics::getAllelesFstats(const PolymorphismMultiGContainer & pmgc, unsigned int locus_position, const set<unsigned int> & groups) throw (Exception) {
	map<unsigned int, MultilocusGenotypeStatistics::VarComp> vc = getVarianceComponents(pmgc, locus_position, groups);
	map<unsigned int, MultilocusGenotypeStatistics::Fstats> f_stats;
	for (map<unsigned int, MultilocusGenotypeStatistics::VarComp>::iterator it = vc.begin() ; it != vc.end() ; it++) {
        double abc = it->second.a + it->second.b + it->second.c;
        double bc  = it->second.b + it->second.c;

		if (abc == 0)
        {
                f_stats[it->first].Fit = NAN;
                f_stats[it->first].Fst = NAN;
        }
        {
                f_stats[it->first].Fit = 1. - it->second.c / abc;
        		f_stats[it->first].Fst = it->second.a / abc;
        }
        if (bc == 0) f_stats[it->first].Fis = NAN;
		else f_stats[it->first].Fis = 1. - it->second.c / bc;
	}
	return f_stats;
}

map<unsigned int, double> MultilocusGenotypeStatistics::getAllelesFit(const PolymorphismMultiGContainer & pmgc, unsigned int locus_position, const set<unsigned int> & groups) throw (Exception) {
	map<unsigned int, MultilocusGenotypeStatistics::VarComp> values = getVarianceComponents(pmgc, locus_position, groups);
	map<unsigned int, double> Fit;
	for (map<unsigned int, MultilocusGenotypeStatistics::VarComp>::iterator it = values.begin() ; it != values.end() ; it++) {
		Fit[it->first] = it->second.a + it->second.b + it->second.c;
		if (Fit[it->first] == 0.)
			throw ZeroDivisionException("MultilocusGenotypeStatistics::getAllelesFit.");
		Fit[it->first] = 1. - it->second.c / Fit[it->first];
	}
	return Fit;
}

map<unsigned int, double> MultilocusGenotypeStatistics::getAllelesFst(const PolymorphismMultiGContainer & pmgc, unsigned int locus_position, const set<unsigned int> & groups) throw (Exception) {
	if (groups.size() <= 1)
		throw BadIntegerException("MultilocusGenotypeStatistics::getAllelesFst: groups must be >= 2.", groups.size());
	map<unsigned int, MultilocusGenotypeStatistics::VarComp> values = getVarianceComponents(pmgc, locus_position, groups);
	map<unsigned int, double> Fst;
	for (map<unsigned int, MultilocusGenotypeStatistics::VarComp>::iterator it = values.begin() ; it != values.end() ; it++) {
		Fst[it->first] = it->second.a + it->second.b + it->second.c;
		if (Fst[it->first] == 0.)
			throw ZeroDivisionException("MultilocusGenotypeStatistics::getAllelesFst.");
		Fst[it->first] = it->second.a / Fst[it->first];
	}
	return Fst;
}



map<unsigned int, double> MultilocusGenotypeStatistics::getAllelesFis(const PolymorphismMultiGContainer & pmgc, unsigned int locus_position, const set<unsigned int> & groups) throw (Exception) {
	map<unsigned int, MultilocusGenotypeStatistics::VarComp> values = getVarianceComponents(pmgc, locus_position, groups);
	map<unsigned int, double> Fis;
	for (map<unsigned int, MultilocusGenotypeStatistics::VarComp>::iterator it = values.begin() ; it != values.end() ; it++) {
		Fis[it->first] = it->second.b + it->second.c;
		if (Fis[it->first] == 0.)
			throw ZeroDivisionException("MultilocusGenotypeStatistics::getAllelesFis.");
		Fis[it->first] = 1. - it->second.c / Fis[it->first];
	}
	return Fis;
}

map<unsigned int, MultilocusGenotypeStatistics::VarComp> MultilocusGenotypeStatistics::getVarianceComponents(const PolymorphismMultiGContainer & pmgc, unsigned int locus_position, const set<unsigned int> & groups) throw (ZeroDivisionException) {
	map<unsigned int, MultilocusGenotypeStatistics::VarComp> values;
	// Base values computation
	double nbar = 0.;
	double nc = 0.;
	vector<unsigned int> ids = getAllelesIdsForGroups(pmgc, locus_position, groups);
	map<unsigned int, double> pbar;
	map<unsigned int, double> s2;
	map<unsigned int, double> hbar;
	for (unsigned int i = 0 ; i < ids.size() ; i++) {
		pbar[ids[i]] = 0.;
		s2[ids[i]] = 0.;
		hbar[ids[i]] = 0.;
	}
	unsigned int r = groups.size();
	for (set<unsigned int>::iterator set_it = groups.begin() ; set_it != groups.end() ; set_it++)
    {
        int i  = (*set_it);
		double ni = (double) pmgc.getLocusGroupSize( i, locus_position);
		set<unsigned int> group_id;
		group_id.insert( i );
		map<unsigned int, double> pi = getAllelesFrqForGroups(pmgc, locus_position, group_id);
		map<unsigned int, double> hi = getHeterozygousFrqForGroups(pmgc, locus_position, group_id);
		nbar += ni;
		if (r > 1)
			nc += ni * ni;

		for (map<unsigned int, double>::iterator it = pi.begin() ; it != pi.end() ; it++)
			pbar[it->first] += ni * it->second;
		for (map<unsigned int, double>::iterator it = hi.begin() ; it != hi.end() ; it++)
			hbar[it->first] += ni * it->second;

		group_id.clear();	//khalid
	}
	nbar = nbar / (double) r;
	if (nbar <= 1)
		throw ZeroDivisionException("MultilocusGenotypeStatistics::getVarianceComponents.");
	if (r > 1)
		nc = (((double) r * nbar) - (nc / ((double) r * nbar))) / ((double) r - 1.);
	for (map<unsigned int, double>::iterator it = pbar.begin() ; it != pbar.end() ; it++)
		it->second = it->second / ((double) r * nbar);
	for (map<unsigned int, double>::iterator it = hbar.begin() ; it != hbar.end() ; it++)
		it->second = it->second / ((double) r * nbar);

	for (set<unsigned int>::iterator set_it = groups.begin() ; set_it != groups.end() ; set_it++)
    {
        int i  = (*set_it);
		double ni = (double) pmgc.getLocusGroupSize( i, locus_position);
		set<unsigned int> group_id;
		group_id.insert( i );
		map<unsigned int, double> pi = getAllelesFrqForGroups(pmgc, locus_position, group_id);
		for (unsigned int j = 0 ; j < ids.size() ; j++)
			pi[ids[j]];
		for (map<unsigned int, double>::iterator it = pi.begin() ; it != pi.end() ; it++)
			s2[it->first] += ni * (it->second - pbar[it->first]) * (it->second - pbar[it->first]);
		group_id.clear();	//khalid
	}
	for (map<unsigned int, double>::iterator it = s2.begin() ; it != s2.end() ; it++)
		it->second = it->second / (((double) r - 1.) * nbar);

	// a, b, c computation
	for (unsigned int i = 0 ; i < ids.size() ; i++) {
		values[ids[i]];
		}

	for (map<unsigned int, MultilocusGenotypeStatistics::VarComp>::iterator it = values.begin() ; it != values.end() ; it++)
    {
		it->second.a = (nbar / nc) * (s2[it->first] - ((1. / (nbar - 1.)) * ((pbar[it->first] * (1. - pbar[it->first])) - (s2[it->first] * ((double) r - 1.) / (double) r) - ((1. / 4.) * hbar[it->first]))));
		it->second.b = (nbar / (nbar - 1.)) * ((pbar[it->first] * (1. - pbar[it->first])) - (s2[it->first] * ((double) r - 1.) / (double) r) - ((((2. * nbar) - 1.) / (4. * nbar)) * hbar[it->first]));
		it->second.c = hbar[it->first] / 2.;
	}
	return values;
}

//khalid 9 juin 2006
double MultilocusGenotypeStatistics::getWCMultilocusFst(const PolymorphismMultiGContainer & pmgc, vector<unsigned int> locus_positions, const set<unsigned int> & groups) throw (Exception)
{
    double A, B, C;
    A=B=C = 0.0;
    for (int i=0; i < locus_positions.size(); i++)
    {
     //limiter les calculs pour les locus polymorphes pour ces groupes
     vector<unsigned int> ids = getAllelesIdsForGroups(pmgc, i, groups);
     if (ids.size() >= 2)
     {
      map<unsigned int, MultilocusGenotypeStatistics::VarComp> values = getVarianceComponents(pmgc, locus_positions[i], groups);
      for (map<unsigned int, MultilocusGenotypeStatistics::VarComp>::iterator it = values.begin() ; it != values.end() ; it++)
      {
		A += it->second.a ;
        B += it->second.b ;
        C += it->second.c ;
      }
     }
    }
    if ((A+B+C) == 0)  throw ZeroDivisionException("MultilocusGenotypeStatistics::getWCMultilocusFst.");
    return A/(A+B+C);
}
//khalid juin 2006
double MultilocusGenotypeStatistics::getWCMultilocusFis(const PolymorphismMultiGContainer & pmgc, vector<unsigned int> locus_positions, const set<unsigned int> & groups) throw (Exception)
{
    double A, B, C;
    A=B=C = 0.0;
    for (int i=0; i < locus_positions.size(); i++)
    {
     //limiter les calculs pour les locus polymorphes pour ces groupes
     vector<unsigned int> ids = getAllelesIdsForGroups(pmgc, i, groups);
     if (ids.size() >= 2)
     {
      map<unsigned int, MultilocusGenotypeStatistics::VarComp> values = getVarianceComponents(pmgc, locus_positions[i], groups);
      for (map<unsigned int, MultilocusGenotypeStatistics::VarComp>::iterator it = values.begin() ; it != values.end() ; it++)
      {
		A += it->second.a ;
        B += it->second.b ;
        C += it->second.c ;
      }
     }
    }
    if ((B+C) == 0)  throw ZeroDivisionException("MultilocusGenotypeStatistics::getWCMultilocusFis.");
    return 1.0 - C/(B+C);
}

//khalid 19/06/2006
MultilocusGenotypeStatistics::PermResults MultilocusGenotypeStatistics::getWCMultilocusFstAndPerm(const PolymorphismMultiGContainer & pmgc, vector<unsigned int> locus_positions,set<unsigned int> groups, int nb_perm) throw (Exception)
{
    //extraire un  PolymorphismMultiGContainer avec uniquement les groupes dans groups
    PolymorphismMultiGContainer sub_pmgc =  PolymorphismMultiGContainerTools::extractGroups(pmgc, groups);
    double nb_sup = 0.0;
    double nb_inf = 0.0;
    PermResults results;
    results.Statistic =  getWCMultilocusFst(sub_pmgc, locus_positions, groups);
    if (nb_perm > 0)
    {

     for (int  i = 0; i < nb_perm; i++)
     {
      PolymorphismMultiGContainer permuted_pmgc = PolymorphismMultiGContainerTools::permutMultiG( sub_pmgc);
      double Fst_perm =  getWCMultilocusFst(permuted_pmgc, locus_positions, groups);
      //cout << Fst_perm << endl;
      if (Fst_perm > results.Statistic) nb_sup++;
      if (Fst_perm < results.Statistic) nb_inf++;
     }

     nb_sup /= (double) nb_perm;
     nb_inf /= (double) nb_perm;
    }

    results.Percent_sup = nb_sup;
    results.Percent_inf = nb_inf;
    return results;
}

MultilocusGenotypeStatistics::PermResults MultilocusGenotypeStatistics::getWCMultilocusFisAndPerm(const PolymorphismMultiGContainer & pmgc, vector<unsigned int> locus_positions,set<unsigned int> groups, int nb_perm) throw (Exception)
{
    //extraire un  PolymorphismMultiGContainer avec uniquement les groupes dans groups
    PolymorphismMultiGContainer sub_pmgc =  PolymorphismMultiGContainerTools::extractGroups(pmgc, groups);
    double nb_sup = 0.0;
    double nb_inf = 0.0;
    PermResults results;
    results.Statistic =  getWCMultilocusFis(sub_pmgc, locus_positions, groups);
    if (nb_perm > 0)
    {

     for (int  i = 0; i < nb_perm; i++)
     {
      PolymorphismMultiGContainer permuted_pmgc = PolymorphismMultiGContainerTools::permutIntraGroupAlleles(sub_pmgc, groups);
      double Fis_perm =  getWCMultilocusFis(permuted_pmgc, locus_positions, groups);

      if (Fis_perm > results.Statistic) nb_sup++;
      if (Fis_perm < results.Statistic) nb_inf++;
     }

     nb_sup /= (double) nb_perm;
     nb_inf /= (double) nb_perm;
    }

    results.Percent_sup = nb_sup;
    results.Percent_inf = nb_inf;
    return results;
}

//khalid 14 juin 2006
double MultilocusGenotypeStatistics::getRHMultilocusFst(const PolymorphismMultiGContainer & pmgc, vector<unsigned int> locus_positions, const set<unsigned int> & groups) throw (Exception)
{
    double Au, Bu, Cu;
    double RH      = 0.0;
    int nb_alleles = 0;
    int total_alleles   = 0;

    for (int i=0; i < locus_positions.size(); i++)
    {
     //limiter les calculs pour les locus polymorphes pour ces groupes
     vector<unsigned int> ids = getAllelesIdsForGroups(pmgc, i, groups);
     if (ids.size() >= 2)
     {
      nb_alleles = 0;
      //les frequences allèliques moyennes
      map < unsigned int , double > P = MultilocusGenotypeStatistics::getAllelesFrqForGroups (pmgc ,locus_positions[i] , groups);
      //les composantes de la variance selon W&C
      map<unsigned int, MultilocusGenotypeStatistics::VarComp> values = getVarianceComponents(pmgc, locus_positions[i], groups);
      for (map<unsigned int, MultilocusGenotypeStatistics::VarComp>::iterator it = values.begin() ; it != values.end() ; it++)
      {

		Au = it->second.a ;
        Bu = it->second.b ;
        Cu = it->second.c ;
        if ((Au + Bu + Cu) != 0 )
        {
                double Pu = P[it->first]; //it->first contient le n° de l'allèle
                RH += (1- Pu) * Au / (Au + Bu + Cu);
                nb_alleles ++;
        }
      }
      total_alleles += (nb_alleles - 1) ;
     }
    }
    if (total_alleles == 0)  throw ZeroDivisionException("MultilocusGenotypeStatistics::getRHMultilocusFst.");
    return  RH/ double(total_alleles);
}

//khalid 13 juin 2006
DistanceMatrix * MultilocusGenotypeStatistics::getDistanceMatrix(const PolymorphismMultiGContainer & pmgc, vector<unsigned int> locus_positions, const set<unsigned int> & groups, string distance_methode) throw (Exception)
{//retourne une matrice de distance pour chaque paire de l'ensemble de groupes

    vector< string > names = pmgc.getAllGroupsNames();
    vector< unsigned int> grp_ids_vect;
    for (set<unsigned int>::iterator i=groups.begin(); i != groups.end(); i++)
    {
         //names.push_back(TextTools::toString(*i));
         grp_ids_vect.push_back(*i);
    }

    DistanceMatrix *_dist = new DistanceMatrix(names);
    for (int i=0; i < groups.size(); i++) (* _dist)(i, i) = 0;

    set<unsigned int> pairwise_grp;
    //for each pair in the groups set

    for ( unsigned int j =0; j < groups.size () - 1; j ++)
    {
     for ( unsigned int k=j +1; k < groups.size (); k ++)
     {
      double distance = 0;
	  if (distance_methode ==  "nei72" )
   	  distance = MultilocusGenotypeStatistics::getDnei72( pmgc,locus_positions, grp_ids_vect[j] , grp_ids_vect[k] );
	  else if  (distance_methode == "nei78" )
	  distance = MultilocusGenotypeStatistics::getDnei78( pmgc,locus_positions, grp_ids_vect[j] , grp_ids_vect[k] );
	  else if (distance_methode == "WC" ) // Fst multilocus selon W&C
	  {
       pairwise_grp.insert(grp_ids_vect[j] );
       pairwise_grp.insert(grp_ids_vect[k] );
	   distance = MultilocusGenotypeStatistics::getWCMultilocusFst( pmgc, locus_positions, pairwise_grp) ;
	   pairwise_grp.clear();
      }
      else if (distance_methode == "RH" ) // Fst multilocus selon ponderation Robertson & Hill
	  {
       pairwise_grp.insert(grp_ids_vect[j] );
       pairwise_grp.insert(grp_ids_vect[k] );
	   distance = MultilocusGenotypeStatistics::getRHMultilocusFst( pmgc, locus_positions, pairwise_grp) ;
	   pairwise_grp.clear();
      }
      else if (distance_methode == "Nm" ) // Nm déduit des Fst multilocus selon W&C modèle en îles Fst = 1/(1+4Nm)
	  {
       pairwise_grp.insert(grp_ids_vect[j] );
       pairwise_grp.insert(grp_ids_vect[k] );
	   distance = MultilocusGenotypeStatistics::getWCMultilocusFst( pmgc, locus_positions, pairwise_grp) ;
	   if (distance != 0)
       	   distance = 0.25 * (1- distance)/distance;
       else distance = NAN;
	   pairwise_grp.clear();
      }
      else if (distance_methode == "D" ) // D=-ln(1-Fst) of Reynolds, Weir and Cockerham, 1983
	  {
       pairwise_grp.insert(grp_ids_vect[j] );
       pairwise_grp.insert(grp_ids_vect[k] );
	   distance = MultilocusGenotypeStatistics::getWCMultilocusFst( pmgc, locus_positions, pairwise_grp) ;
	   if (distance != 1)
       	   distance =  -log(1- distance);
       else distance = NAN;
	   pairwise_grp.clear();
      }
      else if (distance_methode == "Rousset" ) // Calcul de Fst/(1-Fst). Rousset F. 1997
	  {
       pairwise_grp.insert(grp_ids_vect[j] );
       pairwise_grp.insert(grp_ids_vect[k] );
	   distance = MultilocusGenotypeStatistics::getWCMultilocusFst( pmgc, locus_positions, pairwise_grp) ;
	   if (distance != 1)
       	   distance = distance/(1- distance);
       else distance = NAN;
	   pairwise_grp.clear();
      }



      (* _dist)(k, j) =  distance;
      (* _dist)(j, k) =  distance;
     }//for k
    }//for j

    return _dist;
}
