/*
 * File MultilocusGenotypeStatistics.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Tuesday August 03 2004
 *
 * Copyright (C) 2004 Sylvain Gaillard and the
 *                    PopLib Development Core Team
 *
 * This file is part of PopLib.
 *
 * PopLib is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * PopLib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PopLib; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include "MultilocusGenotypeStatistics.h"

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
			if (! pmgc.getMultilocusGenotype(i)->isMonolocusGenotypeMissing(locus_position) && find(groups.begin(), groups.end(), pmgc.getGroupId(i)) != groups.end()) {
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
			if (! pmgc.getMultilocusGenotype(i)->isMonolocusGenotypeMissing(locus_position) && find(groups.begin(), groups.end(), pmgc.getGroupId(i)) != groups.end())
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
			if (! pmgc.getMultilocusGenotype(i)->isMonolocusGenotypeMissing(locus_position) && find(groups.begin(), groups.end(), pmgc.getGroupId(i)) != groups.end())
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
			if (! pmgc.getMultilocusGenotype(i)->isMonolocusGenotypeMissing(locus_position) && find(groups.begin(), groups.end(), pmgc.getGroupId(i)) != groups.end()) {
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
			if (! pmgc.getMultilocusGenotype(i)->isMonolocusGenotypeMissing(locus_position) && find(groups.begin(), groups.end(), pmgc.getGroupId(i)) != groups.end()) {
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

double MultilocusGenotypeStatistics::getDnei72(const PolymorphismMultiGContainer & pmgc, unsigned int locus_position, unsigned int grp1, unsigned int grp2) throw (Exception) {
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
	try {
		allele_ids = getAllelesIdsForGroups(pmgc, locus_position, groups_id);
		allele_frq1 = getAllelesFrqForGroups(pmgc, locus_position, group1_id);
		allele_frq2 = getAllelesFrqForGroups(pmgc, locus_position, group2_id);
	}
	catch (Exception & e) {
		throw e;
	}
	for (unsigned int i = 0 ; i < allele_ids.size() ; i++) {
		map<unsigned int, double>::iterator it1 = allele_frq1.find(i);
		map<unsigned int, double>::iterator it2 = allele_frq2.find(i);
		double tmp_frq1 = (it1 != allele_frq1.end()) ? it1->second : 0;
		double tmp_frq2 = (it2 != allele_frq2.end()) ? it2->second : 0;
		Jx += tmp_frq1 * tmp_frq1;
		Jy += tmp_frq2 * tmp_frq2;
		Jxy += tmp_frq1 * tmp_frq2;
	}
	return -log(Jxy / sqrt(Jx * Jy));
}

double MultilocusGenotypeStatistics::getDnei78(const PolymorphismMultiGContainer & pmgc, unsigned int locus_position, unsigned int grp1, unsigned int grp2) throw (Exception) {
	map<unsigned int, double> allele_frq1, allele_frq2;
	vector<unsigned int> allele_ids;
	set<unsigned int> group1_id;
	set<unsigned int> group2_id;
	set<unsigned int> groups_id;
	double Jx = 0.;
	double Jy = 0.;
	double Jxy = 0.;
	unsigned int nx, ny;
	group1_id.insert(grp1);
	group2_id.insert(grp2);
	groups_id.insert(grp1);
	groups_id.insert(grp2);
	try {
		allele_ids = getAllelesIdsForGroups(pmgc, locus_position, groups_id);
		allele_frq1 = getAllelesFrqForGroups(pmgc, locus_position, group1_id);
		allele_frq2 = getAllelesFrqForGroups(pmgc, locus_position, group2_id);
		nx = countBiAllelicForGroups(pmgc, locus_position, group1_id);
		ny = countBiAllelicForGroups(pmgc, locus_position, group2_id);
	}
	catch (Exception & e) {
		throw e;
	}
	for (unsigned int i = 0 ; i < allele_ids.size() ; i++) {
		map<unsigned int, double>::iterator it1 = allele_frq1.find(i);
		map<unsigned int, double>::iterator it2 = allele_frq2.find(i);
		double tmp_frq1 = (it1 != allele_frq1.end()) ? it1->second : 0;
		double tmp_frq2 = (it2 != allele_frq2.end()) ? it2->second : 0;
		Jx += tmp_frq1 * tmp_frq1;
		Jy += tmp_frq2 * tmp_frq2;
		Jxy += tmp_frq1 * tmp_frq2;
	}
	return -log(Jxy / sqrt((((2 * nx * Jx) - 1)/((2 * nx) -1)) * (((2 * ny * Jy) - 1)/((2 * ny) -1))));
}

map<unsigned int, double> MultilocusGenotypeStatistics::getWCFit(const PolymorphismMultiGContainer & pmgc, unsigned int locus_position, const set<unsigned int> & groups) throw (Exception) {
	map<unsigned int, MultilocusGenotypeStatistics::FstatBases> values = getVarianceComponents(pmgc, locus_position, groups);
	map<unsigned int, double> Fit;
	for (map<unsigned int, MultilocusGenotypeStatistics::FstatBases>::iterator it = values.begin() ; it != values.end() ; it++) {
		Fit[it->first] = it->second.a + it->second.b + it->second.c;
		if (Fit[it->first] == 0.)
			throw ZeroDivisionException("MultilocusGenotypeStatistics::getWCFit.");
		Fit[it->first] = 1. - it->second.c / Fit[it->first];
	}
	return Fit;
}

map<unsigned int, double> MultilocusGenotypeStatistics::getWCFst(const PolymorphismMultiGContainer & pmgc, unsigned int locus_position, const set<unsigned int> & groups) throw (Exception) {
	if (groups.size() <= 1)
		throw BadIntegerException("MultilocusGenotypeStatistics::getWCFst: groups must be >= 2.", groups.size());
	map<unsigned int, MultilocusGenotypeStatistics::FstatBases> values = getVarianceComponents(pmgc, locus_position, groups);
	map<unsigned int, double> Fst;
	for (map<unsigned int, MultilocusGenotypeStatistics::FstatBases>::iterator it = values.begin() ; it != values.end() ; it++) {
		Fst[it->first] = it->second.a + it->second.b + it->second.c;
		if (Fst[it->first] == 0.)
			throw ZeroDivisionException("MultilocusGenotypeStatistics::getWCFst.");
		Fst[it->first] = it->second.a / Fst[it->first];
	}
	return Fst;
}

map<unsigned int, double> MultilocusGenotypeStatistics::getWCFis(const PolymorphismMultiGContainer & pmgc, unsigned int locus_position, const set<unsigned int> & groups) throw (Exception) {
	map<unsigned int, MultilocusGenotypeStatistics::FstatBases> values = getVarianceComponents(pmgc, locus_position, groups);
	map<unsigned int, double> Fis;
	for (map<unsigned int, MultilocusGenotypeStatistics::FstatBases>::iterator it = values.begin() ; it != values.end() ; it++) {
		Fis[it->first] = it->second.b + it->second.c;
		if (Fis[it->first] == 0.)
			throw ZeroDivisionException("MultilocusGenotypeStatistics::getWCFis.");
		Fis[it->first] = 1. - it->second.c / Fis[it->first];
	}
	return Fis;
}

map<unsigned int, MultilocusGenotypeStatistics::FstatBases> MultilocusGenotypeStatistics::getVarianceComponents(const PolymorphismMultiGContainer & pmgc, unsigned int locus_position, const set<unsigned int> & groups) throw (ZeroDivisionException) {
	map<unsigned int, MultilocusGenotypeStatistics::FstatBases> values;
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
	for (unsigned int i = 0 ; i < r ; i++) {
		double ni = (double) pmgc.getLocusGroupSize(i, locus_position);
		set<unsigned int> group_id;
		group_id.insert(* groups.begin() + i);
		map<unsigned int, double> pi = getAllelesFrqForGroups(pmgc, locus_position, group_id);
		map<unsigned int, double> hi = getHeterozygousFrqForGroups(pmgc, locus_position, group_id);
		nbar += ni;
		if (r > 1)
			nc += ni * ni;
		for (map<unsigned int, double>::iterator it = pi.begin() ; it != pi.end() ; it++)
			pbar[it->first] += ni * it->second;
		for (map<unsigned int, double>::iterator it = hi.begin() ; it != hi.end() ; it++)
			hbar[it->first] += ni * it->second;
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
	for (unsigned int i = 0 ; i < r ; i++) {
		double ni = (double) pmgc.getLocusGroupSize(i, locus_position);
		set<unsigned int> group_id;
		group_id.insert(* groups.begin() + i);
		map<unsigned int, double> pi = getAllelesFrqForGroups(pmgc, locus_position, group_id);
		for (unsigned int i = 0 ; i < ids.size() ; i++)
			pi[ids[i]];
		for (map<unsigned int, double>::iterator it = pi.begin() ; it != pi.end() ; it++)
			s2[it->first] += ni * (it->second - pbar[it->first]) * (it->second - pbar[it->first]);
	}
	for (map<unsigned int, double>::iterator it = s2.begin() ; it != s2.end() ; it++)
		it->second = it->second / (((double) r - 1.) * nbar);
	
	// a, b, c computation
	for (unsigned int i = 0 ; i < ids.size() ; i++) {
		values[ids[i]];
	}
	for (map<unsigned int, MultilocusGenotypeStatistics::FstatBases>::iterator it = values.begin() ; it != values.end() ; it++) {
		it->second.a = (nbar / nc) * (s2[it->first] - ((1. / (nbar - 1.)) * ((pbar[it->first] * (1. - pbar[it->first])) - (s2[it->first] * ((double) r - 1.) / (double) r) - ((1. / 4.) * hbar[it->first]))));
		it->second.b = (nbar / (nbar - 1.)) * ((pbar[it->first] * (1. - pbar[it->first])) - (s2[it->first] * ((double) r - 1.) / (double) r) - ((((2. * nbar) - 1.) / (4. * nbar)) * hbar[it->first]));
		it->second.c = hbar[it->first] / 2.;
	}
	return values;
}
