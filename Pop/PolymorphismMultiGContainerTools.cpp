/*
 * File PolymorphismMultiGContainerTools.cpp
 * Author : Sylvain Gailard <yragael2001@yahoo.fr>
 * Last modification : Thursday September 30 2004
 *
 * Copyright (C) 2004 Sylvain Gaillard and the
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

#include "PolymorphismMultiGContainerTools.h"

PolymorphismMultiGContainerTools::~PolymorphismMultiGContainerTools() {}

PolymorphismMultiGContainer PolymorphismMultiGContainerTools::permutMultiG(const PolymorphismMultiGContainer & pmgc) {
	PolymorphismMultiGContainer permuted_pmgc(pmgc);
	vector<unsigned int> groups;
	for (unsigned int i = 0 ; i < permuted_pmgc.size() ; i++)
		groups.push_back(permuted_pmgc.getGroupId(i));
	groups = RandomTools::getSample(groups, groups.size());
	for (unsigned int i = 0 ; i < permuted_pmgc.size() ; i++)
		permuted_pmgc.setGroupId(i, groups[i]);
	return permuted_pmgc;
}

PolymorphismMultiGContainer PolymorphismMultiGContainerTools::permutMonoG(const PolymorphismMultiGContainer & pmgc, const set<unsigned int> & groups) {
	PolymorphismMultiGContainer permuted_pmgc;
	unsigned int loc_num = pmgc.getNumberOfLoci();
	vector<vector<const MonolocusGenotype *> > mono_gens;
	mono_gens.resize(loc_num);
	// Get all the MonolocusGenotypes to permut
	for (unsigned int i = 0 ; i < pmgc.size() ; i++) {
		if (groups.find(pmgc.getGroupId(i)) != groups.end()) {
			for (unsigned int j = 0 ; j < loc_num ; j++)
				mono_gens[j].push_back(pmgc.getMultilocusGenotype(i)->getMonolocusGenotype(j));
		}
	}
	// Permut the MonolocusGenotypes
	for (unsigned int i = 0 ; i < loc_num ; i++)
		mono_gens[i] = RandomTools::getSample(mono_gens[i], mono_gens[i].size());
	// Build the new PolymorphismMultiGContainer
	unsigned int k = 0;
	for (unsigned int i = 0 ; i < pmgc.size() ; i++) {
		if (groups.find(pmgc.getGroupId(i)) != groups.end()) {
			MultilocusGenotype tmp_mg(loc_num);
			for (unsigned int j = 0 ; j < loc_num ; j++) {
				if (mono_gens[j][k] != NULL)
					tmp_mg.setMonolocusGenotype(j, * (mono_gens[j][k]));
			}
			permuted_pmgc.addMultilocusGenotype(tmp_mg, pmgc.getGroupId(i));
			k++;
		}
		else {
			permuted_pmgc.addMultilocusGenotype(* (pmgc.getMultilocusGenotype(i)), pmgc.getGroupId(i));
		}
	}
	return permuted_pmgc;
}

PolymorphismMultiGContainer PolymorphismMultiGContainerTools::permutAlleles(const PolymorphismMultiGContainer & pmgc, const set<unsigned int> & groups) {
	PolymorphismMultiGContainer permuted_pmgc;
	unsigned int loc_num = pmgc.getNumberOfLoci();
	vector<vector<unsigned int> > alleles;
	alleles.resize(loc_num);
	// Get all the alleles to permut
	for (unsigned int i = 0 ; i < pmgc.size() ; i++) {
		if (groups.find(pmgc.getGroupId(i)) != groups.end()) {
			for (unsigned int j = 0 ; j < loc_num ; j++)
				if (pmgc.getMultilocusGenotype(i)->getMonolocusGenotype(j) != NULL)
					for (unsigned int k = 0 ; k < pmgc.getMultilocusGenotype(i)->getMonolocusGenotype(j)->getAlleleIndex().size() ; k++)
						alleles[j].push_back(pmgc.getMultilocusGenotype(i)->getMonolocusGenotype(j)->getAlleleIndex()[k]);
		}
	}
	// Permut the alleles
	for (unsigned int i = 0 ; i < loc_num ; i++)
		alleles[i] = RandomTools::getSample(alleles[i], alleles[i].size());
	// Build the new PolymorphismMultiGContainer
	vector<unsigned int> k(loc_num,0);
	for (unsigned int i = 0 ; i < pmgc.size() ; i++) {
		if (groups.find(pmgc.getGroupId(i)) != groups.end()) {
			MultilocusGenotype tmp_mg(loc_num);
			for (unsigned int j = 0 ; j < loc_num ; j++) {
				if (pmgc.getMultilocusGenotype(i)->getMonolocusGenotype(j) != NULL) {
					if (pmgc.getMultilocusGenotype(i)->getMonolocusGenotype(j)->getAlleleIndex().size() == 1)
						tmp_mg.setMonolocusGenotype(j, MonoAlleleMonolocusGenotype(alleles[j][k[j]++]));
					if (pmgc.getMultilocusGenotype(i)->getMonolocusGenotype(j)->getAlleleIndex().size() == 2)
						tmp_mg.setMonolocusGenotype(j, BiAlleleMonolocusGenotype(alleles[j][k[j]++], alleles[j][k[j]++]));
				}
			}
			permuted_pmgc.addMultilocusGenotype(tmp_mg, pmgc.getGroupId(i));
		}
		else {
			permuted_pmgc.addMultilocusGenotype(* (pmgc.getMultilocusGenotype(i)), pmgc.getGroupId(i));
		}
	}
	return permuted_pmgc;
}
