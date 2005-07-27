/*
 * File PolymorphismMultiGContainerTools.cpp
 * Author : Sylvain Gailard <yragael2001@yahoo.fr>
 * Last modification : Thursday September 30 2004
 *
*/
/*
Copyright or © or Copr. CNRS, (November 17, 2004)


This software is a computer program whose purpose is to provide classes
for sequences analysis.

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
