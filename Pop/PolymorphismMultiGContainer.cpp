/*
 * File PolymorphismMultiGContainer.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Friday July 23 2004
 */

#include "PolymorphismMultiGContainer.h"

//** Constructors : **********************************************************/
PolymorphismMultiGContainer::PolymorphismMultiGContainer() {}

PolymorphismMultiGContainer::PolymorphismMultiGContainer(const PolymorphismMultiGContainer & pmgc) {
	for (unsigned int i = 0 ; i < pmgc.size() ; i++) {
		_multilocusGenotypes.push_back(new MultilocusGenotype(* pmgc.getMultilocusGenotype(i)));
		_groups.push_back(pmgc.getGroupId(i));
	}
}

//** Destructor : ************************************************************/
PolymorphismMultiGContainer::~PolymorphismMultiGContainer() {
	clear();
}

//** Other methodes : ********************************************************/
PolymorphismMultiGContainer & PolymorphismMultiGContainer::operator= (const PolymorphismMultiGContainer & pmgc) {
	for (unsigned int i = 0 ; i < pmgc.size() ; i++) {
		_multilocusGenotypes.push_back(new MultilocusGenotype(* pmgc.getMultilocusGenotype(i)));
		_groups.push_back(pmgc.getGroupId(i));
	}
	return * this;
}

Clonable * PolymorphismMultiGContainer::clone() const {
	return new PolymorphismMultiGContainer(* this);
}

void PolymorphismMultiGContainer::addMultilocusGenotype(const MultilocusGenotype & mg, unsigned int group) {
	_multilocusGenotypes.push_back(new MultilocusGenotype(mg));
	_groups.push_back(group);
}

const MultilocusGenotype * PolymorphismMultiGContainer::getMultilocusGenotype(unsigned int position) const throw (IndexOutOfBoundsException) {
	if (position >= size())
		throw IndexOutOfBoundsException("PolymorphismMultiGContainer::getMultilocusGenotype: position out of bounds.", position, 0, size() - 1);
	return _multilocusGenotypes[position];
}

MultilocusGenotype * PolymorphismMultiGContainer::removeMultilocusGenotype(unsigned int position) throw (IndexOutOfBoundsException) {
	if (position >= size())
		throw IndexOutOfBoundsException("PolymorphismMultiGContainer::removeMultilocusGenotype: position out of bounds.", position, 0, size() - 1);
	MultilocusGenotype * tmp_mg = _multilocusGenotypes[position];
	_multilocusGenotypes.erase(_multilocusGenotypes.begin() + position);
	_groups.erase(_groups.begin() + position);
	return tmp_mg;
}

void PolymorphismMultiGContainer::deleteMultilocusGenotype(unsigned int position) throw (IndexOutOfBoundsException) {
	if (position >= size())
		throw IndexOutOfBoundsException("PolymorphismMultiGContainer::deleteMultilocusGenotype: position out of bounds.", position, 0, size() - 1);
	delete _multilocusGenotypes[position];
	_multilocusGenotypes.erase(_multilocusGenotypes.begin() + position);
	_groups.erase(_groups.begin() + position);
}

unsigned int PolymorphismMultiGContainer::getGroupId(unsigned int position) const throw (IndexOutOfBoundsException) {
	if (position >= size())
		throw IndexOutOfBoundsException("PolymorphismMultiGContainer::getGroup: position out of bounds.", position, 0, size() - 1);
	return _groups[position];
}

set<unsigned int> PolymorphismMultiGContainer::getAllGroupsIds() const {
	set<unsigned int> groups_ids;
	for (unsigned int i = 0 ; i < size() ; i++)
		groups_ids.insert(_groups[i]);
	return groups_ids;
}

bool PolymorphismMultiGContainer::groupExists(unsigned int group) const {
	for (unsigned int i = 0 ; i < size() ; i++)
		if (_groups[i] == group)
			return true;
	return false;
}

unsigned int PolymorphismMultiGContainer::getNumberOfGroups() const {
	return getAllGroupsIds().size();
}

unsigned int PolymorphismMultiGContainer::getGroupSize(unsigned int group) const {
	unsigned int counter = 0;
	for (unsigned int i = 0 ; i < size() ; i++)
		if (_groups[i] == group)
			counter++;
	return counter;
}

unsigned int PolymorphismMultiGContainer::getLocusGroupSize(unsigned int group, unsigned int locus_position) const {
	unsigned int counter = 0;
	for (unsigned int i = 0 ; i < size() ; i++) {
		try {
			if (_groups[i] == group && ! _multilocusGenotypes[i]->isMonolocusGenotypeMissing(locus_position))
				counter++;
		}
		catch (IndexOutOfBoundsException & ioobe) {
			throw IndexOutOfBoundsException("PolymorphismMultiGContainer::getGroupSize: locus_position out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
		}
	}
	return counter;
}

unsigned int PolymorphismMultiGContainer::size() const {
	return _multilocusGenotypes.size();
}

void PolymorphismMultiGContainer::clear() {
	for (unsigned int i = 0 ; i < _multilocusGenotypes.size() ; i++)
		delete _multilocusGenotypes[i];
	_multilocusGenotypes.clear();
	_groups.clear();
}

vector<unsigned int> PolymorphismMultiGContainer::getAllelesIdsForAllGroups(unsigned int locus_position) const throw (IndexOutOfBoundsException) {
	set<unsigned int> groups_ids = getAllGroupsIds();
	try {
		return getAllelesIdsForGroups(locus_position, groups_ids);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("PolymorphismMultiGContainer::getAllelesIdsForAllGroups: locus_position out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
}

vector<unsigned int> PolymorphismMultiGContainer::getAllelesIdsForGroups(unsigned int locus_position, const set<unsigned int> & groups) const throw (IndexOutOfBoundsException) {
	map<unsigned int, unsigned int> tmp_alleles;
	try {
		tmp_alleles = getAllelesMapForGroups(locus_position, groups);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("PolymorphismMultiGContainer::getAllelesIdsForGroups: locus_position out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
	return MapTools::getKeys(tmp_alleles);
}

unsigned int PolymorphismMultiGContainer::countGametesForGroups(unsigned int locus_position, const set<unsigned int> & groups) const throw (IndexOutOfBoundsException) {
	map<unsigned int, unsigned int> allele_count;
	unsigned int nb_tot_allele = 0;
	try {
		allele_count = getAllelesMapForGroups(locus_position, groups);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("PolymorphismMultiGContainer::countGametesForGroups: locus_position out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
	vector<unsigned int> counter = MapTools::getValues(allele_count);
	for (unsigned int i = 0 ; i < counter.size() ; i++)
		nb_tot_allele += counter[i];
	return nb_tot_allele;
}

map<unsigned int, unsigned int> PolymorphismMultiGContainer::getAllelesMapForAllGroups(unsigned int locus_position) const throw (IndexOutOfBoundsException) {
	set<unsigned int> groups_ids = getAllGroupsIds();
	try {
		return getAllelesMapForGroups(locus_position, groups_ids);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("PolymorphismMultiGContainer::getAllelesMapForAllGroups: locus_position out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
}

map<unsigned int, unsigned int> PolymorphismMultiGContainer::getAllelesMapForOneGroup(unsigned int locus_position, unsigned int group) const throw (Exception) {
	if (! groupExists(group))
		throw GroupNotFoundException("PolymorphismMultiGContainer::getAllelesMapForOneGroup: group not found.", group);
	set<unsigned int> group_id;
	group_id.insert(group);
	try {
		return getAllelesMapForGroups(locus_position, group_id);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("PolymorphismMultiGContainer::getAllelesMapForOneGroup: locus_position out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
}

map<unsigned int, unsigned int> PolymorphismMultiGContainer::getAllelesMapForGroups(unsigned int locus_position, const set<unsigned int> & groups) const throw (IndexOutOfBoundsException) {
	map<unsigned int, unsigned int> alleles_count;
	for (unsigned int i = 0 ; i < size() ; i++) {
		try {
			if (! _multilocusGenotypes[i]->isMonolocusGenotypeMissing(locus_position) && find(groups.begin(), groups.end(), _groups[i]) != groups.end()) {
				vector<unsigned int> tmp_alleles = _multilocusGenotypes[i]->getMonolocusGenotype(locus_position)->getAlleleIndex();
				for (unsigned int j = 0 ; j < tmp_alleles.size() ; j++) alleles_count[tmp_alleles[j]]++;
			}
		}
		catch (IndexOutOfBoundsException & ioobe) {
			throw IndexOutOfBoundsException("PolymorphismMultiGContainer::getAllelesMapForGroups: locus_position out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
		}
	}
	return alleles_count;
}

map<unsigned int, double> PolymorphismMultiGContainer::getAllelesFrqForAllGroups(unsigned int locus_position) const throw (Exception) {
	set<unsigned int> groups_ids = getAllGroupsIds();
	try {
		return getAllelesFrqForGroups(locus_position, groups_ids);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("PolymorphismMultiGContainer::getAllelesFrqForAllGroups: locus_position out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
	catch (ZeroDivisionException & zde) {
		throw ZeroDivisionException("PolymorphismMultiGContainer::getAllelesFrqForAllGroups.");
	}
}

map<unsigned int, double> PolymorphismMultiGContainer::getAllelesFrqForOneGroup(unsigned int locus_position, unsigned int group) const throw (Exception) {
	if (! groupExists(group))
		throw GroupNotFoundException("PolymorphismMultiGContainer::getAllelesFrqForOneGroup: group not found.", group);
	set<unsigned int> group_id;
	group_id.insert(group);
	try {
		return getAllelesFrqForGroups(locus_position, group_id);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("PolymorphismMultiGContainer::getAllelesFrqForOneGroup: locus_position out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
	catch (ZeroDivisionException & zde) {
		throw ZeroDivisionException("PolymorphismMultiGContainer::getAllelesFrqForOneGroup.");
	}
}

map<unsigned int, double> PolymorphismMultiGContainer::getAllelesFrqForGroups(unsigned int locus_position, const set<unsigned int> & groups) const throw (Exception) {
	map<unsigned int, double> alleles_frq;
	unsigned int nb_tot_allele = 0;
	map<unsigned int, unsigned int> tmp_alleles;
	try {
		tmp_alleles = getAllelesMapForGroups(locus_position, groups);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("PolymorphismMultiGContainer::getAllelesFrqForGroups: locus_position out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
	vector<unsigned int> counter = MapTools::getValues(tmp_alleles);
	for (unsigned int i = 0 ; i < counter.size() ; i++)
		nb_tot_allele += counter[i];
	if (nb_tot_allele == 0)
		throw ZeroDivisionException("PolymorphismMultiGContainer::getAllelesFrqForGroups.");
	for (map<unsigned int, unsigned int>::iterator it = tmp_alleles.begin() ; it != tmp_alleles.end() ; it++)
		alleles_frq[it->first] = it->second / nb_tot_allele;
	return alleles_frq;
}

unsigned int PolymorphismMultiGContainer::countNonMissingForAllGroups(unsigned int locus_position) const throw (IndexOutOfBoundsException) {
	set<unsigned int> groups = getAllGroupsIds();
	try {
		return countNonMissingForGroups(locus_position, groups);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("PolymorphismMultiGContainer::countNonMissingForAllGroups: locus_position out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
}

unsigned int PolymorphismMultiGContainer::countNonMissingForOneGroup(unsigned int locus_position, unsigned int group) const throw (Exception) {
	if (! groupExists(group))
		throw GroupNotFoundException("PolymorphismMultiGContainer::countNonMissing: group not found.", group);
	set<unsigned int> group_id;
	group_id.insert(group);
	try {
		return countNonMissingForGroups(locus_position, group_id);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("PolymorphismMultiGContainer::countNonMissingForOneGroup: locus_position out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
}

unsigned int PolymorphismMultiGContainer::countNonMissingForGroups(unsigned int locus_position, const set<unsigned int> & groups) const throw (IndexOutOfBoundsException) {
	unsigned int counter = 0;
	for (unsigned int i = 0 ; i < size() ; i++) {
		try {
			if (! _multilocusGenotypes[i]->isMonolocusGenotypeMissing(locus_position) && find(groups.begin(), groups.end(), _groups[i]) != groups.end())
				counter++;
		}
		catch (IndexOutOfBoundsException & ioobe) {
			throw IndexOutOfBoundsException("PolymorphismMultiGContainer::countNonMissing: locus_position out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
		}
	}
	return counter;
}

unsigned int PolymorphismMultiGContainer::countBiAllelicForAllGroups(unsigned int locus_position) const throw (IndexOutOfBoundsException) {
	set<unsigned int> groups_ids = getAllGroupsIds();
	try {
		return countBiAllelicForGroups(locus_position, groups_ids);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("PolymorphismMultiGContainer::countBiAllelicForAllGroups: locus_position out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
}

unsigned int PolymorphismMultiGContainer::countBiAllelicForOneGroup(unsigned int locus_position, unsigned int group) const throw (Exception) {
	if (! groupExists(group))
		throw GroupNotFoundException("PolymorphismMultiGContainer::countBiAllelic: group not found.", group);
	set<unsigned int> group_id;
	group_id.insert(group);
	try {
		return countBiAllelicForGroups(locus_position, group_id);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("PolymorphismMultiGContainer::countBiAllelicForOneGroup: locus_position out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
}

unsigned int PolymorphismMultiGContainer::countBiAllelicForGroups(unsigned int locus_position, const set<unsigned int> & groups) const throw (IndexOutOfBoundsException) {
	unsigned int counter = 0;
	for (unsigned int i = 0 ; i < size() ; i++) {
		try {
			if (! _multilocusGenotypes[i]->isMonolocusGenotypeMissing(locus_position) && find(groups.begin(), groups.end(), _groups[i]) != groups.end())
				if ((_multilocusGenotypes[i]->getMonolocusGenotype(locus_position)->getAlleleIndex()).size() == 2)
					counter++;
		}
		catch (IndexOutOfBoundsException & ioobe) {
			throw IndexOutOfBoundsException("PolymorphismMultiGContainer::countBiAllelic: locus_position out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
		}
	}
	return counter;
}

map<unsigned int, unsigned int> PolymorphismMultiGContainer::countHeterozygousForAllGroups(unsigned int locus_position) const throw (IndexOutOfBoundsException) {
	set<unsigned int> groups_ids = getAllGroupsIds();
	try {
		return countHeterozygousForGroups(locus_position, groups_ids);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("PolymorphismMultiGContainer::countHeterozygousForAllGroups: locus_position out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
}

map<unsigned int, unsigned int> PolymorphismMultiGContainer::countHeterozygousForOneGroup(unsigned int locus_position, unsigned int group) const throw (Exception) {
	if (! groupExists(group))
		throw GroupNotFoundException("PolymorphismMultiGContainer::countHeterozygousForOneGroup: group not found.", group);
	set<unsigned int> group_id;
	group_id.insert(group);
	try {
		return countHeterozygousForGroups(locus_position, group_id);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("PolymorphismMultiGContainer::countHeterozygousForOneGroup: locus_position out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
}

map<unsigned int, unsigned int> PolymorphismMultiGContainer::countHeterozygousForGroups(unsigned int locus_position, const set<unsigned int> & groups) const throw (IndexOutOfBoundsException) {
	map<unsigned int, unsigned int> counter;
	for (unsigned int i = 0 ; i < size() ; i++) {
		try {
			if (! _multilocusGenotypes[i]->isMonolocusGenotypeMissing(locus_position) && find(groups.begin(), groups.end(), _groups[i]) != groups.end()) {
				const MonolocusGenotype * tmp_mg = _multilocusGenotypes[i]->getMonolocusGenotype(locus_position);
				if ((tmp_mg->getAlleleIndex()).size() == 2) {
					if (! dynamic_cast<const BiAlleleMonolocusGenotype *>(tmp_mg)->isHomozygous()) {
						vector<unsigned int> tmp_alleles = tmp_mg->getAlleleIndex();
						for (unsigned int j = 0 ; j < tmp_alleles.size() ; j++) counter[tmp_alleles[j]]++;
					}
				}
			}
		}
		catch (IndexOutOfBoundsException & ioobe) {
			throw IndexOutOfBoundsException("PolymorphismMultiGContainer::countHeterozygous: locus_position out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
		}
	}
	return counter;
}

map<unsigned int, double> PolymorphismMultiGContainer::getHeterozygousFrqForAllGroups(unsigned int locus_position) const throw (Exception) {
	set<unsigned int> groups_ids = getAllGroupsIds();
	try {
		return getHeterozygousFrqForGroups(locus_position, groups_ids);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("PolymorphismMultiGContainer::getHeterozygousFrqForAllGroups: locus_position out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
	catch (ZeroDivisionException & zde) {
		throw ZeroDivisionException("PolymorphismMultiGContainer::getHeterozygousFrqForAllGroups.");
	}
}

map<unsigned int, double> PolymorphismMultiGContainer::getHeterozygousFrqForOneGroup(unsigned int locus_position, unsigned int group) const throw (Exception) {
	if (! groupExists(group))
		throw GroupNotFoundException("PolymorphismMultiGContainer::getHeterozygousFrqForOneGroup: group not found.", group);
	set<unsigned int> group_id;
	group_id.insert(group);
	try {
		return getHeterozygousFrqForGroups(locus_position, group_id);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("PolymorphismMultiGContainer::getHeterozygousFrqForOneGroup: locus_position out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
	catch (ZeroDivisionException & zde) {
		throw ZeroDivisionException("PolymorphismMultiGContainer::getHeterozygousFrqForOneGroup.");
	}
}

map<unsigned int, double> PolymorphismMultiGContainer::getHeterozygousFrqForGroups(unsigned int locus_position, const set<unsigned int> & groups) const throw (Exception) {
	map<unsigned int, double> freq;
	unsigned int counter = 0;
	for (unsigned int i = 0 ; i < size() ; i++) {
		try {
			if (! _multilocusGenotypes[i]->isMonolocusGenotypeMissing(locus_position) && find(groups.begin(), groups.end(), _groups[i]) != groups.end()) {
				const MonolocusGenotype * tmp_mg = _multilocusGenotypes[i]->getMonolocusGenotype(locus_position);
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
			throw IndexOutOfBoundsException("PolymorphismMultiGContainer::getHeterozygousFrqForGroups: locus_position out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
		}
	}
	if (counter == 0)
		throw ZeroDivisionException("PolymorphismMultiGContainer::getHeterozygousFrqForGroups.");
	for (map<unsigned int, double>::iterator i = freq.begin() ; i != freq.end() ; i++)
		i->second = i->second / counter;
	return freq;
}

double PolymorphismMultiGContainer::getHobsForGroups(unsigned int locus_position, const set<unsigned int> & groups) const throw (Exception) {
	map<unsigned int, double> heterozygous_frq;
	double frq = 0.;
	try {
		heterozygous_frq = getHeterozygousFrqForGroups(locus_position, groups);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("PolymorphismMultiGContainer::getHobsForGroups: locus_position out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
	catch (ZeroDivisionException & zde) {
		throw ZeroDivisionException("PolymorphismMultiGContainer::getHobsForGroups.");
	}
	for (map<unsigned int, double>::iterator it = heterozygous_frq.begin() ; it != heterozygous_frq.end() ; it++)
		frq += it->second;
	return frq / heterozygous_frq.size();
}

double PolymorphismMultiGContainer::getHexpForGroups(unsigned int locus_position, const set<unsigned int> & groups) const throw (Exception) {
	map<unsigned int, double> allele_frq;
	double frqsqr = 0.;
	try {
		allele_frq = getAllelesFrqForGroups(locus_position, groups);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("PolymorphismMultiGContainer::getHexpForGroups: locus_position out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
	catch (ZeroDivisionException & zde) {
		throw ZeroDivisionException("PolymorphismMultiGContainer::getHexpForGroups.");
	}
	for (map<unsigned int, double>::iterator it = allele_frq.begin() ; it != allele_frq.end() ; it++)
		frqsqr += it->second * it->second;
	return (1 - frqsqr);
}

double PolymorphismMultiGContainer::getHnbForGroups(unsigned int locus_position, const set<unsigned int> & groups) const throw (Exception) {
	unsigned int nb_alleles;
	double Hexp;
	try {
		nb_alleles = countGametesForGroups(locus_position, groups);
		Hexp = getHexpForGroups(locus_position, groups);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("PolymorphismMultiGContainer::getHnbForGroups: locus_position out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
	catch (ZeroDivisionException & zde) {
		throw ZeroDivisionException("PolymorphismMultiGContainer::getHnbForGroups.");
	}
	return (2 * nb_alleles * Hexp  / ((2 * nb_alleles) - 1));
}

double PolymorphismMultiGContainer::getDnei72(unsigned int locus_position, unsigned int grp1, unsigned int grp2) const throw (Exception) {
	map<unsigned int, double> allele_frq1, allele_frq2;
	vector<unsigned int> allele_ids;
	set<unsigned int> groups_id;
	double Jx = 0.;
	double Jy = 0.;
	double Jxy = 0.;
	groups_id.insert(grp1);
	groups_id.insert(grp2);
	try {
		allele_ids = getAllelesIdsForGroups(locus_position, groups_id);
		allele_frq1 = getAllelesFrqForOneGroup(locus_position, grp1);
		allele_frq2 = getAllelesFrqForOneGroup(locus_position, grp2);
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

double PolymorphismMultiGContainer::getDnei78(unsigned int locus_position, unsigned int grp1, unsigned int grp2) const throw (Exception) {
	map<unsigned int, double> allele_frq1, allele_frq2;
	vector<unsigned int> allele_ids;
	set<unsigned int> groups_id;
	double Jx = 0.;
	double Jy = 0.;
	double Jxy = 0.;
	unsigned int nx, ny;
	groups_id.insert(grp1);
	groups_id.insert(grp2);
	try {
		allele_ids = getAllelesIdsForGroups(locus_position, groups_id);
		allele_frq1 = getAllelesFrqForOneGroup(locus_position, grp1);
		allele_frq2 = getAllelesFrqForOneGroup(locus_position, grp2);
		nx = countBiAllelicForOneGroup(locus_position, grp1);
		ny = countBiAllelicForOneGroup(locus_position, grp2);
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
