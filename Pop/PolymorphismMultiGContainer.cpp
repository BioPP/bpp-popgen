/*
 * File PolymorphismMultiGContainer.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Tuesday July 20 2004
 */

#include "PolymorphismMultiGContainer.h"

//** Constructors : **********************************************************/
PolymorphismMultiGContainer::PolymorphismMultiGContainer() {}

PolymorphismMultiGContainer::PolymorphismMultiGContainer(const PolymorphismMultiGContainer & pmgc) {
	for (unsigned int i = 0 ; i < pmgc.size() ; i++) {
		_multilocusGenotypes.push_back(new MultilocusGenotype(* getMultilocusGenotype(i)));
		_groups.push_back(pmgc.getGroup(i));
	}
}

//** Destructor : ************************************************************/
PolymorphismMultiGContainer::~PolymorphismMultiGContainer() {
	clear();
}

//** Other methodes : ********************************************************/
PolymorphismMultiGContainer & PolymorphismMultiGContainer::operator= (const PolymorphismMultiGContainer & pmgc) {
	for (unsigned int i = 0 ; i < pmgc.size() ; i++) {
		_multilocusGenotypes.push_back(new MultilocusGenotype(* getMultilocusGenotype(i)));
		_groups.push_back(pmgc.getGroup(i));
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

unsigned int PolymorphismMultiGContainer::getGroup(unsigned int position) const throw (IndexOutOfBoundsException) {
	if (position >= size())
		throw IndexOutOfBoundsException("PolymorphismMultiGContainer::getGroup: position out of bounds.", position, 0, size() - 1);
	return _groups[position];
}

bool PolymorphismMultiGContainer::groupExists(unsigned int group) const {
	for (unsigned int i = 0 ; i < size() ; i++)
		if (_groups[i] == group)
			return true;
	return false;
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

map<unsigned int, unsigned int> PolymorphismMultiGContainer::countAlleles(unsigned int locus_position) const throw (IndexOutOfBoundsException) {
	map<unsigned int, unsigned int> counter;
	for (unsigned int i = 0 ; i < size() ; i++) {
		try {
			if (! _multilocusGenotypes[i]->isMonolocusGenotypeMissing(locus_position)) {
				vector<unsigned int> tmp_alleles = _multilocusGenotypes[i]->getMonolocusGenotype(locus_position)->getAlleleIndex();
				for (unsigned int j = 0 ; j < tmp_alleles.size() ; j++) counter[tmp_alleles[j]]++;
			}
		}
		catch (IndexOutOfBoundsException & ioobe) {
			throw IndexOutOfBoundsException("PolymorphismMultiGContainer::countAlleles: locus_position out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
		}
	}
	return counter;
}

map<unsigned int, unsigned int> PolymorphismMultiGContainer::countAlleles(unsigned int locus_position, unsigned int group) const throw (Exception) {
	if (! groupExists(group))
		throw GroupNotFoundException("PolymorphismMultiGContainer::countAlleles: group not found.", group);
	map<unsigned int, unsigned int> counter;
	for (unsigned int i = 0 ; i < size() ; i++) {
		try {
			if (! _multilocusGenotypes[i]->isMonolocusGenotypeMissing(locus_position) && _groups[i] == group) {
				vector<unsigned int> tmp_alleles = _multilocusGenotypes[i]->getMonolocusGenotype(locus_position)->getAlleleIndex();
				for (unsigned int j = 0 ; j < tmp_alleles.size() ; j++) counter[tmp_alleles[j]]++;
			}
		}
		catch (IndexOutOfBoundsException & ioobe) {
			throw IndexOutOfBoundsException("PolymorphismMultiGContainer::countAlleles: locus_position out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
		}
	}
	return counter;
}

unsigned int PolymorphismMultiGContainer::countNonMissing(unsigned int locus_position) const throw (IndexOutOfBoundsException) {
	unsigned int counter = 0;
	for (unsigned int i = 0 ; i < size() ; i++) {
		try {
			if (! _multilocusGenotypes[i]->isMonolocusGenotypeMissing(locus_position))
				counter++;
		}
		catch (IndexOutOfBoundsException & ioobe) {
			throw IndexOutOfBoundsException("PolymorphismMultiGContainer::countNonMissing: locus_position out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
		}
	}
	return counter;
}

unsigned int PolymorphismMultiGContainer::countNonMissing(unsigned int locus_position, unsigned int group) const throw (Exception) {
	if (! groupExists(group))
		throw GroupNotFoundException("PolymorphismMultiGContainer::countNonMissing: group not found.", group);
	unsigned int counter = 0;
	for (unsigned int i = 0 ; i < size() ; i++) {
		try {
			if (! _multilocusGenotypes[i]->isMonolocusGenotypeMissing(locus_position) && _groups[i] == group)
				counter++;
		}
		catch (IndexOutOfBoundsException & ioobe) {
			throw IndexOutOfBoundsException("PolymorphismMultiGContainer::countNonMissing: locus_position out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
		}
	}
	return counter;
}

unsigned int PolymorphismMultiGContainer::countBiAllelic(unsigned int locus_position) const throw (IndexOutOfBoundsException) {
	unsigned int counter = 0;
	for (unsigned int i = 0 ; i < size() ; i++) {
		try {
			if (! _multilocusGenotypes[i]->isMonolocusGenotypeMissing(locus_position))
				if ((_multilocusGenotypes[i]->getMonolocusGenotype(locus_position)->getAlleleIndex()).size() == 2)
					counter++;
		}
		catch (IndexOutOfBoundsException & ioobe) {
			throw IndexOutOfBoundsException("PolymorphismMultiGContainer::countBiAllelic: locus_position out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
		}
	}
	return counter;
}

unsigned int PolymorphismMultiGContainer::countBiAllelic(unsigned int locus_position, unsigned int group) const throw (Exception) {
	if (! groupExists(group))
		throw GroupNotFoundException("PolymorphismMultiGContainer::countBiAllelic: group not found.", group);
	unsigned int counter = 0;
	for (unsigned int i = 0 ; i < size() ; i++) {
		try {
			if (! _multilocusGenotypes[i]->isMonolocusGenotypeMissing(locus_position) && _groups[i] == group)
				if ((_multilocusGenotypes[i]->getMonolocusGenotype(locus_position)->getAlleleIndex()).size() == 2)
					counter++;
		}
		catch (IndexOutOfBoundsException & ioobe) {
			throw IndexOutOfBoundsException("PolymorphismMultiGContainer::countBiAllelic: locus_position out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
		}
	}
	return counter;
}

map<unsigned int, unsigned int> PolymorphismMultiGContainer::countHeterozygous(unsigned int locus_position) const throw (IndexOutOfBoundsException) {
	map<unsigned int, unsigned int> counter;
	for (unsigned int i = 0 ; i < size() ; i++) {
		try {
			if (! _multilocusGenotypes[i]->isMonolocusGenotypeMissing(locus_position)) {
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

map<unsigned int, unsigned int> PolymorphismMultiGContainer::countHeterozygous(unsigned int locus_position, unsigned int group) const throw (Exception) {
	if (! groupExists(group))
		throw GroupNotFoundException("PolymorphismMultiGContainer::countHeterozygous: group not found.", group);
	map<unsigned int, unsigned int> counter;
	for (unsigned int i = 0 ; i < size() ; i++) {
		try {
			if (! _multilocusGenotypes[i]->isMonolocusGenotypeMissing(locus_position) && _groups[i] == group) {
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
