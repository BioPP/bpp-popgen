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

unsigned int PolymorphismMultiGContainer::size() const {
	return _multilocusGenotypes.size();
}

void PolymorphismMultiGContainer::clear() {
	for (unsigned int i = 0 ; i < _multilocusGenotypes.size() ; i++)
		delete _multilocusGenotypes[i];
	_multilocusGenotypes.clear();
	_groups.clear();
}
