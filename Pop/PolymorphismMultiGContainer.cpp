/*
 * File PolymorphismMultiGContainer.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Tuesday September 28 2004
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

bool PolymorphismMultiGContainer::isAligned() const {
	unsigned int value = 0;
	for (unsigned int i = 0 ; i < size() ; i++) {
		if (i == 0)
			value = _multilocusGenotypes[i]->size();
		else
			if (_multilocusGenotypes[i]->size() != value)
				return false;
	}
	return true;
}

unsigned int PolymorphismMultiGContainer::getNumberOfLoci() const throw (Exception) {
	if (!isAligned())
		throw Exception("MultilocusGenotypes are not aligned.");
	if (size() < 1)
		return 0;
	return _multilocusGenotypes[0]->size();
}

unsigned int PolymorphismMultiGContainer::getGroupId(unsigned int position) const throw (IndexOutOfBoundsException) {
	if (position >= size())
		throw IndexOutOfBoundsException("PolymorphismMultiGContainer::getGroupId: position out of bounds.", position, 0, size() - 1);
	return _groups[position];
}

void PolymorphismMultiGContainer::setGroupId(unsigned int position, unsigned int group_id) throw (IndexOutOfBoundsException) {
	if (position >= size())
		throw IndexOutOfBoundsException("PolymorphismMultiGContainer::setGroupId: position out of bounds.", position, 0, size() - 1);
	_groups[position] = group_id;
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
