/*
 * File: PolymorphismSequenceContainer.h
 * Authors: Eric Bazin <bazin@univ-montp2.fr>
 *          Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Wednesday August 04 2004
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

// from PolyLib
#include "PolymorphismSequenceContainer.h"

//** Class constructor: *******************************************************/
PolymorphismSequenceContainer::PolymorphismSequenceContainer(const Alphabet *alpha) : VectorSiteContainer(alpha) {}

PolymorphismSequenceContainer::PolymorphismSequenceContainer(unsigned int size, const Alphabet *alpha) : VectorSiteContainer(size, alpha) {
	_count.resize(size);
	_ingroup.resize(size);
	_group.resize(size);
}

PolymorphismSequenceContainer::PolymorphismSequenceContainer(const OrderedSequenceContainer & sc) : VectorSiteContainer(sc) {
	for (unsigned int i = 0 ; i < sc.getNumberOfSequences() ; i++) {
		_ingroup.push_back(true);
		_count.push_back(1);
		_group.push_back(1);
	}
}

PolymorphismSequenceContainer::PolymorphismSequenceContainer(const SiteContainer & sc) : VectorSiteContainer(sc) {
	for (unsigned int i = 0 ; i < sc.getNumberOfSequences() ; i++) {
		_ingroup.push_back(true);
		_count.push_back(1);
		_group.push_back(1);
	}
}

PolymorphismSequenceContainer::PolymorphismSequenceContainer(const PolymorphismSequenceContainer & psc) : VectorSiteContainer(psc) {
	unsigned int nbSeq = psc.getNumberOfSequences();
	_count.resize(nbSeq);
	_ingroup.resize(nbSeq);
	_group.resize(nbSeq);
	for(unsigned int i = 0; i < nbSeq; i++) {
		_count[i] = psc.getSequenceCount(i);
		_ingroup[i] = psc.isIngroupMember(i);
		_group[i] = psc.getGroupId(i);
	}
}

PolymorphismSequenceContainer & PolymorphismSequenceContainer::operator= (const PolymorphismSequenceContainer & psc) {
	// Setting up alphabet
	_alphabet = psc.getAlphabet();
	// Setting up general comments
	setGeneralComments(psc.getGeneralComments());
	// Setting up the sequences names
	_names.resize(psc.getNumberOfSequences());
	setSequencesNames(psc.getSequencesNames(), true);
	// Setting up the sites
	for(unsigned int i = 0; i < psc.getNumberOfSites(); i++) {
		addSite(* psc.getSite(i), true);
	}
	// Setting up the sequences comments, numbers and ingroup state
	unsigned int nbSeq = psc.getNumberOfSequences();
	_comments.resize(nbSeq);
	_count.resize(nbSeq);
	_ingroup.resize(nbSeq);
	_group.resize(nbSeq);
	for(unsigned int i = 0; i < nbSeq; i++) {
		_comments[i] = new Comments(psc.getComments(i));
		_count[i] = psc.getSequenceCount(i);
		_ingroup[i] = psc.isIngroupMember(i);
		_group[i] = psc.getGroupId(i);
	}
	Sequence * s = NULL;
	_sequences = vector<Sequence *>(nbSeq, s);
	return * this;
}

//** Class destructor: *******************************************************/
PolymorphismSequenceContainer::~PolymorphismSequenceContainer() {
	clear();
}

Clonable * PolymorphismSequenceContainer::clone() const {
	return dynamic_cast<Clonable *>(dynamic_cast<OrderedSequenceContainer *>
			(dynamic_cast<SiteContainer *>
			 (new PolymorphismSequenceContainer(*this))));
}

//** Other methodes: *********************************************************/

Sequence * PolymorphismSequenceContainer::removeSequence(unsigned int index) throw (IndexOutOfBoundsException) {
	if (index >= getNumberOfSequences())
		throw IndexOutOfBoundsException("PolymorphismSequenceContainer::removeSequence: index out of bounds.", index, 0, getNumberOfSequences());
	_count.erase(_count.begin() + index);
	_ingroup.erase(_ingroup.begin() + index);
	return VectorSiteContainer::removeSequence(index);
}

Sequence * PolymorphismSequenceContainer::removeSequence(const string &name) throw (SequenceNotFoundException) {
	try {
		return removeSequence(getSequencePosition(name));
	}
	catch (SequenceNotFoundException & snfe) {
		throw SequenceNotFoundException("PolymorphismSequenceContainer::removeSequence.", name);
	}
}

void PolymorphismSequenceContainer::deleteSequence(unsigned int index) throw (IndexOutOfBoundsException) {
	try {
		delete removeSequence(index);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("PolymorphismSequenceContainer::deleteSequence.", index, 0, getNumberOfSequences());
	}
}

void PolymorphismSequenceContainer::deleteSequence(const string &name) throw (SequenceNotFoundException) {
	try {
		delete removeSequence(name);
	}
	catch (SequenceNotFoundException & snfe) {
		throw SequenceNotFoundException("PolymorphismSequenceContainer::deleteSequence.", name);
	}
}

void PolymorphismSequenceContainer::addSequence(const Sequence &sequence, unsigned int effectif, bool checkNames) throw (Exception) {
	try {
		VectorSiteContainer::addSequence(sequence, checkNames);
	}
	catch (Exception & e) {
		throw e;
	}
	_count.push_back(effectif);
	_ingroup.push_back(true);
	_group.push_back(0);
}

void PolymorphismSequenceContainer::clear() {
	VectorSiteContainer::clear();
	_count.clear();
	_ingroup.clear();
	_group.clear();
}

unsigned int PolymorphismSequenceContainer::getGroupId(unsigned int index) const throw (IndexOutOfBoundsException) {
	if (index >= getNumberOfSequences())
		throw IndexOutOfBoundsException("PolymorphismSequenceContainer::getGroupId: index out of bounds.", index, 0, getNumberOfSequences());
	return _group[index];
}

unsigned int PolymorphismSequenceContainer::getGroupId(const string &name) const throw (SequenceNotFoundException) {
	try {
		return _group[getSequencePosition(name)];
	}
	catch (SequenceNotFoundException & snfe) {
		throw SequenceNotFoundException("PolymorphismSequenceContainer::getGroupId.", name);
	}
}

set<unsigned int> PolymorphismSequenceContainer::getAllGroupsIds() const {
	set<unsigned int> grp_ids;
	for (unsigned int i = 0 ; i < _group.size() ; i++)
		grp_ids.insert(_group[i]);
	return grp_ids;
}

void PolymorphismSequenceContainer::setGroupId(unsigned int index, unsigned int group_id) throw (IndexOutOfBoundsException) {
	if (index >= getNumberOfSequences())
		throw IndexOutOfBoundsException("PolymorphismSequenceContainer::setGroupId: index out of bounds.", index, 0, getNumberOfSequences());
	_group[index] = group_id;
}

void PolymorphismSequenceContainer::setGroupId(const string &name, unsigned int group_id) throw (SequenceNotFoundException) {
	try {
		_group[getSequencePosition(name)] = group_id;
	}
	catch (SequenceNotFoundException & snfe) {
		throw SequenceNotFoundException("PolymorphismSequenceContainer::setGroupId.", name);
	}
}

unsigned int PolymorphismSequenceContainer::getNumberOfGroups() const {
	return getAllGroupsIds().size();
}

bool PolymorphismSequenceContainer::isIngroupMember(unsigned int index) const throw (IndexOutOfBoundsException) {
	if (index >= getNumberOfSequences())
		throw IndexOutOfBoundsException("PolymorphismSequenceContainer::isIngroupMember: index out of bounds.", index, 0, getNumberOfSequences());
	return _ingroup[index];
}

bool PolymorphismSequenceContainer::isIngroupMember(const string &name) const throw (SequenceNotFoundException) {
	try {
		return _ingroup[getSequencePosition(name)];
	}
	catch (SequenceNotFoundException & snfe) {
		throw SequenceNotFoundException("PolymorphismSequenceContainer::isIngroupMember.", name);
	}
}

void PolymorphismSequenceContainer::setAsIngroupMember(unsigned int index) throw (IndexOutOfBoundsException) {
	if (index >= getNumberOfSequences())
		throw IndexOutOfBoundsException("PolymorphismSequenceContainer::setAsIngroupMember.", index, 0, getNumberOfSequences());
	_ingroup[index] = true;
}

void PolymorphismSequenceContainer::setAsIngroupMember(const string &name) throw (SequenceNotFoundException) {
	try {
		unsigned int seqPos = getSequencePosition(name);
		_ingroup[seqPos] = true;
	}
	catch (SequenceNotFoundException & snfe) {
		throw SequenceNotFoundException("PolymorphismSequenceContainer::setAsIngroupMember.", name);
	}
}

void PolymorphismSequenceContainer::setAsOutgroupMember(unsigned int index) throw (IndexOutOfBoundsException) {
	if (index >= getNumberOfSequences())
		throw IndexOutOfBoundsException("PolymorphismSequenceContainer::setAsOutgroupMember.", index, 0, getNumberOfSequences());
	_ingroup[index] = false;
}

void PolymorphismSequenceContainer::setAsOutgroupMember(const string &name) throw (SequenceNotFoundException) {
	try {
		unsigned int seqPos = getSequencePosition(name);
		_ingroup[seqPos] = false;
	}
	catch (SequenceNotFoundException & snfe) {
		throw SequenceNotFoundException("PolymorphismSequenceContainer::setAsOutgroupMember.", name);
	}
}

void PolymorphismSequenceContainer::setSequenceCount(unsigned int index, unsigned int count) throw (Exception) {
	if (index >= getNumberOfSequences())
		throw IndexOutOfBoundsException("PolymorphismSequenceContainer::setSequenceCount.", index, 0, getNumberOfSequences());
	if (count < 1)
		throw BadIntegerException("PolymorphismSequenceContainer::setSequenceCount: count can't be < 1.", count);
	_count[index] = count;
}

void PolymorphismSequenceContainer::setSequenceCount(const string &name, unsigned int count) throw (Exception) {
	try {
		setSequenceCount(getSequencePosition(name), count);
	}
	catch (BadIntegerException & bie) {
		throw bie;
	}
	catch (SequenceNotFoundException & snfe) {
		throw SequenceNotFoundException("PolymorphismSequenceContainer::setSequenceCount.", name);
	}
}

void PolymorphismSequenceContainer::incrementSequenceCount(unsigned int index) throw (IndexOutOfBoundsException) {
	if (index >= getNumberOfSequences())
		throw IndexOutOfBoundsException("PolymorphismSequenceContainer::incrementSequenceCount.", index, 0, getNumberOfSequences());
	_count[index]++;
}

void PolymorphismSequenceContainer::incrementSequenceCount(const string &name) throw (SequenceNotFoundException) {
	try {
		incrementSequenceCount(getSequencePosition(name));
	}
	catch (SequenceNotFoundException & snfe) {
		throw SequenceNotFoundException("PolymorphismSequenceContainer::incrementSequenceCount.", name);
	}
}

void PolymorphismSequenceContainer::decrementSequenceCount(unsigned int index) throw (Exception) {
	if (index >= getNumberOfSequences())
		throw IndexOutOfBoundsException("PolymorphismSequenceContainer::decrementSequenceCount.", index, 0, getNumberOfSequences());
	if (_count[index]-1 < 1)
		throw BadIntegerException("PolymorphismSequenceContainer::decrementSequenceCount: count can't be < 1.", _count[index]-1);
	_count[index]--;
}

void PolymorphismSequenceContainer::decrementSequenceCount(const string &name) throw (Exception) {
	try {
		decrementSequenceCount(getSequencePosition(name));
	}
	catch (BadIntegerException & bie) {
		throw bie;
	}
	catch (SequenceNotFoundException & snfe) {
		throw SequenceNotFoundException("PolymorphismSequenceContainer::decrementSequenceCount.", name);
	}
}

unsigned int PolymorphismSequenceContainer::getSequenceCount(unsigned int index) const throw (IndexOutOfBoundsException) {
	if (index >= getNumberOfSequences())
		throw IndexOutOfBoundsException("PolymorphismSequenceContainer::getSequenceCount.", index, 0, getNumberOfSequences());
	return _count[index];
}

unsigned int PolymorphismSequenceContainer::getSequenceCount(const string &name) const throw (SequenceNotFoundException) {
	try {
		return getSequenceCount(getSequencePosition(name));
	}
	catch (SequenceNotFoundException & snfe) {
		throw SequenceNotFoundException("PolymorphismSequenceContainer::getSequenceCount.", name);
	}
}

unsigned int PolymorphismSequenceContainer::getPhase(const string &setName) const throw (Exception) {
	try {
		unsigned int phase;
		unsigned int index;
		Comments maseFileHeader = getGeneralComments();
		for(unsigned int i = 0; i < maseFileHeader.size(); i++) {
			string current = maseFileHeader[i];
			
			index = current.find("# of regions");
			if(index < current.npos) {
				StringTokenizer st(string(current.begin() + index + 12 , current.end()), " \t\n\f\r=;");
				unsigned int numberOfSegments = TextTools::toInt(st.nextToken());
				//cout << "Number of regions: " << st.nextToken() << endl;
				string name;
				while(st.hasMoreToken()) {
					name = st.nextToken();
					//cout << "Name of regions: " << name << endl;
				}
				if(name == setName) {
					return phase;
				}
			}
			
			index = current.find("/codon_start");
			if(index < current.npos) {
				StringTokenizer st(string(current.begin() + index + 12, current.end()), " \t\n\f\r=;");
				phase = TextTools::toInt(st.nextToken());
			}
		}
		cout << "No phase for " << setName << endl;
	}
	catch (...) {
	}
}
