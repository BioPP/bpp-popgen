/*
 * File Group.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Thursday July 29 2004
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

#include "Group.h"

//** Class constructors: ******************************************************/
Group::Group(unsigned int group_id) {
	setGroupId(group_id);
}

Group::Group(const Group & group) {
	setGroupId(group.getGroupId());
	for (unsigned int i = 0 ; i < group.getNumberOfIndividuals() ; i++)
		addIndividual(* (group.getIndividualAtPosition(i)));
}

Group::Group(const Group & group, unsigned int group_id) {
	setGroupId(group_id);
	for (unsigned int i = 0 ; i < group.getNumberOfIndividuals() ; i++)
		addIndividual(* (group.getIndividualAtPosition(i)));
}

//** Class destructor: ********************************************************/
Group::~Group () {}

//** Other methodes: **********************************************************/
Group & Group::operator= (const Group & group) {
	setGroupId(group.getGroupId());
	for (unsigned int i = 0 ; i < group.getNumberOfIndividuals() ; i++)
		addIndividual(* (group.getIndividualAtPosition(i)));
	return * this;
}

void Group::setGroupId(unsigned int group_id) {
	_id = group_id;
}

unsigned int Group::getGroupId() const {
	return _id;
}

void Group::addIndividual(const Individual & ind) throw (BadIdentifierException) {
	try {
		getIndividualPosition(ind.getId());
		throw BadIdentifierException("Group::addIndividual: individual id already used.", ind.getId());
	}
	catch (BadIdentifierException & bie) {}
	_individuals.push_back(new Individual(ind));
}

void Group::addEmptyIndividual(const string & individual_id) throw (BadIdentifierException) {
	for (unsigned int i = 0 ; i < getNumberOfIndividuals() ; i++)
		if (_individuals[i]->getId() == individual_id)
			throw BadIdentifierException("Group::addEmptyIndividual: individual_id already in use.", individual_id);
	_individuals.push_back(new Individual(individual_id));
}

unsigned int Group::getIndividualPosition(const string individual_id) const throw (IndividualNotFoundException) {
	for (unsigned int i = 0 ; i < getNumberOfIndividuals() ; i++)
		if (_individuals[i]->getId() == individual_id)
			return i;
	throw IndividualNotFoundException("Group::getIndividualPosition: individual_id not found.", individual_id);
}

Individual * Group::removeIndividualById(const string individual_id) throw (IndividualNotFoundException) {
	try {
		unsigned int indPos = getIndividualPosition(individual_id);
		Individual * ind = _individuals[indPos];
		_individuals.erase(_individuals.begin() + indPos);
		return ind;
	}
	catch (IndividualNotFoundException & infe) {
		throw IndividualNotFoundException("Group::removeIndividualById: individual_id not found.", individual_id);
	}
}

Individual * Group::removeIndividualAtPosition(unsigned int individual_position) throw (IndexOutOfBoundsException) {
	if (individual_position >= _individuals.size())
		throw IndexOutOfBoundsException("Group::removeIndividualAtPosition.", individual_position, 0, _individuals.size());
	Individual * ind = _individuals[individual_position];
	_individuals.erase(_individuals.begin() + individual_position);
	return ind;
}

void Group::deleteIndividualById(const string individual_id) throw (IndividualNotFoundException) {
	try {
		Individual * ind = removeIndividualById(individual_id);
		delete ind;
	}
	catch (IndividualNotFoundException & infe) {
		throw IndividualNotFoundException("Group::deleteIndividualById: individual_id not found.", individual_id);
	}
}

void Group::deleteIndividualAtPosition(unsigned int individual_position) throw (IndexOutOfBoundsException) {
	try {
		Individual * ind = removeIndividualAtPosition(individual_position);
		delete ind;
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("Group::deleteIndividualAtPosition.", individual_position, 0, getNumberOfIndividuals());
	}
}

void Group::clear() {
	for (unsigned int i = 0 ; i < _individuals.size() ; i++)
		delete(_individuals[i]);
	_individuals.clear();
}

const Individual * Group::getIndividualById(const string individual_id) const {
	for (unsigned int i = 0 ; i < _individuals.size() ; i++) {
		if (_individuals[i]->getId() == individual_id)
			return getIndividualAtPosition(i);
	}
	return NULL;
}

const Individual * Group::getIndividualAtPosition(unsigned int individual_position) const
throw (IndexOutOfBoundsException) {
	if (individual_position >= _individuals.size())
		throw IndexOutOfBoundsException("Group::getIndividualAtPosition: individual_position out of bounds.", individual_position, 0, _individuals.size());
	return _individuals[individual_position];
}

unsigned int Group::getNumberOfIndividuals() const {
	return _individuals.size();
}

unsigned int Group::getMaxNumberOfSequences() const {
	unsigned int maxnum = 0;
	for (unsigned int i = 0 ; i < getNumberOfIndividuals() ; i++) {
		vector<unsigned int> seqpos = _individuals[i]->getSequencesPositions();
		for (unsigned int j = 0 ; j < seqpos.size() ; j++)
			if (maxnum < seqpos[j]) maxnum = seqpos[j];
	}
	return maxnum + 1;
}

//-- Dealing with individual's properties -----------------
void Group::setIndividualSexAtPosition(unsigned int individual_position, const unsigned short sex) throw (IndexOutOfBoundsException) {
	if (individual_position >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::setIndividualSexAtPosition: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
	_individuals[individual_position]->setSex(sex);
}

unsigned short Group::getIndividualSexAtPosition(unsigned int individual_position) const throw (IndexOutOfBoundsException) {
	if (individual_position >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::getIndividualSexAtPosition: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
	return _individuals[individual_position]->getSex();
}

void Group::setIndividualDateAtPosition(unsigned int individual_position, const Date & date) throw (IndexOutOfBoundsException) {
	if (individual_position >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::setIndividualDateAtPosition: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
	_individuals[individual_position]->setDate(date);
}

const Date * Group::getIndividualDateAtPosition(unsigned int individual_position) const throw (Exception) {
	if (individual_position >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::getIndividualDateAtPosition: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
	try {
		return _individuals[individual_position]->getDate();
	}
	catch (NullPointerException & npe) {
		throw NullPointerException("Group::getIndividualDateAtPosition: individual has no date.");
	}
}

void Group::setIndividualCoordAtPosition(unsigned int individual_position, const Coord<double> & coord) throw (IndexOutOfBoundsException) {
	if (individual_position >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::setIndividualCoordAtPosition: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
	_individuals[individual_position]->setCoord(coord);
}

const Coord<double> * Group::getIndividualCoordAtPosition(unsigned int individual_position) const throw (Exception) {
	if (individual_position >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::getIndividualCoordAtPosition: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
	try {
		return _individuals[individual_position]->getCoord();
	}
	catch (NullPointerException & npe) {
		throw NullPointerException("Group::getIndividualCoordAtPosition: individual has no coordinates.");
	}
}

void Group::setIndividualLocalityAtPosition(unsigned int individual_position, const Locality<double> * locality) throw (IndexOutOfBoundsException) {
	if (individual_position >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::setIndividualLocalityAtPosition: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
	_individuals[individual_position]->setLocality(locality);
}

const Locality<double> * Group::getIndividualLocalityAtPosition(unsigned int individual_position) const throw (Exception) {
	if (individual_position >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::getIndividualLocalityAtPosition: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
	try {
		return _individuals[individual_position]->getLocality();
	}
	catch (NullPointerException & npe) {
		throw NullPointerException("Group::getIndividualLocalityAtPosition: individuals has no locality.");
	}
}

void Group::addIndividualSequenceAtPosition(unsigned int individual_position, unsigned int sequence_position, const Sequence & sequence) throw (Exception) {
	if (individual_position >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::addIndividualSequenceAtPosition: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
	try {
		_individuals[individual_position]->addSequence(sequence_position, sequence);
	}
	catch (AlphabetMismatchException & ame) {
		throw AlphabetMismatchException("Group::addIndividualSequenceAtPosition: sequence's alphabet doesn't match.", ame.getAlphabets()[0], ame.getAlphabets()[1]);
	}
	catch (BadIdentifierException & bie) {
		throw BadIdentifierException("Group::addIndividualSequenceAtPosition: sequence's name already in use.", bie.getIdentifier());
	}
	catch (BadIntegerException & bie) {
		throw BadIntegerException("Group::addIndividualSequenceAtPosition: sequence_position already in use.", bie.getBadInteger());
	}
}

const Sequence * Group::getIndividualSequenceByName(unsigned int individual_position, const string & sequence_name) const throw (Exception) {
	if (individual_position >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::getIndividualSequenceByName: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
	try {
		return _individuals[individual_position]->getSequenceByName(sequence_name);
	}
	catch (NullPointerException & npe) {
		throw NullPointerException("Group::getIndividualSequenceByName: no sequence data in individual.");
	}
	catch (SequenceNotFoundException & snfe) {
		throw SequenceNotFoundException("Group::getIndividualSequenceByName: sequence_name not found.", snfe.getSequenceId());
	}
}

const Sequence * Group::getIndividualSequenceAtPosition(unsigned int individual_position, unsigned int sequence_position) const throw (Exception) {
	if (individual_position >= getNumberOfIndividuals())
		 throw IndexOutOfBoundsException("Group::getIndividualAtPosition: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
	try {
		return _individuals[individual_position]->getSequenceAtPosition(sequence_position);
	}
	catch (NullPointerException & npe) {
		throw NullPointerException("Group::getIndividualSequenceAtPosition: no sequence data in individual.");
	}
	catch (SequenceNotFoundException & snfe) {
		throw SequenceNotFoundException("Group::getIndividualSequenceAtPosition: sequence_position not found.", snfe.getSequenceId());
	}
}

void Group::deleteIndividualSequenceByName(unsigned int individual_position, const string & sequence_name) throw (Exception) {
	if (individual_position >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::deleteIndividualSequenceByName: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
	try {
		_individuals[individual_position]->deleteSequenceByName(sequence_name);
	}
	catch (NullPointerException & npe) {
		throw NullPointerException("Group::deleteSequenceByName: no sequence data in individual.");
	}
	catch (SequenceNotFoundException & snfe) {
		throw SequenceNotFoundException("Group::deleteSequenceByName: sequence_name not found.", snfe.getSequenceId());
	}
}

void Group::deleteIndividualSequenceAtPosition(unsigned int individual_position, unsigned int sequence_position) throw (Exception) {
	if (individual_position >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::deleteIndividualSequenceAtPosition: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
	try {
		_individuals[individual_position]->deleteSequenceAtPosition(sequence_position);
	}
	catch (NullPointerException & npe) {
		throw NullPointerException("Group::deleteSequenceAtPosition: no sequence data in individual.");
	}
	catch (SequenceNotFoundException & snfe) {
		throw SequenceNotFoundException("Group::deleteSequenceAtPosition: sequence_position not found.", snfe.getSequenceId());
	}
}

bool Group::hasIndividualSequences(unsigned int individual_position) const throw (IndexOutOfBoundsException) {
	if (individual_position >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::hasIndividualSequences: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
	return _individuals[individual_position]->hasSequences();
}

vector<string> Group::getIndividualSequencesNames(unsigned int individual_position) const throw (Exception) {
	if (individual_position >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::getIndividualSequencesNames: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
	try {
		return _individuals[individual_position]->getSequencesNames();
	}
	catch (NullPointerException & npe) {
		throw NullPointerException("Group::getSequencesNames: no sequence data in individual.");
	}
}

unsigned int Group::getIndividualSequencePosition(unsigned int individual_position, const string & sequence_name) const throw (Exception) {
	if (individual_position >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::getIndividualSequencePosition: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
	try {
		return _individuals[individual_position]->getSequencePosition(sequence_name);
	}
	catch (NullPointerException & npe) {
		throw NullPointerException("Group::getSequencePosition: no sequence data in individual.");
	}
	catch (SequenceNotFoundException & snfe) {
		throw SequenceNotFoundException("Group::getSequencePosition: sequence_name not found.", snfe.getSequenceId());
	}
}

unsigned int Group::getIndividualNumberOfSequences(unsigned int individual_position) const throw (Exception) {
	if (individual_position >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::getIndividualNumberOfSequences: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
	try {
		return _individuals[individual_position]->getNumberOfSequences();
	}
	catch (NullPointerException & npe) {
		throw NullPointerException("Group::getIndividualNumberOfSequences: no sequence data in individual.");
	}
}

void Group::setIndividualSequences(unsigned int individual_position, const MapSequenceContainer & msc) throw (IndexOutOfBoundsException) {
	if (individual_position >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::setIndividualSequences: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
	_individuals[individual_position]->setSequences(msc);
}

void Group::setIndividualGenotype(unsigned int individual_position, const MultilocusGenotype & genotype) throw (IndexOutOfBoundsException) {
	if (individual_position >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::setIndividualGenotype: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
	_individuals[individual_position]->setGenotype(genotype);
}

void Group::initIndividualGenotype(unsigned int individual_position, unsigned int loci_number) throw (Exception) {
	if (individual_position >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::initIndividualGenotype: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
	try {
		_individuals[individual_position]->initGenotype(loci_number);
	}
	catch (BadIntegerException & bie) {
		throw BadIntegerException("Group::initIndividualGenotype: loci_number must be > 0.", bie.getBadInteger());
	}
	catch (Exception) {
		throw Exception("Group::initIndividualGenotype: individual already has a genotype.");
	}
}

void Group::deleteIndividualGenotype(unsigned int individual_position) throw (IndexOutOfBoundsException) {
	if (individual_position >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::deleteIndividualGenotype: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
	_individuals[individual_position]->deleteGenotype();
}

bool Group::hasIndividualGenotype(unsigned int individual_position) const throw (IndexOutOfBoundsException) {
	if (individual_position >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::hasIndividualGenotype: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
	return _individuals[individual_position]->hasGenotype();
}

void Group::setIndividualMonolocusGenotype(unsigned int individual_position, unsigned int locus_position, const MonolocusGenotype & monogen) throw (Exception) {
	if (individual_position >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::setIndividualMonolocusGenotype: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
	try {
		_individuals[individual_position]->setMonolocusGenotype(locus_position, monogen);
	}
	catch (NullPointerException & npe) {
		throw NullPointerException("Group::setIndividualMonolocusGenotype: individual has no genotype.");
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("Group::setIndividualMonolocusGenotype: locus_position excedes the number of locus.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
}

void Group::setIndividualMonolocusGenotypeByAlleleKey(unsigned int individual_position, unsigned int locus_position, const vector<unsigned int> allele_keys) throw (Exception) {
	if (individual_position >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::setIndividualMonolocusGenotypeByAlleleKey: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
	try {
		_individuals[individual_position]->setMonolocusGenotypeByAlleleKey(locus_position, allele_keys);
	}
	catch (NullPointerException & npe) {
		throw NullPointerException("Group::setIndividualMonolocusGenotypeByAlleleKey: individual has no genotype.");
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("Group::setIndividualMonolocusGenotypeByAlleleKey: locus_position excedes the number of locus.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
	catch (Exception) {
		throw Exception("Group::setIndividualMonolocusGenotypeByAlleleKey: no key in allele_keys.");
	}
}

void Group::setIndividualMonolocusGenotypeByAlleleId(unsigned int individual_position, unsigned int locus_position, const vector<string> allele_id, const LocusInfo & locus_info) throw (Exception) {
	if (individual_position >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::setIndividualMonolocusGenotypeByAlleleId: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
	try {
		_individuals[individual_position]->setMonolocusGenotypeByAlleleId(locus_position, allele_id, locus_info);
	}
	catch (NullPointerException & npe) {
		throw NullPointerException("Group::setIndividualMonolocusGenotypeByAlleleId: individual has no genotype.");
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("Group::setIndividualMonolocusGenotypeByAlleleId: locus_position excedes the number of locus.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
	catch (AlleleNotFoundException & anfe) {
		throw AlleleNotFoundException("Group::setIndividualMonolocusGenotypeByAlleleId: id not found.", anfe.getIdentifier());
	}
}

const MonolocusGenotype *  Group::getIndividualMonolocusGenotype(unsigned int individual_position, unsigned int locus_position) const throw (Exception) {
	if (individual_position >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::getIndividualMonolocusGenotype: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
	try {
		return _individuals[individual_position]->getMonolocusGenotype(locus_position);
	}
	catch (NullPointerException & npe) {
		throw NullPointerException("Group::getIndividualMonolocusGenotype: individual has no genotype.");
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("Group::getIndividualMonolocusGenotype: locus_position excedes the number of locus.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
}

bool Group::hasSequenceData() const {
	for (unsigned int i = 0 ; i < getNumberOfIndividuals() ; i++)
		if (hasIndividualSequences(i)) return true;
	return false;
}

const Alphabet * Group::getAlphabet() const throw (NullPointerException) {
	for (unsigned int i = 0 ; i < getNumberOfIndividuals() ; i++)
		if (hasIndividualSequences(i))
			return _individuals[i]->getSequenceAlphabet();
	throw NullPointerException("Group::getAlphabet: individual has no sequence data.");
}

unsigned int Group::getGroupSizeForLocus(unsigned int locus_position) const {
	unsigned int count = 0;
	for (unsigned int i = 0 ; i < _individuals.size() ; i++)
		if (_individuals[i]->hasGenotype() && !_individuals[i]->getGenotype()->isMonolocusGenotypeMissing(locus_position))
			count++;
	return count;
}

unsigned int Group::getGroupSizeForSequence(unsigned int sequence_position) const {
	unsigned int count = 0;
	for (unsigned int i = 0 ; i < _individuals.size() ; i++)
		if (_individuals[i]->hasSequences()) {
			try {
				_individuals[i]->getSequenceAtPosition(sequence_position);
				count++;
			}
			catch (...) {}
		}
	return count;
}
