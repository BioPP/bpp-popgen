/*
 * File Group.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Monday July 05 2004
 */

#include "Group.h"

//** Class constructor: *******************************************************/
Group::Group() {}

//** Class destructor: ********************************************************/
Group::~Group () {}

//** Other methodes: **********************************************************/
void Group::setGroupId(const string group_id) {
	_id = group_id;
}

string Group::getGroupId() const {
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

Individual * Group::removeIndividualByIndex(unsigned int index) throw (IndexOutOfBoundsException) {
	if (index >= _individuals.size())
		throw IndexOutOfBoundsException("Group::removeIndividualByIndex.", index, 0, _individuals.size());
	Individual * ind = _individuals[index];
	_individuals.erase(_individuals.begin() + index);
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

void Group::deleteIndividualByIndex(unsigned int index) throw (IndexOutOfBoundsException) {
	try {
		Individual * ind = removeIndividualByIndex(index);
		delete ind;
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("Group::deleteIndividualByIndex.", index, 0, getNumberOfIndividuals());
	}
}

void Group::clear() {
	for (unsigned int i = 0 ; i < _individuals.size() ; i++)
		delete(_individuals[i]);
	_individuals.clear();
}

void Group::append(const Group & group) {
	for (unsigned int i = 0 ; i < group.getNumberOfIndividuals() ; i++)
		addIndividual(*(group.getIndividualByIndex(i)));
}

const Individual * Group::getIndividualById(const string individual_id) const
throw (IndividualNotFoundException) {
	for (unsigned int i = 0 ; i < _individuals.size() ; i++) {
		if (_individuals[i]->getId() == individual_id)
			return getIndividualByIndex(i);
	}
	throw IndividualNotFoundException("Group::getIndividualById: individual_id not found.", individual_id);
}

const Individual * Group::getIndividualByIndex(unsigned int index) const
throw (IndexOutOfBoundsException) {
	if (index >= _individuals.size())
		throw IndexOutOfBoundsException("Group::getIndividualByIndex: index out of bounds.", index, 0, _individuals.size());
	return _individuals[index];
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
void Group::setIndividualSexByIndex(unsigned int individual_index, const unsigned short sex) throw (IndexOutOfBoundsException) {
	if (individual_index >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::setIndividualSexByIndex: individual_index out of bounds.", individual_index, 0, getNumberOfIndividuals());
	_individuals[individual_index]->setSex(sex);
}

unsigned short Group::getIndividualSexByIndex(unsigned int individual_index) const throw (IndexOutOfBoundsException) {
	if (individual_index >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::getIndividualSexByIndex: individual_index out of bounds.", individual_index, 0, getNumberOfIndividuals());
	return _individuals[individual_index]->getSex();
}

void Group::setIndividualDateByIndex(unsigned int individual_index, const Date & date) throw (IndexOutOfBoundsException) {
	if (individual_index >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::setIndividualDateByIndex: individual_index out of bounds.", individual_index, 0, getNumberOfIndividuals());
	_individuals[individual_index]->setDate(date);
}

const Date * Group::getIndividualDateByIndex(unsigned int individual_index) const throw (Exception) {
	if (individual_index >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::getIndividualDateByIndex: individual_index out of bounds.", individual_index, 0, getNumberOfIndividuals());
	try {
		return _individuals[individual_index]->getDate();
	}
	catch (NullPointerException & npe) {
		throw NullPointerException("Group::getIndividualDateByIndex: individual has no date.");
	}
}

void Group::setIndividualCoordByIndex(unsigned int individual_index, const Coord<double> & coord) throw (IndexOutOfBoundsException) {
	if (individual_index >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::setIndividualCoordByIndex: individual_index out of bounds.", individual_index, 0, getNumberOfIndividuals());
	_individuals[individual_index]->setCoord(coord);
}

const Coord<double> * Group::getIndividualCoordByIndex(unsigned int individual_index) const throw (Exception) {
	if (individual_index >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::getIndividualCoordByIndex: individual_index out of bounds.", individual_index, 0, getNumberOfIndividuals());
	try {
		return _individuals[individual_index]->getCoord();
	}
	catch (NullPointerException & npe) {
		throw NullPointerException("Group::getIndividualCoordByIndex: individual has no coordinates.");
	}
}

void Group::setIndividualLocalityByIndex(unsigned int individual_index, const Locality<double> * locality) throw (IndexOutOfBoundsException) {
	if (individual_index >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::setIndividualLocalityByIndex: individual_index out of bounds.", individual_index, 0, getNumberOfIndividuals());
	_individuals[individual_index]->setLocality(locality);
}

const Locality<double> * Group::getIndividualLocalityByIndex(unsigned int individual_index) const throw (Exception) {
	if (individual_index >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::getIndividualLocalityByIndex: individual_index out of bounds.", individual_index, 0, getNumberOfIndividuals());
	try {
		return _individuals[individual_index]->getLocality();
	}
	catch (NullPointerException & npe) {
		throw NullPointerException("Group::getIndividualLocalityByIndex: individuals has no locality.");
	}
}

void Group::addIndividualSequenceByIndex(unsigned int individual_index, unsigned int sequence_index, const Sequence & sequence) throw (Exception) {
	if (individual_index >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::addIndividualSequenceByIndex: individual_index out of bounds.", individual_index, 0, getNumberOfIndividuals());
	try {
		_individuals[individual_index]->addSequence(sequence_index, sequence);
	}
	catch (AlphabetMismatchException & ame) {
		throw AlphabetMismatchException("Group::addIndividualSequenceByIndex: sequence's alphabet doesn't match.", ame.getAlphabets()[0], ame.getAlphabets()[1]);
	}
	catch (BadIdentifierException & bie) {
		throw BadIdentifierException("Group::addIndividualSequenceByIndex: sequence's name already in use.", bie.getIdentifier());
	}
	catch (BadIntegerException & bie) {
		throw BadIntegerException("Group::addIndividualSequenceByIndex: sequence_index already in use.", bie.getBadInteger());
	}
}

const Sequence * Group::getIndividualSequenceByName(unsigned int individual_index, const string & sequence_name) const throw (Exception) {
	if (individual_index >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::getIndividualSequenceByName: individual_index out of bounds.", individual_index, 0, getNumberOfIndividuals());
	try {
		return _individuals[individual_index]->getSequenceByName(sequence_name);
	}
	catch (NullPointerException & npe) {
		throw NullPointerException("Group::getIndividualSequenceByName: no sequence data in individual.");
	}
	catch (SequenceNotFoundException & snfe) {
		throw SequenceNotFoundException("Group::getIndividualSequenceByName: sequence_name not found.", snfe.getSequenceId());
	}
}

const Sequence * Group::getIndividualSequenceByIndex(unsigned int individual_index, unsigned int sequence_index) const throw (Exception) {
	if (individual_index >= getNumberOfIndividuals())
		 throw IndexOutOfBoundsException("Group::getIndividualByIndex: individual_index out of bounds.", individual_index, 0, getNumberOfIndividuals());
	try {
		return _individuals[individual_index]->getSequenceByIndex(sequence_index);
	}
	catch (NullPointerException & npe) {
		throw NullPointerException("Group::getIndividualSequenceByIndex: no sequence data in individual.");
	}
	catch (SequenceNotFoundException & snfe) {
		throw SequenceNotFoundException("Group::getIndividualSequenceByIndex: sequence_index not found.", snfe.getSequenceId());
	}
}

void Group::deleteIndividualSequenceByName(unsigned int individual_index, const string & sequence_name) throw (Exception) {
	if (individual_index >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::deleteIndividualSequenceByName: individual_index out of bounds.", individual_index, 0, getNumberOfIndividuals());
	try {
		_individuals[individual_index]->deleteSequenceByName(sequence_name);
	}
	catch (NullPointerException & npe) {
		throw NullPointerException("Group::deleteSequenceByName: no sequence data in individual.");
	}
	catch (SequenceNotFoundException & snfe) {
		throw SequenceNotFoundException("Group::deleteSequenceByName: sequence_name not found.", snfe.getSequenceId());
	}
}

void Group::deleteIndividualSequenceByIndex(unsigned int individual_index, unsigned int sequence_index) throw (Exception) {
	if (individual_index >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::deleteIndividualSequenceByIndex: individual_index out of bounds.", individual_index, 0, getNumberOfIndividuals());
	try {
		_individuals[individual_index]->deleteSequenceByIndex(sequence_index);
	}
	catch (NullPointerException & npe) {
		throw NullPointerException("Group::deleteSequenceByIndex: no sequence data in individual.");
	}
	catch (SequenceNotFoundException & snfe) {
		throw SequenceNotFoundException("Group::deleteSequenceByIndex: sequence_index not found.", snfe.getSequenceId());
	}
}

bool Group::hasIndividualSequences(unsigned int individual_index) const throw (IndexOutOfBoundsException) {
	if (individual_index >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::hasIndividualSequences: individual_index out of bounds.", individual_index, 0, getNumberOfIndividuals());
	return _individuals[individual_index]->hasSequences();
}

vector<string> Group::getIndividualSequencesNames(unsigned int individual_index) const throw (Exception) {
	if (individual_index >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::getIndividualSequencesNames: individual_index out of bounds.", individual_index, 0, getNumberOfIndividuals());
	try {
		return _individuals[individual_index]->getSequencesNames();
	}
	catch (NullPointerException & npe) {
		throw NullPointerException("Group::getSequencesNames: no sequence data in individual.");
	}
}

unsigned int Group::getIndividualSequencePosition(unsigned int individual_index, const string & sequence_name) const throw (Exception) {
	if (individual_index >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::getIndividualSequencePosition: individual_index out of bounds.", individual_index, 0, getNumberOfIndividuals());
	try {
		return _individuals[individual_index]->getSequencePosition(sequence_name);
	}
	catch (NullPointerException & npe) {
		throw NullPointerException("Group::getSequencePosition: no sequence data in individual.");
	}
	catch (SequenceNotFoundException & snfe) {
		throw SequenceNotFoundException("Group::getSequencePosition: sequence_name not found.", snfe.getSequenceId());
	}
}

unsigned int Group::getIndividualNumberOfSequences(unsigned int individual_index) const throw (Exception) {
	if (individual_index >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::getIndividualNumberOfSequences: individual_index out of bounds.", individual_index, 0, getNumberOfIndividuals());
	try {
		return _individuals[individual_index]->getNumberOfSequences();
	}
	catch (NullPointerException & npe) {
		throw NullPointerException("Group::getIndividualNumberOfSequences: no sequence data in individual.");
	}
}

void Group::setIndividualSequences(unsigned int individual_index, const MapSequenceContainer & msc) throw (IndexOutOfBoundsException) {
	if (individual_index >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::setIndividualSequences: individual_index out of bounds.", individual_index, 0, getNumberOfIndividuals());
	_individuals[individual_index]->setSequences(msc);
}

void Group::setIndividualGenotype(unsigned int individual_index, const MultilocusGenotype & genotype) throw (IndexOutOfBoundsException) {
	if (individual_index >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::setIndividualGenotype: individual_index out of bounds.", individual_index, 0, getNumberOfIndividuals());
	_individuals[individual_index]->setGenotype(genotype);
}

void Group::initIndividualGenotype(unsigned int individual_index, unsigned int loci_number) throw (Exception) {
	if (individual_index >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::initIndividualGenotype: individual_index out of bounds.", individual_index, 0, getNumberOfIndividuals());
	try {
		_individuals[individual_index]->initGenotype(loci_number);
	}
	catch (BadIntegerException & bie) {
		throw BadIntegerException("Group::initIndividualGenotype: loci_number must be > 0.", bie.getBadInteger());
	}
	catch (Exception) {
		throw Exception("Group::initIndividualGenotype: individual already has a genotype.");
	}
}

void Group::deleteIndividualGenotype(unsigned int individual_index) throw (IndexOutOfBoundsException) {
	if (individual_index >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::deleteIndividualGenotype: individual_index out of bounds.", individual_index, 0, getNumberOfIndividuals());
	_individuals[individual_index]->deleteGenotype();
}

bool Group::hasIndividualGenotype(unsigned int individual_index) const throw (IndexOutOfBoundsException) {
	if (individual_index >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::hasIndividualGenotype: individual_index out of bounds.", individual_index, 0, getNumberOfIndividuals());
	_individuals[individual_index]->hasGenotype();
}

void Group::setIndividualMonolocusGenotype(unsigned int individual_index, unsigned int locus_index, const MonolocusGenotype & monogen) throw (Exception) {
	if (individual_index >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::setIndividualMonolocusGenotype: individual_index out of bounds.", individual_index, 0, getNumberOfIndividuals());
	try {
		_individuals[individual_index]->setMonolocusGenotype(locus_index, monogen);
	}
	catch (NullPointerException & npe) {
		throw NullPointerException("Group::setIndividualMonolocusGenotype: individual has no genotype.");
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("Group::setIndividualMonolocusGenotype: locus_index excedes the number of locus.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
}

void Group::setIndividualMonolocusGenotypeByAlleleKey(unsigned int individual_index, unsigned int locus_index, const vector<unsigned int> allele_keys) throw (Exception) {
	if (individual_index >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::setIndividualMonolocusGenotypeByAlleleKey: individual_index out of bounds.", individual_index, 0, getNumberOfIndividuals());
	try {
		_individuals[individual_index]->setMonolocusGenotypeByAlleleKey(locus_index, allele_keys);
	}
	catch (NullPointerException & npe) {
		throw NullPointerException("Group::setIndividualMonolocusGenotypeByAlleleKey: individual has no genotype.");
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("Group::setIndividualMonolocusGenotypeByAlleleKey: locus_index excedes the number of locus.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
	catch (Exception) {
		throw Exception("Group::setIndividualMonolocusGenotypeByAlleleKey: no key in allele_keys.");
	}
}

void Group::setIndividualMonolocusGenotypeByAlleleId(unsigned int individual_index, unsigned int locus_index, const vector<string> allele_id, const LocusInfo & locus_info) throw (Exception) {
	if (individual_index >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::setIndividualMonolocusGenotypeByAlleleId: individual_index out of bounds.", individual_index, 0, getNumberOfIndividuals());
	try {
		_individuals[individual_index]->setMonolocusGenotypeByAlleleId(locus_index, allele_id, locus_info);
	}
	catch (NullPointerException & npe) {
		throw NullPointerException("Group::setIndividualMonolocusGenotypeByAlleleId: individual has no genotype.");
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("Group::setIndividualMonolocusGenotypeByAlleleId: locus_index excedes the number of locus.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
	catch (AlleleNotFoundException & anfe) {
		throw AlleleNotFoundException("Group::setIndividualMonolocusGenotypeByAlleleId: id not found.", anfe.getIdentifier());
	}
}

const MonolocusGenotype *  Group::getIndividualMonolocusGenotype(unsigned int individual_index, unsigned int locus_index) const throw (Exception) {
	if (individual_index >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::getIndividualMonolocusGenotype: individual_index out of bounds.", individual_index, 0, getNumberOfIndividuals());
	try {
		return _individuals[individual_index]->getMonolocusGenotype(locus_index);
	}
	catch (NullPointerException & npe) {
		throw NullPointerException("Group::getIndividualMonolocusGenotype: individual has no genotype.");
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("Group::getIndividualMonolocusGenotype: locus_index excedes the number of locus.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
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

bool Group::hasAllelicData() const {
	for (unsigned int i = 0 ; i < getNumberOfIndividuals() ; i++)
		if (hasIndividualGenotype(i)) return true;
	return false;
}
