/*
 * File Group.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Monday June 21 2004
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
	catch (NullPointerException npe) {
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
	catch (NullPointerException npe) {
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
	catch (NullPointerException npe) {
		throw NullPointerException("Group::getIndividualLocalityByIndex: individuals has no locality.");
	}
}

vector<string> Group::getIndividualSequencesKeys(unsigned int individual_index) const throw (IndexOutOfBoundsException) {
	if (individual_index >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::getIndividualSequencesKeys: individual_index out of bounds.", individual_index, 0, getNumberOfIndividuals());
	return _individuals[individual_index]->getSequencesKeys();
}

void Group::addSequenceToIndividualByIndex(unsigned int individual_index, const string & seq_set, const Sequence & sequence) throw (Exception) {
	if (individual_index >= getNumberOfIndividuals())
		throw IndexOutOfBoundsException("Group::addSequenceToIndividualByIndex: individual_index out of bounds.", individual_index, 0, getNumberOfIndividuals());
	try {
		_individuals[individual_index]->addSequence(seq_set, sequence);
	}
	catch (AlphabetMismatchException ame) {
		throw AlphabetMismatchException("Group::addSequenceToIndividualByIndex: sequence's alphabet doesn't match.", ame.getAlphabets()[0], ame.getAlphabets()[1]);
	}
	catch (BadIdentifierException bie) {
		throw BadIdentifierException("Group::addSequenceToIndividualByIndex: sequence's name already in use.", bie.getIdentifier());
	}
}

//-- Dealing with sequence containers ---------------------
VectorSequenceContainer * Group::getVectorSequenceContainer(const string & seqset_id) const  throw (Exception) {
	if (_individuals.size() == 0)
		throw Exception("Group::getVectorSequenceContainer: this group is empty.");
	VectorSequenceContainer * vsc = new VectorSequenceContainer(_individuals[0]->getVectorSequenceContainer(seqset_id)->getAlphabet());
	for (unsigned int i = 0 ; i < _individuals.size() ; i++) {
		try {
			SequenceContainerTools::append<VectorSequenceContainer>(*vsc, *(_individuals[i]->getVectorSequenceContainer(seqset_id)));
		}
		catch (Exception & e) {
			throw e;
		}
	}
}

VectorSiteContainer * Group::getVectorSiteContainer(const string & seqset_id) const  throw (Exception) {
	if (_individuals.size() == 0)
		throw Exception("Group::getVectorSequenceContainer: this group is empty.");
	VectorSiteContainer * vsc = new VectorSiteContainer(_individuals[0]->getVectorSequenceContainer(seqset_id)->getAlphabet());
	for (unsigned int i = 0 ; i < _individuals.size() ; i++) {
		try {
			SequenceContainerTools::append<VectorSiteContainer>(*vsc, *(_individuals[i]->getVectorSequenceContainer(seqset_id)));
		}
		catch (Exception & e) {
			throw e;
		}
	}
}
