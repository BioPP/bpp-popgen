/*
 * File Group.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Monday June 14 2004
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

void Group::addIndividual(const Individual & ind) {
	_individuals.push_back(new Individual(ind));
}

Individual * Group::removeIndividualById(const string individual_id) throw (BadIdentifierException) {
	for (vector<Individual *>::iterator it = _individuals.begin() ; it != _individuals.end() ; it++) {
		if ((*it)->getId() == individual_id) {
			Individual * ind = *it;
			_individuals.erase(it);
			return ind;
		}
	}
	throw BadIdentifierException("Group::removeIndividualById: individual_id not found.", individual_id);
}

Individual * Group::removeIndividualByIndex(unsigned int index) throw (IndexOutOfBoundsException) {
	if (index >= _individuals.size())
		throw IndexOutOfBoundsException("Group::removeIndividualByIndex: indes out of bounds.", index, 0, _individuals.size());
	vector<Individual *>::iterator it = _individuals.begin()+index;
	Individual * ind = *it;
	_individuals.erase(it);
	return ind;
}

void Group::deleteIndividualById(const string individual_id) throw (BadIdentifierException) {
	try {
		Individual * ind = removeIndividualById(individual_id);
		delete ind;
	}
	catch (BadIdentifierException bie) {
		throw bie;
	}
}

void Group::deleteIndividualByIndex(unsigned int index) throw (IndexOutOfBoundsException) {
	try {
		Individual * ind = removeIndividualByIndex(index);
		delete ind;
	}
	catch (IndexOutOfBoundsException ioobe) {
		throw ioobe;
	}
}

void Group::clear() {
	for (unsigned int i = 0 ; i < _individuals.size() ; i++)
		delete(_individuals[i]);
	_individuals.clear();
}

const Individual * Group::getIndividualById(const string individual_id) const
throw (BadIdentifierException) {
	for (unsigned int i = 0 ; i < _individuals.size() ; i++) {
		if (_individuals[i]->getId() == individual_id)
			return getIndividualByIndex(i);
	}
	throw BadIdentifierException("Group::getIndividualById: individual_id not found.", individual_id);
}

const Individual * Group::getIndividualByIndex(unsigned int index) const
throw (IndexOutOfBoundsException) {
	if (index >= _individuals.size())
		throw IndexOutOfBoundsException("Group::getIndividualByIndex: index out of bounds.", index, 0, _individuals.size());
	return _individuals[index];
}

int Group::getNumberOfIndividuals() {
	return _individuals.size();
}

VectorSequenceContainer * Group::getVectorSequenceContainer(const string & seqset_id) const  throw (Exception) {
	if (_individuals.size() == 0)
		throw Exception("Group::getVectorSequenceContainer: this group is empty.");
	VectorSequenceContainer * vsc = new VectorSequenceContainer(_individuals[0]->getVectorSequenceContainer(seqset_id)->getAlphabet());
	for (unsigned int i = 0 ; i < _individuals.size() ; i++) {
		try {
			SequenceContainerTools::append<VectorSequenceContainer>(*vsc, *(_individuals[i]->getVectorSequenceContainer(seqset_id)));
		}
		catch (Exception e) {
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
		catch (Exception e) {
			throw e;
		}
	}
}
