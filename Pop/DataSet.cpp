/*
 * File DataSet.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Monday June 21 2004
 */

#include "DataSet.h"

//** Class constructor: *******************************************************/
DataSet::DataSet() {
	_analyzedLoci = NULL;
}

//** Class destructor: *******************************************************/
DataSet::~DataSet() {
	if (getNumberOfGroups() > 0)
		for (unsigned int i = 0 ; i < getNumberOfGroups() ; i++)
			delete _groups[i];
	if (_analyzedLoci != NULL) delete _analyzedLoci;
	if (getNumberOfLocalities() > 0)
		for (unsigned int i = 0 ; i < getNumberOfLocalities() ; i++)
			delete _localities[i];
}

//** Other methodes: *********************************************************/
// Dealing with file ---------------------------------------
void DataSet::readFile(const string & path) {}

void DataSet::writeFile(const string & path) {}

// Dealing with Localities ---------------------------------
void DataSet::addLocality(Locality<double> & locality) throw (BadIdentifierException) {
	for (unsigned int i = 0 ; i < _localities.size() ; i++)
		if (_localities[i]->getName() == locality.getName())
			throw BadIdentifierException("DataSet::addLocality: locality name already in use.", locality.getName());
	_localities.push_back(new Locality<double>(locality));
}

unsigned int DataSet::getLocalityPosition(const string & name) const throw (LocalityNotFoundException) {
	for (unsigned int i = 0 ; i < _localities.size() ; i++)
		if (_localities[i]->getName() == name)
			return i;
	throw LocalityNotFoundException("DataSet::getLocalityPosition: Locality not found.", name);
}

const Locality<double> * DataSet::getLocalityByIndex(unsigned int index) const throw (IndexOutOfBoundsException) {
	if (index >= _localities.size()) throw IndexOutOfBoundsException("DataSet::getLocalityByIndex: index out of bounds.", index, 0, _localities.size());
	return _localities[index];
}

const Locality<double> * DataSet::getLocalityByName(const string & name) const throw (LocalityNotFoundException) {
	try {
		return (getLocalityByIndex(getLocalityPosition(name)));
	}
	catch (LocalityNotFoundException & lnfe) {
		throw LocalityNotFoundException("DataSet::getLocalityByName: Locality not found.", name);
	}
}

void DataSet::deleteLocalityByIndex(unsigned int index) throw (IndexOutOfBoundsException) {
	if (index >= _localities.size()) throw IndexOutOfBoundsException("DataSet::deleteLocalityByIndex: index out of bounds.", index, 0, _localities.size());
	delete _localities[index];
	_localities.erase(_localities.begin() + index);
}	

void DataSet::deleteLocalityByName(const string & name) throw (LocalityNotFoundException) {
	try {
		deleteLocalityByIndex(getLocalityPosition(name));
	}
	catch (LocalityNotFoundException & lnfe) {
		throw LocalityNotFoundException("DataSet::deleteLocalityByName: Locality not found.", name);
	}
}

unsigned int DataSet::getNumberOfLocalities() const {
	return _localities.size();
}

// Dealing with groups -------------------------------------
void DataSet::addGroup(const Group & group) {
	_groups.push_back(new Group(group));
}

const Group * DataSet::getGroup(unsigned int index) const throw (IndexOutOfBoundsException) {
	if (index >= _groups.size())
		throw IndexOutOfBoundsException("DataSet::getGroup.", index, 0, _groups.size());
	return _groups[index];
}

void DataSet::deleteGroup(unsigned int index) throw (IndexOutOfBoundsException) {
	if (index >= _groups.size())
		throw IndexOutOfBoundsException("DataSet::deleteGroup.", index, 0, _groups.size());
	delete _groups[index];
	_groups.erase(_groups.begin() + index);
}

unsigned int DataSet::getNumberOfGroups() const {
	return _groups.size();
}

void DataSet::mergeGroups(vector<unsigned int> & groups) throw (IndexOutOfBoundsException) {
	for (unsigned int i = 0 ; i < groups.size() ; i++)
		if (groups[i] >= getNumberOfGroups())
			throw IndexOutOfBoundsException("DataSet::mergeGroups.", groups[i], 0, getNumberOfGroups());
	sort(groups.begin(), groups.end());
	for (unsigned int i = 1 ; i < groups.size() ; i++) {
		_groups[groups[0]]->append(*getGroup(groups[i]));
		deleteGroup(groups[i]-i+1);
	}
}

void DataSet::addIndividualToGroup(unsigned int group, const Individual & individual) throw (Exception) {
	if (group >= getNumberOfGroups())
		throw IndexOutOfBoundsException("DataSet::addIndividualToGroup: group out of bounds.", group, 0, getNumberOfGroups());
	try {
		_groups[group]->addIndividual(individual);
	}
	catch (BadIdentifierException & bie) {
		throw BadIdentifierException("DataSet::addIndividualToGroup: individual's id already in use in this group.", bie.getIdentifier());
	}
}

unsigned int DataSet::getIndividualPositionFromGroup(const string & individual_id, unsigned int group) const throw (Exception) {
	if (group >= getNumberOfGroups())
		throw IndexOutOfBoundsException("DataSet::getIndividualPositionFromGroup: group out of bounds.", group, 0, getNumberOfGroups());
	try {
		return _groups[group]->getIndividualPosition(individual_id);
	}
	catch (IndividualNotFoundException infe) {
		throw IndividualNotFoundException("DataSet::getIndividualPositionFromGroup: individual_id not found.", infe.getIdentifier());
	}
}

const Individual * DataSet::getIndividualByIndexFromGroup(unsigned int individual_index, unsigned int group) const throw (IndexOutOfBoundsException) {
	if (group >= getNumberOfGroups())
		throw IndexOutOfBoundsException("DataSet::getIndividualByIndexFromGroup: group out of bounds.", group, 0, getNumberOfGroups());
	try {
		return _groups[group]->getIndividualByIndex(individual_index);
	}
	catch (IndexOutOfBoundsException ioobe) {
		throw IndexOutOfBoundsException("DataSet::getIndividualByIndexFromGroup: individual_index out of bouds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
}

const Individual * DataSet::getIndividualByIdFromGroup(const string & individual_id, unsigned int group) const throw (Exception) {
	if (group >= getNumberOfGroups())
		throw IndexOutOfBoundsException("DataSet::getIndividualByIdFromGroup: group out of bounds.", group, 0, getNumberOfGroups());
	try {
		return _groups[group]->getIndividualById(individual_id);
	}
	catch (IndividualNotFoundException infe) {
		throw IndividualNotFoundException("DataSet::getIndividualByIdFromGroup: individual_id not found.", infe.getIdentifier());
	}
}

void DataSet::deleteIndividualByIndexFromGroup(unsigned int individual_index, unsigned int group) throw (IndexOutOfBoundsException) {
	if (group >= getNumberOfGroups())
		throw IndexOutOfBoundsException("DataSet::deleteIndividualByIndexFromGroup: group out of bounds.", group, 0, getNumberOfGroups());
	try {
		_groups[group]->deleteIndividualByIndex(individual_index);
	}
	catch (IndexOutOfBoundsException ioobe) {
		throw IndexOutOfBoundsException("DataSet::deleteIndividualByIndexFromGroup: individual_index out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
}

void DataSet::deleteIndividualByIdFromGroup(const string & individual_id, unsigned int group) throw (Exception) {
	if (group >= getNumberOfGroups())
		throw IndexOutOfBoundsException("DataSet::deleteIndividualByIdFromGroup: group out of bounds.", group, 0, getNumberOfGroups());
	try {
		_groups[group]->deleteIndividualById(individual_id);
	}
	catch (IndividualNotFoundException infe) {
		throw IndividualNotFoundException("DataSet::deleteIndividualByIdFromGroup: individual_id not found.", infe.getIdentifier());
	}
}

// Dealing with AnalyzedLoci -------------------------------
void DataSet::setAnalyzedLoci(AnalyzedLoci & analyzeLoci) throw (Exception) {
	if (_analyzedLoci != NULL) {
		try {
			deleteAnalyzedLoci();
		}
		catch (Exception & e) {
			throw Exception ("DataSet::setAnalyzedLoci: at least one individual has a genotype of the actual AnalyzedLoci.");
		}
	_analyzedLoci = new AnalyzedLoci(analyzeLoci);
	}
}

void DataSet::initAnalyzedLoci(unsigned int number_of_loci) throw (Exception) {
	if (_analyzedLoci != NULL)
		throw Exception("DataSet::initAnalyzedLoci: _analyzedLoci already initialyzed.");
	_analyzedLoci = new AnalyzedLoci(number_of_loci);
}

void DataSet::deleteAnalyzedLoci() throw (Exception) {
	if (_analyzedLoci != NULL) {
		for (unsigned int i = 0 ; i < getNumberOfGroups() ; i++)
			for (unsigned int j = 0 ; j < getGroup(i)->getNumberOfIndividuals() ; j++)
				if (getGroup(i)->getIndividualByIndex(j)->hasGenotype())
					throw Exception("DataSet::deleteAnalyzedLoci: at least one individual has a genotype.");
		delete _analyzedLoci;
	}
}

void DataSet::setLocusInfo(unsigned int locus_index, const LocusInfo & locus) throw (Exception) {
	if (_analyzedLoci == NULL)
		throw NullPointerException("DataSet::setLocusInfo: there's no AnalyzedLoci to setup.");
	for (unsigned int i = 0 ; i < getNumberOfGroups() ; i++)
		for (unsigned int j = 0 ; j < getGroup(i)->getNumberOfIndividuals() ; j++)
			if (getGroup(i)->getIndividualByIndex(j)->hasGenotype())
				throw Exception("DataSet::setLocusInfo: at least one individual has a genotype.");
	try {
		_analyzedLoci->setLocusInfo(locus_index, locus);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("DataSet::setLocusInfo: locus_index out of bounds.", locus_index, 0, _analyzedLoci->getNumberOfLoci());
	}
}

const LocusInfo * DataSet::getLocusInfoByName(const string & locus_name) const throw (Exception) {
	if (_analyzedLoci == NULL)
		throw NullPointerException("DataSet::getLocusInfoByName: there's no AnalyzedLoci.");
	try {
		return _analyzedLoci->getLocusInfoByName(locus_name);
	}
	catch (LocusNotFoundException & lnfe) {
		throw LocusNotFoundException("DataSet::getLocusInfoByName: locus_name not found", locus_name);
	}
}

const LocusInfo * DataSet::getLocusInfoByIndex(unsigned int locus_index) const throw (Exception) {
	if (_analyzedLoci == NULL)
		throw NullPointerException("DataSet::getLocusInfoByIndex: there's no AnalyzedLoci.");
	try {
		return _analyzedLoci->getLocusInfoByIndex(locus_index);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("DataSet::getLocusInfoByIndex: locus_index out of bounds.", locus_index, ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
	catch (NullPointerException & npe) {
		throw NullPointerException("DataSet::getLocusInfoByIndex: no locus defined here");
	}
}

void DataSet::addAlleleInfoByLocusName(const string & locus_name, const AlleleInfo & allele) throw (Exception) {
	if (_analyzedLoci == NULL)
		throw NullPointerException("DataSet::addAlleleInfoByLocusName: there's no AnalyzedLoci.");
	try {
		_analyzedLoci->addAlleleInfoByLocusName(locus_name, allele);
	}
	catch (LocusNotFoundException & lnfe) {
		throw LocusNotFoundException("DataSet::addAlleleInfoByLocusName: locus_name not found.", lnfe.getIdentifier());
	}
	catch (BadIdentifierException & bie) {
		throw BadIdentifierException("DataSet::addAlleleInfoByLocusName: allele's id already in use.", bie.getIdentifier());
	}
}

void DataSet::addAlleleInfoByLocusIndex(unsigned int locus_index, const AlleleInfo & allele) throw (Exception) {
	if (_analyzedLoci == NULL)
		throw NullPointerException("DataSet::addAlleleInfoByLocusIndex: there's no AnalyzedLoci.");
	try {
		_analyzedLoci->addAlleleInfoByLocusIndex(locus_index, allele);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("DataSet::addAlleleInfoByLocusIndex: locus_index out of bounds.", locus_index, ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
	catch (BadIdentifierException & bie) {
		throw BadIdentifierException("DataSet::addAlleleInfoByLocusIndex: allele'e id already in use.", bie.getIdentifier());
	}
}

unsigned int DataSet::getNumberOfLoci() const throw (NullPointerException) {
	if (_analyzedLoci == NULL)
		throw NullPointerException("DataSet::getNumberOfLoci: there's no AnalyzedLoci.");
	return _analyzedLoci->getNumberOfLoci();
}

unsigned int DataSet::getPloidyByLocusName(const string & locus_name) const throw (Exception) {
	if (_analyzedLoci == NULL)
		throw NullPointerException("DataSet::getPloidyByLocusName: there's no AnalyzedLoci.");
	try {
		return _analyzedLoci->getPloidyByLocusName(locus_name);
	}
	catch (LocusNotFoundException & lnfe) {
		throw LocusNotFoundException("DataSet::getPloidyByLocusName: locus_name not found.", lnfe.getIdentifier());
	}
}

unsigned int DataSet::getPloidyByLocusIndex(unsigned int locus_index) const throw (Exception) {
	if (_analyzedLoci == NULL)
		throw NullPointerException("DataSet::getPloidyByLocusIndex: there's no AnalyzedLoci.");
	try {
		return _analyzedLoci->getPloidyByLocusIndex(locus_index);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("DataSet::getPloidyByLocusIndex: locus_index out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
}
