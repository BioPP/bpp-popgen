/*
 * File DataSet.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Friday June 25 2004
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

bool DataSet::hasLocality() const {
	return (getNumberOfLocalities() > 0);
}

// Dealing with groups -------------------------------------
void DataSet::addGroup(const Group & group) {
	_groups.push_back(new Group(group));
}

void DataSet::addEmptyGroup() {
	_groups.push_back(new Group());
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

// Dealing with individuals -------------------------------
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

void DataSet::addEmptyIndividualToGroup(unsigned int group, const string & individual_id) throw (Exception) {
	if (group >= getNumberOfGroups())
		throw IndexOutOfBoundsException("DataSet::addEmptyIndividual: group out of bounds.", group, 0, getNumberOfGroups());
	try {
		_groups[group]->addEmptyIndividual(individual_id);
	}
	catch (BadIdentifierException & bie) {
		throw BadIdentifierException("DataSet::addEmptyIndividual: individual_id already in use.", bie.getIdentifier());
	}
}

unsigned int DataSet::getNumberOfIndividualsInGroup(unsigned int group_index) const throw (IndexOutOfBoundsException) {
	if (group_index >= getNumberOfGroups())
		throw IndexOutOfBoundsException("DataSet::getNumberOfIndividualsInGroup: group_index out of bounds.", group_index, 0, getNumberOfGroups());
	return _groups[group_index]->getNumberOfIndividuals();
}

unsigned int DataSet::getIndividualPositionInGroup(unsigned int group_index, const string & individual_id) const throw (Exception) {
	if (group_index >= getNumberOfGroups())
		throw IndexOutOfBoundsException("DataSet::getIndividualPositionFromGroup: group_index out of bounds.", group_index, 0, getNumberOfGroups());
	try {
		return _groups[group_index]->getIndividualPosition(individual_id);
	}
	catch (IndividualNotFoundException infe) {
		throw IndividualNotFoundException("DataSet::getIndividualPositionFromGroup: individual_id not found.", infe.getIdentifier());
	}
}

const Individual * DataSet::getIndividualByIndexFromGroup(unsigned int group_index, unsigned int individual_index) const throw (IndexOutOfBoundsException) {
	if (group_index >= getNumberOfGroups())
		throw IndexOutOfBoundsException("DataSet::getIndividualByIndexFromGroup: group_index out of bounds.", group_index, 0, getNumberOfGroups());
	try {
		return _groups[group_index]->getIndividualByIndex(individual_index);
	}
	catch (IndexOutOfBoundsException ioobe) {
		throw IndexOutOfBoundsException("DataSet::getIndividualByIndexFromGroup: individual_index out of bouds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
}

const Individual * DataSet::getIndividualByIdFromGroup(unsigned int group_index, const string & individual_id) const throw (Exception) {
	if (group_index >= getNumberOfGroups())
		throw IndexOutOfBoundsException("DataSet::getIndividualByIdFromGroup: group_index out of bounds.", group_index, 0, getNumberOfGroups());
	try {
		return _groups[group_index]->getIndividualById(individual_id);
	}
	catch (IndividualNotFoundException infe) {
		throw IndividualNotFoundException("DataSet::getIndividualByIdFromGroup: individual_id not found.", infe.getIdentifier());
	}
}

void DataSet::deleteIndividualByIndexFromGroup(unsigned int group_index, unsigned int individual_index) throw (IndexOutOfBoundsException) {
	if (group_index >= getNumberOfGroups())
		throw IndexOutOfBoundsException("DataSet::deleteIndividualByIndexFromGroup: group_index out of bounds.", group_index, 0, getNumberOfGroups());
	try {
		_groups[group_index]->deleteIndividualByIndex(individual_index);
	}
	catch (IndexOutOfBoundsException ioobe) {
		throw IndexOutOfBoundsException("DataSet::deleteIndividualByIndexFromGroup: individual_index out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
}

void DataSet::deleteIndividualByIdFromGroup(unsigned int group_index, const string & individual_id) throw (Exception) {
	if (group_index >= getNumberOfGroups())
		throw IndexOutOfBoundsException("DataSet::deleteIndividualByIdFromGroup: group_index out of bounds.", group_index, 0, getNumberOfGroups());
	try {
		_groups[group_index]->deleteIndividualById(individual_id);
	}
	catch (IndividualNotFoundException infe) {
		throw IndividualNotFoundException("DataSet::deleteIndividualByIdFromGroup: individual_id not found.", infe.getIdentifier());
	}
}

void DataSet::setIndividualSexInGroup(unsigned int group_index, unsigned int individual_index, const unsigned short sex) throw (IndexOutOfBoundsException) {
	if (group_index >= getNumberOfGroups())
		throw IndexOutOfBoundsException("DataSet::setIndividualSexInGroup: group_index out of bounds.", group_index, 0, getNumberOfGroups());
	try {
		_groups[group_index]->setIndividualSexByIndex(individual_index, sex);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("DataSet::setIndividualSexInGroup: individual_index out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
}

unsigned short DataSet::getIndividualSexInGroup(unsigned int group_index, unsigned int individual_index) const throw (IndexOutOfBoundsException) {
	if (group_index >= getNumberOfGroups())
		throw IndexOutOfBoundsException("DataSet::getIndividualSexInGroup: group_index out of bounds.", group_index, 0, getNumberOfGroups());
	try {
		return _groups[group_index]->getIndividualSexByIndex(individual_index);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("DataSet::getIndividualSexInGroup: individual_index out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
}

void DataSet::setIndividualDateInGroup(unsigned int group_index, unsigned int individual_index, const Date & date) throw (IndexOutOfBoundsException) {
	if (group_index >= getNumberOfGroups())
		throw IndexOutOfBoundsException("DataSet::setIndividualDateInGroup: group_index out of bounds.", group_index, 0, getNumberOfGroups());
	try {
		_groups[group_index]->setIndividualDateByIndex(individual_index, date);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("DataSet::setIndividualDateInGroup: individual_index out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
}

const Date * DataSet::getIndividualDateInGroup(unsigned int group_index, unsigned int individual_index) const throw (Exception) {
	if (group_index >= getNumberOfGroups())
		throw IndexOutOfBoundsException("DataSet::getIndividualDateInGroup: group_index out of bounds.", group_index, 0, getNumberOfGroups());
	try {
		return _groups[group_index]->getIndividualDateByIndex(individual_index);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("DataSet::getIndividualDateInGroup: individual_index out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
	catch (NullPointerException) {
		throw NullPointerException("DataSet::getIndividualDateInGroup: individual has no date.");
	}
}

void DataSet::setIndividualCoordInGroup(unsigned int group_index, unsigned int individual_index, const Coord<double> & coord) throw (IndexOutOfBoundsException) {
	if (group_index >= getNumberOfGroups())
		throw IndexOutOfBoundsException("DataSet::setIndividualCoordInGroup: group_index out of bounds.", group_index, 0, getNumberOfGroups());
	try {
		_groups[group_index]->setIndividualCoordByIndex(individual_index, coord);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("DataSet::setIndividualCoordInGroup: individual_index out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
}

const Coord<double> * DataSet::getIndividualCoordInGroup(unsigned int group_index, unsigned int individual_index) const throw (Exception) {
	if (group_index >= getNumberOfGroups())
		throw IndexOutOfBoundsException("DataSet::getIndividualCoordInGroup: group_index out of bounds.", group_index, 0, getNumberOfGroups());
	try {
		return _groups[group_index]->getIndividualCoordByIndex(individual_index);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("DataSet::getIndividualCoordByIndex: individual_index out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
	catch (NullPointerException) {
		throw NullPointerException("DataSet::getIndividualCoordInGroup: individual has no coordinate.");
	}
}

void DataSet::setIndividualLocalityInGroupByName(unsigned int group_index, unsigned int individual_index, const string & locality_name) throw (Exception) {
	if (group_index >= getNumberOfGroups())
		throw IndexOutOfBoundsException("DataSet::setIndividualLocalityInGroup: group_index out of bounds.", group_index, 0, getNumberOfGroups());
	try {
		_groups[group_index]->setIndividualLocalityByIndex(individual_index, getLocalityByName(locality_name));
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("DataSet::setIndividualLocalityInGroup: individual_index out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
	catch (LocalityNotFoundException & lnfe) {
		throw LocalityNotFoundException("DataSet::setIndividualLocalityInGroup: locality_name not found.", lnfe.getIdentifier());
	}
}

const Locality<double> * DataSet::getIndividualLocalityInGroup(unsigned int group_index, unsigned int individual_index) const throw (Exception) {
	if (group_index >= getNumberOfGroups())
		throw IndexOutOfBoundsException("DataSet::getIndividualLocalityInGroup: group_index out of bounds.", group_index, 0, getNumberOfGroups());
	try {
		return _groups[group_index]->getIndividualLocalityByIndex(individual_index);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("DataSet::getIndividualLocalityInGroup: individual_index out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
	catch (NullPointerException) {
		throw NullPointerException("DataSet::getIndividualLocalityInGroup: individual has no locality.");
	}
}

void DataSet::addIndividualSequenceInGroup(unsigned int group_index, unsigned int individual_index, unsigned int sequence_index, const Sequence & sequence) throw (Exception) {
	if (group_index >= getNumberOfGroups())
		throw IndexOutOfBoundsException("DataSet::addIndividualSequenceInGroup: group_index out of bounds.", group_index, 0, getNumberOfGroups());
	try {
		_groups[group_index]->addIndividualSequenceByIndex(individual_index, sequence_index, sequence);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("DataSet::addIndividualSequenceInGroup: individual_index out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
	catch (AlphabetMismatchException & ame) {
		throw AlphabetMismatchException("DataSet::addIndividualSequenceInGroup: sequence's alphabet doesn't match.", ame.getAlphabets()[0], ame.getAlphabets()[1]);
	}
	catch (BadIdentifierException & bie) {
		throw BadIdentifierException("DataSet::addIndividualSequenceInGroup: sequence's name already in use.", bie.getIdentifier());
	}
	catch (BadIntegerException & bie) {
		throw BadIntegerException("DataSet::addIndividualSequenceInGroup: sequence_index already in use.", bie.getBadInteger());
	}
}

const Sequence * DataSet::getIndividualSequenceByNameInGroup(unsigned int group_index, unsigned int individual_index, const string & sequence_name) const throw (Exception) {
	if (group_index >= getNumberOfGroups())
		throw IndexOutOfBoundsException("DataSet::getIndividualSequenceByNameInGroup: group_index out of bounds.", group_index, 0, getNumberOfGroups());
	try {
		return _groups[group_index]->getIndividualSequenceByName(individual_index, sequence_name);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("DataSet::getIndividualSequenceByNameInGroup: individual_index out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
	catch (NullPointerException) {
		throw NullPointerException("DataSet::getIndividualSequenceByNameInGroup: individual has no sequences.");
	}
	catch (SequenceNotFoundException & snfe) {
		throw SequenceNotFoundException("DataSet::getIndividualSequenceByNameInGroup: sequence_name not found.", snfe.getSequenceId());
	}
}

const Sequence * DataSet::getIndividualSequenceByIndexInGroup(unsigned int group_index, unsigned int individual_index, unsigned int sequence_index) const throw (Exception) {
	if (group_index >= getNumberOfGroups())
		throw IndexOutOfBoundsException("DataSet::getIndividualSequenceByIndexInGroup: group_index out of bounds.", group_index, 0, getNumberOfGroups());
	try {
		return _groups[group_index]->getIndividualSequenceByIndex(individual_index, sequence_index);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		if (string(ioobe.what()).find("individual_index") < string(ioobe.what()).size())
			throw IndexOutOfBoundsException("DataSet::getIndividualSequenceByIndexInGroup: individual_index out of bounds.", ioobe.getBadInteger(),ioobe.getBounds()[0], ioobe.getBounds()[1]);
		//if (string(ioobe.what()).find("sequence_index") < string(ioobe.what()).size())
		else
			throw IndexOutOfBoundsException("DataSet::getIndividualSequenceByIndexInGroup: sequence_index out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
	catch (NullPointerException) {
		throw NullPointerException("DataSet::getIndividualSequenceByIndexInGroup: individual has no sequences.");
	}
}

void DataSet::deleteIndividualSequenceByNameInGroup(unsigned int group_index, unsigned int individual_index, const string & sequence_name) throw (Exception) {
	if (group_index >= getNumberOfGroups())
		throw IndexOutOfBoundsException("DataSet::deleteIndividualSequenceByNameInGroup: group_index out of bounds.", group_index, 0, getNumberOfGroups());
	try {
		_groups[group_index]->deleteIndividualSequenceByName(individual_index, sequence_name);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("DataSet::deleteIndividualSequenceByNameInGroup: individual_index out of bounds.", ioobe.getBadInteger(),ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
	catch (NullPointerException) {
		throw NullPointerException("DataSet::deleteIndividualSequenceByNameInGroup: individual has no sequences.");
	}
	catch (SequenceNotFoundException & snfe) {
		throw SequenceNotFoundException("DataSet::deleteIndividualSequenceByNameInGroup: sequence_name not found.", snfe.getSequenceId());
	}
}

void DataSet::deleteIndividualSequenceByIndexInGroup(unsigned int group_index, unsigned int individual_index, unsigned int sequence_index) throw (Exception) {
	if (group_index >= getNumberOfGroups())
		throw IndexOutOfBoundsException("DataSet::deleteIndividualSequenceByIndexInGroup: group_index out of bounds.", group_index, 0, getNumberOfGroups());
	try {
		_groups[group_index]->deleteIndividualSequenceByIndex(individual_index, sequence_index);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		if (string(ioobe.what()).find("individual_index") < string(ioobe.what()).size())
			throw IndexOutOfBoundsException("DataSet::deleteIndividualSequenceByIndexInGroup: individual_index out of bounds.", ioobe.getBadInteger(),ioobe.getBounds()[0], ioobe.getBounds()[1]);
		//if (string(ioobe.what()).find("sequence_index") < string(ioobe.what()).size())
		else
			throw IndexOutOfBoundsException("DataSet::deleteIndividualSequenceByIndexInGroup: sequence_index out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
	catch (NullPointerException) {
		throw NullPointerException("DataSet::deleteIndividualSequenceByIndexInGroup: individual has no sequences.");
	}
}

vector<string> DataSet::getIndividualSequencesNamesInGroup(unsigned int group_index, unsigned int individual_index) const throw (Exception) {
	if (group_index >= getNumberOfGroups())
		throw IndexOutOfBoundsException("DataSet::getIndividualSequencesNamesInGroup: group_index out of bounds.", group_index, 0, getNumberOfGroups());
	try {
		return _groups[group_index]->getIndividualSequencesNames(individual_index);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("DataSet::getIndividualSequencesNamesInGroup: individual_index out of bounds.", ioobe.getBadInteger(),ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
	catch (NullPointerException) {
		throw NullPointerException("DataSet::getIndividualSequencesNamesInGroup: individual has no sequences.");
	}
}

unsigned int DataSet::getIndividualSequencePositionInGroup(unsigned int group_index, unsigned int individual_index, const string & sequence_name) const throw (Exception) {
	if (group_index >= getNumberOfGroups())
		throw IndexOutOfBoundsException("DataSet::getIndividualSequencePositionInGroup: group_index out of bounds.", group_index, 0, getNumberOfGroups());
	try {
		return _groups[group_index]->getIndividualSequencePosition(individual_index, sequence_name);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("DataSet::getIndividualSequencePositionInGroup: individual_index out of bounds.", ioobe.getBadInteger(),ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
	catch (NullPointerException) {
		throw NullPointerException("DataSet::getIndividualSequencePositionInGroup: individual has no sequences.");
	}
	catch (SequenceNotFoundException & snfe) {
		throw SequenceNotFoundException("DataSet::getIndividualSequencePositionInGroup: sequence_name not found.", snfe.getSequenceId());
	}
}

unsigned int DataSet::getIndividualNumberOfSequencesInGroup(unsigned int group_index, unsigned int individual_index) const throw (Exception) {
	if (group_index >= getNumberOfGroups())
		throw IndexOutOfBoundsException("DataSet::getIndividualNumberOfSequencesInGroup: group_index out of bounds.", group_index, 0, getNumberOfGroups());
	try {
		return _groups[group_index]->getIndividualNumberOfSequences(individual_index);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("DataSet::getIndividualNumberOfSequencesInGroup: individual_index out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
	catch (NullPointerException) {
		throw NullPointerException("DataSet::getIndividualNumberOfSequencesInGroup: individual has no sequences.");
	}
}

void DataSet::addIndividualGenotypeInGroup(unsigned int group_index, unsigned int individual_index, const Genotype & genotype) throw (Exception) {
	if (group_index >= getNumberOfGroups())
		throw IndexOutOfBoundsException("DataSet::getIndividualNumberOfSequencesInGroup: group_index out of bounds.", group_index, 0, getNumberOfGroups());
	try {
		_groups[group_index]->addIndividualGenotype(individual_index, genotype);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("DataSet::addIndividualGenotypeInGroup: individual_index out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
	catch (Exception) {
		throw Exception("DataSet::addIndividualGenotypeInGroup: individual already has a genotype.");
	}
}

void DataSet::initIndividualGenotypeInGroup(unsigned int group_index, unsigned int individual_index, const AnalyzedLoci * analyzed_loci) throw (Exception) {
	if (group_index >= getNumberOfGroups())
		throw IndexOutOfBoundsException("DataSet::initIndividualGenotypeInGroup: group_index out of bounds.", group_index, 0, getNumberOfGroups());
	try {
		_groups[group_index]->initIndividualGenotype(individual_index, analyzed_loci);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("DataSet::initIndividualGenotypeInGroup: individual_index out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
	catch (NullPointerException) {
		throw NullPointerException("DataSet::initIndividualGenotypeInGroup: analyzed_loci is NULL.");
	}
	catch (Exception) {
		throw Exception("DataSet::initIndividualGenotypeInGroup: individual already has a genotype.");
	}
}

void DataSet::deleteIndividualGenotypeInGroup(unsigned int group_index, unsigned int individual_index) throw (IndexOutOfBoundsException) {
	if (group_index >= getNumberOfGroups())
		throw IndexOutOfBoundsException("DataSet::deleteIndividualGenotypeInGroup: group_index out of bounds.", group_index, 0, getNumberOfGroups());
	try {
		_groups[group_index]->deleteIndividualGenotype(individual_index);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("DataSet::deleteIndividualGenotypeInGroup: individual_index out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
}

void DataSet::setIndividualMonolocusGenotypeInGroup(unsigned int group_index, unsigned int individual_index, unsigned int locus_index, const MonolocusGenotype & monogen) throw (Exception) {
	if (group_index >= getNumberOfGroups())
		throw IndexOutOfBoundsException("DataSet::setIndividualMonolocusGenotypeInGroup: group_index out of bounds.", group_index, 0, getNumberOfGroups());
	try {
		_groups[group_index]->setIndividualMonolocusGenotype(individual_index, locus_index, monogen);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		if (string(ioobe.what()).find("individual_index") < string(ioobe.what()).size())
			throw IndexOutOfBoundsException("DataSet::setIndividualMonolocusGenotypeInGroup: individual_index out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
		//if (string(ioobe.what()).find("locus_index") < string(ioobe.what()).size())
		else
			throw IndexOutOfBoundsException("DataSet::setIndividualMonolocusGenotypeInGroup: locus_index out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
	catch (NullPointerException) {
		throw NullPointerException("DataSet::setIndividualMonolocusGenotypeInGroup: individual has no genotype.");
	}
}

void DataSet::setIndividualMonolocusGenotypeByAlleleKeyInGroup(unsigned int group_index, unsigned int individual_index, unsigned int locus_index, const vector<unsigned int> allele_keys) throw (Exception) {
	if (group_index >= getNumberOfGroups())
		throw IndexOutOfBoundsException("DataSet::setIndividualMonolocusGenotypeByAlleleKeyInGroup: group_index out of bounds.", group_index, 0, getNumberOfGroups());
	try {
		_groups[group_index]->setIndividualMonolocusGenotypeByAlleleKey(individual_index, locus_index, allele_keys);
	}
  catch (IndexOutOfBoundsException & ioobe) {
		if (string(ioobe.what()).find("individual_index") < string(ioobe.what()).size())
			throw IndexOutOfBoundsException("DataSet::setIndividualMonolocusGenotypeByAlleleKeyInGroup: individual_index out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
		//if (string(ioobe.what()).find("locus_index") < string(ioobe.what()).size())
		else    
			throw IndexOutOfBoundsException("DataSet::setIndividualMonolocusGenotypeByAlleleKeyInGroup: locus_index out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
	catch (NullPointerException) {
		throw NullPointerException("DataSet::setIndividualMonolocusGenotypeByAlleleKeyInGroup: individual has no genotype.");
	}
	catch (Exception) {
		throw Exception("DataSet::setIndividualMonolocusGenotypeByAlleleKeyInGroup: allele_keys.size() doesn't match ploidy.");
	}
}

void DataSet::setIndividualMonolocusGenotypeByAlleleIdInGroup(unsigned int group_index, unsigned int individual_index, unsigned int locus_index, const vector<unsigned int> allele_id) throw (Exception) {
	if (group_index >= getNumberOfGroups())
		throw IndexOutOfBoundsException("DataSet::setIndividualMonolocusGenotypeByAlleleIdInGroup: group_index out of bounds.", group_index, 0, getNumberOfGroups());
	try {
		_groups[group_index]->setIndividualMonolocusGenotypeByAlleleId(individual_index, locus_index, allele_id);
	}
  catch (IndexOutOfBoundsException & ioobe) {
		if (string(ioobe.what()).find("individual_index") < string(ioobe.what()).size())
			throw IndexOutOfBoundsException("DataSet::setIndividualMonolocusGenotypeByAlleleIdInGroup: individual_index out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
		//if (string(ioobe.what()).find("locus_index") < string(ioobe.what()).size())
		else    
			throw IndexOutOfBoundsException("DataSet::setIndividualMonolocusGenotypeByAlleleIdInGroup: locus_index out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
	catch (NullPointerException) {
		throw NullPointerException("DataSet::setIndividualMonolocusGenotypeByAlleleIdInGroup: individual has no genotype.");
	}
	catch (Exception) {
		throw Exception("DataSet::setIndividualMonolocusGenotypeByAlleleIdInGroup: allele_id.size() doesn't match ploidy.");
	}
}

const MonolocusGenotype * DataSet::getIndividualMonolocusGenotypeInGroup(unsigned int group_index, unsigned int individual_index, unsigned int locus_index) const throw (Exception) {
	if (group_index >= getNumberOfGroups())
		throw IndexOutOfBoundsException("DataSet::getIndividualMonolocusGenotypeInGroup: group_index out of bounds.", group_index, 0, getNumberOfGroups());
	try {
		return _groups[group_index]->getIndividualMonolocusGenotype(individual_index, locus_index);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		if (string(ioobe.what()).find("individual_index") < string(ioobe.what()).size())
			throw IndexOutOfBoundsException("DataSet::getIndividualMonolocusGenotypeInGroup: individual_index out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
		//if (string(ioobe.what()).find("locus_index") < string(ioobe.what()).size())
		else
			throw IndexOutOfBoundsException("DataSet::getIndividualMonolocusGenotypeInGroup: locus_index out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
	catch (NullPointerException) {
		throw NullPointerException("DataSet::getIndividualMonolocusGenotypeInGroup: individual has no genotype.");
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

// General tests ------------------------------------------
bool DataSet::hasSequenceData() const {
	for (unsigned int i = 0 ; i < getNumberOfGroups() ; i++)
		if (_groups[i]->hasSequenceData()) return true;
	return false;
}

const Alphabet * DataSet::getAlphabet() const throw (NullPointerException) {
	for (unsigned int i = 0 ; i < getNumberOfGroups() ; i++)
		if (_groups[i]->hasSequenceData())
			return _groups[i]->getAlphabet();
	throw NullPointerException("DataSet::getAlphabet: no sequence data.");
}

bool DataSet::hasAlleleicData() const {
	for (unsigned int i = 0 ; i < getNumberOfGroups() ; i++)
		if (_groups[i]->hasAllelicData()) return true;
	return false;
}

unsigned int DataSet::getNumberOfSequenceSets() const {
	unsigned int nbss = 0;
	unsigned int tmp;
	for (unsigned int i = 0 ; i < getNumberOfGroups() ; i++) {
		tmp = _groups[i]->getMaxNumberOfSequences();
		if (nbss < tmp) nbss = tmp;
	}
	return nbss;
}
