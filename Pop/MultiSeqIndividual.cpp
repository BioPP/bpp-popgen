/*
 * File MultiSeqIndividual.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Tuesday August 03 2004
 *
*/
/*
Copyright or © or Copr. CNRS, (November 17, 2004)


This software is a computer program whose purpose is to provide classes
for sequences analysis.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/
#include "MultiSeqIndividual.h"

//** Class constructor: *******************************************************/
MultiSeqIndividual::MultiSeqIndividual() {
	_id = "";
	_sex = 0;
	_date = NULL;
	_coord = NULL;
	_locality = NULL;
	_genotype = NULL;
}

MultiSeqIndividual::MultiSeqIndividual(const string & id) {
	_id = id;
	_sex = 0;
	_date = NULL;
	_coord = NULL;
	_locality = NULL;
	_genotype = NULL;
}

MultiSeqIndividual::MultiSeqIndividual(const string & id,
                       const Date & date,
                       const Coord<double> & coord,
                       Locality<double> * locality,
                       const unsigned short sex) {
	_id = id;
	_sex = sex;
	_date = new Date(date);
	_coord = new Coord<double>(coord);
	_locality = locality;
}

MultiSeqIndividual::MultiSeqIndividual(const MultiSeqIndividual &ind) {
	setId(ind.getId());
	setSex(ind.getSex());
	try {
		setDate(* ind.getDate());
	}
	catch (NullPointerException) {
		_date = NULL;
	}
	try {
		setCoord(* ind.getCoord());
	}
	catch (NullPointerException) {
		_coord = NULL;
	}
	try {
		setLocality(ind.getLocality());
	}
	catch (NullPointerException) {
		_locality = NULL;
	}
	if (ind.hasSequences()) {
		vector<string> keys = ind.getSequencesKeys();
		for (unsigned int i = 0 ; i < keys.size() ; i++)
			_sequences[keys[i]] = new VectorSequenceContainer(* const_cast<const VectorSequenceContainer *>(ind.getVectorSequenceContainer(keys[i])));
	}
	this->_genotype = ind.hasGenotype() ? new MultilocusGenotype(* ind.getGenotype()) : NULL;
}

//** Class destructor: *******************************************************/
MultiSeqIndividual::~MultiSeqIndividual () {
	delete this->_date;
	delete this->_coord;
}

//** Other methodes: *********************************************************/
MultiSeqIndividual & MultiSeqIndividual::operator= (const MultiSeqIndividual & ind) {
	setId(ind.getId());
	setSex(ind.getSex());
	try {
		setDate(* ind.getDate());
	}
	catch (NullPointerException) {
		_date = NULL;
	}
	try {
		setCoord(* ind.getCoord());
	}
	catch (NullPointerException) {
		_coord = NULL;
	}
	try {
		setLocality(ind.getLocality());
	}
	catch (NullPointerException) {
		_locality = NULL;
	}
	if (ind.hasSequences()) {
		vector<string> keys = ind.getSequencesKeys();
		for (unsigned int i = 0 ; i < keys.size() ; i++)
			_sequences[keys[i]] = new VectorSequenceContainer(* const_cast<const VectorSequenceContainer *>(ind.getVectorSequenceContainer(keys[i])));
	}
	this->_genotype = ind.hasGenotype() ? new MultilocusGenotype(* ind.getGenotype()) : NULL;
	return * this;
}

// Id
void MultiSeqIndividual::setId(const string id) {
	_id = id;
}

string MultiSeqIndividual::getId() const {
	return _id;
}

// Sex
void MultiSeqIndividual::setSex(const unsigned short sex) {
	_sex = sex;
}

unsigned short MultiSeqIndividual::getSex() const {
	return _sex;
}

// Date
void MultiSeqIndividual::setDate(const Date & date) {
	if (!hasDate()) {
		_date = new Date(date);
	}
	else if (* _date != date) {
		delete _date;
		_date = new Date(date);
	}
}

const Date * MultiSeqIndividual::getDate() const throw (NullPointerException) {
	if (hasDate())
		return new Date(* _date);
	else
		throw(NullPointerException("MultiSeqIndividual::getDate: no date associated to this individual."));
}

bool MultiSeqIndividual::hasDate() const {
	return _date != NULL;
}

// Coord
void MultiSeqIndividual::setCoord(const Coord<double> & coord) {
	if (!hasCoord()) {
		_coord = new Coord<double>(coord);
	}
	else if	(* _coord != coord) {
		delete _coord;
		_coord = new Coord<double>(coord);
	}
}

void MultiSeqIndividual::setCoord(const double x, const double y) {
	if (!hasCoord()) {
		_coord = new Coord<double>(x, y);
	}
	else if (this->getX() != x || this->getY() != y) {
		delete _coord;
		_coord = new Coord<double>(x, y);
	}
}

const Coord<double> * MultiSeqIndividual::getCoord() const throw(NullPointerException) {
	if (hasCoord())
		return new Coord<double>(* _coord);
	else
		throw(NullPointerException("MultiSeqIndividual::getCoord: no coord associated to this individual."));
}

bool MultiSeqIndividual::hasCoord() const {
	return _coord != NULL;
}

void MultiSeqIndividual::setX(const double x) throw(NullPointerException) {
	if (hasCoord())
		_coord->setX(x);
	else
		throw(NullPointerException("MultiSeqIndividual::setX: no coord associated to this individual."));
}

void MultiSeqIndividual::setY(const double y) throw(NullPointerException) {
	if (hasCoord())
		_coord->setY(y);
	else
		throw(NullPointerException("MultiSeqIndividual::setY: no coord associated to this individual."));
}

double MultiSeqIndividual::getX() const throw(NullPointerException) {
	if (hasCoord())
		return _coord->getX();
	else
		throw(NullPointerException("MultiSeqIndividual::getX: no coord associated to this individual."));
}

double MultiSeqIndividual::getY() const throw(NullPointerException) {
	if (hasCoord())
		return _coord->getY();
	else
		throw(NullPointerException("MultiSeqIndividual::getY: no coord associated to this individual."));
}

// Locality
void MultiSeqIndividual::setLocality(const Locality<double> * locality) {
	_locality = locality;
}

const Locality<double> * MultiSeqIndividual::getLocality() const  throw (NullPointerException) {
	if (hasLocality())
		return _locality;
	else
		throw(NullPointerException("MultiSeqIndividual::getLocality: no locality associated to this individual."));
}

bool MultiSeqIndividual::hasLocality() const {
	return _locality != NULL;
}

// Sequences
const VectorSequenceContainer * MultiSeqIndividual::getVectorSequenceContainer(const string & id) const throw (Exception) {
	map<string, VectorSequenceContainer *>::const_iterator it;
	it = _sequences.find(id);
	// Test existence of id in the map.
	if (it == _sequences.end()) {
		string mes = "MultiSeqIndividual::getSequence: sequence set not found (" + id
			+ ").";
		throw(Exception(mes));
	}
	return const_cast<const VectorSequenceContainer *>(it->second);
}

void MultiSeqIndividual::addSequence(const string & id, const Sequence & sequence)
throw (Exception) {
	try {
		_sequences[id]->addSequence(sequence);
	}
	catch (AlphabetMismatchException & ame)
	{
		throw(AlphabetMismatchException("MultiSeqIndividual::addSequence: alphabets don't match.", ame.getAlphabets()[0], ame.getAlphabets()[1]));
	}
	catch (Exception & e)
	{
		throw(BadIdentifierException("MultiSeqIndividual::addSequence: sequence's name already in use.", sequence.getName()));
	}
}

const Sequence * MultiSeqIndividual::getSequence(const string & id, const string & name)
const throw(Exception){
	map<string, VectorSequenceContainer *>::const_iterator it;
	it = _sequences.find(id);
	// Test existence of id in the map.
	if (it == _sequences.end()) {
		string mes = "MultiSeqIndividual::getSequence: sequence set not found (" + id
			+ ").";
		throw(Exception(mes));
	}
	try {
		return const_cast<const VectorSequenceContainer *>(it->second)->getSequence(name);
	}
	catch (SequenceNotFoundException & snfe) {
		throw(snfe);
	}
}

const Sequence * MultiSeqIndividual::getSequence(const string & id, unsigned int i)
const throw(Exception) {
	map<string, VectorSequenceContainer *>::const_iterator it;
	it = _sequences.find(id);
	// Test existence of id in the map.
	if (it == _sequences.end()) {
		string mes = "MultiSeqIndividual::getSequence: sequence set not found (" + id
			+ ").";
		throw(Exception(mes));
	}
	try {
		return const_cast<const VectorSequenceContainer *>(it->second)->getSequence(i);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw(ioobe);
	}
}

vector<string> MultiSeqIndividual::getSequencesKeys() const {
	vector<string> keys;
	map<string, VectorSequenceContainer *>::const_iterator it;
	for (it = _sequences.begin() ; it != _sequences.end() ; it++)
		keys.push_back(it->first);
	return keys;
}

bool MultiSeqIndividual::hasSequences() const {
	return _sequences.size() != 0;
}

unsigned int MultiSeqIndividual::getNumberOfSequenceSet() const {
	return _sequences.size();
}

unsigned int MultiSeqIndividual::getNumberOfSequences(const string & id) const
	throw (Exception) {
	map<string, VectorSequenceContainer *>::const_iterator it;
	it = _sequences.find(id);
	// Test existence of id in the map.
	if (it == _sequences.end()) {
		string mes = "MultiSeqIndividual::getSequence: sequence set not found (" + id
			+ ").";
		throw(Exception(mes));
	}
	
	return const_cast<const VectorSequenceContainer *>(it->second)->getNumberOfSequences();
}

// MultilocusGenotype

void MultiSeqIndividual::addGenotype(const MultilocusGenotype & genotype) {
	_genotype = new MultilocusGenotype(genotype);
}

const MultilocusGenotype * MultiSeqIndividual::getGenotype() const throw (NullPointerException) {
		return _genotype;
}

bool MultiSeqIndividual::hasGenotype() const {
	return _genotype != NULL;
}
