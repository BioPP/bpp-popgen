/*
 * File Individual.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Monday May 24 2004
 */

#include "Individual.h"

//** Class constructor: *******************************************************/
Individual::Individual() {
	_id = "";
	_sex = 0;
	_date = NULL;
	_coord = NULL;
	_locality = NULL;
}

Individual::Individual(const string id,
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

Individual::Individual(const Individual &ind) {
	this->_id = ind.getId();
	this->_sex = ind.getSex();
	this->_date = ind.hasDate() ? new Date(* ind.getDate()) : NULL;
	this->_coord = ind.hasCoord() ? new Coord<double>(* ind.getCoord()) : NULL;
	this->_locality = ind.getLocality();
}

//** Class destructor: *******************************************************/
Individual::~Individual () {
	delete this->_date;
	delete this->_coord;
}

//** Other methodes: *********************************************************/
Clonable * Individual::clone() const {
	return new Individual(* this);
}

Individual & Individual::operator= (const Individual & ind) {
	this->_id = ind.getId();
	this->_sex = ind.getSex();
	this->_date = ind.hasDate() ? new Date(* ind.getDate()) : NULL;
	this->_coord = ind.hasCoord() ? new Coord<double>(* ind.getCoord()) : NULL;
	this->_locality = ind.getLocality();
	return * this;
}

// Id
void Individual::setId(const string id) {
	_id = id;
}

string Individual::getId() const {
	return _id;
}

// Sex
void Individual::setSex(const unsigned short sex) {
	_sex = sex;
}

unsigned short Individual::getSex() const {
	return _sex;
}

// Date
void Individual::setDate(const Date & date) {
	if (!hasDate()) {
		_date = new Date(date);
	}
	else if (* _date != date) {
		delete _date;
		_date = new Date(date);
	}
}

Date * Individual::getDate() const throw(NullPointerException) {
	if (hasDate())
		return new Date(* _date);
	else
		throw(NullPointerException("Individual::getDate: no date associated to this individual."));
}

bool Individual::hasDate() const {
	return _date != NULL;
}

// Coord
void Individual::setCoord(const Coord<double> & coord) {
	if (!hasCoord()) {
		_coord = new Coord<double>(coord);
	}
	else if	(* _coord != coord) {
		delete _coord;
		_coord = new Coord<double>(coord);
	}
}

void Individual::setCoord(const double x, const double y) {
	if (!hasCoord()) {
		_coord = new Coord<double>(x, y);
	}
	else if (this->getX() != x || this->getY() != y) {
		delete _coord;
		_coord = new Coord<double>(x, y);
	}
}

Coord<double> * Individual::getCoord() const throw(NullPointerException) {
	if (hasCoord())
		return new Coord<double>(* _coord);
	else
		throw(NullPointerException("Individual::getCoord: no coord associated to this individual."));
}

bool Individual::hasCoord() const {
	return _coord != NULL;
}

void Individual::setX(const double x) throw(NullPointerException) {
	if (hasCoord())
		_coord->setX(x);
	else
		throw(NullPointerException("Individual::setX: no coord associated to this individual."));
}

void Individual::setY(const double y) throw(NullPointerException) {
	if (hasCoord())
		_coord->setY(y);
	else
		throw(NullPointerException("Individual::setY: no coord associated to this individual."));
}

double Individual::getX() const throw(NullPointerException) {
	if (hasCoord())
		return _coord->getX();
	else
		throw(NullPointerException("Individual::getX: no coord associated to this individual."));
}

double Individual::getY() const throw(NullPointerException) {
	if (hasCoord())
		return _coord->getY();
	else
		throw(NullPointerException("Individual::getY: no coord associated to this individual."));
}

// Locality
void Individual::setLocality(Locality<double> * locality) {
	_locality = locality;
}

Locality<double> * Individual::getLocality() const {
	return _locality;
}

bool Individual::hasLocality() const {
	return _locality != NULL;
}

// Sequences
void Individual::addSequence(const string & id, const Sequence & sequence)
throw (Exception) {
	try {
		_sequences[id]->addSequence(sequence);
	}
	catch (AlphabetMismatchException ame)
	{
		throw(ame);
	}
	catch (Exception e)
	{
		throw(e);
	}
}

const Sequence * Individual::getSequence(const string & id, const string & name)
const throw(Exception){
	map<string, VectorSequenceContainer *>::const_iterator it;
	it = _sequences.find(id);
	// Test existence of id in the map.
	if (it == _sequences.end()) {
		string mes = "Individual::getSequence: sequence set not found (" + id
			+ ").";
		throw(Exception(mes));
	}
	try {
		return const_cast<const VectorSequenceContainer *>(it->second)->getSequence(name);
	}
	catch (SequenceNotFoundException snfe) {
		throw(snfe);
	}
}

const Sequence * Individual::getSequence(const string & id, unsigned int i)
const throw(Exception) {
	map<string, VectorSequenceContainer *>::const_iterator it;
	it = _sequences.find(id);
	// Test existence of id in the map.
	if (it == _sequences.end()) {
		string mes = "Individual::getSequence: sequence set not found (" + id
			+ ").";
		throw(Exception(mes));
	}
	try {
		return const_cast<const VectorSequenceContainer *>(it->second)->getSequence(i);
	}
	catch (IndexOutOfBoundsException ioobe) {
		throw(ioobe);
	}
}

vector<string> Individual::getSequencesKeys() const {
	vector<string> keys;
	map<string, VectorSequenceContainer *>::const_iterator it;
	for (it = _sequences.begin() ; it != _sequences.end() ; it++)
		keys.push_back(it->first);
	return keys;
}
