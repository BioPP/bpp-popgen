/*
 * File MultiSeqIndividual.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Thursday July 29 2004
 *
 * Copyright (C) 2004 Sylvain Gaillard and the
 *                    PopLib Development Core Team
 *
 * This file is part of PopLib.
 *
 * PopLib is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * PopLib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PopLib; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
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
	this->_genotype = ind.hasGenotype() ? new Genotype(* ind.getGenotype()) : NULL;
}

//** Class destructor: *******************************************************/
MultiSeqIndividual::~MultiSeqIndividual () {
	delete this->_date;
	delete this->_coord;
}

//** Other methodes: *********************************************************/
Clonable * MultiSeqIndividual::clone() const {
	return new MultiSeqIndividual(* this);
}

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
	this->_genotype = ind.hasGenotype() ? new Genotype(* ind.getGenotype()) : NULL;
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

// Genotype

void MultiSeqIndividual::addGenotype(const Genotype & genotype) {
	_genotype = new Genotype(genotype);
}

const Genotype * MultiSeqIndividual::getGenotype() const throw (NullPointerException) {
		return _genotype;
}

bool MultiSeqIndividual::hasGenotype() const {
	return _genotype != NULL;
}
