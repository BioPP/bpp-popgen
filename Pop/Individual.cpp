/*
 * File Individual.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Wednesday April 28 2004
 */

#include "Individual.h"

//** Class constructor: *******************************************************/
Individual::Individual() {
	this->_sex = 0;
	this->_date = NULL;
	this->_coord = NULL;
	this->_locality = NULL;
}

Individual::Individual(const Date & date, const Coord<double> & coord, const Locality<double> & locality, const unsigned short sex) {
	this->_sex = sex;
	this->_date = new Date(date);
	this->_coord = new Coord<double>(coord);
	this->_locality = new Locality<double>(locality);
}

Individual::Individual(const Individual &ind) {
	this->_sex = ind.getSex();
	this->_date = ind.hasDate() ? new Date(* ind.getDate()) : NULL;
	this->_coord = ind.hasCoord() ? new Coord<double>(* ind.getCoord()) : NULL;
	this->_locality = ind.hasLocality() ? new Locality<double>(* ind.getLocality()) : NULL;
}

//** Class destructor: *******************************************************/
Individual::~Individual () {
	delete this->_date;
	delete this->_coord;
	delete this->_locality;
}

//** Other methodes: *********************************************************/
Clonable * Individual::clone() const {
	return new Individual(* this);
}

Individual & Individual::operator= (const Individual & ind) {
	this->_sex = ind.getSex();
	this->_date = ind.hasDate() ? new Date(* ind.getDate()) : NULL;
	this->_coord = ind.hasCoord() ? new Coord<double>(* ind.getCoord()) : NULL;
	this->_locality = ind.hasLocality() ? new Locality<double>(* ind.getLocality()) : NULL;
	return * this;
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
void Individual::setLocality(const Locality<double> & locality) {
	if (!hasLocality())
		_locality = new Locality<double>(locality);
	else if (* _locality != locality) {
		delete _locality;
		_locality = new Locality<double>(locality);
	}
}

void Individual::setLocality(const string name, const Coord<double> & coord) {
	if (!hasLocality())
		_locality = new Locality<double>(name, coord);
	else if (_locality->getName() != name || * (_locality->getCoord()) != coord) {
		delete _locality;
		_locality = new Locality<double>(name, coord);
	}
}

void Individual::setLocality(const string name, const double x, const double y) {
	if (!hasLocality())
		_locality = new Locality<double>(name, x, y);
	else if (_locality->getName() != name || _locality->getX() != x || _locality->getY() != y) {
		delete _locality;
		_locality = new Locality<double>(name, x, y);
	}
}

Locality<double> * Individual::getLocality() const throw(NullPointerException) {
	if (hasLocality())
		return new Locality<double>(* _locality);
	else
		throw(NullPointerException("Individual::getLocality: no locality associated to this individual."));
}

bool Individual::hasLocality() const {
	return _locality != NULL;
}
