/*
 * File Individual.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Tuesday April 27 2004
 */

#include "Individual.h"

//** Class constructor: *******************************************************/
Individual::Individual() {
	this->_date = NULL;
	this->_coord = NULL;
	this->_locality = NULL;
}

Individual::Individual(const Individual &ind) {
	this->_sex = ind.getSex();
	this->_date = ind.getDate();
	this->_coord = ind.getCoord();
	this->_locality = ind.getLocality();
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
	this->_date = ind.getDate();
	this->_coord = ind.getCoord();
	this->_locality = ind.getLocality();
	return * this;
}

void Individual::setSex(const unsigned short sex) {
	_sex = sex;
}

unsigned short Individual::getSex() const {
	return _sex;
}

void Individual::setDate(const Date & date) {
	if (_date == NULL) {
		_date = newDate(date);
	}
	else if (* _date != date) {
		delete _date;
		_date = new Date(date);
	}
}

Date Individual::getDate() const throw(NullPointerException) {
	if (_date != NULL)
		return new date(_date);
	else
		throw(NullPointerException("Date not set."));
}

void Individual::setCoord(const Coord<double> & coord) {
	if (_coord == NULL) {
		_coord = new Coord<double>(coord);
	}
	else if	(* _coord != coord) {
		delete _coord;
		_coord = new Coord<double>(coord);
	}
}

void Individual::setCoord(const double x, const double y) {
	if (_coord == NULL) {
		_coord = new Coord<double>(x, y);
	}
	else if (this->getX() != x || this->getY() != y) {
		delete _coord;
		_coord = new Coord<double>(x, y);
	}
}

Coord<double> Individual::getCoord() const throw(NullPointerException) {
	if (_coord != NULL)
		return new Coord<double>(_coord);
	else
		throw(NullPointerException("Coord not set."));
}

void Individual::setX(const double x) throw(NullPointerException) {
	if (_coord != NULL)
		_coord->setX(x);
	else
		throw(NullPointerException("Coord not initialized."));
}

void Individual::setY(const double y) throw(NullPointerException) {
	if (_coord != NULL)
		_coord->setY(y);
	else
		throw(NullPointerException("Coord not initialized."));
}

double Individual::getX() const throw(NullPointerException) {
	if (_coord != NULL)
		return _coord->getX();
	else
		throw(NullPointerException("X not set."));
}

double Individual::getY() const throw(NullPointerException) {
	if (_coord != NULL)
		return _coord->getY();
	else
		throw(NullPointerException("Y not set."));
}

void Individual::setLocality(const Locality<double> & locality) {
	if (_locality == NULL)
		_locality = new Locality<double>(locality);
	else if (* _locality != locality) {
		delete _locality;
		_locality = new Locality<double>(locality);
	}
}

void Individual::setLocality(const string name, const Coord<double> & coord) {
	if (_locality == NULL)
		_locality = new Locality<double>(name, coord);
	else if (_locality->getName() != name || _locality->getCoord() != coord) {
		delete _locality;
		_locality = new Locality<double>(name, coord);
	}
}

void Individual::setLocality(const string name, const double x, const double y) {
	if (_locality == NULL)
		_locality = new Locality<double>(name, x, y);
	else if (_locality->getName() != name || _locality->getX() != x || _locality->getY() != y) {
		delete _locality;
		_locality = new Locality<double>(name, x, y);
	}
}

Locality<double> Individual::getLocality() const throw(NullPointerException) {
	if (_locality != NULL)
		return new Locality<double>(_locality);
	else
		throw(NullPointerException("Locality  not set."));
}
