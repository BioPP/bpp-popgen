/*
 * File Individual.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Monday April 05 2004
 */

#include "Individual.h"

//** Class constructor: *******************************************************/
Individual::Individual() {
	this->_date = NULL;
	this->_coord = NULL;
}

Individual::Individual(const Individual &ind) {
	this->_sex = ind.getSex();
	this->_date = ind.getDate();
	this->_coord = ind.getCoord();
}

//** Class destructor: *******************************************************/
Individual::~Individual () {}

//** Other methodes: *********************************************************/
Clonable * Individual::clone() const {
	return new Individual(* this);
}

Individual & Individual::operator= (const Individual & ind) {
	this->_sex = ind.getSex();
	this->_date = ind.getDate();
	this->_coord = ind.getCoord();
	return * this;
}

void Individual::setSex(const unsigned short sex) {
	_sex = sex;
}

unsigned short Individual::getSex() const {
	return _sex;
}

void Individual::setDate(const Date & date) {
	delete _date;
	_date = new Date(date.getDay(), date.getMonth(), date.getYear());
}

Date * Individual::getDate() const {
	return _date;
}

void Individual::setCoord(const Coord<double> & coord) {
	delete _coord;
	_coord = new Coord<double>(coord.getX(), coord.getY());
}

Coord<double> * Individual::getCoord() const {
	return _coord;
}

void Individual::setX(const double x) {
	if (_coord != NULL)
		_coord->setX(x);
}

void Individual::setY(const double y) {
	if (_coord != NULL)
		_coord->setY(y);
}

double Individual::getX() const {
	if (_coord != NULL)
		return _coord->getX();
	else
		return 0;
}

double Individual::getY() const {
	if (_coord != NULL)
		return _coord->getY();
	else
		return 0;
}
