/*
 * File Date.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Wednesday April 14 2004
 */

// From Utils
#include <Utils/TextTools.h>

//From Local
#include "Date.h"

//** Class constructor: *******************************************************/
Date::Date(const int day, const int month, const int year) {
	_day = day;
	_month = month;
	_year = year;
}
//** Class destructor: ********************************************************/
Date::~Date() {}
//** Other methodes: **********************************************************/
Date & Date::operator= (const Date & date) {
	this->_day = date.getDay();
	this->_month = date.getMonth();
	this->_year = date.getYear();
	return * this;
}

void Date::setDate(const int day, const int month, const int year) {
	_day = day;
	_month = month;
	_year = year;
}

string Date::getDate() const {
	string date;
	date = TextTools::toString(_day) + TextTools::toString(_month) + TextTools::toString(_year);
	return date;
}

int Date::getYear() const {
	return _year;
}

int Date::getMonth() const {
	return _month;
}

int Date::getDay() const {
	return _day;
}
