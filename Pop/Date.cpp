/*
 * File Date.cpp
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

// From Utils
#include <Utils/TextTools.h>

//From Local
#include "Date.h"

//** Class constructor: *******************************************************/
Date::Date(const int day, const int month, const int year) throw(BadIntegerException) {
	if (day >= 1 && day <= 31)
		_day = day;
	else
		throw(BadIntegerException("Date::Date: day must be in [1;31].", day));
	if (month >= 1 && month <= 12)
		_month = month;
	else
		throw(BadIntegerException("Date::Date: month must be in [1;12].", month));
	_year = year;
}

Date::Date(const Date & date) {
	this->_day = date.getDay();
	this->_month = date.getMonth();
	this->_year = date.getYear();
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

void Date::setDate(const int day, const int month, const int year) throw(BadIntegerException) {
	if (day >= 1 && day <= 31)
		_day = day;
	else
		throw(BadIntegerException("Date::Date: day must be in [1;31].", day));
	if (month >= 1 && month <= 12)
		_month = month;
	else
		throw(BadIntegerException("Date::Date: month must be in [1;12].", month));
	_year = year;
}

void Date::setYear(const int year) {
	_year = year;
}

void Date::setMonth(const int month) throw(BadIntegerException) {
	if (month >= 1 && month <= 12)
		_month = month;
	else
		throw(BadIntegerException("Date::Date: month must be in [1;12].", month));
}

void Date::setDay(const int day) throw(BadIntegerException) {
	if (day >= 1 && day <= 31)
		_day = day;
	else
		throw(BadIntegerException("Date::Date: day must be in [1;31].", day));
}

Date * Date::getDate() const {
	return new Date(* this);
}

string Date::getDateStr() const {
	string date, uDay="", uMonth="";
	if (_day < 10) uDay="0";
	if (_month < 10) uMonth="0";
	date = uDay + TextTools::toString(_day) + uMonth + TextTools::toString(_month) + TextTools::toString(_year);
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

bool Date::operator==(const Date & date) const {
	if (_day == date.getDay() && _month == date.getMonth() && _year == date.getYear())
		return true;
	else
		return false;
}

bool Date::operator<(const Date & date) const {
	if (_year < date.getYear() || (_month < date.getMonth() && _year == date.getYear()) || (_day < date.getDay() && _month == date.getMonth() && _year == date.getYear()))
		return true;
	else
		return false;
}
