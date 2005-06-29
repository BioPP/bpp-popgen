/*
 * File Locality.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Thursday July 29 2004
 *
 * Copyright (C) 2004 Sylvain Gaillard and the
 *                    PopGenLib Development Core Team
 *
 * This file is part of PopGenLib.
 *
 * PopGenLib is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * PopGenLib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PopGenLib; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

// Secured inclusion of header's file
#ifndef _LOCALITY_H_
#define _LOCALITY_H_

#include "Coord.h"
#include <string>
using namespace std;
/**
 * @brief The Locality class.
 *
 * This is a class derivated from the Coord class.
 * It's a Coord with a name.
 */
template <class T> class Locality : public Coord<T> {
	public: // Constructors and destructor
		/**
		 * @brief Build a new locality with name and coordinates.
		 *
		 * @param name The name of the locality.
		 * @param x The longitude.
		 * @param y The latitude.
		 */
		Locality(const string name, const T x=0, const T y=0);

		/**
		 * @brief Build a new locality with name and coordinates.
		 *
		 * @param name The name of the locality.
		 * @param coord The coordinates of the locality.
		 */
		Locality(const string name, const Coord<T> coord);

		/**
		 * @brief The Locality copy constructor.
		 */
		Locality(const Locality<T> & locality);

		/**
		 * @brief Destroy a locality.
		 */
		~Locality();

	public: // Methodes
		/**
		 * @brief Implements the Clonable interface.
		 */
		Clonable * clone() const;

		/**
		 * @brief The Locality copy operator.
		 *
		 * @return A ref toward the assigned Locality.
		 */
		Locality & operator= (const Locality & locality);

		/**
		 * @brief The == operator.
		 *
		 * returns true if both name and coordinates are identical between the two Locality objects.
		 */
		virtual bool operator== (const Locality<T> & locality) const;

		/**
		 * @brief The != operator.
		 */
		virtual bool operator!= (const Locality<T> & locality) const;

		/**
		 * @brief Set the name of the locality.
		 */
		void setName(const string name);

		/**
		 * @brief Get the name of the locality.
		 */
		string getName() const;

	protected:
		string _name;
};

//** Class constructor: *******************************************************/
template <class T> Locality<T>::Locality(const string name, const T x, const T y) {
	_name = name;
	Coord<T>::_x = x;
	Coord<T>::_y = y;
}

template <class T> Locality<T>::Locality(const string name, const Coord<T> coord) {
	_name = name;
	Coord<T>::_x = coord.getX();
	Coord<T>::_y = coord.getY();
}

template <class T> Locality<T>::Locality(const Locality<T> & locality) {
	this->_x = locality.getX();
	this->_y = locality.getY();
	this->_name = locality.getName();
}

//** Class destructor: *******************************************************/
template <class T> Locality<T>::~Locality() {}

//** Clonable interface: *****************************************************/
template <class T> Clonable * Locality<T>::clone() const {
	return new Locality<T>(* this);
}

//** Copy operator: **********************************************************/
template <class T> Locality<T> & Locality<T>::operator= (const Locality<T> & locality) {
	this->_x = locality.getX();
	this->_y = locality.getY();
	this->_name = locality.getName();
	return * this;
}

//** Comparison operators: ***************************************************/
template <class T> bool Locality<T>::operator== (const Locality<T> & locality) const {
	return (Coord<T>::_x == locality.getX() && Coord<T>::_y == locality.getY() && _name == locality.getName());
}

template <class T> bool Locality<T>::operator!= (const Locality<T> & locality) const {
	return !(locality == *this);
}

//** Assignation opperators: *************************************************/
template <class T> void Locality<T>::setName(string name) {
	_name = name;
}

//** Consultation opperators: ************************************************/
template <class T> string Locality<T>::getName() const {
	return _name;
}
#endif // _LOCALITY_H_
