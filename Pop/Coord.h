/*
 * File Coord.h
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

// Secured inclusion of header's file
#ifndef _COORD_H_
#define _COORD_H_

// From Utils
#include<Utils/Clonable.h>
/**
 * @brief The Coord class.
 *
 * This is a simple class designed to store the coordinates of a point.
 * The type of the two coordinates is defined as a template.
 */
template <class T> class Coord : public Clonable {
	public: // Constructors and destructor :

		/**
		 * @brief Build a new Coord from two values.
		 * 
		 * The two values are set to 0 if no parametre is given to the constructor.
		 *
		 * @param x The longitude or abscissa.
		 * @param y The latitude or ordinate.
		 */
		Coord(const T x=0, const T y=0);
		
		/**
		 * @brief The Coord copy constructor.
		 */
		Coord(const Coord<T> & coord);
		
		/**
		 * @brief Destroy the Coord object.
		 */
		~Coord();

	public: // Methodes

		/**
		 * @brief Implement the Clonable interface.
		 */
		Clonable * clone() const;

		/**
		 * @brief The Coord copy operator.
		 *
		 * @return A ref toward the assigned Coord.
		 */
		Coord & operator= (const Coord<T> & coord);
		
		/**
		 * @brief Set the two values.
		 */
		void setCoord(const T x, const T y);

		/**
		 * @brief Get the coordinates.
		 *
		 * @return A pointer toward a Coord object with the same coordinates.
		 */
		Coord<T> * getCoord() const;

		/**
		 * @brief Set only the longitude.
		 */
		void setX(const T x);

		/**
		 * @brief Set only the latitude.
		 */
		void setY(const T y);

		/**
		 * @brief Get the longitude.
		 */
		T getX() const;

		/**
		 * @brief Get the latitude.
		 */
		T getY() const;

		/**
		 * @brief Compares two Coord objets.
		 *
		 * Return true if the coordinates of the 2 Coords are equals.
		 */
		bool hasSameCoordsAs(const Coord<T> & coord) const;

		/**
		 * @brief The == operator.
		 *
		 * Return true if the coordinates of the 2 Coords are equals.
		 * Does the same as the asSameCoords() methode.
		 */
		virtual bool operator== (const Coord<T> & coord) const;

		/**
		 * @brief The != operator.
		 */
		virtual bool operator!= (const Coord<T> & coord) const;

	protected:
		T _x;
		T _y;
};

//** Class constructor: *******************************************************/
template <class T> Coord<T>::Coord(const T x, const T y) {
	_x = x;
	_y = y;
}

template <class T> Coord<T>::Coord(const Coord<T> & coord) {
	this->_x = coord.getX();
	this->_y = coord.getY();
}
//** Class destructor: *******************************************************/
template <class T> Coord<T>::~Coord() {}

//** Clonable interface: *****************************************************/
template <class T> Clonable * Coord<T>::clone() const {
	return new Coord<T>(* this);
}

//** Copy operator: **********************************************************/
template <class T> Coord<T> & Coord<T>::operator= (const Coord<T> & coord) {
	this->_x = coord.getX();
	this->_y = coord.getY();
	return * this;
}

//** Assignation opperators: *************************************************/
template <class T> void Coord<T>::setCoord(const T x, const T y) {
	_x = x;
	_y = y;
}

template <class T> void Coord<T>::setX(const T x) {
	_x = x;
}

template <class T> void Coord<T>::setY(const T y) {
	_y = y;
}

//** Consultation opperators: ************************************************/
template <class T> Coord<T> * Coord<T>::getCoord() const {
	return new Coord<T>(_x, _y);
}

template <class T> T Coord<T>::getX() const {
	return _x;
}

template <class T> T Coord<T>::getY() const {
	return _y;
}

//** Comparison operators: ***************************************************/
template <class T> bool Coord<T>::operator== (const Coord<T> & coord) const {
	if (_x == coord.getX() && _y == coord.getY())
		return true;
	else
		return false;
}

template <class T> bool Coord<T>::operator!= (const Coord<T> & coord) const {
	return !(coord == *this);
}

template <class T> bool Coord<T>::hasSameCoordsAs(const Coord<T> & coord) const {
	if (_x == coord.getX() && _y == coord.getY())
		return true;
	else
		return false;
}

#endif // _COORD_H_
