/*
 * File Coord.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Thursday July 29 2004
 *
*/
/*
Copyright or © or Copr. CNRS, (November 17, 2004)


This software is a computer program whose purpose is to provide classes
for population genetics analysis.

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
