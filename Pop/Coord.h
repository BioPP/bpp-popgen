/*
 * File Coord.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Wednesday May 05 2004
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
