/*
 * File Coord.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Monday April 05 2004
 */

// Secured inclusion of header's file
#ifndef _COORD_H_
#define _COORD_H_

/**
 * @brief The Coord class.
 *
 * This is a simple class designed to store the coord of a point.
 * The type of the two coord is defined as a template.
 */
template <class T> class Coord {
	public: // Constructors and destructor :

		/**
		 * @brief Build a new Coord from two values.
		 * 
		 * The two values are set to 0 if no parametre is given to the constructor.
		 *
		 * @param x The longitude or abscice.
		 * @param y The latitude or ordonne.
		 */
		Coord(const T x=0, const T y=0);

		/**
		 * @brief Destroy the Coord object.
		 */
		~Coord();

		/**
		 * @brief The Coord copy operator.
		 *
		 * @return A ref toward the assigned Coord.
		 */
		Coord & operator= (const Coord & coord);
		
		/**
		 * @brief Set the two values.
		 */
		void setCoord(const T x, const T y);

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

	protected:
		T _x;
		T _y;
};

//** Class constructor: *******************************************************/
template <class T> Coord<T>::Coord(const T x, const T y) {
	_x = x;
	_y = y;
}

//** Class destructor: *******************************************************/
template <class T> Coord<T>::~Coord() {}

//** Copy operator: **********************************************************/
template <class T> Coord<T> & Coord<T>::operator= (const Coord<T> & coord) {
	this->_x = coord.getX();
	this->_y = coord.getY();
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
template <class T> T Coord<T>::getX() const {
	return _x;
}

template <class T> T Coord<T>::getY() const {
	return _y;
}
#endif // _COORD_H_
