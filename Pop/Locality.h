/*
 * File Locality.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Monday April 05 2004
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
template <class T> class Locality : public Coord <T> {
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
		 * @brief Destroy a locality.
		 */
		~Locality();

		/**
		 * @brief Set the name of the locality.
		 */
		void setName(const string name);

		/**
		 * @brief Get the name of the locality.
		 */
		const string getName();

	protected:
		string _name;
};

//** Class constructor: *******************************************************/
template <class T> Locality<T>::Locality(const string name, const T x, const T y) {
	_name = name;
	_x = x;
	_y = y;
}

template <class T> Locality<T>::Locality(const string name, const Coord<T> coord) {
	_name = name;
	_x = coord.getX();
	_y = coord.getY();
}

//** Class destructor: *******************************************************/
template <class T> Locality<T>::~Locality() {}

//** Assignation opperators: *************************************************/
template <class T> void Locality<T>::setName(string name) {
	_name = name;
}

//** Consultation opperators: ************************************************/
template <class T> const string Locality<T>::getName() {
	return _name;
}
#endif // _LOCALITY_H_
