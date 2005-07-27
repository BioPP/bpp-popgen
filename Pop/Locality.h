/*
 * File Locality.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Thursday July 29 2004
 *
*/
/*
Copyright or © or Copr. CNRS, (November 17, 2004)


This software is a computer program whose purpose is to provide classes
for sequences analysis.

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
