/*
 * File CoordsTools.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Wednesday May 05 2004
 */

// Secured inclusion of header's file
#ifndef _COORDS_TOOLS_H_
#define _COORDS_TOOLS_H_

// From local
#include "Coord.h"

/**
 * @brief Some functions to deal with coordinates.
 */

template <class T> class CoordsTools {
	public:
		CoordsTools();
		~CoordsTools();

		/**
		 * @brief Get the distance between two Coord objects.
		 *
		 * @param co1 A Coord object.
		 * @param co2 An other Coord object.
		 */
		static double getDistanceBetween(const Coord<T> co1, const Coord<T> co2) const;
};

template <class T> CoordsTools<T>::CoordsTools() {}

template <class T> CoordsTools<T>::~CoordsTools() {}

template <class T> double CoordsTools<T>::getDistanceBetween(const Coord<T> co1, const Coord<T> co2) const {
	T base, height;
	base = co1.getX() - co2.getX();
	height = co1.getY() - co2.getY();
	base = base * base;
	height = height * height;
	return sqrt((double)base + (double)height);
}

#endif // _COORDS_TOOLS_H_
