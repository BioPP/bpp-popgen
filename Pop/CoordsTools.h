/*
 * File CoordsTools.h
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
