/*
 * File PolymorphismMultiGContainerTools.h
 * Author : Sylvain Gailard <yragael2001@yahoo.fr>
 * Last modification : Friday September 24 2004
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

//Secured inclusion of header's file
#ifndef _POLYMORPHISMMULTIGCONTAINERTOOLS_H_
#define _POLYMORPHISMMULTIGCONTAINERTOOLS_H_

// From the STL
#include <vector>
using namespace std;

//From the PolyLib library
#include "PolymorphismMultiGContainer.h"

//From the NumCalc library
#include <NumCalc/RandomTools.h>

class PolymorphismMultiGContainerTools {
	public:
		virtual ~PolymorphismMultiGContainerTools();
		
	public:
		static PolymorphismMultiGContainer permutMultiG(const PolymorphismMultiGContainer & pmgc);

};

#endif // _POLYMORPHISMMULTIGCONTAINERTOOLS_H_
