/*
 * File IODataSet.h
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
#ifndef _IODATASET_H_
#define _IODATASET_H_

#include "DataSet.h"

// From STL
#include <iostream>
#include <fstream>
using namespace std;

/**
 * @brief Interface for input/ouput with DataSet.
 *
 * IODataSet is a virtual class.
 * This is an interface to declare commune methodes for in/out action on DataSet.
 */
class IODataSet {
	public: // Class destructor
		virtual ~IODataSet();

	public:
		/**
		 * @brief Get the format's name.
		 */
		virtual const string getFormatName() = 0;

		/**
		 * @brief Get a description of the format.
		 */
		virtual const string getFormatDescription() = 0;
};

#endif // _IODATASET_H_
