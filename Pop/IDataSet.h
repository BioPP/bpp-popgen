/*
 * File IDataSet.h
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
#ifndef _IDATASET_H_
#define _IDATASET_H_

#include "IODataSet.h"

// From Utils
#include <Utils/Exceptions.h>

/**
 * @brief The IDataSet interface.
 */
class IDataSet : public IODataSet {
	public: // Class destructor
		virtual ~IDataSet();

	public:
		/**
		 * @brief Read a DataSet on istream.
		 */
		virtual void read(istream & is, DataSet & data_set) throw (Exception) = 0;

		/**
		 * @brief Read a DataSet from a text file.
		 */
		virtual void read(const string & path, DataSet & data_set) throw (Exception) = 0;

		/**
		 * @brief Read istream and return a DataSet.
		 */
		virtual DataSet * read(istream & is) throw (Exception) = 0;

		/**
		 * @brief Read a text file and return a DataSet.
		 */
		virtual DataSet * read(const string & path) throw (Exception) = 0;
};

#endif // _IDATASET_H_
