/*
 * File AbstractIDataSet.h
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
#ifndef _ABSTRACTIDATASET_H_
#define _ABSTRACTIDATASET_H_

#include "IDataSet.h"

// From Utils
#include <Utils/Exceptions.h>

class AbstractIDataSet : public IDataSet {
	public: // Class destructor
		virtual ~AbstractIDataSet();

	public:
		/**
		 * @name The IDataSet interface.
		 * @{
		 */
		virtual void read(istream & is, DataSet & data_set) throw (Exception) = 0;

		virtual void read(const string & path, DataSet & data_set) throw (Exception);

		virtual DataSet * read(istream & is) throw (Exception);

		virtual DataSet * read(const string & path) throw (Exception);
		/**
		 * @}
		 */
};

#endif // _ABSTRACTIDATASET_H_
