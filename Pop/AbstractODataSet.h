/*
 * File AbstractODataSet.h
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
#ifndef _ABSTRACTODATASET_H_
#define _ABSTRACTODATASET_H_

#include "ODataSet.h"

class AbstractODataSet : public ODataSet {
	public: // Calss destructor
		virtual ~AbstractODataSet();

	public:
		/**
		 * @name The ODataSet interface.
		 * @{
		 */
		virtual void write(ostream & os, const DataSet & data_set) const throw (Exception) = 0;
		virtual void write(const string & path, const DataSet & data_set, bool overwrite) const throw (Exception);
		/**
		 * @}
		 */
};

#endif // _ABSTRACTODATASET_H_
