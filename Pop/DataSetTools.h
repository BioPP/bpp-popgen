/*
 * File DataSetTools.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Wednesday August 04 2004
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
#ifndef _DATASETTOOLS_H_
#define _DATASETTOOLS_H_

// From STL
#include <set>
using namespace std;

// From Utils
#include <Utils/TextTools.h>
#include <Utils/Exceptions.h>

// From SeqLib
#include <Seq/OrderedSequenceContainer.h>

// From local PopLib
#include "DataSet.h"
#include "PolymorphismSequenceContainer.h"

/**
 * @brief A set of tools for DataSet.
 */
class DataSetTools {
	public:
		/**
		 * @brief General method to build a DataSet from an OrderedSequenceContainer.
		 */
		static DataSet * buildDataSet(const OrderedSequenceContainer & osc) throw (Exception);

		/**
		 * @brief Specific methode to build a DataSet from a PolymorphismSequenceContainer.
		 */
		static DataSet * buildDataSet(const PolymorphismSequenceContainer & psc) throw (Exception);
};

#endif // _DATASETTOOLS_H_
