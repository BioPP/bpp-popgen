/*
 * File AbstractODataSet.cpp
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

#include "AbstractODataSet.h"

// From STL
#include <fstream>
using namespace std;

AbstractODataSet::~AbstractODataSet() {}

void AbstractODataSet::write(const string & path, const DataSet & data_set, bool overwrite) const throw (Exception) {
	ofstream output(path.c_str(), overwrite ? (ios::out) : (ios::out|ios::app));
	write(output, data_set);
	output.close();
}
