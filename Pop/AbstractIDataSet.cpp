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

#include "AbstractIDataSet.h"

// From STL
#include <fstream>
using namespace std;

AbstractIDataSet::~AbstractIDataSet() {}

void AbstractIDataSet::read(const string & path, DataSet & data_set) throw (Exception) {
	ifstream input(path.c_str(), ios::in);
	read(input, data_set);
	input.close();
}

DataSet * AbstractIDataSet::read(istream & is) throw (Exception) {
	DataSet * data_set = new DataSet();
	read(is, * data_set);
	return data_set;
}

DataSet * AbstractIDataSet::read(const string & path) throw (Exception) {
	DataSet * data_set = new DataSet();
	read(path, * data_set);
	return data_set;
}
