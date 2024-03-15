// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "AbstractODataSet.h"

using namespace bpp;

// From STL
#include <fstream>

using namespace std;

AbstractODataSet::~AbstractODataSet() {}

void AbstractODataSet::write(const string& path, const DataSet& data_set, bool overwrite) const
{
  ofstream output(path.c_str(), overwrite ? (ios::out) : (ios::out | ios::app));
  write(output, data_set);
  output.close();
}
