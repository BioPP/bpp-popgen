// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "AbstractIDataSet.h"

using namespace bpp;

// From STL
#include <fstream>

using namespace std;

AbstractIDataSet::~AbstractIDataSet() {}

void AbstractIDataSet::read(const std::string& path, DataSet& data_set)
{
  ifstream input(path.c_str(), ios::in);
  read(input, data_set);
  input.close();
}

DataSet* AbstractIDataSet::read(std::istream& is)
{
  DataSet* data_set = new DataSet();
  read(is, *data_set);
  return data_set;
}

DataSet* AbstractIDataSet::read(const std::string& path)
{
  DataSet* data_set = new DataSet();
  read(path, *data_set);
  return data_set;
}
