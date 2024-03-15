// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "DarwinDon.h"

#include <Bpp/Io/OutputStream.h>

using namespace bpp;
using namespace std;

DarwinDon::DarwinDon() {}

DarwinDon::~DarwinDon() {}

void DarwinDon::write(ostream& os, const DataSet& dataset) const
{
  if (!os)
    throw IOException("DarwinDon::write: fail to open stream.");
  StlOutputStreamWrapper out(&os);
  (out << "@DARwin 5.0 - DON").endLine();
  size_t ind_nbr = 0;
  for (size_t i = 0; i < dataset.getNumberOfGroups(); i++)
  {
    ind_nbr += dataset.getNumberOfIndividualsInGroup(i);
  }
  vector<string> header;
  header.push_back("N°");
  header.push_back("Name");
  (out << ind_nbr << "\t" << header.size() - 1).endLine();
  VectorTools::print(header, out, "\t");
  // size_t ind_index = 0;
  for (size_t i = 0; i < dataset.getNumberOfGroups(); i++)
  {
    size_t ind_nbr_ig = dataset.getNumberOfIndividualsInGroup(i);
    for (size_t j = 0; j < ind_nbr_ig; j++)
    {
      (out << j + (i * ind_nbr_ig) + 1 << "\t" << dataset.getIndividualAtPositionFromGroup(i, j).getId()).endLine();
    }
  }
}

void DarwinDon::write(const string& path, const DataSet& dataset, bool overwrite) const
{
  AbstractODataSet::write(path, dataset, overwrite);
}
