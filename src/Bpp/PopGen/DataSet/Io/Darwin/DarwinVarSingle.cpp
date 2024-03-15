// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "DarwinVarSingle.h"

using namespace bpp;
using namespace std;

DarwinVarSingle::DarwinVarSingle(size_t missingData) : missingData_(missingData) {}

DarwinVarSingle::~DarwinVarSingle() {}

void DarwinVarSingle::write(ostream& os, const DataSet& dataset) const
{
  if (!os)
    throw IOException("DarwinVarSingle::write: fail to open stream.");
  StlOutputStreamWrapper out(&os);
  (out << "@DARwin 5.0 - SINGLE").endLine();
  size_t ind_nbr = 0;
  for (size_t i = 0; i < dataset.getNumberOfGroups(); i++)
  {
    ind_nbr += dataset.getNumberOfIndividualsInGroup(i);
  }
  vector<string> header;
  header.push_back("Unit");
  for (size_t i = 0; i < dataset.getNumberOfLoci(); ++i)
  {
    const LocusInfo& li = dataset.getLocusInfoAtPosition(i);
    for (size_t j = 0; j < li.getNumberOfAlleles(); ++j)
    {
      header.push_back(li.getName() + "." + li.getAlleleInfoByKey(j).getId());
    }
  }
  size_t var_nbr = header.size() - 1;
  // header.push_back("Name");
  (out << ind_nbr << "\t" << var_nbr).endLine();
  VectorTools::print(header, out, "\t");
  // size_t ind_index = 0;
  const auto& al = dataset.analyzedLoci();
  for (size_t i = 0; i < dataset.getNumberOfGroups(); ++i)
  {
    size_t ind_nbr_ig = dataset.getNumberOfIndividualsInGroup(i);
    for (size_t j = 0; j < ind_nbr_ig; ++j)
    {
      vector<size_t> var;
      const auto& geno = dataset.getIndividualAtPositionFromGroup(i, j).getGenotype();
      for (size_t k = 0; k < geno.size(); ++k)
      {
        const auto& mg = geno.monolocusGenotype(k);
        if (geno.isMonolocusGenotypeMissing(k))
        {
          for (size_t l = 0; l < al.getNumberOfAlleles()[k]; ++l)
          {
            var.push_back(missingData_);
          }
        }
        else
        {
          for (size_t l = 0; l < al.getNumberOfAlleles()[k]; ++l)
          {
            size_t flag = 0;
            if (VectorTools::contains(mg.getAlleleIndex(), l))
              flag = 1;
            var.push_back(flag);
          }
        }
        // var.push_back((mg->getAlleleIndex()).size());
      }
      (out << j + (i * ind_nbr_ig) + 1 << "\t" << VectorTools::paste(var, "\t")).endLine();
    }
  }
}

void DarwinVarSingle::write(const string& path, const DataSet& dataset, bool overwrite) const
{
  AbstractODataSet::write(path, dataset, overwrite);
}
