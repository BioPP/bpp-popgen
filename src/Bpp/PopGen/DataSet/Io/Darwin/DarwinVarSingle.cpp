//
// File DarwinVarSingle.cpp
// Authors : Sylvain Gaillard
// Last modification : April 7, 2008
//

/*
   Copyright or © or Copr. Bio++ Development Team, (April 7, 2008)

   This software is a computer program whose purpose is to provide classes
   for population genetics analysis.

   This software is governed by the CeCILL  license under French law and
   abiding by the rules of distribution of free software.  You can  use,
   modify and/ or redistribute the software under the terms of the CeCILL
   license as circulated by CEA, CNRS and INRIA at the following URL
   "http://www.cecill.info".

   As a counterpart to the access to the source code and  rights to copy,
   modify and redistribute granted by the license, users are provided only
   with a limited warranty  and the software's author,  the holder of the
   economic rights,  and the successive licensors  have only  limited
   liability.

   In this respect, the user's attention is drawn to the risks associated
   with loading,  using,  modifying and/or developing or reproducing the
   software by the user in light of its specific status of free software,
   that may mean  that it is complicated to manipulate,  and  that  also
   therefore means  that it is reserved for developers  and  experienced
   professionals having in-depth computer knowledge. Users are therefore
   encouraged to load and test the software's suitability as regards their
   requirements in conditions enabling the security of their systems and/or
   data to be ensured and,  more generally, to use and operate it in the
   same conditions as regards security.

   The fact that you are presently reading this means that you have had
   knowledge of the CeCILL license and that you accept its terms.
 */

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
