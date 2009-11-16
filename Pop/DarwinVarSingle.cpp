//
// File DarwinVarSingle.cpp
// Authors : Sylvain Gaillard
// Last modification : April 7, 2008
//

/*
   Copyright or © or Copr. CNRS, (April 7, 2008)

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

DarwinVarSingle::DarwinVarSingle() {}

DarwinVarSingle::~DarwinVarSingle() {}

const string DarwinVarSingle::getFormatName()
{
  return "Darwin .var single data";
}

const string DarwinVarSingle::getFormatDescription()
{
  return "Darwin .var file store data for each marker in each individual (1 variable per allele).";
}

void DarwinVarSingle::write(ostream & os, const DataSet & data_set) const throw (Exception)
{
  if (!os)
    throw IOException("DarwinVarSingle::write: fail to open stream.");
  os << "@DARwin 5.0 - SINGLE" << endl;
  unsigned int ind_nbr = 0;
  for (unsigned int i = 0 ; i < data_set.getNumberOfGroups() ; i++)
    ind_nbr += data_set.getNumberOfIndividualsInGroup(i);
  vector<string> header;
  header.push_back("Unit");
  for (unsigned int i = 0 ; i < data_set.getNumberOfLoci() ; i++) {
    const LocusInfo * li = data_set.getLocusInfoAtPosition(i);
    for (unsigned int j = 0 ; j < li->getNumberOfAlleles() ; j++) {
      header.push_back(li->getName() + "." + li->getAlleleInfoByKey(j)->getId());
    }
  }
  unsigned int var_nbr = header.size() - 1;
  //header.push_back("Name");
  os << ind_nbr << "\t" << var_nbr << endl;
  VectorTools::print(header, os, "\t");
  //unsigned int ind_index = 0;
  for (unsigned int i = 0 ; i < data_set.getNumberOfGroups() ; i++) {
    unsigned int ind_nbr_ig = data_set.getNumberOfIndividualsInGroup(i);
    for (unsigned int j = 0 ; j < ind_nbr_ig ; j++) {
      vector<unsigned int> var;
      const MultilocusGenotype * geno = data_set.getIndividualAtPositionFromGroup(i, j)->getGenotype();
      for (unsigned int k = 0 ; k < geno->size() ; k++) {
        const MonolocusGenotype * mg = geno->getMonolocusGenotype(k);
        const AnalyzedLoci * al = data_set.getAnalyzedLoci();
        if (geno->isMonolocusGenotypeMissing(k)) {
          for (unsigned int l = 0 ; l < al->getNumberOfAlleles()[k] ; l++)
            var.push_back(0);
        }
        else {
          for (unsigned int l = 0 ; l < al->getNumberOfAlleles()[k] ; l++) {
            unsigned int flag = 0;
            if (VectorTools::contains(mg->getAlleleIndex(), l))
              flag = 1;
            var.push_back(flag);
          }
        }
        //var.push_back((mg->getAlleleIndex()).size());
      }
      os << j + (i * ind_nbr_ig) + 1 << "\t" << VectorTools::paste(var, "\t") << endl;
    }
  }
}

void DarwinVarSingle::write(const string & path, const DataSet & data_set, bool overwrite) const throw (Exception)
{
  AbstractODataSet::write(path, data_set, overwrite);
}

