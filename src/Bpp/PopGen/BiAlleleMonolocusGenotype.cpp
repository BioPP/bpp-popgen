//
// File BiAlleleMonolocusGenotype.cpp
// Author : Sylvain Gaillard
// Last modification : Thursday July 29 2004
//

/*
   Copyright or © or Copr. CNRS, (November 17, 2004)

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

#include "BiAlleleMonolocusGenotype.h"

using namespace bpp;
using namespace std;

// ** Class constructor: *******************************************************/

BiAlleleMonolocusGenotype::BiAlleleMonolocusGenotype(
  size_t first_allele_index,
  size_t second_allele_index) : allele_index_(vector<size_t>(2))
{
  allele_index_[0] = first_allele_index;
  allele_index_[1] = second_allele_index;
}

BiAlleleMonolocusGenotype::BiAlleleMonolocusGenotype(std::vector<size_t> allele_index) throw (BadIntegerException) : allele_index_(vector<size_t>(2))
{
  if (allele_index.size() != 2)
    throw BadIntegerException("BiAlleleMonolocusGenotype::BiAlleleMonolocusGenotype: allele_index must contain two values.", static_cast<int>(allele_index.size()));
  allele_index_[0] = allele_index[0];
  allele_index_[1] = allele_index[1];
}

BiAlleleMonolocusGenotype::BiAlleleMonolocusGenotype(const BiAlleleMonolocusGenotype& bmg) : allele_index_(vector<size_t>(2))
{
  for (size_t i = 0; i < 2; i++)
  {
    allele_index_[i] = bmg.getAlleleIndex()[i];
  }
}

// ** Class destructor: ********************************************************/

BiAlleleMonolocusGenotype::~BiAlleleMonolocusGenotype()
{
  allele_index_.clear();
}

// ** Other methodes: **********************************************************/

BiAlleleMonolocusGenotype& BiAlleleMonolocusGenotype::operator=(const BiAlleleMonolocusGenotype& bmg)
{
  for (size_t i = 0; i < 2; i++)
  {
    allele_index_.push_back(bmg.getAlleleIndex()[i]);
  }
  return *this;
}

bool BiAlleleMonolocusGenotype::operator==(const BiAlleleMonolocusGenotype& bmg) const
{
  return (allele_index_[0] == bmg.getAlleleIndex()[0] && allele_index_[1] == bmg.getAlleleIndex()[1])
         || (allele_index_[0] == bmg.getAlleleIndex()[1] && allele_index_[1] == bmg.getAlleleIndex()[0]);
}

size_t BiAlleleMonolocusGenotype::getFirstAlleleIndex() const
{
  return allele_index_[0];
}

size_t BiAlleleMonolocusGenotype::getSecondAlleleIndex() const
{
  return allele_index_[1];
}

bool BiAlleleMonolocusGenotype::isHomozygous() const
{
  return allele_index_[0] == allele_index_[1];
}

std::vector<size_t> BiAlleleMonolocusGenotype::getAlleleIndex() const
{
  return allele_index_;
}

BiAlleleMonolocusGenotype* BiAlleleMonolocusGenotype::clone() const
{
  return new BiAlleleMonolocusGenotype(*this);
}

