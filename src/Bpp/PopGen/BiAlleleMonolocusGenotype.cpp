// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "BiAlleleMonolocusGenotype.h"

using namespace bpp;
using namespace std;

// ** Class constructor: *******************************************************/

BiAlleleMonolocusGenotype::BiAlleleMonolocusGenotype(
    size_t firstAlleleIndex,
    size_t secondAlleleIndex) : alleleIndex_(vector<size_t>(2))
{
  alleleIndex_[0] = firstAlleleIndex;
  alleleIndex_[1] = secondAlleleIndex;
}

BiAlleleMonolocusGenotype::BiAlleleMonolocusGenotype(vector<size_t> alleleIndex) :
  alleleIndex_(2)
{
  if (alleleIndex.size() != 2)
    throw BadSizeException("BiAlleleMonolocusGenotype::BiAlleleMonolocusGenotype: allele_index must contain two values.", alleleIndex.size(), 2);
  alleleIndex_[0] = alleleIndex[0];
  alleleIndex_[1] = alleleIndex[1];
}

BiAlleleMonolocusGenotype::BiAlleleMonolocusGenotype(const BiAlleleMonolocusGenotype& bmg) :
  alleleIndex_(2)
{
  for (size_t i = 0; i < 2; ++i)
  {
    alleleIndex_[i] = bmg.getAlleleIndex()[i];
  }
}

// ** Class destructor: ********************************************************/

BiAlleleMonolocusGenotype::~BiAlleleMonolocusGenotype()
{
  alleleIndex_.clear();
}

// ** Other methodes: **********************************************************/

BiAlleleMonolocusGenotype& BiAlleleMonolocusGenotype::operator=(const BiAlleleMonolocusGenotype& bmg)
{
  for (size_t i = 0; i < 2; ++i)
  {
    alleleIndex_.push_back(bmg.getAlleleIndex()[i]);
  }
  return *this;
}

bool BiAlleleMonolocusGenotype::operator==(const BiAlleleMonolocusGenotype& bmg) const
{
  return (alleleIndex_[0] == bmg.getAlleleIndex()[0] && alleleIndex_[1] == bmg.getAlleleIndex()[1])
         || (alleleIndex_[0] == bmg.getAlleleIndex()[1] && alleleIndex_[1] == bmg.getAlleleIndex()[0]);
}
