// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "MultiAlleleMonolocusGenotype.h"

using namespace bpp;
using namespace std;

// ** Class constructor: *******************************************************/

MultiAlleleMonolocusGenotype::MultiAlleleMonolocusGenotype(std::vector<size_t> allele_index) : alleleIndex_(vector<size_t>(allele_index.size()))
{
  for (size_t i = 0; i < allele_index.size(); ++i)
  {
    alleleIndex_[i] = allele_index[i];
  }
}

MultiAlleleMonolocusGenotype::MultiAlleleMonolocusGenotype(const MultiAlleleMonolocusGenotype& mmg) : alleleIndex_(vector<size_t>(mmg.alleleIndex_.size()))
{
  for (size_t i = 0; i < mmg.getAlleleIndex().size(); ++i)
  {
    alleleIndex_[i] = mmg.getAlleleIndex()[i];
  }
}

// ** Class destructor: ********************************************************/

MultiAlleleMonolocusGenotype::~MultiAlleleMonolocusGenotype()
{
  alleleIndex_.clear();
}

// ** Other methodes: **********************************************************/

MultiAlleleMonolocusGenotype& MultiAlleleMonolocusGenotype::operator=(const MultiAlleleMonolocusGenotype& mmg)
{
  for (size_t i = 0; i < mmg.getAlleleIndex().size(); ++i)
  {
    alleleIndex_.push_back(mmg.getAlleleIndex()[i]);
  }
  return *this;
}

bool MultiAlleleMonolocusGenotype::operator==(const MultiAlleleMonolocusGenotype& mmg) const
{
  return (alleleIndex_[0] == mmg.getAlleleIndex()[0] && alleleIndex_[1] == mmg.getAlleleIndex()[1])
         || (alleleIndex_[0] == mmg.getAlleleIndex()[1] && alleleIndex_[1] == mmg.getAlleleIndex()[0]);
}

bool MultiAlleleMonolocusGenotype::isHomozygous() const
{
  for (size_t i = 1; i < alleleIndex_.size(); ++i)
  {
    if (alleleIndex_[i - 1] != alleleIndex_[i])
      return false;
  }
  return true;
}

