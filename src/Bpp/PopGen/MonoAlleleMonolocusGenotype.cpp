// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "MonoAlleleMonolocusGenotype.h"

using namespace bpp;
using namespace std;

// ** Class constructor: *******************************************************/

MonoAlleleMonolocusGenotype::MonoAlleleMonolocusGenotype(size_t allele_index) :
  alleleIndex_(allele_index) {}

MonoAlleleMonolocusGenotype::MonoAlleleMonolocusGenotype(std::vector<size_t> allele_index) : alleleIndex_(0)
{
  if (allele_index.size() != 1)
    throw BadSizeException("MonoAlleleMonolocusGenotype::MonoAlleleMonolocusGenotype: allele_index must contain one value.", allele_index.size(), 1);
  alleleIndex_ = allele_index[0];
}

MonoAlleleMonolocusGenotype::MonoAlleleMonolocusGenotype(const MonoAlleleMonolocusGenotype& mmg) :
  alleleIndex_(mmg.getAlleleIndex()[0]) {}

// ** Other methods: **********************************************************/

MonoAlleleMonolocusGenotype& MonoAlleleMonolocusGenotype::operator=(const MonoAlleleMonolocusGenotype& mmg)
{
  alleleIndex_ = mmg.getAlleleIndex()[0];
  return *this;
}

bool MonoAlleleMonolocusGenotype::operator==(const MonoAlleleMonolocusGenotype& mmg) const
{
  return alleleIndex_ == mmg.getAlleleIndex()[0];
}

std::vector<size_t> MonoAlleleMonolocusGenotype::getAlleleIndex() const
{
  vector<size_t> index;
  index.push_back(alleleIndex_);
  return index;
}
