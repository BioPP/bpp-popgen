// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Text/TextTools.h>

#include "LocusInfo.h"
#include "GeneralExceptions.h"

using namespace bpp;
using namespace std;

unsigned int LocusInfo::HAPLODIPLOID = 0;
unsigned int LocusInfo::HAPLOID = 1;
unsigned int LocusInfo::DIPLOID = 2;
unsigned int LocusInfo::UNKNOWN = 9999;

// ** Other methods: *********************************************************/

// AlleleInfos
void LocusInfo::addAlleleInfo(const AlleleInfo& allele)
{
  // Check if the allele id is not already in use
  for (const auto& existingAllele : alleles_)
  {
    if (existingAllele->getId() == allele.getId())
      throw BadIdentifierException("LocusInfo::addAlleleInfo: Id already in use.", allele.getId());
  }
  alleles_.push_back(unique_ptr<AlleleInfo>(allele.clone()));
}

const AlleleInfo& LocusInfo::getAlleleInfoById(const std::string& id) const
{
  for (unsigned int i = 0; i < alleles_.size(); i++)
  {
    if (alleles_[i]->getId() == id)
      return *alleles_[i];
  }
  throw AlleleNotFoundException("LocusInfo::getAlleleInfoById: AlleleInfo id unknown.", id);
}

const AlleleInfo& LocusInfo::getAlleleInfoByKey(size_t key) const
{
  if (key >= alleles_.size())
    throw IndexOutOfBoundsException("LocusInfo::getAlleleInfoByKey: key out of bounds.", key, 0, alleles_.size());
  return *(alleles_[key]);
}

unsigned int LocusInfo::getAlleleInfoKey(const std::string& id) const
{
  for (unsigned int i = 0; i < alleles_.size(); ++i)
  {
    if (alleles_[i]->getId() == id)
      return i;
  }
  throw AlleleNotFoundException("LocusInfo::getAlleleInfoKey: AlleleInfo id not found.", id);
}

