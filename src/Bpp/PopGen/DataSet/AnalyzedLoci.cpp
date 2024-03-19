// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "AnalyzedLoci.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

void AnalyzedLoci::setLocusInfo(
    size_t locusPosition,
    const LocusInfo& locus)
{
  if (locusPosition < loci_.size())
    loci_[locusPosition].reset(locus.clone());
  else
    throw IndexOutOfBoundsException("AnalyzedLoci::setLocusInfo: locus_position out of bounds",
          locusPosition, 0, loci_.size());
}

/******************************************************************************/

size_t AnalyzedLoci::getLocusInfoPosition(
    const std::string& locusName) const
{
  for (size_t i = 0; i < loci_.size(); ++i)
  {
    if (loci_[i] && loci_[i]->getName() == locusName)
      return i;
  }
  throw BadIdentifierException("AnalyzedLoci::getLocusInfoPosition: locus not found.", locusName);
}

/******************************************************************************/

const LocusInfo& AnalyzedLoci::getLocusInfoByName(
    const std::string& locusName) const
{
  for (const auto& locus : loci_)
  {
    if (locus && locus->getName() == locusName)
      return *locus;
  }
  throw BadIdentifierException("AnalyzedLoci::getLocusInfo: locus not found.",
        locusName);
}

/******************************************************************************/

const LocusInfo& AnalyzedLoci::getLocusInfoAtPosition(
    size_t locusPosition) const
{
  if (locusPosition >= loci_.size())
    throw IndexOutOfBoundsException("AnalyzedLoci::getLocusInfoAtPosition: locus_position out of bounds.", locusPosition, 0, loci_.size());
  if (loci_[locusPosition])
    return *loci_[locusPosition];
  else
    throw NullPointerException("AnalyzedLoci::getLocusInfo: no locus defined here.");
}

/******************************************************************************/

// AlleleInfo
void AnalyzedLoci::addAlleleInfoByLocusName(const std::string& locusName,
    const AlleleInfo& allele)
{
  bool locusFound = false;
  for (auto& locus : loci_)
  {
    if (locus->getName() == locusName)
    {
      locusFound = true;
      try
      {
        locus->addAlleleInfo(allele);
      }
      catch (BadIdentifierException& bie)
      {
        throw BadIdentifierException("AnalyzedLoci::addAlleleInfoByLocusName: allele id already in use.", bie.getIdentifier());
      }
    }
  }
  if (!locusFound)
    throw LocusNotFoundException("AnalyzedLoci::addAlleleInfoByLocusName: locus_name not found.",
          locusName);
}

/******************************************************************************/

void AnalyzedLoci::addAlleleInfoByLocusPosition(size_t locusPosition,
    const AlleleInfo& allele)
{
  if (locusPosition < loci_.size())
  {
    try
    {
      loci_[locusPosition]->addAlleleInfo(allele);
    }
    catch (BadIdentifierException& bie)
    {
      throw BadIdentifierException("AnalyzedLoci::addAlleleInfoByLocusPosition: allele id is already in use.", bie.getIdentifier());
    }
  }
  else
    throw IndexOutOfBoundsException("AnalyzedLoci::addAlleleInfoByLocusPosition: locus_position out of bounds.",
          locusPosition, 0, loci_.size());
}

/******************************************************************************/

std::vector<size_t> AnalyzedLoci::getNumberOfAlleles() const
{
  vector<size_t> alleleCount;
  for (const auto& locus: loci_)
  {
    alleleCount.push_back(locus->getNumberOfAlleles());
  }
  return alleleCount;
}

/******************************************************************************/

unsigned int AnalyzedLoci::getPloidyByLocusName(const std::string& locusName) const
{
  for (const auto& locus : loci_)
  {
    if (locus && locus->getName() == locusName)
      return locus->getPloidy();
  }
  throw LocusNotFoundException("AnalyzedLoci::getLocusInfo: locus_name not found.",
        locusName);
}

/******************************************************************************/

unsigned int AnalyzedLoci::getPloidyByLocusPosition(size_t locusPosition) const
{
  if (locusPosition >= loci_.size())
    throw IndexOutOfBoundsException("AnalyzedLoci::getPloidyByLocusPosition: locus_position out of bounds.", locusPosition, 0, loci_.size());
  return loci_[locusPosition]->getPloidy();
}

/******************************************************************************/
