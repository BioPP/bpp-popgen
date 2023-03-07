//
// File AnalyzedLoci.cpp
// Author : Sylvain Gaillard
// Last modification : Thursday July 29 2004
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

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

